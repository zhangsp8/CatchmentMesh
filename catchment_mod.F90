MODULE catchment_mod

   USE hydro_data_mod
   USE river_lake_mod, only : binfo

   IMPLICIT NONE

   type catinfo_typ
      integer (kind=4), allocatable :: cid    (:)  ! catchment id
      integer (kind=4), allocatable :: lid    (:)  ! lake id
      integer (kind=4), allocatable :: nswe (:,:)  ! boundaries: north, south, west, east
      type(catinfo_typ), pointer :: next
   END type 

CONTAINS

   SUBROUTINE get_catchment (catsize)

      USE task_mod
      USE cvt_mod
      IMPLICIT NONE

      real (kind=4), intent(in) :: catsize

      integer (kind=4) :: icatdsp, iproc, cdsp, iwork, mesg(2), isrc, nsend, nrecv, ndone
      integer (kind=4) :: catnum, ncat, ncat_work, nexpd
      integer (kind=4) :: icat, i, j, i_up, j_up, idir, iblk, jblk, iblk0, jblk0
      integer (kind=1) :: dir_this
      
      integer (kind=4) :: nrivp, irivp, irivseg, tail, head
      integer (kind=4), allocatable :: river (:,:)

      integer (kind=4) :: nplist, iplist, jplist, npall
      integer (kind=4), allocatable :: pixellist (:,:)

      real (kind=4) :: upa_up(8), elv0, elv
      real (kind=4), allocatable :: elvdata(:)

      type(catinfo_typ), target  :: networkinfo
      type(catinfo_typ), pointer :: thislake, nextlake
      integer (kind=4), allocatable :: cat_id   (:)
      integer (kind=4), allocatable :: lake_id  (:)
      integer (kind=4), allocatable :: bndinfo(:,:)
      
      integer (kind=4) :: south, north, west, east, cdthis
      integer (kind=4), allocatable :: catdown(:)
      
      ! lake
      logical :: is_lake
      integer (kind=4) :: nshorecat, nc, ic, np, ip, jp
      integer (kind=4) :: outlet(2), lakeid, nplake, iplake, icatlast, i_dn, j_dn, npshore
      real    (kind=4) :: shorearea
      integer (kind=4), allocatable :: lakepixels(:,:), lakeindex (:)
      integer (kind=4), allocatable :: shoreline (:,:), shoreindex(:), thiscat(:,:) 
      real    (kind=4), allocatable :: upashore(:), hndshore(:)


      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) THEN

         write(*,'(/3A/)') 'Step 2 ', trim(binfo), ': Finding all catchments and lakes ...'

         ncat = thisinfo%nrivseg

         CALL mpi_bcast (thisinfo%icatdsp, 1, MPI_INTEGER, 0, p_comm_work, p_err)

         cdsp = 0
         DO iproc = 1, p_nwork
            ncat_work = ncat / p_nwork
            IF (iproc <= mod(ncat,p_nwork)) THEN
               ncat_work = ncat_work + 1
            ENDIF

            CALL mpi_send (ncat_work, 1, MPI_INTEGER4, iproc, iproc, p_comm_work, p_err) 

            head = cdsp
            IF (ncat_work > 0) THEN
               irivseg = 0
               DO WHILE (head < thisinfo%nrivpix)
                  head = head + 1
                  DO WHILE (head < thisinfo%nrivpix)
                     IF (thisinfo%riv_pix(3,head+1) == thisinfo%riv_pix(3,head)) THEN
                        head = head + 1
                     ELSE
                        EXIT
                     ENDIF
                  ENDDO

                  irivseg = irivseg + 1
                  IF (irivseg == ncat_work) EXIT
               ENDDO

               nsend = head - cdsp
               CALL mpi_send (nsend, 1, MPI_INTEGER4, iproc, iproc, p_comm_work, p_err) 
               CALL mpi_send (thisinfo%riv_pix(:,cdsp+1:head), nsend*3, MPI_INTEGER4, &
                  iproc, iproc, p_comm_work, p_err) 
            ENDIF

            cdsp = head
         ENDDO

         ! for lake
         ndone = 0
         DO WHILE (.true.)
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, 0, p_comm_work, p_stat, p_err)
            iwork = mesg(1)
            nexpd = mesg(2)
            IF (nexpd > 0) THEN
               CALL mpi_send (ncat, 1, MPI_INTEGER, iwork, 0, p_comm_work, p_err) 
               ncat = ncat + nexpd
            ELSE
               ndone = ndone + 1
               IF (ndone == p_nwork) EXIT
            ENDIF
         ENDDO

         CALL mpi_barrier (p_comm_work, p_err)

         thisinfo%ntotalcat = ncat
         allocate (thisinfo%lake_id    (ncat))
         allocate (thisinfo%bsn_bnds (4,ncat))

         DO iproc = 1, p_nwork
            CALL mpi_recv (nrecv, 1, MPI_INTEGER, iproc, 0, p_comm_work, p_stat, p_err)
            IF (nrecv > 0) THEN
               allocate (cat_id   (nrecv))
               allocate (lake_id  (nrecv))
               allocate (bndinfo(4,nrecv))

               CALL mpi_recv (cat_id,  nrecv,   MPI_INTEGER, iproc, 1, p_comm_work, p_stat, p_err)
               CALL mpi_recv (lake_id, nrecv,   MPI_INTEGER, iproc, 1, p_comm_work, p_stat, p_err)
               CALL mpi_recv (bndinfo, 4*nrecv, MPI_INTEGER, iproc, 1, p_comm_work, p_stat, p_err)

               cat_id = cat_id - thisinfo%icatdsp

               DO icat = 1, nrecv
                  thisinfo%lake_id (cat_id(icat))   = lake_id(icat)
                  thisinfo%bsn_bnds(:,cat_id(icat)) = bndinfo(:,icat)
               ENDDO

               deallocate(cat_id)
               deallocate(lake_id)
               deallocate(bndinfo)
            ENDIF
         ENDDO

         CALL excute_data_task (t_exit)

      ENDIF

      IF (p_is_data) THEN

         CALL data_daemon ()

      ENDIF

      IF (p_is_work) THEN

         CALL sync_window ()

         CALL mpi_bcast (icatdsp, 1, MPI_INTEGER, 0, p_comm_work, p_err)

         CALL mpi_recv (ncat, 1, MPI_INTEGER, 0, p_iam_work, p_comm_work, p_stat, p_err)

         IF (ncat > 0) THEN

            CALL mpi_recv (nrivp, 1, MPI_INTEGER, 0, p_iam_work, p_comm_work, p_stat, p_err)
            allocate (river (3,nrivp))
            CALL mpi_recv (river, nrivp*3, MPI_INTEGER, 0, p_iam_work, p_comm_work, p_stat, p_err)

            allocate (pixellist (2,100000))
            allocate (lakepixels(2,100000))
            allocate (elvdata     (100000))

            allocate (networkinfo%cid    (ncat))
            allocate (networkinfo%lid    (ncat))
            allocate (networkinfo%nswe (4,ncat))

            networkinfo%nswe(:,:) = 0

            thislake => networkinfo
            thislake%next => null()

            irivseg = 0
            head = 0
            DO WHILE (head < nrivp)

               irivseg = irivseg + 1
               head = head + 1
               
               is_lake = (get_lake(river(1,head),river(2,head)) > 0)

               IF (.not. is_lake) THEN
               
                  tail = head
                  DO WHILE (head < nrivp)
                     IF (river(3,head+1) == river(3,head)) THEN
                        head = head + 1
                     ELSE
                        EXIT
                     ENDIF
                  ENDDO

                  npall = 0
                  DO irivp = head, tail, -1

                     i = river(1,irivp)
                     j = river(2,irivp)

                     elv0 = get_elv(i,j)

                     nplist = 0
                     CALL append_plist3 (pixellist, nplist, i, j, elvdata, 0._4)

                     iplist = 1
                     DO WHILE (iplist <= nplist)
                        DO idir = 1, 8
                           CALL nextij (pixellist(1,iplist), pixellist(2,iplist), &
                              ishftc(int(-128,1),idir), i_up, j_up)

                           IF (((pixellist(1,iplist) /= i_up) .or. (pixellist(2,iplist) /= j_up)) &
                              .and. within_region(i_up,j_up)) THEN
                              IF (get_dir(i_up,j_up) == ishftc(int(8,1),idir)) THEN
                                 IF (get_icat(i_up,j_up) == 0) THEN
                                    elv = get_elv(i_up, j_up)
                                    CALL append_plist3 (pixellist, nplist, i_up, j_up, elvdata, elv-elv0)
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDDO

                        iplist = iplist + 1
                     ENDDO                     

                     npall = npall + nplist

                     catnum = river(3,head)
                     CALL update_catchment_pixels     (catnum, nplist, pixellist, elvdata)
                     CALL update_catchment_boundaries (nplist, pixellist, irivseg, networkinfo)
                     
                     networkinfo%cid(irivseg) = catnum
                     networkinfo%lid(irivseg) = 0

                  ENDDO
                     
                  write(*,100) trim(binfo), p_iam_glb, networkinfo%cid(irivseg), networkinfo%nswe(:,irivseg), &
                     irivseg, ncat, npall, head-tail+1
                  100 format('(S2) Catchment ',A,': On', I4, ' ID ', I7, ' (',I6,',',I6,',',I6,',',I6,'),', &
                     ' (', I5, '/', I5, ' done), ', I10, ' Pixels, ', I3, ' Rivers') 

               ELSE
                  ! lake
                  
                  ! 1. lake area.
                  outlet = river(1:2,head)

                  lakeid = get_lake(outlet(1), outlet(2))

                  nplake = 0
                  CALL append_plist (lakepixels, nplake, outlet(1), outlet(2), check_exist = .false.)

                  iplake = 1
                  DO WHILE (iplake <= nplake)
                     DO idir = 1, 8
                        CALL nextij (lakepixels(1,iplake), lakepixels(2,iplake), &
                           ishftc(int(-128,1),idir), i_up, j_up)

                        IF (((lakepixels(1,iplake) /= i_up) .or. (lakepixels(2,iplake) /= j_up)) &
                           .and. within_region(i_up,j_up)) THEN
                           IF (get_dir(i_up,j_up) == ishftc(int(8,1),idir)) THEN
                              IF (get_lake(i_up,j_up) == lakeid) THEN
                                 CALL append_plist3 (lakepixels, nplake, i_up, j_up, elvdata, 0._4)
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDDO

                     iplake = iplake + 1
                  ENDDO                     
                     
                  catnum = river(3,head)
                  CALL update_catchment_pixels     (catnum, nplake, lakepixels, elvdata)
                  CALL update_catchment_boundaries (nplake, lakepixels, irivseg, networkinfo)
                     
                  networkinfo%cid(irivseg) = catnum
                  networkinfo%lid(irivseg) = lakeid
                  
                  write(*,101) trim(binfo), p_iam_glb, networkinfo%cid(irivseg), networkinfo%nswe(:,irivseg), &
                     irivseg, ncat, nplake
                  101 format('(S2) Lake      ',A,': On', I4, ' ID ', I7, ' (',I6,',',I6,',',I6,',',I6,'),', &
                     ' (', I5, '/', I5, ' done), ', I10, ' Pixels') 

                  ! 2. catchments along shore lines.
                  npshore = 0
                  DO iplake = 1, nplake
                     DO idir = 1, 8
                        CALL nextij (lakepixels(1,iplake), lakepixels(2,iplake), &
                           ishftc(int(-128,1),idir), i_up, j_up)

                        IF (((lakepixels(1,iplake) /= i_up) .or. (lakepixels(2,iplake) /= j_up)) &
                           .and. within_region(i_up,j_up)) THEN
                           IF (get_dir(i_up,j_up) == ishftc(int(8,1),idir)) THEN
                              IF (get_lake(i_up,j_up) /= lakeid) THEN
                                 IF (get_upa(i_up,j_up) < catsize) THEN
                                    CALL append_plist (shoreline, npshore, i_up, j_up, check_exist = .false.)
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
                     
                  nshorecat = 0

                  IF (npshore > 0) THEN

                     allocate(upashore(npshore))
                     allocate(hndshore(npshore))

                     DO ip = 1, npshore
                        upashore(ip) = get_upa(shoreline(1,ip),shoreline(2,ip))

                        i_up = shoreline(1,ip)
                        j_up = shoreline(2,ip)
                        CALL nextij (i_up,j_up, get_dir(i_up,j_up), i,j)
                        hndshore(ip) = get_elv(i_up,j_up) - get_elv(i,j)
                     ENDDO

                     allocate(shoreindex(npshore))
                     shoreindex(:) = -1

                     DO WHILE (any(shoreindex == -1))

                        nshorecat = nshorecat + 1
                        shorearea = 0

                        ip = minloc(hndshore, mask = shoreindex == -1, dim=1)                        
                        shoreindex(ip) = 0

                        DO WHILE (any(shoreindex == 0))
                           
                           ip = minloc(hndshore, mask = shoreindex == 0, dim=1) 
                           shorearea = shorearea + upashore(ip)

                           IF (shorearea >= catsize) THEN
                              WHERE (shoreindex == 0) 
                                 shoreindex = -1
                              END WHERE
                              
                              EXIT
                              
                           ENDIF
                           
                           shoreindex(ip) = nshorecat

                           DO jp = 1, npshore
                              IF (shoreindex(jp) == -1) THEN
                                 IF (abs(shoreline(1,ip)-shoreline(1,jp)) <= 1) THEN
                                    IF (shoreline(2,ip) == mglb) THEN
                                       IF ((shoreline(2,jp) == 1) .or. (shoreline(2,jp) == mglb-1)) THEN
                                          shoreindex(jp) = 0
                                       ENDIF
                                    ELSEIF (shoreline(2,ip) == 1) THEN
                                       IF ((shoreline(2,jp) == mglb) .or. (shoreline(2,jp) == 2)) THEN
                                          shoreindex(jp) = 0
                                       ENDIF
                                    ELSEIF (abs(shoreline(2,ip)-shoreline(2,jp)) <= 1) THEN
                                       shoreindex(jp) = 0
                                    ENDIF
                                 ENDIF
                              ENDIF
                           ENDDO

                        ENDDO

                     ENDDO

                     deallocate (upashore)
                     deallocate (hndshore)

                  ENDIF

                  IF (nshorecat > 0) THEN
                  
                     ! inquire catchment number displacement
                     mesg = (/p_iam_work, nshorecat/)
                     CALL mpi_send (mesg,     2, MPI_INTEGER, 0, 0, p_comm_work, p_err) 
                     CALL mpi_recv (icatlast, 1, MPI_INTEGER, 0, 0, p_comm_work, p_stat, p_err)
                     
                     allocate (thislake%next)
                     thislake => thislake%next
                     thislake%next => null()

                     allocate (thislake%cid    (nshorecat))
                     allocate (thislake%lid    (nshorecat))
                     allocate (thislake%nswe (4,nshorecat))

                     thislake%nswe(:,:) = 0

                     DO ic = 1, nshorecat

                        np = count(shoreindex == ic)
                        allocate(thiscat (2,np))
                        thiscat(1,:) = pack(shoreline(1,1:npshore), shoreindex == ic)
                        thiscat(2,:) = pack(shoreline(2,1:npshore), shoreindex == ic)

                        nplist = 0
                        iplist = 1
                        DO ip = 1, np

                           i_up = thiscat(1,ip)
                           j_up = thiscat(2,ip)

                           CALL nextij (i_up,j_up, get_dir(i_up,j_up), i,j)
                           elv0 = get_elv(i,j)
                                          
                           elv = get_elv(i_up, j_up)
                           CALL append_plist3 (pixellist, nplist, i_up, j_up, elvdata, elv-elv0)

                           DO WHILE (iplist <= nplist)

                              i = pixellist(1,iplist); j = pixellist(2,iplist)

                              DO idir = 1, 8
                                 CALL nextij (i, j, ishftc(int(-128,1),idir), i_up, j_up)
                                 IF (((i /= i_up) .or. (j /= j_up)) .and. within_region(i_up,j_up)) THEN
                                    IF (get_dir(i_up,j_up) == ishftc(int(8,1),idir)) THEN
                                       IF (get_icat(i_up,j_up) == 0) THEN
                                          elv = get_elv(i_up, j_up)
                                          CALL append_plist3 (pixellist, nplist, i_up, j_up, &
                                             elvdata, elv-elv0)
                                       ENDIF
                                    ENDIF
                                 ENDIF
                              ENDDO

                              iplist = iplist + 1

                           ENDDO                     
                        ENDDO

                        catnum = icatdsp + icatlast + ic

                        CALL update_catchment_pixels     (catnum, nplist, pixellist, elvdata)
                        CALL update_catchment_boundaries (nplist, pixellist, ic, thislake)

                        thislake%cid(ic) = catnum
                        thislake%lid(ic) = -1

                        write(*,102) trim(binfo), p_iam_glb, thislake%cid(ic), thislake%nswe(:,ic), ic, nshorecat, nplist
                        102 format('(S2) LakeCat   ',A,': On', I4, ' ID ', I7, ' (',I6,',',I6,',',I6,',', &
                           I6,'),', ' (', I5, '/', I5, ' done), ', I10, ' Pixels') 

                        deallocate (thiscat)

                     ENDDO
                           
                     deallocate (shoreline )
                     deallocate (shoreindex)

                  ENDIF

               ENDIF

            ENDDO

            deallocate (river)
            deallocate (pixellist)
            deallocate (lakepixels)
            deallocate (elvdata)

         ENDIF

         mesg = (/p_iam_work, 0/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, 0, 0, p_comm_work, p_err) 
         
         CALL mpi_barrier (p_comm_work, p_err)

         IF (ncat > 0) THEN
            thislake => networkinfo%next
            DO WHILE (associated(thislake))
               ncat = ncat + size(thislake%cid)
               thislake => thislake%next
            ENDDO
         ENDIF
         
         CALL mpi_send (ncat, 1, MPI_INTEGER4, 0, 0, p_comm_work, p_err) 

         IF (ncat > 0) THEN
           
            allocate (cat_id    (ncat))
            allocate (lake_id   (ncat))
            allocate (bndinfo (4,ncat))

            nc = size(networkinfo%cid)
           
            cat_id   (1:nc) = networkinfo%cid
            lake_id  (1:nc) = networkinfo%lid
            bndinfo(:,1:nc) = networkinfo%nswe

            ic = nc
            thislake => networkinfo%next
            DO WHILE (associated(thislake))
               nc = size(thislake%cid)
               
               cat_id   (ic+1:ic+nc) = thislake%cid
               lake_id  (ic+1:ic+nc) = thislake%lid
               bndinfo(:,ic+1:ic+nc) = thislake%nswe

               ic = ic + nc
               thislake => thislake%next
            ENDDO

            CALL mpi_send (cat_id,  ncat,   MPI_INTEGER4, 0, 1, p_comm_work, p_err) 
            CALL mpi_send (lake_id, ncat,   MPI_INTEGER4, 0, 1, p_comm_work, p_err) 
            CALL mpi_send (bndinfo, 4*ncat, MPI_INTEGER4, 0, 1, p_comm_work, p_err) 

         
            deallocate (networkinfo%cid )
            deallocate (networkinfo%lid )
            deallocate (networkinfo%nswe)

            thislake => networkinfo%next
            DO WHILE (associated(thislake))
               nextlake => thislake%next

               deallocate(thislake%cid )
               deallocate(thislake%lid )
               deallocate(thislake%nswe)
               deallocate(thislake)

               thislake => nextlake
            ENDDO

         ENDIF

         
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) THEN

         allocate (thisinfo%bsn_downstream(thisinfo%ntotalcat))

         DO iproc = 1, p_nwork
            CALL mpi_recv (nrecv, 1, MPI_INTEGER, iproc, 0, p_comm_work, p_stat, p_err)
            IF (nrecv > 0) THEN
               allocate(cat_id (ncat))
               allocate(catdown(ncat))

               CALL mpi_recv (cat_id,  nrecv, MPI_INTEGER, iproc, 1, p_comm_work, p_stat, p_err)
               CALL mpi_recv (catdown, nrecv, MPI_INTEGER, iproc, 1, p_comm_work, p_stat, p_err)

               cat_id = cat_id - thisinfo%icatdsp

               DO icat = 1, nrecv
                  thisinfo%bsn_downstream(cat_id(icat)) = catdown(icat)
               ENDDO

               deallocate(cat_id)
               deallocate(catdown)
            ENDIF
         ENDDO
         
         CALL excute_data_task (t_exit)

      ENDIF

      IF (p_is_data) THEN
         CALL data_daemon ()
      ENDIF

      IF (p_is_work) THEN 

         CALL sync_window ()

         IF (ncat > 0) THEN

            allocate (catdown (ncat))
            catdown(:) = -9999

            DO icat = 1, ncat
      
               catnum = cat_id   (icat)
               north  = bndinfo(1,icat)
               south  = bndinfo(2,icat)
               west   = bndinfo(3,icat)
               east   = bndinfo(4,icat)
                        
               ip = 0
               DO i = north, south
                  j = west
                  DO WHILE (.true.)
                     IF (get_icat(i,j) == catnum) THEN
                        dir_this = get_dir (i, j)
                        ! 0: river mouth; -1: inland depression; 
                        IF ((dir_this /= 0) .and. (dir_this /= -1)) THEN
                           CALL nextij (i, j, dir_this, i_dn, j_dn)
                           IF (within_region(i_dn,j_dn)) THEN
                              cdthis = get_icat(i_dn,j_dn)
                              IF (cdthis == 0) cdthis = -3 ! downstream is land, but not in the region.
                           ELSE
                              cdthis = -3
                           ENDIF

                           IF (cdthis /= catnum) THEN
                              IF (catdown(icat) /= -9999) THEN
                                 IF (cdthis /= catdown(icat)) THEN
                                    write(*,'(A,I7,A,I7,A,I7,A,I7,A)') &
                                       'Warning: more than ONE downstream catchment : ',  &
                                       cat_id(icat), ',', lake_id(icat), '(', catdown(icat), '->', &
                                       cdthis, ')'
                                 ENDIF
                              ELSE
                                 catdown(icat) = cdthis
                              ENDIF
                           ENDIF
                        ELSE
                           catdown(icat) = dir_this
                        ENDIF
                     ENDIF

                     IF (j == east) THEN
                        EXIT
                     ELSE
                        j = mod(j,mglb) + 1
                     ENDIF
                  ENDDO
               ENDDO
                           
               write(*,103) trim(binfo), p_iam_glb, cat_id(icat), catdown(icat)
               103 format('(S2) Downstream ', A, ': On', I4, ' From ', I7, ' to ', I7)
            ENDDO
            
         ENDIF

         CALL mpi_send (ncat, 1, MPI_INTEGER4, 0, 0, p_comm_work, p_err) 
         IF (ncat > 0) THEN
            CALL mpi_send (cat_id,  ncat, MPI_INTEGER4, 0, 1, p_comm_work, p_err) 
            CALL mpi_send (catdown, ncat, MPI_INTEGER4, 0, 1, p_comm_work, p_err) 
         ENDIF

         IF (allocated(cat_id )) deallocate(cat_id )
         IF (allocated(bndinfo)) deallocate(bndinfo)
         IF (allocated(lake_id)) deallocate(lake_id)
         IF (allocated(catdown)) deallocate(catdown)

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)

   END SUBROUTINE get_catchment 


   SUBROUTINE update_catchment_pixels (catnum, np, pixellist, elvdata)

      USE task_mod
      IMPLICIT NONE

      integer (kind = 4), intent(in) :: catnum, np
      integer (kind = 4), intent(in) :: pixellist(:,:)
      real    (kind = 4), intent(in) :: elvdata  (:)

      ! Local
      integer :: ip, jp, iblk, jblk, iblk0, jblk0
      integer, allocatable :: rmark (:)

      integer :: scnt, idata
      integer (kind = 4), allocatable :: ibuf (:,:)
      real    (kind = 4), allocatable :: rbuf (:)

      allocate (rmark (np))
      rmark = 0

      DO ip = 1, np
         IF (rmark(ip) == 0) THEN
            rmark(ip) = ip
            iblk = (pixellist(1,ip)-1)/nbox + 1
            jblk = (pixellist(2,ip)-1)/mbox + 1
            DO jp = ip+1, np
               iblk0 = (pixellist(1,jp)-1)/nbox + 1
               jblk0 = (pixellist(2,jp)-1)/mbox + 1
               IF ((iblk0 == iblk) .and. (jblk0 == jblk)) THEN
                  rmark(jp) = ip
               ENDIF
            ENDDO

            scnt = count(rmark(ip:np) == ip)
            allocate (ibuf (2,scnt))
            allocate (rbuf (scnt))

            scnt = 0
            DO jp = ip, np
               IF (rmark(jp) == ip) THEN
                  scnt = scnt + 1
                  ibuf(:,scnt) = pixellist(:,jp)
                  rbuf(scnt) = elvdata(jp)
               ENDIF
            ENDDO

            CALL excute_data_task (t_update_icatch)
            CALL mpi_send (scnt,   1, MPI_INTEGER4, p_data_address, t_update_icatch, p_comm_glb, p_err) 
            CALL mpi_send (iblk,   1, MPI_INTEGER4, p_data_address, t_update_icatch, p_comm_glb, p_err) 
            CALL mpi_send (jblk,   1, MPI_INTEGER4, p_data_address, t_update_icatch, p_comm_glb, p_err) 
            CALL mpi_send (catnum, 1, MPI_INTEGER4, p_data_address, t_update_icatch, p_comm_glb, p_err) 

            CALL mpi_send (ibuf, scnt*2, MPI_INTEGER, p_data_address, t_update_icatch, p_comm_glb, p_err) 
            CALL mpi_send (rbuf, scnt,   MPI_REAL4,   p_data_address, t_update_icatch, p_comm_glb, p_err) 

            deallocate (ibuf)
            deallocate (rbuf)
         ENDIF
      ENDDO

      deallocate (rmark)

   END SUBROUTINE update_catchment_pixels

   SUBROUTINE update_catchment_boundaries (np, pixellist, icatloc, catinfo)

      IMPLICIT NONE

      integer :: np, icatloc
      integer :: pixellist(:,:)
      type(catinfo_typ) :: catinfo

      integer :: ipixel

      IF (any(catinfo%nswe(1:4,icatloc) == 0)) THEN
         catinfo%nswe(1:4,icatloc) = (/pixellist(1,1), pixellist(1,1), pixellist(2,1), pixellist(2,1)/)
      ENDIF

      DO ipixel = 1, np
         catinfo%nswe(1,icatloc) = min     (catinfo%nswe(1,icatloc), pixellist(1,ipixel))
         catinfo%nswe(2,icatloc) = max     (catinfo%nswe(2,icatloc), pixellist(1,ipixel))
         catinfo%nswe(3,icatloc) = min_west(catinfo%nswe(3,icatloc), pixellist(2,ipixel))
         catinfo%nswe(4,icatloc) = max_east(catinfo%nswe(4,icatloc), pixellist(2,ipixel))
      ENDDO

   END SUBROUTINE update_catchment_boundaries


   SUBROUTINE get_basin_neighbour (icat, imin, jmin, np, mp, catch, &
         nnb, nbindex, lenborder)

      USE hydro_data_mod
      USE task_mod

      IMPLICIT NONE

      integer (kind = 4), intent(in) :: icat 
      integer (kind = 4), intent(in) :: imin, jmin
      integer (kind = 4), intent(in) :: np, mp
      integer (kind = 4), intent(in) :: catch (np,mp)    

      integer (kind = 4), intent(out) :: nnb
      integer (kind = 4), allocatable, intent(out) :: nbindex  (:)
      real    (kind = 4), allocatable, intent(out) :: lenborder(:)

      ! Local Variable
      integer :: ithis, jthis, inext, jnext, igthis, jgthis, ignext, jgnext
      real (kind = 4) :: blen

      nnb = 0
      allocate (nbindex  (np*mp))
      allocate (lenborder(np*mp))
      lenborder(:) = 0

      DO ithis = 1, np
         DO jthis = 1, mp
            IF (catch(ithis,jthis) == icat) THEN

               DO inext = ithis-1, ithis+1
                  DO jnext = jthis-1, jthis+1
                     IF ((inext >= 1) .and. (inext <= np) .and. (jnext >= 1) .and. (jnext <= mp) &
                        .and. ((inext /= ithis) .or. (jnext /= jthis))) THEN
                        IF (catch(inext,jnext) /= icat) THEN

                           igthis = ithis + imin - 1
                           ignext = inext + imin - 1
                           jgthis = jthis + jmin - 1;  IF (jgthis > mglb) jgthis = jgthis - mglb
                           jgnext = jnext + jmin - 1;  IF (jgnext > mglb) jgnext = jgnext - mglb

                           blen = dist_between (igthis, jgthis, ignext, jgnext)

                           IF (abs(inext-ithis)+abs(jnext-jthis) == 1) THEN
                              blen = 0.5 * blen
                           ELSE
                              blen = 0.25 * blen
                           ENDIF
                              
                           CALL acc_border(icat, nnb, nbindex, catch(inext,jnext), lenborder, blen)

                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO

            ENDIF
         ENDDO
      ENDDO

      IF (nnb == 0) THEN
         write(*,*) 'Neighbours: ', icat, ' has no neighbours.'
      ELSE
         write(*,*) 'Neighbours: ', icat, ' has ', nnb, 'neighbours', sum(lenborder(1:nnb))/nnb, &
            minval(lenborder(1:nnb)), maxval(lenborder(1:nnb))
      ENDIF
      
   END SUBROUTINE get_basin_neighbour


   SUBROUTINE acc_border (icat, nnb, nbindex, jcat, borderlen, blen)

      IMPLICIT NONE

      integer (kind = 4), intent(in) :: icat 

      integer (kind = 4), intent(inout) :: nnb
      integer (kind = 4), intent(inout) :: nbindex(:)
      integer (kind = 4), intent(in) :: jcat 

      real (kind = 4), intent(inout) :: borderlen(:)
      real (kind = 4), intent(in) :: blen

      ! Local Variables
      integer :: inb

      IF ((jcat <= 0) .and. (jcat /= -9)) RETURN

      IF (jcat == icat) THEN
         write(*,*) 'Warning: border finding error.'
      ENDIF

      IF (nnb > 0) THEN
         inb = findloc(nbindex(1:nnb),jcat,dim=1)
         IF (inb <= 0) THEN
            nnb = nnb + 1
            nbindex(nnb) = jcat
            inb = nnb
         ENDIF
      ELSE
         nnb = 1
         nbindex(nnb) = jcat
         inb = nnb
      ENDIF

      borderlen(inb) = borderlen(inb) + blen

   END SUBROUTINE acc_border

END MODULE catchment_mod
