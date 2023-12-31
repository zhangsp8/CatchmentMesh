MODULE catchment_mod

   use hydro_data_mod

   IMPLICIT NONE

   type catinfo_typ
      integer (kind=4), allocatable :: cid    (:)  ! catchment id
      integer (kind=4), allocatable :: lid    (:)  ! lake id
      integer (kind=4), allocatable :: nswe (:,:)  ! boundaries: north, south, west, east
      type(catinfo_typ), pointer :: next
   END type 

CONTAINS

   subroutine get_catchment (catsize)

      use task_mod
      USE cvt_mod
      implicit none

      real (kind=4), intent(in) :: catsize

      integer (kind=4) :: iproc, cdsp, iwork, mesg(2), isrc, nsend, nrecv, ndone
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
      
      integer (kind=4) :: south, north, west, east
      integer (kind=4), allocatable :: catdown(:)
      
      ! lake
      logical :: is_lake
      integer (kind=4) :: nsmallcat, nshorecat, nlakeall
      integer (kind=4) :: ilakecat, ismallcat, nc, ic, np, ip
      integer (kind=4) :: outlet(2), lakeid, nplake, iplake
      real    (kind=4) :: shorearea
      integer (kind=4) :: icatdsp, i_dn, j_dn, npshore
      integer (kind=4), allocatable :: lakepixels  (:,:), lakeindex(:), shoreindex(:)
      integer (kind=4), allocatable :: smalllist   (:,:) 
      integer (kind=4), allocatable :: shoreline   (:,:), samples(:,:), thiscat(:,:) 
      real    (kind=4), allocatable :: upalake(:), wt(:)
      logical, allocatable :: shoremask(:)
      logical :: nextp

      if (p_is_master) then

         write(*,'(/A/)') 'Step 2 : Finding all catchments and lakes ...'

         ncat = riv_info_pix (3,nrivpix)

         allocate (riv_info_stt(2,ncat))
         allocate (riv_info_end(2,ncat))

         cdsp = 0
         do iproc = 1, p_nwork
            ncat_work = ncat / p_nwork
            if (iproc <= mod(ncat,p_nwork)) then
               ncat_work = ncat_work + 1
            end if

            call mpi_send (ncat_work, 1, MPI_INTEGER4, iproc, iproc, p_comm_work, p_err) 

            head = cdsp
            if (ncat_work > 0) then
               irivseg = 0
               do while (head < nrivpix)
                  head = head + 1
                  riv_info_stt(:,riv_info_pix(3,head)) = riv_info_pix(1:2,head)
                  do while (head < nrivpix)
                     if (riv_info_pix(3,head+1) == riv_info_pix(3,head)) then
                        head = head + 1
                     else
                        exit
                     end if
                  end do
                  riv_info_end(:,riv_info_pix(3,head)) = riv_info_pix(1:2,head)

                  irivseg = irivseg + 1
                  if (irivseg == ncat_work) exit
               end do

               nsend = head - cdsp
               call mpi_send (nsend, 1, MPI_INTEGER4, iproc, iproc, p_comm_work, p_err) 
               call mpi_send (riv_info_pix(:,cdsp+1:head), nsend*3, MPI_INTEGER4, &
                  iproc, iproc, p_comm_work, p_err) 
            end if

            cdsp = head
         end do

         ! for lake
         ndone = 0
         DO WHILE (.true.)
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, 0, p_comm_work, p_stat, p_err)
            iwork = mesg(1)
            nexpd = mesg(2)
            IF (nexpd > 0) THEN
               call mpi_send (ncat, 1, MPI_INTEGER, iwork, 0, p_comm_work, p_err) 
               ncat = ncat + nexpd
            ELSE
               ndone = ndone + 1
               IF (ndone == p_nwork) EXIT
            ENDIF
         ENDDO

         CALL mpi_barrier (p_comm_work, p_err)

         ntotalcat = ncat
         allocate (lake_info_id    (ntotalcat))
         allocate (bsn_info_bnds (4,ntotalcat))

         DO iproc = 1, p_nwork
            CALL mpi_recv (nrecv, 1, MPI_INTEGER, iproc, 0, p_comm_work, p_stat, p_err)
            IF (nrecv > 0) THEN
               allocate (cat_id   (nrecv))
               allocate (lake_id  (nrecv))
               allocate (bndinfo(4,nrecv))

               CALL mpi_recv (cat_id,  nrecv,   MPI_INTEGER, iproc, 1, p_comm_work, p_stat, p_err)
               CALL mpi_recv (lake_id, nrecv,   MPI_INTEGER, iproc, 1, p_comm_work, p_stat, p_err)
               CALL mpi_recv (bndinfo, 4*nrecv, MPI_INTEGER, iproc, 1, p_comm_work, p_stat, p_err)

               DO icat = 1, nrecv
                  lake_info_id (cat_id(icat))   = lake_id(icat)
                  bsn_info_bnds(:,cat_id(icat)) = bndinfo(:,icat)
               ENDDO

               deallocate(cat_id)
               deallocate(lake_id)
               deallocate(bndinfo)
            ENDIF
         ENDDO

         call excute_data_task (t_exit)

      end if

      if (p_is_data) then

         call data_daemon ()

      end if

      if (p_is_work) then

         CALL sync_window ()

         call mpi_recv (ncat, 1, MPI_INTEGER, 0, p_iam_work, p_comm_work, p_stat, p_err)

         if (ncat > 0) then

            call mpi_recv (nrivp, 1, MPI_INTEGER, 0, p_iam_work, p_comm_work, p_stat, p_err)
            allocate (river (3,nrivp))
            call mpi_recv (river, nrivp*3, MPI_INTEGER, 0, p_iam_work, p_comm_work, p_stat, p_err)

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
            do while (head < nrivp)

               irivseg = irivseg + 1
               head = head + 1
               
               is_lake = (get_lake(river(1,head),river(2,head)) > 0)

               IF (.not. is_lake) THEN
               
                  tail = head
                  do while (head < nrivp)
                     if (river(3,head+1) == river(3,head)) then
                        head = head + 1
                     else
                        exit
                     end if
                  end do

                  npall = 0
                  do irivp = head, tail, -1

                     i = river(1,irivp)
                     j = river(2,irivp)

                     elv0 = get_elv(i,j)

                     nplist = 0
                     call append_plist3 (pixellist, nplist, i, j, elvdata, 0._4)

                     iplist = 1
                     do while (iplist <= nplist)
                        do idir = 1, 8
                           call nextij (pixellist(1,iplist), pixellist(2,iplist), &
                              ishftc(int(-128,1),idir), i_up, j_up)

                           if (is_feasible_step(pixellist(1,iplist),pixellist(2,iplist),i_up,j_up)) then
                              if (get_dir(i_up,j_up) == ishftc(int(8,1),idir)) then
                                 if (get_icat(i_up,j_up) == 0) then
                                    elv = get_elv(i_up, j_up)
                                    call append_plist3 (pixellist, nplist, i_up, j_up, elvdata, elv-elv0)
                                 end if
                              end if
                           end if
                        end do

                        iplist = iplist + 1
                     end do                     

                     npall = npall + nplist

                     catnum = river(3,head)
                     CALL update_catchment_pixels     (catnum, nplist, pixellist, elvdata)
                     CALL update_catchment_boundaries (nplist, pixellist, irivseg, networkinfo)
                     
                     networkinfo%cid(irivseg) = catnum
                     networkinfo%lid(irivseg) = 0

                  ENDDO
                     
                  write(*,100) p_iam_glb, networkinfo%cid(irivseg), networkinfo%nswe(:,irivseg), &
                     irivseg, ncat, npall, head-tail+1
                  100 format('(S2) Catchment : On', I4, ' ID ', I7, ' (',I6,',',I6,',',I6,',',I6,'),', &
                     ' (', I5, '/', I5, ' done), ', I10, ' Pixels, ', I3, ' Rivers') 

               ELSE
                  ! lake
                  
                  ! 1. lake area.
                  outlet = river(1:2,head)

                  lakeid = get_lake(outlet(1), outlet(2))

                  nplake = 0
                  CALL append_plist (lakepixels, nplake, outlet(1), outlet(2))

                  iplake = 1
                  DO WHILE (iplake <= nplake)
                     do idir = 1, 8
                        call nextij (lakepixels(1,iplake), lakepixels(2,iplake), &
                           ishftc(int(-128,1),idir), i_up, j_up)

                        if (is_feasible_step(lakepixels(1,iplake),lakepixels(2,iplake),i_up,j_up)) then
                           if (get_dir(i_up,j_up) == ishftc(int(8,1),idir)) then
                              if (get_lake(i_up,j_up) == lakeid) then
                                 call append_plist3 (lakepixels, nplake, i_up, j_up, elvdata, 0._4)
                              end if
                           end if
                        end if
                     end do

                     iplake = iplake + 1
                  end do                     
                     
                  catnum = river(3,head)
                  CALL update_catchment_pixels     (catnum, nplake, lakepixels, elvdata)
                  CALL update_catchment_boundaries (nplake, lakepixels, irivseg, networkinfo)
                     
                  networkinfo%cid(irivseg) = catnum
                  networkinfo%lid(irivseg) = lakeid
                  
                  write(*,101) p_iam_glb, networkinfo%cid(irivseg), networkinfo%nswe(:,irivseg), &
                     irivseg, ncat, nplake
                  101 format('(S2) Lake      : On', I4, ' ID ', I7, ' (',I6,',',I6,',',I6,',',I6,'),', &
                     ' (', I5, '/', I5, ' done), ', I10, ' Pixels') 

                  ! 2. small upstream catchment
                  allocate (upalake (nplake))
                  nsmallcat = 0
                  DO iplake = 1, nplake
                     upa_up(:) = 0
                     do idir = 1, 8
                        call nextij (lakepixels(1,iplake), lakepixels(2,iplake), &
                           ishftc(int(-128,1),idir), i_up, j_up)

                        if (is_feasible_step(lakepixels(1,iplake),lakepixels(2,iplake),i_up,j_up)) then
                           if (get_dir(i_up,j_up) == ishftc(int(8,1),idir)) then
                              if (get_lake(i_up,j_up) /= lakeid) then
                                 IF (get_upa(i_up,j_up) < catsize) THEN
                                    upa_up(idir) = get_upa(i_up,j_up)
                                 ENDIF
                              end if
                           end if
                        end if
                     end do

                     DO WHILE (sum(upa_up) > catsize)
                        idir = maxloc(upa_up, dim = 1)
                        call nextij (lakepixels(1,iplake), lakepixels(2,iplake), &
                           ishftc(int(-128,1),idir), i_up, j_up)
                        call append_plist (smalllist, nsmallcat, i_up, j_up)
                        upa_up(idir) = 0.
                     ENDDO
                                    
                     upalake(iplake) = sum(upa_up)
                  ENDDO

                  ! 3. catchments along shore lines.
                  shorearea = sum(upalake, mask = upalake > 0)
                  IF (shorearea > 0) THEN
                     nshorecat = max(min(nint(shorearea/catsize), count(upalake > 0)), 1)
                  ELSE
                     nshorecat = 0
                  ENDIF

                  nlakeall = nsmallcat + nshorecat
                  
                  IF (nlakeall > 0) THEN
                  
                     ! inquire catchment number displacement
                     mesg = (/p_iam_work, nlakeall/)
                     CALL mpi_send (mesg,    2, MPI_INTEGER, 0, 0, p_comm_work, p_err) 
                     CALL mpi_recv (icatdsp, 1, MPI_INTEGER, 0, 0, p_comm_work, p_stat, p_err)
                     
                     allocate (thislake%next)
                     thislake => thislake%next
                     thislake%next => null()

                     ilakecat = 0
                     allocate (thislake%cid    (nlakeall))
                     allocate (thislake%lid    (nlakeall))
                     allocate (thislake%nswe (4,nlakeall))

                     thislake%nswe(:,:) = 0

                     DO ismallcat = 1, nsmallcat

                        outlet = smalllist(:,ismallcat)

                        call nextij (outlet(1), outlet(2), get_dir(outlet(1),outlet(2)), i_dn, j_dn)

                        elv0 = get_elv(i_dn,j_dn)
                        elv  = get_elv(outlet(1),outlet(2))

                        nplist = 0
                        call append_plist3 (pixellist, nplist, outlet(1), outlet(2), elvdata, elv-elv0)

                        iplist = 1
                        do while (iplist <= nplist)
                           do idir = 1, 8
                              call nextij (pixellist(1,iplist), pixellist(2,iplist), &
                                 ishftc(int(-128,1),idir), i_up, j_up)

                              if (is_feasible_step(pixellist(1,iplist),pixellist(2,iplist),i_up,j_up)) then
                                 if (get_dir(i_up,j_up) == ishftc(int(8,1),idir)) then
                                    elv = get_elv(i_up, j_up)
                                    call append_plist3 (pixellist, nplist, i_up, j_up, elvdata, elv-elv0)
                                 end if
                              end if
                           end do

                           iplist = iplist + 1
                        end do                     

                        ilakecat = ilakecat + 1
                        catnum = icatdsp + ilakecat 

                        CALL update_catchment_pixels     (catnum, nplist, pixellist, elvdata)
                        CALL update_catchment_boundaries (nplist, pixellist, ilakecat, thislake)

                        thislake%cid(ilakecat) = catnum
                        thislake%lid(ilakecat) = -1
                  
                        write(*,102) p_iam_glb, thislake%cid(ilakecat), thislake%nswe(:,ilakecat), &
                           ilakecat, nlakeall, nplist
                        102 format('(S2) LakeCat   : On', I4, ' ID ', I7, ' (',I6,',',I6,',',I6,',',I6,'),', &
                           ' (', I5, '/', I5, ' done), ', I10, ' Pixels') 

                     ENDDO

                     CALL sync_window ()

                     IF (nshorecat > 0) THEN

                        allocate (shoremask (nplake)) 
                        shoremask = upalake > 0

                        npshore   = count(shoremask)
                        shorearea = sum(upalake, mask = shoremask)
                        nc = max(min(nint(shorearea/catsize), npshore), 1)

                        allocate (shoreline (2,npshore))
                        shoreline(1,:) = pack(lakepixels(1,1:nplake), shoremask)
                        shoreline(2,:) = pack(lakepixels(2,1:nplake), shoremask)

                        allocate (samples (2,npshore))
                        samples = shoreline

                        west = shoreline(2,1)
                        east = shoreline(2,1)
                        DO ip = 2, npshore
                           west = min_west(west, shoreline(2,ip))
                           east = max_east(east, shoreline(2,ip))
                        ENDDO
                        IF (east < west) THEN
                           DO ip = 1, npshore
                              IF (shoreline(2,ip) < west) THEN
                                 samples(2,ip) = samples(2,ip) + mglb
                              ENDIF
                           ENDDO
                        ENDIF

                        allocate (wt (npshore))
                        wt = pack(upalake, mask = shoremask)

                        allocate (shoreindex (npshore))

                        CALL cvt (0, npshore, samples, nc, shoreindex, weight = wt)

                        DO ic = 1, nc

                           np = count(shoreindex == ic)
                           allocate(thiscat (2,np))
                           thiscat(1,:) = pack(shoreline(1,:), shoreindex == ic)
                           thiscat(2,:) = pack(shoreline(2,:), shoreindex == ic)
                           
                           nplist = 0
                           iplist = 1
                           DO ip = 1, np

                              elv0 = get_elv(thiscat(1,ip),thiscat(2,ip))
                              nextp = .true.

                              do while ((iplist <= nplist) .or. nextp)
                                 IF (nextp) THEN
                                    i = thiscat(1,ip); j = thiscat(2,ip)
                                 ELSE
                                    i = pixellist(1,iplist); j = pixellist(2,iplist)
                                 ENDIF

                                 do idir = 1, 8
                                    call nextij (i, j, ishftc(int(-128,1),idir), i_up, j_up)
                              
                                    if (is_feasible_step(i,j,i_up,j_up)) then
                                       if (get_dir(i_up,j_up) == ishftc(int(8,1),idir)) then
                                          IF (get_icat(i_up,j_up) == 0) THEN
                                             elv = get_elv(i_up, j_up)
                                             call append_plist3 (pixellist, nplist, i_up, j_up, &
                                                elvdata, elv-elv0)
                                          ENDIF
                                       end if
                                    end if
                                 end do

                                 IF (nextp) THEN
                                    nextp = .false.
                                 ELSE
                                    iplist = iplist + 1
                                 ENDIF
                              end do                     
                           ENDDO
                        
                           ilakecat = ilakecat + 1
                           catnum = icatdsp + ilakecat
                           
                           CALL update_catchment_pixels     (catnum, nplist, pixellist, elvdata)
                           CALL update_catchment_boundaries (nplist, pixellist, ilakecat, thislake)

                           thislake%cid(ilakecat) = catnum
                           thislake%lid(ilakecat) = -1

                           write(*,103) p_iam_glb, thislake%cid(ilakecat), thislake%nswe(:,ilakecat), &
                              ilakecat, nlakeall, nplist
                           103 format('(S2) LakeCat   : On', I4, ' ID ', I7, ' (',I6,',',I6,',',I6,',', &
                              I6,'),', ' (', I5, '/', I5, ' done), ', I10, ' Pixels') 

                           deallocate (thiscat)

                        ENDDO
                           
                        deallocate (shoremask)
                        deallocate (shoreline)
                        deallocate (samples)
                        deallocate (wt)
                        deallocate (shoreindex)

                     ENDIF
                  ENDIF

                  deallocate (upalake)

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

         ncat = size(networkinfo%cid)
         thislake => networkinfo%next
         DO WHILE (associated(thislake))
            ncat = ncat + size(thislake%cid)
            thislake => thislake%next
         ENDDO
         
         call mpi_send (ncat, 1, MPI_INTEGER4, 0, 0, p_comm_work, p_err) 

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

            call mpi_send (cat_id,  ncat,   MPI_INTEGER4, 0, 1, p_comm_work, p_err) 
            call mpi_send (lake_id, ncat,   MPI_INTEGER4, 0, 1, p_comm_work, p_err) 
            call mpi_send (bndinfo, 4*ncat, MPI_INTEGER4, 0, 1, p_comm_work, p_err) 

         ENDIF
         
         IF (ncat > 0) THEN
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

      call mpi_barrier (p_comm_glb, p_err)

      if (p_is_master) then

         allocate (bsn_info_downstream(ntotalcat))

         DO iproc = 1, p_nwork
            CALL mpi_recv (nrecv, 1, MPI_INTEGER, iproc, 0, p_comm_work, p_stat, p_err)
            IF (nrecv > 0) THEN
               allocate(cat_id (ncat))
               allocate(catdown(ncat))

               CALL mpi_recv (cat_id,  nrecv, MPI_INTEGER, iproc, 1, p_comm_work, p_stat, p_err)
               CALL mpi_recv (catdown, nrecv, MPI_INTEGER, iproc, 1, p_comm_work, p_stat, p_err)

               DO icat = 1, nrecv
                  bsn_info_downstream(cat_id(icat)) = catdown(icat)
               ENDDO

               deallocate(cat_id)
               deallocate(catdown)
            ENDIF
         ENDDO
         
         call excute_data_task (t_exit)

      ENDIF

      IF (p_is_data) then
         call data_daemon ()
      ENDIF

      if (p_is_work) THEN 

         CALL sync_window ()

         if (ncat > 0) then

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
                        ! 0: river mouth; -1: inland depression; -9: ocean;
                        IF ((dir_this /= 0) .and. (dir_this /= -1) .and. (dir_this /= -9)) THEN
                           CALL nextij (i, j, dir_this, i_dn, j_dn)
                           IF (is_feasible_step(i,j,i_dn,j_dn)) then
                              IF (get_icat(i_dn,j_dn) /= catnum) THEN
                                 IF (catdown(icat) /= -9999) THEN
                                    IF (get_icat(i_dn,j_dn) /= catdown(icat)) THEN
                                       write(*,'(A,I7,A,I7,A,I7,A,I7,A)') &
                                          'Warning: more than ONE downstream catchment : ',  &
                                          cat_id(icat), ',', lake_id(icat), '(', catdown(icat), '->', &
                                          get_icat(i_dn,j_dn), ')'
                                    ENDIF
                                 ELSE
                                    catdown(icat) = get_icat(i_dn,j_dn)
                                 ENDIF
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
                           
               write(*,104) p_iam_glb, cat_id(icat), catdown(icat)
               104 format('(S2) Downstream : On', I4, ' From ', I7, ' to ', I7)
            ENDDO
            
         ENDIF

         call mpi_send (ncat, 1, MPI_INTEGER4, 0, 0, p_comm_work, p_err) 
         IF (ncat > 0) THEN
            call mpi_send (cat_id,  ncat, MPI_INTEGER4, 0, 1, p_comm_work, p_err) 
            call mpi_send (catdown, ncat, MPI_INTEGER4, 0, 1, p_comm_work, p_err) 
         ENDIF

         IF (allocated(cat_id )) deallocate(cat_id )
         IF (allocated(bndinfo)) deallocate(bndinfo)
         IF (allocated(lake_id)) deallocate(lake_id)
         IF (allocated(catdown)) deallocate(catdown)

      ENDIF

      call mpi_barrier (p_comm_glb, p_err)

   end subroutine get_catchment 


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

      do ip = 1, np
         if (rmark(ip) == 0) then
            rmark(ip) = ip
            iblk = (pixellist(1,ip)-1)/nbox + 1
            jblk = (pixellist(2,ip)-1)/mbox + 1
            do jp = ip+1, np
               iblk0 = (pixellist(1,jp)-1)/nbox + 1
               jblk0 = (pixellist(2,jp)-1)/mbox + 1
               if ((iblk0 == iblk) .and. (jblk0 == jblk)) then
                  rmark(jp) = ip
               end if
            end do

            scnt = count(rmark(ip:np) == ip)
            allocate (ibuf (2,scnt))
            allocate (rbuf (scnt))

            scnt = 0
            do jp = ip, np
               if (rmark(jp) == ip) then
                  scnt = scnt + 1
                  ibuf(:,scnt) = pixellist(:,jp)
                  rbuf(scnt) = elvdata(jp)
               end if
            end do

            idata = bkid(iblk,jblk)
            call excute_data_task (t_update_icatch, idata)
            call mpi_send (scnt,   1, MPI_INTEGER4, idata, t_update_icatch, p_comm_glb, p_err) 
            call mpi_send (iblk,   1, MPI_INTEGER4, idata, t_update_icatch, p_comm_glb, p_err) 
            call mpi_send (jblk,   1, MPI_INTEGER4, idata, t_update_icatch, p_comm_glb, p_err) 
            call mpi_send (catnum, 1, MPI_INTEGER4, idata, t_update_icatch, p_comm_glb, p_err) 

            call mpi_send (ibuf, scnt*2, MPI_INTEGER, idata, t_update_icatch, p_comm_glb, p_err) 
            call mpi_send (rbuf, scnt,   MPI_REAL4,   idata, t_update_icatch, p_comm_glb, p_err) 

            deallocate (ibuf)
            deallocate (rbuf)
         end if
      end do

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
      REAL    (kind = 4), allocatable, intent(out) :: lenborder(:)

      ! Local Variable
      INTEGER :: npxl, istart, jstart, ithis, jthis, inext, jnext, idir
      INTEGER :: igthis, jgthis, ignext, jgnext
      REAL*4  :: blen
      INTEGER*1 :: sdir

      IF (any(catch(1,:) == icat)) THEN
         istart = 1
         jstart = findloc(catch(1,:), icat, dim=1)
      ELSE
         istart = 2
         jstart = findloc(catch(2,:), icat, dim=1)
      ENDIF

      npxl = count(catch == icat)

      nnb = 0

      IF (npxl > 1) THEN

         ithis = istart
         jthis = jstart

         allocate (nbindex  (np*mp-npxl))
         allocate (lenborder(np*mp-npxl))

         lenborder(:) = 0

         sdir = int(64,1)

         DO WHILE (.true.)
            do idir = 1, 8
               call nextij (ithis, jthis, ishftc(sdir,idir), inext, jnext, is_local=.true.)
               IF ((inext >= 1) .and. (inext <= np) .and. (jnext >= 1) .and. (jnext <= mp)) THEN
                  IF (catch(inext,jnext) == icat) THEN
                     exit
                  ENDIF
               ENDIF
            ENDDO

            igthis = ithis + imin - 1
            ignext = inext + imin - 1
            jgthis = jthis + jmin - 1
            IF (jgthis > mglb) jgthis = jgthis - mglb
            jgnext = jnext + jmin - 1
            IF (jgnext > mglb) jgnext = jgnext - mglb

            blen = dist_between (igthis, jgthis, ignext, jgnext)

            IF (abs(inext-ithis)+abs(jnext-jthis) == 2) THEN
               IF     ((inext == ithis+1) .and. (jnext == jthis+1)) THEN
                  CALL acc_border(icat, nnb, nbindex, catch(ithis,jnext), lenborder, blen)
               ELSEIF ((inext == ithis+1) .and. (jnext == jthis-1)) THEN
                  CALL acc_border(icat, nnb, nbindex, catch(inext,jthis), lenborder, blen)
               ELSEIF ((inext == ithis-1) .and. (jnext == jthis-1)) THEN
                  CALL acc_border(icat, nnb, nbindex, catch(ithis,jnext), lenborder, blen)
               ELSEIF ((inext == ithis-1) .and. (jnext == jthis+1)) THEN
                  CALL acc_border(icat, nnb, nbindex, catch(inext,jthis), lenborder, blen)
               ENDIF
            ELSEIF (abs(inext-ithis)+abs(jnext-jthis) == 1) THEN
               IF ((inext == ithis+1) .and. (jthis < mp)) THEN
                  CALL acc_border(icat, nnb, nbindex, catch(ithis,jthis+1), lenborder, blen/2.0_4)
                  CALL acc_border(icat, nnb, nbindex, catch(inext,jthis+1), lenborder, blen/2.0_4)
               ELSEIF ((inext == ithis-1) .and. (jthis > 1)) THEN
                  CALL acc_border(icat, nnb, nbindex, catch(ithis,jthis-1), lenborder, blen/2.0_4)
                  CALL acc_border(icat, nnb, nbindex, catch(inext,jthis-1), lenborder, blen/2.0_4)
               ELSEIF ((jnext == jthis+1) .and. (ithis > 1)) THEN
                  CALL acc_border(icat, nnb, nbindex, catch(ithis-1,jthis), lenborder, blen/2.0_4)
                  CALL acc_border(icat, nnb, nbindex, catch(ithis-1,jnext), lenborder, blen/2.0_4)
               ELSEIF ((jnext == jthis-1) .and. (ithis < np)) THEN
                  CALL acc_border(icat, nnb, nbindex, catch(ithis+1,jthis), lenborder, blen/2.0_4)
                  CALL acc_border(icat, nnb, nbindex, catch(ithis+1,jnext), lenborder, blen/2.0_4)
               ENDIF
            ENDIF

            IF ((inext == istart) .and. (jnext == jstart)) THEN
               exit
            ELSE
               ithis = inext
               jthis = jnext
               sdir = ishftc(sdir,mod(idir+4,8)+1)
            ENDIF

         ENDDO

      ENDIF

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

      REAL (kind = 4), intent(inout) :: borderlen(:)
      REAL (kind = 4), intent(in) :: blen

      ! Local Variables
      INTEGER :: inb

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
