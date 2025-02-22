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

      integer (kind=4) :: catnum, ncat, icat, i, j, i_up, j_up, idir, iblk, jblk
      integer (kind=1) :: dir_this
      
      integer (kind=4) :: irivp, irivseg, tail, head

      integer (kind=4) :: nplist, iplist, jplist, npall
      integer (kind=4), allocatable :: pixellist (:,:)

      real (kind=4) :: upa_up(8), elv0, elv
      real (kind=4), allocatable :: elvdata(:)

      type(catinfo_typ), target  :: networkinfo
      type(catinfo_typ), pointer :: thislake, nextlake
      
      integer (kind=4) :: south, north, west, east, cdthis
      
      ! lake
      logical :: is_lake
      integer (kind=4) :: nshorecat, nc, ic, np, ip, jp
      integer (kind=4) :: outlet(2), lakeid, nplake, iplake, i_dn, j_dn, npshore
      real    (kind=4) :: shorearea
      integer (kind=4), allocatable :: lakepixels(:,:), lakeindex (:)
      integer (kind=4), allocatable :: shoreline (:,:), shoreindex(:), thiscat(:,:) 
      real    (kind=4), allocatable :: upashore(:), hndshore(:)


      IF (p_is_master) THEN

         write(*,'(/3A/)') 'Step 2 ', trim(binfo), ': Finding all catchments and lakes ...'

         ncat = thisinfo%nrivseg

         IF (ncat > 0) THEN

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
            DO WHILE (head < thisinfo%nrivpix)

               irivseg = irivseg + 1
               head = head + 1
               
               is_lake = (get_lake(thisinfo%riv_pix(1,head),thisinfo%riv_pix(2,head)) > 0)

               IF (.not. is_lake) THEN
               
                  tail = head
                  DO WHILE (head < thisinfo%nrivpix)
                     IF (thisinfo%riv_pix(3,head+1) == thisinfo%riv_pix(3,head)) THEN
                        head = head + 1
                     ELSE
                        EXIT
                     ENDIF
                  ENDDO

                  npall = 0
                  DO irivp = head, tail, -1

                     i = thisinfo%riv_pix(1,irivp)
                     j = thisinfo%riv_pix(2,irivp)

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

                     catnum = thisinfo%riv_pix(3,head)
                     CALL update_catchment_pixels     (catnum, nplist, pixellist, elvdata)
                     CALL update_catchment_boundaries (nplist, pixellist, irivseg, networkinfo)
                     
                     networkinfo%cid(irivseg) = catnum
                     networkinfo%lid(irivseg) = 0

                  ENDDO
                     
                  write(*,100) trim(binfo), networkinfo%cid(irivseg), networkinfo%nswe(:,irivseg), &
                     irivseg, thisinfo%nrivseg, npall, head-tail+1
                  100 format('(S2) Catchment ',A,': ID ', I7, ' (',I6,',',I6,',',I6,',',I6,'),', &
                     ' (', I5, '/', I5, ' done), ', I10, ' Pixels, ', I5, ' Rivers') 

               ELSE
                  ! lake
                  
                  ! 1. lake area.
                  outlet = thisinfo%riv_pix(1:2,head)

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
                     
                  catnum = thisinfo%riv_pix(3,head)
                  CALL update_catchment_pixels     (catnum, nplake, lakepixels, elvdata)
                  CALL update_catchment_boundaries (nplake, lakepixels, irivseg, networkinfo)
                     
                  networkinfo%cid(irivseg) = catnum
                  networkinfo%lid(irivseg) = lakeid
                  
                  write(*,101) trim(binfo), networkinfo%cid(irivseg), networkinfo%nswe(:,irivseg), &
                     irivseg, thisinfo%nrivseg, nplake
                  101 format('(S2) Lake      ',A,': ID ', I7, ' (',I6,',',I6,',',I6,',',I6,'),', &
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

                        catnum = thisinfo%icatdsp + ncat + ic
                        CALL update_catchment_pixels     (catnum, nplist, pixellist, elvdata)
                        CALL update_catchment_boundaries (nplist, pixellist, ic, thislake)

                        thislake%cid(ic) = catnum
                        thislake%lid(ic) = -1
                     
                        write(*,102) trim(binfo), thislake%cid(ic), thislake%nswe(:,ic), ic, nshorecat, nplist
                        102 format('(S2) LakeCat   ',A,': ID ', I7, ' (',I6,',',I6,',',I6,',', &
                           I6,'),', ' (', I5, '/', I5, ' done), ', I10, ' Pixels') 

                        deallocate (thiscat)

                     ENDDO
                        
                     ncat = ncat + nshorecat
                           
                     deallocate (shoreline )
                     deallocate (shoreindex)

                  ENDIF

               ENDIF

            ENDDO

            deallocate (pixellist)
            deallocate (lakepixels)
            deallocate (elvdata)


            allocate (thisinfo%bsn_index  (ncat))
            allocate (thisinfo%lake_id    (ncat))
            allocate (thisinfo%bsn_nswe (4,ncat))
            
            thisinfo%bsn_index(1:ncat) = (/(ic, ic = 1, ncat)/) + thisinfo%icatdsp

            nc = size(networkinfo%cid)
            thisinfo%lake_id   (1:nc) = networkinfo%lid
            thisinfo%bsn_nswe(:,1:nc) = networkinfo%nswe

            ic = nc
            thislake => networkinfo%next
            DO WHILE (associated(thislake))
               nc = size(thislake%cid)
               
               thisinfo%lake_id   (ic+1:ic+nc) = thislake%lid
               thisinfo%bsn_nswe(:,ic+1:ic+nc) = thislake%nswe

               ic = ic + nc
               thislake => thislake%next
            ENDDO

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


            allocate (thisinfo%bsn_downstream (ncat))
            thisinfo%bsn_downstream(:) = -9999

            DO icat = 1, ncat
      
               catnum = icat + thisinfo%icatdsp
               north  = thisinfo%bsn_nswe(1,icat)
               south  = thisinfo%bsn_nswe(2,icat)
               west   = thisinfo%bsn_nswe(3,icat)
               east   = thisinfo%bsn_nswe(4,icat)

               IF ((south-1)/nbox - (north-1)/nbox > 1) THEN
                  
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
                                 IF (thisinfo%bsn_downstream(icat) /= -9999) THEN
                                    IF (cdthis /= thisinfo%bsn_downstream(icat)) THEN
                                       write(*,'(A,I7,A,I7,A,I7,A,I7,A)') &
                                          'Warning: more than ONE downstream catchment : ',  &
                                          catnum, ',', thisinfo%lake_id(icat), &
                                          '(', thisinfo%bsn_downstream(icat), '->', cdthis, ')'
                                    ENDIF
                                 ELSE
                                    thisinfo%bsn_downstream(icat) = cdthis
                                 ENDIF
                              ENDIF
                           ELSE
                              thisinfo%bsn_downstream(icat) = dir_this
                           ENDIF
                        ENDIF

                        IF (j == east) THEN
                           EXIT
                        ELSE
                           j = mod(j,mglb) + 1
                        ENDIF
                     ENDDO
                  ENDDO
                        
               ELSE

                  j = west
                  DO WHILE (.true.)
                     DO i = north, south
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
                                 IF (thisinfo%bsn_downstream(icat) /= -9999) THEN
                                    IF (cdthis /= thisinfo%bsn_downstream(icat)) THEN
                                       write(*,'(A,I7,A,I7,A,I7,A,I7,A)') &
                                          'Warning: more than ONE downstream catchment : ',  &
                                          catnum, ',', thisinfo%lake_id(icat), &
                                          '(', thisinfo%bsn_downstream(icat), '->', cdthis, ')'
                                    ENDIF
                                 ELSE
                                    thisinfo%bsn_downstream(icat) = cdthis
                                 ENDIF
                              ENDIF
                           ELSE
                              thisinfo%bsn_downstream(icat) = dir_this
                           ENDIF
                        ENDIF
                     ENDDO

                     IF (j == east) THEN
                        EXIT
                     ELSE
                        j = mod(j,mglb) + 1
                     ENDIF
                  ENDDO

               ENDIF
                           
               write(*,103) trim(binfo), catnum, thisinfo%bsn_downstream(icat)
               103 format('(S2) Downstream ', A, ': From ', I7, ' to ', I7)
            ENDDO
            
         ENDIF
            
         thisinfo%ntotalcat = ncat
                  
      ENDIF

   END SUBROUTINE get_catchment 


   SUBROUTINE update_catchment_pixels (catnum, np, pixellist, elvdata)

      IMPLICIT NONE

      integer (kind = 4), intent(in) :: catnum, np
      integer (kind = 4), intent(in) :: pixellist(:,:)
      real    (kind = 4), intent(in) :: elvdata  (:)

      ! Local
      integer :: ip, iblk, jblk, i, j

      DO ip = 1, np
         iblk = (pixellist(1,ip)-1)/nbox + 1
         jblk = (pixellist(2,ip)-1)/mbox + 1
         
         IF (.not. blks(iblk,jblk)%ready)  CALL readin_blocks(iblk, jblk)

         i = pixellist(1,ip) - blks(iblk,jblk)%idsp
         j = pixellist(2,ip) - blks(iblk,jblk)%jdsp
         IF ((blks(iblk,jblk)%icat (i,j) /= 0) &
            .and. (blks(iblk,jblk)%icat (i,j) /= catnum)) THEN
            write(*,*) 'mismatch while update catnum', &
               blks(iblk,jblk)%icat(i,j), catnum, pixellist(:,ip)
            STOP
         ENDIF
         blks(iblk,jblk)%icat (i,j) = catnum
         blks(iblk,jblk)%hnd  (i,j) = elvdata(ip)
      ENDDO

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

END MODULE catchment_mod
