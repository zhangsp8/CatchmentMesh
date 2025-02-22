MODULE information_mod

CONTAINS

   SUBROUTINE get_information (ntotalall, maxhunum, maxnnb)

      USE task_mod
      USE utils_mod
      USE hydro_data_mod
      USE catchment_mod

      IMPLICIT NONE

      integer (kind=4), intent(in)  :: ntotalall
      integer (kind=4), intent(in)  :: maxhunum
      integer (kind=4), intent(out) :: maxnnb

      ! Local Variables
      real    (kind=8), allocatable :: longitude(:)    
      real    (kind=8), allocatable :: latitude (:)    

      integer (kind=4), allocatable :: catch (:,:)    
      integer (kind=1), allocatable :: dir   (:,:)    
      real    (kind=4), allocatable :: hnd   (:,:)
      real    (kind=4), allocatable :: elv   (:,:)
      integer (kind=4), allocatable :: hunit (:,:)    
      logical, allocatable :: hmask (:,:)

      integer (kind=4) :: numblocks, ntotalcat
      integer (kind=4) :: catnum, lakeid, imin, imax, jmin, jmax, icat, jcat

      integer (kind=4) :: np, mp
      integer (kind=4) :: i, j, inext, jnext
      integer (kind=4) :: ilat, jlon, ilat0, jlon0, ilat1, jlon1

      integer (kind=4), allocatable :: indxhu (:)    
      real    (kind=4), allocatable :: areahu (:)
      integer (kind=4), allocatable :: npxlhu (:)    
      real    (kind=4), allocatable :: handhu (:)
      real    (kind=4), allocatable :: elvahu (:)
      integer (kind=4), allocatable :: nexthu (:)    
      real    (kind=4), allocatable :: plenhu (:)
      real    (kind=4), allocatable :: lfachu (:)
      logical, allocatable :: mask (:)

      real    (kind=4) :: elvacat

      integer (kind=4), allocatable :: route (:,:)    
      real    (kind=4), allocatable :: hdnxt (:,:)
      real    (kind=4), allocatable :: plen2 (:,:)

      integer (kind=4) :: nu, maxnum, unum, nusend, nurecv
      integer (kind=4) :: nnode, inode
      real    (kind=4) :: dist, dx1, dx2, dy1, dy2

      integer (kind=4) :: iwork, zero, ndone, mesg(5)
      integer (kind=4), allocatable :: datai(:)
      real    (kind=4), allocatable :: datar(:)

      integer (kind = 4) :: nnb
      integer (kind = 4), allocatable :: nbindex  (:)
      real    (kind = 4), allocatable :: lenborder(:)

      type basin_neighbour_type
         integer :: nnb
         integer, allocatable :: nbr_index(:)
         real*4 , allocatable :: lenborder(:)
      END type basin_neighbour_type

      type(basin_neighbour_type), allocatable :: bsn_nbr(:)

      integer :: inb, i_in_j
      integer, allocatable :: idxtmp(:)
      real*4 , allocatable :: lentmp(:)

      integer :: errorcode

      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) THEN

         write(*,'(/A/)') 'Step 4 : Get network information ...'

         numblocks = 0
         thisinfo => allinfo
         DO WHILE (associated(thisinfo))
            numblocks = numblocks + 1
            thisinfo => thisinfo%next
         ENDDO
   
         CALL mpi_bcast (numblocks, 1, MPI_INTEGER, p_master_address, p_comm_glb, p_err)

         maxnnb = 0

         thisinfo => allinfo
         DO WHILE (.true.)

            ntotalcat = thisinfo%ntotalcat

            allocate (thisinfo%hru_indx (maxhunum,ntotalcat))
            allocate (thisinfo%hru_area (maxhunum,ntotalcat))
            allocate (thisinfo%hru_next (maxhunum,ntotalcat))
            allocate (thisinfo%hru_hand (maxhunum,ntotalcat))
            allocate (thisinfo%hru_elva (maxhunum,ntotalcat))
            allocate (thisinfo%hru_plen (maxhunum,ntotalcat))
            allocate (thisinfo%hru_lfac (maxhunum,ntotalcat))

            allocate (thisinfo%bsn_num_hru (ntotalcat))
            allocate (thisinfo%bsn_elva    (ntotalcat))

            thisinfo%hru_indx = -1
            thisinfo%hru_area = 0
            thisinfo%hru_next = -1
            thisinfo%hru_hand = -1
            thisinfo%hru_elva = 0
            thisinfo%hru_plen = -1
            thisinfo%hru_lfac = 0.

            allocate (bsn_nbr (ntotalcat))

            icat = 1
            ndone = 0
            DO WHILE (.true.)

               CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, MPI_ANY_SOURCE, 0, p_comm_glb, p_stat, p_err)

               iwork  = mesg(1)
               catnum = mesg(2)
               IF (catnum > 0) THEN

                  jcat = catnum - thisinfo%icatdsp

                  lakeid = thisinfo%lake_id(jcat)

                  CALL mpi_recv (nurecv, 1, MPI_INTEGER, iwork, 1, p_comm_glb, p_stat, p_err)

                  thisinfo%bsn_num_hru(jcat) = nurecv

                  IF (lakeid <= 0) THEN
                     CALL mpi_recv (thisinfo%hru_indx(1:nurecv,jcat), nurecv, MPI_INTEGER, iwork, 1, p_comm_glb, p_stat, p_err)
                     CALL mpi_recv (thisinfo%hru_area(1:nurecv,jcat), nurecv, MPI_REAL4,   iwork, 1, p_comm_glb, p_stat, p_err)
                     CALL mpi_recv (thisinfo%hru_next(1:nurecv,jcat), nurecv, MPI_INTEGER, iwork, 1, p_comm_glb, p_stat, p_err)
                     CALL mpi_recv (thisinfo%hru_hand(1:nurecv,jcat), nurecv, MPI_REAL4,   iwork, 1, p_comm_glb, p_stat, p_err)
                     CALL mpi_recv (thisinfo%hru_elva(1:nurecv,jcat), nurecv, MPI_REAL4,   iwork, 1, p_comm_glb, p_stat, p_err)
                     CALL mpi_recv (thisinfo%hru_plen(1:nurecv,jcat), nurecv, MPI_REAL4,   iwork, 1, p_comm_glb, p_stat, p_err)
                     CALL mpi_recv (thisinfo%hru_lfac(1:nurecv,jcat), nurecv, MPI_REAL4,   iwork, 1, p_comm_glb, p_stat, p_err)
                  ENDIF

                  CALL mpi_recv (thisinfo%bsn_elva(jcat), 1, MPI_REAL4, iwork, 1, p_comm_glb, p_stat, p_err)

                  CALL mpi_recv (nnb, 1, MPI_INTEGER, iwork, 1, p_comm_glb, p_stat, p_err)
                  bsn_nbr(jcat)%nnb = nnb
                  IF (nnb > 0) THEN
                     maxnnb = max(nnb, maxnnb)
                     allocate (bsn_nbr(jcat)%nbr_index(nnb))
                     allocate (bsn_nbr(jcat)%lenborder(nnb))
                     CALL mpi_recv (bsn_nbr(jcat)%nbr_index, nnb, MPI_INTEGER, iwork, 1, p_comm_glb, p_stat, p_err)
                     CALL mpi_recv (bsn_nbr(jcat)%lenborder, nnb, MPI_REAL4  , iwork, 1, p_comm_glb, p_stat, p_err)
                  ENDIF

               ENDIF

               IF (icat <= ntotalcat) THEN

                  catnum = thisinfo%icatdsp+icat

                  CALL mpi_send (catnum, 1, MPI_INTEGER, iwork, 1, p_comm_glb, p_err) 
               
                  imin = thisinfo%bsn_nswe(1,icat) 
                  imax = thisinfo%bsn_nswe(2,icat) 
                  jmin = thisinfo%bsn_nswe(3,icat) 
                  jmax = thisinfo%bsn_nswe(4,icat) 
               
                  imin = max(imin-1, inorth)
                  imax = min(imax+1, isouth)

                  IF ((jwest == 1) .and. (jeast == mglb)) THEN
                     jmin = jmin-1
                     jmax = jmax+1
                  ELSE
                     IF (jmin /= jwest) jmin = jmin-1
                     IF (jmax /= jeast) jmax = jmax+1
                  ENDIF
                  IF (jmin == 0)   jmin = mglb
                  IF (jmax > mglb) jmax = 1

                  np = imax - imin + 1
                  mp = jmax - jmin + 1
                  IF (mp < 0) mp = mp + mglb

                  allocate (catch (np,mp))
                  allocate (dir   (np,mp))
                  allocate (hnd   (np,mp))
                  allocate (elv   (np,mp))
                  allocate (hunit (np,mp))

                  CALL aggregate_data (imin, imax, jmin, jmax, np, mp, &
                     icat = catch, dir = dir, hnd = hnd, elv = elv, hunit = hunit)

                  mesg(1:5) = (/thisinfo%lake_id (icat), imin, jmin, np, mp/)
                  CALL mpi_send (mesg(1:5), 5, MPI_INTEGER,  iwork, 1, p_comm_glb, p_err) 
                  CALL mpi_send (catch, np*mp, MPI_INTEGER,  iwork, 1, p_comm_glb, p_err) 
                  CALL mpi_send (dir,   np*mp, MPI_INTEGER1, iwork, 1, p_comm_glb, p_err) 
                  CALL mpi_send (hnd,   np*mp, MPI_REAL4,    iwork, 1, p_comm_glb, p_err) 
                  CALL mpi_send (elv,   np*mp, MPI_REAL4,    iwork, 1, p_comm_glb, p_err) 
                  CALL mpi_send (hunit, np*mp, MPI_INTEGER,  iwork, 1, p_comm_glb, p_err) 

                  write(*,100) catnum, ntotalall
                  100 format('(S4) Catchment information, ID (', I10, '/', I10, ') in progress.') 
                  
                  deallocate (catch)
                  deallocate (dir  )
                  deallocate (hnd  )
                  deallocate (elv  )
                  deallocate (hunit)

                  icat = icat + 1
               ELSE
                  zero = 0
                  CALL mpi_send (zero, 1, MPI_INTEGER, iwork, 1, p_comm_glb, p_err) 
                  ndone = ndone + 1
               ENDIF

               IF (ndone == p_nwork) EXIT
            ENDDO

            allocate (thisinfo%bsn_num_nbr (ntotalcat))
            allocate (thisinfo%bsn_idx_nbr (maxnnb,ntotalcat))
            allocate (thisinfo%bsn_len_bdr (maxnnb,ntotalcat))

            thisinfo%bsn_idx_nbr(:,:) = -1
            thisinfo%bsn_len_bdr(:,:) = 0

            DO icat = 1, ntotalcat
               thisinfo%bsn_num_nbr(icat) = bsn_nbr(icat)%nnb 
               IF (bsn_nbr(icat)%nnb > 0) THEN
                  thisinfo%bsn_idx_nbr(1:bsn_nbr(icat)%nnb,icat) = bsn_nbr(icat)%nbr_index 
                  thisinfo%bsn_len_bdr(1:bsn_nbr(icat)%nnb,icat) = bsn_nbr(icat)%lenborder 
               ENDIF
            ENDDO

            DO icat = 1, ntotalcat
               IF (bsn_nbr(icat)%nnb > 0) THEN
                  deallocate (bsn_nbr(icat)%nbr_index)
                  deallocate (bsn_nbr(icat)%lenborder)
               ENDIF
            ENDDO
            deallocate (bsn_nbr)

            IF (trim(storage_type) == 'block') THEN
               CALL flush_blocks(.false.)
            ENDIF

            CALL mpi_barrier (p_comm_glb, p_err)
         
            IF (associated(thisinfo%next)) THEN
               thisinfo => thisinfo%next
            ELSE
               EXIT
            ENDIF

         ENDDO

      ELSEIF (p_is_work) THEN

         CALL mpi_bcast (numblocks, 1, MPI_INTEGER, p_master_address, p_comm_glb, p_err)

         DO WHILE (numblocks > 0)

            mesg(1:2) = (/p_iam_glb, -1/)
            CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_master_address, 0, p_comm_glb, p_err) 

            DO WHILE (.true.)

               CALL mpi_recv (catnum, 1, MPI_INTEGER, p_master_address, 1, p_comm_glb, p_stat, p_err)

               IF (catnum == 0) EXIT

               CALL mpi_recv (mesg(1:5), 5, MPI_INTEGER, p_master_address, 1, p_comm_glb, p_stat, p_err)

               lakeid = mesg(1)
               imin   = mesg(2)
               jmin   = mesg(3)
               np     = mesg(4)
               mp     = mesg(5)

               allocate (catch (np,mp))
               allocate (dir   (np,mp))
               allocate (hnd   (np,mp))
               allocate (elv   (np,mp))
               allocate (hunit (np,mp))
               
               CALL mpi_recv (catch, np*mp, MPI_INTEGER,  p_master_address, 1, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (dir,   np*mp, MPI_INTEGER1, p_master_address, 1, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (hnd,   np*mp, MPI_REAL4,    p_master_address, 1, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (elv,   np*mp, MPI_REAL4,    p_master_address, 1, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (hunit, np*mp, MPI_INTEGER,  p_master_address, 1, p_comm_glb, p_stat, p_err)

               allocate (hmask (np,mp))
               hmask = (catch == catnum) 
               maxnum = maxval(hunit,mask=hmask)

               allocate (npxlhu (0:maxnum));  npxlhu(:) = 0
               DO j = 1, mp
                  DO i = 1, np
                     IF (hmask(i,j)) THEN
                        npxlhu(hunit(i,j)) = npxlhu(hunit(i,j)) + 1
                     ENDIF
                  ENDDO
               ENDDO

               nu = count(npxlhu > 0)
               deallocate (npxlhu)

               IF (lakeid <= 0) THEN

                  allocate (areahu (0:maxnum))
                  allocate (npxlhu (0:maxnum))
                  allocate (handhu (0:maxnum))
                  allocate (elvahu (0:maxnum))
                  allocate (nexthu (0:maxnum))
                  allocate (plenhu (0:maxnum))
                  allocate (lfachu (0:maxnum))

                  areahu(:) = 0
                  npxlhu(:) = 0
                  handhu(:) = 0
                  elvahu(:) = 0
                  nexthu(:) = -1
                  plenhu(:) = 0 
                  lfachu(:) = 0 

                  allocate (route (2,np*mp))
                  allocate (hdnxt (np,mp))
                  allocate (plen2 (np,mp))
                  plen2 = -1.

                  DO j = 1, mp
                     DO i = 1, np
                        IF (hmask(i,j)) THEN

                           unum = hunit(i,j)

                           npxlhu(unum) = npxlhu(unum) + 1

                           jlon = j + jmin - 1
                           IF (jlon > mglb) jlon = jlon - mglb
                           areahu(unum) = areahu(unum) + get_area (i+imin-1,jlon)

                           elvahu(unum) = elvahu(unum) + elv(i,j)
                           handhu(unum) = handhu(unum) + hnd(i,j)

                           IF (unum > 0) THEN

                              IF (plen2(i,j) == -1) THEN
                                 nnode = 1
                                 route(:,nnode) = (/i,j/)

                                 DO WHILE (.true.)
                                    CALL nextij (route(1,nnode),route(2,nnode),&
                                       dir(route(1,nnode),route(2,nnode)),inext,jnext)
                                    IF (hmask(inext,jnext) .and. (plen2(inext,jnext) == -1) &
                                       .and. (hunit(inext,jnext) == unum))THEN
                                       nnode = nnode + 1
                                       route(:,nnode) = (/inext,jnext/)
                                    ELSE
                                       EXIT
                                    ENDIF
                                 ENDDO

                                 DO inode = nnode, 1, -1
                                    ilat0 = route(1,inode) + imin - 1
                                    jlon0 = route(2,inode) + jmin - 1
                                    IF (jlon0 > mglb) jlon0 = jlon0 - mglb

                                    ilat1 = inext + imin - 1
                                    jlon1 = jnext + jmin - 1
                                    IF (jlon1 > mglb) jlon1 = jlon1 - mglb

                                    dist = dist_between (ilat0, jlon0, ilat1, jlon1)

                                    IF ((hunit(inext,jnext) /= unum) .or. (.not. hmask(inext,jnext))) THEN
                                       plen2(route(1,inode),route(2,inode)) = dist
                                    ELSE
                                       plen2(route(1,inode),route(2,inode)) = plen2(inext,jnext) + dist
                                    ENDIF

                                    inext = route(1,inode)
                                    jnext = route(2,inode)
                                 ENDDO
                              ENDIF

                              plenhu(unum) = plenhu(unum) + plen2(i,j)

                              CALL nextij (i, j, dir(i,j), inext, jnext)
                              IF ((hunit(inext,jnext) /= unum) .or. (.not. hmask(inext,jnext))) THEN
                                 IF (hmask(inext,jnext)) THEN
                                    IF ((nexthu(unum) /= -1) .and. (hunit(inext,jnext) /= nexthu(unum))) THEN
                                       write(*,*) 'Warning: hillslope more than one downstreams!', &
                                       catnum, unum, nexthu(unum), catch(inext,jnext), hunit(inext,jnext)
                                    ELSE
                                       nexthu(unum) = hunit(inext,jnext)
                                    ENDIF
                                 ELSE
                                    IF ((nexthu(unum) /= -1) .and. (nexthu(unum) /= -2)) THEN
                                       write(*,*) 'Warning: LakeCat more than one downstreams!', &
                                       catnum, unum, nexthu(unum), catch(inext,jnext), hunit(inext,jnext)
                                    ELSE
                                       nexthu(unum) = -2
                                    ENDIF
                                 ENDIF

                                 IF ((dir(i,j) == 1) .or. (dir(i,j) == 16)) THEN 
                                    lfachu(unum) = lfachu(unum) + dlon(i+imin-1)
                                 ELSEIF ((dir(i,j) == 4) .or. (dir(i,j) == 64)) THEN
                                    lfachu(unum) = lfachu(unum) + dlat(i+imin-1)
                                 ELSE
                                    lfachu(unum) = lfachu(unum) + sqrt(dlon(i+imin-1)**2 + dlat(i+imin-1)**2)
                                 ENDIF

                              ENDIF
                           ENDIF

                        ENDIF
                     ENDDO
                  ENDDO

                  WHERE (npxlhu > 0)
                     handhu = handhu / npxlhu
                     elvahu = elvahu / npxlhu
                     plenhu = plenhu / npxlhu
                  endwhere 

                  IF (count(nexthu == -2) > 1) THEN
                     write(*,*) 'Warning: more than one lowest hydro unit in ', catnum, imin, imax, jmin, jmax
                     CALL mpi_abort (p_comm_glb, errorcode, p_err)
                  ENDIF

               ENDIF

               IF (count(hmask) <= 0) write(*,*) catnum, lakeid, imin, imax, jmin, jmax, catch

               elvacat = sum(elv, mask = hmask) / count(hmask)

               CALL get_basin_neighbour (catnum, imin, jmin, np, mp, catch, &
                  nnb, nbindex, lenborder)

               mesg(1:2) = (/p_iam_glb, catnum/)
               CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, 0, 0, p_comm_glb, p_err) 

               CALL mpi_send (nu, 1, MPI_INTEGER, 0, 1, p_comm_glb, p_err) 

               IF (lakeid <= 0) THEN

                  allocate (mask (0:maxnum))

                  mask = npxlhu > 0
                  nusend = count(mask)

                  allocate (datar (nusend))
                  allocate (datai (nusend))

                  datai(1:nusend) = pack( (/(i,i=0,maxnum)/), mask)
                  CALL mpi_send (datai(1:nusend), nusend, MPI_INTEGER, p_master_address, 1, p_comm_glb, p_err) 

                  datar(1:nusend) = pack(areahu, mask)
                  CALL mpi_send (datar(1:nusend), nusend, MPI_REAL4, p_master_address, 1, p_comm_glb, p_err) 

                  datai(1:nusend) = pack(nexthu, mask)
                  CALL mpi_send (datai(1:nusend), nusend, MPI_INTEGER, p_master_address, 1, p_comm_glb, p_err) 

                  datar(1:nusend) = pack(handhu, mask)
                  CALL mpi_send (datar(1:nusend), nusend, MPI_REAL4, p_master_address, 1, p_comm_glb, p_err) 

                  datar(1:nusend) = pack(elvahu, mask)
                  CALL mpi_send (datar(1:nusend), nusend, MPI_REAL4, p_master_address, 1, p_comm_glb, p_err) 

                  datar(1:nusend) = pack(plenhu, mask)
                  CALL mpi_send (datar(1:nusend), nusend, MPI_REAL4, p_master_address, 1, p_comm_glb, p_err) 

                  datar(1:nusend) = pack(lfachu, mask)
                  CALL mpi_send (datar(1:nusend), nusend, MPI_REAL4, p_master_address, 1, p_comm_glb, p_err) 

               ENDIF

               CALL mpi_send (elvacat, 1, MPI_REAL4, p_master_address, 1, p_comm_glb, p_err) 

               CALL mpi_send (nnb, 1, MPI_INTEGER, p_master_address, 1, p_comm_glb, p_err) 
               IF (nnb > 0) THEN
                  CALL mpi_send (nbindex  (1:nnb), nnb, MPI_INTEGER, p_master_address, 1, p_comm_glb, p_err) 
                  CALL mpi_send (lenborder(1:nnb), nnb, MPI_REAL4  , p_master_address, 1, p_comm_glb, p_err) 
               ENDIF

               IF (allocated (catch)    ) deallocate (catch)
               IF (allocated (dir  )    ) deallocate (dir  )
               IF (allocated (hnd  )    ) deallocate (hnd  )
               IF (allocated (elv  )    ) deallocate (elv  )
               IF (allocated (hunit)    ) deallocate (hunit)
               IF (allocated (hmask)    ) deallocate (hmask)

               IF (allocated (areahu)   ) deallocate (areahu)
               IF (allocated (npxlhu)   ) deallocate (npxlhu)
               IF (allocated (handhu)   ) deallocate (handhu)
               IF (allocated (elvahu)   ) deallocate (elvahu)
               IF (allocated (nexthu)   ) deallocate (nexthu)
               IF (allocated (plenhu)   ) deallocate (plenhu)
               IF (allocated (lfachu)   ) deallocate (lfachu)
               IF (allocated (mask)     ) deallocate (mask)
               IF (allocated (datar)    ) deallocate (datar)
               IF (allocated (datai)    ) deallocate (datai)

               IF (allocated (route)    ) deallocate (route)
               IF (allocated (hdnxt)    ) deallocate (hdnxt)
               IF (allocated (plen2)    ) deallocate (plen2)

            ENDDO

            CALL mpi_barrier (p_comm_glb, p_err)

            numblocks = numblocks - 1

         ENDDO

      ENDIF

   END SUBROUTINE get_information

END MODULE information_mod
