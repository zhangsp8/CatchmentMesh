MODULE information_mod

CONTAINS

   subroutine get_information (maxhunum, maxnnb)

      use task_mod
      USE utils_mod
      use hydro_data_mod
      USE catchment_mod

      implicit none

      integer (kind=4), intent(in)  :: maxhunum
      integer (kind=4), intent(out) :: maxnnb

      ! Local Variables
      REAL    (kind=8), allocatable :: longitude(:)    
      REAL    (kind=8), allocatable :: latitude (:)    

      integer (kind=4), allocatable :: catch (:,:)    
      integer (kind=1), allocatable :: dir   (:,:)    
      real    (kind=4), allocatable :: hnd   (:,:)
      real    (kind=4), allocatable :: elv   (:,:)
      integer (kind=4), allocatable :: hunit (:,:)    
      logical, allocatable :: hmask (:,:)

      integer (kind=4) :: catnum, lakeid, imin, imax, jmin, jmax, icat

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
      LOGICAL, allocatable :: mask (:)

      real    (kind=4) :: elvacat

      integer (kind=4), allocatable :: route (:,:)    
      real    (kind=4), allocatable :: hdnxt (:,:)
      real    (kind=4), allocatable :: plen2 (:,:)

      integer (kind=4) :: nu, unum, nusend, nurecv
      integer (kind=4) :: nnode, inode
      real    (kind=4) :: dist, dx1, dx2, dy1, dy2

      integer (kind=4) :: iwork, zero, ndone, mesg(5)
      integer (kind=4), allocatable :: datai(:)
      real    (kind=4), allocatable :: datar(:)

      INTEGER (kind = 4) :: nnb
      integer (kind = 4), allocatable :: nbindex  (:)
      REAL    (kind = 4), allocatable :: lenborder(:)

      TYPE basin_neighbour_type
         INTEGER :: nnb
         INTEGER, allocatable :: nbr_index(:)
         REAL*4 , allocatable :: lenborder(:)
      END TYPE basin_neighbour_type

      TYPE(basin_neighbour_type), allocatable :: bsn_nbr(:)

      INTEGER :: inb, jcat, i_in_j
      INTEGER, allocatable :: idxtmp(:)
      REAL*4 , allocatable :: lentmp(:)

      integer :: fid, errorcode

      if (p_is_master) then

         write(*,'(/A/)') 'Step 4 : Get network information ...'

         allocate (hru_info_indx (maxhunum,ntotalcat))
         allocate (hru_info_area (maxhunum,ntotalcat))
         allocate (hru_info_next (maxhunum,ntotalcat))
         allocate (hru_info_hand (maxhunum,ntotalcat))
         allocate (hru_info_elva (maxhunum,ntotalcat))
         allocate (hru_info_plen (maxhunum,ntotalcat))
         allocate (hru_info_lfac (maxhunum,ntotalcat))

         allocate (bsn_info_elva (ntotalcat))

         hru_info_indx = -1
         hru_info_area = 0
         hru_info_next = -1
         hru_info_hand = -1
         hru_info_elva = 0
         hru_info_plen = -1
         hru_info_lfac = 0.

         allocate (indxhu (maxhunum))
         allocate (areahu (maxhunum))
         allocate (nexthu (maxhunum))
         allocate (handhu (maxhunum))
         allocate (elvahu (maxhunum))
         allocate (plenhu (maxhunum))
         allocate (lfachu (maxhunum))

         maxnnb = 0
         allocate (bsn_nbr (ntotalcat))

         open (unit=fid, file='mismatch.txt', status='replace')

         icat = 1
         ndone = 0
         DO WHILE (.true.)
            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, 0, p_comm_work, p_stat, p_err)

            iwork  = mesg(1)
            catnum = mesg(2)
            IF (catnum > 0) THEN
               lakeid = lake_info_id(catnum)
               IF (lakeid <= 0) THEN
                  CALL mpi_recv (nurecv, 1, MPI_INTEGER, iwork, 1, p_comm_work, p_stat, p_err)
                  CALL mpi_recv (indxhu(1:nurecv), nurecv, MPI_INTEGER, iwork, 1, p_comm_work, p_stat, p_err)
                  CALL mpi_recv (areahu(1:nurecv), nurecv, MPI_REAL4,   iwork, 1, p_comm_work, p_stat, p_err)
                  CALL mpi_recv (nexthu(1:nurecv), nurecv, MPI_INTEGER, iwork, 1, p_comm_work, p_stat, p_err)
                  CALL mpi_recv (handhu(1:nurecv), nurecv, MPI_REAL4,   iwork, 1, p_comm_work, p_stat, p_err)
                  CALL mpi_recv (elvahu(1:nurecv), nurecv, MPI_REAL4,   iwork, 1, p_comm_work, p_stat, p_err)
                  CALL mpi_recv (plenhu(1:nurecv), nurecv, MPI_REAL4,   iwork, 1, p_comm_work, p_stat, p_err)
                  CALL mpi_recv (lfachu(1:nurecv), nurecv, MPI_REAL4,   iwork, 1, p_comm_work, p_stat, p_err)
               ENDIF
                  
               CALL mpi_recv (bsn_info_elva(catnum), 1, MPI_REAL4, iwork, 1, p_comm_work, p_stat, p_err)

               CALL mpi_recv (nnb, 1, MPI_INTEGER, iwork, 1, p_comm_work, p_stat, p_err)
               bsn_nbr(catnum)%nnb = nnb
               IF (nnb > 0) THEN
                  maxnnb = max(nnb, maxnnb)
                  allocate (bsn_nbr(catnum)%nbr_index(nnb))
                  allocate (bsn_nbr(catnum)%lenborder(nnb))
                  CALL mpi_recv (bsn_nbr(catnum)%nbr_index, nnb, MPI_INTEGER, iwork, 1, p_comm_work, p_stat, p_err)
                  CALL mpi_recv (bsn_nbr(catnum)%lenborder, nnb, MPI_REAL4  , iwork, 1, p_comm_work, p_stat, p_err)
               ENDIF

               IF (lakeid <= 0) THEN
                  hru_info_indx(1:nurecv,catnum) = indxhu(1:nurecv)
                  hru_info_area(1:nurecv,catnum) = areahu(1:nurecv) 
                  hru_info_next(1:nurecv,catnum) = nexthu(1:nurecv)
                  hru_info_hand(1:nurecv,catnum) = handhu(1:nurecv)
                  hru_info_elva(1:nurecv,catnum) = elvahu(1:nurecv)
                  hru_info_plen(1:nurecv,catnum) = plenhu(1:nurecv)
                  hru_info_lfac(1:nurecv,catnum) = lfachu(1:nurecv)
               ENDIF

            ENDIF

            IF (icat <= ntotalcat) THEN

               call mpi_send (icat, 1, MPI_INTEGER, iwork, 2, p_comm_work, p_err) 
               call mpi_send (lake_info_id (icat),   1, MPI_INTEGER, iwork, 2, p_comm_work, p_err) 
               call mpi_send (bsn_info_bnds(:,icat), 4, MPI_INTEGER, iwork, 2, p_comm_work, p_err) 

               write(*,100) icat, ntotalcat
               100 format('(S4) Catchment information, ID (', I10, '/', I10, ') in progress.') 

               icat = icat + 1
            ELSE
               zero = 0
               call mpi_send (zero, 1, MPI_INTEGER, iwork, 2, p_comm_work, p_err) 
               ndone = ndone + 1
            ENDIF

            IF (ndone == p_nwork) exit
         ENDDO
         
         allocate (bsn_info_num_hru (ntotalcat))
         DO icat = 1, ntotalcat
            bsn_info_num_hru(icat) = count(hru_info_indx(:,catnum) >= 0)
         ENDDO

         DO icat = 1, ntotalcat
            DO inb = 1, bsn_nbr(icat)%nnb

               jcat = bsn_nbr(icat)%nbr_index(inb)

               IF (jcat > 0) THEN
                  IF (bsn_nbr(jcat)%nnb > 0) THEN

                     i_in_j = findloc(bsn_nbr(jcat)%nbr_index, icat, dim=1)

                     IF (i_in_j <= 0) THEN
                        allocate(idxtmp(bsn_nbr(jcat)%nnb))
                        allocate(lentmp(bsn_nbr(jcat)%nnb))
                        idxtmp = bsn_nbr(jcat)%nbr_index
                        lentmp = bsn_nbr(jcat)%lenborder

                        deallocate(bsn_nbr(jcat)%nbr_index)
                        deallocate(bsn_nbr(jcat)%lenborder)

                        bsn_nbr(jcat)%nnb = bsn_nbr(jcat)%nnb + 1
                        allocate(bsn_nbr(jcat)%nbr_index (bsn_nbr(jcat)%nnb))
                        allocate(bsn_nbr(jcat)%lenborder (bsn_nbr(jcat)%nnb))

                        bsn_nbr(jcat)%nbr_index(1:bsn_nbr(jcat)%nnb-1) = idxtmp 
                        bsn_nbr(jcat)%lenborder(1:bsn_nbr(jcat)%nnb-1) = lentmp 
                        bsn_nbr(jcat)%nbr_index(bsn_nbr(jcat)%nnb)   = icat
                        bsn_nbr(jcat)%lenborder(bsn_nbr(jcat)%nnb)   = 0 

                        i_in_j = bsn_nbr(jcat)%nnb
                        maxnnb = max(maxnnb, bsn_nbr(jcat)%nnb)

                        deallocate(idxtmp)
                        deallocate(lentmp)
                     ENDIF
                  ELSE
                     bsn_nbr(jcat)%nnb = 1
                     allocate (bsn_nbr(jcat)%nbr_index(1))
                     allocate (bsn_nbr(jcat)%lenborder(1))
                     i_in_j = 1
                     bsn_nbr(jcat)%nbr_index(1) = icat
                     bsn_nbr(jcat)%lenborder(1) = 0
                  ENDIF

                  write(*,'(A,I7,A,I7,A,E20.4,A,E20.4,A)') '(S4) Basin Pair: (', icat, ',', jcat, ')', &
                     bsn_nbr(icat)%lenborder(inb), '(->)', bsn_nbr(jcat)%lenborder(i_in_j), '(<-)'

                  IF (abs(bsn_nbr(icat)%lenborder(inb)-bsn_nbr(jcat)%lenborder(i_in_j)) &
                     > 0.1 * max(bsn_nbr(icat)%lenborder(inb),bsn_nbr(jcat)%lenborder(i_in_j))) THEN
                     write(fid,'(A,I7,A,I7,A,E20.4,A,E20.4,A)') '(S4) Border mismatch between : (', &
                        icat, ',', jcat, ')', &
                        bsn_nbr(icat)%lenborder(inb), '(->)', bsn_nbr(jcat)%lenborder(i_in_j), '(<-)'
                  ENDIF

                  bsn_nbr(icat)%lenborder(inb) = &
                     max(bsn_nbr(icat)%lenborder(inb), bsn_nbr(jcat)%lenborder(i_in_j))

                  bsn_nbr(jcat)%lenborder(i_in_j) = bsn_nbr(icat)%lenborder(inb)

               ELSE
                  write(*,'(A,I7,A,I7,A,E20.4,A)') '(S4) Basin to Ocean : (', icat, ',', jcat, ')', &
                     bsn_nbr(icat)%lenborder(inb), '(->)'
               ENDIF

            ENDDO
         ENDDO

         allocate (bsn_info_num_nbr (ntotalcat))
         allocate (bsn_info_idx_nbr (maxnnb,ntotalcat))
         allocate (bsn_info_len_bdr (maxnnb,ntotalcat))

         bsn_info_idx_nbr(:,:) = -1
         bsn_info_len_bdr(:,:) = 0

         DO icat = 1, ntotalcat
            DO inb = 1, bsn_nbr(icat)%nnb
               bsn_info_num_nbr(icat) = bsn_nbr(icat)%nnb 
               bsn_info_idx_nbr(1:bsn_nbr(icat)%nnb,icat) = bsn_nbr(icat)%nbr_index 
               bsn_info_len_bdr(1:bsn_nbr(icat)%nnb,icat) = bsn_nbr(icat)%lenborder 
            ENDDO
         ENDDO

         deallocate (areahu)
         deallocate (handhu)
         deallocate (elvahu)
         deallocate (nexthu)
         deallocate (plenhu)
         deallocate (lfachu)

         close(fid)

         call mpi_barrier (p_comm_work, p_err)
         call excute_data_task (t_exit)

      elseif (p_is_data) then

         call data_daemon ()

      elseif (p_is_work) then

         CALL sync_window ()

         mesg(1:2) = (/p_iam_work, -1/)
         call mpi_send (mesg(1:2), 2, MPI_INTEGER, 0, 0, p_comm_work, p_err) 

         DO WHILE (.true.)

            CALL mpi_recv (catnum, 1, MPI_INTEGER, 0, 2, p_comm_work, p_stat, p_err)

            IF (catnum == 0) exit

            CALL mpi_recv (lakeid, 1, MPI_INTEGER, 0, 2, p_comm_work, p_stat, p_err)

            CALL mpi_recv (mesg(1:4), 4, MPI_INTEGER, 0, 2, p_comm_work, p_stat, p_err)
            imin = mesg(1) 
            imax = mesg(2) 
            jmin = mesg(3) 
            jmax = mesg(4) 

            imin = max(imin-1, inorth)
            imax = min(imax+1, isouth)

            jmin = jmin-1
            jmax = jmax+1
            IF (jwest == 0) THEN
               IF (jmin == 0) jmin = mglb
               IF (jmax > mglb) jmax = 1
            ELSE
               IF (jmin == jwest-1) jmin = jwest
               IF (jmax == jeast+1) jmax = jeast
            ENDIF

            np = imax - imin + 1
            mp = jmax - jmin + 1
            if (mp < 0) mp = mp + mglb

            allocate (latitude (np))
            allocate (longitude(mp))
            allocate (catch (np,mp))
            allocate (dir   (np,mp))
            allocate (hnd   (np,mp))
            allocate (elv   (np,mp))
            allocate (hunit (np,mp))

            call aggregate_data (imin, imax, jmin, jmax, & 
               np, mp, longitude, latitude, &
               icat = catch, dir = dir, hnd = hnd, elv = elv, hunit = hunit)

            allocate (hmask (np,mp))
            hmask = (catch == catnum) 

            IF (lakeid <= 0) THEN

               nu = maxval(hunit, mask = hmask) + 1

               allocate (areahu (0:nu-1))
               allocate (npxlhu (0:nu-1))
               allocate (handhu (0:nu-1))
               allocate (elvahu (0:nu-1))
               allocate (nexthu (0:nu-1))
               allocate (plenhu (0:nu-1))
               allocate (lfachu (0:nu-1))

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

               do j = 1, mp
                  do i = 1, np
                     IF (hmask(i,j)) THEN

                        unum = hunit(i,j)

                        npxlhu(unum) = npxlhu(unum) + 1

                        jlon = j + jmin - 1
                        if (jlon > mglb) jlon = jlon - mglb
                        areahu(unum) = areahu(unum) + get_area (i+imin-1,jlon)

                        elvahu(unum) = elvahu(unum) + elv(i,j)
                        handhu(unum) = handhu(unum) + hnd(i,j)

                        IF (unum > 0) THEN

                           IF (plen2(i,j) == -1) THEN
                              nnode = 1
                              route(:,nnode) = (/i,j/)

                              do while (.true.)
                                 call nextij (route(1,nnode),route(2,nnode),&
                                    dir(route(1,nnode),route(2,nnode)),inext,jnext)
                                 IF (hmask(inext,jnext) .and. (plen2(inext,jnext) == -1) &
                                    .and. (hunit(inext,jnext) == unum))THEN
                                    nnode = nnode + 1
                                    route(:,nnode) = (/inext,jnext/)
                                 else
                                    exit
                                 ENDIF
                              ENDDO

                              do inode = nnode, 1, -1
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
                              end do
                           ENDIF

                           plenhu(unum) = plenhu(unum) + plen2(i,j)

                           CALL nextij (i, j, dir(i,j), inext, jnext)
                           IF ((hunit(inext,jnext) /= unum) .or. (.not. hmask(inext,jnext))) THEN
                              IF (hmask(inext,jnext)) THEN
                                 IF ((nexthu(unum) /= -1) .and. (hunit(inext,jnext) /= nexthu(unum))) THEN
                                    write(*,*) 'Warning: hillslope more than one downstreams!', &
                                    catnum, unum, nexthu(unum), hunit(inext,jnext)
                                 ELSE
                                    nexthu(unum) = hunit(inext,jnext)
                                 ENDIF
                              ELSE
                                 IF ((nexthu(unum) /= -1) .and. (nexthu(unum) /= -2)) THEN
                                    write(*,*) 'Warning: LakeCat more than one downstreams!', &
                                    catnum, unum, nexthu(unum), hunit(inext,jnext)
                                 ELSE
                                    nexthu(unum) = -2
                                 ENDIF
                              ENDIF

                              IF ((dir(i,j) == 1) .or. (dir(i,j) == 16)) THEN 
                                 lfachu(unum) = lfachu(unum) + dlon(i+imin-1)
                              elseif ((dir(i,j) == 4) .or. (dir(i,j) == 64)) THEN
                                 lfachu(unum) = lfachu(unum) + dlat(i+imin-1)
                              ELSE
                                 lfachu(unum) = lfachu(unum) + sqrt(dlon(i+imin-1)**2 + dlat(i+imin-1)**2)
                              ENDIF

                           ENDIF
                        ENDIF

                     ENDIF
                  ENDDO
               ENDDO

               where (npxlhu > 0)
                  handhu = handhu / npxlhu
                  elvahu = elvahu / npxlhu
                  plenhu = plenhu / npxlhu
               endwhere 

               IF (count(nexthu == -2) > 1) THEN
                  write(*,*) 'Error: more than one lowest hydro unit in ', catnum, imin, imax, jmin, jmax
                  CALL mpi_abort (p_comm_glb, errorcode, p_err)
               ENDIF

            ENDIF

            IF (count(hmask) <= 0) write(*,*) catnum, lakeid, imin, imax, jmin, jmax

            elvacat = sum(elv, mask = hmask) / count(hmask)

            CALL get_basin_neighbour (catnum, imin, jmin, np, mp, catch, &
               nnb, nbindex, lenborder)

            mesg(1:2) = (/p_iam_work, catnum/)
            call mpi_send (mesg(1:2), 2, MPI_INTEGER, 0, 0, p_comm_work, p_err) 

            IF (lakeid <= 0) THEN

               allocate (mask (0:nu-1))

               mask = npxlhu > 0
               nusend = count(mask)
               call mpi_send (nusend, 1, MPI_INTEGER, 0, 1, p_comm_work, p_err) 

               allocate (datar (nusend))
               allocate (datai (nusend))

               datai(1:nusend) = pack( (/(i,i=0,nu-1)/), mask)
               call mpi_send (datai(1:nusend), nusend, MPI_INTEGER, 0, 1, p_comm_work, p_err) 

               datar(1:nusend) = pack(areahu, mask)
               call mpi_send (datar(1:nusend), nusend, MPI_REAL4, 0, 1, p_comm_work, p_err) 

               datai(1:nusend) = pack(nexthu, mask)
               call mpi_send (datai(1:nusend), nusend, MPI_INTEGER, 0, 1, p_comm_work, p_err) 

               datar(1:nusend) = pack(handhu, mask)
               call mpi_send (datar(1:nusend), nusend, MPI_REAL4, 0, 1, p_comm_work, p_err) 

               datar(1:nusend) = pack(elvahu, mask)
               call mpi_send (datar(1:nusend), nusend, MPI_REAL4, 0, 1, p_comm_work, p_err) 

               datar(1:nusend) = pack(plenhu, mask)
               call mpi_send (datar(1:nusend), nusend, MPI_REAL4, 0, 1, p_comm_work, p_err) 

               datar(1:nusend) = pack(lfachu, mask)
               call mpi_send (datar(1:nusend), nusend, MPI_REAL4, 0, 1, p_comm_work, p_err) 

            ENDIF
               
            call mpi_send (elvacat, 1, MPI_REAL4, 0, 1, p_comm_work, p_err) 

            call mpi_send (nnb, 1, MPI_INTEGER, 0, 1, p_comm_work, p_err) 
            IF (nnb > 0) THEN
               call mpi_send (nbindex  (1:nnb), nnb, MPI_INTEGER, 0, 1, p_comm_work, p_err) 
               call mpi_send (lenborder(1:nnb), nnb, MPI_REAL4  , 0, 1, p_comm_work, p_err) 
            ENDIF

            IF (allocated (longitude)) deallocate (longitude)
            IF (allocated (latitude )) deallocate (latitude )
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

         end do

         call mpi_barrier (p_comm_work, p_err)

      end if

   end subroutine get_information

END MODULE information_mod
