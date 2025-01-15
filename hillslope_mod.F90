MODULE hillslope_mod

   USE river_lake_mod, only: binfo

CONTAINS

   SUBROUTINE get_hillslope_hydrounits (catsize, lakecellsize, nlev_max, maxhunum)

      USE task_mod
      USE utils_mod
      USE hydro_data_mod

      IMPLICIT NONE

      real    (kind=4), intent(in)  :: catsize
      real    (kind=4), intent(in)  :: lakecellsize
      integer (kind=4), intent(in)  :: nlev_max

      integer (kind=4), intent(inout) :: maxhunum

      real    (kind=4) :: levsize

      real    (kind=8), allocatable :: longitude(:)    
      real    (kind=8), allocatable :: latitude (:)    

      integer (kind=4), allocatable :: catch (:,:)    
      integer (kind=1), allocatable :: dir   (:,:)    
      real    (kind=4), allocatable :: hnd   (:,:)
      real    (kind=4), allocatable :: area  (:,:)
      integer (kind=4), allocatable :: hunit (:,:)    
      logical, allocatable :: hmask (:,:)

      integer (kind=4) :: icat
      integer (kind=4) :: catnum, lakeid, imin, imax, jmin, jmax
      integer (kind=4) :: bnds(6)
      logical :: end_of_data

      integer (kind=4) :: np, mp, dsize 
      integer (kind=4) :: i, j, iblk, jblk, iblk1, jblk1
      integer (kind=4) :: i0, i1, j0, j1, il0, il1, jl0, jl1
      integer (kind=4) :: iwork, ndone
      integer (kind=4) :: jlon
      integer (kind=4) :: mesg(4), zero

      integer (kind=4), allocatable :: hulist(:)
      integer (kind=4) :: iloc


      ! Dividing catchment into units.
      IF (p_is_master) THEN

         write(*,'(/3A/)') 'Step 3 ',trim(binfo),': Finding all hillslope units and lake elements ...'

         icat = 1
         ndone = 0

         DO WHILE (.true.)

            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, 0, p_comm_work, p_stat, p_err)

            iwork = mesg(1)

            IF (icat <= thisinfo%ntotalcat) THEN

               catnum = icat + thisinfo%icatdsp
               CALL mpi_send (catnum, 1, MPI_INTEGER, iwork, 1, p_comm_work, p_err) 
               CALL mpi_send (thisinfo%lake_id   (icat), 1, MPI_INTEGER, iwork, 1, p_comm_work, p_err) 
               CALL mpi_send (thisinfo%bsn_bnds(:,icat), 4, MPI_INTEGER, iwork, 1, p_comm_work, p_err) 

               write(*,100) trim(binfo), catnum, thisinfo%bsn_bnds(:,icat), icat, thisinfo%ntotalcat
               100 format('(S3) Hillslopes and Lake Elements ',A,': ', I7, ' (',I6,',',I6,',',I6,',',I6,'),', &
                  '(', I8, '/', I8, ') in progress.') 

               icat = icat + 1
            ELSE
               zero = 0
               CALL mpi_send (zero, 1, MPI_INTEGER, iwork, 1, p_comm_work, p_err) 
               ndone = ndone + 1
            ENDIF

            maxhunum = max(mesg(2), maxhunum)

            IF (ndone == p_nwork) EXIT
         ENDDO

         CALL mpi_barrier (p_comm_work, p_err)
         CALL excute_data_task (t_exit)

      ELSEIF (p_is_data) THEN

         CALL data_daemon ()

      ELSEIF (p_is_work) THEN

         ! CALL sync_window ()

         levsize = catsize / nlev_max

         mesg(1:2) = (/p_iam_work, 0/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, 0, 0, p_comm_work, p_err) 

         DO WHILE (.true.)

            CALL mpi_recv (catnum, 1, MPI_INTEGER, 0, 1, p_comm_work, p_stat, p_err)

            IF (catnum == 0) EXIT

            CALL mpi_recv (lakeid,    1, MPI_INTEGER, 0, 1, p_comm_work, p_stat, p_err)
            CALL mpi_recv (mesg(1:4), 4, MPI_INTEGER, 0, 1, p_comm_work, p_stat, p_err)

            imin = mesg(1) 
            imax = mesg(2) 
            jmin = mesg(3) 
            jmax = mesg(4) 

            np = imax - imin + 1
            mp = jmax - jmin + 1
            IF (mp < 0) mp = mp + mglb
            dsize = np * mp

            allocate (longitude(mp))
            allocate (latitude (np))
            allocate (catch (np,mp))
            allocate (dir   (np,mp))
            allocate (hnd   (np,mp))

            CALL aggregate_data (imin, imax, jmin, jmax, & 
               np, mp, longitude, latitude, &
               icat = catch, dir = dir, hnd = hnd)

            allocate (area (np,mp))
            DO j = 1, mp
               DO i = 1, np
                  jlon = j + jmin - 1
                  IF (jlon > mglb) jlon = jlon - mglb

                  area (i,j) = get_area (i+imin-1,jlon)
               ENDDO
            ENDDO

            allocate (hmask (np,mp))
            allocate (hunit (np,mp))

            hmask = (catch == catnum) 
            IF (any(hmask)) THEN
               IF (lakeid <= 0) THEN
                  CALL divide_hillslope_into_hydrounits ( lakeid, &
                     np, mp, dir, hnd, area, hmask, levsize, hunit) 
               ELSE
                  CALL divide_lake (lakeid, np, mp, area, hmask, lakecellsize, hunit)
               ENDIF
            ENDIF

            iblk = (imin-1)/nbox + 1
            jblk = (jmin-1)/mbox + 1
            DO WHILE (.true.)
               CALL block_iterator (imin, imax, jmin, jmax, iblk, jblk, &
                  i0, i1, j0, j1, il0, il1, jl0, jl1, &
                  iblk1, jblk1, end_of_data)
                  
               bnds  = (/iblk,jblk,i0,i1,j0,j1/)

               dsize = (j1-j0+1)*(i1-i0+1)

               CALL excute_data_task (t_update_hunit)
               CALL mpi_send (catnum, 1, MPI_INTEGER, p_data_address, t_update_hunit, p_comm_glb, p_err) 
               CALL mpi_send (bnds,   6, MPI_INTEGER, p_data_address, t_update_hunit, p_comm_glb, p_err) 
               CALL mpi_send (hunit(il0:il1,jl0:jl1), dsize, MPI_INTEGER, p_data_address, &
                  t_update_hunit, p_comm_glb, p_err)

               IF (end_of_data) THEN
                  EXIT
               ELSE
                  iblk = iblk1; jblk = jblk1
               ENDIF
            ENDDO

            maxhunum = 0
            
            IF (lakeid <= 0) THEN
               allocate (hulist(np*mp))
               DO j = 1, mp
                  DO i = 1, np
                     IF (catch(i,j) == catnum) THEN
                        CALL insert_into_sorted_list1 (hunit(i,j), maxhunum, hulist, iloc)
                     ENDIF
                  ENDDO
               ENDDO
               deallocate(hulist)
            ENDIF

            mesg(1:2) = (/p_iam_work, maxhunum/)
            CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, 0, 0, p_comm_work, p_err) 

            deallocate (longitude)
            deallocate (latitude )
            deallocate (catch)
            deallocate (dir  )
            deallocate (hnd  )
            deallocate (area )
            deallocate (hmask)
            deallocate (hunit)

         ENDDO

         CALL mpi_barrier (p_comm_work, p_err)

      ENDIF

   END SUBROUTINE get_hillslope_hydrounits


   SUBROUTINE divide_hillslope_into_hydrounits ( lakeid, &
         np, mp, dir, hnd, area, hmask, levsize, &
         hunit)

      USE utils_mod
      USE hydro_data_mod

      ! Flow direction is prepared in 1-byte SIGNED integer (int8) and is defined as:
      !    1: east, 2: southeast, 4: south, 8: southwest, 16: west, 32: northwest, 64: north. 128: northeast
      !    0: river mouth, -1: inland depression, -9: undefined (ocean)
      ! IF a flow direction file is opened as UNSIGNED integer, undefined=247 and inland depression=255

      IMPLICIT NONE

      integer :: lakeid

      integer (kind = 4), intent(in) :: np, mp

      integer (kind = 1), intent(in) :: dir (np,mp)  ! flow direction
      real    (kind = 4), intent(in) :: hnd (np,mp)  ! height above nearest drainage
      real    (kind = 4), intent(in) :: area(np,mp)  ! area

      logical, intent(in) :: hmask (np,mp)

      real    (kind = 4), intent(in) :: levsize  ! max level size

      ! hillslope units
      integer (kind = 4), intent(inout) :: hunit (np,mp)    

      integer (kind = 4) :: nunit, iunit, junit, unit_next, unext

      integer (kind = 4) :: npxl, ipxl
      integer (kind = 4) :: i, j, inext, jnext, iloc, ij(2)
      logical :: is_new

      real    (kind = 4), allocatable :: h1d (:)
      integer (kind = 4), allocatable :: order1d   (:)

      integer (kind = 4), allocatable :: ij_sorted (:,:)

      logical :: msk   (np,mp)
      integer (kind = 4) :: order2d (np,mp)
      integer (kind = 4) :: hsgrp   (np,mp)    

      real (kind = 4) :: harea, area_this
      real (kind = 4), allocatable :: area_unit(:)

      logical :: found
      integer (kind = 4) :: fstloc, inb, jnb

      integer (kind = 4) :: hdown (np,mp)    
      integer (kind = 4) :: inode, nnode
      integer (kind = 4), allocatable :: route (:, :)

      logical (kind = 4), allocatable :: is_neighbour(:,:), small_units (:)
      integer (kind = 4) :: unitneighbour
      real    (kind = 4) :: areaneighbour

      integer (kind = 4), allocatable :: iunit_next (:)
      integer (kind = 4) :: nisland

      IF (lakeid == 0) THEN 
         ! lakeid = 0 refers to catchments with river pixels (hnd = 0)
         msk = hmask .and. (hnd > 0)
      ELSE
         ! lakeid < 0 refers to lake upstream areas without river pixels (hnd = 0)
         msk = hmask 
      ENDIF

      WHERE (hmask) 
         hunit = 0
      END WHERE
      WHERE (msk)   
         hunit = -1
      END WHERE

      harea = sum(area, mask=msk)

      IF (harea < levsize * 1.2) THEN
         WHERE (msk) 
            hunit = 1
         END WHERE
      ELSE

         ! divide catchment into hillslope levels
         npxl = count(msk)

         allocate (h1d     (npxl))
         allocate (order1d (npxl))

         allocate (ij_sorted (2,npxl))
         allocate (route     (2,npxl))

         h1d = pack(hnd, msk)

         order1d = (/ (ipxl, ipxl = 1, npxl) /)

         CALL quicksort (npxl, h1d, order1d)
         
         order1d(order1d) = (/ (ipxl, ipxl = 1, npxl) /)
         order2d = unpack(order1d, msk, -1)

         DO j = 1, mp
            DO i = 1, np
               IF (msk(i,j)) THEN
                  ij_sorted(:,order2d(i,j)) = (/i,j/)
               ENDIF
            ENDDO
         ENDDO

         IF (lakeid < 0) THEN
            area_this = 0
            DO j = 1, mp
               DO i = 1, np
                  CALL nextij (i,j, dir(i,j),inext,jnext, is_local = .true.)
                  IF ((inext >= 1) .and. (inext <= np) .and. (jnext >= 1) .and. (jnext <= mp)) THEN 
                     IF (.not. hmask(inext,jnext)) THEN
                        hunit(i,j) = 1
                        area_this = area_this + area(i,j)
                     ENDIF
                  ELSE
                     hunit(i,j) = 1
                     area_this = area_this + area(i,j)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

         nunit = 0
         fstloc = 1
         DO WHILE (any(hmask .and. (hunit == -1)))

            nunit = nunit + 1
         
            ! find first pixel.
            DO ipxl = fstloc, npxl
               i = ij_sorted(1,ipxl)
               j = ij_sorted(2,ipxl)
               IF (hunit(i,j) == -1) THEN
                  CALL nextij (i,j, dir(i,j),inext,jnext, is_local = .true.)
                  IF (hunit(inext,jnext) /= -1) THEN
                     EXIT
                  ENDIF
               ENDIF
            ENDDO

            hunit(i,j) = nunit

            IF (.not. ((lakeid < 0) .and. (nunit == 1))) THEN
               area_this = area(i,j)
               unit_next = hunit(inext,jnext)
            ELSE
               area_this = area_this + area(i,j)
               unit_next = -2
            ENDIF

            DO WHILE (area_this < levsize)
               
               i = ij_sorted(1,ipxl)
               j = ij_sorted(2,ipxl)

               found = .false.
               DO jnb = max(j-1,1), min(j+1,mp)
                  DO inb = max(i-1,1), min(i+1,np)
                     IF ((inb /= i) .or. (jnb /= j)) THEN
                        IF (hmask(inb,jnb) .and. (hunit(inb,jnb) == -1)) THEN

                           CALL nextij (inb,jnb, dir(inb,jnb),inext,jnext, is_local = .true.)

                           IF ((hunit(inext,jnext) == nunit) .or. (hunit(inext,jnext) == unit_next)) THEN
                              IF (order2d(inb,jnb) < ipxl) THEN
                                 found = .true.
                                 ipxl = order2d(inb,jnb)
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO

               IF (.not. found) THEN

                  ipxl = ipxl + 1

                  DO WHILE (ipxl <= npxl)

                     i = ij_sorted(1,ipxl)
                     j = ij_sorted(2,ipxl)

                     IF (hunit(i,j) == -1) THEN

                        CALL nextij (i,j, dir(i,j),inext,jnext, is_local = .true.)

                        IF (hunit(inext,jnext) == nunit) THEN
                           found = .true.
                        ELSEIF (hunit(inext,jnext) == unit_next) THEN
                           DO jnb = max(j-1,1), min(j+1,mp)
                              DO inb = max(i-1,1), min(i+1,np)
                                 IF ((inb /= i) .or. (jnb /= j)) THEN
                                    IF (hmask(inb,jnb) .and. (hunit(inb,jnb) == nunit)) THEN
                                       found = .true.
                                       EXIT
                                    ENDIF
                                 ENDIF
                              ENDDO
                              IF (found) EXIT
                           ENDDO
                        ENDIF

                     ENDIF

                     IF (found) THEN
                        EXIT
                     ELSE
                        ipxl = ipxl + 1
                     ENDIF 
                  ENDDO

               ENDIF

               IF (found) THEN
                  i = ij_sorted(1,ipxl)
                  j = ij_sorted(2,ipxl)
                  hunit(i,j) = nunit
                  area_this = area_this + area(i,j)
               ELSE
                  EXIT
               ENDIF

            ENDDO
               
            DO WHILE ((fstloc < npxl) .and. (hunit(ij_sorted(1,fstloc),ij_sorted(2,fstloc)) /= -1))
               fstloc = fstloc + 1
            ENDDO

         ENDDO

         allocate (area_unit(1:nunit))
         area_unit(:) = 0
         DO j = 1, mp
            DO i = 1, np
               IF (hmask(i,j) .and. (hunit(i,j) > 0)) THEN
                  area_unit(hunit(i,j)) = area_unit(hunit(i,j)) + area(i,j)
               ENDIF
            ENDDO
         ENDDO

         allocate (iunit_next (nunit))
         iunit_next = -1
         hdown = -9
         DO j = 1, mp
            DO i = 1, np
               IF (hmask(i,j) .and. (hunit(i,j) > 0) .and. (hdown(i,j) == -9)) THEN

                  nnode = 1
                  route(:,nnode) = (/i,j/)

                  DO WHILE (.true.)
                     CALL nextij (route(1,nnode),route(2,nnode),&
                        dir(route(1,nnode),route(2,nnode)),inext,jnext, is_local = .true.)

                     IF ((inext >= 1) .and. (inext <= np) .and. (jnext >= 1) .and. (jnext <= mp)) THEN
                        IF (hmask(inext,jnext) .and. (hunit(inext,jnext) == hunit(i,j)) &
                           .and. (hdown(inext,jnext) == -9)) THEN
                           nnode = nnode + 1
                           route(:,nnode) = (/inext,jnext/)
                        ELSE
                           IF (.not. hmask(inext,jnext)) THEN
                              hdown(route(1,nnode),route(2,nnode)) = 0
                           ELSEIF (hunit(inext,jnext) /= hunit(i,j)) THEN
                              hdown(route(1,nnode),route(2,nnode)) = hunit(inext,jnext)
                           ELSE
                              hdown(route(1,nnode),route(2,nnode)) = hdown(inext,jnext)
                              IF (iunit_next(hunit(i,j)) /= hdown(inext,jnext)) THEN
                                 write(*,*) 'Warning: more than one downstreams!'
                              ENDIF
                           ENDIF

                           EXIT
                        ENDIF
                     ELSE
                        hdown(route(1,nnode),route(2,nnode)) = 0
                        EXIT
                     ENDIF
                  ENDDO

                  iunit_next(hunit(i,j)) = hdown(route(1,nnode),route(2,nnode))
                  DO inode = nnode-1, 1, -1
                     hdown(route(1,inode),route(2,inode)) = hdown(route(1,inode+1), route(2,inode+1))
                  ENDDO

               ENDIF
            ENDDO
         ENDDO

         allocate (is_neighbour (nunit,nunit))
         is_neighbour(:,:) = .false.
         DO j = 1, mp
            DO i = 1, np
               IF (hmask(i,j) .and. (hunit(i,j) > 0)) THEN
                  DO jnext = max(1,j-1), min(mp,j+1)
                     DO inext = max(1,i-1), min(np,i+1)
                        IF (hmask(inext,jnext) .and. (hunit(inext,jnext) > 0) &
                           .and. (hunit(inext,jnext) /= hunit(i,j))) THEN
                           is_neighbour(hunit(i,j),hunit(inext,jnext)) = .true.
                           is_neighbour(hunit(inext,jnext),hunit(i,j)) = .true.
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         allocate (small_units (nunit))

         small_units = ((area_unit < levsize/2.0) .and. (area_unit > 0))
         DO WHILE (any(small_units))
            iunit = minloc(area_unit, mask=small_units, dim=1)
            areaneighbour = levsize * 1.2
            unitneighbour = 0
            DO junit = 1, nunit
               IF ((area_unit(junit) > 0) .and. (junit /= iunit) &
                  .and. (iunit_next(iunit) == iunit_next(junit)) &
                  .and. (is_neighbour(iunit,junit)) &
                  .and. (area_unit(junit) + area_unit(iunit) < levsize * 1.2) &
                  .and. (area_unit(junit) < areaneighbour)) THEN
                  unitneighbour = junit
                  areaneighbour = area_unit(junit)
               ENDIF
            ENDDO

            IF (unitneighbour > 0) THEN
               IF (iunit > unitneighbour) THEN
                  area_unit(unitneighbour) = area_unit(unitneighbour) + area_unit(iunit)
                  area_unit(iunit) = 0.

                  WHERE (iunit_next == iunit)
                     iunit_next = unitneighbour
                  END WHERE
                  iunit_next(iunit) = -1

                  is_neighbour(unitneighbour,:) = is_neighbour(unitneighbour,:) .or. is_neighbour(iunit,:)
                  is_neighbour(:,unitneighbour) = is_neighbour(:,unitneighbour) .or. is_neighbour(:,iunit)

                  WHERE (hmask .and. (hunit == iunit))
                     hunit = unitneighbour
                  END WHERE
               ELSE
                  area_unit(iunit) = area_unit(iunit) + area_unit(unitneighbour)
                  area_unit(unitneighbour) = 0.

                  WHERE (iunit_next == unitneighbour)
                     iunit_next = iunit
                  END WHERE
                  iunit_next(unitneighbour) = -1

                  is_neighbour(iunit,:) = is_neighbour(iunit,:) .or. is_neighbour(unitneighbour,:)
                  is_neighbour(:,iunit) = is_neighbour(:,iunit) .or. is_neighbour(:,unitneighbour)

                  WHERE (hmask .and. (hunit == unitneighbour))
                     hunit = iunit
                  END WHERE
               ENDIF

               small_units = ((area_unit < levsize/2.0) .and. (area_unit > 0))
            ELSE
               small_units(iunit) = .false.
            ENDIF
         ENDDO

         DO iunit = nunit, 1, -1
            IF ((iunit_next(iunit) > 0) .and. all(iunit_next /= iunit) &
               .and. (area_unit(iunit) < levsize * 0.1)) THEN
               junit = iunit_next(iunit)
               area_unit(junit) = area_unit(junit) + area_unit(iunit)
               area_unit(iunit) = 0.
               iunit_next(iunit) = -1
               WHERE (hmask .and. (hunit == iunit))
                  hunit = junit
               END WHERE
            ENDIF
         ENDDO

         iunit = 1
         DO WHILE (iunit <= nunit)
            IF (count(iunit_next == iunit) == 1) THEN
               junit = findloc(iunit_next, iunit, dim=1)
               IF (area_unit(iunit) + area_unit(junit) < levsize * 1.2) THEN
                  area_unit(iunit) = area_unit(iunit) + area_unit(junit)
                  area_unit(junit) = 0.
                  iunit_next(junit) = -1
                  WHERE (iunit_next == junit)
                     iunit_next = iunit
                  END WHERE
                  WHERE (hmask .and. (hunit == junit))
                     hunit = iunit
                  END WHERE
               ELSE
                  iunit = iunit + 1
               ENDIF
            ELSE
               iunit = iunit + 1
            ENDIF
         ENDDO

         junit = 0
         DO iunit = 1, nunit
            IF (area_unit(iunit) > 0) THEN
               junit = junit + 1
               WHERE (hmask .and. (hunit == iunit)) 
                  hunit = junit
               END WHERE
               WHERE (iunit_next == iunit)
                  iunit_next = junit
               END WHERE
               iunit_next(junit) = iunit_next(iunit)
               area_unit (junit) = area_unit (iunit)
            ENDIF
         ENDDO

         nunit = junit

         nisland = 0
         DO iunit = 1, nunit
            IF ((iunit_next(iunit) == 0) &
               .and. all(iunit_next(1:nunit) /= iunit) &
               .and. (area_unit(iunit) < levsize * 0.1)) THEN
               nisland = nisland + 1
               WHERE (hmask .and. (hunit == iunit)) 
                  hunit = -1
               END WHERE
            ENDIF
         ENDDO

         IF (nisland > 0) THEN
            junit = 0
            DO iunit = 1, nunit
               IF (.not. ((iunit_next(iunit) == 0) &
                  .and. all(iunit_next(1:nunit) /= iunit) &
                  .and. (area_unit(iunit) < levsize * 0.1))) THEN
                  junit = junit + 1
                  WHERE (hmask .and. (hunit == iunit)) 
                     hunit = junit
                  END WHERE
               ENDIF
            ENDDO

            WHERE (hmask .and. (hunit > 0  )) 
               hunit = hunit + 1
            END WHERE
            WHERE (hmask .and. (hunit == -1)) 
               hunit = 1
            END WHERE
         ENDIF

         deallocate (iunit_next)
         deallocate (area_unit )

         deallocate (route)

         deallocate (h1d    ) 
         deallocate (order1d)

         deallocate (ij_sorted)

         deallocate (is_neighbour)
         deallocate (small_units )

      ENDIF 

   END SUBROUTINE divide_hillslope_into_hydrounits

   SUBROUTINE divide_lake (lakeid, np, mp, area, hmask, lakecellsize, hunit)

      USE cvt_mod
      IMPLICIT NONE

      integer :: lakeid

      integer(kind=4), intent(in) :: np, mp
      real   (kind=4), intent(in) :: area (np,mp), lakecellsize
      logical,         intent(in) :: hmask(np,mp)

      integer(kind=4), intent(inout) :: hunit(np,mp)    

      ! Local Variables
      integer(kind=4) :: nplake, nlakeelm, iplake, i, j
      real   (kind=4) :: lakearea
      integer(kind=4) :: hunit_(np,mp)
      integer(kind=4), allocatable :: samplepoints(:,:), lakeindex(:)


      nplake   = count(hmask)
      lakearea = sum(area, mask = hmask)
      nlakeelm = max(nint(lakearea/lakecellsize), 1)

      allocate (samplepoints(2,nplake))

      iplake = 0
      DO j = 1, mp
         DO i = 1, np
            IF (hmask(i,j)) THEN
               iplake = iplake + 1
               samplepoints(:,iplake) = (/i, j/)
            ENDIF
         ENDDO
      ENDDO

      allocate (lakeindex (nplake))

      CALL cvt (lakeid, nplake, samplepoints, nlakeelm, lakeindex)

      hunit_ = hunit
      hunit  = unpack(lakeindex, mask = hmask, field = hunit_) 

      deallocate (samplepoints)
      deallocate (lakeindex)

   END SUBROUTINE divide_lake

END MODULE hillslope_mod
