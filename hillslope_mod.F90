MODULE hillslope_mod

CONTAINS

   subroutine get_hillslope_hydrounits (catsize, lakecellsize, nlev_max, maxhunum)

      use task_mod
      USE utils_mod
      use hydro_data_mod

      implicit none

      real    (kind=4), intent(in)  :: catsize
      real    (kind=4), intent(in)  :: lakecellsize
      integer (kind=4), intent(in)  :: nlev_max
      integer (kind=4), intent(out) :: maxhunum

      real    (kind=4) :: levsize

      REAL    (kind=8), allocatable :: longitude(:)    
      REAL    (kind=8), allocatable :: latitude (:)    

      integer (kind=4), allocatable :: catch (:,:)    
      integer (kind=1), allocatable :: dir   (:,:)    
      real    (kind=4), allocatable :: hnd   (:,:)
      real    (kind=4), allocatable :: area  (:,:)
      integer (kind=4), allocatable :: hunit (:,:)    
      logical, allocatable :: hmask (:,:)

      integer (kind=4) :: icat, jcat
      integer (kind=4) :: catnum, lakeid, imin, imax, jmin, jmax
      integer (kind=4) :: bnds(6)
      logical :: end_of_data

      integer (kind=4) :: np, mp, dsize 
      integer (kind=4) :: i, j, iblk, jblk, iblk1, jblk1
      integer (kind=4) :: i0, i1, j0, j1, il0, il1, jl0, jl1
      integer (kind=4) :: idata, iwork, ndone
      integer (kind=4) :: jlon
      integer (kind=4) :: mesg(4), zero

      integer (kind=4), allocatable :: hulist(:)
      integer (kind=4) :: iloc

      LOGICAL, allocatable :: catdone(:)

      ! Dividing catchment into units.
      if (p_is_master) then

         write(*,'(/A/)') 'Step 3 : Finding all hillslope units and lake elements ...'

         maxhunum = 0

         icat = 1
         ndone = 0

         allocate(catdone(ntotalcat))
         catdone = .false.
         jcat = 0

         DO WHILE (.true.)

            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, 0, p_comm_work, p_stat, p_err)

            iwork = mesg(1)

            IF (icat <= ntotalcat) THEN

               jcat = mod(jcat+1234, ntotalcat) + 1
               DO WHILE (catdone(jcat))
                  jcat = mod(jcat+1,ntotalcat)+1
               ENDDO
               catdone(jcat) = .true.

               call mpi_send (jcat, 1, MPI_INTEGER, iwork, 1, p_comm_work, p_err) 
               call mpi_send (lake_info_id   (jcat), 1, MPI_INTEGER, iwork, 1, p_comm_work, p_err) 
               call mpi_send (bsn_info_bnds(:,jcat), 4, MPI_INTEGER, iwork, 1, p_comm_work, p_err) 

               write(*,100) jcat, bsn_info_bnds(:,jcat), icat, ntotalcat
               100 format('(S3) Hillslopes and Lake Elements: ', I7, ' (',I6,',',I6,',',I6,',',I6,'),', &
                  '(', I8, '/', I8, ') in progress.') 

               icat = icat + 1
            ELSE
               zero = 0
               call mpi_send (zero, 1, MPI_INTEGER, iwork, 1, p_comm_work, p_err) 
               ndone = ndone + 1
            ENDIF

            maxhunum = max(mesg(2), maxhunum)

            IF (ndone == p_nwork) exit
         ENDDO

         call mpi_barrier (p_comm_work, p_err)
         call excute_data_task (t_exit)

      elseif (p_is_data) then

         call data_daemon ()

      elseif (p_is_work) then

         CALL sync_window ()

         levsize = catsize / nlev_max

         mesg(1:2) = (/p_iam_work, 0/)
         call mpi_send (mesg(1:2), 2, MPI_INTEGER, 0, 0, p_comm_work, p_err) 

         DO WHILE (.true.)

            CALL mpi_recv (catnum, 1, MPI_INTEGER, 0, 1, p_comm_work, p_stat, p_err)

            IF (catnum == 0) exit

            CALL mpi_recv (lakeid,    1, MPI_INTEGER, 0, 1, p_comm_work, p_stat, p_err)
            CALL mpi_recv (mesg(1:4), 4, MPI_INTEGER, 0, 1, p_comm_work, p_stat, p_err)

            imin = mesg(1) 
            imax = mesg(2) 
            jmin = mesg(3) 
            jmax = mesg(4) 

            np = imax - imin + 1
            mp = jmax - jmin + 1
            if (mp < 0) mp = mp + mglb
            dsize = np * mp

            allocate (longitude(mp))
            allocate (latitude (np))
            allocate (catch (np,mp))
            allocate (dir   (np,mp))
            allocate (hnd   (np,mp))

            call aggregate_data (imin, imax, jmin, jmax, & 
               np, mp, longitude, latitude, &
               icat = catch, dir = dir, hnd = hnd)

            allocate (area (np,mp))
            do j = 1, mp
               do i = 1, np
                  jlon = j + jmin - 1
                  if (jlon > mglb) jlon = jlon - mglb

                  area (i,j) = get_area (i+imin-1,jlon)
               end do
            end do

            allocate (hmask (np,mp))
            allocate (hunit (np,mp))

            hmask = (catch == catnum) 
            if (any(hmask)) then
               IF (lakeid <= 0) THEN
                  call divide_hillslope_into_hydrounits ( lakeid, &
                     np, mp, dir, hnd, area, hmask, levsize, hunit) 
               ELSE
                  CALL divide_lake (lakeid, np, mp, area, hmask, lakecellsize, hunit)
               ENDIF
            end if

            iblk = (imin-1)/nbox + 1
            jblk = (jmin-1)/mbox + 1
            do while (.true.)
               call block_iterator (imin, imax, jmin, jmax, iblk, jblk, &
                  i0, i1, j0, j1, il0, il1, jl0, jl1, &
                  iblk1, jblk1, end_of_data)

               idata = bkid(iblk,jblk)

               if (idata /= -1) then
                  bnds  = (/iblk,jblk,i0,i1,j0,j1/)
                  dsize = (j1-j0+1)*(i1-i0+1)

                  call excute_data_task (t_update_hunit, idata)
                  call mpi_send (catnum, 1, MPI_INTEGER, idata, t_update_hunit, p_comm_glb, p_err) 
                  call mpi_send (bnds,   6, MPI_INTEGER, idata, t_update_hunit, p_comm_glb, p_err) 
                  call mpi_send (hunit(il0:il1,jl0:jl1), dsize, MPI_INTEGER, idata, &
                     t_update_hunit, p_comm_glb, p_err)
               end if

               if (end_of_data) then
                  exit
               else
                  iblk = iblk1; jblk = jblk1
               end if
            end do

            maxhunum = 0
            
            IF (lakeid <= 0) THEN
               allocate (hulist(np*mp))
               do j = 1, mp
                  do i = 1, np
                     IF (catch(i,j) == catnum) THEN
                        CALL insert_into_sorted_list1 (hunit(i,j), maxhunum, hulist, iloc)
                     ENDIF
                  ENDDO
               ENDDO
               deallocate(hulist)
            ENDIF

            mesg(1:2) = (/p_iam_work, maxhunum/)
            call mpi_send (mesg(1:2), 2, MPI_INTEGER, 0, 0, p_comm_work, p_err) 

            deallocate (longitude)
            deallocate (latitude )
            deallocate (catch)
            deallocate (dir  )
            deallocate (hnd  )
            deallocate (area )
            deallocate (hmask)
            deallocate (hunit)

         end do

         call mpi_barrier (p_comm_work, p_err)

      end if

   end subroutine get_hillslope_hydrounits


   subroutine divide_hillslope_into_hydrounits ( lakeid, &
         np, mp, dir, hnd, area, hmask, levsize, &
         hunit)

      use utils_mod
      use hydro_data_mod

      ! Flow direction is prepared in 1-byte SIGNED integer (int8) and is defined as:
      !    1: east, 2: southeast, 4: south, 8: southwest, 16: west, 32: northwest, 64: north. 128: northeast
      !    0: river mouth, -1: inland depression, -9: undefined (ocean)
      ! If a flow direction file is opened as UNSIGNED integer, undefined=247 and inland depression=255

      implicit none

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

      integer (kind = 4) :: nborder
      integer (kind = 4), allocatable :: blist(:,:)

      integer (kind = 4) :: hdown (np,mp)    
      integer (kind = 4) :: inode, nnode
      integer (kind = 4), allocatable :: route (:, :)

      LOGICAL (kind = 4), allocatable :: is_neighbour(:,:), small_units (:)
      INTEGER (kind = 4) :: unitneighbour
      REAL    (kind = 4) :: areaneighbour

      integer (kind = 4), allocatable :: iunit_next (:)
      integer (kind = 4) :: nisland

      IF (lakeid == 0) THEN 
         ! lakeid = 0 refers to catchments with river pixels (hnd = 0)
         msk = hmask .and. (hnd > 0)
      ELSE
         ! lakeid < 0 refers to lake upstream areas without river pixels (hnd = 0)
         msk = hmask 
      ENDIF

      where (hmask) 
         hunit = 0
      END where
      where (msk)   
         hunit = -1
      END where

      harea = sum(area, mask=msk)

      if (harea < levsize * 1.2) then
         where (msk) 
            hunit = 1
         END where
      else

         ! divide catchment into hillslope levels
         npxl = count(msk)

         allocate (h1d     (npxl))
         allocate (order1d (npxl))

         allocate (ij_sorted (2,npxl))
         allocate (blist     (2,npxl))
         allocate (route     (2,npxl))

         h1d = pack(hnd, msk)

         order1d = (/ (ipxl, ipxl = 1, npxl) /)

         call quicksort (npxl, h1d, order1d)

         order1d(order1d) = (/ (ipxl, ipxl = 1, npxl) /)
         order2d = unpack(order1d, msk, -1)

         do j = 1, mp
            do i = 1, np
               if (msk(i,j)) then
                  ij_sorted(:,order2d(i,j)) = (/i,j/)
               end if
            end do
         end do

         IF (lakeid < 0) THEN
            area_this = 0
            do j = 1, mp
               do i = 1, np
                  call nextij (i,j, dir(i,j),inext,jnext, is_local = .true.)
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
         do while (any(hmask .and. (hunit == -1)))

            nunit = nunit + 1

            do ipxl = 1, npxl
               i = ij_sorted(1,ipxl)
               j = ij_sorted(2,ipxl)
               if (hunit(i,j) == -1) exit
            ENDDO

            hunit(i,j) = nunit

            IF (.not. ((lakeid < 0) .and. (nunit == 1))) THEN
               area_this = area(i,j)
            ELSE
               area_this = area_this + area(i,j)
            ENDIF

            CALL find_drainage_neighbour (np, mp, dir, hunit, hmask, i, j, route, nnode, unit_next)

            nborder = 0
            CALL find_nearest_pixels (np, mp, i, j, hmask, hunit, order2d, order1d, nborder, blist)

            DO WHILE (area_this < levsize)

               ipxl = 1
               DO WHILE (ipxl <= nborder)
                  IF (hunit(blist(1,ipxl),blist(2,ipxl)) == -1) THEN
                     CALL find_drainage_neighbour (np, mp, dir, hunit, hmask, &
                        blist(1,ipxl), blist(2,ipxl), route, nnode, unext)
                     IF ((unext == unit_next) .or. (unext == nunit)) THEN
                        DO inode = 1, nnode
                           i = route(1,inode)
                           j = route(2,inode)
                           IF (hunit(i,j) == -1) THEN
                              hunit(i,j) = nunit
                              area_this = area_this + area(i,j)
                              CALL find_nearest_pixels (np, mp, i, j, &
                                 hmask, hunit, order2d, order1d, nborder, blist)
                           ENDIF
                        ENDDO
                        exit
                     ELSE
                        ipxl = ipxl + 1
                     ENDIF
                  ELSE
                     ipxl = ipxl + 1
                  ENDIF
               ENDDO

               IF (ipxl > nborder) THEN
                  exit
               ENDIF
            ENDDO

         ENDDO

         allocate (area_unit(1:nunit))
         area_unit(:) = 0
         do j = 1, mp
            do i = 1, np
               IF (hmask(i,j) .and. (hunit(i,j) > 0)) THEN
                  area_unit(hunit(i,j)) = area_unit(hunit(i,j)) + area(i,j)
               ENDIF
            ENDDO
         ENDDO

         allocate (iunit_next (nunit))
         iunit_next = -1
         hdown = -9
         do j = 1, mp
            do i = 1, np
               if (hmask(i,j) .and. (hunit(i,j) > 0) .and. (hdown(i,j) == -9)) then

                  nnode = 1
                  route(:,nnode) = (/i,j/)

                  do while (.true.)
                     call nextij (route(1,nnode),route(2,nnode),&
                        dir(route(1,nnode),route(2,nnode)),inext,jnext, is_local = .true.)

                     IF ((inext >= 1) .and. (inext <= np) .and. (jnext >= 1) .and. (jnext <= mp)) THEN
                        if (hmask(inext,jnext) .and. (hunit(inext,jnext) == hunit(i,j)) &
                           .and. (hdown(inext,jnext) == -9)) then
                           nnode = nnode + 1
                           route(:,nnode) = (/inext,jnext/)
                        else
                           if (.not. hmask(inext,jnext)) then
                              hdown(route(1,nnode),route(2,nnode)) = 0
                           elseif (hunit(inext,jnext) /= hunit(i,j)) then
                              hdown(route(1,nnode),route(2,nnode)) = hunit(inext,jnext)
                           else
                              hdown(route(1,nnode),route(2,nnode)) = hdown(inext,jnext)
                              if (iunit_next(hunit(i,j)) /= hdown(inext,jnext)) then
                                 write(*,*) 'Warning: more than one downstreams!'
                              end if
                           end if

                           exit
                        end if
                     ELSE
                        hdown(route(1,nnode),route(2,nnode)) = 0
                        EXIT
                     ENDIF
                  end do

                  iunit_next(hunit(i,j)) = hdown(route(1,nnode),route(2,nnode))
                  do inode = nnode-1, 1, -1
                     hdown(route(1,inode),route(2,inode)) = hdown(route(1,inode+1), route(2,inode+1))
                  end do

               end if
            end do
         end do

         allocate (is_neighbour (nunit,nunit))
         is_neighbour(:,:) = .false.
         do j = 1, mp
            do i = 1, np
               if (hmask(i,j) .and. (hunit(i,j) > 0)) then
                  do jnext = max(1,j-1), min(mp,j+1)
                     do inext = max(1,i-1), min(np,i+1)
                        if (hmask(inext,jnext) .and. (hunit(inext,jnext) > 0) &
                           .and. (hunit(inext,jnext) /= hunit(i,j))) then
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
                  .and. (area_unit(junit) < areaneighbour)) then
                  unitneighbour = junit
                  areaneighbour = area_unit(junit)
               ENDIF
            ENDDO

            IF (unitneighbour > 0) THEN
               IF (iunit > unitneighbour) THEN
                  area_unit(unitneighbour) = area_unit(unitneighbour) + area_unit(iunit)
                  area_unit(iunit) = 0.

                  where (iunit_next == iunit)
                     iunit_next = unitneighbour
                  END where
                  iunit_next(iunit) = -1

                  is_neighbour(unitneighbour,:) = is_neighbour(unitneighbour,:) .or. is_neighbour(iunit,:)
                  is_neighbour(:,unitneighbour) = is_neighbour(:,unitneighbour) .or. is_neighbour(:,iunit)

                  where (hmask .and. (hunit == iunit))
                     hunit = unitneighbour
                  END where
               ELSE
                  area_unit(iunit) = area_unit(iunit) + area_unit(unitneighbour)
                  area_unit(unitneighbour) = 0.

                  where (iunit_next == unitneighbour)
                     iunit_next = iunit
                  END where
                  iunit_next(unitneighbour) = -1

                  is_neighbour(iunit,:) = is_neighbour(iunit,:) .or. is_neighbour(unitneighbour,:)
                  is_neighbour(:,iunit) = is_neighbour(:,iunit) .or. is_neighbour(:,unitneighbour)

                  where (hmask .and. (hunit == unitneighbour))
                     hunit = iunit
                  END where
               ENDIF

               small_units = ((area_unit < levsize/2.0) .and. (area_unit > 0))
            ELSE
               small_units(iunit) = .false.
            ENDIF
         ENDDO

         do iunit = nunit, 1, -1
            if ((iunit_next(iunit) > 0) .and. all(iunit_next /= iunit) &
               .and. (area_unit(iunit) < levsize * 0.1)) THEN
               junit = iunit_next(iunit)
               area_unit(junit) = area_unit(junit) + area_unit(iunit)
               area_unit(iunit) = 0.
               iunit_next(iunit) = -1
               where (hmask .and. (hunit == iunit))
                  hunit = junit
               END where
            end if
         end do

         iunit = 1
         DO WHILE (iunit <= nunit)
            IF (count(iunit_next == iunit) == 1) THEN
               junit = findloc(iunit_next, iunit, dim=1)
               IF (area_unit(iunit) + area_unit(junit) < levsize * 1.2) THEN
                  area_unit(iunit) = area_unit(iunit) + area_unit(junit)
                  area_unit(junit) = 0.
                  iunit_next(junit) = -1
                  where (iunit_next == junit)
                     iunit_next = iunit
                  END where
                  where (hmask .and. (hunit == junit))
                     hunit = iunit
                  END where
               ELSE
                  iunit = iunit + 1
               ENDIF
            ELSE
               iunit = iunit + 1
            ENDIF
         ENDDO

         junit = 0
         do iunit = 1, nunit
            if (area_unit(iunit) > 0) then
               junit = junit + 1
               where (hmask .and. (hunit == iunit)) 
                  hunit = junit
               END where
               where (iunit_next == iunit)
                  iunit_next = junit
               END where
               iunit_next(junit) = iunit_next(iunit)
               area_unit (junit) = area_unit (iunit)
            end if
         end do

         nunit = junit

         nisland = 0
         do iunit = 1, nunit
            if ((iunit_next(iunit) == 0) &
               .and. all(iunit_next(1:nunit) /= iunit) &
               .and. (area_unit(iunit) < levsize * 0.1)) then
               nisland = nisland + 1
               where (hmask .and. (hunit == iunit)) 
                  hunit = -1
               END where
            end if
         end do

         if (nisland > 0) then
            junit = 0
            do iunit = 1, nunit
               if (.not. ((iunit_next(iunit) == 0) &
                  .and. all(iunit_next(1:nunit) /= iunit) &
                  .and. (area_unit(iunit) < levsize * 0.1))) then
                  junit = junit + 1
                  where (hmask .and. (hunit == iunit)) 
                     hunit = junit
                  END where
               end if
            ENDDO

            where (hmask .and. (hunit > 0  )) 
               hunit = hunit + 1
            END where
            where (hmask .and. (hunit == -1)) 
               hunit = 1
            END where
         end if

         deallocate (iunit_next)
         deallocate (area_unit )

         deallocate (route)

         deallocate (h1d    ) 
         deallocate (order1d)

         deallocate (ij_sorted)
         deallocate (blist)

         deallocate (is_neighbour)
         deallocate (small_units )

      end if 

   end subroutine divide_hillslope_into_hydrounits

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

   SUBROUTINE find_nearest_pixels (np, mp, i, j, hm, hu, o2d, o1d, nb, blist)

      USE utils_mod
      IMPLICIT NONE

      integer (kind = 4), intent(in) :: np, mp
      integer (kind = 4), intent(in) :: i, j
      logical, intent(in) :: hm (:,:)
      integer (kind = 4), intent(in) :: hu (:,:)
      integer (kind = 4), intent(in) :: o2d(:,:)

      integer (kind = 4), intent(inout) :: o1d(:)
      integer (kind = 4), intent(inout) :: nb
      integer (kind = 4), intent(inout) :: blist(:,:)

      ! Local Variables
      integer (kind = 4) :: inext, jnext, iloc
      LOGICAL :: is_new

      do jnext = max(j-1,1), min(j+1,mp)
         do inext = max(i-1,1), min(i+1,np)
            if (hm(inext,jnext) .and. (hu(inext,jnext) == -1)) then
               call insert_into_sorted_list1 (o2d(inext,jnext), nb, o1d, iloc, is_new)

               if (is_new) then
                  if (iloc < nb) then
                     blist(:,iloc+1:nb) = blist(:,iloc:nb-1)
                  end if

                  blist(:,iloc) = (/inext,jnext/)
               end if
            end if

         end do
      end do

   END SUBROUTINE find_nearest_pixels

   SUBROUTINE find_drainage_neighbour (np, mp, dir, hu, hmsk, i, j, route, nnode, unext)

      use utils_mod
      use hydro_data_mod
      IMPLICIT NONE

      integer :: np, mp
      integer (kind = 1), intent(in) :: dir (:,:)  ! flow direction
      integer (kind = 4), intent(in) :: hu  (:,:)    
      logical,            intent(in) :: hmsk(:,:)    
      integer (kind = 4), intent(in) :: i, j

      integer (kind = 4), intent(inout) :: route(:,:)
      integer (kind = 4), intent(out)   :: nnode
      integer (kind = 4), intent(out)   :: unext

      ! Local Variables
      integer (kind = 4) :: inext, jnext

      nnode = 1
      route(:,nnode) = (/i,j/)
      do while (.true.)
         call nextij (route(1,nnode),route(2,nnode),&
            dir(route(1,nnode),route(2,nnode)),inext,jnext, is_local = .true.)

         IF ((inext >= 1) .and. (inext <= np) .and. (jnext >= 1) .and. (jnext <= mp)) THEN 
            if ((hu(inext,jnext) == -1) .and. hmsk(inext,jnext)) THEN
               nnode = nnode + 1
               route(:,nnode) = (/inext,jnext/)
            else
               IF (.not. hmsk(inext,jnext)) THEN
                  unext = -1
               ELSE
                  unext = hu(inext,jnext)
               ENDIF
               exit
            end if
         ELSE
            unext = -1
            EXIT
         ENDIF
      end do

   END SUBROUTINE find_drainage_neighbour 

END MODULE hillslope_mod
