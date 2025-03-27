MODULE hillslope_mod

   USE river_lake_mod, only: binfo

CONTAINS

   SUBROUTINE get_hillslope_hydrounits (catsize, lakecellsize, nlev_max, nhrumax)

      USE task_mod
      USE utils_mod
      USE hydro_data_mod

      IMPLICIT NONE

      real    (kind=4), intent(in)    :: catsize
      real    (kind=4), intent(in)    :: lakecellsize
      integer (kind=4), intent(in)    :: nlev_max
      integer (kind=4), intent(inout) :: nhrumax

      ! local
      integer (kind=4), allocatable :: catch (:,:)    
      integer (kind=1), allocatable :: dir   (:,:)    
      real    (kind=4), allocatable :: hnd   (:,:)
      real    (kind=4), allocatable :: elv   (:,:)
      real    (kind=4), allocatable :: area  (:,:)
      integer (kind=4), allocatable :: hunit (:,:)    
      logical, allocatable :: hmask (:,:)
      
      integer (kind=4), allocatable :: indxhu (:)    
      real    (kind=4), allocatable :: areahu (:)
      integer (kind=4), allocatable :: npxlhu (:)    
      real    (kind=4), allocatable :: handhu (:)
      real    (kind=4), allocatable :: elvahu (:)
      integer (kind=4), allocatable :: nexthu (:)    
      real    (kind=4), allocatable :: plenhu (:)
      real    (kind=4), allocatable :: lfachu (:)
      
      real    (kind=4), allocatable :: datar (:)
      integer (kind=4), allocatable :: datai (:)    
      
      integer (kind=4), allocatable :: route (:,:)    
      real    (kind=4), allocatable :: hdnxt (:,:)
      real    (kind=4), allocatable :: plen2 (:,:)

      real    (kind=4) :: levsize, elvacat, dist
      integer (kind=4) :: ntotalcat, icat, jcat, nhru, maxnum, unum
      integer (kind=4) :: catnum, lakeid, imin, imax, jmin, jmax, ilat0, ilat1, jlon0, jlon1
      integer (kind=4) :: np, mp, ip, jp, i, j, iblk, jblk, jlon, nnode, inode, inext, jnext
      integer (kind=4) :: iwork, mesg(5), zero, errorcode
      logical, allocatable :: mask(:)
   
      type hru_typ
         ! hydrounit information
         integer :: nhru
         integer (kind=4), allocatable :: indx (:)    
         real    (kind=4), allocatable :: area (:)
         real    (kind=4), allocatable :: hand (:)
         real    (kind=4), allocatable :: elva (:)
         integer (kind=4), allocatable :: next (:)    
         real    (kind=4), allocatable :: plen (:)
         real    (kind=4), allocatable :: lfac (:)
      END type hru_typ
      
      type(hru_typ), allocatable :: hruall (:)


      ! Dividing catchment into units.
      IF (p_is_master) THEN

         write(*,'(/3A/)') 'Step 3 ',trim(binfo),': Finding all hillslope units and lake elements ...'
            
         ntotalcat = thisinfo%ntotalcat
            
         allocate (hruall (ntotalcat))
         allocate (thisinfo%bsn_elva (ntotalcat))

         DO icat = 1, thisinfo%ntotalcat + p_nwork
            
            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, MPI_ANY_SOURCE, 0, p_comm_glb, p_stat, p_err)

            iwork  = mesg(1)
            catnum = mesg(2)

            IF (catnum > 0) THEN

               jcat = catnum - thisinfo%icatdsp

               imin = thisinfo%bsn_nswe(1,jcat) 
               imax = thisinfo%bsn_nswe(2,jcat) 
               jmin = thisinfo%bsn_nswe(3,jcat) 
               jmax = thisinfo%bsn_nswe(4,jcat) 
               
               CALL enlarge_nswe (imin, imax, jmin, jmax)

               np = imax - imin + 1
               mp = jmax - jmin + 1
               IF (mp < 0) mp = mp + mglb
            
               allocate (hunit (np,mp))

               CALL mpi_recv (hunit, np*mp, MPI_INTEGER, iwork, 2, p_comm_glb, p_stat, p_err)

               CALL mpi_recv (nhru, 1, MPI_INTEGER, iwork, 2, p_comm_glb, p_stat, p_err)

               hruall(jcat)%nhru = nhru
                  
               IF (thisinfo%lake_id(jcat) <= 0) THEN
               
                  allocate (hruall(jcat)%indx (nhru));  hruall(jcat)%indx = -1
                  allocate (hruall(jcat)%area (nhru));  hruall(jcat)%area = 0
                  allocate (hruall(jcat)%next (nhru));  hruall(jcat)%next = -1
                  allocate (hruall(jcat)%hand (nhru));  hruall(jcat)%hand = -1
                  allocate (hruall(jcat)%elva (nhru));  hruall(jcat)%elva = 0
                  allocate (hruall(jcat)%plen (nhru));  hruall(jcat)%plen = -1
                  allocate (hruall(jcat)%lfac (nhru));  hruall(jcat)%lfac = 0.

                  CALL mpi_recv (hruall(jcat)%indx, nhru, MPI_INTEGER, iwork, 2, p_comm_glb, p_stat, p_err)
                  CALL mpi_recv (hruall(jcat)%area, nhru, MPI_REAL4,   iwork, 2, p_comm_glb, p_stat, p_err)
                  CALL mpi_recv (hruall(jcat)%next, nhru, MPI_INTEGER, iwork, 2, p_comm_glb, p_stat, p_err)
                  CALL mpi_recv (hruall(jcat)%hand, nhru, MPI_REAL4,   iwork, 2, p_comm_glb, p_stat, p_err)
                  CALL mpi_recv (hruall(jcat)%elva, nhru, MPI_REAL4,   iwork, 2, p_comm_glb, p_stat, p_err)
                  CALL mpi_recv (hruall(jcat)%plen, nhru, MPI_REAL4,   iwork, 2, p_comm_glb, p_stat, p_err)
                  CALL mpi_recv (hruall(jcat)%lfac, nhru, MPI_REAL4,   iwork, 2, p_comm_glb, p_stat, p_err)

               ENDIF

               CALL mpi_recv (thisinfo%bsn_elva(jcat), 1, MPI_REAL4, iwork, 2, p_comm_glb, p_stat, p_err)

               DO ip = 1, np
                  DO jp = 1, mp
                     iblk = (imin+ip-2)/nbox + 1
                     jblk = mod(jmin+jp-2,mglb)/mbox + 1
                     i = mod(imin-1+ip-1,nbox) + 1
                     j = mod(jmin-1+jp-1,mbox) + 1
                     
                     IF (blks(iblk,jblk)%icat(i,j) == catnum) THEN

                        blks(iblk,jblk)%hunit(i,j) = hunit(ip,jp)

                     ENDIF
                  ENDDO
               ENDDO

               deallocate (hunit)

            ENDIF

            IF (icat <= thisinfo%ntotalcat) THEN
               
               imin = thisinfo%bsn_nswe(1,icat) 
               imax = thisinfo%bsn_nswe(2,icat) 
               jmin = thisinfo%bsn_nswe(3,icat) 
               jmax = thisinfo%bsn_nswe(4,icat) 

               CALL enlarge_nswe (imin, imax, jmin, jmax)

               np = imax - imin + 1
               mp = jmax - jmin + 1
               IF (mp < 0) mp = mp + mglb

               allocate (catch (np,mp))
               allocate (dir   (np,mp))
               allocate (hnd   (np,mp))
               allocate (elv   (np,mp))

               CALL aggregate_data (imin, imax, jmin, jmax, np, mp, &
                  icat = catch, dir = dir, hnd = hnd, elv = elv)
               
               catnum = icat + thisinfo%icatdsp
               CALL mpi_send (catnum, 1, MPI_INTEGER, iwork, 1, p_comm_glb, p_err) 
               
               mesg(1:5) = (/thisinfo%lake_id(icat), imin, jmin, np, mp/)
               CALL mpi_send (mesg(1:5), 5, MPI_INTEGER,  iwork, 1, p_comm_glb, p_err) 
               CALL mpi_send (catch, np*mp, MPI_INTEGER,  iwork, 1, p_comm_glb, p_err) 
               CALL mpi_send (dir,   np*mp, MPI_INTEGER1, iwork, 1, p_comm_glb, p_err) 
               CALL mpi_send (hnd,   np*mp, MPI_REAL4,    iwork, 1, p_comm_glb, p_err) 
               CALL mpi_send (elv,   np*mp, MPI_REAL4,    iwork, 1, p_comm_glb, p_err) 

               write(*,100) trim(binfo), catnum, thisinfo%bsn_nswe(:,icat), icat, thisinfo%ntotalcat
               100 format('(S3) Hillslopes and Lake Elements ',A,': ', I7, ' (',I6,',',I6,',',I6,',',I6,'),', &
                  '(', I8, '/', I8, ') in progress.') 

               deallocate (catch)
               deallocate (dir  )
               deallocate (hnd  )
               deallocate (elv  )

            ELSE
               zero = 0
               CALL mpi_send (zero, 1, MPI_INTEGER, iwork, 1, p_comm_glb, p_err) 
            ENDIF

         ENDDO
         
         allocate (thisinfo%bsn_num_hru (thisinfo%ntotalcat))
         DO icat = 1, thisinfo%ntotalcat
            thisinfo%bsn_num_hru(icat) = hruall(icat)%nhru
            IF (thisinfo%lake_id(icat) <= 0) THEN
               nhrumax = max(nhrumax, hruall(icat)%nhru)
            ENDIF
         ENDDO

         IF (nhrumax > 0) THEN

            allocate (thisinfo%hru_indx (nhrumax,thisinfo%ntotalcat)); thisinfo%hru_indx(:,:) = -1
            allocate (thisinfo%hru_area (nhrumax,thisinfo%ntotalcat)); thisinfo%hru_area(:,:) = 0.
            allocate (thisinfo%hru_hand (nhrumax,thisinfo%ntotalcat)); thisinfo%hru_hand(:,:) = 0.
            allocate (thisinfo%hru_elva (nhrumax,thisinfo%ntotalcat)); thisinfo%hru_elva(:,:) = 0.
            allocate (thisinfo%hru_next (nhrumax,thisinfo%ntotalcat)); thisinfo%hru_next(:,:) = -1
            allocate (thisinfo%hru_plen (nhrumax,thisinfo%ntotalcat)); thisinfo%hru_plen(:,:) = 0.
            allocate (thisinfo%hru_lfac (nhrumax,thisinfo%ntotalcat)); thisinfo%hru_lfac(:,:) = 0.
         
            DO icat = 1, thisinfo%ntotalcat
               IF (thisinfo%lake_id(icat) <= 0) THEN
                  nhru = hruall(icat)%nhru
                  thisinfo%hru_indx (1:nhru,icat) = hruall(icat)%indx
                  thisinfo%hru_area (1:nhru,icat) = hruall(icat)%area
                  thisinfo%hru_hand (1:nhru,icat) = hruall(icat)%hand
                  thisinfo%hru_elva (1:nhru,icat) = hruall(icat)%elva
                  thisinfo%hru_next (1:nhru,icat) = hruall(icat)%next
                  thisinfo%hru_plen (1:nhru,icat) = hruall(icat)%plen
                  thisinfo%hru_lfac (1:nhru,icat) = hruall(icat)%lfac
               ENDIF
            ENDDO
         ENDIF

         DO icat = 1, thisinfo%ntotalcat
            IF (allocated(hruall(icat)%indx)) deallocate(hruall(icat)%indx)
            IF (allocated(hruall(icat)%area)) deallocate(hruall(icat)%area)
            IF (allocated(hruall(icat)%next)) deallocate(hruall(icat)%next)
            IF (allocated(hruall(icat)%hand)) deallocate(hruall(icat)%hand)
            IF (allocated(hruall(icat)%elva)) deallocate(hruall(icat)%elva)
            IF (allocated(hruall(icat)%plen)) deallocate(hruall(icat)%plen)
            IF (allocated(hruall(icat)%lfac)) deallocate(hruall(icat)%lfac)
         ENDDO

         IF (allocated(hruall)) deallocate(hruall)

      ELSEIF (p_is_work) THEN

         levsize = catsize / nlev_max

         mesg(1:2) = (/p_iam_glb, 0/)
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
            
            CALL mpi_recv (catch, np*mp, MPI_INTEGER,  p_master_address, 1, p_comm_glb, p_stat, p_err)
            CALL mpi_recv (dir,   np*mp, MPI_INTEGER1, p_master_address, 1, p_comm_glb, p_stat, p_err)
            CALL mpi_recv (hnd,   np*mp, MPI_REAL4,    p_master_address, 1, p_comm_glb, p_stat, p_err)
            CALL mpi_recv (elv,   np*mp, MPI_REAL4,    p_master_address, 1, p_comm_glb, p_stat, p_err)

            allocate (area  (np,mp))
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

            maxnum = maxval(hunit,mask=hmask)

            allocate (npxlhu (0:maxnum));  npxlhu(:) = 0
            DO j = 1, mp
               DO i = 1, np
                  IF (hmask(i,j)) THEN
                     npxlhu(hunit(i,j)) = npxlhu(hunit(i,j)) + 1
                  ENDIF
               ENDDO
            ENDDO

            nhru = count(npxlhu > 0)

            IF (lakeid <= 0) THEN

               allocate (areahu (0:maxnum));  areahu(:) = 0.
               allocate (handhu (0:maxnum));  handhu(:) = 0.
               allocate (elvahu (0:maxnum));  elvahu(:) = 0.
               allocate (nexthu (0:maxnum));  nexthu(:) = -1
               allocate (plenhu (0:maxnum));  plenhu(:) = 0.
               allocate (lfachu (0:maxnum));  lfachu(:) = 0. 

               npxlhu(:) = 0

               allocate (route (2,np*mp))
               allocate (hdnxt (np,mp))
               allocate (plen2 (np,mp));  plen2 = -1.

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

            mesg(1:2) = (/p_iam_glb, catnum/)
            CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_master_address, 0, p_comm_glb, p_err) 

            CALL mpi_send (hunit, np*mp, MPI_INTEGER, p_master_address, 2, p_comm_glb, p_err) 

            CALL mpi_send (nhru, 1, MPI_INTEGER, p_master_address, 2, p_comm_glb, p_err) 

            IF (lakeid <= 0) THEN

               allocate (mask (0:maxnum))

               mask = npxlhu > 0

               allocate (datar (nhru))
               allocate (datai (nhru))

               datai(1:nhru) = pack( (/(i,i=0,maxnum)/), mask)
               CALL mpi_send (datai(1:nhru), nhru, MPI_INTEGER, p_master_address, 2, p_comm_glb, p_err) 

               datar(1:nhru) = pack(areahu, mask)
               CALL mpi_send (datar(1:nhru), nhru, MPI_REAL4, p_master_address, 2, p_comm_glb, p_err) 

               datai(1:nhru) = pack(nexthu, mask)
               CALL mpi_send (datai(1:nhru), nhru, MPI_INTEGER, p_master_address, 2, p_comm_glb, p_err) 

               datar(1:nhru) = pack(handhu, mask)
               CALL mpi_send (datar(1:nhru), nhru, MPI_REAL4, p_master_address, 2, p_comm_glb, p_err) 

               datar(1:nhru) = pack(elvahu, mask)
               CALL mpi_send (datar(1:nhru), nhru, MPI_REAL4, p_master_address, 2, p_comm_glb, p_err) 

               datar(1:nhru) = pack(plenhu, mask)
               CALL mpi_send (datar(1:nhru), nhru, MPI_REAL4, p_master_address, 2, p_comm_glb, p_err) 

               datar(1:nhru) = pack(lfachu, mask)
               CALL mpi_send (datar(1:nhru), nhru, MPI_REAL4, p_master_address, 2, p_comm_glb, p_err) 

            ENDIF

            CALL mpi_send (elvacat, 1, MPI_REAL4, p_master_address, 2, p_comm_glb, p_err) 

            IF (allocated (catch) ) deallocate (catch)
            IF (allocated (dir  ) ) deallocate (dir  )
            IF (allocated (hnd  ) ) deallocate (hnd  )
            IF (allocated (elv  ) ) deallocate (elv  )
            IF (allocated (area ) ) deallocate (area )
            IF (allocated (hunit) ) deallocate (hunit)
            IF (allocated (hmask) ) deallocate (hmask)

            IF (allocated (areahu)) deallocate (areahu)
            IF (allocated (npxlhu)) deallocate (npxlhu)
            IF (allocated (handhu)) deallocate (handhu)
            IF (allocated (elvahu)) deallocate (elvahu)
            IF (allocated (nexthu)) deallocate (nexthu)
            IF (allocated (plenhu)) deallocate (plenhu)
            IF (allocated (lfachu)) deallocate (lfachu)
            IF (allocated (mask)  ) deallocate (mask)
            IF (allocated (datar) ) deallocate (datar)
            IF (allocated (datai) ) deallocate (datai)

            IF (allocated (route) ) deallocate (route)
            IF (allocated (hdnxt) ) deallocate (hdnxt)
            IF (allocated (plen2) ) deallocate (plen2)

         ENDDO

      ENDIF
         
      CALL mpi_barrier (p_comm_glb, p_err)

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
      integer (kind = 4) :: i, j, inext, jnext, ilast, jlast, iloc, ij(2)
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
            
         IF (.not. any(hmask .and. (hunit == -1))) THEN
            RETURN
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

            nnode = 1
            route(:,1) = (/i, j/)

            DO WHILE (area_this < levsize)
               
               found = .false.

               ipxl = order2d(route(1,1),route(2,1))

               DO inode = 1, nnode

                  i = route(1,inode)
                  j = route(2,inode)

                  ipxl = min(order2d(i,j), ipxl) 

                  DO jnb = max(j-1,1), min(j+1,mp)
                     DO inb = max(i-1,1), min(i+1,np)
                        IF ((inb /= i) .or. (jnb /= j)) THEN
                           IF (hmask(inb,jnb) .and. (hunit(inb,jnb) == -1) .and. (order2d(inb,jnb) < ipxl)) THEN

                              inext = inb
                              jnext = jnb
                              DO WHILE (hunit(inext,jnext) == -1)
                                 ilast = inext
                                 jlast = jnext
                                 CALL nextij (ilast,jlast, dir(ilast,jlast),inext,jnext, is_local = .true.)
                              ENDDO

                              IF ((hunit(inext,jnext) == nunit) .or. (hunit(inext,jnext) == unit_next)) THEN
                                 found = .true.
                                 ipxl = min(order2d(inb,jnb), ipxl)
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDDO
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

                  nnode = 0
                  DO WHILE (hunit(i,j) == -1)
                     nnode = nnode + 1
                     route(:,nnode) = (/i, j/)
                     hunit(i,j) = nunit
                     area_this = area_this + area(i,j)
                        
                     CALL nextij (i,j, dir(i,j),inext,jnext, is_local = .true.)
                     i = inext
                     j = jnext
                  ENDDO

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
