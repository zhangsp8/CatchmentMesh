MODULE cvt_mod

CONTAINS

   SUBROUTINE cvt (areaid, point_num, points, cell_num, cell_id, weight)

      !*****************************************************************************80
      !
      !  CVT computes a Centroidal Voronoi Tessellation.
      !
      !  Reference:
      !    Qiang Du, Vance Faber, Max Gunzburger,
      !    Centroidal Voronoi Tessellations: Applications and Algorithms,
      !    SIAM Review, Volume 41, 1999, pages 637-676.
      !
      !  Original Author: Shupeng Zhang, June 2023.
      !
      !*****************************************************************************80

      implicit none

      integer, parameter :: rk = 8

      integer areaid
      integer point_num
      integer points (2,point_num) ! Input:  grid points in the area
      integer cell_num             ! Input:  number of cells
      integer cell_id(point_num)   ! Output: cell id where the grid point is located
      
      real*4, optional :: weight (point_num)   ! Input: weight

      ! Local Variables
      real ( kind = rk ), allocatable :: dist (:)

      real ( kind = rk ), allocatable :: sample (:,:)
      real ( kind = rk ), allocatable :: r1 (:,:)
      real ( kind = rk ), allocatable :: r2 (:,:)
      
      real ( kind = rk ), allocatable :: xdst2 (:,:)

      real ( kind = rk ), allocatable :: wts   (:)
      real ( kind = rk ), allocatable :: sumwt (:)

      integer,       parameter :: it_max = 100
      real(kind=rk), parameter :: tol = 1.e-2

      real ( kind = rk ) it_diff, d0, d04, dd
      integer del, it_num, js, jr, icell, jcell, ip, i, j

      allocate (r1 (2, cell_num))

      allocate (sample (2, point_num))
      sample = real(points, kind = rk)

      ! del = point_num / cell_num
      ! DO icell = 1, cell_num
      !    r1(:,icell)  = sample(:,del*(icell-1)+1)
      ! ENDDO

      allocate (dist (point_num))
      CALL random_number(dist)
      DO jr = 1, cell_num
         i = maxloc(dist, dim=1)
         r1(:,jr) = sample(:,i)
         IF (jr == 1) THEN
            DO js = 1, point_num
               dist(js) = (sample(1,js)-r1(1,1))**2 + (sample(2,js)-r1(2,1))**2
            ENDDO
         ELSE
            DO js = 1, point_num
               dist(js) = min((sample(1,js)-r1(1,jr))**2 + (sample(2,js)-r1(2,jr))**2, dist(js))
            ENDDO
         ENDIF
      ENDDO
      deallocate (dist)

      allocate (wts (point_num))
      IF (present(weight)) THEN
         wts = weight
      ELSE
         wts(:) = 1.0
      ENDIF

      allocate (r2 (2, cell_num))
      allocate (xdst2 (cell_num,cell_num))
      allocate (sumwt (cell_num))
      
      DO i = 1, cell_num
         DO j = i+1, cell_num
            xdst2(i,j) = (r1(1,i) - r1(1,j))**2 + (r1(2,i) - r1(2,j))**2
            xdst2(j,i) = xdst2(i,j)
         ENDDO
      ENDDO

      IF (point_num > 10000) THEN 
!$OMP PARALLEL DO NUM_THREADS(10) PRIVATE(js,d0,d04,jr,dd)
         do js = 1, point_num
            cell_id(js) = 1
            d0 = (r1(1,cell_id(js)) - sample(1,js))**2 + (r1(2,cell_id(js)) - sample(2,js))**2
            d04 = d0 * 4.0
            DO jr = 1, cell_num
               IF (jr /= cell_id(js)) THEN
                  IF (d04 > xdst2(cell_id(js),jr)) THEN
                     dd = (r1(1,jr) - sample(1,js))**2 + (r1(2,jr) - sample(2,js))**2
                     IF (dd < d0) THEN
                        cell_id(js) = jr
                        d0 = dd
                        d04 = d0 * 4.0
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ELSE
         do js = 1, point_num
            cell_id(js) = 1
            d0 = (r1(1,cell_id(js)) - sample(1,js))**2 + (r1(2,cell_id(js)) - sample(2,js))**2
            d04 = d0 * 4.0
            DO jr = 1, cell_num
               IF (jr /= cell_id(js)) THEN
                  IF (d04 > xdst2(cell_id(js),jr)) THEN
                     dd = (r1(1,jr) - sample(1,js))**2 + (r1(2,jr) - sample(2,js))**2
                     IF (dd < d0) THEN
                        cell_id(js) = jr
                        d0 = dd
                        d04 = d0 * 4.0
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      it_num = 0
      !  Carry out the iteration.
      do while ( it_num < it_max )

         it_num = it_num + 1

         !  Estimate the new centroids.
         r2(1:2,1:cell_num) = 0.0D+00
         sumwt (1:cell_num) = 0.0D+00
         do js = 1, point_num
            r2(1:2,cell_id(js)) = r2(1:2,cell_id(js)) + sample(1:2,js) * wts(js)
            sumwt(cell_id(js)) = sumwt(cell_id(js)) + wts(js)
         end do

         DO WHILE (any(sumwt == 0))
            icell = findloc(sumwt, 0, dim=1)
            jcell = maxloc(sumwt, dim=1)

            ip = findloc(cell_id, jcell, dim=1)
            cell_id(ip) = icell

            r2(1:2,icell) = r2(1:2,icell) + sample(1:2,ip) * wts(ip)
            r2(1:2,jcell) = r2(1:2,jcell) - sample(1:2,ip) * wts(ip)
            sumwt(icell) = wts(ip)
            sumwt(jcell) = sumwt(jcell) - wts(ip)
         ENDDO

         do jr = 1, cell_num
            r2(1:2,jr) = r2(1:2,jr) / real ( sumwt(jr), kind = rk )
         end do

         !  Determine the sum of the distances between the old generators 
         !  and the estimated centroids.
         it_diff = sqrt(sum((r2(1:2,1)-r1(1:2,1))**2))
         do jr = 2, cell_num
            it_diff = max(it_diff, sqrt(sum((r2(1:2,jr)-r1(1:2,jr))**2)))
         end do
      
         !  Replace the generators by the centroids.
         r1 = r2
         DO i = 1, cell_num
            DO j = i+1, cell_num
               xdst2(i,j) = (r1(1,i) - r1(1,j))**2 + (r1(2,i) - r1(2,j))**2
               xdst2(j,i) = xdst2(i,j)
            ENDDO
         ENDDO

         ! Update clusters
         IF (point_num > 10000) THEN 
!$OMP PARALLEL DO NUM_THREADS(10) PRIVATE(js,d0,d04,jr,dd)
            do js = 1, point_num
               d0 = (r1(1,cell_id(js)) - sample(1,js))**2 + (r1(2,cell_id(js)) - sample(2,js))**2
               d04 = d0 * 4.0
               DO jr = 1, cell_num
                  IF (jr /= cell_id(js)) THEN
                     IF (d04 > xdst2(cell_id(js),jr)) THEN
                        dd = (r1(1,jr) - sample(1,js))**2 + (r1(2,jr) - sample(2,js))**2
                        IF (dd < d0) THEN
                           cell_id(js) = jr
                           d0 = dd
                           d04 = d0 * 4.0
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ELSE
            do js = 1, point_num
               d0 = (r1(1,cell_id(js)) - sample(1,js))**2 + (r1(2,cell_id(js)) - sample(2,js))**2
               d04 = d0 * 4.0
               DO jr = 1, cell_num
                  IF (jr /= cell_id(js)) THEN
                     IF (d04 > xdst2(cell_id(js),jr)) THEN
                        dd = (r1(1,jr) - sample(1,js))**2 + (r1(2,jr) - sample(2,js))**2
                        IF (dd < d0) THEN
                           cell_id(js) = jr
                           d0 = dd
                           d04 = d0 * 4.0
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
         
         write(*,'(A,I0,A,I0,A,I0,A,I0,A,ES15.4)') 'CVT: ID ', areaid, &
            ', Points ', point_num, ', Cells ', cell_num, ' Iter ', it_num, ', Diff ', it_diff

         IF (it_diff < tol) EXIT

      end do

      deallocate (sample)
      deallocate (r1)
      deallocate (r2)
      deallocate (sumwt)
      deallocate (wts)
      deallocate (xdst2)

   END SUBROUTINE

END MODULE
