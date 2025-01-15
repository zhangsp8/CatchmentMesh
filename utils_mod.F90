MODULE utils_mod

   interface quicksort
      MODULE procedure quicksort_single
      MODULE procedure quicksort_double
      MODULE procedure quicksort_int32
   END interface quicksort


CONTAINS

   !--------------------------------------------------
   SUBROUTINE insert_into_sorted_list1 (x, n, list, iloc, is_new_out)

      IMPLICIT NONE

      integer, intent(in) :: x
      integer, intent(inout) :: n
      integer, intent(inout) :: list(:)
      integer, intent(out)   :: iloc
      logical, intent(out), optional :: is_new_out

      logical :: is_new

      ! Local variables
      integer :: iright, ileft

      IF (n == 0) THEN
         iloc = 1
         is_new = .true.
      ELSEIF (x <= list(1)) THEN
         iloc = 1
         is_new = (x /= list(1))
      ELSEIF (x > list(n)) THEN
         iloc = n + 1
         is_new = .true. 
      ELSEIF (x == list(n)) THEN
         iloc = n
         is_new = .false. 
      ELSE
         iright = 1
         ileft  = n

         DO WHILE (.true.)
            IF (ileft - iright > 1) THEN
               iloc = (iright + ileft) / 2
               IF (x > list(iloc)) THEN
                  iright = iloc
               ELSEIF (x < list(iloc)) THEN
                  ileft = iloc
               ELSE
                  is_new = .false.
                  EXIT 
               ENDIF
            ELSE
               iloc = ileft
               is_new = .true.
               EXIT
            ENDIF
         ENDDO
      ENDIF

      IF (is_new) THEN
         IF (iloc <= n) THEN
            list(iloc+1:n+1) = list(iloc:n)
         ENDIF

         list(iloc) = x
         n = n + 1
      ENDIF

      IF (present(is_new_out)) THEN
         is_new_out = is_new
      ENDIF

   END SUBROUTINE insert_into_sorted_list1

   !--------------------------------------------------
   SUBROUTINE insert_into_sorted_list2 (x, y, n, xlist, ylist, iloc, is_new_out)

      IMPLICIT NONE

      integer, intent(in) :: x, y
      integer, intent(inout) :: n
      integer, intent(inout) :: xlist(:), ylist(:)
      integer, intent(out)   :: iloc
      logical, intent(out), optional :: is_new_out

      ! Local variables
      logical :: is_new
      integer :: ileft, iright

      IF (n == 0) THEN
         iloc = 1
         is_new = .true.
      ELSEIF ((y < ylist(1)) .or. ((y == ylist(1)) .and. (x <= xlist(1)))) THEN
         iloc = 1
         is_new = (x /= xlist(1)) .or. (y /= ylist(1))
      ELSEIF ((y > ylist(n)) .or. ((y == ylist(n)) .and. (x > xlist(n)))) THEN
         iloc = n + 1
         is_new = .true.
      ELSEIF ((x == xlist(n)) .and. (y == ylist(n))) THEN
         iloc = n
         is_new = .false.
      ELSE
         ileft  = 1
         iright = n

         DO WHILE (.true.)
            IF (iright - ileft > 1) THEN
               iloc = (ileft + iright) / 2
               IF ((y > ylist(iloc)) .or. ((y == ylist(iloc)) .and. (x > xlist(iloc)))) THEN
                  ileft = iloc
               ELSEIF ((y < ylist(iloc)) .or. ((y == ylist(iloc)) .and. (x < xlist(iloc)))) THEN
                  iright = iloc
               ELSE
                  is_new = .false.
                  EXIT
               ENDIF
            ELSE
               iloc = iright
               is_new = .true.
               EXIT
            ENDIF
         ENDDO
      ENDIF

      IF (is_new) THEN
         IF (iloc <= n) THEN
            xlist(iloc+1:n+1) = xlist(iloc:n)
            ylist(iloc+1:n+1) = ylist(iloc:n)
         ENDIF

         xlist(iloc) = x
         ylist(iloc) = y
         n = n + 1
      ENDIF

      IF (present(is_new_out)) THEN
         is_new_out = is_new
      ENDIF

   END SUBROUTINE insert_into_sorted_list2

   !--------------------------------------------------
   FUNCTION find_in_sorted_list1 (x, n, list) result(iloc)

      IMPLICIT NONE

      integer :: iloc

      integer, intent(in) :: x
      integer, intent(in) :: n
      integer, intent(in) :: list (n)

      ! Local variables
      integer :: i, ileft, iright

      iloc = 0
      IF (n > 0) THEN
         IF ((x >= list(1)) .and. (x <= list(n))) THEN
            IF (x == list(1)) THEN
               iloc = 1
            ELSEIF (x == list(n)) THEN
               iloc = n
            ELSE
               ileft  = 1
               iright = n

               DO WHILE (iright - ileft > 1)
                  i = (ileft + iright) / 2
                  IF (x == list(i)) THEN
                     iloc = i
                     EXIT
                  ELSEIF (x > list(i)) THEN
                     ileft = i
                  ELSEIF (x < list(i)) THEN
                     iright = i
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDIF

   END FUNCTION find_in_sorted_list1

   !--------------------------------------------------
   FUNCTION find_in_sorted_list2 (x, y, n, xlist, ylist) result(iloc)

      IMPLICIT NONE

      integer :: iloc

      integer, intent(in) :: x, y
      integer, intent(in) :: n
      integer, intent(in) :: xlist(:), ylist(:)

      ! Local variables
      integer :: i, ileft, iright

      iloc = 0
      IF (n < 1) RETURN

      IF ((y < ylist(1)) .or. ((y == ylist(1)) .and. (x < xlist(1)))) THEN
         iloc = 0
      ELSEIF ((y > ylist(n)) .or. ((y == ylist(n)) .and. (x > xlist(n)))) THEN
         iloc = 0
      ELSEIF ((x == xlist(1)) .and. (y == ylist(1))) THEN
         iloc = 1
      ELSEIF ((x == xlist(n)) .and. (y == ylist(n))) THEN
         iloc = n
      ELSE
         ileft  = 1
         iright = n

         DO WHILE (.true.)
            IF (iright - ileft > 1) THEN
               i = (ileft + iright) / 2
               IF ((y == ylist(i)) .and. (x == xlist(i))) THEN
                  iloc = i
                  EXIT
               ELSEIF ((y > ylist(i)) .or. ((y == ylist(i)) .and. (x > xlist(i)))) THEN
                  ileft = i
               ELSEIF ((y < ylist(i)) .or. ((y == ylist(i)) .and. (x < xlist(i)))) THEN
                  iright = i
               ENDIF
            ELSE
               iloc = 0
               EXIT
            ENDIF
         ENDDO
      ENDIF

   END FUNCTION find_in_sorted_list2
   
   !-----------------------------------------------------
   SUBROUTINE find_floor_lat (y, n, lat, iloc, yinc)

      USE precision
      IMPLICIT NONE

      real(r8), intent(in) :: y
      integer,  intent(in) :: n
      real(r8), intent(in) :: lat (n)
      
      integer, intent(out) :: iloc
      integer, intent(out), optional :: yinc 

      ! Local variables
      integer :: i, iright, ileft

      IF (lat(1) < lat(n))  THEN
         IF (present(yinc))  yinc =  1
         IF (y <= lat(1)) THEN
            iloc = 1
         ELSEIF (y >= lat(n)) THEN
            iloc = n
         ELSE
            ileft = 1;  iright = n

            DO WHILE (iright - ileft > 1)
               i = (iright + ileft) / 2
               IF (y >= lat(i)) THEN
                  ileft = i
               ELSE
                  iright = i
               ENDIF
            ENDDO

            iloc = ileft
         ENDIF
      ELSE
         IF (present(yinc))  yinc = -1
         IF (y >= lat(1)) THEN
            iloc = 1
         ELSEIF (y <= lat(n)) THEN
            iloc = n
         ELSE
            ileft = 1;  iright = n

            DO WHILE (iright - ileft > 1)
               i = (iright + ileft) / 2
               IF (y >= lat(i)) THEN
                  iright = i
               ELSE
                  ileft = i
               ENDIF
            ENDDO

            iloc = iright
         ENDIF
      ENDIF

   END SUBROUTINE find_floor_lat

   !-----------------------------------------------------
   SUBROUTINE find_ceiling_lat (y, n, lat, iloc, yinc)

      USE precision
      IMPLICIT NONE

      real(r8), intent(in) :: y
      integer,  intent(in) :: n
      real(r8), intent(in) :: lat (n)
      
      integer, intent(out) :: iloc
      integer, intent(out), optional :: yinc 

      ! Local variables
      integer :: i, iright, ileft

      IF (lat(1) < lat(n))  THEN
         IF (present(yinc))  yinc = 1
         IF (y <= lat(1)) THEN
            iloc = 1
         ELSEIF (y >= lat(n)) THEN
            iloc = n
         ELSE
            ileft = 1;  iright = n

            DO WHILE (iright - ileft > 1)
               i = (iright + ileft) / 2
               IF (y > lat(i)) THEN
                  ileft = i
               ELSE
                  iright = i
               ENDIF
            ENDDO

            iloc = ileft
         ENDIF
      ELSE
         IF (present(yinc))  yinc = -1
         IF (y >= lat(1)) THEN
            iloc = 1
         ELSEIF (y <= lat(n)) THEN
            iloc = n
         ELSE
            ileft = 1;  iright = n

            DO WHILE (iright - ileft > 1)
               i = (iright + ileft) / 2
               IF (y > lat(i)) THEN
                  iright = i
               ELSE
                  ileft = i
               ENDIF
            ENDDO

            iloc = iright
         ENDIF
      ENDIF

   END SUBROUTINE find_ceiling_lat

   !-----------------------------------------------------
   SUBROUTINE find_floor_lon (x, n, lon, iloc)

      USE precision
      IMPLICIT NONE

      real(r8), intent(in) :: x
      integer,  intent(in) :: n
      real(r8), intent(in) :: lon (n)
      
      integer, intent(out) :: iloc

      ! Local variables
      integer :: i, iright, ileft

      IF (lon(1) >= lon(n)) THEN  
         ! for the cast 0~180 -> -180~0
         IF ((x >= lon(n)) .and. (x < lon(1))) THEN
            iloc = n; RETURN
         ENDIF
      ELSE   
         ! for the CASE -180~180
         IF ((x < lon(1)) .or. (x >= lon(n))) THEN
            iloc = n; RETURN
         ENDIF 
      ENDIF

      ileft = 1; iright = n
      DO WHILE (iright - ileft > 1)
         i = (iright + ileft)/2
         IF (lon(iright) > lon(i)) THEN
            IF (x >= lon(i)) THEN
               ileft = i
            ELSE
               iright = i
            ENDIF
         ELSE
            IF ((x < lon(iright)) .or. (x >= lon(i))) THEN
               ileft = i
            ELSE
               iright = i
            ENDIF 
         ENDIF
      ENDDO
      
      iloc = ileft

   END SUBROUTINE find_floor_lon

   !-----------------------------------------
   FUNCTION max_lon (lon1, lon2) result (lon)
      USE precision
      IMPLICIT NONE

      real(r8) :: lon
      real(r8), intent(in) :: lon1, lon2

      real(r8) :: lonr, lonl

      lonr = max(lon1,lon2);   lonl = min(lon1,lon2)
      IF (lonr - lonl > 180.0) THEN
         lon = lonl
      ELSE
         lon = lonr
      ENDIF

   END FUNCTION max_lon

   !-----------------------------------------
   FUNCTION min_lon (lon1, lon2) result (lon)
      USE precision
      IMPLICIT NONE

      real(r8) :: lon
      real(r8), intent(in) :: lon1, lon2

      real(r8) :: lonr, lonl

      lonr = max(lon1,lon2);   lonl = min(lon1,lon2)
      IF (lonr - lonl > 180.0) THEN
         lon = lonr
      ELSE
         lon = lonl
      ENDIF

   END FUNCTION min_lon

   !-----------------------------------------------------
   SUBROUTINE find_ceiling_lon (x, n, lon, iloc)

      USE precision
      IMPLICIT NONE

      real(r8), intent(in) :: x
      integer,  intent(in) :: n
      real(r8), intent(in) :: lon (n)
      
      integer, intent(out) :: iloc

      ! Local variables
      integer :: i, iright, ileft

      IF (lon(1) >= lon(n)) THEN  
         ! for the cast 0~180 -> -180~0
         IF ((x > lon(n)) .and. (x <= lon(1))) THEN
            iloc = 1; RETURN
         ENDIF
      ELSE   
         ! for the CASE -180~180
         IF ((x <= lon(1)) .or. (x > lon(n))) THEN
            iloc = 1; RETURN
         ENDIF 
      ENDIF

      ileft = 1; iright = n
      DO WHILE (iright - ileft > 1)
         i = (iright + ileft)/2
         IF (lon(iright) > lon(i)) THEN
            IF (x > lon(i)) THEN
               ileft = i
            ELSE
               iright = i
            ENDIF
         ELSE
            IF ((x <= lon(iright)) .or. (x > lon(i))) THEN
               ileft = i
            ELSE
               iright = i
            ENDIF 
         ENDIF
      ENDDO
      
      iloc = iright

   END SUBROUTINE find_ceiling_lon


   !-----------------------------------------------------
   recursive SUBROUTINE quicksort_single (nA, A, order)

      USE precision
      IMPLICIT NONE

      integer, intent(in) :: nA
      real(4), intent(inout) :: A     (nA)
      integer, intent(inout) :: order (nA)

      ! Local variables
      integer :: left, right
      real(4) :: pivot
      integer :: marker
      real(4) :: rtemp 
      integer :: itemp

      IF (nA > 1) THEN

         pivot = A (nA/2)
         left  = 0
         right = nA + 1

         DO WHILE (left < right)
            right = right - 1
            DO WHILE (A(right) > pivot)
               right = right - 1
            ENDDO

            left = left + 1
            DO WHILE (A(left) < pivot)
               left = left + 1
            ENDDO

            IF (left < right) THEN
               rtemp    = A(left)
               A(left)  = A(right)
               A(right) = rtemp

               itemp        = order(left)
               order(left)  = order(right)
               order(right) = itemp

            ENDIF
         ENDDO

         IF (left == right) THEN
            marker = left + 1
         ELSE
            marker = left
         ENDIF

         CALL quicksort_single (marker-1,    A(:marker-1), order(:marker-1))
         CALL quicksort_single (nA-marker+1, A(marker:),   order(marker:)  )

      ENDIF

   END SUBROUTINE quicksort_single

   !-----------------------------------------------------
   recursive SUBROUTINE quicksort_double (nA, A, order)

      USE precision
      IMPLICIT NONE

      integer , intent(in) :: nA
      real(r8), intent(inout) :: A     (nA)
      integer , intent(inout) :: order (nA)

      ! Local variables
      integer  :: left, right
      real(r8) :: pivot
      integer  :: marker
      real(r8) :: rtemp 
      integer  :: itemp

      IF (nA > 1) THEN

         pivot = A (nA/2)
         left  = 0
         right = nA + 1

         DO WHILE (left < right)
            right = right - 1
            DO WHILE (A(right) > pivot)
               right = right - 1
            ENDDO

            left = left + 1
            DO WHILE (A(left) < pivot)
               left = left + 1
            ENDDO

            IF (left < right) THEN
               rtemp    = A(left)
               A(left)  = A(right)
               A(right) = rtemp

               itemp        = order(left)
               order(left)  = order(right)
               order(right) = itemp

            ENDIF
         ENDDO

         IF (left == right) THEN
            marker = left + 1
         ELSE
            marker = left
         ENDIF

         CALL quicksort_double (marker-1,    A(:marker-1), order(:marker-1))
         CALL quicksort_double (nA-marker+1, A(marker:),   order(marker:)  )

      ENDIF

   END SUBROUTINE quicksort_double

   !-----------------------------------------------------
   recursive SUBROUTINE quicksort_int32 (nA, A, order)

      USE precision
      IMPLICIT NONE

      integer, intent(in) :: nA
      integer, intent(inout) :: A     (nA)
      integer, intent(inout) :: order (nA)

      ! Local variables
      integer  :: left, right
      real(r8) :: pivot
      integer  :: marker
      integer  :: itemp

      IF (nA > 1) THEN

         pivot = A (nA/2)
         left  = 0
         right = nA + 1

         DO WHILE (left < right)
            right = right - 1
            DO WHILE (A(right) > pivot)
               right = right - 1
            ENDDO

            left = left + 1
            DO WHILE (A(left) < pivot)
               left = left + 1
            ENDDO

            IF (left < right) THEN
               itemp    = A(left)
               A(left)  = A(right)
               A(right) = itemp

               itemp        = order(left)
               order(left)  = order(right)
               order(right) = itemp

            ENDIF
         ENDDO

         IF (left == right) THEN
            marker = left + 1
         ELSE
            marker = left
         ENDIF

         CALL quicksort_int32 (marker-1,    A(:marker-1), order(:marker-1))
         CALL quicksort_int32 (nA-marker+1, A(marker:),   order(marker:)  )

      ENDIF

   END SUBROUTINE quicksort_int32

   !-----------------------------------------------------
   recursive FUNCTION quickselect (nA, A, k) result(selected)

      USE precision
      IMPLICIT NONE

      real(r8) :: selected

      integer , intent(in)    :: nA
      real(r8), intent(inout) :: A (nA)
      integer,  intent(in)    :: k

      ! Local variables
      integer  :: left, right
      real(r8) :: pivot
      integer  :: marker
      real(r8) :: rtemp 

      IF (nA > 1) THEN

         pivot = A (nA/2)
         left  = 0
         right = nA + 1

         DO WHILE (left < right)
            right = right - 1
            DO WHILE (A(right) > pivot)
               right = right - 1
            ENDDO

            left = left + 1
            DO WHILE (A(left) < pivot)
               left = left + 1
            ENDDO

            IF (left < right) THEN
               rtemp    = A(left)
               A(left)  = A(right)
               A(right) = rtemp
            ENDIF
         ENDDO

         IF (left == right) THEN
            marker = left + 1
         ELSE
            marker = left
         ENDIF

         IF (k < marker) THEN
            selected = quickselect (marker-1, A(:marker-1), k)
         ELSE
            selected = quickselect (nA-marker+1, A(marker:), k-marker+1)
         ENDIF

      ELSE
         selected = A(1) 
      ENDIF

   END FUNCTION quickselect
   
   
   ! ------------------------
   FUNCTION median(x, n) result(mval)

      USE precision
      IMPLICIT NONE

      real(r8) :: mval

      integer,  intent(in) :: n
      real(r8), intent(in) :: x(n)

      real(r8), allocatable :: xtemp(:)

      allocate (xtemp(n))
      xtemp = x
      
      IF (mod(n,2) == 0) THEN
         mval = (quickselect(n,xtemp,n/2) + quickselect(n,xtemp,n/2+1)) / 2.0_r8
      ELSE
         mval = quickselect(n,xtemp,n/2+1)
      ENDIF

      deallocate (xtemp)

   END FUNCTION median

   
   !-----------------------------------------------------
   FUNCTION areaquad (lat0, lat1, lon0, lon1) result(area)

      USE precision
      USE MathConstants, only : deg2rad
      IMPLICIT NONE

      real(r8), parameter :: re = 6.37122e3 ! kilometer 
      real(r8) :: area

      real(r8), intent(in) :: lat0, lat1, lon0, lon1

      real(r8) :: dx, dy

      IF (lon1 < lon0) THEN
         dx = (lon1 + 360 -lon0) * deg2rad
      ELSE
         dx = (lon1 - lon0) * deg2rad
      ENDIF

      dy = sin(lat1 * deg2rad) - sin(lat0 * deg2rad)

      area = dx * dy * re * re

   END FUNCTION areaquad

END MODULE utils_mod
