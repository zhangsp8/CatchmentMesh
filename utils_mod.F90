module utils_mod

   interface quicksort
      module procedure quicksort_single
      module procedure quicksort_double
      module procedure quicksort_int32
   end interface quicksort


contains

   !--------------------------------------------------
   subroutine insert_into_sorted_list1 (x, n, list, iloc, is_new_out)

      implicit none

      integer, intent(in) :: x
      integer, intent(inout) :: n
      integer, intent(inout) :: list(:)
      integer, intent(out)   :: iloc
      logical, intent(out), optional :: is_new_out

      logical :: is_new

      ! Local variables
      integer :: iright, ileft

      if (n == 0) then
         iloc = 1
         is_new = .true.
      elseif (x <= list(1)) then
         iloc = 1
         is_new = (x /= list(1))
      elseif (x > list(n)) then
         iloc = n + 1
         is_new = .true. 
      elseif (x == list(n)) then
         iloc = n
         is_new = .false. 
      else
         iright = 1
         ileft  = n

         do while (.true.)
            if (ileft - iright > 1) then
               iloc = (iright + ileft) / 2
               if (x > list(iloc)) then
                  iright = iloc
               elseif (x < list(iloc)) then
                  ileft = iloc
               else
                  is_new = .false.
                  exit 
               end if
            else
               iloc = ileft
               is_new = .true.
               exit
            end if
         end do
      end if

      if (is_new) then
         if (iloc <= n) then
            list(iloc+1:n+1) = list(iloc:n)
         end if

         list(iloc) = x
         n = n + 1
      end if

      if (present(is_new_out)) then
         is_new_out = is_new
      end if

   end subroutine insert_into_sorted_list1

   !--------------------------------------------------
   SUBROUTINE insert_into_sorted_list2 (x, y, n, xlist, ylist, iloc, is_new_out)

      IMPLICIT NONE

      INTEGER, intent(in) :: x, y
      INTEGER, intent(inout) :: n
      INTEGER, intent(inout) :: xlist(:), ylist(:)
      INTEGER, intent(out)   :: iloc
      LOGICAL, intent(out), optional :: is_new_out

      ! Local variables
      LOGICAL :: is_new
      INTEGER :: ileft, iright

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
                  exit
               ENDIF
            ELSE
               iloc = iright
               is_new = .true.
               exit
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

      INTEGER :: iloc

      INTEGER, intent(in) :: x
      INTEGER, intent(in) :: n
      INTEGER, intent(in) :: list (n)

      ! Local variables
      INTEGER :: i, ileft, iright

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
                     exit
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

      INTEGER :: iloc

      INTEGER, intent(in) :: x, y
      INTEGER, intent(in) :: n
      INTEGER, intent(in) :: xlist(:), ylist(:)

      ! Local variables
      INTEGER :: i, ileft, iright

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

         DO while (.true.)
            IF (iright - ileft > 1) THEN
               i = (ileft + iright) / 2
               IF ((y == ylist(i)) .and. (x == xlist(i))) THEN
                  iloc = i
                  exit
               ELSEIF ((y > ylist(i)) .or. ((y == ylist(i)) .and. (x > xlist(i)))) THEN
                  ileft = i
               ELSEIF ((y < ylist(i)) .or. ((y == ylist(i)) .and. (x < xlist(i)))) THEN
                  iright = i
               ENDIF
            ELSE
               iloc = 0
               exit
            ENDIF
         ENDDO
      ENDIF

   END FUNCTION find_in_sorted_list2
   
   !-----------------------------------------------------
   subroutine find_floor_lat (y, n, lat, iloc, yinc)

      use precision
      implicit none

      real(r8), intent(in) :: y
      integer,  intent(in) :: n
      real(r8), intent(in) :: lat (n)
      
      integer, intent(out) :: iloc
      integer, intent(out), optional :: yinc 

      ! Local variables
      integer :: i, iright, ileft

      if (lat(1) < lat(n))  then
         if (present(yinc))  yinc =  1
         if (y <= lat(1)) then
            iloc = 1
         elseif (y >= lat(n)) then
            iloc = n
         else
            ileft = 1;  iright = n

            do while (iright - ileft > 1)
               i = (iright + ileft) / 2
               if (y >= lat(i)) then
                  ileft = i
               else
                  iright = i
               end if
            end do

            iloc = ileft
         end if
      else
         if (present(yinc))  yinc = -1
         if (y >= lat(1)) then
            iloc = 1
         elseif (y <= lat(n)) then
            iloc = n
         else
            ileft = 1;  iright = n

            do while (iright - ileft > 1)
               i = (iright + ileft) / 2
               if (y >= lat(i)) then
                  iright = i
               else
                  ileft = i
               end if
            end do

            iloc = iright
         end if
      end if

   end subroutine find_floor_lat

   !-----------------------------------------------------
   subroutine find_ceiling_lat (y, n, lat, iloc, yinc)

      use precision
      implicit none

      real(r8), intent(in) :: y
      integer,  intent(in) :: n
      real(r8), intent(in) :: lat (n)
      
      integer, intent(out) :: iloc
      integer, intent(out), optional :: yinc 

      ! Local variables
      integer :: i, iright, ileft

      if (lat(1) < lat(n))  then
         if (present(yinc))  yinc = 1
         if (y <= lat(1)) then
            iloc = 1
         elseif (y >= lat(n)) then
            iloc = n
         else
            ileft = 1;  iright = n

            do while (iright - ileft > 1)
               i = (iright + ileft) / 2
               if (y > lat(i)) then
                  ileft = i
               else
                  iright = i
               end if
            end do

            iloc = ileft
         end if
      else
         if (present(yinc))  yinc = -1
         if (y >= lat(1)) then
            iloc = 1
         elseif (y <= lat(n)) then
            iloc = n
         else
            ileft = 1;  iright = n

            do while (iright - ileft > 1)
               i = (iright + ileft) / 2
               if (y > lat(i)) then
                  iright = i
               else
                  ileft = i
               end if
            end do

            iloc = iright
         end if
      end if

   end subroutine find_ceiling_lat

   !-----------------------------------------------------
   subroutine find_floor_lon (x, n, lon, iloc)

      use precision
      implicit none

      real(r8), intent(in) :: x
      integer,  intent(in) :: n
      real(r8), intent(in) :: lon (n)
      
      integer, intent(out) :: iloc

      ! Local variables
      integer :: i, iright, ileft

      if (lon(1) >= lon(n)) then  
         ! for the cast 0~180 -> -180~0
         if ((x >= lon(n)) .and. (x < lon(1))) then
            iloc = n; return
         end if
      else   
         ! for the case -180~180
         if ((x < lon(1)) .or. (x >= lon(n))) then
            iloc = n; return
         end if 
      end if

      ileft = 1; iright = n
      do while (iright - ileft > 1)
         i = (iright + ileft)/2
         if (lon(iright) > lon(i)) then
            if (x >= lon(i)) then
               ileft = i
            else
               iright = i
            end if
         else
            if ((x < lon(iright)) .or. (x >= lon(i))) then
               ileft = i
            else
               iright = i
            end if 
         end if
      end do
      
      iloc = ileft

   end subroutine find_floor_lon

   !-----------------------------------------
   function max_lon (lon1, lon2) result (lon)
      use precision
      implicit none

      real(r8) :: lon
      real(r8), intent(in) :: lon1, lon2

      real(r8) :: lonr, lonl

      lonr = max(lon1,lon2);   lonl = min(lon1,lon2)
      if (lonr - lonl > 180.0) then
         lon = lonl
      else
         lon = lonr
      end if

   end function max_lon

   !-----------------------------------------
   function min_lon (lon1, lon2) result (lon)
      use precision
      implicit none

      real(r8) :: lon
      real(r8), intent(in) :: lon1, lon2

      real(r8) :: lonr, lonl

      lonr = max(lon1,lon2);   lonl = min(lon1,lon2)
      if (lonr - lonl > 180.0) then
         lon = lonr
      else
         lon = lonl
      end if

   end function min_lon

   !-----------------------------------------------------
   subroutine find_ceiling_lon (x, n, lon, iloc)

      use precision
      implicit none

      real(r8), intent(in) :: x
      integer,  intent(in) :: n
      real(r8), intent(in) :: lon (n)
      
      integer, intent(out) :: iloc

      ! Local variables
      integer :: i, iright, ileft

      if (lon(1) >= lon(n)) then  
         ! for the cast 0~180 -> -180~0
         if ((x > lon(n)) .and. (x <= lon(1))) then
            iloc = 1; return
         end if
      else   
         ! for the case -180~180
         if ((x <= lon(1)) .or. (x > lon(n))) then
            iloc = 1; return
         end if 
      end if

      ileft = 1; iright = n
      do while (iright - ileft > 1)
         i = (iright + ileft)/2
         if (lon(iright) > lon(i)) then
            if (x > lon(i)) then
               ileft = i
            else
               iright = i
            end if
         else
            if ((x <= lon(iright)) .or. (x > lon(i))) then
               ileft = i
            else
               iright = i
            end if 
         end if
      end do
      
      iloc = iright

   end subroutine find_ceiling_lon


   !-----------------------------------------------------
   recursive subroutine quicksort_single (nA, A, order)

      use precision
      implicit none

      integer, intent(in) :: nA
      real(4), intent(inout) :: A     (nA)
      integer, intent(inout) :: order (nA)

      ! Local variables
      integer :: left, right
      real(4) :: pivot
      integer :: marker
      real(4) :: rtemp 
      integer :: itemp

      if (nA > 1) then

         pivot = A (nA/2)
         left  = 0
         right = nA + 1

         do while (left < right)
            right = right - 1
            do while (A(right) > pivot)
               right = right - 1
            end do

            left = left + 1
            do while (A(left) < pivot)
               left = left + 1
            end do

            if (left < right) then
               rtemp    = A(left)
               A(left)  = A(right)
               A(right) = rtemp

               itemp        = order(left)
               order(left)  = order(right)
               order(right) = itemp

            end if
         end do

         if (left == right) then
            marker = left + 1
         else
            marker = left
         end if

         call quicksort_single (marker-1,    A(:marker-1), order(:marker-1))
         call quicksort_single (nA-marker+1, A(marker:),   order(marker:)  )

      end if

   end subroutine quicksort_single

   !-----------------------------------------------------
   recursive subroutine quicksort_double (nA, A, order)

      use precision
      implicit none

      integer , intent(in) :: nA
      real(r8), intent(inout) :: A     (nA)
      integer , intent(inout) :: order (nA)

      ! Local variables
      integer  :: left, right
      real(r8) :: pivot
      integer  :: marker
      real(r8) :: rtemp 
      integer  :: itemp

      if (nA > 1) then

         pivot = A (nA/2)
         left  = 0
         right = nA + 1

         do while (left < right)
            right = right - 1
            do while (A(right) > pivot)
               right = right - 1
            end do

            left = left + 1
            do while (A(left) < pivot)
               left = left + 1
            end do

            if (left < right) then
               rtemp    = A(left)
               A(left)  = A(right)
               A(right) = rtemp

               itemp        = order(left)
               order(left)  = order(right)
               order(right) = itemp

            end if
         end do

         if (left == right) then
            marker = left + 1
         else
            marker = left
         end if

         call quicksort_double (marker-1,    A(:marker-1), order(:marker-1))
         call quicksort_double (nA-marker+1, A(marker:),   order(marker:)  )

      end if

   end subroutine quicksort_double

   !-----------------------------------------------------
   recursive subroutine quicksort_int32 (nA, A, order)

      use precision
      implicit none

      integer, intent(in) :: nA
      integer, intent(inout) :: A     (nA)
      integer, intent(inout) :: order (nA)

      ! Local variables
      integer  :: left, right
      real(r8) :: pivot
      integer  :: marker
      integer  :: itemp

      if (nA > 1) then

         pivot = A (nA/2)
         left  = 0
         right = nA + 1

         do while (left < right)
            right = right - 1
            do while (A(right) > pivot)
               right = right - 1
            end do

            left = left + 1
            do while (A(left) < pivot)
               left = left + 1
            end do

            if (left < right) then
               itemp    = A(left)
               A(left)  = A(right)
               A(right) = itemp

               itemp        = order(left)
               order(left)  = order(right)
               order(right) = itemp

            end if
         end do

         if (left == right) then
            marker = left + 1
         else
            marker = left
         end if

         call quicksort_int32 (marker-1,    A(:marker-1), order(:marker-1))
         call quicksort_int32 (nA-marker+1, A(marker:),   order(marker:)  )

      end if

   end subroutine quicksort_int32

   !-----------------------------------------------------
   recursive function quickselect (nA, A, k) result(selected)

      use precision
      implicit none

      real(r8) :: selected

      integer , intent(in)    :: nA
      real(r8), intent(inout) :: A (nA)
      integer,  intent(in)    :: k

      ! Local variables
      integer  :: left, right
      real(r8) :: pivot
      integer  :: marker
      real(r8) :: rtemp 

      if (nA > 1) then

         pivot = A (nA/2)
         left  = 0
         right = nA + 1

         do while (left < right)
            right = right - 1
            do while (A(right) > pivot)
               right = right - 1
            end do

            left = left + 1
            do while (A(left) < pivot)
               left = left + 1
            end do

            if (left < right) then
               rtemp    = A(left)
               A(left)  = A(right)
               A(right) = rtemp
            end if
         end do

         if (left == right) then
            marker = left + 1
         else
            marker = left
         end if

         if (k < marker) then
            selected = quickselect (marker-1, A(:marker-1), k)
         else
            selected = quickselect (nA-marker+1, A(marker:), k-marker+1)
         end if

      else
         selected = A(1) 
      end if

   end function quickselect
   
   
   ! ------------------------
   function median(x, n) result(mval)

      use precision
      implicit none

      real(r8) :: mval

      integer,  intent(in) :: n
      real(r8), intent(in) :: x(n)

      real(r8), allocatable :: xtemp(:)

      allocate (xtemp(n))
      xtemp = x
      
      if (mod(n,2) == 0) then
         mval = (quickselect(n,xtemp,n/2) + quickselect(n,xtemp,n/2+1)) / 2.0_r8
      else
         mval = quickselect(n,xtemp,n/2+1)
      end if

      deallocate (xtemp)

   end function median

   
   !-----------------------------------------------------
   function areaquad (lat0, lat1, lon0, lon1) result(area)

      use precision
      use MathConstants, only : deg2rad
      implicit none

      real(r8), parameter :: re = 6.37122e3 ! kilometer 
      real(r8) :: area

      real(r8), intent(in) :: lat0, lat1, lon0, lon1

      real(r8) :: dx, dy

      if (lon1 < lon0) then
         dx = (lon1 + 360 -lon0) * deg2rad
      else
         dx = (lon1 - lon0) * deg2rad
      end if

      dy = sin(lat1 * deg2rad) - sin(lat0 * deg2rad)

      area = dx * dy * re * re

   end function areaquad

end module utils_mod
