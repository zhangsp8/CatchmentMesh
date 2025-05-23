MODULE precision
!-------------------------------------------------------------------------------
! Purpose:
!	Define the precision to USE for floating point and integer operations
!	throughout the model.
!-------------------------------------------------------------------------------
  integer, parameter :: r4 = selected_real_kind(5)
  integer, parameter :: r8 = selected_real_kind(12)
  integer, parameter :: i8 = selected_int_kind(13)
END MODULE precision
