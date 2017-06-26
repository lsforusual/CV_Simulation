subroutine I1I2 (i, N, m, i1, i2)
! Computes i1 and i2 for the stretch of m points along N, so that
! they centre i where possible, but make an asymmetric pair where not.
  implicit none
  integer :: i, N, m, i1, i2

  integer :: mm

  mm = (m+1) / 2
  i1 = MAX(0, i-mm+1);  i2 = i1 + m - 1  ! Bulk
  if (i2 > N) then;  i2 = N;  i1 = N - m + 1;  endif
end subroutine I1I2
