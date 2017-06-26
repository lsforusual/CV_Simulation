function U_TRAP (x, u, N)
! Integrates the set of function points u(0:N) at positions
! x(0:N), by simple trapezium integration. The x-intervals
! need not be equal.
  use STUFF;   implicit none
  integer        :: N
  real(kind=dbl) :: U_TRAP, x(0:N), u(0:N), integ

  integer        :: i

  integ = 0
  do i = 0, N-1
    integ = integ + (u(i)+u(i+1))/2 * (x(i+1)-x(i))
  enddo
  U_TRAP = integ
end function U_TRAP
