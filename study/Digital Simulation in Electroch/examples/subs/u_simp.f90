function U_SIMP (x, u, N)
! Integrates the set of function points u(0:N) at positions
! x(0:N). The integration is done by fitting a parabola to
! all sets of three consecutive points and integrating the
! two panels. For all but the very first and last panels,
! there are always two estimates; these are averaged.
  use STUFF;   implicit none
  integer        :: N
  real(kind=dbl) :: U_SIMP, x(0:N), u(0:N), int

  integer        :: i
  real(kind=dbl) :: int1, int2, x1, x2, x3, u1, u2, u3, &
                    h1, h2, a0, a1, a2

  x1 = x(0);  x2 = x(1);  x3 = x(2) ! First two panels
  h1 = x2 - x1;  h2 = x3 - x2
  u1 = u(0);  u2 = u(1);  u3 = u(2)
  a2 = (h2*(u1-u2) + h1*(u3-u2)) / (h1*h2*(h1+h2))
  a1 = (u2-u1) / h1 - a2*(x1+x2)
  a0 = u2 - a1*x2 - a2*x2**2
  int1 =   a0*x2 + a1/2*x2**2 + a2/3*x2**3 &
         - a0*x1 - a1/2*x1**2 - a2/3*x1**3
  int2 =   a0*x3 + a1/2*x3**2 + a2/3*x3**3 &
         - a0*x2 - a1/2*x2**2 - a2/3*x2**3
  int = int1 + int2/2
  do i = 2, N-1                     ! All others now
    x1 = x(i-1);  x2 = x(i);  x3 = x(i+1)
    h1 = x2 - x1;  h2 = x3 - x2
    u1 = u(i-1);  u2 = u(i);  u3 = u(i+1)
    a2 = (h2*(u1-u2) + h1*(u3-u2)) / (h1*h2*(h1+h2))
    a1 = (u2-u1) / h1 - a2*(x1+x2)
    a0 = u2 - a1*x2 - a2*x2**2
    int1 =   a0*x2 + a1/2*x2**2 + a2/3*x2**3 &
           - a0*x1 - a1/2*x1**2 - a2/3*x1**3
    int2 =   a0*x3 + a1/2*x3**2 + a2/3*x3**3 &
           - a0*x2 - a1/2*x2**2 - a2/3*x2**3
    int = int + (int1+int2)/2
  enddo
  int = int + int2/2
  U_SIMP = int
end function U_SIMP
