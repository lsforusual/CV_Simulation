subroutine MINMAX (y1, y2, y3, y0, x0)
! From the 3 y-values, assumed to lie at x = (-1, 0 and +1),
! a minmax Y0 value is computed, by fitting a parabola to them.
! Y0 is returned together with the x value, X0, at Y0.
  use STUFF;  implicit none
  real(kind=dbl) :: y1, y2, y3, y0, x0

  real(kind=dbl) :: a0, a1, a2

  a0 = y2
  a2 = (y1 - 2*y2  + y3) / 2
  a1 = a0 + a2 - y1
  x0 = -a1 / 2 / a2
  y0 = a0 + a1*x0 + a2*x0**2
end subroutine MINMAX
