function G0FORN (u, x, n)
! Computes G0 for points u(0:n-1) at positions x(0:n-1), arbitrarily spaced.
  use STUFF;   implicit none
  integer   :: n
  real(dbl) :: G0FORN, u(0:n-1), x(0:n-1)

  real(dbl) :: deriv, w(15) ! W are the weights, dummies

  call FORNBERG (1, x(0), x, u, n, deriv, w)
  G0FORN = deriv
end function G0FORN
