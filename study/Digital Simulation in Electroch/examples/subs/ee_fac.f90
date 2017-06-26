function EE_FAC (H1, N, lim, ierr)
! Binary search for the gamma value that fits (H1, N, lim),
! i.e. such that the series sum of h(i), i = 1, N is LIM,
! and h(i) = h(i-1) * gamma.
! If the function value exceeds 2, or the number of iterations
! hits the maximum allowed, the error condition ierr = 1 is
! returned, else ierr = 0. The calling program must check for this.
  use STUFF;  implicit none
  integer        :: N, ierr
  real(kind=dbl) :: EE_FAC, H1, lim

  integer,parameter        :: nitmax = 30
  real(kind=dbl),parameter :: eps = 1.0E-08_dbl
  integer        :: i
  real(kind=dbl) :: a, b, f, gamma

  a = 1;  b = 2     ! The starting limits for gamma
  do i = 1, nitmax
    gamma = (a + b) / 2
    f = (gamma**N - 1)/(gamma-1) - lim / H1 ! Func value
    if (f <= 0) then
       a = gamma
    else
       b = gamma
    endif
    if (ABS(b-a) <= eps .OR. i == nitmax .OR. gamma > 1.999) exit
  enddo
  gamma = (a + b) / 2  ! A final correction
  EE_FAC = gamma
  ierr = 0
  if (gamma > 1.999 .OR. i >= nitmax) ierr = 1 ! Bad gamma flag
end function EE_FAC
