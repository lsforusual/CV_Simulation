function SV_FAC (H1, N, lim)
! Binary search for the alpha value that fits (H1, N, lim),
! i.e. such that the series sum of h(i), i = 1, N is LIM,
! H(i) = H(i-1) * (1 + alpha*H(i-1)/H(1)), the S&V sequence.
  use STUFF;  implicit none
  integer        :: N
  real(kind=dbl) :: SV_FAC, H1, H, X, lim

  integer,parameter        :: nitmax = 30
  real(kind=dbl),parameter :: eps = 1.0E-08_dbl
  integer        :: i, j
  real(kind=dbl) :: a, b, f, alpha

  a = 0;  b = 1     ! The starting range for alpha (a safe bet)
  do j = 1, nitmax
    alpha = (a + b) / 2
    X = H1;  H = H1
    do i = 2, N
      H = H * (1 + alpha*H/H1)
      X = X + H
    enddo
    f = X - lim ! Func value
    if (f <= 0) then
       a = alpha
    else
       b = alpha
    endif
    if (ABS(b-a) <= eps) exit
    if (j == nitmax) STOP " Not possible to find alpha. Exitus."
  enddo
  alpha = (a + b) / 2  ! A final correction
  SV_FAC = alpha
end function SV_FAC
