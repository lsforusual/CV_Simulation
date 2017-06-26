function G0BETA (i, n)
! Returns the beta parameter at i given n of the n-point
! approximation), i.e. for the form y'_0(n) for n: 2..7.
! If i or n is outside the respective range, zero is returned.
  use STUFF;  implicit none
  real(kind=dbl) :: G0BETA
  integer        :: i, n

  integer        :: numer(0:6,2:7) = &
                      (/ -1,   1,    0,   0,    0,   0,   0,  &
                         -3,   4,   -1,   0,    0,   0,   0,  &
                        -11,  18,   -9,   2,    0,   0,   0,  &
                        -25,  48,  -36,  16,   -3,   0,   0,  &
                       -137, 300, -300, 200,  -75,  12,   0,  &
                       -147, 360, -450, 400, -225,  72, -10   /)
  integer        :: denom(2:7) = (/  1, 2, 6, 12, 60, 60 /)
  save           :: numer, denom

  if (i>=0.and.i<=6  .and. n>=2.and.n<=7) then
    G0BETA = REAL(numer(i,n),dbl) / REAL(denom(n),dbl)
  else
    G0BETA = 0
  endif
end function G0BETA
