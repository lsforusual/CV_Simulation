function UMDE_REF_FLUX (T)
! Computes the UMDE flux at T by interpolating in the table of known
! exact values for a range of T. If T lies outside the table's range,
! MO_I is called instead. This avoids the slight inaccuracies of MO_I
! near the cross-over point.
  use STUFF;   implicit none

  real(kind=dbl) :: UMDE_REF_FLUX, T

  integer,parameter :: N = 30
  integer        :: index
  real(kind=dbl) :: err_short, err_long, w(5), MO_I
  real(kind=dbl) :: Ttab(N) =(/0.010,  0.020,  0.030,  0.040, &
                               0.050,  0.070,  0.100,  0.150, &
                               0.200,  0.250,  0.275,  0.300, &
                               0.325,  0.350,  0.400,  0.450, &
                               0.500,  0.550,  0.600,  0.650, &
                               0.700,  0.800,  0.900,  1.000, &
                               1.100,  1.200,  1.400,  1.600, &
                               1.800, 2.00 /)
  real(kind=dbl) :: Itab(N) =(/5.2377, 3.9481, 3.3793, 3.0416, &
                               2.8120, 2.5125, 2.2477, 2.0020, &
                               1.8577, 1.7606, 1.7226, 1.6896, &
                               1.6607, 1.6351, 1.5915, 1.5556, &
                               1.5255, 1.4997, 1.4774, 1.4578, &
                               1.4404, 1.4108, 1.3864, 1.3659, &
                               1.3482, 1.3330, 1.3076, 1.2873, &
                               1.2705, 1.2564 /)
  real(kind=dbl) :: Ttablog(N), Tlog, sd

  Ttablog = Log(Ttab); Tlog = LOG(T)
  if (T < Ttab(1)) then
    UMDE_REF_FLUX = MO_I (T)
  else if (T > Ttab(N)) then
    UMDE_REF_FLUX = MO_I (T)
  else
    call FIND_INDEX (T, Ttab, 30, 5, index)
    call FORNBERG (0, Tlog, Ttablog(index),Itab(index), 5, UMDE_REF_FLUX, w)
  endif

CONTAINS
subroutine FIND_INDEX (T, Ttab, N, Npt, index)
! TO find the index best for the group of npt points
! along the table of values TTABS, for interpolation.
! It will be central to the group, or close.
  use STUFF;   implicit none

  integer        :: N, Npt, index
  real(kind=dbl) :: T, Ttab(N)

  integer        :: i

  if (ABS(T-Ttab(1)) < small) then
    index = 1
  else if (ABS(T-Ttab(N)) < small) then
    index = N - Npt + 1
  else
    do i = 1, N-1
      index = i;  if (T>Ttab(i) .and. T<= Ttab(i+1)) exit
    enddo
    index = index - Npt/2
    index = MAX (index,1);  index = MIN (index,N-Npt+1)
  endif
end subroutine FIND_INDEX
end function UMDE_REF_FLUX
