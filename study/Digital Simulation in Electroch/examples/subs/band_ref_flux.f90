function BAND_REF_FLUX (T)
! Computes the UMDE flux at T by interpolating in the table of known
! exact values for a range of T. If T lies outside the table's range,
! MO_I is called instead. This avoids the slight inaccuracies of MO_I
! near the cross-over point.
  use STUFF;   implicit none
  real(kind=dbl) :: T, BAND_REF_FLUX

  integer,parameter :: N = 20
  integer        :: index
  real(kind=dbl) :: err_short, err_long, w(5), AOKI_BAND_I
  real(kind=dbl) :: Ttab(N) = &
                    (/   0.5_dbl,    0.7_dbl,    1.0_dbl,    1.5_dbl, &
                         2.0_dbl,    3.0_dbl,    5.0_dbl,    7.0_dbl, &
                        10.0_dbl,   15.0_dbl,   20.0_dbl,   30.0_dbl, &
                        50.0_dbl,   70.0_dbl,  100.0_dbl,  150.0_dbl, &
                       200.0_dbl,  300.0_dbl,  500.0_dbl,  700.0_dbl /)
  real(kind=dbl) :: Itab(N) = &
                    (/ 2.5948_dbl, 2.3462_dbl, 2.1226_dbl, 1.9092_dbl, &
                       1.7792_dbl, 1.6205_dbl, 1.4534_dbl, 1.3593_dbl, &
                       1.2711_dbl, 1.1826_dbl, 1.1265_dbl, 1.0553_dbl, &
                       0.9767_dbl, 0.9306_dbl, 0.8861_dbl, 0.8402_dbl, &
                       0.8103_dbl, 0.7714_dbl, 0.7273_dbl, 0.7008_dbl /)
  real(kind=dbl) :: Ttablog(N), Tlog, sd

  Ttablog = Log(Ttab); Tlog = LOG(T)
  if (T < Ttab(1)) then
    BAND_REF_FLUX = AOKI_BAND_I (T)
  else if (T > Ttab(N)) then
    BAND_REF_FLUX = AOKI_BAND_I (T)
  else
    call FIND_INDEX (T, Ttab, N, 5, index) ! Interpolate, using FORNBERG
    call FORNBERG (0, Tlog, Ttablog(index),Itab(index), 5, BAND_REF_FLUX, w)
  endif

CONTAINS
subroutine FIND_INDEX (T, Ttab, N, Npt, index)
! Finds the index best suited for the group of npt points
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
end function BAND_REF_FLUX
