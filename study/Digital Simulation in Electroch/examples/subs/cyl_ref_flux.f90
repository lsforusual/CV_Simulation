function CYL_REF_FLUX (T)
! To evaluate, by interpolation, the F(theta) function of Jaeger and Clark.
! We multiply it by pi^2/4 to let it conform to the cylinder current function
! of Szabo et al (1987).
  use STUFF;   implicit none
  integer,parameter :: N=30
  integer           :: index
  real(kind=dbl)    :: CYL_REF_FLUX, T, Ttablog(N), Tlog, w(5)
  real(kind=dbl)    :: Ttab(1:N) = (/ &
       0.01_dbl,  0.02_dbl,  0.03_dbl,  0.05_dbl,  0.07_dbl,           &
       0.10_dbl,  0.15_dbl,  0.2_dbl,   0.3_dbl,   0.5_dbl,   0.7_dbl, &
       1.0_dbl,   1.5_dbl,   2.0_dbl,   3.0_dbl,   5.0_dbl,   7.0_dbl, &
      10.0_dbl,  15.0_dbl,  20.0_dbl,  30.0_dbl,  50.0_dbl,  70.0_dbl, &
     100.0_dbl, 150.0_dbl, 200.0_dbl, 300.0_dbl, 500.0_dbl, 700.0_dbl, &
     1000.0_dbl /)
  real(kind=dbl)    :: JCtab(1:N) = (/ &
     15.122_dbl, 11.033_dbl, 9.218_dbl, 7.394_dbl, 6.421_dbl,            &
      5.549_dbl,  4.726_dbl, 4.232_dbl, 3.643_dbl, 3.044_dbl, 2.720_dbl, &
      2.427_dbl,  2.147_dbl, 1.975_dbl, 1.767_dbl, 1.550_dbl, 1.429_dbl, &
      1.317_dbl,  1.207_dbl, 1.138_dbl, 1.052_dbl, 0.958_dbl, 0.904_dbl, &
      0.853_dbl,  0.802_dbl, 0.768_dbl, 0.725_dbl, 0.676_dbl, 0.648_dbl, &
      0.620_dbl /)

!print '("  i      T       JC       pi^2/4*JC")'
!do i = 1, N
!print '(i4, f8.2, 3f10.3)', i, Ttab(i), JCtab(i), pi**2/4*JCtab(i)
!enddo
Ttablog = Log(Ttab); Tlog = LOG(T)
  if (T < Ttab(1)) then
    CYL_REF_FLUX = 0
  else if (T > Ttab(N)) then
    CYL_REF_FLUX = 0
  else
    call FIND_INDEX (T, Ttab, 30, 5, index)
    call FORNBERG (0, Tlog, Ttablog(index), JCtab(index), 5, CYL_REF_FLUX, w)
    CYL_REF_FLUX = CYL_REF_FLUX * pi**2 / 4
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
end function CYL_REF_FLUX
