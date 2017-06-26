program LSV_CN
! LSV simulation of a reversible reaction on an exponentially
! expanding grid of points, using Crank-Nicolson.
! The (dimensionless) current G is written out (along with p)
! into a file for plotting, and the interpolated peak value
! and its potential are reported.
! Data input:
! 1. p(start), p(stop), nT, X(1), N, np
! 2. output data file name
! The p's are dimensionless potential units
! nT is the number of steps per p-unit
! X(1) is the first X position away from the electrode
! N is the number of internal concentration points
! np is the number of points used for the G-approximation.
  use STUFF;   implicit none
  integer        :: nT, nstep, N, iT, i, np, jmax(1), imax
  real(kind=dbl),allocatable :: &
                    X(:), cA(:), cB(:), a1(:), a3(:), ad(:), &
                    bdA(:), bdB(:), G(:), p(:)
  real(kind=dbl) :: H1, gamma, dT, a2, Xlim, pstart, pstop,   &
                    dp, Gmax, xmax, pmax, EE_FAC, GU

  call DATA_IN (pstart, pstop, nT, H1, N, np, dT, Xlim, nstep)
  ALLOCATE (X(0:N+1), cA(0:N+1), cB(0:N+1), a1(N), a3(N),    &
            ad(N), bdA(N), bdB(N), G(nstep), p(0:nstep))
  cA = 1;  cB = 0              ! Initialisation; only A present
  gamma = EE_FAC (H1, N, Xlim)
  print '(" gamma  =", f12.5, " (found by iteration)")', gamma
  x(0) = 0;  x(1) = H1                            ! x-values:
  do i = 2, N+1;  x(i) = H1 * (gamma**i-1) / (gamma-1);  enddo
  call PRECALC (X, N, H1, gamma, dT, a1, a2, a3, ad) ! coeffs
  write (1,'("# gamma =",f8.5)') gamma  ! Last 2 header lines
  write (1,'("#   p       G(dir),    G(trans)")')
  dp = - dT;    p(0) = pstart
  do iT = 1, nstep
    p(iT) = p(iT-1) + dp
    call STEP (cA, cB, N, p(iT), np, a1, a2, a3, ad, bdA, bdB)
    G(iT) = GU (cA, x, np)
    write (1,'(3f10.6)') p(iT), G(iT)
  enddo
  jmax = MAXLOC (G);  imax = jmax(1);  Gmax = G(imax)
  call MINMAX (G(imax-1), G(imax), G(imax+1), Gmax, xmax)
  pmax = p(imax) + xmax*dp
  print '(//" G-peak of", f10.4, " found at p =", f10.4)', &
                          Gmax,                    pmax
end program LSV_CN

subroutine DATA_IN (pstart,pstop,nT, H1, N, np, dT, Xlim, nstep)
! Reads in input data, does some calculation derived from it,
! echoes some of this to the terminal and the plotting file.
  use STUFF;  implicit none
  character(len=40) :: filename
  integer           :: nT, N, np, nstep
  real(kind=dbl)    :: pstart, pstop, H1, dT, Xlim

  print '(" pstart, pstop, nT, X(1), N, np?")'
  read *,   pstart, pstop, nT, H1,   N, np
  print '(" Output file name?")';  read (*,'(a)') filename
  OPEN (unit=1, file=filename, status='new', action='write')

  dT = 1.0_dbl / nT
  Xlim = 6 * SQRT(ABS(pstop-pstart))
  nstep = NINT(nT * ABS(pstop-pstart)) ! Total no. of steps

  print '(" LSV_CN_UN:")'               ! Data echo etc
  print '(" nT     =", i6, " per p-unit")', nT
  print '(" pstart =", f10.3)', pstart
  print '(" pstop  =", f10.3)', pstop
  print '(" X(1)   =", f12.5)', H1
  print '(" N      =", i6)',    N
  print '(" Xlim   =", f10.3)', Xlim
  print '(" G, using", i6, " points.")', np

! Header (comment) lines in plot file:
  write (1,'("# LSV_CN_UN: nT/p-unit, np =", 2i5)') nT, np
  write (1,'("# pstart, pstop =", 2f10.3)') pstart, pstop
  write (1,'("# No. of conc. nodes between boundaries:",i6)') N
  write (1,'("# X(1), Xlim =", f10.4, f10.2)') H1, Xlim
end subroutine DATA_IN


subroutine PRECALC (X, N, H1, gamma, dT, a1, a2, a3, ad)
! To precalculate the constants.
  use STUFF;  implicit none
  integer        :: N
  real(kind=dbl) :: X(0:N+1), a1(N),a2,a3(N), H1, gamma, &
                    dT, ad(N)
                    
  integer        :: i
  real(kind=dbl) :: alpha1, alpha2, alpha3

  a2 = 1 / gamma                          ! The a's:
  do i = 1, N
    alpha1 =   2 / (x(i)-x(i-1)) / (x(i+1)-x(i-1))
    alpha2 = - 2 / (x(i)-x(i-1)) / (x(i+1)-x(i)  )
    alpha3 =   2 / (x(i+1)-x(i)) / (x(i+1)-x(i-1))
    a1(i) =   (alpha2 - 2/dT) / alpha1
    a3(i) = - (alpha2 + 2/dT) / alpha1
  enddo
  ad(N) = a1(N)                         ! a' (same for A & B):
  do i = N-1, 1, -1
    ad(i) = a1(i) - a2/ad(i+1)
  enddo
end subroutine PRECALC

subroutine STEP (cA, cB, N, p, np, a1, a2, a3, ad, bdA, bdB)
! Solves for the next cA and cB at potential p, using implicit
! (Nernstian) boundary conditions and unequal intervals, and
! the precalculated a-coefficients.
  use STUFF;  implicit none
  integer        :: N, np
  real(kind=dbl) :: cA(0:N+1), cB(0:N+1), p, a1(N), a2, a3(N),&
                    ad(N), bdA(N), bdB(N)

  integer        :: i
  real(kind=dbl) :: bAi, bBi, Pee, Q, R, G0FUNC
  real(kind=dbl),allocatable :: uA(:), uB(:), vA(:), vB(:)

! Generating the 2 sets of b':
  bAi = - cA(N-1) + a3(N)*cA(N) - a2*cA(N+1)
  bdA(N) = bAi - a2*cA(N+1)
  bBi = - cB(N-1) + a3(N)*cB(N) - a2*cB(N+1)
  bdB(N) = bBi - a2*cB(N+1)
  do i = N-1, 1, -1
    bAi = - cA(i-1) + a3(i)*cA(i) - a2*cA(i+1)
    bdA(i) = bAi - a2*bdA(i+1)/ad(i+1)
    bBi = - cB(i-1) + a3(i)*cB(i) - a2*cB(i+1)
    bdB(i) = bBi - a2*bdB(i+1)/ad(i+1)
  enddo

! The u's and v's; u(0) = 0 and v(0) = 1.
  ALLOCATE (uA(0:np-1), vA(0:np-1), uB(0:np-1), vB(0:np-1))
  uA(0) = 0;  uB(0) = 0;     vA(0) = 1; vB(0) = 1
  do i = 1, np-1
      uA(i) = (bdA(i) - uA(i-1)) / ad(i)
      vA(i) = -vA(i-1) / ad(i)
      uB(i) = (bdB(i) - uB(i-1)) / ad(i)
      vB(i) = -vB(i-1) / ad(i)
  enddo
  Pee = G0FUNC (vA, np, 1.0_dbl);  Q = G0FUNC (vB, np, 1.0_dbl)
  R = -G0FUNC (uA+uB, np, 1.0_dbl)
  cA(0) = R / (Pee + EXP(-p)*Q)
  cB(0) = EXP(-p) * cA(0)
  DEALLOCATE (uA, vA, uB, vB)
! Forward again, replacing all c with c' values:
  do i = 1, N
      cA(i) = (bdA(i) - cA(i-1)) / ad (i)
      cB(i) = (bdB(i) - cB(i-1)) / ad(i)
  enddo
end subroutine STEP
