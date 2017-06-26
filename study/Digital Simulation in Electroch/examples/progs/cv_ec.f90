program CV_EC
! CV simulation of a reversible reaction on an exponentially
! expanding grid of points, using Crank-Nicolson. EC mechamism.
! The (dimensionless) current G is written out (along with p)
! into a file for plotting, and the interpolated peak value
! and its potential are reported.
! Data input:  1. p(start), p(stop), nT, X(1), N, K
!              2. output data file name
! The p's are dimensionless potential units
! nT is the number of steps per p-unit
! N is the number of internal concentration points
  use STUFF;   implicit none
  integer        :: nT, nstep, N, iT, i, jmax(1), imax, isweep
  real(kind=dbl),allocatable :: &
                    X(:), cA(:), cB(:), a1(:), a3(:), ak(:),  &
                    adA(:), adB(:), bdA(:), bdB(:), G(:)
  real(kind=dbl) :: H1, gamma, dT, a2, Xlim, pstart,prev, K, &
                    p,dp, Gmax, xmax, pmax, G0FORN, EE_FAC

  call DATA_IN (pstart, prev, nT, H1, N, K, dT, Xlim, nstep)
  ALLOCATE (X(0:N+1), cA(0:N+1),cB(0:N+1), a1(N),a3(N),ak(N), &
            adA(N), adB(N), bdA(N), bdB(N), G(2*nstep))
  cA = 1;  cB = 0              ! Initialisation; only A present
  gamma = EE_FAC (H1, N, Xlim)
  print '(" gamma  =", f12.5, " (found by iteration)")', gamma

  x(0) = 0;  x(1) = H1                            ! x-values:
  do i = 2, N+1;  x(i) = H1 * (gamma**i-1) / (gamma-1);  enddo

  call PRECALC (X, N, H1, gamma, dT, K, a1, a2, a3, ak, &
                adA, adB)                              ! coeffs
  write (4,'("# gamma =",f8.5)') gamma  ! Last 2 header lines
  write (4,'("#   p       G(dir),    G(trans)")')

  dp = - dT;    p = pstart;  i = 0
  do isweep = 1, 2
    do iT = 1, nstep
      p = p + dp
      call STEP (cA, cB, N, p, a1, a2, a3, ak, &
                 adA, adB, bdA, bdB)
      i = i + 1
      G(i) = G0FORN (CA, X, 4) ! 4-pt current approx.
      write (4,'(3f10.6)') p, G(i)
    enddo
    dp = - dp ! Reverse sweep:
  enddo
end program CV_EC


subroutine DATA_IN (pstart,prev, nT, H1, N, K, dT, Xlim, nstep)
! Reads in input data, does some calculation derived from it,
! echoes some of this to the terminal and the plotting file.
  use STUFF;  implicit none
  integer           :: nT, N, nstep
  real(kind=dbl)    :: pstart, prev, H1, dT, K, Xlim

  print '(" pstart, prev, nT, X(1), N, K?")'
  read *,   pstart, prev, nT, H1,   N, K
  call FILSPC (4)

  dT = 1.0_dbl / nT
  Xlim = 6 * SQRT(2*ABS(prev-pstart))
  nstep = NINT(nT * ABS(prev-pstart)) ! Total no. steps/sweep

  print '(" CV_EC:")'               ! Data echo etc
  print '(" nT     =", i6, " per p-unit")', nT
  print '(" pstart =", f10.3)', pstart
  print '(" prev  =", f10.3)', prev
  print '(" X(1)   =", f12.5)', H1
  print '(" N      =", i6)',    N
  print '(" hcr K  =", f10.3)', K
  print '(" Xlim   =", f10.3)', Xlim

! Header (comment) lines in plot file:
  write (1,'("# CV_EC: nT/p-unit, K =", i5,f8.2)') nT, K
  write (1,'("# pstart, prev =", 2f10.3)') pstart, prev
  write (1,'("# No. of conc. nodes between boundaries:",i6)') N
  write (1,'("# X(1), Xlim =", f10.4, f10.2)') H1, Xlim
end subroutine DATA_IN


subroutine PRECALC (X, N, H1, gamma, dT, K, &
                    a1, a2, a3, ak, adA, adB)
! To precalculate the constants.
  use STUFF;  implicit none
  integer        :: N
  real(kind=dbl) :: X(0:N+1), a1(N), a2, a3(N), H1, gamma, &
                    dT, K, ak(N), adA(N), adB(N)
                    
  integer        :: i
  real(kind=dbl) :: alpha1, alpha2, alpha3

  a2 = 1 / gamma                          ! The a's:
  do i = 1, N
    alpha1 =   2 / (x(i)-x(i-1)) / (x(i+1)-x(i-1))
    alpha2 = - 2 / (x(i)-x(i-1)) / (x(i+1)-x(i)  )
    alpha3 =   2 / (x(i+1)-x(i)) / (x(i+1)-x(i-1))
    a1(i) =   (alpha2 - 2/dT) / alpha1
    a3(i) = - (alpha2 + 2/dT) / alpha1
    ak(i) = K / alpha1
  enddo
  adA(N) = a1(N); adB(N) = a1(N) - ak(N) ! a' for A & B):
  do i = N-1, 1, -1
    adA(i) = a1(i) - a2/adA(i+1)
    adB(i) = a1(i)-ak(i) - a2/adB(i+1)
  enddo
end subroutine PRECALC


subroutine STEP (cA,cB, N, p, a1,a2,a3, ak, adA,adB, bdA,bdB)
! Solves for the next cA and cB at potential p, using implicit
! (Nernstian) boundary conditions and unequal intervals, and
! the precalculated a-coefficients.
  use STUFF;  implicit none
  integer        :: N
  real(kind=dbl) :: cA(0:N+1), cB(0:N+1), p, a1(N), a2, a3(N),&
                    ak(N), adA(N), adB(N), bdA(N), bdB(N)

  integer        :: i
  real(kind=dbl) :: bAi, bBi, Pee(2,2), Q(2), C0(2)

! Generating the 2 sets of b':
  bAi = - cA(N-1) + a3(N)*cA(N) - a2*cA(N+1)
  bdA(N) = bAi - a2*cA(N+1)
  bBi = - cB(N-1) + (a3(N)+ak(N))*cB(N) - a2*cB(N+1)
  bdB(N) = bBi - a2*cB(N+1)
  do i = N-1, 1, -1
    bAi = - cA(i-1) + a3(i)*cA(i) - a2*cA(i+1)
    bdA(i) = bAi - a2*bdA(i+1)/adA(i+1)
    bBi = - cB(i-1) + (a3(i)+ak(i))*cB(i) - a2*cB(i+1)
    bdB(i) = bBi - a2*bdB(i+1)/adB(i+1)
  enddo

! The boundary conditions. 2-pt makes it simple, no U-V.
  Pee (1,1) = 1;             Pee (1,2) = -EXP(p)
  Pee (2,1) = 1 + 1/adA(1);  Pee (2,2) = 1 + 1/adB(1)
  Q(1) = 0;  Q(2) = bdA(1) / adA(1) + bdB(1) / adB(1)
  call MATINV (Pee, 2)    ! 2*2 mat inversion
  C0(:) = MATMUL (Pee, Q)
  CA(0) = C0(1);  CB(0) = C0(2)
  
! Forward again, replacing all c with c' values:
  do i = 1, N
      cA(i) = (bdA(i) - cA(i-1)) / adA(i)
      cB(i) = (bdB(i) - cB(i-1)) / adB(i)
  enddo
end subroutine STEP
