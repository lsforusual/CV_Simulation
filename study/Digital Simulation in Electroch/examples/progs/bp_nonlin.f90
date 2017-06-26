program BP_NONLIN
! CN simulation of the Birk & Perone 2nd-order mechanism with
! Cottrell-like conditions (C0 = 0). The nonlinear term, C^2,
! is left as is, and the Newton procedure is used.
! An expanding grid of points is used.
! The current G is written out with its error for plotting.
! Data input:  1. nT, X(1), N, np, K
!              2. output data file name
! X(1) (H1) is the first X position away from the electrode
! N is the number of internal concentration points
! np is the number of points used for the G-approximation.
! K is the dimensionless 2nd-order rate constant.
  use STUFF;   implicit none
  integer        :: nT, N, iT, i, np, inext, Npear, nit
  real(kind=dbl),allocatable :: &
                    X(:), C(:), D(:), dD(:), a1(:), a3(:),    &
                    ak(:), ad(:), bd(:), F(:)
  real(kind=dbl) :: H1, gamma, dT, dTP, K, a2, Xlim, bpa(10), &
                    rat, G, G_anal, err, Cbnew, T,            &
                    G0FORN, EE_FAC, GBPFUNC

  call DATA_IN (nT, H1, N, np, dT, K)
  Xlim = 6
  ALLOCATE (X(0:N+1), C(0:N+1), D(0:N+1), dD(0:N),            &
            a1(N), ak(N), a3(N), ad(N), bd(N), F(N))
  C = 1; C(0) = 0              ! Initialisation; only A present
  gamma = EE_FAC (H1, N, Xlim)
  bpa(1) = 4/pi - 1
  bpa(2:10) = (/ 0.08327_dbl, 0.02893_dbl, 0.01162_dbl, &
                 0.00540_dbl, 0.00286_dbl, 0.00169_dbl, &
                 0.00108_dbl, 0.00074_dbl, 0.00053_dbl /)
  print '(" gamma  =", f12.5, " (found by iteration)")', gamma
  x(0) = 0;  x(1) = H1                            ! x-values:
  do i = 2, N+1;  x(i) = H1 * (gamma**i-1) / (gamma-1);  enddo

  write (4,'("# gamma =",f8.5)') gamma   ! Last 2 header lines
  Npear = NINT (dT/H1**2) ! Max. lambda sets no. Pearson steps
  print '(i6," Pearson substeps in the first interval.")',Npear
  print '(6x, "T", 9x, "G", 5x, "LOG|rel.err|", 12x, "No. iter")'
  dTP = dT / Npear        ! Sub-dT for Pearson steps
  call COEFFS (X, N, gamma, K, dTP, a1, ak, a2, a3)
  do i = 1, Npear
    T = i * dTP;  Cbnew = 1 / (1 + K*T)
    call STEP (C, D, dD, N, K, gamma, dTP, Cbnew,  &
                  a1, ak, a2, a3, F, ad, bd, nit)
  enddo
  G = G0FORN (c, x, np); G_anal = GBPFUNC (K, dT, bpa)
  err = LOG10(ABS(G-G_anal)/G_anal)
  print '(2f10.3, f10.2, i20)', dT, G, err, nit
  write (4,'(3f10.6)') dT, G, err
  inext = 2
  call COEFFS (X, N, gamma, K, dT, a1, ak, a2, a3)
  do iT = 2, nT
    T = i * dT;  Cbnew = 1 / (1 + K*T)
    call STEP (C, D, dD, N, K, gamma, dT, Cbnew,   &
                 a1, ak, a2, a3, F, ad, bd, nit)
    G = G0FORN (c, x, np); G_anal = GBPFUNC (K, iT*dT, bpa)
    err = LOG10(ABS(G-G_anal)/G_anal)
    if (iT == inext) then
      print '(2f10.3, f10.2, i20)', iT*dT, G, err, nit
      inext = 2 * inext;  if (inext > nT) inext = nT
    endif
    write (4,'(3f10.6)') iT*dT, G, err
  enddo
end program BP_NONLIN

subroutine DATA_IN (nT, H1, N, np, dT, K)
! Reads in input data, does some calculation derived from it,
! echoes some of this to the terminal and the plotting file.
  use STUFF;  implicit none
  integer           :: nT, N, np
  real(kind=dbl)    :: H1, dT, K

  print '(" nT, X(1), N, np, K?")'
  read *,   nT, H1,   N, np, K
  call FILSPC (4)

  dT = 1.0_dbl / nT
! Data echo etc
  print'(" BP_NONLIN (Birk & Perone potential-step):")'
  print '(" nT     =", i6)', nT
  print '(" X(1)   =", f12.5)', H1
  print '(" N      =", i6)',    N
  print '(" G, using", i6, " points.")', np
  print '(" K      =", f9.2)', K

! Header (comment) lines in plot file:
  write (1,'("# BP_LIN: nT, K", i5, f8.2)') nT, K
  write (1,'("# No. of conc. nodes between boundaries:",i6)') N
  write (1,'("# X(1) =", f10.4, f10.2)') H1
end subroutine DATA_IN

subroutine COEFFS (X, N, gamma, K, dT, a1, ak, a2, a3)
! To precalculate the constants.
  use STUFF;  implicit none
  integer        :: N
  real(kind=dbl) :: X(0:N+1), gamma, K, dT, a1(N), &
                    a2, ak(N), a3(N)

  integer        :: i
  real(kind=dbl) :: alpha1, alpha2, alpha3

  a2 = 1 / gamma                          ! The a's:
  do i = 1, N
    alpha1 =   2 / (x(i)-x(i-1)) / (x(i+1)-x(i-1))
    alpha2 = - 2 / (x(i)-x(i-1)) / (x(i+1)-x(i)  )
    alpha3 =   2 / (x(i+1)-x(i)) / (x(i+1)-x(i-1))
    a1(i)  =   (alpha2 - 2/dT) / alpha1
    ak(i)  =   - K / alpha1
    a3(i)  =   (alpha2 + 2/dT) / alpha1
  enddo
end subroutine COEFFS

subroutine STEP (C, D, dD, N, K, gamma, dT, Cbnew,  &
                   a1, ak, a2, a3, F, ad, bd, nit)
! Solves for the next C by Newton iteration.
  use STUFF;  implicit none
  integer        :: N, nit
  real(kind=dbl) :: C(0:N+1), D(0:N+1), dD(0:N), K, gamma, &
                    dT, Cbnew, a1(N), ak(N), a2, a3(N),    &
                    ad(N), bd(N), F(N)

  integer        :: i
  real(kind=dbl) :: A_D, B, norm, eps=1.0e-08_dbl

  D = C;  D(N+1) = Cbnew;  nit = 0
  do
    nit = nit + 1
    do i = 1, N       ! Generate F = AD - B
      A_D = D(i-1) + a1(i)*D(i) + ak(i)*D(i)**2 + a2*D(i+1)
      B  = -C(i-1) - a3(i)*C(i) - ak(i)*C(i)**2 - a2*C(i+1)
      F(i) = A_D - B
    enddo
!   Generating a' and b' for the Jacobian:
    ad(N) = a1(N) + 2*ak(N)*D(N)
    bd(N) = - F(N)
    do i = N-1, 1, -1
      ad(i) = a1(i) + 2*ak(i)*D(i) - a2/ad(i+1)
      bd(i) = -F(i) - a2*bd(i+1)/ad(i+1)
    enddo

!   Forward again, now calculating the d's:
    dD(0) = 0
    do i = 1, N
      dD(i) = (bd(i) - dD(i-1)) / ad(i)
    enddo
    D(1:N) = D(1:N) + dD(1:N) ! Correction vector added
    norm = SQRT (SUM (dD(1:N)**2) / N) ! Correction vector norm
    if (norm < eps) exit
  enddo
  C(1:N+1) = D(1:N+1)
end subroutine STEP


function GBPFUNC (K, T, a)
! To evaluate the current at T, given K, from the
! formula in Britz + Kastening 1974.
  use STUFF;  implicit none
  real(kind=dbl) :: GBPFUNC, K, T, a(10)

  integer        :: i
  real(kind=dbl) :: theta, theta_term, term, sum, ratio

  theta = K * T;  theta_term = theta / ( 1 + theta)
  sum = 1;  term = theta_term
  do i = 1, 10
    sum = sum + a(i) * term
    term = term * theta_term ! Next power please
  enddo
  ratio = (1 + theta) / sum
  GBPFUNC = 1 / SQRT(pi*T) / ratio
end function GBPFUNC
