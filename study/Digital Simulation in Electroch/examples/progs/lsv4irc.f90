program LSV4IRC
! To simulate a reversible LSV sweep, including uncompensated
! resistance and capacity effects. Unequal intervals and a
! 4-pt approximation in space, and LMG-2 are used. The
! nonlinear boundary condition system is solved by Newton
! iteration. Vector U contains the 6 boundary values. Because
! of the extrapolation, it must be represented twice.
! Data input:
! 1. output data file name
! 2. H1, N (first interval, total no. of points in X)
! 3. p1, p2, nper (start & finish p, no. pts per p-unit)
! 4. gamC, rho (dimensionless capacity and R_u)
  use STUFF;  implicit none
  integer                    :: N, np, nper, ip, i, nitmax, nit
  real(kind=dbl)             :: dT, dp, H1, Xlim, GA, GB, GC, &
                                p1, p2, p, gamC, rho, gamma,  &
                                EE_FAC, GU
  real(kind=dbl),allocatable :: X(:), CA(:),CB(:),DA(:),DB(:),&
                                a1(:), a2(:), a3(:), a4(:),   &
                                b1(:), b2(:), b3(:), b4(:),   &
                                ad(:), bdA(:), bdB(:),        &
                                U1(:), U2(:)

  call DATA_IN (H1, N, p1, p2, nper, gamC, rho)
  dT = 1.0_dbl / nper
  ALLOCATE (CA(0:N+1), CB(0:N+1), DA(0:N+1), DB(0:N+1),  &
            X(0:N+2), a1(N), a2(N), a3(N), a4(N), b1(N), &
            b2(N), b3(N), b4(N), ad(N), bdA(N), bdB(N),  &
            U1(6), U2(6))
  CA = 1; CB = 0; DA = 1;  DB = 0            ! C's initialised
  X(0) = 0;   X(1) = H1;  Xlim = 6*SQRT(ABS(p2-p1)) ! X-values:
  gamma = EE_FAC (X(1), N, Xlim)
  do i = 2, N+2;  X(i) = H1 * (gamma**i-1)/(gamma-1); enddo
  np = ABS(p2 - p1) * nper
  call HEADERS (N, np, nper, p1, p2, H1, gamma, gamC, rho)

  call COEFFS (X, N, dT,   a1, a2, a3, a4)
  call COEFFS (X, N, dT/2, b1, b2, b3, b4)
  dp = -dT;  p = p1 - dp ! First step goes to start potential
  U1(1) = 1;  U1(2) = 0;  U1(3:5) = 0;  U1(6) = p
  U2 = U1
  nitmax = 0
  do ip = 1, np                                ! The sweep
    p = p + dp/2
    DA = CA;  DB = CB    ! Copy the present concs into the D's
    call STEP (DA, DB, X, N, b1, b2, b3, b4, ad, bdA, bdB, &
               p, dT/2, gamC, rho, U2, nit)        ! half step
    p = p + dp/2
    call STEP (DA, DB, X, N, b1, b2, b3, b4, ad, bdA, bdB, &
               p, dT/2, gamC, rho, U2, nit)        ! & another
    call STEP (CA, CB, X, N, a1, a2, a3, a4, ad, bdA, bdB, &
               p, dT, gamC, rho, U1, nit)          ! Whole step
    CA = 2*DA - CA;   CB = 2*DB - CB        ! The extrapolation
    U1 = 2*U2 - U1          ! The boundary values extrapolated
    write (1,'(5f10.4)') p, U1(6), U1(3), U1(5), U1(3)+U1(5)
! writing out           pnom, p,    GA,    GC,    GA + GC
    nitmax = MAX (nitmax, nit)
  enddo
  print '(" Max no. iter =", i6)', nitmax
  DEALLOCATE (CA, CB, DA, DB, X, a1, a2, a3, a4, b1, b2, b3, &
                                     b4, ad, bdA, bdB, U1, U2)
end program LSV4IRC


subroutine DATA_IN (H1, N, p1, p2, nper, gamC, rho)
  use STUFF;   implicit none
  integer        :: N, nper
  real(kind=dbl) :: H1, p1, p2, gamC, rho

  character(len=40) :: filename

  print '(" Output data file name?")'; read '(a40)', filename
  OPEN (unit=1, file=TRIM(filename), status='new')  ! Data file
  print '(" X(1), N?")';      read *, H1, N
  print '(" p1, p2, nper?")'; read *, p1, p2, nper
  print '(" gamC, rho?")';    read *, gamC, rho
end subroutine DATA_IN


subroutine HEADERS (N, np, nper, p1, p2, H1, gamma, gamC, rho)
  use STUFF;   implicit none
  integer        :: N, np, nper
  real(kind=dbl) :: p1, p2, H1, gamma, gamC, rho

  print '(" LSV, 4-pt 2nd der, with iR & C-dl.")'
  print '(" p1, p2   =", 2f9.2)', p1, p2
  print '(" rho      =", f9.2)', rho
  print '(" gam_c    =", f9.2)', gamc
  print '(" n/p unit =", i6)', nper
  print '(" X(1)     =", f10.3)', H1
  print '(" N        =", i6, " pts along X.")', N
  print '(" gamma    =", f12.5, " (found by iteration)")', gamma
  write (1,'("# LSV sim., 4-pt LMG-2, with iR & C).")')
  write (1,'("# X(1), N =", f8.4, i6)') H1, N
  write (1,'("# gamma =",f8.5)') gamma  ! Last 2 header lines
  write (1,'("# p1, p2, nper =", 2f8.2, i6)') p1, p2, nper
  write (1,'("# dimensionless gam_c, rho =", 2f8.3)') gamC, rho
  write (1,'("#", 3x, "p(nom)", 4x,"p(act)", 6x, "GA", 8x, &
             & "GC", 5x, "G(total)")')
end subroutine HEADERS


subroutine COEFFS (X, N, dT, a1, a2, a3, a4)
! To precalculate the constants.
  use STUFF;  implicit none
  integer        :: N
  real(kind=dbl) :: X(0:N+1), K, a1(N), a2(N),  &
                    a3(N), a4(N), dT
  integer        :: i
  real(kind=dbl) :: dumval, dummy(4), alfa(4)

  do i = 1, N ! The four-point formulas:
    call U_DERIV (dummy(1:4), X(i-1:i+2), 2, 4, 2, dumval, alfa)
    a1(i) =   (alfa(2) - 1/dT) / alfa(1)
    a2(i) = alfa(3) / alfa(1)
    a3(i) = alfa(4) / alfa(1)
    a4(i) = -1 / dT / alfa(1)
  enddo
end subroutine COEFFS


subroutine STEP (CA, CB, X, N, a1, a2, a3, a4, ad, bdA, bdB, &
                                     p, dT, gamC, rho, U, nit)
! To solve the BI system, by the backwards/forwards
! scheme. Reversible boundary condition at potential p, iR & C
  use STUFF; implicit none
  integer        :: N, nit
  real(kind=dbl) :: CA(0:N+1), CB(0:N+1), X(0:N+2), a1(N), &
                    a2(N), a3(N), a4(N), ad(N), bdA(N),    &
                    bdB(N), p, dT, gamC, rho, U(6)
  integer        :: i
  real(kind=dbl) :: bi, dummy(3)

! Backwards from C(N) to generate a' and b' values recursively
! for A and B separately. The ad's are common but the bd's not.
  ad(N) = a1(N)
  bi = a4(N) * CA(N); bdA(N) = bi - (a2(N)+a3(N))*CA(N+1)
  bi = a4(N) * CB(N); bdB(N) = bi - (a2(N)+a3(N))*CB(N+1)
  ad(N-1) = a1(N-1) - a2(N-1)/ad(N)
  bi = a4(N-1) * CA(N-1)
  bdA(N-1) = bi - a3(N-1)*CA(N+1) - a2(N)*bdA(N)/ad(N)
  bi = a4(N-1) * CB(N-1)
  bdB(N-1) = bi - a3(N-1)*CB(N+1) - a2(N)*bdB(N)/ad(N)
  do i = N-2, 1, -1
    ad(i) = a1(i) - (a2(i) - a3(i)/ad(i+2)) / ad(i+1)
    bi = a4(i)*CA(i)
    bdA(i) = bi - a2(i) * bdA(i+1)/ad(i+1) &
               - a3(i) * (bdA(i+2) - bdA(i+1)/ad(i+1)) / ad(i+2)
    bi = a4(i)*CB(i)
    bdB(i) = bi - a2(i) * bdB(i+1)/ad(i+1) &
               - a3(i) * (bdB(i+2) - bdB(i+1)/ad(i+1)) /ad(i+2)
  enddo
! C_0's etc, all packed in U:
  call BOUND_VALS (CA, CB, X, N, ad, bdA, bdB, p, dT, gamC, &
                                                 rho, U, nit)
  CA(0) = U(1);  CB(0) = U(2)

! Forward again, replacing all C with C' values:
  do i = 1, N
    CA(i) = (bdA(i) - CA(i-1)) / ad(i)
    CB(i) = (bdB(i) - CB(i-1)) / ad(i)
  enddo
end subroutine STEP


subroutine BOUND_VALS (CA, CB, X, N, ad, bdA, bdB, pnom, dT,  &
                                             gamC, rho, U, nit)
! To solve for the boundary unknowns, given the other current
! parameters. Pnom is the desired, sweep, potential. Vector U
! contains present values and is updated to the new ones, in
! order: CA(0), CB(0), GA, GB, GC, p.
! Currents are approximated by the 3-point formula, using GU.
  use STUFF;   implicit none
  integer        :: N, nit
  real(kind=dbl) :: CA(0:N+1), CB(0:N+1), X(0:N+2), ad(N),  &
                    bdA(N), bdB(N), pnom, dT, GA, GC, gamC, &
                    rho, U(6)

  integer        :: i, nitmax=20, ip(6), ier
  real(kind=dbl) :: uA(0:2), uB(0:2), vA(0:2), vB(0:2), p,    &
                    VeeA, VeeB, YUA, YUB, J(6,6), F(6), GU

  uA(0) = 0;  uB(0) = 0;  vA(0) = 1;  vB(0) = 1
  do i = 1, 2 ! 3-pt approx on unequal grid
    uA(i) = (bdA(i) - uA(i-1)) / ad(i)
    uB(i) = (bdB(i) - uB(i-1)) / ad(i)
    vA(i) = - vA(i-1) / ad(i);      vB(i) = - vB(i-1) / ad(i)
  enddo

  GA = U(3);  GC = U(5) ! Must keep old G values for use in F
  VeeA = - GU (vA, X, 3);  VeeB = - GU (vB, X, 3)
  YUA  =   GU(uA, X, 3);   YUB  =   GU (uB, X, 3)
  do i = 1, nitmax  ! The Newton iteration:
    nit = i
    J = 0;  F = 0
    J(1,1) = 1;      J(1,2) = -EXP(U(6));  J(1,6) = J(1,2)*U(2)
    J(2,1) = VeeA;   J(2,3) = 1
    J(3,2) = VeeB;   J(3,4) = 1
    J(4,3) = 1;      J(4,4) = 1
    J(5,3) = gamC*rho/dT; J(5,5) = 1 + gamC*rho/dT
    J(6,3) = -rho;        J(6,5) = -rho;            J(6,6) = 1

    F(1) = U(1) - EXP(u(6))*U(2)
    F(2) = VeeA*U(1) + U(3) - YUA
    F(3) = VeeB*U(2) + U(4) - YUB
    F(4) = U(3) + U(4)
    F(5) = gamC*rho*U(3)/dT +(1+gamC*rho/dT)*U(5) -gamC*(1 + rho*(GA+GC)/dT)
    F(6) = - rho*U(3) - rho*U(5) + U(6) - pnom
    F = - F ! Newton formula is J.U = -F
    call DEC (6, 6, J, IP, ier)
    if (ier /= 0) STOP "Trouble with DEC, abort"
    call SOL (6, 6, J, F, ip)
    U = U + F  ! F is now the changes vector, is added to U
    if (SQRT(SUM(F**2)) < small) exit
  enddo
  if (nit == nitmax) STOP "No convergence for BC!"
end subroutine BOUND_VALS
