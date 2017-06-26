module PARAMS
  use STUFF;   implicit none
  integer,parameter     :: MaxNewt=4
  integer               :: N, NU, NT
  real(dbl),parameter   :: TOL=1.0E-06_dbl
  real(dbl),allocatable :: U1(:,:), U2(:,:), U3(:,:), deltaU(:,:), &
                           E(:), X(:), RHS(:,:), F(:,:), alpha(:,:), beta(:,:)
  real(dbl)             :: L, dXmin, dT, dP, dC, theta, zR, zC, dE, rho
end module PARAMS

program RPC
! Simulating the reduction of R, charged with zR, to neutral P, with
! only the counter ion C, charged with zC, present, chronoamperometry
! at dimensionless potential theta, reversible redox reaction.
! This is the hypothetical system used in both G.Myla1999 and Bien2004b.
! Solution of the tridiagonal block Jacobian system is by block-Thomas.
  use STUFF;   use PARAMS;   implicit none

  read *, dP, dC
  read *, zR, zC
  read *, theta 
  read *, N, dXmin
  read *, dT, NT
  L = SQRT(MAX(1.0_dbl,dP,dC)) * 6 * SQRT(dT*NT)
  print '(" RPC: R->P at plane electrode,  block Thomas solution.")'
  print '(" L             =", f9.2)', L
  print '(" z_R, z_C      =", 2f7.0)', zR, zC
  print '(" d_P, d_C      =", 3f10.3)', dP, dC
  print '(" theta         =", f10.3)', theta ! As defined by G.Bien2002
  print '(" dT, NT, Tmax  =", f11.4, i6, f10.3)', dT, NT, dT*NT
  print '(" dX(min)       =", es12.1)', dXmin
  print '(" N, NU         =", 2i6)', N, 4*(N+1)
  dE = 1 * dC * (zR-zC) / (zR-zC*dC)
  rho = 1 / SQRT(dE*dP) * (1-zR/zC)
  print '(/" rho           =", f10.3)', rho
  print '(" TOL(norm)     =", es12.1)', TOL
  call RUN ()
end program RPC

subroutine RUN ()
! Plain exp-exp intervals are used; beta are the spatial first derivative
! coefficients, alpha the second-. The ENC is applied.
! The three concentrations for R, P and C plus the potential are
! stored as U, in the order [R, P, C, psi...], and we have past, present
! and unknowns U1, U2, U3, resp., for BDF, started with one BI step.
! The current at Tmax is compared with the Cottrell current, for excess
! electrolyte. The ratios have known values from the cited papers.
! The integrated charge balance is printed at the end and should be zero.
  use STUFF;   use PARAMS;   implicit none
  integer   :: ier, NNewt, Newtmax, i, iT, inext
  real(dbl) :: curr, analcurr, Cottcurr, integ1, integ2, integ3, ratio, &
               G0FORN, U_SIMP

  NU = 4 * (N+1)
  ALLOCATE (U1(4,0:N), U2(4,0:N), U3(4,0:N), E(0:N), deltaU(4,0:N), &
            RHS(4,0:N), F(4,0:N), alpha(3,0:N), beta(3,0:N), X(0:N))

! Initialisation
  U1(1,:) = 1         ! Bulk conc for R
  U1(2,:) = 0         !  "    !   for P
  U1(3,:) = -zC / zR  !  "    !   for C
  U1(4,:) = 0         ! Bulk value of psi
  U2 = U1;   U3 = U1
  call PRECOMP ()
  print '(4x,"iT", 8x, "T", 10x, "curr", 4x, "curr/analcurr", &
          & 4x, "curr/Cott-curr", 4x, "NNewt")'
  Newtmax = 0;  inext = 1
  do iT = 1, NT
    call ONE_STEP ()
    curr = zR * (G0FORN (U3(1,0:4), x(0:4), 5) &
                 + zR * U3(1,0) * G0FORN(U3(4,0:4), x(0:4), 5))
    analcurr = zR * rho / (1 + rho*theta) / SQRT(pi*iT*dT)
    Cottcurr = zR / SQRT(pi*iT*dT) / (1 + theta)
    if (iT == inext .OR. iT == NT) then
      print '(i6, 4f12.3, i16)', &
              iT, iT*dT, curr, curr/analcurr, curr/Cottcurr, NNewt-1
      inext = inext * 2
    endif
    U1 = U2;  U2 = U3
  enddo
  ratio = SQRT(dP) * rho * (1 + theta) / ( 1 + rho*theta)
  print '(/" Expected:", 20x, f12.3, f12.3)', curr/analcurr, ratio
  print '(/" Max Newton steps:", i6)', Newtmax
  integ1 = U_SIMP (X, U3(1,0:N), N)
  integ3 = U_SIMP (X, U3(3,0:N), N)
  print '( " Charge balance at end          :", f10.3)', &
            zR*integ1 + zC*integ3
  DEALLOCATE (U1, U2, U3, F, deltaU, X, alpha, beta)

CONTAINS

subroutine ONE_STEP ()
  call RHS_CALC (iT)
  NNewt = 0
  do
    NNewt = NNewt + 1
    call F_CALC (iT)
    F = F - RHS
    if (SQRT(SUM(F**2)) < TOL) exit
    call SOLVE (iT)
    U3 = U3 + deltaU
    call E_CALC (U3, X, E, N)
    if (SQRT(SUM(deltaU**2)) < TOL) exit
    if (NNewt == maxNewt) exit
    Newtmax = MAX(Newtmax, NNewt)
  enddo
end subroutine ONE_STEP
end subroutine RUN


subroutine PRECOMP ()
! Precomputes the X values, and the derivative coefficients alpha & beta.
  use STUFF;   use PARAMS;   implicit none

  integer   :: i, k
  real(dbl) :: gamma, EE_FAC

! The X values:
  gamma = EE_FAC(dXmin, N, L)
  print '(" Gamma         =", f10.3)', gamma
  X(0) = 0
  do i = 1, N
    x(i) = dXmin * (gamma**i-1) / (gamma-1)
  enddo
! The coeffs:
  call FORN(1, 3, X(0), X(0:2), beta(1:3,0)) ! zero d*/dX at X = 0
  alpha = 0 ! alpha(1:3,0) is not used but must be there for BLOCKS_CALC
  do i = 1, N-1
    call FORN(1, 3, X(i), X(i-1:i+1), beta(1:3,i))
    call FORN(2, 3, X(i), X(i-1:i+1), alpha(1:3,i))
  enddo
end subroutine PRECOMP

subroutine RHS_CALC (iT)
  use STUFF;   use PARAMS;   implicit none
  integer :: iT

  integer :: i

  RHS = 0
  if (iT == 1) then                               ! BI step
    do i = 1, N-1
      RHS(1,i) = - U2(1,i) / dT
      RHS(2,i) = - U2(2,i) / (dP*dT)
      RHS(3,i) = - U2(3,i) / (dC*dT)
      RHS(4,i) = 0
    enddo
  else                                            ! BDF step
    do i = 1, N-1
      RHS(1,i) = (U1(1,i) - 4*U2(1,i)) / (2*dT)
      RHS(2,i) = (U1(2,i) - 4*U2(2,i)) / (2*dP*dT)
      RHS(3,i) = (U1(3,i) - 4*U2(3,i)) / (2*dC*dT)
      RHS(4,i) = 0
    enddo
  endif
  RHS(1,N) = 1
  RHS(2,N) = 0
  RHS(3,N) = -zC / zR
  RHS(4,N) = 0
end subroutine RHS_CALC

subroutine F_CALC (iT)
! Computes the current F
  use STUFF;   use PARAMS;   implicit none
  integer :: iT

  integer   :: i
  real(dbl) :: R_X, P_X, C_X, psi_X, R_XX, P_XX, C_XX, psi_XX

  R_X   = SUM(beta(1:3,0)*U3(1,0:2))
  P_X   = SUM(beta(1:3,0)*U3(2,0:2))
  C_X   = SUM(beta(1:3,0)*U3(3,0:2))
  psi_X = SUM(beta(1:3,0)*U3(4,0:2))
  F(1,0) = U3(1,0) - theta*U3(2,0)             ! Nernst
  F(2,0) = R_X + zR*U3(1,0)*psi_X + dP*P_X     ! R,P, Flux equality
  F(3,0) = C_X + zC*U3(3,0)*psi_X              ! zero C flux
  F(4,0) = zR * U3(1,0)  + zC * U3(3,0)        ! ENC
  do i = 1, N-1
    psi_X  = SUM( beta(1:3,i) * U3(4,i-1:i+1))
    psi_XX = SUM(alpha(1:3,i) * U3(4,i-1:i+1))
!   R:
    R_X  = SUM( beta(1:3,i) * U3(1,i-1:i+1))
    R_XX = SUM(alpha(1:3,i) * U3(1,i-1:i+1))
    F(1,i) =  R_XX - U3(1,i)/dT  +  zR * (U3(1,i)*psi_XX + R_X*psi_X)
      if (iT > 1) F(1,i) = F(1,i) - U3(1,i)/(2*dT)
!   P:
    P_XX = SUM(alpha(1:3,i) * U3(2,i-1:i+1))
    F(2,i) =  P_XX - U3(2,i)/(dP*dT)
      if (iT > 1) F(2,i) = F(2,i) - U3(2,i)/(2*dP*dT)  ! BDF
!   C:
    C_X  = SUM( beta(1:3,i) * U3(3,i-1:i+1))
    C_XX = SUM(alpha(1:3,i) * U3(3,i-1:i+1))
    F(3,i) =  C_XX - U3(3,i)/(dC*dT)  + zC * (U3(3,i)*psi_XX + C_X*psi_X)
      if (iT > 1) F(3,i) = F(3,i) - U3(3,i)/(2*dC*dT)  ! BDF
!   ENC for psi:
    F(4,i) = zR * U3(1,i)  + zC * U3(3,i)
  enddo
! Bulk values, i = N:
  F(1,N) = U3(1,N)
  F(2,N) = U3(2,N)
  F(3,N) = U3(3,N)
  F(4,N) = 0
end subroutine F_CALC

subroutine SOLVE (iT)
! Solves the block system, using block-Thomas for the tridiagonal
! block system.
  use STUFF;   use PARAMS;   implicit none
  integer :: iT

  integer   :: i, row1, col1, sp, ier, IP(12)
  real(dbl) :: Li(4,4), Mi(4,4), Ri(4,4), V(4,4,N), W(4,N), mat(4,4), &
               inv(4,4), vec(4), L0(4,4), M0(4,4), R0(4,4), L1(4,4),  &
               M1(4,4), R1(4,4), B0(4), B1(4)

! Last blocks. We know deltaU_N = 0. This provides P_N-1 and Q_N-1:
  call BLOCKS_CALC (Li, Mi, Ri, N-1, iT) 
  call MAT_INV (Mi(1:4,1:4), 4, inv(1:4,1:4))
  V(1:4,1:4,N-1) = - MATMUL(inv(1:4,1:4), Li(1:4,1:4))
  W(1:4,N-1)     =   MATMUL (inv(1:4,1:4), -F(1:4,N-1))
! All other P and Q up to row zero:
  do i = N-2, 1, -1
    call BLOCKS_CALC (Li, Mi, Ri, i, iT)
    mat(1:4,1:4) = Mi(1:4,1:4) + MATMUL(Ri(1:4,1:4), V(1:4,1:4,i+1))
    call MAT_INV (mat(1:4,1:4), 4, inv(1:4,1:4))
    vec(1:4) = F(1:4,i) + MATMUL(Ri(1:4,1:4), W(1:4,i+1))
    V(1:4,1:4,i) = - MATMUL(inv(1:4,1:4), Li(1:4,1:4))
    W(1:4,i) = - MATMUL(inv(1:4,1:4), vec(1:4))
  enddo
! Now to solve for U(0):
  call BLOCKS_CALC (L0, M0, R0, 0, iT)
  M0 = M0 + MATMUL(R0, V(:,:,2));  B0 = -F(:,0) - MATMUL(R0, W(:,2))
  L0 = L0 + MATMUL(M0, V(:,:,1));  B0 = B0 - MATMUL(M0, W(:,1))
  call MAT_INV (L0, 4, inv)
  deltaU(1:4,0) = MATMUL(inv, B0)
! Complete the solution for all i > 0:
  do i = 1, N-1
    deltaU(1:4,i) = MATMUL(V(:,:,i), deltaU(1:4,i-1)) + W(1:4,i)
  enddo
end subroutine SOLVE

subroutine BLOCKS_CALC (Li, Mi, Ri, i, iT)
! To compute the Jacobian block at index i
  use STUFF;   use PARAMS;   implicit none
  integer   :: i, iT
  real(dbl) :: Li(4,4), Mi(4,4), Ri(4,4)

  real(dbl) :: R_X, C_X, psi_X, psi_XX,     &
               al1, al2, al3, be1, be2, be3

  Li = 0;  Mi = 0;  Ri = 0
  al1 = alpha(1,i);  al2 = alpha(2,i);  al3 = alpha(3,i)
  be1 =  beta(1,i);  be2 =  beta(2,i);  be3 =  beta(3,i)
  if (i == 0) then
!   L0:
    psi_x  = SUM(beta(1:3,0) * U3(4,0:2))
    Li(1,1) = 1;   Li(1,2) = -theta
    Li(2,1) = be1 + zR*psi_X;  Li(2,2) = dP * be1;  Li(2,4) = zR*be1*U3(1,0)
    Li(3,3) = be1 + zC*psi_X;                       Li(3,4) = zC*be1*U3(3,0)
    Li(4,1) = zR;  Li(4,3) = zC
!   C0:
    Mi(2,1) = be2;  Mi(2,2) = dP * be2;  Mi(2,4) = zR * be2 * U3(1,0)
    Mi(3,3) = be2;                       Mi(3,4) = zC * be2 * U3(3,0)
!   R0:
    Ri(2,1) = be3;  Ri(2,2) = dP * be3;  Ri(2,4) = zR * be3 * U3(1,0)
    Ri(3,3) = be3;                       Ri(3,4) = zC * be3 * U3(3,0)
  else
    R_X    = SUM( beta(1:3,i) * U3(1,i-1:i+1))
    C_X    = SUM( beta(1:3,i) * U3(3,i-1:i+1))
    psi_X  = SUM( beta(1:3,i) * U3(4,i-1:i+1))
    psi_XX = SUM(alpha(1:3,i) * U3(4,i-1:i+1))
!   Left block:
    Li(1,1) = al1 + zR*be1*psi_X;  Li(1,4) = zR*(al1*U3(1,i) + be1*R_X)
    Li(2,2) = al1
    Li(3,3) = al1 + zC*be1*psi_X;  Li(3,4) = zC*(al1*U3(3,i) + be1*C_X)
!   Center block:
    Mi(1,1) = al2 - 1/dT + zR*(psi_XX + be2*psi_X)
      if (iT > 1) Mi(1,1) = Mi(1,1) - 1/(2*dT)       ! BDF
      Mi(1,4) = zR*(al2*U3(1,i) + be2*R_X)
    Mi(2,2) = al2 - 1/(dP*dT)
      if (iT > 1) Mi(2,2) = Mi(2,2) - 1/(2*dP*dT)
    Mi(3,3) = al2 - 1/(dC*dT)  + zC*(psi_XX + be2*psi_X)
      if (iT > 1) Mi(3,3) = Mi(3,3) - 1/(2*dC*dT)    ! BDF
      Mi(3,4) = zC*(al2*U3(3,i) + be2*C_X)
    Mi(4,1)  = zR;   Mi(4,3) = zC
!   Right block:
    Ri(1,1) = al3 + zR*be3*psi_X;  Ri(1,4) = zR * (al3*U3(1,i) + be3*R_X)
    Ri(2,2) = al3
    Ri(3,3) = al3 + zC*be3*psi_X;  Ri(3,4) = zC * (al3*U3(3,i) + be3*C_X)
  endif
! No block at i = N needed, this is done locally in SOLVE.
end subroutine BLOCKS_CALC

subroutine E_CALC (U, X, E, N)
! Computes the potential field by differentiating the potentials.
! 5-point approximations are used.
  use STUFF;   implicit none
  integer   :: N
  real(dbl) :: U(4,0:N), X(0:N), E(0:N)

  integer   :: i
  real(dbl) :: w(5)

  call FORNBERG (1, X(0), X(0:4), U(4,0:4), 5, E(0), w(1:5)) ! i = 0
  call FORNBERG (1, X(1), X(0:4), U(4,0:4), 5, E(1), w(1:5)) ! i = 1
  do i = 2, N-2
    call FORNBERG (1, X(i), X(i-2:i+2), U(4,i-2:i+2), 5, E(i), w(1:5))
  enddo
  call FORNBERG (1, X(N-1), X(N-4:N), U(4,N-4:N), 5, E(N-1), w(1:5)) ! N-1
  call FORNBERG (1, X(N),   X(N-4:N), U(4,N-4:N), 5, E(N),   w(1:5)) ! N
end subroutine E_CALC
