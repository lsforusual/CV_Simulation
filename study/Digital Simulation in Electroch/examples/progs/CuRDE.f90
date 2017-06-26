module PARAMS
  use STUFF;   implicit none 
  integer,parameter     :: MaxNewt=30
  integer               :: N, NU
  real(dbl),parameter   :: TOL=1.0E-07, Kconst = 2.13593_dbl
  real(dbl),allocatable :: U(:,:), E(:), X(:), F(:,:), deltaU(:,:),  &
                           alpha(:,:), beta(:,:), d(:), z(:), Kd(:,:)
  real(dbl)             :: dXmin, L, S
end module PARAMS

program Block_Cudep_RDE
! Simulates Cu deposition at an RDE in the presence of H2SO4 at
! S times the concentration of CuSO4.
! The steady state  current should approach (normalised) unity for
! large concentrations of H2SO4 and 1.885 for zero H2SO4.
! Block-tri is used, with ENC even at X = 0.
  use STUFF;   use PARAMS;   implicit none

  ALLOCATE (d(3), z(3))
  call FILSPC (4)
  d(1) = 1
  read *, d(2), d(3), S
  read *, N, dXmin
  print '(" ----St-st. Cu Dep. at RDE. June 2015")'
  print '(" A, B, C = [Cu++], [H+], [SO4--], resp.")'
  print '(" S             =", f11.4)', S
  print '(" Init. A, B, C =", 3f10.3)', 1.0, 2*S, 1+S
  L = 5 * SQRT(MAXVAL(d))
  print '(" L             =", f9.2)', L
  print '(" d_A, d_B, d_C =", 3f10.3)', d(1:3)
  print '(" H2SO4 conc    =", f10.3)', S
  print '(" dX(min)       =", es12.1)', dXmin
  print '(" N, NU         =", 2i6)', N, 4*(N+1)
  print '(" TOL(norm)     =", es12.1)', TOL
  print '(" Newman''s sqrt(r) =", f10.3)', SQRT(S/(1+S))
  call CPUNUL
  call RUN ()
  call CPUOUT ("Final")
end program Block_Cudep_RDE


subroutine RUN ()
! Exp-exp intervals are used; beta are  the spatial first derivative
! coefficients, alpha the second-. ENC is applied.
! The three concentrations for Cu++, H+ and SO4-- are stored as
! U(i,:), i = 1, 2, 3, respectively, and psi as U(4,*). The relative
! diffusion coeffs are d(1:3) and the charges z(1:3). The parameter
! K / d * X^2 is stored as Kd(3,1:N). These arrays make the code for
! F and Jac more compact and clearer.
  use STUFF;   use PARAMS;   implicit none

  integer   :: ier, NNewt, Newtmax, i
  real(dbl) :: curr, integ(3), G0FORN, USIMP

  NU = 4 * (N+1)
  ALLOCATE (U(4,0:N), deltaU(4,0:N), E(0:N), F(4,0:N), &
            alpha(3,N), beta(3,0:N), X(0:N), Kd(3,N))

  U(1,:) = 1;  U(1,0) = 0;   z(1) =  2   ! Cu++
  U(2,:) = 2*S;              z(2) =  1   ! H+
  U(3,:) = 1 + S;            z(3) = -2   ! SO4--
  U(4,:) = 0                             ! psi
  call PRECOMP ()
  print *
  Newtmax = 0
  call STEADY ()
  call  E_CALC (U, X, E, N)
  curr = G0FORN (U(1,0:4), x(0:4), 5)
  print '(/" Current =", f10.4)', curr
  print '(/" Max Newton steps:", i6)', Newtmax
  write (4,'("# Profiles at steady state, S =", f8.3)') S
  write (4,'(8x,"X",12x,"Cu++",11x,"H+",10X,"SO4--",11x,"Psi", "sqrt(r)")')
  do i = 0, N
    write (4,'(6es15.6)') X(i), U(1:4,i), SQRT(S/(1+S))
  enddo
  integ(1) = USIMP (X, U(1,0:N), N)
  integ(2) = USIMP (X, U(2,0:N), N)
  integ(3) = USIMP (X, U(3,0:N), N)
  print '(" Final charge balance:", f12.5)', SUM(z(1:3)*integ(1:3))
  DEALLOCATE (U, deltaU, F, X, alpha, beta, Kd, d, z)

CONTAINS
subroutine STEADY ()
integer :: i

  NNewt = 0
  do
    NNewt = NNewt + 1
    call F_CALC ()
    call SOLVE ()
    U = U + deltaU
    if (SQRT(SUM(deltaU**2)) < TOL) exit
    if (NNewt == maxNewt) exit
    Newtmax = MAX(Newtmax, NNewt)
  enddo
end subroutine STEADY
end subroutine RUN

subroutine PRECOMP ()
! Precomputes the X values, and the derivative coefficients alpha & beta.
  use STUFF;   use PARAMS;   implicit none

  integer   :: i, k
  real(dbl) :: gamma, EE_FAC, deriv

! The X values, and K' at same time:
  gamma = EE_FAC (dXmin, N, L)
  print '(" Gamma         =", f10.4)', gamma
  X(0) = 0
  do i = 1, N
    X(i) = dXmin * (gamma**i-1)/(gamma-1)
    do k = 1, 3
      Kd(k,i) = Kconst * X(i)**2 / d(k)
    enddo
  enddo
! The coeffs:
  call FORN(1, 3, X(0), X(0:2), beta(1:3,0)) ! zero d*/dX at X = 0
  do i = 1, N-1
    call FORN(1, 3, X(i), X(i-1:i+1), beta(1:3,i))
    call FORN(2, 3, X(i), X(i-1:i+1), alpha(1:3,i))
  enddo
  call FORN(1, 3, X(N), X(N-2:N), beta(1:3,N)) ! zero d*/dX at X = L
end subroutine PRECOMP


subroutine F_CALC ()
! Computes the current F, RHS subtracted (mostly zeroes).
  use STUFF;   use PARAMS;   implicit none

  integer   :: i, sp
  real(dbl) :: U_X,  U_XX, psi_X, psi_XX

  F = 0 ! Default

  psi_X  = SUM(beta(1:3,0)  * U(4,0:2))
  F(2,0) = SUM (beta(1:3,0) * U(2,0:2)) + z(2)*U(2,0)*psi_X
  F(3,0) = SUM (beta(1:3,0) * U(3,0:2)) + z(3)*U(3,0)*psi_X
  F(4,0) = SUM(z(1:3) * U(1:3,0))
  do i = 1, N-1
    psi_X  = SUM(beta(1:3,i)  * U(4,i-1:i+1))
    psi_XX = SUM(alpha(1:3,i) * U(4,i-1:i+1))
    do sp = 1, 3 ! Three ionic species
      U_X  = SUM(beta(1:3,i)  * U(sp,i-1:i+1))
      U_XX = SUM(alpha(1:3,i) * U(sp,i-1:i+1))
      F(sp,i) = U_XX + Kd(sp,i)*U_X + z(sp)*(U(sp,i)*psi_XX + U_X*psi_X)
    enddo
    F(4,i) = SUM(z(1:3)*U(1:3,i))
  enddo
  F(1:4,N) = 0              ! Bulk values
end subroutine F_CALC

subroutine SOLVE ()
! Solves the block system, using block-Thomas for the tridiagonal
! block system of 4 x 4 blocks.
  use STUFF;   use PARAMS;   implicit none
  integer :: iT

  integer   :: i, row1, col1, sp, ier, IP(12)
  real(dbl) :: L0(4,4), M0(4,4), R0(4,4), Li(4,4), Mi(4,4), Ri(4,4), &
               V(4,4,N), W(4,N), mat(4,4), inv(4,4), vec(4)

! Blocks at N-1. We know deltaU_N = 0. This provides P_N-1 and Q_N-1:
  call BLOCKS_CALC (Li, Mi, Ri, N-1) 
  call MAT_INV (Mi(1:4,1:4), 4, inv(1:4,1:4))
  V(1:4,1:4,N-1) = - MATMUL(inv(1:4,1:4), Li(1:4,1:4))
  W(1:4,N-1)     =   MATMUL (inv(1:4,1:4), -F(1:4,N-1))
! All other V and W up to row zero:
  do i = N-2, 1, -1
    call BLOCKS_CALC (Li, Mi, Ri, i)
    mat(1:4,1:4) = Mi(1:4,1:4) + MATMUL(Ri(1:4,1:4), V(1:4,1:4,i+1))
    call MAT_INV (mat(1:4,1:4), 4, inv(1:4,1:4))
    vec(1:4) = F(1:4,i) + MATMUL(Ri(1:4,1:4), W(1:4,i+1))
    V(1:4,1:4,i) = - MATMUL(inv(1:4,1:4), Li(1:4,1:4))
    W(1:4,i) = - MATMUL(inv(1:4,1:4), vec(1:4))
  enddo
! Now to solve for U(0), using substitutions for deltaU(:,1)
! and deltaU(:,2), i.e. V(:,:,1) & W(:,1) and V(:,:,2) & W(:,2).
  call BLOCKS_CALC (L0, M0, R0, 0)
  mat(:,:) =  L0(:,:)                                       &
             + MATMUL (M0(:,:), V(:,:,1))                   &
             + MATMUL (R0(:,:), MATMUL(V(:,:,2), V(:,:,1)))
  vec(:) = - F(:,0)                                     &
           - MATMUL (M0(:,:), W(:,1))                   &
           - MATMUL (R0(:,:), MATMUL(V(:,:,2), W(:,1))) &
           - MATMUL (R0(:,:), W(:,2))
  call MAT_INV (mat, 4, inv)
  deltaU(:,0) = MATMUL(inv, vec)
! Complete the solution for all i > 0:
  do i = 1, N-1
    deltaU(:,i) = MATMUL(V(:,:,i), deltaU(:,i-1)) + W(:,i)
  enddo
end subroutine SOLVE


subroutine BLOCKS_CALC (Li, Mi, Ri, i)
! To compute the Jacobian block at index i
  use STUFF;   use PARAMS;   implicit none
  integer   :: i
  real(dbl) :: Li(4,4), Mi(4,4), Ri(4,4)

  integer   :: sp
  real(dbl) :: psi_X, psi_XX, U_X, al1, al2, al3, be1, be2, be3

  Li = 0;  Mi = 0;  Ri = 0
  if (i == 0) then
    psi_X = SUM(beta(1:3,0)*U(4,0:2))
    Li(1,1) = 1
    Li(2,2) = beta(1,0) + z(2)*psi_X;  Li(2,4) = beta(1,0) * z(2) * U(2,0)
    Li(3,3) = beta(1,0) + z(3)*psi_X;  Li(3,4) = beta(1,0) * z(3) * U(3,0)
    Li(4,1:3) = z(1:3)
    Mi(2,2) = beta(2,0);  Mi(2,4) = beta(2,0) * z(2) * U(2,0)
    Mi(3,3) = beta(2,0);  Mi(3,4) = beta(2,0) * z(3) * U(3,0)
    Ri(2,2) = beta(3,0);  Ri(2,4) = beta(3,0) * z(2) * U(2,0)
    Ri(3,3) = beta(3,0);  Ri(3,4) = beta(3,0) * z(3) * U(3,0)
  else if (i < N) then
    psi_X  = SUM(beta(1:3,i)  * U(4,i-1:i+1))
    psi_XX = SUM(alpha(1:3,i) * U(4,i-1:i+1))
    al1 = alpha(1,i);  al2 = alpha(2,i);  al3 = alpha(3,i)
    be1 = beta(1,i);   be2 = beta(2,i);   be3 = beta(3,i)
!   Left block:
    do sp = 1, 3
      U_X  = SUM(beta(1:3,i) * U(sp,i-1:i+1))
      Li(sp,sp) = al1 + Kd(sp,i)*be1 + z(sp)*be1*psi_X
      Li(sp,4) = z(sp) * (al1*U(sp,i) + be1*U_X)
    enddo
!   Centre block:
    do sp = 1, 3
      U_X   = SUM(beta(1:3,i) * U(sp,i-1:i+1))
      Mi(sp,sp) = al2 + Kd(sp,i)*be2 + z(sp)*(psi_XX + be2*psi_X)
      Mi(sp,4) = z(sp) * (al2*U(sp,i) + be2*U_X)
    enddo
    Mi(4,1:3) = z(1:3)
!   Right block:
    do sp = 1, 3
      U_X   = SUM(beta(1:3,i) * U(sp,i-1:i+1))
      Ri(sp,sp) = al3 + Kd(sp,i)*be3 + z(sp)*be3*psi_X
      Ri(sp,4) = z(sp) * (al3*U(sp,i) + be3*U_X)
    enddo
  else ! (i = N)
    STOP " No blocks for i = N, abort."
  endif
end subroutine BLOCKS_CALC


subroutine E_CALC (U, X, E, N)
! Computes the potential field by differentiating the potentials.
! 5-point approximations are used, and FORNBERG.
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
