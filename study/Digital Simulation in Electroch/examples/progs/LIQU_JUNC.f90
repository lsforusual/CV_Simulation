module PARAMS
  use STUFF;   implicit none
  integer,parameter     :: MaxNewt=1999
  integer               :: M, NT, N
  real(dbl),parameter   :: TOL=1.0E-05_dbl
  real(dbl)             :: CL, CR, L, d_A, d_B, dT, dXmin, f_A, f_B
end module PARAMS

program LIQU_JUNC
! To simulate a double-ended chamber with two different concentrations
! CL and CR of A+B-, meeting at a sharp interface at x = 0, t = 0.
! We input M, the no. of points in each chamber, so that the total
! number of points between -L <= X <= +L is 2M+1. The end points have fixed
! values for A and B and psi at the LH end (zero) but at the RH end its 1st
! derivative is zero. The ends are all included as unknowns for convenience.
! Profiles are written out at the end, T=Tmax.
  use STUFF;   use PARAMS;   implicit none
  real(dbl) :: Hend

  call FILSPC (4)
  read *, CL, CR
  read *, d_A, d_B
  read *, dT, NT
  read *, M, dXmin, L
  print '(/" Liquid junuction block-diag. ABE. BI/BDF3")'
  print '(" CL, CR   =", 2f10.3)', CL, CR
  print '(" d_A, d_B =", 2f11.4)', d_A, d_B
  print '(" dT       =", f10.3)', dT
  print '(" NT       =", i6)', NT
  print '(" Tmax     =", f9.2)', NT * dT
  print '(" L        =", f9.2)', L
  print '(" dX(min)  =", es12.1)', dXmin
  print '(" M, N     =", 3i6)', M, 2*M+1
  print '(" TOL(norm)=", es12.1)', TOL
  Hend = -(d_A-d_B)/(d_A+d_B) * LOG(CR/CL)
  f_A = d_A * dT;   f_B = d_B * dT
  print '(" Henderson:", f11.4, " units")', Hend
  call RUN ()
end program LIQU_JUNC

subroutine RUN ()
  use STUFF;   use PARAMS;   implicit none
  character(len=30)     :: filename
  integer               :: iT, ier, nextout, Newtmax, maxNewtmax
  real(dbl),allocatable :: A1(:), B1(:), E1(:), A2(:), B2(:), E2(:),  &
                           A3(:), B3(:), E3(:), psi(:), X(:), U(:,:), &
                           F(:,:), RHS(:,:), alpha(:,:), beta(:,:),   &
                           v(:,:,:), w(:,:), deltaU(:,:)
  real(dbl)             :: Li(3,3), Ci(3,3), Ri(3,3)

  N = 2*M+1;   ! Known end points are included as "unknowns"
  ALLOCATE (A1(N), B1(N), E1(N), A2(N), B2(N), E2(N), A3(N), B3(N), E3(N), &
            psi(N), U(3,N), F(3,N), RHS(3,N), v(3,3,N), w(3,N),            &
            deltaU(3,N), alpha(3,N), beta(3,N), X(N))

  A1(1:M) = CL;   A1(M+1) = (Cl+CR)/2;   A1(M+2:N) = CR ! Initial values
  B1(1:M) = CL;   B1(M+1) = (Cl+CR)/2;   B1(M+2:N) = CR
  A2 = A1; B2 = B1; A3 = A1;  B3 = B1
  E1 = 0;  E2 = E1;  E3 = E1
  call PRECOMP (X, beta, alpha, M, N, L, dXmin) ! X, alpha & beta arrays.
  nextout = 1
  maxNewtmax = 0
  print '(4x,"iT", 8x,"T", 10x, "psi(N)", 5x, "NNewt")'
  do iT = 1, NT
    call ONE_STEP ()
    call MAYBE_OUTPUT ()
    A1 = A2;  B1 = B2;  E1 = E2
    A2 = A3;  B2 = B3;  E2 = E3
    maxNewtmax = MAX (maxNewtmax, Newtmax)
  enddo
  print '(" Max Newtons needed:", i4)', maxNewtmax
  DEALLOCATE (A1, B1, E1, A2, B2, E2, A3, B3, E3, psi, U, F, RHS, &
            v, w, deltaU, alpha, beta, X)

CONTAINS
  subroutine ONE_STEP ()
    use STUFF;   use PARAMS;   implicit none

    integer :: NNewt

    call RHS_CALC (A1, B1, A2, B2, RHS, iT)
    NNewt = 0;  Newtmax = 0
    do
      NNewt = NNewt + 1
      call ABE_INTO_U (A3, B3, E3, U, N)
      call F_CALC (A3, B3, E3, alpha, beta, F, iT)
      F = F - RHS
      if (SQRT(SUM(F**2)) < TOL) exit
      call THOMAS (A3, B3, E3, F, alpha, beta, deltaU, iT)
      U = U + deltaU
      if (SQRT(SUM(deltaU**2)) < TOL) exit
      if (NNewt == maxNewt) exit
      Newtmax = MAX(Newtmax, NNewt)
      call U_INTO_ABE (U, A3, B3, E3, N)
    enddo
  end subroutine ONE_STEP

  subroutine MAYBE_OUTPUT ()
    integer :: i

    if (iT == NT) then
      call E_INTEG (X, E3, psi, N)
      write (4,'("# iT, T =", i6, f10.3)') iT, dT*iT
      do i = 1, N
        write (4,'(6es15.6)') &
                   X(i), A3(i), B3(i), E3(i), psi(i), A3(i)-B3(i)
      enddo
      CLOSE (4)
    endif
    if (iT == nextout .OR. iT == NT) then
      call E_INTEG (X, E3, psi, N)
      print '(i8, 2f12.6, i8)', iT, iT*dT, psi(N), Newtmax
      nextout = 2 * nextout
    endif
  end subroutine MAYBE_OUTPUT
end subroutine RUN

subroutine ABE_INTO_U (A, B, E, U, N)
! Deals the three arrays A, B, E, into U in the order
! A_1, B_1, E_1, A_2, B_2, E_2, ...
  use STUFF;   implicit none
  integer   :: N
  real(dbl) :: A(N), B(N), E(N), U(3,N)

  integer   :: i

  do i = 1, N
    U(1,i) = A(i)
    U(2,i) = B(i)
    U(3,i) = E(i)
  enddo
end subroutine ABE_INTO_U

subroutine U_INTO_ABE (U, A, B, E, N)
! Picks the three arrays A, B, E, out of U.
  use STUFF;   implicit none
  integer   :: N
  real(dbl) :: U(3,N), A(N), B(N), E(N)

  integer   :: i

  do i = 1, N
    A(i) = U(1,i)
    B(i) = U(2,i)
    E(i) = U(3,i)
  enddo
end subroutine U_INTO_ABE

subroutine RHS_CALC (A1, B1, A2, B2, RHS, iT)
! Computes the current RHS.
  use STUFF;   use PARAMS;   implicit none
  integer   :: iT
  real(dbl) :: A1(N), B1(N), A2(N), B2(N), RHS(3,N)

  integer   :: i

  RHS(1:2,1) = CL;  RHS(3,1) = 0
  if (iT == 1) then  ! First BI step
    do i = 2, N-1
      RHS(1,i) = -A2(i) / f_A
      RHS(2,i) = -B2(i) / f_B
      RHS(3,i) = 0
    enddo
  else               ! Subsequent BDF steps
    do i = 2, N-1
      RHS(1,i) = (A1(i) -4*A2(i)) / (2*f_A)
      RHS(2,i) = (B1(i) -4*B2(i)) / (2*f_B)
      RHS(3,i) = 0
    enddo
  endif
  RHS(1:2,N) = CR;  RHS(3,N) = 0
end subroutine RHS_CALC

subroutine F_CALC (A, B, E, alpha, beta, F, iT)
! Computes the current F function values.
  use STUFF;   use PARAMS;   implicit none
  integer   :: iT
  real(dbl) :: A(N), B(N), E(N), alpha(3,N), beta(3,N), F(3,N)

  integer   :: i
  real(dbl) :: DelA, Del2A, DelB, Del2B, DelE

  F(1:2,1) = CL;  F(3,1) = 0
  do i = 2, N-1
    DelA  = SUM(beta (1:3,i) * A(i-1:i+1))
    Del2A = SUM(alpha(1:3,i) * A(i-1:i+1))
    DelB  = SUM(beta (1:3,i) * B(i-1:i+1))
    Del2B = SUM(alpha(1:3,i) * B(i-1:i+1))
    DelE  = SUM(beta (1:3,i) * E(i-1:i+1))
    F(1,i) = Del2A - A(i)/f_A - DelA*E(i) - A(i)*DelE ! For A
      if (iT > 1) F(1,i) = F(1,i) - 0.5_dbl*A(i)/f_A
    F(2,i) = Del2B - B(i)/f_B + DelB*E(i) + B(i)*DelE ! For B
      if (iT > 1) F(2,i) = F(2,i) - 0.5_dbl*B(i)/f_B
    F(3,i) = A(i) - B(i) - DelE                       ! For E
  enddo
  F(1,N) = A(N);   F(2,N) = B(N);  F(3,N) = 0
end subroutine F_CALC

subroutine THOMAS (A, B, E, F, alpha, beta, deltaU, iT)
  use STUFF;   use PARAMS;   implicit none
  integer   :: iT
  real(dbl) :: A(N), B(N), E(N), F(3,N), deltaU(3,N), &
               alpha(3,N), beta(3,N)

  integer   :: i
  real(dbl) :: Li(3,3), Ci(3,3), Ri(3,3), mat(3,3), inv(3,3), vec(3), &
               v(3,3,N), w(3,N)

 ! Start; deltaU(1) = 0.
  call BLOCKS_CALC (A, B, E, alpha, beta, Li, Ci, Ri, iT, 2)
  call MAT_INV (Ci, 3, inv)
  v(:,:,2) = MATMUL (inv(:,:), -Ri(:,:))
  w(:,2)   = MATMUL (inv(:,:), -F(:,2))
  do i = 3, N-1
    call BLOCKS_CALC (A, B, E, alpha, beta, Li, Ci, Ri, iT, i)
    mat(1:3,1:3) = MATMUL(Li(:,:), v(:,:,i-1)) + Ci(:,:)
    call MAT_INV (mat, 3, inv)
    v(:,:,i) = MATMUL (inv(:,:), -Ri(:,:))
    vec(:)   = MATMUL (Li(:,:), w(:,i-1)) + F(:,i)
    w(:,i)   = MATMUL (inv(:,:), -vec(:))
  enddo
! Solving for delta U, by mat * delta U = vec
  deltaU(:,N) = 0
  do i = N-1, 2, -1
    deltaU(:,i) = MATMUL (v(:,:,i), deltaU(:,i+1)) + w(:,i)
  enddo
end subroutine THOMAS

subroutine BLOCKS_CALC (A, B, E, alpha, beta, Li, Ci, Ri, iT, i)
! Computes the current 3x3 blocks at index i.
  use STUFF;   use PARAMS;   implicit none
  integer   :: iT, i
  real(dbl) :: A(N), B(N), E(N), alpha(3,N), beta(3,N), &
               Li(3,3), Ci(3,3), Ri(3,3)

  real(dbl) :: DelA, DelB, DelE

  DelA    = SUM(beta(1:3,i)  * A(i-1:i+1))
  DelB    = SUM(beta(1:3,i)  * B(i-1:i+1))
  DelE  = SUM(beta(1:3,i)  * E(i-1:i+1))
  ! Left block for i-1:
  Li(1,1) = alpha(1,i) - beta(1,i)*E(i)            ! A row
  Li(1,2) = 0
  Li(1,3) = -beta(1,i)*A(i)
  Li(2,1) = 0                                      ! B row
  Li(2,2) = alpha(1,i) + beta(1,i)*E(i)
  Li(2,3) =  beta(1,i)*B(i)
  Li(3,1) = 0;  Li(3,2) = 0;  Li(3,3) = -beta(1,i) ! Whole E row
  ! Centre block, for i:
  Ci(1,1) = alpha(2,i) - 1/f_A - beta(2,i)*E(i) - delE ! A row
    if (iT > 1) Ci(1,1) = Ci(1,1) - 0.5_dbl/f_A        ! BDF
  Ci(1,2) = 0
  Ci(1,3) = - delA - beta(2,i)*A(i)
  Ci(2,1) = 0                                        ! B row
  Ci(2,2) = alpha(2,i) - 1/f_B + beta(2,i)*E(i) + delE
    if (iT > 1) Ci(2,2) = Ci(2,2) - 0.5_dbl/f_B           ! BDF
  Ci(2,3) = delB + beta(2,i)*B(i) 
  Ci(3,1) = 1;  Ci(3,2) = -1;   Ci(3,3) = -beta(2,i) ! Whole E row
  ! Right block for i+1:
  Ri(1,1) = alpha(3,i) - beta(3,i)*E(i)              ! A row
  Ri(1,2) = 0
  Ri(1,3) = -beta(3,i)*A(i)
  Ri(2,1) = 0                                        ! B row
  Ri(2,2) = alpha(3,i) + beta(3,i)*E(i)
  Ri(2,3) = beta(3,i)*B(i)
  Ri(3,1) = 0;  Ri(3,2) = 0;  Ri(3,3) = -beta(3,i)   ! Whole E row
end subroutine BLOCKS_CALC

subroutine PRECOMP (X, beta, alpha, M, N, L, dXmin)
! Precomputes the X values, and the derivative coefficients alpha & beta.
  use STUFF;   implicit none
  integer   :: M, N
  real(dbl) :: X(N), beta(3,N), alpha(3,N), L, dXmin

  integer   :: i
  real(dbl) :: gamma
  real(dbl),allocatable :: Y(:)

  ALLOCATE (Y(0:M))
! The X values:
  call SIGMOID_EXPANSION (dXmin, M, L, gamma, Y)
  print '(" Gamma      =", f10.3)', 1+gamma
  print '(" Largest interval in X =", f10.3)', Y(M) - Y(M-1)
  X(M+1) = 0
  do i = 1, M
    X(M+1+i) = Y(i)
    X(M+1-i) = - Y(i)
  enddo
  X(1) = -L;  X(N) = L ! Just to make sure.
  DEALLOCATE (Y)
! The coeffs:
  do i = 2, N-1
    call FORN(1, 3, X(i), X(i-1:i+1), beta(1:3,i))
    call FORN(2, 3, X(i), X(i-1:i+1), alpha(1:3,i))
  enddo
  call FORN(1, 3, X(N), X(N-2:N), beta(1:3,N)) ! RH pt for E, 1st deriv = 0
end subroutine PRECOMP

subroutine E_INTEG (X, E, psi, N)
! Integrates the E profile to produce psi, the potential, relative
! to the LH end of the scale, which is arbitrarily set to zero.
! Simple trapezium integration is used.
  use STUFF;   implicit none
  integer   :: N
  real(dbl) :: X(N), E(N), psi(N)

  integer   :: i
  real(kind=dbl) :: USIMP

  psi(1) = 0
  psi(2) = (E(2) + E(1)) / 2 * (X(2) - X(1)) ! Trapeze here
  do i = 3, N
    psi(i) = USIMP (X, E, i)
  enddo
  psi = -psi
end subroutine E_INTEG
