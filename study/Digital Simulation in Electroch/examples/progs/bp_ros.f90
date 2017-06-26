program BP_ROS
! Rosenbrock method used on a Birk-Perone simulation.
! The Rosenbrock orders to choose are, i.e.
!  ORDER = 2: 2nd order ROS2 as in Lang's book, App. C
!  ORDER = 3: 3rd order ROWDA3 as in Roch1988 or Bien1999.
! Equal intervals are used, with lambda as input.
  use STUFF;   use ROSTUFF;   implicit none
  integer        :: order, N, np, nT, inext, i, iT
  real(kind=dbl) :: lambda, K, H, dT, T, G, G_anal,logerr,   &
                    bpa(10), G0FUNC, GBPFUNC
  real(kind=dbl),allocatable :: C(:), F(:), B(:), Sel(:),    &
                                k1(:),k2(:),k3(:), ad(:),bd(:)
! The variables:
! F receives dT*F(C+ ...)
! B is the RHS at any stage
! Sel is the selection diagonal, [0,1...1,0], needed for B
! ad, bd are the a' and b' coeffs for the solution.

  print '(" Orders:  2 = ROS2;  3 = ROWDA3")'
  print '(" nT, lambda, K, np, order?")'
  read *, nT, lambda, K, np, order
  call FILSPC (4)
  dT = 1.0_dbl / nT;   H = SQRT (dT/lambda)
  N = NINT (6 / H)
  ALLOCATE (C(0:N+1), F(0:N+1), B(0:N+1), Sel(0:N+1),    &
            k1(0:N+1), k2(0:N+1), k3(0:N+1), ad(N), bd(N))
  call ROCOEFFS (order)
  bpa(1) = 4/pi - 1 ! The coefficients for the series solution:
  bpa(2:10) = (/ 0.08327_dbl, 0.02893_dbl, 0.01162_dbl, &
                 0.00540_dbl, 0.00286_dbl, 0.00169_dbl, &
                 0.00108_dbl, 0.00074_dbl, 0.00053_dbl /)
  Sel = 1;  Sel(0) = 0;  Sel(N+1) = 0 ! Selection diagonal set
  C = 1;  C(0) = 0
  print '(" B & P by Rosenbrock order", i2, ".")', order
  print '(" nT       =", i6)', nT
  print '(" lambda   =", f8.2)', lambda
  print '(" K        =", f8.2)', K
  print '(" H        =", f10.4)', H
  print '(" N        =", i6)', N
  print '(" np(curr) =", i6)', np
  print *
  print '("   iT", 5x, "T",  9x, "G", 4x, "log|err|")'
  inext = 1
  do iT = 1, nT
      T = (iT-1) * dT
      call STEP (C, N, K, dT, H, lambda, T, F, B, Sel, &
                          ad, bd, k1, k2, k3, order)
      if (iT == inext .OR. iT == nT) then
      T = iT * dT
      G = G0FUNC(C, np, H)
      G_anal = GBPFUNC (K, T, bpa)
      logerr = LOG10(ABS(G - G_anal)/G_anal)
      print '(i6, f8.3, f10.3, f8.2)', iT, T, G, logerr
      write (4,'(3f10.6)') iT*dT, G, logerr
      inext = 2 * inext
    endif
  enddo
end program BP_ROS


subroutine STEP (C, N, K, dT, H, lambda, T, F, B, Sel, &
                            ad, bd, k1, k2, k3, order)
! Takes one Rosenbrock step.
  use STUFF;   use ROSTUFF;   implicit none

  integer        :: N, order
  real(kind=dbl) :: C(0:N+1), K, T, dT, H, lambda, F(0:N+1), &
                    B(0:N+1), Sel(0:N+1), ad(N), bd(N),      &
                    k1(0:N+1), k2(0:N+1), k3(0:N+1)
  integer        :: i

  select case (order)
    case (2) ! ROS2
      call EVAL_B (K,dT,H,lambda, T, 1, k1,k2, C, F, Sel, N, B)
      call SOLVE (K, H, lambda, N, ad, bd, B, k1)      ! --> k1
      call EVAL_B (K,dT,H,lambda, T, 2, k1,k2, C, F, Sel, N, B)
      call SOLVE (K, H, lambda, N, ad, bd, B, k2)      ! --> k2
      C(0:N+1) = C(0:N+1) + m1*k1(0:N+1) + m2*k2(0:N+1)
    case (3) ! ROWDA3
      call EVAL_B (K,dT,H,lambda, T, 1, k1,k2, C, F, Sel, N, B)
      call SOLVE (K, H, lambda, N, ad, bd, B, k1)      ! --> k1
      call EVAL_B (K,dT,H,lambda, T, 2, k1,k2, C, F, Sel, N, B)
      call SOLVE (K, H, lambda, N, ad, bd, B, k2)      ! --> k2
      call EVAL_B (K,dT,H,lambda, T, 3, k1,k2, C, F, Sel, N, B)
      call SOLVE (K, H, lambda, N, ad, bd, B, k3)      ! --> k3
      C(0:N+1) = C(0:N+1) &
                 + m1*k1(0:N+1) + m2*k2(0:N+1) + m3*k3(0:N+1)
    case default; STOP ' Disallowed Rosenbrock order.'
  end select
end subroutine STEP


subroutine EVAL_B (K, dT, H, lambda, T, s, k1, k2, C, &
                                                 F, Sel, N, B)
! Evaluates the current B vector.
  use STUFF;   use ROSTUFF;   implicit none
  integer        :: s, N
  real(kind=dbl) :: K, dT,H,lambda, T, k1(0:N+1), k2(0:N+1),   &
                    C(0:N+1), B(0:N+1), F(0:N+1), Sel(0:N+1)

  real(kind=dbl) :: Tee, Ft

! First augment C (according to the stage s), temporarily
! storing it in B:
  select case (s)
    case (1);  B(0:N+1) = C(0:N+1);  Tee = T
    case (2);  B(0:N+1) = C(0:N+1) + a21*k1(0:N+1)
               Tee = t + alpha2 * dT
    case (3);  B(0:N+1) = C(0:N+1) + a31*k1(0:N+1) &
                                   + a32*k2(0:N+1)
               Tee = T + alpha3*dT
  end select

  call EVAL_F (K, H, Tee, B, F, N)         ! -dT*F/lambda
  Ft = dT * H**2* K/(1+K*T)**2 /lambda   ! Ft/lambda

  select case (s) ! Now make B, adding some other bits
    case (1)
        B(0:N) = F(0:N)
        B(N+1) = F(N+1) - gamma1*Ft
    case (2)
         B(0:N) = F(0:N) - c21*Sel(0:N+1)*k1(0:N+1)/lambda
         B(N+1) = F(N+1) - gamma2*Ft
    case (3);  B(0:N) = F(0:N) &
                        - c31*Sel(0:N+1)*k1(0:N+1)/lambda &
                        - c32*Sel(0:N+1)*k2(0:N+1)/lambda
               B(N+1) = F(N+1) - gamma3*Ft
  end select

end subroutine EVAL_B


subroutine EVAL_F (K, H, T, C, F, N)
! Evaluates - dT * F(C) / lambda
  use STUFF;   use ROSTUFF;   implicit none
  integer        :: N
  real(kind=dbl) :: K, T, H, C(0:N+1), F(0:N+1)

  integer        :: i

  F(0) = - H**2 * C(0)
  do i = 1, N
    F(i) = - C(i-1) + 2*C(i) + K*H**2*C(i)**2 - C(i+1)
  enddo
  F(N+1) = H**2 * (-C(N+1) + 1/(1 + K*T))

end subroutine EVAL_F

subroutine SOLVE (K, H, lambda, N, ad, bd, B, ki)
! Solves the system, using Thomas
  use STUFF;   use ROSTUFF;   implicit none
  integer        :: N, np
  real(kind=dbl) :: K, H, lambda, ad(N), bd(N), B(0:N+1), &
                    ki(0:N+1)

  integer        :: i
  real(kind=dbl) :: ai

  ki(N+1) = B(N+1) / H**2
!  ai = - 1/lambda/gamma  - 2 - K*H**2*ki(N)
  ai = - 1/lambda/gamma  - 2 - 2*K*H**2*ki(N)
  ad(N) = ai;  bd(N) = B(N) - ki(N+1)
  do i = N-1, 1, -1
!    ai = - 1/lambda/gamma  - 2 - K*H**2*ki(i)
    ai = - 1/lambda/gamma  - 2 - 2*K*H**2*ki(i)
    ad(i) = ai - 1/ad(i+1);  bd(i) = B(i) - bd(i+1)/ad(i+1)
  enddo
  ki(0) = 0
  do i = 1, N
    ki(i) = (bd(i) - ki(i-1)) / ad(i)
  enddo
end subroutine SOLVE


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
