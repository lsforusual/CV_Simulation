program CV_CAT
! Block-Thomas method for CV of the simple catalytic (EC')
! mechanism A + ne- --> B;  B --> A (rate constant K,
! dimensionless). Species A and B have initial uniform concs
! of 1 and 0 resp. A plot file of the current G is output.
! The simple 2-point G-approximation is used, also for the
! boundary conditions, justified by the small X(1).
! Data input:
! 1. p-start, p-reverse, N(per p-unit), K
! 2. H1, N      (first, X(1), and N, no.of X points)
! 3. Plotting file name
  use STUFF;  implicit none
  integer        ::  N, nT, iscan, iT, Ntot, i
  real(kind=dbl) :: K, dT, H1, gamma, pstart, prev, a2,  p, dp
  real(kind=dbl),allocatable ::                              &
                    CA(:), CB(:), X(:),                      &
                    a3(:), ak(:), Adinv(:,:,:), Bd(:,:)

  call DATIN (pstart, prev, nT, N, K, H1, gamma)
  ALLOCATE (CA(0:N+1), CB(0:N+1), X(0:N+1),      &
            a3(N), ak(N), Adinv(2,2,N), Bd(2,N))
  X(0) = 0; X(1) = H1
  do i = 2, N+1;   X(i) = H1 * (gamma**i-1)/(gamma-1);   enddo
  CA = 1;  CB = 0        ! Initial concs at p-start.
  dT = 1.0_dbl / nT
  write (4,'("# CV_CAT; H1, N =", f8.3, i4)') X(1), N
  write (4,'("# pstart, prev, nT/P-unit:", 2f6.1, i3)') &
                pstart, prev, nT
  write (4,'("# K =", f10.4)') K
  write (4,'("#    p      G")')
  call PRECALC (X, N, a2, a3, ak, Adinv, gamma, K, dT)
  Ntot = NINT (nT * ABS(pstart-prev)) ! No. points per sweep
  p = pstart;  dp = -dT
  do iscan = 1, 2
    do iT = 1, Ntot
      p = p + dp
      call BACKSC(CA, CB, N, a2, a3, ak, Adinv, Bd)   ! B' scan
      call CZERO (CA, CB, N, p, Adinv, Bd)    ! Boundary values
      call FORWSC (CA, CB, N, Adinv, Bd)           ! Final scan
      write (4,'(2f12.5)') p, (CA(1)-CA(0)) / H1
    enddo
    dp = - dp ! Reverse sweep direction
  enddo
end program CV_CAT

subroutine DATIN (pstart, prev, nT, N, K, H1, gamma)
  use STUFF;   implicit none
  integer           :: nT, N
  real(kind=dbl)    :: pstart, prev, K, H1, gamma

  real(kind=dbl)    :: dT, H0, Xlim, EE_FAC

  print '(" p(start), p(rev), no/unit p, K? ")'
  read *, pstart, prev, nT, K
  print '(" H1, N? ")';   read *, H1, N
  call FILSPC (4)

  dT = 1.0_dbl / nT
  Xlim = 6 * SQRT (2*ABS(pstart-prev))  ! Outer X limit and N
  gamma = EE_FAC (H1, N, Xlim)      !      ... set gamma.

  print '(" Program CV_CAT"/)'
  print '(" p_start =", f10.4)', pstart
  print '(" p_rev   =", f10.4)', prev
  print '(" nT      =", i5, " per unit p")', nT
  print '(" N       =", i5)', N
  print '(" H1      =", f10.4)', H1
  print '(" gamma   =", f10.4)', gamma
  print '(" K       =", f10.4)', K

end subroutine DATIN

subroutine PRECALC (X, N, a2, a3, ak, Adinv, gamma, K, dT)
! Precomputes the coeffs and the A' inverse matrices.
  use STUFF;  implicit none
  integer        :: N
  real(kind=dbl) :: X(0:N+1), a2, a3(N), ak(N), &
                    Adinv(2,2,N), K, dT, gamma

  integer        :: i
  real(kind=dbl) :: alpha1, alpha2, alpha3
  real(kind=dbl),allocatable :: a1(:)

  ALLOCATE (a1(N)) ! a1 is not needed other than here.
  a2 = 1 / gamma                          ! The a's:
  do i = 1, N
    alpha1 =   2 / (x(i)-x(i-1)) / (x(i+1)-x(i-1))
    alpha2 = - 2 / (x(i)-x(i-1)) / (x(i+1)-x(i)  )
    alpha3 =   2 / (x(i+1)-x(i)) / (x(i+1)-x(i-1))
    a1(i) =   (alpha2 - 2/dT) / alpha1
    a3(i) = - (alpha2 + 2/dT) / alpha1
    ak(i) = K / alpha1
  enddo
! The pre-inverted A' matrices, Adinv. First the new A' must be
! produced. It is put directly into A'inv, and then inverted.
  Adinv(1,1,N) = a1(N);  Adinv(1,2,N) =  ak(N) ! Outer A matrix
  Adinv(2,1,N) = 0;      Adinv(2,2,N) = a1(N) - ak(N)
  call MATINV (Adinv(:,:,N), 2)
  do i = N-1, 1, -1                            ! The rest
    Adinv(1,1,i) = a1(i);  Adinv(1,2,i) =  ak(i) ! A(i)
    Adinv(2,1,i) = 0;      Adinv(2,2,i) = a1(i) - ak(i)
    Adinv(:,:,i) = Adinv(:,:,i) - a2 * Adinv(:,:,i+1) ! A'(i)
    call MATINV (Adinv(:,:,i), 2)  ! That's it.
  enddo
  DEALLOCATE (a1)
end subroutine PRECALC


subroutine BACKSC (CA, CB, N, a2, a3, ak, Adinv, Bd)
! Performs the backwards, RL scan, computing the array
! of B' vectors.
  use STUFF;  implicit   none
  integer        :: N
  real(kind=dbl) :: CA(0:N+1), CB(0:N+1), a2, a3(N), ak(N), &
              &     Adinv(2,2,N), Bd(2,N)

  integer        :: i
  real(kind=dbl) :: C(2), B(2)

! The array of A' inverses are already precalculated.
! The recursive expression for B' is:
! B'(i) = B(i) - a2*A'inv*B'(i+1)

! At i = N, things are special.
  C(1) = CA(N+1); C(2) = CB(N+1)     ! vec C is the bulk values
  B(1) = -CA(N-1) + a3(N)*CA(N) - ak(N)*CB(N) - a2*CA(N+1)
  B(2) = -CB(N-1)     + (a3(N) + ak(N))*CB(N) - a2*CB(N+1)
  Bd(:,N) = B(:) - a2*C(:)
! Now the rest:
  do i = N-1, 1, -1
    B(1) = -CA(i-1) + a3(i)*CA(i) - ak(i)*CB(i) - a2*CA(i+1)
    B(2) = -CB(i-1)     + (a3(i) + ak(i))*CB(i) - a2*CB(i+1)
    Bd(:,i) = B(:) - a2 * MATMUL(Adinv(:,:,i+1),Bd(:,i+1))
  enddo
end subroutine BACKSC

subroutine CZERO (CA, CB, N, p, Adinv, Bd)
! This calculates the new boundary, C'0 values. The simple
! 2-point flux approximation is used.
  use STUFF;  implicit none
  integer        :: N
  real(kind=dbl) :: CA(0:N+1),CB(0:N+1), p,Adinv(2,2,N),Bd(2,N)

  real(kind=dbl) :: Pee(2,2), Q(2), M0(2,2), M1(2,2), C0(2)

! Reversible system, thus the 2 boundary conditions,
! Nernst equation and flux equality, in M0 and M1,
! M0*C0 + M1*C1 = 0:
  M0(1,1) = 1;  M0(1,2) = -EXP(p);  M1(1,1) =  0;  M1(1,2) =  0
  M0(2,1) = 1;  M0(2,2) =    1   ;  M1(2,1) = -1;  M1(2,2) = -1
! Now P and Q (for the equ P*C0 = Q) (but p is taken already)
  Pee(:,:) = M0(:,:) - MATMUL(M1(:,:),Adinv(:,:,1))
  Q(:) = - MATMUL (M1(:,:), MATMUL (Adinv(:,:,1), Bd(:,1)))
  call MATINV (Pee, 2)                  ! P is now inverted
  C0(:) = MATMUL(Pee(:,:), Q(:))
  CA(0) = C0(1);  CB(0) = C0(2)         ! Final boundary values
end subroutine CZERO


subroutine FORWSC (CA, CB, N, Adinv, Bd)
! Does the final forward scan, calculating the new concs C'.
  use STUFF;  implicit none
  integer        :: N
  real(kind=dbl) :: CA(0:N+1), CB(0:N+1), Adinv(2,2,N), Bd(2,N)

  integer        :: i
  real(kind=dbl) :: C(2)

! At each i (1..N), we solve the equation
!   C'(i) = A'inv(i) * [B'(i) - C'(i-1)].
  C(1) = CA(0); C(2) = CB(0)
  do i = 1, N
     C(1:2) = MATMUL(Adinv(1:2,1:2,i), Bd(1:2,i)-C(1:2))
     CA(i) = C(1);  CB(i) = C(2)
  enddo
end subroutine FORWSC
