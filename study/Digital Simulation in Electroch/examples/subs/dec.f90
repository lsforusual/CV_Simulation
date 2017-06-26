SUBROUTINE DEC (nmax,n,A,IP,IER)
! Minimally onverted to f90 format and adapted for use with STUFF, db, 2001.
  use STUFF;   implicit none
  INTEGER        :: nmax,n,IP(*), IER, NM1, K, KP1, M, I, J
  real(kind=dbl) :: A(nmax,nmax), T
!----------------------------------------------------------------------
!   MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
!   INPUT..
!      N = ORDER OF MATRIX.
!      NDIM = DECLARED DIMENSION OF ARRAY A . Alas, no more.
!      A = MATRIX TO BE TRIANGULARIZED.
!   OUTPUT..
!      A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U .
!      A(I,J), I.GT.J = MULTIPLIERS - LOWER TRIANGULAR FACTOR, I - L.
!      IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
!      IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR 0.
!      IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
!            SINGULAR AT STAGE K.
!   USE SOL TO OBTAIN SOLUTION OF LINEAR SYSTEM.
!   DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
!   IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
!
!   REFERENCE..
!
!      C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
!      C.A.C.M. 15(1972), P. 274.
!----------------------------------------------------------------------
   IER=0;   IP(N)=1
   IF (N > 1) then
     NM1 = N-1
     DO K = 1, NM1
       KP1 = K+1;   M = K
       DO I = KP1,N
         IF(ABS(A(I,K)) > ABS(A(M,K))) M=I
       enddo
       IP(K)=M;   T=A(M,K)
       IF(M /= K) then
         IP(N)=-IP(N)
         A(M,K)=A(K,K);   A(K,K)=T
       endif
       IF (T == 0.0) then
         ier = k;  ip(n) = 0;  exit
       endif
       T = 1.0_dbl / T
       DO I = KP1, N
         A(I,K)=-A(I,K)*T
       enddo
       DO J=KP1,N
         T=A(M,J)
         A(M,J)=A(K,J);   A(K,J)=T
         IF (T == 0.0) cycle
         DO I = KP1, N
           A(I,J) = A(I,J) + A(I,K)*T
         enddo
       enddo
     enddo
   endif
   if (n==1 .and. a(1,1)==0.0) then
     ier = 1;  ip(1) = 0
   endif
!*---------------------END OF SUBROUTINE DEC----------------------------
END SUBROUTINE DEC
