SUBROUTINE SOL (Nmax,N,A,B,IP)
! Minimally converted to f90 format and adapted for use with STUFF, db, 2001.
      use STUFF;  implicit none
      INTEGER        :: nmax, N, IP(*), NM1, K, KP1, M, I, KB, KM1
      real(kind=dbl) :: A(nmax,nmax), B(*), T
!-----------------------------------------------------------------------
!   SOLUTION OF LINEAR SYSTEM, A*X = B .
!   INPUT..
!     N = ORDER OF MATRIX.
!     NDIM = DECLARED DIMENSION OF ARRAY A . No no.
!     A = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
!     B = RIGHT HAND SIDE VECTOR.
!     IP = PIVOT VECTOR OBTAINED FROM DEC.
!   DO NOT USE IF DEC HAS SET IER .NE. 0.
!   OUTPUT..
!     B = SOLUTION VECTOR, X .
!-----------------------------------------------------------------------
  IF (N > 1) then
    NM1=N-1
    DO K=1,NM1
      KP1=K+1
      M=IP(K)
      T=B(M)
      B(M)=B(K)
      B(K)=T
      DO I=KP1,N
        B(I)=B(I)+A(I,K)*T
      enddo
    enddo
    DO KB=1,NM1
      KM1=N-KB
      K=KM1+1
      B(K)=B(K)/A(K,K)
      T=-B(K)
      DO I=1,KM1
        B(I)=B(I)+A(I,K)*T
      enddo
    ENDDO
  endif
  B(1)=B(1)/A(1,1)
!-------------------------END OF SUBROUTINE SOL -----------------------
END SUBROUTINE SOL
