subroutine MATINV (A, n)
! Uses LUD (as in Press et al) to invert matrix A in place;
! that is, A is changed on return. If this is not wanted, use
! the parallel routine MAT_INV, which in any case is the more
! efficient.
! The matrix must have dimensions (n,n); use a section, e.g.
! A(1:n,1:n) in the call, if needed.
  use STUFF;    implicit none
  integer        :: n
  real(kind=dbl) :: A(n,n)

  integer                       :: i, j, jerr
  integer,allocatable           :: indx(:)
  real(kind=dbl),allocatable    :: y(:,:)
  real(kind=dbl)                :: det, Inv(2,2)

  if (n == 2) then ! Do it "by hand":
    det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
    Inv(1,1) =  A(2,2) / det;  Inv(1,2) = -A(1,2) / det
    Inv(2,1) = -A(2,1) / det;  Inv(2,2) =  A(1,1) / det
    A(1:2,1:2) = Inv(1:2,1:2)
  else
    ALLOCATE (indx(n), y(n,n))
    y = 0       ! Set up unit matrix y:
    do i=1,n
      Y(i,i)=1
    enddo
    call DEC (n, n, A, indx, jerr) ! Now LUD A
    if (jerr == 0) then
      do j=1,n
        call SOL (n, n, A, y(1,j), indx) ! jth column of inverse
      enddo
      A(1:n,1:n) = y(1:n,1:n) !Overwrite orig. A with its inverse
    else
      STOP ' MATINV: Something wrong.'
    endif
    DEALLOCATE (indx, y)
  endif
end subroutine MATINV
