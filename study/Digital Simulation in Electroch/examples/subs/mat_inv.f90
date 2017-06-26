subroutine MAT_INV (A, n, V)
! Uses LUD (as in Press et al) to invert matrix A to V, except
! for matrices up to 4 x 4, which are done "manually".
! The matrix must have dimensions (n,n); use a section, e.g.
! A(1:n,1:n) in the call, if needed.
! A is not changed, the result goes into V.
  use STUFF;    implicit none
  integer        :: n
  real(kind=dbl) :: A(n,n), V(n,n)

  integer                       :: i, j, jerr
  integer,allocatable           :: indx(:)
  real(kind=dbl),allocatable    :: y(:,:)
  real(kind=dbl)                :: det, Inv(2,2)

  if (n == 2) then ! Do it "by hand":
    det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
    Inv(1,1) =  A(2,2) / det;  Inv(1,2) = -A(1,2) / det
    Inv(2,1) = -A(2,1) / det;  Inv(2,2) =  A(1,1) / det
    V(1:2,1:2) = Inv(1:2,1:2)
  else if (n == 3) then
    call MAT_INV_3X3 (A, V)
  else if (n == 4) then
    call MAT_INV_4X4 (A, V)
  else
    ALLOCATE (indx(n), y(n,n))
    y = 0       ! Set up unit matrix y:
    do i=1,n
      Y(i,i)=1
    enddo
    V(1:N,1:N) = A(1:N,1:N) ! Copy for the LUD
    call DEC (n, n, V, indx, jerr) ! Now LUD A
    if (jerr == 0) then
      do j=1,n
        call SOL (n, n, V, y(1,j), indx) ! jth column of inverse
      enddo
      V(1:n,1:n) = y(1:n,1:n)
    else
      STOP ' MATINV: Something wrong.'
    endif
    DEALLOCATE (indx, y)
  endif
end subroutine MAT_INV

subroutine MAT_INV_3X3 (A, V)
  use STUFF;   implicit none
  integer   :: n
  real(dbl) :: A(3,3), V(3,3)

  real(dbl) :: det, DET3x3

  det =  DET3x3 (A)

  V(1,1) =  (A(2,2)*A(3,3) - A(3,2)*A(2,3)) / det
  V(1,2) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3)) / det
  V(1,3) =  (A(1,2)*A(2,3) - A(2,2)*A(1,3)) / det
  V(2,1) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3)) / det
  V(2,2) =  (A(1,1)*A(3,3) - A(3,1)*A(1,3)) / det
  V(2,3) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3)) / det
  V(3,1) =  (A(2,1)*A(3,2) - A(3,1)*A(2,2)) / det
  V(3,2) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2)) / det
  V(3,3) =  (A(1,1)*A(2,2) - A(2,1)*A(1,2)) / det
end subroutine MAT_INV_3X3

subroutine MAT_INV_4X4 (A, V)
  use STUFF;   implicit none
  integer   :: n
  real(dbl) :: A(4,4), T(4,4), V(4,4)

  integer   :: isign, jsign, i, j
  real(dbl) :: det, subT(3,3), DET3X3, DET4X4

  det =  DET4X4 (A)
  T = TRANSPOSE (A)
  isign = 1
  do i = 1, 4
    jsign =  isign
    do j = 1, 4
      call EXTRACT3X3 (T, i, j, subT)
      V(i,j) = jsign * DET3X3(subT)
      jsign = -jsign
    enddo
    isign = -isign
  enddo
  V = V / det
end subroutine MAT_INV_4X4


function DET3X3 (A)
! Computes the determinant of a 3 x 3 matrix.
  use STUFF;   implicit none
  real(dbl) :: A(3,3), DET3X3

  DET3X3 =   A(1,1) * (A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
           - A(1,2) * (A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
           + A(1,3) * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
end function DET3X3

function DET4X4 (A)
! Computes the determinant of a 4 x 4 matrix.
  use STUFF;   implicit none
  real(dbl) :: A(4,4), DET4X4

  integer   :: j, sign
  real(dbl) :: subA(3,3), sum, DET3X3

  sum = 0
  sign = 1
  do j = 1, 4
    call EXTRACT3X3 (A, 1, j, subA)
    sum = sum  + sign * A(1,j) * DET3X3(subA)
    sign = - sign
  enddo
  DET4X4 = sum
end function DET4X4

subroutine EXTRACT3X3 (A, row, col, subA)
! Extracts the submatrix that lacks the the ith row and the jth column.
  use STUFF;   implicit none
  integer   :: row, col
  real(dbl) :: A(4,4), subA(3,3)

  integer   :: i, j, newi, newj

  newi = 0
  do i = 1, 4          ! Select from row k
    if (i == row) cycle  ! But not if it is row to be skipped
    newi = newi + 1
    newj = 0
    do j = 1, 4          ! Now select jth element
      if (j == col) cycle  ! but not if it is col to be skipped
      newj = newj + 1
      subA(newi,newj) = A(i,j)
    enddo
  enddo
end subroutine EXTRACT3X3
