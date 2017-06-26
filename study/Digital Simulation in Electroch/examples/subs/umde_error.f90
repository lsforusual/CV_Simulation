subroutine UMDE_ERROR (T, I, I_exact, error)
! From a given (T,I) pair, computes the error in I, by comparison
! with the current from the function UMDE_FLUX.
  use STUFF;   implicit none
  real(kind=dbl) :: T, I, I_exact, error

  real(kind=dbl) :: UMDE_REF_FLUX

  I_exact = UMDE_REF_FLUX (T)
  error = I - I_exact
end subroutine UMDE_ERROR
