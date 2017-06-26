subroutine BAND_ERROR (T, I, I_exact, error)
! From a given (T,I) pair, computes the error in I, by comparison
! with the current from the function BAND_REF_FLUX.
  use STUFF;   implicit none
  real(kind=dbl) :: T, I, I_exact, error

  real(kind=dbl) :: BAND_REF_FLUX

  I_exact = BAND_REF_FLUX (T)
  error = I - I_exact
end subroutine BAND_ERROR
