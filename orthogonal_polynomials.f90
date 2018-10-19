!=====================================================================!
!                                                                     !
!=====================================================================!

module orthogonal_polynomials

  use iso_fortran_env, only: dp => REAL64

  implicit none
  
  type :: opoly

     ! order of the polynomial
     integer :: order

   contains

     procedure :: evalf   ! evaluate the polynomial
     procedure :: evaldf  ! evaluate the derivative of polynomial
     procedure :: evald2f ! evaluate the second derivative of polynomial
     
  end type opoly

contains
  
  pure subroutine evalf(this, x, y)
    
    class(opoly), intent(in) :: this
    real(dp), intent(in)     :: x
    real(dp), intent(out)    :: y

    if (this % order .eq. 0) y = 1.0_dp
    if (this % order .eq. 1) y = x
    if (this % order .eq. 2) y = x*x - 1.0_dp
    if (this % order .eq. 3) y = x*x*x - 3.0_dp*x

  end subroutine evalf
  
  pure subroutine evaldf(this, x, dy)

    class(opoly), intent(in) :: this
    real(dp), intent(in)     :: x
    real(dp), intent(out)    :: dy

    if (this % order .eq. 0) dy = 0.0_dp
    if (this % order .eq. 1) dy = 1.0_dp
    if (this % order .eq. 2) dy = 2.0_dp*x
    if (this % order .eq. 3) dy = 3.0_dp*(x*x - 1.0_dp)

  end subroutine evaldf
  
  pure subroutine evald2f(this, x, d2y)

    class(opoly), intent(in) :: this
    real(dp), intent(in)     :: x
    real(dp), intent(out)    :: d2y

    if (this % order .eq. 0) d2y = 0.0_dp
    if (this % order .eq. 1) d2y = 0.0_dp
    if (this % order .eq. 2) d2y = 2.0_dp
    if (this % order .eq. 3) d2y = 6.0_dp*x

  end subroutine evald2f
  
end module orthogonal_polynomials

program main

  use orthogonal_polynomials, only: opoly

  type(opoly) :: hermite
  real(8) :: x, y
  
  x = 2.0d0
  hermite % order = 0
  call hermite % evalf(x, y)

  print *, y
end program main
