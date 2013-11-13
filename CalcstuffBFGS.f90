subroutine CalcstuffBFGS(X,ndvart,fobj,dfdD,fct)
  implicit none

  integer  :: ndvart,fct
  double precision :: X(ndvart),fobj,dfdD(ndvart),x3
  double precision ::  rho, L, sigmay, pi, p, E, Fs  

  call get_f(ndvart,fct,x,fobj)
  call get_df(ndvart,fct,x,dfDD)

  return
end subroutine CalcstuffBFGS
