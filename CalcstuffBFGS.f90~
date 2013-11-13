subroutine CalcstuffBFGS(X,ndvart,fobj,dfdD,fct)
  use dimpce,only:fctindx
  implicit none

  integer  :: ndvart,fct
  double precision :: X(ndvart),fobj,dfdD(ndvart),x3
  double precision ::  rho, L, sigmay, pi, p, E, Fs  

  fctindx=fct

  call get_f(ndvart,12,x,fobj)
  call get_df(ndvart,12,x,dfDD)


  return
end subroutine CalcstuffBFGS
