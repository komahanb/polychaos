!================================== Objective functions

subroutine threebarf(fctindx,dat,x,f)
  implicit none

  integer,intent(in) ::fctindx
  real*8,intent(in)  ::dat(20)
  real*8,intent(in)  ::x(6)
  real*8::E,H,gamma,theta,P,phi(3),L(3),Pu,Pv
  real*8::tensile_sigma1_max,tensile_sigma2_max,tensile_sigma3_max
  real*8::comp_sigma1_max,comp_sigma2_max,comp_sigma3_max,max_u_disp,max_v_disp
  real*8::Fs,sigma(3),u,v
  real*8,intent(out)::f

  ! Read/store needed params

  H=dat(1)
  E=dat(2)
  gamma=dat(3)
  theta=dat(4)
  P=dat(5)

  tensile_sigma1_max=dat(6)      
  tensile_sigma2_max=dat(7)
  tensile_sigma3_max=dat(8)

  comp_sigma1_max=dat(9)      
  comp_sigma2_max=dat(10)
  comp_sigma3_max=dat(11)

  max_u_disp=dat(12)
  max_v_disp=dat(13)

  Fs=dat(14)

  phi(1)=x(4)
  phi(2)=x(5)
  phi(3)=x(6)

  L(1)=H/sin(phi(1))
  L(2)=H/sin(phi(2))
  L(3)=H/sin(phi(3))


  pu=P*cos(theta)
  pv=P*sin(theta)

  
!  if(fctindx.gt.0) then

  f=0.0 ! initialize

  if (fctindx.eq.0) then
     F= x(1)*gamma*L(1) + x(2)*gamma*L(2) +  x(3)*gamma*L(3)
  else if (fctindx.eq.1) then
     F = (sigma(1) - tensile_sigma1_max)/tensile_sigma1_max   !tensile 1
  else if (fctindx.eq.2) then
     F = (sigma(2) - tensile_sigma2_max)/tensile_sigma2_max   !tensile 2
  else if (fctindx.eq.3) then 
     F = (sigma(3) - tensile_sigma3_max)/tensile_sigma3_max   !tensile 3
  else if (fctindx.eq.4) then
     F = (-1.0*sigma(1)/comp_sigma1_max) -1.0                 !compressive 1
  else if (fctindx.eq.5) then
     F = (-1.0*sigma(2)/comp_sigma2_max) -1.0                 !compressive 2
  else if (fctindx.eq.6) then
     F = (-1.0*sigma(3) / comp_sigma3_max) -1.0               !compressive 3
  else if (fctindx.eq.7) then
     F = (u -max_u_disp)/max_u_disp
  else if (fctindx.eq.8) then
     F = (v -max_v_disp)/max_v_disp
  else
     stop'Wrong Fctindx for 6d threebar'
  end if

  return
end subroutine threebarf

!===================================== Gradients


subroutine threebardf(fctindx,dat,x,df)
  implicit none

  integer,intent(in) ::fctindx
  real*8,intent(in) ::dat(20)
  real*8,intent(in)  ::x(6)
  real*8,intent(out)::df(6)
  real*8::E,H,gamma,theta,P,phi(3),L(3),Pu,Pv
  real*8::tensile_sigma1_max,tensile_sigma2_max,tensile_sigma3_max
  real*8::comp_sigma1_max,comp_sigma2_max,comp_sigma3_max,max_u_disp,max_v_disp
  real*8::Fs,sigma(3),u,v

  ! Read/store needed params

  H=dat(1)
  E=dat(2)
  gamma=dat(3)
  theta=dat(4)
  P=dat(5)

  tensile_sigma1_max=dat(6)      
  tensile_sigma2_max=dat(7)
  tensile_sigma3_max=dat(8)

  comp_sigma1_max=dat(9)      
  comp_sigma2_max=dat(10)
  comp_sigma3_max=dat(11)

  max_u_disp=dat(12)
  max_v_disp=dat(13)

  Fs=dat(14)

  phi(1)=x(4)
  phi(2)=x(5)
  phi(3)=x(6)

  L(1)=H/sin(phi(1))
  L(2)=H/sin(phi(2))
  L(3)=H/sin(phi(3))

  pu=P*cos(theta)
  pv=P*sin(theta)

  df=0.0 !initialize 


     df(1)=(H*gamma)/sin(phi(1))

     df(2)=(H*gamma)/sin(phi(2))

     df(3)=(H*gamma)/sin(phi(3))

     df(4)=-(X(1)*H*gamma*cos(phi(1)))/sin(phi(1))**2.0

     df(5)=-(X(2)*H*gamma*cos(phi(2)))/sin(phi(2))**2.0

     df(6)=-(X(3)*H*gamma*cos(phi(3)))/sin(phi(3))**2.0

  return
end subroutine threebardf

!===================================== Hessian

!!$
!!$subroutine threebard2f()
!!$  implicit none
!!$
!!$! Just kidding
!!$end subroutine threebard2f

