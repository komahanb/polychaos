subroutine evalcostf(stat,dim,fct,x,fv,gv,hv) !hpcb to be implemented
  use dimpce
  use omp_lib

  implicit none

  integer ::stat,fct,dim
  real*8 ::x(DIM),fv,gv(DIM),hv(DIM,DIM)

  !CFD
  integer:: flag
  real*8:: v(dim),time,temp

  !=======================================================
  ! Evaluate Function value and gradients and also hessian
  !======================================================

  if (fct.le.8) then

     call get_f(dim,fct,x(1:DIM),fv)

     if (stat.gt.0) then
        call get_df(dim,fct,x(1:DIM),gv)
     end if

     if (stat.gt.1) then
        call get_dff(dim,fct,x(1:DIM),hv)
     end if

  else if (Fct.eq.10) then

     ! Flow solver needs input 1) Angle of attack in rads and 2) Mach number

     !print *,'Now using coarser mesh'

     if (stat.eq.0) flag=0
     if (stat.eq.1) flag=1
     if (stat.eq.2) flag=3

     CALL chdir('lowfid') ! Comment when using fine mesh

     call omp_set_num_threads(omp_get_max_threads())

     !!     if (x(2).gt.0.75d0 .and. x(2).lt.1.15d0 ) then

     !!     call Eulersolve(x,dim,0,fv,gv,hv,0,v,fctindx)
     !        call Eulersolve(x,dim,0,fv,gv,hv,0,v,fctindx)       

     !!   print *,' >> Avoiding gradients in M [0.75,1.15] '

     !!        gv=0.0d0
     !!        hv=0.0d0

     !!     else

     call Eulersolve(x,dim,0,fv,gv,hv,flag,v,fctindx)

     !!     end if

     CALL chdir('../') !Comment when using fine mesh

     !% temp program
!!$
!!$!!     if (x(2).gt.0.75d0 .and. x(2).lt.1.15d0 ) flag=0
!!$
!!$    CALL chdir('lowfid')
!!$     call omp_set_num_threads(omp_get_max_threads())
!!$
!!$     do alpha=2,2
!!$        if (alpha.eq.1)   x(1)=0.0174532925d0 ! 1 deg
!!$        if (alpha.eq.2)   x(1)=0.0698131701d0 ! 4 deg
!!$
!!$        do fctindx=4,4,4
!!$
!!$           if (fctindx.eq.0) then !Drag
!!$              if (alpha.eq.1) then
!!$                 open(unit=99,file='dragalpha1')
!!$              else
!!$                 open(unit=99,file='dragalpha4')
!!$              end if
!!$
!!$           else if (fctindx.eq.4) then !Lift
!!$
!!$              if (alpha.eq.1) then
!!$                 open(unit=99,file='liftalpha1')
!!$              else
!!$                 open(unit=99,file='liftalpha4')
!!$              end if
!!$
!!$           end if
!!$
!!$         
!!$           write(99,*) 'variables= "x","f","g1","g2"'
!!$           write(99,*) 'zone i=',500,' datapacking=point'
!!$
!!$           do temp=0.5,1.5,0.002
!!$              x(2)=temp
!!$              print *,x,dim, v,fctindx 
!!$!              stop
!!$
!!$              call Eulersolve(x,dim,0,fv,gv,hv,flag,v,fctindx)
     !!fv=sin(x(1)+x(2))
!!$              write(99,'(4f30.16)')x(2),fv,gv(1),gv(2)
!!$           end do
!!$           close(99)
!!$           
!!$        end do
!!$     end do
!!$  
!!$     CALL chdir('../')      
!!$     stop     
!!$     

  end if
  return
end subroutine evalcostf


subroutine get_f(dim,fct,x,f)
  use dimpce
  implicit none
  
    integer :: fct,dim,k
    real*8, intent(in)  :: x(dim)
    real*8, intent(out) :: f

    real*8 :: rho, L, sigmay, pi, Fs, p, E !Truss design parameters

    f=0.0d0
 
    if (fct.eq.1) then     

       f=0.0
       do k=1,DIM
          f=f+x(k)
       end do
       f=cos(f)


    else if (fct.eq.2) then

       f=1.0
       do k=1,DIM
          f=f+x(k)**2
       end do
       f=1.0/f

    else if (fct.eq.3) then

       f=0.0
       do k=1,DIM
          f=f+x(k)**2
       end do
       f=f

    else if (fct.eq.4) then

       f=0.0
       do k=1,DIM
          f=f+x(k)
       end do
       f=exp(f)

    else if (fct.eq.5) then

       f=0.0
       do k=1,DIM
          f=f+x(k)**3
       end do
       f=f

    else if (fct.eq.6) then

       f=0.0
       do k=1,DIM-1
          f=f+100.0*(x(k+1)-x(k)**2)**2+(1-x(k))**2
       end do




    else if (fct.eq.7) then


       f=0.0
       do k=1,DIM        
          f=f+sin(3.0*x(k)-1.5)+cos(3.0*x(k)-3.0)
       end do
       f=f



    else if (fct.eq.8) then

       rho=0.2836
       sigmay=36260.0
       p=25000.0
       L=5.0
       E=30e6
       pi=4.0*atan(1.0)

       Fs=1.0

       if (fctindx.eq.0) then
          !---- OBJECTIVE FUNCTION
          f = rho*x(1)*L+rho*x(2)*sqrt(L**2+x(3)**2)
       else if (fctindx.eq.1) then
          !---- INEQUALITY CONSTRAINTS
          f = p*Fs*sqrt(L**2+x(3)**2) / (x(2)*x(3)*sigmay) - 1.0
       else if (fctindx.eq.2) then  
          f = p*Fs*L / (x(1)*x(3)*sigmay) - 1.0
       else if (fctindx.eq.3) then
          f = 4.0*p*Fs*L**3 / (x(1)**2*x(3)*E*pi) - 1.0
       end if

    else

       stop'Unsupported function number'
       
    end if

  end subroutine get_f
  
  subroutine get_df(dim,fct,x,df)
    use dimpce
    implicit none

    integer :: fct,dim,k
    real*8, intent(in)  :: x(dim)
    real*8:: fac
    real*8, intent(out) :: df(DIM)

    real*8 :: rho, L, sigmay, pi, Fs, p, E

    if (fct.eq.1) then

       fac=0.0
       do k=1,DIM
          fac=fac+x(k)
       end do

       fac=-sin(fac)

       do k=1,DIM
          df(k)=fac
       end do

    else if (fct.eq.2) then 

       fac=1.0
       do k=1,DIM
          fac=fac+x(k)**2
       end do
       fac=1.0/fac**2

       do k=1,DIM
          df(k)=-2.0*x(k)*fac
       end do

    else if (fct.eq.3) then  

       do k=1,DIM
          df(k)=2.0*x(k)
       end do

    else if (fct.eq.4) then  

       !f=exp(x+y)
       fac=0.0
       do k=1,DIM
          fac=fac+x(k)
       end do
       fac=exp(fac)
       do k=1,DIM
          df(k)=fac
       end do

    else if (fct.eq.5) then  

       !f=(x**2+y**2)

       do k=1,DIM
          df(k)=3.0*x(k)**2
       end do

 else if (fct.eq.6) then

  
       df(1)=-200.0*(x(2)-x(1)**2)*2.0*x(1)-2.d0*(1-x(1))
       do k=2,DIM-1
          df(k)=200.0*(x(k)-x(k-1)**2)-200.0*(x(k+1)-x(k)**2)*2.0*x(k)-2.0*(1-x(k))
       end do
       df(DIM)=200.0*(x(DIM)-x(DIM-1)**2)

    else if (fct.eq.7) then

       do k=1,DIM
          df(k)=3.0*cos(3.0*x(k)-1.5)-3.0*sin(3.0*x(k)-3.0d0)
       end do

   
    else if (fct.eq.8)then
       
       rho=0.2836
       sigmay=36260.0
       p=25000.0
       L=5.0
       E=30e6
       pi=4.0*atan(1.0)

       Fs=1.0

       if (fctindx.eq.0) then
          !---- OBJECTIVE FUNCTION
          df(1) = rho*L
          df(2) = rho*sqrt(L**2+x(3)**2)
          df(3) = rho*x(2)*x(3) / sqrt(L**2+x(3)**2)
       else if (fctindx.eq.1) then
!---- INEQUALITY CONSTRAINTS
          df(1) = 0.0
          df(2) =-p*Fs*sqrt(L**2+x(3)**2) / (x(2)**2*x(3)*sigmay)
          df(3) =-p*Fs*L**2 /sqrt(L**2/x(3)**2+1.0)/ (x(2)*x(3)**3*sigmay)
       else if (fctindx.eq.2) then  
          df(1) =-p*Fs*L / (x(1)**2*x(3)*sigmay)
          df(2) = 0.0
          df(3) =-p*Fs*L / (x(1)*x(3)**2*sigmay)
       else if (fctindx.eq.3) then
          df(1) =-8.0*p*Fs*L**3 / (pi*E*x(1)**3*x(3))
          df(2) = 0.0
          df(3) =-4.0*p*Fs*L**3 / (pi*E*x(1)**2*x(3)**2)
       end if

    else
       
          stop'Unsupported function number'

    end if

  end subroutine get_df

!+++++++++++++++++++++++++++++++++++++++++++


  subroutine get_dff(DIM,fct,x,d2f)
    use dimpce
    implicit none

    integer :: DIM,fct,j,k
    real*8 :: x(DIM),d2f(DIM,DIM),fac,A,omeg

    real*8 :: rho, L, sigmay, pi, Fs, p, E


    if (fct.eq.1) then !Cosine
       !f=cos(x+y)

       fac=0.0
       do k=1,DIM
          fac=fac+x(k)
       end do
       
       fac=-cos(fac)

       do j=1,DIM
          do k=1,DIM
             d2f(j,k)=fac
          end do
       end do
       
        
    else if (fct.eq.2) then !Runge
       
       !f=1.0/(1.0+x**2+y**2)

       fac=1.0
       do k=1,DIM
          fac=fac+x(k)**2
       end do
       fac=1.0/fac
       
       do j=1,DIM
          do k=1,DIM
             if (j.ne.k) then
                d2f(j,k)=8.0*x(j)*x(k)*fac**3
             else
                d2f(j,k)=(-2.0+8.0*x(k)**2*fac)*fac**2
             end if
          end do
       end do
  
    else  if (fct.eq.3) then !quadratic
       

       
       d2f(:,:)=0.0d0
       do k=1,DIM
          d2f(k,k)=2.0d0
       end do
       
    else if (fct.eq.4) then !exponential
       
       !f=cos(x+y)
       
       fac=0.0
       do k=1,DIM
          fac=fac+x(k)
       end do
       
       fac=exp(fac)

       do j=1,DIM
          do k=1,DIM
             d2f(j,k)=fac
          end do
       end do
      
     
    else if (fct.eq.5) then !Cubic OK
       
       !f=x**3+y**3
       
       d2f(:,:)=0.0
       do k=1,DIM
          d2f(k,k)=6.0*x(k)
       end do 


    else if (fct.eq.6) then  !Rosenbrock

       !f=(1.0-x)**2 + 100.0*(y-x**2)**2

       d2f(:,:)=0.0

       d2f(1,1)=(-400.0*x(2)+1200.0*x(1)**2+2.0)
       do k=2,DIM-1
          d2f(k,k)=200.0-400.0*x(k+1)+1200.0*x(k)**2+2.0
       end do
       d2f(DIM,DIM)=200.0

       do k=2,DIM
          d2f(k,k-1)=-400.0*x(k-1)
          d2f(k-1,k)=d2f(k,k-1)
       end do

    else if (fct.eq.7) then ! 6:sin(3x-1.5)+cos(3y-3)
       
       d2f(:,:)=0.0
       do j=1,dim
          do k=1,dim
             d2f(j,k)= 3.0*cos(3.0*x(k) - 1.5d0) - 3.0*sin(3.0*x(k) - 3.0)
          end do
       end do

    else if (fct.eq.8) then

 
       rho=0.2836
       sigmay=36260.0
       p=25000.0
       L=5.0
       E=30e6
       pi=4.0*atan(1.0)

       Fs=1.0

       d2f(:,:)=0.0d0
       if (fctindx.eq.0) then
          !---- OBJECTIVE FUNCTION
          d2f(1,1) = 0.0d0
          d2f(2,1) = 0.0d0
          d2f(3,1) = rho*x(2)*x(3) / sqrt(L**2+x(3)**2)

       else if (fctindx.eq.1) then
!---- INEQUALITY CONSTRAINTS
!!$          df(1) = 0.0
!!$          df(2) =-p*Fs*sqrt(L**2+x(3)**2) / (x(2)**2*x(3)*sigmay)
!!$          df(3) =-p*Fs*L**2 /sqrt(L**2/x(3)**2+1.0)/ (x(2)*x(3)**3*sigmay)
       else if (fctindx.eq.2) then  
!!$          df(1) =-p*Fs*L / (x(1)**2*x(3)*sigmay)
!!$          df(2) = 0.0
!!$          df(3) =-p*Fs*L / (x(1)*x(3)**2*sigmay)
       else if (fctindx.eq.3) then
!!$          df(1) =-8.0*p*Fs*L**3 / (pi*E*x(1)**3*x(3))
!!$          df(2) = 0.0
!!$          df(3) =-4.0*p*Fs*L**3 / (pi*E*x(1)**2*x(3)**2)
       end if

    else

       print *,'Warning',fct
       stop'Invalid fct for gradients'

    end if
      
  end subroutine get_dff
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine gather_costfn(dim,fct,stat,npts,RN,fpcb,gpcb,hpcb)
    use dimpce
    implicit none
    include "collsub.h"

    integer,intent(in)::dim,fct,stat,npts
    integer::i,j,k,ii
    real*8,intent(in)::RN(DIM,MAXPTS)
    real*8::x(dim)
    real*8::fv,gv(dim),hv(dim,dim)

    real*8,intent(out)::fpcb(MAXPTS),gpcb(DIM,MAXPTS),hpcb(dim,dim,maxpts)
    !====================================================
    ! Polynomial chaos sample points to file
    !=====================================================

    if (fct.ne.10) then
       open(unit=44,file='output/outputsamplePC')

    else if (fct.eq.10) then

       if (fctindx.eq.0) then
          open(unit=44,file='output/tecsamp00.dat')
       else 
          open(unit=44,file='output/tecsamp04.dat')
       end if

    end if

    if (dim.eq.2) then
       write(44,*) 'variables= "x","y","f"'
       write(44,*) 'zone i=',npts,' datapacking=point'

    else if (Dim.eq.1) then
       write(44,*) 'variables= "x","f"'
       write(44,*) 'zone i=',npts,' datapacking=point'
    end if

    do j=1,npts

       x(1:DIM)=RN(1:DIM,j)
       !                 print *,'x:',x
       call evalcostf(stat,dim,fct,x,fv,gv,hv)
       fpcb(j)=fv                   !Store the cost function value at appropriate location
       if (stat.gt.0) then
          gpcb(1:dim,j)=gv(1:dim)   !Store them appropriately
       end if

       if (stat.gt.1) then
          hpcb(1:dim,1:dim,j)=hv(1:dim,1:dim)
       end if

       if (fct.eq.10) then
          x(1)=x(1)*180.0/4.0/atan(1.0)   !RADINS TO degree
       end if

       write (44,*) (x(ii),ii=1,dim),fv
!       write (filenum,*) (x(ii),ii=1,dim),fv

    end do ! loop over all points

    close(44)

    return
  end subroutine gather_costfn
