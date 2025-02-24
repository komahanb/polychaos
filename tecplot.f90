subroutine tecplot(ndim,dimpc,ipar,par,fct,npts,xcof) 
  use dimpce
  implicit none

  include "collsub.h"
  !#include 'mpif.h'

  integer ::ndim,npts
  integer :: ifac
  integer :: i,j,k,fct

  double precision :: yhat
  double precision :: fex
  double precision, dimension(ndim) :: x

  real*8::xcof(MAXTRM),xloc(nDIM)

  integer :: dimpc,ipar(MAXVAR)
  real*8:: par(nDIM,MAXPAR)
  integer, dimension(ndim) :: itp
  real*8::yhatprime(ndim),yhatprimetmp
  real*8::yhatdbleprime(ndim,ndim),yhatdbleprimetmp

  call get_ifac(NDIM,ifac)

  if (fct.ne.20) then
     !=========================================
     ! Exact function with output to file
     !=========================================

     open(unit=43,file='output/outputex')
     if (ndim.eq.2) then
        write(43,*) 'variables= "x","y","f"'
        write(43,*) 'zone i=',ifac,' j=',ifac,' datapacking=point'
     else if (nDim.eq.1) then
        write(43,*) 'variables= "x","f"'
        write(43,*) 'zone i=',ifac,' datapacking=point'
     end if


     !=====================================================
     ! Evaluate PC in the domain
     !=====================================================

     open(unit=45,file='output/outputPC')
     if (ndim.eq.2) then
        write(45,*) 'variables= "x","y","fPC"'
        write(45,*) 'zone i=',ifac,' j=',ifac,' datapacking=point'
     else if (nDim.eq.1) then
        write(45,*) 'variables= "x","fPC"'
        write(45,*) 'zone i=',ifac,' datapacking=point'
     end if

     
     !====================================================
     ! Polynomial chaos sample points to file
     !=====================================================
     if (fct.ne.20) then
        open(unit=44,file='output/outputsamplePC')

     else

        if (fctindx.eq.0) then
           open(unit=44,file='output/tecsamp00.dat')
        else 
           open(unit=44,file='output/tecsamp04.dat')
        end if


     end if

     if (ndim.eq.2) then
        write(44,*) 'variables= "x","y","f"'
        write(44,*) 'zone i=',npts,' datapacking=point'

     else if (nDim.eq.1) then
        write(44,*) 'variables= "x","f"'
        write(44,*) 'zone i=',npts,' datapacking=point'
     end if

     open(unit=55,file='output/outputsamplePCtmp')
     read(55,*)
     read(55,*)

     do j=1,npts
        read(55,*)(X(i),i=1,nDIM) !Read from temp file.
        call get_f(ndim,fct,x,fex)
        if (fct.eq.20) then
           x(1)=x(1)*180.0/4.0/atan(1.0)   !RADINS TO degree
        end if

        write (44,*) (X(i),i=1,nDIM),fex
     enddo

     close(44) ! final version samples

     close(55) ! temp version samples


     do i=1,ifac**ndim

        ! Initialize for safety

        fex=0.0d0 !exact value
        x=0.0d0 !location
        yhat=0.0d0 !PC value

        call make_cartesian_higherD(i,ndim,ifac,itp,x)!get the x(DIM)

        x(1:nDIM)=par(1:nDIM,1)+(par(1:nDIM,2)-par(1:nDIM,1))*x(1:nDIM) !scaling to domain size

        call evalPC(ndim,dimPC,ipar,xcof,x,yhat,yhatprime,yhatdbleprime) ! Evaluates PC at x and returns yhat 

        write (45,*) (x(k),k=1,nDIM),yhat !PC to file

        if(fct.ne.20) then

           call get_f(ndim,fct,x,fex)

           !!      x(1) = x(1) * 180.0 / 4.0 /atan(1.0)   !RADINS TO degree  
           write(43,*) (x(k),k=1,nDIM),fex !Exact to file

        end if

     end do


     close(43)

     close(45)


  end if !fct.ne.20

  return
end subroutine tecplot

!=================================================


subroutine get_ifac(NDIM,ifac)

  implicit none

  integer :: NDIM

  integer :: ifac

  if (NDIM.eq.1) then

     ifac=1001

  else if (ndim.eq.2) then

     ifac=101

  else

     ifac = int( exp(log(2.d6)/dble(ndim)))
     if(ifac.le.2)ifac = 2

  end if

end subroutine get_ifac
