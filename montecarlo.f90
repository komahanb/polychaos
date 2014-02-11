subroutine montecarlo(stat,fct,NDIM,dimpc,nterms,npts,ipar,xcof)
  use dimpce
  implicit none
  include "collsub.h"
  include 'mpif.h'

  integer,intent(in)::stat,fct,ndim,dimpc,nterms,npts
  integer,intent(in)::ipar(MAXVAR)
  real*8,intent(in) ::xcof
  integer::ndimtmp

  real*8:: fv,gv(nDIM),hv(nDIM,nDIM)

  integer :: idec,is,ie,id,istat,ierr

  integer :: mreg(MAXDAT,NDIM)
  integer :: npdf

  real*8  :: DS(2,NDIM),scal
  real*8  :: fpcb(MAXPTS),coll(MAXPTS,MAXVAR),PL(nDIM,0:MAXTRM),DPL(nDIM,0:MAXTRM),xcoftmp
  integer :: NMCS,nmcstmp,seed,i,j,k,readMcsamples,kk,jj

  real*8 :: xavgtmp(ndim),xvartmp(ndim),xstdtmp(ndim)

  character*60 :: histname
  character*2  :: fctindxnumber

  character*60 :: filename

  real*8::average(ndim)

  double precision :: ymin,ymax,yminglb,ymaxglb,yhmin,yhmax

  integer :: ict,ictglb,ihst
  double precision :: MCm,MCmglb,MCd,MCdglb,hglb,width,pdf
  double precision, dimension(ndim):: xmin,xmax,xminglb,xmaxglb,MCmprime,MCmprimeglb,MCdprime,MCdprimeglb

  real*8,dimension(ndim,ndim)::MCmdbleprimeglb,MCmdbleprime,MCddbleprimeglb,MCddbleprime

  double precision :: st,pst,dx,p1,p2,xmintmp,xmaxtmp,dinvnorm
  double precision :: yhat,RMSE,EI,ran
  double precision, dimension(ndim)      :: x,df,x0,Ddiff,sigma,xbar,v
  double precision, dimension(ndim,ndim) :: d2f
  double precision :: f,muy1,muy2,sigmay1,sigmay2,fobjlin,fobjquad,Javg(3),Jvar(3),freal,Jstd(3)

  ! FOR MPI

  double precision, allocatable, dimension(:)   :: MNCf
  double precision, allocatable, dimension(:,:) :: MNCx
  double precision, allocatable, dimension(:,:) :: HST,HSTglb

  integer::fctindxtmp
  real*8::yhatprime(ndim),yhatprimetmp
  real*8::yhatdbleprime(ndim,ndim),yhatdbleprimetmp

  integer :: expensive
  
 
  if(id_proc.eq.0) then

!!$     ! Scaling to kriging domain
!!$
!!$     do k=1,ndim
!!$        xi(k) = (xi(k)-DS(1,k))/(DS(2,k)-DS(1,k))
!!$     end do
!!$

     write(filenum,'(1x,a)')'===================================='
     write(filenum,'(1x,a)')'  Cross Validation of the model     '
     write(filenum,'(1x,a)')'===================================='

     write(filenum,*)
     write(filenum,*)

     write(filenum,'(1x,a,2e15.5)')'>> Evaluating PC at  :', xi
     call evalPC(ndim,dimPC,ipar,xcof,xi,yhat,yhatprime,yhatdbleprime) ! Evaluates PC at x and returns yhat, yhatprime, yhatdbleprime
     write(filenum,'(1x,a,e15.5)')'>> PCE  Value f(xi)   :', yhat

     write(filenum,*)
     write(filenum,*)

     write(filenum,'(1x,a,2e15.5)')'>> Evaluating Real Fn at  :', xi

     call evalcostf(stat,ndim,fct,xi,freal,df,d2f) 

     write(filenum,'(1x,a,e15.5)')'>> Real     Value f(xi)   :', freal

     write(filenum,*)
     write(filenum,*)

     cverror=abs(yhat-freal)

     write(filenum,'(1x,a,e15.5)')'>> CV Error       :', cverror

     write(filenum,*)


  end if

  call MPI_BCAST(cverror,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  return
end subroutine montecarlo
!!$!================================================
!!$subroutine find_st(pst,dx,st)
!!$  implicit none
!!$  ! find st which ensures the probability within [-dx:dx] is pst
!!$  double precision, intent(in)  :: pst,dx
!!$  double precision, intent(out) :: st
!!$  integer :: iter
!!$  double precision :: pout,s,ds,vv,vp,dv,sini
!!$
!!$  if(pst.le.0.d0.or.pst.ge.1.d0)stop'pst in find_st'
!!$  pout = (1.d0-pst)/2.d0
!!$  sini = 1.d0
!!$190 continue
!!$  s    = sini
!!$  ds   = sini*1.d-3
!!$  iter = 0
!!$200 continue
!!$  iter = iter + 1
!!$  call CDF(-1.d0*dx,0.d0,s,   vv)
!!$  if(dabs(vv-pout).le.1.d-10)go to 210
!!$  call CDF(-1.d0*dx,0.d0,s+ds,vp)
!!$  dv = (vp-vv)/ds
!!$  !       write(filenum,'(5e15.5)')s,vv,pout,vv-pout,dv
!!$  if(dv.eq.0.d0)stop'dv = 0.d0 in find_st'
!!$  if(iter.ge.100)stop'iter>100 in find_st'
!!$  s = s - (vv-pout)/dv
!!$  if(s.le.0.d0)then
!!$     sini = sini * 0.1d0
!!$     go to 190
!!$  end if
!!$  go to 200
!!$210 continue
!!$  !       write(filenum,'(4e15.5)')s,vv,pout,vv-pout
!!$  st = s
!!$
!!$end subroutine find_st
!!$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$subroutine find_x(mode,ran,xc,st,xout)
!!$  implicit none
!!$  ! find xout which ensures the CDF at xout is ran
!!$  ! mode=0 : analytical CDF (fast but less robust?)
!!$  ! mode=1 : numerical  CDF (time comsuming, but robust)
!!$  integer, intent(in) :: mode
!!$  double precision, intent(in)  :: ran,xc,st
!!$  double precision, intent(out) :: xout
!!$  integer :: iter 
!!$  double precision :: x,vv,dv,vp,dx
!!$
!!$  x    = xc
!!$  dx   = 1.d-4
!!$  iter = 0
!!$200 continue
!!$  iter = iter + 1
!!$  if(mode.eq.0)then
!!$     call CDF(x,xc,st,vv)
!!$  else
!!$     call CDF_Numerical(x,xc,st,vv)
!!$  end if
!!$  if(dabs(vv-ran).le.1.d-10)go to 210
!!$  if(mode.eq.0)then
!!$     call CDF(x+dx,xc,st,vp)
!!$  else
!!$     call CDF_Numerical(x+dx,xc,st,vp)
!!$  end if
!!$  dv = (vp-vv)/dx
!!$  !       write(filenum,'(4e15.5)')x,vv,vv-ran,dv
!!$  if(dv.eq.0.d0)stop'dv=0 in find_x'
!!$  if(iter.ge.100)stop'iter>100 in find_x'
!!$  x = x - (vv-ran)/dv
!!$  go to 200
!!$210 continue
!!$  !       write(filenum,'(3e15.5)')x,vv,vv-ran
!!$  xout = x
!!$
!!$end subroutine find_x
!!$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$subroutine CDF(xin,xc,st,vout)
!!$  implicit none
!!$  double precision, intent(in)  :: xin,xc,st
!!$  double precision, intent(out) :: vout
!!$  double precision :: vtmp
!!$  !       vout = 0.5d0 * (1.d0 + erf( (xin-xc)/(st*dsqrt(2.d0)) ))
!!$  call ERF_MINE1( (xin-xc)/(st*dsqrt(2.d0)), vtmp )
!!$  vout = 0.5d0 * (1.d0 + vtmp)
!!$end subroutine CDF
!!$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$subroutine DCDF(xin,xc,st,dvout)
!!$  implicit none
!!$  double precision, intent(in)  :: xin,xc,st
!!$  double precision, intent(out) :: dvout
!!$  double precision :: dvtmp
!!$  call DERF_MINE( (xin-xc)/(st*dsqrt(2.d0)), dvtmp )
!!$  dvout = 0.5d0*dvtmp/(st*dsqrt(2.d0))
!!$end subroutine DCDF
!!$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$subroutine ERF_MINE1(xin,yout)
!!$  implicit none
!!$  double precision, intent(in)  :: xin
!!$  double precision, intent(out) :: yout
!!$  integer :: i,k,n
!!$  double precision :: vsum,kai
!!$  ! n is the order of Taylor
!!$  ! Maybe accurate within the range of [-4:4] with n=100
!!$  n = 100
!!$  vsum = 0.d0
!!$  do i=0,n
!!$     kai = 1.d0
!!$     do k=1,i
!!$        kai = kai * (-1.d0) * xin**2 / dble(k)
!!$     end do
!!$     vsum = vsum + kai*xin/dble(2*i+1)
!!$  end do
!!$  yout = vsum*2.d0/(dsqrt(3.141592653589793238d0))
!!$
!!$  if(yout.gt.1.d0)write(*,'(a,2e15.5)')'*ERF>1 ',xin,yout-1.d0
!!$end subroutine ERF_MINE1
!!$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$subroutine DERF_MINE(xin,dyout)
!!$  implicit none
!!$  double precision, intent(in)  :: xin
!!$  double precision, intent(out) :: dyout
!!$  double precision :: vsum
!!$
!!$  vsum  = exp(-1.d0*xin**2)
!!$  dyout = vsum*2.d0/(dsqrt(3.141592653589793238d0))
!!$
!!$end subroutine DERF_MINE
!!$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$subroutine CDF_Numerical(xin,xc,st,cdf)
!!$  implicit none
!!$  double precision, intent(in)  :: xin,xc,st
!!$  double precision, intent(out) :: cdf
!!$  double precision :: vtmp
!!$  integer :: i,num
!!$  double precision :: xs,xe,dx,x1,x2,pdf1,pdf2
!!$  if(xin.lt.xc)then
!!$     cdf = 0.d0
!!$     xs  = xin -2.d0
!!$     xe  = xin
!!$  else if(xin.ge.xc)then
!!$     cdf = 0.5d0
!!$     xs  = xc
!!$     xe  = xin
!!$  end if
!!$  num = 1001
!!$  dx  = (xe-xs)/dble(num-1)
!!$  do i=1,num-1
!!$     x1 = xs + dble(i-1)*dx
!!$     x2 = xs + dble(i  )*dx
!!$     call normal_dist(x1,xc,st,pdf1)
!!$     call normal_dist(x2,xc,st,pdf2)
!!$     cdf = cdf + (pdf1+pdf2)*dx*0.5d0
!!$  end do
!!$end subroutine CDF_Numerical
!!$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$subroutine normal_dist(xin,xc,st,y)
!!$  implicit none
!!$  double precision, intent(in)  :: xin,xc,st
!!$  double precision, intent(out) :: y
!!$  double precision :: pi
!!$
!!$  pi = 4.d0*datan(1.d0)
!!$  y = exp(-1.d0*(xin-xc)**2/2.d0/st/st)/dsqrt(2.d0*pi*st**2)
!!$
!!$end subroutine normal_dist
!!$
!!$!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
