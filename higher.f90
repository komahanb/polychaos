subroutine RMSE_Higher(stat,ndim,fct,npts,dimPC,ipar,par,xcof)

  use dimpce
  !  use timer_mod

  implicit none

  include "collsub.h"
  include 'mpif.h'


  integer ::ndim
  integer :: ifac,ierr
  integer :: i,j,k,istat,stat,fct,npts,iii,jjj,kk,jj
  integer :: idec,is,ie,imode,id
  integer :: ict,ictglb,idy,idyy,ide,idee
  integer, dimension(ndim) :: itp
  double precision :: yhat,RMSE,EI,maxerror,maxerrorglb
  double precision :: yprimehat,RMSEhat,maxerrorhat,maxerrorhatglb, errorhat,errorhatglb
  double precision :: ydbleprimehat,RMSEdblehat,maxerrordblehat,maxerrordblehatglb, errordblehat,errordblehatglb

  double precision, dimension(ndim) :: x,xe,xy,df,v

  real*8::xcof(MAXTRM),xloc(nDIM),par(nDIM,MAXPAR)
  integer :: dimpc,ipar(MAXVAR)

  double precision:: error, yhatmin,yhatglb, fv,gv(nDIM),hv(NDIM,NDIM),errorglb,diffloc,distloc
  character*60 :: filename
  character(len=13) :: Cexac
  double precision,dimension(51,51,4) :: TEC2
  !CFD
  real*8::yhatprime(ndim),yhatprimetmp
  real*8::yhatdbleprime(ndim,ndim),yhatdbleprimetmp

  call get_ifac(NDIM,ifac)

  if (ndim.eq.5)  ifac=10

  if(ndim.eq.3) ifac=30

  if (ndim.eq.2) then
     if (fct.eq.20) then
        ifac = 51 ! I have database of only  51*51 for CFD
     else 
        ifac=101
     end if
  end if




  if (id_proc.eq.0) then

     call getfilename(ndim,fct,dimpc,stat,0,filename)

     outfile(:)=filename(:)

     if (dyncyccnt.eq.1)  then

        open(51,file=filename,form='formatted',status='unknown')!,position='append')
     else

        open(51,file=filename,form='formatted',status='old',position='append')

     end if

     if (dyncyccnt.eq.1)   write(51,'(a)') 'dimpc    npts    diff       maxdiff graddiff hessdiff  Maxdisc Meandisc'



     write(filenum,*)
     write(filenum,*)'     >> RMSE on ',ifac**ndim,' points . . .'
     write(filenum,*)

  end if !master

  call MPI_Barrier(MPI_COMM_WORLD,ierr) ! All wait until master gets the coeffs

  !===============================
  ! CFD TEST CASE SPECIAL CODING
  !===============================

  if (id_proc.eq.0.and.ndim.eq.2.and.fct.eq.20) then

     write(filenum,*) '>>  Reading from the existing CFD database for RMSE'

     if (fctindx.eq.4) then
        open(15,file='output/tecex10Lift',form='formatted',status='old')

     else if (fctindx.eq.0) then

        open(15,file='output/tecex10Drag',form='formatted',status='old')

     else 

        stop'Wrong fctindx-->Check it'

     end if

     read(15,*)
     read(15,*)
     read(15,*)
     read(15,*)((TEC2(iii,jjj,1),iii=1,ifac),jjj=1,ifac) ! is in degrees
     read(15,*)((TEC2(iii,jjj,2),iii=1,ifac),jjj=1,ifac) ! mach number
     read(15,*)((TEC2(iii,jjj,3),iii=1,ifac),jjj=1,ifac) ! lift/drag
     close(15)                            

     TEC2(:,:,1) = TEC2(:,:,1) *4.0d0*atan(1.0)/180.0  !Deg to rad

     ! Find RMSE between Exact database and PC Surrogate

     error=0.0d0 !Initialize

     do i=1,ifac
        do j=1,ifac

           fv=  TEC2(i,j,3)

           x(1)=TEC2(i,j,1)
           x(2)=TEC2(i,j,2)

           call evalPC(ndim,dimPC,ipar,xcof,x,yhat,yhatprime,yhatdbleprime) ! Evaluates PC at x and returns yhat 

           tec2(i,j,4)=yhat 

           error = error + (yhat-fv)**2
          !         print *, x,fv,yhat,error
           if (abs(fv-yhat).gt.maxerror) maxerror=abs(fv-yhat)

        end do
     end do

     error = dsqrt(error/dble(ifac*ifac)) !Normalize

     write(filenum,*)
     write(filenum,'(6x,a,e20.10)') '>> RMSE compared to Database = ',error
     write(51,'(2i8,3e15.8)') dimpc,npts,error,maxerror !write to file



     ! Tecplot output of PC surrogate


     TEC2(:,:,1) = TEC2(:,:,1) *180.0 / 4.0 /atan(1.0)   !RADINS TO degree

     if (fctindx.eq.4) then
        open(110,file='output/tecex104.dat',form='formatted',status='unknown')
     else if (fctindx.eq.0) then
        open(110,file='output/tecex100.dat',form='formatted',status='unknown')
     end if

     write(110,'(a)')'TITLE = " "'
     write(110,'(a)')'VARIABLES = "x" "y" "f"'
     write(110,'(3(a,i5),a)')'ZONE T="Kriging", I=',ifac,', J=',ifac,', K=',1,', F=BLOCK'
     write(110,*)((TEC2(i,j,1),i=1,ifac),j=1,ifac)
     write(110,*)((TEC2(i,j,2),i=1,ifac),j=1,ifac)
     write(110,*)((TEC2(i,j,3),i=1,ifac),j=1,ifac)
     close(110)

     if (fctindx.eq.4) then
        open(110,file='output/tecpc104.dat',form='formatted',status='unknown')
     else if (fctindx.eq.0) then
        open(110,file='output/tecpc100.dat',form='formatted',status='unknown')
     end if

     write(110,'(a)')'TITLE = " "'
     write(110,'(a)')'VARIABLES = "x" "y" "f"'
     write(110,'(3(a,i5),a)')'ZONE T="Kriging", I=',ifac,', J=',ifac,', K=',1,', F=BLOCK'
     write(110,*)((TEC2(i,j,1),i=1,ifac),j=1,ifac)
     write(110,*)((TEC2(i,j,2),i=1,ifac),j=1,ifac)
     write(110,*)((TEC2(i,j,4),i=1,ifac),j=1,ifac)
     close(110)

     return

  end if ! CFD



  do k=1,1!nfunc-1

     idec = dble(ifac**ndim)/dble(num_proc)
     is   = idec*id_proc + 1
     ie   = idec*(id_proc+1)
     if(id_proc.eq.num_proc-1)ie = ifac**ndim
     write(filenum,'(11x,a,2i10,a,i3)') '>> [',is,ie,' ] for Processor',id_proc

     ict     = 0
     error    = 0.d0
     errorhat=0.0d0
     errordblehat=0.0d0

     yhatmin = 1.d10
     maxerror= 0.d0

     do i=is,ie

        ! Initialize for safety
        !exact value
        fv=0.0d0 
        gv=0.0d0 
        hv=0.0d0

        x=0.0d0 !location

        !PC value
        yhat=0.0d0
        yhatprime=0.0d0
        yhatdbleprime=0.0d0

        call make_cartesian_higherD(i,ndim,ifac,itp,x) !get the x(DIM)
!        x(1)=0.0d0
!        x(2)=0.0d0
!        x(3)=0.0d0

     
        x(1:nDIM)=par(1:nDIM,1)+(par(1:nDIM,2)-par(1:nDIM,1))*x(1:nDIM) !scaling x to domain size
      
        call evalPC(ndim,dimPC,ipar,xcof,x,yhat,yhatprime,yhatdbleprime) ! Evaluates PC at x and returns yhat,yhatprime,yhatdbleprime

!!$        if(yhat.le.yhatmin)then
!!$           yhatmin = yhat
!!$           xy(:)   = x(:)
!!$        end if

        call evalcostf(2,ndim,fct,x,fv,gv,hv) !  to force computing grad and Hess norm

        !   call get_f(ndim,fct,x,fv)

       !        print *,'x:',X,'Exact:',fv,'PC:',yhat
        !      print*,TEC2(i,1),tec2(i,2)
!stop
        ! Error in Cost function
        error = error + (yhat-fv)**2
        if (abs(fv-yhat).gt.maxerror) maxerror=abs(fv-yhat)

        ! Error in gradients
        do kk=1,ndim
           errorhat=errorhat+(yhatprime(kk)-gv(kk))**2
        end do

        ! Error in Hessian
        do jj=1,ndim
           do kk=1,ndim
              errordblehat=errordblehat+(yhatdbleprime(kk,jj)-hv(kk,jj))**2
           end do
        end do
!        if (i.eq.125000) print *,'x:',X,'Exact:',fv,'PC:',yhat
     end do ! main loop for Cartesian (i)
!stop

     ! Informtion Sharing

     call MPI_Barrier(MPI_COMM_WORLD,ierr) ! All wait until master gets the coeffs

     call MPI_ALLREDUCE(error,errorglb,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     errorglb = dsqrt(errorglb/dble(ifac**ndim)) 

     call MPI_ALLREDUCE(errorhat,errorhatglb,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     errorhatglb = dsqrt(errorhatglb/dble(ifac**ndim)) 

     call MPI_ALLREDUCE(errordblehat,errordblehatglb,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     errordblehatglb = dsqrt(errordblehatglb/dble(ifac**ndim)) 

     call MPI_ALLREDUCE(maxerror,maxerrorglb,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)

!!$
!!$     call MPI_ALLREDUCE(yhatmin,yhatglb,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
!!$     idy = -1
!!$     if(yhatmin.eq.yhatglb) idy = id_proc
!!$     call MPI_ALLREDUCE(idy,idyy,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
!!$
!!$     if(idyy.lt.0.or.idyy.ge.num_proc) stop'idyy'
!!$     call MPI_BCAST(xy(1),ndim,MPI_DOUBLE_PRECISION,idyy,MPI_COMM_WORLD,ierr)
     
     if(id_proc.eq.0)then

        write(filenum,*)
        write(filenum,'(6x,a,e20.10)') '>> RMSE compared to Analytical function = ',errorglb
        write(filenum,*)

        ! write to errornorm file
        write(51,'(2i8,6e15.8)') dimpc,npts,errorglb,maxerrorglb, errorhatglb, errordblehatglb,difflocmax,diffloc2 !, ttotal !write to file

        rmsemat(runnum,dyncyccnt,1)=npts
        rmsemat(runnum,dyncyccnt,2)=errorglb
        outfile=filename

     end if

  end do ! function loop (k)

  close(51) !close the file

  return
end subroutine RMSE_Higher

!=================================================



subroutine make_cartesian_higherD(inp,ndim,ifac,itp,x)
  implicit none
  integer, intent(in) :: inp,ndim,ifac
  integer, dimension(ndim), intent(out) :: itp
  double precision, dimension(ndim), intent(out) :: x
  integer :: j,l,indx

  j  = inp

  do 100 l=1,ndim

     indx = int(j/ifac**(ndim-l))
     if(mod(j,ifac**(ndim-l)).eq.0) indx = indx - 1
     if(indx.eq.0.or.indx.eq.ifac-1)then
     else if(indx.ge.1.and.indx.le.ifac-2)then
     else
        stop'invalid indx'
     end if
     itp(l) = indx
     x(l)   = dble(indx)/dble(ifac-1)
     j      = j - indx*(ifac**(ndim-l))
     if (x(l).lt.0.d0.or.x(l).gt.1.d0) stop'x(l)>1 or x(l)<0'

100  continue

   end subroutine make_cartesian_higherD
