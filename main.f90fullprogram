program main
  use dimpce
  !  use timer_mod
  implicit none

  INCLUDE "collsub.h"
  include 'mpif.h'

  integer :: DIM
  parameter (DIM=6)

  !indices
  integer :: i,j,k,ii,jj,kk,fuct

  !Monte Carlo
  integer::nmcs,readMCsamples

  !flags
  integer :: stat,makesamples,solver,ierr

  !PC vitals
  integer::DIMPC,numpc,npts,nterms
  real*8 :: coll(MAXPTS,DIM),par(DIM,MAXPAR)
  real*8 :: RN(DIM,MAXPTS),DS(2,DIM)
  real*8 :: fpcb(MAXPTS),gpcb(DIM,MAXPTS),xcof(MAXTRM),hpcb(dim,dim,maxpts)
  integer :: ipar(MAXVAR)


  !RandomNumber 
  integer:: seed

  !LU Decomposition
  integer :: indx(MAXDAT)
  real*8 :: W(MAXDAT),V(MAXDAT,MAXDAT),rhsF(MAXTRM),z(MAXDAT),tmp(MAXTRM)


  !Orthogonal Polynomial related
  integer::mreg(MAXDAT,DIM)
  real*8 :: PL(DIM,0:MAXTRM),DPL(DIM,0:MAXTRM),xmat(MAXDAT,MAXDAT)


  !Name related declarations
  character*60 :: filename
  character*2 :: dimnumber,fctnumber

  ! Cost func evals
  integer::fct
  real*8:: fv,gv(DIM),hv(DIM,DIM)
  real*8 :: x(DIM)

  !!Dynamic samples

  real*8::xvec(maxdat,dim)
  integer::Randomini,maxorderwant
  integer::index,evalfunction
  integer::nptsold,ntermsold
  integer::  nptstoaddpercyc
  
  call MPI_START 
 
  mainprog=.true.

  !Settings
  
  filenum=6 ! 6 for screen, any other number for fort.x
  
   if(id_proc.eq.0)  then
     write(filenum,*)
     write(filenum,*) '======================================================='
     write(filenum,*) 'Non-Intrusive Polynomial Chaos (Stochastic Collocation)'
     write(filenum,*) '======================================================='
     write(filenum,*) 'Author: Komahan Boopathy, University of Dayton, OH'
     write(filenum,*) 'Email : komahan.cool@gmail.com'
     write(filenum,*) '======================================================='
     write(filenum,*)
!     write(filenum,'(x,a,i3)')'>> Number of Processors = ',num_proc
  end if !master

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  !                 Non-parallel region                        !
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

  !============================================================
  ! Initial settings
  !============================================================

  
  makesamples=1 ! 0=read, 1= Make via LHS for building surrogate

  ! Choice of orthogonal basis
  ! 1=Legendre
  ! 2=Hermite

  do j=1,DIM
     ipar(j)=1  
  end do

  casemode=1 !0=RMSE only, 1=Stats+RMSE

  evlfnc=1  ! For montecarlo, should the program evaluate the exact fn (CFD)?

  if (casemode.eq.1) then
     ! This file is read again in montecarlo subroutine. Here it is needed to define the domain size when doing casemode=1(Stats)
     open(10,file='MC.inp',form='formatted',status='old')
!     read(10,*) (xavg(i),i=1,dim)
!     read(10,*) (xstd(i),i=1,dim)     
     read(10,*)
     read(10,*)
     read(10,*) NMCS!,ndimtmp
     read(10,*) !npdf
     read(10,*) readMCsamples
     read(10,*) evlfnc
     close(10)
  end if

  ! Mean setup
  do i =1,dim
  xavg(i)=1.0d0
  end do

  probtype=1
  
  !variance setup 
  if (probtype.eq.1) then
    xstd(1:dim)=0.05
  else if (probtype.eq.2) then
    xstd(1:dim)=xavg(1:dim)*0.05
  else	
    stop"Wrong prob type"
  end if	


  do  dynamics=1,1

     !===============================
     ! Main Program
     !===============================

     DO OS=2,2 ! Ratio of Over Sampling ratio 1 or 2 (2 is recommended)

        do  stat=0,0,1   
           
           !0= Function only
           !1= Function + Gradient
           !2= Function +Gradient +Hessian   

           call solvertype(stat,os,solver)

           !1 : cos(x+y) (Nd)
           !2 : 1.0/(1.0+x**2+y**2)  (Nd)
           !3 : x**2+y**2  (Nd)
           !4 : exp(x+y)  (Nd)
           !5 : x**3+y**3 (Nd)
           !6 : Rosenbrock (Nd)
           !7 : sin(3x-1.5)+cos(3y-3) (Nd)
           !8 : Two bar truss design (3d)
           !9 : Short tubular column (3d)
           !10: Cantilever beam (2d)
           !11: Three bar truss (3d)
           !12: Threebar truss (6d)           
           !20: CFD

           do fct=1,1,1

              fctindx=0 

!!$              if (fuct.eq.1) fct=4
!!$              if (fuct.eq.2) fct =2
!!$              if (fuct.eq.3) fct =6
!!$              if (fuct.eq.4) fct =10


              !Domain size
              if (casemode.eq.0) then !RMSE comparisons only

                 if (fct.lt.7)then ! Analytic functions

                    DO j = 1, DIM
                       par(j,1)=-2.0
                       par(j,2)=2.0
                    END DO

                 else if (Fct.eq.8) then ! Truss design

                    if (fctindx.gt.3.or.dim.ne.3) stop'Wrong Function Index or Dimension for this Truss test case. Please check'
                    !0=obj, 1,2,3= constraints

                    DO j = 1, DIM

                       par(j,1)=1.0d0
                       par(j,2)=3.0d0

                    end do


                 else if (Fct.eq.9) then ! Short Column
                    
                    if (fctindx.gt.2.or.dim.ne.3) stop'Wrong Function Index or Dimension for this Truss test case. Please check'
                    !0=obj, 1,2,3= constraints
                    
                    DO j = 1, DIM ! Figure out a best range
                       
                       par(j,1)=1.0d0
                       par(j,2)=3.0d0
                       
                    end do


                 else if (Fct.eq.8) then ! Cantilever beam problem

                    if (fctindx.gt.3.or.dim.ne.3) stop'Wrong Function Index or Dimension for this Truss test case. Please check'
                    !0=obj, 1,2,3= constraints

                    DO j = 1, DIM

                       par(j,1)=1.0d0
                       par(j,2)=3.0d0

                    end do
                    

                 else if (fct.eq.20) then
                    
                    call cfdparams(par) !? on all threads?
                    
                 end if

              else if (Casemode.eq.1) then ! stats + rmse domain

                 ! statistics--> construct PC surrogate between mean and 3 SD's
                 do i=1,dim
                    par(i,1)=xavg(i)-3.0*xstd(i)                 
                    par(i,2)=xavg(i)+3.0*xstd(i)
                 end do

              else

                 print *, casemode              
                 stop 'Wrong Casemode'

              end if

              !              do fctindx=0,4,4  

              !             if (fctindx.eq.4) call system('cp MCSampCFD00.dat MCSampCFD04.dat')

              dyncyccnt=0

              do DIMPC =5,5 !order 5D requires 3003 terms

                 dyncyccnt=dyncyccnt+1

                 ! Initialize timer
                 !              if (id_proc.eq.0) call TimerInit()
                 !              if (id_proc.eq.0) call TimerStart('PC')

                 ! Get number of terms in the expansion
                 call combination(DIM+DIMPC,DIM,nterms)
                 ! Get number of points based on stat,solver,oversamp ratios

                 call getnpts(solver,stat,dim,nterms,OS,npts) 

                 if (dyncyccnt.eq.1) then

                    if(id_proc.eq.0) call sampdist(stat,DIM,DIMPC,ipar,par,makesamples,nterms,npts,fct,RN)

                 else

                    if (dynamics.eq.1) then 

                       nptstoaddpercyc=abs(npts-nptsold)

                       ! Trick to force a least number of points
!!$                       if (nptstoaddpercyc.le.dim) then
!!$                          npts=npts+dim
!!$                          nptstoaddpercyc=abs(npts-nptsold)
!!$                       end if

                       ! Call the dynamic sample routine


                       if (id_proc.eq.0) then
                          write(filenum,*)
                          write(filenum,*) '================================================='
                          write(filenum,*) '        Dynamic Point Selection                  '
                          write(filenum,*) '================================================='
                          write(filenum,*)
                          write(filenum,*)' >> NPTS to add this cycle:',nptstoaddpercyc 
                       end if

                       call dynsampdist(stat,DIM,DIMPC,ipar,par,makesamples,ntermsold,nterms,nptsold,npts,nptstoaddpercyc,fct,fpcb,gpcb,hpcb,xcof,RN)

                    else ! random points again call the same routine

                       if(id_proc.eq.0)  call sampdist(stat,DIM,DIMPC,ipar,par,makesamples,nterms,npts,fct,RN)

                    end if !dynamics

                 end if ! dyncyccnt

                 call MPI_Barrier(MPI_COMM_WORLD,ierr) ! All wait until master gets here

                 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                 !                 Non-parallel region                        !
                 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

                 if (id_proc.eq.0) then

                    !==================================================
                    !                SETUP MATRIX
                    !==================================================

                    write(filenum,*)
                    write(filenum,*) '================================================='
                    write(filenum,*) '        Matrix Setup & Pseudo-Inverse            '
                    write(filenum,*) '================================================='
                    write(filenum,*)

                    call setupmat(stat,dim,dimpc,npts,ipar,RN,numpc,xmat)

                    write(filenum,*)' >> Number of columns in the matrix (nterms):',nterms 
                    write(filenum,*)' >> Number of rows    in the matrix  (numpc):',numpc
                    if (solver.eq.1)then
                       write(filenum,*)'    >> Solver : LU Decomposition'
                    else if (solver.eq.2) then
                       write(filenum,*)'    >> Solver : SV Decomposition'
                    else
                       stop'Invalid Solver'
                    end if

                    !==============================================================
                    ! Decomposing the matrix to enable to use it for multiple RHS's
                    !==============================================================
                    write(filenum,*)'    >> Finding Pseudo-Inverse of the matrix . . . '

                    call pseudoinv(solver,numpc,nterms,xmat,indx,W,V)

                    write(filenum,*)
                    write(filenum,*) '================================================='
                    write(filenum,*) '                RHS Setup                        '
                    write(filenum,*) '================================================='
                    write(filenum,*)

                    ! Evaluates cost and grad and stores them as arrays fpcb(maxpts),gpcb(dim,maxpts)

                    if (dynamics.eq.0.or.dyncyccnt.eq.1) then
                       call gather_costfn(dim,fct,stat,npts,RN,fpcb,gpcb,hpcb)
                       ! for further dyncyccnts they are automatically computed within dynsamp routine
                    end if

                    !=======================================================
                    ! Setting up RHS
                    !=======================================================

                    write(filenum,*)'>> Cost Function Number :',fct
                    if (stat.eq.0) write(filenum,*)'    >> Setting up the RHS with Function evals . . .'
                    if (stat.eq.1) write(filenum,*)'    >> Setting up the RHS with Func+Grad evals . . .'
                    if (stat.eq.2) write(filenum,*)'    >> Setting up the RHS with Func+Grad+Hess evals . . .'

                    ! Put fpcb,gpcb,hpcb all into one rhsF
                    call setupRHS(stat,dim,npts,fpcb,gpcb,hpcb,numpc,rhsF)

                    write(filenum,*)'    >> Number of RHS entries :',numpc ! Atleast equal to number of coeffients
                    write(filenum,*)
                    do i=1,numpc
                       write(filenum,*)i, rhsF(i)
                    end do

                    write(filenum,*)
                    write(filenum,*) '================================================='
                    write(filenum,*) '            Solve Linear System                  '
                    write(filenum,*) '================================================='
                    write(filenum,*)

                    !=======================================================
                    ! Backsubstitute for any RHS (like a black box approach)
                    !=======================================================

                    write(filenum,*)'>> Solving for the coefficients . . . '
                    call backsub(solver,xmat,numpc,nterms,indx,W,V,rhsF,XCOF)
                    write(filenum,*)'>> Coeffs are succesfully determined . . .'
                    write(filenum,*)'>> Number of coefficients : ',nterms
                    write(filenum,*)

                 end if !master 

                 !-------------------------- Info Sharing -------------------------------!

                 call MPI_Barrier(MPI_COMM_WORLD,ierr) ! All wait until master gets the coeffs
                 call MPI_BCAST(npts,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
                 call MPI_BCAST(numpc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                 call MPI_BCAST(nterms,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                 call MPI_BCAST(xcof,maxtrm,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

                 !-------------------------- End  Sharing -------------------------------!

                 !              if (id_proc.eq.0) call TimerStop('PC') 
                 !              if (id_proc.eq.0) call TimerReport(ttotal)


                 !++++++++++++++++
                 !Parallel Region!
                 !++++++++++++++++

                 !================================================================
                 ! Post processing--> RMSE, Output statistics,Tecplot
                 !================================================================
                 if(id_proc.eq.0) then

                    write(filenum,*)
                    write(filenum,*) '================================================='
                    write(filenum,*) '            MONTE-CARLO on Surrogate             '
                    write(filenum,*) '================================================='
                    write(filenum,*)
                 end if

                 call montecarlo(stat,fct,DIM,dimpc,nterms,npts,ipar,xcof)

                 call MPI_Barrier(MPI_COMM_WORLD,ierr)

                 !=======================================================
                 ! Calculate RMSE
                 !=======================================================

                 if(id_proc.eq.0) then
                    write(filenum,*)
                    write(filenum,*) '================================================='
                    write(filenum,*) '             RMSE on Surrogate                   '
                    write(filenum,*) '================================================='
                    write(filenum,*)
                 end if
                 call MPI_Barrier(MPI_COMM_WORLD,ierr)
                 if(casemode.eq.0)   call RMSE_Higher(stat,dim,fct,npts,dimPC,ipar,par,xcof)
                 call MPI_Barrier(MPI_COMM_WORLD,ierr)

                 !================================================================
                 ! Tecplot output to file
                 !================================================================

                 if(id_proc.eq.0) then
                    if(dim.le.2) then
                       write(filenum,*)
                       write(filenum,*) '================================================='
                       write(filenum,*) '             Tecplot Output                      '
                       write(filenum,*) '================================================='
                       write(filenum,*)
                       write(filenum,*)'>> Writing Tecplot output to file . . .'

                       call tecplot(dim,dimpc,ipar,par,fct,npts,xcof) 
                    end if
                 end if

                 nptsold=npts
                 ntermsold=nterms

              enddo !order

              !           end do !fct indx

           enddo ! fct

        enddo ! F or FG (Stat)

     end do !Oversamp loop

  end do !dynamics loop

  if (id_proc.eq.0) then
     write(filenum,*)
     write(filenum,*)'>> Program terminated successfully'
     write(filenum,*) 
  end if

!!$  if (id_proc.eq.0) then
!!$     print *, fmean,fmeanprime(1:dim)
!!$     print *,fvar,fvarprime(1:dim)
!!$  end if

  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  call stop_all

end program main

!+======================================================================

function PHI(dim,x1,x2,scal)

  implicit none
  integer :: dim,j
  real*8 :: x1(dim),x2(dim),phi,x(dim),xnorm,tmp,scal

  x(:)=x1(:)-x2(:)

  xnorm=0.0
  do j=1,dim
     xnorm=xnorm+x(j)**2
  end do
  xnorm=SQRT(xnorm)/scal

  tmp=max(1.0-xnorm,0.0)
  phi=xnorm**3


  RETURN
end function phi



SUBROUTINE OTHPL(KF,N,X,PL,DPL)
  !
  !       ==========================================================
  !       Purpose: Compute orthogonal polynomials: Hn(x) or Ln(x),
  !                or Pn(x), and their derivatives
  !       Input :  KF --- Function code
  !                       KF=1 for Legendre polynomials Pn(x)
  !                       KF=2 for Hermite polynomial Hn(x)
  !                       KF=5 for Laguerre polynomial Ln(x)
  !                n ---  Order of orthogonal polynomials
  !                x ---  Argument of orthogonal polynomials
  !       Output:  PL(n) --- Hn(x) or Ln(x) or Pn(x)
  !                DPL(n)--- Hn'(x) or Ln'(x) or Pn'(x)
  !       =========================================================
  !
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  DIMENSION PL(0:N),DPL(0:N)

  Y0=1.0D0 
  DY0=0.0D0 
  PL(0)=1.0D0      
  DPL(0)=0.0D0

  IF (KF.EQ.2) THEN
     Y1=X
     DY1=1.0D0
     PL(1)=X
     DPL(1)=1.0D0
  ELSE IF (KF.EQ.5) THEN
     Y1=1.0D0-X
     DY1=-1.0D0
     PL(1)=1.0D0-X
     DPL(1)=-1.0D0
  ELSE IF (KF.EQ.1) THEN
     Y1=X
     PL(1)=X
     DPL(1)=1.0D0
  ELSE
     PRINT *, 'Invalid choice of orthogonal polynomial'
     STOP
  ENDIF

  DO K=2,N

     IF (KF.EQ.2) THEN
        A=1.0D0
        B=0.0D0
        C=(K-1.0D0)
     ELSE IF (KF.EQ.5) THEN
        A=-1.0D0/K
        B=2.0D0+A
        C=1.0D0+A
     ELSE IF (KF.EQ.1) THEN
        A=2.0D0-1.0D0/K
        B=0.D0
        C=1.0D0-1.0D0/K
     ENDIF

     YN=(A*X+B)*Y1-C*Y0
     PL(K)=YN

     IF (KF.GE.2) THEN
        DYN=A*Y1+(A*X+B)*DY1-C*DY0 
        DPL(K)=DYN
        DY0=DY1
        DY1=DYN
     ELSE
        IF (DABS(X).EQ.1.0D0) THEN
           DPL(K)=0.5D0*X**(K+1)*K*(K+1.0D0)
        ELSE
           DPL(K)=K*(Y1-X*YN)/(1.0D0-X*X)
        ENDIF
     ENDIF

     Y0=Y1
     Y1=YN       

  END DO

  RETURN
END SUBROUTINE OTHPL


subroutine i_to_s(intval,s)

  use dimpce, only: OS,fctindx
  implicit none

  integer idig,intval,ipos,ival,i
  character (len=*) s


  s(:) = ' '

  ival = intval

  !  Working from right to left, strip off the digits of the integer
  !  and place them into S(1:len ( s )).
  !
  ipos = len(s)

  do while ( ival /= 0 )

     idig = mod( ival, 10 )
     ival = ival / 10
     s(ipos:ipos) = char(idig + 48 )
     ipos = ipos - 1

  end do
  !
  !  Fill the empties with zeroes.
  !
  do i = 1, ipos
     s(i:i) = '0'
  end do

  return

end subroutine i_to_s

!=================================


subroutine getfilename(dim,fct,dimpc,stat,casemode,filename)
  use dimpce,only:OS,fctindx,dynamics
  implicit none  
  character*2 :: dimnumber,fctnumber,ordnumber,OSnumber,fctindxnumber
  integer ::lenc,fct,dim,stat,casemode,dimpc
  character*60 :: filename

  call i_to_s(fct,fctnumber)
  call i_to_s(dim,dimnumber)
  call i_to_s(dimpc,ordnumber)
  call i_to_s(OS,OSnumber)
  call i_to_s(fctindx,fctindxnumber)


  if (casemode.eq.0) then

     filename='norm/errornormdim'
     lenc=12+5

  else if (Casemode.eq.1) then

     filename='stats/statsdim'
     lenc=8+6

  else if (casemode.eq.55)then  !Storing samples

     filename='samp/sampdim'
     lenc=7+5
  else

     stop'Wrong flag in case mode'  

  end if

  filename(lenc+1:lenc+2)=dimnumber
  lenc=lenc+2

  if (casemode.ne.55)then
     filename(lenc+1:lenc+3)='fct'
     lenc=lenc+3
     filename(lenc+1:lenc+2)=fctnumber
     lenc=lenc+2
     if (fct.eq.20) then
        filename(lenc+1:lenc+1)='i'
        lenc=lenc+1
        filename(lenc+1:lenc+2)=fctindxnumber
        lenc=lenc+2
     end if
  else
     filename(lenc+1:lenc+3)='ord'
     lenc=lenc+3
     filename(lenc+1:lenc+2)=ordnumber
     lenc=lenc+2
  end if


  if (OS.gt.1) then

     filename(lenc+1:lenc+3)='Rat'
     lenc=lenc+3
     filename(lenc+1:lenc+2)=OSnumber
     lenc=lenc+2

  else 


  end if


  if (stat.eq.2) then 
     filename(lenc+1:lenc+3)='FGH'
     lenc=lenc+3    
  else if (stat.eq.1) then
     filename(lenc+1:lenc+2)='FG'
     lenc=lenc+2 
  else if (stat.eq.0) then
     filename(lenc+1:lenc+1)='F'
     lenc=lenc+1 
  else
     print *, 'Wrong value in hstat'
     stop'Wrong value in stat'
  end if
  if  (dynamics.eq.1) then
     filename(lenc+1:lenc+6)='mirdyn'
     lenc=lenc+6
  end if
  return
end subroutine getfilename

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!$
!!$subroutine evalPC(ndim,dimPC,ipar,xcof,x,yhat)
!!$
!!$  implicit none
!!$  include "collsub.h"
!!$
!!$  integer :: NDIM,DIMPC,nterms,mreg(maxdat,ndim) ! multi-index variables
!!$  integer :: k,j,jj,kk
!!$
!!$  integer ::  ipar(MAXVAR)
!!$  double precision ::  PL(NDIM,0:DIMPC),DPL(NDIM,0:DIMPC),DDPL(NDIM,0:dimpc)
!!$  double precision ::  PL1(NDIM,0:DIMPC),DPL1(NDIM,0:DIMPC),DDPL1(NDIM,0:dimpc)
!!$
!!$  double precision :: x(ndim),yhat
!!$  double precision :: xcof(MAXTRM),xcoftmp
!!$  integer::ipartmp
!!$  real*8::  xtmp,PLtmp(0:dimpc),DPLtmp(0:dimpc),DDPLtmp(0:dimpc)
!!$
!!$  call multidx(MAXDAT,NDIM,DIMPC,mreg,nterms) ! get multiindex notation for tensor procduct
!!$
!!$  !Initialize for safety
!!$
!!$  dpltmp=0.0d0
!!$  ddpltmp=0.0d0
!!$  pltmp=0.0d0
!!$
!!$  yhat = 0.0d0
!!$  do kk=1,nterms
!!$     xcoftmp=xcof(kk)
!!$     do jj=1,nDIM 
!!$
!!$        ipartmp=ipar(jj) ! Normal or uniform
!!$ 
!!$        xtmp=x(jj) !location to evalutuate the orthogonal polynomials
!!$
!!$        call ortho(ipartmp,DIMPC,xtmp,PLtmp,dpltmp,ddpltmp) !get values and derivatives
!!$
!!$        ! Store it in the way it is needed
!!$
!!$        PL(jj,0:dimpc)=pltmp(0:dimpc) 
!!$        DPL(jj,0:dimpc)=Dpltmp(0:dimpc)
!!$        DDPL(jj,0:dimpc)=DDpltmp(0:dimpc)
!!$
!!$        xcoftmp=xcoftmp*PL(jj,mreg(kk,jj)) ! the derivatives are not used at all
!!$
!!$     enddo
!!$     yhat = yhat + xcoftmp ! PC value
!!$  end do
!!$
!!$  return
!!$
!!$end subroutine evalPC

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine ortho(dist,dimpc,x,pl,dpl,ddpl)
  implicit none
  integer :: dist,dimpc,i
  real*8:: x,pl(0:dimpc),dpl(0:dimpc),ddpl(0:dimpc)


  if (dist.eq.1) then

     do i=0,dimpc
        call LEGEN(i,X,pl(i),dpl(i),Ddpl(i))        
     end do
     !           call LEGENDRE(dimpc,X,PL,DPL,ddpl) !something wrong in DDPL, DPL

  else if (dist.eq.2) then

     call hermite(dimpc,X,PL,DPL,ddpl,2)

  else
     stop'Not implemented yet'
  end if

end subroutine ortho

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine solvertype(stat,os,solver)
  implicit none

  integer,intent(in)::stat,os
  integer,intent(out)::solver
  if (stat.eq.0) then
     if (OS.eq.1) then
        solver=1
     else
        solver=2
     end if
  else
     solver =2
  end if

  return
end subroutine solvertype
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine cfdparams(par)
  implicit none
  include "collsub.h"
  integer,parameter::dim=2   
  integer::j

  real*8,intent(out)::  par(DIM,MAXPAR)

  double precision :: Initialmach, Finalmach, InitialAOA,FinalAOA

  InitialMach=0.5d0
  FinalMach=1.5d0
  InitialAOA=0.0d0  !degrees
  FinalAOA=5.0d0    !degrees  

  !Domain size
  DO j = 1, DIM

     if (j.eq.1) then
        par(j,1)=InitialAOA*4.0d0*atan(1.0)/180.0   ! in radian
        par(j,2)=FinalAOA*4.0d0*atan(1.0)/180.0   ! in radian
     else
        par(j,1)=InitialMach !initial mach number
        par(j,2)=FinalMach    !final mach number
     end if

  END DO

end subroutine cfdparams
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine pseudoinv(solver,numpc,nterms,xmat,indx,W,V)

  implicit none
  include "collsub.h"

  integer,intent(in)::solver,numpc,nterms
  integer,intent(out)::indx(MAXDAT)
  real*8,intent(in)::xmat(MAXDAT,MAXDAT)
  real*8,intent(out)::W(MAXDAT),V(MAXDAT,MAXDAT)
  real*8::dd


  if (solver.eq.1) then
     call ludcmp(xmat,numpc,MAXDAT,indx,dd)
     W=0.0d0
     V=0.0d0
  else

     call svdcmp(xmat,numpc,nterms,MAXDAT,MAXDAT,W,V)     
     indx=0

  end if

end subroutine pseudoinv

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine backsub(solver,xmat,numpc,nterms,indx,W,V,rhsF,XCOF)
  use dimpce,only:filenum
  implicit none
  include "collsub.h"
  INTEGER::I,J,K

  integer,intent(in)::solver,numpc,nterms
  integer,intent(in)::indx(MAXDAT)
  real*8,intent(in)::xmat(MAXDAT,MAXDAT),rhsf(MAXDAT)
  real*8,intent(in)::W(MAXDAT),V(MAXDAT,MAXDAT)
  real*8,intent(out)::XCOF(MAXDAT)


  if (solver.eq.1) then

     call lubksb(xmat,nterms,MAXDAT,indx,rhsF) !tmp=rhsfunc
     xcof=rhsF

  else

     call svbksb(xmat,W,V,numpc,nterms,MAXDAT,MAXDAT,rhsF,xcof)

  end if

  do j=1,nterms
     write(filenum,*) j,xcof(j)
  enddo

end subroutine backsub

!Unused routines

SUBROUTINE HERM(N,X,Y,DY,D2Y)                                   
  !  *************************************************************
  !  *     COMPUTES THE VALUE OF THE HERMITE POLYNOMIAL OF DEGREE N            
  !  *     AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT               
  !  *     N  = DEGREE OF THE POLYNOMIAL                                       
  !  *     X  = POINT IN WHICH THE COMPUTATION IS PERFORMED                    
  !  *     Y  = VALUE OF THE POLYNOMIAL IN X                                   
  !  *     DY = VALUE OF THE FIRST DERIVATIVE IN X                             
  !  *     D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
  !  *************************************************************
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               

  Y   = 1.D0                                                     
  DY  = 0.D0                                                     
  D2Y = 0.D0                                                     
  IF (N .EQ. 0) RETURN                                              

  Y   = 2.D0*X                                                   
  DY  = 2.D0                                                     
  D2Y = 0.D0                                                     
  IF (N .EQ. 1) RETURN                                              

  YP=1.D0                                                        
  DO 1 K=2,N                                                        
     DK = DFLOAT(K-1)                                               
     YM  = Y                                                        
     Y   = 2.D0*X*Y-2.D0*DK*YP                                      
     YPM = YP                                                       
     YP  = YM                                                       
1    CONTINUE                                                          
     DN  = 2.D0*DFLOAT(N)                                           
     DNN = 2.D0*DFLOAT(N-1)                                         
     DY  = DN*YP                                                    
     D2Y = DN*DNN*YPM                                               

     RETURN                                                            

   END SUBROUTINE HERM

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   SUBROUTINE LEGEN(N,X,Y,DY,D2Y)                                   
     !  **************************************************************          
     !  *   COMPUTES THE VALUE OF THE LEGENDRE POLYNOMIAL OF DEGREE N           
     !  *   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT               
     !  *   N  = DEGREE OF THE POLYNOMIAL                                       
     !  *   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED                    
     !  *   Y  = VALUE OF THE POLYNOMIAL IN X                                   
     !  *   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
     !  *   D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
     !  **************************************************************          
     IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               

     Y   = 1.D0                                                     
     DY  = 0.D0                                                     
     D2Y = 0.D0                                                     
     IF (N .EQ. 0) RETURN                                              

     Y   = X                                                        
     DY  = 1.D0                                                     
     D2Y = 0.D0                                                     
     IF(N .EQ. 1) RETURN                                               

     YP   = 1.D0                                                    
     DYP  = 0.D0                                                    
     D2YP = 0.D0                                                    
     DO 1 I=2,N                                                        
        C1 = DFLOAT(I)                                                 
        C2 = 2.D0*C1-1.D0                                              
        C4 = C1-1.D0                                                   
        YM = Y                                                         
        Y  = (C2*X*Y-C4*YP)/C1                                         
        YP = YM                                                        
        DYM  = DY                                                      
        DY   = (C2*X*DY-C4*DYP+C2*YP)/C1                               
        DYP  = DYM                                                     
        D2YM = D2Y                                                     
        D2Y  = (C2*X*D2Y-C4*D2YP+2.D0*C2*DYP)/C1                       
        D2YP = D2YM                                                    
1       CONTINUE                                                          

        RETURN                                                            
        !  END DO

      End SUBROUTINE LEGEN
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      subroutine evalPC(ndim,dimPC,ipar,xcof,x,yhat,yhatprime,yhatdbleprime)

        implicit none
        include "collsub.h"

        integer :: NDIM,DIMPC,nterms,mreg(maxdat,ndim) ! multi-index variables
        integer :: k,j,jj,kk

        integer ::  ipar(MAXVAR)
        double precision ::  PL(NDIM,0:DIMPC),DPL(NDIM,0:DIMPC),DDPL(NDIM,0:dimpc)
        double precision ::  PL1(NDIM,0:DIMPC),DPL1(NDIM,0:DIMPC),DDPL1(NDIM,0:dimpc)

        double precision :: x(ndim),yhat
        double precision :: xcof(MAXTRM),xcoftmp
        integer::ipartmp
        real*8::  xtmp,PLtmp(0:dimpc),DPLtmp(0:dimpc),DDPLtmp(0:dimpc)

        real*8::yhatprime(ndim),yhatprimetmp
        real*8::yhatdbleprime(ndim,ndim),yhatdbleprimetmp

        call multidx(MAXDAT,NDIM,DIMPC,mreg,nterms) ! get multiindex notation for tensor procduct

        !  print *,xcof(1:nterms)
        dpltmp=0.0
        ddpltmp=0.0
        pltmp=0.0

        yhat = 0.0d0
        do kk=1,nterms
           xcoftmp=xcof(kk)
           do jj=1,nDIM 

              ipartmp=ipar(jj)
              xtmp=x(jj)
              call ortho(ipartmp,DIMPC,xtmp,PLtmp,dpltmp,ddpltmp)


              PL(jj,0:dimpc)=pltmp(0:dimpc)
              DPL(jj,:)=Dpltmp(:)
              DDPL(jj,:)=DDpltmp(:)

              xcoftmp=xcoftmp*PL(jj,mreg(kk,jj))
           enddo

           yhat = yhat + xcoftmp

        end do

        yhatprime=0.0d0 !Initialize
        xcoftmp=0.0d0   !Initialize

        ! Gradient values
        do kk=1,nDIM
           yhatprimetmp=0.0d0
           do j=1,nterms
              xcoftmp=xcof(j) 
              do k=1,nDIM
                 ipartmp=ipar(k)
                 xtmp=x(k)
                 !       print *, X
                 call ortho(ipartmp,DIMPC,xtmp,PLtmp,dpltmp,ddpltmp)
                 PL(k,0:dimpc)=pltmp(0:dimpc)
                 DPL(k,0:dimpc)=Dpltmp(0:dimpc)
                 !        print*,'pl:',pltmp(0:dimpc)
                 !        print*,'dpl:',dpltmp(0:dimpc)

                 !!           DDPL(k,:)=DDpltmp(:)

                 if (k.eq.kk) then
                    xcoftmp = xcoftmp*DPL(k,mreg(j,k))
                 else
                    xcoftmp = xcoftmp*PL(k,mreg(j,k))
                 end if

              end do
              !      print*,'xcof:',xcoftmp
              yhatprimetmp=yhatprimetmp+xcoftmp
           end   do
           yhatprime(kk)=yhatprimetmp
           !     print *,yhatprimetmp
        end do
        !print *, x,yhatprime 

        !stop
        yhatdbleprime=0.0d0 !Initialize
        xcoftmp=0.0d0 !Initialize
        ddpl=0.0d0
        dpl=0.0d0
        pl=0.0d0

        do jj=1,nDIM

           do kk=1,nDIM
              yhatDBLEprimetmp=0.0d0

              do j=1,nterms
                 xcoftmp=xcof(j) 

                 do k=1,nDIM !x1,x2,x3.....

                    ipartmp=ipar(k)
                    xtmp=x(k)
                    call ortho(ipartmp,DIMPC,xtmp,PLtmp,dpltmp,ddpltmp)
                    PL(k,0:dimpc)=pltmp(0:dimpc)
                    DPL(k,:)=Dpltmp(:)
                    DDPL(k,:)=DDpltmp(:)
!!$              print *, PL(k,0:dimpc)
!!$              print *, DPL(k,0:dimpc)
!!$              print *, ddPL(k,0:dimpc)

                    if (kk.eq.jj) then

                       if (k.eq.kk) then
                          XCOFTMP = XCOFTMP*DDPL(k,mreg(j,k))
                       else
                          XCOFTMP= XCOFTMP*PL(k,mreg(j,k))
                       end if

                    else if (kk.ne.jj)then !off diagonal elements

                       if (k.eq.kk.or.k.eq.jj) then
                          XCOFTMP = XCOFTMP*DPL(k,mreg(j,k))
                       else
                          XCOFTMP = XCOFTMP*PL(k,mreg(j,k))         
                       end if

                    end if
                 end do
                 yhatDBLEprimetmp=yhatDBLEprimetmp+xcoftmp

              end   do

              yhatdbleprime(kk,jj)=yhatdbleprimetmp
           end do
        end do

!!$  print *,'x:',x
!!$
!!$  do jj=1,ndim
!!$
!!$    write(filenum,*)(yhatdbleprime(kk,jj),kk=1,ndim)
!!$ end do



        return
      end subroutine evalPC
