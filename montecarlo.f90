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
  !  real*8,intent(out) ::fvar,fmean
  !  double precision, dimension(20),intent(out) :: fmeanprime,fvarprime
  
  !  real*8::fstd!,fmeanprime,fvarprime,fstdprime
  !  real*8 :: xavg(ndim),xvar(ndim),xstd(ndim)
  real*8 :: xavgtmp(ndim),xvartmp(ndim),xstdtmp(ndim)
  ! double precision, dimension(20)   :: fstdprime,xavgt,xstdt  
  
  character*60 :: filename
  !  double precision :: MNCf(MAXNMCS)
  !  double precision :: MNCx(NDIM,MAXNMCS)
  real*8::average(ndim)
  ! double precision :: HST(npdf,MAXNMCS),HSTglb(npdf,MAXNMCS)
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

!  real*8::fvtemp(maxnmcs)

  real*8::yhatprime(ndim),yhatprimetmp
  real*8::yhatdbleprime(ndim,ndim),yhatdbleprimetmp

  if(ndim.gt.20) stop'Max dimensions Exceeded. Change preallocated arrays'

  call multidx(MAXDAT,nDIM,DIMPC,mreg,nterms) ! get multiindex notation for tensor procduct

  open(10,file='MC.inp',form='formatted',status='old')
  read(10,*) !(xavg(i),i=1,ndim)
  read(10,*) !(xstd(i),i=1,ndim)     
  read(10,*)
  read(10,*)
  read(10,*) NMCS!,ndimtmp
  read(10,*) npdf
  read(10,*) readMCsamples
  read(10,*) evlfnc
  close(10)

!  print*, xavg,xstd

  allocate(MNCf(NMCS))
  allocate(MNCx(ndim,NMCS))
  allocate(HST(npdf,3))
  allocate(HSTglb(npdf,3))

  if (id_proc.eq.0) then

     write(filenum,*) 
     write(filenum,'(6x,a,i8)')'>> # of Monte-Carlo Samples = ',NMCS
  end if

  !Initialize variables for safety

  MNCf   = 0.d0
  MNCx   = 0.d0
  HST    = 0.d0
  HSTglb = 0.d0
  ictglb=0

  if(id_proc.eq.0)then
       if (readMcsamples.eq.1) then

        if (fct.lt.20) then
           write(filenum,*)'     >> Reading ',nmcs,' Monte-Carlo samples from MCsamp.dat'
           write(filenum,*)

           open(120,file='MCsamp.dat',form='formatted',status='old')
           read(120,'(2i8)') NMCStmp, ndimtmp
           read(120,*) (xavgtmp(i),i=1,ndim)
           read(120,*) (xstdtmp(i),i=1,ndim)
           do j=1,NMCStmp
              read(120,*) (MNCx(k,j),k=1,ndim)
           end do
           close(120)
          
        end if

        if (fct.eq.20) then
           if (fctindx.eq.0) then
              write(filenum,*)'     >> Reading ',nmcs,' Monte-Carlo samples from MCSampCFD00.dat'
              open(120,file='MCSampCFD00.dat',form='formatted',status='old')
           else if (fctindx.eq.4) then
              write(filenum,*)'     >> Reading ',nmcs,' Monte-Carlo samples from MCSampCFD04.dat'
              open(120,file='MCSampCFD04.dat',form='formatted',status='old')
           end if

           read(120,'(3i8)') NMCStmp,ndimtmp,fctindxtmp
           read(120,*) (xavgtmp(i),i=1,ndim)
           read(120,*) (xstdtmp(i),i=1,ndim)
           do j=1,NMCStmp
              read(120,*) (MNCx(k,j),k=1,ndim)!,fvtemp(j)
           end do
           close(120)
       
        end if

        ! Check if input and available sample data are clean and the same
        do i=1,ndim
           if (xavgtmp(i).ne.xavg(i)) then
              print *, xavg,xavgtmp
              print *, xstd,xstdtmp
           stop'xavgtmp should be equal to xavg when reading MC samples'
        end if
        if (xstdtmp(i).ne.xstd(i)) stop'xstdtmp should be equal to xstd when reading MC samples'
        end do
        if (ndimtmp.ne.ndim) stop'ndimtmp should be equal to ndim when reading MC samples'
        if (nmcstmp.ne.nmcs) stop'NMCStmp should be equal to NMCS when reading MC samples'

        
     else !make samples and write to new file

        call get_seed(seed)
        call latin_random(ndim,NMCS,seed,MNCx)

        ! Normally distribute with given mean and standard deviation
        do j = 1, NMCS
           do k=1,ndim 
              MNCx(k,j)=xavg(k)+dinvnorm(MNCx(k,j))*xstd(k)
           end do
        end do

        if (fct.lt.20) then

           ! write monte carlo samples to file
           write(filenum,*)'     >> Writing ',nmcs,' Monte-Carlo samples to MCsamp.dat'
           write(filenum,*)
           open(120,file='MCsamp.dat',form='formatted',status='unknown')
           write(120,'(2i8)') NMCS,ndim
           write(120,*) (xavg(i),i=1,ndim)
           write(120,*) (xstd(i),i=1,ndim)
           do j=1,NMCS
              write(120,*) (MNCx(k,j),k=1,ndim)
           end do
           close(120)

        else

           if (fctindx.eq.0) then
              write(filenum,*)'     >> Writing ',nmcs,' Monte-Carlo samples to MCSampCFD00.dat'
              open(121,file='MCSampCFD00.dat',form='formatted',status='unknown')
           else if (fctindx.eq.4) then
              write(filenum,*)'     >> Writing ',nmcs,' Monte-Carlo samples to MCSampCFD04.dat'
              open(121,file='MCSampCFD04.dat',form='formatted',status='unknown')
           end if
           
           write(121,'(3i8)') NMCS,ndim,fctindx
           write(121,*) (xavg(i),i=1,ndim)
           write(121,*) (xstd(i),i=1,ndim)
           do j=1,NMCS
              write(121,*) (MNCx(k,j),k=1,ndim)!,fvtemp(j)
           end do
           close(121)

        end if ! which fct
       
     end if ! read or create MC samples
!!$
!!$     do j = 1, NMCS
!!$        do k=1,ndim 
!!$           MNCx(k,j) = (MNCx(k,j)-DS(1,k))/(DS(2,k)-DS(1,k))
!!$        end do
!!$     end do

!!$     average(:)=0.0d0
!!$     ! Map samples to PC space
!!$     do j=1,NMCS
!!$!        do k=1,ndim
!!$           !         MNCx(k,j) =(MNCx(k,j)-DS(1,k))/(DS(2,k)-DS(1,k))
!!$!                           print *,DS(1,k),(DS(2,k)-DS(1,k))
!!$ !       end do
!!$        !stop
!!$        average(:)=average(:)+ MNCx(:,j)
!!$     end do
!!$     average(:)=dble(average(:))/dble(NMCS)
!!$     print *, average
!!$     print *,'AVG:',xavg
!!$     print*,'STD:',xstd
!!$     stop

  end if !id_proc.eq.0

  !  call MPI_BCAST(MNCx(1,1),NMCS*ndim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  ! call MPI_BCAST(MNCx(1,1),NMCS*1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  !  do jj=1,ndim


  call MPI_BCAST(MNCx,NMCS*ndim,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  !end do

  do k=1,1 !ftgt,ftgts

     idec = dble(NMCS)/dble(num_proc)
     is   = idec*id_proc + 1
     ie   = idec*(id_proc+1)
     if(id_proc.eq.num_proc-1)ie = NMCS 
     write(filenum,'(11x,a,2i10,a,i3)') '>> [',is,ie,' ] for Processor',id_proc

     ymin = 1.d10
     ymax =-1.d10
     xmin = 1.d10
     xmax =-1.d10
     MCm  = 0.d0
     MCd  = 0.d0
     MCmprime(:)=0.d0
     MCdprime(:)=0.d0
     MCmdbleprime(:,:)=0.d0
     MCddbleprime(:,:)=0.d0

     ict  = 0

!!$     if (fct.eq.20.and.id_proc.eq.0) then
!!$        open(10,file='MCCFDvalues.dat',form='formatted',status='unknown')
!!$        write(120,*) NMCS, ndim,fctindx
!!$        write(120,*) (xavg(i),i=1,ndim)
!!$        write(120,*) (xstd(i),i=1,ndim)
!!$     end if
     
     do i=is,ie ! Main Loop for MC

        x(:)    = MNCx(:,i)                           

        do j=1,ndim
           xmin(j) = min(xmin(j),x(j))
           xmax(j) = max(xmax(j),x(j))
        end do
!!$
!!$        
!!$        do j=1,nDIM
!!$           scal=DS(2,j)-DS(1,j)
!!$           x(j)=x(j)*scal+DS(1,j)
!!$        end do
        
        call evalPC(ndim,dimPC,ipar,xcof,x,yhat,yhatprime,yhatdbleprime) ! Evaluates PC at x and returns yhat

        !!  call evalcostf(stat,ndim,fct,x,fv,gv,hv,DS)  !Evaluate cost function


!                print *, 'x:',x,'PC:',yhat!,'Ex:',fv ,id_proc  
        
        ict    = ict + 1
        MCm    = MCm + yhat    ! for mean
        MCd    = MCd + yhat**2 ! for variance 
        
      
        do j=1,ndim
           MCmprime(j) = MCmprime(j) + yhatprime(j)  !for derivative of mean
           MCdprime(j) = MCdprime(j) + 2.0*yhat*yhatprime(j)  !for derivative of variance
        end do


        ! For hessian of mean and variance
        
        
        do jj=1,ndim

           do kk=1,ndim
              

              MCmdbleprime(jj,kk)=  MCmdbleprime(jj,kk) + yhatdbleprime(jj,kk) ! Mean
              
              MCddbleprime(jj,kk)=  MCddbleprime(jj,kk) + 2.0*(yhat*yhatdbleprime(jj,kk) + yhatprime(kk)**2)
              
           end do
        end do

        MNCf(i) = yhat
        if(yhat.ge.ymax)then
           ymax = yhat
        end if
        if(yhat.le.ymin)then
           ymin = yhat
        end if

     end do ! main loop for monte carlo
     
     ! Information Sharing
     do id=0,num_proc-1
        is   = idec*id + 1
        ie   = idec*(id+1)
        if(id.eq.num_proc-1)ie = NMCS
        call MPI_BCAST(MNCf(is),ie-is+1,MPI_DOUBLE_PRECISION,id,MPI_COMM_WORLD,ierr)
     end do
     call MPI_ALLREDUCE(ict,ictglb,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MCm,MCmglb,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MCd,MCdglb,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(ymin,yminglb,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(ymax,ymaxglb,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
     do i=1,ndim
        call MPI_ALLREDUCE(xmin(i),xminglb(i),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(xmax(i),xmaxglb(i),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MCmprime(i),MCmprimeglb(i),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

        call MPI_ALLREDUCE(MCdprime(i),MCdprimeglb(i),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

     end do


     ! Hessian
     do i=1,ndim
        do j=1,ndim
           call MPI_ALLREDUCE(MCmdbleprime(i,j),MCmdbleprimeglb(i,j),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
           call MPI_ALLREDUCE(MCddbleprime(i,j),MCddbleprimeglb(i,j),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

        end do
     end do
     


     if(ictglb.ne.NMCS) then
        print *,ictglb,nmcs
        stop'ictglb.ne.NMCS'
     end if

     fmean = MCmglb / dble(ictglb)
     fvar = MCdglb / dble(ictglb) - fmean**2  

     !Derivative
     do i=1,ndim
        fmeanprime(i) = MCmprimeglb(i) / dble(ictglb)
        fvarprime(i) =  MCdprimeglb(i) / dble(ictglb) - 2.0 * fmean * fmeanprime(i)
     end do

     !Hessian
     
     do i=1,ndim
        do j=1,ndim

           fmeandbleprime(i,j)=MCmdbleprimeglb(i,j)/dble(ictglb)
           
           fvardbleprime(i,j)=MCddbleprimeglb(i,j)/dble(ictglb) -  2.0*(fmean*fmeandbleprime(i,j) + fmeanprime(j)**2)

        end do
     end do

     ! Histogram

     yhmin = yminglb
     yhmax = ymaxglb
     width = (yhmax-yhmin)/dble(npdf)
     do i=1,npdf
        HST(i,1) = yhmin + dble(i-1)*width
        HST(i,2) = yhmin + dble(i  )*width
     end do
     HST(npdf,2) = yhmax
     HST(:,3)    = 0.d0
     MCd = 0.d0
     is   = idec*id_proc + 1
     ie   = idec*(id_proc+1)
     if(id_proc.eq.num_proc-1)ie = NMCS
     do 250 i=is,ie
        if(MNCf(i).le.yhmin)then
           HST(1,   3) = HST(1,   3) + 1.d0
        else if(MNCf(i).ge.yhmax)then
           HST(npdf,3) = HST(npdf,3) + 1.d0
        else
           do 260 j=1,npdf
              if(MNCf(i).ge.HST(j,1).and.MNCf(i).lt.HST(j,2))then
                 HST(j,3) = HST(j,3) + 1.d0
                 go to 250
              end if
260           continue
              write(filenum,'(2i8,3e15.5)')id_proc,i,yhmin,MNCf(i),yhmax
              stop'without the border of histogram?'
           end if

250        continue

           do i=1,npdf
              call MPI_ALLREDUCE(HST(i,3),hglb,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!!$              hglb=0.0d0 !remove it
              HSTglb(i,1) = HST(i,1)
              HSTglb(i,2) = HST(i,2)
              HSTglb(i,3) = hglb
           end do
           ! output
           if(id_proc.eq.0)then
              
              fstd=sqrt(fvar)
              fmeanprime(1:ndim) = fmeanprime(1:ndim) !/(DS(2,1:ndim)-DS(1,1:ndim))
              fvarprime= fvarprime(1:ndim) !/(DS(2,1:ndim)-DS(1,1:ndim))
              write(filenum,*)
              write(filenum,'(6x,a,2i4)')'>> Fct,  Function index:',fct,fctindx
              write(filenum,*)
              write(filenum,'(6x,a,6e15.5)')'>> xavg ',xavg(1:ndim)
              write(filenum,'(6x,a,6e15.5)')'>> xstd ',xstd(1:ndim)  
              write(filenum,*)
              write(filenum,'(6x,a,3e15.5)')'>> Mean, Var, SD  = ',fmean,fvar,fstd
              write(filenum,*)
              write(filenum,'(6x,a,20e15.5)')'>> Meanprime = ',fmeanprime(1:ndim)
              write(filenum,'(6x,a,20e15.5)')'>> Varprime = ',fvarprime(1:ndim)
              write(filenum,*)

              write(filenum,'(6x,a,20e15.5)')'>> Meandbleprime = '
              do i=1,ndim
                 write(filenum,'(8x,20e15.5)')fmeandbleprime(1:ndim,i)
              end do

              write(filenum,'(6x,a,20e15.5)')'>> Vardbleprime = '
              do i=1,ndim
                 write(filenum,'(8x,20e15.5)')fvardbleprime(1:ndim,i)
              end do
              write(filenum,*)
              write(filenum,'(6x,a,99f8.3)')'>> Range of X(min) = ',(xminglb(i),i=1,ndim)
              write(filenum,'(6x,a,99f8.3)')'>> Range of X(max) = ',(xmaxglb(i),i=1,ndim)
              write(filenum,'(6x,a,2e15.5)')'>> Ymax/min = ',ymaxglb,yminglb
              write(filenum,'(6x,a,2(e15.5,a))')'>> Output Histogram  [',yhmin,' :',yhmax,' ]'
!              write(filenum,'(6x,a,i8)')'>> ictglb:',ictglb
              write(filenum,*)
              pdf = 0.d0
              open(10,file='HISTG.dat',form='formatted',status='unknown')
              do i=1,npdf
                 pdf = pdf + HSTglb(i,3)/dble(NMCS)
                 write(10,'(4e15.5)')(HSTglb(i,1)+HSTglb(i,2))*0.5d0,HSTglb(i,3),HSTglb(i,3)/dble(NMCS),pdf
              end do
              close(10)
           end if

        end do! dummy loop for function index (k)

        if (id_proc.eq.0) then
           
           if (fct.eq.20) then

              
              if (fctindx.eq.0) then
                 open(10,file='MCCFDvalues00.dat',form='formatted',status='unknown')
              else if (fctindx.eq.4) then
                 open(10,file='MCCFDvalues04.dat',form='formatted',status='unknown')
              end if

           
           end if

           Javg(:)=0.0
           Jvar(:)=0.0
           Jstd(:)=0.0

           do i=1,NMCS

!!$              ! Real function evaluation
!!$              if (fct.lt.20) then
!!$                 call evalcostf(stat,ndim,fct,MNCX(:,i),freal,df,d2f,DS) 
!!$
!!$              else
!!$                 read(10,*) (x(kk),kk=1,ndim),freal               
!!$              end if
!!$              
              ! Real function evaluation

              if (fct.eq.20) then   
                 
                 if (evlfnc.eq.1) then

                    if (i.eq.1) write(filenum,'(a,i4,i4)')'>> Doing CFD Evaluation for MC - Fct, Fctindx',fct,fctindx
                    call evalcostf(stat,ndim,fct,MNCX(:,i),freal,df,d2f,DS) 
                    write(10,*) (MNCx(kk,i),kk=1,ndim),freal 

                 else !read

                    if (i.eq.1) write(filenum,'(a,i4,i4)')'>> Reading CFD database for MC - Fct, Fctindx',fct,fctindx
                    read(10,*) (MNCx(kk,i),kk=1,ndim),freal               
                 end if

              else !other functions not a problem

                 call evalcostf(stat,ndim,fct,MNCX(:,i),freal,df,d2f,DS) 

              end if

          
                 Javg(1)=Javg(1)+freal
                 Jvar(1)=Jvar(1)+freal**2


              if (stat.gt.0) then
    
                 Ddiff(:)=MNCx(:,i)-xbar(:)
                 
                 ! Linear extrapolation        
                 fobjlin=freal
                 do k=1,ndim
                    fobjlin=fobjlin+df(k)*Ddiff(k)
                 end do
                 
                 Javg(2)=Javg(2)+fobjlin
                 Jvar(2)=Jvar(2)+fobjlin**2


                 if (stat.eq.2) then

                    ! Quadratic extrapolation
                    fobjquad=fobjlin
                    do j=1,ndim
                       do k=1,ndim
                          fobjquad=fobjquad+0.5*d2f(j,k)*Ddiff(j)*Ddiff(k)
                       end do
                    end do
                    
                    Javg(3)=Javg(3)+fobjquad
                    Jvar(3)=Jvar(3)+fobjquad**2

                 end if

              end if

           end do

           if (fct.eq.20)  close(10)

           Javg(:)=Javg(:)/real(NMCS)      
           Jvar(:)=Jvar(:)/real(NMCS)-Javg(:)**2
           Jstd(:)=sqrt(Jvar(:))

!!$           if (fct.eq.20) then
!!$              if (xstd(1).eq.0.01) then
!!$                 Javg(1)= 0.25906E+00    
!!$                 Jvar(1)= 0.31653E-01
!!$              else if (xstd(1).eq.0.0075) then
!!$                 Javg(1)= 0.26361E+00    
!!$                 Jvar(1)= 0.18123E-01
!!$              else if (xstd(1).eq.0.005) then
!!$                 Javg(1)= 0.26671E+00    
!!$                 Jvar(1)= 0.86688E-02
!!$              else if (xstd(1).eq.0.0025) then
!!$                 Javg(1)= 0.26819E+00    
!!$                 Jvar(1)= 0.20390E-02
!!$              end if
!!$           end if


           write(filenum,'(6x,a,3e15.5)')'Real : Mean, Var, SD',Javg(1),Jvar(1),Jstd(1)
!           write(filenum,*)
           
           write(filenum,'(6x,a,3e15.5)')'PC   : Mean, Var, SD',fmean,fvar,fstd
 !          write(filenum,*)

           write(filenum,'(6x,a,3e15.5)')'Error: Mean, Var, SD',abs(fmean-Javg(1)),abs(fvar-Jvar(1)),abs(fstd-Jstd(1))

  !         write(filenum,*)

           ! Stats output
           
           call getfilename(ndim,fct,dimpc,stat,1,filename)
           if (dyncyccnt.eq.1)  then
              open(94,file=filename,form='formatted',status='unknown')!,position='append')
           else
              open(94,file=filename,form='formatted',status='old',position='append')
           end if
           if (dyncyccnt.eq.1) then ! Write the title line
              write(94,'(11a)') '      dimpc ',' npts   ','RealAVG   ','       RealVAR   ', '   RealSTD   ','           PcAVG   ', '   PcVAR   ', '             PcSTD  ','        ErrAVG   ','     ErrVAR', '           ErrSTD'

!!              if (stat.eq.1) write(94,'(a131)') 'Nhs      RealAVG        RealVAR         KrigAVG         KrigVAR        ErrorAVG       ErrorVAR      MM1AVG         MM1VAR          LinAVG           LinVAR'
 !!             if (stat.eq.2) write(94,'(a196)') 'Nhs      RealAVG        RealVAR         KrigAVG         KrigVAR        ErrorAVG       ErrorVAR   MM1AVG         MM1VAR          LinAVG           LinVAR          MM2AVG         MM2VAR          QuadAVG        QuadVAR'
           end if

           write(94,'(2i8,9e16.8)') dimpc,npts,Javg(1),Jvar(1),Jstd(1),fmean,fvar,fstd,abs(fmean-Javg(1)),abs(fvar-Jvar(1)),abs(fstd-Jstd(1))
!!$           

!!$           if (hstat.gt.0) then
!!$
!!$              write(filenum,'(6x,a,2e15.5)')'MM1: Mean and Variance',muy1,sigmay1
!!$              write(filenum,*)
!!$
!!$              write(filenum,'(6x,a,2e15.5)')'Lin: Mean and Variance',Javg(2),Jvar(2)
!!$              write(filenum,*)
!!$
!!$              if (hstat.ne.3) write(94,'(i4,8e16.8)') nhs,Javg(1),Jvar(1),fmean,fvar,abs(fmean-Javg(1)),abs(fvar-Jvar(1)),muy1,sigmay1,Javg(2),Jvar(2)
!!$
!!$              if (hstat.eq.3) then
!!$                 write(filenum,'(6x,a,2e15.5)')'MM2: Mean and Variance',muy2,sigmay2
!!$                 write(filenum,*)
!!$                 write(filenum,'(6x,a,2e15.5)')'Quad: Mean and Variance',Javg(3),Jvar(3)
!!$                 write(filenum,*)
!!$                 write(94,'(i4,12e16.8)') nhs,Javg(1),Jvar(1),fmean,fvar,abs(fmean-Javg(1)),abs(fvar-Jvar(1)),muy1,sigmay1,Javg(2),Jvar(2),muy2,sigmay2,Javg(3),Jvar(3)
!!$              end if

           close(94)
           
        end if

        deallocate(MNCf,MNCx) 
        deallocate(HST,HSTglb) 

      end subroutine montecarlo
!================================================
    subroutine find_st(pst,dx,st)
        implicit none
! find st which ensures the probability within [-dx:dx] is pst
        double precision, intent(in)  :: pst,dx
        double precision, intent(out) :: st
        integer :: iter
        double precision :: pout,s,ds,vv,vp,dv,sini

        if(pst.le.0.d0.or.pst.ge.1.d0)stop'pst in find_st'
        pout = (1.d0-pst)/2.d0
        sini = 1.d0
190     continue
        s    = sini
        ds   = sini*1.d-3
        iter = 0
200     continue
        iter = iter + 1
        call CDF(-1.d0*dx,0.d0,s,   vv)
        if(dabs(vv-pout).le.1.d-10)go to 210
        call CDF(-1.d0*dx,0.d0,s+ds,vp)
        dv = (vp-vv)/ds
!       write(filenum,'(5e15.5)')s,vv,pout,vv-pout,dv
        if(dv.eq.0.d0)stop'dv = 0.d0 in find_st'
        if(iter.ge.100)stop'iter>100 in find_st'
        s = s - (vv-pout)/dv
        if(s.le.0.d0)then
          sini = sini * 0.1d0
          go to 190
        end if
        go to 200
210     continue
!       write(filenum,'(4e15.5)')s,vv,pout,vv-pout
        st = s

        end subroutine find_st
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine find_x(mode,ran,xc,st,xout)
        implicit none
! find xout which ensures the CDF at xout is ran
! mode=0 : analytical CDF (fast but less robust?)
! mode=1 : numerical  CDF (time comsuming, but robust)
        integer, intent(in) :: mode
        double precision, intent(in)  :: ran,xc,st
        double precision, intent(out) :: xout
        integer :: iter 
        double precision :: x,vv,dv,vp,dx

        x    = xc
        dx   = 1.d-4
        iter = 0
200     continue
        iter = iter + 1
        if(mode.eq.0)then
          call CDF(x,xc,st,vv)
        else
          call CDF_Numerical(x,xc,st,vv)
        end if
        if(dabs(vv-ran).le.1.d-10)go to 210
        if(mode.eq.0)then
          call CDF(x+dx,xc,st,vp)
        else
          call CDF_Numerical(x+dx,xc,st,vp)
        end if
        dv = (vp-vv)/dx
!       write(filenum,'(4e15.5)')x,vv,vv-ran,dv
        if(dv.eq.0.d0)stop'dv=0 in find_x'
        if(iter.ge.100)stop'iter>100 in find_x'
        x = x - (vv-ran)/dv
        go to 200
210     continue
!       write(filenum,'(3e15.5)')x,vv,vv-ran
        xout = x

        end subroutine find_x
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine CDF(xin,xc,st,vout)
        implicit none
        double precision, intent(in)  :: xin,xc,st
        double precision, intent(out) :: vout
        double precision :: vtmp
!       vout = 0.5d0 * (1.d0 + erf( (xin-xc)/(st*dsqrt(2.d0)) ))
        call ERF_MINE1( (xin-xc)/(st*dsqrt(2.d0)), vtmp )
        vout = 0.5d0 * (1.d0 + vtmp)
        end subroutine CDF
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine DCDF(xin,xc,st,dvout)
        implicit none
        double precision, intent(in)  :: xin,xc,st
        double precision, intent(out) :: dvout
        double precision :: dvtmp
        call DERF_MINE( (xin-xc)/(st*dsqrt(2.d0)), dvtmp )
        dvout = 0.5d0*dvtmp/(st*dsqrt(2.d0))
        end subroutine DCDF
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine ERF_MINE1(xin,yout)
        implicit none
        double precision, intent(in)  :: xin
        double precision, intent(out) :: yout
        integer :: i,k,n
        double precision :: vsum,kai
! n is the order of Taylor
! Maybe accurate within the range of [-4:4] with n=100
        n = 100
        vsum = 0.d0
        do i=0,n
          kai = 1.d0
          do k=1,i
            kai = kai * (-1.d0) * xin**2 / dble(k)
          end do
          vsum = vsum + kai*xin/dble(2*i+1)
        end do
        yout = vsum*2.d0/(dsqrt(3.141592653589793238d0))

        if(yout.gt.1.d0)write(*,'(a,2e15.5)')'*ERF>1 ',xin,yout-1.d0
        end subroutine ERF_MINE1
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine DERF_MINE(xin,dyout)
        implicit none
        double precision, intent(in)  :: xin
        double precision, intent(out) :: dyout
        double precision :: vsum

        vsum  = exp(-1.d0*xin**2)
        dyout = vsum*2.d0/(dsqrt(3.141592653589793238d0))

        end subroutine DERF_MINE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine CDF_Numerical(xin,xc,st,cdf)
        implicit none
        double precision, intent(in)  :: xin,xc,st
        double precision, intent(out) :: cdf
        double precision :: vtmp
        integer :: i,num
        double precision :: xs,xe,dx,x1,x2,pdf1,pdf2
        if(xin.lt.xc)then
          cdf = 0.d0
          xs  = xin -2.d0
          xe  = xin
        else if(xin.ge.xc)then
          cdf = 0.5d0
          xs  = xc
          xe  = xin
        end if
        num = 1001
        dx  = (xe-xs)/dble(num-1)
        do i=1,num-1
           x1 = xs + dble(i-1)*dx
           x2 = xs + dble(i  )*dx
           call normal_dist(x1,xc,st,pdf1)
           call normal_dist(x2,xc,st,pdf2)
           cdf = cdf + (pdf1+pdf2)*dx*0.5d0
        end do
        end subroutine CDF_Numerical
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine normal_dist(xin,xc,st,y)
        implicit none
        double precision, intent(in)  :: xin,xc,st
        double precision, intent(out) :: y
        double precision :: pi

        pi = 4.d0*datan(1.d0)
        y = exp(-1.d0*(xin-xc)**2/2.d0/st/st)/dsqrt(2.d0*pi*st**2)

        end subroutine normal_dist

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
