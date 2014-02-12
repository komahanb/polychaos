subroutine dynsampdist(stat,nDIM,DIMPC,ipar,par,makesamples,ntermsold,nterms,nptsold,npts,nptstoaddpercyc,fct,fpcb,gpcb,hpcb,xcof,RN)
  use dimpce
  implicit none
  include 'collsub.h'
  include 'mpif.h'
  real*8::fv,gv(ndim),hv(ndim,ndim)
  integer ::stat, fct,ierr
  integer ::ndim,dimpc,makesamples
  integer,intent(in)::nterms
  integer,intent(in)::npts,nptsold,nptstoaddpercyc
  integer :: ipar(MAXVAR)
  integer :: seed
  real*8 :: dinvnorm,par(nDIM,MAXPAR),x(NDIM)
  real*8,intent(inout)::RN(ndim,maxout)
  real*8::RNTEMP(NDIM,nptstoaddpercyc)
  character*60::filename

  integer :: i,ii,j,jj,jjj,k,kk,kp,l,NTOEX,NTOEXtmp,triangle_node(ndim+1,100000),triangle_num,triangle_coor_num,NCP,node,knnptr(20000),nseed
  integer :: mode
  double precision :: Dtoextmp(ndim),Dtoex(ndim,100000),dist(100000),minftoex(100000),maxftoex(100000),distmean,ftoextry(2,100000)
  double precision :: diff2,RMSE(100000),RMSEmean,EI,Ddibtmp(ndim,0:1000),Dgdibtmp(ndim,0:1000),fdibtmp(0:1000),gdibtmp(ndim,0:1000),hdibtmp(ndim,ndim,0:1000),diffloctmp,difflocmin,difflocavg, SIGMAmean,distcomp
  double precision, dimension(nptstoaddpercyc) :: f
  double precision, dimension(ndim,nptstoaddpercyc) :: df,Dad,v
  double precision, dimension(ndim,ndim,nptstoaddpercyc) :: d2f
  double precision, DIMENSION(200) :: SIGV
  double precision, DIMENSION(200) :: SIGG
  integer :: Taylororder,  NCPG,idec,is,ie,id,point,kpc,nptstoaddpercyctmp,  nptstoaddpercycorig
  double precision :: BETA, GAMM, SIGMA(100000),distmeanloc,RMSEmeanloc,minftoexloc(10000),maxftoexloc(10000),phi,invphi
  character*60::export
  integer::ntermsold
  real*8,intent(inout):: fpcb(MAXPTS),gpcb(NDIM,MAXPTS),hpcb(NDIM,NDIM,maxpts)
  real*8,intent(in):: xcof(MAXTRM)

  
  real*8::diffloc

  real*8::derivdummy(ndim),dblederivdummy(ndim,ndim)

  integer::npass

  if (lhsdyn) then

  !==========================================
  ! LHS - DYNAMIC TRAINING POINT SELECTION 
  !==========================================

     if (id_proc.eq.0) then

        write(filenum,*) '>> [Picking points dynamically for LHS]' 

      
        call get_seed(nseed)
        call latin_random(ndim,nptstoaddpercyc,nseed,Dtoex) 

        do j=1,ndim
           do i =1,nptstoaddpercyc
              dtoex(j,i)=par(j,1)+(par(j,2)-par(j,1))*DTOEX(j,i) !Scale the test candidates to the domain
           end do
        end do


        do ii=1,nptstoaddpercyc
     
           RN(1:ndim,nptsold+ii)=Dtoex(:,ii)
           x(1:ndim)=Dtoex(1:ndim,ii)

           call evalcostf(stat,ndim,fct,x,fv,gv,hv)
           fpcb(nptsold+ii)=fv                   !Store the cost function value at appropriate location
           if (stat.gt.0) then
              gpcb(1:ndim,nptsold+ii)=gv(1:ndim)   !Store them appropriately
           end if

           if (stat.gt.1) then
              hpcb(1:ndim,1:ndim,nptsold+ii)=hv(1:ndim,1:ndim)
           end if

        end do !! ii loop

     end if

     call MPI_Barrier(MPI_COMM_WORLD,ierr)           
     call MPI_BCAST(RN(:,:),ndim*maxout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  
     call MPI_BCAST(fpcb,maxpts,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
     call MPI_BCAST(gpcb,ndim*maxpts,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
     call MPI_BCAST(hpcb,ndim*ndim*maxpts,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 

     goto 888
     return
  end if


  !==========================================
  ! DYNAMIC TRAINING POINT SELECTION WITH MIR
  !==========================================
   
  !  if (abs(nptsold-npts).lt.nptstoaddpercyc) npts=npts+nptstoaddpercyc
  !  if (abs(nptsold-npts).lt.ndim) npts=npts+nptstoaddpercyc

  call mirtunableparams(fct,stat,ndim,nptsold,ncp,taylororder)
  !  print *,nptsold, ncp

  NTOEX=10201

!  NTOEX=5000*NDIM

  if (NTOEX.gt.25000) NTOEX=25000

  !  if (stat.eq.2)  NTOEX=int(NTOEX*nDIM*nDIM)

  if(id_proc.eq.0)  write(filenum,*) '>> Number of test candidates',NTOEX

  if (id_proc.eq.0) then
     call get_seed(nseed)
     call latin_random(ndim,NTOEX,nseed,Dtoex) 

     do j=1,ndim
        do i =1,NTOEX
           dtoex(j,i)=par(j,1)+(par(j,2)-par(j,1))*DTOEX(j,i) !Scale the test candidates to the domain
        end do
     end do

     !        write(export, '(a,i3.3,a)')'testcand.dat'
     !        open(10,file='./KrigSamples/'//export,form='formatted',status='unknown')
     !        read(10,*)(Dtoex(:,i),i=1,NToex)
     !        close(10)
  end if !master

  ! Information sharing by master with slaves        

  call MPI_BCAST(Dtoex(:,:),100000,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
  !  do kk=1,ndim
  call MPI_BCAST(RN(:,:),ndim*maxout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
  !end do
  !if (id_proc.eq.1)  print*, RN(1:ndim,1:nptsold),id_proc,ierr,nptsold,maxout

  call MPI_BCAST(fpcb,maxpts,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(gpcb,ndim*maxpts,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)           


  !  print *, DTOEX(1:nDIm,1:ntoex),id_proc,ntoex

  idec = dble(NTOEX)/dble(num_proc)
  is   = idec*id_proc + 1
  ie   = idec*(id_proc+1)
  if(id_proc.eq.num_proc-1)ie =NTOEX 
  write(filenum,'(11x,a,2i10,a,i3)') '>> [',is,ie,' ] for Processor',id_proc            

  do k=is,ie !Main loop for test candidates

     ! Still need to make sure points are not collinear in higher dimensions!

     call knn(Dtoex(:,k),RN,knnptr,ndim,nptsold,NCP)!? !Find NCP closest points

     NCPG=0
     do j=0,NCP-1  
        node=knnptr(j+1)
        Ddibtmp(1:ndim,j)=RN(1:ndim,node)
        fdibtmp(j)=fpcb(node)
        !        print *, fdibtmp(j)
        if (stat.gt.0) then
           Dgdibtmp(1:ndim,NCPG)=RN(1:ndim,node) 
           gdibtmp(1:ndim,NCPG)=gpcb(1:ndim,node)
           NCPG=NCPG+1
           !           print *, gdibtmp(1:ndim,NCPG-1)
        end if
     end do

     ! Calculate the best parameters beta and gamma
     sigv=0.d0
     sigg=0.d0
    ! CALL MIR_BETA_GAMMA(1, ndim, NCP, Ddibtmp(:,0:NCP-1), fdibtmp(0:NCP-1), SIGV, NCPG , Dgdibtmp(:,0:NCPG-1), gdibtmp(:,0:NCPG-1), SIGG, Taylororder, 1, dble(1.0), BETA, GAMM, IERR)
    ! if (IERR.ne.0)  stop'Error in MIR Local Surrogate Beta Gamma'

     beta=0.5d0

     if (fct.eq.4) then

        gamm=1.0d0

     else if (fct.eq.2) then

        gamm=5.0d0

     else if (fct.eq.6) then

        gamm=0.5d0

     else
           
        gamm=5
        
     end if

     CALL MIR_EVALUATE(1, ndim, 1, Dtoex(:,k), NCP, Ddibtmp(:,0:NCP-1), fdibtmp(0:NCP-1), SIGV, NCPG , Dgdibtmp(:,0:NCPG-1), gdibtmp(:,0:NCPG-1), SIGG, BETA, GAMM, Taylororder, 1, ftoextry(1,k), SIGMA(k), IERR)

     if (IERR.ne.0)stop'Error in MIR Local Surrogate Beta Gamma'


     call evalPC(ndim,dimPC-1,ipar,xcof,Dtoex(:,k),ftoextry(2,k),derivdummy,dblederivdummy) ! Evaluates PC at x and returns yhat 

     !      print *,k,ftoextry(2,k),ftoextry(1,k)!(ftoextry(2,k)-ftoextry(1,k))!,ftoextry(2,k),ftoextry(1,k),

  end do !NTOEX


  ! Information Sharing --Exchange the function values
  do id=0,num_proc-1
     is   = idec*id + 1
     ie   = idec*(id+1)
     if(id.eq.num_proc-1)ie = NTOEX
     call MPI_BCAST(ftoextry(1,is:ie),ie-is+1,MPI_DOUBLE_PRECISION,id,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(ftoextry(2,is:ie),ie-is+1,MPI_DOUBLE_PRECISION,id,MPI_COMM_WORLD,ierr)
  end do

!!$
!!$  if (id_proc.eq.0) then
!!$  do i=1,npts
!!$     print *, RN(1:NDIM,i)
!!$  end do
!!$
!!$  print*,'after'
!!$end if

  if (id_proc.eq.0) then

     do k =1,NTOEX
        minftoex(k)=ftoextry(1,k)
        maxftoex(k)=ftoextry(1,k)
        if (ftoextry(2,k).gt.maxftoex(k)) then
           maxftoex(k)=ftoextry(2,k)
        else 
           minftoex(k)=ftoextry(2,k)
        end if
     end do

     ! Figure out distance to closest real sample point and the mean of all these distances
     distmean=0.0
     do k=1,NTOEX
        dist(k)=1000000000000.0
        do kk=1,nptsold!nhs
           diff2=0.0
           do jj=1,ndim
              diff2=diff2+(Dtoex(jj,k)-RN(jj,kk))**2
           end do
           if ( diff2.lt.dist(k) ) then !current sample closer than previous samples update the smallest distance
              dist(k)=diff2
           end if
        end do
        dist(k)=SQRT(dist(k))
        distmean= distmean+dist(k)
     end do

     distmean=distmean/real(NTOEX)

     ! Pick test candidate with largest difference in values, but above distcomp distance to nearest neighbours
     
     !    if (NTOEX.lt.10) then
     distcomp=distmean 
     !    else
     !       distcomp=1.1*distmean !0.618
     !    end if

     diffloctmp=0.0d0
     diffloc2=0.0d0
     difflocmax=0.0d0
     do k=1,NTOEX
        diffloctmp=(maxftoex(k)-minftoex(k))**2
        diffloc2=diffloc2+diffloctmp
        difflocmax=max(difflocmax,dsqrt(diffloctmp))
     end do
     diffloc2=diffloc2/dble(ntoex)
     diffloc2=sqrt(diffloc2)

     diffloc=0.0
     do ii=1,nptstoaddpercyc
        
        npass=0
        do while (npass.ne.1)

        kp=0
        diffloctmp=0.0d0

        do k=1,NTOEX
           !! Successful Training point passes the follwing tests:
           !! 1. The local difference between the local and global surrogate should > 
           !! 2. The next training point should be atleast distcomp away from the closest existing sample
           !! 3. RMSE should be greater than RMSE mean of Kriging
           !! 4. SIGMA should be greater than SIGMA mean of MIR 

           if ((maxftoex(k)-minftoex(k)).gt.diffloctmp .and. dist(k).ge.distcomp ) then !.and. SIGMA(k).ge.SIGMAmean .and. RMSE(k).ge.RMSEmean
              diffloctmp=maxftoex(k)-minftoex(k)
              kp=k
           end if
        end do
        diffloc=max(diffloc,diffloctmp)
        
        
        if (kp.eq.0) then ! if no successful candidate



           write (filenum,*) '  >>No passing candidate found . . .'
           write (filenum,*) '  >>Relaxing geometric constraint by 2 percent . . .'
           distcomp=0.98*distcomp

        else 

           npass=npass+1 ! one successful candidate

           write(filenum,*)
           write(filenum,'(a,E10.2,a,i6,a,i4)') '>>Loc diff is',diffloctmp,' for test candidate',kp,' at iteration',dyncyccnt        
           write(filenum,'(a,20F10.2)')' x = ',DTOEX(1:NDIM,kp)


        end if

     end do! while loop execute until a passing candidate is found


!!$
!!$           write (filenum,*) 'Could not find suitable test candidate just take the one with largest difference'
!!$           diffloctmp=0.0
!!$           do k=1,NTOEX
!!$              if ((maxftoex(k)-minftoex(k)).gt.diffloctmp) then
!!$                 diffloctmp=maxftoex(k)-minftoex(k)
!!$                 kp=k
!!$              end if
!!$           end do
!!$

     ! Trick to not consider this point again
     maxftoex(kp)=minftoex(kp)

     ! Update other minimum distances
     do k=1,NTOEX
        if (k.ne.kp) then
           diff2=0.0
           do jj=1,ndim
              diff2=diff2+(Dtoex(jj,k)-Dtoex(jj,kp))**2
           end do
           diff2=SQRT(diff2)
           if ( diff2.lt.dist(k) ) then
              dist(k)=diff2
           end if
        end if
     end do

        RN(1:ndim,nptsold+ii)=Dtoex(:,kp)

        x(1:ndim)=Dtoex(1:ndim,kp)

        call evalcostf(stat,ndim,fct,x,fv,gv,hv)
        fpcb(nptsold+ii)=fv                   !Store the cost function value at appropriate location
        if (stat.gt.0) then
           gpcb(1:ndim,nptsold+ii)=gv(1:ndim)   !Store them appropriately
        end if

        if (stat.gt.1) then
           hpcb(1:ndim,1:ndim,nptsold+ii)=hv(1:ndim,1:ndim)
        end if

        !        if (fct.eq.20) then
        !           x(1)=x(1)*180.0/4.0/atan(1.0)   !RADINS TO degree
        !        end if

!!$
!!$        ! Evaluate desired quantities
!!$        if(nstyle.eq.0) then
!!$           Dad(:,ii)=Dtoex(:,kp)
!!$           call evalfunc(Dtoex(:,kp),ndim,fct,0,hstatad(ii),f(ii),df(:,ii),d2f(:,:,ii),v(:,ii))
!!$        end if

        !INSTEAD STORE INTO RN

     end do !! ii loop

!     difflocmax=diffloc !stores maxdiff to write as output

  END IF! master thread 

  call MPI_BCAST(RN(:,:),ndim*maxout,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  
  call MPI_BCAST(fpcb,maxpts,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
  call MPI_BCAST(gpcb,ndim*maxpts,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
  call MPI_BCAST(hpcb,ndim*ndim*maxpts,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 

  goto 888

  !1 ) Assume that the model is built with a particular order say 3, [x,f(x)] are already in a file named sampleandvalues.dat (for now lets forget about center and borders)

  !2 ) Now for order 4  we want to choose nptstoaddpercyc points to samplepoints.dat in this subroutine and pass it as RNF

  !3 ) Make setting up RHS part of the program, read the previous [x,f(x)] from the file samplepoints.dat and then evaluate the remaining points from evalfunc.f90 and immediately append the values to samplepoints.dat


!!$  if (id_proc.eq.0) then
!!$  ! dont forget to scale and output the RN
!!$
!!$  do i=1,npts
!!$     print *,i,RN(1:NDIM,i)
!!$  end do
!!$  end if
!!$  call MPI_Barrier(MPI_COMM_WORLD,ierr)!slaves wait until master arrives
!!$  call stop_all

888 continue

  if (id_proc.eq.0 ) then

     if (fct.eq.20) then

        if (fctindx.eq.0) then
           open(unit=88,file='output/tecsamp00.dat')
        else 
           open(unit=88,file='output/tecsamp04.dat')
        end if

        if (ndim.eq.2) then
           write(88,*) 'variables= "x","y","f"'
           write(88,*) 'zone i=',npts,' datapacking=point'

        else if (nDim.eq.1) then
           write(88,*) 'variables= "x","f"'
           write(88,*) 'zone i=',npts,' datapacking=point'
        end if

        do j=1,npts

           x(1:nDIM)=RN(1:nDIM,j)

           x(1)=x(1)*180.0/4.0/atan(1.0)   !RADINS TO degree

           write (88,*) (x(ii),ii=1,ndim),fpcb(j)

        end do ! loop over all points
        close(88)
     end if

     if(ndim.le.2) then

        !====================================================

        ! Polynomial chaos sample points to file

        !=====================================================
        if (ndim.eq.2) then

           open(unit=43,file='output/outputsamplePCtmp')

           write(43,*) 'variables= "x","y","fpcb"'
           write(43,*) 'zone i=',npts,' datapacking=point'

           !  write(43,*) 'variables= "x","y"'
           !  write(43,*) 'zone i=',npts,' j=',npts,' datapacking=point'

        else if (nDim.eq.1) then

           write(43,*) 'variables= "x"'
           write(43,*) 'zone i=',npts,' datapacking=point'

        end if

        do j=1,npts
           write (43,*) (RN(i,j),i=1,nDIM),0.0d0
        enddo

        close(43)

     end if
  end if

  call MPI_Barrier(MPI_COMM_WORLD,ierr)!slaves wait until master arrives

  return

end subroutine dynsampdist



subroutine mirtunableparams(fct,stat,ndim,nhs,ncp,taylororder)
  implicit none
  integer,INTENT(IN)::fct,ndim,nhs,stat
  INTEGER,INTENT(OUT)::NCP,TAYLORORDER

  if (nhs.le.50)  then  

     NCP=nhs

  else

     NCP=50
     !    tAYLORORDER=5
  end if

  ! Higher dimensional test functions

  if(ndim.gt.2) then

     !NCP= ceiling(0.25*dble(nhs))
     ! use 25% of the existing data points?

     if (Nhs.lt.50) then
        NCP=nhs
     else
        ncp=50
     end if

  end if

  ! Recommended Taylor order of expansion by Qiqi Wang

  if (stat.eq.0) then

     Taylororder=NCP

  else

     Taylororder=NCP+ndim*NCP

  end if

  return
end subroutine mirtunableparams

!+++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine knn(SC,sample,knnptr,ndim,nhs,NCP)
  ! Subroutine to find NCP closest neighbours from array sample to point SC
  implicit none
  integer :: ndim,nhs,NCP,j,k,knnptr(NCP)
  real*8 ::SC(ndim),sample(ndim,nhs),norm2(nhs),sn

     ! Calculate all distances
     do j=1,nhs
        norm2(j)=0.0
        do k=1,ndim
           norm2(j)=norm2(j)+(sample(k,j)-SC(k))**2
        end do
     end do

     ! Find NCP closest neighbours

     do k=1,NCP
        sn=10000000000.0
        do j=1,nhs
           if (norm2(j).lt.sn) then
              sn=norm2(j)
              knnptr(k)=j
           end if
        end do
        norm2(knnptr(k))=10000000000.0
     end do

   end subroutine knn
