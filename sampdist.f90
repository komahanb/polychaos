subroutine sampdist(stat,DIM,DIMPC,ipar,par,makesamples,nterms,npts,fct,RN)
  use dimpce
  implicit none
  include 'collsub.h'

  integer ::stat,fct, DIM, dimpc,makesamples,j,i

  integer :: ipar(MAXVAR),ierr
  integer :: seed

  integer,intent(in)::nterms
  integer,intent(in)::npts

  real*8,intent(out) :: RN(DIM,MAXPTS)

  real*8::par(DIM,MAXPAR)

  character*60::filename
!  double precision,dimension(2,dim) :: bound

  integer::nptstmp

!  integer::ctr
!  common/global/ctr


  !========================================================
  ! Figure out how many points are required and make them
  !=========================================================

  write(filenum,*)
  write(filenum,*) '================================================='
  write(filenum,*) '         Summary of function requirements        '
  write(filenum,*) '================================================='
  write(filenum,*)


  write(filenum,*)'>> Number of Function evals needed :',npts

  if (Stat.eq.1) then

     write(filenum,*)'>> Number of Gradient evals needed :',npts
     write(filenum,*)'>> Equivalent Function evals needed :', npts+npts*2

  else if (stat.eq.2) then

     write(filenum,*)'>> Number of Gradient evals needed :',npts
     write(filenum,*)'>> Number of Hessian  evals needed :',npts
     write(filenum,*)'>> Equivalent Function evals needed :',npts+npts*2+npts*dim

  end if

  write(filenum,*)

  call getfilename(dim,fct,dimpc,stat,55,filename)

  if (makesamples.eq.1) then !
 
!        nptstmp=1 ! (center of the domain)

        ! Initialize Random number generators with system time 

!        call get_seed(seed)
        ! Generate random numbers and store in RN
 !       call latin_random(dim,npts-1,seed,RN(:,2:npts))

        if (randomflag.eq.1) then

           write(filenum,*)'>> Initial Sample Points by Latin Hypercube'
           write(filenum,*)'>> Including center of the domain'
           
           ! force center of the domain
           do j=1,dim
              RN(j,1)=0.5d0
           end do

           call get_seed(seed)
           call latin_random(dim,npts-1,seed,RN(:,2:npts))       

        else if (randomflag.eq.2) then

           write(filenum,*)'>> Initial Sample Points by NIEDER Sequence'
           call get_seed(seed)
           call nieder(seed,dim,npts,RN(:,1:npts))
           !                 print*,sample(:,1:nhs-nhstmp)
           !                 stop    
        else if (randomflag.eq.3) then

           write(filenum,*)'>> Initial Sample Points by Halton Sequence'
           call halton_real(dim,npts,RN(:,1:npts))
           !                print*,sample(:,1:nhs-nhstmp)
           !                stop
        else if (randomflag.eq.4) then

           write(filenum,*)'>> Initial Sample Points by Hammersley Sequence'
           call hammersley_real(dim,npts,RN(:,1:npts))
           !                 print*,sample(:,1:nhs-nhstmp)
           !                stop

        else if (randomflag.eq.5) then

           write(filenum,*)'>> Initial Sample Points by Sobol Sequence'
           call get_seed(seed)
           call sobol_real(seed,dim,npts,RN(:,1:npts))
           !                 print*,sample(:,1:nhs-nhstmp)
           !                stop

        else if (randomflag.eq.6) then

           write(filenum,*)'>> Initial Sample Points by Faure Sequence'
           call get_seed(seed)
           call faure_real(seed,dim,npts,RN(:,1:npts))
           !                print*,sample(:,1:nhs-nhstmp)
           !                stop

        else !if (randomflag.eq.6) then

           write(filenum,*)'>> Initial Sample Points by Which Sequence'
           print*,"Cool.. Please go ahead and implement"
           stop

        end if


        ! Write those points to a file (sample points can be used again)
        open(30,file=filename,form='formatted',status='unknown')
        write(30,*)(RN(:,i),i=1,npts)
        close(30)

        write(filenum,'(2a)')' >> Training points are written to file :',filename

  else ! read from the file

     write(filenum,'(2a)')' >> Training points are read from file : ',filename

     open(30,file=filename,form='formatted',status='unknown')
     read(30,*)(RN(:,i),i=1,npts)
     close(30)

  end if

  !====================================
  ! Scale it to the domain bounds (par)
  !====================================
  
  do j=1,DIM

     do i = 1,npts
        RN(j,i)=par(j,1)+(par(j,2)-par(j,1))*RN(j,i)
     end do

  end do


!? Do I need this anymore? Tecplot may be affected
  
  if(dim.le.2) then

     !====================================================
     ! Polynomial chaos sample points to file
     !=====================================================
     if (dim.eq.2) then

        open(unit=43,file='output/outputsamplePCtmp')

        write(43,*) 'variables= "x","y","fpcb"'
        write(43,*) 'zone i=',npts,' datapacking=point'

     else if (Dim.eq.1) then

        write(43,*) 'variables= "x"'
        write(43,*) 'zone i=',npts,' datapacking=point'

     end if

     do j=1,npts
        write (43,*) (RN(i,j),i=1,DIM),0.0d0
     enddo

     close(43)
     
  end if

end subroutine sampdist

!+++++++++++++++++++++++++++++++++++++++++++++++++++++


  subroutine getnpts(solver,stat,dim,nterms,oversamp,npts)

    implicit none

    integer,intent(in)  ::solver,stat,nterms,oversamp,dim
    integer,intent(out) ::npts
    
    if (oversamp.eq.0) stop 'Wrong ovesampling factor'    

    if (stat.eq.0) then ! F alone

       if (solver.eq.1) then !LU decomposition

          npts=nterms
          
       else !SV decomposition

          npts=oversamp*nterms

       end if

    else if (stat.eq.1) then ! F and G

       if (solver.eq.1) then !LU decomposition

          npts=ceiling(real(nterms)/real(1+DIM))

       else

          npts=ceiling(real(oversamp*nterms)/real(1+DIM))

       end if

    else if (stat.eq.2) then ! F , G and Hessian

       if (solver.eq.1) then !LU

          npts=ceiling(real(nterms)/(real(DIM)*real(1+DIM)/real(2.0d0)))

       else
          npts=ceiling(real(oversamp*nterms)/real(1+DIM+(DIM*DIM+DIM)/2))

       end if

    else

       Stop'Unsupported Stat'

    end if

    return
  end subroutine getnpts
!+===================================================

!!$     
!!$     if (int((2**DIM)+1).le.npts) then
!!$        
!!$        bound(1,:)=0.0d0!0.05
!!$        bound(2,:)=1.0d0!0.95
!!$        
!!$        ctr=2
!!$        call recurcorn(1,dim,RN,npts,bound)
!!$        
!!$        do j=1,dim
!!$           RN(j,(2**Dim)+1)=0.5d0
!!$        end do
!!$        
!!$        nptstmp=(2**DIM)+1
!!$        
!!$        write(filenum,*)'>> Including center, bounds of the domain',nptstmp
!!$        
!!$
!!$
!!$        if (nptstmp.lt.npts) then
!!$           write(filenum,*)'>> Remaining ',npts-nptstmp,'training points via LHS'
!!$           call get_seed(seed) ! Initialize Random number generators with system time 
!!$           call latin_random(DIM,npts-nptstmp,seed,RN(:,nptstmp+1:npts)) ! Get lattice hypercube random numbers
!!$        end if
!!$        
!!$     else
!!$        

!======================================================================
! Not used

function DINVNORM(p)

  implicit none

  real*8 :: dinvnorm,p,p_low,p_high
  real*8 :: a1,a2,a3,a4,a5,a6
  real*8 :: b1,b2,b3,b4,b5
  real*8 :: c1,c2,c3,c4,c5,c6
  real*8 :: d1,d2,d3,d4
  real*8 :: z,q,r
  a1=-39.6968302866538
  a2=220.946098424521
  a3=-275.928510446969
  a4=138.357751867269
  a5=-30.6647980661472
  a6=2.50662827745924
  b1=-54.4760987982241
  b2=161.585836858041
  b3=-155.698979859887
  b4=66.8013118877197
  b5=-13.2806815528857
  c1=-0.00778489400243029
  c2=-0.322396458041136
  c3=-2.40075827716184
  c4=-2.54973253934373
  c5=4.37466414146497
  c6=2.93816398269878
  d1=0.00778469570904146
  d2=0.32246712907004
  d3=2.445134137143
  d4=3.75440866190742
  p_low=0.02425
  p_high=1-p_low
  if(p.lt.p_low) goto 201
  if(p.ge.p_low) goto 301
201 q=dsqrt(-2*dlog(p))
  z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
  goto 204
301 if((p.ge.p_low).and.(p.le.p_high)) goto 202
  if(p.gt.p_high) goto 302
202 q=p-0.5
  r=q*q
  z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
  goto 204
302 if((p.gt.p_high).and.(p.lt.1)) goto 203
203 q=dsqrt(-2*dlog(1-p))
  z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
204 dinvnorm=z

  return

end function dinvnorm



recursive subroutine recurcorn(dimact,ndim,sample,nhs,bounds)

  implicit none

  common/global/ctr

  integer :: dimact,ndim,nhs,ctr,j,ctrin
  real*8 :: sample(ndim,nhs),bounds(2,ndim)

  ctrin=ctr
  do j=2,1,-1
     sample(dimact,ctr)=bounds(j,dimact)
     if (dimact.gt.1) then
        sample(1:dimact-1,ctr)=sample(1:dimact-1,ctrin)
     end if
     if (dimact.lt.ndim) call recurcorn(dimact+1,ndim,sample,nhs,bounds)
     if (j.ne.1) ctr=ctr+1
  end do
end subroutine recurcorn

subroutine conv_normal(n,array)
 
  INTEGER :: n
  INTEGER :: i
  REAL*8 :: array(n), pi, temp, mean = 0.5, sd = 0.15
 
  pi = 4.0d0*ATAN(1.0d0)
  CALL RANDOM_NUMBER(array) ! Uniform distribution
 
! Now convert to normal distribution
  DO i = 1, n-1, 2
    temp = sd * SQRT(-2.0*LOG(array(i))) * COS(2*pi*array(i+1)) + mean
    array(i+1) = sd * SQRT(-2.0*LOG(array(i))) * SIN(2*pi*array(i+1)) + mean
    array(i) = temp
  END DO
 
! Check mean and standard deviation
  mean = SUM(array)/n
  sd = SQRT(SUM((array - mean)**2)/n)
 
!  WRITE(*, "(A,F8.6)") "Mean = ", mean
!  WRITE(*, "(A,F8.6)") "Standard Deviation = ", sd
  
end subroutine conv_normal
