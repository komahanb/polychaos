subroutine exactoutputfile(jex,kex,par,ipar,fct,coll,dim,dimpc,&
     & nterms,npts,fpcb,rhsF,mreg,fex,fpc)
  use dimpce
  implicit none
  INCLUDE "collsub.h"

  integer :: Jex,Kex,k,j,fct,npts,dim,i,kk,mreg(MAXDAT,DIM),nterms,&
       & dimpc,ipar(MAXVAR),ict
  real*8  :: par(MAXVAR,MAXPAR) 
  real*8  :: fex(Jex,Kex),dxex,dyex!,xex(Jex,Kex),yex(Jex,Kex)
  real*8  :: fpcb(MAXPTS),coll(MAXPTS,MAXVAR),coor(DIM), &
       & PL(DIM,0:MAXTRM),DPL(DIM,0:MAXTRM), &
       & fpc(Jex,Kex),xcoftmp,rhsF(MAXTRM),ddpl(dim,0:maxtrm)
  real*8 :: x(DIM),xin(dim),Xex(DIM,JEx*Kex),xcof(MAXTRM)

  !=========================================

  ! Exact function with output to file

  !=========================================

  open(unit=43,file='output/outputex')

  if (dim.eq.2) then
     write(43,*) 'variables= "x","y","f"'
     write(43,*) 'zone i=',Jex,' j=',Kex,' datapacking=point'
  else if (Dim.eq.1) then
     write(43,*) 'variables= "x","f"'
     write(43,*) 'zone i=',Jex,' datapacking=point'
  end if

  ict = 0
  do i=1,Jex

     do j=1,Kex
        ict = ict + 1
        xin(1) = dble(i-1)/dble(Jex-1)

        DO k=2,DIM
           xin(k) = dble(j-1)/dble(Kex-1)
        END DO

        Xex(1:DIM,ict)=par(1:DIM,1)+(par(1:DIM,2)-par(1:DIM,1))*xin(1:DIM)

        fex(i,j)=f(dim,xex(1:DIM,ict),fct)

        if (dim.le.2)   write(43,*) (xex(k,ict),k=1,DIM),fex(i,j)

     end do
  end do

  close(43)


  !====================================================

  ! Polynomial chaos sample points to file

  !=====================================================

  open(unit=43,file='output/outputsamplePC')
  if (dim.eq.2) then

     write(43,*) 'variables= "x","y","fpcb"'
     write(43,*) 'zone i=',npts,' j=',npts,' datapacking=point'

  else if (Dim.eq.1) then

     write(43,*) 'variables= "x","fpcb"'
     write(43,*) 'zone i=',npts,' datapacking=point'

  end if

  !           write(43,*) 'variables= "x","y","fpcb"'
  !           write(43,*) 'zone i=',npts,' datapacking=point'

  do j=1,npts
     write (43,*) (coll(j,i),i=1,DIM),fpcb(j)
  enddo

  close(43)


  !=====================================================

  ! Evaluate PC in the domain

  !=====================================================

  write(filenum,*)'   >> Evaluating PC in the domain'

  open(unit=43,file='output/outputPC')

!!$  write(43,*) 'variables= "x","y","fpc"'
!!$  write(43,*) 'zone i=',Jex,' j=',Kex,' datapacking=point'
!!$  


  if(dim.eq.2)then
     write(43,*) 'variables= "x","y","fpc"'
     write(43,*) 'zone i=',Jex,' j=',Kex,' datapacking=point'
  else if (Dim.eq.1)then
     write(43,*) 'variables= "x","fpc"'
     write(43,*) 'zone i=',Jex,' datapacking=point'
  end if


   ict = 0
  do j=1,Jex            
     do k=1,Kex
        ict = ict + 1

        fpc(j,k) = 0.0
        do kk=1,nterms
           xcoftmp=xcof(kk)
           do i=1,DIM 
              COOR(1:DIM)=xex(1:DIM,ict)
!              call OTHPL(ipar(i),DIMPC,COOR(i),PL(i,:),DPL(i,:))
!              call ortho(ipar(k),DIMPC,coor(i),PL(i,:),DPL(i,:),ddpl(i,:)) 
              call LEGENDRE(dimpc,coor(i),PL(i,:),DPL(i,:),ddpl(i,:))
              xcoftmp=xcoftmp*PL(i,mreg(kk,i))
           enddo
           fpc(j,k) = fpc(j,k) + xcoftmp
        end do


        write(43,*)(XEX(i,ict),i=1,DIM),fpc(j,k)
     enddo
  enddo
  write(filenum,*)'   >> Completed evaluating PC  in the domain'


contains

  ! Function

  function f(dim,x,fct)

    implicit none

    integer :: fct,dim,k
    real*8  :: x(dim),f


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
       do k=1,DIM        
          f=f+sin(3.0*x(k)-1.5)+cos(3.0*x(k)-3.0)
       end do
       f=f

    end if
  end function f

end subroutine exactoutputfile
