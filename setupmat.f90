subroutine setupmat(stat,dim,dimpc,npts,ipar,RN,numpc,xmat)
  implicit none
  include "collsub.h"

  integer :: stat,dim,dimpc,nterms,i,j,mreg(MAXDAT,DIM),npts,ipar(MAXVAR),kk,k,jj
  integer,intent(out)::numpc
  real*8  :: xmat(MAXDAT,MAXDAT), coll(MAXPTS,DIM), PL(DIM,0:MAXTRM),DPL(DIM,0:MAXTRM),ddpl(dim,0:maxtrm)!,xmatG(MAXDAT,MAXDAT)
  real*8::RN(dim,maxpts)

  !-------------------------  
  ! (2)  Assemble Matrix
  !--------------------------

  do j=1,npts
     coll(j,1:dim)=RN(1:dim,j)
  end do

  call multidx(MAXDAT,DIM,DIMPC,mreg,nterms) ! get multiindex notation for tensor product

!!$
!!$print *,nterms
!!$ do i=1,nterms
!!$   write(filenum,'(i5,a,99i4)')i,':',(mreg(i,j),j=1,DIM)
!!$ end do
!!$  stop


  PL(:,:)=0.0d0
  DPL(:,:)=0.0d0
  ddpl(:,:)=0.0d0

  if (stat.eq.0) then

     xmat(:,:)=1.0d0

     do i=1,npts !number of equations
        do k=1,DIM
           call ortho(ipar(k),DIMPC,coll(i,k),PL(k,:),DPL(k,:),ddpl(k,:)) 
!           call OTHPL(ipar(k),DIMPC,coll(i,k),PL(k,:),DPL(k,:)) 

        end do

        do j=1,nterms !number of terms in each equation
           do k=1,DIM
!              print *,i, mreg(j,k),xmat(i,j), PL(k,mreg(j,k))
              xmat(i,j) = xmat(i,j)*PL(k,mreg(j,k))
              ! !print*, xmat(i,j)
           end do
        end do
     end do

     numpc=npts


  else if (stat.eq.1) then

     xmat(:,:)=1.0d0

     numpc=0
     do i=1,npts
        do k=1,DIM
!           call OTHPL(ipar(k),DIMPC,coll(i,k),PL(k,:),DPL(k,:)) 
         call ortho(ipar(k),DIMPC,coll(i,k),PL(k,:),DPL(k,:),ddpl(k,:)) 
!           call LEGENDRE(dimpc,X,PL(k,:),DPL(k,:),ddpl(k,:))
        end do

        ! Function value

        numpc=numpc+1
        do j=1,nterms
           do k=1,DIM
              xmat(numpc,j) = xmat(numpc,j)*PL(k,mreg(j,k))
!              print *,i, mreg(j,k),xmat(i,j), PL(k,mreg(j,k)),DPL(k,mreg(j,k))
           end do
        end do

        ! Gradient values
        do kk=1,DIM
           numpc=numpc+1
           do j=1,nterms
              do k=1,DIM
                 if (k.eq.kk) then
                    xmat(numpc,j) = xmat(numpc,j)*DPL(k,mreg(j,k))
                 else
                    xmat(numpc,j) = xmat(numpc,j)*PL(k,mreg(j,k))
                 end if
              end do
           end   do
        end do
     end do

  else if (stat.eq.2) then

     xmat(:,:)=1.0d0

     numpc=0
     do i=1,npts
        do k=1,DIM

!           coll(i,1:dim)=0.5d0

           call ortho(ipar(k),DIMPC,coll(i,k),PL(k,:),DPL(k,:),ddpl(k,:)) 
        end do

        ! Function value

        numpc=numpc+1
        do j=1,nterms
           do k=1,DIM
              xmat(numpc,j) = xmat(numpc,j)*PL(k,mreg(j,k))
           end do
        end do

        ! Gradient values
        do kk=1,DIM
           numpc=numpc+1
           do j=1,nterms
              do k=1,DIM
                 if (k.eq.kk) then
                    xmat(numpc,j) = xmat(numpc,j)*DPL(k,mreg(j,k))
                 else
                    xmat(numpc,j) = xmat(numpc,j)*PL(k,mreg(j,k))
                 end if
              end do
           end   do
        end do

        ! Hessian values
        
        do kk=1,DIM

           do jj=1,DIM
             
              if (kk.eq.jj.or.kk.gt.jj)  numpc=numpc+1

              do j=1,nterms

                 do k=1,DIM !x1,x2,x3.....

                    if (kk.eq.jj) then

                       if (k.eq.kk) then
                          xmat(numpc,j) = xmat(numpc,j)*DDPL(k,mreg(j,k))
                       else
                          xmat(numpc,j) = xmat(numpc,j)*PL(k,mreg(j,k))
                       end if

                    else if (kk.gt.jj)then !upper diagonal elements

                       if (k.eq.kk.or.k.eq.jj) then
                          xmat(numpc,j) = xmat(numpc,j)*DPL(k,mreg(j,k))
                       else
                          xmat(numpc,j) = xmat(numpc,j)*PL(k,mreg(j,k))         
                       end if
!!$
!!$                    else !lower diagonal elements
!!$
!!$                     !  numpc=numpc-1
!!$!                       stop'Something is wrong'

                    end if

                 end do
              end   do
           end do
        end do

        
!!$        do ii=1,numpc
!!$           write(filenum,*) (xmat(ii,j),j=1,nterms)
!!$           write(filenum,*)
!!$        end do
!!$
!!$stop

     end do !npts loop



     
  end if

end subroutine setupmat
 
