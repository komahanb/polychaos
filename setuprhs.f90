subroutine setupRHS(stat,dim,npts,fpcb,gpcb,hpcb,numpc,rhsF)

  implicit none
  include 'collsub.h'

  integer,intent(in):: stat,dim,npts
  integer::j,k,i
  integer,intent(out)::numpc
  
  real*8,intent(in) :: fpcb(MAXPTS),gpcb(DIM,MAXPTS),hpcb(dim,dim,maxpts)

  real*8::rhsF(maxdat)
  real*8 :: tmp(maxdat)

  if (stat.eq.0) then
     numpc=0
     do j=1,npts
        numpc=numpc+1
        tmp(j) = fpcb(j)
     enddo
     rhsF(:)=tmp(:)

  else if (stat.eq.1) then

     numpc=0
     do j=1,npts
        numpc=numpc+1
        tmp(numpc) = fpcb(j)

        do k=1,DIM

           numpc=numpc+1
           tmp(numpc) = gpcb(k,j)

        end do
     enddo
     rhsF=tmp

  else if (stat.eq.2) then

     
     numpc=0
     do j=1,npts
        numpc=numpc+1
        tmp(numpc) = fpcb(j)

        do k=1,DIM

           numpc=numpc+1
           tmp(numpc) = gpcb(k,j)

        end do


        do k=1,dim

           do i=1,dim

              if (k.eq.i.or.k.gt.i) numpc=numpc+1

              if (k.eq.i) then !diagonal

                 tmp(numpc) = hpcb(k,i,j)

              else if (k.gt.i) then !upper diagonal

                 tmp(numpc) = hpcb(k,i,j)

              else

                 ! numpc=numpc-1

              end if

           end do
        end do

     enddo

     rhsF=tmp

  else

     stop 'unknown stat'

  end if

end subroutine setupRHS
