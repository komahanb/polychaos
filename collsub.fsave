c$$$   Subroutines to identify sample points and to perform transformations 
c$$$   This file contains all collocation related routines.
c$$$
c$$$
c$$$    Input Variables: 
c$$$      ipar    : Array of distribution type of input variables
c$$$                (e.g., 1 => uniform, 2 => normal, etc.)
c$$$      par     : 2-D array of distribution parameters for each input
c$$$                (e.g., par(1,2) => 2nd distribution parameter of the
c$$$                first input, which has distribution iparno(1)
c$$$      ndim    : Number of input random variables
c$$$      norder  : Order of the polynomial chaos expansion 
c$$$      IPRN    : Print option (IPRN=1 => print the collocation points)
c$$$      ITRNS   : Transform the collocation points into the actual sample
c$$$                values (ITRNS=1 => transform). If not, just leave them 
c$$$                as sample values of the standard random variables.
c$$$      MAXVAR  : Maximum number of inputs (required for physical dimension 
c$$$                of coll array)
c$$$      MAXPTS  : Maximum number of collocation points (required for
c$$$                physical dimension of coll array)
c$$$      MAXPAR  : Maximum number of parameters a distribution can have
c$$$                (e.g., normal distribution has 2 params, etc.)
c$$$    Output Variables:
c$$$      coll    : array of the collocation points coll(point_number,variable)
c$$$      ictr    : total number of collocation points
      

      subroutine colgen (coll,par,ipar,ndim,norder,IPRN,ITRNS,ictr)

      INCLUDE "collsub.h"
     
      dimension coll(MAXPTS,MAXVAR), iparno(MAXDST), par(MAXVAR,MAXPAR)
      dimension ipar(MAXVAR), param(MAXPAR), ocoll(MAXPTS,MAXVAR)

      data iparno /2,2,2,2,1,1,0,0,0,0,
     .           0,0,0,0,0,0,0,0,0,0/

      do i=1,MAXPTS
        do j=1,MAXVAR
          coll(i,j) = 0.0
        enddo
      enddo

      jctr = 2
      
      if (norder.ne.1) goto 20
C 1st order terms
      call trm1(coll,ndim,1.0,jctr)      

C  2nd order terms
 20   if (norder.ne.2) goto 30
      x = 1.732
      call trm1(coll,ndim,x,jctr)
      call trm22(coll,ndim,x,-x,jctr)
      call trm22(coll,ndim,-x,x,jctr)
      call trm22(coll,ndim,x,x,jctr)
      call trm22(coll,ndim,-x,-x,jctr)
      
 30   if (norder.ne.3) goto 40
C  3rd order terms
      x = 2.334
      y = 0.742
      call trm1(coll,ndim,x,jctr)
      call trm1(coll,ndim,y,jctr)
      call trm22(coll,ndim,y,y,jctr)
      call trm22(coll,ndim,-y,-y,jctr)
      call trm22(coll,ndim,-y,y,jctr)
      call trm22(coll,ndim,y,-y,jctr)
      call trm22(coll,ndim,x,y,jctr)
      call trm22(coll,ndim,y,x,jctr)
      call trm22(coll,ndim,x,-y,jctr)
      call trm22(coll,ndim,-y,x,jctr)
      call trm33(coll,ndim,y,-y,y,jctr)
      call trm33(coll,ndim,-y,y,y,jctr)
      call trm33(coll,ndim,y,y,-y,jctr)
      call trm33(coll,ndim,y,-x,y,jctr)
      call trm33(coll,ndim,-y,x,-y,jctr)
      call trm33(coll,ndim,-x,y,y,jctr)
      call trm33(coll,ndim,x,-y,-y,jctr)
      call trm33(coll,ndim,y,y,-x,jctr)
      call trm33(coll,ndim,-y,-y,x,jctr)
      
 40   if (norder.ne.4) goto 50
C  4th order terms
      y = 1.356
      x = 2.860
      call trm1(coll,ndim,x,jctr)
      call trm1(coll,ndim,y,jctr)
      call trm22(coll,ndim,y,y,jctr)
      call trm22(coll,ndim,-y,-y,jctr)
      call trm22(coll,ndim,y,-y,jctr)
      call trm22(coll,ndim,-y,y,jctr)
      call trm22(coll,ndim,x,y,jctr)
      call trm22(coll,ndim,y,x,jctr)
      call trm22(coll,ndim,x,-y,jctr)
      call trm22(coll,ndim,-y,x,jctr)
      call trm33(coll,ndim,y,-y,y,jctr)
      call trm33(coll,ndim,-y,y,y,jctr)
      call trm33(coll,ndim,y,y,-y,jctr)
      call trm33(coll,ndim,y,y,y,jctr)
      call trm33(coll,ndim,-y,-y,-y,jctr)
      call trm33(coll,ndim,y,y,-y,jctr)
      call trm33(coll,ndim,y,-x,y,jctr)
      call trm33(coll,ndim,-y,x,-y,jctr)
      call trm33(coll,ndim,-x,y,y,jctr)
      call trm33(coll,ndim,x,-y,-y,jctr)
      call trm33(coll,ndim,y,y,-x,jctr)
      call trm33(coll,ndim,-y,-y,x,jctr)      
      call trm44(coll,ndim,y,-y,y,-y,jctr)      
      call trm44(coll,ndim,-y,y,-y,y,jctr)            
      call trm44(coll,ndim,-x,y,-y,x,jctr)      
       
 50   if (norder.ne.5) goto 999
C  5th order terms
      call trm1(coll,ndim,3.324,jctr)
      call trm1(coll,ndim,0.617,jctr)
      call trm1(coll,ndim,1.890,jctr)
!      write (filenum,*) 'ERROR... Not Yet Implemented'
      stop'NOT yet implemented '

 999  ictr = jctr - 1

      do i=1,ndim
        do j=1,ictr
          ocoll(j,i) = coll(j,i)
        enddo
      enddo

      if (ITRNS.EQ.1) then
        do i=1,ndim
          do j=1,iparno(ipar(i))
            param(j) = par(i,j)
          enddo
          do j=1,ictr
            coll(j,i) = xitran(coll(j,i), ipar(i), param)
          enddo
        enddo
      endif
        
      if (IPRN.EQ.1) then
        do j=1,ictr
          write (*,21) j, (coll(j,i),i=1,ndim)
        enddo
      endif

 21   format (i7,2X,50(G11.4,X))

      return
      end


C  Transformation from a Gaussian to another distribution type
C  The following transformations are available in this version

C  1  - Uniform (a,b)              UNIFORM
C  2  - Normal (mu, sigma)         NORMAL
C  3  - Lognormal (mu, sigma)      LOGNORMAL
C  4  - Gamma (a,b)                GAMMA
C  5  - Exponential (lambda)       EXPONENTIAL
C  6  - Weibull (a)                WEIBULL
C  7  - Extreme Value              EXTREME

C  Reading the distribution information from a file
C  Read Variable number and Type of distribution (numerical or name)
C  

      function xitran(xi, itrans, param)
      INCLUDE "collsub.h"

      dimension param(MAXPAR)

      sqrt2 = 1.41421356

      if ( (itrans.lt.1).OR.(itrans.gt.MXDNOW) ) then
!        write (filenum,*) 'Error: Unknown distribution type'
        stop'Unknownd distribution'
      endif

      if (itrans.eq.1) xnew = param(1) + (param(2)-param(1))*
     .     (0.5 + 0.5*erf(xi/sqrt2))
      if (itrans.eq.2) xnew = param(1) + param(2)*xi
      if (itrans.eq.3) xnew = exp(param(1) + param(2)*xi)
      if (itrans.eq.4) xnew = param(1)*param(2)*
     .  (xi*sqrt2/(3.0*sqrt(param(1))) + 1 - 2.0/(9.*param(1)))**3.0
      if (itrans.eq.5) xnew = -1.0/param(1) * 
     .     alog(0.5 + 0.5*erf(xi/sqrt2))
      if (itrans.eq.6) xnew = (-1.0*log(0.5 + 0.5*erf(xi/sqrt2)))
     .     **(1.0/param(1))
      if (itrans.eq.7) xnew = -log(-1.0*log(0.5 + 0.5*erf(xi/sqrt2)))

      xitran=xnew
      return
      end
      
      

      subroutine trm1(coll,ndim,xval,jctr)
      INCLUDE "collsub.h"
      dimension coll(MAXPTS,MAXVAR)
      do j=1,ndim
          coll(jctr,j) = xval
          coll(jctr+ndim,j) = -xval
          jctr = jctr + 1
      enddo  
      jctr = jctr + ndim

      return
      end

      subroutine trm22(coll,ndim,xval,yval,jctr)
      INCLUDE "collsub.h"
      dimension coll(MAXPTS,MAXVAR)

C      ixtr = (ndim*(ndim-1))/2
      do j=1,ndim-1
        do k=j+1,ndim
          coll(jctr,j) = xval
          coll(jctr,k) = yval
C          coll(jctr+ixtr,j) = yval
C          coll(jctr+ixtr,k) = xval
          jctr = jctr + 1
        enddo  
      enddo
C      jctr = jctr + ixtr

      return
      end

      subroutine trm33(coll,ndim,xval,yval,zval,jctr)
      INCLUDE "collsub.h"
      dimension coll(MAXPTS,MAXVAR)

C      ixtr = (ndim*(ndim-1)*(ndim-2))/6
      do j=1,ndim-2
        do k=j+1,ndim-1
          do l=k+1,ndim
            coll(jctr,j) = xval
            coll(jctr,k) = yval
            coll(jctr,l) = zval
            jctr = jctr + 1
          enddo  
        enddo
      enddo

      return
      end


      subroutine trm44(coll,ndim,x,y,z,w,jctr)
      INCLUDE "collsub.h"
      dimension coll(MAXPTS,MAXVAR)

C      ixtr = (ndim*(ndim-1)*(ndim-2)*(ndim-3))/2
      do j=1,ndim-3
        do k=j+1,ndim-2
          do l=k+1,ndim-1
            do m=l+1,ndim
              coll(jctr,j) = x
              coll(jctr,k) = y
              coll(jctr,l) = z
              coll(jctr,m) = w
              jctr = jctr + 1
            enddo  
          enddo
        enddo
      enddo

      return
      end

      subroutine trm55(coll,ndim,x,y,z,w,u,jctr)
      INCLUDE "collsub.h" 
      dimension coll(MAXPTS,MAXVAR)

C      ixtr = (ndim*(ndim-1)*(ndim-2)*(ndim-3))/2
      do j=1,ndim-4
        do k=j+1,ndim-3
          do l=k+1,ndim-2
            do m=l+1,ndim-1
              do n=m+1,ndim
                coll(jctr,j) = x
                coll(jctr,k) = y
                coll(jctr,l) = z
                coll(jctr,m) = w
                coll(jctr,n) = u
                jctr = jctr + 1
              enddo  
            enddo
          enddo
        enddo
      enddo
        

      return
      end


      subroutine multidx(maxdim,DIM,DIMPC,mregout,nterms)

      implicit none
      integer :: i,ii,j,k,isum,DIM,DIMPC,nterms,ent,maxdim
      integer :: mreg(maxdim,DIM),mregout(maxdim,DIM)
      
      call combination(DIM+DIMPC,DIM,nterms)    
                       
      mreg(:,:) = 0 

      if (DIMPC.ne.0) then

         do 100 i=1,nterms
            isum = 0
            do j=1,DIM
               isum = isum + mreg(i,j)
            end do
            if(isum.ne.DIMPC)then
               mreg(i+1,:) = mreg(i,:)
               mreg(i+1,1) = mreg(i+1,1) + 1
               go to 100
            else
               do j=1,DIM
                  if(mreg(i,j).ne.0)then
                     if(j.eq.DIM) go to 200
                     mreg(i+1,:) = mreg(i,:)
                     mreg(i+1,j) = 0
                     mreg(i+1,j+1) = mreg(i+1,j+1) + 1
                     go to 100
                  end if
               end do
               stop 'No target j in Make_Mreg'
            end if
 100     continue

!         write(filenum,'(3i6)')DIMPC,nterms,DIM
         stop 'Error in Make_Mreg'
 200     continue
        
         if(i.ne.nterms)then
 !           write(filenum,*) i,nterms,DIMPC
            stop 'i.ne.nterms in Make_Mreg'
         end if

      end if


!     Resort

      mregout(:,:) = 0 
      ent=2
      do ii=1,DIMPC
         do i=2,nterms
         
            isum = 0
            do j=1,DIM
               isum = isum + mreg(i,j)
            end do
            
            if (isum.eq.ii) then
               mregout(ent,:)=mreg(i,:)
               mreg(i,:)=0
               ent=ent+1
            end if
        
         end do
      end do

   
      end subroutine multidx


      subroutine combination(n,m,l)
      implicit none
      ! l = nCm = (n!)/(m!)/((n-m)!)
      integer, intent(in)  :: n,m
      integer, intent(out) :: l
      integer :: i,i1,i2,i3,nbig,npet 

      if(n.le.0.or.m.le.0) stop'n<1 or m<1 in combination'
      if(n.lt.m) stop'n<m in combination'
      

      if(m.ne.n-m)then
         nbig = max(m,n-m)
         npet = min(m,n-m)
      else 
         nbig = m
         npet = n-m
      end if

      i1 = 1
      i2 = 1
      do i=nbig+1,n
         i1 = i1 * i
      end do
      do i=1,npet 
         i2 = i2 * i
      end do
      l = int( dble(i1)/dble(i2) )
  
      end subroutine combination
