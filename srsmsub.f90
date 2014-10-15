!!$c$$$  Subroutines for SRSM: polynomial chaos generation and parameter
!!$c$$$  estimation routines (including collocation and regression)
!!$
!!$
!!$C   Given a vector Xi, this routine computes the values of the terms of
!!$C   polynomial chaos expansions, that could either be exported as an
!!$C   array, or could be used to calculate the sum
!!$C
!!$C   norder is the order of polynomial chaos expansion to be computed
!!$C   nvar is the order of interactive terms desired 
!!$C   for example, nvar = 1 implies that only homogeneous terms are needed
!!$C   nvar = 2 means that terms containing upto two variables are needed
!!$C
!!$C   ictr returns the number of terms in the polynomial chaos expansion
!!$C
!!$C   coeffs contains the coefficients of the polynomial chaos expansion, 
!!$C   if known. Then iopt should be set to 1. If coeffs are not known,
!!$C   this routine helps in formulating the linear equations for the
!!$C   determination of the coefficients of the polynomial chaos expansion.
!!$C   If the coefficients are known, then this can be used to generate 
!!$C   Monte Carlo samples to analyze the distribution.
!!$
!!$C   If the above comments don't convey any meaning, or if they sound 
!!$C   cryptic, please read the documentation that came with this program.
!!$C   If there is no documentation, please contact ... (will fill in details)
!!$C   
!!$
      subroutine chpoly(Xi, n, terms, ictr, norder, nvar)
      include "collsub.h"
  
      dimension Xi(MAXVAR), terms(MAXTRM)
      
C   jctr is the counter for the term in consideration
      jctr = 1

C   "zeroth" order term!
      terms(jctr) = 1.0
      jctr = jctr + 1
      
C   1st order terms
      do i=1,n
        terms(jctr) = Xi(i)
        jctr = jctr + 1
      enddo
      if (norder.lt.2) goto 999

C   2nd order homogeneous terms
      do i=1,n
        terms(jctr) = Xi(i)*Xi(i) - 1
        jctr = jctr + 1
      enddo  

      if (nvar.lt.2) goto 31
      
C   2nd order non-homogeneous terms
      do i=1,n-1
        do j = i+1, n
          terms(jctr) = Xi(i)*Xi(j) 
          jctr = jctr + 1
        enddo  
      enddo

 31   if (norder.lt.3) goto 999

C   3rd order homogeneous terms
      do i=1,n
        terms(jctr) = Xi(i)**3.0 - 3.0*Xi(i)
        jctr = jctr + 1
      enddo  

      if (nvar.lt.2) goto 41

C   3rd order 2-variable terms
      do i=1,n-1
        do j = i+1,n
          terms(jctr) = Xi(i)*Xi(i)*Xi(j) - Xi(j)
          jctr = jctr + 1
          terms(jctr) = Xi(j)*Xi(j)*Xi(i) - Xi(i)
          jctr = jctr + 1
        enddo  
      enddo

      if (nvar.lt.3) goto 41

C   3rd order 3-variable terms
      do i=1,n-2
        do j = i+1,n-1
          do k = j+1,n
            terms(jctr) = Xi(i)*Xi(j)*Xi(k) 
            jctr = jctr + 1
          enddo  
        enddo
      enddo
      
 41   if (norder.lt.4) goto 999

C   4th order homogeneous terms      
      do i=1,n
        terms(jctr) = Xi(i)**4.0 - 6.0*Xi(i)*Xi(i) + 3.0
        jctr = jctr + 1
      enddo  

      if (nvar.lt.2) goto 51

C   4th order 2-variable terms (1st two are asymmetric and next two
C   are symmetric)

      do i=1,n-1
        do j = i+1,n
          terms(jctr) = Xi(i)**3.0*Xi(j) - 3.0*Xi(i)*Xi(j)
          jctr = jctr + 1
          terms(jctr) = Xi(j)**3.0*Xi(i) - 3.0*Xi(i)*Xi(j)
          jctr = jctr + 1
          terms(jctr) = Xi(j)*Xi(j)*Xi(i)*Xi(i) - Xi(i)*Xi(i) -
     #                  Xi(j)*Xi(j) + 1.0
          jctr = jctr + 1
        enddo  
      enddo

      if (nvar.lt.3) goto 51

C   4th order 3-variable terms
      do i=1,n-1
        do j = i+1,n
C   Asymmetric with the third variable, so the third variable = 1,n
          do k = 1,n
            if ((i.ne.k).and.(j.ne.k)) then
              terms(jctr) = Xi(k)*Xi(k)*(Xi(i)*Xi(j) - 1.0)
              jctr = jctr + 1
            endif
          enddo  
        enddo
      enddo

      if (nvar.lt.4) goto 51

C   4th order 4-variable terms
      do i=1,n-3
        do j = i+1,n-2
          do k = j+1,n-1
            do l = k+1,n
              terms(jctr) = Xi(i)*Xi(j)*Xi(k)*Xi(l) 
              jctr = jctr + 1
            enddo
          enddo  
        enddo
      enddo

 51   if (norder.lt.5) goto 999

C   5th order homogeneous terms      
      do i=1,n
        terms(jctr) = Xi(i)**5.0 - 10.0*Xi(i)**3.0 + 15.0*Xi(i)
        jctr = jctr + 1
      enddo  

      if (nvar.lt.2) goto 999

C   5th order 2-variable terms (all are asymmetric)

      do i=1,n-1
        do j = i+1,n
          terms(jctr) = Xi(i)**4.0*Xi(j) - 6.0*Xi(i)*Xi(i)*Xi(j) +
     #                  3*Xi(j)
          jctr = jctr + 1
          terms(jctr) = Xi(j)**4.0*Xi(i) - 6.0*Xi(j)*Xi(j)*Xi(i) +
     #                  3*Xi(i)
          jctr = jctr + 1
          terms(jctr) = (3.0*Xi(i) - 1.0)*(1.0-Xi(j)*Xi(j))*Xi(i)*Xi(i)
          jctr = jctr + 1
          terms(jctr) = (3.0*Xi(j) - 1.0)*(1.0-Xi(i)*Xi(i))*Xi(j)*Xi(j)
          jctr = jctr + 1
        enddo  
      enddo

      if (nvar.lt.3) goto 999

C   5th order 3-variable terms
      do i=1,n-1
        do j = i+1,n
C   Asymmetric with the third variable, so the third variable = 1,n
          do k = 1,n
            if ((i.ne.k).and.(j.ne.k)) then
              terms(jctr) = (Xi(k)**3.0 - 3.0*Xi(k))*Xi(i)*Xi(j)
              jctr = jctr + 1
              terms(jctr) = Xi(k)**3.0 * (1.0 - Xi(i)*Xi(i)) *
     #                      (1.0 - Xi(j)*Xi(j))
              jctr = jctr + 1
            endif
          enddo  
        enddo
      enddo

      if (nvar.lt.4) goto 999

C   5th order 4-variable terms
      do i=1,n-2
        do j = i+1,n-1
          do k = j+1,n
            do l = 1,n
              if ((i.ne.l).and.(j.ne.l).and.(k.ne.l)) then
                terms(jctr) = Xi(l)*Xi(l)*(Xi(i)*Xi(j)*Xi(k) - 1.0)
                jctr = jctr + 1
              endif
            enddo
          enddo  
        enddo
      enddo

      if (nvar.lt.5) goto 999

C   5th order 5-variable terms
      do i=1,n-4
        do j = i+1,n-3
          do k = j+1,n-2
            do l = k+1,n-1
              do m = l+1,n
                terms(jctr) = Xi(i)*Xi(j)*Xi(k)*Xi(l)*Xi(m) 
                jctr = jctr + 1
              enddo
            enddo
          enddo  
        enddo
      enddo

 999  ictr = jctr - 1
c$$$      if (iopt.eq.1) then
c$$$        sum = 0.0
c$$$        do i=1,jctr-1
c$$$          sum = sum + coeffs(i)*terms(i)
c$$$        enddo
c$$$      endif
C      print *, sum, jctr-1

      return
      end


C   Informative comments later...
C   2*MAXTRM is used since we need more runs than the number of terms
C   for the regression method

      subroutine genmat(xin,xout,xcof,ndim,norder,nvar,nouts,nruns,IX)
      include "collsub.h"

      parameter(IDIM=MAXDAT)

      dimension xin(MAXDAT,MAXVAR),xout(MAXDAT,MAXOUT)
      dimension xcof(MAXVAR,MAXTRM)
      dimension indx(MAXDAT), xmat(MAXDAT,MAXDAT)
      dimension xi(MAXVAR), terms(MAXTRM)
      dimension W(MAXDAT), V(MAXDAT,MAXDAT), tmp(MAXDAT), z(MAXDAT)

      do i=1,nruns
        do j=1,ndim
          xi(j) = xin(i,j)
        enddo
        call chpoly(xi,ndim,terms,ictr,norder,nvar)
        do j=1,ictr
          xmat(i,j) = terms(j)
        enddo
      enddo

      if (IX.le.0) then
        call svdcmp(xmat,nruns,ictr,MAXDAT,IDIM,W,V)
        
        do i=1,nouts
          if (IX.eq.-10) then
            do j=1,nruns
              tmp(j) = alog(xout(j,i))
            enddo
          elseif (IX.eq.-20) then
            do j=1,nruns
              tmp(j) = exp(xout(j,i))
            enddo
          else
            do j=1,nruns
              tmp(j) = xout(j,i)
            enddo             
          endif
            
          call svbksb(xmat,W,V,nruns,ictr,MAXDAT,IDIM,tmp,z)
          do j=1,ictr
            xcof(i,j) = z(j)
          enddo
        enddo
      else
        call ludcmp(xmat,ictr,MAXDAT,indx,d)
        do i=1,nouts
          if (IX.eq.10) then
            do j=1,nruns
              tmp(j) = alog(xout(j,i))
            enddo
          elseif (IX.eq.20) then
            do j=1,nruns
              tmp(j) = exp(xout(j,i))
            enddo
          else
            do j=1,nruns
              tmp(j) = xout(j,i)
            enddo             
          endif
        
          call lubksb(xmat,ictr,MAXDAT,indx,tmp)
          do j=1,ictr
            xcof(i,j) = tmp(j)
          enddo
        enddo
      endif

      return
      end

