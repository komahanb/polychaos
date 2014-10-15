      FUNCTION pythag(a,b)
      REAL*8:: a,b,pythag
      REAL*8:: absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1.0d0+(absb/absa)**2)
      else
        if(absb.eq.0.)then
          pythag=0.0d0
        else
          pythag=absb*sqrt(1.0d0+(absa/absb)**2)
        endif
      endif
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software #,V41.04'21v.
