      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
      include 'collsub.h'
      INTEGER m,mp,n,np,NMAX
      REAL b(mp),u(mp,np),v(np,np),w(np),x(np)
      PARAMETER (NMAX=MAXDAT)
      INTEGER i,j,jj
      REAL s,tmp(NMAX)
      do 12 j=1,n
        s=0.
        if(w(j).ne.0.)then
          do 11 i=1,m
            s=s+u(i,j)*b(i)
11        continue
          s=s/w(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0.
        do 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      continue
        x(j)=s
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #,V41.04'21v.
