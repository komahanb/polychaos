      integer MAXVAR,MAXTRM,MAXOUT,MAXDAT,MAXBIN
      parameter (MAXVAR=20)
      parameter (MAXTRM=10000) 
      parameter (MAXOUT=10000)
      parameter (MAXDAT=10000) ! Increase to increase the multiindex
      parameter (MAXBIN=100000)


!   Collocation generation related parameters
      integer MAXDST,MAXPTS,MAXPAR
      parameter (MAXDST=20)
      parameter (MAXPTS=MAXDAT)
      parameter (MAXPAR=20)

!   Currently only seven distributions are considered
      integer MXDNOW
      parameter (MXDNOW=7)

