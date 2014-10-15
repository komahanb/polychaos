      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(n),info
      double precision a(lda,n)
!!$c
!!$c     dgefa factors a double precision matrix by gaussian elimination.
!!$cc
!!$c     on entry
!!$c
!!$c        a       double precision(lda, n)
!!$c                the matrix to be factored.
!!$c
!!$c        lda     integer
!!$c                the leading dimension of the array  a .
!!$c
!!$c        n       integer
!!$c                the order of the matrix  a .
!!$c
!!$c     on return
!!$c
!!$c        a       an upper triangular matrix and the multipliers
!!$c                which were used to obtain it.
!!$c                the factorization can be written  a = l*u  where
!!$c                l  is a product of permutation and unit lower
!!$c                triangular matrices and  u  is upper triangular.
!!$c
!!$c        ipvt    integer(n)
!!$c                an integer vector of pivot indices.
!!$c
!!$c        info    integer
!!$c                = 0  normal value.
!!$c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!!$c                     condition for this subroutine, but it does
!!$c                     indicate that dgesl or dgedi will divide by zero
!!$c                     if called.  use  rcond  in dgeco for a reliable
!!$c                     indication of singularity.
!!$c
!!$c     linpack. this version dated 08/14/78 .
!!$c     cleve moler, university of new mexico, argonne national lab.
!!$c
!!$c     subroutines and functions
!!$c
!!$c     blas daxpy,dscal,idamax
!!$c
!!$c     internal variables
!!$c
      double precision t
      integer idamax,j,k,kp1,l,nm1
!!$c
!!$c
!!$c     gaussian elimination with partial pivoting
!!$c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
!c
!c        find l = pivot index
!c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
!c
!c        zero pivot implies this column already triangularized
!c
         if (a(l,k) .eq. 0.0d0) go to 40
!c
!c           interchange if necessary
!c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
!c
!c           compute multipliers
!c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
!c
!c           row elimination with column indexing
!c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end








      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(n),job
      double precision a(lda,n),b(n)
!!$c
!!$c     dgesl solves the double precision system
!!$c     a * x = b  or  trans(a) * x = b
!!$c     using the factors computed by dgeco or dgefa.
!!$c
!!$c     on entry
!!$c
!!$c        a       double precision(lda, n)
!!$c                the output from dgeco or dgefa.
!!$c
!!$c        lda     integer
!!$c                the leading dimension of the array  a .
!!$c
!!$c        n       integer
!!$c                the order of the matrix  a .
!!$c
!!$c        ipvt    integer(n)
!!$c                the pivot vector from dgeco or dgefa.
!!$c
!!$c        b       double precision(n)
!!$c                the right hand side vector.
!!$c
!!$c        job     integer
!!$c                = 0         to solve  a*x = b ,
!!$c                = nonzero   to solve  trans(a)*x = b  where
!!$c                            trans(a)  is the transpose.
!!$c
!!$c     on return
!!$c
!!$c        b       the solution vector  x .
!!$c
!!$c     error condition
!!$c
!!$c        a division by zero will occur if the input factor contains a
!!$c        zero on the diagonal.  technically this indicates singularity
!!$c        but it is often caused by improper arguments or improper
!!$c        setting of lda .  it will not occur if the subroutines are
!!$c        called correctly and if dgeco has set rcond .gt. 0.0
!!$c        or dgefa has set info .eq. 0 .
!!$c
!!$c     to compute  inverse(a) * c  where  c  is a matrix
!!$c     with  p  columns
!!$c           call dgeco(a,lda,n,ipvt,rcond,z)
!!$c           if (rcond is too small) go to ...
!!$c           do 10 j = 1, p
!!$c              call dgesl(a,lda,n,ipvt,c(1,j),0)
!!$c        10 continue
!!$c
!!$c     linpack. this version dated 08/14/78 .
!!$c     cleve moler, university of new mexico, argonne national lab.
!!$c
!!$c     subroutines and functions
!!$c
!!$c     blas daxpy,ddot
!!$c
!!$c     internal variables
!!$c
      double precision ddot,t
      integer k,kb,l,nm1

      nm1 = n - 1
      if (job .ne. 0) go to 50
!!$c
!!$c        job = 0 , solve  a * x = b
!!$c        first solve  l*y = b
!!$c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
!c
!c        now solve  u*x = y
!c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
!c
!c        job = nonzero, solve  trans(a) * x = b
!c        first solve  trans(u)*y = b
!c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
!c
!c        now solve trans(l)*x = y
!c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end



      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
!*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
!*     ..
!*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!*     ..
!*
!*  Purpose
!*  =======
!*
!*     forms the dot product of two vectors.
!*     uses unrolled loops for increments equal to one.
!*     jack dongarra, linpack, 3/11/78.
!*     modified 12/3/93, array(1) declarations changed to array(*)
!*
!*
!*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MOD
!*     ..
      DDOT = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!*
!*        code for unequal increments or equal increments
!*          not equal to 1
!*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DTEMP = DTEMP + DX(IX)*DY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
!*
!*        code for both increments equal to 1
!*
!*
!*        clean-up loop
!*
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF (N.LT.5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) &
               + DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END




      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
!*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
!*     ..
!*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
!*     ..
!*
!*  Purpose
!*  =======
!*
!*     constant times a vector plus a vector.
!*     uses unrolled loops for increments equal to one.
!*     jack dongarra, linpack, 3/11/78.
!*     modified 12/3/93, array(1) declarations changed to array(*)
!*
!*
!*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MOD
!*     ..
      IF (N.LE.0) RETURN
      IF (DA.EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
!*
!*        code for unequal increments or equal increments
!*          not equal to 1
!*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
!*
!*        code for both increments equal to 1
!*
!*
!*        clean-up loop
!*
   20 M = MOD(N,4)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N.LT.4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
          DY(I) = DY(I) + DA*DX(I)
          DY(I+1) = DY(I+1) + DA*DX(I+1)
          DY(I+2) = DY(I+2) + DA*DX(I+2)
          DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
      END




      SUBROUTINE DSCAL(N,DA,DX,INCX)
!*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
!*     ..
!*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!*     ..
!*
!*  Purpose
!*  =======
!**
!*     scales a vector by a constant.
!*     uses unrolled loops for increment equal to one.
!*     jack dongarra, linpack, 3/11/78.
!*     modified 3/93 to return if incx .le. 0.
!*     modified 12/3/93, array(1) declarations changed to array(*)
!*
!*
!*     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC MOD
!*     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
!*
!*        code for increment not equal to 1
!*
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
          DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
!*
!*        code for increment equal to 1
!*
!*
!*        clean-up loop
!*
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N.LT.5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          DX(I) = DA*DX(I)
          DX(I+1) = DA*DX(I+1)
          DX(I+2) = DA*DX(I+2)
          DX(I+3) = DA*DX(I+3)
          DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END





      INTEGER FUNCTION IDAMAX(N,DX,INCX)
!*     .. Scalar Arguments ..
      INTEGER INCX,N
!*     ..
!*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
!*     ..
!*
!*  Purpose
!*  =======
!*
!*     finds the index of element having max. absolute value.
!*     jack dongarra, linpack, 3/11/78.
!*     modified 3/93 to return if incx .le. 0.
!*     modified 12/3/93, array(1) declarations changed to array(*)
!*
!*
!*     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
!*     ..
!*     .. Intrinsic Functions ..
      INTRINSIC DABS
!*     ..
      IDAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 20
!*
!*        code for increment not equal to 1
!*
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
          IF (DABS(DX(IX)).LE.DMAX) GO TO 5
          IDAMAX = I
          DMAX = DABS(DX(IX))
    5     IX = IX + INCX
   10 CONTINUE
      RETURN
!*
!*        code for increment equal to 1
!*
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
          IF (DABS(DX(I)).LE.DMAX) GO TO 30
          IDAMAX = I
          DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      END
