
SUBROUTINE LEGENDRE(order,X,POLY,DPOLY,ddpoly)

  !       ==========================================================================
  !       Purpose: Compute LEGENDRE polynomials: Pn(x) and their derivatives
  !       Input :  N ---  Order of LEGENDRE POLYNOMIAL
  !                X ---  Argument of LEGENDRE POLYNOMIAL

  !       Output:  POLY(n) --- Pn(x)    -- POLYNOMIAL VALUE
  !                DPOLY(n)--- Pn'(x)   -- DERIVATIVE VAL OF POLYNOMIAL
  !       ==========================================================================

  IMPLICIT NONE

  INTEGER,intent(in) :: order
  integer::K,i,n

  real*8 :: a,b,c

  DOUBLE PRECISION,intent(in) :: X
  real*8               ::PL(0:order+5),DPL(0:ORDER+5),ddpl(0:order+5)
  double precision,intent(out):: POLY(order+1),DPOLY(order+1),DDPOLY(order+1)

  DOUBLE PRECISION ::PL1,DPL1,PL0,DPL0,ddpl1,ddpl0
  
  n=order
  n=n+1

  IF (N.GE.0) THEN

     PL0  = 1.0D0
     DPL0 = 0.0D0
     ddpl0=0.0d0

     PL(0)=PL0
     DPL(0)=DPL0
     ddpl(0)=ddpl0

  END IF

  IF (N.GE.1) THEN    

     PL1  = X 
     DPL1 = 1.0D0  
     ddpl1= 0.0d0

     PL(1)=PL1
     DPL(1)=DPL1
     ddpl(1)=ddpl1

  END IF

  IF (N.GE.2) THEN 

     DO K=2,N+1
        PL(K) =DBLE((2.0d0*k-1)*X*PL(K-1)/k) - DBLE((k-1)*PL(K-2)/k)
     END DO

     DO K=2,N
        DPL(K)= (K+1)*(X*PL(K)-PL(K+1))/(1-X**2)
     END DO

     do k=2,n+1

        !        ddpl(k)=(k+1)( (((k+2)*X**2)+1)*PL(k) -((2*k+5)*x*pl(k+1)) + ((k+2)*pl(k+2)) ) /((X**2)-1)**2


        a=((k+2.0d0)*X**2+1)*pl(k)

        b=-(2.0d0*k+5)*x*pl(k+1)

        c=(k+2.0d0)*pl(k+2)

        ddpl(k)=(a+b+c)*(k+1)/(1-x**2)**2

     end do


  END IF

  do i=1,n
     
     POLY(i)=PL(i-1)
     DPOLY(i)=DPL(i-1)
     DDPOLY(i)=ddpl(i-1)

  end do
! PRINT *,DDPoly(1:N)

!  PRINT *,ddPL(0:N)
!  print *,ddpoly(1:n+1)
  
  !  IF (N.GT.499) STOP 'UNSUPPORTED ORDER IN LEGENDRE POLYNOMIAL. INCREASE THE ALLOCATED MEMORY TO A HIGHER NUMBER IN SUBROUTINE LEGENDRE'

END SUBROUTINE LEGENDRE

SUBROUTINE HERMITE(order,X,POLY,DPOLY,ddpoly,FLAG)

  !       ==========================================================================
  !       Purpose: Compute HERMITE polynomials: Hn(x) or He(x) and their derivatives
  !       Input :  HF --- Function code
  !                       FLAG=1 for PHYSICIST HERMITE
  !                       FLAG=2 for PROBABILIST HERMITE
  !                N ---  Order of HERMITE POLYNOMIAL
  !                X ---  Argument of HERMITE POLYNOMIAL

  !       Output:  POLY(n) --- Hn(x)  or Hen(x) -- POLYNOMIAL VALUE
  !                DPOLY(n)--- Hn'(x) or Hen'(x) --DERIVATIVE VAL OF POLYNOMIAL
  !       ==========================================================================

  IMPLICIT NONE
  !  INTEGER        :: N,K,FLAG
  !  DOUBLE PRECISION :: X,PL(0:N),DPL(0:N),POLY(n+1),DPOLY(n+1)
  !  DOUBLE PRECISION ::PL1,DPL1,PL0,DPL0
  INTEGER        :: order,K,FLAG,i,n

  real*8 :: a,b,c

  DOUBLE PRECISION,intent(in) :: X
  double precision::PL(0:order+5),DPL(0:order+5),ddpl(0:order+5)
  double precision,intent(out):: POLY(order+1),DPOLY(order+1),DDPOLY(order+1)

  DOUBLE PRECISION ::PL1,DPL1,PL0,DPL0,ddpl1,ddpl0

  n=order
  n=n+1

  IF (FLAG.EQ.1) THEN ! PHYSICIST HERMITE not used in our research

     IF (N.GE.0) THEN

        PL0  = 1.0D0
        DPL0 = 0.0D0
        ddpl0= 0.0d0

        PL(0)=PL0
        DPL(0)=DPL0
        ddpl(0)=ddpl0

     END IF

     IF (N.GE.1) THEN    

        PL1  = 2.0D0*X 
        DPL1 = 2.0D0 
        ddpl1 =0.0d0

        PL(1)=PL1
        DPL(1)=DPL1
        ddpl(1)=ddpl1

     END IF

     IF (N.GE.2) THEN 

        DO K=2,N

           PL(K) = 2.0D0*X*PL(K-1)-(2.0D0*(K-1)*PL(K-2))
           DPL(K)= 2.0D0*(K)*PL(K-1)
           ddpl(k)=4.0d0*(k-1)*k*pl(k-2)

        END DO

     END IF

  ELSE IF (FLAG.EQ.2) THEN ! PROBABILIST HERMITEE

     IF (N.GE.0) THEN

        PL0   = 1.0D0
        DPL0  = 0.0D0
        ddpl0 = 0.0d0

        PL(0)=PL0
        DPL(0)=DPL0
        ddpl(0)=ddpl0

     END IF

     IF (N.GE.1) THEN    

        PL1  = X 
        DPL1 = 1.0D0 
        ddpl1 = 0.0d0

        PL(1)=PL1
        DPL(1)=DPL1
        ddpl(1)=ddpl1

     END IF

     IF (N.GE.2) THEN 

        DO K=2,N+3

           PL(K)  = X*PL(K-1)-((K-1)*PL(K-2))
           DPL(K) = (K)*PL(K-1)
           ddpl(k)=(k-1)*k*pl(k-2)
        END DO

     END IF


  ELSE

     STOP'INVALID HERMITE FLAG'

  END IF

  do i=1,n
     
     POLY(i)=PL(i-1)
     DPOLY(i)=DPL(i-1)
     DDPOLY(i)=ddpl(i-1)

  end do

! PRINT *,DDPoly(1:N)

  !  IF (N.GT.499) STOP 'UNSUPPORTED ORDER IN HERMITE POLYNOMIAL. INCREASE THE ALLOCATED MEMORY TO A HIGHER NUMBER IN SUBROUTINE HERMITE'

END SUBROUTINE HERMITE

!!$
!!$SUBROUTINE LEGENDRE(N,X,POLY,DPOLY,ddpoly)
!!$
!!$  !       ==========================================================================
!!$  !       Purpose: Compute LEGENDRE polynomials: Pn(x) and their derivatives
!!$  !       Input :  N ---  Order of LEGENDRE POLYNOMIAL
!!$  !                X ---  Argument of LEGENDRE POLYNOMIAL
!!$
!!$  !       Output:  POLY(n) --- Pn(x)    -- POLYNOMIAL VALUE
!!$  !                DPOLY(n)--- Pn'(x)   -- DERIVATIVE VAL OF POLYNOMIAL
!!$  !       ==========================================================================
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER        :: N,K,FLAG
!!$
!!$  real*8 :: a,b,c
!!$
!!$  DOUBLE PRECISION :: X,PL(0:N),DPL(0:N),ddpl(0:N)
!!$  double precision:: POLY(n+1),DPOLY(n+1),DDPOLY(n+1)
!!$
!!$  DOUBLE PRECISION ::PL1,DPL1,PL0,DPL0,ddpl1,ddpl0
!!$
!!$  IF (N.GE.0) THEN
!!$
!!$     PL0  = 1.0D0
!!$     DPL0 = 0.0D0
!!$     ddpl0=0.0d0
!!$
!!$     PL(0)=PL0
!!$     DPL(0)=DPL0
!!$     ddpl(0)=ddpl0
!!$
!!$  END IF
!!$
!!$  IF (N.GE.1) THEN    
!!$
!!$     PL1  = X 
!!$     DPL1 = 1.0D0  
!!$     ddpl1= 0.0d0
!!$
!!$     PL(1)=PL1
!!$     DPL(1)=DPL1
!!$     ddpl(1)=ddpl1
!!$
!!$  END IF
!!$
!!$  IF (N.GE.2) THEN 
!!$
!!$     DO K=2,N+1
!!$        PL(K) =DBLE((2.0d0*k-1)*X*PL(K-1)/k) - DBLE((k-1)*PL(K-2.0d0)/k)
!!$     END DO
!!$
!!$     DO K=2,N
!!$        DPL(K)= (K+1)*(X*PL(K)-PL(K+1))/(1-X**2)
!!$     END DO
!!$
!!$     do k=2,n
!!$
!!$        !        ddpl(k)=(k+1)( (((k+2)*X**2)+1)*PL(k) -((2*k+5)*x*pl(k+1)) + ((k+2)*pl(k+2)) ) /((X**2)-1)**2
!!$
!!$
!!$        a=((k+2.0d0)*X**2+1)*pl(k)
!!$
!!$        b=-(2.0d0*k+5)*x*pl(k+1)
!!$
!!$        c=(k+2.0d0)*pl(k+2)
!!$
!!$        ddpl(k)=(a+b+c)*(k+1)/(1-x**2)**2
!!$
!!$     end do
!!$
!!$
!!$  END IF
!!$
!!$  POLY(1:N+1)=PL(0:N)
!!$  DPOLY(1:N+1)=DPL(0:N)
!!$  ddpoly(1:n+1)=ddpl(0:N)
!!$
!!$  !PRINT *,PL(0:N)
!!$
!!$  !  IF (N.GT.499) STOP 'UNSUPPORTED ORDER IN LEGENDRE POLYNOMIAL. INCREASE THE ALLOCATED MEMORY TO A HIGHER NUMBER IN SUBROUTINE LEGENDRE'
!!$
!!$END SUBROUTINE LEGENDRE

SUBROUTINE LAGUERRE(N,X,POLY,DPOLY,ddpoly)

  !       ==========================================================================
  !       Purpose: Compute LAGUERRE polynomials: Ln(x) and their derivatives
  !       Input :  N ---  Order of LAGUERRE POLYNOMIAL
  !                X ---  Argument of LAGUERRE POLYNOMIAL

  !       Output:  POLY(n) --- Ln(x)    -- POLYNOMIAL VALUE
  !                DPOLY(n)--- Ln'(x)   -- DERIVATIVE VAL OF POLYNOMIAL
  !       ==========================================================================

  IMPLICIT NONE
  !  INTEGER        :: N,K,FLAG
  !  DOUBLE PRECISION :: X,PL(0:N),DPL(0:N),POLY(n+1),DPOLY(n+1)
  !  DOUBLE PRECISION ::PL1,DPL1,PL0,DPL0

  INTEGER        :: N,K,FLAG

  real*8 :: a,b,c

  DOUBLE PRECISION :: X,PL(0:N),DPL(0:N),ddpl(0:N)
  double precision:: POLY(n+1),DPOLY(n+1),DDPOLY(n+1)

  DOUBLE PRECISION ::PL1,DPL1,PL0,DPL0,ddpl1,ddpl0

  IF (N.GE.0) THEN

     PL0  = 1.0D0
     DPL0 = 0.0D0
     ddpl0=0.0d0

     PL(0)=PL0
     DPL(0)=DPL0
     ddpl(0)=ddpl0

  END IF

  IF (N.GE.1) THEN    

     PL1  = -X+1 
     DPL1 = -1.0D0
     ddpl1 = 0.0d0

     PL(1)=PL1
     DPL(1)=DPL1
     ddpl(1)=ddpl1

  END IF

  IF (N.GE.2) THEN 

     DO K=2,N
        PL(K) = ((2.0d0*k-1-X)*PL(K-1)/K) - ((k-1)*PL(K-2)/K)
        DPL(K)= (K/X)*(PL(K)- PL(K-1))
        ddpl(k)=(pl(k-2))**1
     END DO


  END IF

  POLY(1:N+1)=PL(0:N)
  DPOLY(1:N+1)=DPL(0:N)
  ddpoly(1:n+1)=ddpl(0:N)

  !PRINT *,PL(0:N)

END SUBROUTINE LAGUERRE
