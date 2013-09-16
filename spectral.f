*     VIA CAMPI 213/B, 41100 MODENA, ITALY                         
*     E-MAIL: funaro@unimo.it                                                                      
*****************************************************************       
C                                                                       
      SUBROUTINE GAMMAF(X,GX)                                           
*****************************************************************       
*     COMPUTES THE GAMMA FUNCTION AT A GIVEN POINT                      
*     X = ARGUMENT GREATER THAN 1.E-75 AND SMALLER THAN 57.             
*     GX= VALUE OF GAMMA IN X                                           
*****************************************************************       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION C(11)                                                   
                                                                        
      IF (X .LE. 1.D-75 .OR. X .GE. 57.D0) THEN                         
         WRITE(filenum,*) 'ARGUMENT OUT OF RANGE IN SUBROUTINE GAMMAF'        
         RETURN                                                         
      ENDIF                                                             
                                                                        
         PI  = 3.14159265358979323846D0                                 
         EPS = 1.D-14                                                   
         XX  = X                                                        
         GX  = 1.0D0                                                    
                                                                        
1     IF (DABS(XX-1.D0) .LT. EPS) RETURN                                
      IF (XX .GE. 1.D0) THEN                                            
         XX = XX-1.D0                                                   
         GX = GX*XX                                                     
         GOTO 1                                                         
      ENDIF                                                             
         IND = 0                                                        
      IF (XX .LT. .5D0) THEN                                            
         IND = 1                                                        
         GX  = GX*PI/DSIN(PI*XX)                                        
         XX  = 1.D0-XX                                                  
      ENDIF                                                             
                                                                        
         PR = 1.D0                                                      
         S  = 0.426401432711220868D0                                    
         C(1)  = -0.524741987629368444D0                                
         C(2)  =  0.116154405493589130D0                                
         C(3)  = -0.765978624506602380D-2                               
         C(4)  =  0.899719449391378898D-4                               
         C(5)  = -0.194536980009534621D-7                               
         C(6)  =  0.199382839513630987D-10                              
         C(7)  = -0.204209590209541319D-11                              
         C(8)  =  0.863896817907000175D-13                              
         C(9)  =  0.152237501608472336D-13                              
         C(10) = -0.82572517527771995D-14                               
         C(11) =  0.29973478220522461D-14                               
                                                                        
      DO 2 K=1,11                                                       
         PR = PR*(XX-DFLOAT(K))/(XX+DFLOAT(K-1))                        
         S  = S+C(K)*PR                                                 
2     CONTINUE                                                          
         G  = S*DEXP(1.D0-XX)*(XX+4.5D0)**(XX-.5D0)                     
      IF (IND .EQ. 1) THEN                                              
         GX = GX/G                                                      
      ELSE                                                              
         GX = GX*G                                                      
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE VAJAPO(N,A,B,X,Y,DY,D2Y)                               
************************************************************            
*   COMPUTES THE VALUE OF THE JACOBI POLYNOMIAL OF DEGREE N             
*   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT               
*   N  = DEGREE OF THE POLYNOMIAL                                       
*   A  = PARAMETER >-1                                                  
*   B  = PARAMETER >-1                                                  
*   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED                    
*   Y  = VALUE OF THE POLYNOMIAL IN X                                   
*   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
************************************************************            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
                                                                        
         Y   = 1.D0                                                     
         DY  = 0.D0                                                     
         D2Y = 0.D0                                                     
      IF (N .EQ. 0) RETURN                                              
                                                                        
         AB  = A+B                                                      
         Y   = .5D0*(AB+2.D0)*X+.5D0*(A-B)                              
         DY  = .5D0*(AB+2.D0)                                           
         D2Y = 0.D0                                                     
      IF(N .EQ. 1) RETURN                                               
                                                                        
         YP   = 1.D0                                                    
         DYP  = 0.D0                                                    
         D2YP = 0.D0                                                    
      DO 1 I=2,N                                                        
         DI = DFLOAT(I)                                                 
         C0 = 2.D0*DI+AB                                                
         C1 = 2.D0*DI*(DI+AB)*(C0-2.D0)                                 
         C2 = (C0-1.D0)*(C0-2.D0)*C0                                    
         C3 = (C0-1.D0)*(A-B)*AB                                        
         C4 = 2.D0*(DI+A-1.D0)*C0*(DI+B-1.D0)                           
         YM = Y                                                         
         Y  = ((C2*X+C3)*Y-C4*YP)/C1                                    
         YP = YM                                                        
         DYM  = DY                                                      
         DY   = ((C2*X+C3)*DY-C4*DYP+C2*YP)/C1                          
         DYP  = DYM                                                     
         D2YM = D2Y                                                     
         D2Y  = ((C2*X+C3)*D2Y-C4*D2YP+2.D0*C2*DYP)/C1                  
         D2YP = D2YM                                                    
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE VALEPO(N,X,Y,DY,D2Y)                                   
**************************************************************          
*   COMPUTES THE VALUE OF THE LEGENDRE POLYNOMIAL OF DEGREE N           
*   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT               
*   N  = DEGREE OF THE POLYNOMIAL                                       
*   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED                    
*   Y  = VALUE OF THE POLYNOMIAL IN X                                   
*   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
**************************************************************          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
                                                                        
         Y   = 1.D0                                                     
         DY  = 0.D0                                                     
         D2Y = 0.D0                                                     
      IF (N .EQ. 0) RETURN                                              
                                                                        
         Y   = X                                                        
         DY  = 1.D0                                                     
         D2Y = 0.D0                                                     
      IF(N .EQ. 1) RETURN                                               
                                                                        
         YP   = 1.D0                                                    
         DYP  = 0.D0                                                    
         D2YP = 0.D0                                                    
      DO 1 I=2,N                                                        
         C1 = DFLOAT(I)                                                 
         C2 = 2.D0*C1-1.D0                                              
         C4 = C1-1.D0                                                   
         YM = Y                                                         
         Y  = (C2*X*Y-C4*YP)/C1                                         
         YP = YM                                                        
         DYM  = DY                                                      
         DY   = (C2*X*DY-C4*DYP+C2*YP)/C1                               
         DYP  = DYM                                                     
         D2YM = D2Y                                                     
         D2Y  = (C2*X*D2Y-C4*D2YP+2.D0*C2*DYP)/C1                       
         D2YP = D2YM                                                    
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE VACHPO(N,X,Y,DY,D2Y)                                   
***************************************************************         
*   COMPUTES THE VALUE OF THE CHEBYSHEV POLYNOMIAL OF DEGREE N          
*   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT               
*   N  = DEGREE OF THE POLYNOMIAL                                       
*   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED                    
*   Y  = VALUE OF THE POLYNOMIAL IN X                                   
*   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
***************************************************************         
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
                                                                        
         Y   = 1.D0                                                     
         DY  = 0.D0                                                     
         D2Y = 0.D0                                                     
      IF (N .EQ. 0) RETURN                                              
                                                                        
         Y   = X                                                        
         DY  = 1.D0                                                     
         D2Y = 0.D0                                                     
      IF (N .EQ. 1) RETURN                                              
                                                                        
         YP   = 1.D0                                                    
         DYP  = 0.D0                                                    
         D2YP = 0.D0                                                    
      DO 1 K=2,N                                                        
         YM  = Y                                                        
         Y   = 2.D0*X*Y-YP                                              
         YP  = YM                                                       
         DYM = DY                                                       
         DY  = 2.D0*X*DY+2.D0*YP-DYP                                    
         DYP = DYM                                                      
         D2YM= D2Y                                                      
         D2Y = 2.D0*X*D2Y+4.D0*DYP-D2YP                                 
         D2YP= D2YM                                                     
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE VALAPO(N,A,X,Y,DY,D2Y)                                 
**************************************************************          
*   COMPUTES THE VALUE OF THE LAGUERRE POLYNOMIAL OF DEGREE N           
*   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT               
*   N  = DEGREE OF THE POLYNOMIAL                                       
*   A  = PARAMETER >-1                                                  
*   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED                    
*   Y  = VALUE OF THE POLYNOMIAL IN X                                   
*   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
**************************************************************          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
                                                                        
         Y   = 1.D0                                                     
         DY  = 0.D0                                                     
         D2Y = 0.D0                                                     
      IF (N .EQ. 0) RETURN                                              
                                                                        
         Y   = 1.D0+A-X                                                 
         DY  = -1.D0                                                    
         D2Y = 0.D0                                                     
      IF (N .EQ. 1) RETURN                                              
                                                                        
         YP  = 1.D0                                                     
         DYP = 0.D0                                                     
         D2YP= 0.D0                                                     
      DO 1 K=2,N                                                        
         DK = DFLOAT(K)                                                 
         B1 = (2.D0*DK+A-1.D0-X)/DK                                     
         B2 = (DK+A-1.D0)/DK                                            
         YM = Y                                                         
         Y  = B1*Y-B2*YP                                                
         YP = YM                                                        
         DYM = DY                                                       
         DY  = B1*DY-YP/DK-B2*DYP                                       
         DYP = DYM                                                      
         D2YM= D2Y                                                      
         D2Y = B1*D2Y-2.D0*DYP/DK-B2*D2YP                               
         D2YP= D2YM                                                     
 1    CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE VAHEPO(N,X,Y,DY,D2Y)                                   
*************************************************************           
*   COMPUTES THE VALUE OF THE HERMITE POLYNOMIAL OF DEGREE N            
*   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT               
*   N  = DEGREE OF THE POLYNOMIAL                                       
*   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED                    
*   Y  = VALUE OF THE POLYNOMIAL IN X                                   
*   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
*************************************************************           
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
                                                                        
         Y   = 1.D0                                                     
         DY  = 0.D0                                                     
         D2Y = 0.D0                                                     
      IF (N .EQ. 0) RETURN                                              
                                                                        
         Y   = 2.D0*X                                                   
         DY  = 2.D0                                                     
         D2Y = 0.D0                                                     
      IF (N .EQ. 1) RETURN                                              
                                                                        
         YP=1.D0                                                        
      DO 1 K=2,N                                                        
         DK = DFLOAT(K-1)                                               
         YM  = Y                                                        
         Y   = 2.D0*X*Y-2.D0*DK*YP                                      
         YPM = YP                                                       
         YP  = YM                                                       
1     CONTINUE                                                          
         DN  = 2.D0*DFLOAT(N)                                           
         DNN = 2.D0*DFLOAT(N-1)                                         
         DY  = DN*YP                                                    
         D2Y = DN*DNN*YPM                                               
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE VALASF(N,A,X,Y,DY)                                     
*********************************************************************   
*   COMPUTES THE VALUES OF THE SCALED LAGUERRE FUNCTION OF DEGREE N     
*   AND ITS FIRST DERIVATIVE AT A GIVEN POINT                           
*   N  = DEGREE                                                         
*   A  = PARAMETER >-1                                                  
*   X  = POINT (NON NEGATIVE) IN WHICH THE COMPUTATION IS PERFORMED     
*   Y  = VALUE OF THE FUNCTION IN X                                     
*   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
*********************************************************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
                                                                        
         Y  = 1.D0                                                      
         DY = 0.D0                                                      
      IF (N .EQ. 0) RETURN                                              
                                                                        
         C0 = 1.D0/(4.D0+X)                                             
         C1 = 4.D0*C0/(A+1.D0)                                          
         Y  = C1*(A+1.D0-X)                                             
         DY = -C0*C1*(A+5.D0)                                           
      IF (N .EQ. 1) RETURN                                              
                                                                        
         YP  = 1.D0                                                     
         DYP = 0.D0                                                     
      DO 1 K=2,N                                                        
         DK  = DFLOAT(K)                                                
         DK4 = 4.D0*DK                                                  
         C0  = 1.D0/(DK4+X)                                             
         C1  = DK4+X-2.D0                                               
         C2  = 1.D0/(C1-2.D0)                                           
         C3  = DK4*C0/(DK+A)                                            
         C4  = 2.D0*DK+A-1.D0                                           
         C5  = C4-X                                                     
         C6  = C4+DK4                                                   
         C7  = C2*DFLOAT(4*(K-1)**2)                                    
         DYM = DY                                                       
         DY  = C3*(C5*DY-C0*C6*Y+C7*(2.D0*C0*C1*C2*YP-DYP))             
         DYP = DYM                                                      
         YM  = Y                                                        
         Y   = C3*(C5*Y-C7*YP)                                          
         YP  = YM                                                       
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE VAHESF(N,X,Y,DY)                                       
***********************************************************             
*   COMPUTES THE VALUES OF THE SCALED HERMITE FUNCTION OF               
*   DEGREE N AND ITS FIRST DERIVATIVE AT A GIVEN POINT                  
*   N  = DEGREE                                                         
*   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED                    
*   Y  = VALUE OF THE FUNCTION IN X                                     
*   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
***********************************************************             
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
                                                                        
         Y  = 1.D0                                                      
         DY = 0.D0                                                      
      IF (N .EQ. 0) RETURN                                              
                                                                        
         XX = X*X                                                       
         M  = N/2                                                       
      IF(N .EQ. 2*M) THEN                                               
         A  = -.5D0                                                     
      CALL VALASF(M,A,XX,Y,DY)                                          
         DY = 2.D0*X*DY                                                 
      ELSE                                                              
         A  = .5D0                                                      
      CALL VALASF(M,A,XX,Y,DY)                                          
         DY = Y+2.D0*XX*DY                                              
         Y  = X*Y                                                       
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE ZEJAGA(N,A,B,CS,DZ)                                    
***************************************************************         
*   COMPUTES THE ZEROES OF THE JACOBI POLYNOMIAL OF DEGREE N            
*   N  = THE NUMBER OF ZEROES                                           
*   A  = PARAMETER BETWEEN -1/2 AND 1/2                                 
*   B  = PARAMETER BETWEEN -1/2 AND 1/2                                 
*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N                             
*   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N          
***************************************************************         
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), DZ(1)                                            
      IF (N .EQ .0) RETURN                                              
                                                                        
         AB = A+B                                                       
         CS(1) = (B-A)/(AB+2.D0)                                        
         DZ(1) = .5D0*AB+1.D0                                           
      IF(N .EQ. 1) RETURN                                               
                                                                        
         EPS= 1.E-17                                                    
         PH = 1.57079632679489661923D0                                  
         C  = PH/(2.D0*DFLOAT(N)+AB+1.D0)                               
      DO 1 I=1,N                                                        
         DI  = DFLOAT(I)                                                
         CSX = -DCOS(C*(4.D0*DI+AB-1.D0))                               
      DO 2 IT=1,8                                                       
      CALL VAJAPO(N,A,B,CSX,Y,DY,D2Y)                                   
      IF(DABS(Y) .LT. EPS) GOTO 3                                       
         CSX = CSX-Y/DY                                                 
2     CONTINUE                                                          
3     IF(DABS(CSX) .LT. EPS) CSX=0.D0                                   
         CS(I) = CSX                                                    
         DZ(I) = DY                                                     
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE ZELEGA(N,CS,DZ)                                        
***************************************************************         
*   COMPUTES THE ZEROES OF THE LEGENDRE POLYNOMIAL OF DEGREE N          
*   N  = THE NUMBER OF ZEROES                                           
*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N                             
*   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N          
***************************************************************         
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), DZ(1)                                            
      IF (N .EQ .0) RETURN                                              
                                                                        
         CS(1) = 0.D0                                                   
         DZ(1) = 1.D0                                                   
      IF(N .EQ. 1) RETURN                                               
                                                                        
         N2 = N/2                                                       
         IN = 2*N-4*N2-1                                                
         PH = 1.57079632679489661923D0                                  
         C  = PH/(2.D0*DFLOAT(N)+1.D0)                                  
      DO 1 I=1,N2                                                       
         DI  = DFLOAT(I)                                                
         CSX = DCOS(C*(4.D0*DI-1.D0))                                   
      DO 2 IT=1,8                                                       
      CALL VALEPO(N,CSX,Y,DY,D2Y)                                       
         CSX = CSX-Y/DY                                                 
2     CONTINUE                                                          
         CS(I) = -CSX                                                   
         CS(N-I+1) = CSX                                                
         DZ(I) = DY*DFLOAT(IN)                                          
         DZ(N-I+1) = DY                                                 
1     CONTINUE                                                          
                                                                        
      IF(IN .EQ. -1) RETURN                                             
         CSX = 0.D0                                                     
         CS(N2+1) = CSX                                                 
      CALL VALEPO(N,CSX,Y,DY,D2Y)                                       
         DZ(N2+1) = DY                                                  
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE ZECHGA(N,CS,DZ)                                        
****************************************************************        
*   COMPUTES THE ZEROES OF THE CHEBYSHEV POLYNOMIAL OF DEGREE N         
*   N  = THE NUMBER OF ZEROES                                           
*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N                             
*   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N          
****************************************************************        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), DZ(1)                                            
      IF (N .EQ. 0) RETURN                                              
                                                                        
         CS(1) = 0.D0                                                   
         DZ(1) = 1.D0                                                   
      IF(N .EQ. 1) RETURN                                               
                                                                        
         N2 = N/2                                                       
         IN = 1+4*N2-2*N                                                
         PH = 1.57079632679489661923D0                                  
         DN = DFLOAT(N)                                                 
         C  = PH/DN                                                     
         SI = -1.D0                                                     
      DO 1 I=1,N2                                                       
         DI = DFLOAT(I)                                                 
         CSX = DCOS(C*(2.D0*DI-1.D0))                                   
         CS(I) = -CSX                                                   
         CS(N-I+1) = CSX                                                
         QX = DN/DSQRT(1.D0-CSX*CSX)                                    
         DZ(I) = QX*SI*DFLOAT(IN)                                       
         DZ(N-I+1) = -QX*SI                                             
         SI = -SI                                                       
1     CONTINUE                                                          
                                                                        
      IF(IN .EQ. 1) RETURN                                              
         CS(N2+1) = 0.D0                                                
         N4  = N2/2                                                     
         IN2 = 1+4*N4-2*N2                                              
         DZ(N2+1) = DN*DFLOAT(IN2)                                      
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE ZELAGA(N,A,CS,DZ)                                      
************************************************************************
*   COMPUTES THE ZEROES OF THE LAGUERRE POLYNOMIAL OF DEGREE N          
*   N  = THE NUMBER OF ZEROES                                           
*   A  = PARAMETER >-1                                                  
*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N                             
*   DZ = DERIVATIVES OF THE SCALED FUNCTIONS AT THE ZEROES, DZ(I), I=1,N
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), DZ(1)                                            
      IF (N .EQ. 0) RETURN                                              
                                                                        
         A1 = A+1.D0                                                    
         CS(1) = A1                                                     
         DZ(1) = -4.D0/(A1*(A+5.D0))                                    
      IF (N .EQ. 1) RETURN                                              
                                                                        
         PI = 3.14159265358979323846D0                                  
         DN = DFLOAT(N)                                                 
         C1 = 2.D0*DN+A1                                                
      DO 1 M=1,N                                                        
          DM = DFLOAT(M)                                                
          C2 = 2.D0*(DN+.75D0-DM)*PI/C1                                 
          XN = (C2+PI)/2.D0                                             
      DO 2 IT=1,8                                                       
          XP = XN                                                       
          XN = (DSIN(XP)-XP*DCOS(XP)+C2)/(1.D0-DCOS(XP))                
2     CONTINUE                                                          
          Z  = (DCOS(XN/2.D0))**2                                       
          ZD = 1.D0/(Z-1.D0)                                            
          CSX= 2.D0*C1*Z-((1.25D0*ZD+1.D0)*ZD-1.D0+3.D0*A*A)/(6.D0*C1)  
      DO 3 IT=1,6                                                       
      CALL VALASF(N,A,CSX,Y,DY)                                         
          CSX = CSX-Y/DY                                                
3     CONTINUE                                                          
          CS(M) = CSX                                                   
          DZ(M) = DY                                                    
1     CONTINUE                                                          
                                                                        
4         IN = 0                                                        
      DO 5 M=1,N-1                                                      
      IF(CS(M) .LE. CS(M+1)) GOTO 5                                     
          CSM = CS(M)                                                   
          CS(M) = CS(M+1)                                               
          CS(M+1) = CSM                                                 
          DZM = DZ(M)                                                   
          DZ(M) = DZ(M+1)                                               
          DZ(M+1) = DZM                                                 
          IN = 1                                                        
5     CONTINUE                                                          
      IF(IN .EQ. 1) GOTO 4                                              
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE ZEHEGA(N,CS,DZ)                                        
************************************************************************
*   COMPUTES THE ZEROES OF THE HERMITE POLYNOMIAL OF DEGREE N           
*   N  = THE NUMBER OF ZEROES                                           
*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N                             
*   DZ = DERIVATIVES OF THE SCALED FUNCTIONS AT THE ZEROES, DZ(I), I=1,N
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), DZ(1)                                            
      IF (N .EQ. 0) RETURN 
                                             
         M = N/2                                                        
         CS(M+1) = 0.D0                                                 
         DZ(M+1) = 1.D0                                                 
      IF(N .EQ. 1) RETURN                                               
                                                                        
         IN = 2*N-4*M-1                                                 
      IF(IN .EQ. -1) THEN                                               
         A = -.5D0                                                      
      CALL ZELAGA(M,A,CS,DZ)                                            
      ELSE                                                              
         A = .5D0                                                       
      CALL ZELAGA(M,A,CS,DZ)                                            
      ENDIF                                                             
                                                                        
      DO 1 I=1,M                                                        
         CSX = DSQRT(CS(M-I+1))                                         
         CS(N-I+1) = CSX                                                
      IF(IN .EQ. -1) THEN                                               
         DZ(N-I+1) = 2.D0*CSX*DZ(M-I+1)                                 
      ELSE                                                              
         DZ(N-I+1) = 2.D0*CSX*CSX*DZ(M-I+1)                             
      ENDIF                                                             
1     CONTINUE                                                          
                                                                        
      DO 2 I=1,M                                                        
         CS(I) = -CS(N-I+1)                                             
         DZ(I) = DZ(N-I+1)*DFLOAT(IN)                                   
2     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE ZEJAGL(N,A,B,ET,VN)                                    
********************************************************************    
*   COMPUTES THE NODES RELATIVE TO THE JACOBI GAUSS-LOBATTO FORMULA     
*   N  = ORDER OF THE FORMULA                                           
*   A  = PARAMETER BETWEEN -1/2 AND 1/2                                 
*   B  = PARAMETER BETWEEN -1/2 AND 1/2                                 
*   ET = VECTOR OF THE NODES, ET(I), I=0,N                              
*   VN = VALUES OF THE JACOBI POLYNOMIAL AT THE NODES, VN(I), I=0,N     
********************************************************************    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), VN(0:*)                                        
      IF (N .EQ. 0) RETURN                                              
                                                                        
         ET(0) = -1.D0                                                  
         ET(N) = 1.D0                                                   
         X = -1.D0                                                      
      CALL VAJAPO(N,A,B,X,Y,DY,D2Y)                                     
         VN(0) = Y                                                      
         X = 1.D0                                                       
      CALL VAJAPO(N,A,B,X,Y,DY,D2Y)                                     
         VN(N) = Y                                                      
      IF (N .EQ. 1) RETURN                                              
                                                                        
         EPS= 1.E-17                                                    
         PH = 1.57079632679489661923D0                                  
         AB = A+B                                                       
         DN = DFLOAT(N)                                                 
         C  = PH/(2.D0*DN+AB+1.D0)                                      
         N1 = N-1                                                       
         A1 = A+1.D0                                                    
         B1 = B+1.D0                                                    
      DO 1 I=1,N1                                                       
         DI  = DFLOAT(I)                                                
         ETX = -DCOS(C*(4.D0*DI+AB+1.D0))                               
      DO 2 IT=1,8                                                       
      CALL VAJAPO(N1,A1,B1,ETX,Y,DY,D2Y)                                
      IF(DABS(Y) .LE. EPS) GOTO 3                                       
         ETX = ETX-Y/DY                                                 
2     CONTINUE                                                          
3     IF(DABS(ETX) .LE. EPS) ETX=0.D0                                   
         ET(I) = ETX                                                    
         VN(I) = -.5D0*DY*(1.D0-ETX*ETX)/DN                             
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE ZELEGL(N,ET,VN)                                        
*********************************************************************   
*   COMPUTES THE NODES RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA   
*   N  = ORDER OF THE FORMULA                                           
*   ET = VECTOR OF THE NODES, ET(I), I=0,N                              
*   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N   
*********************************************************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), VN(0:*)                                        
      IF (N .EQ. 0) RETURN                                              
                                                                        
         N2 = (N-1)/2                                                   
         SN = DFLOAT(2*N-4*N2-3)                                        
         ET(0) = -1.D0                                                  
         ET(N) = 1.D0                                                   
         VN(0) = SN                                                     
         VN(N) = 1.D0                                                   
      IF (N .EQ. 1) RETURN                                              
                                                                        
         ET(N2+1) = 0.D0                                                
         X = 0.D0                                                       
      CALL VALEPO(N,X,Y,DY,D2Y)                                         
         VN(N2+1) = Y                                                   
      IF(N .EQ. 2) RETURN                                               
                                                                        
         PI = 3.14159265358979323846D0                                  
         C  = PI/DFLOAT(N)                                              
      DO 1 I=1,N2                                                       
         ETX = DCOS(C*DFLOAT(I))                                        
      DO 2 IT=1,8                                                       
      CALL VALEPO(N,ETX,Y,DY,D2Y)                                       
         ETX = ETX-DY/D2Y                                               
2     CONTINUE                                                          
         ET(I) = -ETX                                                   
         ET(N-I) = ETX                                                  
         VN(I) = Y*SN                                                   
         VN(N-I) = Y                                                    
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE ZECHGL(N,ET)                                           
**********************************************************************  
*   COMPUTES THE NODES RELATIVE TO THE CHEBYSHEV GAUSS-LOBATTO FORMULA  
*   N  = ORDER OF THE FORMULA                                           
*   ET = VECTOR OF THE NODES, ET(I), I=0,N                              
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*)                                                 
      IF (N .EQ. 0) RETURN                                              
                                                                        
         ET(0) = -1.D0                                                  
         ET(N) = 1.D0                                                   
      IF (N .EQ. 1) RETURN                                              
                                                                        
         N2 = (N-1)/2                                                   
         ET(N2+1) = 0.D0                                                
      IF(N .EQ. 2) RETURN                                               
                                                                        
         PI = 3.14159265358979323846D0                                  
         C  = PI/DFLOAT(N)                                              
      DO 1 I=1,N2                                                       
         ETX = DCOS(C*DFLOAT(I))                                        
         ET(I) = -ETX                                                   
         ET(N-I) = ETX                                                  
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE ZELAGR(N,A,ET,VN)                                      
************************************************************************
*   COMPUTES THE NODES RELATIVE TO THE LAGUERRE GAUSS-RADAU FORMULA     
*   N  = ORDER OF THE FORMULA                                           
*   A  = PARAMETER >-1                                                  
*   ET = VECTOR OF THE NODES, ET(I), 0=1,N-1                            
*   VN = SCALED LAGUERRE FUNCTION AT THE NODES, VN(I), I=0,N-1          
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), VN(0:*)                                        
      IF (N .EQ. 0) RETURN                                              
                                                                        
         ET(0) = 0.D0                                                   
         VN(0) = 1.D0                                                   
      IF (N .EQ. 1) RETURN                                              
                                                                        
         A1 = A+1.D0                                                    
         PI = 3.14159265358979323846D0                                  
         N1 = N-1                                                       
         DN = DFLOAT(N1)                                                
         C1 = 2.D0*DN+A1+1.D0                                           
      DO 1 M=1,N1                                                       
          DM = DFLOAT(M)                                                
          C2 = 2.D0*(DN+.75D0-DM)*PI/C1                                 
          XN = (C2+PI)/2.D0                                             
      DO 2 IT=1,8                                                       
          XP = XN                                                       
          XN = (DSIN(XP)-XP*DCOS(XP)+C2)/(1.D0-DCOS(XP))                
2     CONTINUE                                                          
          Z  = (DCOS(XN/2.D0))**2                                       
          ZD = 1.D0/(Z-1.D0)                                            
          ETX= 2.D0*C1*Z-((1.25D0*ZD+1.D0)*ZD-1.D0+3.D0*A1*A1)/(6.D0*C1)
      DO 3 IT=1,6                                                       
      CALL VALASF(N1,A1,ETX,Y,DY)                                       
          ETX = ETX-Y/DY                                                
3     CONTINUE                                                          
          ET(M) = ETX                                                   
      CALL VALASF(N,A,ETX,Y,DY)                                         
          VN(M) = Y                                                     
1     CONTINUE                                                          
                                                                        
4         IN = 0                                                        
      DO 5 M=1,N-2                                                      
      IF(ET(M) .LE. ET(M+1)) GOTO 5                                     
          ETM = ET(M)                                                   
          ET(M) = ET(M+1)                                               
          ET(M+1) = ETM                                                 
          VNM = VN(M)                                                   
          VN(M) = VN(M+1)                                               
          VN(M+1) = VNM                                                 
          IN = 1                                                        
5     CONTINUE                                                          
      IF(IN .EQ. 1) GOTO 4                                              
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE WEJAGA(N,A,B,CS,DZ,WE)                                 
****************************************************************        
*   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS FORMULA           
*   N  = ORDER OF THE FORMULA                                           
*   A  = PARAMETER > -1                                                 
*   B  = PARAMETER > -1                                                 
*   CS = ZEROES OF THE JACOBI POLYNOMIAL, CS(I), I=1,N                  
*   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N          
*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N                            
****************************************************************        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), DZ(1), WE(1)                                     
      IF (N .EQ. 0) RETURN                                              
                                                                        
          AB = A+B+2.D0                                                 
          A2 = A+2.D0                                                   
          B2 = B+2.D0                                                   
      CALL GAMMAF(A2,GA2)                                               
      CALL GAMMAF(B2,GB2)                                               
      CALL GAMMAF(AB,GAB)                                               
          C  = .5D0*(2.D0**AB)*GA2*GB2/GAB                              
      DO 1 M=2,N                                                        
          DM = DFLOAT(M)                                                
          C  = C*(DM+A)*(DM+B)/(DM*(DM+A+B))                            
1     CONTINUE                                                          
                                                                        
      DO 2 I=1,N                                                        
          X  = CS(I)                                                    
          DY = DZ(I)                                                    
          WE(I) = C/((1.D0-X*X)*DY*DY)                                  
2     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE WELEGA(N,CS,DZ,WE)                                     
*****************************************************************       
*   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS FORMULA         
*   N  = ORDER OF THE FORMULA                                           
*   CS = ZEROES OF THE LEGENDRE POLYNOMIAL, CS(I), I=1,N                
*   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N          
*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N                            
*****************************************************************       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), DZ(1), WE(1)                                     
      IF (N .EQ. 0) RETURN                                              
                                                                        
          N2 = N/2                                                      
      DO 1 I=1,N2                                                       
          X  = CS(I)                                                    
          DY = DZ(I)                                                    
          WEX = 2.D0/((1.D0-X*X)*DY*DY)                                 
          WE(I) = WEX                                                   
          WE(N-I+1) = WEX                                               
1     CONTINUE                                                          
                                                                        
      IF(N .EQ. 2*N2) RETURN                                            
          DY = DZ(N2+1)                                                 
          WE(N2+1) = 2.D0/(DY*DY)                                       
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE WECHGA(N,WE)                                           
*****************************************************************       
*   COMPUTES THE WEIGHTS RELATIVE TO THE CHEBYSHEV GAUSS FORMULA        
*   N  = ORDER OF THE FORMULA                                           
*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N                            
*****************************************************************       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION  WE(1)                                                  
      IF (N .EQ. 0) RETURN                                              
                                                                        
         PI = 3.14159265358979323846D0                                  
         C  = PI/DFLOAT(N)                                              
      DO 1 I=1,N                                                        
         WE(I) = C                                                      
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE WELAGA(N,A,CS,WE)                                      
****************************************************************        
*   COMPUTES THE WEIGHTS RELATIVE TO THE LAGUERRE GAUSS FORMULA         
*   N  = ORDER OF THE FORMULA                                           
*   A  = PARAMETER > -1                                                 
*   CS = ZEROES OF THE LAGUERRE POLYNOMIAL, CS(I), I=1,N                
*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N                            
****************************************************************        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), WE(1)                                            
      IF (N .EQ. 0) RETURN                                              
                                                                        
          A1 = A+1.D0                                                   
      CALL GAMMAF(A1,GA1)                                               
          N1 = N+1                                                      
          DN = DFLOAT(N1)                                               
          C = GA1/DN                                                    
      DO 1 M=1,N                                                        
          DM = DFLOAT(M)                                                
          C = C*(DM+A)/(DM+1.D0)                                        
1     CONTINUE                                                          
                                                                        
      DO 2 I=1,N                                                        
          X =CS(I)                                                      
      CALL VALAPO(N1,A,X,Y,DY,D2Y)                                      
          WE(I) = C*X/(Y*Y)                                             
2     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE WEHEGA(N,CS,WE)                                        
****************************************************************        
*   COMPUTES THE WEIGHTS RELATIVE TO THE HERMITE GAUSS FORMULA          
*   N  = ORDER OF THE FORMULA                                           
*   CS = ZEROES OF THE HERMITE POLYNOMIAL, CS(I), I=1,N                 
*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N                            
****************************************************************        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), WE(1)                                            
      IF (N .EQ. 0) RETURN                                              
                                                                        
          PR = 1.77245385090551588D0                                    
          R2 = 1.41421356237309515D0                                    
          C  = PR/DFLOAT(N)                                             
          N2 = N/2                                                      
      DO 1 I=1,N2                                                       
          X  = CS(I)                                                    
          YP = 1.D0                                                     
          Y  = R2*X                                                     
      DO 2 K=2,N-1                                                      
          DK = DFLOAT(K)                                                
          RK = DSQRT(DK)                                                
          QK = DSQRT(DK-1.D0)                                           
          YM = Y                                                        
          Y  = (R2*X*Y-QK*YP)/RK                                        
          YP = YM                                                       
2     CONTINUE                                                          
          WEX = C/(Y*Y)                                                 
          WE(I) = WEX                                                   
          WE(N-I+1) = WEX                                               
1     CONTINUE                                                          
                                                                        
      IF(N .EQ. 2*N2) RETURN                                            
          Y = 1.D0                                                      
      DO 3 K=2,N-1,2                                                    
          DK = DFLOAT(K)                                                
          Y  = Y*DSQRT((DK-1.D0)/DK)                                    
3     CONTINUE                                                          
          WE(N2+1) = C/(Y*Y)                                            
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE WEJAGL(N,A,B,ET,WT)                                    
*********************************************************************   
*   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS-LOBATTO FORMULA   
*   N  = ORDER OF THE FORMULA                                           
*   A  = PARAMETER > -1                                                 
*   B  = PARAMETER > -1                                                 
*   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N                       
*   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N                            
*********************************************************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), WT(0:*)                                        
      IF (N .EQ. 0) RETURN                                              
                                                                        
          A1  = A+1.D0                                                  
          B1  = B+1.D0                                                  
          AB  = A+B                                                     
          AB2 = A1+B1                                                   
      CALL GAMMAF(A1,GA1)                                               
      CALL GAMMAF(B1,GB1)                                               
      CALL GAMMAF(AB2,GAB2)                                             
          C = (2.D0**AB)*GA1*GB1/GAB2                                   
          WT(0) = 2.D0*C*A1/AB2                                         
          WT(N) = 2.D0*C*B1/AB2                                         
      IF (N .EQ. 1) RETURN                                              
                                                                        
          N1 = N-1                                                      
          DN = DFLOAT(N)                                                
          C  = C*(2.D0*DN+AB)/(DN+AB+1.D0)                              
          C1 = C*A1/((B+2.D0)*AB2)                                      
          C2 = C*B1/((A+2.D0)*AB2)                                      
          C3 = .5D0*C*A1*B1                                             
      DO 1 K=1,N-2                                                      
          DK = DFLOAT(K)                                                
          C1 = C1*(DK+A1)*DK/((DK+AB2)*(DK+B+2.D0))                     
          C2 = C2*(DK+B1)*DK/((DK+AB2)*(DK+A+2.D0))                     
          C3 = C3*(DK+A1)*(DK+B1)/((DK+2.D0)*(DK+AB+1.D0))              
1     CONTINUE                                                          
                                                                        
          SU = 0.D0                                                     
      DO 2 M=1,N1                                                       
          SU = SU+ET(M)                                                 
2     CONTINUE                                                          
          WT(0) = C1*(DN-1.D0-SU)                                       
          WT(N) = C2*(DN-1.D0+SU)                                       
      DO 3 I=1,N1                                                       
          X = ET(I)                                                     
      CALL VAJAPO(N,A,B,X,Y,DY,D2Y)                                     
          C4 = -C3/Y                                                    
      CALL VAJAPO(N1,A,B,X,Y,DY,D2Y)                                    
          WT(I) = C4/DY                                                 
3     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE WELEGL(N,ET,VN,WT)                                     
*********************************************************************** 
*   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA 
*   N  = ORDER OF THE FORMULA                                           
*   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N                       
*   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N   
*   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N                            
*********************************************************************** 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), VN(0:*), WT(0:*)                               
      IF (N .EQ. 0) RETURN                                              
                                                                        
          N2 = (N-1)/2                                                  
          DN = DFLOAT(N)                                                
          C  = 2.D0/(DN*(DN+1.D0))                                      
      DO 1 I=0,N2                                                       
          X = ET(I)                                                     
          Y = VN(I)                                                     
          WTX = C/(Y*Y)                                                 
          WT(I) = WTX                                                   
          WT(N-I) = WTX                                                 
1     CONTINUE                                                          
                                                                        
      IF(N-1 .EQ. 2*N2) RETURN                                          
          X = 0.D0                                                      
          Y = VN(N2+1)                                                  
          WT(N2+1) = C/(Y*Y)                                            
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE WECHGL(N,WT)                                           
************************************************************************
*   COMPUTES THE WEIGHTS RELATIVE TO THE CHEBYSHEV GAUSS-LOBATTO FORMULA
*   N  = ORDER OF THE FORMULA                                           
*   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N                            
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION  WT(0:*)                                                
      IF (N .EQ. 0) RETURN                                              
                                                                        
          PI = 3.14159265358979323846D0                                 
          C  = PI/DFLOAT(N)                                             
          C2 = .5D0*C                                                   
          WT(0) =C2                                                     
          WT(N) =C2                                                     
      IF (N .EQ. 1) RETURN                                              
                                                                        
      DO 1 I=1,N-1                                                      
          WT(I) = C                                                     
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE WELAGR(N,A,ET,WT)                                      
*********************************************************************   
*   COMPUTES THE WEIGHTS RELATIVE TO THE LAGUERRE GAUSS-RADAU FORMULA   
*   N  = ORDER OF THE FORMULA                                           
*   A  = PARAMETER > -1                                                 
*   ET = LAGUERRE GAUSS-RADAU NODES, ET(I), I=0,N-1                     
*   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N-1                          
*********************************************************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), WT(0:*)                                        
      IF (N .EQ. 0) RETURN                                              
                                                                        
          A1 = A+1.D0                                                   
      CALL GAMMAF(A1,GA1)                                               
          C1 = GA1                                                      
          WT(0) = C1                                                    
      IF (N .EQ. 1) RETURN                                              
                                                                        
          N1 = N-1                                                      
          C2 = GA1                                                      
      DO 1 K=1,N1                                                       
          DK = DFLOAT(K)                                                
          C1 = C1*DK/(DK+A1)                                            
          C2 = C2*(DK+A)/(DK+1.D0)                                      
1     CONTINUE                                                          
          WT(0) = C1                                                    
      DO 2 I=1,N1                                                       
          X = ET(I)                                                     
      CALL VALAPO(N,A,X,Y,DY,D2Y)                                       
          C3 = C2/Y                                                     
      CALL VALAPO(N1,A,X,Y,DY,D2Y)                                      
          WT(I) = C3/DY                                                 
2     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE WECHCC(N,WK)                                           
*********************************************************************   
*   COMPUTES THE WEIGHTS OF THE CLENSHAW-CURTIS FORMULA OF ORDER 2*N    
*   N  = INTEGER PARAMETER                                              
*   WK = VECTOR OF THE WEIGHTS, WK(I), I=0,2*N                          
*********************************************************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION WK(0:*)                                                 
      IF (N .EQ. 0) RETURN                                              
                                                                        
         PI = 3.14159265358979323846D0                                  
         M  = 2*N                                                       
         DN = DFLOAT(N)                                                 
         DM = DFLOAT(M)                                                 
         C  = 1.D0/(DM*DM-1.D0)                                         
         WK(0) = C                                                      
         WK(M) = C                                                      
         WK(N) = 1.33333333333333333333D0                               
      IF (N .EQ. 1) RETURN                                              
                                                                        
      DO 1 J=1,N-1                                                      
          DJ = DFLOAT(J)                                                
          SU = 1.D0-((-1.D0)**J)*C                                      
      DO 2 K=1,N-1                                                      
          DK = 2.D0*DFLOAT(K)                                           
          SU = SU+2.D0*DCOS(DJ*DK*PI/DM)/(1.D0-DK*DK)                   
2     CONTINUE                                                          
          WK(J) = SU/DN                                                 
          WK(M-J) = SU/DN                                               
1     CONTINUE                                                          
                                                                        
          SU = 1.D0-((-1.D0)**N)*C                                      
      DO 3 K=1,N-1                                                      
          DK = 2.D0*DFLOAT(K)                                           
          SU = SU+2.D0*((-1.D0)**K)/(1.D0-DK*DK)                        
3     CONTINUE                                                          
          WK(N) = SU/DN                                                 
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE INJAGA(N,A,B,CS,DZ,QZ,X,QX)                            
************************************************************************
*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED    
*   BY THE VALUES ATTAINED AT THE ZEROES OF THE JACOBI POLYNOMIAL       
*   N  = THE NUMBER OF ZEROES                                           
*   A  = PARAMETER > -1                                                 
*   B  = PARAMETER > -1                                                 
*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N                             
*   DZ = JACOBI POLYNOMIAL DERIVATIVES AT THE ZEROES, DZ(I), I=1,N      
*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N           
*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED                
*   QX = VALUE OF THE POLYNOMIAL IN X                                   
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), DZ(1), QZ(1)                                     
      IF (N .EQ. 0) RETURN                                              
                                                                        
          EPS = 1.D-14                                                  
      CALL VAJAPO(N,A,B,X,Y,DY,D2Y)                                     
          QX = 0.D0                                                     
      DO 1 J=1,N                                                        
          ED = X-CS(J)                                                  
      IF(DABS(ED) .LT. EPS) THEN                                        
          QX = QZ(J)                                                    
      RETURN                                                            
      ELSE                                                              
          QX = QX+QZ(J)*Y/(DZ(J)*ED)                                    
      ENDIF                                                             
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE INLEGA(N,CS,DZ,QZ,X,QX)                                
************************************************************************
*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED    
*   BY THE VALUES ATTAINED AT THE ZEROES OF THE LEGENDRE POLYNOMIAL     
*   N  = THE NUMBER OF ZEROES                                           
*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N                             
*   DZ = LEGENDRE POLYNOMIAL DERIVATIVES AT THE ZEROES, DZ(I), I=1,N    
*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N           
*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED                
*   QX = VALUE OF THE POLYNOMIAL IN X                                   
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), DZ(1), QZ(1)                                     
      IF (N .EQ. 0) RETURN                                              
                                                                        
          EPS = 1.D-14                                                  
      CALL VALEPO(N,X,Y,DY,D2Y)                                         
          QX = 0.D0                                                     
      DO 1 J=1,N                                                        
          ED = X-CS(J)                                                  
      IF(DABS(ED) .LT. EPS) THEN                                        
          QX = QZ(J)                                                    
      RETURN                                                            
      ELSE                                                              
          QX = QX+QZ(J)*Y/(DZ(J)*ED)                                    
      ENDIF                                                             
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE INCHGA(N,DZ,QZ,X,QX)                                   
************************************************************************
*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED    
*   BY THE VALUES ATTAINED AT THE ZEROES OF THE CHEBYSHEV POLYNOMIAL    
*   N  = THE NUMBER OF ZEROES                                           
*   DZ = CHEBYSHEV POLYNOMIAL DERIVATIVES AT THE ZEROES, DZ(I), I=1,N   
*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N           
*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED                
*   QX = VALUE OF THE POLYNOMIAL IN X                                   
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION  DZ(1), QZ(1)                                           
      IF (N .EQ. 0) RETURN                                              
                                                                        
          EPS = 1.D-14                                                  
          PH = 1.57079632679489661923D0                                 
          DN = DFLOAT(N)                                                
          C  = PH/DN                                                    
      CALL VACHPO(N,X,Y,DY,D2Y)                                         
          QX = 0.D0                                                     
      DO 1 J=1,N                                                        
          ED = X+DCOS(C*(2.D0*DFLOAT(J)-1.D0))                          
      IF(DABS(ED) .LT. EPS) THEN                                        
          QX = QZ(J)                                                    
      RETURN                                                            
      ELSE                                                              
          QX = QX+QZ(J)*Y/(DZ(J)*ED)                                    
      ENDIF                                                             
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE INLAGA(N,A,CS,DZ,QZ,X,QX)                              
*********************************************************************   
*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED    
*   BY THE VALUES ATTAINED AT THE ZEROES OF THE LAGUERRE POLYNOMIAL     
*   N  = THE NUMBER OF ZEROES                                           
*   A  = PARAMETER > -1                                                 
*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N                             
*   DZ = SCALED LAGUERRE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N        
*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N           
*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED                
*   QX = VALUE OF THE POLYNOMIAL IN X                                   
*********************************************************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), DZ(1), QZ(1)                                     
      IF (N .EQ. 0) RETURN                                              
                                                                        
          EPS = 1.D-14                                                  
      CALL VALASF(N,A,X,Y,DY)                                           
          QX = 0.D0                                                     
      DO 1 J=1,N                                                        
          ED = X-CS(J)                                                  
      IF(DABS(ED) .LT. EPS) THEN                                        
          QX = QZ(J)                                                    
      RETURN                                                            
      ELSE                                                              
          PR = 1.D0                                                     
      DO 2 K=1,N                                                        
          DK = 4.D0*DFLOAT(K)                                           
          PR = PR*(DK+X)/(DK+CS(J))                                     
2     CONTINUE                                                          
          QX = QX+QZ(J)*Y*PR/(DZ(J)*ED)                                 
      ENDIF                                                             
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE INHEGA(N,CS,DZ,QZ,X,QX)                                
*********************************************************************   
*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED    
*   BY THE VALUES ATTAINED AT THE ZEROES OF THE HERMITE POLYNOMIAL      
*   N  = THE NUMBER OF ZEROES                                           
*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N                             
*   DZ = SCALED HERMITE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N         
*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N           
*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED                
*   QX = VALUE OF THE POLYNOMIAL IN X                                   
*********************************************************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), DZ(1), QZ(1)                                     
      IF (N .EQ. 0) RETURN                                              
                                                                        
          QX = QZ(1)                                                    
      IF (N .EQ. 1) RETURN                                              
                                                                        
          EPS = 1.D-14                                                  
      CALL VAHESF(N,X,Y,DY)                                             
          QX = 0.D0                                                     
      DO 1 J=1,N                                                        
          ED = X-CS(J)                                                  
      IF(DABS(ED) .LT. EPS) THEN                                        
          QX = QZ(J)                                                    
      RETURN                                                            
      ELSE                                                              
          PR = 1.D0                                                     
      DO 2 K=1,N/2                                                      
          DK = 4.D0*DFLOAT(K)                                           
          PR = PR*(DK+X*X)/(DK+CS(J)**2)                                
2     CONTINUE                                                          
          QX = QX+QZ(J)*Y*PR/(DZ(J)*ED)                                 
      ENDIF                                                             
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE INJAGL(N,A,B,ET,VN,QN,X,QX)                            
*********************************************************************   
*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED    
*   BY THE VALUES ATTAINED AT THE JACOBI GAUSS-LOBATTO NODES            
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   A  = PARAMETER > -1                                                 
*   B  = PARAMETER > -1                                                 
*   ET = VECTOR OF THE NODES, ET(I), I=0,N                              
*   VN = VALUES OF THE JACOBI POLYNOMIAL AT THE NODES, VN(I), I=0,N     
*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N            
*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED                
*   QX = VALUE OF THE POLYNOMIAL IN X                                   
*********************************************************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), VN(0:*), QN(0:*)                               
      IF (N .EQ. 0) RETURN                                              
                                                                        
          EPS = 1.D-14                                                  
      CALL VAJAPO(N,A,B,X,Y,DY,D2Y)                                     
          DN = DFLOAT(N)                                                
          C  = 1.D0/(DN*(DN+A+B+1.D0))                                  
          QX = QN(0)*C*DY*(B+1.D0)*(X-1.D0)/VN(0)                       
          QX = QX+QN(N)*C*DY*(A+1.D0)*(X+1.D0)/VN(N)                    
      IF (N .EQ. 1) RETURN                                              
                                                                        
      DO 1 J=1,N-1                                                      
          ED = X-ET(J)                                                  
      IF(DABS(ED) .LT. EPS) THEN                                        
          QX = QN(J)                                                    
      RETURN                                                            
      ELSE                                                              
          QX = QX+QN(J)*C*DY*(X*X-1.D0)/(VN(J)*ED)                      
      ENDIF                                                             
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE INLEGL(N,ET,VN,QN,X,QX)                                
**********************************************************************  
*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED    
*   BY THE VALUES ATTAINED AT THE LEGENDRE GAUSS-LOBATTO NODES          
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   ET = VECTOR OF THE NODES, ET(I), I=0,N                              
*   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N   
*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N            
*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED                
*   QX = VALUE OF THE POLYNOMIAL IN X                                   
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), VN(0:*), QN(0:*)                               
      IF (N .EQ. 0) RETURN                                              
                                                                        
          EPS = 1.D-14                                                  
      CALL VALEPO(N,X,Y,DY,D2Y)                                         
          DN = DFLOAT(N)                                                
          C  = 1.D0/(DN*(DN+1.D0))                                      
          QX = QN(0)*C*DY*(X-1.D0)/VN(0)                                
          QX = QX+QN(N)*C*DY*(X+1.D0)/VN(N)                             
      IF (N .EQ. 1) RETURN                                              
                                                                        
      DO 1 J=1,N-1                                                      
          ED = X-ET(J)                                                  
      IF(DABS(ED) .LT. EPS) THEN                                        
          QX = QN(J)                                                    
      RETURN                                                            
      ELSE                                                              
          QX = QX+QN(J)*C*DY*(X*X-1.D0)/(VN(J)*ED)                      
      ENDIF                                                             
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE INCHGL(N,QN,X,QX)                                      
*********************************************************************   
*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED    
*   BY THE VALUES ATTAINED AT THE CHEBYSHEV GAUSS-LOBATTO NODES         
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N            
*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED                
*   QX = VALUE OF THE POLYNOMIAL IN X                                   
*********************************************************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION  QN(0:*)                                                
      IF (N .EQ. 0) RETURN                                              
                                                                        
          EPS = 1.D-14                                                  
          PI = 3.14159265358979323846D0                                 
      CALL VACHPO(N,X,Y,DY,D2Y)                                         
          DN = DFLOAT(N)                                                
          SN = DFLOAT(1+4*(N/2)-2*N)                                    
          PN = PI/DN                                                    
          C  = 1.D0/(DN*DN)                                             
          QX = .5D0*SN*QN(0)*C*DY*(X-1.D0)                              
          QX = QX+.5D0*QN(N)*C*DY*(X+1.D0)                              
      IF (N .EQ. 1) RETURN                                              
                                                                        
          SJ = -1.D0                                                    
      DO 1 J=1,N-1                                                      
          ED = X+DCOS(PN*DFLOAT(J))                                     
      IF(DABS(ED) .LT. EPS) THEN                                        
          QX = QN(J)                                                    
      RETURN                                                            
      ELSE                                                              
          QX = QX+QN(J)*SN*SJ*C*DY*(X*X-1.D0)/ED                        
      ENDIF                                                             
          SJ = -SJ                                                      
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE INLAGR(N,A,ET,VN,QN,X,QX)                              
*********************************************************************   
*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED    
*   BY THE VALUES ATTAINED AT THE LAGUERRE GAUSS-RADAU NODES            
*   N  = THE NUMBER OF NODES                                            
*   A  = PARAMETER > -1                                                 
*   ET = VECTOR OF THE NODES, ET(I), I=0,N-1                            
*   VN = SCALED LAGUERRE FUNCTION AT THE NODES, VN(I), I=0,N-1          
*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N-1          
*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED                
*   QX = VALUE OF THE POLYNOMIAL IN X                                   
*********************************************************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), VN(0:*), QN(0:*)                               
      IF (N .EQ. 0) RETURN                                              
                                                                        
          QX = QN(0)                                                    
      IF (N .EQ. 1) RETURN                                              
                                                                        
          EPS = 1.D-14                                                  
          DN = DFLOAT(N)                                                
          C  = -1.D0/DN                                                 
          PR = 1.D0                                                     
          SU = 0.D0                                                     
      DO 1 M=1,N                                                        
          DM = 4.D0*DFLOAT(M)                                           
          PR = PR*(DM+X)/DM                                             
          SU = SU+1.D0/(DM+X)                                           
1     CONTINUE                                                          
      CALL VALASF(N,A,X,Y,DY)                                           
          QX = QN(0)*C*(A+1.D0)*PR*(DY+SU*Y)/VN(0)                      
      DO 2 J=1,N-1                                                      
          ED = X-ET(J)                                                  
      IF(DABS(ED) .LT. EPS) THEN                                        
          QX = QN(J)                                                    
      RETURN                                                            
      ELSE                                                              
          PR = 1.D0                                                     
      DO 3 K=1,N                                                        
          DK = 4.D0*DFLOAT(K)                                           
          PR = PR*(DK+X)/(DK+ET(J))                                     
3     CONTINUE                                                          
          QX = QX+QN(J)*C*PR*(DY+SU*Y)*X/(VN(J)*ED)                     
      ENDIF                                                             
2     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE NOLEGA(N,QZ,WE,QI,QM)                                  
*********************************************************************   
*   COMPUTES THE NORMS OF A POLYNOMIAL DEFINED AT THE LEGENDRE ZEROES   
*   N  = THE NUMBER OF ZEROES                                           
*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N           
*   WE = VECTOR OF THE LEGENDRE GAUSS WEIGHTS, WE(I), I=1,N             
*   QI = INTEGRAL NORM OF THE POLYNOMIAL                                
*   QM = MAXIMUM VALUE OF THE POLYNOMIAL AT THE ZEROES                  
*********************************************************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION QZ(1), WE(1)                                            
      IF (N .EQ. 0) RETURN                                              
                                                                        
          EPS = 1.D-14                                                  
          SU  = 0.D0                                                    
          QM  = 0.D0                                                    
      DO 1 J=1,N                                                        
          Y  = DABS(QZ(J))                                              
      IF(Y .GT. QM) QM=Y                                                
      IF(Y .LT. EPS) GOTO 1                                             
          SU = SU+Y*Y*WE(J)                                             
1     CONTINUE                                                          
          QI = DSQRT(SU)                                                
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE NOCHGA(N,DZ,QZ,WK,QW,QI,QM)                            
**********************************************************************  
*   COMPUTES THE NORMS OF A POLYNOMIAL DEFINED AT THE CHEBYSHEV ZEROES  
*   N  = THE NUMBER OF ZEROES                                           
*   DZ = CHEBYSHEV POLYNOMIAL DERIVATIVES AT THE ZEROES, DZ(I), I=1,N   
*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N           
*   WK = VECTOR OF THE CLENSHAW-CURTIS WEIGHTS, WE(I), I=0,2*N          
*   QW = WEIGHTED INTEGRAL NORM OF THE POLYNOMIAL                       
*   QI = INTEGRAL NORM OF THE POLYNOMIAL                                
*   QM = MAXIMUM VALUE OF THE POLYNOMIAL AT THE ZEROES                  
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION DZ(1), QZ(1), WK(0:*)                                   
      IF (N .EQ. 0) RETURN                                              
                                                                        
          PI = 3.14159265358979323846D0                                 
          DN = DFLOAT(N)                                                
          X  = -1.D0                                                    
      CALL INCHGA(N,DZ,QZ,X,QX)                                         
          S1 = 0.D0                                                     
          S2 = QX*QX*WK(0)                                              
          QM = 0.D0                                                     
      DO 1 J=1,N                                                        
          DJ = DFLOAT(J)                                                
          J2 = 2*J                                                      
          X  = -DCOS(PI*DJ/DN)                                          
          Y  = DABS(QZ(J))                                              
      IF(Y .GT. QM) QM=Y                                                
      CALL INCHGA(N,DZ,QZ,X,QX)                                         
          S1 = S1+Y*Y                                                   
          S2 = S2+Y*Y*WK(J2-1)+QX*QX*WK(J2)                             
1     CONTINUE                                                          
          QW = DSQRT(S1*PI/DN)                                          
          QI = DSQRT(S2)                                                
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE NOJAGL(N,A,B,VN,QN,WT,QW,QS,QM)                        
**********************************************************************  
*   COMPUTES THE NORMS OF A POLYNOMIAL DEFINED AT THE JACOBI GAUSS-     
*   LOBATTO NODES                                                       
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   A  = PARAMETER > -1                                                 
*   B  = PARAMETER > -1                                                 
*   VN = VALUES OF THE JACOBI POLYNOMIAL AT THE NODES, VN(I), I=0,N     
*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N            
*   WT = VECTOR OF THE JACOBI GAUSS-LOBATTO WEIGHTS, WT(I), I=0,N       
*   QW = WEIGHTED INTEGRAL NORM OF THE POLYNOMIAL                       
*   QS = QUADRATURE NORM OF THE POLYNOMIAL                              
*   QM = MAXIMUM VALUE OF THE POLYNOMIAL AT THE NODES                   
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION VN(0:*), QN(0:*), WT(0:*)                               
      IF (N .EQ. 0) RETURN                                              
                                                                        
          EPS = 1.D-14                                                  
          A1  = A+1.D0                                                  
          B1  = B+1.D0                                                  
          AB  = A+B                                                     
          AB2 = AB+2.D0                                                 
          DN = DFLOAT(N)                                                
          C  = ((2.D0)**(AB+1.D0))*(DN+AB+1.D0)/(2.D0*DN+AB+1.D0)       
      CALL GAMMAF(A1,GA1)                                               
      CALL GAMMAF(B1,GB1)                                               
      CALL GAMMAF(AB2,GAB2)                                             
          C  = C*GA1*GB1/GAB2                                           
      DO 1 K=1,N                                                        
          DK = DFLOAT(K)                                                
          C  = C*(DK+A)*(DK+B)/(DK*(DK+AB+1.D0))                        
1     CONTINUE                                                          
                                                                        
          S1 = 0.D0                                                     
          S2 = 0.D0                                                     
          S3 = 0.D0                                                     
          QM = 0.D0                                                     
      DO 2 J=0,N                                                        
          Y1 = QN(J)                                                    
          YM = DABS(Y1)                                                 
      IF(YM .GT. QM) QM=YM                                              
          Y2 = VN(J)                                                    
          S2 = S2+Y1*Y2*WT(J)                                           
          S3 = S3+Y2*Y2*WT(J)                                           
      IF(YM .LT. EPS) GOTO 2                                            
          S1 = S1+Y1*Y1*WT(J)                                           
2     CONTINUE                                                          
          QS = DSQRT(S1)                                                
          QW = DSQRT(S1-(S3-C)*S2*S2/(S3*S3))                           
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE NOLEGL(N,VN,QN,WT,QI,QS,QM)                            
**********************************************************************  
*   COMPUTES THE NORMS OF A POLYNOMIAL DEFINED AT THE LEGENDRE GAUSS-   
*   LOBATTO NODES                                                       
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N   
*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N            
*   WT = VECTOR OF THE LEGENDRE GAUSS-LOBATTO WEIGHTS, WT(I), I=0,N     
*   QW = INTEGRAL NORM OF THE POLYNOMIAL                                
*   QS = QUADRATURE NORM OF THE POLYNOMIAL                              
*   QM = MAXIMUM VALUE OF THE POLYNOMIAL AT THE NODES                   
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION VN(0:*), QN(0:*), WT(0:*)                               
      IF (N .EQ. 0) RETURN                                              
                                                                        
          EPS = 1.D-14                                                  
          DN = DFLOAT(N)                                                
          C  = .5D0*DN*(DN+1.D0)/(2.D0*DN+1.D0)                         
                                                                        
          S1 = 0.D0                                                     
          S2 = 0.D0                                                     
          QM = 0.D0                                                     
      DO 1 J=0,N                                                        
          Y1 = QN(J)                                                    
          YM = DABS(Y1)                                                 
      IF(YM .GT. QM) QM=YM                                              
          Y2 = VN(J)                                                    
          S2 = S2+Y1*Y2*WT(J)                                           
      IF(YM .LT. EPS) GOTO 1                                            
          S1 = S1+Y1*Y1*WT(J)                                           
1     CONTINUE                                                          
          QS = DSQRT(S1)                                                
          QI = DSQRT(DABS(S1-C*S2*S2))                                  
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE NOCHGL(N,QN,WK,QW,QI,QS,QM)                            
**********************************************************************  
*   COMPUTES THE NORMS OF A POLYNOMIAL DEFINED AT THE CHEBYSHEV GAUSS-  
*   LOBATTO NODES                                                       
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N            
*   WK = VECTOR OF THE CLENSHAW-CURTIS WEIGHTS, WE(I), I=0,2*N          
*   QW = WEIGHTED INTEGRAL NORM OF THE POLYNOMIAL                       
*   QI = INTEGRAL NORM OF THE POLYNOMIAL                                
*   QS = QUADRATURE NORM OF THE POLYNOMIAL                              
*   QM = MAXIMUM VALUE OF THE POLYNOMIAL AT THE NODES                   
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION QN(0:*), WK(0:*)                                        
      IF (N .EQ. 0) RETURN                                              
                                                                        
          PI = 3.14159265358979323846D0                                 
          PH = 1.57079632679489661923D0                                 
          DN = DFLOAT(N)                                                
          SN = DFLOAT(1+4*(N/2)-2*N)                                    
          Y  = QN(0)                                                    
          S1 = .5D0*Y*Y                                                 
          S2 = .5D0*Y*SN                                                
          S3 = Y*Y*WK(0)                                                
          QM = DABS(Y)                                                  
                                                                        
          SJ = -1.D0                                                    
      DO 1 J=1,N-1                                                      
          J2 = 2*J                                                      
          DJ = DFLOAT(J2-1)                                             
          X  = -DCOS(PH*DJ/DN)                                          
          Y  = QN(J)                                                    
          YM = DABS(Y)                                                  
      IF(YM .GT. QM) QM=YM                                              
      CALL INCHGL(N,QN,X,QX)                                            
          S1 = S1+Y*Y                                                   
          S2 = S2+Y*SN*SJ                                               
          S3 = S3+QX*QX*WK(J2-1)+Y*Y*WK(J2)                             
          SJ = -SJ                                                      
1     CONTINUE                                                          
          N2 = 2*N                                                      
          DD = DFLOAT(N2-1)                                             
          X  = -DCOS(PH*DD/DN)                                          
          Y  = QN(N)                                                    
          YM = DABS(Y)                                                  
      IF(YM .GT. QM) QM=YM                                              
      CALL INCHGL(N,QN,X,QX)                                            
          S1 = S1+.5D0*Y*Y                                              
          S2 = S2+.5D0*Y                                                
          S3 = S3+QX*QX*WK(N2-1)+Y*Y*WK(N2)                             
                                                                        
          QW = DSQRT(DABS(PI*S1/DN-PH*S2*S2/(DN*DN)))                   
          QI = DSQRT(S3)                                                
          QS = DSQRT(PI*S1/DN)                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE COJAGA(N,A,B,CS,QZ,WE,CO)                              
**********************************************************************  
*   COMPUTES THE JACOBI FOURIER COEFFICIENTS OF A POLYNOMIAL            
*   INDIVIDUATED BY THE VALUES ATTAINED AT THE JACOBI ZEROES            
*   N  = THE NUMBER OF ZEROES                                           
*   A  = PARAMETER >-1                                                  
*   B  = PARAMETER >-1                                                  
*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N                             
*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N           
*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N                            
*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N-1         
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), QZ(1), WE(1), CO(0:*)                            
      IF (N .EQ. 0) RETURN                                              
                                                                        
          A1  = A+1.D0                                                  
          B1  = B+1.D0                                                  
          AB  = A+B                                                     
          AB2 = AB+2.D0                                                 
      CALL GAMMAF(A1,GA1)                                               
      CALL GAMMAF(B1,GB1)                                               
      CALL GAMMAF(AB2,GAB2)                                             
          C  = ((2.D0)**(AB+1.D0))*GA1*GB1/GAB2                         
                                                                        
          SU = 0.D0                                                     
      DO 1 J=1,N                                                        
          SU = SU+QZ(J)*WE(J)                                           
          CO(J-1) = 0.D0                                                
1     CONTINUE                                                          
          CO(0) = SU/C                                                  
      IF (N .EQ. 1) RETURN                                              
                                                                        
      DO 2 J=1,N                                                        
          X  = CS(J)                                                    
          YP = QZ(J)*WE(J)                                              
          Y  = .5D0*YP*(AB2*X+A-B)                                      
      DO 3 K=1,N-1                                                      
          CO(K) = CO(K)+Y                                               
          DK = DFLOAT(K+1)                                              
          CC = 2.D0*DK+AB                                               
          C1 = 2.D0*DK*(DK+AB)*(CC-2.D0)                                
          C2 = (CC-1.D0)*(CC-2.D0)*CC                                   
          C3 = (CC-1.D0)*(A-B)*AB                                       
          C4 = 2.D0*(DK+A-1.D0)*CC*(DK+B-1.D0)                          
          YM = Y                                                        
          Y  = ((C2*X+C3)*Y-C4*YP)/C1                                   
          YP = YM                                                       
3     CONTINUE                                                          
2     CONTINUE                                                          
                                                                        
      DO 4 K=1,N-1                                                      
          DK = DFLOAT(K)                                                
          C  = C*(DK+A)*(DK+B)/DK                                       
          CO(K) = CO(K)*(2.D0*DK+AB+1.D0)/C                             
          C  = C/(DK+AB+1.D0)                                           
4     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE COLEGA(N,CS,QZ,WE,CO)                                  
**********************************************************************  
*   COMPUTES THE LEGENDRE FOURIER COEFFICIENTS OF A POLYNOMIAL          
*   INDIVIDUATED BY THE VALUES ATTAINED AT THE LEGENDRE ZEROES          
*   N  = THE NUMBER OF ZEROES                                           
*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N                             
*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N           
*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N                            
*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N-1         
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), QZ(1), WE(1), CO(0:*)                            
      IF (N .EQ. 0) RETURN                                              
                                                                        
          SU = 0.D0                                                     
      DO 1 J=1,N                                                        
          SU = SU+QZ(J)*WE(J)                                           
          CO(J-1) = 0.D0                                                
1     CONTINUE                                                          
          CO(0) = .5D0*SU                                               
      IF (N .EQ. 1) RETURN                                              
                                                                        
      DO 2 J=1,N                                                        
          X  = CS(J)                                                    
          YP = QZ(J)*WE(J)                                              
          Y  = X*YP                                                     
      DO 3 K=1,N-1                                                      
          CO(K) = CO(K)+Y                                               
          DK = DFLOAT(K+1)                                              
          C1 = 2.D0*DK-1.D0                                             
          C2 = DK-1.D0                                                  
          YM = Y                                                        
          Y  = (C1*X*Y-C2*YP)/DK                                        
          YP = YM                                                       
3     CONTINUE                                                          
2     CONTINUE                                                          
                                                                        
      DO 4 K=1,N-1                                                      
          CO(K) = .5D0*CO(K)*(2.D0*DFLOAT(K)+1.D0)                      
4     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE COCHGA(N,QZ,CO)                                        
**********************************************************************  
*   COMPUTES THE CHEBYSHEV FOURIER COEFFICIENTS OF A POLYNOMIAL         
*   INDIVIDUATED BY THE VALUES ATTAINED AT THE CHEBYSHEV ZEROES         
*   N  = THE NUMBER OF ZEROES                                           
*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N           
*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N-1         
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION  QZ(1), CO(0:*)                                         
      IF (N .EQ. 0) RETURN                                              
                                                                        
          PH = 1.57079632679489661923D0                                 
          DN = DFLOAT(N)                                                
          SU = 0.D0                                                     
      DO 1 J=1,N                                                        
          SU = SU+QZ(J)                                                 
1     CONTINUE                                                          
          CO(0) = SU/DN                                                 
      IF (N .EQ. 1) RETURN                                              
                                                                        
          SK = -2.D0                                                    
      DO 2 K=1,N-1                                                      
          DK = DFLOAT(K)                                                
          SU = 0.D0                                                     
      DO 3 J=1,N                                                        
          DJ = 2.D0*DFLOAT(J)-1.D0                                      
          SU = SU+QZ(J)*DCOS(DK*DJ*PH/DN)                               
3     CONTINUE                                                          
          CO(K) = SK*SU/DN                                              
          SK = -SK                                                      
2     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE COLAGA(N,A,CS,QZ,WE,CO)                                
**********************************************************************  
*   COMPUTES THE LAGUERRE FOURIER COEFFICIENTS OF A POLYNOMIAL          
*   INDIVIDUATED BY THE VALUES ATTAINED AT THE LAGUERRE ZEROES          
*   N  = THE NUMBER OF ZEROES                                           
*   A  = PARAMETER >-1                                                  
*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N                             
*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N           
*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N                            
*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N-1         
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), QZ(1), WE(1), CO(0:*)                            
      IF (N .EQ. 0) RETURN                                              
                                                                        
          A1  = A+1.D0                                                  
      CALL GAMMAF(A1,C)                                                 
          SU = 0.D0                                                     
      DO 1 J=1,N                                                        
          SU = SU+QZ(J)*WE(J)                                           
          CO(J-1) = 0.D0                                                
1     CONTINUE                                                          
          CO(0) = SU/C                                                  
      IF (N .EQ. 1) RETURN                                              
                                                                        
      DO 2 J=1,N                                                        
          X  = CS(J)                                                    
          YP = QZ(J)*WE(J)                                              
          Y  = (A1-X)*YP                                                
      DO 3 K=1,N-1                                                      
          CO(K) = CO(K)+Y                                               
          DK = DFLOAT(K+1)                                              
          B1 = (2.D0*DK+A-1.D0-X)/DK                                    
          B2 = (DK+A-1.D0)/DK                                           
          YM = Y                                                        
          Y  = B1*Y-B2*YP                                               
          YP = YM                                                       
3     CONTINUE                                                          
2     CONTINUE                                                          
                                                                        
      DO 4 K=1,N-1                                                      
          DK = DFLOAT(K)                                                
          C  = C*(DK+A)/DK                                              
          CO(K) = CO(K)/C                                               
4     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE COHEGA(N,CS,QZ,WE,CO)                                  
**********************************************************************  
*   COMPUTES THE HERMITE FOURIER COEFFICIENTS OF A POLYNOMIAL           
*   INDIVIDUATED BY THE VALUES ATTAINED AT THE HERMITE ZEROES           
*   N  = THE NUMBER OF ZEROES                                           
*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N                             
*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N           
*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N                            
*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N-1         
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), QZ(1), WE(1), CO(0:*)                            
      IF (N .EQ. 0) RETURN                                              
                                                                        
          PR = 1.77245385090551588D0                                    
          SU = 0.D0                                                     
      DO 1 J=1,N                                                        
          SU = SU+QZ(J)*WE(J)                                           
          CO(J-1) = 0.D0                                                
1     CONTINUE                                                          
          CO(0) = SU/PR                                                 
      IF (N .EQ. 1) RETURN                                              
                                                                        
      DO 2 J=1,N                                                        
          X  = CS(J)                                                    
          YP = QZ(J)*WE(J)/PR                                           
          Y  = X*YP                                                     
      DO 3 K=1,N-1                                                      
          CO(K) = CO(K)+Y                                               
          DK = DFLOAT(K+1)                                              
          YM = Y                                                        
          Y  = (X*Y-.5D0*YP)/DK                                         
          YP = YM                                                       
3     CONTINUE                                                          
2     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE COJAGL(N,A,B,ET,VN,QN,WT,CO)                           
**********************************************************************  
*   COMPUTES THE JACOBI FOURIER COEFFICIENTS OF A POLYNOMIAL            
*   INDIVIDUATED BY ITS VALUES AT THE JACOBI GAUSS-LOBATTO NODES        
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   A  = PARAMETER >-1                                                  
*   B  = PARAMETER >-1                                                  
*   ET = VECTOR OF THE NODES, ET(I), I=0,N                              
*   VN = VALUES OF THE JACOBI POLYNOMIAL AT THE NODES, VN(I), I=0,N     
*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N            
*   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N                            
*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N           
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), VN(0:*), QN(0:*), WT(0:*), CO(0:*)             
          CO(0) = QN(0)                                                 
      IF (N .EQ. 0) RETURN                                              
                                                                        
          A1  = A+1.D0                                                  
          B1  = B+1.D0                                                  
          AB  = A+B                                                     
          AB2 = AB+2.D0                                                 
          CO(1) = (QN(1)-QN(0))/AB2                                     
          CO(0) = .5D0*(QN(0)+QN(1)-(A-B)*CO(1))                        
      IF (N .EQ. 1) RETURN                                              
                                                                        
      CALL GAMMAF(A1,GA1)                                               
      CALL GAMMAF(B1,GB1)                                               
      CALL GAMMAF(AB2,GAB2)                                             
          C  = ((2.D0)**(AB+1.D0))*GA1*GB1/GAB2                         
          SU = 0.D0                                                     
      DO 1 J=0,N                                                        
          SU = SU+QN(J)*WT(J)                                           
          CO(J) = 0.D0                                                  
1     CONTINUE                                                          
          CO(0) = SU/C                                                  
                                                                        
          CN = 0.D0                                                     
      DO 2 J=0,N                                                        
          X  = ET(J)                                                    
          YP = QN(J)*WT(J)                                              
          Y  = .5D0*YP*(AB2*X+A-B)                                      
          CN = CN+VN(J)*VN(J)*WT(J)                                     
      DO 3 K=1,N                                                        
          CO(K) = CO(K)+Y                                               
          DK = DFLOAT(K+1)                                              
          CC = 2.D0*DK+AB                                               
          C1 = 2.D0*DK*(DK+AB)*(CC-2.D0)                                
          C2 = (CC-1.D0)*(CC-2.D0)*CC                                   
          C3 = (CC-1.D0)*(A-B)*AB                                       
          C4 = 2.D0*(DK+A-1.D0)*CC*(DK+B-1.D0)                          
          YM = Y                                                        
          Y  = ((C2*X+C3)*Y-C4*YP)/C1                                   
          YP = YM                                                       
3     CONTINUE                                                          
2     CONTINUE                                                          
                                                                        
      DO 4 K=1,N-1                                                      
          DK = DFLOAT(K)                                                
          C  = C*(DK+A)*(DK+B)/DK                                       
          CO(K) = CO(K)*(2.D0*DK+AB+1.D0)/C                             
          C  = C/(DK+AB+1.D0)                                           
4     CONTINUE                                                          
          CO(N) = CO(N)/CN                                              
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE COLEGL(N,ET,QN,WT,CO)                                  
**********************************************************************  
*   COMPUTES THE LEGENDRE FOURIER COEFFICIENTS OF A POLYNOMIAL          
*   INDIVIDUATED BY ITS VALUES AT THE LEGENDRE GAUSS-LOBATTO NODES      
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   ET = VECTOR OF THE NODES, ET(I), I=0,N                              
*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N            
*   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N                            
*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N           
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), QN(0:*), WT(0:*), CO(0:*)                      
          CO(0) = QN(0)                                                 
      IF (N .EQ. 0) RETURN                                              
                                                                        
          CO(0) = .5D0*(QN(0)+QN(1))                                    
          CO(1) = .5D0*(QN(1)-QN(0))                                    
      IF (N .EQ. 1) RETURN                                              
                                                                        
          SU = 0.D0                                                     
      DO 1 J=0,N                                                        
          SU = SU+QN(J)*WT(J)                                           
          CO(J) = 0.D0                                                  
1     CONTINUE                                                          
          CO(0) = .5D0*SU                                               
                                                                        
      DO 2 J=0,N                                                        
          X  = ET(J)                                                    
          YP = QN(J)*WT(J)                                              
          Y  = X*YP                                                     
      DO 3 K=1,N                                                        
          CO(K) = CO(K)+Y                                               
          DK = DFLOAT(K+1)                                              
          C1 = 2.D0*DK-1.D0                                             
          C2 = DK-1.D0                                                  
          YM = Y                                                        
          Y  = (C1*X*Y-C2*YP)/DK                                        
          YP = YM                                                       
3     CONTINUE                                                          
2     CONTINUE                                                          
                                                                        
          DN = DFLOAT(N)                                                
          CO(N) = .5D0*DN*CO(N)                                         
      IF (N .EQ. 1) RETURN                                              
                                                                        
      DO 4 K=1,N-1                                                      
          CO(K) = .5D0*CO(K)*(2.D0*DFLOAT(K)+1.D0)                      
4     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE COCHGL(N,QN,CO)                                        
**********************************************************************  
*   COMPUTES THE CHEBYSHEV FOURIER COEFFICIENTS OF A POLYNOMIAL         
*   INDIVIDUATED BY ITS VALUES AT THE CHEBYSHEV GAUSS-LOBATTO NODES     
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N            
*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N           
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION  QN(0:*), CO(0:*)                                       
          CO(0) = QN(0)                                                 
      IF (N .EQ. 0) RETURN                                              
                                                                        
          PI = 3.14159265358979323846D0                                 
          DN = DFLOAT(N)                                                
          DD = DFLOAT(1+4*(N/2)-2*N)                                    
          CO(0) = .5D0*(QN(0)+QN(N))                                    
          CO(N) = .5D0*(QN(0)+DD*QN(N))                                 
      IF (N .EQ. 1) RETURN                                              
                                                                        
          S0 = CO(0)                                                    
          SN = CO(N)                                                    
          SJ = -1.D0                                                    
      DO 1 J=1,N-1                                                      
          S0 = S0+QN(J)                                                 
          SN = SN+QN(J)*SJ                                              
          SJ = -SJ                                                      
1     CONTINUE                                                          
          CO(0) = S0/DN                                                 
          CO(N) = DD*SN/DN                                              
                                                                        
          SK = -1.D0                                                    
      DO 2 K=1,N-1                                                      
          DK = DFLOAT(K)                                                
          SU = .5D0*(QN(0)+QN(N)*SK)                                    
      DO 3 J=1,N-1                                                      
          DJ = DFLOAT(J)                                                
          SU = SU+QN(J)*DCOS(DK*DJ*PI/DN)                               
3     CONTINUE                                                          
          CO(K) = 2.D0*SK*SU/DN                                         
          SK = -SK                                                      
2     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE COLAGR(N,A,ET,QN,WT,CO)                                
**********************************************************************  
*   COMPUTES THE LAGUERRE FOURIER COEFFICIENTS OF A POLYNOMIAL          
*   INDIVIDUATED BY ITS VALUES AT THE LAGUERRE GAUSS-RADAU NODES        
*   N  = THE NUMBER OF NODES                                            
*   A  = PARAMETER >-1                                                  
*   ET = VECTOR OF THE NODES, ET(I), I=0,N-1                            
*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N-1          
*   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N-1                          
*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N-1         
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), QN(0:*), WT(0:*), CO(0:*)                      
      IF (N .EQ. 0) RETURN                                              
                                                                        
          A1  = A+1.D0                                                  
      CALL GAMMAF(A1,C)                                                 
          SU = 0.D0                                                     
      DO 1 J=0,N-1                                                      
          SU = SU+QN(J)*WT(J)                                           
          CO(J) = 0.D0                                                  
1     CONTINUE                                                          
          CO(0) = SU/C                                                  
      IF (N .EQ. 1) RETURN                                              
                                                                        
      DO 2 J=0,N-1                                                      
          X  = ET(J)                                                    
          YP = QN(J)*WT(J)                                              
          Y  = (A1-X)*YP                                                
      DO 3 K=1,N-1                                                      
          CO(K) = CO(K)+Y                                               
          DK = DFLOAT(K+1)                                              
          B1 = (2.D0*DK+A-1.D0-X)/DK                                    
          B2 = (DK+A-1.D0)/DK                                           
          YM = Y                                                        
          Y  = B1*Y-B2*YP                                               
          YP = YM                                                       
3     CONTINUE                                                          
2     CONTINUE                                                          
                                                                        
      DO 4 K=1,N-1                                                      
          DK = DFLOAT(K)                                                
          C  = C*(DK+A)/DK                                              
          CO(K) = CO(K)/C                                               
4     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE PVJAEX(N,A,B,X,CO,Y,DY,D2Y)                            
**********************************************************************  
*   COMPUTES THE VALUE OF A POLYNOMIAL OF DEGREE N AND ITS FIRST AND    
*   SECOND DERIVATIVES BY KNOWING THE JACOBI FOURIER COEFFICIENTS       
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   A  = PARAMETER >-1                                                  
*   B  = PARAMETER >-1                                                  
*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED                
*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N           
*   Y  = VALUE OF THE POLYNOMIAL IN X                                   
*   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CO(0:*)                                                 
          Y   = CO(0)                                                   
          DY  = 0.D0                                                    
          D2Y = 0.D0                                                    
      IF (N .EQ. 0) RETURN                                              
                                                                        
          AB  = A+B                                                     
          P   = .5D0*((AB+2.D0)*X+A-B)                                  
          DP  = .5D0*(AB+2.D0)                                          
          D2P = 0.D0                                                    
          Y   = CO(0)+P*CO(1)                                           
          DY  = DP*CO(1)                                                
          D2Y = 0.D0                                                    
      IF (N .EQ. 1) RETURN                                              
                                                                        
          PP  = 1.D0                                                    
          DPP = 0.D0                                                    
          D2PP = 0.D0                                                   
      DO 1 K=2,N                                                        
          DK = DFLOAT(K)                                                
          CC = 2.D0*DK+AB                                               
          C1 = 2.D0*DK*(DK+AB)*(CC-2.D0)                                
          C2 = (CC-1.D0)*(CC-2.D0)*CC                                   
          C3 = (CC-1.D0)*(A-B)*AB                                       
          C4 = 2.D0*(DK+A-1.D0)*CC*(DK+B-1.D0)                          
          PM = P                                                        
          P  = ((C2*X+C3)*P-C4*PP)/C1                                   
          Y  = Y+P*CO(K)                                                
          PP = PM                                                       
          DPM = DP                                                      
          DP  = ((C2*X+C3)*DP-C4*DPP+C2*PP)/C1                          
          DY  = DY+DP*CO(K)                                             
          DPP  = DPM                                                    
          D2PM = D2P                                                    
          D2P  = ((C2*X+C3)*D2P-C4*D2PP+2.D0*C2*DPP)/C1                 
          D2Y  = D2Y+D2P*CO(K)                                          
          D2PP = D2PM                                                   
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE PVLEEX(N,X,CO,Y,DY,D2Y)                                
**********************************************************************  
*   COMPUTES THE VALUE OF A POLYNOMIAL OF DEGREE N AND ITS FIRST AND    
*   SECOND DERIVATIVES BY KNOWING THE LEGENDRE FOURIER COEFFICIENTS     
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED                
*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N           
*   Y  = VALUE OF THE POLYNOMIAL IN X                                   
*   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CO(0:*)                                                 
          Y   = CO(0)                                                   
          DY  = 0.D0                                                    
          D2Y = 0.D0                                                    
      IF (N .EQ. 0) RETURN                                              
                                                                        
          P   = X                                                       
          DP  = 1.D0                                                    
          D2P = 0.D0                                                    
          Y   = CO(0)+P*CO(1)                                           
          DY  = DP*CO(1)                                                
          D2Y = 0.D0                                                    
      IF (N .EQ. 1) RETURN                                              
                                                                        
          PP  = 1.D0                                                    
          DPP = 0.D0                                                    
          D2PP = 0.D0                                                   
      DO 1 K=2,N                                                        
          DK = DFLOAT(K)                                                
          C2 = 2.D0*DK-1.D0                                             
          C4 = DK-1.D0                                                  
          PM = P                                                        
          P  = (C2*X*P-C4*PP)/DK                                        
          Y  = Y+P*CO(K)                                                
          PP = PM                                                       
          DPM = DP                                                      
          DP  = (C2*X*DP-C4*DPP+C2*PP)/DK                               
          DY  = DY+DP*CO(K)                                             
          DPP  = DPM                                                    
          D2PM = D2P                                                    
          D2P  = (C2*X*D2P-C4*D2PP+2.D0*C2*DPP)/DK                      
          D2Y  = D2Y+D2P*CO(K)                                          
          D2PP = D2PM                                                   
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE PVCHEX(N,X,CO,Y,DY,D2Y)                                
**********************************************************************  
*   COMPUTES THE VALUE OF A POLYNOMIAL OF DEGREE N AND ITS FIRST AND    
*   SECOND DERIVATIVES BY KNOWING THE CHEBYSHEV FOURIER COEFFICIENTS    
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED                
*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N           
*   Y  = VALUE OF THE POLYNOMIAL IN X                                   
*   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CO(0:*)                                                 
          Y   = CO(0)                                                   
          DY  = 0.D0                                                    
          D2Y = 0.D0                                                    
      IF (N .EQ. 0) RETURN                                              
                                                                        
          P   = X                                                       
          DP  = 1.D0                                                    
          D2P = 0.D0                                                    
          Y   = CO(0)+P*CO(1)                                           
          DY  = DP*CO(1)                                                
          D2Y = 0.D0                                                    
      IF (N .EQ. 1) RETURN                                              
                                                                        
          PP  = 1.D0                                                    
          DPP = 0.D0                                                    
          D2PP = 0.D0                                                   
      DO 1 K=2,N                                                        
          PM = P                                                        
          P  = 2.D0*X*P-PP                                              
          Y  = Y+P*CO(K)                                                
          PP = PM                                                       
          DPM = DP                                                      
          DP  = 2.D0*X*DP+2.D0*PP-DPP                                   
          DY  = DY+DP*CO(K)                                             
          DPP  = DPM                                                    
          D2PM = D2P                                                    
          D2P  = 2.D0*X*D2P+4.D0*DPP-D2PP                               
          D2Y  = D2Y+D2P*CO(K)                                          
          D2PP = D2PM                                                   
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE PVLAEX(N,A,X,CO,Y,DY,D2Y)                              
**********************************************************************  
*   COMPUTES THE VALUE OF A POLYNOMIAL OF DEGREE N AND ITS FIRST AND    
*   SECOND DERIVATIVES BY KNOWING THE LAGUERRE FOURIER COEFFICIENTS     
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   A  = PARAMETER >-1                                                  
*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED                
*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N           
*   Y  = VALUE OF THE POLYNOMIAL IN X                                   
*   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CO(0:*)                                                 
          Y   = CO(0)                                                   
          DY  = 0.D0                                                    
          D2Y = 0.D0                                                    
      IF (N .EQ. 0) RETURN                                              
                                                                        
          P   = 1.D0+A-X                                                
          DP  = -1.D0                                                   
          D2P = 0.D0                                                    
          Y   = CO(0)+P*CO(1)                                           
          DY  = DP*CO(1)                                                
          D2Y = 0.D0                                                    
      IF (N .EQ. 1) RETURN                                              
                                                                        
          PP  = 1.D0                                                    
          DPP = 0.D0                                                    
          D2PP = 0.D0                                                   
      DO 1 K=2,N                                                        
          DK = DFLOAT(K)                                                
          B1 = (2.D0*DK+A-1.D0-X)/DK                                    
          B2 = (DK+A-1.D0)/DK                                           
          PM = P                                                        
          P  = B1*P-B2*PP                                               
          Y  = Y+P*CO(K)                                                
          PP = PM                                                       
          DPM = DP                                                      
          DP  = B1*DP-PP/DK-B2*DPP                                      
          DY  = DY+DP*CO(K)                                             
          DPP  = DPM                                                    
          D2PM = D2P                                                    
          D2P  = B1*D2P-2.D0*DPP/DK-B2*D2PP                             
          D2Y  = D2Y+D2P*CO(K)                                          
          D2PP = D2PM                                                   
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE PVHEEX(N,X,CO,Y,DY,D2Y)                                
**********************************************************************  
*   COMPUTES THE VALUE OF A POLYNOMIAL OF DEGREE N AND ITS FIRST AND    
*   SECOND DERIVATIVES BY KNOWING THE HERMITE FOURIER COEFFICIENTS      
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED                
*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N           
*   Y  = VALUE OF THE POLYNOMIAL IN X                                   
*   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
**********************************************************************  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CO(0:*)                                                 
          Y   = CO(0)                                                   
          DY  = 0.D0                                                    
          D2Y = 0.D0                                                    
      IF (N .EQ. 0) RETURN                                              
                                                                        
          P   = 2.D0*X                                                  
          DP  = 2.D0                                                    
          D2P = 0.D0                                                    
          Y   = CO(0)+P*CO(1)                                           
          DY  = DP*CO(1)                                                
          D2Y = 0.D0                                                    
      IF (N .EQ. 1) RETURN                                              
                                                                        
          PP  = 1.D0                                                    
          DPP = 0.D0                                                    
          D2PP = 0.D0                                                   
      DO 1 K=2,N                                                        
          DK = DFLOAT(K)                                                
          PM = P                                                        
          P  = 2.D0*X*P-2.D0*PP*(DK-1.D0)                               
          Y  = Y+P*CO(K)                                                
          DY = DY+2.D0*DK*PM*CO(K)                                      
          D2Y = D2Y+4.D0*DK*(DK-1.D0)*PP*CO(K)                          
          PP  = PM                                                      
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE NOJAEX(N,A,B,CO,QW)                                    
***************************************************************         
*   COMPUTES THE INTEGRAL NORM OF A POLYNOMIAL BY KNOWING THE           
*   FOURIER COEFFICIENTS WITH RESPECT TO THE JACOBI BASIS               
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   A  = PARAMETER >-1                                                  
*   B  = PARAMETER >-1                                                  
*   CO = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N                   
*   QW = WEIGHTED INTEGRAL NORM OF THE POLYNOMIAL                       
***************************************************************         
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CO(0:*)                                                 
                                                                        
          EPS = 1.D-14                                                  
          A1  = A+1.D0                                                  
          B1  = B+1.D0                                                  
          AB  = A+B                                                     
          AB2 = AB+2.D0                                                 
          DN = DFLOAT(N)                                                
      CALL GAMMAF(A1,GA1)                                               
      CALL GAMMAF(B1,GB1)                                               
      CALL GAMMAF(AB2,GAB2)                                             
          C  = ((2.D0)**(AB+1.D0))*GA1*GB1/GAB2                         
          V  = DABS(CO(0))                                              
          QW = V*DSQRT(C)                                               
      IF (N .EQ. 0) RETURN                                              
                                                                        
          SU = 0.D0                                                     
      IF (V .LT. EPS) GOTO 1                                            
          SU = C*V*V                                                    
1     DO 2 K=1,N                                                        
          DK = DFLOAT(K)                                                
          C  = C*(DK+A)*(DK+B)/DK                                       
          V  = DABS(CO(K))                                              
      IF (V .LT. EPS) GOTO 3                                            
          SU = SU+C*V*V/(2.D0*DK+AB+1.D0)                               
3         C  = C/(DK+AB+1.D0)                                           
2     CONTINUE                                                          
          QW = DSQRT(SU)                                                
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE NOLEEX(N,CO,QI)                                        
****************************************************************        
*   COMPUTES THE INTEGRAL NORM OF A POLYNOMIAL BY KNOWING THE           
*   FOURIER COEFFICIENTS WITH RESPECT TO THE LEGENDRE BASIS             
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   CO = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N                   
*   QI = INTEGRAL NORM OF THE POLYNOMIAL                                
****************************************************************        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CO(0:*)                                                 
                                                                        
          EPS = 1.D-14                                                  
          SU = 0.D0                                                     
      DO 1 K=0,N                                                        
          DK = DFLOAT(K)                                                
          V  = DABS(CO(K))                                              
      IF (V .LT. EPS) GOTO 1                                            
          SU = SU+V*V/(2.D0*DK+1.D0)                                    
1     CONTINUE                                                          
          QI = DSQRT(2.D0*SU)                                           
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE NOCHEX(N,CO,QW,QI)                                     
****************************************************************        
*   COMPUTES THE INTEGRAL NORMS OF A POLYNOMIAL BY KNOWING THE          
*   FOURIER COEFFICIENTS WITH RESPECT TO THE CHEBYSHEV BASIS            
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   CO = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N                   
*   QW = WEIGHTED INTEGRAL NORM OF THE POLYNOMIAL                       
*   QI = INTEGRAL NORM OF THE POLYNOMIAL                                
****************************************************************        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CO(0:*)                                                 
                                                                        
          EPS = 1.D-14                                                  
          PR = 1.77245385090551588D0                                    
          R2 = 1.41421356237309515D0                                    
          V  = DABS(CO(0))                                              
          QW = PR*V                                                     
          QI = R2*V                                                     
      IF (N .EQ. 0) RETURN                                              
                                                                        
          SU = 0.D0                                                     
      IF (V .LT. EPS) GOTO 1                                            
          SU = 2.D0*V*V                                                 
1     DO 2 K=1,N                                                        
          V  = DABS(CO(K))                                              
      IF (V .LT. EPS) GOTO 2                                            
          SU = SU+V*V                                                   
2     CONTINUE                                                          
          QW = PR*DSQRT(.5D0*SU)                                        
                                                                        
          SU = 0.D0                                                     
      DO 3 K=0,N,2                                                      
          V  = CO(K)                                                    
      DO 3 M=0,N,2                                                      
          D1 = 1.D0-DFLOAT((K+M)*(K+M))                                 
          D2 = 1.D0-DFLOAT((K-M)*(K-M))                                 
          C  =1.D0/D1+1.D0/D2                                           
          SU = SU+C*V*CO(M)                                             
3     CONTINUE                                                          
      DO 4 K=1,N,2                                                      
          V  = CO(K)                                                    
      DO 4 M=1,N,2                                                      
          D1 = 1.D0-DFLOAT((K+M)*(K+M))                                 
          D2 = 1.D0-DFLOAT((K-M)*(K-M))                                 
          C  =1.D0/D1+1.D0/D2                                           
          SU = SU+C*V*CO(M)                                             
4     CONTINUE                                                          
          QI = DSQRT(SU)                                                
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE NOLAEX(N,A,CO,QW)                                      
*****************************************************************       
*   COMPUTES THE INTEGRAL NORM OF A POLYNOMIAL BY KNOWING THE           
*   FOURIER COEFFICIENTS WITH RESPECT TO THE LAGUERRE BASIS             
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   A  = PARAMETER >-1                                                  
*   CO = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N                   
*   QW = WEIGHTED INTEGRAL NORM OF THE POLYNOMIAL                       
*****************************************************************       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CO(0:*)                                                 
                                                                        
          EPS = 1.D-14                                                  
          A1  = A+1.D0                                                  
      CALL GAMMAF(A1,C)                                                 
          V  = DABS(CO(0))                                              
          QW = V*DSQRT(C)                                               
      IF (N .EQ. 0) RETURN                                              
                                                                        
          SU = 0.D0                                                     
      IF (V .LT. EPS) GOTO 1                                            
          SU = C*V*V                                                    
1     DO 2 K=1,N                                                        
          DK = DFLOAT(K)                                                
          C  = C*(DK+A)/DK                                              
          V  = DABS(CO(K))                                              
      IF (V .LT. EPS) GOTO 2                                            
          SU = SU+C*V*V                                                 
2     CONTINUE                                                          
          QW = DSQRT(SU)                                                
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE NOHEEX(N,CO,QW)                                        
****************************************************************        
*   COMPUTES THE INTEGRAL NORM OF A POLYNOMIAL BY KNOWING THE           
*   FOURIER COEFFICIENTS WITH RESPECT TO THE HERMITE BASIS              
*   N  = THE DEGREE OF THE POLYNOMIAL                                   
*   CO = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N                   
*   QW = WEIGHTED INTEGRAL NORM OF THE POLYNOMIAL                       
****************************************************************        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CO(0:*)                                                 
                                                                        
          EPS = 1.D-14                                                  
          PR = 1.33133536380038953D0                                    
          R2 = 1.41421356237309515D0                                    
          V  = DABS(CO(0))                                              
          QW = V*PR                                                     
      IF (N .EQ. 0) RETURN                                              
                                                                        
          SU = 0.D0                                                     
      IF (V .LT. EPS) GOTO 1                                            
          SU = V*V                                                      
1         C  = 1.D0                                                     
      DO 2 K=1,N                                                        
          DK = DFLOAT(K)                                                
          C  = C*R2*DSQRT(DK)                                           
          V  = DABS(CO(K))                                              
      IF (V .LT. EPS) GOTO 2                                            
          SU = SU+(C*V)*(C*V)                                           
2     CONTINUE                                                          
          QW = PR*DSQRT(SU)                                             
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE COJADE(N,G,CO,CD,CD2)                                  
************************************************************************
*   COMPUTES THE FOURIER COEFFICIENTS OF THE DERIVATIVES OF A POLYNOMIAL
*   FROM ITS FOURIER COEFFICIENTS WITH RESPECT TO THE JACOBI BASIS      
*   N   = THE DEGREE OF THE POLYNOMIAL                                  
*   G   = PARAMETER >-1                                                 
*   CO  = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N                  
*   CD  = COEFFICIENTS OF THE FIRST DERIVATIVE, CD(I), I=0,N            
*   CD2 = COEFFICIENTS OF THE SECOND DERIVATIVE, CD2(I), I=0,N          
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CO(0:*), CD(0:*), CD2(0:*)                              
                                                                        
          CD(N)  = 0.D0                                                 
          CD2(N) = 0.D0                                                 
      IF (N .EQ. 0) RETURN                                              
                                                                        
          CD(0) = (G+1.D0)*CO(1)                                        
          CD2(N-1) = 0.D0                                               
      IF (N .EQ. 1) RETURN                                              
                                                                        
          DN = DFLOAT(N)                                                
          G2 = 2.D0*G                                                   
          CD(N-1) = (2.D0*DN+G2-1.D0)*(DN+G)*CO(N)/(DN+G2)              
      DO 1 K=0,N-2                                                      
          KR = N-K-2                                                    
      IF (KR .NE. 0) THEN                                               
          DK = DFLOAT(KR)                                               
          C1 = (2.D0*DK+G2+1.D0)*(DK+G+1.D0)/(DK+G2+1.D0)               
          C2 = (DK+G+2.D0)/((2.D0*DK+G2+5.D0)*(DK+G2+2.D0))             
          CD(KR)  = C1*(C2*CD(KR+2)+CO(KR+1))                           
          CD2(KR) = C1*(C2*CD2(KR+2)+CD(KR+1))                          
      ELSE                                                              
          CD(0) = .25D0*(G+2.D0)*CD(2)/(G+2.5D0)+(G+1.D0)*CO(1)         
          CD2(0) = .25D0*(G+2.D0)*CD2(2)/(G+2.5D0)+(G+1.D0)*CD(1)       
      ENDIF                                                             
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE COLEDE(N,CO,CD,CD2)                                    
************************************************************************
*   COMPUTES THE FOURIER COEFFICIENTS OF THE DERIVATIVES OF A POLYNOMIAL
*   FROM ITS FOURIER COEFFICIENTS WITH RESPECT TO THE LEGENDRE BASIS    
*   N   = THE DEGREE OF THE POLYNOMIAL                                  
*   CO  = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N                  
*   CD  = COEFFICIENTS OF THE FIRST DERIVATIVE, CD(I), I=0,N            
*   CD2 = COEFFICIENTS OF THE SECOND DERIVATIVE, CD2(I), I=0,N          
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CO(0:*), CD(0:*), CD2(0:*)                              
                                                                        
          CD(N)  = 0.D0                                                 
          CD2(N) = 0.D0                                                 
      IF (N .EQ. 0) RETURN                                              
                                                                        
          DN = DFLOAT(N)                                                
          CD(N-1)  = (2.D0*DN-1.D0)*CO(N)                               
          CD2(N-1) = 0.D0                                               
      IF (N .EQ. 1) RETURN                                              
                                                                        
      DO 1 K=0,N-2                                                      
          KR = N-K-2                                                    
          DK = 2.D0*DFLOAT(KR)+1.D0                                     
          CD(KR)  = DK*(CD(KR+2)/(DK+4.D0)+CO(KR+1))                    
          CD2(KR) = DK*(CD2(KR+2)/(DK+4.D0)+CD(KR+1))                   
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE COCHDE(N,CO,CD,CD2)                                    
************************************************************************
*   COMPUTES THE FOURIER COEFFICIENTS OF THE DERIVATIVES OF A POLYNOMIAL
*   FROM ITS FOURIER COEFFICIENTS WITH RESPECT TO THE CHEBYSHEV BASIS   
*   N   = THE DEGREE OF THE POLYNOMIAL                                  
*   CO  = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N                  
*   CD  = COEFFICIENTS OF THE FIRST DERIVATIVE, CD(I), I=0,N            
*   CD2 = COEFFICIENTS OF THE SECOND DERIVATIVE, CD2(I), I=0,N          
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CO(0:*), CD(0:*), CD2(0:*)                              
                                                                        
          CD(N)  = 0.D0                                                 
          CD2(N) = 0.D0                                                 
      IF (N .EQ. 0) RETURN                                              
                                                                        
          CD(0) = CO(1)                                                 
          CD2(N-1) = 0.D0                                               
      IF (N .EQ. 1) RETURN                                              
                                                                        
          DN = DFLOAT(N)                                                
          CD(N-1) = 2.D0*DN*CO(N)                                       
      DO 1 K=0,N-2                                                      
          KR = N-K-2                                                    
      IF (KR .NE. 0) THEN                                               
          DK = 2.D0*(DFLOAT(KR)+1.D0)                                   
          CD(KR)  = CD(KR+2)+DK*CO(KR+1)                                
          CD2(KR) = CD2(KR+2)+DK*CD(KR+1)                               
      ELSE                                                              
          CD(0)  = .5D0*CD(2)+CO(1)                                     
          CD2(0) = .5D0*CD2(2)+CD(1)                                    
      ENDIF                                                             
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE COLADE(N,CO,CD,CD2)                                    
************************************************************************
*   COMPUTES THE FOURIER COEFFICIENTS OF THE DERIVATIVES OF A POLYNOMIAL
*   FROM ITS FOURIER COEFFICIENTS WITH RESPECT TO THE LAGUERRE BASIS    
*   N   = THE DEGREE OF THE POLYNOMIAL                                  
*   CO  = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N                  
*   CD  = COEFFICIENTS OF THE FIRST DERIVATIVE, CD(I), I=0,N            
*   CD2 = COEFFICIENTS OF THE SECOND DERIVATIVE, CD2(I), I=0,N          
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CO(0:*), CD(0:*), CD2(0:*)                              
                                                                        
          CD(N)  = 0.D0                                                 
          CD2(N) = 0.D0                                                 
      IF (N .EQ. 0) RETURN                                              
                                                                        
          CD(N-1)  = -CO(N)                                             
          CD2(N-1) = 0.D0                                               
      IF (N .EQ. 1) RETURN                                              
                                                                        
      DO 1 K=0,N-2                                                      
          KR = N-K-2                                                    
          CD(KR)  = CD(KR+1)-CO(KR+1)                                   
          CD2(KR) = CD2(KR+2)-CD(KR+1)                                  
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE COHEDE(N,CO,CD,CD2)                                    
************************************************************************
*   COMPUTES THE FOURIER COEFFICIENTS OF THE DERIVATIVES OF A POLYNOMIAL
*   FROM ITS FOURIER COEFFICIENTS WITH RESPECT TO THE HERMITE BASIS     
*   N   = THE DEGREE OF THE POLYNOMIAL                                  
*   CO  = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N                  
*   CD  = COEFFICIENTS OF THE FIRST DERIVATIVE, CD(I), I=0,N            
*   CD2 = COEFFICIENTS OF THE SECOND DERIVATIVE, CD2(I), I=0,N          
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CO(0:*), CD(0:*), CD2(0:*)                              
                                                                        
          CD(N)  = 0.D0                                                 
          CD2(N) = 0.D0                                                 
      IF (N .EQ. 0) RETURN                                              
                                                                        
          DN = DFLOAT(N)                                                
          CD(N-1)  = 2.D0*DN*CO(N)                                      
          CD2(N-1) = 0.D0                                               
      IF (N .EQ. 1) RETURN                                              
                                                                        
      DO 1 K=0,N-2                                                      
          KR = N-K-2                                                    
          DK = 2.D0*DFLOAT(KR)+2.D0                                     
          CD(KR)  = DK*CO(KR+1)                                         
          CD2(KR) = DK*CD(KR+1)                                         
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE DEJAGA(N,A,B,CS,DZ,QZ,DQZ)                             
************************************************************************
*   COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE JACOBI ZEROES FROM   
*   THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS            
*   N   = THE NUMBER OF ZEROES                                          
*   A   = PARAMETER >-1                                                 
*   B   = PARAMETER >-1                                                 
*   CS  = ZEROES OF THE JACOBI POLYNOMIAL, CS(I), I=1,N                 
*   DZ  = JACOBI DERIVATIVES AT THE ZEROES, DZ(I), I=1,N                
*   QZ  = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N          
*   DQZ = DERIVATIVES OF THE POLYNOMIAL AT THE ZEROES, DQZ(I), I=1,N    
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), DZ(1), QZ(1), DQZ(1)                             
      IF (N .EQ. 0) RETURN                                              
                                                                        
      DO 1 I=1,N                                                        
          SU = 0.D0                                                     
          CI = CS(I)                                                    
          DI = DZ(I)                                                    
      DO 2 J=1,N                                                        
      IF (I .NE. J) THEN                                                
          CJ = CS(J)                                                    
          DJ = DZ(J)                                                    
          SU = SU+QZ(J)/(DJ*(CI-CJ))                                    
      ELSE                                                              
          SU = SU+.5D0*QZ(I)*((A+B+2.D0)*CI+A-B)/(DI*(1.D0-CI*CI))      
      ENDIF                                                             
2     CONTINUE                                                          
          DQZ(I) = DI*SU                                                
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE DELAGA(N,A,CS,QZ,DQZ)                                  
************************************************************************
*   COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE LAGUERRE ZEROES FROM 
*   THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS            
*   N   = THE NUMBER OF ZEROES                                          
*   A   = PARAMETER >-1                                                 
*   CS  = ZEROES OF THE LAGUERRE POLYNOMIAL, CS(I), I=1,N               
*   QZ  = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N          
*   DQZ = DERIVATIVES OF THE POLYNOMIAL AT THE ZEROES, DQZ(I), I=1,N    
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), QZ(1), DQZ(1)                                    
      IF (N .EQ. 0) RETURN                                              
                                                                        
      DO 1 J=1,N                                                        
          CJ = CS(J)                                                    
      CALL VALAPO(N,A,CJ,Y,DY,D2Y)                                      
          QZ(J) = QZ(J)/DY                                              
1     CONTINUE                                                          
                                                                        
      DO 2 I=1,N                                                        
          SU = 0.D0                                                     
          CI = CS(I)                                                    
      CALL VALAPO(N,A,CI,Y,DI,D2Y)                                      
      DO 3 J=1,N                                                        
      IF (I .NE. J) THEN                                                
          CJ = CS(J)                                                    
          SU = SU+QZ(J)/(CI-CJ)                                         
      ELSE                                                              
          SU = SU+.5D0*QZ(I)*(CI-A-1.D0)/CI                             
      ENDIF                                                             
3     CONTINUE                                                          
          DQZ(I) = DI*SU                                                
2     CONTINUE                                                          
                                                                        
      DO 4 I=1,N                                                        
          CI = CS(I)                                                    
      CALL VALAPO(N,A,CI,Y,DY,D2Y)                                      
          QZ(I) = DY*QZ(I)                                              
4     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE DEHEGA(N,CS,QZ,DQZ)                                    
************************************************************************
*   COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE HERMITE ZEROES FROM  
*   THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS            
*   N   = THE NUMBER OF ZEROES                                          
*   CS  = ZEROES OF THE HERMITE POLYNOMIAL, CS(I), I=1,N                
*   QZ  = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N          
*   DQZ = DERIVATIVES OF THE POLYNOMIAL AT THE ZEROES, DQZ(I), I=1,N    
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), QZ(1), DQZ(1)                                    
      IF (N .EQ. 0) RETURN                                              
                                                                        
      DO 1 J=1,N                                                        
          CJ = CS(J)                                                    
      CALL VAHEPO(N,CJ,Y,DY,D2Y)                                        
          QZ(J) = QZ(J)/DY                                              
1     CONTINUE                                                          
                                                                        
      DO 2 I=1,N                                                        
          SU = 0.D0                                                     
          CI = CS(I)                                                    
      CALL VAHEPO(N,CI,Y,DI,D2Y)                                        
      DO 3 J=1,N                                                        
      IF (I .NE. J) THEN                                                
          CJ = CS(J)                                                    
          SU = SU+QZ(J)/(CI-CJ)                                         
      ELSE                                                              
          SU = SU+CI*QZ(I)                                              
      ENDIF                                                             
3     CONTINUE                                                          
          DQZ(I) = DI*SU                                                
2     CONTINUE                                                          
                                                                        
      DO 4 I=1,N                                                        
          CI = CS(I)                                                    
      CALL VAHEPO(N,CI,Y,DY,D2Y)                                        
          QZ(I) = DY*QZ(I)                                              
4     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE DEJAGL(N,A,B,ET,VN,QN,DQN)                             
************************************************************************
*   COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE JACOBI GAUSS-LOBATTO 
*   NODES FROM THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS 
*   N   = THE DEGREE OF THE POLYNOMIAL                                  
*   A   = PARAMETER >-1                                                 
*   B   = PARAMETER >-1                                                 
*   ET  = VECTOR OF THE NODES, ET(I), I=0,N                             
*   VN  = VALUES OF THE JACOBI POLYNOMIAL AT THE NODES, VN(I), I=0,N    
*   QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N           
*   DQN = DERIVATIVES OF THE POLYNOMIAL AT THE NODES, DQZ(I), I=0,N     
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), VN(0:*), QN(0:*), DQN(0:*)                     
          DQN(0) = 0.D0                                                 
      IF (N .EQ. 0) RETURN                                              
                                                                        
          A1 = A+1.D0                                                   
          B1 = B+1.D0                                                   
          AB = A+B                                                      
          DN = DFLOAT(N)                                                
          C1 = DN*(DN+AB+1.D0)                                          
          C2 = A1*VN(0)/(B1*VN(N))                                      
          DQN(0) = .5D0*((A-C1)*QN(0)/(B+2.D0)-C2*QN(N))                
          DQN(N) = .5D0*(QN(0)/C2+(C1-B)*QN(N)/(A+2.D0))                
      IF (N .EQ. 1) RETURN                                              
                                                                        
          S1 = DQN(0)                                                   
          S2 = DQN(N)                                                   
          C3 = -VN(0)/B1                                                
          C4 = VN(N)/A1                                                 
      DO 1 J=1,N-1                                                      
          VJ = QN(J)/VN(J)                                              
          EJ = ET(J)                                                    
          S1 = S1+C3*VJ/(1.D0+EJ)                                       
          S2 = S2+C4*VJ/(1.D0-EJ)                                       
1     CONTINUE                                                          
          DQN(0) = S1                                                   
          DQN(N) = S2                                                   
                                                                        
      DO 2 I=1,N-1                                                      
          VI = VN(I)                                                    
          EI = ET(I)                                                    
          SU = B1*QN(0)/((1.D0+EI)*VN(0))-A1*QN(N)/((1.D0-EI)*VN(N))    
      DO 3 J=1,N-1                                                      
      IF (I .NE. J) THEN                                                
          VJ = VN(J)                                                    
          EJ = ET(J)                                                    
          SU = SU+QN(J)/(VJ*(EI-EJ))                                    
      ELSE                                                              
          SU = SU+.5D0*QN(I)*(AB*EI+A-B)/(VI*(1.D0-EI*EI))              
      ENDIF                                                             
3     CONTINUE                                                          
          DQN(I) = VI*SU                                                
2     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE DELEGL(N,ET,VN,QN,DQN)                                 
************************************************************************
*  COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE LEGENDRE GAUSS-LOBATTO
*  NODES FROM THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS  
*   N   = THE DEGREE OF THE POLYNOMIAL                                  
*   ET  = VECTOR OF THE NODES, ET(I), I=0,N                             
*   VN  = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N  
*   QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N           
*   DQN = DERIVATIVES OF THE POLYNOMIAL AT THE NODES, DQZ(I), I=0,N     
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), VN(0:*), QN(0:*), DQN(0:*)                     
          DQN(0) = 0.D0                                                 
      IF (N .EQ. 0) RETURN                                              
                                                                        
      DO 1 I=0,N                                                        
          SU = 0.D0                                                     
          VI = VN(I)                                                    
          EI = ET(I)                                                    
      DO 2 J=0,N                                                        
      IF (I .EQ. J) GOTO 2                                              
          VJ = VN(J)                                                    
          EJ = ET(J)                                                    
          SU = SU+QN(J)/(VJ*(EI-EJ))                                    
2     CONTINUE                                                          
          DQN(I) = VI*SU                                                
1     CONTINUE                                                          
                                                                        
          DN = DFLOAT(N)                                                
          C  = .25D0*DN*(DN+1.D0)                                       
          DQN(0) = DQN(0)-C*QN(0)                                       
          DQN(N) = DQN(N)+C*QN(N)                                       
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE DECHGL(N,ET,QN,DQN)                                    
************************************************************************
* COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE CHEBYSHEV GAUSS-LOBATTO
* NODES FROM THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS   
*   N   = THE DEGREE OF THE POLYNOMIAL                                  
*   ET  = VECTOR OF THE NODES, ET(I), I=0,N                             
*   QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N           
*   DQN = DERIVATIVES OF THE POLYNOMIAL AT THE NODES, DQZ(I), I=0,N     
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), QN(0:*), DQN(0:*)                              
          DQN(0) = 0.D0                                                 
      IF (N .EQ. 0) RETURN                                              
                                                                        
          DN = DFLOAT(N)                                                
          CN = (2.D0*DN*DN+1.D0)/6.D0                                   
          SN = DFLOAT(1+4*(N/2)-2*N)                                    
          DQN(0) = -CN*QN(0)-.5D0*SN*QN(N)                              
          DQN(N) = .5D0*SN*QN(0)+CN*QN(N)                               
      IF (N .EQ. 1) RETURN                                              
                                                                        
          S1 = DQN(0)                                                   
          S2 = DQN(N)                                                   
          SGN = -1.D0                                                   
      DO 1 J=1,N-1                                                      
          EJ = ET(J)                                                    
          QJ = 2.D0*SGN*QN(J)                                           
          S1 = S1-QJ/(1.D0+EJ)                                          
          S2 = S2+QJ*SN/(1.D0-EJ)                                       
          SGN = -SGN                                                    
1     CONTINUE                                                          
          DQN(0) = S1                                                   
          DQN(N) = S2                                                   
                                                                        
          SGNI = -1.D0                                                  
      DO 2 I=1,N-1                                                      
          EI = ET(I)                                                    
          SU = .5D0*SGNI*(QN(0)/(1.D0+EI)-SN*QN(N)/(1.D0-EI))           
          SGNJ = -1.D0                                                  
      DO 3 J=1,N-1                                                      
      IF (I .NE. J) THEN                                                
          EJ = ET(J)                                                    
          SU = SU+SGNI*SGNJ*QN(J)/(EI-EJ)                               
      ELSE                                                              
          SU = SU-.5D0*EI*QN(I)/(1.D0-EI*EI)                            
      ENDIF                                                             
          SGNJ = -SGNJ                                                  
3     CONTINUE                                                          
          DQN(I) = SU                                                   
          SGNI = -SGNI                                                  
2     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE DELAGR(N,A,ET,QN,DQN)                                  
************************************************************************
*   COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE LAGUERRE GAUSS-RADAU 
*   NODES FROM THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS 
*   N   = THE NUMBER OF NODES                                           
*   A   = PARAMETER >-1                                                 
*   ET  = VECTOR OF THE NODES, ET(I), I=0,N-1                           
*   QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N-1         
*   DQN = DERIVATIVES OF THE POLYNOMIAL AT THE NODES, DQZ(I), I=0,N-1   
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), QN(0:*), DQN(0:*)                              
          DN = DFLOAT(N)                                                
          DQN(0) = (1.D0-DN)*QN(0)/(A+2.D0)                             
      IF (N .EQ. 1) RETURN                                              
                                                                        
          A1 = A+1.D0                                                   
          SU = DQN(0)                                                   
          X  = 0.D0                                                     
      CALL VALAPO(N,A,X,Y,DY,D2Y)                                       
          C = Y                                                         
      DO 1 J=1,N-1                                                      
          EJ = ET(J)                                                    
      CALL VALAPO(N,A,EJ,Y,DY,D2Y)                                      
          QN(J) = QN(J)/Y                                               
          SU = SU-C*QN(J)/(A1*EJ)                                       
1     CONTINUE                                                          
          DQN(0) = SU                                                   
                                                                        
      DO 2 I=1,N-1                                                      
          EI = ET(I)                                                    
      CALL VALAPO(N,A,EI,Y,DY,D2Y)                                      
          SU = A1*QN(0)/(C*EI)                                          
      DO 3 J=1,N-1                                                      
      IF (I .NE. J) THEN                                                
          EJ = ET(J)                                                    
          SU = SU+QN(J)/(EI-EJ)                                         
      ELSE                                                              
          SU = SU+.5D0*QN(I)*(EI-A)/EI                                  
      ENDIF                                                             
3     CONTINUE                                                          
          DQN(I) = Y*SU                                                 
2     CONTINUE                                                          
                                                                        
      DO 4 J=1,N-1                                                      
          EJ = ET(J)                                                    
      CALL VALAPO(N,A,EJ,Y,DY,D2Y)                                      
          QN(J) = Y*QN(J)                                               
4     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE DMJAGL(N,NM,A,B,ET,VN,DMA)                             
************************************************************************
*  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE        
*  JACOBI GAUSS-LOBATTO NODES                                           
*  N   = PARAMETER RELATIVE TO THE DIMENSION OF THE MATRIX              
*  NM  = ORDER OF THE MATRIX AS DECLARED IN THE MAIN DIMENSION STATEMENT
*  A   = PARAMETER >-1                                                  
*  B   = PARAMETER >-1                                                  
*  ET  = VECTOR OF THE NODES, ET(I), I=0,N                              
*  VN  = VALUES OF THE JACOBI POLYNOMIAL AT THE NODES, VN(I), I=0,N     
*  DMA = DERIVATIVE MATRIX, DMA(I,J), I=0,N  J=0,N                      
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), VN(0:*), DMA(0:NM,0:*)                         
          DMA(0,0) = 0.D0                                               
      IF (N .EQ. 0) RETURN                                              
                                                                        
          A1 = A+1.D0                                                   
          B1 = B+1.D0                                                   
          AB = A+B                                                      
          DN = DFLOAT(N)                                                
          C1 = DN*(DN+AB+1.D0)                                          
          C2 = A1*VN(0)/(B1*VN(N))                                      
          DMA(0,0) = .5D0*(A-C1)/(B+2.D0)                               
          DMA(N,N) = .5D0*(C1-B)/(A+2.D0)                               
          DMA(0,N) = -.5D0*C2                                           
          DMA(N,0) = .5D0/C2                                            
      IF (N .EQ. 1) RETURN                                              
                                                                        
          C3 = VN(0)/B1                                                 
          C4 = VN(N)/A1                                                 
      DO 1 J=1,N-1                                                      
          VJ = VN(J)                                                    
          EJ = ET(J)                                                    
          DMA(0,J) = -C3/(VJ*(1.D0+EJ))                                 
          DMA(N,J) = C4/(VJ*(1.D0-EJ))                                  
          DMA(J,0) = VJ/(C3*(1.D0+EJ))                                  
          DMA(J,N) = -VJ/(C4*(1.D0-EJ))                                 
1     CONTINUE                                                          
                                                                        
      DO 2 I=1,N-1                                                      
          VI = VN(I)                                                    
          EI = ET(I)                                                    
      DO 3 J=1,N-1                                                      
      IF (I .NE. J) THEN                                                
          VJ = VN(J)                                                    
          EJ = ET(J)                                                    
          DMA(I,J) = VI/(VJ*(EI-EJ))                                    
      ELSE                                                              
          DMA(I,I) = .5D0*(AB*EI+A-B)/(1.D0-EI*EI)                      
      ENDIF                                                             
3     CONTINUE                                                          
2     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE DMLEGL(N,NM,ET,VN,DMA)                                 
************************************************************************
*  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE        
*  LEGENDRE GAUSS-LOBATTO NODES                                         
*  N   = PARAMETER RELATIVE TO THE DIMENSION OF THE MATRIX              
*  NM  = ORDER OF THE MATRIX AS DECLARED IN THE MAIN DIMENSION STATEMENT
*  ET  = VECTOR OF THE NODES, ET(I), I=0,N                              
*  VN  = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N   
*  DMA = DERIVATIVE MATRIX, DMA(I,J), I=0,N  J=0,N                      
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), VN(0:*), DMA(0:NM,0:*)                         
          DMA(0,0) = 0.D0                                               
      IF (N .EQ. 0) RETURN                                              
                                                                        
      DO 1 I=0,N                                                        
          VI = VN(I)                                                    
          EI = ET(I)                                                    
      DO 2 J=0,N                                                        
      IF (I .NE. J) THEN                                                
          VJ = VN(J)                                                    
          EJ = ET(J)                                                    
          DMA(I,J) = VI/(VJ*(EI-EJ))                                    
      ELSE                                                              
          DMA(I,I) = 0.D0                                               
      ENDIF                                                             
2     CONTINUE                                                          
1     CONTINUE                                                          
                                                                        
          DN = DFLOAT(N)                                                
          C  = .25D0*DN*(DN+1.D0)                                       
          DMA(0,0) = -C                                                 
          DMA(N,N) = C                                                  
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE DMCHGL(N,NM,ET,DMA)                                    
************************************************************************
*  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE        
*  CHEBYSHEV GAUSS-LOBATTO NODES                                        
*  N   = PARAMETER RELATIVE TO THE DIMENSION OF THE MATRIX              
*  NM  = ORDER OF THE MATRIX AS DECLARED IN THE MAIN DIMENSION STATEMENT
*  ET  = VECTOR OF THE NODES, ET(I), I=0,N                              
*  DMA = DERIVATIVE MATRIX, DMA(I,J), I=0,N  J=0,N                      
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), DMA(0:NM,0:*)                                  
          DMA(0,0) = 0.D0                                               
      IF (N .EQ. 0) RETURN                                              
                                                                        
          DN = DFLOAT(N)                                                
          CN = (2.D0*DN*DN+1.D0)/6.D0                                   
          SN = DFLOAT(1+4*(N/2)-2*N)                                    
          DMA(0,0) = -CN                                                
          DMA(N,N) = CN                                                 
          DMA(0,N) = -.5D0*SN                                           
          DMA(N,0) = .5D0*SN                                            
      IF (N .EQ. 1) RETURN                                              
                                                                        
          SGN = -1.D0                                                   
      DO 1 J=1,N-1                                                      
          EJ = ET(J)                                                    
          DMA(0,J) = -2.D0*SGN/(1.D0+EJ)                                
          DMA(N,J) = 2.D0*SGN*SN/(1.D0-EJ)                              
          DMA(J,0) = .5D0*SGN/(1.D0+EJ)                                 
          DMA(J,N) = -.5D0*SGN*SN/(1.D0-EJ)                             
          SGN = -SGN                                                    
1     CONTINUE                                                          
                                                                        
          SGNI = -1.D0                                                  
      DO 2 I=1,N-1                                                      
          EI = ET(I)                                                    
          SGNJ = -1.D0                                                  
      DO 3 J=1,N-1                                                      
      IF (I .NE. J) THEN                                                
          EJ = ET(J)                                                    
          DMA(I,J) = SGNI*SGNJ/(EI-EJ)                                  
      ELSE                                                              
          DMA(I,I) = -.5D0*EI/(1.D0-EI*EI)                              
      ENDIF                                                             
          SGNJ = -SGNJ                                                  
3     CONTINUE                                                          
          SGNI = -SGNI                                                  
2     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE DMLAGR(N,NM,A,ET,DMA)                                  
************************************************************************
*  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE        
*  LAGUERRE GAUSS-RADAU NODES                                           
*  N   = DIMENSION OF THE MATRIX                                        
*  NM  = ORDER OF THE MATRIX AS DECLARED IN THE MAIN DIMENSION STATEMENT
*  A   = PARAMETER >-1                                                  
*  ET  = VECTOR OF THE NODES, ET(I), I=0,N-1                            
*  DMA = DERIVATIVE MATRIX, DMA(I,J), I=0,N-1  J=0,N-1                  
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), DMA(0:NM,0:*)                                  
          DN = DFLOAT(N)                                                
          DMA(0,0) = (1.D0-DN)/(A+2.D0)                                 
      IF (N .EQ. 1) RETURN                                              
                                                                        
          A1 = A+1.D0                                                   
          X  = 0.D0                                                     
      CALL VALAPO(N,A,X,Y,DY,D2Y)                                       
          C = Y                                                         
      DO 1 J=1,N-1                                                      
          EJ = ET(J)                                                    
      CALL VALAPO(N,A,EJ,Y,DY,D2Y)                                      
          DMA(0,J) = -C/(A1*EJ*Y)                                       
          DMA(J,0) = A1*Y/(C*EJ)                                        
1     CONTINUE                                                          
                                                                        
      DO 2 I=1,N-1                                                      
          EI = ET(I)                                                    
      CALL VALAPO(N,A,EI,Y,DY,D2Y)                                      
      DO 3 J=1,N-1                                                      
      IF (I .NE. J) THEN                                                
          EJ = ET(J)                                                    
          DMA(I,J) = Y/(EI-EJ)                                          
      ELSE                                                              
          DMA(I,I) = .5D0*(EI-A)/EI                                     
      ENDIF                                                             
3     CONTINUE                                                          
2     CONTINUE                                                          
                                                                        
      DO 4 J=1,N-1                                                      
          EJ = ET(J)                                                    
      CALL VALAPO(N,A,EJ,Y,DY,D2Y)                                      
      DO 5 I=1,N-1                                                      
      IF (I .EQ. J) GOTO 5                                              
          DMA(I,J) = DMA(I,J)/Y                                         
5     CONTINUE                                                          
4     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE FCCHGA(N,QZ,CO)                                        
************************************************************************
*  COMPUTES USING FFT THE CHEBYSHEV FOURIER COEFFICIENTS OF A POLYNOMIAL
*  INDIVIDUATED BY THE VALUES ATTAINED AT THE CHEBYSHEV ZEROES          
*  N   = THE NUMBER OF ZEROES                                           
*  QZ  = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N           
*  CO  = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N-1         
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION QZ(1), CO(0:*)                                          
      IF(N .EQ. 0) RETURN                                               
                                                                        
          PH = 1.57079632679489661923D0                                 
          R2 = 1.41421356237309515D0                                    
          DN = DFLOAT(N)                                                
          CO(0) = QZ(1)/DN                                              
      IF(N .EQ. 1) RETURN                                               
                                                                        
          N2 = N/2                                                      
          C  = 2.D0/DSQRT(DN)                                           
          SN = DFLOAT(1+4*N2-2*N)                                       
          CO(N-N2-1) = QZ(N)                                            
      DO 1 I=1,N2                                                       
          CO(I-1) = QZ(2*I-1)                                           
          CO(N-I) = QZ(2*I)                                             
1     CONTINUE                                                          
                                                                        
      CALL C06EAF(CO,N,IFAIL)                                           
      IF (IFAIL .NE. 0) THEN                                            
      WRITE(filenum,*) 'IFAIL IS NOT ZERO IN SUBROUTINE C06EAF'               
      ENDIF                                                             
                                                                        
          CO(0) = .5D0*C*CO(0)                                          
      IF (2*N2 .EQ. N) THEN                                             
          CO(N2) = C*((-1.D0)**N2)*CO(N2)/R2                            
      ENDIF                                                             
      IF (N .EQ. 2) RETURN                                              
          SM = -1.D0                                                    
      DO 2 M=1,N-N2-1                                                   
          AR = PH*DFLOAT(M)/DN                                          
          CS = DCOS(AR)                                                 
          SI = DSIN(AR)                                                 
          V1 = C*SM*(CO(M)*CS+CO(N-M)*SI)                               
          V2 = C*SM*SN*(CO(M)*SI-CO(N-M)*CS)                            
          CO(M) = V1                                                    
          CO(N-M) = V2                                                  
          SM = -SM                                                      
2     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE FCCHGL(N,QN,CO)                                        
************************************************************************
*  COMPUTES USING FFT THE CHEBYSHEV FOURIER COEFFICIENTS OF A POLYNOMIAL
*  INDIVIDUATED BY ITS VALUES AT THE CHEBYSHEV GAUSS-LOBATTO NODES      
*  N   = THE DEGREE OF THE POLYNOMIAL                                   
*  QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N            
*  CO  = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N           
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION QN(0:*), CO(0:*)                                        
      IF(N .EQ. 0) RETURN                                               
                                                                        
          PI = 3.14159265358979323846D0                                 
          DN = DFLOAT(N)                                                
          N2 = N/2                                                      
          SN = DFLOAT(1+4*N2-2*N)                                       
          S1 = .5D0*(QN(0)+QN(N))                                       
          S2 = .5D0*(SN*QN(0)+QN(N))                                    
          CO(0) = S1                                                    
          CO(N) = S2                                                    
      IF(N .EQ. 1) RETURN                                               
                                                                        
          CO(0) = QN(0)                                                 
      IF (2*N2 .EQ. N) THEN                                             
          CO(N2) = QN(N)                                                
      ENDIF                                                             
      IF(N .EQ. 2) GOTO 2                                               
                                                                        
      DO 1 I=1,N-N2-1                                                   
          I2 = 2*I                                                      
          CO(I) = QN(I2)                                                
          CO(N-I) = QN(I2-1)-QN(I2+1)                                   
1     CONTINUE                                                          
                                                                        
2     CALL C06EBF(CO,N,IFAIL)                                           
      IF (IFAIL .NE. 0) THEN                                            
      WRITE(filenum,*) 'IFAIL IS NOT ZERO IN SUBROUTINE C06EBF'               
      ENDIF                                                             
                                                                        
          SJ = -1.D0                                                    
      DO 3 J=1,N-1                                                      
          S1 = S1+QN(J)                                                 
          S2 = S2+SJ*SN*QN(J)                                           
          SJ = -SJ                                                      
3     CONTINUE                                                          
          CO(0) = S1/DN                                                 
          CO(N) = S2/DN                                                 
          C  = .5D0/DSQRT(DN)                                           
      IF (2*N2 .EQ. N) THEN                                             
          CO(N2) = 2.D0*C*((-1.D0)**N2)*CO(N2)                          
      ENDIF                                                             
      IF(N .EQ. 2) RETURN                                               
                                                                        
          SM = -1.D0                                                    
      DO 4 M=1,N-N2-1                                                   
          AR = PI*DFLOAT(M)/DN                                          
          SI = .5D0/DSIN(AR)                                            
          V1 = CO(M)*(1.D0+SI)+CO(N-M)*(1.D0-SI)                        
          V2 = CO(M)*(1.D0-SI)+CO(N-M)*(1.D0+SI)                        
          CO(M) = C*SM*V1                                               
          CO(N-M) = C*SM*SN*V2                                          
          SM = -SM                                                      
4     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE FVCHGL(N,CO,QN)                                        
*******************************************************************     
*  COMPUTES USING FFT THE VALUES OF A POLYNOMIAL AT THE CHEBYSHEV       
*  GAUSS-LOBATTO NODES FROM ITS CHEBYSHEV FOURIER COEFFICIENTS          
*  N   = THE DEGREE OF THE POLYNOMIAL                                   
*  CO  = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N           
*  QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N            
*******************************************************************     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION QN(0:*), CO(0:*)                                        
      IF(N .EQ. 0) RETURN                                               
                                                                        
          PI = 3.14159265358979323846D0                                 
          DN = DFLOAT(N)                                                
          N2 = N/2                                                      
          SN = DFLOAT(1+4*N2-2*N)                                       
          S1 = CO(0)+SN*CO(N)                                           
          S2 = CO(0)+CO(N)                                              
          QN(0) = S1                                                    
          QN(N) = S2                                                    
      IF(N .EQ. 1) RETURN                                               
                                                                        
          QN(0) = 2.D0*CO(0)                                            
      IF (2*N2 .EQ. N) THEN                                             
          QN(N2) = 2.D0*CO(N)                                           
      ENDIF                                                             
      IF(N .EQ. 2) GOTO 2                                               
                                                                        
      DO 1 I=1,N-N2-1                                                   
          I2 = 2*I                                                      
          QN(N-I) = CO(I2+1)-CO(I2-1)                                   
          QN(I) = CO(I2)                                                
1     CONTINUE                                                          
      IF (2*N2 .NE. N) THEN                                             
          QN(N2+1) = 2.D0*CO(N)-CO(N-2)                                 
      ENDIF                                                             
                                                                        
2     CALL C06EBF(QN,N,IFAIL)                                           
      IF (IFAIL .NE. 0) THEN                                            
      WRITE(filenum,*) 'IFAIL IS NOT ZERO IN SUBROUTINE C06EBF'               
      ENDIF                                                             
                                                                        
          SJ = -1.D0                                                    
      DO 3 J=1,N-1                                                      
          S1 = S1+SJ*CO(J)                                              
          S2 = S2+CO(J)                                                 
          SJ = -SJ                                                      
3     CONTINUE                                                          
          QN(0) = S1                                                    
          QN(N) = S2                                                    
          C  = .25D0*DSQRT(DN)                                          
      IF (2*N2 .EQ. N) THEN                                             
          QN(N2) = 2.D0*C*QN(N2)                                        
      ENDIF                                                             
      IF(N .EQ. 2) RETURN                                               
                                                                        
      DO 4 M=1,N-N2-1                                                   
          AR = PI*DFLOAT(M)/DN                                          
          SI = .5D0/DSIN(AR)                                            
          V1 = QN(M)*(1.D0+SI)+QN(N-M)*(1.D0-SI)                        
          V2 = QN(M)*(1.D0-SI)+QN(N-M)*(1.D0+SI)                        
          QN(M) = C*V1                                                  
          QN(N-M) = C*V2                                                
4     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE FDCHGL(N,QN,CD,DQN)                                    
*********************************************************************   
*  COMPUTES USING FFT THE FOURIER COEFFICIENTS AND THE VALUES AT THE    
*  CHEBYSHEV GAUSS-LOBATTO NODES OF THE DERIVATIVE OF A POLYNOMIAL      
*  FROM THE VALUES ATTAINED BY THE POLYNOMIAL AT THE NODES              
*  N   = THE DEGREE OF THE POLYNOMIAL                                   
*  QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N            
*  CD  = FOURIER COEFFICIENTS OF THE DERIVATIVE, CD(I), I=0,N           
*  DQN = VALUES OF THE DERIVATIVE AT THE NODES, DQN(I), I=0,N           
*********************************************************************   
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION QN(0:*), CD(0:*), DQN(0:*)                              
      IF(N .EQ. 0) RETURN                                               
                                                                        
      CALL FCCHGL(N,QN,DQN)                                             
          CD(N) = 0.D0                                                  
          CD(0) = DQN(1)                                                
      IF(N .EQ. 1) GOTO 2                                               
                                                                        
          DN = DFLOAT(N)                                                
          CD(N-1) = 2.D0*DN*DQN(N)                                      
      DO 1 K=0,N-2                                                      
          KR = N-K-2                                                    
      IF(KR .NE. 0) THEN                                                
          DK = 2.D0*(DFLOAT(KR)+1.D0)                                   
          CD(KR) = CD(KR+2)+DK*DQN(KR+1)                                
      ELSE                                                              
          CD(0) = .5D0*CD(2)+DQN(1)                                     
      ENDIF                                                             
1     CONTINUE                                                          
                                                                        
2     CALL FVCHGL(N,CD,DQN)                                             
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE EDLEGL(N,ET,VN,FU,S1,S2,PA,SO,DSO,D2SO)                
************************************************************************
*  COMPUTES THE APPROXIMATE SOLUTION OF A DIRICHLET PROBLEM BY          
*  COLLOCATION AT THE LEGENDRE GAUSS-LOBATTO NODES                      
*  N   = DEGREE OF THE APPROXIMATING POLYNOMIAL                         
*  ET  = VECTOR OF THE NODES, ET(I), I=0,N                              
*  VN  = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N   
*  FU  = VALUES OF THE RIGHT-HAND SIDE AT THE NODES, FU(I), I=0,N       
*  S1  = DIRICHLET DATUM AT X=-1                                        
*  S2  = DIRICHLET DATUM AT X=1                                         
*  PA  = PARAMETER >0                                                   
*  SO  = VALUES OF THE APPROXIMATED SOLUTION AT THE NODES, SO(I), I=0,N 
*  DSO = DERIVATIVE OF THE SOLUTION AT THE NODES, DSO(I), I=0,N         
*  D2SO= SECOND DERIVATIVE OF THE SOLUTION AT THE NODES, D2SO(I), I=0,N 
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION ET(0:*), VN(0:*), FU(0:*), SO(0:*), DSO(0:*), D2SO(0:*) 
      IF (N .EQ. 0) RETURN                                              
                                                                        
          SO(0) = S1                                                    
          SO(N) = S2                                                    
          C1 = .5D0*(S2-S1)                                             
          DSO(0) = C1                                                   
          DSO(N) = C1                                                   
          D2SO(0) = 0.D0                                                
          D2SO(N) = 0.D0                                                
      IF (N .EQ. 1) RETURN                                              
                                                                        
      DO 1 J=1,N-1                                                      
          SO(J) = .5D0*(S1+S2)+C1*ET(J)                                 
1     CONTINUE                                                          
                                                                        
          TH = .58D0                                                    
          IM = 42+DLOG(1.D0+PA)                                         
      DO 2 IT=1,IM                                                      
      CALL DELEGL(N,ET,VN,SO,DSO)                                       
      CALL DELEGL(N,ET,VN,DSO,D2SO)                                     
      DO 3 I=1,N-1                                                      
          D2SO(I) = D2SO(I)+FU(I)-PA*SO(I)                              
3     CONTINUE                                                          
                                                                        
          D2SO(N) = 0.D0                                                
          C2 = 2.D0/(ET(1)-ET(0))                                       
          BE = C2/(ET(2)-ET(0))                                         
          DSO(1) = D2SO(1)                                              
          C3 = 0.D0                                                     
      IF(N .EQ. 2) GOTO 5                                               
                                                                        
      DO 4 I=1,N-2                                                      
          D1 = ET(I)-ET(I-1)                                            
          D2 = ET(I+1)-ET(I)                                            
          D3 = ET(I+2)-ET(I+1)                                          
          C4 = -BE*C3+2.D0/(D1*D2)+PA                                   
          D2SO(I) = C4                                                  
          BE = 2.D0/(C4*D2*(D2+D3))                                     
          DSO(I+1) = D2SO(I+1)+BE*DSO(I)                                
          C3 = 2.D0/(D2*(D1+D2))                                        
4     CONTINUE                                                          
                                                                        
5         D2SO(N-1) = -BE*C3+C2/(ET(2)-ET(1))+PA                        
          DSO(N) = 0.D0                                                 
      DO 6 I=1,N-1                                                      
          IR = N-I                                                      
          D1 = ET(IR)-ET(IR-1)                                          
          D2 = ET(IR+1)-ET(IR)                                          
          DSO(IR) = (DSO(IR)+2.D0*DSO(IR+1)/(D2*(D1+D2)))/D2SO(IR)      
6     CONTINUE                                                          
                                                                        
      DO 7 I=1,N-1                                                      
          SO(I) = SO(I)+TH*DSO(I)                                       
7     CONTINUE                                                          
                                                                        
2     CONTINUE                                                          
                                                                        
      CALL DELEGL(N,ET,VN,SO,DSO)                                       
      CALL DELEGL(N,ET,VN,DSO,D2SO)                                     
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
      SUBROUTINE DMLEGA(N,NM,CS,DZ,DGA)                                 
************************************************************************
*  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE        
*  LEGENDRE GAUSS NODES                                                 
*  N   = DIMENSION OF THE MATRIX                                        
*  NM  = ORDER OF THE MATRIX AS DECLARED IN THE MAIN DIMENSION STATEMENT
*  CS  = VECTOR OF THE ZEROES, CS(I), I=1,N                             
*  DZ  = LEGENDRE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N               
*  DGA = DERIVATIVE MATRIX, DGA(I,J), I=1,N  J=1,N                      
************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION CS(1), DZ(1), DGA(NM,NM)                                
          DGA(1,1) = 0.D0                                               
      IF (N .LE. 1) RETURN                                              
                                                                        
      DO 1 I=1,N                                                        
          VI = DZ(I)                                                    
          ZI = CS(I)                                                    
      DO 2 J=1,N                                                        
      IF (I .NE. J) THEN                                                
          VJ = DZ(J)                                                    
          ZJ = CS(J)                                                    
          DGA(I,J) = VI/(VJ*(ZI-ZJ))                                    
      ELSE                                                              
          DGA(I,I) = ZI/(1.D0-ZI*ZI)                                    
      ENDIF                                                             
2     CONTINUE                                                          
1     CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               