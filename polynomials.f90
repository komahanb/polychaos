      SUBROUTINE LEGEN(N,X,Y,DY,D2Y)                                   
**************************************************************
*     COMPUTES THE VALUE OF THE LEGENDRE POLYNOMIAL OF DEGREE N           
*     AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT               
*     N  = DEGREE OF THE POLYNOMIAL                                       
*     X  = POINT IN WHICH THE COMPUTATION IS PERFORMED                    
*     Y  = VALUE OF THE POLYNOMIAL IN X                                   
*     DY = VALUE OF THE FIRST DERIVATIVE IN X                             
*     D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
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
 1    CONTINUE                                                          
      
      RETURN                                                            
      END                                                               
C     
      SUBROUTINE HERM(N,X,Y,DY,D2Y)                                   
*************************************************************
*     COMPUTES THE VALUE OF THE HERMITE POLYNOMIAL OF DEGREE N            
*     AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT               
*     N  = DEGREE OF THE POLYNOMIAL                                       
*     X  = POINT IN WHICH THE COMPUTATION IS PERFORMED                    
*     Y  = VALUE OF THE POLYNOMIAL IN X                                   
*     DY = VALUE OF THE FIRST DERIVATIVE IN X                             
*     D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
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
 1    CONTINUE                                                          
      DN  = 2.D0*DFLOAT(N)                                           
      DNN = 2.D0*DFLOAT(N-1)                                         
      DY  = DN*YP                                                    
      D2Y = DN*DNN*YPM                                               
      
      RETURN                                                            
      END                                                               
C     
      SUBROUTINE VALASF(N,A,X,Y,DY)                                     
*********************************************************************
*     COMPUTES THE VALUES OF THE SCALED LAGUERRE FUNCTION OF DEGREE N     
*     AND ITS FIRST DERIVATIVE AT A GIVEN POINT                           
*     N  = DEGREE                                                         
*     A  = PARAMETER >-1                                                  
*     X  = POINT (NON NEGATIVE) IN WHICH THE COMPUTATION IS PERFORMED     
*     Y  = VALUE OF THE FUNCTION IN X                                     
*     DY = VALUE OF THE FIRST DERIVATIVE IN X                             
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
 1    CONTINUE                                                          
      
      RETURN                                                            
      END                                                               
C     
