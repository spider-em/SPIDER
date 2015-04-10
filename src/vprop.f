
C ++********************************************************************
C                                                                      *
C                                                                      *
C                                                                      *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@wadsworth.org                                        *
C=*                                                                    *
C=* SPIDER is free software; you can redistribute it and/or            *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* SPIDER is distributed in the hope that it will be useful,          *
C=* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C=* merchantability or fitness for a particular purpose.  See the GNU  *
C=* General Public License for more details.                           *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program. If not, see <http://www.gnu.org/licenses> *
C=*                                                                    *
C **********************************************************************
C                                                                      *
C VPROP                                                                     *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C
C Double precision version 10/21/96 PP
C======================================================================

      SUBROUTINE VPROP( NDIM, N, W, D, S, KODE ) 
                            
      DIMENSION  W(NDIM,N) , D(N) , S(N) 
                                       
      DATA  SEUIL / 1.0E-7 /
	DOUBLE PRECISION Q,H,T,U,V,P,B
                                                  
      KODE  = 0                                                                 
      CALL  TRIDI (NDIM,N,W,D,S)                                                
        IF (N .EQ. 1)        GO TO  140                                         
        DO   I = 2,N                                                          
      S(I-1)= S(I)                                                              
        ENDDO
      S(N)  = 0.0                                                               
        DO 90  K = 1,N                                                          
      M     = 0                                                                 
   20   DO   J = K,N                                                          
        IF (J .EQ. N)        GO TO   40                                         
      ABJ   = ABS(S(J))                                                         
      EPS   = SEUIL*(ABS(D(J)) + ABS(D(J+1)))                                   
        IF (ABJ .LE. EPS)    GO TO   40                                         
	ENDDO
   40 H     = D(K)                                                              
        IF (J .EQ. K)        GO TO   90                                         
        IF (M .EQ. 30)       GO TO  130                                         
      M     = M + 1                                                             
      Q     = (D(K+1) - H) / (2.0*S(K))                                         
      T     = DSQRT (Q*Q + 1.0)                                                  
      Q     = D(J) - H + S(K) / (Q + DSIGN(T,Q))                                 
      U     = 1.0                                                               
      V     = 1.0                                                               
      H     = 0.0                                                               
      JK    = J - K                                                             
        DO  IJK = 1,JK                                                        
      I     = J - IJK                                                           
      P     = U * S(I)                                                          
      B     = V * S(I)                                                          
        IF (ABS(P).LT.ABS(Q))GO TO   50                                         
      V     = Q / P                                                             
      T     = DSQRT (V*V + 1.0)                                                  
      S(I+1)= P * T                                                             
      U     = 1.0 / T                                                           
      V     = V * U                                                             
                             GO TO   60                                         
   50 U     = P / Q                                                             
      T     = DSQRT (U*U + 1.0)                                                  
      S(I+1)= Q * T                                                             
      V     = 1.0 / T                                                           
      U     = U * V                                                             
   60 Q     = D(I+1) - H                                                        
      T     = (D(I) - Q)*U + 2.0*V*B                                            
      H     = U * T                                                             
      D(I+1)= Q + H                                                             
      Q     = V*T - B                                                           
        DO   L = 1,N                                                          
      P     = W(L,I+1)                                                          
      W(L,I+1)= U*W(L,I) + V*P                                                  
      W(L,I)= V*W(L,I) - U*P                                                    
	ENDDO
	ENDDO
      D(K)  = D(K) - H                                                          
      S(K)  = Q                                                                 
      S(J)  = 0.0                                                               
                             GO TO   20                                         
   90   CONTINUE                                                                
        DO 120  IJ = 2,N                                                        
      I     = IJ - 1                                                            
      L     = I                                                                 
      H     = D(I)                                                              
        DO 100  M = IJ,N                                                        
        IF (D(M) .LE. H)     GO TO  100                                         
      L     = M                                                                 
      H     = D(M)                                                              
  100   CONTINUE                                                                
        IF (L .EQ. I)        GO TO  120                                         
      D(L)  = D(I)                                                              
      D(I)  = H                                                                 
        DO   M = 1,N                                                         
      H     = W(M,I)                                                            
      W(M,I)= W(M,L)                                                            
      W(M,L)= H                                                                 
	ENDDO
  120   CONTINUE                                                                
                             GO TO  140                                         
  130 KODE  = K                                                                 
  140     RETURN                                                                
          END                                                                   

