
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
C                                                                      *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C Double precision version  10/21/96  PP
C
C***********************************************************************

      SUBROUTINE   TRIDI ( NDIM, N, W, D, S )                                   

      DIMENSION  D(N) , S(N) , W(NDIM,N)   
      DOUBLE PRECISION B,C,P,Q,PB
                                     
        IF (N .EQ. 1)        GO TO  130                                         
        DO   I2 = 2,N                                                        
      B     = 0.0                                                               
      C     = 0.0                                                               
      I     = N - I2 + 2                                                        
      K     = I - 1                                                             
        IF (K .LT. 2)        GO TO   20        
        DO   L= 1,K                                                           
      C     = C + ABS(W(I,L))                                                   
        ENDDO
        IF (C .NE. 0.0)      GO TO   30                                         
   20 S(I)  = W(I,K)                                                            
                             GO TO  110                                         
   30   DO   L= 1,K                                                           
      W(I,L)= W(I,L) / C                                                        
      B     = B + W(I,L)*W(I,L)                                                 
	ENDDO
      P     = W(I,K)                                                            
      Q     = -DSIGN (DSQRT(B) , P)                                               
      S(I)  = C * Q                                                             
      B     = B - P*Q                                                           
      W(I,K)= P - Q                                                             
      P     = 0.0                                                               
        DO    M = 1,K                                                         
      W(M,I)= W(I,M) / (B*C)                                                    
      Q     = 0.0                                                               
        DO    L = 1,M                                                         
      Q     = Q + W(M,L)*W(I,L)                                                 
        ENDDO
      M1    = M + 1                                                             
        IF (K .LT. M1)       GO TO   70                                         
        DO    L = M1,K                                                        
      Q     = Q + W(L,M)*W(I,L)                                                 
        ENDDO
   70 S(M)  = Q / B                                                             
      P     = P + S(M)*W(I,M)                                                   
	ENDDO
      PB    = P / (B+B)                                                         
        DO   M = 1,K                                                         
      P     = W(I,M)                                                            
      Q     = S(M) - PB*P                                                       
      S(M)  = Q                                                                 
        DO    L = 1,M                                                         
      W(M,L)= W(M,L) - P*S(L) - Q*W(I,L)                                        
	ENDDO
	ENDDO
        DO   L = 1,K                                                         
      W(I,L)= C * W(I,L)                                                        
        ENDDO
  110 D(I)  = B                                                                 
	ENDDO
  130 S(1)  = 0.0                                                               
      D(1)  = 0.0                                                               
        DO 180  I = 1,N                                                         
      K     = I - 1                                                             
        IF (D(I) .EQ. 0.0)   GO TO  160                                         
        DO   M = 1,K                                                         
      Q     = 0.0                                                               
        DO   L = 1,K                                                         
      Q     = Q + W(I,L)*W(L,M)                                                 
        ENDDO
        DO   L = 1,K                                                         
      W(L,M)= W(L,M) - Q*W(L,I)                                                 
	ENDDO
	ENDDO
  160 D(I)  = W(I,I)                                                            
      W(I,I)= 1.0                                                               
        IF (K .LT. 1)        GO TO  180                                         
        DO   M = 1,K                                                         
      W(I,M)= 0.0                                                               
      W(M,I)= 0.0                                                               
	ENDDO
  180   CONTINUE                                                                
          END

