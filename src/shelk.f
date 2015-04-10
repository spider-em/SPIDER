C **********************************************************************
C SHELK.F
C                                                             07.06.80 *        
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
C
C  PURPOSE:  SORT VECTOR X() IN ASCENDING ORDER                                  *        
C            WARNING: ORIGINAL VECTOR(X) GETS CLOBBERED UPON RETURN         
C            BUT THE ORIGINAL LOCATIONS ARE SAVED IN VECTOR KX()                 *        
C            REFERENCES:
C            (1) J.BOOTHROYD  SHELLSORT ALGORITHM.201 
C                COMM.ACM  VOL.6  (1963),NO.8,PP.445   
C            (2) D.A.SHELL  A HIGH-SPEED SORTING PROCEDURE         
C                COMM.ACM  VOL.2(1959),PP.30-32                                       *        
C                                                                    
C*---------------------------------------------------------------------*        

      SUBROUTINE SHELK(N, X, KX )
                                            
C     I DO NOT KNOW IF SAVE IS NEEDED FEB 99 al
      SAVE

      DIMENSION  X(N) , KX(N) 
                                                  
      DO   J = 1,N                                                          
         KX(J) = J                                                                 
      ENDDO

      I     = 1                                                                 
   20 I     = I + I                                                             
      IF (I .LE. N) GO TO  20                                          
      M     = I - 1                                                             
   30 M     = M / 2                                                             
      IF (M .EQ. 0) GO TO  70                                          
      K     = N - M  
                                                           
      DO 60 J = 1,K                                                          
         JM = J + M                                                             
   40    JM = JM - M                                                            
         IF (JM .LE. 0) GO TO  60                                          
         L = JM + M                                                            
         IF (X(L) .GE. X(JM)) GO TO  60                                          
         PIV   = X(JM)                                                             
         X(JM) = X(L)                                                              
         X(L)  = PIV                                                               
         KPIV  = KX(JM)                                                            
         KX(JM)= KX(L)                                                             
         KX(L) = KPIV                                                              
         GO TO  40 
                                         
   60   CONTINUE                                                                
        GO TO  30 
                                         
   70   RETURN                                                                  
        END   
