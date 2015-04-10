
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

C*--------------------------------------------------------------------*         
C*                                                                    *         
C*      EXHAUSTIVE RANDOM DRAWING OF K INDICES BETWEEN 1 AND N.       *         
C*      STORED IN IDK(*), IDN(N) IS A WORKING ARRAY.                  *         
C*      CALL ... SEN3A. REPLACED BY RAN FUNCTION WITH TIME-DEPENDENT  *         
C*      SEED J.F. 8/20/86. OPTION FOR SEED SPECIFICATION FROM PREVIOUS*
C*	RUN INSTALLED J.F. 11/26/86				      *         
C==06.07.80                                                           *         
C*--------------------------------------------------------------------*         

      SUBROUTINE RETIR(N,IDN,K,IDK)                                       

      DIMENSION  IDN(N) , IDK(K)

      DO J = 1,N                                                          
         IDN(J) = J                                                                 
      ENDDO

      KFIN  = K                                                                 
      IF (K .GT. N) KFIN = N 
                                
      DO  L = 1,KFIN                                                       
          X = N - L + 1                                                         
C         I = X*RAN(IST) + 1.0
	
          CALL RANDOM_NUMBER(VALUE)
          I = X * VALUE + 1.0
                                              
          IDK(L) = IDN(I)                                                            
          IF (I .NE. N) THEN 
             N1 = N - L                                                             
             IF (N1 .GT. 0) THEN                                
                DO   J = I,N1                                                         
                   IDN(J)= IDN(J+1)                                                          
                ENDDO
             ENDIF
          ENDIF
       ENDDO

       RETURN                                                                 
       END 
          

  
