C++*********************************************************************
C
C RADAVI.F
C
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
C  RADAVI.F
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE RADAVI(LUN1,LUN2,NSAM1,NROW1,MAXDIM)

 

         COMMON BUF(1)

         IRMAX = MIN0(NSAM1/2,NROW1/2)
         CALL CRCSE2(LUN1,BUF,NSAM1,NROW1,
     &               BUF(NSAM1+1),BUF(NSAM1+IRMAX+1),IRMAX)
C         ASUM = 0
C         DO I=1,IRMAX
C            ASUM=ASUM+BUF(I+NSAM1)
C         ENDDO
C         AVA = ASUM/FLOAT(IRMAX)
C         PLACE OUTSIDE OF THE CIRCLE (IN THE CORNERS) LAST ELEMNT, NOT THE
C         AVERAGE
         AVA=BUF(IRMAX+NSAM1)
         DO I=1,NROW1
            KI = I-NROW1/2-1
            KI=KI*KI
            DO K=1,NSAM1
               KK  = K-NSAM1/2-1
               R   = SQRT(FLOAT(KI)+FLOAT(KK*KK))+1.0
               IR  = INT(R)
               IF (IR .LT. IRMAX) THEN
                  BUF(K) = BUF(IR+NSAM1)+(BUF(IR+1+NSAM1)-
     &                     BUF(IR+NSAM1))*(R-IR)
               ELSE
                  BUF(K) = AVA
               ENDIF
            ENDDO
            CALL WRTLIN(LUN2,BUF,NSAM1,I)
         ENDDO

         END
 
