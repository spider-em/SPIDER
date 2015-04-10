C++*********************************************************************
C
C CRCSE3.FOR
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
C  CRCSE3(LUN,LUN2,NX,NY,IR)
C
C  PURPOSE: ROTATIONAL AVERAGING INTO A SINGLE LINE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE CRCSE3(LUN,LUN2,NX,NY,NZ,IR,SEC)

         IMPLICIT NONE

         INTEGER    :: LUN,LUN2,NX,NY,NZ,IR
         REAL       :: SEC(IR)

         REAL       :: BUF(NX),SNO(IR)
         REAL       :: R,xd
         INTEGER    :: K,KK,J,KJ,NR1,KI,L,I
     
         SEC = 0.0  ! ARRAY ZERO
         SNO = 0.0  ! ARRAY ZERO
   
         DO K=1,NZ
            KK = K-NZ/2-1

            IF (IABS(KK) .LE. IR-1) THEN
               DO J=1,NY
                  KJ  = J-NY/2-1
                  NR1 = J+(K-1)*NY

                  IF (IABS(KJ) .LE. IR-1) THEN
                     CALL REDLIN(LUN,BUF,NX,NR1)

                     DO I=1,NX
                        KI = I - NX / 2 - 1
                        R  = SQRT(FLOAT(KJ*KJ)  + 
     &                            FLOAT(KI*KI)  +
     &                            FLOAT(KK*KK)) + 1.0
                        L  = R
                        IF (L .LE. IR-1) THEN
                           XD       = R-L
                           SEC(L)   = SEC(L)+BUF(I)*(1.0-XD)
                           SEC(L+1) = SEC(L+1)+BUF(I)*XD
                           SNO(L)   = SNO(L)+1.0-XD
                           SNO(L+1) = SNO(L+1)+XD
                        ENDIF
                     ENDDO             
                  ENDIF
               ENDDO
            ENDIF
         ENDDO        

         DO I=1,IR
            SEC(I) = SEC(I) / AMAX1(1.0,SNO(I))
         ENDDO

         CALL WRTLIN(LUN2,SEC,IR,1)

         END
