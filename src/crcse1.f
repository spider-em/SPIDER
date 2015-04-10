C++*********************************************************************
C
C CRCSE1.F
C              RITHALF ADDED                     NOV 2013 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2013  Health Research Inc.,                         *
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
C  CRCSE1(LUN,LUN2,NX,NY,IR,SEC,HALF)
C
C  PURPOSE: ROTATIONAL AVERAGING INTO A SINGLE LINE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE CRCSE1(LUN,LUN2,NX,NY,IR,SEC,HALF)

         IMPLICIT NONE
         INCLUDE 'CMBLOCK.INC'

         INTEGER           :: LUN,LUN2,NX,NY,IR
         REAL              :: SEC(IR)
         CHARACTER (LEN=1) :: HALF

         REAL              :: BUF(NX),SNO(IR)
         REAL              :: R,XD
         INTEGER           :: J,KJ,I,KI,L,IGO,IEND
         
         SEC = 0.0  ! ARRAY ZERO
         SNO = 0.0  ! ARRAY ZERO
 
         IGO  = 1
         IEND = NX

         IF (HALF == 'L') THEN
           WRITE(NOUT,*) ' USING LEFT-HALF OF IMAGE ONLY'
           IGO  = 1
           IEND = NX / 2

         ELSEIF (HALF == 'R') THEN
           WRITE(NOUT,*) ' USING RIGHT-HALF OF IMAGE ONLY'
           IGO  = NX / 2 + 1
           IEND = NX
         ENDIF

         DO J=1,NY
            KJ = J-NY/2-1

            IF (IABS(KJ) <= (IR-1) )  THEN
               CALL REDLIN(LUN,BUF,NX,J)

               DO I=IGO,IEND
                  KI = I - NX / 2 - 1

C                 FIND RADIUS
                  R  = SQRT(FLOAT(KJ*KJ) + 
     &                      FLOAT(KI*KI)) + 1.0
                  L = R
                  IF (L <= IR-1) THEN
C                    INSIDE CIRCLE
                     XD       = R - L
                     SEC(L)   = SEC(L)   + BUF(I) * (1.0-XD)
                     SEC(L+1) = SEC(L+1) + BUF(I) * XD
                     SNO(L)   = SNO(L)   + 1.0 - XD
                     SNO(L+1) = SNO(L+1) + XD
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

         DO I=1,IR
            SEC(I) = SEC(I) / MAX(1.0,SNO(I))
         ENDDO

         CALL WRTLIN(LUN2,SEC,IR,1)

         END
