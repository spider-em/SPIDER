C++*******************************************************************
C
C  AVERG.F
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
C  AVERG(LUN1,XAVG,J,K,NX,NY,DS,DR,MOVWAY,X)
C      
C  XAVG: RECURSIVE LOCAL AVERAGE OF PICTURE AREA
C  X( ): DR PICTURE LINES STORED  X[J,K]<>X((MOD(J-1,DR))*NX+K)
C  J,K PICTURE COORD
C  NX,NY: PICTURE SIZE
C  DS,DR: LOCAL BOX  SIZE
C  MOVWAY= 1: READ IN DR LINE AND COMPUTE 1ST AVERAGE VALUE FROM SCRATCH
C             IF J,K IS IN A BORDER AREA COMPUTE XAVG FOR NEAREST
C             ACTIVE POINT
C        = 2: MOVE RIGHT ONE PIXEL
C        = 3: MOVE LEFT  ONE PIXEL
C        = 4: MOVE DOWN  ONE PIXEL, READ IN NEXT LINE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*******************************************************************

        SUBROUTINE AVERG(LUN1,XAVG,J,K,NX,NY,DS,DR,MOVWAY,X)

        INTEGER :: DS,DR, HDS,HDR, EDGFST,EDGLST,EDGLFT,EDGRT

        REAL    :: X(*)

        SAVE       ! WHY?? al


        IF (MOVWAY .EQ. 1) THEN       
C          INITIALIZE PARM AND XAVG
           HDR    = DR  / 2
           HDS    = DS  / 2
           EDGFST = HDR + 1
           EDGLST = NY  - HDR
           EDGLFT = HDS + 1
           EDGRT  = NX  - HDS

C          CHECK IF IN BORDER
           XAVG = 0.
           JT   = J
           KT   = K
           IF (J .LT. EDGFST) JT = EDGFST
           IF (J .GT. EDGLST) JT = EDGLST
           IF (K .LT. EDGLFT) KT = EDGLFT
           IF (K .GT. EDGRT)  KT = EDGRT

C          READ IN DR LINES AND COMPUTE AVERAGE FROM SCRATCH
           DO  JJ=JT-HDR,JT+HDR
              LSTART = (MOD(JJ-1,DR)) * NX

              CALL REDLIN(LUN1,X(LSTART+1),NX,JJ)

              DO  KK=KT-HDS, KT+HDS
                 XAVG = XAVG + X(LSTART+KK)
              ENDDO
           ENDDO

           XAVG = XAVG / FLOAT(DR*DS)

        ELSEIF (MOVWAY .EQ. 2) THEN
C          MOVE RT ONE PIXEL, UPDATE XAVG
           JT = J
           IF (K .EQ. EDGLFT .OR. K.LT.EDGLFT .OR. K.GT.EDGRT)RETURN
           IF (J .LT. EDGFST) JT = EDGFST
           IF (J .GT. EDGLST) JT = EDGLST
           CORECT = 0.

           DO  JJ=JT-HDR, JT+HDR
              LSTART = (MOD(JJ-1,DR))*NX
              CORECT = CORECT + X(LSTART+K+HDS) - X(LSTART+K-HDS-1)
           ENDDO

           XAVG = CORECT/FLOAT(DR*DS)+XAVG

        ELSEIF (MOVWAY .EQ. 3) THEN
C          MOVE LEFT ONE PIXEL UPDATE XAVG
           JT = J
           IF (K .LT. EDGLFT .OR. K.GT.EDGRT .OR. K.EQ.EDGRT)RETURN
           IF (J .LT. EDGFST) JT = EDGFST
           IF (J .GT. EDGLST) JT = EDGLST

           CORECT = 0.

           DO  JJ=JT-HDR,JT+HDR
              LSTART = (MOD(JJ-1,DR))*NX
              CORECT = CORECT + X(LSTART+K-HDS) - X(LSTART+K+HDS+1)
           ENDDO

           XAVG = CORECT / FLOAT(DR*DS) + XAVG

        ELSEIF (MOVWAY .EQ. 4) THEN
C          J INCREMENTED ONE
           IF (J.LE.EDGFST .OR. J.GT.EDGLST)RETURN
           KT = K
           IF (K .LT. EDGLFT) KT = EDGLFT
           IF (K .GT. EDGRT)  KT = EDGRT
           CORECT = 0.
           JOLD   = J-HDR-1
           LSTART = (MOD(JOLD-1,DR))*NX

           DO  KK=KT-HDS, KT+HDS
              CORECT = CORECT - X(LSTART+KK)
           ENDDO

           CALL REDLIN(LUN1,X(LSTART+1),NX,J+HDR)

           DO  KK=KT-HDS,KT+HDS
              CORECT = CORECT + X(LSTART+KK)
           ENDDO

           XAVG = XAVG + CORECT / FLOAT(DS*DR)

        ENDIF

        END
