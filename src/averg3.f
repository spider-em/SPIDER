C++*******************************************************************
C
C AVERG3.F
C         X PASSED                              AUG 12 ARDEAN LEITH
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
C  AVERG3(LUN1,XAVG,J,K,L,NX,NY,NZ,DS,DR,DL,MOVWAY,LREAD,X)
C      
C  XAVG: RECURSIVE LOCAL AVERAGE OF 3D PICTURE VOLUME
C  X( ): DL PICTURE SLICES STORED
C  X[J,K,L]<>X(MOD(L-1,DL)*NY+J*NX+K)
C  J,K,L      PICTURE COORD
C  NX,NY,NZ:  PICTURE SIZE
C  DS,DR,DL:  LOCAL BOX  SIZE
C  MOVWAY= 1: READ IN DL SLICES AND COMPUTE 1ST AVERAGE VALUE FROM SCRATCH
C             IF J,K,L IS IN A BORDER AREA COMPUTE XAVG FOR NEAREST
C             ACTIVE POINT
C        = 2: MOVE RIGHT ONE PIXEL
C        = 3: MOVE LEFT  ONE PIXEL
C        = 4: MOVE DOWN  RIGHT ONE PIXEL
C        = 5: MOVE DOWN  LEFT  ONE PIXEL
C        = 6: MOVE UP,   READ IN NEXT SLICE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*******************************************************************

        SUBROUTINE AVERG3(LUN1,XAVG, J,K,L, NX,NY,NZ,
     &                    DS,DR,DL, MOVWAY,LREAD,X)

        INTEGER :: LUN1
        REAL    :: XAVG
        INTEGER :: J,K,L, NX,NY,NZ,DS,DR,DL,MOVWAY,LREAD

        INTEGER :: HDS,HDR,HDL
        INTEGER :: EDGFST,EDGLST,EDGLFT,EDGRT,EDGDW,EDGUP

        REAL    :: X(*)
        SAVE             ! WHY??al

        IF (MOVWAY == 1) THEN        
C          INITIALIZE PARM AND XAVG ---------------------- 1
           HDR    = DR / 2
           HDS    = DS / 2
           HDL    = DL / 2
           EDGFST = HDR + 1
           EDGLST = NY  - HDR
           EDGLFT = HDS + 1
           EDGRT  = NX  - HDS
           EDGDW  = HDL + 1
           EDGUP  = NZ  - HDL

C          CHECK IF IN BORDER
           XAVG = 0.
           JT   = J
           KT   = K
           KL   = L
           IF (J .LT. EDGFST) JT = EDGFST
           IF (J .GT. EDGLST) JT = EDGLST
           IF (K .LT. EDGLFT) KT = EDGLFT
           IF (K .GT. EDGRT)  KT = EDGRT
           IF (L .LT. EDGDW)  KL = EDGDW
           IF (L .GT. EDGUP)  KL = EDGUP

C          READ IN DL SLICES AND COMPUTE AVERAGE FROM SCRATCH
           DO LL=1,DL
              LSTART = (LL-1)* NX * NY

              DO JJ=1,NY
                 NL = MOD(LL-1,DL)  * NY + JJ
                 NI = LSTART + (JJ-1) * NX
                 CALL REDLIN(LUN1,X(NI+1),NX,NL)

                 LREAD = LREAD + 1
               ENDDO
           ENDDO

           DO LL=KL-HDL,KL+HDL
              LSTART = MOD(LL-1,DL) * NX * NY

              DO  JJ=JT-HDR,JT+HDR
                 NI = LSTART + (JJ-1) * NX

                 DO  KK=KT-HDS,KT+HDS
                    XAVG = XAVG + X(NI+KK)
                 ENDDO
              ENDDO
           ENDDO

           XAVG = XAVG / FLOAT(DR*DS*DL)
           RETURN

        ELSEIF (MOVWAY == 2) THEN 
C          MOVE RT ONE PIXEL, UPDATE XAVG ----------------------2
           JT = J
           LT = L
           IF (K .LE. EDGLFT .OR. K.GT. EDGRT) RETURN
           IF (J .LT. EDGFST) JT = EDGFST
           IF (J .GT. EDGLST) JT = EDGLST
           IF (L .LT. EDGDW)  LT = EDGDW
           IF (L .GT. EDGUP)  LT = EDGUP
           CORECT = 0.0
           DO LL=LT-HDL,LT+HDL
              LSTART = MOD(LL-1,DL)*NY*NX

              DO  JJ=JT-HDR,JT+HDR
                 NI     = LSTART + (JJ-1) * NX
                 CORECT = CORECT + X (NI+K+HDS) - X(NI+K-HDS-1)
              ENDDO
           ENDDO
           XAVG = CORECT / FLOAT(DR*DS*DL) + XAVG
           RETURN

        ELSEIF (MOVWAY == 3) THEN
C          MOVE LEFT ONE PIXEL UPDATE XAVG ---------------------- 3
           JT=J
           LT=L
           IF (K .LT. EDGLFT .OR. K .GE. EDGRT) RETURN
           IF (J .LT. EDGFST) JT = EDGFST
           IF (J .GT. EDGLST) JT = EDGLST
           IF (L .LT. EDGDW)  LT = EDGDW
           IF (L .GT. EDGUP)  LT = EDGUP
           CORECT = 0.0
           DO    LL=LT-HDL,LT+HDL
              LSTART = MOD(LL-1,DL)*NY*NX

              DO  JJ=JT-HDR,JT+HDR
                 NI     = LSTART+(JJ-1)*NX
                 CORECT = CORECT + X(NI+K-HDS) - X(NI+K+HDS+1)
              ENDDO
           ENDDO
           XAVG = CORECT / FLOAT(DR*DS*DL) + XAVG
           RETURN

        ELSEIF (MOVWAY == 4) THEN
C          J INCREMENTED ONE ----------------------------------- 4
           IF(J.LE.EDGFST .OR. J.GT.EDGLST) RETURN
           LT = L
           KT = K
           IF(K .LT. EDGLFT) KT = EDGLFT
           IF(K .GT. EDGRT)  KT = EDGRT
           IF(L .LT. EDGDW)  LT = EDGDW
           IF(L .GT. EDGUP)  LT = EDGUP
           CORECT = 0.0
           DO LL=LT-HDL,LT+HDL
              LSTART = MOD(LL-1,DL)*NY*NX

              DO  KK=KT-HDS,KT+HDS
                 NI     = LSTART + KK+(J-HDR-2) * NX
                 NO     = LSTART + KK+(J-1+HDR) * NX
                 CORECT = CORECT + X(NO) - X(NI)
              ENDDO
           ENDDO
           XAVG = XAVG + CORECT / FLOAT(DS*DR*DL)
           RETURN

        ELSEIF (MOVWAY == 5) THEN
C          J DECREMENTED ONE ----------------------------------- 5

           IF(J .LT. EDGFST .OR. J .GE. EDGLST) RETURN
           LT = L
           KT = K
           IF (K .LT. EDGLFT) KT = EDGLFT
           IF (K .GT. EDGRT)  KT = EDGRT
           IF (L .LT. EDGDW)  LT = EDGDW
           IF (L .GT. EDGUP)  LT = EDGUP
           CORECT = 0.0
           DO    LL=LT-HDL,LT+HDL
              LSTART = MOD(LL-1,DL)*NY*NX

              DO  KK=KT-HDS,KT+HDS
                 NI     = LSTART + KK+(J+HDR)   * NX
                 NO     = LSTART + KK+(J-HDR-1) * NX
                 CORECT = CORECT + X(NO) - X(NI)
              ENDDO
           ENDDO
           XAVG = XAVG + CORECT / FLOAT(DS*DR*DL)
           RETURN

        ELSEIF (MOVWAY == 6) THEN
C          L INCREMENTED ONE ----------------------------------- 6

           IF (L .LE. EDGDW .OR. L .GT. EDGUP) RETURN
           JT = J
           IF (J .LT. EDGFST) JT = EDGFST
           IF (J .GT. EDGLST) JT = EDGLST
           KT = K
           IF (K .LT. EDGLFT) KT = EDGLFT
           IF (K .GT. EDGRT)  KT = EDGRT

           CORECT = 0.0
           LOLD   = L - HDL -1
           LSTART = MOD(LOLD-1,DL)* NX * NY

           DO JJ=JT-HDR, JT+HDR
              NI = (JJ-1)*NX + LSTART

              DO  KK = KT-HDS, KT+HDS
                 CORECT = CORECT - X(NI+KK)
              ENDDO
           ENDDO

           DO JJ=1,NY
              NL    = L*NY + JJ
              NI    = LSTART + (JJ-1)*NX

              LREAD = LREAD + 1

              CALL  REDLIN(LUN1,X(NI+1),NX,lread)

              !it = mod(nl,ny) 
              !if (it == 1) write(6,*) ' read line1:',nl

           ENDDO

           DO JJ=JT-HDR, JT+HDR
              NI = (JJ-1) * NX + LSTART

              DO  KK=KT-HDS, KT+HDS
                 CORECT = CORECT + X(NI+KK)
              ENDDO
           ENDDO

           XAVG = XAVG + CORECT / FLOAT(DS*DR*DL)

           RETURN
        ENDIF

        END
