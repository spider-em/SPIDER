
C **********************************************************************
C                                                                      *
C  RPRQ.F                                                              *
C              COMMON PAR REMOVED                DEC 2010 ARDEAN LEITH *
C              RI REMOVED                        DEC 2010 ARDEAN LEITH *
C              PHI.. REMOVED                     FEB 2012 ARDEAN LEITH *
C                                                                      *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C  PURPOSE:         BACKPROJECTS B INTO CUBE. BILINEAR INTERPOLATION   *
C                                                                      *
C  PARAMETERS:      N       DIMENSIONS OF B                  SENT      *
C                   B       IMAGE                            SENT      *
C                   CUBE    VOLUME                           SENT/RET. *
C                   IPCUBE  SCAN LINE INDICES                SENT      *
C                   NN      2ND DIMENSION OF IPCUBE          SENT      *
C                   DM      ROTATION  MATRIX                 SENT      *
C                   LDP                                      SENT      *
C                   LDPNM                                    SENT      *
C                   IRTFLG                                   RET.      *
C                                                                      *
C  IPCUBE: 1 - BEGINNING VOXEL ON THIS LINE                            *
C          2 - ENDING VOXEL ON THIS LINE                               *
C          3 - IX     BEGINNING VOXEL COORDINATES                      *
C          4 - IY                                                      *
C          5 - IZ                                                      *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE RPRQ(N,B,CUBE,IPCUBE,NN,
     &                  DM,LDP,LDPNM,IRTFLG)

        INCLUDE 'CMBLOCK.INC'

        INTEGER           :: N
        REAL              :: B(N,N),CUBE(*)
        INTEGER           :: IPCUBE(5,NN)
        INTEGER           :: NN
        REAL              :: PHI,THETA,PSI
        REAL              :: DM(9)
        INTEGER           :: LDP,LDPNM,IRTFLG

        DOUBLE PRECISION  :: CPHI,SPHI,CTHE,STHE,CPSI,SPSI
        LOGICAL           :: ALLOK
        INTEGER           :: IBAD,JBAD
     
        ALLOK = .TRUE.

        DM1   = DM(1)
        DM4   = DM(4)

c$omp parallel do private(i,j,xb,yb,xbb,ybb,iqx,iqy,dipy),
c$omp&                    shared(allok,ibad,jbad)
        DO I=1,NN
           XB = (IPCUBE(3,I)-LDP) * DM(1) + (IPCUBE(4,I)-LDP)*DM(2) +
     &          (IPCUBE(5,I)-LDP) * DM(3)
           YB = (IPCUBE(3,I)-LDP) * DM(4) + (IPCUBE(4,I)-LDP)*DM(5) +
     &          (IPCUBE(5,I)-LDP) * DM(6)

           DO J=IPCUBE(1,I),IPCUBE(2,I)

              XBB = (J - IPCUBE(1,I)) * DM1 + XB
              IQX = IFIX(XBB + FLOAT(LDPNM))

              YBB = (J - IPCUBE(1,I)) * DM4 + YB
              IQY = IFIX(YBB + FLOAT(LDPNM))

              IF (IQX < 1 .OR. IQX >= N .OR.
     &            IQY < 1 .OR. IQY >= N) THEN
                 ALLOK = .FALSE.
                 IBAD  = IQX
                 JBAD  = IQY
              ELSE
                 DIPY = YBB + LDPNM - IQY

C                EVEN FASTER VERSION
                 CUBE(J) = CUBE(J) + B(IQX,IQY) + 
     &              DIPY * (B(IQX,IQY+1)   - B(IQX,IQY)) +
     &              (XBB + LDPNM - IQX) * (B(IQX+1,IQY) - B(IQX,IQY) + 
     &              DIPY * (B(IQX+1,IQY+1) - B(IQX+1,IQY) -
     &                      B(IQX,IQY+1)   + B(IQX,IQY)))
              ENDIF
           ENDDO
        ENDDO

        IF (ALLOK) THEN
           IRTFLG = 0
        ELSE
           IF (IBAD < 1 .OR. IBAD >= N)  THEN
              CALL ERRT(102,'Outside image, reduce the radius',IBAD)
           ELSE
              CALL ERRT(102,'Outside image, reduce the radius',JBAD)
           ENDIF
           IRTFLG = 1
        ENDIF

        END

C               FASTER VERSION
C               CUBE(J) = CUBE(J)
C    &          +B(IQX,IQY)+DIPY*(B(IQX,IQY+1)-B(IQX,IQY))
C    &          +DIPX*(B(IQX+1,IQY)-B(IQX,IQY)
C    &          +DIPY*(B(IQX+1,IQY+1)-B(IQX+1,IQY)
C    &          -B(IQX,IQY+1)+B(IQX,IQY)))
C
C               ORIGINAL VERSION
C                CUBE(J) = CUBE(J)
C     &                 +(1.0-DIPX)*(1.0-DIPY)*B(MAP(IQX,IQY))
C     &                 +     DIPX *(1.0-DIPY)*B(MAP(IQX+1,IQY))
C     &                 +(1.0-DIPX)*     DIPY *B(MAP(IQX,IQY+1))
C     &                 +     DIPX *     DIPY *B(MAP(IQX+1,IQY+1))

