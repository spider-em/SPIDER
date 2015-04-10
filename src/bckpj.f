C ++********************************************************************
C                                                                      *
C   BCKPJ.F    SPEEDED UP                       FEB 2000 ARDEAN LEITH  *
C              COMMON PAR REMOVED               DEC 2010 ARDEAN LEITH  *
C              FBS                              OCT 2011 G. KISHCHENKO *
C              FBS2                             DEC 2011 ARDEAN LEITH  *
C              RENAMED FROM: BCKCQ              JAN 2012 ARDEAN LEITH  *                                                                   *
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
C  BCKPJ_LIN(CUBE,NNN,DM,B,N,IPCUBE,NN,LDP,LDPNM)                        *
C  BCKPJ_FBS(CUBE,NNN,DM,PROJ,XDER,YDER,XYDER,
C            NXLD,NX,IPCUBE,NN, LDP,LDPNM)
C                                                                      *
C  PURPOSE:    BACKPROJECTS B INTO CUBE.                               *
C                                                                      *
C  PARAMETERS: CUBE                                          SENT/RET. *
C              NNN      DIMENSIONS OF CUBE ARRAY             SENT      *
C              DM       TRANSFORM MATRIX                     SENT      *
C              B        IMAGE (AND IMAGE LOC +1 DIFFERENCE   SENT      *
C              N        DIMENSIONS OF B                      SENT      *
C              IPCUBE   SCAN LINE INDICES                    SENT      *
C              NN       2ND DIMENSION OF IPCUBE              SENT      *
C              LDP                                           SENT      *
C              LDPNM                                         SENT      *
C                                                                      *
C  IPCUBE:     1 - BEGINNING VOXEL ON THIS LINE                        *
C              2 - ENDING VOXEL ON THIS LINE                           *
C              3 - IX     BEGINNING VOXEL COORDINATES                  *
C              4 - IY                                                  *
C              5 - IZ                                                  *
C                                                                      *
C  CONTAINS:   BCKPJ_LIN   BILINEAR INTERP                            *
C              BCKPJ_FBS    FBS INTERP                                 *
C              BCKPJ        BILINEAR INTERP WITH COMMON PAR            *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE BCKPJ_LIN(CUBE,NNN,DM,B,N,IPCUBE,NN, LDP,LDPNM)

        IMPLICIT NONE

        REAL    :: CUBE(NNN)
        INTEGER :: NNN
        REAL    :: DM(9), B(4,N,N)
        INTEGER :: N,NN
        INTEGER :: IPCUBE(5,NN)
        INTEGER :: LDP,LDPNM

        INTEGER :: I,J,JT,IQX,IQY
        REAL    :: XB,YB,XBB,YBB

#ifdef SP_MP
c$omp   parallel do private(i,j,xb,yb,xbb,ybb,iqx,iqy,jt)
        DO I=1,NN

           XB = ((IPCUBE(3,I)-LDP)*DM(1)  + (IPCUBE(4,I)-LDP) * DM(2) +
     &           (IPCUBE(5,I)-LDP)*DM(3)) + LDPNM
           YB = ((IPCUBE(3,I)-LDP)*DM(4)  + (IPCUBE(4,I)-LDP) * DM(5) +
     &           (IPCUBE(5,I)-LDP)*DM(6)) + LDPNM

           DO J=IPCUBE(1,I),IPCUBE(2,I)

	      JT   = J - IPCUBE(1,I)

	      XBB  = XB + JT * DM(1)
	      YBB  = YB + JT * DM(4)

              IQX  = IFIX(XBB)
              IQY  = IFIX(YBB)

              CUBE(J) = CUBE(J) + 
     &                       B(1,IQX,IQY) + (YBB-IQY) * B(2,IQX,IQY) + 
     &          (XBB-IQX) * (B(3,IQX,IQY) + (YBB-IQY) * B(4,IQX,IQY) )

           ENDDO
        ENDDO
#else
        DO I=1,NN
           XB = ((IPCUBE(3,I)-LDP)*DM(1)  +(IPCUBE(4,I)-LDP)*DM(2) +
     &           (IPCUBE(5,I)-LDP)*DM(3)) + LDPNM
           YB = ((IPCUBE(3,I)-LDP)*DM(4)  +(IPCUBE(4,I)-LDP)*DM(5) +
     &           (IPCUBE(5,I)-LDP)*DM(6)) + LDPNM

           DO J=IPCUBE(1,I),IPCUBE(2,I)

              IQX  = IFIX(XB)
              IQY  = IFIX(YB)

              CUBE(J) = CUBE(J) + 
     &                      B(1,IQX,IQY) + (YB-IQY) * B(2,IQX,IQY) + 
     &          (XB-IQX) * (B(3,IQX,IQY) + (YB-IQY) * B(4,IQX,IQY) )

              XB = XB + DM(1)
              YB = YB + DM(4)
           ENDDO
        ENDDO
#endif
        END

C       --------------------- BCKPJ_FBS -------------------------------

        SUBROUTINE BCKPJ_FBS(CUBE,NNN,DM,
     &                       PROJ,XDER,YDER,XYDER,
     &                       NXLD,NX,IPCUBE,NN, LDP,LDPNM)
C       2D FOURIER-BASED SPLINE INTERPOLATION

        IMPLICIT NONE

        REAL,    INTENT(OUT) :: CUBE(NNN)
        INTEGER              :: NNN
        REAL,    INTENT(IN)  :: DM(9)
        REAL,    INTENT(IN)  :: PROJ (NX  ,NX)
        REAL,    INTENT(IN)  :: XDER (NXLD,NX)
        REAL,    INTENT(IN)  :: YDER (NXLD,NX)
        REAL,    INTENT(IN)  :: XYDER(NXLD,NX)
        INTEGER              :: NXLD,NX,NN
        INTEGER, INTENT(IN)  :: IPCUBE(5,NN)
        INTEGER              :: LDP, LDPNM

        REAL    :: XB,YB
        INTEGER :: I, J

        REAL    :: fbs2

c$omp   parallel do private(i,j,xb,yb)
        DO I=1,NN
           XB = ((IPCUBE(3,I)-LDP)*DM(1)  + (IPCUBE(4,I)-LDP)*DM(2) +
     &          ( IPCUBE(5,I)-LDP)*DM(3)) + LDPNM
           YB = ((IPCUBE(3,I)-LDP)*DM(4)  + (IPCUBE(4,I)-LDP)*DM(5) +
     &          ( IPCUBE(5,I)-LDP)*DM(6)) + LDPNM

           DO J=IPCUBE(1,I),IPCUBE(2,I)

              CUBE(J) = CUBE(J) + 
     &                  FBS2(XB,YB, NXLD,NX,NX, PROJ,NX,
     &                       XDER,YDER,XYDER, .TRUE.)

              XB = XB + DM(1)
              YB = YB + DM(4)
           ENDDO
        ENDDO

        END



C       --------------------- NEVER -------------------------------


#ifdef NEVER
              XBB  = (J-IPCUBE(1,I)) * DM1 + XB
              YBB  = (J-IPCUBE(1,I)) * DM4 + YB
              IQX  = IFIX(XBB+FLOAT(LDPNM))
              IQY  = IFIX(YBB+FLOAT(LDPNM))
              DIPX = XBB+LDPNM-IQX
              DIPY = YBB+LDPNM-IQY
              IF (iqx .ne. iqx1 .or. iqy .ne. iqy1 .or.
     &            dipx1 .ne. dipx .or. dipy .ne. dipy1) then
                  write(6,*) 'simailar test failed'
                  write(6,*) iqx,iqy, ' != ',iqx1,iqy1
                  write(6,*) dipx,dipy, ' != ',dipx1,dipy1
                  stop 'simalar test failed'
              endif

#endif

#ifdef NEVER
              XBB  = (J-IPCUBE(1,I)) * DM1 + XB
              YBB  = (J-IPCUBE(1,I)) * DM4 + YB
              IQX  = IFIX(XBB+FLOAT(LDPNM))
              IQY  = IFIX(YBB+FLOAT(LDPNM))
              DIPX = XBB+LDPNM-IQX
              DIPY = YBB+LDPNM-IQY
              IF (iqx .ne. iqx1 .or. iqy .ne. iqy1 .or.
     &            dipx1 .ne. dipx .or. dipy .ne. dipy1) then
                  write(6,*) 'simailar test failed'
                  write(6,*) iqx,iqy, ' != ',iqx1,iqy1
                  write(6,*) dipx,dipy, ' != ',dipx1,dipy1
                  stop 'simalar test failed'
              endif


C             REFORMATTED FASTER VERSION :
              CUBE(J) = CUBE(J) +
     &               B(IQX,  IQY) + 
     &         DIPY*(B(IQX,  IQY+1) - B(IQX,  IQY)) + 
     &         DIPX*(B(IQX+1,IQY)   - B(IQX,  IQY)  +
     &         DIPY*(B(IQX+1,IQY+1) - B(IQX+1,IQY)  -
     &              (B(IQX,  IQY+1) - B(IQX,  IQY))  ))

C             FASTER VERSION :
              CUBE(J)=CUBE(J)
     &        +B(IQX,IQY)+DIPY*(B(IQX,IQY+1)-B(IQX,IQY))
     &        +DIPX*(B(IQX+1,IQY)-B(IQX,IQY)
     &        +DIPY*(B(IQX+1,IQY+1)-B(IQX+1,IQY)
     &        -B(IQX,IQY+1)+B(IQX,IQY)))

C             SLOWER VERSION :
C     &                 +(1.0-DIPX)*(1.0-DIPY)*B(MAP(IQX,IQY))
C     &                 +     DIPX *(1.0-DIPY)*B(MAP(IQX+1,IQY))
C     &                 +(1.0-DIPX)*     DIPY *B(MAP(IQX,IQY+1))
C     &                 +     DIPX *     DIPY *B(MAP(IQX+1,IQY+1))

#endif
