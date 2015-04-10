C ++********************************************************************
C                                                                      *
C   PRJCQ.F     SPEEDED UP                       FEB 2000 ARDEAN LEITH *
C               COMMON PAR REMOVED               DEC 2010 ARDEAN LEITH *                                                                  *
C               PRJCQ_FBS3                       DEC 2011 G KISHCHENKO *
C               PRJCQ_FBS3 REWRITE               DEC 2011 ARDEAN LEITH *
C               PRJCQ_FBS3 REMOVED               APR 2012 G KISHCHENKO *
C               PRJCQ_N3 REMOVED                 APR 2012 G KISHCHENKO *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@wadsworth.org                                        *
C=*                                                                    *
C=* SPIDER is free software; you can reIribute it and/or               *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* SPIDER is Iributed in the hope that it will be useful,             *
C=* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C=* merchantability or fitness for a particular purpose.  See the GNU  *
C=* General Public License for more details.                           *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program. If not, see <http://www.gnu.org/licenses> *
C=*                                                                    *
C **********************************************************************
C                                                                      *
C  PRJCQ_N(CUBE,LTC,DM,PRJ,N,IPCUBE,NN,LDP,LDPNM)
C                                                                      *
C  PURPOSE:  PROJECTS CUBE INTO B                                      * 
C            PROJECTION FROM VOLUME TO PLANE USING 2D BILINEAR 
C            INTERPOLATION SEE ATTIC OR RCS FOR UNSUCCESSFUL 3D 
C            TRILINEAR & FBS ATTEMPTS
C                                                                      *
C  PARAMETERS:  CUBE    INPUT VOLUME                             SENT  *
C               B       OUTPUT PROJECTION IMAGE                  RET.  *
C               NN      CUBE DIMENSIONS                          SENT  *
C               N       IMAGE DIMENSION (SQUARE)                 SENT  *
C                                                                      *
C  CALLERS: bpcg, bprp
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE PRJCQ_N(CUBE,LTC,DM,PRJ,N,IPCUBE,NN,LDP,LDPNM)

        IMPLICIT NONE

        REAL      :: CUBE(LTC)
        INTEGER   :: LTC
        REAL      :: DM(9) 
        REAL      :: PRJ(N,N)
        INTEGER   :: N
        INTEGER   :: IPCUBE(5,NN)
        INTEGER   :: NN,LDP,LDPNM

        INTEGER   :: I,J,IQX,IQY
        REAL      :: DX,DY,DY1M,DX1M,XB,YB

        DO I=1,NN
           XB = (IPCUBE(3,I)-LDP) * DM(1) + (IPCUBE(4,I)-LDP) * DM(2) +
     &          (IPCUBE(5,I)-LDP) * DM(3) + LDPNM
           YB = (IPCUBE(3,I)-LDP) * DM(4) + (IPCUBE(4,I)-LDP) * DM(5) +
     &          (IPCUBE(5,I)-LDP) * DM(6) + LDPNM

           DO J=IPCUBE(1,I),IPCUBE(2,I)
              IQX    = IFIX(XB)
              IQY    = IFIX(YB)

              DX   = (XB - IQX)
              DY   = (YB - IQY) * CUBE(J)

              DX1M = (1.0 - DX)
              DY1M = (CUBE(J) - DY)

              PRJ(IQX,IQY)    = PRJ(IQX,IQY)     + DX1M * DY1M 
              PRJ(IQX+1,IQY)  = PRJ(IQX+1,IQY)   + DX   * DY1M 
              PRJ(IQX,IQY+1)  = PRJ(IQX,IQY+1)   + DX1M * DY         
              PRJ(IQX+1,IQY+1)= PRJ(IQX+1,IQY+1) + DX   * DY         

              XB = XB + DM(1)
              YB = YB + DM(4)
           ENDDO
        ENDDO

        END

C******************************* PRJCQ_FBS3 ****************************

C       PROJECTION USING 3D FBS INTERPOLATION
C       CALLED WITHIN:  parallel do private(l_th),schedule(static)        
C                       PRJOUT IS LOCAL TO EACH THREAD        

        SUBROUTINE  PRJCQ_FBS3(CUBE, DM, 
     &                         PRJOUT, N,NXLD, XYZ,
     &                         X1, Y1, Z1,
     &                         XY2,XZ2,YZ2)

        IMPLICIT NONE

        REAL              :: CUBE(N,N,N)
        REAL              :: DM(9) 
        REAL              :: PRJOUT(N,N)
        INTEGER           :: N,NXLD
        REAL              :: XYZ(NXLD, N,N)
        REAL              :: X1 (NXLD, N,N)
        REAL              :: Y1 (NXLD, N,N)
        REAL              :: Z1 (NXLD, N,N)
        REAL              :: XY2(NXLD, N,N)
        REAL              :: XZ2(NXLD, N,N)
        REAL              :: YZ2(NXLD, N,N)

        INTEGER           :: I,J,K
        INTEGER           :: IOX, IOY, IOZ

	REAL              :: XB, YB, ZB
        REAL              :: DM1,DM2,DM3

        REAL              :: fbs3

        DM1  = DM(1)
        DM2  = DM(2)
        DM3  = DM(3)

        PRJOUT = 0   ! ARRAY ZERO

        DO K=1,N
          DO J=1,N

             XB  = -N/2*DM(1) + (J-N/2-1)*DM(4) +
     &             (K-N/2-1)*DM(7) + N/2+1

             YB  = -N/2*DM(2) + (J-N/2-1)*DM(5) +
     &             (K-N/2-1)*DM(8) + N/2+1

             ZB  = -N/2*DM(3) + (J-N/2-1)*DM(6) +
     &               (K-N/2-1)*DM(9) + N/2+1

             DO I=1,N
                IOX = IFIX(XB)
                IOY = IFIX(YB)
                IOZ = IFIX(ZB)

                IF (IOX >= 1 .AND. IOX <  N  .AND.
     &              IOY >= 1 .AND. IOY <  N  .AND.
     &              IOZ >= 1 .AND. IOZ <  N) THEN

                   PRJOUT(I,J) = PRJOUT(I,J) + 
     &                  FBS3(XB,YB,ZB,
     &                         NXLD, N, N, N,
     &                         CUBE,N,  XYZ,
     &                         X1, Y1, Z1,
     &                         XY2,XZ2,YZ2)
                ENDIF
                XB = XB + DM1
                YB = YB + DM2
                ZB = ZB + DM3
             ENDDO
          ENDDO
        ENDDO

        END


