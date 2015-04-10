
C ++********************************************************************
C                                                                      *
C  BPCM                                                                *
C                                                                      *
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
C                                                                      *
C  BPCM(B,NNNN,NSAM,NROW,LPRJ,CUBE,NX3D,NY3D,NZC,DM,IOPIC,FIRST)                                                            *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE BPCM(B,NNNN,NSAM,NROW,LPRJ,CUBE,NX3D,NY3D,NZC,DM,
     &                   IOPIC,FIRST)

        DIMENSION         B(NNNN,NROW,LPRJ),CUBE(NX3D,NY3D)
        DIMENSION         DM(9,LPRJ)
        COMMON /PAR/      LDPX,LDPY,LDPZ,LDPNMX,LDPNMY,NZ1
        LOGICAL           FIRST

        DO K=1,NZC
           KZ = K-1+NZ1
           DO J=1,NY3D
              IF (FIRST)  THEN
                 DO I=1,NX3D
                    CUBE(I,J) = 0.0
		 ENDDO
              ELSE
                 CALL REDLIN(IOPIC,CUBE(1,J),NX3D,J+(K-1)*NY3D)
              ENDIF
	   ENDDO

           DO LP=1,LPRJ
c$omp         parallel do private(i,j,iqx,iqy,xb,yb,xbb,ybb,dipx,dipy)
              DO    J=1,NY3D
                 XBB =
     &            (1-LDPX)*DM(1,LP)+(J-LDPY)*DM(2,LP)+(KZ-LDPZ)*DM(3,LP)
                 YBB =
     &            (1-LDPX)*DM(4,LP)+(J-LDPY)*DM(5,LP)+(KZ-LDPZ)*DM(6,LP)
                 DO I=1,NX3D

                   XB   = (I-1)*DM(1,LP)+XBB
                   IQX  = IFIX(XB+FLOAT(LDPNMX))
                   IF (IQX.LT.1 .OR. IQX.GE.NSAM)  GOTO  101
                   YB   = (I-1)*DM(4,LP)+YBB
                   IQY  = IFIX(YB+FLOAT(LDPNMY))
                   IF (IQY.LT.1 .OR. IQY.GE.NROW)  GOTO  101
                   DIPX = XB+LDPNMX-IQX
                   DIPY = YB+LDPNMY-IQY

C faster version :
                   CUBE(I,J) = CUBE(I,J)
     &              +B(IQX,IQY,LP)+DIPY*(B(IQX,IQY+1,LP)-B(IQX,IQY,LP))
     &              +DIPX*(B(IQX+1,IQY,LP)-B(IQX,IQY,LP)
     &              +DIPY*(B(IQX+1,IQY+1,LP)-B(IQX+1,IQY,LP)
     &              -B(IQX,IQY+1,LP)+B(IQX,IQY,LP)))

C orig. version :
c                  CUBE(I,J) = CUBE(I,J)
c     &                 +(1.0-DIPX)*(1.0-DIPY)*B(MAP(IQX,IQY))
c     &                 +     DIPX *(1.0-DIPY)*B(MAP(IQX+1,IQY))
c     &                 +(1.0-DIPX)*     DIPY *B(MAP(IQX,IQY+1))
c     &                 +     DIPX *     DIPY *B(MAP(IQX+1,IQY+1))

101                CONTINUE
	  	 ENDDO
	    ENDDO
	  ENDDO

          DO J=1,NY3D
             CALL WRTLIN(IOPIC,CUBE(1,J),NX3D,J+(K-1)*NY3D)
          ENDDO
	ENDDO
        END
