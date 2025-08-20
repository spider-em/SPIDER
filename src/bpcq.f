
C ++********************************************************************
C                                                                      *
C  BPCQ.F                                                                   *
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
C  BPCQ(B,NNNN,NSAM,NROW,CUBE,NX3D,NY3D,NZC,DM)                        *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE BPCQ(B,NNNN,NSAM,NROW,CUBE,NX3D,NY3D,NZC,DM)

        INCLUDE 'PAR.INC'
C	     PAR includes INTEGER LDPX,LDPY,LDPZ,LDPNMX,LDPNMY,NZ1,LDP,NM,LDPNM
	
        DIMENSION    B(NNNN,NROW),CUBE(NX3D,NY3D,NZC)
        DIMENSION    DM(9)

c$omp   parallel do private(k,j,i,kz,iqx,iqy,dipx,dipy,xb,yb,xbb,ybb)
        DO K=1,NZC
           KZ=K-1+NZ1
           DO J=1,NY3D
              XBB = (1-LDPX)*DM(1)+(J-LDPY)*DM(2)+(KZ-LDPZ)*DM(3)
              YBB = (1-LDPX)*DM(4)+(J-LDPY)*DM(5)+(KZ-LDPZ)*DM(6)
              DO I=1,NX3D
                 XB  = (I-1)*DM(1)+XBB
                 IQX = IFIX(XB+FLOAT(LDPNMX))
                 IF (IQX.LT.1 .OR. IQX.GE.NSAM)  GOTO  101
                 YB   = (I-1)*DM(4)+YBB
                 IQY  = IFIX(YB+FLOAT(LDPNMY))
                 IF (IQY.LT.1 .OR. IQY.GE.NROW)  GOTO  101
                 DIPX = XB+LDPNMX-IQX
                 DIPY = YB+LDPNMY-IQY

C FASTER VERSION :
                 CUBE(I,J,K) = CUBE(I,J,K)
     &             +B(IQX,IQY)+DIPY*(B(IQX,IQY+1)-B(IQX,IQY))
     &             +DIPX*(B(IQX+1,IQY)-B(IQX,IQY)
     &             +DIPY*(B(IQX+1,IQY+1)-B(IQX+1,IQY)
     &             -B(IQX,IQY+1)+B(IQX,IQY)))

C ORIG VERSION :
C                CUBE(I,J,K) = CUBE(I,J,K)
C     &                 +(1.0-DIPX)*(1.0-DIPY)*B(MAP(IQX,IQY))
C     &                 +     DIPX *(1.0-DIPY)*B(MAP(IQX+1,IQY))
C     &                 +(1.0-DIPX)*     DIPY *B(MAP(IQX,IQY+1))
C     &                 +     DIPX *     DIPY *B(MAP(IQX+1,IQY+1))

101           CONTINUE
	      ENDDO
	   ENDDO
	ENDDO

        END
