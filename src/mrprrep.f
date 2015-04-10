
C ++********************************************************************
C 
C  MRPRREP                                                                     *
C                  LONG FILE NAMES                 JAN 89 al
C                  MAXNAM                          JUL 14 ARDEAN LEITH
C                                              
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014 Health Research Inc.,                          *
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
C  MRPRREP(LUN3,LUNP,MAXDIM,IER)                              
C                                                                      
C  PURPOSE:                                                            
C       COMPUTES THE PROJECTION OF A 3-D ARRAY IN ARBITRARY
C       DIRECTION, WITH PARTIAL EXPONENTIAL ATTENUATION. 
C       M.R. 5/87 
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE MRPRREP(LUN3,LUNP,MAXDIM,IER)

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

	COMMON  ADUM(256),BUF(1)

	REAL    DPO(3),DP(3),XP(2),X(3),DM(3,3)
	INTEGER IROW,INDP,MAXDIM,MA,IMA

        CHARACTER(LEN=MAXNAM) :: FILNAM
        CHARACTER *1          :: NULL

	DATA PI/3.14159265/

        NULL = CHAR(0)

        MAXIM  = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN3,'O',IFORM,NSAM,NROW,NSLICE,
     &                   MAXIM,'THREED',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
	IF (IFORM .NE. 3) GOTO 9100

	CALL RDPRMI(NSAMP,NDUM,NOT_USED,'PROJECTION SAMPLE DIM.')
	NROWP=NSAMP*(FLOAT(NROW)/FLOAT(NSAM))
	IF (NSAMP .EQ. 0) GOTO 9200
	MA=FLOAT(NSAMP)*FLOAT(NROWP)+0.5
C       INITIALIZE BUFFER WITH 1.
	DO  IMA = 1, MA
          BUF(IMA) = 1.
	ENDDO

        MAXIM  = 0
        IFORM  = 1
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNP,'O',IFORM,NSAMP,NROWP,1,
     &                   MAXIM,'OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9300

	IER = 0
	CALL RDPRM(PHI,NOT_USED,'AZIMUTH ANGLE (PHI)')
	CALL RDPRM(THETA,NOT_USED,'TILT ANGLE (THETA)')
	PHIR=PHI*PI/180.
	THETAR=THETA*PI/180.
        CALL RDPRM2(FFF,GGG,NOT_USED,'SCALE FACTS IN EXP,SUM')
        CALL RDPRM2(THRESH,SCON,NOT_USED,
     &       'THRESHOLD VALUE (BELOW),ADD. CONST.')

C       FOLLOWING ARE INCREMENTS THAT WERE IN THE ORIGINAL PROGRAM
C       ON THE SCALE OF THE OBJECT.
	DX=1.
	DY=1.
	DZ=1.
	DXP=1.
	DYP=1.

	DP(1)=FLOAT(NSAMP)/2.
	DP(2)=FLOAT(NROWP)/2.
	DPO(1)=FLOAT(NSAM)/2.
	DPO(2)=FLOAT(NROW)/2.
	DPO(3)=FLOAT(NSLICE)/2.

C CALCULATE THE DIRECTION OF PROJECTION AND PROJECTION PLANE.

	CPHI=COS(PHIR)
	SPHI=SIN(PHIR)
	CTHE=COS(THETAR)
	STHE=SIN(THETAR)

	DM(1,1)=CPHI*CTHE
	DM(2,1)=SPHI*CTHE
	DM(3,1)=-STHE
	DM(1,2)=-SPHI
	DM(2,2)=CPHI
	DM(3,2)=0.0
	DM(1,3)=CPHI*STHE
	DM(2,3)=SPHI*STHE
	DM(3,3)=CTHE
	WRITE(NOUT,200)(DM(I,3),I=1,3)
  200	FORMAT(' ** DIRECTION OF PROJECTION : ',3F8.3)

C COMPUTATION OF THE PROJECTION POINT

        IROW = 0                                 
	X(3) =-DZ
	DO  I1= 1, NSLICE
	X(3)=X(3)+DZ
	X(2)=-DY

	DO   I2 = 1, NROW
	X(2)=X(2)+DY
	X(1)=-DX
C	IROW=FLOAT(I1-1)*FLOAT(NROW)+I2+0.5
        IROW = IROW + 1                             
	CALL REDLIN(LUN3, BUF, NSAM, IROW)          
C	READ(LUN3'IROW+1)(BUF(K),K=1,NSAM)

	DO 400 I3 = 1, NSAM
	X(1) = X(1)+DX
	XP(1)=0.0
	XP(2)=0.0

	DO  I4 = 1, 3
	XP(1)=XP(1)+(X(I4)-DPO(I4))*DM(I4,1)
  	XP(2)=XP(2)+(X(I4)-DPO(I4))*DM(I4,2)
	ENDDO
	XP(1)=XP(1)+DP(1)
	XP(2)=XP(2)+DP(2)

	IF (XP(1) .LT. 0.0  .OR.  XP(1) .GT. NSAMP-DXP  .OR.
     &      XP(2) .LT. 0.0  .OR.  XP(2) .GT. NROWP-DYP)  GO TO 400

	ZWX =XP(1)
C	ZWX =XP(1)/DXP
	ZWY =XP(2)
C	ZWY =XP(2)/DYP
	IPX =IFIX(ZWX)
	IPY =IFIX(ZWY)
	DIPX=ZWX-FLOAT(IPX)
	DIPY=ZWY-FLOAT(IPY)

	W=BUF(I3)
C THRESHOLD 3D
        IF(W.LE.THRESH) GOTO 400
C ADD CONSTANT
        W=W+SCON
        INDP = NSAMP
        INDP = INDP*IPY + IPX + NSAM + 1
        AAA=EXP(-(1.0-DIPX)*(1.0-DIPY)*W*FFF)
        BBB=+(1.0-DIPX)*(1.0-DIPY)*W*GGG
	BUF(INDP)=BUF(INDP)*AAA+BBB 
        AAA=EXP(-DIPX*(1.0-DIPY)*W*FFF)
        BBB=+DIPX*(1.0-DIPY)*W*GGG
	BUF(INDP+1)=BUF(INDP+1)*AAA+BBB
	INDP=INDP+NSAMP
        AAA=EXP(-(1.0-DIPX)*DIPY*W*FFF)
        BBB=+(1.0-DIPX)*DIPY*W*GGG
	BUF(INDP)=BUF(INDP)*AAA+BBB
        AAA=EXP(-DIPX*DIPY*W*FFF)
        BBB=+DIPX*DIPY*W*GGG
	BUF(INDP+1)=BUF(INDP+1)*AAA+BBB
  400	CONTINUE
  	ENDDO
	ENDDO

C NOW STORE PROJECTION ARRAY INTO FILE
	DO  I = 1, NROWP
        INDP = NSAMP
        INDP = INDP*(I-1) + NSAM + 1
	CALL WRTLIN(LUNP, BUF(INDP), NSAMP, I)
	ENDDO

	CLOSE(UNIT=LUN3)
	CLOSE(UNIT=LUNP)
	RETURN

9100	CALL ERRT(2,'PROJ3 ',NE)
	RETURN

9200	CALL ERRT(31,'PROJ3 ',NE)
	RETURN

9300	CALL ERRT(4,'PROJ3 ',NE)
	RETURN
	END
