
C++*********************************************************************
C
C MRREFL.F          FILENAMES LENGTHENED          JAN 89 al
C                   OPFILEC                       FEB 03 ARDEAN LEITH
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
C       MRREFL
C
C--*******************************************************************

	SUBROUTINE MRREFL

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        CHARACTER(LEN=MAXNAM)   ::  FNAM,OUTP

        COMMON BUF(500)

        COMMON /COMMUN/ FNAM,OUTP

	LUN1 = 11
	LUN2 = 12

        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FNAM,LUN1,'O',IFORM,NSAM,NROW,NSLICE,
     &                   MAXIM,'INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
 
        MAXIM  = 0
	NSLICE = 1
        CALL OPFILEC(LUN1,.TRUE.,OUTP,LUN2,'U',IFORM,NSAM,NROW,NSLICE,
     &                   MAXIM,' ',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CALL RDPRM(FCT,NOT_USED,'INTERPOLATION FACTOR OF INPUT')
  	CALL RDPRM2(BETA,FBET,NOT_USED,'ILLUMIN. ANGLE, INTENSITY')
	CALL RDPRM(DIF,NOT_USED,'MIX FACTOR')
	IF (FCT.EQ.0) FCT=1.

	PI   = 3.142
	BETA = BETA/180*PI
	SB   = SIN(BETA)
	CB   = COS(BETA)
C       ERROR: THIS WILL NOT WORK ON SHORT IMAGES SHOULD USE REDHED al
	CALL REDLIN(LUN1,BUF,NSAM,0)
	FMIN = BUF(8)
	FMAX = BUF(7)
	CALL NORM3(LUN1,NSAM,NROW,NSLICE,FMAX,FMIN,AV)
	DO L=1,NROW-1
	   CALL REDLIN(LUN1,BUF,NSAM,L)
	   CALL REDLIN(LUN1,BUF(NSAM+1),NSAM,L+1)
	   NNEW = 2*NSAM
           BUF(NNEW+NSAM)=0.
	
           DO I=1,NSAM-1
              BUF(NNEW+I)=0.
              I1=I
              I2=I+1
              I3=I+NSAM
              A1=BUF(I1)-FMIN
              A2=BUF(I2)-FMIN
              A3=BUF(I3)-FMIN
              IF(A1.LE.0.001) GOTO 2
              IF(A2.LE.0.001) GOTO 2
              IF(A3.LE.0.001) GOTO 2
              AD=DIF*(A1+A2+A3)/(3*(FMAX-FMIN))
              XX=(A1-A2)*FCT
              YY=(A1-A3)*FCT
C             INCIDENT SECOND LIGHT VECTOR:
              XS=CB
              YS=0
              ZS=SB
              XS2=XS*XS
              ZS2=ZS*ZS
              SABS=SQRT(XS2+ZS2)

              X2=XX*XX
              Y2=YY*YY
              CABS=SQRT(X2+Y2+1)
              CS=1./CABS
              SI=CS
              R=ABS(A1)
C             SPECULAR REFLECTION:
C             N*L:      
              XS1=XX*XS+ZS
              RX=2*XX*XS1-XS
              RY=2*YY*XS1-YS
              RZ=2*1*XS1-ZS
C             INTENSITY IN Z-COMPONENT, WHICH IS THE ONE ONE SEES:
              IF(RZ.LT.0) THEN
                 RIZ=0
              ELSE
                 RIZ=RZ/SQRT(RX*RX+RY*RY+RZ*RZ)
              ENDIF         
      
              BUF(NNEW+I)=SI-AD+DIF+RIZ*FBET
2	      CONTINUE
           ENDDO
           BUF(NNEW+NSAM)=0.
           CALL WRTLIN(LUN2,BUF(NNEW+1),NSAM,L)
        ENDDO

	DO  I=1,NSAM
           BUF(I)=0.0
        ENDDO
  
        CALL WRTLIN(LUN2,BUF(1),NSAM,NROW)
	CLOSE(LUN1)
	CLOSE(LUN2)

	RETURN
	END

