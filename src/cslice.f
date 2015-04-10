
C++*********************************************************************
C
C  CSLICE.F         FILE NAMES LENGTHENED                 ARDEAN LEITH
C                   USED OPFILE                    NOV 00 ARDEAN LEITH
C
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
C
C  CSLICE:  SELECT CENTRAL SLICE OF A 3-D IMAGE WITH ARBITRARY
C           AZIMUTH AND INCLINATION.
C
C--*********************************************************************

	SUBROUTINE CSLICE

	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'

        COMMON /IOBUF/ A0(NBUFSIZ)
	COMMON  BUF(1)
 
        CHARACTER(LEN=MAXNAM) :: FILNAM,FILOUT
        
        INTEGER, PARAMETER    :: LUNI = 17
        INTEGER, PARAMETER    :: LUNO = 16
        REAL, PARAMETER       :: PI   = 3.14159

        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNI,'O',IFORM,NX,NY,NZ,
     &                   MAXIM,'INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

	IF (IFORM .NE. 3) THEN
           CALL ERRT(101,'OPERATION ONLY WORKS ON VOLUMES',NDUM)
           RETURN
        ENDIF

        FMININ = FMIN
	NY3    = NY*NZ
	S3     = SQRT(3.)
	MAXREC = NY*NZ
	NX2    = NX*S3+0.5
	NY2    = NY*S3+0.5
	NXH    = NX2/2+1
	NYH    = NY2/2+1

        MAXIM  = 0
	IFORM  = 1
        CALL OPFILEC(0,.TRUE.,FILOUT,LUNO,'U',IFORM,NX2,NY2,1,
     &                   MAXIM,'OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

	CALL RDPRM1S(PHI,NOT_USED,'AZIMUTH',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

	CALL RDPRM1S(THETA,NOT_USED,'INCLINATION',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

	PHI   = PHI   * PI/180.
	THETA = THETA * PI/180.

 	KXM   = NX/2+1
 	KYM   = NY/2+1
 	KZM   = -9999   !NZ/2+1

	CALL RDPRI3S(KXM,KYM,KZM,NOT_USED,
     &               'X,Y, & Z POSITION',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000
 
        IF (KZM == -9999) THEN
C           NOT ENTERED, MAYBE LEGACY PROCEDURE
 	    KZM   = NZ/2+1
	    CALL RDPRI1S(KZM,NDUM,NOT_USED,'Z POSITION',IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9000
        ENDIF

	SPHI   = SIN(PHI)
	CPHI   = COS(PHI)
	STHETA = SIN(THETA)
	CTHETA = COS(THETA)
	XFACT  = CTHETA*SPHI
	YFACT  = CTHETA*CPHI

C       DIMENSIONS OF 2-D FILE ARE NX2 BY NY2
	DO IY=1,NY2
           DO IX=1,NX2
              X = KXM + ((FLOAT(IX-NXH))*CPHI)+((FLOAT(IY-NYH))*XFACT)
              Y = KYM - ((FLOAT(IX-NXH))*SPHI)+((FLOAT(IY-NYH))*YFACT)
              Z = KZM + (FLOAT(IY-NYH))*STHETA

C             DETERMINE THE 8 SURROUNDING COEFFICIENTS
              KXBOT = X
              KXTOP = KXBOT+1
              XDEC  = X -(FLOAT(KXBOT))
              XREM  = 1. -XDEC

              KYBOT = Y
              KYTOP = KYBOT+1
              YDEC  = Y -(FLOAT(KYBOT))
              YREM  = 1. -YDEC

              KZBOT = Z
              KZTOP = KZBOT+1
              ZDEC  = Z -(FLOAT(KZBOT))
              ZREM  = 1. -ZDEC

C             CHECK IF COORDINATES ARE INSIDE THE VOLUME;  
C             CONTINUE IF THEY ARE, OTHERWISE, SET = 0.
              IF (KXTOP .LE. NX .AND. KYTOP .LE. NY .AND.
     &            KZTOP .LE. NZ .AND.
     &            KXBOT .GE. 1  .AND. KYBOT .GE. 1 .AND. KZBOT .GE. 1)
     &            THEN
                 CONTINUE
              ELSE
                 BUF(IX) = FMININ
                 CYCLE
              ENDIF

              IREC1 = (KZBOT-1) * NY+KYBOT
              IF (IREC1 > MAXREC) THEN
                 CALL ERRT(102,RECORD NUMBER OUT OF VOLUME,IREC1)
                 GOTO 9000
              ENDIF

              CALL REDLIN(LUNI,A0,NX,IREC1)

              PT1 = A0(KXBOT)
              PT2 = A0(KXTOP)

              IREC2 = (KZBOT-1)*NY+KYTOP
              IF (IREC2 > MAXREC)  THEN
                 CALL ERRT(102,RECORD NUMBER OUT OF VOLUME,IREC2)
                 GOTO 9000
              ENDIF

              CALL REDLIN(LUNI,A0,NX,IREC2)
              PT3 = A0(KXBOT)
              PT4 = A0(KXTOP)

              IREC3 = (KZTOP-1)*NY+KYBOT
              IF (IREC3 > MAXREC)  THEN
                 CALL ERRT(102,RECORD NUMBER OUT OF VOLUME,IREC3)
                 GOTO 9000
              ENDIF

              CALL REDLIN(LUNI,A0,NX,IREC3)
              PT5 = A0(KXBOT)
              PT6 = A0(KXTOP)

              IREC4 = (KZTOP-1)*NY+KYTOP
              IF (IREC4 > MAXREC)  THEN
                 CALL ERRT(102,RECORD NUMBER OUT OF VOLUME,IREC4)
                 GOTO 9000
              ENDIF

              CALL REDLIN(LUNI,A0,NX,IREC4)
              PT7 = A0(KXBOT)
              PT8 = A0(KXTOP)

C             WRITE(4,8888)PT1,PT2,PT3,PT4,PT5,PT6,PT7,PT8
C8888         FORMAT(' PTS 1-8 = ',8F7.3)

C             INTERPOLATE
	      BUF(IX) = ZREM*(XDEC*(YREM*PT2+YDEC*PT4) +
     &           XREM*(YREM*PT1+YDEC*PT3))+
     &           ZDEC*(XDEC*(YREM*PT6+YDEC*PT8)+
     &           XREM*(YREM*PT5+YDEC*PT7))

	   ENDDO
	   CALL WRTLIN(LUNO,BUF,NX2,IY)

	ENDDO

9000	CLOSE(LUNO)
	CLOSE(LUNI)

	END
