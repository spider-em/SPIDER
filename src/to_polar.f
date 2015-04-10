C **********************************************************************
C             OPFILEC                              FEB 03 ARDEAN LEITH
C             COSMETIC & ERROR TRAP REWRITE        DEC 12 ARDEAN LEITH
C             NUMBER OF RADII BUG                  AUG 13 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2013  Health Research Inc.,                         *
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
C PURPOSE:  CREATE POLAR REPRESENTATION OF IMAGE WITH CIRCUMFERANCE 
C           (RING) ON THE X DIMENSION
C
C **********************************************************************

        SUBROUTINE TO_POLAR

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        CHARACTER(LEN=MAXNAM) :: FILNAM

        REAL, ALLOCATABLE     :: BUF(:,:)
        REAL, ALLOCATABLE     :: OUT(:)
        CHARACTER(LEN=1)      :: ASK

        INTEGER               :: MAXIM,ITYPE,NX,NY,NZ,IRTFLG,IER,NLET
        INTEGER               :: MR,NR,NOT_USED,ITEMP,NA,NXP,NYP
        INTEGER               :: IXC,IYC,J,I
        REAL                  :: PI,DFI,FI,XS,YS
        real                  :: quadri

        INTEGER, PARAMETER    :: LUNIN  = 20
        INTEGER, PARAMETER    :: LUNOUT = 21

        CHARACTER(LEN=1)      :: NULL = CHAR(0)

C       OPEN INPUT FILE
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNIN,'O',ITYPE,NX,NY,
     &              NZ,MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (NZ > 1)  THEN
           CALL ERRT(101,'DOES NOT WORK ON VOLUMES',IER)
           RETURN
        ENDIF

C       NAME FOR OUTPUT FILE
        CALL FILERD(FILNAM,NLET,NULL,'OUTPUT',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        CALL  RDPRMI(MR,NR,NOT_USED,'INNER and OUTER RADII')
        IF (IRTFLG .NE. 0) GOTO 9999
 
        ITEMP = MIN(((NX-1)/2),((NY-1)/2))

        IF (MR  < 0 ) THEN
           CALL ERRT(102,'INVALID INNER RADIUS',MR)
           GOTO 9999
        ELSEIF (NR > ITEMP) THEN
           CALL ERRT(102,'INVALID OUTER RADIUS',NR)
           GOTO 9999
        ENDIF

        CALL  RDPRMC(ASK,NA,.TRUE.,
     &              'FULL OR HALF CIRCLE (F/H)',NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999


	PI  = 4*DATAN(1.0D0)

	NXP = INT(2*PI*NR)
	IF (ASK == 'H')  NXP = NXP / 2

	NYP = NR - MR + 1

C       OPEN OUTPUT FILE
        MAXIM = 0
	ITYPE = 1
        CALL OPFILEC(0,.FALSE.,FILNAM,LUNOUT,'N',ITYPE,NXP,NYP,
     &               NZ,MAXIM,'POLAR',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       FIND TOTAL NUMBER OF RINGS
        ALLOCATE (BUF(NX,NY),
     &            OUT(NXP), 
     &            STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'TOPOLAR; BUF, OUT',NX*NY + NXP)
           GOTO 9999
        ENDIF

C       INPUT IMAGE
        CALL REDVOL(LUNIN,NX,NY, 1,1, BUF,IRTFLG)

	IXC = NX/2+1         ! CENTER
	IYC = NY/2+1

        DFI = PI / NXP
	IF (ASK == 'F')  DFI = 2 * PI/NXP

C       DEPENDENT ON CIRCULAR CLOSURE IN QUADRI!
	DO J=MR,NR
	  DO I=1,NXP
	     FI     = (I-1) * DFI
	     XS     = COS(FI) * J
	     YS     = SIN(FI) * J
	     OUT(I) = QUADRI(XS+IXC,YS+IYC,NX,NY,BUF) * SQRT(REAL(J))
	  ENDDO

	  CALL WRTLIN(LUNOUT,OUT,NXP,J-MR+1)
	ENDDO

9999    IF (ALLOCATED(BUF)) DEALLOCATE (BUF)
        IF (ALLOCATED(OUT)) DEALLOCATE (OUT)

        CLOSE(LUNIN)
        CLOSE(LUNOUT) 
 
        END

