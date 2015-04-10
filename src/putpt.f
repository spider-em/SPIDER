
C ++********************************************************************
C
C PUTPT.FOR        LONG FILNAMES                   JAN 89 al
C                  SUPERCEDES SECTION OF PICKPT           jf
C                  MAXNAM                          JUL 14 ARDEAN LEITH
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C  PUTPT(LUN2,NDOC,NSAM,NROW)
C
C  PURPOSE:   CREATE CCF/PEAK MAP BY SUPERPOSING, ONTO AN IMAGE OR CCF, 
C	      PIXELS IN PEAK POSITIONS ACCORDING TO DOCUMENT FILE    
C
C  PARAMETERS:    LUN2          LOGICAL UNIT NUMBER OF INPUT FILE
C	          NDOC          LOGICAL UNIT NUMBER OF DOCUMENT FILE
C	          NSAM,NROW     DIMENSIONS OF INPUT FILE
C
C--*********************************************************************
 
	SUBROUTINE PUTPT(LUN2,NDOC,NSAM,NROW,NSLICE)
 
	PARAMETER (MAXPEAK=9999)
	COMMON     NPEAK(MAXPEAK),BUF(1)

	INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=MAXNAM):: DOCNAM

	DIMENSION     PLIST(10)
        CHARACTER     NULL

        NULL  = CHAR(0)

        FMAXO = FMAX
        FMINO = FMIN
        IF (IMAMI .NE. 1) 
     &     CALL NORM3(LUN2,NSAM,NROW,NSLICE,FMAXO,FMINO,AVDO)

	CALL FILERD(DOCNAM,NLETD,NULL,'DOCUMENT',IRTFLG)
	IF (IRTFLG .NE. 0) RETURN

	CALL RDPRMI(NCOLX,NCOLY,NOT_USED,'X-COL, Y-COL')

C       NCOLP IS COLUMN OF DOC FILE CONTAINING PEAK HEIGHT (INTENSITY)
        CALL RDPRMI(NCOLP,NDUM,NOT_USED,'PEAK COLUMN')

	CALL RDPRM2(XFACT,YFACT,NOT_USED,'X-FACTOR, Y-FACTOR')
	IF (XFACT.EQ.0.0) XFACT=1.0
	IF (YFACT.EQ.0.0) YFACT=1.0

	NVALU = MAX0(NCOLX,NCOLY)
        NVALU = MAX0(NVALU,NCOLP)
        NUM   = MAXPEAK

        CALL RDPRAI(NPEAK,MAXPEAK,NUM,0,9999,'.KEY NUMBERS',
     &              NULL,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        NOPEN = 0
	WRITE(NOUT,*) ' NUMBER OF LOCATIONS SET: ' ,NUM

	DO  I=1,NUM
C          COORDS FROM DOCUMENT FILE

           CALL UNSAV(DOCNAM,NOPEN,NDOC,NPEAK(I),PLIST,NVALU,LERR,1)
           IF (LERR.NE.0) GOTO 9300
           NOPEN=1

           IF (NCOLX.EQ.0) THEN
             IXCOR = NPEAK(I) * XFACT
           ELSE
             IXCOR = PLIST(NCOLX)*XFACT+0.5
           ENDIF

           IYCOR = PLIST(NCOLY)*YFACT+0.5
           IF (YFACT.LT.0.0) IYCOR=NROW+IYCOR

           IF((IXCOR+1.GT.NSAM .OR. IXCOR-1.LE.0) .OR.
     &        (IYCOR+1.GT.NROW .OR. IYCOR-1.LE.0)) THEN
              WRITE(NOUT,721) NPEAK(I)
721           FORMAT(' *** PEAK NO.',I4,' OUT OF LIMITS, CONTINUING.')

           ELSE
             CALL REDLIN(LUN2,BUF,NSAM,IYCOR)
             BUF(IXCOR+1) = FMAXO
             BUF(IXCOR-1) = FMAXO
       
             HEIGHT = PLIST(NCOLP)
             BUF(IXCOR)  = HEIGHT

             CALL WRTLIN(LUN2,BUF,NSAM,IYCOR)
             CALL REDLIN(LUN2,BUF,NSAM,IYCOR-1)
             BUF(IXCOR-1) = FMAXO
             BUF(IXCOR)   = FMAXO
             BUF(IXCOR+1) = FMAXO
             CALL WRTLIN(LUN2,BUF,NSAM,IYCOR-1)
             CALL REDLIN(LUN2,BUF,NSAM,IYCOR+1)
             BUF(IXCOR-1) = FMAXO
             BUF(IXCOR)   = FMAXO
             BUF(IXCOR+1) = FMAXO
             CALL WRTLIN(LUN2,BUF,NSAM,IYCOR+1)
           ENDIF
	ENDDO

C       ZERO THE STATISTICS ON FILE
        CALL SETPRM(LUN2,NSAM,NROW,0.0,0.0,0.0,'U')
        
9300	CLOSE(LUN2)
	RETURN
	END
