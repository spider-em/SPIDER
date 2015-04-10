

C++*********************************************************************
C
C PUTPT2.F                         
C                  USED RDPRM3S                    AUG 99 ARDEAN LEITH
C                  MAXNAM                          JUL 14 ARDEAN LEITH
C
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
C PUTPT2(LUN2,NDOC,NSAM,NROW)
C
C PURPOSE:  SUPERPOSE ONTO AN IMAGE, 
C	    PIXELS IN LOCATIONS READ FROM  DOCUMENT FILE    
C            
C PARAMETERS: LUN2	        LOGICAL UNIT NUMBER OF I/O FILE
C	      NDOC	        LOGICAL UNIT NUMBER OF DOCUMENT FILE
C	      NSAM,NROW,NSLICE  DIMENSIONS OF INPUT FILE
C
C--*********************************************************************
 
	SUBROUTINE PUTPT2(LUN2,NDOC,NSAM,NROW,NSLICE)

C       DOC FIL CAN ONLY CONTAIN 99999 KEYS NOW
	PARAMETER (MAXPEAK=99999)
	COMMON NPEAK(MAXPEAK),BUF(4096)

	INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=MAXNAM):: DOCNAM

	DIMENSION     PLIST(10)
        CHARACTER     NULL

#ifdef USE_MPI
        include 'mpif.h'
        INTEGER MYPID, ICOMM, MPIERR
        ICOMM   = MPI_COMM_WORLD
        MPIERR = 0
        CALL MPI_COMM_RANK(ICOMM, MYPID, MPIERR)
#else
        MYPID = -1
#endif

        NULL  = CHAR(0)

	CALL FILERD(DOCNAM,NLETD,NULL,'DOCUMENT',IRTFLG)
	IF (IRTFLG .NE. 0) RETURN

        IF (IFORM .EQ. 3) THEN
           CALL RDPRI3S(NCOLX,NCOLY,NCOLZ,NOT_USED,
     &               'X-COL, Y-COL, Z-COL',IRTFLG)
	   IF (IRTFLG .NE. 0) RETURN
        ELSE
	   CALL RDPRMI(NCOLX,NCOLY,NOT_USED,'X-COL, Y-COL')
           NCOLZ = 0
        ENDIF

C       NCOLP IS COLUMN OF DOC FILE CONTAINING PEAK HEIGHT (INTENSITY)
        CALL RDPRM2(COLP,HEIGHT,NOT_USED,
     &     'INTENSITY COLUMN (& INTENSITY IF NOT IN COL. OF DOC FILE)')
        NCOLP = 0
        IF (COLP .GT. 0) THEN
           NCOLP = COLP
        ELSEIF (COLP .LT. 0) THEN
           HEIGHT = -COLP
        ENDIF

        IF (IFORM .EQ. 3) THEN
           CALL RDPRM3S(XFACT,YFACT,ZFACT,NOT_USED,
     &          'X-FACTOR, Y-FACTOR, Z-FACTOR',IRTFLG)
 	   IF (IRTFLG .NE. 0) RETURN

           CALL RDPRM3S(XOFF,YOFF,ZOFF,NOT_USED,
     &                 'X-OFFSET, Y-OFFSET, Z-OFFSET',IRTFLG)
 	   IF (IRTFLG .NE. 0) RETURN

        ELSE
	   CALL RDPRM2(XFACT,YFACT,NOT_USED,'X-FACTOR, Y-FACTOR')
	   CALL RDPRM2(XOFF,YOFF,NOT_USED,'X-OFFSET, Y-OFFSET')
           ZOFF = 0.0
        ENDIF

	IF (XFACT .EQ. 0.0) XFACT = 1.0
	IF (YFACT .EQ. 0.0) YFACT = 1.0
	IF (ZFACT .EQ. 0.0) ZFACT = 1.0

	NVALU = MAX0(NCOLX,NCOLY,NCOLZ)
        NVALU = MAX0(NVALU,NCOLP)
        NUM   = MAXPEAK

        CALL RDPRAI(NPEAK,MAXPEAK,NUM,0,99999,'KEY NUMBERS',
     &              NULL,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        NOPEN  = 0
        NUMSET = 0
	DO  I=1,NUM
C          GET COORDS FROM DOCUMENT FILE

           CALL UNSAV(DOCNAM,NOPEN,NDOC,NPEAK(I),PLIST,NVALU,LERR,1)
           IF (LERR .NE. 0) GOTO 9300
           NOPEN = 1

           IF (NCOLX.EQ.0) THEN
              IXCOR = NPEAK(I) * XFACT - XOFF
           ELSE
              IXCOR = PLIST(NCOLX) * XFACT + 0.5 - XOFF
           ENDIF

           IYCOR = PLIST(NCOLY) * YFACT + 0.5 - SIGN(1.,YFACT) * YOFF
           IF (YFACT .LT. 0.0) IYCOR = NROW + IYCOR

           IF (IFORM .EQ. 3) THEN
              IZCOR = PLIST(NCOLZ) * ZFACT + 0.5 - SIGN(1.,ZFACT) * ZOFF
           ELSE
              IZCOR = 1
           ENDIF

           IF ((IXCOR .GT. NSAM   .OR. IXCOR .LE. 0) .OR.
     &         (IYCOR .GT. NROW   .OR. IYCOR .LE. 0) .OR.
     &         (IZCOR .GT. NSLICE .OR. IZCOR .LE. 0)) THEN
               IF (MYPID .LE. 0) 
     &         WRITE(NOUT,721) NPEAK(I),IXCOR,IYCOR,IZCOR
721            FORMAT(' *** LOCATION: ',I4,' : (',I5,',',I5,',',I5,
     &                ') OUTSIDE IMAGE, CONTINUING')
           ELSE
              IREC = (IZCOR -1) * NROW + IYCOR
              CALL REDLIN(LUN2,BUF,NSAM,IREC)
              IF (NCOLP .GE. 1)  HEIGHT = PLIST(NCOLP)

              BUF(IXCOR) = HEIGHT
              CALL WRTLIN(LUN2,BUF,NSAM,IREC)
              NUMSET = NUMSET + 1
           ENDIF
        ENDDO

9300    IF (MYPID .LE. 0) WRITE(NOUT,90) NUMSET
 90     FORMAT(' NUMBER OF LOCATIONS SET: ',I5)

        IF (NUMSET .GT. 0) THEN
C          SET FMIN, FMAX AS UNDETERMINED
           CALL SETPRM(LUN2,NSAM,NROW,0.0,0.0,0.0,'U')
        ENDIF
        
        CLOSE(LUN2)
	RETURN
	END
