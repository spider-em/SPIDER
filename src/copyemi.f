
C ++********************************************************************
C
C COPYEMI                     MODIFIED FROM COPYCCP4 FEB 02 ArDean Leith         
C                             NPIX8                  DEC 08 ArDean Leith
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
C COPYEMI(LUNSPI,LUNNEW)
C                                                                      
C PURPOSE: CONVERTS EMI IMAGES TO SPIDER IMAGES
C
C NOTES: DATA IN EMI FILE:
C	 IMODE   0 : IMAGE STORED AS INTEGER*1  
C                1 : IMAGE STORED AS INTEGER*2
C                2 : IMAGE STORED AS REALS
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE COPYEMI(LUNSPI,LUNEMI)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
 
        INTEGER * 1     IDATA(4)
        REAL *4         DATA
        EQUIVALENCE(DATA,IDATA)

        INTEGER * 2     I2BUF(1000)
        COMMON          BUF(NBUFSIZ),I2BUF

        CHARACTER(LEN=MAXNAM) :: FILNAM,EMIFILE
	CHARACTER (LEN=1)     :: NULL
        LOGICAL               :: FLIP,FOLD
        INTEGER * 8           :: NPIX8

        NULL = CHAR(0)

C       OPEN EMI FILE AS DIRECT ACCESS, UNFORMATTED, RECL= 2000 BYTES
        LENREC = 2000
        CALL OPAUXFILE(.TRUE.,EMIFILE,DATEXC,LUNEMI,LENREC,'O',
     &                       'EMI INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
 
        READ(LUNEMI,REC=1,IOSTAT=IERR) I2BUF
        IF (IERR .NE. 0) THEN
           WRITE(NOUT,*) '*** ERROR: (',IERR,') READING EMI HEADER'
           CALL ERRT(100,'COPYEMI',NE)
           GOTO 9999
        ENDIF

	CALL GETHEDEMI(I2BUF,NSAM,NROW,NSLICE,IMODE,
     &                 FLIP,FOLD,IOFFSET,IRTFLG)

C       OPEN SPIDER OUTPUT FILE	
        IFORM  = 1
        IF (NSLICE .GT. 1) IFORM = 3
        MAXIM  = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNSPI,'U',IFORM,NSAM,NROW,NSLICE,
     &             MAXIM,'SPIDER OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999


C       EXTRACT DATA FROM EMI FILE AFTER HEADER & PUT IN SPIDER FILE

C       NPIX8 = TOTAL NUMBER OF PIXELS IN EMI FILE
        NPIX8 = NSAM * NROW 
        NPIX8 = NPIX8 * NSLICE

C       CLOSE EMI FILE
        CLOSE(LUNEMI)

C       REOPEN EMI FILE AS NSAM*IMODE/8 BYTE, DIRECT ACCESS, UNFORMATTED
        LENOPEN = NSAM * (IMODE / 8)
        CALL OPAUXFILE(.FALSE.,EMIFILE,NULL,LUNEMI,LENOPEN,'O',
     &                    ' ',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (IMODE .EQ. 8) THEN
C          8 BIT INTEGER EMI INPUT FILE

           CALL RAW8TOSPI(LUNEMI,LUNSPI,NSAM,NPIX8,IOFFSET,.TRUE.,
     &                   LENOPEN,BUF,IRTFLG)

        ELSEIF (IMODE .EQ. 16) THEN
C          16 BIT INTEGER EMI FILE (HEADER LENGTH DIVISABLE BY TWO)
           CALL RAW16TOSPI(LUNEMI,LUNSPI,NSAM,NPIX8,IOFFSET,FLIP,
     &                   FOLD,LENOPEN,BUF,IRTFLG)

        ELSEIF (IMODE .EQ. 32) THEN
C          32 BIT FOATING POINT IMAGE
           NFLIP   = -1
           IOFFSET = 802
           FOLD    = .FALSE.
           FLIP    = .TRUE.
           CALL RAW32TOSPI(LUNEMI,LUNSPI,NSAM,NPIX8,IOFFSET,
     &                    NFLIP,LENOPEN,BUF,IRTFLG)

C          SLOWER ALTERNATIVE (NEED TO OPEN WITH LENOPEN + 1)
C          CALL RAWTOSPI(LUNEMI,LUNSPI,NSAM,NPIX8,IOFFSET,
C    &                    IMODE,FOLD,FLIP,IDUM,IRTFLG)

        ELSE
           CALL ERRT(102,'CAN NOT CONVERT EMI MODE',IMODE)
        ENDIF


9999    CLOSE(LUNSPI)
        CLOSE(LUNEMI)

        RETURN
        END

C ------------------------ FOLDNFLIP  -------------------------------

      SUBROUTINE FOLDNFLIP(I1IN,I1OUT)

C     FOLDS NEGATIVES AND ASSIGNS 11IN TO I1OUT AND FLIPS BYTES 
C     WITHIN WORD

      INTEGER * 1    I1IN(2),I1OUT(2)

      I1OUT(1) = I1IN(2)
      IF (I1OUT(1) .LT. 0) I1OUT(1) = 256 + I1OUT(1)

      I1OUT(2) = I1IN(1)
      IF (I1OUT(2) .LT. 0) I1OUT(2) = 256 + I1OUT(2)

      RETURN
      END


C ++********************************************************************
C                                                                      *
C  GETHEDEMI                                                           *            *
C                                                                      *
C **********************************************************************
C **********************************************************************
C                                                                      *
C  GETHEDEMI(HEADBUF,NSAM,NROW,NSLICE,IMODE,DMIN,DMAX,
C             DMEAN,RMS,NSYMBT,FLIP,MACHST,IRTFLG)
C                                                                      *
C  PURPOSE:     DECODE EMISPEC HEADER                      *                                          *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C MAP/IMAGE HEADER FORMAT
C 
C 	LENGTH = 1024 BYTES, ORGANIZED AS 56 LONG WORDS 
C 
C  1	NSAM            # OF COLUMNS 
C  2	NROW		# OF ROWS
C  3	NSLICE          # OF SECTIONS 	
C  4	IMODE		 32	IMAGE : 32-BIT REALS			*
C 
C 		DATA RECORDS FOLLOW.
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	SUBROUTINE GETHEDEMI(I2BUF,NSAM,NROW,NSLICE,
     &                       IMODE,FLIP,FOLD,IOFFSET,IRTFLG)

        INCLUDE 'CMBLOCK.INC'

        INTEGER * 2               I2BUF(*),I2FLIP
        LOGICAL ::                BIGENDARCH,BIGENDED,BIGENDFILE
        LOGICAL ::                FLIP,FOLD

C       GET CURRENT ARCHITECTURE ENDED-NESS
        BIGENDARCH = BIGENDED(0)
 
C       GET CURRENT FILE ENDED-NESS FLAG

        IMODE    = 32
        FOLD    = .FALSE.
        FLIP    = .TRUE.
        IOFFSET =  802

        CALL FOLDNFLIP(I2BUF(396),I2FLIP)
        NSAM = I2FLIP

        CALL FOLDNFLIP(I2BUF(398),I2FLIP)
        NROW = I2FLIP

        NSLICE  = 1
	 
C       WRITE OUT HEADER INFORMATION
        IF (VERBOSE) THEN    	
	   WRITE(NOUT,1000) NSAM,NROW,NSLICE,IMODE,IOFFSET

1000	   FORMAT(/
     &     7X,'Number of columns, rows, sections ........ ',3I6/
     &     7X,'Mode ..................................... ',I6/
     &     7X,'Header offset ............................ ',3I6//)


        ENDIF

	RETURN
        END

#ifdef NEVER
        DO I=1,NHEADER/2
            CALL FOLDNFLIP(I2BUF(I),I2FLIP(i))
        END DO

        DO I = 1,NHEADER/2
           IF (I2FLIP(I) .EQ. 1024) 
     &         WRITE(6,*) 'I2FLIP(',I,'): ',I2FLIP(I),I2FLIP(I+1),
     &                     I2FLIP(I+2),I2FLIP(I+3)
           IF (I2FLIP(I) .EQ. 1026) 
     &         WRITE(6,*) 'I2FLIP(',I,'): ',I2FLIP(I),I2FLIP(I+1),
     &                     I2FLIP(I+2),I2FLIP(I+3)
        ENDDO
#endif       


C
