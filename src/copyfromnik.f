C ++********************************************************************
C                                                                      *
C COPYFROMNIK              CREATED               JAN 2005 ARDEAN LEITH                                            *
C                          READING IBOFF         FEB 2006 ARDEAN LEITH
C                          NPIX8                 DEC 2008 ARDEAN LEITH
C                          SWABTAG               OCT 2009 ARDEAN LEITH
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
C COPYFROMNIK(LUNIN,LUNOUT,MAXMEM)                                                                     *
C 
C PURPOSE: CONVERT NIKON COOLSCAN TIFF  IMAGE TO SPIDER IMAGE
C                                                                     
C PARAMETER:                                                                      *
C       COMMON TIFF TAGS 
C	TIFFTAG_IMAGEWIDTH		256
C	TIFFTAG_IMAGELENGTH		257
C	TIFFTAG_BITSPERSAMPLE		258
C	TIFFTAG_COMPRESSION		259
C	TIFFTAG_PHOTOMETRIC		262
C	TIFFTAG_STRIPOFFSETS		273
C	TIFFTAG_ROWSPERSTRIP		278
C	TIFFTAG_STRIPBYTECOUNTS		279
C	TIFFTAG_RESOLUTIONUNIT		296
C
C od -u -N200 RSp614.tif
C *********************************************************************

        SUBROUTINE COPYFROMNIK(LUNIN,LUNOUT,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        COMMON   BUFRAW(NBUFSIZ)

        CHARACTER(LEN=MAXNAM)         :: FILOLD,FILNEW
        CHARACTER(LEN=1)              :: NULL,CTEMP
        CHARACTER(LEN=4)              :: C4

        INTEGER, PARAMETER            :: NLISTMAX = 16000
        CHARACTER(LEN=NLISTMAX)       :: CTEXT
        INTEGER*2                     :: ITEXT(NLISTMAX/2)
        INTEGER*2                     :: I2VAL
        INTEGER*4                     :: I4VAL
        LOGICAL                       :: SWAB,FOLD,SWABTAG
        INTEGER *8                    :: NPIX8

        NULL = CHAR(0)

22      CALL OPAUXFILE(.TRUE.,FILOLD,NULL,LUNIN,0,
     &                 'O','NIKON TIFF',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
       
        READ(LUNIN,80,IOSTAT=IRTFLG) CTEXT
80      FORMAT(A)
        IF (IRTFLG .NE. 0) RETURN

        WRITE(NOUT,*) ' '

        IMAGIC = TRANSFER(CTEXT(3:4),I2VAL)
        SWAB   = .FALSE.

        IF (IMAGIC .NE. 42) THEN
C          SWAP BYTES FIRST
           WRITE(NOUT,*) CTEXT
 
            WRITE(NOUT,*) 
     &         ' Swapped bytes trying to match TIFF magic number'

            DO I = 1,NLISTMAX,2
               CTEMP          = CTEXT(I:I)
               CTEXT(I:I)     = CTEXT(I+1:I+1)
               CTEXT(I+1:I+1) = CTEMP
            ENDDO
            SWAB = .TRUE.
        ENDIF

        IMAGIC = TRANSFER(CTEXT(3:4),I2VAL)
        WRITE(NOUT,90) IMAGIC
90      FORMAT('  MAGIC NUMBER (42):    ',I7) 

        IF (IMAGIC .NE. 42) THEN
            CALL ERRT(102,'BAD MAGIC NUMBER',IMAGIC)
            RETURN
         ENDIF

c       INIK = INDEX(CTEXT,'Nikon')
c       WRITE(NOUT,*) ' inik: ',inik
C       WRITE(NOUT,*) ' got input: ',ctext
       
        WRITE(NOUT,91) CTEXT(1:2)
91      FORMAT('  ENDEDNESS:                 ',A) 
              
        IBOFFSET = TRANSFER(CTEXT(5:8),I2VAL)
        WRITE(NOUT,92) IBOFFSET
92      FORMAT('  TIFF TAG DIR. OFFSET: ',I7) 
        WRITE(NOUT,*) ' AT ALBANY THIS OFFSET USED TO BE: 799' 
        !iboffset = 799 ! this offset used with Raj's scanner

        IF (IBOFFSET .GE. NLISTMAX) THEN
           CALL ERRT(102,'NLISTMAX OVERFLOW, IBOFFSET',IBOFFSET)
           RETURN
        ENDIF

        NDIR = TRANSFER(CTEXT(IBOFFSET+1:IBOFFSET+2),I2VAL)
        WRITE(NOUT,93) NDIR
93      FORMAT('  TIFF TAGS:            ',I7,/) 

        NUNUSED = 0

        DO IDIR = 1,NDIR
           IOFF = (IDIR - 1) * 12 + IBOFFSET + 3

           IF ((IOFF + 12) .GT. NLISTMAX) THEN
              WRITE(NOUT,*) ' *** SKIPPING SOME TAGS TO AVOID OVERFLOW' 
              WRITE(NOUT,*) ' *** CHECK YOUR RESULTS' 
C              CALL ERRT(102,'NLISTMAX OVERFLOW, IOFF',IOFF)
              EXIT
           ENDIF

           ITAG   = TRANSFER(CTEXT(IOFF  :IOFF+1), I2VAL)
           ITYPE  = TRANSFER(CTEXT(IOFF+2:IOFF+3), I2VAL)

           C4(1:2) = CTEXT(IOFF+6:IOFF+7)
           C4(3:4) = CTEXT(IOFF+4:IOFF+5)
           ICOUNT  = TRANSFER(C4, I4VAL)

           SWABTAG = .FALSE.
           IF (ICOUNT .GT. 65535) THEN
              C4(3:4) = CTEXT(IOFF+6:IOFF+7)
              C4(1:2) = CTEXT(IOFF+4:IOFF+5)
              ICOUNT  = TRANSFER(C4, I4VAL)
              SWABTAG = .TRUE.
           ENDIF

           IF (SWABTAG) THEN
              C4(3:4) = CTEXT(IOFF+10:IOFF+11)
              C4(1:2) = CTEXT(IOFF+8:IOFF+9)
           ELSE
              C4(1:2) = CTEXT(IOFF+10:IOFF+11)
              C4(3:4) = CTEXT(IOFF+8:IOFF+9)
           ENDIF
           IVAL    = TRANSFER(C4,I4VAL)

c          WRITE(NOUT,*)' '
           WRITE(NOUT,98)ITAG,ITYPE,ICOUNT,IVAL
98         FORMAT('  TAG: ',I7,' TYPE:',I3,'  LENGTH: ',I6,
     &            '   VALUE:',I10, i8) 
 
           IREL = IOFF + 8

           IF (ITAG .EQ. 256) THEN
C             NSAM = TRANSFER(CTEXT(IOFF+8:IOFF+9),  I2VAL)
              NSAM = IVAL

c              WRITE(NOUT,94) NSAM
94            FORMAT('  NSAM:                 ',I7) 

           ELSEIF (ITAG .EQ. 257) THEN
C             NROW = TRANSFER(CTEXT(IOFF+8:IOFF+9),  I2VAL)
              NROW = IVAL

c              WRITE(NOUT,95) NROW
95            FORMAT('  NROW:                 ',I7) 

           ELSEIF (ITAG .EQ. 273) THEN
              IHEAD = IVAL

c              WRITE(NOUT,96) IHEAD
96            FORMAT('  HEADER BYTES:         ',I7) 

           ELSE
c             WRITE(NOUT,*) ' UNUSED TAG:      ',ITAG
              NUNUSED = NUNUSED + 1
           ENDIF

        ENDDO

        WRITE(NOUT,*) ' '
        IF (SWABTAG) WRITE(NOUT,*)' SWAPPED BYTES IN TAG CONTENTS'

        WRITE(NOUT,97) NUNUSED
97      FORMAT('  UNUSED TAGS:  ',I7,/) 

        WRITE(NOUT,99) NSAM,NROW,IHEAD
99      FORMAT('  NSAM: ',I7,'  NROW: ',I7,'     HEADER BYTES: ',I7,/) 
      
C       OPEN NEW OUTPUT SPIDER IMAGE FILE WITH SPECIFIED SIZE
        IFORM  = 1
        NSLICE = 1
        MAXIM = 0 
        CALL OPFILEC(0,.TRUE.,FILNEW,LUNOUT,'U',IFORM,NSAM,NROW,NSLICE,
     &              MAXIM,'NEW SPIDER IMAGE',.FALSE.,IRTFLG)
       IF (IRTFLG .EQ. -1) GOTO 22
       IF (IRTFLG .NE. 0) GOTO 9999

        CLOSE(LUNIN)
C       REOPEN RAW FILE AS DIRECT ACCESS, UNFORMATTED REC. 
        LENOPEN = NSAM * 2
        CALL OPAUXFILE(.FALSE.,FILOLD,CHAR(0),LUNIN,LENOPEN,'O',
     &                   ' ',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        NPIX8 = NSAM * NROW 
        NPIX8 = NPIX8 * NSLICE

        FOLD  = .TRUE.

c       swab  = .true.

C       SWAB=.TRUE.  FOLD=.TRUE. IS BEST

        CALL RAW16TOSPI(LUNIN,LUNOUT,NSAM,NPIX8,IHEAD,SWAB,
     &                   FOLD,LENOPEN,BUFRAW,IRTFLG)

9999    CLOSE(LUNIN)
        CLOSE(LUNOUT)

        RETURN
        END
