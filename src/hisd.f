
C++*********************************************************************
C
C HISD.F
C          ADDED NBINS                           OCT 2006 ARDEAN LEITH
C          'HIS D                                OCT 2006 ARDEAN LEITH
C          'HIS' ONLY NOW                        FEB 2014 ARDEAN LEITH
C          WANTDOC DECLARED                      JUN 2015 ARDEAN LEITH
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2015  Health Research Inc.,                         *
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
C    HISD(NDOC)
C
C    PURPOSE:    COMPUTE HISTOGRAM FROM ONE COLUMN OF
C                DOCUMENT FILE, DISPLAY HISTOGRAM ON 
C                LINE PRINTER OR IN DOC FILE.
C
C    NOTES:      FCHAR = HD      - FULL  HISTOGRAM TO RESULTS FILE
C    NOTES:      FCHAR = HD D    - FULL  HISTOGRAM TO DOC FILE
C    NOTES:      FCHAR = HD RD   - RANGE HISTOGRAM TO DOC FILE
C    NOTES:      FCHAR = HD R    - RANGE HISTOGRAM TO RESULTS FILE
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

	SUBROUTINE HISD(LUNDOC)

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'F90ALLOC.INC'

        REAL, POINTER         :: PBUF(:,:)
        CHARACTER(LEN=MAXNAM) :: DOCFIL
        CHARACTER(LEN=MAXNAM) :: COMMENT
        REAL, ALLOCATABLE     :: FREQ(:)
	REAL                  :: PLIST(4)
        LOGICAL               :: ADDEXT,ASKNAM,ISOLD,APPEND
        LOGICAL               :: MESSAGE,NEWFILE
        LOGICAL               :: WANTDOC

	IF (FCHAR(4:4) == 'D' ) THEN
	  WRITE(NOUT,*) " OBSOLETE OPERATION, USE 'HD' NOW"
        ENDIF

        MAXX = 0
        MAXY = 0
        CALL GETDOCDAT('INPUT DOCUMENT',.TRUE.,DOCFIL,NDOC,
     &                 .TRUE.,MAXX, MAXY,PBUF,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

	CALL RDPRI1S(ICOL,NOT_USED,
     &              'REGISTER (COLUMN) USED FOR HISTOGRAM',IRTFLG)

        IF (ICOL <= 0 .OR. ICOL >= MAXX) THEN
           CALL ERRT(100,'HISD',NE)
           GOTO 999
        ENDIF

        NBINS = 128
	CALL RDPRI1S(NBINS,NOT_USED,'NUMBER OF BINS',IRTFLG)
        IF (NBINS <= 0) THEN
           CALL ERRT(100,'MUST HAVE AT LEAST 1 BIN',NE)
           GOTO 999
        ENDIF

        ALLOCATE(FREQ(NBINS),STAT=IRTFLG)
        IF (IRTFLG. NE. 0) THEN 
           CALL ERRT(46,'HISD, NBINS',NBINS)
           GOTO 999
        ENDIF

C       IF ONE OF THE FCHAR LETTERS IS 'R' ASK USER FOR THE RANGE
	IF (FCHAR(4:4) == 'R' .OR. FCHAR(5:5) == 'R') THEN
	  CALL RDPRM2(HMIN,HMAX,NOT_USED,'HISTOGRAM RANGE MIN & MAX')

	ELSE
C         FIND MIN & MAX FOR HISTOGRAM BY READING DOC FILE 
	  HMAX  = -1.0E23
	  HMIN  = -HMAX

          DO IKEY = 1,MAXY
             IF (PBUF(1, IKEY) > 0) THEN
C               KEY EXISTS
                BVAL = PBUF( ICOL + 1, IKEY)
	        HMAX = MAX(HMAX,BVAL)
	        HMIN = MIN(HMIN,BVAL)
             ENDIF
          ENDDO
        ENDIF

C       ZERO THE HISTOGRAM
        FREQ = 0.0 

C       CALCULATE HISTOGRAM HERE
        HDIFF  = HMAX - HMIN
        FF     = (NBINS - 1.0) / HDIFF
        BSIZ   = 1.0 / FF
        FNUMEL = 0
        FKEYS  = 0.0
        Y      = 0.0
        HAV    = 0.0
        HAV2   = 0.0

        DO IKEY = 1,MAXY
           IF (PBUF( 1, IKEY) .GT. 0) THEN
C            KEY EXISTS
             BVAL  = PBUF( ICOL + 1, IKEY)
             JBIN  = INT((BVAL - HMIN) * FF) + 1.5
             FKEYS = FKEYS + 1

             IF (JBIN >=1 .AND. JBIN <= NBINS)  THEN
C               WITHIN HISTOGRAM RANGE
                FREQ(JBIN) = FREQ(JBIN) + 1.0
                FNUMEL     = FNUMEL + 1

                HAV        = HAV  + BVAL 
                HAV2       = HAV2 + DBLE(BVAL) * DBLE(BVAL)
             ENDIF
           ENDIF
        ENDDO

        ADDEXT   = .TRUE.
        ASKNAM   = .TRUE.
        ISOLD    = .FALSE.
        APPEND   = .FALSE. 
        MESSAGE  = .TRUE. 

        CALL OPENDOC(DOCFIL,ADDEXT,NLET,LUNDOC,LUNDOCT,ASKNAM,
     &           'OUTPUT DOC FILE (OR * IF NONE)',ISOLD,APPEND,MESSAGE,
     &            NEWFILE,IRTFLG)

        WANTDOC = (IRTFLG == 0)

        IF (WANTDOC) THEN
C          OUTPUT TO DOC FILE

C                    123456789 123456789 123456789 123456789 123456789 
           COMMENT = ' BIN   BOUNDARY     OCCUPANCY'
           CALL LUNDOCPUTCOM(LUNDOCT,COMMENT,IRTFLG)

           DO IBIN=1,NBINS
              PLIST(1) = HMIN + (BSIZ *(IBIN - 1))
              PLIST(2) = FREQ(IBIN)
              CALL LUNDOCWRTDAT(LUNDOCT,IBIN,PLIST,2,IRTFLG)
           ENDDO
           CLOSE(LUNDOC)

        ELSE
C          OUTPUT IN RESULTS FILE EVEN IF INTERACTIVE INPUT
           WRITE(NDAT,*) ' '
           CALL GRAPHS(NDAT,FREQ,NBINS,1,1,1.0,IRTFLG)

           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(100,'GRAPHING',NE)
              GOTO 999
           ENDIF
        ENDIF

C       FIND MAXIMUM FREQ.  IN HISTOGRAM (HISMAX) & LOCATION
        HISMAX = FREQ(1)
        MAXBIN = 1
        DO  IBIN = 2,NBINS
           IF (FREQ(IBIN) >= HISMAX) THEN
              HISMAX = FREQ(IBIN)
              MAXBIN = IBIN
           ENDIF
        ENDDO

C       CONVERT LOCATION MAXBIN TO CORRESPONDING IMAGE INTENSITY (MODE)
        IF (MAXBIN == 1 .OR. MAXBIN == NBINS) THEN
           BMODE  = MAXBIN - 1
        ELSE
           YM1    = FREQ(MAXBIN-1)
           Y1     = FREQ(MAXBIN+1)
           BMODE  = FLOAT(MAXBIN-1) +  
     &                (YM1-Y1) *.5 / (YM1+Y1-2.* HISMAX)
        ENDIF
        HMODE  = HMIN + BMODE * BSIZ

        DTOP  = HAV2 - HAV * HAV / DBLE(FNUMEL)

        IF (DTOP < 0.0) THEN
C          SQRT OF NEGATIVE NUMBER
           WRITE(NOUT,*) '*** HISD, SQRT(',DTOP,') IMPOSSIBLE' 
           CALL ERRT(100,'HIST',NE)
           GOTO 999

        ELSEIF (FNUMEL == 1) THEN
           CALL ERRT(101,'CAN NOT DETERMINE S.D. (ONLY ONE VALUE)',NE)
           GOTO 999
        ENDIF

        HSIG   = DSQRT(DTOP / DBLE(FNUMEL - 1))
        HAV    = HAV  / DBLE(FNUMEL)
        FNBINS = NBINS

        WRITE(NOUT,*) ' '
        WRITE(NOUT,90) HMIN,HMAX,FKEYS,FNUMEL,
     &                 FNBINS,BSIZ,HAV,HMODE,HSIG

90      FORMAT('  HISTOGRAM RANGE:   ',G11.4,'   .........     ',G11.4,/,
     &         '  FILE KEYS:         ',G11.4,'   HIST. KEYS:   ',G11.4,/,
     &         '  NO. OF BINS:       ',G11.4,'   BIN SIZE:     ',G11.4,/,
     &         '  HIST. MEAN:        ',G11.4,'   HIST. MODE:   ',G11.4,/,
     &         '  HIST. S.D.:        ',G11.4)

        WRITE(NOUT,*) ' '

        IF (NOUT .NE. NDAT) THEN
C           ECHO TO RESULTS FILE
            WRITE(NDAT,*) ' '
            WRITE(NDAT,90) HMIN,HMAX,FKEYS,FNUMEL,
     &                     FNBINS,BSIZ,HAV,HMODE,HSIG
            WRITE(NDAT,*) ' '
        ENDIF

999     IF (ASSOCIATED(PBUF)) DEALLOCATE(PBUF)
        IF (ALLOCATED(FREQ))  DEALLOCATE(FREQ)

        END

