
C++*********************************************************************
C
C  COPYE.F   -- CREATED JULY 17 1989 al                     
C               MAXNAM                             JUL 14 ARDEAN LEITH
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
C    COPYE(LUN1,LUN2,NSAM,NROW,NSLICE)
C
C    PURPOSE:   COPIES A SPIDER IMAGE FILE INTO AN EDITABLE (ASCII
C               TEXT) IMAGE FILE
C
C    PARAMETERS:   LUN       INPUT FILE UNIT NUMBER             (SENT)
C                  LUN2      OUTPUT FILE UNIT NUMBER (OPENED)   (SENT)
C                  NSAM      SAMPLES PER LINE IN IMAGE          (SENT)
C                  NROW      NO. OF ROWS IN IMAGE               (SENT)
C                  NSLICE    NO. OF SLICES IN IMAGE             (SENT)
C--*********************************************************************

	SUBROUTINE COPYE(LUN1,LUN2,NSAM,NROW,NSLICE)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        COMMON /IOBUF/  BUF(NBUFSIZ)

        CHARACTER(LEN=MAXNAM) :: FILNAM

	CHARACTER *60   FMT,FORMIN
        CHARACTER *15   FMTDEF
        CHARACTER(LEN=1) :: NULL = CHAR(0)

        DATA FMTDEF /'(6(1X,1PG12.4))'/

        CALL OPAUXFILE(.TRUE.,FILNAM,DATEXC,LUN2,0,'N',
     &                 'OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

20      CALL RDPRMC(FORMIN,NC,.TRUE.,
     &     'FORMAT DESCRIPTION (OR <CR> FOR 6(1X,G12.4))',
     &     NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9998

	IF (NC .GT. 58) THEN
          CALL ERRT(101,' *** FORMAT LENGTH MUST BE <59 CHARACTERS.',NE)
	  GOTO 20

        ELSEIF (NC .EQ. 0) THEN
C         USE DEFAULT FORMAT
          FMT(1:60) = FMTDEF

        ELSE
C         SET OUTPUT FORMAT
          FMT(1:60) = '(' // FORMIN(1:58) // ')'
        ENDIF

        IRECT = NROW * NSLICE

        WRITE(LUN2,*) NSAM,NROW,NSLICE

	DO  I=1,IRECT
          CALL REDLIN(LUN1,BUF,NSAM,I)
          WRITE(LUN2,FMT) (BUF(K),K=1,NSAM)
	ENDDO

9998    CLOSE(LUN2)
9999    CLOSE(LUN1)

	RETURN
        END

