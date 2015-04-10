
C++*********************************************************************
C
C  COPYU.F    -- CREATED MAR 89
C                REWRITTEN                           FEB 99 ARDEAN LEITH
C                NATIVE BYTE_ORDER                   JUL 09 ARDEAN LEITH
C
C **********************************************************************
C *  AUTHOR: A. LEITH                                                      *
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
C   COPYU(LUNO,LUNN,NSAM,NROW,NSLICE)
C
C   PURPOSE:      CONVERT SPIDER IMAGE FILES TO 8 BIT RAW
C
C   CALLED BY:    UTIL2
C
C--*********************************************************************

        SUBROUTINE COPYU(LUNO,LUNN,NSAM,NROW,NSLICE)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER (LEN=MAXNAM) :: FILNEW
        CHARACTER (LEN=1)      :: NULL,ANS
        LOGICAL                :: HEADER,NORMAL

        NULL = CHAR(0)

        CALL RDPRMC(ANS,NA,.TRUE.,
     &    'DO YOU WANT TO NORMALIZE THE OUTPUT TO 0-255? (Y/N)',
     &	   NULL,IRT)
        NORMAL = (ANS .NE. 'N' .AND. ANS .NE. 'n') 

        CALL RDPRMC(ANS,NA,.TRUE.,
     &     'DO YOU WANT TO KEEP THE HEADER? (N/Y)', NULL,IRT)
        HEADER =  (ANS .EQ. 'Y' .OR. ANS .EQ. 'y') 
        IF (HEADER) THEN
          WRITE(NOUT,*)
     &         'SORRY: THIS HEADER OPTION NO LONGER AVAILABLE'
          HEADER = .FALSE.
        ENDIF

        IF (IMAMI .EQ. 0) THEN
C          MUST NORMALIZE INPUT IMAGE FIRST
           CALL NORM3(LUNO,NSAM,NROW,NSLICE,FMAX,FMIN,AV)
        ENDIF

        IPAD    = 0
        NSAMRE  = MOD(NSAM,4)
        CALL RDPRMC(ANS,NC,.TRUE.,
     &      'PAD TO INCREMENT OF 4 IF NECESSARY? (N/Y)',NULL,IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 9001

        IF ((NSAMRE .NE. 0).AND.(ANS .EQ. 'Y' .OR. ANS .EQ. 'y'))THEN
C          WANT TO PAD
           IPAD = 4 - NSAMRE
        ENDIF

C       OUTPUT FILE HAS 1 BYTE RECORDS, DIRECT FORMATTED RECORDS
        LENREC = -1
        CALL OPAUXFILE(.TRUE.,FILNEW,DATEXC,LUNN,LENREC,'N',
     &                  'EIGHT BIT RAW',.TRUE.,IRTFLG)

        CALL UNIXTOUNIX8(LUNO,LUNN,0,255,FMIN,FMAX,NSAM,NROW,NSLICE,
     &                   NORMAL,IPAD,IRTFLG)

C       CLOSE INPUT AND OUTPUT FILES
        CLOSE(LUNN)
9001    CLOSE(LUNO)

        RETURN
        END
