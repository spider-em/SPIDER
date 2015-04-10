
C++*********************************************************************
C
C SETVAL.F           REMOVED FROM UTIL1           OCT  88 ARDEAN LEITH
C                    LONG FILE NAMES              FEB  89 ARDEAN LEITH
C                    CAN OPEN STACK WITHOUT @ NOW SEPT 98 ARDEAN LEITH
C                    OPFILEC                      FEB  03 ARDEAN LEITH
C                    RDPRAF REMOVED               DEC  05 ARDEAN LEITH 
C                    RDPRM3S                      AUG  13 ARDEAN LEITH 
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
C  SETVAL(LUN1)
C
C  PURPOSE:   SET PARTICULAR LABEL LOCATIONS TO SPECIFIED VALUES
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE SETVAL(LUN1,NX1,NY1,NZ1)

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 

        CHARACTER(LEN=MAXNAM) :: FILNAM
        CHARACTER(LEN=1)      :: SET,ANS

        INTEGER, PARAMETER    :: NMAX = 200
        INTEGER               :: ILIST(NMAX) 
        REAL                  :: FLIST(NMAX) 
        REAL                  :: FARRAY(4)

        REAL, PARAMETER       :: PI   = 3.14159
        CHARACTER(LEN=1)      :: NULL = CHAR(0)
        INTEGER, PARAMETER    :: LUN2 = 21


10      CALL RDPRMC(SET,NCHAR,.TRUE.,
     &    '(A)NGLES, (BUF), (P)COPY, (C)LEAR, OR (F)IX',NULL,IRTFLG)

        IF (SET == 'A') THEN
C          CHANGE THE TILT ANGLES ---------------------------------- A

           CALL RDPRI1S(ILPOS,NOT_USED,
     &                  'ANGLE SET 1, 2, OR 3',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          GET EULER ANGLES
           PHI   = 0.0
           THETA = 0.0
           PSI   = HUGE(PSI)

           CALL RDPRM3S(PHI,THETA,PSI,NOT_USED,
     &                 'PHI, THETA, & PSI',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (PSI == HUGE(PSI)) THEN
              PSI = 0.0
              CALL  RDPRM1S(PSI,NOT_USED,'PSI',IRTFLG)
              IF (IRTFLG .NE. 0) RETURN
           ENDIF

           FARRAY(1) = 1.0
           FARRAY(2) = PHI
           FARRAY(3) = THETA
           FARRAY(4) = PSI

           IF (ILPOS == 1) THEN
             LOCATION = 14
             NANG     = 4
             CALL SETLAB(LUN1,NX1,DUM,LOCATION,NANG,FARRAY,'U',IRTFLG)

           ELSEIF (ILPOS == 2)THEN
             LOCATION = 30
             NANG     = 4
             CALL SETLAB(LUN1,NX1,DUM,LOCATION,NANG,FARRAY,'U',IRTFLG)

           ELSEIF (ILPOS == 3)THEN
             FARRAY(1) = 2
             LOCATION  = 30
             NANG      = 1
             CALL SETLAB(LUN1,NX1,DUM,LOCATION,NANG,FARRAY,'U',IRTFLG)
             LOCATION  = 34
             NANG      = 3
             CALL SETLAB(LUN1,NX1,DUM,LOCATION,NANG,
     &                   FARRAY(2),'U',IRTFLG)
           ELSE
              CALL ERRT(102,'INCORRECT NUMBER FOR ANGLE SET',ILPOS)
              GOTO 9000
           ENDIF

        ELSEIF (SET == 'B') THEN
C          CHANGE A PARTICULAR BUFFER LOCATION IN THE FILE LABEL ---- B

1601       NVAL1 = NMAX
           CALL RDPRAI(ILIST,NMAX,NVAL1,1,1024,
     &          'NUMBER(S) OF HEADER LOCATION TO BE CHANGED',
     &          NULL,IRTFLG)
           IF (IRTFLG == -1) GOTO 10

           NMAX2 = NVAL1
           CALL RDPRA(
     &              'NEW VALUE FOR EACH HEADER LOCATION CHANGED',
     &               NMAX2,0,.FALSE.,FLIST,NVAL2,IRTFLG)
           IF (IRTFLG == -1) GOTO 1601

           IF (NVAL2 .NE. NVAL1) THEN
               CALL ERRT(102,'INCORRECT NUMBER OF VALUES',NVAL2)
               GOTO 9000
           ENDIF
              
           IF (IFORM == 8 .OR. IFORM == 11) THEN
C             FOR 8 BIT FILES
              LENREC = NX1 / 4
              IF ((LENREC * 4) .LT. NX1) LENREC = LENREC + 1
           ELSE
C             FOR NORMAL 32 BIT SPIDER FILES
              LENREC = NX1
           ENDIF

           DO I = 1,NVAL1
             CALL SETLAB(LUN1,LENREC,DUM,ILIST(I),1,FLIST(I),'U',IRTFLG)
           ENDDO

       ELSEIF (SET == 'P') THEN
C         COPY ANGLES FROM HEADER OF INPUT FILE TO  OUTPUT FILE ----- P
C         OUTPUT FILE PREEXISTS, IT IS NOT CREATED HERE.

          MAXIM = 0
          CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'O',ITYPE,
     &            NX2,NY2,NZ2,MAXIM,'OUTPUT',.TRUE.,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

          CALL COPYANGLES(LUN1,LUN2,NX1,NX2) 
                 
        ELSEIF (SET == 'C') THEN
C          CLEAR ---------------------------------------------------- C

           SIG   = -1.0
           CALL SETPRM(LUN1,NX1,NY1,0.0,0.0,0.0,'U')

        ELSEIF (SET == 'F') THEN
C          SET FMIN, FMAX, AV, AND S.D.

C          INPUT PHI AND THETA  TO BE STORED IN FMAX AND FMIN
           CALL RDPRM(FMAX,NOT_USED,'IMAGE MAXIMUM')

           CALL RDPRM(FMIN,NOT_USED,'IMAGE MINIMUM')

           CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &           'AVERAGE AND STANDARD DEVIATION AVAILABLE (Y/N)',
     &           NULL,IRTFLG)

           IF (ANS == 'N') THEN
               CALL ERRT(101,'MUST PROVIDE AVERAGE AND S.D. NOW',NE)
               GOTO 9000
           ENDIF

           CALL RDPRM2(AV,SIG,NOT_USED, 
     &               'AVERAGE & STANDARD DEVIATION')

           IF (IFORM < 0) THEN
C             FOURIER
              FMAX = FMAX * PI / 180.0
              FMIN = FMIN * PI / 180.0
           ENDIF
           CALL SETPRM(LUN1,NX1,NY1,FMAX,FMIN,AV,'U')

        ELSE
C          UNKNOWN OPTION ------------------------------------------- ?
           CALL ERRT(101,'UNKNOWN OPTION',NE)
           GOTO 9000
        ENDIF

9000    CONTINUE
        CLOSE(LUN1)
        CLOSE(LUN2)

        END
