
C***********************************************************************
C
C SETVAL.F           REMOVED FROM UTIL1           OCT  88 ArDean Leith
C                    LONG FILE NAMES              FEB  89 ArDean Leith
C                    CAN OPEN STACK WITHOUT @ NOW SEPT 98 ArDean Leith
C                    OPFILEC                      FEB  03 ArDean Leith
C                    RDPRAF REMOVED               DEC  05 ArDean Leith 
C                    RDPRM3S                      AUG  13 ArDean Leith 
C                    MRC SUPPORT                  OCT  19 ArDean Leith 
C                    CORIGIN                      NOV  19 ArDean Leith 
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2019  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C  SETVAL(LUN1,NX1,NY1,NZ1)
C
C  PURPOSE:    SET PARTICULAR HEADER LABEL  LOCATIONS 
C
C  CALL TREE:  UTIL1 --> SETVAL --> SETLAB_R 
C                               --> SETLAB_R_MRC
C                               --> LUNSETHAND_MRC
C  FOR MRC FILES:
C     14-16 CELLB           CELL ANGLES IN DEGREES                          
C     20    DMIN            MINIMUM DENSITY VALUE
C     21    DMAX            MAXIMUM DENSITY VALUE  
C     22    DMEAN           MEAN    DENSITY VALUE  
C     25-26 EXTRA           EXTRA, USER DEFINED STORAGE SPACE
C     29-49 EXTRA           EXTRA, USER DEFINED STORAGE SPACE
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE SETVAL(LUN1,NX1,NY1,NZ1)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 

        INTEGER               :: LUN1,NX1,NY1,NZ1

        CHARACTER(LEN=MAXNAM) :: FILNAM
        CHARACTER(LEN=1)      :: SET,ANS

        INTEGER, PARAMETER    :: NMAX = 200
        INTEGER               :: ILIST(NMAX) 
        REAL                  :: FLIST(NMAX) 
        REAL                  :: FARRAY(7)
        LOGICAL               :: IS_MRC,IS_MRC2
        INTEGER               :: ILPOS,I,LENREC,NVAL2,NMAX2
        INTEGER               :: MAXIM,ITYPE,NVAL1,NANG,IDUM,IRTFLG
        INTEGER               :: NOT_USED,NCHAR,IGO

        REAL                  :: PHI,THETA,PSI
        INTEGER               :: NX2,NY2,NZ2

        REAL, PARAMETER       :: PI   = 3.14159
        CHARACTER(LEN=2)      :: CORIG
        CHARACTER(LEN=1)      :: CHANDED
        CHARACTER(LEN=4)      :: CAXIS
        CHARACTER(LEN=1)      :: NULL = CHAR(0)
        INTEGER, PARAMETER    :: LUN2 = 21

        CALL LUNGETIS_MRC(LUN1,IS_MRC,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

10      IF (IS_MRC) THEN
          CALL RDPRMC(SET,NCHAR,.TRUE.,
     &    '(A)NGLES, (L)OCS, (D)ATA-ORIGIN, (P)COPY, (C)LEAR, OR (F)IX',
     &     NULL,IRTFLG)
        ELSE

          CALL RDPRMC(SET,NCHAR,.TRUE.,
     &      '(A)NGLES, (L)OCS, (P)COPY, (C)LEAR, OR (F)IX',NULL,IRTFLG)
        ENDIF

        IF (IRTFLG .NE. 0) RETURN


        IF (SET == 'A') THEN
C          CHANGE THE TILT ANGLES ---------------------------------- A

           ILPOS = 1
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

           FARRAY(1) = 3.0
           FARRAY(2) = PHI
           FARRAY(3) = THETA
           FARRAY(4) = PSI

           IF (IS_MRC) THEN
C             MRC FILE
              IF (ILPOS == 1) THEN
                IGO = 42

                CALL SETLAB_R_MRC(LUN1,IGO,4,FARRAY,.FALSE.,IRTFLG)

              ELSEIF (ILPOS == 2) THEN 
                IGO       = 46
                FARRAY(1) = 6.0
                CALL SETLAB_R_MRC(LUN1,42, 1,FARRAY,   .FALSE.,IRTFLG)
                CALL SETLAB_R_MRC(LUN1,IGO,3,FARRAY(2),.FALSE.,IRTFLG)

              ELSE  
                CALL ERRT(101,'ONLY 2 ANGLE SETS IN MRC FILES',IDUM)
                IRTFLG = 1
              ENDIF
              GOTO 9000

           ELSE
C            SPIDER FILE
             IF (ILPOS == 1) THEN
                IGO       = 14
                NANG      = 4
                CALL SETLAB_R(LUN1,IGO,NANG,FARRAY,.FALSE.,IRTFLG)

             ELSEIF (ILPOS == 2) THEN
                IGO       = 30
                NANG      = 4
                CALL SETLAB_R(LUN1,IGO,NANG,FARRAY,.FALSE.,IRTFLG)

             ELSEIF (ILPOS == 3)THEN
                FARRAY(1) = 2
                IGO       = 30
                NANG      = 1
                CALL SETLAB_R(LUN1,IGO,NANG,FARRAY,.FALSE.,IRTFLG)
                IGO       = 34
                NANG      = 3
                CALL SETLAB_R(LUN1,IGO,NANG,FARRAY(2),.FALSE.,IRTFLG)

             ELSE
           
                CALL ERRT(102,'INCORRECT NUMBER FOR ANGLE SET',ILPOS)
                IRTFLG =1
             ENDIF
             GOTO 9000
          ENDIF

        ELSEIF (SET == 'D') THEN
C          CHANGE DATA ORIGIN & HANDEDNESS  ------------------------ D

           IF (.NOT. IS_MRC) THEN
              CALL ERRT(101,'OPERATION ONLY FOR MRC FILES',IDUM)
              IRTFLG = 1
              GOTO 9000
           ENDIF
 
           CORIG = 'LL'
           CALL RDPRMC(CORIG,NCHAR,.TRUE.,
     &      'ORIGIN  AT UPPER LEFT (UL) or LOWER LEFT(LL)', NULL,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

           CHANDED = 'L'
           CALL RDPRMC(CHANDED,NCHAR,.TRUE., 
     &       'HANDEDNESS (L or R)', NULL,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

C          SETS  AXIS ORIGIN LOCATION & VOLUME HANDEDNESS (CHARACTERS)
C          THIS IS NOT A STANDARD MRC DEFINED HEADER POSITION!!

           CAXIS = CORIG // ' ' // CHANDED
           CALL LUNSETHAND_MRC(LUN1,CAXIS,IRTFLG)
           !write(3,*) 'in setval - caxis:',caxis

C          WRITE ALTERED HEADER BACK IN THE MRC FILE
           CALL LUNWRTHED_MRC(LUN1,IRTFLG)
           GOTO 9000


        ELSEIF (SET == 'B' .OR. SET == 'L') THEN
C          CHANGE A PARTICULAR BUFFER LOCATION IN THE FILE LABEL ---- B

1601       NVAL1 = NMAX
           CALL RDPRAI(ILIST,NMAX,NVAL1,1,1024,
     &          'NUMBER(S) OF HEADER LOCATION TO BE CHANGED',
     &          NULL,IRTFLG)
           IF (IRTFLG == -1) GOTO 10

           NMAX2 = NVAL1
           CALL RDPRA('NEW VALUE FOR EACH HEADER LOCATION CHANGED',
     &               NMAX2,0,.FALSE.,FLIST,NVAL2,IRTFLG)
           IF (IRTFLG == -1) GOTO 1601

           IF (NVAL2 .NE. NVAL1) THEN
               CALL ERRT(102,'INCORRECT NUMBER OF VALUES',NVAL2)
               GOTO 9000
           ENDIF
              
           IF (IFORM == 8 .OR. IFORM == 11) THEN
C             FOR 8 BIT FILES
              LENREC = NX1 / 4
              IF ((LENREC * 4) < NX1) LENREC = LENREC + 1
           ELSE
C             FOR NORMAL 32 BIT SPIDER FILES
              LENREC = NX1
           ENDIF

           DO I = 1,NVAL1            ! MRC OK
             CALL SETLAB_R(LUN1,ILIST(I),1,FLIST(I),.FALSE.,IRTFLG)
           ENDDO


       ELSEIF (SET == 'P') THEN
C         COPY ANGLES FROM HEADER OF INPUT FILE TO  OUTPUT FILE ----- P
C         OUTPUT FILE MUST EXIST, IT IS NOT CREATED HERE.

          MAXIM = 0
          CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'O',ITYPE,
     &            NX2,NY2,NZ2,MAXIM,'OUTPUT',.TRUE.,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9000

          CALL LUNGETIS_MRC(LUN2,IS_MRC2,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9000

C          IF (IS_MRC .NE. IS_MRC2) THEN !gfort error,  Sept 2025 al
          IF (IS_MRC .NEQV. IS_MRC2) THEN 
             CALL ERRT(101,'INPUT/OUTPUT MUST BE SAME MRC/SPIDER',IDUM)
             IRTFLG =1
             GOTO 9000
          ENDIF

          IF (IS_MRC) THEN
C            MRC FILE
             IGO = 43      ! SET ONE
             CALL GETLAB_R_MRC(LUN1,IGO,3,FARRAY,IRTFLG)
             CALL SETLAB_R_MRC(LUN2,IGO,3,FARRAY,.FALSE.,IRTFLG)
             IF (IRTFLG .NE. 0) GOTO 9000
              
             IGO = 46      ! SET TWOO
             CALL GETLAB_R_MRC(LUN1,IGO,3,FARRAY,IRTFLG)
             CALL SETLAB_R_MRC(LUN2,IGO,3,FARRAY,.FALSE.,IRTFLG)
             IF (IRTFLG .NE. 0) GOTO 9000 
 
             FARRAY(1) = 6
             CALL SETLAB_R_MRC(LUN2,42,1,FARRAY,.FALSE.,IRTFLG)

         ELSE
C            SPIDER FILE
             IGO = 14     ! SET ONE
             CALL GETLAB(  LUN1,IGO,4,FARRAY,IRTFLG)
             CALL SETLAB_R(LUN2,IGO,4,FARRAY,.FALSE.,IRTFLG)
       
             IGO = 30     ! SET TWO & THREE
             CALL GETLAB(  LUN1,IGO,7,FARRAY,IRTFLG)
             CALL SETLAB_R(LUN2,IGO,7,FARRAY,.FALSE.,IRTFLG)
          
             IGO = 51     ! SET FOUR & FIVE
             CALL GETLAB(LUN1,  IGO,6,FARRAY,IRTFLG)
             CALL SETLAB_R(LUN2,IGO,6,FARRAY,.FALSE.,IRTFLG)
          ENDIF


        ELSEIF (SET == 'C') THEN
C          CLEAR FMIN... STATISTICS  ------------------------------- C

           SIG   = -1.0
           CALL SETPRMB(LUN1, 0.0,0.0, 0.0,SIG)    ! MRC OK


        ELSEIF (SET == 'F') THEN
C          SET FMIN, FMAX, AV, AND S.D. ---------------------------- F

           CALL RDPRM1S(FMAX,NOT_USED,'IMAGE MAXIMUM',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

           CALL RDPRM1S(FMIN,NOT_USED,'IMAGE MINIMUM',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

           CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &           'AVERAGE AND STANDARD DEVIATION AVAILABLE (Y/N)',
     &           NULL,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

           IF (ANS == 'N') THEN
               CALL ERRT(101,'MUST PROVIDE AVERAGE AND S.D. NOW',IDUM)
               GOTO 9000
           ENDIF

           CALL RDPRM2S(AV,SIG,NOT_USED,
     &                  'AVERAGE & STANDARD DEVIATION',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

           IF (IFORM < 0) THEN
C             FOURIER
              FMAX = FMAX * PI / 180.0
              FMIN = FMIN * PI / 180.0
           ENDIF

           CALL SETPRMB(LUN1, FMAX,FMIN, AV,SIG)    ! MRC OK

        ELSE
C          UNKNOWN OPTION ------------------------------------------- ?
           CALL ERRT(101,'UNKNOWN OPTION',IDUM)
           GOTO 9000
        ENDIF

9000    CONTINUE
        CLOSE(LUN1)
        CLOSE(LUN2)

        END
