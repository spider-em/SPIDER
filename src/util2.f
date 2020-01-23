C++*********************************************************************
C
C UTIL2.F                   NEW
C                           CHANGED:         Mahieddine Ladjadj 4/23/93
C                           MODIFIED              Jing Su       8/24/93
C                           REWRITTEN             ArDean Leith  8/30/96
C                           ADDED "PP S"          ArDean Leith  4/08/98
C                           STACKS IN "AD"        ArDean Leith 01/11/99
C                           REWROTE "AD"..        ArDean Leith 04/03/99
C                           USED RDPRM3S          ArDean Leith 08/05/99
C                           USED ALLOCATE         ArDean Leith 01/04/01
C                           'PA' & 'IN' 3D        ArDean Leith 03/01/01
C                           'CE VAR'              ArDean Leith 05/01/01
C                           'CE RAN'              ArDean Leith 05/02/01
C                           'CE MAX'              ArDean Leith 05/03/01
C                           'CE MIN'              ArDean Leith 05/03/01
C                           'CE LAP'              ArDean Leith 05/03/01
C                           'CE SOB'              ArDean Leith 05/03/01
C                           'CE PRE'              ArDean Leith 05/03/01
C                           'CE TOP'              ArDean Leith 05/04/01
C                           'CE RIDGE'            ArDean Leith 05/08/01
C                           'CE HURST'            ArDean Leith 05/08/01
C                           'CE HARALICK'         ArDean Leith 05/16/01
C                           NORM3 IN CE           ArDean Leith 04/02/02
C                           'CE LAHE'             ArDean Leith 04/10/02
C                           'CE AD'               ArDean Leith 04/18/02
C                           'CE OR'               ArDean Leith 04/18/02
C                           'AR SCA'              ArDean Leith 09/11/02
C                           'AR SCA' NORM3        ArDean Leith 10/04/02
C                           STACKS SUPPORT        ArDean Leith 10/04/02 
C                           'CE L' REMOVED        ArDean Leith 11/19/02 
C                           'WI' x,Y,Z            ArDean Leith 12/02/02  
C                           OPFILEC               ArDean Leith  3/18/03  
C                           'AD F'                ArDean Leith  3/24/03    
C                           'AD S'                ArDean Leith  4/21/03    
C                           'DIV'                 ArDean Leith  5/30/03    
C                           'SQRT'                ArDean Leith  5/30/03   
C                           RDPRM3S BUG           ArDean Leith  9/05/03   
C                           GPRP                  ArDean Leith  9/08/03
C                           MPI                   Chao Yang    10/31/03
C                           USEBORDER             ArDean Leith 11/21/03
C                           'PD' LOCATION         ArDean Leith  5/14/04
C                           MPI REMOVED           Chao Yang    11/19/04
C                           'CE WA'               ArDean Leith 11/19/04
C                           'CE ME'               ArDean Leith 06/22/05
C                           'BL' AV BUG           ArDean Leith 03/30/06
C                           'WI' 1 SLICE BUG      ArDean Leith 10/19/06
C                           NSLICEW = -999999999  ArDean Leith 03/05/07
C                           SETPRMB PARAMETERS    ArDean Leith 05/19/09
C                           'IP' SETPRMS          ArDean Leith 11/23/10
C                           REMOVED 'NK'          ArDean Leith 11/24/10
C                           'IP' CUTOUT           ArDean Leith 11/23/10
C                           'IP FT'               ArDean Leith  5/23/11
C                           'IP FP'               ArDean Leith  5/23/11
C                           'MM' FOURIER OK       ArDean Leith  7/23/11
C                           'PD' PROMPT           ArDean Leith 11/03/11
C                           'AS AV' PROMPT        ArDean Leith 03/13/12
C                           'AD L', 'MU L'...     ArDean Leith 11/01/12
C                           'MA SOFT'             ArDean Leith 04/25/13
C                           'IN S' BUG            ArDean Leith 01/10/14
C                           'CE' UNCOCUMENTED     ArDean Leith 02/10/14
C                           'MA' FMIN             ArDean Leith 12/10/14
C                           'DIS' ADDED           ArDean Leith  1/05/15
C                           PUTLIN PARAMETERS     ArDean Leith  6/13/18
C                           COSMETIC              ArDean Leith 10/10/19
C                           'MD MRC'              ArDean Leith 11/20/19
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
C  UTIL2(MAXDIM)
C
C  PARAMETERS:      MAXDIM         MAX LENGTH FOR UNLABELED COMMON
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE UTIL2(MAXDIM)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        INTEGER                    :: MAXDIM

        COMMON BUF(1)  ! LEGACY SIZE SHOULD BE UPDATED!! al

        CHARACTER(LEN=MAXNAM)      :: FILNAM1,FILNAM2,FILNAM3,FILNAM

        INTEGER, PARAMETER         :: NFUNC=33
        CHARACTER(LEN=2)           :: FUNC(NFUNC)

        REAL                       :: FWA(3)
        REAL, ALLOCATABLE          :: Q(:)
        INTEGER                    :: IORDER(3),ISIZE(3)
        INTEGER                    :: IRTFLG
        CHARACTER(LEN=MAXNAM)      :: EXPR
        CHARACTER(LEN=1)           :: NULL = CHAR(0)
        CHARACTER(LEN=1)           :: MODE
        CHARACTER(LEN=2)           :: ANS
        LOGICAL                    :: NORMIT,USEBORDER

        INTEGER, PARAMETER         :: LUN1  = 21
        INTEGER, PARAMETER         :: LUN2  = 22
        INTEGER, PARAMETER         :: LUN3  = 23
        INTEGER, PARAMETER         :: LUN4  = 70


        DATA FUNC/'AD', 'BL', 'CP', 'IN', 'IP',
     &            'MU', 'PA', 'PD', 'SH', 'SQ',
     &            'SU', 'WI', 'CE', 'AR', 'MR',
     &            'DF', 'MA', 'WV', 'PP', 'SZ', 
     &            'WU', 'MM', 'CM', 'PV', 'NK', 
     &            'AS', 'MN', 'TH', 'GP', 'RP',
     &            'MX', '30', '31'/ 

        MAXIM  = 0
        MAXIM2 = 0
        IRTFLG = 0

        IF (FCHAR(1:4) == 'SQRT') GOTO 21 ! SQRT
        IF (FCHAR(1:2) == '12')   GOTO 6  ! DIV
        IF (FCHAR(1:2) == '30')   GOTO 66 ! DISP
        IF (FCHAR(1:2) == '31')   GOTO 67 ! MRC

        DO  IFUNC = 1, NFUNC
          IF (FCHAR(1:2) == FUNC(IFUNC)) THEN
            GOTO ( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
     &            11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
     &            21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31 ) IFUNC
          ENDIF
        ENDDO
C       FUNCTION NOT FOUND IN UTIL2
        RETURN 


C       OPERATION ----------- DISP = 30 ------------------------ 'DISP' 
66      CALL DISP()              
        GOTO 9000

C       OPERATION ----------- DISP = 31 ------------------------ 'MRC' 
67      IF (FCHAR(4:4) == 'H')  THEN
C          INQUIRE MRC HEADER CONTENTS
           CALL LISTHEDMRC(LUN1,IRTFLG)

        ELSE IF (FCHAR(4:4) == 'M') THEN
C          HANDLE 'MRC MD'
           CALL SETMODE_MRC()
        ENDIF
        GOTO 9000

C       OPERATION  ADD ------------------------------------------- 'AD' 
C       AD         ADD IMAGES           

1       IF (FCHAR(4:4) == 'S') THEN
C          ADD SERIES OF IMAGES, FASTER WITH LESS MEMORY ALLOCATED
           CALL ADS(LUN1,LUN2,LUN3)
           
        ELSE
           SIGN = +1.0
           IF ( INDEX(FCHAR(4:) ,'F') > 0) SIGN = 1000.0
           IF ( INDEX(FCHAR(4:) ,'R') > 0) SIGN = 2000.0

           IF ( INDEX(FCHAR(4:) ,'2') > 0) THEN
C             IMAGE LIST SUPPORTED  'AD 2', 'AD 2F', 'AD 2R'
              CALL UTIL2SUPL('FIRST  INPUT FILE NAME OR TEMPLATE~',
     &                       'SECOND INPUT FILE NAME OR TEMPLATE~',
     &                       'OUTPUT FILE NAME OR TEMPLATE~',
     &                        SIGN)
           ELSE
              CALL UTIL2SUP('FIRST INPUT',
     &                      'NEXT INPUT','SUMMED OUTPUT',
     &                      LUN1,LUN2,LUN3, SIGN)
           ENDIF
        ENDIF
        GOTO 9000

C       OPERATION ------   BLANK --------------------------------- 'BL' 

2       IFORM   = 3
        NSAM2   = 0
        NSLICE2 = -9999
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'U',IFORM,
     &               NSAM2,NROW2,NSLICE2,
     &               MAXIM,'BLANK OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        CALL RDPRMC(ANS,NA,.TRUE.,'AVERAGE? (Y/N)',NULL,IRTFLG)
        IF (ANS(:1) == 'Y') THEN
           CALL OPFILEC(0,.TRUE.,FILNAM2,LUN1,'O',IDUM,
     &                  NSAMR,NROWR,NDUM,
     &                 MAXIM2,'REFERENCE',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000
           IF (IMAMI .NE. 1) THEN
              CALL NORM3(LUN1,NSAMR,NROWR,NDUM,FMAXR,FMINR,AVR)
              B = AVR
           ELSE
              B = AV
           ENDIF
           CLOSE(LUN1)
        ELSE
C           INPUT BACKGROUND VALUE
            CALL RDPRM(B,NOT_USED,'BACKGROUND')
        ENDIF
        CALL BLANK(LUN2,NSAM2,NROW2,NSLICE2,B)
        GOTO 9000

C       OPERATION ----------- COPY  ------------------------------ 'CP' 
3       CALL COPY1(MAXDIM)              
        GOTO 9000

C       OPERATION -------- INSERT -------------------------------- 'IN' 
C               IN      : INSERT
C               IN S    : INSERT AND CONTRAST STRETCH

4       CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'SMALL INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000
        AV1    = AV
        IF (FCHAR(4:4) == 'S' .AND. IMAMI .NE. 1) THEN
           CALL NORM3(LUN1,NSAM1,NROW1,NSLICE1,FMAX1,FMIN1,AV1)
        ELSE
           FMIN1 = FMIN
           FMAX1 = FMAX
        ENDIF

        CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'O',IFORM,
     &               NSAM2,NROW2,NSLICE2,
     &               MAXIM2,'LARGE INPUT (OVERWRITTEN!)',
     &              .TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        IF (FCHAR(4:4) == 'S' .AND. IMAMI .NE. 1) THEN
           CALL NORM3(LUN2,NSAM2,NROW2,NSLICE2,FMAX2,FMIN2,AV2)
        ENDIF

        !write(6,*) ' fmin, fmin1: ', fmin,fmin1
        !write(6,*) ' fmax, fmax1: ', fmax,fmax1

        IF (IFORM < 0) THEN
C          WRONG LARGE INPUT FORMAT
           IER = 2
           GOTO 9900
        ENDIF

        NSAMS   = 1
        NROWS   = 1
        NSLICES = 1
        CALL RDPRI3S(NSAMS,NROWS,NSLICES,NOT_USED,
     &     'TOP LEFT COORDINATES',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        IN = 1
        CALL PATCH(LUN1,LUN2,NSAM1,NROW1,NSLICE1,
     &             NSAM2,NROW2,NSLICE2,
     &             NSAMS,NROWS,NSLICES, 
     &             IN,AV1,FCHAR(4:4),FMIN1,FMAX1,
     &             .FALSE.)

C       SET UNDETERMINED SATATISTICS FLAG
        CALL SETPRMB(LUN2, 0.0,0.0, 0.0,0.0)
        GOTO 9000


C       OPERATION  INTERPOLATE  ---------------------------------- 'IP'
C       IP         BILINEAR INTERPOLATION 
C       IP T       TRIANGULAR INTERPOLATION
C       IP FT      FOURIER INTERPOLATION
C       IP FP      FOURIER BASED POLYNOMIAL INTERPOLATION

5       CALL INTERPS(MAXDIM)
        GOTO 9000


C       OPERATION  MULTIPLY -------------------------------------- 'MU' 
C       MU            MULTIPLY REAL OR COMPLEX FILES 
C       MU D OR DIV   DIVIDE REAL FILES
C       MU M          MULTIPLY FIRST COMPLEX FILE BY THE SECOND CONJUGATED.
C       MU O          MULTIPLY WITH ARITHMETIC OR

6       SIGN = +2.0
        IF ( INDEX(FCHAR(4:) ,'D') > 0 .OR. 
     &            FCHAR(1:2) == '12')   THEN
C          'MU D' OR 'DIV'   DIVIDE REAL FILES
           IF ( INDEX(FCHAR(4:) ,'2') > 0) THEN
C             IMAGE LIST SUPPORTED
              CALL UTIL2SUPL('INPUT FILE NAME OR TEMPLATE~',
     &                       'DIVISOR FILE NAME OR TEMPLATE~',
     &                       'OUTPUT FILE NAME OR TEMPLATE~',
     &                       SIGN)
           ELSE
              CALL UTIL2SUP('INPUT','DIVISOR','OUTPUT',
     &                       LUN1,LUN2,LUN3, SIGN)
           ENDIF

        ELSEIF (INDEX (FCHAR(4:),'O') > 0) THEN
C          'MU O', 'MU 2O'          ARITHMETIC OR
           IF ( INDEX(FCHAR(4:) ,'2') > 0) THEN
              CALL UTIL2SUPL('INPUT FILE NAME OR TEMPLATE~',
     &                      'SECOND INPUT FILE NAME OR TEMPLATE~',
     &                      'OUTPUT FILE NAME OR TEMPLATE~',
     &                      SIGN)
           ELSE
              CALL UTIL2SUP('INPUT','SECOND INPUT','OUTPUT',
     &               LUN1,LUN2,LUN3, SIGN)
           ENDIF
        ELSE
C          MULTIPLICATION
           IF ( INDEX(FCHAR(4:) ,'2') > 0) THEN
C             IMAGE LIST SUPPORTED
              CALL UTIL2SUPL('INPUT FILE NAME OR TEMPLATE~',
     &                      'MULTIPLIER FILE NAME OR TEMPLATE~',
     &                      'OUTPUT FILE NAME OR TEMPLATE~',
     &                       SIGN)
           ELSE
              CALL UTIL2SUP('INPUT','MULTIPLIER','OUTPUT',
     &                       LUN1,LUN2,LUN3, SIGN)
           ENDIF
        ENDIF
        RETURN

C       OPERATION PATCH (ADDS INPUT TO ORIG) --------------------- 'PA' 

C       OPEN FIRST INPUT FILE
7       CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM1,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'SMALL INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000
        FMIN1 = FMIN
        FMAX1 = FMAX

C       OPEN SECOND INPUT FILE
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'O',IFORM2,
     &             NSAM2,NROW2,NSLICE2,
     &             MAXIM2,'LARGE INPUT (OVERWRITTEN!)',.TRUE.,IRTFLG)
        IF (IRTFLG == -1) GOTO 7
        IF (IRTFLG .NE. 0) GOTO 9000

        IF (IFORM <= 0) THEN
           IER = 2
           GOTO 9900
        ENDIF 
   
        NSAMS   = 1
        NROWS   = 1
        NSLICES = 1
        CALL RDPRI3S(NSAMS,NROWS,NSLICES,NOT_USED,
     &     'TOP LEFT COORDINATES',IRTFLG)

        CALL PATCH(LUN1,LUN2,NSAM1,NROW1,NSLICE1, NSAM2,NROW2,
     &             NSLICE2, NSAMS,NROWS,NSLICES, 0,0,FCHAR(4:4),
     &             FMIN1,FMAX1,.FALSE.)

C       SET UNDETERMINED STATISTICS FLAG
        CALL SETPRMB(LUN2, 0.0,0.0, 0.0,0.0)
        GOTO 9000                            


C       OPERATION   PAD  ----------------------------------------- 'PD' 

C       EMBED A PICT. OR VOL. IN A LARGER EMPTY ARRAY.

C       OPEN INPUT FILE
8       CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &              MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        AV1    = AV
        IMAMI1 = IMAMI
        FMIN1  = FMIN
        FMAX1  = FMAX

C       OPEN THE OUTPUT FILE
        NSAM2   = 0
        NSLICE2 = 0
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &               NSAM2,NROW2,NSLICE2,
     &               MAXIM2,'OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9900

        IN = 2
        CALL RDPRMC(ANS,NC,.TRUE.,
     &  'PAD TYPE= AVERAGE, SET, BORDER, MIN, OR CIRCULAR (A/S/B/M/C)',
     &       NULL,IRTFLG)
       !&   'AVERAGE? (Y/N), (B)ORDER,(M)INIMUM,(C)IRCULAR'

        IF (ANS(1:1) == 'A') THEN  ! OBSOLETE INPUT RESPONSES
           ANS(1:1) = 'Y'
        ELSEIF (ANS(1:1) == 'S') THEN
           ANS(1:1) = 'N'
        ENDIF

        IF (NC >= 2 .AND. ANS(2:2).EQ. 'C') IN = 3
        IF (ANS(:1) == 'Y' .OR. ANS(:1) == 'M') THEN
           IF (IMAMI1 .NE. 1) THEN
              CALL NORM3(LUN1,NSAM1,NROW1,NSLICE1,FMAX1,FMIN1,AV1)
           ENDIF

           IF (ANS(:1) == 'M') THEN
             B = FMIN1
           ELSE
             B = AV1
           ENDIF
        ELSE
           IF (ANS(:1) .NE. 'B') THEN
              CALL RDPRM(B,NOT_USED,'BACKGROUND')
           ENDIF
        ENDIF

        NSAMS  = 1
        NROWS  = 1
        IF (NSLICE2 <= 1) THEN
C          PAD INTO A INTO AN IMAGE
           NSLICS = 1
           CALL RDPRIS(NSAMS,NROWS,NOT_USED,
     &              'TOP LEFT COORDINATES',IRTFLG)
        ELSE
C          PAD INTO A VOLUME
           NSLICS = -1000000
           CALL RDPRI3S(NSAMS,NROWS,NSLICS,NOT_USED,
     &              'TOP LEFT COORDINATES',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

           IF (NSLICS <= -1000000) THEN
C             NEED TO INQUIRE AS TO NSLICS
              NSLICS = 1
              CALL RDPRI1S(NSLICS,NOT_USED,
     &                    'TOP Z COORDINATE',IRTFLG)
           ENDIF
        ENDIF
        IF (IRTFLG .NE. 0) GOTO 9000

C       A 'B' IS A SIGNAL FOR BORDERING PATCH
        USEBORDER = (ANS(:1) == 'B')

        CALL PATCH(LUN1,LUN2,NSAM1,NROW1,NSLICE1,
     &             NSAM2,NROW2,NSLICE2,
     &             NSAMS,NROWS,NSLICS, IN,B,FCHAR(4:4),
     &             FMIN1,FMAX1,USEBORDER)

C       SET UNDETERMINED STATISTICS FLAG
        CALL SETPRMB(LUN2, 0.0,0.0, 0.0,0.0)
        GOTO 9000

C       OPERATION  ------------------------------------------------ 'SH' 
C               SH       SHIFT
C               SH F     SHIFT USING FOURIER INTERPOLATION

C       OPEN FIRST INPUT FILE
9       CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NSAM1,NROW1,NSLICE1,MAXIM,
     &               'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

C       OPEN OUTPUT FILE
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &               NSAM1,NROW1,NSLICE1,MAXIM2,
     &               'OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        SAMS  = 0
        ROWS  = 0
        SLICS = HUGE(SLICS)

        IF (IFORM == 3) THEN
C          SHIFT IN 3-D
           CALL RDPRM3S(SAMS,ROWS,SLICS,NOT_USED,
     &              'SHIFT COMPONENTS IN X, Y, & Z',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000
           TVAL = HUGE(SLICS)
           IF (SLICS == TVAL) THEN
              CALL RDPRM1S(SLICS,NOT_USED,
     &                    'SHIFT COMPONENT IN Z',IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9000
           ENDIF

           IF (SAMS  == FLOAT(IFIX(SAMS)) .AND.
     &         ROWS  == FLOAT(IFIX(ROWS)) .AND.
     &         SLICS == FLOAT(IFIX(SLICS)))  THEN

C              INTEGER SHIFT
               NSAMS  = SAMS
               NROWS  = ROWS
               NSLICS = SLICS
               IF (2*NSAM1 > MAXDIM) THEN
                  IER = 6
                  GOTO 9900
               ENDIF

               CALL SHIFT3(LUN1,LUN2,NSAM1,NROW1,NSLICE1,
     &                     NSAMS,NROWS,NSLICS)
          ELSE
            IF (FCHAR(4:5) == 'F')  THEN
               NNNN    = NSAM1+2-MOD(NSAM1,2)
               MEMWANT = NNNN*NROW1*NSLICE1

               ALLOCATE(Q(MEMWANT),STAT=IRTFLG)
               IF (IRTFLG .NE. 0)  THEN
                  IER = 6
                  GOTO 9900
               ENDIF

               DO J = 1, NROW1*NSLICE1
                  CALL REDLIN(LUN1,Q(1 + (J-1)*NNNN),NSAM1,J)
               ENDDO
               INS = +1
               CALL  FMRS_3(Q(1),NSAM1,NROW1,NSLICE1,INS)
               IF (INS == 0)  THEN
                  IER = 38
                  GOTO 9900
               ENDIF
               CALL SHIFT_PF(Q(1),NNNN/2,NSAM1,NROW1,NSLICE1,
     &                       SAMS,ROWS,SLICS)
               INS=-1
               CALL FMRS_3(Q(1),NSAM1,NROW1,NSLICE1,INS)
               DO J = 1, NROW1*NSLICE1
                  CALL WRTLIN(LUN2,Q(1 + (J-1)*NNNN),NSAM1,J)
               ENDDO

               IF (ALLOCATED(Q)) DEALLOCATE(Q)
            ELSE
                IF (5*NSAM1 > MAXDIM) THEN
                IER = 6
                GOTO 9900
             ENDIF
            CALL  SHIFT_3D(LUN1,LUN2,BUF(1),BUF(4*NSAM1+1),
     &                  NSAM1,NROW1,NSLICE1,SAMS,ROWS,SLICS)
            ENDIF
          ENDIF

        ELSEIF (IFORM == 1)  THEN
           CALL RDPRM2S(SAMS,ROWS,NOT_USED,
     &                 'SHIFT COMPONENTS IN X & Y',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

         IF (REAL(IFIX(SAMS)).EQ.SAMS .AND.
     &       REAL(IFIX(ROWS)).EQ.ROWS) THEN
C          INTEGER SHIFT
           NSAMS = SAMS
           NROWS = ROWS
           IF (2*NSAM1 > MAXDIM) THEN
              IER = 6
              GOTO 9900
           ENDIF
           CALL  SHIFT2(LUN1,LUN2,NSAM1,NROW1,NSAMS,NROWS)

         ELSE
            IF (FCHAR(4:5).EQ.'F')  THEN
               NNNN    = NSAM1+2-MOD(NSAM1,2)
               MEMWANT = NNNN*NROW1
               ALLOCATE(Q(MEMWANT),STAT=IRTFLG)

               IF (IRTFLG .NE. 0)  THEN
                  IER = 6
                  GOTO 9900
               ENDIF

               DO  J = 1, NROW1
                  CALL REDLIN(LUN1,Q(1 + (J-1)*NNNN),NSAM1,J)
               ENDDO
               INS = +1
               CALL  FMRS_2(Q(1),NSAM1,NROW1,INS)
               IF (INS == 0)  THEN
                  IER = 38
                  GOTO 9900
               ENDIF
               NSLICE1 = 1
               SLICS   = 0.0
               CALL SHIFT_PF(Q(1),NNNN/2,NSAM1,NROW1,NSLICE1,
     &                        SAMS,ROWS,SLICS)
               INS=-1
               CALL FMRS_2(Q(1),NSAM1,NROW1,INS)
               DO J = 1, NROW1
                   CALL WRTLIN(LUN2,Q(1 + (J-1)*NNNN),NSAM1,J)
               ENDDO

               IF (ALLOCATED(Q)) DEALLOCATE(Q)

            ELSE
C              BILINEAR INTERPOLATION
               IF (6*NSAM1 > MAXDIM) THEN
                  IER = 6
                  GOTO 9900
               ENDIF
               NNROWS = 1
               NNROWE = NROW1
               NNROWK = 1
               CALL SHIFTR(LUN1,LUN2,NSAM1,NROW1,NNROWS,NNROWE,
     &                  NNROWK,SAMS,ROWS)
            ENDIF
          ENDIF
        ELSE
           CALL ERRT(2,'SH',NE)
        ENDIF   
        GOTO 9000

C       OPERATION  ----------------------------------------------- 'SQ'
C       SQUARES INPUT POINT BY POINT 

10      CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

C       OPEN OUTPUT FILE (SAME SIZE AS INPUT)
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN3,'U',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM2,'OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        CALL IMSQ(.FALSE.,LUN1,LUN3,ITYPE,
     &             NSAM1,NROW1,NSLICE1,IRTFLG)
        GOTO 9000

C       OPERATION ------------------------------------------------- 'SU' 
C       SUBTRACTS ONE OR MORE IMAGE FROM THE FIRST IMAGE

11      SIGN = -1.0
        IF ( INDEX(FCHAR(4:) ,'2') > 0) THEN
C          IMAGE LIST SUPPORTED
           CALL UTIL2SUPL(
     &            'INPUT FILE NAME OR TEMPLATE (E.G. STK@****)~',
     &            'SUBTRACTED FILE NAME OR TEMPLATE~',
     &            'OUTPUT FILE NAME OR TEMPLATE~',
     &            SIGN)
        ELSE
           CALL UTIL2SUP('INPUT','SUBTRACTED','OUTPUT',
     &                   LUN1,LUN2,LUN3, SIGN)
        ENDIF
        RETURN


C       OPERATION  ------------------------------------------------ 'WI' 
C       WI        WINDOW 
C       WI B      WINDOW USING SPECIFIED BACKGROUND 

C       OPEN INPUT FILE
12      CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &              MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

C       OPEN OUTPUT FILE
        NSAM2   = 0
        NSLICE2 = 0
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &              NSAM2,NROW2,NSLICE2,
     &              MAXIM2,'OUTPUT',.FALSE.,IRTFLG)

        NSAMW   = 1
        NROWW   = 1
        NSLICEW = -999999999    ! SOMEDAY THIS WILL CAUSE A PROBLEM,
C                               ! BUT IT WILL BE AFTER I'M DEAD
        CALL RDPRI3S(NSAMW,NROWW,NSLICEW,NOT_USED,
     &              'TOP LEFT COORDINATES',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        IF (NSLICE1 > 1 .AND. NSLICEW == -999999999) THEN
C          FOR LEGACY INPUT OF X, Y ONLY, THEN Z ON NEXT LINE 
           NSLICEW = 1
           CALL RDPRI1S(NSLICEW,NOT_USED,'TOP SLICE',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000
        ELSEIF (NSLICE1 <= 1) THEN
           NSLICEW = 1
        ENDIF

        IF ((NSAMW               > NSAM1)   .OR.
     &      (NROWW               > NROW1)   .OR.
     &      (NSLICEW             > NSLICE1) .OR.
     &      ((NSAMW + NSAM2)     <= 1) .OR. 
     &      ((NROWW + NROW2)     <= 1) .OR. 
     &      ((NSLICEW + NSLICE2) <= 1)) THEN
            WRITE(NOUT,*)' WARNING: NO INPUT PIXELS WITHIN OUTPUT IMAGE'
        ENDIF

        BACK = 0.0
        IF (FCHAR(4:4) == 'B')
     &     CALL RDPRM(BACK,NOT_USED,'BACKGROUND')

        CALL WINDOW(LUN1,LUN2, NSAM1,NROW1,NSLICE1,
     &                         NSAMW,NROWW,NSLICEW, 
     &                         NSAM2,NROW2,NSLICE2,BACK)
        GOTO 9000


C       OPERATION ------------------------------------------------ 'CE' 
C       CE        CONTRAST ENHANCEMENT 
C       CE FIT       FIT HISTOGRAM
C       CE OD        FIT HISTOGRAM FOR OD MICROGRAPHS.
C       CE GNC       USING GRADUATED NON-CONVEX RESTORATION 
C       CE MED       USING MEDIAN FILTERING 
C       CE VAR       USING VARIANCE FILTERING 
C       CE VS        USING VARIANCE SMOOTHING FILTERING 
C       CE G?        USING GRADIENT FILTER 
C       CE RAN       USING RANGE FILTER 
C       CE MAX       USING MAX FILTER 
C       CE MIN       USING MIN FILTER 
C       CE LAP       USING LAPLACIAN FILTER 
C       CE SOBEL     USING SOBEL FILTER 
C       CE PREWITT   USING PREWITT FILTER 
C       CE RIDGE     RIDGE FOLLOWER 
C       CE HURST     USING HURST FILTER 
C       CE HARALICK  USING HARALICK FILTER 
C       CE LAHE      USING LOCAL AREA HISTOGRAM FILTER 
C       CE AD        USING ANISOTROPIC DIFFUSION FILTER 
C       CE OR        USING LOCAL ORIENTATION 
C       CE ME        USING MAXIMUM ENTROPY THRESHOLD 
 
13      IF (FCHAR(4:6) == 'FIT')  THEN
C          FITS HISTOGRAM OF IMAGE FILE TO THE HISTOGRAM OF REF. FILE.
           CALL  HISTE(MAXDIM)
           GOTO  9911

        ELSEIF (FCHAR(4:5) == 'OD')  THEN
C          ADJUST OPTICAL DENSITIES
           CALL  HISTOD
           GOTO  9911

        ENDIF

C       OPEN INPUT FILE, SOME OPERATIONS CAN TAKE WHOLE STACKS
        IF (FCHAR(4:5) == 'AD') MAXIM = -1
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        IF (IMAMI.NE.1) 
     &      CALL NORM3(LUN1,NSAM1,NROW1,NSLICE1,FMAX,FMIN,AV)
        FMIN1 = FMIN
        FMAX1 = FMAX
        SIG1  = SIG

C       OPEN OUTPUT FILE
        IF (MAXIM > 0) MAXIM2 = MAXIM
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &              NSAM1,NROW1,NSLICE1,
     &               MAXIM2,'OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        IF (FCHAR(4:6) == 'GNC') THEN
C          GRADUATED NON CONVEX RESTORATION
           CALL GNC(LUN1,LUN2,NSAM1,NROW1)

        ELSEIF (FCHAR(4:5) == 'AD')  THEN
C          ANISO DIFFUSION (CAN HANDLE WHOLE STACKS) 
           CALL ANISO(LUN1,LUN2,NSAM1,NROW1,NSLICE1,MAXIM,IRTFLG)

        ELSEIF (FCHAR(4:5) == 'OR')  THEN

C          OPEN SECOND OUTPUT FILE
           CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN3,'U',ITYPE,
     &             NSAM1,NROW1,NSLICE1,
     &             MAXIM2,'CONFIDENCE OUTPUT',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

           CALL ORIENT(LUN1,LUN2,LUN3,NSAM1,NROW1,NSLICE1,IRTFLG)

        ELSEIF (FCHAR(4:6) == 'MED') THEN
C          MEDIAN FILTER
           CALL MEDIAN(LUN1,LUN2,NSAM1,NROW1,NSLICE1)

        ELSEIF (FCHAR(4:5) == 'HA') THEN
C          HARALICK TEXTURE FILTER
           CALL FILTER_HAR(LUN1,LUN2,NSAM1,NROW1,NSLICE1,
     &                     FMIN1,FMAX1)

        ELSEIF (FCHAR(4:6) == 'LAH' ) THEN
C          LAHE --> 'LH'
           CALL FILTER(LUN1,LUN2,NSAM1,NROW1,NSLICE1,
     &                 MAXIM,'LH',FMIN1,FMAX1,SIG1)

        ELSEIF (FCHAR(4:5) == 'ST') THEN
C          CONTRAST STRETCHING
           CALL ENHANC(FILNAM2,LUN1,LUN2,NSAM1,NROW1,NSLICE1)

        ELSEIF (FCHAR(4:5) == 'WA')  THEN
C          WATERSHED 
           CALL WATERSHED(LUN1,LUN2,NSAM1,NROW1,NSLICE1,FMIN1)

        ELSEIF (FCHAR(4:5) == 'HI') THEN
C          HISTOGRAM EQUALIZATION
           CALL EHIST(FILNAM,LUN1,LUN2,NSAM1,NROW1,NSLICE1)

        ELSEIF (FCHAR(4:6) == 'MET') THEN
C          MAX. ENTROPY THRESHOLDING
           CALL MEHIST(LUN1,LUN2,NSAM1,NROW1,NSLICE1,FMIN1,FMAX1)

        ELSEIF (FCHAR(4:4) == 'V' .OR. FCHAR(4:4) == 'G' .OR.
     &          FCHAR(4:4) == 'R' .OR. FCHAR(4:4) == 'M' .OR.
     &          FCHAR(4:4) == 'L' .OR. FCHAR(4:4) == 'S' .OR.
     &          FCHAR(4:4) == 'P' .OR. FCHAR(4:4) == 'T' .OR.
     &          FCHAR(4:4) == 'F' .OR. FCHAR(4:5) == 'HU') THEN
C          VARIANCE, GRADIENT, RANGE, MAX, MIN, LAPLACIAN, SOBEL, 
C          PREWITT, TOP-HAT, FREI-CHEN, RIDGE, HURST
           CALL FILTER(LUN1,LUN2,NSAM1,NROW1,NSLICE1,
     &                 MAXIM,FCHAR(4:5),
     &                 FMIN1,FMAX1,SIG1)

        ELSE
C          HANDLE OLD PLAIN 'CE' (SHOULD NOT BE USED NOW)
           CALL RDPRMC(MODE,NC,.TRUE.,
     &        'STRETCH, HISTOGRAM EQUALIZE, OR LOCAL? (S/H/L)', 
     &        NULL,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

           IF (MODE == 'H') THEN
C             SAME AS 'CE HI'  HISTOGRAM EQUALIZATION
              WRITE(NOUT,*) ' OBSOLETE OPERATION, USE: <CE HI> INSTEAD'
              CALL EHIST(FILNAM,LUN1,LUN2,NSAM1,NROW1,NSLICE1)

           ELSEIF (MODE == 'L') THEN
C             A TILED LOCAL HISTOGRAM (A POOR ENHANCEMENT IDEA!!)
              CALL LOCAL(LUN1,LUN2,NSAM1,NROW1,NSLICE1)

           ELSE
C             SAME AS 'CE ST' CONTRAST STRETCH
              WRITE(NOUT,*) ' OBSOLETE OPERATION, USE: <CE ST> INSTEAD'
              CALL ENHANC(FILNAM2,LUN1,LUN2,NSAM1,NROW1,NSLICE1)
           ENDIF
        ENDIF
        GOTO 9000


C       OPERATION  ----------------------------------------------  'AR' 
C       AR        ARITHMETIC OPERATION

14      CONTINUE

C       OPEN INPUT FILE,  STACK OK FOR PLAIN 'AR'
        IF (FCHAR(4:5) .NE. 'IF') MAXIM = -1
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

C       OPEN OUTPUT FILE
        IF (MAXIM > 0) MAXIM2 = MAXIM
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM2,'OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        IF (FCHAR(4:5) == 'IF') THEN
           CALL ARITHL(LUN1,LUN2,NSAM1,NROW1,NSLICE1)
           GOTO 9000

        ELSE IF (FCHAR(4:5) == 'SC') THEN
C          'AR SCA'
           CALL RDPRM2S(FLOW,FHI,NOT_USED,
     &                  'NEW IMAGE MIN. & MAX.',IRTFLG)
        ELSE
C          PLAIN 'AR' 
           IRTFLG = -999     ! NO UPPERCASE
           CALL RDPRMC(EXPR,NLET,.TRUE.,
     &                 'FORMULA: P2=',NULL,IRTFLG)
        ENDIF
        IF (IRTFLG .NE. 0) GOTO 9000

        NORMIT = (FCHAR(4:5) == 'SC')
        IMGNUM = -3   

        !write(6,*)' In util2 -imgnum,nslice,maxim:',imgnum,nslice1,maxim

        DO WHILE (IMGNUM < MAXIM) 
           !write(6,*) ' In util2 - imgnum,maxim 2:',imgnum,maxim
           CALL GETSTACK(LUN1,LUN2,IMGNUM,MAXIM,
     &                   VERBOSE,.FALSE.,FDUM,NORMIT,IRTFLG)
           !write(6,*) ' In util2 - irtflg:',irtflg
           IF (IRTFLG .NE. 0) GOTO 9000

           IF (FCHAR(4:5) == 'SC') THEN
              CALL ARITHSCA(LUN1,LUN2,NSAM1,NROW1,NSLICE1,
     &                      FMIN,FMAX,FLOW,FHI)
           ELSE
              CALL ARITH(LUN1,LUN2,
     &                   NSAM1,NROW1,NSLICE1,EXPR(1:NLET))
           ENDIF
           !write(6,*) ' In util2 -:',imgnum,nslice1,maxim

        ENDDO

        GOTO 9000

C       ------------------- MIRROR SYMMETRY ---------------------  'MR'

C       OPEN INPUT FILE
15      CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &              MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

C       OPEN THE OUTPUT FILE
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM2,'OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        CALL MIRROR(LUN1,LUN2,NSAM1,NROW1,NSLICE1)
        GOTO 9000

C       DENSITY FOLDOVER ----------------------------------------- 'DF' 

16      CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

C       OPEN OUTPUT FILE
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM2,'OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        CALL DENOV(LUN1,LUN2,NSAM1,NROW1,NSLICE1)
        GOTO 9000

C       OPERATION  ----------------------------------------------- 'MA' 
C       MA       MASK
C       MA X     LATERAL MASKING IN SAMPLE DIRECTION
C       MA       MASK
C       MA SOFT  RADIAL SOFT MASKING 

17      CONTINUE
        IF (FCHAR(4:4) == 'S') THEN
           CALL SOFTMASK()
           GOTO 9000
        ENDIF


C       OPEN INPUT FILE
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000
        AV1   = AV
        FMIN1 = FMIN
        IF (IMAMI .NE. 1) 
     &      CALL NORM3(LUN1,NSAM1,NROW1,NSLICE1,FMAX1,FMIN1,AV1)

C       OPEN THE OUTPUT FILE
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM2,'OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        CALL MASK(LUN1,LUN2,NSAM1,NROW1,NSLICE1,AV1,FMIN1)
        GOTO 9000


C       OPERATION ------------------------------------------------ 'WV' 
C       WV        WINDOW AVERAGING
C       WV S      WINDOW AVERAGING -- SEQUENTIAL DOCUMENT SEARCH
C       WV P      WINDOW AVERAGING OVER PATCHES       

C       TRAP FOR WINDOW AVERAGING WITH PATCH OPTION
18      IF ( FCHAR(4:4) == 'P') THEN
          CALL WINAVE2(LUN1,LUN2,LUN3)
          GOTO 9000
        END IF

C       OPEN INPUT FILE
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

C       OPEN THE OUTPUT FILE
        NSAM2   = 0
        NSLICE2 = 0
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &               NSAM2,NROW2,NSLICE2,
     &               MAXIM2,'OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           IER = 4
           GOTO 9900
        ENDIF

        CALL WINAVE(LUN1,LUN2,LUN3,NSAM1,NROW1,NSAM2,NROW2)
        GOTO 9000

C       OPERATION  ----------------------------------------------- 'PP' 
C       PP P    PUT POINT 
C       PP L    PUT POINT SPECIFIED INTENSITY
C       PP LL   PUT LINE

19      CONTINUE
        IF (FCHAR(4:5) == 'LL') THEN
C          PUT CONTINUOUS LINE IN IMAGE FILE
           CALL PUTLIN(LUN1,LUN2)


        ELSEIF (FCHAR(4:4) == 'S') THEN
C          CONVERT FUNCTION TO A SURFACE IN A VOLUME

6400       ISIZE(1) = 32
           ISIZE(2) = 32
           ISIZE(3) = 32
           CALL RDPRI3S(ISIZE(1),ISIZE(2),ISIZE(3),NOT_USED,
     &                  'X, Y, & Z SIZES',IRTFLG)
           IF (IRTFLG == -1) GOTO 9000
          
6401       FWA(1)   = 1.0
           FWA(2)   = 1.0
           FWA(3)   = 1.0
           CALL RDPRM3S(FWA(1),FWA(2),FWA(3),NOT_USED,
     &           'NUMBER OF REPEATS IN X, Y, & Z',IRTFLG)
           IF (IRTFLG == -1) GOTO 6400

6402       IORDER(1) = 1
           IORDER(2) = 2
           IORDER(3) = 3
           CALL RDPRI3S(IORDER(1),IORDER(2),IORDER(3),
     &         NOT_USED,'ORDER FOR X, Y, & Z',IRTFLG)
           IF (IRTFLG == -1) GOTO 6401
           IF ((IORDER(1) < 1 .OR. IORDER(1) > 3) .OR.
     &         (IORDER(2) < 1 .OR. IORDER(2) > 3) .OR.
     &         (IORDER(3) < 1 .OR. IORDER(3) > 3)) THEN
              CALL ERRT(101,'IMPROPER ORDER (MUST BE 1...3)',NDUM)
              GOTO 6402
           ENDIF

6403       CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &         'SURFACE (G,D,D2)',NULL,IRTFLG)
           IF (IRTFLG == -1) GOTO 6402

           IFORM  = 3
           NSAM   = ISIZE(IORDER(1))
           NROW   = ISIZE(IORDER(2))
           NSLICE = ISIZE(IORDER(3))
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'U',ITYPE,
     &                  NSAM,NROW,NSLICE,
     &                  MAXIM,'OUTPUT',.FALSE.,IRTFLG)
           IF (IRTFLG == -1) GOTO 6401

           TINY = 1.0E-2
           CALL SURFTOVOL(LUN1,ISIZE,FWA,IORDER,ANS,TINY)
    
        ELSE
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &                  NSAM1,NROW1,NSLICE1,
     &                  MAXIM,'INPUT',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000
           IF (NSAM1 > MAXDIM) THEN
C             TOO LONG LINE LENGTH
              IER = 9
              GOTO 9900
           ENDIF

           IF (FCHAR(4:4) == 'P') THEN
              CALL PUTPT(LUN1,LUN2,NSAM1,NROW1,NSLICE1)

           ELSEIF (FCHAR(4:4) == 'L') THEN
              CALL PUTPT2(LUN1,LUN2,NSAM1,NROW1,NSLICE1)

           ELSEIF (FCHAR(4:4) == 'V') THEN
C             CONVERT HEIGHT FIELD TO A BINARY VOLUME

C             NORMALIZE FILE IF NECESSARY
              FMIN1 = FMIN
              FMAX1 = FMAX
              IF (IMAMI .NE. 1) 
     &            CALL NORM3(LUN1,NSAM1,NROW1,NSLICE1,
     &                       FMAX1,FMIN1,AV)

              WRITE(NOUT,90) FMIN,FMAX
90            FORMAT(' IMAGE DEPTH RANGE: ',G11.3,' ... ',G11.3)

6511            CALL RDPRI1S(NSLICE2,NOT_USED,
     &          'NUMBER OF SLICES IN OUTPUT VOLUME',IRTFLG)
              IF (IRTFLG == -1) GOTO 19

              IFORM = 3
              CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &                     NSAM1,NROW1,NSLICE2, 
     &                     MAXIM2,'OUTPUT',.FALSE.,IRTFLG)
              IF (IRTFLG == -1) GOTO 6511

              CALL IMTOVOL(LUN1,NSAM1,NROW1,NSLICE2,
     &                     LUN2,FMIN1,FMAX1, MAXDIM)
 
           ELSE
C             INPUT FROM TERMINAL
              CALL PUTPT1(LUN1,NSAM1,NROW1,NSLICE1)

           ENDIF
        ENDIF
        GOTO 9000

C       OPERATION ------------------------------------------------ 'SZ' 
C       SZ      SHEARS AN IMAGE BY OFSETTING EACH ROW   

C       OPEN FIRST INPUT FILE
20      CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        CALL SQUEEZ(LUN1,LUN2,NSAM1,NROW1,NSLICE1,IERR)
        GOTO 9000

C       OPERATION   SQUARE ROOT (WURZEL) ------------------------- 'WU' 
C       TAKES THE SQUARE ROOT OF AN IMAGE.

C       OPEN INPUT FILE
21      CALL OPFILEC(0,.TRUE.,FILNAM1,LUN1,'O',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000
        IF (FMIN < 0) THEN
           CALL ERRT(101,'*** SQ. ROOT OF NEGATIVE NUMBER AVOIDED',IE)
           GOTO 9000
        ENDIF

C       OPEN OUTPUT FILE
        CALL OPFILEC(LUN1,.TRUE.,FILNAM3,LUN3,'U',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        CALL IMSQ(.TRUE.,LUN1,LUN3,ITYPE,
     &             NSAM1,NROW1,NSLICE1,IRTFLG)
        GOTO 9000

C       OPERATION    MASK MULTIPLICATION ------------------------- 'MM' 
C       MM           MASK MULT.
C       MM C         MULT. CONTINUOUS
C       CAN HANDLE FOURIER x REAL

22      CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM1,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'MASK REFERENCE',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000
        FILNAM2 = FILNAM

C       OPEN IMAGE INPUT FILE
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'O',IFORM2,
     &               NSAM2,NROW2,NSLICE2,
     &               MAXIM2,'IMAGE (OVERWRITTEN!)',.TRUE.,IRTFLG)
        IF (IRTFLG == -1) GOTO 22
        IF (IRTFLG .NE. 0) GOTO 9000

C       MASKMU EXTENDED TO 3-D (JMC)
        CALL MASKMU(LUN1,LUN2,NSAM1,NROW1,NSLICE1)

C       SET UNDETERMINED SATATISTICS FLAG 
        CALL SETPRMB(LUN2, 0.0,0.0, 0.0,0.0)
        GOTO 9000

C       OPERATION ------------------------------------------------ 'CM'

C       OPEN OUTPUT FILE
23      CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'U',ITYPE,256,256,1,
     &             MAXIM,'OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000
    
C       NO MANUAL CHAPTER EITHER!!
C       CALL CLUMAP(LUN1,LUN2,256) !!!! commented out??????
        GOTO 9000

C       OPERATION  ----------------------------------------------- 'PV' 
24      CALL ERRT(100,'PV OPERATION NO LONGER SUPPORTED',NE)
CCC     CALL PDPVAX(LUN1,LUN2,FCHAR(4:4))
        GOTO 9000

C       OPERATION  ----------------------------------------------- 'NK' 
C       NK   SHRINK
                
C       POORLY IMPLEMENTED SHRINK OPERATION REMOVED NOV 2010 al
25      CALL ERRT(101,
     &   'OPERATION FAILS ON CONVEX SURFACES, USE <SE> INSTEAD',NE)  
C       CALL SHR(LUN1,LUN2,NSAM1,NSAM2)
        GOTO 9000

C       OPERATION ------------------------------------------------ 'AS'
C       AS F      AVG. AND COMPUTE STATISTICS AND Q FACTOR MAP 
C       AS DC     AVG. WITH VAR.  (NOW "AS R")        
C       AS        AVG. WITH VAR.  
C       AS R      AVG. WITH VAR.  BETTER       
C       AS A      AVG. WITH OPTIONAL VAR. & SUBSETS  (NEW 2012)

C       ADDS 2 OR MORE (UP TO 500) PICTURES TOGETHER, COMPUTING
C          SEVERAL MEASURES OF VARIANCE (PER POINT, TOTAL, ETC.)
                
26      IF (FCHAR(4:4) == 'F') THEN
C          AVG AND COMPUTE STATISTICS AND Q FACTOR MAP 
           CALL QFACT(LUN1,LUN2,LUN3)

        ELSE IF (FCHAR(4:5) == 'DC') THEN
           WRITE(NOUT,*) ' OBSOLETE OPERATION, RENAMED: <AS R>'
           FCHAR(4:5) = 'R '
           CALL ADDS(LUN1,LUN2,LUN3,LUN4,MAXDIM)

        ELSE IF (FCHAR(4:5) == 'AD') THEN
C          ADD IMAGE TO PRECOMPUTED AVERAGE WITH OFF IN INPUT REGISTER  
           WRITE(NOUT,*) ' OBSOLETE OPERATION'
C          AVG WITH OPTIONAL VARIANCE & SUBSETS  
           CALL ADDS(LUN1,LUN2,LUN3,LUN4,MAXDIM)

        ELSE IF (FCHAR(4:4) == 'S') THEN
C          AVG WITH OPTIONAL VARIANCE & SUBSETS  
           CALL ADDS_N()

        ELSE IF (FCHAR(4:4) == 'R') THEN
C          AVG WITH VARIANCE COMPUTED BETTER       
           CALL ADDS(LUN1,LUN2,LUN3,LUN4,MAXDIM)

        ELSE
C          AVG WITH VARIANCE COMPUTED POOR = OBSOLETE
           WRITE(NOUT,*)
     &         ' GIVES POOR VARIANCES, USE <AS AV> or <AS R> INSEEAD'
           CALL ADDS(LUN1,LUN2,LUN3,LUN4,MAXDIM)
        ENDIF
        GOTO 9000

C       OPERATION  ----------------------------------------------- 'MN' 
C       MN        MONTAGE   
C       MN S      MONTAGE WITH INDIVIDUAL SCALING  

27      CALL MONT(MAXDIM)
        GOTO 9000


C       OPERATION ------------------------------------------------ 'TH'
 
C       TH        THRESHOLD    
C       TH F      THRESHOLD--FIXUP CONSTANT 
C       TH M      THRESHOLD AND CREATE MASK 
C       TH C      THRESHOLD--CHANGE A VALUE 
   
C       REPLACES ALL VALUES WITHIN AN IMAGE OR VOLUME BEYOND A
C       SPECIFIED THRESHOLD BY THE THRESHOLD VALUE.

C       OPEN INPUT FILE
28      CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000
        IF (NSAM1 > MAXDIM) THEN
          IER = 9
          GOTO 9900
        END IF

        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM2,'OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        IF (FCHAR(4:4) == 'C') THEN
           CALL CHANGEVAL(LUN1,LUN2,NSAM1,NROW1,NSLICE1)
        ELSE
           CALL THRESH(LUN1,LUN2,NSAM1,NROW1,NSLICE1)
        ENDIF
        GOTO 9000

C       OPERATION  GP -- (GET PIXEL VALUE)----------------------- 'GP' 
29      CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        CALL GPRP(LUN1,NSAM1,NROW1,NSLICE1,FCHAR)
        GOTO 9000

C       OPERATION  RP -- (REPLACE PIXEL) ------------------------ 'RP' 
C       RP        REPLACE PIXEL         
30      CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        CALL GPRP(LUN1,NSAM1,NROW1,NSLICE1,FCHAR)
        GOTO 9000

C       OPERATION    MAXIMUM ------------------------------------ 'MX' 
C       COMPARES CORRESPONDING PIXELS OF 2 REAL IMAGES AND WRITES
C       THE MAXIMUM PIXEL VALUE AT THE CORRESPONDING POSITION OF THE
C       OUTPUT FILE

31      CALL OPFILEC(0,.TRUE.,FILNAM1,LUN1,'O',ITYPE,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'FIRST INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000
        IF (ITYPE < 0) THEN
           IER = 2
           GOTO 9900
        ENDIF     

C       OPEN IMAGE INPUT FILE
        CALL OPFILEC(0,.TRUE.,FILNAM2,LUN2,'O',ITYPE,
     &               NSAM2,NROW2,NSLICE2,
     &               MAXIM2,'SECOND INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000
        IF (ITYPE < 0) THEN
           IER = 2
           GOTO 9900
        ENDIF
        
        IF ((NSAM1   .NE. NSAM2) .OR.
     &      (NROW1   .NE. NROW2) .OR. 
     &      (NSLICE1 .NE. NSLICE2) ) THEN 
           IER = 1
           GOTO 9900     
        ENDIF

        IF ((MAXIM .NE. -2) .AND. (MAXIM2 .NE. -2) ) THEN
           IER = 2
           GOTO 9900             
        ENDIF
                        
        NSAM3   = NSAM1
        NROW3   = NROW1
        NSLICE3 = NSLICE1
        MAXIM   = 0
        CALL OPFILEC(0,.TRUE.,FILNAM3,LUN3,'U',ITYPE,
     &               NSAM3,NROW3,NSLICE3,
     &               MAXIM,'OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        CALL MX(LUN1,LUN2,LUN3,NSAM1,NROW1,NSLICE1)
        GOTO 9000



9900    CALL ERRT(IER,'UTIL2 ',NE)


9000    CLOSE(LUN1)
        CLOSE(LUN2)
        CLOSE(LUN3)

9911    RETURN
        END



