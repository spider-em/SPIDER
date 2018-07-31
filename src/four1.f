C++*********************************************************************
C
C  FOUR1.F                                        08/22/96
C             OPFILEC                             FEB  03 ARDEAN LEITH
C             MPI                                 OCT  03 CHAO YANG
C             REMOVED UNDOCUMENTED OLD 'FD R'     MAY  14 ARDEAN LEITH
C             ADDED 'FQ Q'                        NOV  14 ARDEAN LEITH
C             ADDED 'FSC MA'                      APR  16 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2016  Health Research Inc.,                         *
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
C FOUR1(MAXMEM)
C
C CALLS: FQ, FT, FF, FL, FP, EF, PW, RF, CF, GF, RD, FD, 16=FSC, 17=FRC
C
C--*********************************************************************

        SUBROUTINE FOUR1(MAXMEM)

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 

        INTEGER, PARAMETER    :: NFUNC=14
        CHARACTER(LEN=2)      :: FUNC(NFUNC)

        CHARACTER(LEN=MAXNAM) :: FILNAM,FILNAM2
        REAL                  :: VALUES(6)
        LOGICAL               :: FSCOP

        INTEGER, PARAMETER    :: LUN1 =21
        INTEGER, PARAMETER    :: LUN2 =22
        INTEGER, PARAMETER    :: LUNF =23
        INTEGER, PARAMETER    :: LUN3 =24

        DATA FUNC/'FQ','FT', 'FF', 'FL', 'FP', 
     &            'EF','PW', 'RF', 'CF', 'GF', 
     &            'RD','FD','16','17'/

        MAXIM  = 0
        MAXIM2 = 0
        IRTFLG = 0

C                 FQ, FT, FF, FL, FP, EF, PW, RF, CF, GF, RD, FD, 
C                 16 = FSC, 17=FRC

        DO IFUNC = 1,NFUNC
            IF (FCHAR(1:2) == FUNC(IFUNC)) THEN
              GOTO ( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,8,8), (IFUNC)
            ENDIF
        ENDDO

C       OPERATION NOT IN FOUR1, RETURN TO SPIDER
        RETURN

C       ---------------- QUICK FILTERING ------------------------- 'FQ'

C       NEW INCORE FQ OPERATION
1       CALL FOUR_FQ
        RETURN

C       ---------------- FOURIER TRANSFORM ----------------------- 'FT'
2       IF (FCHAR(4:4) == 'R') THEN
           CALL FTR
        ELSE
           CALL FOUR1C
        ENDIF
        RETURN

C       ---------------- FOURIER FILTER -------------------------- 'FF'

C       APPLIES FILTERS TO 2-D OR 3-D FOURIER TRANSFORMS.
3       IF (FCHAR(4:7) == 'PLOT') THEN
           CALL FILTPLOT(MAXMEM)
           RETURN
        ENDIF

        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,
     &               NX,NY,NZ,
     &               MAXIM,'INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0)  RETURN

        IF (IFORM .NE. -11 .AND. IFORM.NE. -12 .AND.
     &      IFORM .NE. -21 .AND. IFORM.NE. -22) THEN
           CALL  ERRT(101,'OPERATION INCONSISTENT WITH DATA FORMAT',NE)
           GOTO 9001
        ENDIF

        IF (FCHAR(4:4) == 'S') THEN
           CALL  ERRT(41,'FF S',NE)
C          CALL FSHADO(LUN1,NX,NY)

        ELSEIF (FCHAR(4:4) == 'L' .OR. FCHAR(4:4) == 'B') THEN
           CALL  ERRT(41,'FF L/B',NE)
C          CALL FILTB(LUN1,NX,NY)

        ELSE
           NXO = NX-MOD(-IFORM,10)
           CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',IFORM,
     &             NX,NY,NZ,
     &             MAXIM,'OUTPUT',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0)  THEN
              CLOSE(LUN1)
              CALL ERRT(4,'FF',NE)
              RETURN
           ENDIF
           CALL FFILTS(LUN1,LUN2,NX,NY,NZ,NXO)
        ENDIF
        GOTO 9000

C       ---------------- FOURIER LISTING ------------------------- 'FL'
C       LISTS MODULI AND PHASES OF 2-D FOURIER TRANSFORMATION.
4       CALL  ERRT(101,'OBSOLETE OPERATION',NE)
        RETURN


C       ---------------- FOURIER INTERPOLATION ------------------- 'FP'
5       CALL FOUR1A_FP
        RETURN


c       --------- EXTRACT FOURIER -------------------------------- 'EF'
C       EXTRACTS CENTRAL SECTION FROM 3-D FOURIER UNDER ARBITRARY ANGLES.
6       CALL ERRT(101,'OBSOLETE OPERATION',NE)
        RETURN


C ---------------------- POWER SPECTRUM -------------------------- 'PW'
7       CALL FOUR1B
        RETURN

C       -------------- R-FACTOR ---------------------------------- 'RF'
C       -------------- FSC      ---------------------------------- 'FSC'
C       -------------- FRC      ---------------------------------- 'FRC'
C       COMPUTES MEASURES OF PROXIMITY BETWEEN 2 GIVEN TRANSFORMS
8       IF (FCHAR(4:6)     == '3SN')  THEN
           CALL SSNR3

        ELSEIF (FCHAR(4:6) == '3NN')  THEN
           CALL SSNR3DNN

        ELSEIF (FCHAR(1:2) == '16' .AND. FCHAR(4:4) == 'N') THEN
           CALL PR3D_NEW()    ! OPERATION: 'FSC NEW' 

        ELSEIF (FCHAR(1:2) == '16' .AND. FCHAR(4:4) == 'M') THEN
           CALL PR3D(.TRUE.)    ! OPERATION: 'FSC MA' (undocumented)

        ELSEIF (FCHAR(4:4) == '3' .OR. FCHAR(1:2) == '16') THEN
           FSCOP = (FCHAR(1:2) == '16')
           CALL PR3D(FSCOP)    ! OPERATION: 'FSC' OR 'RF 3'

        ELSEIF (FCHAR(4:5) == 'SN') THEN
           CALL SSNRB

        ELSE
           FSCOP = (FCHAR(1:2) == '17')  ! OPERATION: 'FRC'
           CALL RFACTSDO(FSCOP)

        ENDIF
        RETURN

C       ---------------- CONSTRUCT FOURIER ----------------------- 'CF'
C       CONSTRUCT FOURIER FILE FROM AMPLITUDES & PHASES OF REFLECTIONS.
9       CALL  ERRT(101,'OBSOLETE OPERATION',NE)
        RETURN

C       ---------------- GENERAL FILTER -------------------------- 'GF'
C       FOR QUASI-OPTICAL FOURIER FILTRATION 
10      CALL  ERRT(101,'OBSOLETE OPERATION',NE)
        RETURN

C       ---------------- REDUCE TRANSFORM ------------------------ 'RD'
C       GENERATES REDUCED FOURIER TRANSFORM FROM MASKED FOURIER
11      CALL  ERRT(101,'OBSOLETE OPERATION',NE)
        RETURN


C       FILTER ACCORDING TO A DOCUMENT FILE --------------------- 'FD'

12      IF (FCHAR(4:4) == 'R') THEN
           CALL  ERRT(101,'REMOVED UNDOCUMENTED OPERATION IN 2014',NE)
C          CALL  RADWEIGHT

        ELSE
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NX,NY,NZ,
     &               MAXIM,'INPUT',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           MAXIM = 0
           CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',IFORM,
     &               NX,NY,NZ, MAXIM,'OUTPUT',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           CALL  FILTDOC(LUN1,LUN2,NX,NY,NZ,IFORM)
        ENDIF

9000    CLOSE(LUN2)
9001    CLOSE(LUN1)

        END
