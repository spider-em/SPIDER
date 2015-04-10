C++*********************************************************************
C
C CORR1.F       REWRITTEN                          AUG 96 PP
C               ADDED 'CC H'                       MAR 02 ARDEAN LEITH
C               OPFILEC                            FEB 03 ARDEAN LEITH
C               'CC MS' bug                        OCT 03 ARDEAN LEITH
C               FMRS USED, UNUSED ALLOCS REMOVED   FEB 08 ARDEAN LEITH
C               PHASE REMOVED                      FEB 08 ARDEAN LEITH
C               MOD PGI COMPILER BUG               FEB 08 ARDEAN LEITH
C               INPLACE PARAMETERS                 APR 09 ARDEAN LEITH
C               RECOVERED 'CC H'                   JUL 11 ARDEAN LEITH
C               NRMS --> NORMVALSP                 DEC 11 ARDEAN LEITH
C               NPROMPTS, ==                       MAR 15 ARDEAN LEITH
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
C  CORR1
C
C  'CC N' CALL TREE:
C  CORR1 ---> READV
C        ---> NORMVALSP
C        ---> FMRS 
C                                        
C        ---> READV
C        ---> NORMVALSP
C        ---> FMRS 
C
C        ---> CCRS
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE CORR1

        USE TYPE_KINDS
        INTEGER(KIND=I_8)    :: IPLAN = 0     !STRUCTURE POINTER 

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

        CHARACTER(LEN=MAXNAM) :: FILNAM1,FILNAM2,FILNAMM

        COMPLEX, ALLOCATABLE  :: QK1(:,:,:)
        COMPLEX, ALLOCATABLE  :: QK2(:,:,:)
        COMPLEX, ALLOCATABLE  :: QKB(:)

        CHARACTER(LEN=1)      :: NULL = CHAR(0)

        PARAMETER (NFUNC=3)
        CHARACTER(LEN=2)      :: FUNC(NFUNC)

        LOGICAL               :: ACASE, BOTH_INC
        LOGICAL               :: SPIDER_SIGN
        LOGICAL               :: SPIDER_SCALE
	DOUBLE PRECISION      :: DAV,DSIG

        INTEGER, PARAMETER    :: LUN1 = 21
        INTEGER, PARAMETER    :: LUN2 = 22
        INTEGER, PARAMETER    :: LUN3 = 23

        DATA            FUNC/'AC',  'CC',  'CN'/

        IRTFLG = 0
 
C       DETERMINE IFUNC
        DO  IFUNC = 1,NFUNC
          IF (FCHAR(1:2) == FUNC(IFUNC)(1:2)) THEN
            GOTO 1111
          ENDIF
        ENDDO
C       OPERATION NOT HERE, RETURN TO CALLER
        RETURN 

1111    CONTINUE
C       CATCH EXCEPTIONS
C       ---------------------------------------------------------- CC C
        IF (FCHAR(1:4) == 'CC C')  THEN
           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FILNAM1,LUN1,'O',IFORM,
     &                NX1,NY1,NZ1,MAXIM,'FIRST INPUT',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FILNAM2,LUN2,'O',IFORM,
     &                NX2,NY2,NZ2,MAXIM,'SECOND INPUT',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              GOTO 998

           ELSEIF (NX1 .NE. NX2 .OR.
     &             NY1 .NE. NY2 .OR.
     &             NZ1 .NE. NZ2) THEN
              CALL ERRT(1,'CORR1',NE)
              GOTO 998
           ENDIF
           
           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FILNAMM,LUN3,'O',IFORM,
     &                NXM,NYM,NZM,MAXIM,'MASK',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0)  GOTO 998

           IF (NX1 .NE. NXM .OR. 
     &         NY1 .NE. NYM .OR.
     &         NZ1 .NE. NZM) THEN
               CALL ERRT(1,'CCC',NE)
               GOTO 998
           ENDIF
        
           CALL CCC(LUN1,FILNAM1,NX1,NY1,NZ1,
     &              LUN2,FILNAM2,NX2,NY2,NZ2,
     &              LUN3,FILNAMM)
          
           GOTO 998

C           ----------------------------------------------------- CC MS
        ELSEIF (FCHAR(1:5) == 'CC MS') THEN
C          CROSS CORRELATION - MASKED AND NORMALIZED 
           CALL MCCF
           RETURN

C          ------------------------------------------------------ AC MS
        ELSEIF (FCHAR(1:5) == 'AC MS')  THEN
          IF (FCHAR(6:6) == 'S')  THEN
C            AUTO CORRELATION - MASKED AND NORMALIZED 
              CALL MACF('S')
           ELSE
C             AUTO CORRELATION - MASKED AND NORMALIZED 
              CALL MACF(' ')
           ENDIF
           RETURN

C          ------------------------------------------------------ CC P
        ELSEIF (FCHAR(1:4) == 'CC P')  THEN
C          THIS APPEARS TO BE UNUSED ELSEWHERE IN SPIDER??
           CALL POLAR_CC 
           RETURN
        ENDIF


C       --------------------------------------------------------- CC H
C       --------------------------------------------------------- CC N 
C       --------------------------------------------------------- CC 
C       --------------------------------------------------------- AC 
C       OPEN FIRST INPUT FILE, FOURIER INPUT ALLOWED
        MAXIM = 0 
        CALL OPFILEC(0,.TRUE.,FILNAM1,LUN1,'O',IFORM,
     &               NX1,NY1,NZ1,MAXIM,'INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (FCHAR(4:4) == 'N' .AND. IFORM < 0) THEN
C          NO STATISTICS IN FILE, CAN NOT NORMALIZE OUTPUT
           CALL ERRT(101, 
     &        'CAN NOT NORMALIZE OUTPUT - FILE LACKS STATISTICS.',NE)
           GOTO 998

        ELSEIF (FCHAR(4:4) == 'N')  THEN
C          FOR 'CC N'
           IMAMI1 = IMAMI
           SIG1   = SIG
        ENDIF
        IFORM1 = IFORM

C       CALCULATE DIMENSIONS
        IF (IFORM1 > 0)  THEN
C          REAL SPACE INPUT
           LS1    = NX1 + 2 - MOD(NX1,2)
           LREC1  = NX1
           IFORM3 = IFORM1
           NX3  = NX1
        ELSE
C          FOURIER SPACE INPUT
           LS1   = NX1
           LREC1 = NX1
           NX1 = NX1 - MOD(-IFORM1,10)
           IF (IFORM1 > -20)  THEN
              IFORM3 = 1
           ELSE
              IFORM3 = 3
           ENDIF
           NX3 = NX1
        ENDIF

        IF (FCHAR(1:2) == 'AC') THEN
C          AUTO CORRELATION WANTED, NO SECOND FILE
           ACASE  = .TRUE.
           IFORM2 = IFORM1
        ELSE 
C          GET NAME FOR SECOND INPUT FILE
           CALL FILERD(FILNAM2,NLET,NULL,'REFERENCE',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 998

           IF (FILNAM1 == FILNAM2) THEN
C             FILENAMES ARE SAME, AUTOCORRELATION WANTED
              ACASE  = .TRUE.
              IFORM2 = IFORM1
           ELSE
C             OPEN SECOND INPUT FILE, FOURIER INPUT ALLOWED 
              MAXIM = 0
              CALL OPFILEC(0,.FALSE.,FILNAM2,LUN2,'O',IFORM,
     &                NX2,NY2,NZ2,MAXIM,NULL,.TRUE.,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 998

              IFORM2 = IFORM
              ACASE  = .FALSE.

              IF (FCHAR(4:4) == 'N' .AND. IFORM < 0) THEN
C                NO STATISTICS IN FILE, CAN NOT NORMALIZE OUTPUT
                 CALL ERRT(101, 
     &           'CAN NOT NORMALIZE OUTPUT - FILE LACKS STATISTICS.',I)
                 GOTO 998

              ELSEIF (FCHAR(4:4) == 'N')  THEN
C                FOR 'CC N'
                 IMAMI2 = IMAMI
                 SIG2   = SIG
              ENDIF
C             CALCULATE DIMENSIONS
              IF (IFORM2 > 0)  THEN
C                REAL SPACE INPUT
                 LS2   = NX2 + 2 - MOD(NX2,2)
              ELSE
C                FOURIER SPACE INPUT
                 LS2   = NX2
                 NX2 = NX2 - MOD(-IFORM2,10)
              ENDIF

C             CHECK THAT DIMENSIONS ARE THE SAME FOR BOTH FILES
              IF (NX1 .NE. NX2 .OR. 
     &            NY1 .NE. NY2 .OR.
     &            NZ1 .NE. NZ2) THEN
C                 ERROR. IMAGES DO NOT HAVE SAME DIMENSIONS
                  CALL ERRT(1,'CORR1',NE)
                  GOTO 998
              ENDIF
           ENDIF
        ENDIF

C       -------------------------------------------------- FIRST INPUT
        ALLOCATE(QK1(LS1/2,NY1,NZ1), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           MWANT = LS1 / 2 * NY1 * NZ1 
           CALL ERRT(46,'CORR1, QK1',MWANT)
           GOTO 996
        ENDIF

        CALL READV(LUN1,QK1,LS1,NY1,LREC1,NY1,NZ1)

        IF (FCHAR(4:4) == 'N' .AND. IMAMI1 .NE. 1) THEN
           !CALL NRMS(QK1,LS1,NX1,NY1,NZ1,SIG1)
           CALL NORMVALSP(QK1,NX1,NY1,NZ1,
     &                        LS1,  NY1,NZ1, 
     &                        DAV,DSIG, .TRUE.)
           SIG1 = DSIG
        ENDIF

        IF (IFORM1 > 0)  THEN
C          REAL SPACE INPUT, TRANSFORM IT TO FOURIER
           INV          = +1
           SPIDER_SIGN  = .TRUE.
           SPIDER_SCALE = .TRUE.
           CALL FMRS(QK1, NX1,NY1,NZ1,IPLAN, 
     &               SPIDER_SIGN,SPIDER_SCALE, INV,IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(101,'FOURIER TRANSFORM FAILED',NE)
              GOTO 996
           ENDIF
        ENDIF

        IF (FCHAR(4:4) == 'N')  THEN
           QK1(1,1,1) = (0.0,0.0)
        ENDIF

C       ------------------------------------------------- SECOND INPUT
        IF (.NOT. ACASE .AND. IFORM2 > 0)  THEN
C          CROSS-CORRELATION WITH REAL IMAGES, USE IN-CORE
           BOTH_INC = .TRUE.

           ALLOCATE(QK2(LS2/2,NY2,NZ2), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              MWANT = LS2 / 2 * NY2 * NZ2 
              CALL ERRT(46,'CORR1, QK2',MWANT)
              GOTO 998
           ENDIF

           CALL READV(LUN2,QK2,LS2,NY2,NX2,NY2,NZ2)

           IF (FCHAR(4:4) == 'N' .AND. IMAMI2 .NE. 1) THEN
              !CALL NRMS(QK2,LS2,NX2,NY2,NZ2,SIG2)
              CALL NORMVALSP(QK2,NX2,NY2,NZ2,
     &                           LS2,  NY2,NZ2, 
     &                           DAV,DSIG, .TRUE.)
              SIG2 = DSIG
           ENDIF
              
C          REAL SPACE INPUT, TRANSFORM IT TO FOURIER
           INV          = +1
           SPIDER_SIGN  = .TRUE.
           SPIDER_SCALE = .TRUE.
           CALL FMRS(QK2, NX2,NY2,NZ2, IPLAN, 
     &               SPIDER_SIGN,SPIDER_SCALE, INV,IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(101,FOURIER TRANSFORM FAILED,NE)
              GOTO 996
           ENDIF
        ELSE
           BOTH_INC = .FALSE.

           MWANT = NX1 + 2 - MOD(NX1,2) / 2 
           ALLOCATE (QKB(MWANT), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN 
              CALL ERRT(46,'CORR1, QKB',MWANT)
              GOTO 996
           ENDIF
        ENDIF

C       OPEN OUTPUT FILE
        IFORM = IFORM3
        MAXIM = 0
        CALL OPFILEC(LUN1,.TRUE.,FILNAM1,LUN3,'U',IFORM3,
     &               NX3,NY1,NZ1,MAXIM,'OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 996

C       -------------------------------------------------- CORRELATION

        LS = NX1 + 2 - MOD(NX1,2)

        IF (ACASE)  THEN
C          'AC' AUTO CORRELATION ------------------------------- 'AC'
           IF (NZ1 .LE. 1)  THEN
C            IMAGE INPUT
             IF (FCHAR(4:4) == 'S' .OR. FCHAR(4:5)== 'NS')  THEN
                CALL ACRS_2S(QK1, LS,NX1,NY1)
             ELSE            
                CALL ACRS_2(QK1,  LS,NX1,NY1)
             ENDIF
           ELSE
C            VOLUME INPUT
             IF (FCHAR(4:4) == 'S' .OR. FCHAR(4:5) == 'NS')  THEN
                CALL ACRS_3S(QK1, LS,NX1,NY1,NZ1)
             ELSE            
                CALL ACRS_3(QK1,  LS,NX1,NY1,NZ1)
             ENDIF
           ENDIF
           SIG2 = SIG1

        ELSEIF (FCHAR(4:4) == 'H')  THEN
C          'CC' PHASE CROSS CORRELATION --------------------- 'CC H'

           IF (BOTH_INC)  THEN
C             BOTH FILES ARE AVAILABLE IN BUFFERS

              IF (NZ1 > 1) THEN
                 CALL CCRS_PH_3(QK1,QK2, LS,NX1,NY1,NZ1,IRTFLG)
              ELSE
                 CALL CCRS_PH_2(QK1,QK2, LS,NX1,NY1,IRTFLG)
              ENDIF
              IF (IRTFLG .NE. 0) GOTO 996

           ELSE 
C             MUST READ SECOND FILE
              IF (NZ1 > 1) THEN
                 CALL CCRD_PH_3(LUN2,QK1,QKB, LS,NX1,NY1,NZ1)
              ELSE
                 CALL CCRD_PH_2(LUN2,QK1,QKB, LS,NX1,NY1)
              ENDIF
           ENDIF

        ELSEIF (FCHAR(1:2) == 'CC')  THEN
C          'CC' CROSS CORRELATION ------------------------------ 'CC'

           IF (BOTH_INC)  THEN
C             BOTH FILES ARE AVAILABLE IN BUFFERS
              SPIDER_SIGN  = .TRUE.
              SPIDER_SCALE = .TRUE.
              CALL CCRS(QK1,QK2, LS,NX1,NY1,NZ1,
     &                  SPIDER_SIGN,SPIDER_SCALE, IRTFLG)
              IF (IRTFLG .NE. 0) THEN
                 CALL ERRT(101,'CCRS FAILED',NE)
                 GOTO 996
              ENDIF

           ELSE 
C             MUST READ SECOND FILE
              IF (NZ1 > 1) THEN
                 CALL CCRD_3(LUN2,QK1,QKB, LS,NX1,NY1,NZ1)
              ELSE
                 CALL CCRD_2(LUN2,QK1,QKB, LS,NX1,NY1)
              ENDIF
           ENDIF

        ELSEIF (FCHAR(1:2) == 'CN')  THEN
C          'CN' CONVOLUTION, (NOT CORRELATION)
           IF (BOTH_INC)  THEN
C             BOTH FILES ARE AVAILABLE IN BUFFERS
              IF (NZ1 > 1) THEN
                 CALL CNRS_3(QK1,QK2, LS,NX1,NY1,NZ1)
              ELSE
                 CALL CNRS_2(QK1,QK2, LS,NX1,NY1)
              ENDIF
           ELSE
C             MUST READ SECOND FILE
              IF (NZ1 > 1) THEN
                 CALL CNRD_3(LUN2,QK1,QKB, LS,NX1,NY1,NZ1)
              ELSE
                 CALL CNRD_2(LUN2,QK1,QKB, LS,NX1,NY1)
              ENDIF
           ENDIF
        ENDIF

C       ------------------------------------------------------ OUTPUT
        IF (FCHAR(4:4) == 'N') THEN
C          NORMALIZE  HERE
           FAN = 1.0 / (NX1*NY1*FLOAT(NZ1)-1.0) / SIG1 / SIG2
           QK1 = QK1 * FAN
        ENDIF

C       THIS ONLY WRITES FIRST LS1 VALUES FROM EACH ROW
        CALL WRITEV(LUN3,QK1,LS1,NY1,NX1,NY1,NZ1)

996     IF (ALLOCATED(QKB)) DEALLOCATE (QKB)
        IF (ALLOCATED(QK1)) DEALLOCATE(QK1)
        IF (ALLOCATED(QK2)) DEALLOCATE(QK2)       

998     CLOSE(LUN3)
        CLOSE(LUN2)
        CLOSE(LUN1)

        RETURN
        END
