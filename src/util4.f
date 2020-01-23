C++*********************************************************************
C
C  UTIL4.F            ADDED IQ                      SEP 97 ArDean Leith
C                     ADDED IQ SYNC                 JUN 99 ArDean Leith
C                     ADDED NEG                     JUN 99 ArDean Leith
C                     NEG BUG                       FEB 01 ArDean Leith
C                     'AP MQ I'                     APR 01 ArDean Leith
C                     ADDED 'AP RQ'                 OCT 01 Haixiao Gao
C                     ADDED 'IQ W'                  MAR 02 ArDean Leith
C                     ADDED 'IQ PAR'                JUN 02 ArDean Leith
C                     ADDED 'IQ GONE'               AUG 02 ArDean Leith
C                     ADDED 'MS' VOLUMES            AUG 02 ArDean Leith
C                     ADDED 'MS I'                  JAN 03 ArDean Leith
C                     OPFILEC                       FEB 03 ArDean Leith
C                     REMOVED 'AP MR'               APR 03 ArDean Leith
C                     USED APMASTER                 AUG 03 ArDean Leith
C                     MPI                           FEB 04 Chao Yang
C                     ADDED 'IQ PID'                JAN 05 ArDean Leith
C                     ADDED 'IQ R'                  NOV 05 ArDean Leith
C                     'MS IF' IFORM BUG             FEB 07 ArDean Leith
C                     'AP C'                        JUN 08 ArDean Leith
C                     REMOVED VAR3* NO MAN          MAY 09 ArDean Leith
C                     'AP SCC' OUT OF APMASTER      AUG 08 ArDean Leith
C                     'IQ DI'                       JAN 10 ArDean Leith
C                     INQUIRESYNC(.FALSE.,REMOVE)   AUG 10 ArDean Leith
C                     'AP TOOL'                     OCT 10 ArDean Leith
C                     REPLACED VAR3Q,VAR 3R         MAY 11 ArDean Leith
C                     NX, 'AP FOU' 'H' mode         MAY 12 ArDean Leith
C                     REMOVED 'HF'                  DEC 12 ArDean Leith
C                     INQUIRECOMP                   APR 13 ArDean Leith
C                     'IQ PAR' PROMPT               JUN 15 ArDean Leith
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
C   UTIL4    DRIVER FOR CERTAIN ROUTINES
C
C--*********************************************************************

        SUBROUTINE UTIL4(MAXDIM)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER (LEN=MAXNAM) :: CID,CORRECT,FILNAM
        CHARACTER (LEN=1)      :: NULL,MODE
        LOGICAL                :: MAKEREFFILE,USEREFFILE,FLIP,FOLD

        LOGICAL                :: GETANS,UPPER,WANTSUB,SAYPRMT,SAYANS
        LOGICAL                :: STRIP,ENDATSEMI,REMOVE

        INTEGER                :: getpid  ! SYSTEM INTRINSIC FUNCTION

        INTEGER, PARAMETER     :: LUN  = 20
        INTEGER, PARAMETER     :: LUN1 = 21

        NULL   = CHAR(0)
        MAXIM1 = 0
        IRTFLG = 0

        IF (FCHAR(1:2) == 'AP')  THEN
C          OPERATION AP ------------------------------------------- AP

           IF (FCHAR(4:5) == 'RH' .OR.
     &         FCHAR(4:5) == 'MH' .OR.
     &         FCHAR(4:5) == 'NH')  THEN
              CALL ERRT(101,'OBSOLETE OPERATION, USE <AP REF>',NDUM)
              RETURN

           ELSEIF (FCHAR(4:5) == 'MQ' .OR. 
     &             FCHAR(4:5) == 'NQ' .OR. 
     &             FCHAR(4:5) == 'RQ') THEN
              CALL ERRT(101,'OBSOLETE OPERATION, USE <AP SH>',NDUM)
              RETURN

           ELSEIF (FCHAR(4:5) == 'ORN' .OR.
     &             FCHAR(4:5) == 'ORM') THEN
              CALL ERRT(101,'OBSOLETE OPERATION, USE <OR SH>',NDUM)
              RETURN

           ELSEIF (FCHAR(4:5) == 'CA')  THEN
              WRITE(NOUT,91)
91            FORMAT('  OBSOLETE OPERATION, USE: <AP C>',/)
              CALL HALI(.FALSE.)

           ELSEIF (FCHAR(4:5) == 'CM')  THEN
              WRITE(NOUT,91)
              CALL HALI(.FALSE.)

           ELSEIF (FCHAR(4:4) == 'C')   THEN
              CALL HALI(.TRUE.)

           ELSEIF (FCHAR(4:5) == 'RA')  THEN
              CALL FALB

           ELSEIF (FCHAR(4:5) == 'SA')  THEN
              CALL SAQB

           ELSEIF (FCHAR(4:5) == 'SR')  THEN
              CALL GALI

           ELSEIF (FCHAR(4:5) == 'MS')  THEN
              CALL MULTISHIFT

           ELSEIF (FCHAR(4:6) == 'SCC') THEN
C             2D & 3D PADDED, CROSS CORRELATION MULTI-REF SHIFT ALIGNMENT
              CALL APSCC()

           ELSEIF (FCHAR(4:5) == 'TO')  THEN
              MODE = 'F'
              CALL APMASTER_TOOL(MODE,FCHAR(4:5))

           ELSEIF (FCHAR(4:5) == 'FO')  THEN
              MODE = 'H'
              CALL APMASTER(MODE,FCHAR(4:))

           ELSE
              MODE = 'F'
              CALL APMASTER(MODE,FCHAR(4:))

          ENDIF


        ELSEIF(FCHAR(1:2) == 'HF')  THEN
C          OBSOLETE OPERATION REMOVED DEC 2012 --------------------- HF
           CALL ERRT(101,'OBSOLETE OPERATION REMOVED DEC 2012',NE)
           !CALL HF

        ELSEIF(FCHAR(1:2) == 'AT')  THEN
C          AUTOMATIC PARTICLE PICKING.------------------------------ AT
           IF (FCHAR(4:5) == 'IT')  THEN
              CALL DOCS1(MAXDIM)

           ELSEIF(FCHAR(4:5) == 'PK')  THEN
             CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &            NX,NY,NZ,MAXIM1,'INPUT',.FALSE.,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

              CALL ATPK(LUN1,NX,NY,NZ)
              CLOSE(LUN1)

           ELSEIF(FCHAR(4:5) == 'MC')  THEN
              CALL ATMC

           ELSEIF(FCHAR(4:5) == 'SA')  THEN
              CALL ATSA(MAXDIM)

           ELSEIF(FCHAR(4:5) == 'WN')  THEN
              CALL ATWN(MAXDIM)
           ENDIF


        ELSEIF(FCHAR(1:4) == 'MS I')  THEN
C          MAKE INLINE OR FILE BASED INDEXED STACK --------------- MS I
C          MAKE AN INLINE OR FILE BASED INDEXED FOURIER STACK----- MS IF

           MAXIM  = -1
           NX   = 0
           NZ = 0
           IFORM  = 0
           IF (FCHAR(5:5) == 'F') IFORM = -1
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN,'N',IFORM,NX,NY,NZ,
     &                   MAXIM,'NEW INDEXED STACK',.FALSE.,IRTFLG)
           CLOSE(LUN)


        ELSEIF(FCHAR(1:2) == 'MS')  THEN
C          MAKE AN INLINE STACK ---------------------------------- MS
C          MAKE AN INLINE FOURIER STACK -------------------------- MS F

C          SOLICIT FILE NAME
           IF (FCHAR(4:4) == 'F')  THEN
              CALL FILERD(FILNAM,NLET,NULL,
     &                   'NEW INLINE FOURIER STACK',IRTFLG)
              IFORM = -1
           ELSE
              CALL FILERD(FILNAM,NLET,NULL,'NEW INLINE STACK',IRTFLG)
              IFORM  = 0
           ENDIF
           IF (IRTFLG .NE. 0) RETURN

           IF (FILNAM(1:1) .NE. '_') THEN
              CALL ERRT(101,'NOT AN INLINE FILE',NE)
              RETURN
           ENDIF
           IF (NLET .LT. MAXNAM) FILNAM(NLET+1:) = CHAR(0)

           MAXIM  = 1
           NX   = 0
           NZ = 0
           CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'N',IFORM,NX,NY,NZ,
     &                   MAXIM,' ',.FALSE.,IRTFLG)


        ELSEIF(FCHAR(1:2) == 'NE')  THEN
C          NEGATE/INVERT AN IMAGE --------------------------------- NE

C          OPEN INPUT FILE
           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN,'O',IFORM,NX,NY,NZ,
     &             MAXIM,'INPUT',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          NEED FMAX & AV BELOW
           FMAXVAL = FMAX
           AVVAL   = AV
           IF (IMAMI==0) CALL NORM3(LUN,NX,NY,NZ,
     &                                FMAXVAL,FMINVAL,AVVAL)

C          OPEN OUTPUT FILE
           MAXIM = 0
           CALL OPFILEC(LUN,.TRUE.,FILNAM,LUN1,'U',IFORM,
     &                NX,NY,NZ,MAXIM,'OUTPUT',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CLOSE(LUN)
              RETURN
           ENDIF

C          FOR NEG A  (FCHAR ONLY HAS FIRST TWO LETTERS!!)
           IF (FCHAR(4:4) .NE. 'A')    THEN
C             NEGATE THEN ADD ORIGINAL FMAX TO EACH VALUE
              CALL NEGATE(LUN,LUN1,NX,NY,NZ,FMAXVAL)
           ELSE
C             NEGATE AROUND AVERAGE VALUE
              CALL NEGATI(LUN,LUN1,NX,NY,NZ,AVVAL)
           ENDIF
           CLOSE(LUN)
           CLOSE(LUN1)
           RETURN


        ELSEIF(FCHAR(1:2) == 'IQ')  THEN
C          INQUIRE SOMETHING -------------------------------------- IQ

           IF (FCHAR(4:5) == 'FI')  THEN
C             SEE IF FILE EXISTS
              CALL INQUIREIF('FI')

           ELSE IF (FCHAR(4:4) == 'D')  THEN
C             INQUIRE IF DIR EXISTS
              CALL INQUIREIF('DIR')

           ELSE IF (FCHAR(4:5) == 'SY')  THEN
C             WAIT TILL FILE EXISTS
              REMOVE =  (FCHAR(9:9) == 'D')
              CALL INQUIRESYNC(.FALSE.,REMOVE)

           ELSE IF (FCHAR(4:4) == 'G')  THEN
C             WAIT TILL FILE GONE
              CALL INQUIRESYNC(.TRUE.,.FALSE.)

           ELSE IF (FCHAR(4:4) == 'R')  THEN
C             CHECK ON REGISTER VARIABLE CONTENTS
              CALL INQUIREREG(.TRUE.,.TRUE.,IRTFLG)

           ELSE IF (FCHAR(4:4) == 'A')  THEN
C             CHECK ON ALLOCABLE MEMORY -------------------------- IQ A
              CALL RDPRM1S(GSTART,NOT_USED,'MEMORY DESIRED',IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

              CALL INQUIREALLOC(GSTART,IMBYTES,.TRUE.,IRTFLG)
              RGOT = IMBYTES
              CALL REG_SET_NSEL(1,1,RGOT,0.0, 0.0, 0.0, 0.0, IRTFLG)

           ELSE IF (FCHAR(4:4) == 'W')  THEN
C             CHECK ON MACHINE ARCHITECTURE
              CALL INQUIREARCH(LUN,FLIP,FOLD,IRTFLG)


           ELSE IF (FCHAR(4:4) == 'C')  THEN
C             CHECK ON COMPILATION
              CALL INQUIRECOMP(IVAL)

           ELSEIF (FCHAR(4:6) == 'PID') THEN
C             INQUIRE PROCESS ID  ------------------------------ IQ PID

#if defined (SP_GFORTRAN)
              IPID = getpid()
#else
              IPID = getpid(IPID)
#endif
              WRITE(NOUT,92) IPID
92            FORMAT('  Current process id: ',I9,/)

              CALL REG_GET_USED(NSEL_USED)
              IF (NSEL_USED .GT. 0) THEN
C                OUTPUT TO SPIDER'S REGISTERS/REAL VARIABLES
                 FPID = IPID
                 CALL REG_SET_NSEL(1,1,FPID,FPID,FPID,FPID,FPID,IRTFLG)
              ENDIF

           ELSEIF (FCHAR(4:4) == 'P') THEN
C             TEST OF PARAMETER SUBSTITUTION MECHANISM --------- IQ PAR

C             DO NOT UPPERCASE THE INPUT LINE, DO NOT STRIP AFTER ;
              GETANS    = .TRUE.
              UPPER     = .FALSE.
              WANTSUB   = .TRUE.
              SAYPRMT   = .TRUE.
              SAYANS    = .TRUE.
              ENDATSEMI = .TRUE.
              STRIP     = .TRUE.

              CALL RDPR('VARIABLE NAME (WITH [])',NCHAR,CID,GETANS,
     &             UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI,STRIP,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

              IRTFLG = -999
              CALL RDPRMC(CORRECT,NLET2,.TRUE.,'CORRECT VALUE',
     &                 NULL,IRTFLG)

              IF (IRTFLG == 0 .AND. 
     &            CID(1:NCHAR) .NE. CORRECT(1:NLET2)) THEN
                 WRITE(NOUT,90) CID(1:NCHAR), CORRECT(1:NLET2)
90               FORMAT(' *** GOT: ',A,'  SHOULD BE: ',A)
                 CALL ERRT(101,'SYMBOL SUBSTITUTION INCORRECT',NE)
               ENDIF

           ELSE
C             UNKNOWN OPTION
              CALL ERRT(101,'UNKNOWN OPTION',NE)
           ENDIF


        ELSEIF(FCHAR(1:2) == 'VA')  THEN
C          VARIANCE CALCULATION ----------------------------------- VA

           IF (FCHAR(4:4) == 'F')  THEN
C             FOURIER SPACE VARIANCE CALCULATION ----------------- VA F
              CALL VARF

           ELSEIF (FCHAR(4:5) == '3R')  THEN
C             VARIANCE CALCULATION ------------------------------- VA 3R
              CALL VAR3R

           ELSEIF (FCHAR(4:5) == '3Q')  THEN
C             VARIANCE CALCULATION ------------------------------- VA 3Q
              CALL VAR3D('Q')

           ELSE
              CALL ERRT(101,'UNDOCUMENTED, BUGGY OPERATION REMOVED',NE)
	   ENDIF


        ELSEIF (FCHAR(1:2) == 'SN')  THEN
C          SNR FROM FSC -------------------------------------------- SN

           IF (FCHAR(4:5) == 'RB')  THEN
C             APPROXIMATE SNR BY BUTTERWORTH FILTER ------------- SN RB
              CALL SNRB

           ELSEIF (FCHAR(4:5) == 'RF')  THEN
C             CREATE BUTTERWORTH FILTER  ------------------------ SN RF
              CALL SNRF
	   ENDIF
        ENDIF

        END









