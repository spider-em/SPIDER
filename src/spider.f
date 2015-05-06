
C ++********************************************************************
C                                                                      *
C SPIDER : (System for Processing Image Data in                        *
C           Electron microscopy and Related fields)                    *
C                                                                      *
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
C                                                                      *
C  MAIN SUBROUTINE OF SPIDER IMAGE PROCESSING SYSTEM.                  *
C  UPDATE VERSION STATEMENT (MARKED BY CHERE) WHEN VERSION IS CHANGED! *
C                                                                      *
C  LUN ASSIGNMENTS: LUN     INTERACTIVE   IN PROC.         CONNECTS    *
C                   NLOGP        1           1                LOG      *
C                   NLOG         1           0                         *
C                   LUNSPIRE     2           2               SPIRE     *
C                   NDAT         3           3              RESULTS    *
C                   NDOC         4           4             A DOCFILE   *
C                   NSTDINP      5           5               STDIN     *
C                   NIN          5          5/1                        *
C                   NOUT         6           3                         *
C                   NSTDOUTP     6           6              STDOUT     *
C                                                                      *
C                   LUNTEXT    103                         SYMPAR TXT  *
C                   LUNDO      301                         LUNDO FILE  *
C                              200...200+MAXICDOCS          SAVDOCQ    *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

#ifdef  HAS_IEEE
#ifndef USE_MPI
#ifndef __APPLE__
        USE, intrinsic :: ieee_exceptions
#endif
#endif
#endif

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

C       @@@@@@@@@@@@@@@@@@@@@ PARAMETER  INITIALIZATION @@@@@@@@@@@@@@@@

C       ONE PLUS THE MAXIMUM NUMBER OF REGISTERS PER KEY IN 'UD IC'
        INTEGER, PARAMETER :: MAXREG=7     

C       MAXIMUM NUMBER OF KEYS IN 'UD IC'
        INTEGER, PARAMETER :: MAXKEY=9999  

C       MAXIMUM NUMBER OF REGISTER ARGUMENTS SENT TO A PROCEDURE
        INTEGER, PARAMETER :: NPARG=24    

C       MAXIMUM NUMBER OF NESTED 'IF's
        INTEGER, PARAMETER :: MAXIF=20      
    
C       COMMON BLOCK SPACE RESERVATION USED AT ALBANY, IS NOW ONLY 5 MB. 
C       (MOST ROUTINES NOW USE RUN-TIME ALLOCATION OF MEMORY WHICH 
C       IS INDEPENDENT OF THE COMMON BLOCK AND MAXDI ASSIGNMENT.
C       WE USUALLY HAVE > 2 GB RAM AVAILABLE ON ALBANY MACHINES
 
        INTEGER, PARAMETER  :: MAXDI = 5000000
        INTEGER             :: PLINEGO(MAXDI/5)
        CHARACTER           :: PDATA(4*4*MAXDI/5)
        COMMON   PLINEGO,PDATA

C       @@@@@@@@@@@@@@@@@@@  DECLARATIONS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        INTEGER   :: IDOSTK(7,MAXPRC),NARGSREC(MAXPRC),LOOPSV(8,MAXPRC)
        INTEGER   :: IARGSENT(NPARG,MAXPRC),IARGSREC(NPARG,MAXPRC)
        INTEGER   :: IFSV(MAXPRC)
        INTEGER   :: NUML(3)
        INTEGER   :: getpid
        INTEGER   :: NEWLOOP

        LOGICAL   :: USEELSE(MAXIF,MAXPRC)
        LOGICAL   :: JUMP,EX,ISDIGI,ISCHAR,DELETIT,LISTIT
        LOGICAL   :: RES_TO_TERM,GLOBAL,ISATAT

        CHARACTER(LEN=MAXNAM) :: PNAME
        CHARACTER(LEN=MAXNAM) :: PROCFL(MAXPRCNAM)
        CHARACTER(LEN=MAXNAM) :: RESULT,LOG,SPIRE_FILE

        CHARACTER(LEN=160)    :: MESG,PLINE,ARG4,ARGNOW,FCHARNOBLANK
        CHARACTER(LEN=40)     :: CVERS
        CHARACTER(LEN=12)     :: CDATT
        CHARACTER(LEN=8)      :: ZEIT
        CHARACTER(LEN=7)      :: RESULM
        CHARACTER(LEN=5)      :: LABEL 
        CHARACTER(LEN=4)      :: CXNUM,CREG
        CHARACTER(LEN=3)      :: LOGM
        CHARACTER(LEN=2)      :: NQ12
        CHARACTER(LEN=1)      :: NULL,RESPONSE

C       DBUF = COMMON TEMPORARY BUFFER FOR DOCUMENT FILE(S) 
C       NO LONGER USED FOR 'INCORE' DOCUMENT FILES 
        COMMON /DOC_BUF/ DBUF(MAXREG,MAXKEY,2)

        COMMON /DRIV1_COM/ T1,LOOPREG,CXNUM

C       MAKE SURE FIRST SIZING OF IOBUF IS LARGE
        COMMON /IOBUF/   BUFIO(NBUFSIZ)

C       MAKE SURE FIRST SIZING OF COMMUN IS LARGE
        INTEGER, PARAMETER :: NCOMSIZ = 8000
        COMMON /COMMUN/ BUFC(NCOMSIZ)

        COMMON /LUNDOECHO/ LUNDONOW,NDOLINE,NINSAVEOF

        INTEGER                    :: ITIME(8)

C       RANDOM NUMBER GENERATOR SEED
#if defined (SP_GFORTRAN)
        INTEGER,ALLOCATABLE        :: ISEED(:)       
#else
        INTEGER                    :: ISEED(2)       
#endif

C       FOR LOCAL VARIABLE HANDLING 
        INTEGER, DIMENSION(MAXPRC) :: IPSTACK,IPNUMSTACK,IPARNUM
        COMMON /QSTR_STUFF1/ ISTOPR,NSTDOUTP,NSTDINP,IWHERE,IPSTACK,
     &                       IPNUMSTACK,IPARNUM

        COMMON /PROC_STUFF/ NUMPROCNOW

C       SIZING OF GOTLAB FOR CHECKING DUPLICATE LABELS
        INTEGER, PARAMETER       :: MAXNUMLAB = 50
        INTEGER                  :: LABGOT(MAXNUMLAB)

C       LOGICAL UNIT NUMBERS DEFINED HERE
        INTEGER, PARAMETER       :: LUNDO    = 300   
        INTEGER, PARAMETER       :: LUNSPIRE = 2   
        INTEGER, PARAMETER       :: LUNTEXT  = 103   

        INTEGER  :: omp_get_stack_size
        INTEGER  :: isiz1,isiz2

C       @@@@@@@@@@@@@@@@@@@@@@@@@@ DATA STATEMENTS @@@@@@@@@@@@@@@@@@@@
C       @@@@@@@@@@@@@@@@@@@@@@ VERSION INITIALIZATION @@@@@@@@@@@@@@@@@

CHERE               123456789 123456789 123456789 1234567890 
        DATA CVERS/'VERSION:  UNIX  22.11 ISSUED:  5/05/2015'/

        DATA RESULM/'results'/
        DATA LOGM/'LOG'/
        DATA RES_TO_TERM/.FALSE./

C       SOME DO LOOP PARAMETERS
        DATA IDOTOP,IFLEVEL/1,0/

#ifdef USE_MPI
        INCLUDE 'mpif.h'
        DOUBLE PRECISION :: TT0, TT1   
        LOGICAL          :: ONLYONE_RED,ONLYONE_WRT
        COMMON /COMM_MPI/ONLYONE_RED,ONLYONE_WRT

C       @@@@@@@@@@@@@@@@@@@@@@@@@@@  CODE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        ONLYONE_RED = .TRUE.
        ONLYONE_WRT = .TRUE.

        CALL MPI_INIT(MPIERR)
        ICOMM = MPI_COMM_WORLD
        CALL MPI_COMM_RANK(ICOMM, MYPID,  MPIERR)
        CALL MPI_COMM_SIZE(ICOMM, NPROCS, MPIERR)
#ifdef MPI_DEBUG
        IF (MYPID == 0)WRITE(6,*) ' NPROCS, ICOMM: ', NPROCS,ICOMM
        TT0 = MPI_WTIME()
#endif
#else
C       NOT USING MPI
        MYPID = -1
#endif

#ifdef SP_MP
C       NEEDED BY PGI 2013.10 COMPILER

#ifndef USE_MPI
#ifndef __APPLE__
        isiz1 = omp_get_stack_size()
        CALL omp_set_stack_size(65536)
        isiz2 = omp_get_stack_size()
        !write(6,*) ' OMP Stack size: ',isiz1,' -->',isiz2

#endif
#endif
#endif

C       INITIALIZE COMMON BLOCK LUNS
C       NLOGP IS FOR LOG FILE,  NLOG = NLOGP WHEN LOG IS IN USE
        NLOGP    = 1  
        NLOG     = NLOGP
        NECHO    = 0        ! COUNTER FOR LOG FILE OUTPUT

C       NDAT IS FOR RESULTS FILE,  NOUT=NDAT=3 WHEN RESULTS FILE IN USE
        NDAT      = 3

        NSTDINP   = 5
        NIN       = NSTDINP
        NSTDOUTP  = 6
        NOUT      = NSTDOUTP

C       INITIALIZE SOME COMMON BLOCK DATA ELEMENTS (SEE: SETMODE.F)
        ISTOP          = 1
        ISTOPR         = 1
        IBCNT          = 0
        NLOOP          = 0       ! TOTAL NUMBER OF LOOP ITERATIONS
        ILOOP          = 1       ! CURRENT LOOP ITERATION
        IFOUND         = 1
        NTRACE         = 0
        VERBOSE        = .TRUE.
        SILENT         = .FALSE.
        LEGACYPAR      = .FALSE. ! () IN DO NO LONGER! DEC 2010
        USE_FBS_INTERP = .FALSE. ! NEW JUL 2011
        NECHO          = 0       ! COUNTER FOR LOG FILE OUTPUT
        MAXDIM         = MAXDI   ! SET SIZE OF COMMON BUFFER
        NUMFFTWTH      = 0       ! NUMBER OF FFTW3 THREADS
        NULL           = CHAR(0)
        NQ12           = CHAR(34) // CHAR(39)   ! QUOTES

C       SET ISEED  TO INITIAL "RANDOM" VALUE USING CLOCK
#if defined (SP_GFORTRAN)
        CALL RANDOM_SEED(SIZE = N)
        ALLOCATE(ISEED(N))
        CALL SYSTEM_CLOCK(COUNT = ICLOCK)
        ISEED = ICLOCK + 37 * (/ (I - 1, I = 1, N) /)
        CALL RANDOM_SEED(PUT = ISEED)
        DEALLOCATE(ISEED)
#else
        CALL DATE_AND_TIME(VALUES=ITIME)       ! GET CURRENT TIME
        ISEED(1) = ITIME(4) * (360000*ITIME(5) + 6000*ITIME(6) + 
     &             100*ITIME(7) + ITIME(8))
        IF (ISEED(1) == 0 .AND. ISEED(2) == 0) THEN
C          KLUDGE TO PREVENT ERROR ON SOME INTEL PROCESSORS
           write(0,*) ' Using default random number seed'
           CALL RANDOM_SEED()
        ELSE
           CALL RANDOM_SEED(PUT = ISEED)
        ENDIF
#endif

C       TIME TM IS ASSUMED AT BEGINNING OF RUN.  
        T1 = SECNDS(0.0)

C       INITIALIZE DO_LOOP STACK VARIABLES
        DO I=1,5
           IDOSTK(I,1) = 1
        ENDDO
        IDOSTK(4,1) = 0
        IDOSTK(6,1) = 0
        IDOSTK(7,1) = 1
        LOOPREG     = 0    ! REG. FOR LOOP COUNTER
        IFSV(1)     = 0

C       INITIAL MODE IS INTERACTIVE.
        COPT = 'I'

C       INITIALIZE PROCFL & NUMPRC USED TO TRACK OF PROCEDURE LISTING:
C       PROCEDURES LOADED & LISTED ON 1'ST OCCURRENCE ONLY
        DO NUMPROCNOWT=1,MAXPRCNAM
           PROCFL(NUMPROCNOWT)(:MAXNAM) = ' '
        ENDDO
        NUMPRC     = 1
        NUMPROCNOW = 1

C       PUT THE STARTING PROCEDURE FILE NAME ON THE STACK
        IPNUMSTACK(ISTOP) = 1
        PROCFL(1)(:MAXNAM) = 'INTERACTIVE'

C       INITIALIZE REGISTERS FOR MAIN BANK
        CALL REG_INIT(1,IRTFLG)
 
C       INITIALIZE GLOBAL SYMBOL STORAGE
        CALL SETSYMPAR(NULL,RESPONSE,.FALSE.,IRTFLG)
C       INITIALIZE SYMBOL INFO FOR INITIAL PROCEDURE = 1
        CALL SETSYMPAR(NULL,RESPONSE,.TRUE.,IRTFLG)
      
C       GET THE DATE AND TIME
        CALL DATE_2K(CDATT)
        CALL MYTIME(CTIM)

C       GET DATEXC, RESULTS FILE VERSION, FIRST OP & REGISTER SETTING
C       AND USE_SPIRE SETTING
        CALL INITUNIX(NUMARG,FCHAR,NALPH,CXNUM,MESG)

        IF (MYPID <= 0) THEN
C       PRINT OUT HEADING WITH VERSION AND RELEASE DATES
        WRITE(NOUT,*)' '
        WRITE(NOUT,9090)
        WRITE(NOUT,9091)
        WRITE(NOUT,9092)
        WRITE(NOUT,9093)CVERS
        WRITE(NOUT,9094)CDATT(1:11),CTIM

 9090   FORMAT('  \\__`O O''__/        SPIDER -- COPYRIGHT')
 9091   FORMAT('  ,__xXXXx___        HEALTH RESEARCH INC., ALBANY, NY.') 
 9092   FORMAT('   __xXXXx__')
 9093   FORMAT('  /  /xxx\\  \\        ',A)
 9094   FORMAT('    /     \\          DATE:     ',A,'    AT  ',A,/)

        WRITE(NOUT,9097)
 9097   FORMAT('  If SPIDER is useful, please cite:',/,
     &        '  Frank J, Radermacher M, Penczek P, Zhu J, Li Y,',
     &         ' Ladjadj M,  Leith A.',/,
     &       '  SPIDER and WEB: Processing and visualization of images',
     &         ' in 3D electron ',/,
     &        '  microscopy and related fields.  J. Struct. Biol.',
     &         ' 1996; 116: 190-199.')
        WRITE(NOUT,*) ' '

        ENDIF

        IF (NUMARG <= 0) THEN
C          GET THE PROJECT AND DATA EXTENSION FROM USER
           NLOG       = 0
           FCHAR(1:3) = 'NC' // CHAR(0)
           CALL DRIV1(MAXDIM)
           NLOG       = NLOGP
        ENDIF

C       CREATE NAME FOR LOG FILE AND OPEN THE LOG FILE
        CALL FILNAMANDEXT(LOGM,PRJEXC,LOG,NLET,.FALSE.,IER)
        IF (IER .NE. 0)
     &     STOP '*** UNABLE TO CONSTRUCT LOG FILE NAME ***'

        OPEN(NLOG,FILE=LOG,STATUS='UNKNOWN',IOSTAT=IER)
        IF (IER .NE. 0) 
     &     STOP '*** UNABLE TO OPEN LOG FILE ***'

C       CREATE NAME FOR THE RESULTS FILE
        CALL FILNAMANDEXT(RESULM,PRJEXC,RESULT,NRESUL,.FALSE.,IER)
        IF (IER .NE. 0) 
     &     STOP '*** UNABLE TO CONSTRUCT RESULTS FILE NAME ***'

C       INCREMENT THE RESULTS FILE VERSION IF EXISTING
        IF (MYPID <= 0) THEN
#ifndef SP_NO_VERSION
           CALL NEXTVERSION(RESULT(1:11),RESULT,NDAT,CXNUM)
#endif

           INQUIRE(FILE=RESULT,EXIST=EX)
           IF (EX) THEN
              OPEN(NDAT,FILE=RESULT,STATUS='OLD',POSITION='APPEND',
     &             IOSTAT=IER)
           ELSE
              OPEN(NDAT,FILE=RESULT,STATUS='REPLACE',IOSTAT=IER)
           ENDIF
           IF (IER .NE. 0) STOP '*** UNABLE TO OPEN RESULTS FILE ***'

C          PRINT OUT HEADING WITH VERSION AND RELEASE DATES
           WRITE(NDAT,9090)
           WRITE(NDAT,9091)
           WRITE(NDAT,9092)
           WRITE(NDAT,9093)CVERS
           WRITE(NDAT,9094)CDATT(1:11),CTIM

           WRITE(NDAT,9096) PRJEXC(1:3),DATEXC(1:3)
9096       FORMAT(/'  PROJECT EXTENSION: ',A3,
     &             '   DATA EXTENSION: ',A3,/)

           WRITE(NOUT,9095)MESG
9095       FORMAT('  Running: ',A)
C          FLUSH RESULTS FILE TO ENSURE THAT IT IS CREATED NOW
           CALL FLUSHRESULTS

           IF (USE_SPIRE) THEN
C             SPIRE IN USE, OPEN SPIRE OUTPUT FILE
              SPIRE_FILE = 'spireout' // RESULT(8:) 
              OPEN(LUNSPIRE,FILE=SPIRE_FILE,STATUS='REPLACE',IOSTAT=IER)

              MESG = ' Running: ' // MESG  
              CALL SPIREOUT(MESG,IRTFLG)     ! Running: SPIDER executable
              CALL SPIREOUT(CVERS,IRTFLG)    ! Version
              MESG = 'DATE:    ' // CDATT(1:11) // '    AT  ' // CTIM  
              CALL SPIREOUT(MESG,IRTFLG)     ! Date and time
              CALL SPIREOUT(RESULT,IRTFLG)   ! Results file name
#ifdef SP_NT
              CALL ERRT(101,'PROCESS ID NOT AVAILABLE IN WINDOWS',NE)
#else
#ifdef SP_GFORTRAN
              IPID =  getpid()
#else
              IPID =  getpid(IPID)
#endif
              WRITE(MESG,9098) IPID
9098          FORMAT('  Current process id: ',I9)
              CALL SPIREOUT(MESG,IRTFLG)     ! Process id
#endif
              CALL FLUSHFILE(LUNSPIRE)
           ENDIF

        ENDIF

C       SKIP OPERATION INPUT IF FIRST OPERATION ON COMMAND LINE
        IF (NUMARG > 1) GOTO 5300
        GOTO 5000

C       @@@@@@@@@@@@@@@@@@@@@@@@@@@@ OPERATION @@@@@@@@@@@@@@@@@@@@@@@@

5000    CONTINUE

C       SEE IF WE MUST COPY LINE TO INTERACTIVE DO LOOP FILE IN RDPR
        LUNDONOW = 0
        IF (COPT == 'I' .AND. NLOOP > 0 .AND. ILOOP == 1) THEN 
           LUNDONOW = LUNDO
           IF (IDOTOP >= 2 ) THEN
              LUNDONOW = LUNDO
              DO I = 1,IDOTOP-1
                 IF (IDOSTK(1,I) > 1) LUNDONOW = 0
              ENDDO
           ELSEIF (IDOTOP == 1) THEN
              IF (NLOOP <= 0) LUNDONOW = 0
           ENDIF
        ENDIF  

C       GET THE NEXT OPERATION, DOES NOT CONVERT INPUT TO UPPERCASE
        SILENT  = .FALSE.
        FCHAR   = NULL
        CALL RDPROP('OPERATION',FCHAR,NALPH,IRTFLG)

C       IF OPERATION IS A COMMENT OR NULL, IGNORE IT.
5100    IF (NALPH  < 1 .OR. 
     &      FCHAR(:1) == ';' .OR. 
     &      FCHAR(:1) == '!') GOTO 5000 

5300    IF (NTRACE > 0 .AND. MYPID <= 0) THEN
           IF (IABSLP .NE. 0 .AND. LOOPREG .NE. 0) THEN
              CALL REG_GET_NAME(LOOPREG,MESG,NLET,IRTFLG)
              WRITE(NSTDOUTP,9039) MESG(1:NLET),IABSLP
9039          FORMAT('  LOOP INDEX (',A,') = ',I8)
           ENDIF
        ENDIF

        !write(6,*) ' idotop: ',idotop, ' at: ',fchar(1:nalph)

C       IF THE FIRST CHARACTER IS '@', GOTO PROCEDURE EVALUATION
        IF (FCHAR(:1) == '@') GOTO 5600

C       IF VARIABLE ASSIGNMENT', SET VARIABLE
        !write(6,*) 'fchar: ',fchar(1:nalph),':', nalph

        ILBRAK = INDEX(FCHAR(1:NALPH),'[')
        IF (ILBRAK > 0) THEN
           IEQ =  INDEX(FCHAR(1:NALPH),'=')
           IF (IEQ > ILBRAK) THEN

C             EQUAL SIGN AFTER BRACKET, CHECK FOR GLOBAL
              GLOBAL = (FCHAR(1:2) == 'GL' .OR. FCHAR(1:2) == 'gl'
     &             .OR. FCHAR(1:2) == 'Gl' .OR. FCHAR(1:2) == 'gL')
 
C             CREATE AND ASSIGN SYMBOLIC (STRING) VARIABLE
              CALL EQU_SYMPAR(FCHAR(1:NALPH),GLOBAL,IRTFLG)
              IF (IRTFLG == 0) GOTO 5000   ! HAVE SET THE VARIABLE
              IF (IRTFLG == 1) GOTO 5000   ! ERROR IN SETTING

C             REGISTER ASSIGNMENTS RETURN: IRTFLG = 2
              IF (IRTFLG == 2) THEN
C                REGISTER ASSIGNMENT, SUBSTITUTE FOR ALL SYM. VAR.
                 CALL SUBSYMPAR(FCHAR(1:NALPH),FCHAR,NALPH,0,IRTFLG)
              ENDIF
           ENDIF
        ENDIF
        
C       IF THE FIRST OR SECOND CHARACTER IS NEITHER A 
C       LETTER NOR A DIGIT, CONSIDER OPERATION AN EXPRESSION
        IF(((.NOT. ISCHAR(FCHAR(1:1)))   .AND.
     &      (.NOT. ISDIGI(FCHAR(1:1))))  .OR.
     &     ((.NOT. ISCHAR(FCHAR(2:2)))   .AND.
     &      (.NOT. ISDIGI(FCHAR(2:2))))) GOTO 6800

C       IF THE FIRST THREE CHARACTERS ARE LETTERS AND THE FORTH IS '('
C       THEN IT MUST BE AN ON-LINE FUNCTION CALL. GOTO EXPRESSION EVAL.
C       FCHAR(5:5) ALLOWS MIS-TYPING SQRT(...) FOR SQR(...)

        IF ((ISCHAR(FCHAR(1:1))) .AND. (ISCHAR(FCHAR(2:2))) .AND. 
     &      (ISCHAR(FCHAR(3:3))) .AND.
     &      (FCHAR(4:4)=='('   .OR. FCHAR(5:5)=='(' )) GOTO 6800

C       IF THE OPERATION STARTS WITH A [] , EVALUATE EXPRESSION.
        IF (FCHAR(1:1) == '[') GOTO 6800

C       CHAR FOLLOWED BY 2 DIGITS IS OLD STYLE BATCH (B01) CALL
        IF (ISCHAR(FCHAR(1:1)) .AND.
     &      ISDIGI(FCHAR(2:2)) .AND. ISDIGI(FCHAR(3:3)) .AND.
     &      NALPH == 3) GOTO 5600

C       OK TO TRANSLATE OPERATION STRING TO UPPER CASE NOW
        CALL SSUPCAS_NOVAR(FCHAR)

C       IF A LABEL 'LB<DIGIT> IS FOUND, AND A DO-LOOP IS IN EFFECT ...
        IF (FCHAR(1:2) == 'LB' .AND. NLOOP >= 1) GOTO 8800
C       IF A ENDDO IS FOUND, AND A DO-LOOP IS IN EFFECT ...
        IF (FCHAR(1:5) == 'ENDDO'  .AND. NLOOP >= 1) GOTO 8800
        IF (FCHAR(1:6) == 'END DO' .AND. NLOOP >= 1) GOTO 8800

C       IF A LABEL IS FOUND, AND NO DO-LOOP IS IN EFFECT ...
        IF (FCHAR(1:2) == 'LB' .AND. NLOOP <= 0) GOTO 5000

C       RESET FFTW3 CACHE 
        CALL FMRS_DEPLAN(IRTFLG)

C       TSWITCH IS MAIN SELECTION PROGRAM FOR OPERATIONS OUTSIDE MAIN
        CALL TSWITCH(IWHICH,ICOM,MAXDIM,IRTFLG)
 
C       IF OPERATION FOUND OUTSIDE OF SPIDER MAIN, GET NEXT OPERATION
        IF (IRTFLG == 0) GOTO 5000
         
C       OPERATION IS NOT IN OUTSIDE MENU. SEARCH LOCAL SUBMENU FOR
C       SPECIFIC LOOPING, IF, EN, ETC OPERATIONS
C       HANDLE LONG LOCAL 'OPERATIONS'

        IF (FCHAR(1:6) == 'ELSEIF' .OR.
     &      FCHAR(1:7) == 'ELSE IF')  GOTO 10795
        IF (FCHAR(1:2) == 'IF'  .OR.
     &      FCHAR(1:2) == 'GO')       GOTO 10800
        IF (FCHAR(1:4) == 'ELSE')     GOTO 10798
        IF (FCHAR(1:5) == 'ENDIF')    GOTO 10799
        IF (FCHAR(1:5) == 'CYCLE')    GOTO 10796
        IF (FCHAR(1:4) == 'EXIT')     GOTO 10797

C       HANDLE SHORT LOCAL OPS: 'EN','DO','LB','RE','MD','OF','IQ VER'
        SELECT CASE(FCHAR(1:2))
          CASE ('EN')
             GOTO 8400
          CASE ('DO')
             GOTO 8600
          CASE ('LB')
             GOTO 8800
          CASE ('RE')
             GOTO 10000
          CASE ('MD')
             GOTO 8500
          CASE ('IQ')
             GOTO 8300
        END SELECT

C       ANY REMAINING OPERATION ASSUMED TO BE ARITHMETIC EXPRESSION
6800    CALL ARASQ(FCHAR,NALPH,GLOBAL,IFLAG)
        IF (IFLAG .NE. 0)THEN
C          EXPRESSION IS NO GOOD - IF IN PROCEDURE, TERMINATES
           MESG = 'UNDEFINED OPERATION OR EXPRESSION: ' // FCHAR
           CALL ERRT(101,MESG,NE)
        ENDIF
        GOTO 5000

C       COMMON TERMINATE ON ERROR SEQUENCE IN PROCEDURE -------------
9999    CALL ENDIT('SHOULD NEVER GET HERE',.FALSE.,RESULT)
        STOP '**** FATAL ERROR'

C       @@@@@@@@@@@@@@@@@@@@@ START PROCEDURE @@@@@@@@@@@@@@@@@@@@@@@@@

5600    CONTINUE

C       SET FLAG FOR PROCEDURE/INTERACTIVE MODE WHEN PROCEDURE WAS CALLED
        FROMBATCH = (COPT == 'B')
        NCHAR     = NALPH
        ISATAT    = (FCHAR(1:2) == '@@')

C       COPY ARGUMENT LIST IF PRESENT -------------------------------
C       IF REGISTERS ARE FOUND, THEN NARGREG > 0
        NCHAR    = NALPH
        ILEFPAR = INDEX(FCHAR(1:NCHAR),'(')
        IF (ILEFPAR > 0 .AND. ISATAT) THEN
           CALL ERRT(101,
     &           'ARGUMENT TRANSFER NOT ALLOWED WITH @@PROCEDURE',NE)
           GOTO 5000

        ELSEIF (ILEFPAR > 0) THEN
C          SAVE REGISTER LIST SENT TO NEW PROC. IN IARGSENT

           CALL REG_GET_SEL(ISTOP,FCHAR(ILEFPAR:NCHAR),.TRUE.,.TRUE.,
     &                      NARGREG,IRTFLG)
           CALL REG_GET_SELS(IARGSENT(1,ISTOP+1),NPARG,NREG,IRTFLG)
C          CUT OFF THE ARGUMENT STRING
           NCHAR = ILEFPAR - 1
        ENDIF

C       TRY TO FIND PROCEDURE FILE LISTED IN THE OPERATION
        CALL GETPROCFILE(FCHAR,NCHAR,PNAME,NPNAME,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 5000

C       INCREMENT PROCEDURE STACK LEVEL COUNTER
        ISTOP  = ISTOP + 1 
        ISTOPR = ISTOP
        IF (ISATAT) ISTOPR = ISTOPR - 1

        IF (ISTOP > MAXPRC) THEN
C          LIMIT IS MAXPRC PROCEDURES !! LET USER KNOW 
           CALL ERRT(102,'PROCEDURE NESTING LEVEL EXCEEDED',MAXPRC)
        ENDIF

C       SAVE DO-LOOP INFO FOR CALLING PROCEDURE IN PROCEDURE'S STACK
        LOOPSV(1,ISTOP) = ILOOP     ! CURRENT LOOP ITERATION
        LOOPSV(2,ISTOP) = IABSLP    ! CURRENT LOOP COUNTER
        LOOPSV(3,ISTOP) = LOOPREG   ! REG. FOR LOOP COUNTER
        LOOPSV(4,ISTOP) = NLOOP     ! NUMBER OF ITERATIONS FOR LOOP
        LOOPSV(6,ISTOP) = IDOTOP    ! LOOP STACK POINTER
        LOOPSV(7,ISTOP) = LBNO      ! LABEL/NDOLINE FOR END OF CURRENT LOOP
        LOOPSV(8,ISTOP) = LOOPINC   ! LOOP COUNTER INCREMENT

C       SAVE THE IFLEVEL VALUES IN THE CALLED PROCEDURE STACK
        IFSV(ISTOP) = IFLEVEL
        IFLEVEL     = 0

C       INITIALIZE VARIABLE INFO FOR NEW PROCEDURE
        CALL SETSYMPAR(NULL,RESPONSE,.TRUE.,IRTFLG)

C       RE-SET DO-LOOP INFO INSIDE CALLED PROCEDURE
        ILOOP   = 1      ! CURRENT LOOP ITERATION
        IABSLP  = 0      ! CURRENT LOOP COUNTER
        NLOOP   = 0      ! NUMBER OF ITERATIONS FOR LOOP
        LOOPREG = 0      ! REG. FOR LOOP COUNTER
        LOOPINC = 1      ! LOOP COUNTER INCREMENT

C       SAVE ALL THE DO-LOOP VALUES IN THE CALLED PROCEDURE STACK
        IDOTOP           = IDOTOP + 1 ! LOOP STACK POINTER
        IDOSTK(1,IDOTOP) = ILOOP      ! CURRENT LOOP ITERATION
        IDOSTK(2,IDOTOP) = IABSLP     ! CURRENT LOOP COUNTER
        IDOSTK(3,IDOTOP) = LOOPREG    ! REG. FOR LOOP COUNTER
        IDOSTK(4,IDOTOP) = NLOOP      ! NUMBER OF ITERATIONS FOR LOOP
        IDOSTK(5,IDOTOP) = ISTOP      ! PROCEDURE STACK LEVEL
        IDOSTK(6,IDOTOP) = LBNO       ! LABEL/NDOLINE FOR END OF CURRENT LOO
        IDOSTK(7,IDOTOP) = LOOPINC    ! LOOP COUNTER INCREMENT

C       PUT IBCNT ON STACK, IF IT ISN'T FIRST CALL, STACK OFFSET
        IPSTACK(ISTOP) = IBCNT

        IF (COPT == 'I') THEN
C         WRITE IT, SINCE NOT ECHOED IN RESULTS FILE IN INTERACTIVE MODE
          IF (MYPID <= 0) THEN
             WRITE(NDAT, 6320)
6320         FORMAT(/,' .OPERATION:')
             WRITE(NDAT, 6340) FCHAR(1:80)
6340         FORMAT(5X,A)
          ENDIF

C         WE'RE NOW IN PROCEDURE MODE, WRITE OUT HEADING.
          COPT = 'B'
        ENDIF

        IF (MYPID <= 0) THEN
           IF (VERBOSE) WRITE(NDAT, *) ' '
           WRITE(NDAT, 6380) PNAME(1:NPNAME)
6380       FORMAT('  -- START OF: ',A,'    --')
        ENDIF 

C       INITIALIZE NEW BANK OF REGISTERS FOR NEW PROCEDURE
        CALL REG_INIT(ISTOPR,IRTFLG)

        IF (NUMARG >= 4) THEN
C          SET INITIAL REGISTERS IN THIS PROCEDURE FROM COMMAND LINE
C          !write(6,*) 'Number of intial arguments: ',numarg
           CALL getarg(4,ARG4)
           NLETA = lnblnk(ARG4)

           DO NARG = 5,NUMARG
C             CONCATENATE FOLLOWING ARGUMENTS (IF ANY) ONTO ARG4
              CALL getarg(NARG,ARGNOW)
              NLETN = lnblnk(ARGNOW)
              ARG4  = ARG4(1:NLETA) // ' ' // ARGNOW(1:NLETN) 
              NLETA = NLETA + NLETN + 1
           ENDDO

           !write(6,*) 'intial arguments: ',ARG4(1:NLETA)
           IFIRST = 1
           DO WHILE (IFIRST < NLETA) 
              CALL GETNEXTTOKEN(ARG4,IFIRST,IGO,IEND)
              IF (IGO <= 0) EXIT
              ARGNOW = ARG4(IGO:IEND)
              IF (ARGNOW(1:1) .NE. '[') THEN
C                NO [] AROUND VARIABLE NAME, ADD IT
                 NTOEQ  = INDEX(ARGNOW,'=') - 2
                 IF (NTOEQ < 0)THEN
C                   EXPRESSION IS NO GOOD - IF BATCH, TERMINATES
                    MESG = 'INVALID INITIAL REGISTER: ' // ARGNOW
                    CALL ERRT(101,MESG,NE)
                    EXIT
                 ENDIF
                 IF (ARG4(IGO:IGO) == 'X' .OR. 
     &               ARG4(IGO:IGO) == 'x') THEN
                    ARGNOW = '[_' // ARG4(IGO+1:IGO+NTOEQ) // ']' //
     &                              ARG4(IGO+NTOEQ+1:IEND)
                 ELSE    
                    ARGNOW = '[' // ARG4(IGO:IGO+NTOEQ) // ']' //
     &                              ARG4(IGO+NTOEQ+1:IEND)
                 ENDIF                  
              ENDIF 
              NLETN = lnblnk(ARGNOW)
              !write(6,*) ' current argument: ',ARGNOW(1:NLETN)

C             SET THE REGISTER, GLOBAL = .FALSE.
              CALL ARASQ(ARGNOW,NLETN,.FALSE.,IFLAG)
              IF (IFLAG .NE. 0)THEN
C                EXPRESSION IS NO GOOD - IF BATCH, TERMINATES
                 MESG = 'INVALID INITIAL REGISTER: '//ARGNOW(:NLETA)
                 CALL ERRT(101,MESG,NE)
                 EXIT
              ENDIF
              IFIRST = IEND + 1
           ENDDO
           NUMARG = 0
        ENDIF 

C       SET ALL THE LUNS & CURRENT PROCEDURE LINE COUNTER
        NIN   = 1
        NOUT  = NDAT
        IBCNT = 0

C       CLOSE THE LOG FILE IF INTERACTIVE,
        IF (NLOG .NE. 0) CLOSE(NLOG)
        NLOG = 0

C       SEE IF PROCEDURE ALREADY LOADED, IF NOT LOAD & LIST
        LISTIT = .TRUE.
        DO NUMPROCNOW = 1, NUMPRC 
           IF (PNAME == PROCFL(NUMPROCNOW)(:MAXNAM)) THEN
C            DON'T NEED TO LOAD PROC FILE OR LIST INTO RESULTS FILE
             LISTIT = .FALSE.
             EXIT
          ENDIF
        ENDDO

        IF (LISTIT) THEN
C          MUST LOAD & LIST

C          OPEN NEW PROCEDURE FILE.
           OPEN(NIN,FILE=PNAME,STATUS='UNKNOWN')

           NUMPRC = NUMPRC + 1
           IF (NUMPRC > MAXPRCNAM) THEN
C             TOO MANY PROCEDURES TO STORE
              CALL ERRT(102,'NO. OF PROCEDURES EXCEEDS MAXPRCNAM',
     &                  MAXPRCNAM)
              GOTO 5000
           ENDIF
           NUMPROCNOW = NUMPRC
           PROCFL(NUMPROCNOW)(:MAXNAM) = PNAME

C          READ IN PROCEDURE LINES & LIST IN RESULTS FILE
           NLINES   = 0
           NCHARS   = 1
           NUMLABS  = 0

           IF (VERBOSE .AND. MYPID <= 0) WRITE(NDAT,*) ' '

           DO 
              READ(NIN,3950,IOSTAT=IERR) PLINE
3950          FORMAT(A)
              IF (IERR .NE. 0) EXIT

              NLINES = NLINES + 1
              NCHAR  = lnblnk(PLINE)

              IF (VERBOSE .and. MYPID <= 0) 
     &           WRITE(NDAT,3960) NLINES, PLINE(:NCHAR)
3960          FORMAT(3X,I4,4X,A)

C             FIND MESG = PLINE WITHOUT WHITESPACE
              CALL SHRINK(PLINE,MESG,NCHARM)

              IF (NLINES == 1 .AND. NCHAR > 0) THEN
C                FIRST LINE. MODERNIZE ANY OLD STYLE ARGUMENTS 
                 IF (NCHARM > 0) THEN
C                   FIRST LINE IS NOT ALL WHITESPACE 

C                   CHECK FOR OLD STYLE ARGUMENTS []
                    CALL CHARINSIDE(MESG,'[',']',.FALSE.,.FALSE.,
     &                              ILEFBRAK,IRITBRAK,NINBRAK)

                    IF (ILEFBRAK == 1 .AND. IRITBRAK == NCHARM) THEN
C                      MODERNIZE OLD STYLE ARGUMENTS 
                       MESG(ILEFBRAK:ILEFBRAK) = '('
                       MESG(IRITBRAK:IRITBRAK) = ')'
                       IF (MYPID <= 0) WRITE(NOUT,*) 
     &                ' *** PLEASE CONVERT PROCEDURE ARGUMENTS TO: (..)'
                    ENDIF   ! (ILEFBRAK == 1 .... 

                    IF (MESG(1:1) == '(') THEN
C                      STRIP WHITESPACE FROM FIRST LINE
                       PLINE = MESG
                       NCHAR = NCHARM 
                    ENDIF   ! (MESG(1:1) == '(')
                 ENDIF      ! (NCHARM > 0) 
              ENDIF         ! (NLINES == 1) 

C             CHECK FOR DUPLICATE LABEL
              IF ((MESG(1:1) == 'L' .OR. MESG(1:1) == 'l') .AND.
     &            (MESG(2:2) == 'B' .OR. MESG(2:2) == 'b') .AND.
     &            ISDIGI(MESG(3:3)) ) THEN

C                OPERATION IS A GOTO LABEL
                 CALL GETLBNO(PLINE,ILBNO,IRTFLG)
                 IF (NUMLABS > 0) THEN
C                   CHECK THRU EXISTING LABELS IN LIST
                    DO I = 1,NUMLABS
                       IF (LABGOT(I) == ILBNO) THEN
C                         ALREADY HAVE THIS LABEL
                          CALL ERRT(102,'DUPLICATE LABEL',ILBNO)
                       ENDIF
                    ENDDO
                 ENDIF
                
C                WRITE(NOUT,*)' PUSHING LABEL: ',ILBNO
                 NUMLABS = NUMLABS + 1
                 IF (NUMLABS <= MAXNUMLAB) THEN
                    LABGOT(NUMLABS) = ILBNO
                 ELSE
                    WRITE(NOUT,*)' --- WARNING TOO MANY LABELS, ',
     &                           ' DUPLICATE LABEL CHECKING ABANDONED'
                    NUMLABS = MAXNUMLAB
                 ENDIF
              ENDIF ! END OF: IF ((MESG(1:1) == 'L' ......

              PLINEGO(NLINES) = NCHARS
              DO I = 1,NCHAR
                 PDATA(NCHARS+I-1) = PLINE(I:I)
              ENDDO
              NCHARS = NCHARS + NCHAR
           ENDDO

           IF (VERBOSE .AND. MYPID <= 0) WRITE(NDAT,*) ' '

           CLOSE(NIN)
C          SAVE PROCEDURE LINES FOR FUTURE USE
           CALL PROC_SET(NUMPRC,NCHARS,NLINES,PLINEGO(1),
     &                   PDATA(1),IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 10000
        ENDIF

C       PUT NEW PROCEDURE NUMBER ON THE STACK
        IPNUMSTACK(ISTOP) = NUMPROCNOW
        IPARNUM(ISTOP)    = 0   !NUMBER FOR OLD IMPLIED <n> SYM. PAR.
        NARGSREC(ISTOP)   = 0

C       CHECK FIRST LINE OF PROCEDURE FOR ARGUMENT TRANSFER
        IBCNT = IBCNT + 1
        CALL PROC_GETPLINE(IBCNT,NUMPROCNOW,PLINE,NUMCHR,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(102,'COULD NOT READ PROC. LINE',IBCNT)
           GOTO 10000
        ENDIF

        IF (PLINE(1:1) == '(') THEN
C          SUBSTITUTE FOR ANY ARGUMENT TRANSFER INTO PROCEDURE
C          @@PROC HAS NONE
           CALL RDPR(PLINE,NCHAR,PLINE,.FALSE.,.FALSE.,.TRUE.,
     &               .FALSE.,.FALSE.,.TRUE.,.TRUE.,IRTFLG)

           IF (PLINE(2:2) == '[') THEN
C             REGISTER ARGUMENT TRANSFER TO PROCEDURE

C             GET REGISTER LIST IN NSELREG
              CALL REG_GET_SEL(ISTOP,PLINE(1:NCHAR),.TRUE.,.FALSE.,
     &                         NARGREG1,IRTFLG)

C             SAVE REGISTER LIST RECEIVED IN IARGSREC
              CALL REG_GET_SELS(IARGSREC(1,ISTOP),NPARG,NREG,IRTFLG)

              IF (NARGREG .NE. NARGREG1) THEN
                WRITE(NOUT,*) 
     &          '*** REGISTER ARGUMENTS SENT TO PROCEDURE:',NARGREG,
     &                        ' <> ARGUMENTS IN PROCEDURE:',NARGREG1
                CALL ERRT(102, 
     &            'WRONG NUMBER OF REGISTERS SENT TO PROCEDURE',NARGREG)
             ENDIF
C            STORE NUMBER OF CURRENT REGISTERS SENT TO THIS PROC.
             NARGSREC(ISTOP) = NARGREG1

C            UPDATE NEW PROCEDURES REGISTERS FROM CALLER'S REGISTERS
             CALL REG_LIST_COPY(NARGREG,
     &                         IARGSENT(1,ISTOP),IARGSREC(1,ISTOP))

           ELSE
             WRITE(NOUT,*) 
     &       '*** REGISTER ARGUMENTS SENT TO PROCEDURE:',NARGREG,
     &       ' BUT NO REGISTER ARGUMENTS IN PROCEDURE'
             CALL ERRT(102, 
     &          'WRONG NUMBER OF REGISTERS SENT TO PROCEDURE',NARGREG)
           ENDIF

        ELSEIF (PLINE(1:1) .NE. '(') THEN

C         NO ARGUMENT TRANSFER, WILL WANT TO REREAD THIS LINE
          IBCNT = IBCNT - 1
        ENDIF

        GOTO 5000


C       @@@@@@@@@@@@@@@@@@ RETURN FROM PROCEDURE @@@@@@@@@@@@@@@@@@@@@@

C       RETURN FROM PROCEDURE--------------------------------------- RE

C       ISTOP TELLS US HOW DEEPLY NESTED WE ARE.  POINTS TO CURRENT
C       TOP OF STACK, WHICH IS CURRENT PROCEDURE FILE.

10000   IF (COPT == 'I') THEN
           CALL ERRT(101,'OPERATION NOT ALLOWED IN INTERACTIVE MODE',N)
           GOTO 5000
        ELSEIF (ISTOP <= 1) THEN
C          TOO MANY RETURNS, HALT IN ERRT
           CALL ERRT(101,'TOO MANY PROCEDURE RETURNS GIVEN',NE)
           GOTO 9999
        ENDIF

C       TAKE PROCEDURE OFF THE STACK, AND CLOSE IT UP

        IF (NARGSREC(ISTOP) > 0) THEN
C          SAVE REGISTER VALUES SPECIFIED IN RECEIVED ARGUMENT LIST
           CALL REG_LIST_COPY(NARGSREC(ISTOP),
     &                        IARGSREC(1,ISTOP),IARGSENT(1,ISTOP))
        ENDIF

C       RETRIEVE DO-LOOP INFO FROM LOWER PROCEDURE LEVEL
        ILOOP    = LOOPSV(1,ISTOP)    ! CURRENT LOOP ITERATION
        IABSLP   = LOOPSV(2,ISTOP)    ! CURRENT LOOP COUNTER
        LOOPREG  = LOOPSV(3,ISTOP)    ! REG. FOR LOOP COUNTER
        NLOOP    = LOOPSV(4,ISTOP)    ! NUMBER OF ITERATIONS FOR LOOP
        IDOTOP   = LOOPSV(6,ISTOP)    ! LOOP STACK POINTER
        LBNO     = LOOPSV(7,ISTOP)    ! CURRENT LOOP LABEL NO.
        LOOPINC  = LOOPSV(8,ISTOP)    ! LOOP COUNTER INCREMENT

C       RETRIEVE IFLEVEL FROM LOWER PROCEDURE LEVEL
        IFLEVEL          = IFSV(ISTOP)

C       RETRIEVE PROCEDURE INFO FROM LOWER PROCEDURE LEVEL
        IBCNT     = IPSTACK(ISTOP)

        ISTOP     = ISTOP - 1
        FROMBATCH = ISTOP > 2
        ISTOPR    = ISTOP             ! REGISTER STACK

C       SET CURRENT LOOPREG IN THIS PROCEDURE.
        IF (LOOPREG > 0)
     &     CALL REG_SET_BYNUM(LOOPREG,REAL(IABSLP),IRTFLG)

C       SIGNAL END OF CURRENT PROCEDURE 
        IF (MYPID <= 0) THEN
            IF (VERBOSE) WRITE(NDAT,*) ' '
            WRITE(NDAT,10080) PNAME(1:NPNAME)
10080       FORMAT('  -- END OF: ',A,'  --')
            IF (VERBOSE) WRITE(NDAT,*) ' '
        ENDIF

C       RETRIEVE CALLER INFO FROM LOWER PROCEDURE LEVEL 
        NUMPROCNOW  = IPNUMSTACK(ISTOP)

        IF (ISTOP > 1) THEN
C          NEW PROCEDURE NAME IS NOW PNAME AT TOP OF STACK.
           PNAME(1:MAXNAM) = PROCFL(NUMPROCNOW)

C          IBCNT OFFSETS INPUT FOR SOLICITATIONS THAT WERE DONE BY
C          PROCEDURE PRIOR TO CALLING A CHILD PROCEDURE
           GOTO 5000
        ENDIF

C       INTERACTIVE MODE.  PUT ALL THE LUN'S BACK, AND
C       REOPEN THE LOG FILE, SINCE IT WAS CLOSED FOR PROCEDURE, 

        COPT  = 'I'
        IBCNT = 0
        NIN   = 5
        NOUT  = NSTDOUTP
        NLOG  = NLOGP

C       USE APPEND FOR LOG FILE, SINCE WE WANT TO ADD ON TO FROM BEFORE
        OPEN(NLOG,FILE=LOG,STATUS='OLD',
     &        ACCESS='SEQUENTIAL',POSITION='APPEND')

C       IF NOT 1ST TIME THROUGH THIS OLD LOOP, READ FROM DO-LOOP FILE
        IF (COPT == 'I' .AND. ILOOP > 1)  NIN = LUNDO

        GOTO 5000

C@@@@@@@@@@@@@@@@@@@ OTHER LOCAL OPERATIONS @@@@@@@@@@@@@@@@@@@@@@@@@@@


C     IQ VAR      ------------------------------------------------- IQ
C            123456789 123456789 123456789 1234567890 
C     CVERS/'VERSION:  UNIX  20.07 ISSUED:  1/30/2012'/ 
8300  READ(CVERS(17:21),*) FVERS
      WRITE(NOUT,8301) FVERS
8301  FORMAT('  SPIDER VERSION: ',F5.2,/)
      CALL REG_SET_NSEL(1,1,FVERS, 0.0,0.0,0.0,0.0,IRTFLG)   
      GOTO 5000


C     SET OPTIONS ------------------------------------------------- MD
8500  CALL SETMODE(RES_TO_TERM)
      GOTO 5000


C     END SPIDER. ------------------------------------------------- EN
C     CLOSE RESULTS & LOG FILE

8400  DELETIT = FCHAR(1:4) == 'EN D'
      CALL SPIREOUT('**** SPIDER NORMAL STOP ****',IRTFLG)
      IF (MYPID <= 0) CALL ENDIT(' COMPLETED',DELETIT,RESULT)

#ifdef USE_MPI
#ifdef MPI_DEBUG
      TT1 = MPI_WTIME()
      IF (MYPID == 0) WRITE(6,8405) TT1-TT0
8405  FORMAT(' TOTAL TIME = ', 1PE11.3)  
#endif
      CALL MPI_FINALIZE(IRC)
#else
#ifdef  HAS_IEEE
#ifndef __APPLE__
C     DO NOT REPORT IEEE INEXACT ....
      call ieee_set_flag( ieee_inexact,       .FALSE. )
      call ieee_set_flag( ieee_denorm,        .FALSE. )
      call ieee_set_flag( ieee_invalid,       .FALSE. )
      call ieee_set_flag( ieee_overflow,      .FALSE. ) ! MAYBE SHUD??
      call ieee_set_flag( ieee_underflow,     .FALSE. ) ! MAYBE SHUD??
      call ieee_set_flag( ieee_divide_by_zero,.FALSE. ) ! MAYBE SHUD??
#endif
#endif
#endif

      STOP ' **** SPIDER NORMAL STOP ****'


C START OF DO LOOP ------------------------------------------------ DO
C ILOOP COUNTS THE NUMBER OF TIMES WE'VE BEEN THRU THE LOOP.
C NLOOP IS THE NUMBER OF ITERATIONS FOR THE LOOP (IF < 1, NOT IN LOOP).
C ILOOP AND NLOOP ARE NEEDED IN THE RDPR* INPUT ROUTINES.
C IABSLP IS THE ACTUAL VALUE OF THE CURRENT ITERATION.

C LBNO HOLDS THE # FROM "DO LB#".
C THIS IS COMPARED WHEN A LABEL "LB#" IS ENCOUNTERED TO SEE IF THE 
C DO-LOOP SHOULD BE UPDATED OR IF THE LABEL SHOULD BE IGNORED.

8600    IF (COPT == 'I' .AND. NLOOP <= 0) THEN
C          NOT IN A LOOP, MUST OPEN NEW DOLOOP SCRATCH FILE
           OPEN(LUNDO,STATUS='SCRATCH',ACCESS='SEQUENTIAL',
     &          FORM='FORMATTED',IOSTAT=IOERR)
           IF (IOERR .NE. 0) THEN
              CALL ERRT(101,'UNABLE TO OPEN TEMP. DO-LOOP FILE',NE)
              GOTO 5000
           ENDIF
C          COPY FIRST DO LINE TO INTERACTIVE DO LOOP FILE
           WRITE(LUNDO,*) FCHAR(1:NALPH)
           NINSAVE = NIN
           NDOLINE = 1
        ENDIF
        ILOOP = 1

C       PRESERVE FCHAR,REMOVE BLANKS FROM FCHAR, ADDS BLANKS AT END 
        FCHARNOBLANK(1:160) = FCHAR(1:160)  
        CALL SHRINKQ(FCHAR,80,FCHAR,NLET)

        IF (NLET == 2) THEN
C          BARE DO LOOP, APPEND SOME DUMMY LOOP PARAMETERS
C                       123456789012345678901234567890
           FCHAR(1:) = 'DO[__DumIndx_]=1,999999999'
           NLET=26           
        ENDIF

        ILOCLIM = INDEX(FCHAR,'=') + 1
        IF (.NOT. ISCHAR(FCHAR(ILOCLIM-2:ILOCLIM-2))) THEN
C          NEW STYLE LOOP INDEX [var], GET REGISTER NUMBER FROM TOKEN
           CALL REG_FIND_IREG('LOC',FCHAR(1:ILOCLIM-2),IDUM,
     &                        LOOPREG,IERR)
           IF (IERR .NE. 0) THEN
               CALL ERRT(101,'CAN NOT PARSE DO LOOP',NDUM)
               GOTO 5000
           ENDIF
        ELSE
C          OLD STYLE LOOP INDEX (K,I,...)
C          CONVERT DO LOOP VARIABLE TO REGISTER (E.G. A --> _A)

           CREG = '[_' // FCHAR(ILOCLIM-2:ILOCLIM-2) // ']'
           CALL REG_FIND_IREG('LOC',CREG,IDUM,LOOPREG,IERR)
           !!CALL REG_GET_VAR(0,CREG,.TRUE.,VALUE,LOOPREG,IENDVAR,IERR)
        ENDIF

C       PARSE OUT LIMITS OF DO LOOP,  (LIMITS CAN BE IN REGISTERS)
C       LOWER LIMIT PUT IN IABSLP, REPITITIONS PUT IN NLOOP
        NC = NLET - ILOCLIM + 1 
        CALL CHKSTR(FCHAR(ILOCLIM:NLET),NC,'I',NUML,FDUM,3,NVAL,IRTFLG)
        IF (IRTFLG .NE. 0 .OR. NVAL < 2) THEN
           CALL ERRT(101,'CAN NOT PARSE DO LOOP',NDUM)
           GOTO 5000
        ENDIF
        IABSLP   = NUML(1)

C       SET LOOP COUNTER INCREMENT
        LOOPINC  = 1
        IF (NVAL > 2)  LOOPINC = NUML(3)
        NLOOP    = (NUML(2) - IABSLP) / LOOPINC + 1
        
C       PUT LOOP START IN REGISTER FOR LOOP COUNTER
        IF (LOOPREG > 0) THEN
           CALL REG_SET_BYNUM(LOOPREG,REAL(IABSLP),IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 5000
        ENDIF

C       PUT IABSLP IN REGISTER [_0] ALSO
        CALL REG_SET(0,REAL(IABSLP),NULL,IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 5000

C       ADD THIS NEW DO LOOP TO THE TOP OF THE LOOP STACK
        IDOTOP = IDOTOP + 1      ! LOOP STACK POINTER

        !write(6,*) ' set idotop: ',idotop,' for:',fchar(1:nalph)

        IF (IDOTOP > MAXPRC) THEN
C          NESTING LEVEL IS MAXPRC, ALWAYS HALTS IN ERRT
           WRITE(NOUT,6171) MAXPRC
6171       FORMAT(' *** LOOP NESTING LEVEL (',I3,') EXCEEDED')
           CALL ERRT(101,'LOOP NESTING LEVEL EXCEEDED',NE)
           GOTO 9999
        ENDIF

C       GET ENDING LABEL NUMBER FOR THIS LOOP
        CALL GETLBNO(FCHAR,LBNO,IRTFLG)
        IF (LBNO < 0) THEN
C          MAKE NEG. LB# FOR PROCEDURE LINE 
           LBNO = -IBCNT
           IF (COPT == 'I') LBNO = -NDOLINE
        ENDIF

C       SAVE ALL THE VALUES IN THE STACKS
        IDOSTK(1,IDOTOP) = ILOOP      ! CURRENT LOOP ITERATION
        IDOSTK(2,IDOTOP) = IABSLP     ! CURRENT LOOP COUNTER
        IDOSTK(3,IDOTOP) = LOOPREG    ! REG. FOR LOOP COUNTER
        IDOSTK(4,IDOTOP) = NLOOP      ! NUMBER OF ITERATIONS FOR LOOP
        IDOSTK(5,IDOTOP) = ISTOP      ! PROCEDURE STACK LEVEL
        IDOSTK(6,IDOTOP) = LBNO       ! LABEL FOR END OF CURRENT LOOP
        IDOSTK(7,IDOTOP) = LOOPINC    ! LOOP INCREMENT

        IF (VERBOSE .and. MYPID <= 0) THEN 
           WRITE(NOUT, 8860) FCHARNOBLANK,IABSLP
        ENDIF
        GOTO 5000

C END OF DO LOOP -------------------------------------------- LB,ENDDO
C       IF # IN LB# IS THE SAME AS CURRENT DO-LOOP, CONTINUE AS
C       USUAL AT THE END OF A DO-LOOP.  OTHERWISE, IGNORE THE LABEL
C       AND GO BACK TO OPERATION AND READ IN THE NEXT LINE.
C       GET LABEL NUMBER FROM THE LINE

8800    IF (FCHAR(1:2) == 'LB') THEN
           CALL GETLBNO(FCHAR,ILBNO,IRTFLG)

C          IF NUMBER IN LB## IS NOT SAME AS CURRENT DO-LOOP, IGNORE 
C          THIS LABEL, GO BACK TO INPUT AND READ IN THE NEXT LINE.
           IF (ILBNO .NE. IDOSTK(6,IDOTOP)) GOTO 5000
        ELSE
           ILBNO = LBNO
        ENDIF

C       MAKE SURE THE LABEL BEING HUNTED FOR IS IN CURRENT PROC.
        IF (ISTOP .NE. IDOSTK(5,IDOTOP)) GOTO 5000

C       NORMAL ENDING OF DO-LOOP
C       INCREASE THE NUMBER OF TIMES WE'VE BEEN THRU THE LOOP AS WELL
C       AS THE ACTUAL ITERATION VALUE OF THE LOOP
        ILOOP  = ILOOP  + 1
        IABSLP = IABSLP + LOOPINC

C       UPDATE REGISTER FOR CURRENT LOOP COUNT X0
        CALL REG_SET(0,REAL(IABSLP),NULL,IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 5000

        IF (ILOOP <= NLOOP) THEN
C         WE'RE NOT DONE WITH THIS LOOP YET
          IF (COPT == 'I' .AND. ILOOP == 2) THEN
C            USE DOLOOP SCRATCH FILE FOR INPUT NOW
             NIN = LUNDO
          ENDIF
          IBCNT = 0
          IF (LOOPREG > 0)
     &          CALL REG_SET_BYNUM(LOOPREG,REAL(IABSLP),IRTFLG)

C         SAVE CHANGED VALUES IN THE CURRENT DO-LOOP STACK
          IDOSTK(1,IDOTOP)    = ILOOP
          IDOSTK(2,IDOTOP)    = IABSLP

C         FIND START OF THE DOLOOP WE'RE WORKING ON, AND CONTINUE
          CALL SEARCHQ(ILBNO,IER)
          IF (IER .NE. 0) THEN
C            SHOULD HALT IN ERRT
             CALL ERRT(102,'END-OF-FILE IN DO-LOOP SEARCH LINE',ILBNO)
             GOTO 9999
          ENDIF

          IF (VERBOSE .AND. MYPID <= 0) 
     &       WRITE(NOUT, 8860) FCHAR(1:20),IABSLP
8860      FORMAT(5X,A20,'    / ',I8)
 
          GOTO 5000
        ENDIF

C       DONE WITH THIS LOOP, RESET VALUES
C       POP THE DO-LOOP STACK, GET ALL THE LOOP VALUES BACK

        IDOTOP = IDOTOP - 1      ! LOOP STACK POINTER

        IF (IDOTOP < 1 .OR. IDOTOP > MAXPRC) THEN
C          NESTING LEVEL IS MAXPRC.  IF EXCEEDED, HALTS IN ERRT
           CALL ERRT(101,'PGM ERROR, LOOP NESTING LEVEL EXCEEDED',NE)
           GOTO 9999
        ENDIF
        ILOOP       = IDOSTK(1,IDOTOP)    ! CURRENT LOOP ITERATION
        IABSLP      = IDOSTK(2,IDOTOP)    ! CURRENT LOOP COUNTER
        LOOPREG     = IDOSTK(3,IDOTOP)    ! REGISTER FOR LOOP COUNTER
        NLOOP       = IDOSTK(4,IDOTOP)    ! NUMBER OF ITERATIONS FOR LOOP
        LBNO        = IDOSTK(6,IDOTOP)    ! LOOP LABEL NUMBER
        LOOPINC     = IDOSTK(7,IDOTOP)    ! LOOP COUNT INCREMENT
        !write(6,*) ' down idotop: ',idotop,lbno, ' for: ',fchar(1:nalph)
 
C       PUT CURRENT IABSLP IN X0 & LOOPREG
        CALL REG_SET(0,REAL(IABSLP),NULL,IRTFLG)
        IF (LOOPREG > 0) THEN
           CALL REG_SET_BYNUM(LOOPREG,REAL(IABSLP),IRTFLG)
        ENDIF

        IF (COPT == 'I' .AND. ILOOP == 1) THEN
C          FIRST TIME THROUGH THIS LOOP, READ FROM TERMINAL NOW
           NIN = NINSAVE

           IF (IDOTOP >= 2 ) THEN
              DO I = 1,IDOTOP-1
C                NOT 1ST TIME THROUGH A HIGHER LOOP, READ DO-LOOP FILE
                 IF (IDOSTK(1,I) > 1) NIN = LUNDO
              ENDDO
           ENDIF
        ENDIF

C       DELETE DOLOOP SCRATCH FILE IF NOT IN NESTED INTERACTIVE DO-LOOP
        IF (COPT == 'I' .AND. IDOTOP == 1) THEN
           CLOSE(LUNDO)
           NDOLINE = 0
        ENDIF

C       GO TO THE TOP OF THE FILE, RESET IBCNT, WHICH TELLS US
C       HOW MANY LINES WE'VE READ IN, AND UPDATE THE STACK VALUES 
        GOTO 5000

C LOGICAL ELSEIF ----------------------------------------------- ELSEIF
10795   IF (IFLEVEL >= 1 .AND. (.NOT. USEELSE(ISTOP,IFLEVEL))) THEN
C          DO NOT NEED TO PROCESS THESE OPERATIONS, SKIP THEM
C          KEEP READING INPUT LINES TILL CORRESPONDING ENDIF FOUND
           CALL FINDENDIF('ENDIF',IFLEVEL,IRTFLG)
           GOTO 5000
        ENDIF
C       CONTINUES INTO 10800

C LOGICAL IF -- OR GOTO---------------------------------------- IF/GOTO
10800   CALL LOGIFQ(FCHAR,LABEL,JUMP,IER)

        IF (IER .NE. 0) THEN
C          ERROR DETECTED, CAN NOT JUMP, (ERRT CALLED IN LOGIFQ)
           GOTO 5000

        ELSEIF (LABEL(1:5) == 'ENDDO' .AND. JUMP) THEN
C          WANT IMMEDIATE EXIT FROM CURRENT LOOP
           GOTO 10797

        ELSEIF (LABEL(1:5) == 'CYCLE' .AND. JUMP) THEN
C             WANT IMMEDIATE CYCLING OF CURRENT LOOP
              GOTO 10796

        ELSEIF (LABEL(1:4) .NE. 'ELSE' .AND. .NOT. JUMP) THEN
C          NO ERROR AND DO NOT WANT TO JUMP, CONTINUE AS NORMAL
           GOTO 5000

        ELSEIF (LABEL(1:1) == ' ') THEN
C          IF FIRST LABEL ELEMENT BLANK, LOGIFQ JUST SETS AN ARITHMETIC 
C          EXPRESSION. NO NEED TO JUMP TO ANY LABEL, CONTINUE AS NORMAL
           GOTO 5000

        ELSEIF (FCHAR(1:6) == 'ELSEIF') THEN
C          THIS IS AN: ELSEIF....THEN....ELSE OPERATION
           USEELSE(ISTOP,IFLEVEL) = JUMP
           IF (JUMP) THEN
C             'IF' IS FALSE, JUMP TO CORRESPONDING ELSE OR ENDIF
C             DECREMENTS IFLEVEL ALSO
              CALL FINDENDIF('ELSE',IFLEVEL,IRTFLG)
           ENDIF
           GOTO 5000
        ENDIF

        IF (LABEL(1:4) == 'ELSE') THEN
C          THIS IS AN: IF....THEN....ELSE OPERATION
           IFLEVEL = IFLEVEL + 1

           IF (IFLEVEL <= 0) THEN
C             USE-ELSE UNDERFLOW, WILL HALT IN PROCEDURE MODE IN ERRT
              CALL ERRT(101,' 2 IN PROGRAM IFLEVEL <= 0',NE)
              GOTO 5000
           ELSEIF (IFLEVEL > MAXPRC) THEN
C             USEELSE OVERFLOW, WILL HALT IN PROCEDURE MODE IN ERRT
              CALL ERRT(101,'IF..ELSE NESTING LEVEL EXCEEDED',NE)
              GOTO 5000
           ENDIF

           USEELSE(ISTOP,IFLEVEL) = JUMP
           IF (JUMP) THEN
C             'IF' IS FALSE, JUMP TO CORRESPONDING ELSE OR ENDIF
C             DECREMENTS IFLEVEL ALSO
              CALL FINDENDIF('ELSE',IFLEVEL,IRTFLG)
           ENDIF
           GOTO 5000
        ENDIF

        IF (.NOT. ISDIGI(LABEL(4:4))) LABEL(4:4) = ' '
C       MAKE SECOND DIGIT OF LB# A BLANK IN THIS CASE

C       CONTINUE INPUT LINES TILL LABEL IS FOUND. NO EFFECT ON IFLEVEL
        CALL FINDLBQ(LABEL,IDOTOP,NLOOP,IDOSTK,NEWLOOP,IFLEVEL,IRT)
C       ERRT CALLED IN FINDLBQ NOW, SO IT SHOULD HALT THERE

        IF (NEWLOOP < 0) THEN
C          HAVE PASSED CURRENT DO-LOOP END, MUST POP DO-LOOP STACK
           IF (COPT == 'I') THEN
              CALL ERRT(13,'SPIDER',NE)
              GOTO 5000
           ENDIF

           IDOTOP      = IDOTOP + NEWLOOP   ! LOOP STACK POINTER

           IF (IDOTOP > MAXPRC .OR. IDOTOP <= 0) THEN
C             NESTING LEVEL IS MAXPRC.  IF EXCEEDED, HALTS IN ERRT
              CALL ERRT(101,'PGM ERROR 2, LOOP NESTING LEVEL 2',NE)
              GOTO 9999
           ENDIF
           ILOOP       = IDOSTK(1,IDOTOP)   ! CURRENT LOOP ITERATION
           IABSLP      = IDOSTK(2,IDOTOP)   ! CURRENT LOOP COUNTER
           LOOPREG     = IDOSTK(3,IDOTOP)   ! REG. FOR LOOP COUNTER
           NLOOP       = IDOSTK(4,IDOTOP)   ! NUMBER OF ITERATIONS FOR LOOP
           LBNO        = IDOSTK(6,IDOTOP)   ! LOOP LABEL NUMBER
           LOOPINC     = IDOSTK(7,IDOTOP)   ! LOOP COUNTER INCREMENT
 
C          PUT CURRENT IABSLP IN X0 & LOOPREG
           CALL REG_SET(0,REAL(IABSLP),NULL,IRTFLG)
           IF (LOOPREG > 0) THEN
              CALL REG_SET_BYNUM(LOOPREG,REAL(IABSLP),IRTFLG)
           ENDIF
        ENDIF

C       PROCESS THE OPERATION RETURNED FROM FINDLBQ
        GOTO 5100

C       LOOP CYCLE ---------------------------------------------- CYCLE
C       CONTINUE INPUT LINES TILL ENDDO IS FOUND. NO EFFECT ON IFLEVEL

10796   CALL FINDLBQ('ENDDO',IDOTOP,NLOOP,IDOSTK,NEWLOOP,IFLEVEL,IRT)

C       PROCESS THE 'ENDDO' OPERATION RETURNED FROM FINDLBQ
        GOTO 5100

C       LOOP EXIT ------------------------------------------------ EXIT
C       CONTINUE INPUT LINES TILL ENDDO IS FOUND. NO EFFECT ON IFLEVEL

10797   CALL FINDLBQ('ENDDO',IDOTOP,NLOOP,IDOSTK,NEWLOOP,IFLEVEL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 8400 

C       MAKE THIS LAST ITERATION OF THE CURRENT LOOP
        ILOOP = NLOOP
        !write(6,*) ' iloop,newloop: ',iloop,newloop

C       PROCESS THE 'ENDDO' OPERATION RETURNED FROM FINDLBQ
        GOTO 5100


C LOGICAL ELSE -------------------------------------------------- ELSE
10798   CONTINUE
        IF (IFLEVEL <= 0) THEN
           CALL ERRT(101,' IN PROGRAM IFLEVEL <= 0',NE)
        ELSEIF (.NOT. USEELSE(ISTOP,IFLEVEL)) THEN
C          DO NOT NEED TO PROCESS THESE OPERATIONS, SKIP THEM
C          KEEP READING INPUT LINES TILL CORRESPONDING ENDIF FOUND
           CALL FINDENDIF('ENDIF',IFLEVEL,IRTFLG)

        ELSE IF (USEELSE(ISTOP,IFLEVEL)) THEN
C          MUST USE THIS ELSE CLAUSE
           CONTINUE
        ENDIF
        GOTO 5000

C LOGICAL ENDIF ------------------------------------------------ ENDIF
10799   CONTINUE
        IFLEVEL = IFLEVEL - 1
        GOTO 5000
      
        END



C **************************************************************************
C
C    SHRINK(INSTR,OUTSTR,LENOUT)

C    PURPOSE:  SHRINK STRING BY IGNORING ALL NON-PRINTING CHARACTERS
C
C    PARAMETERS: INSTR      INPUT STRING TO BE SHRANK         SENT
C                OUTSTR     OUPUT SHRUNKEN STRING             RET.
C                LENOUT     LENGTH OF SHRUNKEN STRING         RET.
C
C *************************************************************************

        SUBROUTINE SHRINK(INSTR,OUTSTR,LENOUT)

        CHARACTER *(*) INSTR,OUTSTR

        LENIN  = LEN(INSTR)
        LENMAX = LEN(OUTSTR)

        LENOUT = 0
        DO I=1,LENIN
           IF (INSTR(I:I) >= '!' .AND. INSTR <=  '~') THEN
              LENOUT = LENOUT + 1
              OUTSTR(LENOUT:LENOUT) = INSTR(I:I)
           ENDIF
        ENDDO

        IF (LENOUT. LT. LENMAX) THEN
C          PUT BLANKS AT END OF OUTSTR
           OUTSTR(LENOUT+1:LENMAX) = ' '
        ENDIF

        RETURN
        END

C **************************************************************************
C
C    EQU_SYMPAR(LINE,GLOBAL,IRTFLG)
  
C    PURPOSE: CREATE A SYMBOLIC (STRING) VARIABLE FROM COMMAND LINE
C             ASSIGNMENT
C
C    PARAMETERS:   LINE       INPUT STRING 
C                  GLOBAL     FLAG FOR GLOBAL VARIABLE
C                  IRTFLG     ERROR FLAG (2 IF NOT FOR STRING VAR)
C
C *************************************************************************

      SUBROUTINE EQU_SYMPAR(LINE,GLOBAL,IRTFLG)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER (LEN=*)   :: LINE
      CHARACTER (LEN=160) :: SYMPARID,SYMPARVAL,MSG
      CHARACTER(LEN=1)    :: NQ1,NQ2,CIGO,CEND
      LOGICAL             :: GLOBAL,LOCAL,ISREGVAR

C     FOR ISTOP 
      INTEGER, DIMENSION(MAXPRC) :: IPSTACK,IPNUMSTACK,IPARNUM
      COMMON /QSTR_STUFF1/ ISTOP,ITI,ITIN,IWHERE,IPSTACK,
     &                     IPNUMSTACK,IPARNUM

      CALL SET_MPI(ICOMM,MYPID,MPIERR)

      NQ1       = CHAR(39)   ! '
      NQ2       = CHAR(34)   ! "

      IRTFLG    = 1
      LENT      = lnblnkn(LINE)

C     LOCATE FIRST VARIABLE NAME IN LINE (SHOULD'VE ALEADY BEEN TESTED)

      CALL CHARINSIDE(LINE,'[',']',.TRUE.,.FALSE.,IP1,IP2,NCHARID)
      IF (NCHARID <= 0) THEN
         LENT = lnblnkn(LINE)
         MSG = 'NO VARIABLE NAME ([NAME]) IN: ' // LINE(1:LENT)
         CALL ERRT(101,MSG,NE)
         RETURN
      ENDIF

C     GET  VARIABLE NAME FROM LINE 
      SYMPARID = '<' // LINE(IP1:IP2) // '>'  
      NCHARI   =  IP2 - IP1 + 3    
      !write(6,*) ' Got symparid:',symparid(:nchari),':',nchari

      CALL REG_FIND(0,SYMPARID(:NCHARI),VALDUM,IREG,IRTFLG)
      !write(6,*) ' Queried reg var:', SYMPARID(1:NCHARI), ireg,valdum

C     SEE IF TRYING TO SET A REGISTER VARIABLE
      IF (IREG > 0 .OR. SYMPARID(2:2) == '_' ) THEN
         IRTFLG = 2
         RETURN
      ENDIF 

C     GET ASSIGNED VARIABLE VALUE FROM LINE, MAY BE AN EXPRESSION
      NEQ    = 0         ! = SIGN COUNTER
      NCHARV = 0         ! = CHARACTER COUNTER
      IFIRST = IP2 + 2   ! = STARTING LOCATION
      DO 
         CALL GETNEXTTOKEN2(LINE,IFIRST,IGO,IEND) 
         IF (IEND <= 0) EXIT     ! NO MORE TOKENS

C        FOUND A TOKEN, IT SHOULD BE A: =, SYM. VARIABLE, OR QUOTED STRING
         IFIRST = IEND + 1         ! NEXT START FOR TOKEN SEARCH
         CIGO   = LINE(IGO:IGO)
         CEND   = LINE(IEND:IEND)

         IF (CIGO == '=' ) THEN
C           TOKEN IS AN EQUAL SIGN             
            IF (NEQ > 0) THEN
               MSG = 'EXTRA = SYMBOL IN: '// LINE
               CALL ERRT(101,MSG,NE)
               RETURN
            ENDIF
            NEQ = NEQ + 1 

         ELSEIF ((CIGO == NQ1 .AND. CEND == NQ1) .OR.
     &           (CIGO == NQ2 .AND. CEND == NQ2)) THEN
C           TOKEN IS QUOTED TEXT STRING             
            SYMPARVAL(NCHARV+1:) = LINE(IGO+1:IEND-1)
            NCHARV               = NCHARV + (IEND - IGO - 1)

         ELSEIF (CIGO == '[' .AND.CEND == ']') THEN

            !write(6,*)' Calling issympar:',line(igo:iend),':'
C           TOKEN IS A [] VARIABLE, IS IT A SYM. STRING VARIABLE?
            CALL ISSYMPAR(LINE(IGO:IEND),-1,ICVAR,IRTFLG)
            !write(6,*)' Issympar:',line(igo:iend),':',ICVAR,IRTFLG

            IF (ICVAR <= 0 .OR. IRTFLG .NE. 0) THEN
C              RIGHT SIDE NOT A SYMVAR, MAY BE REG. ASSIGNMENT INSTEAD?
               !write(6,*)' Not sympar:',line(igo:iend),':',ICVAR,IRTFLG
           
               IRTFLG = 2
               RETURN
            ENDIF

C           UNQUOTED SYM [] VARIABLE, SUBSTITUTE VALUE FOR IT 
            !write(6,*)' Sub:',line(igo:iend),':',symparval(ncharv+1:),':'
            CALL SYMPAR_SUB(LINE(IGO:IEND),SYMPARVAL(NCHARV+1:),
     &                      NCHARS,ISTOP,.TRUE.,IRTFLG)
            !write(6,*)' Sub, nchars,irtflg:',nchars,irtflg

            NCHARV = NCHARV + NCHARS
         ENDIF
      ENDDO

      IF (NCHARV < 1) THEN
C        MAY BE A REGISTER ASSIGNMENT INSTEAD?
         IRTFLG = 2
         RETURN
      ENDIF
       
      !write(6,*) ' Symparval:',symparval(:ncharv),':',ncharv

      LOCAL = .NOT. GLOBAL
      CALL SETSYMPAR(SYMPARID(:NCHARI),SYMPARVAL(:NCHARV),LOCAL,IRTFLG)

      IRTFLG = 0
      END




C       -------------  GETPROCFILE ----------------------- GETPROCFILE

C       GET PROCEDURE FILE NAME FROM STR LINE

        SUBROUTINE GETPROCFILE(STR,NCHAR,PNAME,NPNAME,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=*)      :: STR,PNAME
        CHARACTER(LEN=MAXNAM) :: PNAMEM
        CHARACTER(LEN=160)    :: MESG
        LOGICAL               :: EX

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C       FIND FIRST CHARACTER IN FILENAME
        NFSTRT = INDEX(STR(1:NCHAR),'@',BACK=.TRUE.) + 1

C       SUBSTITUTE FOR ALL SYM. VAR. IN STR
        !write(6,*) ' Checking: ',STR(NFSTRT:NCHAR)
        CALL SUBSYMPAR(STR(NFSTRT:NCHAR),PNAMEM,NCT,0,IRTFLG)

C       WILL STOP IN FILNAMANDEXT ON ERRT
        CALL FILNAMANDEXT(PNAMEM(:NCT),PRJEXC,PNAME,NPNAME,.TRUE.,IER)

C       TRY TO FIND PROCEDURE IN USER'S DIRECTORY & PROJECT EXTENSION

        !write(6,*) ' Checking: ',pname(NFSTRT:NCHAR)
        IF (MYPID <= 0) INQUIRE(FILE=PNAME,EXIST=EX)
        CALL BCAST_MPI('SPIDER','EX',EX,1,'L',ICOMM)

        IF (.NOT. EX) THEN
C         PROCEDURE FILE DOESN'T EXIST IN USER'S DIRECTORY
          IF (MYPID <= 0)  WRITE(NOUT,90) PNAME 
90        FORMAT(' NO LOCAL PROCEDURE FILE: ',A) 

C         TRY AGAIN UNDER 'PROC:*.spi' IN PROC DIR.
C         7/10/88 PROC IS ENV. VAR. FOR DIR. WHERE *.spi FILES ARE al

          CALL MYGETENV('SPPROC_DIR',PNAME,NCHART,
     &                 'dir-for-proc-files',IER)
          IF (IER .NE. 0) CALL ERRT(101,'NO ENVIRONMENT VARIABLE',NE)
          NCHARTN = LNBLNKN(PNAMEM)
          PNAME   = PNAME(:NCHART) // PNAMEM(:NCHARTN) // '.spi'
          NPNAME  = NCHART + NCHARTN + 4

          IF (MYPID <= 0) INQUIRE(FILE=PNAME,EXIST=EX)
          CALL BCAST_MPI('SPIDER','EX',EX,1,'L',ICOMM)

          IF (.NOT. EX) THEN
C            THE *.spi FILE DOES NOT EXIST. NOTIFY USER
             MESG = 'PROCEDURE FILE.spi DOES NOT EXIST: '//
     &              PNAME(:NPNAME) 
             CALL ERRT(101,MESG,NE)
             IRTFLG = 1
             RETURN
          ENDIF
        ENDIF

        IRTFLG = 0
        END




#ifdef MPI_DEBUG

        SUBROUTINE PI3F

        DOUBLE PRECISION, PARAMETER ::
     &                    PI25DT = 3.141592653589793238462643D0

        DOUBLE PRECISION :: DMYPI, PI, H, SUM, X, F, A
        INTEGER          :: RC
        ! FUNCTION TO INTEGRATE
        F(A) = 4.D0 / (1.D0 + A*A)

        INCLUDE 'mpif.h'

        ICOMM = MPI_COMM_WORLD
        CALL MPI_COMM_RANK(ICOMM, MYPID,  MPIERR)
        CALL MPI_COMM_SIZE(ICOMM, NPROCS, MPIERR)

        IF(MYPID == 0)WRITE(6,*) ' NPROCS = ', NPROCS,' icomm:',icomm
        SIZETYPE   = 1
        SUMTYPE    = 2
 
        DO j= 35,-1,-7
            IF ( MYPID == 0 ) THEN
               N = j
c              CALL RDPRI1S( N,notused,
c     &           'ENTER THE NUMBER OF INTERVALS: (0 QUITS)',irtflg)
c 99           FORMAT(I10)

              WRITE(6,*) 'INTERVALS: ',N
            ENDIF
      
c              WRITE(6,*) ' bcasting: ',N
           CALL MPI_BCAST(N,1,MPI_INTEGER,0,ICOMM,MPIERR)

           ! CHECK FOR QUIT SIGNAL
           IF ( N <= 0 ) EXIT

           ! CALCULATE THE INTERVAL SIZE
           H = 1.0D0 / N
 
           SUM  = 0.0D0
           DO I = MYPID+1, N, NPROCS
              X   = H * (DBLE(I) - 0.5D0)
              SUM = SUM + F(X)
           ENDDO
           DMYPI = H * SUM

C          COLLECT THE PARTIAL SUMS
c              WRITE(6,*) ' reducing: ',N
           CALL MPI_REDUCE(DMYPI,PI,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, 
     &                  ICOMM,MPIERR)

C          NODE 0 PRINTS THE ANSWER.
           IF (MYPID == 0) THEN
              WRITE(6, 97) PI, ABS(PI - PI25DT)
 97           FORMAT('  Pi is approximately: ', F10.8,
     &               '  Error is: ', F10.8)
           ENDIF
        ENDDO

        END

#endif
