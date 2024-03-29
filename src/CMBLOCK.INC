
C  LEGACY COMMON BLOCKS USED IN MANY PROGRAMS FOR COMMUNICATION

C  CMBLOCK.INC  CHANGED 9/21/79,  1/1/81, 9/27/82, 11/16/87,   4/9/91,
C                       10/2/91, 2/23/92, 3/10/92, 12/28/93, 12/30/93,
C                       17/4/95, 3/4/98,  9/14/00,  4/17/03,  6/24/05,
C                       6/24/06, 12/11/07,6/17/11,  4/02/12, 11/18/19

C...    BATCH    COMMON BLOCK FOR BATCH MODE OPERATIONS ---------------

        COMMON    /BATCH/IBCNT,NLOOP,IABSLP,ILOOP
        INTEGER  ::      IBCNT,NLOOP,IABSLP,ILOOP

C       IBCNT    LINE COUNTER IN BATCH MODE
C       NLOOP    TOTAL NO. OF LOOPS TO BE EXECUTED IN CURRENT DO LOOP
C       IABSLP   EXACT CURRENT DO LOOP COUNTER
C       ILOOP    OFFSET OF CURRENT DO LOOP COUNTER FROM STARTING 
C                LOOP COUNT + 1

        COMMON     /BATCH2/ COPT
        CHARACTER(LEN=1) :: COPT

C       COPT     = 'I' FOR INTERACTIVE MODE, = 'B' FOR BATCH MODE

C...    UNITS    COMMON BLOCK FOR UNIT NOS. ---------------------------

        COMMON /UNITS/LUNC,NIN,NOUT,NECHO,IFOUND,NLOG,NDAT
        INTEGER    :: LUNC,NIN,NOUT,NECHO,IFOUND,NLOG,NDAT

C       LUNC     (CURRENTLY UNUSED)
C       NIN      UNIT NO. FOR COMMAND INPUT STREAM
C       NOUT     UNIT NO. FOR OUTPUT 
C                (IN BATCH MODE IT IS ATTACHED TO RESULTS FILE)
C       NECHO    (CURRENTLY UNUSED) 
C       IFOUND   CONTROL VARIABLE FOR PRINTING OF MESSAGES BY OPEN
C       NLOG     UNIT NO. FOR LOG FILE
C       NDAT     UNIT NO. FOR RESULTS FILE

C...    IPRTT    COMMON BLOCK FOR INPUT/OUTPUT CONTROL ----------------
        LOGICAL  :: SILENT,VERBOSE,USE_SPIRE
        INTEGER  :: NUMFFTWTH,NTRACE,NALPH
        COMMON /IPRTT/NUMFFTWTH,NTRACE,NALPH,VERBOSE,USE_SPIRE,SILENT

        LOGICAL     :: LEGACYPAR,USE_FBP_INTERP,USE_FBS_INTERP,
     &                 USE_LONGCOL,IN_PARALLEL
        COMMON /IPRTT2/LEGACYPAR,USE_FBP_INTERP,USE_FBS_INTERP,
     &                 USE_LONGCOL,IN_PARALLEL

C       NTRACE     ENABLE/DISABLE TRACING ON TERMINAL DURING BATCH MODE
C       NALPH      NO. OF CHARACTERS IN THE MOST RECENT OPERATION 
C       VERBOSE    USER SET CONTROL VARIABLE FOR PRINTING OF INFO
C       USE_SPIRE  USING SPIRE 
C       SILENT     OPERATION SET CONTROL VARIABLE FOR PRINTING OF INFO
C       LEGACYPAR  NEEDS () AROUND SIMPLE LIST INSIDE LOOPED INPUT
C       NUMFFTWTH  NUMBER OF FFTW3 THREADS

C...    MASTER   COMMON BLOCK FOR LABEL OF MOST RECENT FILE OPENED ----

        COMMON /MASTER/ NSAMC,NROWC,IREC,NHISTREC,IFORM,IMAMI,
     &                  FMAX,FMIN,AV,SIG,IHIST,IPROT

        INTEGER  ::     NSAMC,NROWC,IREC,NHISTREC,IFORM,IMAMI,IHIST
        INTEGER  ::     IPROT
        REAL     ::     FMAX,FMIN,AV,SIG

C       NSAMC    NO. OF SAMPLES (COLUMNS) PER ROW
C       NROWC    NO. OF ROWS PER SLICE
C       IREC     TOTAL NO. OF RECORDS IN THE FILE
C       NHISTREC (CURRENTLY UNUSED) 
C       IFORM    FILE TYPE SPECIFIER
C       IMAMI    MAXIMUM/MINIMUM FLAG. SET AT 0 WHEN FILE IS CREATED,
C                SET AT 1 WHEN THE MAXIMUM AND MINIMUM HAVE BEEN 
C                COMPUTED AND STORED INTO THE HEADER RECORD.
C       FMAX     MAXIMUM VALUE
C       FMIN     MINIMUM VALUE
C       AV       AVERAGE VALUE
C       SIG      STANDARD DEVIATION, A -1 INDICATES SIG NOT CALCULATED  
C       IHIST    (CURRENTLY UNUSED)
C       IPROT    (CURRENTLY UNUSED)


C...    MASTER1  COMMON BLOCK FOR LABEL INFORMATION OF MOST RECENT FILE

        COMMON /MASTER1/      DATEXC,PRJEXC,CDAT,CTIM,CTIT
        CHARACTER(LEN=4)   :: DATEXC,PRJEXC
        CHARACTER(LEN=10)  :: CDAT
        CHARACTER(LEN=8)   :: CTIM
        CHARACTER(LEN=160) :: CTIT

C       DATEXC   CHARACTER VARIABLE CONTAINING DATA EXTENSION CODE
C       PRJEXC   CHARACTER VARIABLE CONTAINING PROJECT EXTENSION CODE
C       CDAT     CHARACTER VARIABLE CONTAINING DATE OF FILE
C                CDAT IS UNUSED NOW AS LENGTH SHOULD BE 12 NOW (Nov 99)
C       CTIM     CHARACTER VARIABLE CONTAINING TIME OF FILE
C       CTIT     CHARACTER VARIABLE CONTAINING TITLE OF FILE


C...    PROC     COMMON BLOCK FOR PROCEDURE OPERATION -----------------
C                FROMBATCH USED IN errt.f, rdpra.f, setmode.f, 
C                                  spider.f, sympar.f
        LOGICAL     :: FROMBATCH
        COMMON  /PROC/ FROMBATCH

C...    FUNCTION   COMMON BLOCK FOR CURRENT OPERATION -----------------

C       FCHAR IS USED TO READ OPERATION CODES ETC. FROM INPUT STREAM.
        CHARACTER(LEN=160) :: FCHAR
        COMMON      /FUNCTION/FCHAR

C...    MRC_ORIGIN   COMMON BLOCK FOR CURRENT DEFAULT MRC ORIGIN -----

        CHARACTER(LEN=4)  :: MRC_AXIS
        COMMON   /MRC_ORIGIN/MRC_AXIS



