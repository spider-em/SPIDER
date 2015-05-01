C++*********************************************************************
C
C    APSH_PS.F      USED CMLIMIT                  AUG 00 ARDEAN LEITH
C                   ADDED REF_CIRC FILE           APR 01 ARDEAN LEITH
C                   NORMASS -> NORMAS             OCT 01 ARDEAN LEITH
C                   PROMPTS                       JAN 02 ARDEAN LEITH
C                   OPFILEC                       FEB 03 ARDEAN LEITH
C                   APMASTER REWRITE              AUG 03 ARDEAN LEITH
C                   USED AP_OUT                   FEB 04 ARDEAN LEITH
C                   NORMASS USED FOR ALTIX        JUN 04 ARDEAN LEITH
C                   INTERPOLATION ON EDGE BUG     AUG 04 ARDEAN LEITH
C                   CROSRNG_E SPEEDS UP           AUG 04 ARDEAN LEITH
C                   AP_END CALL HAS PARLIST       OCT 04 ARDEAN LEITH
C                   PEAKV = 1                     JAN 05 ARDEAN LEITH
C                   RANGNEWT OPT 64 BUG FIXED     MAY 05 ARDEAN LEITH
C                   DISCARD MIRROR ...            JUN 06 ARDEAN LEITH
C                   AP_STAT CALL                  JAN 07 ARDEAN LEITH
C                   REF-RINGS FILE ADDED          MAR 08 ARDEAN LEITH
C                   CROSRNG_E REWRITE             MAR 08 ARDEAN LEITH
C                   IREFLIST VAR. NAME            MAR 08 ARDEAN LEITH
C                   APRINGS_ONE_NEW               MAR 08 ARDEAN LEITH
C                   AP_END CALL CHANGED           NOV 08 ARDEAN LEITH
C                   AP_STAT_ADD                   NOV 08 ARDEAN LEITH
C                   NPROJ BUG                     AUG 09 ARDEAN LEITH
C                   ISHRANGEX                     FEB 10 ARDEAN LEITH
C                   LENTT, TT REMOVED             JUN 10 ARDEAN LEITH
C                   CROSRNG_2                     JUN 10 ARDEAN LEITH
C                   AP_STAT NBORDER               OCT 10 ARDEAN LEITH
C                   RENAMED FROM MRQLI_PS         JAN 11 ARDEAN LEITH
C                   AP_GETDATA                    DEC 11 ARDEAN LEITH
C                   ROTFIRST                      DEC 11 ARDEAN LEITH
C                   FBS_WANTED                    JAN 12 ARDEAN LEITH
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C APSH_PS 
C  
C PURPOSE: FIND ROTATIONAL AND SHIFT PARAMETERS TO ALIGN A SERIES OF
C           REFERENCE IMAGES WITH SAMPLE IMAGES
C           USED IF:  CIRCREF_IN_CORE .AND. NUMEXP >= NUMTH 
C
C PARAMETERS:
C       IREFLIST            LIST OF REF. IMAGE FILE NUMBERS      (SENT)
C       NUMREF              NO. OF IMAGES                        (SENT)
C       IEXPLIST            LIST OF EXP. IMAGE FILE NUMBERS      (SENT)
C       NUMEXP              NO. OF IMAGES                        (SENT)
C       NX,NY           ACTUAL (NOT-WINDOWED) IMAGE SIZE     (SENT)
C       ISHRANGEX,ISHRANGEY ALLOWED SHIFT RANGE                  (SENT)
C       ISTEP               SHIFT STEP WITHIN RANGE              (SENT)
C       NRING                                                    (SENT)
C       LCIRC                                                    (SENT)
C       NUMR                                                     (SENT)
C       CIRCREF                                                  (SENT)
C       CIRCREF_IN_CORE                                          (SENT)
C       MODE                                                     (SENT)
C       REFANGDOC           REF. ANGLES FILE NAME                (SENT)
C       EXPANGDOC           EXP. ANGLES FILE NAME                (SENT)
C       SCRFILE                                                  (SENT)
C       FFTW_PLANS                                               (SENT)
C       REFPAT              REF. IMAGE SERIES FILE TEMPLATE      (SENT)
C       EXPPAT              EXP. IMAGE SERIES FILE TEMPLATE      (SENT)
C       RANGE                                                    (SENT)
C       ROTFIRST            USE RTSQ ON EXP INPUT IMAGES         (SENT)
C       CKMIRROR            LOGICAL FLAG TO CHECK MIRRORING      (SENT)
C       CTYPE                                                    (SENT)
C       LUNDOC                                                   (SENT)
C
C
C NOTE:  IF USING APSH_PS, MOST MEMORY DEMAND APPEARS TO BE DEPENDENT 
C        ON LCIRC & NUMREF.  LCIRC IS THE TOTAL LENGTH OF THE ARRAY
C        THAT HOLDS THE CIRCULAR RINGS, SO IT IS DEPENDENT ON
C        NUMBER OF RINGS AND THEIR RADIUS.  NUMREF IS NUMBER OF REFERENCE
C        IMAGES. BIGGEST ARRAY ALLOCATED IS: CIRCREF(LCIRC,NUMREF)
C        ANOTHER SMALL ALLOCATED ARRAY IS: A(NX,NY,NUMTH)
C        FOR 83 IMAGES of 125x125 WITH RINGS AT 5...47 SIZE=45MB
C        ARRAYS ONLY APPEAR TO TAKE: 3.6MB?
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************


C ---------------------  NON-MPI CODE --------------------------------
#ifndef USE_MPI

       SUBROUTINE APSH_PS(IREFLIST,NUMREF,IEXPLIST,NUMEXP, 
     &               NX,NY,ISHRANGEX,ISHRANGEY,ISTEP,
     &               NRING,LCIRC,NUMR,CIRCREF,CIRCREF_IN_CORE,
     &               MODE, REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS,
     &               REFPAT,EXPPAT,RANGE,ROTFIRST,
     &               CKMIRROR,CTYPE,LUNDOC,FBS_WANTED)

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

	INTEGER                :: IREFLIST(NUMREF) 
	INTEGER                :: NUMREF 
	INTEGER                :: IEXPLIST(NUMEXP) 
	INTEGER                :: NUMEXP 
        INTEGER                :: NX,NY,ISHRANGEX,ISHRANGEY,ISTEP
        INTEGER                :: NRING,LCIRC
        INTEGER                :: NUMR(3,NRING)
	REAL                   :: CIRCREF(LCIRC,NUMREF)
	LOGICAL                :: CIRCREF_IN_CORE
        CHARACTER (LEN=1)      :: MODE
        CHARACTER (LEN=*)      :: REFANGDOC,EXPANGDOC
        CHARACTER (LEN=*)      :: SCRFILE
        INTEGER *8             :: FFTW_PLANS(*)
        CHARACTER (LEN=*)      :: REFPAT,EXPPAT 
        REAL                   :: RANGE
	LOGICAL                :: ROTFIRST,CKMIRROR
        CHARACTER (LEN=*)      :: CTYPE
        INTEGER                :: LUNDOC
        LOGICAL                :: FBS_WANTED


C       AUTOMATIC ARRAYS
	REAL                   :: ANGOUT(3)

C       ALLOCATED ARRAYS
	REAL, ALLOCATABLE      :: EXPBUF(:,:,:)
	INTEGER, ALLOCATABLE   :: NPROJA(:) 
	REAL, ALLOCATABLE      :: DLIST(:,:) 
	REAL, ALLOCATABLE      :: REFDIR(:,:),EXPDIR(:,:) 
	REAL, ALLOCATABLE      :: ANGREF(:,:),ANGEXP(:,:)
	REAL, ALLOCATABLE      :: TMPBUF(:,:)

        INTEGER, PARAMETER     :: NLISTMAX = 15
        REAL                   :: PARLIST(NLISTMAX)

        REAL                   :: ADUM
        LOGICAL                :: ANGINHEADER  
	LOGICAL                :: GOTREFANG,LIMITRANGE
	LOGICAL                :: MIRRORNEW

        INTEGER                :: NBORDER = 0      ! # BORDER PIXELS
        INTEGER                :: NSUBPIX = 0      ! # SUBPIX PIXELS

        INTEGER, PARAMETER     :: LUNPIC  = 77
        INTEGER, PARAMETER     :: LUNANG   = 78
        INTEGER, PARAMETER     :: LUNRING = 50

        LOGICAL, PARAMETER     :: MPIBCAST = .FALSE.
	REAL, PARAMETER        :: QUADPI     = 3.14159265358979323846
	REAL, PARAMETER        :: DGR_TO_RAD = (QUADPI/180)


        MYPID = -1
 
C       SET TYPE OF OUTPUT DOC FILES WANTED
        NWANTOUT = 7
        IF (CTYPE(1:2) == 'SH') NWANTOUT = 15

C       INITIALIZE CCROT STATISTICS COUNTERS
        ANGDIFTHR   = 0.0
        CALL  AP_STAT_ADD(-1,CCROT,ANGDIF,ANGDIFTHR,CCROTLAS,
     &                   CCROTAVG,IBIGANGDIF,ANGDIFAVG,IMPROVCCROT,
     &                   CCROTIMPROV,IWORSECCROT,CCROTWORSE)

        LIMITRANGE  = (RANGE > 0.0)
        RANGECOS    = COS(RANGE * DGR_TO_RAD)  
        !write(6,*) ' range,rangecos:',range,rangecos

C       FIND NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)
        CALL FLUSHRESULTS()

        !write(6,*) ' aprings_new called '
C       READ REFERENCE IMAGES INTO REFERENCE RINGS (CIRCREF) ARRAY 
        CALL APRINGS_NEW(IREFLIST,NUMREF, NX,NY,
     &               NRING,LCIRC,NUMR,MODE,FFTW_PLANS,
     &               REFPAT,LUNPIC,CIRCREF,CIRCREF_IN_CORE,
     &               LUNRING,SCRFILE,IRTFLG)

	ALLOCATE(EXPBUF(NX,NY,NUMTH), 
     &           DLIST(5,NUMTH), 
     &           NPROJA(NUMTH), 
     &           STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           MWANT = NX*NY*NUMTH + 6*NUMTH  
           CALL  ERRT(46,'APSH_PS; EXPBUF & DLIST',MWANT)
           RETURN
        ENDIF 
        DLIST = 0.0    ! ZEROS WHOLE DLIST

        GOTREFANG = .FALSE.
        IF (LIMITRANGE .OR. CTYPE(1:2) == 'SH') THEN
C          REFANGLES FILE FOR RESTRICTED ANGULAR SEARCH  OR 'SH'
	   ALLOCATE(REFDIR(3,NUMREF),
     &              ANGREF(3,NUMREF), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
               MWANT = 6*NUMREF  
               CALL ERRT(46,'APSH_PS; REFDIR, ANGREF',MWANT)
               GOTO 9999
           ENDIF 
C          READ REF. ANGLES INTO ANGREF
C          CONVERT REF. ANGLES TO UNITARY DIRECTIONAL VECTORS (REFDIR).
	   CALL AP_GETANGAS(IREFLIST,NUMREF,0,REFANGDOC,REFPAT,
     &                      LUNPIC,LUNANG,3,ANGREF,GOTREFANG,NGOTREF,
     &                      .TRUE.,REFDIR,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
        ENDIF

        NGOTPAR = 0
        IF ((CTYPE(1:2) == 'SH' .AND. EXPANGDOC .NE. CHAR(0))) THEN
	   ALLOCATE(ANGEXP(8,NUMEXP), 
     &              EXPDIR(3,NUMEXP), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
              MWANT = 11*NUMEXP 
              CALL ERRT(46,'APSH_PS; ANGEXP....',MWANT)
              GOTO 9999
           ENDIF 

C          READ EXP. ANGLES INTO ANGEXP
C          CONVERT EXP. ANGLES TO UNITARY DIRECTIONAL VECTORS(EXPDIR).
	   CALL AP_GETANGAS(IEXPLIST,NUMEXP,0,EXPANGDOC,EXPPAT,
     &                      LUNPIC,LUNANG,8,ANGEXP,GOTEXPANG,NGOTPAR,
     &                      .TRUE.,EXPDIR,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

       ELSE
C          DUMMY ALLOCATE TO AVOID BUS ERROR ON SOME SYSTEMS
	   ALLOCATE(ANGEXP(8,1), 
     &              EXPDIR(3,1), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
              MWANT = 11  
              CALL ERRT(46,'APSH_PS; ANGEXP....',MWANT)
              GOTO 9999
           ENDIF 
        ENDIF

        anginheader = .FALSE. ! unfinished !!!!!!!!!!!
        IF (ROTFIRST) THEN
	   ALLOCATE(TMPBUF(NX,NY), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'APSH_PS; TMPBUF',NX*NY)
              GOTO 9999
           ENDIF
           anginheader = (expangdoc == '*') ! unfinished !!!!!!!!!!!
           IF (FBS_WANTED) THEN
            WRITE(NOUT,*)' ALIGNING INPUT IMAGES WITH FBS INTERPOLATION'
           ELSE
           WRITE(NOUT,*)' ALIGNING INPUT IMAGES WITH QUAD INTERPOLATION'
           ENDIF
        ENDIF 

C       LOOP OVER ALL SETS OF EXPERIMENTAL (SAMPLE) IMAGES -----------
 	DO IEXPT=1,NUMEXP,NUMTH

C          LOAD EXP. IMAGE DATA FOR THIS SET OF IMAGES
           IEND = MIN(NUMEXP,IEXPT+NUMTH-1)
           IF (ROTFIRST) THEN
C             WANT TO ROTATE/SHIFT EXP IMAGES WHEN READING THEM
	      CALL AP_GETDATA_RTSQ(IEXPLIST,NUMEXP, 
     &                    NX,NY, NX,NY,0.0,
     &                    NUMTH,EXPPAT,LUNPIC, IEXPT,IEND,
     &                    ANGINHEADER, ANGEXP, 
     &                    MPIBCAST, TMPBUF,EXPBUF,
     &                    .FALSE., ADUM, ADUM, FBS_WANTED,IRTFLG)
           ELSE
	      CALL AP_GETDATA(IEXPLIST,NUMEXP,
     &                    NX,NY, NX,NY,0.0,
     &                    NUMTH,EXPPAT,LUNPIC, IEXPT,IEND,
     &                    .TRUE., EXPBUF,
     &                    .FALSE.,ADUM,ADUM, IRTFLG)
           ENDIF
           IF (IRTFLG .NE. 0) GOTO 9999

C          NUMTH INPUT IMAGES READY TO BE ALIGNED
c$omp      parallel do private(iexp,it)
	   DO IEXP=IEXPT,MIN(NUMEXP,IEXPT+NUMTH-1)
              IT = IEXP - IEXPT + 1
	      CALL APRQ2D(EXPBUF(1,1,IT),CIRCREF,NUMR,
     &	            NX,NY,ISHRANGEX,ISHRANGEY,ISTEP,
     &	            LCIRC,NRING,NUMREF,MODE,
     &              REFDIR,EXPDIR(1,IEXP),RANGECOS,
     &              DLIST(1,IT),DLIST(2,IT),
     &              DLIST(3,IT),DLIST(4,IT),
     &              DLIST(5,IT),NPROJA(IT),
     &              CKMIRROR,LIMITRANGE,FFTW_PLANS)
	   ENDDO
c$omp      end parallel do

C          OUTPUT FROM DLIST POSITION 
C          1 - NUMBER OF MOST SIMILAR REFERENCE PROJECTION.
C          2 - NOT-NORMALIZED CORRELATION COEFFICIENT.
C          3 - PSI ANGLE. (IN=PLANE ROTATION) RANGNEW
C          4 - SX SHIFT
C          5 - SY SHIFT
    
           DO IEXP=IEXPT,MIN(NUMEXP,IEXPT+NUMTH-1)
              IT     = IEXP - IEXPT + 1
              IMGEXP = IEXPLIST(IEXP)

C             DLIST(1,IT) IS LIST NUMBER OF MOST SIMILAR REF. IMAGE 
C                 (<0 IF MIRRORED, 0 IF NO SIMILAR IMAGE )

              IREF = INT(DLIST(1,IT))
              IF (IREF .LT. 0) THEN
C                MIRRORED REFERENCE IMAGE
                 IMGREF = IREFLIST(-IREF)

C                IREFT IS FOR REFDIR INDEX
                 IREFT     = -IREF
                 MIRRORNEW = .TRUE.

              ELSEIF (IREF == 0) THEN
C                NO NEARBY REFERENCE IMAGE
                 IMGREF = 0

C                IREFT IS FOR REFDIR INDEX
                 IREFT  = 1
                 MIRRORNEW = .FALSE.

              ELSE
                 IMGREF = IREFLIST(IREF)
C                IREFT IS FOR REFDIR INDEX
                 IREFT     = IREF
                 MIRRORNEW = .FALSE.
              ENDIF
 
              CCROT     = DLIST(2,IT)
              RANGNEW   = DLIST(3,IT)
              XSHNEW    = DLIST(4,IT)
              YSHNEW    = DLIST(5,IT)
              NPROJ     = NPROJA(IT)
              PEAKV     = 1.0

              CALL AP_END(IEXP,IMGEXP,IMGREF,
     &                ANGREF(1,IREFT),REFDIR(1,IREFT),
     &                ANGEXP(1,IEXP), EXPDIR(1,IEXP),ISHRANGEX,
     &                GOTREFANG, NGOTPAR, CCROT,PEAKV,
     &                RANGNEW,XSHNEW,YSHNEW, MIRRORNEW,REFPAT,
     &                NPROJ, CTYPE, LUNDOC,PARLIST)

              CALL AP_STAT_ADD(NGOTPAR,CCROT,PARLIST(10),
     &                       ANGDIFTHR,ANGEXP(8,IEXP),
     &                       CCROTAVG,IBIGANGDIF,ANGDIFAVG,IMPROVCCROT,
     &                       CCROTIMPROV,IWORSECCROT,CCROTWORSE)
	   ENDDO
	ENDDO

        IF (LUNDOC > 0) THEN
C         SAVE CCROT & ANGULAR DISPLACEMENT STATISTICS
          CALL AP_STAT(NUMEXP,ANGDIFTHR,IBIGANGDIF,
     &                 ANGDIFAVG, CCROTAVG,
     &                 IMPROVCCROT,CCROTIMPROV,
     &                 IWORSECCROT,CCROTWORSE,
     &                 NBORDER,NSUBPIX,LUNDOC)
        ENDIF

9999    CONTINUE

C       DEALLOCATE  ARRAYS
        IF (ALLOCATED(DLIST))      DEALLOCATE(DLIST)
        IF (ALLOCATED(NPROJA))     DEALLOCATE(NPROJA)
	IF (ALLOCATED(EXPBUF))     DEALLOCATE(EXPBUF)
	IF (ALLOCATED(REFDIR))     DEALLOCATE(REFDIR)
	IF (ALLOCATED(ANGEXP))     DEALLOCATE(ANGEXP)
	IF (ALLOCATED(ANGREF))     DEALLOCATE(ANGREF)

	END


#else

C------------------------  MPI SPECIFIC SUBROUTINE --------------------

       SUBROUTINE APSH_PS(IREFLIST,NUMREF,IEXPLIST,NUMEXP, 
     &               NX,NY,ISHRANGEX,ISHRANGEY,ISTEP,
     &               NRING,LCIRC,NUMR,CIRCREF,CIRCREF_IN_CORE,
     &               MODE, REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS,
     &               REFPAT,EXPPAT,RANGE,ROTFIRST,
     &               CKMIRROR,CTYPE,LUNDOC,FBS_WANTED)

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

	INTEGER               :: IREFLIST(NUMREF) 
	INTEGER               :: NUMREF 
	INTEGER               :: IEXPLIST(NUMEXP) 
	INTEGER               :: NUMEXP 
        INTEGER               :: NX,NY,ISHRANGEX,ISHRANGEY,ISTEP
        INTEGER               :: NRING,LCIRC
        INTEGER               :: NUMR(3,NRING)
	REAL                  :: CIRCREF(LCIRC,NUMREF)
	LOGICAL               :: CIRCREF_IN_CORE
        CHARACTER (LEN=1)     :: MODE
        CHARACTER (LEN=*)     :: REFANGDOC,EXPANGDOC
        CHARACTER (LEN=*)     :: SCRFILE
        INTEGER *8            :: FFTW_PLANS(*)
        CHARACTER (LEN=*)     :: REFPAT,EXPPAT 
        REAL                  :: RANGE
	LOGICAL               :: ROTFIRST,CKMIRROR
        CHARACTER (LEN=*)     :: CTYPE
        INTEGER               :: LUNDOC
        LOGICAL               :: FBS_WANTED

        LOGICAL               :: MIRRORNEW
        LOGICAL               :: GOTREFANG,LIMITRANGE

        INTEGER, PARAMETER    :: NLISTMAX = 15
        REAL                  :: PARLIST(NLISTMAX)
        LOGICAL               :: ANGINHEADER  

        INCLUDE 'mpif.h'

        INTEGER               :: MYPID, COMM, MPIERR, NPROCS
        INTEGER               :: NEXPLOC, NREM
        INTEGER               :: IPROC, ISAM, JROW, JLOC, ILOC
        INTEGER               :: TAG,TI, IMIT,IGLB, IBEG,IEND
        INTEGER               :: ISTAT(MPI_STATUS_SIZE)

C       AUTOMATIC ARRAYS
        REAL                  :: ANGOUT(3)

C       ALLOCATED ARRAYS
        REAL,    ALLOCATABLE  :: DLIST (:,:) 
        REAL,    ALLOCATABLE  :: REFDIR(:,:),EXPDIR(:,:) 
        REAL,    ALLOCATABLE  :: ANGREF(:,:),ANGEXP(:,:)
	REAL,    ALLOCATABLE  :: TMPBUF(:,:)
        REAL,    ALLOCATABLE  :: ALOC(:,:,:)
        REAL,    ALLOCATABLE  :: EXPBUF(:,:,:)
        REAL,    ALLOCATABLE  :: DLISTLOC(:,:)
        INTEGER, ALLOCATABLE  :: NPROJA(:)
        REAL,    ALLOCATABLE  :: PARTAB(:,:),PARTABLOC(:,:)
        INTEGER, ALLOCATABLE  :: NBASE(:), PSIZE(:)

#ifdef MPI_DEBUG
        DOUBLE PRECISION      :: TCOM1, TCOM2
#endif

        PARAMETER (QUADPI = 3.14159265358979323846)
        PARAMETER (DGR_TO_RAD =   (QUADPI/180))

        INTEGER, PARAMETER     :: LUNPIC  = 77
        INTEGER, PARAMETER     :: LUNANG  = 78
        INTEGER, PARAMETER     :: LUNRING = 50

        LOGICAL, PARAMETER     :: MPIBCAST = .FALSE.

        COMM = MPI_COMM_WORLD
        CALL MPI_COMM_RANK(COMM, MYPID , MPIERR)
        CALL MPI_COMM_SIZE(COMM, NPROCS, MPIERR)
 
C       SET TYPE OF OUTPUT DOC FILES WANTED
        NWANTOUT = 7
        IF (CTYPE(1:2) == 'SH') NWANTOUT = 15

C       INITIALIZE CCROT STATISTICS COUNTERS
        ANGDIFTHR   = 0.0

C       INITIALIZE CCROT STATISTICS
        CALL AP_STAT_ADD(-1,CCROT,ANGDIF,ANGDIFTHR,CCROTLAS,
     &                   CCROTAVG,IBIGANGDIF,ANGDIFAVG,IMPROVCCROT,
     &                   CCROTIMPROV,IWORSECCROT,CCROTWORSE)

        LIMITRANGE  = (RANGE > 0.0) 
        MAXRIN      = NUMR(3,NRING)
        RANGECOS    = COS(RANGE*DGR_TO_RAD)

C       FIND NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)
        CALL FLUSHRESULTS()

C       READ REFERENCE IMAGES INTO REFERENCE RINGS (CIRCREF) ARRAY 
        CALL APRINGS_NEW(IREFLIST,NUMREF,  NX,NY, 
     &               NRING,LCIRC,NUMR,MODE,FFTW_PLANS,
     &               REFPAT,LUNPIC,CIRCREF,CIRCREF_IN_CORE,
     &               LUNRING,SCRFILE,IRTFLG)
 
        GOTREFANG = .FALSE.
        IF (LIMITRANGE .OR.  CTYPE(1:2) == 'SH') THEN
C          REFANGLES FILE FOR RESTRICTED ANGULAR SEARCH  OR 'SH'
	   ALLOCATE(REFDIR(3,NUMREF),
     &              ANGREF(3,NUMREF), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
               MWANT = 6*NUMREF  
               CALL ERRT(46,'APSH_PS; REFDIR, ANGREF',MWANT)
               RETURN
           ENDIF 

C          READ REF. ANGLES INTO ANGREF FROM HEADER OR DOC FILE
	   CALL AP_GETANGAS(IREFLIST,NUMREF,0,REFANGDOC,REFPAT,
     &                      LUNPIC,LUNANG,3,ANGREF,GOTREFANG,NGOTREF,
     &                      .TRUE.,REFDIR,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
        ENDIF

        NGOTPAR = 0
        IF (CTYPE(1:2) == 'RQ' .OR.
     &     (CTYPE(1:2) == 'SH' .AND. 
     &     EXPANGDOC .NE. CHAR(0))) THEN

           ALLOCATE(ANGEXP(8,NUMEXP), 
     &              EXPDIR(3,NUMEXP), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              MWANT = 11*NUMEXP 
              CALL ERRT(46,'APSH_PS; ANGEXP & EXPDIR',MWANT)
              RETURN
           ENDIF 

C          READ EXP. ANGLES INTO ANGEXP
C          CONVERT EXP. ANGLES TO UNITARY DIRECTIONAL VECTORS(EXPDIR).
           CALL AP_GETANGAS(IEXPLIST,NUMEXP,0,EXPANGDOC,EXPPAT,
     &                      LUNPIC,LUNANG,8,ANGEXP,NGOTPAR,GOTEXPANG,
     &                      .TRUE.,EXPDIR,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

       ELSE
C          DUMMY ALLOCATE TO AVOID BUS ERROR ON SOME SYSTEMS
           ALLOCATE(ANGEXP(8,1), 
     &              EXPDIR(3,1), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              MWANT = 11*1 
              CALL ERRT(46,'APSH_PS; EXPDIR',MWANT)
              RETURN
           ENDIF 
       ENDIF

C       === PARTITION AND DISTRIBUTE EXPERIMENTAL IMAGES ===

        ALLOCATE(PSIZE(NPROCS),
     &           NBASE(NPROCS), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'APSH_PS; PSIZE & NBASE.',2*NPROCS)
           RETURN
        ENDIF

        CALL SETPART(NUMEXP, PSIZE, NBASE)
        NEXPLOC = PSIZE(MYPID+1)

        ALLOCATE(ALOC  (NX,NY,NEXPLOC),
     &           EXPBUF(NX,NY,PSIZE(1)), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MWANT = NX*NY*NEXPLOC + NX*NY*PSIZE(1)
           CALL ERRT(46,'APSH_PS; ALOC...',MWANT)
           RETURN
        ENDIF
        ALOC = 0.0

#ifdef MPI_DEBUG
        WRITE(6,111) NBASE(MYPID+1), MYPID
        CALL FLUSHFILE(6)
 111    FORMAT(' APSH_PS: NBASE = ', I5, ' MYPID = ', I5)
        TCOM0 = MPI_WTIME()
#endif

        ANGINHEADER = .FALSE. ! unfinished !!!!!!!!!!!
        IF (ROTFIRST) THEN
	   ALLOCATE(TMPBUF(NX,NY), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'APSH_PS; TMPBUF',NX*NY)
              GOTO 9999
           ENDIF
           anginheader = (expangdoc == '*') ! unfinished !!!!!!!!!!!
           IF (FBS_WANTED .AND. MYPID == 0) THEN
            WRITE(NOUT,*)' ALIGNING INPUT IMAGES WITH FBS INTERPOLATION'
           ELSE
           WRITE(NOUT,*)' ALIGNING INPUT IMAGES WITH QUAD INTERPOLATION'
           ENDIF
        ENDIF 


C       PROCESS 0 READS IMAGES & DISTRIBUTES TO DIFFERENT PROCESSORS. ===
        DO IPROC = 1, NPROCS
           NLOC = PSIZE(IPROC)

C          === CALCULATE THE GLOBAL INDEX ===

           IBEG = NBASE(IPROC) + 1
           IEND = NBASE(IPROC) + NLOC 

C          ALTHOUGH THE FOLLOWING IS CALLED BY ALL PROCESSORS,
C          ONLY ONE PROCESSOR (MYPID=0) READS IMAGES INTO EXPBUF
C          IT DOES NOT BROADCAST EXPBUF!

           IF (ROTFIRST) THEN
C             WANT TO ROTATE/SHIFT EXP IMAGES WHEN READING THEM
	      CALL AP_GETDATA_RTSQ(IEXPLIST,NUMEXP, 
     &                    NX,NY, NX,NY,0.0,
     &                    NUMTH,EXPPAT,LUNPIC, IBEG,IEND,
     &                    ANGINHEADER, ANGEXP, 
     &                    MPIBCAST, TMPBUF,EXPBUF,
     &                    .FALSE., ADUM, ADUM, FBS_WANTED,IRTFLG)
           ELSE
              CALL AP_GETDATA(IEXPLIST,NUMEXP,
     &                    NX,NY, NX,NY,0.0,
     &                    NUMTH,EXPPAT,LUNPIC, IBEG,IEND,
     &                    .FALSE., EXPBUF,
     &                    .FALSE., ADUM,ADUM,IRTFLG)
           ENDIF
           IF (IRTFLG .NE. 0) STOP

           IF (IPROC > 1) THEN
              IF (MYPID == 0) THEN
#ifdef MPI_DEBUG
                 WRITE(6,222) IPROC-1
                 CALL FLUSHFILE(6)
 222             FORMAT(' APSH_PS: SENDING TO PID = ', I3)
#endif
                 CALL MPI_SEND(EXPBUF   , NX*NY*NLOC, MPI_REAL,
     &                         IPROC-1, IPROC-1       , COMM    ,
     &                         MPIERR)

              ELSEIF (MYPID == IPROC-1) THEN

C                === SLAVES RECEIVE LOCAL PIECES ===

                 CALL MPI_RECV(ALOC , NX*NY*NLOC, MPI_REAL,
     &                         0    , MPI_ANY_TAG   , COMM    ,
     &                         ISTAT, MPIERR)
                 IF (MPIERR .NE. 0) THEN
                     WRITE(6,*) ' RECV FAILED'
                     STOP
                 ENDIF
#ifdef MPI_DEBUG
                 WRITE(6,223) MYPID
                 CALL FLUSHFILE(6)
 223             FORMAT(' APSH_PS: RECEIVED BY MYPID = ', I3)
#endif
              ENDIF
           ELSEIF (MYPID == 0) THEN  

C             === SIMPLY COPY FROM EXPBUF TO ALOC ===
              DO JLOC = 1, NLOC
                 DO ISAM = 1, NX
                    DO JROW = 1, NY
                       ALOC(ISAM,JROW,JLOC) = EXPBUF(ISAM,JROW,JLOC)
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
        ENDDO

        IF (ALLOCATED(EXPBUF)) DEALLOCATE(EXPBUF)
        IF (ALLOCATED(TMPBUF)) DEALLOCATE(TMPBUF)

        ALLOCATE(PARTAB(15,NUMEXP), 
     &           PARTABLOC(15,NEXPLOC),
     &           DLISTLOC(5,NEXPLOC),
     &           NPROJA(NEXPLOC),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MWANT = 15*NUMEXP + 15*NEXPLOC + 6*NEXPLOC
           CALL ERRT(46,'APSH_PS; PARTAB..',MWANT)
           RETURN
        ENDIF

#ifdef MPI_DEBUG
        TCOM1 = MPI_WTIME()
        IF (MYPID == 0) WRITE(6,440) TCOM1-TCOM0
 440    FORMAT(' APSH_PS: DATA DIST TIME = ', 1PE11.3)
        WRITE(6,444) MYPID
 444    FORMAT(' APSH_PS: CALLING APRQ2D.., MYPID = ', I3)
#endif

        PARTAB     = 0.0
        PARTABLOC = 0.0

        DO  IT=1,NEXPLOC
           IGLB = NBASE(MYPID+1) + IT
           CALL APRQ2D(ALOC(1,1,IT),CIRCREF,NUMR,
     &              NX,NY,ISHRANGEX,ISHRANGEY,ISTEP,
     &              LCIRC,NRING,NUMREF,MODE,
     &              REFDIR,EXPDIR(1,IGLB),RANGECOS,
     &              DLISTLOC(1,IT),DLISTLOC(2,IT),
     &              DLISTLOC(3,IT),DLISTLOC(4,IT),
     &              DLISTLOC(5,IT),NPROJA(IT),
     &              CKMIRROR,LIMITRANGE, FFTW_PLANS)

C          OUTPUT IN DLISTLOC
C          1 - NUMBER OF THE MOST SIMILAR REFERENCE PROJECTION.
C          2 - NOT-NORMALIZED CORRELATION COEFFICIENT.
C          3 - PSI ANGLE. (IN=PLANE ROTATION)
C          4 - SX SHIFT
C          5 - SY SHIFT
 
           IMGEXP = IEXPLIST(IGLB)
 
C          DLISTLOC(1,IT) IS LIST NUMBER OF MOST SIMILAR REF. IMAGE
C              (<0 IF MIRRORED, 0 IF NONE )
 
           IREF = INT(DLISTLOC(1,IT))
           IF (IREF < 0) THEN
C             MIRRORED REFERENCE IMAGE
              IMGREF = IREFLIST(-IREF)
 
C             IREFT IS FOR REFDIR INDEX
              IREFT     = -IREF
              MIRRORNEW = .TRUE.
 
           ELSEIF (IREF == 0) THEN
C             NO NEARBY REFERENCE IMAGE
              IMGREF = 0
 
C             IREFT IS FOR REFDIR INDEX
              IREFT  = 1
              MIRRORNEW = .FALSE.
 
           ELSE
              IMGREF = IREFLIST(IREF)
C             IREFT IS FOR REFDIR INDEX
              IREFT  = IREF
              MIRRORNEW = .FALSE.
           ENDIF
 
           CCROT     = DLISTLOC(2,IT)
           RANGNEW   = DLISTLOC(3,IT)
           XSHNEW    = DLISTLOC(4,IT)
           YSHNEW    = DLISTLOC(5,IT)
           NPROJ     = NPROJA(IT)
           PEAKV     = 1.0
 
           CALL AP_END(IGLB,IMGEXP,IMGREF,
     &             ANGREF(1,IREFT),REFDIR(1,IREFT),
     &             ANGEXP(1,IGLB), EXPDIR(1,IGLB),ISHRANGEX,
     &             GOTREFANG, NGOTPAR, CCROT,PEAKV,
     &             RANGNEW,XSHNEW,YSHNEW,MIRRORNEW, REFPAT,
     &             NPROJ, CTYPE,  LUNDOC, PARTABLOC(1,IT))
        ENDDO
C
        DO IPROC = 1, NPROCS
           PSIZE(IPROC) = 15 * PSIZE(IPROC)
           NBASE(IPROC) = 15 * NBASE(IPROC)
        ENDDO
        CALL MPI_ALLGATHERV(PARTABLOC, PSIZE(MYPID+1),
     &                      MPI_REAL, PARTAB, PSIZE, NBASE,
     &                      MPI_REAL, COMM  , MPIERR)

C       WRITE HAS BEEN SYNCHRONIZED WITHIN LUNDOCWRTDAT

        IF (NWANTOUT > 7) THEN
           IF (LUNDOC > 0) THEN
C             SAVE IN ALIGNMENT DOC FILE

C             <,<,<, MIR-REF#,IMG#,INPLANE<, SX,SY,NPROJ,DIF,CCROT,INPLANE<,SX,SY

              DO IT = 1, NUMEXP
                 CALL LUNDOCWRTDAT(LUNDOC,IT,PARTAB(1,IT),
     &                             NWANTOUT,IRTFLG)
              ENDDO
           ENDIF
        ELSE
           IF (LUNDOC >. 0) THEN
C             SAVE IN ALIGNMENT DOC FILE
C             MIR-REF#, CCROT, INPLANE<, SX,SY, IMG#, < DIFF
              DO IT = 1, NUMEXP
                 CALL LUNDOCWRTDAT(LUNDOC,IT,PARTAB(1,IT),
     &                             NWANTOUT,IRTFLG)
              ENDDO
           ENDIF
        ENDIF
9999    CONTINUE


C       DEALLOCATE  ARRAYS
        IF (ALLOCATED(ALOC))       DEALLOCATE(ALOC)
        IF (ALLOCATED(PSIZE))      DEALLOCATE(PSIZE)
        IF (ALLOCATED(NBASE))      DEALLOCATE(NBASE)
        IF (ALLOCATED(NPROJA))     DEALLOCATE(NPROJA)
        IF (ALLOCATED(DLISTLOC))   DEALLOCATE(DLISTLOC)
        IF (ALLOCATED(PARTAB))     DEALLOCATE(PARTAB)
        IF (ALLOCATED(PARTABLOC))  DEALLOCATE(PARTABLOC)
        IF (ALLOCATED(REFDIR))     DEALLOCATE(REFDIR)
        IF (ALLOCATED(ANGEXP))     DEALLOCATE(ANGEXP)
        IF (ALLOCATED(ANGREF))     DEALLOCATE(ANGREF)

        END
C ------------------------- END OF MPI CODE -----------------------
#endif






C+**********************************************************************
C
C APRQ2D.F
C 
C  PARAMETERS:
C                DIREF    NUMBER OF  MOST SIMILAR REF. PROJ.  (OUTPUT)
C                            (NEGATIVE IF MIRRORED)
C                CCROT    CORR COEFF.                         (OUTPUT)
C                RANGNEW  INPLANE ANGLE                       (OUTPUT)
C                XSHSUM   SHIFT                               (OUTPUT)
C                YSHSUM   SHIFT                               (OUTPUT)
C                NPROJ                                        (OUTPUT)
C
C-**********************************************************************

	SUBROUTINE APRQ2D(EXPBUF,CIRCREF,NUMR,
     &	             NX,NY,ISHRANGEX,ISHRANGEY,ISTEP,
     &	             LCIRC,NRING,NUMREF,MODE,
     &               REFDIR,EXPDIR,RANGECOS,
     &               DIREF,CCROT,RANGNEW,XSHSUM,YSHSUM,NPROJ,
     &               CKMIRROR,LIMITRANGE,FFTW_PLANS)

C       NOTE: RUNS WITHIN OMP PARALLEL SECTION OF CODE IF NOT UNDER MPI!

        INCLUDE 'MAKE_CLOSE_LIST.INC'  

	REAL                      :: EXPBUF(NX,NY)
        REAL                      :: CIRCREF(LCIRC,NUMREF)
        INTEGER                   :: NUMR(3,NRING) 

	DOUBLE PRECISION          :: FITP(-1:1,-1:1)
        CHARACTER (LEN=1)         :: MODE
	REAL                      :: REFDIR(3,NUMREF), EXPDIR(3)
        INTEGER *8                :: FFTW_PLANS(*)

C       AUTOMATIC ARRAYS
	DOUBLE PRECISION          :: FIT(-ISTEP:ISTEP,-ISTEP:ISTEP)
	REAL                      :: ROTMP(-ISTEP:ISTEP,-ISTEP:ISTEP)

C       ALLOCATABLE ARRAYS
        INTEGER, POINTER          :: LCG(:)
        REAL, ALLOCATABLE         :: CIRCEXP(:)

	DOUBLE PRECISION          :: CCROTD,PEAK,CCROTD_INTERP
	DOUBLE PRECISION          :: CCOA
        LOGICAL                   :: CKMIRROR,LIMITRANGE
        LOGICAL                   :: MIRRORED
        LOGICAL                   :: ISMIRRORED,USE_UN,USE_MIR
        REAL                      :: ADUM

        LOGICAL, PARAMETER        :: USE_OMP = .FALSE.

        REAL, PARAMETER           :: QUADPI = 3.1415926535
        REAL, PARAMETER           :: DGR_TO_RAD = (QUADPI/180)


        PEAK   = 0.0
        WR     = 0.0    ! DUMMY VALUE FLAG FOR APRINGS CALL

        MAXRIN = NUMR(3,NRING)   ! ACTUAL LENGTH OF LONGEST RING
#ifdef SP_LIBFFTW3
        MAXRIN = NUMR(3,NRING) - 2
#endif

	ALLOCATE(CIRCEXP(LCIRC),STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           CALL  ERRT(46,'APRQ2D; CIRCEXP',LCIRC)
           RETURN
        ENDIF

        IEND  = NUMREF

C       write(6,*) 'numref,range: ',numref,rangecos,'expdir: ',expdir

C       IF LIMITRANGE, MAKE LIST OF CLOSE REF. IMAGES, RETURNS: NPROJ
        CALL MAKE_CLOSE_LIST(NUMREF,LIMITRANGE,
     &                       REFDIR,EXPDIR,
     &                       RANGECOS, CKMIRROR, 
     &                       LCG, NPROJ, IRTFLG)

	IF (NPROJ .LE. 0) THEN
C          THERE IS NO REFERENCE WITHIN SEARCH RANGE
           XSHSUM  = 0
	   YSHSUM  = 0
           DIREF   = 0
           RANGNEW = 0
           CCROT   = -1.0 
           GOTO 9999	
        ENDIF
        IEND = NPROJ

	CCROTD = -HUGE(CCROTD)

c       GO THROUGH CENTERS FOR SHIFT ALIGNMENT ------------------------
	DO JT=-ISHRANGEY,ISHRANGEY,ISTEP
	   CNR2 = NY / 2 + 1 + JT

	   DO IT=-ISHRANGEX,ISHRANGEX,ISTEP
	      CNS2 = NX / 2 + 1 + IT

C             NORMALIZE IMAGE VALUES UNDER THE MASK OVER VARIANCE RANGE
C             INTERPOLATE TO POLAR COORD., CREATE FOURIER OF: A_CIRC
C             NO WEIGHTING OF RINGS
              !write(6,*) ' aprings one:,',jt,it,cns2,cnr2

	      CALL APRINGS_ONE_NEW(NX,NY, CNS2,CNR2, EXPBUF,.FALSE.,
     &                             MODE,NUMR,NRING,LCIRC,WR,FFTW_PLANS,
     &                             CIRCEXP,IRTFLG)

C             COMPARE EXP. IMAGE WITH ALL REFERENCE IMAGES ------------
	      DO IRR=1,IEND
                 IR = IRR
                 IF (LIMITRANGE) IR = ABS(LCG(IRR))
                
                 IF (CKMIRROR .AND. LIMITRANGE) THEN
C                   ONLY SEARCH EITHER MIRRORED OR NON-MIRRORED
                    USE_UN  = (LCG(IRR) .GE. 0)
                    USE_MIR = (LCG(IRR) .LT. 0)
                 ELSE
C                   SEARCH BOTH MIRRORED & NON-MIRRORED IF CHKMIR
                    USE_UN  = .TRUE.
                    USE_MIR = CKMIRROR
                 ENDIF

	         CALL CROSRNG_2(CIRCREF(1,IR),CIRCEXP,
     &                          LCIRC,NRING, MAXRIN,NUMR, 
     &                          USE_OMP,FFTW_PLANS(1),
     &                          USE_UN,USE_MIR,
     &                          ISMIRRORED,CCOA,RANGO)

                 IF (CCOA .GE. CCROTD)  THEN
C                   BETTER MATCH WITH THIS POSITION 
	            CCROTD  = CCOA
	            IBE     = IR
	            ISX     = IT
	            ISY     = JT
	            RANGNEW = ANG_N(RANGO,MODE,MAXRIN)
	            IDIS    = IR
                    IF (ISMIRRORED) IDIS = -IR
	         ENDIF

	      ENDDO ! END OF:  DO IRR=1,IEND ---------------------- REFS.
	   ENDDO    ! END OF:  DO IT=-ISHRANGEX,ISHRANGEX,ISTEP
	ENDDO       ! END OF:  DO JT=-ISHRANGEY,ISHRANGEY,ISTEP

        SX       = ISX              ! BEST SHIFTS
        SY       = ISY
        CCROT    = CCROTD           ! BEST ROTATION

#ifdef NEVER
        write(6,921) idis,isx,isy,ccrotd,rangnew
921     format(' 1 ',i5,' (',i3,',',i3,'): ',f12.4,' ',2f8.2,f6.1)
#endif
       
C       CHECK LOCATIONS WITHIN ISHRANGE AROUND MAX  ------------------
        DIREF  = IDIS
C       WHEN (IDIS .LT. 0)  CHECK MIRRORED ONLY
        MIRRORED = (IDIS .LT. 0)

C       DO NOT INTERPOLATE FOR POINT ON THE EDGE
        IF (IABS(ISX).NE.ISHRANGEX .AND. IABS(ISY).NE.ISHRANGEY) THEN
C          NOT ON BOUNDARY, HAVE TO FIND NEIGHBOURING VALUES

	   FIT(0,0)   = CCROTD
           ROTMP(0,0) = RANGNEW

           DO JT=-ISTEP,ISTEP
	      CNR2 = NY / 2 + 1 + JT + ISY
	      DO IT=-ISTEP,ISTEP
	         CNS2 = NX / 2 + 1 + IT + ISX

   	         IF (IT.NE.0 .OR. JT.NE.0) THEN
C                   NORMALIZE IMAGE VALUES UNDER THE MASK OVER VARIANCE RANGE
C                   INTERPOLATE INTO POLAR COORDINATES
C                   CREATES FOURIER OF: CIRCEXP

                    !!write(6,*) ' aprings max:,',jt,it,cns2,cnr2
 
	            CALL APRINGS_ONE_NEW(NX,NY, CNS2,CNR2, 
     &                           EXPBUF,.FALSE.,
     &                           MODE,NUMR,NRING,LCIRC,WR,FFTW_PLANS,
     &                           CIRCEXP,IRTFLG)

	            CALL CROSRNG_2(CIRCREF(1,IBE),CIRCEXP,
     &                          LCIRC,NRING, MAXRIN,NUMR, 
     &                          USE_OMP,FFTW_PLANS(1),
     &                          .NOT. MIRRORED,MIRRORED,
     &                          ISMIRRORED,FIT(IT,JT),ROTMP(IT,JT))

C                   RECORD BEST ANGLE IN ROTMP
                    ROTMP(IT,JT) = ANG_N(ROTMP(IT,JT),MODE,MAXRIN)

                 ENDIF  ! END OF:  IF (IT.NE.0 .OR. JT.NE.0)
	      ENDDO     ! END OF:  DO IT=-ISTEP,ISTEP
           ENDDO        ! END OF:  DO JT=-ISTEP,ISTEP

#ifdef NEVER
          fit(it,jt) = fit(it,jt) 
          write(6,961) it,jt, fit(it,jt),cns2,cnr2, rotmp(it,jt)
961       format(' fit(',i2,',',i2,') : ',f12.4,' ',2f8.2,f6.1)
#endif


C          FIND THE MAXIMUM CC ANGLE WITHIN +/-ISTEP
C          MAXIMUM CANNOT BE ON THE EDGE, I.E., IT,JT/=ISTEP

	   AFIT     = FIT(0,0)
	   JTMA     = 0
	   ITMA     = 0
           RANGNEWT = ROTMP(0,0)

	   IF (ISTEP > 1)  THEN
	      DO JT=-ISTEP+1,ISTEP-1
	         DO IT=-ISTEP+1,ISTEP-1
	            IF (FIT(IT,JT) > AFIT)  THEN
	               AFIT     = FIT(IT,JT)
	               RANGNEWT = ROTMP(IT,JT) !compiler bug on OPT64
	               ITMA     = IT
	               JTMA     = JT
	            ENDIF
	         ENDDO  ! END OF:  DO IT=-ISTEP,ISTEP
              ENDDO     ! END OF:  DO JT=-ISTEP,ISTEP
           ENDIF

C          TEMP VARIABLE OVERCOMES COMPILER BUG ON OPT 64 PGI 6.0
           RANGNEW = RANGNEWT
           CCROTD  = AFIT
        
C          COPY VALUES AROUND THE PEAK.
           DO JT=-1,1
              DO IT=-1,1
                 FITP(IT,JT) = FIT(ITMA+IT,JTMA+JT)
c                write(6,910) it,jt,fitp(it,jt)
910              format(' fitp(',i5,',',i5,') : ',f14.4)
              ENDDO
           ENDDO

C          UPDATE INTEGRAL LOCATION OF PEAK
	   ISX    = ISX + ITMA
	   ISY    = ISY + JTMA
           SX     = ISX
           SY     = ISY

#ifdef NEVER
	   cnr2 = NY / 2 + 1 + isy
           cns2 = NX / 2 + 1 + isx
           write(6,905) idis,isx,isy, ccrotd,rangnew, cns2,cnr2,sx,sy
905        format(' 2 ',i5,' (',i3,',',i3,'): ',f12.4,' ',5f8.2)
#endif

C          SUB-PIXEL INTERPOLATION -----------------------------------
           CALL PARABLD(FITP,SSX,SSY,PEAK)

	   IF (ABS(SSX) .LT. 1.0 .AND. ABS(SSY) .LT. 1.0)  THEN
C             NOT ON EDGE OF 3x3 INTERPOLATION AREA

	      CNS2 = NX/2+1 + SX + SSX
	      CNR2 = NY/2+1 + SY + SSY

C             NORMALIZE IMAGE VALUES UNDER THE MASK OVER VARIANCE RANGE
C             INTERPOLATE INTO POLAR COORDINATES, 
C             CREATE FOURIER OF: CIRCEXP
              !write(6,*) ' aprings sub:,',sx,sy,cns2,cnr2
	      CALL APRINGS_ONE_NEW(NX,NY, CNS2,CNR2, EXPBUF,.FALSE.,
     &                           MODE,NUMR,NRING,LCIRC,WR,FFTW_PLANS,
     &                           CIRCEXP,IRTFLG)

	      CALL CROSRNG_2(CIRCREF(1,IBE),CIRCEXP,
     &                       LCIRC,NRING, MAXRIN,NUMR, 
     &                       USE_OMP,FFTW_PLANS(1),
     &                      .NOT. MIRRORED,MIRRORED,
     &                       ISMIRRORED,CCROTD_INTERP,RANGNEW_INTERP)

#ifdef NEVER
	      rt1 = ang_n(rangnew_interp,mode,maxrin)
              write(6,904) idis,isx,isy, ccrotd_interp,rt1,
     &                     cns2,cnr2, sx+ssx,sy+ssy
904           format(' 3 ',i5,' (',i3,',',i3,'): ',f12.4,' ',5f8.2)
#endif

              IF (CCROTD_INTERP > CCROTD) THEN
C                USE SUB-PIXEL LOCATION
                 CCROTD  = CCROTD_INTERP
	         RANGNEW = ANG_N(RANGNEW_INTERP,MODE,MAXRIN)
	         SX      = SX + SSX 
	         SY      = SY + SSY 
              ENDIF
           ENDIF
        ENDIF

        CCROT = CCROTD
	SX    = -SX
	SY    = -SY
 
C       HAVE TO CHANGE ORDER OF SHIFT & ROTATION.
C       IN THIS PROGRAM IMAGE IS SHIFTED FIRST, ROTATED SECOND.
C       IN 'RT SQ' IT IS ROTATION FIRST, SHIFT SECOND.
C       THIS CODE CORRESPONDS TO 'SA P'.
	CO     =  COS(RANGNEW * DGR_TO_RAD)
	SO     = -SIN(RANGNEW * DGR_TO_RAD)
	XSHSUM = SX*CO - SY*SO
	YSHSUM = SX*SO + SY*CO
C       ALMOST ZERO IS LIKELY TO BE ZERO
        IF (ABS(XSHSUM)  .LT. 0.08) XSHSUM  = 0.0
        IF (ABS(YSHSUM)  .LT. 0.08) YSHSUM  = 0.0
        IF (ABS(RANGNEW) .LT. 0.08) RANGNEW = 0.0

#ifdef NEVER
        cns2 = NX / 2 + 1 - sx
        cnr2 = NY / 2 + 1 - sy
        write(6,906) idis,isx,isy, ccrotd,rangnew,
     &               cns2,cnr2, xshsum,yshsum
906     format(' 4 ',i5,' (',i3,',',i3,'): ',f12.4,' ',5f8.2)

        write(6,*) ' ------------------------------------- '
        write(6,*) '  '
#endif

9999    IF (ASSOCIATED(LCG))     DEALLOCATE(LCG)
        IF (ALLOCATED(CIRCEXP))  DEALLOCATE(CIRCEXP)
        NULLIFY(LCG)

	END

