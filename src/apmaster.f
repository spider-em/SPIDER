
C++*********************************************************************
C
C APMASTER.F        CAN OUTPUT TO REGISTERS NOW   MAY 01 ARDEAN LEITH
C                   CAN GET ANGLES FROM HEADER    JUN 01 ARDEAN LEITH
C                   NORMASS -> NORMAS             OCT 01 ARDEAN LEITH
C                   SAVDN1 + SAVD BUG             JAN 02 ARDEAN LEITH
C                   PROMPTS                       JAN 02 ARDEAN LEITH
C                   UNSAV LOOP IMPROVED           SEP 02 ARDEAN LEITH
C                   ADDED ANG. DIFFERENCE         OCT 02 ARDEAN LEITH
C                   OPFILEC                       FEB 03 ARDEAN LEITH
C                   MERGED WITH DSFR & DSGRS      AUG 03 ARDEAN LEITH
C                   MERGED WITH DSFS              AUG 03 ARDEAN LEITH
C                   MERGED WITH AP_SH1            SEP 03 ARDEAN LEITH
C                   HEADER OUTPUT CHANGED         OCT 03 ARDEAN LEITH
C                   MPI OUTPUT CHANGED            FEB 04 Chao Yang
C                   'AP SH'                       FEB 04 ARDEAN LEITH
C                   'DOC FILE HEADERS'            APR 04 ARDEAN LEITH
C                   OR REF FILE                   JUN 04 ARDEAN LEITH
C                   PSI,THE,PHI                   JUN 04 ARDEAN LEITH
C                   REF_RINGS AUTO CREATION       JAN 05 ARDEAN LEITH
C                   ANG DIFF. THRESHOLD           FEB 05 ARDEAN LEITH
C                   'OR' HAD OUTPUT FILE BUG      AUG 05 ARDEAN LEITH
C                   'AP RQN' MIRRORED BUG         DEC 05 ARDEAN LEITH
C                   'AP SCC' INCORPORATED         FEB 08 ARDEAN LEITH
C                   APRINGS_INIT_PLANS            APR 08 ARDEAN LEITH
C                   OBSOLETE OPERATION MSG.       MAY 08 ARDEAN LEITH
C                   'OR NQ, OR MQ --> OR SH'      JUN 08 ARDEAN LEITH
C                   APRINGS RAYS                  JUN 08 ARDEAN LEITH
C                   FFTW3_KILLPLANS               JAN 09 ARDEAN LEITH
C                   'AP SH' CIRCREF ALLOC MSG.    AUG 09 ARDEAN LEITH
C                   MOVED 'AP SCC' OUT            AUG 09 ARDEAN LEITH
C                   ISHRANGEX                     FEB 10 ARDEAN LEITH
C                   CUDA, APSH_SS PARAMETERS      APR 10 ARDEAN LEITH
C                   REMOVED OBSOLETE OPERATIONS   JUN 10 ARDEAN LEITH
C                   'AP REF' REGISTER ONLY BUG    AUG 10 ARDEAN LEITH
C                   'AP OR' NO CALLING MSG        OCT 10 ARDEAN LEITH
C                   MRQLI & DSGR RENAME           JAN 11 ARDEAN LEITH
C                   AP_ RENAME                    FEB 11 ARDEAN LEITH
C                   WEIGHT = (YN == 'Y' ..        APR 11 ARDEAN LEITH
C                   ALLOCATABLE NPLANS            APR 11 ARDEAN LEITH
C                   APRINGS_INIT_PLANS PARAMS     JUN 11 ARDEAN LEITH
C                   AP FOU                        AUG 11 ARDEAN LEITH
C                   CKMIRROR PARSING              AUG 11 ARDEAN LEITH
C                   RAY1,RAY2                     NOV 11 ARDEAN LEITH
C                   ROTFIRST                      DEC 11 ARDEAN LEITH
C                   FBS_WANTED                    JAN 12 ARDEAN LEITH
C                   AP SHC PARAM. SUMMARY         FEB 12 ARDEAN LEITH
C                   AP FOU PATM                   JUN 12 ARDEAN LEITH
C                   DENOISE, ROTFIRST=FBS         SEP 12 ARDEAN LEITH
C                   RING LIMIT TRAP               FEB 13 ARDEAN LEITH
C
C **********************************************************************
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* Authors: J. Frank & A. Leith                                       *
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
C=* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
C=* General Public License for more details.                           *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program. If not, see <http://www.gnu.org/licenses> *
C=*                                                                    *
C **********************************************************************
C
C    APMASTER(MODE,CTYPE)                                  
C
C    MASTER IO AND INITIALIZATION ROUTINE FOR MOST 'AP ..' OPERATIONS
C  
C    'AP AL'   -- ??
C    'OR SH'
C    'AP I'    -- CREATE RINGS FILE ONLY
C    'AP MI'   -- 
C    'AP REF'  -- APRE_P or APRE_PM
C    'AP REFT' -- FORCES: APREF_PM
C    'AP REFF' -- FORCES: NON-INCORE EVEN IF SIZE IS OK
C    'AP REFB' -- FORCES: DOC FILE OUTPUT EVEN IF HAS OPERATION REGISTERS
C    'AP SH'   -- APSH_SS or APSH_PS
C    'AP SHC'  -- APSH_PSC; COEFF, NON-TRANSFORMED RINGS, CPLX VAR.
C    'AP SHF'  -- NON-INCORE EVEN IF SIZE IS OK
C    'AP SHT'  -- FORCES: APSH_SS
C    'AP SHG'  -- GPU (MAY NOT BE LINKED)
C    'AP I'    -- CREATE RINGS FILE ONLY
C    'AP FOU'  -- NEW ALGORITHM
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE APMASTER(MODE,CTYPE)

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

	CHARACTER(LEN=1)       :: MODE
	CHARACTER(LEN=*)       :: CTYPE

        INTEGER, ALLOCATABLE   :: IMGLST(:)
        INTEGER, ALLOCATABLE   :: NUMR(:,:)
        REAL,    ALLOCATABLE   :: CIRCREF(:,:)

#ifndef SP_32
        INTEGER *8             :: IASK8,IOK
        INTEGER *8,ALLOCATABLE :: FFTW_PLANS(:)
#else
        INTEGER                :: IASK8,IOK
        INTEGER,   ALLOCATABLE :: FFTW_PLANS(:)
#endif

        CHARACTER (LEN=MAXNAM) :: ASK,SCRFILE,FILNAM,REFANGDOC
        CHARACTER (LEN=MAXNAM) :: REFPAT,EXPPAT,EXPANGDOC,OUTANG

	CHARACTER(LEN=1)       :: NULL,ANS,YN,CTEMP
	CHARACTER(LEN=80)      :: PROMPT,MSG
	CHARACTER(LEN=220)     :: COMMEN
        LOGICAL                :: CIRCREF_IN_CORE,CKMIRROR
        LOGICAL                :: WINDOW,NEWFILE,WEIGHT,GPU,WANTDOC
        LOGICAL                :: ROTFIRST,GOTMIR,FBS_WANTED 
        LOGICAL                :: DENOISE,GOTRTSH 
        REAL                   :: VALUES(6)

        INTEGER, PARAMETER     :: LUNREF  = 50
        INTEGER, PARAMETER     :: LUNEXP  = 51
        INTEGER, PARAMETER     :: LUNRING = 52
	!USED IN CALLED ROUTINE
        INTEGER, PARAMETER     :: INPIC   = 77 
        INTEGER, PARAMETER     :: INANG   = 78 
        INTEGER, PARAMETER     :: NDOC    = 55 
        INTEGER, PARAMETER     :: NSCF    = 50 

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID 

        NULL   = CHAR(0)
        NILMAX = NIMAX

#ifdef USE_MPI
        IF (CTYPE(1:3) == 'SHC') THEN
           CALL ERRT(101,"SUB-OPERATION NOT ON MPI, USE 'AP SH'",NDUM)
           RETURN
        ENDIF
#endif

        IF (CTYPE(1:2) == 'OR') THEN
           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,REFPAT,LUNREF,'O',ITYPE,NX,NY,
     &               NZ,MAXIM,'REFERENCE',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           NUMREF    = 1
           INUMBR(1) = 0
        ELSE
           NUMREF = 0
C          ASK FOR TEMPLATE AND NUMBERS FOR REFERENCE IMAGES
	   CALL FILELIST(.TRUE.,LUNREF,REFPAT,NLET,INUMBR,NILMAX,NUMREF,
     &           'TEMPLATE FOR REFERENCE IMAGES',IRTFLG)
	   IF (IRTFLG .NE. 0) RETURN

           IF (MYPID <= 0) WRITE(NOUT,2001) NUMREF
2001       FORMAT('  Number of reference images: ',I7)
        ENDIF

C       NUMREF - TOTAL NUMBER OF REF. IMAGES
        IF (NUMREF <= 0)  THEN
           WRITE(NOUT,*)
     &     ' *** ERROR: Operation requires file template and number(s)'
           WRITE(6,*)
     &     '*** ERROR: Operation requires file template and number(s)'
           CALL ERRT(101,'No reference images',IDUM)
           GOTO 9999
        ENDIF

C       GET FIRST REFERENCE IMAGE TO DETERMINE DIMENSIONS
        IF (CTYPE(1:2) == 'OR') THEN
           FILNAM = REFPAT
        ELSE
           NLET = 0
           CALL FILGET(REFPAT,FILNAM,NLET,INUMBR(1),INTFLG)
        ENDIF

        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FILNAM,LUNREF,'O',IFORM,NX,NY,NZ,
     &               MAXIM,' ',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9999
        CLOSE(LUNREF)

        ISHRANGE = 1
        IVAL2    = 999999
        IVAL3    = 999999
        IF (CTYPE(1:2) == 'SH'  .OR.
     &      CTYPE(1:2) == 'OR') THEN
	   CALL RDPRI3S(ISHRANGE,IVAL2,IVAL3,NOT_USED,
     &     'TRANSLATION SEARCH RANGE IN X, IN Y (OPTIONAL), STEP SIZE',
     &      IRTFLG)
           IF (IRTFLG .NE. 0)  GOTO 9999
	   ISHRANGE   = MAX(ISHRANGE,0)       
	   ISHRANGEX  = ISHRANGE           ! _ & 5
           IF (IVAL2 == 999999) THEN       ! 5
              ISHRANGEY = ISHRANGEX
              ISTEP     = 1
           ELSEIF (IVAL3 == 999999) THEN   ! 5,6
              ISHRANGEY = ISHRANGEX
              ISTEP     = IVAL2
           ELSE                            ! 5,3,1
              ISHRANGEY = MAX(0,IVAL2)
              ISTEP     = IVAL3
           ENDIF
	   ISTEP = MAX(ISTEP,1)

C          CHECK SEARCH RANGE AND STEP SIZE.
	   IF (ISHRANGEX  > NX/2-2)  THEN
	      CALL ERRT(102,'X SEARCH MUST BE LESS THAN',NX/2-2)
              GOTO 9999
	   ELSEIF (ISHRANGEY  > NX/2-2)  THEN
	      CALL ERRT(102,'Y SEARCH MUST BE LESS THAN',NX/2-2)
              GOTO 9999
	   ELSEIF (ISHRANGEX > 0 .AND.
     &             MOD(ISHRANGEX,ISTEP) .NE. 0) THEN
	      CALL ERRT(102,'X SEARCH RANGE MUST BE DIVISIBLE BY',ISTEP)
              GOTO 9999
	   ELSEIF (ISHRANGEY > 0 .AND.
     &        MOD(ISHRANGEY,ISTEP) .NE. 0) THEN
	      CALL ERRT(102,'Y SEARCH RANGE MUST BE DIVISIBLE BY',ISTEP)
              GOTO 9999
	   ENDIF

        ELSEIF ( CTYPE(1:3) == 'FOU' ) THEN
	   CALL RDPRIS(ISHRANGE,IVAL2,NOT_USED,
     &     'TRANSLATION SEARCH RANGE IN X AND IN Y (OPTIONAL)',IRTFLG)
           IF (IRTFLG .NE. 0)  GOTO 9999
	   ISHRANGE   = MAX(ISHRANGE,1)       
	   ISHRANGEX  = ISHRANGE            ! _ & 5
           IF (IVAL2 == 999999) THEN        ! 5
              ISHRANGEY = ISHRANGEX
           ELSE                             ! 5,3,1
              ISHRANGEY = MAX(1,IVAL2)
           ENDIF
	   ISTEP = 1

C          CHECK SEARCH RANGE.
	   IF (ISHRANGEX  > NX/2-2)  THEN
	      CALL ERRT(102,'X SEARCH MUST BE LESS THAN',NX/2-2)
              GOTO 9999
	   ELSEIF (ISHRANGEY  > NX/2-2)  THEN
	      CALL ERRT(102,'Y SEARCH MUST BE LESS THAN',NX/2-2)
              GOTO 9999
	   ENDIF

        ELSEIF ( CTYPE(1:3) == 'REF' ) THEN
           CALL RDPRI1S(ISHRANGE,NOT_USED,
     &         'TRANSLATION SEARCH RANGE (ZERO FOR NONE)' ,IRTFLG)
           IF (IRTFLG .NE. 0)  GOTO 9999
        ENDIF

        IRAY   = 1
        RAY1   = 0               ! FIRST RAY ANGLE
        RAY2   = 360             ! LAST  RAY ANGLE

        MR     = 5
        MAXRAD = MIN( (NX - (NX/2+1)),(NY - (NY/2+1)) ) 
        NR     = MAXRAD - 2      ! DEFAULT VALUE

        !write(6,*) 'center:' , (NX/2+1),(NY/2+1)
        !write(6,*) 'nx,ny,maxrad:' ,nx,ny,maxrad 

	IF (CTYPE(1:2) == 'OR' .OR. 
     &      CTYPE(1:3) == 'FOU' ) THEN
           ISKIP = 1
           CALL RDPRIS(MR,NR,NOT_USED,'FIRST & LAST RING',IRTFLG)
           IF (IRTFLG .NE. 0)  GOTO 9999
 
        ELSE
           ISKIP     = 0
           ISKIP     = 0
           ISKIP     = 0

           VALUES(1) = MR
           VALUES(2) = NR
           VALUES(3) = ISKIP
           VALUES(4) = IRAY
           VALUES(5) = RAY1
           VALUES(6) = RAY2

           CALL RDPRA('FIRST, LAST RING, RING STEP & RAY STEP',
     &                6,0,.TRUE.,VALUES,NGOT,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (NGOT  > 0) THEN
C             COPY THE RETURNED VALUES 
              IF (NGOT >= 1) MR    = VALUES(1)
              IF (NGOT >= 2) NR    = VALUES(2)
              IF (NGOT >= 3) ISKIP = VALUES(3)
              IF (NGOT >= 4) IRAY  = VALUES(4)
              IF (NGOT >= 5) RAY1  = VALUES(5)
              IF (NGOT >= 6) RAY2  = VALUES(6)
           ENDIF
  
           IF ( CTYPE(1:2) .NE. 'SH' .AND. ISKIP <= 0) THEN
              ISKIP = 1
              CALL RDPRI1S(ISKIP,NOT_USED,'RING STEP',IRTFLG)
           ENDIF
        ENDIF
        IF (IRTFLG .NE. 0) GOTO 9999
        ISKIP = MAX(1,ISKIP)
        IF (IRAY .NE. 1 .AND. IRAY .NE. 2 .AND. 
     &      IRAY .NE. 4 .AND. IRAY .NE. 8) THEN
	   CALL ERRT(101,'RAY STEP MUST BE 1,2,4, OR 8',NE)
           GOTO 9999
        ENDIF

	IF (MR <= 0) THEN
	   CALL ERRT(101,'FIRST RING MUST BE > 0',NE)
	   GOTO 9999

	ELSEIF (NR < MR)  THEN 
	   CALL ERRT(102,'LAST RING MUST BE > ',MR)
	   GOTO 9999

	ELSEIF (CTYPE(1:2) == 'SH') THEN
C          CHECK SEARCH RANGE AND STEP SIZE TOGETHER
	   IF ((ISHRANGEX+NR)  >= MAXRAD)  THEN
	      CALL ERRT(102,
     &          'LAST RING + X TRANSLATION MUST BE LESS THAN',MAXRAD)
              GOTO 9999
	   ELSEIF ((ISHRANGEY+NR)  >= MAXRAD)  THEN
	      CALL ERRT(102,
     &          'LAST RING + Y TRANSLATION MUST BE LESS THAN',MAXRAD)
              GOTO 9999
	   ENDIF
	ELSEIF (NR >= MAXRAD)  THEN 
	   CALL ERRT(102,'LAST RING MUST BE LESS THAN',MAXRAD)
	   GOTO 9999
        ENDIF        

        REFANGDOC = NULL
        IF (CTYPE(1:3) == 'REF' .OR. 
     &      CTYPE(1:3) == 'FOU' .OR.
     &      CTYPE(1:2) == 'SH') THEN
C          GET NAME OF REFERENCE IMAGES ANGLES DOCUMENT FILE
           CALL FILERD(REFANGDOC,NREFA,NULL,
     &		'REFERENCE IMAGES ANGLES DOCUMENT',IRTFLG)
C          FILERD WILL RETURN IRTFLG=-1 IF "*" !!!!
           IF (NREFA <= 0) REFANGDOC = NULL
         ENDIF

C        FIND NUMBER OF REFERENCE-RINGS IN USE
         NRING=0
         DO I=MR,NR,ISKIP
            NRING = NRING + 1
	 ENDDO

         ALLOCATE(NUMR(3,NRING),STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'APMASTER; NUMR',3*NRING)
            GOTO 9999
         ENDIF

C        INITIALIZE NUMR ARRAY WITH RING RADII
         NRING = 0
         DO I=MR,NR,ISKIP
            NRING         = NRING + 1
            NUMR(1,NRING) = I
	 ENDDO

C        CALCULATES NUMR & LCIRC, EACH RING HAS FFT PAD OF 2 FLOATS
         CALL ALPRBS_Q_NEW(NUMR,NRING,LCIRC,MODE,IRAY)

C        FIND NUMBER OF OMP THREADS
         CALL GETTHREADS(NUMTH)

         NPLANS = 42 + 2   ! I THOUGHT 14 WAS ENUFF!!
         ALLOCATE(FFTW_PLANS(NPLANS),STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'APMASTER; FFTW_PLANS',NPLANS)
            GOTO 9999
         ENDIF

         IASK8 = (LCIRC * NUMREF)*4
         CALL BIGALLOC(IASK8,IOK,.FALSE.,.FALSE.,IRTFLG)

         GPU = (CTYPE(1:4) == 'SH G')
         IF (GPU) THEN
C           GPU MUST USE IN_CORE
            CIRCREF_IN_CORE = .TRUE.
            
         ELSEIF (CTYPE(1:4) == 'REFF' .OR. 
     &           CTYPE(1:3) == 'SHF'  .OR. 
     &           CTYPE(1:1) == 'I' ) THEN
C           INITIATE NON-INCORE EVEN IF SIZE IS OK

            ALLOCATE(CIRCREF(LCIRC,NUMTH),STAT=IRTFLG)
            IF (IRTFLG .NE. 0) THEN
               IF (MYPID <= 0) WRITE(NOUT,92) LCIRC,NUMTH
92             FORMAT ('  CAN NOT ALLOCATE: CIRCREF(',I8,' X ',I8,')') 
	       CALL ERRT(46,'APMASTER; CIRCREF',LCIRC*NUMTH)
	       GOTO 9999
            ENDIF
            IF (MYPID <= 0) WRITE(NOUT,93) LCIRC,NUMTH
93          FORMAT ('  DISK BASED RINGS FILE, ALLOCATED: CIRCREF(',I8,
     &                 ' X ',I8,')') 

            CIRCREF_IN_CORE = .FALSE.
         ELSE

            ALLOCATE(CIRCREF(LCIRC,NUMREF),STAT=IRTFLG)

            NTOT = LCIRC * NUMREF
            IF (IRTFLG == 0) THEN 
C              CIRCREF ALLOCATION SUCCEEDED 
               CIRCREF_IN_CORE = .TRUE.
               IF (MYPID <= 0) WRITE(NOUT,91) LCIRC,NUMREF,NTOT
91             FORMAT('  ALLOCATED: CIRCREF(',I8,' X ',I8,'): ',I10) 

	    ELSE
               CIRCREF_IN_CORE = .FALSE.
               IF (MYPID <= 0) WRITE(NOUT,90) LCIRC,NUMREF,NTOT
90             FORMAT('  CAN NOT ALLOCATE: CIRCREF(',I8,' X ',I8,'): ',
     &            I10,'  WILL USE REFERENCE-RINGS FILE.',/,
     &                '  MAY BE VERY SLOW! ',  
     &                'ADVISE YOU USE FEWER REFERENCES, INSTEAD.',/) 
 
C              GWP - HAVE TO FIX THE ALLOCATION HERE FOR DEC UNIX
               IF (ALLOCATED(CIRCREF)) DEALLOCATE(CIRCREF)
	       ALLOCATE(CIRCREF(LCIRC,NUMTH),STAT=IRTFLG)
	       IF (IRTFLG .NE. 0) THEN
	          CALL  ERRT(46,'APMASTER; CIRCREF',LCIRC*NUMTH)
	          GOTO 9999
	       ENDIF
            ENDIF
         ENDIF

         IF ((CTYPE(1:2) .EQ. 'AL') .OR. 
     &       (CTYPE(1:2) .NE. 'SH'  .AND.
     &        CTYPE(1:3) .NE. 'FOU' .AND.
     &        CTYPE(1:2) .NE. 'OR') .OR.
     &       (CTYPE(1:2) .EQ. 'SH'  .AND. .NOT. CIRCREF_IN_CORE) .OR.
     &       (CTYPE(1:3) .EQ. 'FOU' .AND. .NOT. CIRCREF_IN_CORE)) THEN
C
C           ~9 IS TO ACCEPT EXTENSION IF FILE IS NAMED
            CALL FILERD(ASK,NA,NULL,'REFERENCE-RINGS~9',IRTFLG)
            IF (IRTFLG .NE. 0)  GOTO 9999

            SCRFILE = ASK
            IF (ASK(1:NA) == 'W') THEN
               CALL ERRT(101,
     &         'OBSOLETE, USE: <AP I> TO CREATE REFERENCE-RINGS FILE',N)
               GOTO 9999
            ELSEIF (NA <= 3 .AND. ASK(1:1) == 'N' .AND. 
     &              .NOT. CIRCREF_IN_CORE) THEN
               CALL ERRT(101,
     &         'OBSOLETE, USE: <AP I> TO CREATE REFERENCE-RINGS FILE',N)
               GOTO 9999
            ELSEIF (NA <= 3 .AND. ASK(1:1) == 'Y') THEN
               SCRFILE = 'SCRATCH.file'
               IF (MYPID <= 0) WRITE(NOUT,*) 
     &               'OBSOLETE, GIVE NAME FOR REFERENCE-RINGS FILE'
            ENDIF 
         ELSE
            SCRFILE = CHAR(0) 
         ENDIF


         IF (CTYPE(1:1) == 'I' .OR. 
     &       CTYPE(1:2) == 'MI') THEN
C           ----------------- 'I' ------------------------ APRINGS_FILL
C          CREATE REFERENCE RINGS FILE FOR OUTPUT
           NSL = 1
           CALL OPFILEC(0,.FALSE.,SCRFILE,LUNRING,'B',IFORM,
     &               LCIRC,NUMREF,NSL,MAXIM,' ', .FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9990

           CALL RDPRMC(YN,NLET,.TRUE.,'WEIGHT THE RINGS? (Y/N)',
     &                 NULL,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9990
           WEIGHT = (YN == 'Y' .OR. YN == 'y')

C          INITIALIZE FFTW3 PLANS FOR USE WITHIN OMP || SECTIONS
           CALL APRINGS_INIT_PLANS(NUMR,NRING,
     &                          FFTW_PLANS,NPLANS,NX,NY,IRTFLG)

           CALL APRINGS_FILL_NEW(INUMBR,NUMREF,
     &                  NX,NY,NUMTH,
     &                  NRING,LCIRC,NUMR,MODE,FFTW_PLANS, 
     &                  REFPAT,LUNREF,
     &                  CIRCREF,NUMTH,LUNRING,WEIGHT,.FALSE.,IRTFLG)

           CLOSE(LUNRING)
           GOTO 9989
        ENDIF

399     IF (CTYPE(1:2) == 'OR') THEN
C           GET NAME OF SINGLE IMAGE TO BE ALIGNED
            MAXIM = 0
            CALL OPFILEC(0,.TRUE.,EXPPAT,LUNEXP,'O',IFORM,
     &                  NX,NY,NZ,MAXIM,'EXPERIMENTAL IMAGE',
     &                  .FALSE.,IRTFLG)
	    IF (IRTFLG .NE. 0) GOTO 9999

	    ALLOCATE(IMGLST(1),STAT=IRTFLG)
	    IF (IRTFLG .NE. 0) THEN
                CALL ERRT(46,'APMASTER; IMGLST',1)
                GOTO 9999
            ENDIF

            IMGLST(1) = 0
            NUMEXP    = 1
        ELSE
C           GET LIST OF EXPERIMENTAL IMAGES TO BE ALIGNED
	    ALLOCATE(IMGLST(NILMAX),STAT=IRTFLG)
	    IF (IRTFLG .NE. 0) THEN
               CALL ERRT(46,'APMASTER; IMGLST',NILMAX)
               GOTO 9999
            ENDIF

	    CALL FILELIST(.TRUE.,LUNEXP,EXPPAT,NLEP,
     &         IMGLST,NILMAX,NUMEXP,
     &         'TEMPLATE FOR IMAGE SERIES TO BE ALIGNED',IRTFLG)
	    IF (IRTFLG .NE. 0) GOTO 9999

            IF (MYPID <= 0) WRITE(NOUT,2002) NUMEXP
2002        FORMAT('  Number of experimental images: ',I6)
        ENDIF   


        EXPANGDOC = NULL
        IF (CTYPE(1:2) == 'SH'  .OR. 
     &      CTYPE(1:3) == 'REF' .OR.
     &      CTYPE(1:3) == 'FOU') THEN

C          GET NAME OF DOC FILE HOLDING EXPERIMENTAL IMAGES ANGLES
           CALL FILERD(EXPANGDOC,NEXPA,NULL,
     &                 'EXPERIMENTAL IMAGES ALIGNMENT DOCUMENT',IRTFLG)
           IF (NEXPA == 0) EXPANGDOC = NULL
        ENDIF

        RANGE     = 0.0
        ANGDIFTHR = 0.0
        IF (CTYPE(1:3) == 'REF' .OR.
     &      CTYPE(1:3) == 'FOU' .OR.
     &     (CTYPE(1:2) == 'SH'  .AND. .NOT. GPU) ) THEN

           CALL RDPRM2S(RANGE,ANGDIFTHR,NOT_USED,
     &      'RANGE OF PROJECTION ANGLE SEARCH & ANGLE CHANGE THRESHOLD',
     &       IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (RANGE  > 0.0 .AND. EXPANGDOC == NULL) THEN
             CALL ERRT(101,
     &       'MUST SPECIFY EXPERIMENTAL IMAGES ALIGNMENT DOCUMENT FILE',
     &       IDUM)
             RANGE = 0.0
           ENDIF
           RANGE = MIN(RANGE,180.0)  ! 3D MAX. DIFFERENCE POSSIBLE
       ENDIF

       DENOISE    = .FALSE.
       ROTFIRST   = .FALSE.
       CKMIRROR   = .TRUE.
       FBS_WANTED = .FALSE.

       IF ((CTYPE(1:3) == 'REF') .OR.
     &      CTYPE(1:3) == 'FOU'  .OR.
     &      CTYPE(1:3) == 'ORS'  .OR.
     &     (CTYPE(1:2) == 'SH'   .AND. .NOT. GPU)) THEN

          IF (CTYPE(1:3) == 'SHC') THEN
C                      123456789 123456789 1234567890
             PROMPT = 'CHECK MIRRORED POSITIONS?, ' //
     &                'SHIFT/ROTATE INPUT?, '       //
     &                'DENOISE? (Y/N)' 
             CALL RDPRMC(ASK,NA,.TRUE.,PROMPT(1:62),NULL,IRTFLG)
             IF (IRTFLG .NE. 0) GOTO 9999

             GOTMIR     = .FALSE.
             GOTRTSH    = .FALSE.
             CKMIRROR   = .FALSE.
             ROTFIRST   = .FALSE.
             FBS_WANTED = .FALSE.
             DENOISE    = .FALSE.

             DO I = 1,NA
                CTEMP = ASK(I:I)
                IF (CTEMP == '1' .OR. 
     &              CTEMP == 'Y' .OR. 
     &              CTEMP == 'Q'.OR.
     &              CTEMP == 'F') THEN
                   IF (.NOT. GOTMIR) THEN
                      CKMIRROR   = .TRUE.
                      GOTMIR     = .TRUE.
        !write(6,*) '1',ctemp,CKMIRROR,gotrtsh,fbs_wanted,denoise

                   ELSEIF (.NOT. GOTRTSH) THEN 
                      GOTRTSH    = .TRUE.
                      ROTFIRST   = .TRUE.
                      IF (CTEMP == 'F' ) FBS_WANTED = .TRUE. 
        !write(6,*) '2',ctemp,CKMIRROR,gotrtsh,fbs_wanted,denoise

                    ELSE  
                      IF (CTEMP == 'Y') DENOISE = .TRUE. 
        !write(6,*) '3',ctemp,CKMIRROR,gotrtsh,fbs_wanted,denoise
                      EXIT
                   ENDIF

                ELSEIF (CTEMP == '0' .OR. 
     &                  CTEMP == 'N') THEN
                   IF (.NOT. GOTMIR) THEN
                      CKMIRROR = .FALSE.
                      GOTMIR   = .TRUE.
                   ELSEIF (.NOT. GOTRTSH) THEN 
                      ROTFIRST = .FALSE.
                      GOTRTSH  = .TRUE.
                   ELSE  
                      DENOISE  = .FALSE.
                      EXIT
                   ENDIF
                ENDIF
             ENDDO
 
          ELSEIF ((CTYPE(1:3) == 'REF') .OR.
     &             CTYPE(1:3) == 'FOU'  .OR.
     &             CTYPE(1:2) == 'SH') THEN
C                      123456789 123456789 1234567890
             PROMPT = 'CHECK MIRRORED POSITIONS?, ' //
     &                'SHIFT/ROTATE INPUT? (Y/N)' 
             CALL RDPRMC(ASK,NA,.TRUE.,PROMPT(1:52),NULL,IRTFLG)
             IF (IRTFLG .NE. 0) GOTO 9999

             CKMIRROR = .FALSE.
             GOTMIR   = .FALSE.
             GOTRTSH  = .FALSE.
             DO I = 1,NA
                CTEMP = ASK(I:I)
                IF (CTEMP == '1' .OR. 
     &              CTEMP == 'Y' .OR. 
     &              CTEMP == 'Q'.OR.
     &              CTEMP == 'F') THEN
                   IF (.NOT. GOTMIR) THEN
                      CKMIRROR   = .TRUE.
                      GOTMIR     = .TRUE.
                   ELSEIF (.NOT. GOTRTSH) THEN 
                      GOTRTSH    = .TRUE.
                      ROTFIRST   = .TRUE.
                      IF (CTEMP == 'Q') FBS_WANTED = .FALSE. 
                   ENDIF

                ELSEIF (CTEMP == '0' .OR. 
     &                  CTEMP == 'N') THEN
                   IF (.NOT. GOTMIR) THEN
                      CKMIRROR = .FALSE.
                      GOTMIR   = .TRUE.
                   ELSEIF (.NOT. GOTRTSH) THEN 
                      ROTFIRST = .FALSE.
                      GOTRTSH  = .TRUE.
                   ENDIF
                ENDIF
             ENDDO
 
          ELSE

             CALL RDPRMC(ASK,NA,.TRUE.,
     &          'CHECK MIRRORED POSITIONS? (Y/N)', NULL,IRTFLG)
             IF (IRTFLG .NE. 0) GOTO 9999

             DO I = 1,NA
                CTEMP = ASK(I:I)
                IF (CTEMP == '1' .OR. CTEMP == 'Y') THEN
                   CKMIRROR = .TRUE.
                   EXIT
                ELSEIF (CTEMP == '0' .OR. CTEMP == 'N') THEN
                   CKMIRROR = .FALSE.
                   EXIT
                ENDIF
             ENDDO
          ENDIF

       ELSEIF (GPU) THEN
C          MIRRORING TOO COMPLEX LOGIC FOR GPU
           CKMIRROR = .FALSE.
       ENDIF

C      GET NAME FOR OUTPUT DOC FILE
       CALL REG_GET_USED(NSEL_USED)
       LENC = lnblnkn(CTYPE)

       WANTDOC = .TRUE.
       IF  (FCHAR(1:2) == 'OR') WANTDOC = .FALSE.
       IF  (NSEL_USED  > 0)       WANTDOC = .FALSE.
       IF ((NSEL_USED  > 0) .AND.
     &     (LENC >= 4 .AND. CTYPE(4:4) == 'D') .OR.
     &     (LENC >= 5 .AND. CTYPE(5:5) == 'D')) WANTDOC = .TRUE.
       !write(6,*) 'wantdoc: ',wantdoc,lenc,ctype

        IF (.NOT. WANTDOC) THEN
C        NO OUTPUT FILE WANTED
         OUTANG  = NULL
         NOUTANG = 0
       ELSE
C        OPEN OUTPUT DOC FILE (FOR APPENDING)
         NOUTANG = NDOC
         CALL OPENDOC(OUTANG,.TRUE.,NLET,NDOC,NOUTANG,.TRUE.,
     &           'OUTPUT ALIGNMENT DOCUMENT',.FALSE.,.TRUE.,.TRUE.,
     &            NEWFILE,IRTFLG)
         IF (IRTFLG == -1) THEN
C           DO NOT WANT OUTPUT DOC FILE
            NOUTANG = 0
         ELSEIF (IRTFLG .NE. 0) THEN
            RETURN
         ELSE
C           WANT OUTPUT DOC FILE
            IF (CTYPE(1:2) == 'SH'  .OR. 
     &          CTYPE(1:3) == 'REF' .OR.
     &          CTYPE(1:3) == 'FOU' ) THEN
               IF (USE_LONGCOL) THEN    ! FROM CMBLOCK.INC
               COMMEN ='                 '                           //
     &         ' PSI,          THE,          PHI,         REF#,     '//
     &         '   EXP#,    CUM. {INPLANE,      SX,          SY},  ' //
     &         '        NPROJ,         DIFF,        CCROT,  '        //
     &         '     INPLANE,          SX,            SY,       MIR-CC'


               ELSE
               COMMEN =' KEY       '              //
     &         'PSI,    THE,    PHI,   REF#,    ' //
     &         'EXP#,  CUM.{ROT,   SX,    SY},  ' //
     &         'NPROJ,   DIFF,      CCROT,    '   //
     &         'ROT,     SX,     SY,   MIR-CC'
               ENDIF

            ELSE
               COMMEN = '      ' //
     &         'MIR-REF#,     CCROT,     INPLANE <,      SX,         '//
     &         'SY,           IMG#,       < DIFF'
            ENDIF
            CALL LUNDOCPUTCOM(NOUTANG,COMMEN,IRTFLG)
         ENDIF
       ENDIF


C        -------PARAMETER INPUT FINISHED, CALCULATE NOW ---------


C       INITIALIZE FFTW3 PLANS FOR USE WITHIN OMP || SECTIONS
        CALL APRINGS_INIT_PLANS(NUMR,NRING,
     &                          FFTW_PLANS,NPLANS,NX,NY,IRTFLG)

        IF (CTYPE(1:3) == 'FOU') THEN
C          --------------------'FOU'   ------------------------ 'AP FOU'
 
          IF (MYPID <= 0) WRITE(NOUT,*) 
     &           ' Calling: APFOU_PATM FOR: ',CTYPE(1:4),' -----------'
          CALL FLUSHFILE(6)

          CALL APFOU_PATM(INUMBR,NUMREF, IMGLST,NUMEXP, 
     &                NX,NY, RANGE,ANGDIFTHR,
     &                NRING,LCIRC,NUMR,CIRCREF,CIRCREF_IN_CORE,
     &                REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS,
     &                REFPAT,EXPPAT, CKMIRROR,CTYPE,
     &                ROTFIRST,ISHRANGEX,ISHRANGEY,
     &                NOUTANG,FBS_WANTED)

         ELSEIF (CTYPE(1:3) == 'REF') THEN
C           --------------------'REF'   ----------------------- 'AP REF'
 
            IF ((CIRCREF_IN_CORE .AND. 
     &           NUMTH > 1 .AND. NUMEXP > NUMTH) .OR. 
     &           CTYPE(3:3) == 'T' .OR. CTYPE(4:4) == 'T') THEN

               IF (MYPID <= 0) WRITE(NOUT,*) 
     &           ' Calling: APREF_PM FOR: ',CTYPE(1:3),' -----------'
               CALL FLUSHFILE(6)

               CALL APREF_PM(INUMBR,NUMREF,IMGLST,NUMEXP,
     &                NX,NY, ANGDIFTHR,RANGE,
     &                NRING,LCIRC,NUMR,CIRCREF,
     &                MODE,REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS,
     &                REFPAT,EXPPAT,
     &                CKMIRROR,CTYPE,ROTFIRST,ISHRANGE,
     &                NOUTANG,FBS_WANTED)
            ELSE 
               IF (MYPID <= 0) WRITE(NOUT,*) 
     &            ' Calling: APREF_P FOR: ',CTYPE(1:3),' -----------'
               CALL FLUSHFILE(6)


               CALL APREF_P(INUMBR,NUMREF,IMGLST,NUMEXP,
     &               NX,NY,RANGE,ANGDIFTHR,
     &               NRING,LCIRC,NUMR,CIRCREF,CIRCREF_IN_CORE,
     &               MODE,REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS,
     &               REFPAT,EXPPAT,
     &               CKMIRROR,CTYPE,ROTFIRST,ISHRANGE,
     &               NOUTANG,FBS_WANTED)
           ENDIF 

        ELSEIF (CTYPE(1:2) == 'SH' .OR.
     &          CTYPE(1:2) == 'AL' .OR.
     &          CTYPE(1:2) == 'OR') THEN
C        ---- ' SH', 'ORS',  ---------------------------------- 'AP_SH'

#ifdef MPI_DEBUG 
              T0 = MPI_WTIME()
#endif

	   IF (GPU) THEN
C             FOR GPU ACCELERATED, USES COEF. TO SPEED UP APRINGS
 
              IF (MYPID <= 0) WRITE(NOUT,*)
     &             ' Calling: APSH_CUDA FOR: ',CTYPE(1:3),' ----------'
              IF (MYPID <= 0) WRITE(6,*)
     &             ' Calling: APSH_CUDA FOR: ',CTYPE(1:3),' ----------'
              CALL FLUSHFILE(6)

              CALL APSH_CUDA(INUMBR,NUMREF, IMGLST,NUMEXP, 
     &               NX,NY,ISHRANGEX,ISHRANGEY,ISTEP,
     &               NUMR,NRING,
     &               MODE, REFANGDOC,EXPANGDOC,FFTW_PLANS,
     &               REFPAT,EXPPAT,NOUTANG)

           ELSEIF (CIRCREF_IN_CORE .AND. (CTYPE(3:3) == 'C')) THEN
C             FOR MP, LARGE NUMBER OF IMAGES TO BE ALIGNED, OR SP.
C             USES COEF. TO SPEED UP APRINGS          
              IF (MYPID <= 0) WRITE(NOUT,*)
     &           ' Calling: APSH_PSC FOR: ',CTYPE(1:5),' ----------'

C             ECHO SEARCH PARAMETERS
              NSHIFTSX = (2*ISHRANGEX + 1) 
              NSHIFTSY = (2*ISHRANGEY + 1)
              NSHIFTS  = NSHIFTSX *  NSHIFTSY 

              MSG = ' Aligning'
              NA  = 9 
              IF (ROTFIRST) THEN
                 IF (FBS_WANTED) THEN
                    MSG = MSG(1: NA) //  ' FBS interpolated, rotated'
                    NA = NA + 26
                 ELSE
                    MSG = MSG(1: NA) //  ' QUAD interpolated, rotated'
                    NA = NA + 27
                 ENDIF
              ELSE
                 MSG = MSG(1: NA)    //  ' unrotated'
                 NA = NA + 10

              ENDIF
              IF (DENOISE) THEN
                 MSG = MSG(1: NA)    // ', denoised'
                 NA = NA + 10
              ENDIF 
              MSG = MSG(1: NA)       // ' exp images'
              NA = NA + 11
              IF (MYPID <= 0) WRITE(NOUT,*) MSG(1:NA)

              IF (MYPID <= 0) WRITE(NOUT,981)-ISHRANGEX,ISHRANGEX,
     &                                       -ISHRANGEY,ISHRANGEY,
     &                                        ISTEP    !,NSHIFTS
981           FORMAT('  Shift search X range:      ',I3,' ...',I3,
     &               '     Y range:',I3,'...',I3, 
     &               '    Step size:',I3)
              IF (MYPID <= 0) WRITE(NOUT,982)MR,NR,ISKIP
982           FORMAT('  Ring search radius range:  ',I3,' ...',I4,
     &               '    Skip increment:',I3)
              IF (MYPID <= 0) WRITE(NOUT,983)RAY1,RAY2
983           FORMAT('  Rotational search range:',F6.1,' ...',F6.1,/)
              CALL FLUSHRESULTS

              CALL APSH_PSC(INUMBR,NUMREF,IMGLST,NUMEXP, 
     &               NX,NY,ISHRANGEX,ISHRANGEY,ISTEP,
     &               NRING,LCIRC,NUMR,CIRCREF,CIRCREF_IN_CORE,
     &               MODE, REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS,
     &               REFPAT,EXPPAT,RANGE,ROTFIRST,DENOISE,
     &               CKMIRROR,CTYPE,NOUTANG,RAY1,RAY2,FBS_WANTED)

	   ELSEIF (CIRCREF_IN_CORE      .AND. 
     &             NUMEXP >= NUMTH      .AND. 
     &             CTYPE(3:3) .NE. 'T'  .AND. 
     &             CTYPE(4:4) .NE. 'T'  .AND.
     &             CTYPE(1:2) .NE. 'OR') THEN
C             FOR MP, LARGE NUMBER OF IMAGES TO BE ALIGNED, OR SP.

              IF (MYPID <= 0) WRITE(NOUT,*)
     &           ' Calling: APSH_PS FOR: ',CTYPE(1:3),' ----------'
              CALL FLUSHRESULTS

C             ECHO SEARCH PARAMETERS
              NSHIFTSX = (2*ISHRANGEX + 1) 
              NSHIFTSY = (2*ISHRANGEY + 1)
              NSHIFTS  = NSHIFTSX *  NSHIFTSY 
              IF (MYPID <= 0) WRITE(NOUT,*)
              IF (MYPID <= 0) WRITE(NOUT,981)-ISHRANGEX,ISHRANGEX,
     &                                       -ISHRANGEY,ISHRANGEY,
     &                                       ISTEP    
              IF (MYPID <= 0) WRITE(NOUT,982)MR,NR,ISKIP
              IF (MYPID <= 0) WRITE(NOUT,983)RAY1,RAY2

              CALL APSH_PS(INUMBR,NUMREF,IMGLST,NUMEXP, 
     &               NX,NY,ISHRANGEX,ISHRANGEY,ISTEP,
     &               NRING,LCIRC,NUMR,CIRCREF,CIRCREF_IN_CORE,
     &               MODE, REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS,
     &               REFPAT,EXPPAT,RANGE,ROTFIRST,
     &               CKMIRROR,CTYPE,NOUTANG,FBS_WANTED)

	   ELSE
C             USE DIFFERENT STRATEGY FOR SMALL NUMBER OF SAMPLE IMAGES  
C             OR FOR CIRCREF FROM RINGS FILE  
C             TO MAKE MP EFFICIENT. ALSO FOR 'OR' OPERATIONS
              IF (CTYPE(1:2) .NE. 'OR' .AND. MYPID <= 0) WRITE(NOUT,*) 
     &           ' Calling: APSH_SS FOR: ',CTYPE(1:3),' ----------'
              CALL FLUSHFILE(6)

              CALL APSH_SS(INUMBR,NUMREF,IMGLST,NUMEXP, 
     &               NX,NY,ISHRANGEX,ISHRANGEY,ISTEP,
     &               NRING,LCIRC,NUMR,CIRCREF,CIRCREF_IN_CORE,
     &               MODE, REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS, 
     &               REFPAT,EXPPAT,RANGE,ROTFIRST,
     &               CKMIRROR,CTYPE,NOUTANG,FBS_WANTED)
 	   ENDIF

#ifdef MPI_DEBUG
           T1 = MPI_WTIME()
           T1 = T1 - T0
           IF (MYPID == 0)  WRITE(6, 222) T1
 222       FORMAT('  AP TIME: ', 1PE11.3)
#endif
         ENDIF

C        UNALLOCATE FFTW3 PLANS (TO REMOVE MEMORY LEAK)
9989     CALL FFTW3_KILLPLANS(FFTW_PLANS,NPLANS,IRTFLG)

9990     IF (MYPID <= 0 .AND. VERBOSE) WRITE (NOUT,2600)
2600     FORMAT (/' ',12('-'),' END OF COMPUTATION ',12('-')/)

9999     IF (ALLOCATED(IMGLST))     DEALLOCATE(IMGLST)
         IF (ALLOCATED(NUMR))       DEALLOCATE(NUMR)
         IF (ALLOCATED(CIRCREF))    DEALLOCATE(CIRCREF)
         IF (ALLOCATED(FFTW_PLANS)) DEALLOCATE(FFTW_PLANS)

         CLOSE(NDOC)

#ifdef USE_MPI_NEVER
         write(0,*) ' apmaster; at final barrier: ',mypid
         CALL MPI_BARRIER(ICOMM,MPIERR)
         write(0,*) ' apmaster; after final barrier: ',mypid
#endif

         END

