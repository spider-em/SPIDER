C++*********************************************************************
C                                                                      *
C APRINGS_NEW   USED CMLIMIT                       AUG 00 ARDEAN LEITH *
C               ADDED REF_CIRC FILE                APR 01 ARDEAN LEITH *
C               NORMASS -> NORMAS                  OCT 01 ARDEAN LEITH *
C               PROMPTS                            JAN 02 ARDEAN LEITH *
C               AUTO REF RINGS FILE                DEC 04 ARDEAN LEITH *
C               SPIDER REF RINGS FILE              FEB 05 ARDEAN LEITH *
C               INULL  = lnblnkn(SCRFILE)          APR 05 ARDEAN LEITH *
C               REWRITE FOR SPEED                  MAR 08 ARDEAN LEITH *
C               BCAST_MPI                          NOV 08 ARDEAN LEITH *
C               OUTPUT TO NOUT                     AUG 09 ARDEAN LEITH *
C               APRINGS_SATU                       OCT 10 ARDEAN LEITH *
C               ALRQ_MS_FBS                        AUG 11 G KISHCHENKO *
C               ALRQ_MS_QUAD                       AUG 11 ARDEAN LEITH *
C               TYPET                              MAY 12 ARDEAN LEITH
C **********************************************************************
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2011  Health Research Inc.,                         *
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
C  APRINGS_NEW(ILIST,NUMREF,NX,NY,
C              NRING,LCIRC,NUMR,MODE,FFTW_PLANS, 
C              REFPAT,LUNREF,CIRCREF,CIRCREF_IN_CORE,
C              LUNRING,SCRFILE,IRTFLG)
C  
C  PARAMETERS: ILIST       LIST OF IMAGE FILE NUMBERS          (INPUT)
C              NIMREF      NO. OF IMAGES                       (INPUT)
C              NX,NY   IMAGE DIMENSIONS                    (INPUT)
C              REFPAT      IMAGE SERIES FILE TEMPLATE          (INPUT)
C              LUNREF      IMAGE FILE IO UNIT                  (INPUT)
C              CIRCREF     OUTPUT ARRAY                        (OUTPUT)
C
C              CIRCREF_IN_CORE   NO OUTPUT ARRAY FLAG          (OUTPUT)
C              LUNRING      REF-RINGS FILE IO UNIT             (INPUT)
C              SCRFILE      REF-RINGS FILE                     (INPUT)
C              IRTFLG       ERROR FLAG                         (OUTPUT)
C
C NOTE: MOST MEMORY DEMAND DEPENDENT ON LCIRC & NUMREF.  LCIRC IS THE 
C       TOTAL LENGTH OF ARRAY THAT HOLDS THE CIRCULAR RINGS, SO IT IS 
C       DEPENDENT ON NUMBER OF RINGS AND THEIR RADIUS. ARRAY ALLOCATED 
C       IS: CIRCREF(LCIRC,NUMREF) ANOTHER SMALL ALLOCATED ARRAY IS: 
C       A(NX,NY,NUMTH)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

cpgi$g opt=3

C       -------------------- APRINGS_NEW ----------------------------

        SUBROUTINE APRINGS_NEW(ILIST,NUMREF, 
     &                     NX,NY,
     &                     NRING,LCIRC,NUMR,MODE,FFTW_PLANS, 
     &                     REFPAT,LUNREF,CIRCREF,CIRCREF_IN_CORE,
     &                     LUNRING,SCRFILE, IRTFLG)

        IMPLICIT NONE
	INCLUDE 'CMBLOCK.INC' 
	INCLUDE 'CMLIMIT.INC' 

        INTEGER                  :: ILIST(NUMREF)
        INTEGER                  :: NUMREF,NX,NY,NRING,LCIRC
        INTEGER                  :: NUMR(3,NRING) 
	CHARACTER(LEN=1)         :: MODE
        INTEGER*8                :: FFTW_PLANS(*) !POINTERS TO STRUCTURES
        CHARACTER (LEN=*)        :: REFPAT
        INTEGER                  :: LUNREF 
        REAL                     :: CIRCREF(LCIRC,*)
        LOGICAL                  :: CIRCREF_IN_CORE
        INTEGER                  :: LUNRING 
        CHARACTER (LEN=*)        :: SCRFILE
        INTEGER                  :: IRTFLG 

        CHARACTER (LEN=MAXNAM)   :: FILNAM
        LOGICAL                  :: USEREFFILE
        LOGICAL                  :: ISOPEN
        LOGICAL                  :: WINDOW,WEIGHT

        INTEGER                  :: ICOMM,MYPID,MPIERR
        INTEGER                  :: ILOCAT,INULL,NLET,NDUM,LUNOP
        INTEGER                  :: NUMTH,LCIRCT,NUMREFT,NSLICE,MAXIM
        INTEGER                  :: IMGNUM,NE,IREF
        REAL                     :: DUM
        LOGICAL                  :: INLNED
        CHARACTER(LEN=3)         :: TYPET = 'FI'
        INTEGER                  :: lnblnkn,lnblnk

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID 

C       SEE IF THIS FILE EXISTS, (RETURNS EX, ISOPEN, LUNOP)

C       SEE IF SCRATCH "FILE" EXISTS (MAY BE INCORE FILE)
        ILOCAT = INDEX(SCRFILE,'@')
        INULL  = lnblnkn(SCRFILE)
        IF (INULL .GT. 0 .AND. SCRFILE(1:1) .NE. '_' .AND. 
     &      ILOCAT .EQ. 0) THEN
C          ADD EXTENSION TO PHYSICAL FILENAME
           CALL FILNAMANDEXT(SCRFILE,DATEXC,FILNAM,NLET,.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ELSE
           FILNAM = SCRFILE
           NLET   = lnblnk(FILNAM)
        ENDIF
#ifdef USE_MPI
        IF (MYPID .LE. 0) THEN
           INQUIRE(FILE=FILNAM,EXIST=USEREFFILE,OPENED=ISOPEN,
     &             NUMBER=LUNOP,IOSTAT=IRTFLG)
        ENDIF
        CALL BCAST_MPI('APRINGS','USEREFFILE',USEREFFILE,1, 'L',ICOMM)
#else
        USEREFFILE = .FALSE. 
        IF (INULL .GT. 0) THEN
           CALL INQUIREIF1(LUNRING,FILNAM,TYPET,USEREFFILE,ISOPEN,
     &                     LUNOP,INLNED,IMGNUM,IRTFLG)
        ENDIF
#endif

        IF (USEREFFILE .AND. MYPID .LE. 0) THEN 
           WRITE(NOUT,*) ' Using existing reference rings file: ',
     &                     FILNAM(1:NLET)
        ELSEIF (INULL > 0  .AND. MYPID .LE. 0) THEN 
           WRITE(NOUT,*) ' No existing reference rings file: ',
     &                     FILNAM(1:NLET)
        ENDIF

        WEIGHT = .TRUE.   ! REFERENCES ARE WEIGHTED

C       NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)

        IF (CIRCREF_IN_CORE  .AND. .NOT. USEREFFILE) THEN
C          CALCULATE REF. RINGS DATA AND FILL CIRCREF ARRAY WITH IT

           CALL APRINGS_FILL_NEW(ILIST,NUMREF,  NX,NY,NUMTH,
     &                      NRING,LCIRC,NUMR,MODE,FFTW_PLANS, 
     &                      REFPAT,LUNREF,
     &                      CIRCREF,NUMREF,0,WEIGHT,
     &                      IRTFLG)

          !IF (MYPID .LE. 0) WRITE(NOUT,*) ' Created incore reference rings'

        ELSEIF (CIRCREF_IN_CORE .AND. USEREFFILE) THEN
C          READ EXISTING REF RINGS FILE AND FILL CIRCREF ARRAY WITH ITS DATA

C          OPEN EXISTING REFERENCE RINGS FILE
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,SCRFILE,LUNRING,'O',IFORM,
     &                 LCIRCT,NUMREFT,NSLICE,MAXIM,' ', .TRUE.,IRTFLG)
	   IF (IRTFLG .NE. 0) RETURN

           IF (LCIRCT .NE. LCIRC .OR. NUMREFT .NE. NUMREF) THEN
               CALL ERRT(101,'REF. RINGS FILE HAS WRONG SIZE',NE)
               IRTFLG = 1
               RETURN
           ENDIF

C          FILL CIRCREF WITH EXISTING RINGS DATA FROM FILE AND RETURN.
           DO IREF=1,NUMREF
              CALL REDLIN(LUNRING,CIRCREF(1,IREF),LCIRC,IREF) 
           ENDDO

           IF (MYPID .LE. 0) THEN
             WRITE(NOUT,*) ' Loaded reference rings file incore'
           ENDIF

        ELSEIF (.NOT. CIRCREF_IN_CORE .AND. USEREFFILE) THEN
C          OPEN EXISTING REF. RINGS FILE TO READ REF RINGS DATA
 
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,SCRFILE,LUNRING,'O',IFORM,
     &                 LCIRCT,NUMREFT,NSLICE,MAXIM,' ', .FALSE.,IRTFLG)
	   IF (IRTFLG .NE. 0) RETURN

           IF (LCIRCT .NE. LCIRC .OR. NUMREFT .NE. NUMREF) THEN
               CALL ERRT(101,'REF. RINGS FILE HAS DIFFERENT SIZE',NE)
               IRTFLG = 1
               RETURN
           ENDIF

        ENDIF

        END


C       --------------------- APRINGS_FILL_NEW --------------------------

        SUBROUTINE APRINGS_FILL_NEW(ILIST,NUMREF,
     &                        NX,NY,NUMTH,
     &                        NRING,LCIRC,NUMR,MODE,FFTW_PLANS, 
     &                        REFPAT,LUNREF,
     &                        CIRCREF,ICORE,LUNRING,WEIGHT,
     &                        IRTFLG)

        IMPLICIT NONE
	INCLUDE 'CMBLOCK.INC' 
	INCLUDE 'CMLIMIT.INC' 

        INTEGER              :: ILIST(NUMREF)
        INTEGER              :: NUMREF,NX,NY,NUMTH,NRING,LCIRC
        INTEGER              :: NUMR(3,NRING)
        CHARACTER(LEN=1)     :: MODE
        INTEGER*8            :: FFTW_PLANS(*) !POINTERS TO STRUCTURES
        CHARACTER (LEN=*)    :: REFPAT 
        INTEGER              :: LUNREF
        REAL                 :: CIRCREF(LCIRC,ICORE)
        INTEGER              :: ICORE,LUNRING
        LOGICAL              :: WEIGHT
        INTEGER              :: IRTFLG

C       AUTOMATIC ARRAYS
        REAL                 :: WR(NRING)

C       ALLOCATABLE ARRAYS
	REAL, ALLOCATABLE    :: BUFIMG(:,:,:)

        INTEGER              :: NXP,NYP,MAXRIN,MWANT,IGO,IEND
        INTEGER              :: IMI,IT,IPT,IV
        INTEGER              :: ICOMM,MYPID,MPIERR 
        REAL                 :: CNS2,CNR2,RDUM
        LOGICAL, PARAMETER   :: MPIBCAST   = .TRUE.
        LOGICAL, PARAMETER   :: WANTSTATS  = .FALSE.

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID 
        IRTFLG = 0

        !IF (MYPID .LE. 0)  WRITE(NOUT,'(A,i5)') 
        !&               '  Filling reference rings with threads:',NUMTH

        IF (WEIGHT) THEN
           MAXRIN = NUMR(3,NRING) - 2

C          RINGWE RETURNS WR WEIGHTS
	   CALL RINGWE_NEW(WR,NUMR,NRING,MAXRIN)
           IF (MODE .EQ. 'H') WR = WR * 0.5
        ELSE
           WR(1) = 0.0
        ENDIF

        MWANT = NX*NY*NUMTH
        ALLOCATE(BUFIMG(NX,NY,NUMTH),  STAT=IRTFLG)

        IF (IRTFLG .NE. 0) THEN
          CALL ERRT(46,'BUFIMG...',MWANT)
          RETURN
        ENDIF 

C       CALCULATE CENTER FOR APRINGS
	CNS2 = NX / 2 + 1
	CNR2 = NY / 2 + 1

        !write(6,*) ' numref,numth: ', numref,numth

C       PREPARE CIRCULAR RINGS DATA FOR ALL REFERENCE IMAGES
C       LOOP OVER ALL REF. IMAGES, CHUNKS OF NUMTH
        DO IGO=1,NUMREF,NUMTH  
           IEND = MIN(IGO+NUMTH-1,NUMREF)

C          LOAD SET OF REFERENCE IMAGES (THERE ARE NUMTH IMAGES IN SET)
           CALL AP_GETDATA(ILIST,NUMREF, 
     &                    NX,NY, NX,NY,0.0,
     &                    NUMTH,REFPAT,LUNREF, IGO,IEND,
     &                    MPIBCAST, BUFIMG, 
     &                    WANTSTATS,RDUM,RDUM,
     &                    IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

c$omp      parallel do private(imi,ipt,it)
	   DO IMI=IGO,IEND      
              IT  = IMI - IGO + 1           ! BUFIMG INDEX
              IPT = IMI                     ! CIRCREF INDEX
              IF (LUNRING .GT. 0) IPT = IT  

C             CONVERT BUFIMG TO POLAR RINGS, FFT, & WEIGHT THE RINGS
 	      CALL APRINGS_ONE_NEW(NX,NY, CNS2,CNR2, 
     &                 BUFIMG(1,1,IT),
     &                .FALSE.,MODE,NUMR,NRING,LCIRC,WR,FFTW_PLANS,
     &                 CIRCREF(1,IPT),IRTFLG)

c              write(6,*) 'circref(1,',ipt,'):',CIRCREF(1,IPT),it
           ENDDO

           IF (LUNRING .GT. 0) THEN
C             SAVE CIRCREF IN FILE OPENED ON LUNRING
	      DO IMI=IGO,IEND
                 IV = IMI - IGO + 1 
                 CALL WRTLIN(LUNRING,CIRCREF(1,IV),LCIRC,IMI)
              ENDDO
           ENDIF
        ENDDO

C       DEALLOCATE ARRAY
9999    IF (ALLOCATED(BUFIMG))   DEALLOCATE(BUFIMG)

	END


C       --------------------- APRINGS_ONE_NEW --------------------------

	SUBROUTINE APRINGS_ONE_NEW(NX,NY, CNS2,CNR2, XIM, USE_OMP,
     &                         MODE,NUMR,NRING,LCIRC,WR,FFTW_PLANS,
     &                         CIRC,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'

        INTEGER,   INTENT(IN)   :: NX,NY
        REAL,      INTENT(IN)   :: CNS2,CNR2

        REAL,      INTENT(INOUT):: XIM(NX,NY)
        LOGICAL,   INTENT(IN)   :: USE_OMP
        CHARACTER(LEN=1)        :: MODE
        INTEGER,   INTENT(IN)   :: NUMR(3,NRING)
        INTEGER,   INTENT(IN)   :: NRING,LCIRC
        REAL,      INTENT(OUT)  :: WR(*)
        INTEGER*8, INTENT(IN)   :: FFTW_PLANS(*)!POINTERS TO STRUCTURES
        REAL,      INTENT(OUT)  :: CIRC(LCIRC)
        INTEGER,   INTENT(OUT)  :: IRTFLG

        LOGICAL                 :: WEIGHTIT
        INTEGER                 :: NSB,NSE,NRB,NRE,NXLD,MAXRING
        INTEGER                 :: I,J,IGOM1,NRAYS,IRING
        DOUBLE PRECISION        :: AVG,SIGINV
        LOGICAL, PARAMETER      :: NEWFFTW     = .TRUE.
        LOGICAL, PARAMETER      :: SPIDER_SIGN = .FALSE.

        integer :: nend

C       FIND PARAMETERS TO NORMALIZE UNDER THE MASK,  
C       TRIED DOING THIS COMPLETELY ON THE
C       POLAR RINGS BUT IT GIVES SOME DIFFERENT REF. CHOICES. al

C       CALCULATE DIMENSIONS FOR NORMALIZING MASK
	NSB  = -CNS2
	NSE  =  NSB + NX - 1
	NRB  = -CNR2
	NRE  =  NRB + NY - 1

        CALL NORMASC(XIM, NSB,NSE,NRB,NRE, NUMR,NUMR(1,NRING),
     &                  AVG,SIGINV,USE_OMP)
        !write(6,*) 'AVG,SIGINV,wr(1):',AVG,SIGINV,wr(1)

C       INTERPOLATE INTO POLAR COORDINATES & APPLY NORMALIZATION
C       CREATING CIRC (RADIAL IMAGE CIRCLES) FOR THIS IMAGE POSITION
        IF (USE_FBS_INTERP) THEN

           !write(6,*)' Using fbs interp on image:',NX,' x ',NY

 	   CALL ALRQ_MS_FBS(XIM,NX,NY, CNS2,CNR2,
     &               NUMR,CIRC,LCIRC, NRING,MODE, NEWFFTW,USE_OMP, 
     &               AVG,SIGINV)
        ELSE
           !write(6,*)' Using quadratic interp on image',NX,' x ',NY

  	   CALL ALRQ_MS_QUAD(XIM,NX,NY, CNS2,CNR2,
     &               NUMR,CIRC,LCIRC, NRING,MODE, NEWFFTW,USE_OMP, 
     &               AVG,SIGINV)
        ENDIF

        !nend = numr(3,1) - 2     ! Length of ring
        !call chkreal('real ring 1', circ,lcirc, nend,1, 1)

 
        WEIGHTIT = (WR(1) > 0.0)

        IF (WR(1) < 0.0) THEN
C          RETURN SUM SQ WEIGHT FOR EACH CIRC RING

           MAXRING         = NUMR(1,NRING) 
           WR(1:MAXRING+1) = 0.0

           DO I=1,NRING

	      IRING = NUMR(1,I)      ! RING NUMBER
	      IGOM1 = NUMR(2,I) - 1  ! START LOC. WITHIN CIRC -1
              NRAYS = NUMR(3,I) - 2  ! NUMBER OF RAYS = POINTS ON RING   

C             NORMALIZATION DETERMINATION
              DO J = 1,NRAYS
                 WR(IRING) = WR(IRING) + CIRC(IGOM1 + J)**2
              ENDDO

              WR(MAXRING + 1) = WR(MAXRING + 1) + WR(IRING)
 
	   ENDDO
        ENDIF

C       FOURIER TRANSFORM CIRC RINGS
        CALL FRNGS_NEWT(CIRC,LCIRC, NUMR,NRING, SPIDER_SIGN,
     &                  FFTW_PLANS, USE_OMP)

        IF (WEIGHTIT) THEN
C          WEIGHT TRANSFORMED CIRC RINGS  USING  FACTORS FROM: WR
           CALL APPLYWS_NEW(CIRC,LCIRC,NUMR,WR,NRING)
        ENDIF

        IRTFLG = 0
        END


C       ---------------  APRINGS_INIT_PLANS -----------------------------

        SUBROUTINE APRINGS_INIT_PLANS(NUMR,NRING,
     &                                FFTW_PLANS,NPLANS,
     &                                NX,NY,IRTFLG)

C       INITIALIZE REV. FFTW3 PLAN(S) FOR USE WITHIN OMP ||

	INCLUDE 'CMBLOCK.INC' 

        INTEGER, INTENT(IN)    :: NUMR(3,NRING)
        INTEGER, INTENT(IN)    :: NRING
        INTEGER, INTENT(IN)    :: NPLANS
        INTEGER, INTENT(IN)    :: NX,NY
        INTEGER, INTENT(OUT)   :: IRTFLG

C       FFTW_PLANS IS AN ARRAY OF POINTERS TO STRUCTURES 
        INTEGER*8, INTENT(OUT) :: FFTW_PLANS(NPLANS)

        FFTW_PLANS = 0
        IRTFLG     = 1
        NUMTH      = 1

        LENO = -1
        DO I=1,NRING
           LEN = NUMR(3,I) - 2     ! LENGTH OF RING

           IF (LEN > LENO) THEN
C             CREATE PLAN FOR THIS NEW RING LENGTH

              INDX = LOG2(LEN) - 1

              IF (INDX < 2 ) THEN
                 write(NOUT,'(A,3I6,I14)')' plan i,len,indx:',i,len,indx
                 CALL ERRT(102,'APRINGS_INIT_PLANS; BAD INDX',INDX)
                 RETURN
              ELSEIF (INDX > NPLANS) THEN
                 write(NOUT,'(A,3I6,I14)')' PLAN I,LEN,INDX:',i,len,indx
                 WRITE(NOUT,*) ' CREATING FFTW_PLANS(',INDX,')'
                 CALL ERRT(102,'APRINGS_INIT_PLANS; OVERFLOW',NPLANS)
                 RETURN
              ENDIF

 	      CALL FFTW3_MAKEPLAN(LEN,1,1, NUMTH,
     &                            FFTW_PLANS(INDX),+1,IRTFLG)

              IF (IRTFLG .NE. 0) RETURN
              !write(6,'(a,3i6,i14)')'plan i,len,indx:',i,len,indx,fftw_plans(indx)

              LENO = LEN
           ENDIF
        ENDDO

C       INITIALIZE REV. FFTW3 PLAN FOR USE WITHIN OMP ||
C       CAN USE INDX = 1 FOR REVERSE PLAN AS IT IS NEVER NEEDED ELSEWHERE.

        MAXRIN = NUMR(3,NRING) - 2
 	CALL FFTW3_MAKEPLAN(MAXRIN,1,1, NUMTH,
     &                      FFTW_PLANS(1),-1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (USE_FBS_INTERP .AND. NX > 0 .AND. NY > 0) THEN
C          CREATE CACHED PLANS FOR FBS ALSO
           CALL FMRS_PLAN(.FALSE.,FDUM,NX,NY,1, NUMTH,+1,IRTFLG)
           CALL FMRS_PLAN(.FALSE.,FDUM,NX,NY,1, NUMTH,-1,IRTFLG)
        ENDIF

        IRTFLG = 0

        END




C **********************************************************************
C
C  APRINGS_SATU 
C  
C  PURPOSE: CONVERT IMAGE--> POLAR RINGS .

C  SOME PARAMETERS:
C       IMGBUF              IMAGE SIZE                        (INPUT)
C       NX,NY           ACTUAL (NOT-WINDOWED) IMAGE SIZE  (INPUT)
C       COEFFS              COEF.                             (INPUT)
C       IXY                 COEF.  LOCS.                      (INPUT)
C       CIRCT               SCRATCH AREA                      (OUTPUT)
C       CIRCIMG             IMAGE IN POLAR RINGS              (OUTPUT)
C       IRTFLG              ERROR FLAG                        (OUTPUT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

cpgi$g opt=3

        SUBROUTINE APRINGS_SATU(IMGBUF,NX,NY,CNS2,CNR2, 
     &                          MODE,NUMR,NRING,LCIRC,USE_OMP,
     &                          WR,FFTW_PLANS, NLOCS,NRAYSC,
     &                          COEFFS,IXY, USECOEF,TRANS,CPLX,
     &                          CIRCT,  CIRCIMG, IRTFLG)

C       NOTE: RUNS WITHIN OMP PARALLEL SECTION OF CODE!

        IMPLICIT NONE

	REAL                   :: IMGBUF(NX,NY)
        INTEGER                :: NX,NY
        REAL                   :: CNS2,CNR2
        CHARACTER (LEN=1)      :: MODE
        INTEGER                :: NRING,LCIRC
        INTEGER                :: NUMR(3,NRING) 
        LOGICAL                :: USE_OMP
        REAL                   :: WR(1)          ! WEIGHTS
        INTEGER *8             :: FFTW_PLANS(*)  ! STRUCTURE POINTERS
        INTEGER                :: NLOCS(2,*)     ! IF TRANS ONLY
        INTEGER                :: NRAYSC
        REAL                   :: COEFFS(6,LCIRC)! IF USECOEF ONLY
        INTEGER                :: IXY(2,LCIRC)   ! IF USECOEF ONLY

        LOGICAL                :: USECOEF        ! CREATE COEF
        LOGICAL                :: TRANS          ! TRANSFORMED RAYS
        LOGICAL                :: CPLX           ! COMPLEX CROSRNG
        REAL                   :: CIRCT(LCIRC)   ! WORK AREA IF TRANS

        REAL,   INTENT(OUT)    :: CIRCIMG(LCIRC)
        INTEGER,INTENT(OUT)    :: IRTFLG

        INTEGER                :: MAXRIN

        IRTFLG = 0
        MAXRIN = NUMR(3,NRING) - 2 ! ACTUAL LENGTH OF LONGEST RING

C       CONVERT IMAGE TO POLAR RINGS ----------------

        IF (TRANS) THEN
C          NOT USING COEFF ARRAY,  TRANSFORMED RINGS

           CALL APRINGS_TRANS_ONE(IMGBUF,NX,NY,  CNS2,CNR2,
     &                            NUMR,NRING, NLOCS,NRAYSC,
     &                            MODE,USE_OMP, WR,FFTW_PLANS,
     &                            CIRCT,LCIRC, CIRCIMG,LCIRC/2)

        ELSEIF (.NOT. USECOEF) THEN
C          NOT USING COEFF ARRAY, NON-TRANSFORMED RINGS

           CALL APRINGS_ONE_NEW(NX,NY, CNS2,CNR2, 
     &                          IMGBUF,USE_OMP,
     &                          MODE,NUMR,NRING,LCIRC,WR,FFTW_PLANS,
     &                          CIRCIMG,IRTFLG)

        ELSEIF (USECOEF) THEN
C          USING COEFF ARRAY, EITHER TRANSFORMED RINGS OR NOT

           CALL APRINGS_ONE_COEF(IMGBUF, NX,NY, CNS2,CNR2, 
     &                           NUMR,NRING, NLOCS,NRAYSC,
     &                           COEFFS,IXY,
     &                           USE_OMP, WR, FFTW_PLANS, TRANS,
     &                           CIRCT,LCIRC, CIRCIMG, IRTFLG)
        ENDIF

        END














 




