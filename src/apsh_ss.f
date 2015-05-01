C++*********************************************************************
C
C    APSH_SS.F      USED CMLIMIT                  AUG 00 ARDEAN LEITH
C                   ADDED REF_CIRC FILE           APR 01 ARDEAN LEITH
C                   NORMASS -> NORMAS             OCT 01 ARDEAN LEITH
C                   PROMPTS                       JAN 02 ARDEAN LEITH
C                   OPFILEC                       FEB 03 ARDEAN LEITH
C                   APMASTER REWRITE              AUG 03 ARDEAN LEITH
C                   'OR MQ' SUPPORT               SEP 03 ARDEAN LEITH
C                   AP_OUT USAGE                  FEB 04 ARDEAN LEITH
C                   'MQ' NOT READ EXP. ANGLES     APR 04 ARDEAN LEITH
C                   BAD PEAK IF INTERP ON EDGE    AUG 04 ARDEAN LEITH
C                   CROSRNG_E SPEEDS UP           AUG 04 ARDEAN LEITH
C                   AP_END CALL HAS PARLIST       OCT 04 ARDEAN LEITH
C                   LIMITRANGE BUG                OCT 04 ARDEAN LEITH
C                   PEAKV = 1                     JAN 05 ARDEAN LEITH
C                   DISCARD MIRROR ...            JUN 06 ARDEAN LEITH
C                   AP_STAT CALL                  JAN 07 ARDEAN LEITH
C                   USE FFTW3 IN APRINGS          MAR 08 ARDEAN LEITH
C                   APRINGS_ONE_NEW               APR 08 ARDEAN LEITH
C                   AP_END CALL ALTERED           NOV 08 ARDEAN LEITH
C                   AP_STAT_ADD                   NOV 08 ARDEAN LEITH
C                   CIRCREF NOT INCORE            AUG 09 ARDEAN LEITH
C                   ISHRANGEX                     FEB 10 ARDEAN LEITH
C                   ANGREF NOT IN OUTPUT          MAY 10 ARDEAN LEITH
C                   CROSRNG NO TT                 JUN 10 ARDEAN LEITH
C                   CROSRNG_2, REMOVED IA_64      JUN 10 ARDEAN LEITH
C                   AP_STAT NBORDER               OCT 10 ARDEAN LEITH
C                   MAKE_CLOSE_LIST, GETANGAS     FEB 11 ARDEAN LEITH
C                   AP_GETDATA USED               NOV 11 ARDEAN LEITH  
C                   ROTFIRST                      DEC 11 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
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
C  APSH_SS(IREFLIST,NUMREF,IEXPLIST,NUMEXP, 
C           NX,NY,NR,ISHRANGEX,ISHRANGEY,ISTEP,
C           NRING,LCIRC,NUMR,CIRCREF,CIRCREF_IN_CORE,
C           MODE,REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS,
C           REFPAT,EXPPAT,RANGE,ROTFIRST,CKMIRROR,
C           CTYPE,LUNDOC,FBS_WANTED)
C
C  PURPOSE: FIND ROTATIONAL AND SHIFT PARAMETERS TO ALIGN A SERIES OF
C           REFERENCE IMAGES WITH SERIES OF SAMPLE IMAGES
C
C           VERSION FOR MP AND A SMALL NUMBER OF IMAGES TO BE ALIGNED. 
C           NOT FOR MPI USE.
C
C  OPERATIONS:  'AP SH', 'AP OR'
C
C PARAMETERS:
C       IREFLIST            LIST OF REF. IMAGE FILE NUMBERS   (INPUT)
C       NUMREF              NO. OF REF. IMAGES                (INPUT)
C       IEXPLIST            LIST OF EXP. IMAGE FILE NUMBERS   (INPUT)
C       NUMEXP              NO. OF EXP IMAGES                 (INPUT)
C       REFANGDOC           REF. ANGLES FILE NAME             (INPUT)
C       EXPANGDOC           EXP. ANGLES FILE NAME             (INPUT)
C       REFPAT              REF. IMAGE SERIES FILE TEMPLATE   (INPUT)
C       EXPPAT              EXP. IMAGE SERIES FILE TEMPLATE   (INPUT)
C       ROTFIRST            USE RTSQ ON EXP INPUT IMAGES      (INPUT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

       SUBROUTINE APSH_SS(IREFLIST,NUMREF,IEXPLIST,NUMEXP, 
     &               NX,NY,ISHRANGEX,ISHRANGEY,ISTEP,
     &               NRING,LCIRC,NUMR,CIRCREF,CIRCREF_IN_CORE,
     &               MODE,REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS,
     &               REFPAT,EXPPAT,RANGE,ROTFIRST,
     &               CKMIRROR,CTYPE,LUNDOC,FBS_WANTED)

        INCLUDE 'MAKE_CLOSE_LIST.INC'  
	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

	INTEGER               :: IREFLIST(NUMREF)
	INTEGER               :: IEXPLIST(NUMEXP)
	INTEGER               :: NUMREF,NUMEXP,NX,NY,ISHRANGEX,ISHRANGEY
	INTEGER               :: ISTEP,NRING,LCIRC
	INTEGER               :: NUMR(3,NRING)
	REAL                  :: CIRCREF(LCIRC,NUMREF)
	LOGICAL               :: CIRCREF_IN_CORE
	LOGICAL               :: MIRRORNEW
	LOGICAL               :: GOTREFANG
	LOGICAL               :: LIMITRANGE
	CHARACTER (LEN=1)     :: MODE,NULL
        CHARACTER (LEN=*)     :: REFANGDOC,EXPANGDOC
        CHARACTER (LEN=*)     :: SCRFILE
        CHARACTER (LEN=*)     :: REFPAT,EXPPAT
        CHARACTER (LEN=*)     :: CTYPE
	LOGICAL               :: ROTFIRST,CKMIRROR
        INTEGER               :: LUNDOC
        LOGICAL               :: FBS_WANTED

	DOUBLE PRECISION      :: FITP(-1:1,-1:1)

	DOUBLE PRECISION      :: CCROTD,PEAK,CCROTD_INTERP
        INTEGER *8            :: FFTW_PLANS(*)

C       AUTOMATIC ARRAYS
	DOUBLE PRECISION      :: CCOA(NUMREF,   -ISHRANGEX:ISHRANGEX,
     &                                          -ISHRANGEY:ISHRANGEY)
	REAL                  :: RANGOA(NUMREF, -ISHRANGEX:ISHRANGEX,
     &                                          -ISHRANGEY:ISHRANGEY)
	LOGICAL               :: ISMIRRORED(NUMREF, 
     &                                          -ISHRANGEX:ISHRANGEX,
     &                                          -ISHRANGEY:ISHRANGEY)

	DOUBLE PRECISION      :: FIT(     -ISTEP:ISTEP,-ISTEP:ISTEP)
	DIMENSION             :: ROTMP(   -ISTEP:ISTEP,-ISTEP:ISTEP)
        LOGICAL               :: ISMIRDUM(-ISTEP:ISTEP,-ISTEP:ISTEP)

        INTEGER, POINTER      :: LCG(:)
	REAL                  :: ANGOUT(3)
	REAL                  :: EXPDIR(3)
        LOGICAL               :: USE_UN,USE_MIR
        LOGICAL               :: ANGINHEADER  

C       ALLOCATED ARRAYS
	REAL, ALLOCATABLE     :: CIRCEXP(:,:,:)
	REAL, ALLOCATABLE     :: A(:,:)
	REAL, ALLOCATABLE     :: REFDIR(:,:) 
	REAL, ALLOCATABLE     :: ANGREF(:,:),ANGEXP(:,:)
	REAL, ALLOCATABLE     :: TMPBUF(:,:)
 
        LOGICAL, PARAMETER    :: MPIBCAST = .TRUE.
        LOGICAL, PARAMETER    :: USE_OMP  = .FALSE.
        INTEGER, PARAMETER    :: NLISTMAX = 15
        REAL                  :: PARLIST(NLISTMAX)

        REAL, PARAMETER       :: QUADPI = 3.1415926535
        REAL, PARAMETER       :: DGR_TO_RAD = (QUADPI/180)

        INTEGER               :: NBORDER = 0       ! # BORDER PIXELS
        INTEGER               :: NSUBPIX = 0       ! # SUBPIX PIXELS

        INTEGER, PARAMETER    :: INPIC   = 77
        INTEGER, PARAMETER    :: INANG   = 78
        INTEGER, PARAMETER    :: LUNRING = 50

        NULL = CHAR(0)

C       INITIALIZE CCROT STATISTICS COUNTERS
        ANGDIFTHR   = 0.0
        CALL  AP_STAT_ADD(-1,CCROT,ANGDIF,ANGDIFTHR,CCROTLAS,
     &                  CCROTAVG,IBIGANGDIF,ANGDIFAVG,IMPROVCCROT,
     &                  CCROTIMPROV,IWORSECCROT,CCROTWORSE)

C       FLAG FOR RESTRICTED PROJECTION RANGE
        LIMITRANGE = RANGE .GT. 0 .AND. RANGE .LT. 360

        MAXRIN = NUMR(3,NRING)
#ifdef SP_LIBFFTW3
        MAXRIN = NUMR(3,NRING) -2  ! ACTUAL LENGTH OF LONGEST RING
#endif

	RANGECOS = COS(RANGE*DGR_TO_RAD)
        WR       = 0.0    ! DUMMY VALUE FLAG FOR APRINGS CALL

C       FIND NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)
        CALL FLUSHRESULTS()

C       READ REFERENCE IMAGES INTO REFERENCE RINGS (CIRCREF) ARRAY OR
C       CREATE REFERENCE RINGS FILE FOR LATER READING 

        CALL APRINGS_NEW(IREFLIST,NUMREF, NX,NY,
     &               NRING,LCIRC,NUMR,MODE,FFTW_PLANS,
     &               REFPAT,INPIC,CIRCREF,CIRCREF_IN_CORE,
     &               LUNRING,SCRFILE,IRTFLG)

        IF (CIRCREF_IN_CORE) THEN
           WRITE(NOUT,91)NUMTH
91         FORMAT('  Ref. rings in core,  Threads: ',I4)
        ELSE
           WRITE(NOUT,92)NUMTH
92         FORMAT('  Ref. rings not in core,  Threads: ',I4)
        ENDIF

	NSISX = MAX(ISHRANGEX/ISTEP, ISTEP)
	NSISY = MAX(ISHRANGEY/ISTEP, ISTEP)
	ALLOCATE(CIRCEXP(LCIRC,-NSISX:NSISX, -NSISY:NSISY),
     &           A(NX,NY), 
     &           STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           NSISTX = 2 * NSISX + 1
           NSISTY = 2 * NSISY + 1
           MWANT = LCIRC*NSISTX*NSISTY + NX*NY 
           CALL  ERRT(46,'APSH_SS; CIRCEXP,....',MWANT)
           GOTO 9999
        ENDIF 

        NULLIFY(LCG)             ! INTEL COMPILER REQUIRES THIS
        NUMREFLCG = NUMREF
        IEND      = NUMREF
        NGOTPAR   = 0
        GOTREFANG = .FALSE.

        IF (LIMITRANGE .OR. CTYPE(1:2) .EQ. 'SH') THEN
C          REFANGLES FILE FOR RESTRICTED ANGULAR SEARCH  OR 'SH'
	   ALLOCATE(ANGREF(3,NUMREF), 
     &              REFDIR(3,NUMREF),STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
               MWANT = 6*NUMREF  
               CALL ERRT(46,'APSH_SS; ANGREF..',MWANT)
               RETURN
           ENDIF 

C          READ REF. ANGLES INTO ANGREF FROM REFANGDOC OR HEADER
C          CONVERT REF. ANGLES TO UNITARY DIRECTIONAL VECTORS (REFDIR).
	   CALL AP_GETANGAS(IREFLIST,NUMREF,0,REFANGDOC,REFPAT,
     &                     INPIC,INANG,3,ANGREF,GOTREFANG,NGOTREF,
     &                     .TRUE.,REFDIR,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
        ENDIF
   
        IF (EXPANGDOC .NE. NULL) THEN
	   ALLOCATE(ANGEXP(8,NUMEXP), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'APSH_SS; ANGEXP',8*NUMEXP)
              RETURN
           ENDIF 

C          LOAD EXP. ANGLES & ALIGNMENT PARAMETERS (ANGEXP) 
C          FROM DOC. FILE (EXPANGDOC) OR IMAGE FILE (REFPAT) HEAD
C          THIS RETURNS NGOTPAR
	   CALL AP_GETANGAS(IEXPLIST,NUMEXP,0,EXPANGDOC,EXPPAT,
     &                     INPIC,INANG,8,ANGEXP,GOTEXPANG,NGOTPAR,
     &                     .FALSE.,EXPDIR,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
        ENDIF

        ANGINHEADER = .FALSE. ! unfinished !!!!!!!!!!!
        IF (ROTFIRST) THEN
	   ALLOCATE(TMPBUF(NX,NY), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'APSH_SS; TMPBUF',NX*NY)
              GOTO 9999
           ENDIF
           anginheader = (expangdoc == '*') ! unfinished !!!!!!!!!!!
           IF (FBS_WANTED) THEN
            WRITE(NOUT,*)' ALIGNING INPUT IMAGES WITH FBS INTERPOLATION'
           ELSE
           WRITE(NOUT,*)' ALIGNING INPUT IMAGES WITH QUAD INTERPOLATION'
           ENDIF
        ENDIF 

C       LOOP OVER EXP. IMAGES TO BE ALIGNED
	DO IEXP=1,NUMEXP
           IMGEXP = IEXPLIST(IEXP)

           IF (LIMITRANGE .AND. EXPANGDOC .NE. NULL) THEN
C             CONVERT EXP. ANGLE TO UNITARY DIRECTIONAL VECTORS (EXPDIR).
	      CALL AP_GETSATA(ANGEXP(1,IEXP),EXPDIR,8,1,IRTFLG)

C             MAKE LIST OF NEARBY REFERENCE IMAGES, RETURNS: NUMREFLCG
              CALL MAKE_CLOSE_LIST(NUMREF,LIMITRANGE,
     &                             REFDIR,EXPDIR,
     &                             RANGECOS, CKMIRROR, 
     &                             LCG, NUMREFLCG, IRTFLG)
              !print *,'numreflcg:',numreflcg
              !print *,'list:',lcg(1:numreflcg)
              IEND = NUMREFLCG

              IF (NUMREFLCG .LE. 0) THEN
C                REPORT THAT THERE ARE NO NEARBY REFERENCE IMAGES
                 IMGREF    = 0
                 PEAKV     = 0.0
                 CCROT     = 0.0
                 XSHNEW    = 0.0
                 YSHNEW    = 0.0
                 MIRRORNEW = .FALSE.

                 CALL AP_END(IEXP,IMGEXP,IMGREF,
     &                ANGREF,REFDIR,
     &                ANGEXP(1,IEXP),EXPDIR,ISHRANGEX,
     &                GOTREFANG, NGOTPAR, CCROT,PEAKV,
     &                RANGNEW,XSHNEW,YSHNEW,MIRRORNEW,REFPAT,
     &                NUMREFLCG, CTYPE, LUNDOC,PARLIST)

                 CYCLE
              ENDIF
           ENDIF

C          LOAD EXP. IMAGE DATA FOR THIS IMAGE
           IMITT = IEXP
           IF (CTYPE(1:2) == 'OR') IMITT = 0

           IF (ROTFIRST) THEN
C             WANT TO ROTATE/SHIFT EXP IMAGE WHEN READING THEM
	      CALL AP_GETDATA_RTSQ(IEXPLIST,NUMEXP, 
     &                    NX,NY, NX,NY,0.0,
     &                    NUMTH,EXPPAT,INPIC, IMITT,IMITT,
     &                    ANGINHEADER, ANGEXP(1,IEXP), 
     &                    MPIBCAST, TMPBUF, A,
     &                    .TRUE., AVI,SIGI, FBS_WANTED,IRTFLG)
           ELSE
	      CALL AP_GETDATA(IEXPLIST,NUMEXP,
     &                     NX,NY, NX,NY, 0.0,
     &                     NUMTH,EXPPAT,INPIC,IMITT,IMITT,
     &                     MPIBCAST,A, 
     &                     .TRUE.,AVI,SIGI, IRTFLG)
           ENDIF

           IF (IRTFLG .NE. 0) GOTO 9999

           IF ( CIRCREF_IN_CORE) THEN
C             USE CIRCREF ARRAY FOR REFERENCE RINGS
              !write(6,*) ' incore, old parallel '

C             LOOP OVER ALL SHIFTS IN ||
c$omp         parallel do private(jt,cnr2,it,cns2,imil,imi,use_un,
c$omp&              use_mir)
             
              DO JT=-ISHRANGEY,ISHRANGEY,ISTEP
C                LOOP OVER SHIFTED CENTERS IN Y
                 CNR2 = NY/2+1+JT

                 DO IT=-ISHRANGEX,ISHRANGEX,ISTEP
C                   LOOP OVER SHIFTED CENTERS IN X
                    CNS2 = NX/2+1+IT

C                   NORMALIZE IMAGE VALUES UNDER THE MASK OVER VARIANCE 
C                   RANGE INTERPOLATE TO POLAR COORDINATES, CREATE 
C                   FOURIER OF: CIRCEXP.  NO WEIGHTING OF RINGS
	            CALL APRINGS_ONE_NEW(NX,NY, CNS2,CNR2, 
     &                           A,.FALSE.,
     &                           MODE,NUMR,NRING,LCIRC,WR,FFTW_PLANS,
     &                           CIRCEXP(1,IT/ISTEP,JT/ISTEP),IRTFLG)

                    DO IMIL=1,IEND
C                      LOOP OVER REFERENCE IMAGES
                       IMI = IMIL
                       IF (LIMITRANGE) IMI = ABS(LCG(IMIL))

                       IF (CKMIRROR .AND. LIMITRANGE) THEN
C                         ONLY SEARCH EITHER MIRRORED OR NON-MIRRORED
                          USE_UN  = (LCG(IMIL) .GE. 0)
                          USE_MIR = (LCG(IMIL) .LT. 0)
                       ELSE
C                         SEARCH BOTH MIRRORED & NON-MIRRORED IF CHKMIR
                          USE_UN  = .TRUE.
                          USE_MIR = CKMIRROR
                       ENDIF

C                      CHECK EITHER MIRRORED/NON-MIRRORED POSITIONS 
                       CALL CROSRNG_2(CIRCREF(1,IMI),
     &                            CIRCEXP(1,IT/ISTEP,JT/ISTEP),
     &                            LCIRC,NRING, MAXRIN,NUMR,
     &                            USE_OMP,FFTW_PLANS(1),
     &                            USE_UN,USE_MIR,   
     &                            ISMIRRORED(IMI,IT,JT),
     &                            CCOA(IMI,IT,JT),RANGOA(IMI,IT,JT))
	            ENDDO  ! END OF: DO IMIL=1,IEND
                 ENDDO     ! END OF: DO IT=-ISHRANGEX,ISHRANGEX,ISTEP
              ENDDO        ! END OF: DO JT=-ISHRANGEY,ISHRANGEY,ISTEP

C$OMP         END PARALLEL DO
C             END OF THE OMP PARALLEL SECTION
           ELSE
              !write(6,*) ' not incore,  new parallel  '
C             USE REFERENCE RINGS FILE (MIGHT BE AN INCORE FILE?) -----

              DO JT=-ISHRANGEY,ISHRANGEY,ISTEP ! LOOP OVER ALL SHIFTS IN ||
C                LOOP OVER SHIFTED CENTERS IN Y
                 CNR2 = NY/2+1+JT

                 DO IT=-ISHRANGEX,ISHRANGEX,ISTEP
C                   LOOP OVER SHIFTED CENTERS IN X
                    CNS2 = NX/2+1+IT

C                   NORMALIZE IMAGE VALUES UNDER THE MASK OVER VARIANCE 
C                   RANGE INTERPOLATE TO POLAR COORDINATES, CREATE 
C                   FOURIER OF: CIRCEXP.  NO WEIGHTING OF RINGS
	            CALL APRINGS_ONE_NEW(NX,NY,  CNS2,CNR2, 
     &                           A,.FALSE.,
     &                           MODE,NUMR,NRING,LCIRC,WR,FFTW_PLANS,
     &                           CIRCEXP(1,IT/ISTEP,JT/ISTEP),IRTFLG)

c$omp               parallel do private(imil,imi,ithread,use_un,use_mir)
c$omp&              schedule(static,1)
                    DO IMIL=1,IEND
C                      LOOP OVER ALL REFERENCE IMAGES IN ||

                       IMI = IMIL
                       IF (LIMITRANGE) IMI = ABS(LCG(IMIL))
C                      FIND THREAD NUMBER 
                       ITHREAD = MOD((IMIL-1),NUMTH) + 1

C                      FILL CIRCREF FROM REFERENCE RINGS FILE
c$omp                  critical
                       CALL REDLIN(LUNRING,CIRCREF(1,ITHREAD),LCIRC,IMI)
c$omp                  end critical

                       IF (CKMIRROR .AND. LIMITRANGE) THEN
C                         ONLY SEARCH EITHER MIRRORED OR NON-MIRRORED
                          USE_UN  = (LCG(IMIL) .GE. 0)
                          USE_MIR = (LCG(IMIL) .LT. 0)
                       ELSE
C                         SEARCH BOTH MIRRORED & NON-MIRRORED IF CHKMIR
                          USE_UN  = .TRUE.
                          USE_MIR = CKMIRROR
                       ENDIF

C                      CHECK EITHER MIRRORED/NON-MIRRORED POSITIONS 
                       CALL CROSRNG_2(CIRCREF(1,ITHREAD),
     &                          CIRCEXP(1,IT/ISTEP,JT/ISTEP),
     &                          LCIRC,NRING, MAXRIN,NUMR,
     &                          USE_OMP,FFTW_PLANS(1),
     &                          USE_UN,USE_MIR,  ISMIRRORED(IMI,IT,JT),
     &                          CCOA(IMI,IT,JT), RANGOA(IMI,IT,JT))
	            ENDDO  ! END OF: DO IMIL=1,IEND
                 ENDDO     ! END OF: DO IT=-ISHRANGEX,ISHRANGEY,ISTEP
              ENDDO        ! END OF: DO JT=-ISHRANGEX,ISHRANGEY,ISTEP
           ENDIF           ! END OF: ELSE/IF ( CIRCREF_IN_CORE) THEN


C          LOCATE BEST CCROT MATCH FROM ALL THE VALUES ACCUMULATED ABOVE
           CCROTD  = -1.0D23

           DO JT=-ISHRANGEY,ISHRANGEY,ISTEP
              DO IT=-ISHRANGEX,ISHRANGEX,ISTEP

                 DO IRR=1,IEND
C                   LOOP OVER REFERENCE IMAGES
                    IR = IRR
                    IF (LIMITRANGE) IR = ABS(LCG(IRR))

                    IF (CCOA(IR,IT,JT) .GE. CCROTD)  THEN
C                      BETTER POSITION    
	               CCROTD    = CCOA(IR,IT,JT)
                       IREF      = IR
                       ISX       = IT
                       ISY       = JT
                       RANGNEW   = ANG_N(RANGOA(IR,IT,JT),MODE,MAXRIN)
                       MIRRORNEW = ISMIRRORED(IR,IT,JT)
                    ENDIF
                 ENDDO     ! END OF: DO IRR=1,IEND
               ENDDO        ! END OF: DO IT=-ISHRANGEX,ISHRANGEX,ISTEP
            ENDDO           ! END OF: DO JT=-ISHRANGEY,ISHRANGEY,ISTEP
C          write(6,*) 'iref: ',iref,  '   ccrotd: ',ccrotd

	   SX      = ISX
           SY      = ISY
	   CCROT   = CCROTD
           IMGREF  = IREFLIST(IREF)

           IF (CIRCREF_IN_CORE) THEN
              IREFT = IREF
           ELSE
C             FILL CIRCREF (1,1) FROM REFERENCE RINGS FILE
              IREFT = 1
              CALL REDLIN(LUNRING,CIRCREF(1,IREFT),LCIRC,IREF)
           ENDIF

#ifdef DEBUG
           write(6,921) imgref,isx,isy,ccrotd,rangnew
921        format(' 1 ',i5,' (',i3,',',i3,'): ',f14.4,' ',2f8.2,f6.1)
#endif

C          CHECK LOCATIONS WITHIN ISHRANGE AROUND MAX  ------------------
   
	   IF (ABS(ISX) .NE. ISHRANGEX .AND. 
     &         ABS(ISY) .NE. ISHRANGEY)  THEN
C             NOT ON BOUNDARY, HAVE TO FIND NEIGHBOURING VALUES

	      FIT(0,0)   = CCROTD
	      ROTMP(0,0) = RANGNEW

c$omp         parallel do private(jt,it,cnr2,cns2)
	      DO JT=-ISTEP,ISTEP
	         DO IT=-ISTEP,ISTEP
	            CNR2 = NY / 2 + 1 + JT + ISY
	            IF (IT.NE.0 .OR. JT.NE.0) THEN
	               CNS2 = NX / 2 + 1 + IT + ISX

C                      NORMALIZE IMAGE VALUES UNDER THE MASK OVER VARIANCE RANGE
C                      INTERPOLATE INTO POLAR COORDINATES
C                      CREATES FOURIER OF: CIRCEXP

	               CALL APRINGS_ONE_NEW(NX,NY, CNS2,CNR2, 
     &                           A,.FALSE.,
     &                           MODE,NUMR,NRING,LCIRC,WR,FFTW_PLANS,
     &                           CIRCEXP(1,IT,JT),IRTFLG)

C                      ONLY CHECK THE MIR/UN-MIR RETURNED ABOVE
                       CALL CROSRNG_2(CIRCREF(1,IREFT),
     &                          CIRCEXP(1,IT,JT),
     &                          LCIRC,NRING, MAXRIN,NUMR,
     &                          USE_OMP,FFTW_PLANS(1),
     &                          .NOT. MIRRORNEW,MIRRORNEW,  
     &                          ISMIRDUM(IT,JT),
     &                          FIT(IT,JT),ROTMP(IT,JT))

C                      RECORD BEST ANGLE IN ROTMP
	               ROTMP(IT,JT) = ANG_N(ROTMP(IT,JT),MODE,MAXRIN)
	            ENDIF
	         ENDDO      ! END OF: DO IT=-ISTEP,ISTEP
	      ENDDO         ! END OF: DO JT=-ISTEP,ISTEP
c$omp         end parallel do 
	         
C             FIND THE MAXIMUM CC ANGLE WITHIN +/-ISTEP
C             MAXIMUM CANNOT BE ON THE EDGE, I.E., IT,JT/=ISTEP

	      AFIT     = FIT(0,0)
	      JTMA     = 0
	      ITMA     = 0
              RANGNEWT = ROTMP(0,0)

	      IF (ISTEP .GT. 1) THEN
	         DO JT=-ISTEP+1,ISTEP-1
	            DO IT=-ISTEP+1,ISTEP-1
	               IF (FIT(IT,JT) .GT. AFIT)  THEN
	                  AFIT     = FIT(IT,JT)
	                  RANGNEWT = ROTMP(IT,JT) !compiler bug on OPT64
	                  ITMA     = IT
	                  JTMA     = JT
	               ENDIF
	            ENDDO
	         ENDDO
	      ENDIF    ! END OF: IF (ISTEP .GT. 1)
C             write(6,*) ((fit(i,j),i=-1,1),j=-1,1)

C             TEMP VARIABLE OVERCOMES COMPILER BUG ON OPT 64 PGI 6.0
              RANGNEW = RANGNEWT

C             COPY VALUES AROUND THE PEAK.
	      DO JT=-1,1
	         DO IT=-1,1
	            FITP(IT,JT) = FIT(ITMA+IT,JTMA+JT)
c                   write(6,910) it,jt,fitp(it,jt)
910                 format(' fitp(',i5,',',i5,') : ',f14.4)
	         ENDDO
	      ENDDO

C             UPDATE INTEGRAL LOCATION OF THE PEAK
              IF (AFIT > CCROTD) NBORDER = NBORDER + 1
	      CCROTD  = AFIT
	      ISX     = ISX + ITMA
	      ISY     = ISY + JTMA
              SX      = ISX
              SY      = ISY

#ifdef DEBUG
	      cnr2 = NY / 2 + 1 + isy
              cns2 = NX / 2 + 1 + isx
              write(6,905)imgref,isx,isy,ccrotd,rangnew,cns2,cnr2,sx,sy
905           format(' 2 ',i5,' (',i3,',',i3,'): ',f12.4,' ',5f8.2)
#endif

C             SUB-PIXEL INTERPOLATION ------------------------------

C             FIND PEAK BY FITTING PARABOLA TO 3x3 NEIGHBORHOOD
	      CALL PARABLD(FITP,SSX,SSY,PEAK)

	      IF (ABS(SSX) .LT. 1.0 .AND. ABS(SSY) .LT. 1.0)  THEN
C                INTERPOLATED LOCATION IS NOT ON BOUNDARY

	         CNS2 = NX/2+1 + SX + SSX
	         CNR2 = NY/2+1 + SY + SSY

C                NORMALIZE IMAGE VALUES UNDER MASK OVER VARIANCE RANGE
C                INTERPOLATE INTO POLAR COORD., CREATE FFT OF: CIRCEXP
C                CAN NOT USE: APRINGS_ONE_COEF AS NOT INTEGRAL SHIFT

	         CALL APRINGS_ONE_NEW(NX,NY, CNS2,CNR2, A,.FALSE.,
     &                           MODE,NUMR,NRING,LCIRC,WR,FFTW_PLANS,
     &                           CIRCEXP,IRTFLG)
 
C                ONLY CHECK THE MIR/UN-MIR RETURNED ABOVE
                 CALL CROSRNG_2(CIRCREF(1,IREFT), CIRCEXP,
     &                          LCIRC,NRING, MAXRIN,NUMR,
     &                          USE_OMP,FFTW_PLANS(1),
     &                          .NOT. MIRRORNEW,MIRRORNEW,  
     &                          ISMIRDUM(0,0),
     &                          CCROTD_INTERP,RANGNEW_INTERP)

#ifdef DEBUG
	         rt1 = ang_n(rangnew_interp,mode,maxrin)
                 write(6,904) imgref,isx,isy, ccrotd_interp,rt1,
     &                     cns2,cnr2, sx+ssx,sy+ssy
904              format(' 3 ',i5,' (',i3,',',i3,'): ',f12.4,' ',5f8.2)
#endif

                 IF (CCROTD_INTERP .GT. CCROTD) THEN
C                   USE SUB-PIXEL LOCATION
                    NSUBPIX = NSUBPIX + 1
                    CCROTD  = CCROTD_INTERP
	            RANGNEW = ANG_N(RANGNEW_INTERP,MODE,MAXRIN)
	            SX      = SX + SSX 
	            SY      = SY + SSY 
                 ENDIF
	      ENDIF  ! END OF: IF (ABS(SX) .LT. 1.0 .....
           ENDIF     ! END OF: IF ( ABS(ISX) .NE. ISHRANGE ......
	
           CCROT = CCROTD
	   SX    = -SX
	   SY    = -SY

C          HAVE TO CHANGE ORDER OF SHIFT & ROTATION.
C          IN THIS PROGRAM IMAGE IS SHIFTED FIRST, ROTATED SECOND.
C          IN 'RT SQ' IT IS ROTATION FIRST, SHIFT SECOND.
C          THIS CODE CORRESPONDS TO 'SA P'.
	   CO        =  COS(RANGNEW * DGR_TO_RAD)
	   SO        = -SIN(RANGNEW * DGR_TO_RAD)
	   XSHNEW    = SX*CO - SY*SO
	   YSHNEW    = SX*SO + SY*CO
           PEAKV     = 1.0

           !write(6,*) 'iref, angref:', iref,angref(1:3,iref)

           CALL AP_END(IEXP,IMGEXP,IMGREF,
     &                ANGREF(1,IREF),REFDIR(1,IREF),
     &                ANGEXP(1,IEXP),EXPDIR,ISHRANGEX,
     &                GOTREFANG, NGOTPAR, CCROT,PEAKV,
     &                RANGNEW,XSHNEW,YSHNEW,MIRRORNEW,REFPAT,
     &                NUMREFLCG, CTYPE, LUNDOC,PARLIST)

           CALL AP_STAT_ADD(NGOTPAR,CCROT,PARLIST(10),
     &                     ANGDIFTHR,ANGEXP(8,IEXP),
     &                     CCROTAVG,IBIGANGDIF,ANGDIFAVG,IMPROVCCROT,
     &                     CCROTIMPROV,IWORSECCROT,CCROTWORSE)

#ifdef DEBUG
           WRITE(6,96) 'HOST EXP:', IMGEXP,
     &                 '   REF:',   IMGREF,
     &                 '   XSH:',   ISX,
     &                 '   YSH:',   ISY,
     &                 '   RAY:',   IRAY,
     &                 '   CC:',    CCROTD
 96        FORMAT ('  ',A,I6, A,I5,  A,I3, A,I3, A,I5, A,F11.2)
           cns2 = NX / 2 + 1 - sx
           cnr2 = NY / 2 + 1 - sy
           write(6,906) imgref,isx,isy, ccrotd,rangnew,
     &                  cns2,cnr2, xshsum,yshsum
906        format(' 4 ',i5,' (',i3,',',i3,'): ',f12.4,' ',5f8.2)

           write(6,*) ' ------------------------------------- '
           write(6,*) '  '
#endif
	ENDDO
	
        CALL REG_GET_USED(NSEL_USED)
        IF (NSEL_USED .GT. 0 .AND. CTYPE(1:2) .EQ. 'OR') THEN
C           OUTPUT TO REGISTER NOT TO DOC FILE (FOR 'OR')
            DMR = 0
            IF (MIRRORNEW) DMR = 1
            CALL REG_SET_NSEL(1,5,RANGNEW,XSHNEW,YSHNEW,
     &                        DMR,CCROT,IRTFLG)
        ENDIF

        IF (LUNDOC .GT. 0 .AND. NUMEXP .GT. 1) THEN
C          SAVE CCROT & ANGULAR DISPLACEMENT STATISTICS
           CALL AP_STAT(NUMEXP,ANGDIFTHR,IBIGANGDIF,
     &                  ANGDIFAVG, CCROTAVG,
     &                  IMPROVCCROT,CCROTIMPROV,
     &                  IWORSECCROT,CCROTWORSE,
     &                  NBORDER,NSUBPIX,LUNDOC)
        ENDIF



9999    IF (ALLOCATED(CIRCEXP))  DEALLOCATE(CIRCEXP)
	IF (ALLOCATED(A))        DEALLOCATE(A)
	IF (ALLOCATED(ANGREF))   DEALLOCATE(ANGREF)
	IF (ALLOCATED(ANGEXP))   DEALLOCATE(ANGEXP)
	IF (ALLOCATED(REFDIR))   DEALLOCATE(REFDIR)
        IF (ASSOCIATED(LCG))     DEALLOCATE(LCG)
        NULLIFY(LCG)

	END

