
C++*********************************************************************
C
C APREF_P.F
C                   CROSRNG_E SPEEDS UP           AUG 04 ARDEAN LEITH
C                   AP_END CALL HAS PARLIST       OCT 04 ARDEAN LEITH
C                   SPIDER REF_RINGS FILE         FEB 05 ARDEAN LEITH
C                   CCROT STATISTICS              FEB 05 ARDEAN LEITH
C                   NON-INCORE PARALLEL           FEB 05 ARDEAN LEITH
C                   APRINGS_ONE                   MAR 05 ARDEAN LEITH
C                   FFTW_PLANS                    MAR 08 ARDEAN LEITH
C                   APRINGS_ONE_NEW               MAR 08 ARDEAN LEITH
C                   MPI PARTAB BUG FIXED          OCT 08 ARDEAN LEITH
C                   MPI LOCAL READ NOW            NOV 08 ARDEAN LEITH
C                   AP_STAT_ADD                   NOV 08 ARDEAN LEITH
C                   CROSRNG_2, TT REMOVED         JUN 10 ARDEAN LEITH
C                   AP_STAT NBORDER               OCT 10 ARDEAN LEITH
C                   RENAMED FROM DSGR_P           JAN 11 ARDEAN LEITH
C                   MAKE_CLOSE_LIST, GETANGAS     FEB 11 ARDEAN LEITH
C                   AP_GETDATA USED               NOV 11 ARDEAN LEITH  
C                   FBS_WANTED                    JAN 12 ARDEAN LEITH
C
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
C  APREF_P(IREFLIST,NUMREF,IEXPLIST,NUMEXP,
C         NSAM,NROW,NR,RANGE,
C         NRING,LCIRC,NUMR,CIRCREF,CIRCREF_IN_CORE,
C         MODE,REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS,
C         REFPAT,EXPPAT,CKMIRROR,CTYPE,ISHRANGE,LUNDOC)
C
C  PURPOSE: FIND ROTATIONAL AND SHIFT PARAMETERS TO ALIGN A SERIES OF
C           REFERENCE IMAGES WITH SAMPLE IMAGES
C
C  ORDER OF PROCESSING:
C 1 - CONVERT EACH REFERENCE IMAGE TO RINGS, FFT ALL THE RINGS.
C     HIGHEST FREQUENCY FOR ALL THE RINGS EXCEPT LAST ARE DIVIDED 
C     BY 2., STORE RINGS IN CIRCREF. APPLY WEIGHTS TO THE RINGS, 
C 2 - CONVERT EACH INPUT IMAGE TO RINGS, FFT ALL THE RINGS.
C 3 - COMPARE EACH INPUT IMAGE WITH ALL THE REFERENCE IMAGES
C     USING CROSRNG_**. 
C     SINCE Y WERE PRE-WEIGHTED THE RESULT IS ALREADY CORRECT.
C
C PARAMETERS:
C       IREFLIST            LIST OF REF. IMAGE FILE NUMBERS   (INPUT)
C       NUMREF              NO. OF IMAGES                     (INPUT)
C       IEXPLIST            LIST OF EXP. IMAGE FILE NUMBERS   (INPUT)
C       NUMEXP              NO. OF IMAGES                     (INPUT)
C       REFANGDOC           REF. ANGLES FILE NAME             (INPUT)
C       EXPANGDOC           EXP. ANGLES FILE NAME             (INPUT)
C       REFPAT              REF. IMAGE SERIES FILE TEMPLATE   (INPUT)
C       EXPPAT              EXP. IMAGE SERIES FILE TEMPLATE   (INPUT)
C
C  OPERATIONS:  'AP REF', 'AP RD', 'AP RN'
C
C--*********************************************************************

#ifndef USE_MPI
C ---------------------  NON-MPI CODE --------------------------------

         SUBROUTINE APREF_P(IREFLIST,NUMREF, IEXPLIST,NUMEXP,
     &          NSAM,NROW,RANGE,ANGDIFTHR,
     &          NRING,LCIRC,NUMR,CIRCREF,CIRCREF_IN_CORE,
     &          MODE,REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS,
     &          REFPAT,EXPPAT,CKMIRROR,CTYPE,
     &          ROTFIRST,ISHRANGE,LUNDOC,FBS_WANTED)

        INCLUDE 'MAKE_CLOSE_LIST.INC'  
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

	INTEGER                      :: IREFLIST(NUMREF) 
	INTEGER                      :: NUMREF 
	INTEGER                      :: IEXPLIST(NUMEXP) 
	INTEGER                      :: NUMEXP 
	INTEGER                      :: NSAM,NROW
        REAL                         :: RANGE,ANGDIFTHR 
	INTEGER                      :: NRING,LCIRC 
        INTEGER                      :: NUMR(3,NRING)
	REAL                         :: CIRCREF(LCIRC,NUMREF)
	LOGICAL                      :: CIRCREF_IN_CORE
	CHARACTER (LEN=1)            :: MODE 
        CHARACTER (LEN=*)            :: REFANGDOC
        CHARACTER (LEN=*)            :: EXPANGDOC
        CHARACTER (LEN=*)            :: SCRFILE
        INTEGER *8                   :: FFTW_PLANS(*)
        CHARACTER (LEN=*)            :: REFPAT,EXPPAT 
	LOGICAL                      :: CKMIRROR
        CHARACTER (LEN=*)            :: CTYPE 
        LOGICAL                      :: ROTFIRST,FBS_WANTED
        INTEGER                      :: ISHRANGE,LUNDOC

        CHARACTER (LEN=74)           :: COMMENT 
	DOUBLE PRECISION             :: CCROTD
	CHARACTER (LEN=1)            :: NULL
        LOGICAL                      :: MIRRORNEW
	LOGICAL                      :: GOTREFANG,GOTEXPANG
	LOGICAL                      :: ONLYMIRROR
        LOGICAL                      :: LIMITRANGE
        LOGICAL                      :: MIRRORED
        LOGICAL                      :: USE_UN,USE_MIR
        LOGICAL                      :: ANGINHEADER

C       ALLOCATABLE ARRAYS
	REAL, ALLOCATABLE            :: EXPBUF(:) 
        REAL, ALLOCATABLE            :: CIRCEXP(:)
	REAL, ALLOCATABLE            :: TMPBUF(:,:)

C       AUTOMATIC ARRAYS
 	DOUBLE PRECISION             :: TOTMIN(NUMREF)
 	LOGICAL                      :: ISMIRRORED(NUMREF)
 	REAL                         :: ROTANGT(NUMREF) 
	REAL                         :: DLIST(6)
	REAL                         :: ANGEXP(8,NUMEXP)
	REAL                         :: ANGREF(3,NUMREF)
	REAL                         :: REFDIR(3,NUMREF)
	REAL                         :: EXPDIR(3)
	INTEGER, POINTER             :: LCG(:)
        INTEGER, PARAMETER           :: NLISTMAX = 15
        REAL                         :: PARLIST(NLISTMAX) 

        INTEGER                      :: NSAID = 0

        LOGICAL, PARAMETER           :: USE_OMP = .FALSE.

        REAL, PARAMETER              :: QUADPI     = 3.1415926535
        REAL, PARAMETER              :: DGR_TO_RAD = (QUADPI/180)

        INTEGER                      :: NBORDER = 0       ! UNUSED
        INTEGER                      :: NSUBPIX = 0       ! UNUSED

        INTEGER, PARAMETER           :: INPIC   = 77
        INTEGER, PARAMETER           :: LUNANG  = 78
        INTEGER, PARAMETER           :: LUNRING = 50

        MYPID = -1                 ! NOT USING MPI
        NULL  = CHAR(0)

        MAXRIN = NUMR(3,NRING)     ! ACTUAL LENGTH OF LONGEST RING
#ifdef SP_LIBFFTW3
        MAXRIN = NUMR(3,NRING) - 2
#endif

	RANGECOS  = COS(RANGE*DGR_TO_RAD)

C       FIND NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)

        IF (MODE .EQ. 'H')  THEN
           DIVAS = 180.0
        ELSE
           DIVAS = 360.0
        ENDIF

C       MAKE EXPBUF AT LEAST NSAM*NROW FOR USE BY APSHIFT
        MWANTX = NSAM * NROW
	ALLOCATE(EXPBUF(MWANTX), CIRCEXP(LCIRC),STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'APREF_P; EXPBUF, CIRCEXP',MWANTX+LCIRC)
           GOTO 9999
        ENDIF

        ANGINHEADER = .FALSE. ! unfinished !!!!!!!!!!!
        IF (ROTFIRST) THEN
	   ALLOCATE(TMPBUF(NSAM,NROW), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'APREF_P; TMPBUF',NSAM*NROW)
              GOTO 9999
           ENDIF
           !anginheader = (rotfirst .and. expangdoc .eq. '*') ! unfinished !!!!!!!!!!!
           IF (FBS_WANTED) THEN
            WRITE(NOUT,*)' ALIGNING INPUT IMAGES WITH FBS INTERPOLATION'
           ELSE
           WRITE(NOUT,*)' ALIGNING INPUT IMAGES WITH QUAD INTERPOLATION'
           ENDIF
        ENDIF 

C       LOAD REF. PROJ. ANGLES FROM DOC. FILE  OR IMAGE HEAD (IF WANTED)
C       CONVERT REF. ANGLES TO UNITARY DIRECTIONAL VECTORS (REFDIR).
	CALL AP_GETANGAS(IREFLIST,NUMREF,0,REFANGDOC,REFPAT,
     &                   INPIC,LUNANG,3,ANGREF,GOTREFANG,NGOTREF,
     &                   .TRUE.,REFDIR,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       READ REFERENCE IMAGES INTO REFERENCE RINGS (CIRCREF) ARRAY OR
C       CREATE REFERENCE RINGS FILE FOR LATER READING 
        CALL APRINGS_NEW(IREFLIST,NUMREF, NSAM,NROW,
     &               NRING,LCIRC,NUMR,MODE,FFTW_PLANS, 
     &               REFPAT,INPIC,CIRCREF,CIRCREF_IN_CORE,
     &               LUNRING,SCRFILE,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (NSAID .LE. 0) THEN
           IF (CIRCREF_IN_CORE) THEN
              WRITE(NOUT,91)NUMTH
91            FORMAT('  Ref. rings in core,  Threads: ',I4)
           ELSE
              WRITE(NOUT,92)NUMTH
92            FORMAT('  Ref. rings not in core,  Threads: ',I4)
           ENDIF
           NSAID = NSAID + 1
        ENDIF

C       LOAD EXP. ANGLES & ALIGNMENT PARAM. FROM DOC. FILE OR 
C       HEADER (IF WANTED) 
        CALL AP_GETANGAS(IEXPLIST,NUMEXP,0,EXPANGDOC,EXPPAT,
     &                   INPIC,LUNANG,8,ANGEXP,GOTEXPANG,NGOTPAR,
     &                   .FALSE.,FDUM,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       INITIALIZE CCROT STATISTICS
        CALL  AP_STAT_ADD(-1,CCROT,ANGDIF,ANGDIFTHR,CCROTLAS,
     &                   CCROTAVG,IBIGANGDIF,ANGDIFAVG,IMPROVCCROT,
     &                   CCROTIMPROV,IWORSECCROT,CCROTWORSE)

C       CALCULATE DIMENSIONS FOR APRINGS
        CNS2 = NSAM / 2 + 1
        CNR2 = NROW / 2 + 1

        DO IEXP=1,NUMEXP
C         LOOP OVER ALL EXPERIMENTAL (SAMPLE) IMAGES ------------------
    
C         CONVERT EXP. ANGLE TO UNITARY DIRECTIONAL VECTORS (EXPDIR).
	  CALL AP_GETSATA(ANGEXP(1,IEXP),EXPDIR,8,1,IRTFLG)

C         LOAD CURRENT (SINGLE) EXPERIMENTAL IMAGE INTO ARRAY EXPBUF
          IF (ROTFIRST) THEN
C            WANT TO ROTATE/SHIFT EXP IMAGES WHEN READING THEM
	     CALL AP_GETDATA_RTSQ(IEXPLIST,NUMEXP,
     &                        NSAM,NROW, NSAM,NROW, 0.0,
     &                        1,EXPPAT,INPIC, IEXP,IEXP,
     &                        ANGINHEADER, ANGEXP(1,IEXP), 
     &                        .TRUE., TMPBUF,EXPBUF,
     &                        .TRUE.,AVI,SIGI, FBS_WANTED,IRTFLG)
          ELSE
 	     CALL AP_GETDATA(IEXPLIST,NUMEXP,
     &                       NSAM,NROW, NSAM,NROW, 0.0,
     &                       1,EXPPAT,INPIC, IEXP,IEXP,
     &                       .TRUE., EXPBUF,
     &                       .TRUE.,AVI,SIGI, IRTFLG)
          ENDIF
          IF (IRTFLG .NE. 0) GOTO 9999

C         EXTRACT EXP. IMAGE POLAR COORD. RINGS, NORMALIZE & FFT THEM
	  CALL APRINGS_ONE_NEW(NSAM,NROW,  CNS2,CNR2, EXPBUF,.FALSE.,
     &                         MODE,NUMR,NRING,LCIRC, 0.0,FFTW_PLANS,
     &                         CIRCEXP,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9999

C         DETERMINE WHICH REF IMAGES ARE TO BE COMPARED
          IEND       = NUMREF
          CCROTD     = -1.0D20
          LIMITRANGE = (RANGECOS < 1.0)
          NULLIFY(LCG)  ! FOR INTEL COMPILER

C         IF LIMITRANGE, LIST NEARBY REF IMAGES, RETURNS: NUMREFLCG
          CALL MAKE_CLOSE_LIST(NUMREF,LIMITRANGE,
     &                         REFDIR,EXPDIR,
     &                         RANGECOS, .TRUE., 
     &                         LCG, NUMREFLCG, IRTFLG)

          IF (NUMREFLCG <= 0) THEN
C            NO REF. IMAGE WITHIN COMPARISON ANGLE
             IREF = 0 
             GOTO 1000
          ENDIF

          IEND = NUMREFLCG

C         REF. IMAGES FOUND WITHIN COMPARISION RANGE

          IF (CIRCREF_IN_CORE) THEN
C            USE CIRCREF FOR REFERENCE RINGS

c$omp        parallel do private(imil,imi,use_un,use_mir)
             DO IMIL=1,IEND
                IMI = IMIL
	        IF (LIMITRANGE) IMI = ABS(LCG(IMIL))

                IF (CKMIRROR .AND. LIMITRANGE) THEN
C                   ONLY SEARCH EITHER MIRRORED OR NON-MIRRORED
                    USE_UN  = (LCG(IMIL) >= 0)
                    USE_MIR = (LCG(IMIL) <  0)
                ELSE
C                   SEARCH BOTH MIRRORED & NON-MIRRORED IF CHKMIR
                    USE_UN  = .TRUE.
                    USE_MIR = CKMIRROR
                ENDIF

C               CHECK MIRRORED/ NON-MIRRORED POSITION 
                CALL CROSRNG_2(CIRCREF(1,IMI),CIRCEXP,
     &                         LCIRC,NRING, MAXRIN,NUMR,
     &                         .FALSE.,FFTW_PLANS(1),
     &                         USE_UN,USE_MIR,
     &                         ISMIRRORED(IMI),
     &                         TOTMIN(IMI),ROTANGT(IMI))
	     ENDDO
c$omp        end parallel do
         ELSE
C            USE REFERENCE RINGS FILE (MIGHT BE AN INCORE FILE)
             
c$omp        parallel do private(imil,imi,ithread,use_un,use_mir)
c$omp&       schedule(static,1)
             DO IMIL=1,IEND     !LOOP OVER ALL REFERENCE IMAGES
                IMI = IMIL
	        IF (LIMITRANGE) IMI = ABS(LCG(IMIL))

C               FIND THREAD NUMBER 
                ITHREAD = MOD((IMIL-1),NUMTH) + 1
 
C               FILL CIRCREF FROM REFERENCE RINGS FILE
c$omp           critical
                CALL REDLIN(LUNRING,CIRCREF(1,ITHREAD),LCIRC,IMI)
c$omp           end critical

                IF (CKMIRROR .AND. LIMITRANGE) THEN
C                   ONLY SEARCH EITHER MIRRORED OR NON-MIRRORED
                    USE_UN  = (LCG(IMIL) .GE. 0)
                    USE_MIR = (LCG(IMIL) .LT. 0)
                ELSE
C                   SEARCH BOTH MIRRORED & NON-MIRRORED IF CHKMIR
                    USE_UN  = .TRUE.
                    USE_MIR = CKMIRROR
                ENDIF

C               CHECK EITHER MIRRORED/ NON-MIRRORED POSITION 
                CALL CROSRNG_2(CIRCREF(1,ITHREAD),CIRCEXP,
     &                          LCIRC,NRING, MAXRIN,NUMR,
     &                          .FALSE.,FFTW_PLANS(1),
     &                          USE_UN,USE_MIR,
     &                          ISMIRRORED(IMI),
     &                          TOTMIN(IMI),ROTANGT(IMI))
	     ENDDO
comp         end parallel do
          ENDIF

C         LOOP OVER ALL RELEVANT REF. IMAGES -------------------------
          IREF   = 0
          CCROTD = -1.0D20

          DO IMIL=1,IEND
             IMI = IMIL
	     IF (LIMITRANGE) IMI = ABS(LCG(IMIL))

             IF (TOTMIN(IMI) >= CCROTD) THEN
C               GOOD MATCH WITH TOTMIN (MIRRORED OR NOT)  POSITION 
                CCROTD    = TOTMIN(IMI)
                RANGNEW   = ROTANGT(IMI)
                MIRRORNEW = ISMIRRORED(IMI)
                IREF      = IMI
	     ENDIF
          ENDDO   ! END OF: DO IMIL=1,IEND -------------------------


1000      RANGNEW  = (RANGNEW-1) / MAXRIN * DIVAS
          CCROT    = CCROTD
          IMGEXP   = IEXPLIST(IEXP)
          PEAKV    = 0.0
          XSHNEW   = 0.0
          YSHNEW   = 0.0

          IF (IREF <= 0) THEN
C             NO NEARBY REFERENCE IMAGE
              IMGREF = 0
C             IREFT IS FOR REFDIR INDEX
              IREFT  = 1
          ELSE
              IMGREF = IREFLIST(IREF)
C             IREFT IS FOR REFDIR INDEX
              IREFT  = IREF
          ENDIF
         
          IF (IMGREF > 0 .AND. ISHRANGE > 0) THEN
C            DETERMINE SHIFT PARAMETERS, NO NEED TO RELOAD EXP. IMAGE  
             NSAMP = 2*NSAM+2
             NROWP = 2*NROW

             CALL APSHIFT(INPIC, REFPAT,IMGREF,
     &                  NSAM,NROW, NSAMP,NROWP,
     &                  EXPBUF,AVI,SIGI, ISHRANGE,
     &                  RANGNEW,XSHNEW,YSHNEW,MIRRORNEW,PEAKV,IRTFLG)
             IF (IRTFLG .NE. 0) RETURN
          ENDIF

C         AP_END WRITES ALIGNMENT PARAMETERS TO DOC FILE 
          NPROJ = NUMREF
          IF (LIMITRANGE) NPROJ = NUMREFLCG 
 
          CALL AP_END(IEXP,IMGEXP,IMGREF, 
     &         ANGREF(1,IREFT),REFDIR(1,IREFT),
     &         ANGEXP(1,IEXP),EXPDIR,ISHRANGE,
     &         GOTREFANG, NGOTPAR, CCROT,PEAKV,
     &         RANGNEW,XSHNEW,YSHNEW,MIRRORNEW,REFPAT,
     &         NPROJ,CTYPE, LUNDOC,PARLIST)

C         WRITE DATA TO IMAGE HEADER 
          CALL AP_END_HEAD(IMGEXP,EXPPAT,INPIC,PARLIST,8,IRTFLG)

          CALL  AP_STAT_ADD(NGOTPAR,CCROT,PARLIST(10),
     &                    ANGDIFTHR,ANGEXP(8,IEXP),
     &                    CCROTAVG,IBIGANGDIF,ANGDIFAVG,IMPROVCCROT,
     &                    CCROTIMPROV,IWORSECCROT,CCROTWORSE)

       ENDDO  ! END OF: DO IEXP=1,NUMEXP  LOOP OVER ALL EXP. IMAGES ----

      IF (NUMEXP > 1) THEN
C         SAVE CCROT & ANGULAR DISPLACEMENT STATISTICS
          CALL AP_STAT(NUMEXP,ANGDIFTHR,IBIGANGDIF,
     &                 ANGDIFAVG, CCROTAVG,
     &                 IMPROVCCROT,CCROTIMPROV,
     &                 IWORSECCROT,CCROTWORSE,
     &                 NBORDER,NSUBPIX,LUNDOC)
       ENDIF

9999  CLOSE(LUNANG)
      IF (.NOT. CIRCREF_IN_CORE) CLOSE(LUNRING)

      IF (ALLOCATED(EXPBUF)) DEALLOCATE(EXPBUF)
      IF (ALLOCATED(CIRCEXP))DEALLOCATE(CIRCEXP)
      IF (ALLOCATED(TMPBUF)) DEALLOCATE(TMPBUF)
      IF (ASSOCIATED(LCG))   DEALLOCATE(LCG)
      NULLIFY(LCG)

      END





#else

C------------------------  MPI SPECIFIC SUBROUTINE --------------------

         SUBROUTINE APREF_P(IREFLIST,NUMREF,IEXPLIST,NUMEXP,
     &              NSAM,NROW,RANGE,ANGDIFTHR,
     &              NRING,LCIRC,NUMR,CIRCREF,CIRCREF_IN_CORE,
     &              MODE,REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS,
     &              REFPAT,EXPPAT,CKMIRROR,CTYPE,
     &              ROTFIRST,ISHRANGE,LUNDOC,FBS_WANTED)

        INCLUDE 'MAKE_CLOSE_LIST.INC'  
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

	INTEGER                          :: IREFLIST(NUMREF) 
	INTEGER                          :: IEXPLIST(NUMEXP) 
        INTEGER                          :: NUMR(3,NRING)
	REAL                             :: CIRCREF(LCIRC,NUMREF)
	CHARACTER (LEN=1)                :: MODE 
        CHARACTER (LEN=*)                :: REFANGDOC
        CHARACTER (LEN=*)                :: EXPANGDOC
        CHARACTER (LEN=*)                :: SCRFILE
        CHARACTER (LEN=*)                :: REFPAT,EXPPAT 
        CHARACTER (LEN=*)                :: CTYPE 
	LOGICAL                          :: ROTFIRST
        CHARACTER (LEN=74)               :: COMMENT 
        LOGICAL                          :: FBS_WANTED

	DOUBLE PRECISION                 :: CCROTD
	CHARACTER(LEN=1)                 :: NULL = CHAR(0)
	LOGICAL                          :: CIRCREF_IN_CORE
	LOGICAL                          :: CKMIRROR
        LOGICAL                          :: MIRRORNEW
	LOGICAL                          :: GOTREFANG,GOTEXPANG
	LOGICAL                          :: ONLYMIRROR
        LOGICAL                          :: LIMITRANGE
        LOGICAL                          :: MIRRORED
        LOGICAL                          :: USE_UN,USE_MIR
        INTEGER *8                       :: FFTW_PLANS(*)

C       ALLOCATABLE ARRAYS
	REAL, ALLOCATABLE                :: EXPBUF(:) 
 	REAL, ALLOCATABLE                :: CIRCEXP(:) 

C       AUTOMATIC ARRAYS
 	DOUBLE PRECISION                 :: TOTMIN(NUMREF)
 	LOGICAL                          :: ISMIRRORED(NUMREF)

 	REAL                             :: ROTANGT(NUMREF) 

	REAL                             :: DLIST(6)
	REAL                             :: EXPDIR(3)
	REAL                             :: ANGEXP(8,NUMEXP)
	REAL                             :: ANGREF(3,NUMREF)
	REAL                             :: REFDIR(3,NUMREF)
	INTEGER, POINTER                 :: LCG(:)
        INTEGER                          :: NSAID = 0
        INTEGER, PARAMETER               :: NLISTMAX = 15
        REAL                             :: PARLIST(NLISTMAX)

        REAL, PARAMETER                  :: QUADPI     = 3.141592653589
        REAL, PARAMETER                  :: DGR_TO_RAD = (QUADPI/180)

        INTEGER, PARAMETER               :: INPIC   = 77
        INTEGER, PARAMETER               :: LUNANG  = 78
        INTEGER, PARAMETER               :: LUNRING = 50

        INCLUDE 'mpif.h'
        INTEGER, ALLOCATABLE             :: PSIZE(:)
        INTEGER, ALLOCATABLE             :: NBASE(:)
        REAL, ALLOCATABLE                :: PARTAB(:,:),PARTABLOC(:,:)

#ifdef MPI_DEBUG
        DOUBLE PRECISION                 :: TCOM0, TCOM1
#endif

        LOGICAL                          :: ONLYONE_RED,ONLYONE_WRT
        COMMON /COMM_MPI/ONLYONE_RED,ONLYONE_WRT

        ICOMM = MPI_COMM_WORLD
        CALL MPI_COMM_RANK(ICOMM, MYPID, MPIERR)
        CALL MPI_COMM_SIZE(ICOMM, NPROCS, MPIERR)

        NR     = NUMR(1,NRING)     ! OUTER RING NUMBER
        MAXRIN = NUMR(3,NRING)     ! ACTUAL LENGTH OF LONGEST RING
#ifdef SP_LIBFFTW3
        MAXRIN = NUMR(3,NRING) - 2
#endif
        IF (ROTFIRST) THEN
           CALL ERRT(101,'ROTFIRST NOT SUPPORTED UNDER MPI',IDUM)
           STOP
        ENDIF

	RANGECOS  = COS(RANGE*DGR_TO_RAD)

C       FIND NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)

        IF (MODE .EQ. 'H')  THEN
           DIVAS = 180.0
        ELSE
           DIVAS = 360.0
        ENDIF

	ALLOCATE(CIRCEXP(LCIRC),STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'APREF_P; CIRCEXP',LCIRC)
           GOTO 9999
        ENDIF

C       LOAD REF. PROJ. ANGLES (ANGREF) FROM DOC. FILE (REFANGDOC) OR
C       REF. IMAGE FILE (REFPAT) HEAD (IF WANTED)
C       CONVERT REF. ANGLES TO UNITARY DIRECTIONAL VECTORS (REFDIR).
	CALL AP_GETANGAS(IREFLIST,NUMREF,0,REFANGDOC,REFPAT,
     &                  INPIC,LUNANG,3,ANGREF,GOTREFANG,NGOTREF,
     &                  .TRUE.,REFDIR,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       READ REFERENCE IMAGES INTO REFERENCE RINGS (CIRCREF) ARRAY OR
C       CREATE REFERENCE RINGS FILE FOR LATER READING 
        CALL APRINGS_NEW(IREFLIST,NUMREF, NSAM,NROW,
     &               NRING,LCIRC,NUMR,MODE,FFTW_PLANS, 
     &               REFPAT,INPIC,CIRCREF,CIRCREF_IN_CORE,
     &               LUNRING,SCRFILE,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (NSAID .LE. 0) THEN
           IF (CIRCREF_IN_CORE) THEN
              IF (MYPID .LE. 0) WRITE(NOUT,91)NUMTH
91            FORMAT('  Ref. rings in core,  Threads: ',I4)
           ELSE
              IF (MYPID .LE. 0) WRITE(NOUT,92)NUMTH
92            FORMAT('  Ref. rings not in core,  Threads: ',I4)
           ENDIF
           NSAID = NSAID + 1
        ENDIF

C       LOAD EXP. PROJ. ANGLES & ALIGNMENT PARAMETERS 
C       FROM DOC. FILE OR IMAGE HEADER (IF WANTED)
        CALL AP_GETANGAS(IEXPLIST,NUMEXP,0,EXPANGDOC,EXPPAT,
     &                   INPIC,LUNANG,8,ANGEXP,GOTEXPANG,NGOTPAR,
     &                   .FALSE.,EXPDIR,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       INITIALIZE CCROT STATISTICS
        CALL  AP_STAT_ADD(-1,CCROT,ANGDIF,ANGDIFTHR,CCROTLAS,
     &                   CCROTAVG,IBIGANGDIF,ANGDIFAVG,IMPROVCCROT,
     &                   CCROTIMPROV,IWORSECCROT,CCROTWORSE)

C       DISTRIBUTE EXPERIMENTAL IMAGES 
        ALLOCATE(PSIZE(NPROCS), NBASE(NPROCS), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'APREF_P; NBASE',2*NPROCS)
           GOTO 9999
        ENDIF

C       FILL PSIZE WITH PARTITION LIMITS, PSIZE(1) = NUMEXP / NPROCS
        CALL SETPART(NUMEXP, PSIZE, NBASE)

#ifdef MPI_DEBUG
        WRITE(6,111) NBASE(MYPID+1), MYPID
 111    FORMAT('  APREF_P; NBASE(MYPID+1): ', I5, ' MYPID: ', I5)
        CALL FLUSHFILE(6)
        CALL MPI_BARRIER(ICOMM,MPIERR)
#endif

        NLOC = PSIZE(MYPID+1)
        ALLOCATE(EXPBUF(NSAM*NROW),
     &           PARTAB(15,NUMEXP), 
     &           PARTABLOC(15,NLOC),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MWANT = NSAM*NROW + 15*(NUMEXP+NLOC) 
           CALL ERRT(46,'APREF_P: PARTAB...',MWANT)
           GOTO 9999
        ENDIF

C       ZERO THE ARRAYS
        PARTAB    = 0.0
        PARTABLOC = 0.0

C       CALCULATE DIMENSIONS FOR APRINGS
        CNS2 = NSAM / 2 + 1
        CNR2 = NROW / 2 + 1

C       FOLLOWING LOOP PERFORMED BY ALL PROCESSORS ON LOCAL DATA

c       write(0,*)' APREF_p; barrier 0,nloc,nprocs: ',nloc,nprocs,mypid
        CALL MPI_BARRIER(ICOMM,MPIERR)

        ONLYONE_RED = .FALSE.  ! DO NOT BROADCAST LOADED DATA
        ONLYONE_WRT = .FALSE.  ! USED TO be NEEDED FOR NORM3 IN AP_GETDATS
                               ! BUT AP_GETDATA DOES NOT WRITE ANYWAY!!

        DO JLOC = 1, NLOC   ! ---------------------------------------
          INUM   = NBASE(MYPID+1) + JLOC
          IMGEXP = IEXPLIST(INUM)
          !write(0,*)' APREF_p;  jloc,inum: ',jloc,inum,imgexp,mypid
         
C         LOAD CURRENT EXPERIMENTAL IMAGE INTO ARRAY EXPBUF
	  CALL AP_GETDATA(IEXPLIST,NUMEXP,
     &                    NSAM,NROW, NSAM,NROW,0.0,
     &                    1,EXPPAT,INPIC, INUM,INUM,
     &                    ONLYONE_RED,EXPBUF,
     &                    .TRUE.,AVI,SIGI, IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9999

C         CONVERT EXP. ANGLE TO UNITARY DIRECTIONAL VECTORS (EXPDIR).
          CALL AP_GETSATA(ANGEXP(1,INUM),EXPDIR,8,1,IRTFLG)

C         EXTRACT EXP. IMAGE POLAR COORD. RINGS, NORMALIZE & FFT THEM
	  CALL APRINGS_ONE_NEW(NSAM,NROW, CNS2,CNR2, EXPBUF,
     &                   .FALSE.,MODE,NUMR,NRING,LCIRC, 0.0,FFTW_PLANS,
     &                   CIRCEXP,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9999

C         DETERMINE WHICH REF IMAGES ARE TO BE COMPARED
          IEND       = NUMREF
          CCROTD     = -1.0D20
          LIMITRANGE = RANGECOS .LT. 1.0
          NULLIFY(LCG)  ! FOR INTEL COMPILER
 
C         IF LIMITRANGE, MAKE LIST OF NEARBY REF IMAGES, RETURNS: NUMREFLCG
          CALL MAKE_CLOSE_LIST(NUMREF,LIMITRANGE,
     &                         REFDIR,EXPDIR,
     &                         RANGECOS, .TRUE., 
     &                         LCG, NUMREFLCG, IRTFLG)

          IF (NUMREFLCG <= 0) THEN
C            NO REF. IMAGE WITHIN COMPARISON ANGLE
             IREF = 0 
             GOTO 1000
          ENDIF

          IEND = NUMREFLCG

C         REF. IMAGES FOUND WITHIN COMPARISION RANGE

          IF (CIRCREF_IN_CORE) THEN
C            USE CIRCREF FOR REFERENCE RINGS

             DO IMIL=1,IEND
                IMI = IMIL
	        IF (LIMITRANGE) IMI = ABS(LCG(IMIL))

                IF (CKMIRROR .AND. LIMITRANGE) THEN
C                   ONLY SEARCH EITHER MIRRORED OR NON-MIRRORED
                    USE_UN  = (LCG(IMIL) .GE. 0)
                    USE_MIR = (LCG(IMIL) .LT. 0)
                ELSE
C                   SEARCH BOTH MIRRORED & NON-MIRRORED IF CHKMIR
                    USE_UN  = .TRUE.
                    USE_MIR = CKMIRROR
                ENDIF

C               CHECK MIRRORED / NON-MIRRORED POSITION 
                CALL CROSRNG_2(CIRCREF(1,IMI),CIRCEXP,
     &                         LCIRC,NRING, MAXRIN,NUMR,
     &                         .FALSE.,FFTW_PLANS(1),
     &                         USE_UN,USE_MIR,
     &                         ISMIRRORED(IMI),
     &                         TOTMIN(IMI),ROTANGT(IMI))

	     ENDDO  ! END OF:   DO IMIL=1,IEND 
          ELSE
C            USE REFERENCE RINGS FILE (MIGHT BE AN INCORE FILE)

             DO IMIL=1,IEND     !LOOP OVER ALL REFERENCE IMAGES
                IMI = IMIL
	        IF (LIMITRANGE) IMI = ABS(LCG(IMIL))

C               FIND THREAD NUMBER 
                ITHREAD = MOD((IMIL-1),NUMTH) + 1
 
C               FILL CIRCREF FROM REFERENCE RINGS FILE
                CALL REDLIN(LUNRING,CIRCREF(1,ITHREAD),LCIRC,IMI)

                IF (CKMIRROR .AND. LIMITRANGE) THEN
C                   ONLY SEARCH EITHER MIRRORED OR NON-MIRRORED
                    USE_UN  = (LCG(IMIL) .GE. 0)
                    USE_MIR = (LCG(IMIL) .LT. 0)
                ELSE
C                   SEARCH BOTH MIRRORED & NON-MIRRORED IF CHKMIR
                    USE_UN  = .TRUE.
                    USE_MIR = CKMIRROR
                ENDIF

C               CHECK EITHER MIRRORED/ NON-MIRRORED POSITION 
                CALL CROSRNG_2(CIRCREF(1,ITHREAD),CIRCEXP,
     &                         LCIRC,NRING, MAXRIN,NUMR,
     &                         .FALSE.,FFTW_PLANS(1),
     &                         USE_UN,USE_MIR,
     &                         ISMIRRORED(IMI),
     &                         TOTMIN(IMI),ROTANGT(IMI))

	     ENDDO     ! END OF:   DO IMIL=1,IEND     
          ENDIF

C         LOOP OVER ALL RELEVANT REF. IMAGES
          IREF   = 0
          CCROTD = -1.0D20

          DO IMIL=1,IEND
             IMI = IMIL
	     IF (LIMITRANGE) IMI = ABS(LCG(IMIL))

             IF (TOTMIN(IMI) .GE. CCROTD) THEN
C               GOOD MATCH WITH TOTMIN (MIRRORED OR NOT)  POSITION 
                CCROTD    = TOTMIN(IMI)
                RANGNEW   = ROTANGT(IMI)
                MIRRORNEW = ISMIRRORED(IMI)
                IREF      = IMI
	     ENDIF
          ENDDO   ! END OF: DO IMIL=1,IEND


1000      RANGNEW  = (RANGNEW-1) / MAXRIN * DIVAS
          CCROT    = CCROTD
          PEAKV    = 0.0
          XSHNEW   = 0.0
          YSHNEW   = 0.0

          IF (IREF .LE. 0) THEN
C             NO NEARBY REFERENCE IMAGE
              IMGREF = 0
C             IREFT IS FOR REFDIR INDEX
              IREFT  = 1
          ELSE
              IMGREF = IREFLIST(IREF)
C             IREFT IS FOR REFDIR INDEX
              IREFT  = IREF
          ENDIF

c         write(0,*)' APREF_p; above apshiftn,imgexp: ',imgexp,mypid
          IF (IMGREF .GT. 0 .AND. ISHRANGE .GT. 0) THEN
C            DETERMINE SHIFT PARAMETERS FOR EXP IMAGE IN EXPBUF 
             NSAMP = 2*NSAM+2
             NROWP = 2*NROW

             CALL APSHIFT(INPIC, REFPAT,IMGREF,
     &                   NSAM,NROW, NSAMP,NROWP,
     &                   EXPBUF,AVI,SIGI, ISHRANGE,
     &                   RANGNEW,XSHNEW,YSHNEW,MIRRORNEW,PEAKV,IRTFLG)
             IF (IRTFLG .NE. 0) RETURN
          ENDIF

          NPROJ = NUMREF
          IF (LIMITRANGE) NPROJ = NUMREFLCG   

C         AP_END PUTS ALIGNMENT PARAMETERS FOR THIS IMAGE IN PARTABLOC

          CALL AP_END(IMGEXP,IMGEXP,IMGREF, 
     &         ANGREF(1,IREFT),REFDIR(1,IREFT),
     &         ANGEXP(1,INUM),EXPDIR,ISHRANGE,
     &         GOTREFANG, NGOTPAR, CCROT,PEAKV,
     &         RANGNEW,XSHNEW,YSHNEW,MIRRORNEW,REFPAT,
     &         NPROJ, CTYPE, LUNDOC,PARTABLOC(1,JLOC))

          CALL AP_STAT_ADD(NGOTPAR,CCROT,PARTABLOC(10,JLOC),
     &                    ANGDIFTHR,ANGEXP(8,INUM),
     &                    CCROTAVG,IBIGANGDIF,ANGDIFAVG,IMPROVCCROT,
     &                    CCROTIMPROV,IWORSECCROT,CCROTWORSE)

        ENDDO     ! END OF:  DO JLOC = 1, NLOC ----------------------

c       write(0,*) ' APREF_p, barrier 2, mypid: ' ,mypid

        CALL MPI_BARRIER(ICOMM,MPIERR)
        ONLYONE_RED = .TRUE.  ! SET BCAST FLAG BACK TO DEFAULT
        ONLYONE_WRT = .TRUE.  ! SET BCAST FLAG BACK TO DEFAULT

C       COLLECT ALL ALIGNMENT PARAMETERS INTO PARTAB
   
        DO IPROC = 1, NPROCS
           PSIZE(IPROC) = 15 * PSIZE(IPROC)
           NBASE(IPROC) = 15 * NBASE(IPROC)
        ENDDO
        CALL MPI_ALLGATHERV(PARTABLOC, PSIZE(MYPID+1),
     &                      MPI_REAL, PARTAB, PSIZE, NBASE,
     &                      MPI_REAL, ICOMM,  MPIERR)

C       WRITE IS SYNCHRONIZED WITHIN LUNDOCWRTDAT
        IF (LUNDOC .GT. 0) THEN
C          SAVE IN ALIGNMENT DOC FILE
C          <,<,<, MIR-REF#,IMG#,INPLANE<, SX,SY,NPROJ,
C                <DIF,CCROT,INPLANE<,SX,SY
           NWANT = 15 
           DO IT = 1, NUMEXP
              CALL LUNDOCWRTDAT(LUNDOC,IT,PARTAB(1,IT),NWANT,IRTFLG)
           ENDDO
        ENDIF

        IF (NUMEXP .GT. 1) THEN   ! I SUSPECT SUMS ARE WRONG FOR MPI!
C         SAVE CCROT & ANGULAR DISPLACEMENT STATISTICS
          CALL AP_STAT(NUMEXP,ANGDIFTHR,IBIGANGDIF,
     &                 ANGDIFAVG, CCROTAVG,
     &                 IMPROVCCROT,CCROTIMPROV,
     &                 IWORSECCROT,CCROTWORSE,
     &                 NBORDER,NSUBPIX,LUNDOC)
        ENDIF

9999   IF (ALLOCATED(PSIZE))     DEALLOCATE(PSIZE)
       IF (ALLOCATED(NBASE))     DEALLOCATE(NBASE)
       IF (ALLOCATED(PARTABLOC)) DEALLOCATE(PARTABLOC)
       IF (ALLOCATED(PARTAB))    DEALLOCATE(PARTAB)
       IF (ALLOCATED(EXPBUF))    DEALLOCATE(EXPBUF)
       IF (ALLOCATED(CIRCEXP))   DEALLOCATE(CIRCEXP)
       IF (ASSOCIATED(LCG))      DEALLOCATE(LCG)
       NULLIFY(LCG)

       CLOSE(LUNANG)
       IF (.NOT. CIRCREF_IN_CORE) CLOSE(LUNRING)

       END

#endif
C ------------------------- END OF MPI CODE -----------------------









