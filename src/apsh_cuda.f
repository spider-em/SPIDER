C++*********************************************************************
C                                                                      *
C    APSH_CUDA.F   ADAPTED FROM APSH_PS           APR 10 ARDEAN LEITH  *
C                  AP_STAT NBORDER                OCT 10 ARDEAN LEITH
C                  APSH not AP_SH                                                             *
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
C  APSH_CUDA 
C  
C  PURPOSE: FIND ROTATIONAL AND SHIFT PARAMETERS TO ALIGN A SERIES OF
C           REFERENCE IMAGES WITH SAMPLE IMAGES
C           USED IF:  GPU USED 
C
C PARAMETERS:
C       IREFLIST            LIST OF REF. IMAGE FILE NUMBERS   (INPUT)
C       NUMREF              NO. OF IMAGES                     (INPUT)
C       IEXPLIST            LIST OF EXP. IMAGE FILE NUMBERS   (INPUT)
C       NUMEXP              NO. OF IMAGES                     (INPUT)
C       NSAM,NROW           ACTUAL (NOT-WINDOWED) IMAGE SIZE  (INPUT)
C       REFANGDOC           REF. ANGLES FILE NAME             (INPUT)
C       EXPANGDOC           EXP. ANGLES FILE NAME             (INPUT)
C       REFPAT              REF. IMAGE SERIES FILE TEMPLATE   (INPUT)
C       EXPPAT              EXP. IMAGE SERIES FILE TEMPLATE   (INPUT)
C       LUNDOC              OUTPUT UNIT                       (INPUT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

#ifdef SP_CUDA
C$PRAGMA C(chkrings_cuda)
C$PRAGMA C(chkrings_cuda_alloc)
#endif

       SUBROUTINE APSH_CUDA(IREFLIST,NUMREF, IEXPLIST,NUMEXP, 
     &               NSAM,NROW, ISHRANGEX,ISHRANGEY,ISTEP,
     &               NUMR,NRING,
     &               MODE, REFANGDOC,EXPANGDOC,FFTW_PLANS,
     &               REFPAT,EXPPAT,LUNDOC)


#ifdef SP_CUDA
      USE CUDAFOR              ! FOR PINNED ALLOCATE?
      USE CHKRINGS_CUDA_INFO
#endif

      IMPLICIT NONE

      INCLUDE 'CMLIMIT.INC'
      INCLUDE 'CMBLOCK.INC'

      INTEGER                    :: IREFLIST(NUMREF) 
      INTEGER                    :: IEXPLIST(NUMEXP) 
      INTEGER                    :: NUMR(3,NRING)

      CHARACTER(LEN=1)           :: MODE,CDUM
      CHARACTER(LEN=*)           :: REFANGDOC,EXPANGDOC
      CHARACTER(LEN=*)           :: REFPAT,EXPPAT 

      INTEGER *8                 :: FFTW_PLANS(*) !POINTERS TO STRUCTURES 
      INTEGER *8                 :: I8,MEMCUDA,NUMVAR

C     AUTOMATIC ARRAYS
      REAL                       :: ANGOUT(3)
      REAL                       :: WR(NRING)

C     ALLOCATED ARRAYS
      REAL, ALLOCATABLE          :: DLIST(:,:) 
      REAL, ALLOCATABLE          :: REFDIR(:,:),EXPDIR(:,:) 
      REAL, ALLOCATABLE          :: ANGREF(:,:),ANGEXP(:,:)
      REAL, ALLOCATABLE          :: FSHV(:,:)

      INTEGER, ALLOCATABLE       :: ISHV(:,:)
      INTEGER, ALLOCATABLE       :: NLOCS(:,:)
      INTEGER, ALLOCATABLE       :: MAX_POS(:)

      !COMPLEX,PINNED,ALLOCATABLE :: CIRCREF(:,:)
      !COMPLEX,PINNED,ALLOCATABLE :: CIRCEXP(:,:)
      !COMPLEX,PINNED,ALLOCATABLE :: QC(:,:,:)

      COMPLEX, ALLOCATABLE       :: CIRCREF(:,:) 
      COMPLEX, ALLOCATABLE       :: CIRCEXP(:,:)
      COMPLEX, ALLOCATABLE       :: QC(:,:,:)

      REAL,    ALLOCATABLE       :: COEFFS(:,:)
      INTEGER, ALLOCATABLE       :: IXY(:,:)
 
      LOGICAL                    :: LPINFLG
      LOGICAL                    :: GOTREFANG,WEIGHT,GOTEXPANG
      LOGICAL                    :: MIRRORNEW

      REAL                       :: ANGDIFAVG,CCROTIMPROV,CCROTWORSE
      REAL                       :: TDUM,WRT
      REAL                       :: ANGDIFTHR,ANGDIF,SCALE
      REAL                       :: CCROTLAS,CCROTAVG,IWORSECCROT,FDUM
      REAL                       :: CCROT,SX,SY,CO,SO,XSHSUM,YSHSUM
      REAL                       :: PEAKV,RANGNEW

      INTEGER                    :: NUMREF,NUMEXP,NRING,NUMTH,NLOAD
      INTEGER                    :: NSAM,NROW,ISHRANGEX,ISHRANGEY
      INTEGER                    :: ISTEP,IBIGANGDIF,NEXPpL
      INTEGER                    :: LUNDOC,MYPID,NWANTOUT,IMPROVCCROT
      INTEGER                    :: NSAMS,NROWS,NRAYS,NRAYSC 
      INTEGER                    :: IRTFLG,LCIRCC,LCIRCR,NSHIFTS,MAXL
      INTEGER                    :: JT,IT,MAXRIN,MWANT,NCAN,NCANREF
      INTEGER                    :: NCANEXP 
      INTEGER                    :: NSHEXP,NREFGOT 
      INTEGER                    :: IM,ISHIFTM,IREFM,IRAYM,NDUM,I
      INTEGER                    :: IEXP,IREF,ISH,NREFLEFT
      INTEGER                    :: IEXPGO,NEXPDONE,NREFDONE
      INTEGER                    :: NGOTREF,NGOTPAR
      INTEGER                    :: IEXPSET,JEXP,ILOC,MAX_POS_SIZE
      INTEGER                    :: NSHIFTST
     

      REAL                       :: TLOAD = 0.0, TCALC = 0.0
      REAL                       :: TFFT  = 0.0, TBLAS = 0.0
      REAL                       :: TSUB  = 0.0, TALL  = 0.0

      REAL, PARAMETER            :: QUADPI = 3.14159265358979323846
      REAL, PARAMETER            :: DGR_TO_RAD = (QUADPI/180)

      INTEGER, PARAMETER         :: LUNIMG   = 77
      INTEGER, PARAMETER         :: INANG    = 78
      INTEGER, PARAMETER         :: LUNRING  = 0 !!!!!
      INTEGER, PARAMETER         :: NLISTMAX = 15

      INTEGER                    :: NBORDER = 0       ! UNUSED
      INTEGER                    :: NSUBPIX = 0       ! UNUSED

      REAL                       :: ANG_N      ! FUNCTION

#ifndef SP_CUDA
      CALL  ERRT(101,'NOT COMPILED FOR CUDA',NDUM)
      RETURN

      END

#else

      MYPID = -1
 
C     SET TYPE OF OUTPUT DOC FILES WANTED
      NWANTOUT = 15

C     FIND NUMBER OF OMP THREADS  (UNUSED FOR NOW)
      CALL GETTHREADS(NUMTH)      
      CALL FLUSHRESULTS()

C     INITIALIZE CCROT STATISTICS COUNTERS
      ANGDIFTHR   = 0.0
      CALL AP_STAT_ADD(-1,CCROT,ANGDIF,ANGDIFTHR,CCROTLAS,
     &                   CCROTAVG,IBIGANGDIF,ANGDIFAVG,IMPROVCCROT,
     &                   CCROTIMPROV,IWORSECCROT,CCROTWORSE)

C     NEED REFANGLES FILE FOR  'SH'
      ALLOCATE(REFDIR(3,NUMREF),
     &         ANGREF(3,NUMREF), STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
          CALL ERRT(46,'REFDIR, ANGREF',6*NUMREF)
          GOTO 9999
      ENDIF
 
C     READ REF. ANGLES INTO ANGREF
C     CONVERT REF. ANGLES TO UNITARY DIRECTIONAL VECTORS (REFDIR).
      CALL AP_GETANGAS(IREFLIST,NUMREF,0,REFANGDOC,REFPAT,
     &                LUNIMG,INANG,3,ANGREF,GOTREFANG,NGOTREF,
     &                .TRUE.,REFDIR,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

      NGOTPAR = 0
      IF (EXPANGDOC .NE. CHAR(0)) THEN

         ALLOCATE(ANGEXP(8,NUMEXP), 
     &            EXPDIR(3,NUMEXP), STAT=IRTFLG)
	 IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'APSH_CUDA; ANGEXP....',11*NUMEXP)
            RETURN
         ENDIF 

C        READ EXP. ANGLES INTO ANGEXP
C        CONVERT EXP. ANGLES TO UNITARY DIRECTIONAL VECTORS(EXPDIR).
	 CALL AP_GETANGAS(IEXPLIST,NUMEXP,0,EXPANGDOC,EXPPAT,
     &                    LUNIMG,INANG,8,ANGEXP,GOTEXPANG,NGOTPAR,
     &                    .TRUE.,EXPDIR,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999
      ELSE
C        DUMMY ALLOCATE TO AVOID BUS ERROR ON SOME SYSTEMS
	 ALLOCATE(ANGEXP(8,1), 
     &           EXPDIR(3,1), STAT=IRTFLG)
	 IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'APSH_CUDA; ANGEXP....',11)
            GOTO 9999
         ENDIF 
      ENDIF

      NRAYS  = NUMR(3,NRING)         ! MAXIMUM RING LENGTH (REAL)
      NRAYSC = NUMR(3,NRING) / 2     ! # OF RAYS   (FOURIER)

      ALLOCATE(NLOCS(2,NRAYSC+1),STAT=IRTFLG) ! ADD AN EXTRA LINE
      IF (IRTFLG .NE. 0) THEN
         CALL ERRT(46,'APSH_CUDA; NLOCS',2*NRAYSC+2)   
         GOTO 9999
      ENDIF

C     SET # OF POINTS ON EACH RAY AND RAY STARTING INDEX IN CIRC
      CALL APRINGS_TRANS_LOCS(NUMR,NRING, NLOCS,NRAYSC)

      LCIRCC  = NLOCS(2,NRAYSC+1)-1 ! TOTAL RINGS LENGTH   (CMPLX)
      LCIRCR  = LCIRCC * 2          ! TOTAL RINGS LENGTH   (REAL)

      WRITE(6,97)'# OF RINGS:',NRING
      WRITE(6,97)'# PTS ON A RING = # RAYS:',NRAYS
      WRITE(6,97)'RAYORDER  RINGS LEN:',LCIRCR
      WRITE(6,97)'RAYORDER  RINGS LEN (CMPLX):',LCIRCC

 97   FORMAT('  ',A,T32,I10)
 98   FORMAT('  ',A,T32,I10,  '    ',A,I5)
 99   FORMAT('  ',A,T32,F10.3,'    ',A,I5)

C     RINGWE RETURNS WR WEIGHTS
      MAXRIN = NUMR(3,NRING) - 2
      CALL RINGWE_NEW(WR,NUMR,NRING,MAXRIN)
      IF (MODE .EQ. 'H') WR = WR * 0.5

      NSHIFTS = ((ISHRANGEX/ISTEP)*2 + 1) * ((ISHRANGEY/ISTEP)*2 + 1)
      ALLOCATE(ISHV(2,NSHIFTS), 
     &         DLIST(NLISTMAX,NUMEXP),
     &         COEFFS(6,LCIRCR),
     &            IXY(2,LCIRCR),      STAT=IRTFLG)
      IF (IRTFLG.NE.0) THEN
         MWANT = 4*NSHIFTS + +NLISTMAX*NUMEXP + 6*LCIRCR + 2*LCIRCR
         CALL ERRT(46,'ISHV, ...',MWANT)
         GOTO 9999
      ENDIF 

      NSHIFTS = 0                    ! LAZY NSHIFTS COUNTER
      DO JT=-ISHRANGEY,ISHRANGEY,ISTEP
         DO IT=-ISHRANGEX,ISHRANGEX,ISTEP
            NSHIFTS       = NSHIFTS + 1
            ISHV(1,NSHIFTS) = IT
            ISHV(2,NSHIFTS) = JT
        ENDDO
      ENDDO
      WRITE(6,97)'# SHIFTS:',NSHIFTS 

      DLIST   = 0                      ! FOR CCROT TEST
      COEFFS  = 0.0
      IXY     = -100                   ! FOR HANDLING CIRC PADS
      COEFFS  = 0.0

      SCALE   = 1.0 / FLOAT(NRAYS-2)   ! FFT SCALING
      WRT     = 0.0                    ! USED WHEN NO IMAGE WEIGHTING

C     INITIALIZE CUDA, GET #CIRCS THAT FIT IN GPU ------------- CUDA
      CALL CUDA_INIT(LCIRCC,NRAYSC,NUMEXP,NUMREF,NCAN,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

      WRITE(6,97)    'MAX # RINGS (FOR: #REF=#EXP):',NCAN
      WRITE(NOUT,97) 'MAX # RINGS (FOR GPU: #REF=#EXP):',NCAN

C     LOAD NLOCS INTO GPU
      !write(6,98)'loading nlocs, nraysc:', nraysc
      CALL CHKRINGS_CUDA_LOAD(LCIRCC,
     &                       .FALSE.,CIRCEXP,NDUM, 
     &                       .FALSE.,CIRCREF,NDUM, 
     &                       .TRUE.,NLOCS,NRAYSC,
     &                        TLOAD, IRTFLG)

        !!!!!!!!!!!!!!!! QUERY
      MEMCUDA      = 1610350592       ! 1,610,350,592 SHOULD USE QUERY !!!!!
      NUMVAR       = FLOAT(MEMCUDA) / 8.0
      MAX_POS_SIZE = 0
      PEAKV        = 1.0

C     LOOP OVER ALL REF IMAGE GROUPS
      NREFDONE  = 0
 
      DO WHILE(NREFDONE < NUMREF)         ! LOOP OVER ALL REF. IMAGES
         NREFLEFT  = NUMREF - NREFDONE
     
         IF ((NREFLEFT + NSHIFTS) <= NCAN) THEN
C           CAN FIT ALL REFs * >= 1 SET OF SHIFTED EXPs
            NCANREF = NREFLEFT                   ! # REF IMAGES TO LOAD
            NEXPpL  = (NCAN-NREFLEFT)/NSHIFTS    ! # EXP/LOOP (INT. DIVISION
            NEXPpL  = MIN(NEXPpl,NUMEXP)
            NCANEXP = NEXPpL * NSHIFTS           ! # SHIFTED EXP IMAGES TO LOAD

         ELSEIF (NSHIFTS <= ((2*NCAN)/3)) THEN
C           CAN FIT  >= ONE SET OF SHIFTED EXPs AND SOME REFS
            NEXPpL  = (2*NCAN/3) / NSHIFTS     ! # EXP / LOOP (INT. DIVISION)
            NEXPpL  = MIN(NEXPpl,NUMEXP)
            NCANEXP = NSHIFTS * NEXPpL         ! # SHIFTED EXP IMAGES TO LOAD

            !I8 = 8*(NCAN*LCIRCC + NCANEXP*LCIRCC + NRAYSC*NCAN*NCANEXP)
            NCANREF = (FLOAT(NUMVAR) - FLOAT(NCANEXP*LCIRCC)) / 
     &              FLOAT(NRAYSC*NCANEXP)   ! # REF IMAGES TO LOAD
         ELSEIF (NSHIFTS < (0.9 * FLOAT(NCAN))) THEN
C           CAN FIT ONE SET OF SHIFTED EXPs AND A FEW REFS
            NEXPpL  = 1                        ! # EXP / LOOP (INT. DIVISION)
            NEXPpL  = MIN(NEXPpl,NUMEXP)
            NCANEXP = NSHIFTS                  ! # SHIFTED EXP IMAGES TO LOAD
            NCANREF = (FLOAT(NUMVAR) - FLOAT(NCANEXP*LCIRCC)) / 
     &                 FLOAT(NRAYSC*NCANEXP)   ! # REF IMAGES TO LOAD
         ELSE
C           CAN NOT FIT ONE SET OF SHIFTED EXPs AND A FEW REFS
            MWANT = ((NCAN * 2)/3)+1
            CALL ERRT(102,'DECREASE # SHIFTED POSITIONS TO <',MWANT)
            GOTO 9999
         ENDIF

         WRITE(6,97)'# REFs THAT FIT:',     NCANREF
         WRITE(6,97)'# EXPs/ LOOP:',        NEXPpL
         WRITE(6,97)'# SHIFTED EXPs USED:', NCANEXP

         NCANREF = MIN(NCANREF,NREFLEFT)  
         WRITE(6,97)'# REFs USED:',         NCANREF

         WRITE(6,97)'GPU MEMORY :',         MEMCUDA
         I8 = (8*NCANEXP*LCIRCC) + 
     &        (8*NCANREF*LCIRCC) +
     &        (8*NRAYSC*NCANREF*NCANEXP)    ! GLPU MEM. USED IN BYTES
         WRITE(6,97)'MEMORY USED:',         I8

         IF (NEXPpL > MAX_POS_SIZE) THEN
            IF (ALLOCATED(MAX_POS)) DEALLOCATE(MAX_POS)
            MAX_POS_SIZE = NEXPpL
            ALLOCATE(MAX_POS(MAX_POS_SIZE),  STAT=IRTFLG)
            IF (IRTFLG .NE. 0) THEN
               CALL ERRT(46,'MAX_POS_SIZE',MAX_POS_SIZE)
               GOTO 9999
            ENDIF
         ENDIF

         !write(6,*) ' Allocating circref(',lcircc,',',ncanref,')'
         !write(6,*) ' Allocating circexp(',lcircc,',',ncanexp,')'
         MWANT = 2*LCIRCC*NCANEXP + 
     &           2*LCIRCC*NCANREF + 
     &           2*NRAYSC*NCANREF*NCANEXP
         ALLOCATE(CIRCREF(LCIRCC,NCANREF),  
     &            CIRCEXP(LCIRCC,NCANEXP),  
     &            QC(NRAYSC,NCANREF,NCANEXP),   STAT=IRTFLG)
!     &           PINNED=LPINFLG, STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'CIRCREF, CIRCEXP & QC',MWANT)
            GOTO 9999
         !ELSEIF (.NOT. LPINFLG) THEN
             !WRITE(NOUT,98)'CIRCREF, CIRCEXP & QC not pinned, Size: ',MWANT
         ENDIF
 
C        ALLOCATE SPACE ON GPU FOR CEDEV, CRDEV, & QCDEV     
         CALL CHKRINGS_CUDA_ALLOC(.TRUE.,LCIRCC,
     &                             NCANREF,NCANEXP,NRAYSC,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999

C        FILL CIRCREF WITH: NCANREF  WEIGHTED IMAGE
C        READ REF. IMAGES INTO TRANSPOSED REF. RINGS ARRAY (CIRCREF) 
C        CREATES & SAVES COEFFS FOR LATER USE 

         !write(6,98)'filling circref:',nrefdone+1,'...',nrefdone+numref-1
         WEIGHT = .TRUE.    ! REF. IMAGES HAVE WEIGHTED FFT'S
         CALL APRINGS_FILL_COEFT(IREFLIST(NREFDONE+1),NCANREF,
     &                           REFPAT,LUNIMG,
     &                           NSAM,NROW,
     &                           0,0,1,
     &                           NUMR,NRING,  NLOCS,NRAYSC,
     &                           COEFFS,IXY,  FFTW_PLANS,
     &                           CIRCREF,LCIRCR,NCANREF,
     &                           1, NREFGOT,
     &                           LUNRING,WEIGHT,MODE,IRTFLG)

         IF (IRTFLG .NE. 0) GOTO 9999
         IF (NREFGOT .NE. NCANREF) THEN
            CALL ERRT(102,'GOT TOO FEW IMAGES',NREFGOT)
            GOTO 9999
         ENDIF
         !call chkcmplx('Ref: 1',  circref,lcircc, 32006,1000, 0)

         NREFDONE = NREFDONE + NREFGOT

C        LOAD CIRCREF INTO GPU
         !write(6,98)'LOADING CIRCREF, IMAGES', ncanref
         CALL CHKRINGS_CUDA_LOAD(LCIRCC,
     &                           .FALSE.,CIRCEXP,NDUM, 
     &                           .TRUE., CIRCREF,NCANREF, 
     &                           .FALSE.,NLOCS,NDUM,
     &                           TLOAD, IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999

         WRITE(6,*)' '

C        READ EXP. IMAGES INTO CIRCEXP RINGS ARRAY ------------- EXP
         WEIGHT = .FALSE.                ! EXP. IMAGES HAVE UNWEIGHTED FFT'S
         NSHEXP   = 0                    ! # OF SHIFTED  EXP. IMAGES
         NEXPDONE = 0                    ! # OF EXP IMAGES FINISHED
         SCALE    = 1.0 / FLOAT(NRAYS-2) ! FFT SCALING FACTOR
         MIRRORNEW = .FALSE.             ! NOT IMPLEMENTED (NO SPACE)

         DO IEXPGO = 1,NUMEXP,NEXPpL   ! LOOP OVER ALL EXP. IMAGE SETS

C           LOAD EXP IMAGE(S), SHIFT, CONVERT TO POLAR RINGS, NORMALIZE,  
C           FFT, WEIGHT, & TRANSPOSE THE RINGS

            NLOAD = MIN(NEXPpl,NUMEXP-IEXPGO+1)  ! # EXP IMAGES TO LOAD

            !write(6,98)'Filling CIRCEXP:',iexpgo,'..',iexplist(iexpgo+nload-1)
            CALL APRINGS_FILL_COEFT(IEXPLIST(IEXPGO),NLOAD,
     &                              EXPPAT,LUNIMG,
     &                              NSAM,NROW,
     &                              ISHRANGEX,ISHRANGEY,ISTEP,
     &                              NUMR,NRING,  NLOCS,NRAYSC,
     &                              COEFFS,IXY,  FFTW_PLANS,
     &                              CIRCEXP,LCIRCR,NCANEXP,
     &                              1, NSHEXP,
     &                              LUNRING,WEIGHT, MODE, IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999
            !call chkring(1,.true.,circexp(1,1),lcircr,numr,nring,nlocs,nraysc)
            !call chkray (1,.true.,circexp(1,1),lcircr,numr,nring,nlocs,nraysc)

C           LOAD NSHEXP CIRCEXP RINGS INTO GPU
            !write(6,98)'Loading CIRCEXP:',iexpgo,'Shifted Exps:',nshexp
            CALL CHKRINGS_CUDA_LOAD(LCIRCC,
     &                             .TRUE., CIRCEXP,NSHEXP, 
     &                             .FALSE.,CIRCREF,NCANREF,
     &                             .FALSE.,NLOCS,NRAYSC,
     &                             TLOAD,IRTFLG)

C           CARRY OUT CC SEARCH HERE
            !write(6,98)'Calling chkrings', nshexp,' ',NSHEXP
            CALL CHKRINGS_CUDA(CIRCREF,CIRCEXP,QC, 
     &                         LCIRCC, NCANREF,NCANEXP, IEXPGO,
     &                         NSHEXP,NLOAD,
     &                         NRING, NRAYS,NRAYSC,
     &                         TCALC,TFFT,TBLAS, 
     &                         MAX_POS,NUMEXP, ISHV)

            NSHIFTST = NSHEXP / NLOAD            ! # SHIFTS / EXP IMAGE  
            DO IEXPSET = 1,NLOAD                 ! LOOP OVER EXP SETS
               IEXP = IEXPLIST(IEXPGO)+IEXPSET-1 ! ACTUAL EXP #
               ILOC = (IEXPSET-1) * NSHIFTS + 1

               CALL CHKCUDAMAXA(NRAYS,NCANREF,NSHIFTST, QC(1,1,ILOC), 
     &                 NLOAD,MAX_POS(IEXPSET), SCALE,
     &                 IEXP, ISHIFTM,IREFM,IRAYM, CCROT,.TRUE.,
     &                 ISHV)

	       RANGNEW   = ANG_N(FLOAT(IRAYM),MODE,MAXRIN)

C              HAVE TO CHANGE ORDER OF SHIFT & ROTATION.
C              IN THIS PROGRAM IMAGE IS SHIFTED FIRST, ROTATED SECOND.
C              IN 'RT SQ' IT IS ROTATION FIRST, SHIFT SECOND.
C              THIS CODE CORRESPONDS TO 'SA P'.
	       SX     = -ISHV(1,ISHIFTM)  
	       SY     = -ISHV(2,ISHIFTM) 
	       CO     =  COS(RANGNEW * DGR_TO_RAD)
	       SO     = -SIN(RANGNEW * DGR_TO_RAD)
	       XSHSUM = SX*CO - SY*SO
	       YSHSUM = SX*SO + SY*CO

               IF (CCROT > DLIST(11,IEXP)) THEN
C                 THIS IS A BETTER REFERENCE (MAY HAVE >1 REF SET)
C                 RETURNS DLIST()
                  CALL AP_END(IEXP,IEXP,IREFM,
     &                  ANGREF(1,IREFM),REFDIR,
     &                  ANGEXP(1,IEXP),EXPDIR,ISHRANGEX,
     &                  GOTREFANG, NGOTPAR, CCROT,PEAKV,
     &                  RANGNEW,XSHSUM,YSHSUM, 
     &                  MIRRORNEW,CDUM,
     &                  NUMREF,'SH', 0,DLIST(2,IEXP))
               ENDIF

            ENDDO

            NEXPDONE = NEXPDONE + NLOAD    ! TOTAL # OF EXP IMAGES 

         ENDDO                        ! END OF: IEXPGO = 1,NUMEXP,NEXPpL

C        DEALLOCATE GPU DEVICE MEMORY
         CALL CHKRINGS_CUDA_ALLOC(.FALSE.,NDUM,NDUM,NDUM,NDUM, IRTFLG)

C        DEALLOCATE CPU LOOP LEVEL ARRAYS
         IF (ALLOCATED(CIRCREF))  DEALLOCATE(CIRCREF)
         IF (ALLOCATED(CIRCEXP))  DEALLOCATE(CIRCEXP)
         IF (ALLOCATED(QC))       DEALLOCATE(QC)

      ENDDO  ! END OF: DO WHILE ()   LOOP OVER ALL REF. IMAGES

      IF (LUNDOC .GT. 0) THEN
         DO IEXP = 1,NUMEXP     ! LOOP OVER ALL EXP IMAGES
C           SAVE IN ALIGNMENT DOC FILE
C           <,<,<, REF#,IMG#,INPLANE<, SX,SY,NPROJ, <DIF,CCROT,IN
            CALL LUNDOCWRTDAT(LUNDOC,IEXP,DLIST(2,IEXP),
     &                        NLISTMAX,IRTFLG)

            CALL AP_STAT_ADD(NGOTPAR,CCROT,DLIST(11,IEXP),
     &                       ANGDIFTHR,ANGEXP(8,IEXP),
     &                       CCROTAVG,IBIGANGDIF,ANGDIFAVG,IMPROVCCROT,
     &                       CCROTIMPROV,IWORSECCROT,CCROTWORSE)
         ENDDO

C        SAVE CCROT & ANGULAR DISPLACEMENT STATISTICS
         CALL AP_STAT(NUMEXP,ANGDIFTHR,IBIGANGDIF,
     &               ANGDIFAVG, CCROTAVG,
     &               IMPROVCCROT,CCROTIMPROV,
     &               IWORSECCROT,CCROTWORSE,
     &               NBORDER,NSUBPIX,LUNDOC)
      ENDIF

      WRITE(6,99) ' ' 
      WRITE(6,99) 'CCROT AVG:',CCROTAVG 


C     OUTPUT TIMING STATISTICS
      TALL  = TLOAD + TCALC + TFFT + TBLAS
      WRITE(6,*)  '  '
      WRITE(6,99) 'CUDA LOAD:',TLOAD,' secs.'
      WRITE(6,99) 'CUDA CALC:',TCALC,' secs.'
      WRITE(6,99) 'CUDA FFT: ',TFFT, ' secs.'
      WRITE(6,99) 'CUDA BLAS:',TBLAS,' secs.'
      WRITE(6,99) 'CUDA ALL :',TALL, ' secs.'
 
9999  CONTINUE

C     DEALLOCATE  REMAINING ARRAYS
      IF (ALLOCATED(CIRCREF))    DEALLOCATE(CIRCREF)
      IF (ALLOCATED(CIRCEXP))    DEALLOCATE(CIRCEXP)
      IF (ALLOCATED(QC))         DEALLOCATE(QC)
      IF (ALLOCATED(MAX_POS))    DEALLOCATE(MAX_POS)

      IF (ALLOCATED(DLIST))      DEALLOCATE(DLIST)
      IF (ALLOCATED(REFDIR))     DEALLOCATE(REFDIR)
      IF (ALLOCATED(ANGEXP))     DEALLOCATE(ANGEXP)
      IF (ALLOCATED(ANGREF))     DEALLOCATE(ANGREF)

      IF (ALLOCATED(ISHV))       DEALLOCATE(ISHV)
      IF (ALLOCATED(NLOCS))      DEALLOCATE(NLOCS)
      IF (ALLOCATED(COEFFS))     DEALLOCATE(COEFFS)
      IF (ALLOCATED(IXY))        DEALLOCATE(IXY)
      END

#endif











