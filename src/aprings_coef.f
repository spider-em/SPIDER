
C++*********************************************************************
C
C APRINGS_*_COEF    USED CMLIMIT                  AUG 00 ARDEAN LEITH
C                   ADDED REF_CIRC FILE           APR 01 ARDEAN LEITH
C                   NORMASS -> NORMAS             OCT 01 ARDEAN LEITH
C                   PROMPTS                       JAN 02 ARDEAN LEITH
C                   AUTO REF RINGS FILE           DEC 04 ARDEAN LEITH
C                   SPIDER REF RINGS FILE         FEB 05 ARDEAN LEITH
C                   INULL  = lnblnkn(SCRFILE)     APR 05 ARDEAN LEITH
C                   REWRITE FOR SPEED             MAR 08 ARDEAN LEITH
C                   BCAST_MPI                     NOV 08 ARDEAN LEITH
C                   OUTPUT TO NOUT                AUG 09 ARDEAN LEITH
C                   ADAPTED FOR COEFFS            MAY 10 ARDEAN LEITH
C                   ADAPTED FOR TRANS             JUN 10 ARDEAN LEITH
C                   AP_GETDATA                    DEC 11 ARDEAN LEITH
C                   TYPET                         MAY 12 ARDEAN LEITH
C                   OUT OF BOUNDS TRAP ENLARGED   FEB 2013 ARDEAN LEITH
C **********************************************************************
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* Authors: J. Frank & A. Leith                                        *
C=* Copyright 1985-2013  Health Research Inc.                          *
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
C=*                                                                    *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program. If not, see <http://www.gnu.org/licenses> *
C=*                                                                    *
C **********************************************************************
C
C  APRINGS_NEW_COEF(ILIST,NUMIMG, NX,NY, 
C                    NRING,LCIRC,NUMR, MODE,FFTW_PLANS, 
C                    FILPAT,LUNIMG, CIRC,CIRC_IN_CORE,
C                    LUNRING,SCRFILE,WEIGHT,TRANS,IRTFLG)
C  
C  PARAMETERS: ILIST          LIST OF IMAGE FILE NUMBERS        (INPUT)
C              NIMREF         NO. OF IMAGES                     (INPUT)
C              NX,NY      IMAGE DIMENSIONS                  (INPUT)
C              FILPAT         IMAGE SERIES FILE TEMPLATE        (INPUT)
C              LUNIMG         IMAGE FILE IO UNIT                (INPUT)
C              CIRC           OUTPUT ARRAY                      (OUTPUT)
C
C              CIRC_IN_CORE   NO OUTPUT ARRAY FLAG              (OUTPUT)
C              LUNRING        RINGS FILE IO UNIT                (INPUT)
C              SCRFILE        RINGS FILE                        (INPUT)
C              WEIGHT         WEIGHT THE FFT'D RINGS            (INPUT)
C              WEIGHT         TRANSFORM RING RAY ORDER          (INPUT)
C              IRTFLG         ERROR FLAG                        (OUTPUT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

cpgi$g opt=O3

        SUBROUTINE APRINGS_NEW_COEF(ILIST,NUMIMG,  NX,NY,
     &                         NRING,LCIRC,NUMR, NLOCS,NRAYSC,
     &                         COEFFS,IXY,
     &                         MODE,FFTW_PLANS,
     &                         FILPAT,LUNIMG,CIRC,CIRC_IN_CORE,
     &                         LUNRING,SCRFILE, WEIGHT,TRANS,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        INTEGER,   INTENT(IN)    :: ILIST(NUMIMG) 
        INTEGER,   INTENT(IN)    :: NUMIMG,NX,NY 
        INTEGER,   INTENT(IN)    :: NUMR(3,NRING)
        INTEGER,   INTENT(IN)    :: NLOCS(2,NRAYSC+1)
        INTEGER,   INTENT(IN)    :: NRING,LCIRC,NRAYSC 
        REAL,      INTENT(OUT)   :: COEFFS(6,LCIRC)
        INTEGER,   INTENT(OUT)   :: IXY(2,LCIRC)
        CHARACTER(LEN=1)         :: MODE
        INTEGER*8, INTENT(IN)    :: FFTW_PLANS(*) ! STRUCTURE POINTERS
        CHARACTER(LEN=*)         :: FILPAT,SCRFILE
        INTEGER,   INTENT(IN)    :: LUNIMG,LUNRING   ! CONSTANTS
        REAL,      INTENT(OUT)   :: CIRC(LCIRC,NUMIMG)
        LOGICAL,   INTENT(IN)    :: CIRC_IN_CORE
        INTEGER,   INTENT(OUT)   :: IRTFLG 
        LOGICAL,   INTENT(IN)    :: WEIGHT     ! REF WEIGHTED, EXP NOT!
        LOGICAL,   INTENT(IN)    :: TRANS      ! TRANSFORMED RINGS

        CHARACTER(LEN=3)         :: TYPET = 'FI'
        CHARACTER(LEN=MAXNAM)    :: FILNAM
        LOGICAL                  :: ISOPEN
        LOGICAL                  :: USERINGFILE

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C       SEE IF THIS FILE EXISTS, (RETURNS EX, ISOPEN, LUNOP)

C       SEE IF SCRATCH "FILE" EXISTS (MAY BE INCORE FILE)
        ILOCAT = INDEX(SCRFILE,'@')
        INULL  = lnblnkn(SCRFILE)
        IF (INULL .GT. 0 .AND. SCRFILE(1:1) .NE. '_' .AND.
     &      ILOCAT == 0) THEN
C          ADD EXTENSION TO PHYSICAL FILENAME
           CALL FILNAMANDEXT(SCRFILE,DATEXC,FILNAM,NLET,.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ELSE
           FILNAM = SCRFILE
           NLET   = lnblnk(FILNAM)
        ENDIF
#ifdef USE_MPI
        IF (MYPID <= 0) THEN
           INQUIRE(FILE=FILNAM,EXIST=USERINGFILE,OPENED=ISOPEN,
     &             NUMBER=LUNOP,IOSTAT=IRTFLG)
        ENDIF
        CALL BCAST_MPI('APRINGS','USERINGFILE',USERINGFILE,1, 'L',ICOMM)
#else
        USERINGFILE = .FALSE.
        IF (INULL .GT. 0) THEN
           CALL INQUIREIF1(LUNRING,FILNAM,TYPET,USERINGFILE,ISOPEN,
     &                     LUNOP,INLNED,IMGNUM,IRTFLG)
        ENDIF
#endif

        IF (USERINGFILE .AND. MYPID <= 0) THEN
           WRITE(NOUT,*) ' Using existing rings file: ',FILNAM(1:NLET)
        ELSEIF (INULL > 0 .AND. MYPID <= 0) THEN
           WRITE(NOUT,*) ' No existing rings file: ',FILNAM(1:NLET)
        ENDIF

        IF (CIRC_IN_CORE  .AND. .NOT. USERINGFILE) THEN !-------------
C          CALCULATE RINGS DATA AND FILL CIRC ARRAY WITH IT

C          CALCULATE RINGS DATA COEFFS (HAS OMP PARALLEL LOOP)
           CALL APRINGS_FIND_COEF(NX,NY, LCIRC, NUMR,NRING, 
     &                            MODE, COEFFS,IXY, IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          USE COEFFS TO CALCULATE RINGS AND FILL CIRC ARRAY
C          HAS OMP PARALLEL LOOPS
           CALL APRINGS_FILL_COEF(ILIST,NUMIMG, FILPAT,LUNIMG,
     &                            NX,NY, 
     &                            0,0,1,
     &                            NUMR,NRING, NLOCS,NRAYSC, 
     &                            COEFFS,IXY, FFTW_PLANS,
     &                            CIRC,LCIRC,NUMIMG,
     &                            1, NCIRCGOT,
     &                            0,WEIGHT, MODE, TRANS,IRTFLG)
   
           IF (MYPID <= 0) WRITE(NOUT,*) ' Created incore rings file'

        ELSEIF (CIRC_IN_CORE .AND. USERINGFILE) THEN !----------------
C          READ EXISTING RINGS FILE AND FILL CIRC ARRAY WITH ITS DATA

C          OPEN AN EXISTING RINGS FILE
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,SCRFILE,LUNRING,'O',IFORM,
     &                 LCIRCT,NUMIMGT,NSLICE,MAXIM,' ', .TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (LCIRCT .NE. LCIRC .OR. NUMIMGT .NE. NUMIMG) THEN
               CALL ERRT(101,'RINGS FILE HAS WRONG SIZE',NE)
               IRTFLG = 1
               RETURN
           ENDIF

C          FILL CIRC WITH EXISTING RINGS DATA FROM FILE AND RETURN.
           DO IREF=1,NUMIMG
              CALL REDLIN(LUNRING,CIRC(1,IREF),LCIRC,IREF)
           ENDDO

           IF (MYPID <= 0)  WRITE(NOUT,*)' Loaded rings file incore'

        ELSEIF (.NOT. CIRC_IN_CORE .AND. USERINGFILE) THEN !-----------
C          OPEN AN EXISTING RINGS FILE SO THAT CAN READ RINGS DATA LATER
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,SCRFILE,LUNRING,'O',IFORM,
     &                 LCIRCT,NUMIMGT,NSLICE,MAXIM,' ', .FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (LCIRCT .NE. LCIRC .OR. NUMIMGT .NE. NUMIMG) THEN
              CALL ERRT(101,'EXISTING RINGS FILE HAS WRONG SIZE',NE)
              IRTFLG = 1
              RETURN
           ENDIF

        ELSEIF (.NOT. CIRC_IN_CORE .AND. .NOT. USERINGFILE) THEN !-----
C          CALCULATE RINGS DATA AND FILL RINGS FILE WITH DATA

C          CREATE RINGS OUTPUT FILE
           NSL = 1
           CALL OPFILEC(0,.FALSE.,SCRFILE,LUNRING,'B',IFORM,
     &                  LCIRC,NUMIMG,NSL,MAXIM,' ', .FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          CALCULATE RINGS DATA COEFFS (HAS OMP PARALLEL LOOP)
           CALL APRINGS_FIND_COEF(NX,NY, LCIRC, NUMR,NRING, 
     &                            MODE, COEFFS,IXY, IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          FILL CIRC RINGS FILE WITH RINGS DATA
C          HAS OMP PARALLEL SECTION
           CALL APRINGS_FILL_COEF(ILIST,NUMIMG, FILPAT,LUNIMG,
     &                            NX,NY, 
     &                            0,0,1,
     &                            NUMR,NRING,  NLOCS,NRAYSC, 
     &                            COEFFS,IXY,  FFTW_PLANS,
     &                            CIRC,LCIRC,NUMIMG,
     &                            1, NCIRCGOT,
     &                            LUNRING,WEIGHT, MODE,TRANS,IRTFLG)

           IF (MYPID <= 0) WRITE (NOUT,92) SCRFILE
92         FORMAT ('  Filled rings file: ',A )

        ENDIF

        END

C       --------------------- APRINGS_FILL_COEF --------------------------


        SUBROUTINE APRINGS_FILL_COEF(IMGLIST,NUMIMG, IMGPAT,LUNIMG,
     &                               NX,NY,
     &                               ISHRANGEX,ISHRANGEY,ISTEP,
     &                               NUMR,NRING,  NLOCS,NRAYSC,
     &                               COEFFS,IXY,  FFTW_PLANS,
     &                               CIRC,LCIRC,NCIRC,
     &                               IGOCIRC, NCIRCGOT,
     &                               LUNRING,WEIGHT, MODE, TRANS,IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        INTEGER,  INTENT(IN)    :: IMGLIST(NUMIMG)
        INTEGER,  INTENT(IN)    :: NUMIMG 
        CHARACTER(LEN=*)        :: IMGPAT
        INTEGER,  INTENT(IN)    :: LUNIMG, NX,NY
        INTEGER,  INTENT(IN)    :: ISHRANGEX,ISHRANGEY,ISTEP
        INTEGER,  INTENT(IN)    :: NUMR(3,NRING)
        INTEGER,  INTENT(IN)    :: NRING
        INTEGER,  INTENT(IN)    :: NLOCS(2,NRAYSC+1)
        INTEGER,  INTENT(IN)    :: NRAYSC
        REAL,     INTENT(OUT)   :: COEFFS(6,LCIRC)
        INTEGER,  INTENT(OUT)   :: IXY(2,LCIRC)
        INTEGER*8,INTENT(IN)    :: FFTW_PLANS(*)  ! STRUCTURE POINTERS
        REAL,     INTENT(OUT)   :: CIRC(LCIRC,NCIRC)
        INTEGER,  INTENT(IN)    :: LCIRC,NCIRC
        INTEGER,  INTENT(IN)    :: IGOCIRC
        INTEGER,  INTENT(INOUT) :: NCIRCGOT
        INTEGER,  INTENT(IN)    :: LUNRING
        LOGICAL,  INTENT(IN)    :: WEIGHT
        CHARACTER(LEN=1)        :: MODE
        LOGICAL,  INTENT(IN)    :: TRANS      ! TRANSFORMED RINGS/RAYS
        INTEGER,  INTENT(OUT)   :: IRTFLG

        REAL, ALLOCATABLE       :: XIM(:,:),CIRCO(:)
        INTEGER                 :: ICOMM,MYPID,MPIERR,MAXRIN
        INTEGER                 :: ISHX,ISHY,IMI,MWANT
        REAL                    :: CNS2,CNR2,CNS2C,CNR2C
        REAL                    :: WR(NRING)
        REAL                    :: ADUM

        LOGICAL, PARAMETER      :: USE_OMP  = .TRUE.

        CALL SET_MPI(ICOMM,MYPID,MPIERR)   ! SETS ICOMM AND MYPID

        IF (WEIGHT) THEN
           MAXRIN = NUMR(3,NRING) - 2   ! FOR FFTW3

C          RINGWE RETURNS WR WEIGHTS
	   CALL RINGWE_NEW(WR,NUMR,NRING,MAXRIN)
           IF (MODE == 'H') WR = WR * 0.5
        ELSE
           WR(1) = 0.0
        ENDIF

C       FILL COEFFS ARRAY WITH RINGS DEFINITIONS  (HAS OMP || LOOP)
        CALL APRINGS_FIND_COEF(NX,NY, LCIRC, 
     &                         NUMR,NRING, MODE, COEFFS,IXY, IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
  
        IF (TRANS) THEN
           MWANT = NX*NY+LCIRC 
           ALLOCATE(XIM(NX,NY),
     &              CIRCO(LCIRC), STAT=IRTFLG)
        ELSE
           MWANT = NX*NY + 1
           ALLOCATE(XIM(NX,NY), CIRCO(1),STAT=IRTFLG)
        ENDIF

        IF (IRTFLG.NE.0) THEN
            CALL ERRT(46,'XIM & CIRCO',MWANT)
            RETURN
        ENDIF

C       CALCULATE CENTERS
	CNS2C = (NX / 2 + 1) 
	CNR2C = (NY / 2 + 1)

        NCIRCGOT = IGOCIRC -1

C       PREPARE CIRCULAR RINGS DATA FOR ALL IMAGES FROM IMGLIST
        DO IMI=1,NUMIMG            ! LOOP OVER ALL IMAGES
 
C          LOAD ONE IMAGE (IMI) INTO XIM
           CALL AP_GETDATA(IMGLIST,NUMIMG, 
     &                    NX,NY, NX,NY,0.0,
     &                    1,IMGPAT,LUNIMG, IMI,IMI,
     &                    .TRUE., XIM, 
     &                    .FALSE.,ADUM,ADUM,  IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          LOOP OVER ANY SHIFTS (CAN LOAD A SET OF SHIFTS)

           DO ISHY=-ISHRANGEY,ISHRANGEY,ISTEP
C             LOOP OVER SHIFTED CENTERS IN Y
 	      CNR2 = CNR2C + ISHY
 
              DO ISHX=-ISHRANGEX,ISHRANGEX,ISTEP
C                LOOP OVER SHIFTED CENTERS IN X
	         CNS2     = CNS2C + ISHX 

                 NCIRCGOT = NCIRCGOT + 1      ! CURRENT CIRC # TO FILL

C                CONVERT XIM TO POLAR RINGS, NORMALIZE UNDER MASK,
C                FFT, & WEIGHT THE RINGS
                 CALL APRINGS_ONE_COEF(XIM, NX,NY, CNS2,CNR2, 
     &                       NUMR,NRING, NLOCS,NRAYSC,
     &                       COEFFS,IXY,
     &                       USE_OMP, WR, FFTW_PLANS, TRANS,
     &                       CIRCO,LCIRC, CIRC(1,NCIRCGOT),IRTFLG)
                 IF (IRTFLG .NE. 0) GOTO 9999

                 IF (LUNRING .GT. 0) THEN
C                   SAVE CIRC IN FILE OPENED ON LUNRING
                    CALL WRTLIN(LUNRING,CIRC(1,NCIRCGOT),
     &                          LCIRC,NCIRCGOT)
                 ENDIF
              ENDDO      ! END OF: ISHX=-ISHRANGEX....
           ENDDO         ! END OF: ISHY=-ISHRANGEY....
        ENDDO            ! END OF: DO IMI=1,NUMIMG


9999    CONTINUE
C       DEALLOCATE ARRAY
        IF (ALLOCATED(XIM))   DEALLOCATE(XIM)
        IF (ALLOCATED(CIRCO)) DEALLOCATE(CIRCO)

        END


C       --------------------- APRINGS_ONE_COEF --------------------------

        SUBROUTINE APRINGS_ONE_COEF(XIM, NX,NY, CNS2,CNR2, 
     &                       NUMR,NRING, NLOCS,NRAYSC,
     &                       COEFFS,IXY,
     &                       USE_OMP, WR, FFTW_PLANS, TRANS,
     &                       CIRCT,LCIRC, CIRC, IRTFLG)

        IMPLICIT NONE

        REAL,      INTENT(IN)    :: XIM(NX,NY)
        INTEGER,   INTENT(IN)    :: NX,NY
        REAL,      INTENT(IN)    :: CNS2,CNR2
        INTEGER,   INTENT(IN)    :: LCIRC        ! CIRCS LENGTH (REAL)
        INTEGER,   INTENT(IN)    :: NUMR(3,NRING)
        INTEGER,   INTENT(IN)    :: NRING
        INTEGER,   INTENT(IN)    :: NLOCS(2,NRAYSC+1)
        INTEGER,   INTENT(IN)    :: NRAYSC
        REAL,      INTENT(IN)    :: COEFFS(6,LCIRC)
        INTEGER,   INTENT(IN)    :: IXY(2,LCIRC)
        LOGICAL,   INTENT(IN)    :: USE_OMP
        REAL,      INTENT(IN)    :: WR(NRING)
        INTEGER*8, INTENT(IN)    :: FFTW_PLANS(*) ! STRUCTURE POINTERS
        LOGICAL,   INTENT(IN)    :: TRANS         ! REFORM RINGS/RAYS
        REAL,      INTENT(OUT)   :: CIRCT(LCIRC)  ! TEMP WORK ARRAY
        REAL,      INTENT(OUT)   :: CIRC(LCIRC)
        INTEGER,   INTENT(OUT)   :: IRTFLG
 
        DOUBLE PRECISION         :: AVT,VRINV
        LOGICAL, PARAMETER       :: SPIDER_SIGN = .FALSE.

        INTEGER                  :: ISHX,ISHY
        INTEGER                  :: NSB,NSE,NRB,NRE   ,I

C       CNS2 AND CNR2 ARE PREDEFINED CENTERS
C       CALCULATE DIMENSIONS FOR NORMALIZING MASK
	NSB  = -CNS2
	NSE  =  NSB + NX - 1
	NRB  = -CNR2
	NRE  =  NRB + NY - 1

C       FIND PARAMETERS TO NORMALIZE UNDER THE MASK,  
C       GET PARAMETERS AVT & VRINV TO NORMALIZE UNDER CIRCULAR MASK
        CALL NORMASC(XIM, NSB,NSE,NRB,NRE, NUMR(1,1),NUMR(1,NRING),
     &               AVT,VRINV, USE_OMP)

        ISHX = CNS2 - (NX / 2 + 1)
        ISHY = CNR2 - (NY / 2 + 1)

        !write(6,'f7.2,2i4') (coeffs(1,i), ixy(1:2,I),i=1,30)

        IF (TRANS) THEN
C          INTERPOLATE IMAGE INTO POLAR COORDINATES & APPLY NORMALIZATION
C          CREATING: CIRC  (RADIAL RINGS) FOR THIS IMAGE POSITION

           CALL APRINGS_USE_COEF(XIM,NX,NY, ISHX,ISHY,
     &                  COEFFS,IXY, AVT,VRINV, CIRCT,LCIRC, USE_OMP)
           !call chkreal('real circt',circt,lcirc, 32006,1000, 1)

C          FOURIER TRANSFORM WHOLE SET OF CIRCT RINGS
           CALL FRNGS_NEWT(CIRCT,LCIRC, NUMR,NRING, SPIDER_SIGN,
     &                     FFTW_PLANS, USE_OMP)
           !call chkcmplx('fftd coeft',circt,lcirc, 32006,1000, 1000)
           !call chkray(1,.false., circt,lcirc, numr,nring, nlocs,nraysc)

C          REFORM RINGS BY RAY ORDER AND APPLY OPTIONAL WEIGHTING
           CALL APRINGS_TRANS_REFORM(NUMR,NRING, NLOCS,NRAYSC,  
     &                               CIRCT,LCIRC/2, CIRC,LCIRC/2, WR)
         ELSE
C          INTERPOLATE IMAGE INTO POLAR COORDINATES & APPLY NORMALIZATION
C          CREATING: CIRC (RADIAL RINGS) FOR THIS IMAGE POSITION

           CALL APRINGS_USE_COEF(XIM,NX,NY, ISHX,ISHY,
     &                      COEFFS,IXY, AVT,VRINV, CIRC,LCIRC, USE_OMP)
           !call chkreal('real circ',circ,lcirc, 32006,1000, 1)

C          FOURIER TRANSFORM WHOLE SET OF CIRC RINGS
           CALL FRNGS_NEWT(CIRC, LCIRC, NUMR,NRING, SPIDER_SIGN,
     &                     FFTW_PLANS, USE_OMP)
           !call chkcmplx('fftd coef',circ,lcirc, 32006,1000, 1000)
           !call chkray(1,.false., circ,lcirc, numr,nring, nlocs,nraysc)

           IF (WR(1) > 0.0) THEN
C             WEIGHT TRANSFORMED CIRC RINGS USING  FACTORS FROM: WR
              CALL APPLYWS_NEW(CIRC,LCIRC, NUMR,WR,NRING)
           ENDIF
        ENDIF

        !write(6,*)' --- transformed wr coeft ------------'  
        !call chkcmplx('wr coef',circ,lcirc, 32006,2000, 17)  
        !call chkray(1,.true.,   circ,lcirc, numr,nring, nlocs,nraysc) 
        !call chkring(113,.true.,circ,lcirc, numr,nring, nlocs,nraysc) 

        IRTFLG = 0

        END

C     ------------------------- APRINGS_USE_COEF ------------------------

       SUBROUTINE APRINGS_USE_COEF(XIM,NX,NY,  IXSH,IYSH,
     &                COEFFS,IXY,   AVT,VRINV,      CIRC,LCIRC, USE_OMP)

C      PURPOSE: GIVEN VALUES FOR: GRIDDED LOCATIONS: F0,F1,F2,F3,F4,FC
C              AND INTERPOLATION COEFF'S: C0,C1,C2,C3,C4,CC
C              RETURNS INTERPOLATED VALUE FOR POINT WHERE COEFFS DERIVED
C          F4
C       F2 F0 F1
C          F3 FC

        IMPLICIT NONE

        INTEGER, INTENT(IN)    :: NX,NY,LCIRC
        REAL,    INTENT(IN)    :: XIM(NX,NY)
        INTEGER, INTENT(IN)    :: IXSH,IYSH
        LOGICAL, INTENT(IN)    :: USE_OMP

        REAL,    INTENT(IN)    :: COEFFS(6,LCIRC)
        INTEGER, INTENT(IN)    :: IXY(2,LCIRC)
        REAL,    INTENT(INOUT) :: CIRC(LCIRC) 
        DOUBLE PRECISION       :: AVT,VRINV

        INTEGER                :: I,IX,IY 

        IF (USE_OMP) THEN
C          FILL THE RINGS WITH ACTUAL IMAGE VALUES

c$omp      parallel do private(i,ix,iy)
           DO  I=1,LCIRC  ! LOOP OVER ALL POINTS ON THE RING SET

              IF (IXY(1,I) == -100 ) THEN
                 CIRC(I) = 0.0   ! THIS IS UNUSED FFT PAD LOCATION
                 CYCLE
              ENDIF

              IX      = IXY(1,I) + IXSH 
              IY      = IXY(2,I) + IYSH

              CIRC(I) = (XIM(IX,  IY)   * COEFFS(1,I) +
     &                   XIM(IX+1,IY)   * COEFFS(2,I) +
     &                   XIM(IX-1,IY)   * COEFFS(3,I) +
     &                   XIM(IX,  IY-1) * COEFFS(4,I) +
     &                   XIM(IX,  IY+1) * COEFFS(5,I) +
     &                   XIM(IX+1,IY+1) * COEFFS(6,I) - AVT) * VRINV
	   ENDDO
        ELSE
C          FILL THE RINGS WITH ACTUAL IMAGE VALUES
           DO  I=1,LCIRC  ! LOOP OVER ALL POINTS ON THE RING SET

              IF (IXY(1,I) == -100 ) THEN
                 CIRC(I) = 0.0   ! THIS IS UNUSED FFT PAD LOCATION
                 CYCLE
              ENDIF

              IX      = IXY(1,I) + IXSH 
              IY      = IXY(2,I) + IYSH

              CIRC(I) = (XIM(IX,  IY)   * COEFFS(1,I) +
     &                   XIM(IX+1,IY)   * COEFFS(2,I) +
     &                   XIM(IX-1,IY)   * COEFFS(3,I) +
     &                   XIM(IX,  IY-1) * COEFFS(4,I) +
     &                   XIM(IX,  IY+1) * COEFFS(5,I) +
     &                   XIM(IX+1,IY+1) * COEFFS(6,I) - AVT) * VRINV
	   ENDDO
        ENDIF

        END

C++***************************APRINGS_FIND_COEF ************************
C
C APRINGS_FIND_COEF.F   MODIFIED FROM ALRQ_MS    MAY 2010 ARDEAN LEITH
C PURPOSE: CREATES ARRAY OF COEFFICIENTS USED TO CREATE POLAR CIRCULAR
C          RINGS. FFT WILL BE USED ON RINGS.

        SUBROUTINE APRINGS_FIND_COEF(NX,NY, LCIRC, 
     &                               NUMR,NRING, MODE, 
     &                               COEFFS,IXY, IRTFLG)

        IMPLICIT NONE

        INTEGER,     INTENT(IN)  :: NX,NY,NRING,LCIRC
        INTEGER,     INTENT(IN)  :: NUMR(3,NRING)
        CHARACTER*1, INTENT(IN)  :: MODE
        REAL,        INTENT(OUT) :: COEFFS(6,LCIRC)
        INTEGER,     INTENT(OUT) :: IXY(2,LCIRC)
        INTEGER,     INTENT(OUT) :: IRTFLG

        REAL                     :: CNS2,CNR2
        REAL                     :: X, Y 
        DOUBLE PRECISION         :: PI,DFI
        INTEGER                  :: IT,INR,IGO,NVAL,LT,NE
        INTEGER                  :: LTIGO,LTLTIGO,LTLTLTIGO,NSIM,JT
        REAL                     :: YQ,FI

        INCLUDE 'CMBLOCK.INC'


C       CNS2 AND CNR2 ARE PREDEFINED CENTERS
	CNS2  = (NX / 2 + 1) 
	CNR2  = (NY / 2 + 1)

        IRTFLG = 0

        PI     = 2 * DATAN(1.0D0)

C       FIND COEFFICIENTS FOR ALL THE RINGS OVER THE IMAGE

c$omp   parallel do private(it,inr,yq,igo,nval,lt,ltigo,ltltigo,
c$omp&                      ltltltigo,nsim,dfi,x,y,jt,fi)

        DO  IT=1,NRING 

           INR  = NUMR(1,IT)        ! RADIUS OF THE CURRENT RING
           YQ   = INR               ! FLOATING POINT RADIUS
           IGO  = NUMR(2,IT)        ! STARTING LOCATION FOR RING 

C          THE ACTUAL, POWER-OF-TWO LENGTH IS NUMR(3,I)-2, ADDITIONAL
C          TWO LOCATIONS ARE ONLY FOR THE NEW FFT.
           NVAL = NUMR(3,IT) - 2    ! LENGTH OF THIS RING

           IF (MODE == 'H')   THEN
              LT = NVAL / 2
           ELSEIF (MODE=='F') THEN
              LT = NVAL / 4
           ENDIF

           LTIGO        = LT + IGO
           LTLTIGO      = LT + LT + IGO
           LTLTLTIGO    = LT + LT + LT + IGO

           NSIM         = LT - 1
           DFI          = PI / (NSIM+1)

C          TO AVOID SLOW BOUNDARY TESTS IN QUADRI_FAST, PUT THEM HERE
           X            = CNS2
           Y            = INR + CNR2
           IF (X < 2.0 .OR. X > (FLOAT(NX)-1.0) .OR. 
     &         Y < 2.0 .OR. Y > (FLOAT(NY)-1.0) ) THEN
               WRITE(NOUT,*) 'For image size1: ',NX,NY,INR
               WRITE(NOUT,90) X,Y
90             FORMAT('  FOR LOCATION: (',F7.1,',',F7.1,')')
               CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
C              RETURN NOT POSSIBLE WITH INTEL PARALLEL COMPILER
               IRTFLG = 1
#ifndef __INTEL_COMPILER
#ifndef SP_GFORTRAN
#ifndef SP_IFC
               EXIT
#endif
#endif
#endif
          ENDIF

           CALL QUADRI_COEFFS(X,Y, COEFFS(1,IGO),IXY(1,IGO))

C          TO AVOID SLOW BOUNDARY TESTS IN QUADRI_FAST, PUT THEM HERE
           X = INR + CNS2
           Y =     + CNR2
           IF (X < 2.0 .OR. X > (FLOAT(NX)-1.0) .OR. 
     &         Y < 2.0 .OR. Y > (FLOAT(NY)-1.0) ) THEN
               WRITE(NOUT,*) 'For image size2: ',NX,NY,INR
               WRITE(NOUT,90) X,Y
               CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
C              RETURN NOT POSSIBLE WITH INTEL PARALLEL COMPILER
               IRTFLG = 1
#ifndef __INTEL_COMPILER
#ifndef SP_GFORTRAN
#ifndef SP_IFC
               EXIT
#endif
#endif
#endif
           ENDIF

           CALL QUADRI_COEFFS(X,Y, COEFFS(1,LTIGO),IXY(1,LTIGO))

           IF (MODE == 'F')  THEN
C             FILL OTHER HALF OF RING

C             TO AVOID SLOW BOUNDARY TESTS IN QUADRI_FAST, PUT THEM HERE
              X = 0.0  + CNS2
              Y = -INR + CNR2
              IF (X < 2.0 .OR. X > (FLOAT(NX)-1.0) .OR. 
     &            Y < 2.0 .OR. Y > (FLOAT(NY)-1.0) ) THEN
                  WRITE(NOUT,*) 'For image size3: ',NX,NY,INR
                  WRITE(NOUT,90) X,Y
                  CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
C                 RETURN NOT POSSIBLE WITH INTEL PARALLEL COMPILER
                  IRTFLG = 1
#ifndef __INTEL_COMPILER
#ifndef SP_GFORTRAN
#ifndef SP_IFC
               EXIT
#endif
#endif
#endif
              ENDIF

              CALL QUADRI_COEFFS(X, Y, COEFFS(1,LTLTIGO),IXY(1,LTLTIGO))

C             TO AVOID SLOW BOUNDARY TESTS IN QUADRI_FAST, PUT THEM HERE
              X = -INR + CNS2
              Y =  0.0 + CNR2
              IF (X < 2.0 .OR. X > (FLOAT(NX)-1.0) .OR. 
     &            Y < 2.0 .OR. Y > (FLOAT(NY)-1.0) ) THEN
                  WRITE(NOUT,*) 'For image size4: ',NX,NY,INR
                  WRITE(NOUT,90) X,Y
                  CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
C                 RETURN NOT POSSIBLE WITH INTEL PARALLEL COMPILER
                  IRTFLG = 1
#ifndef __INTEL_COMPILER
#ifndef SP_GFORTRAN
#ifndef SP_IFC
               EXIT
#endif
#endif
#endif
              ENDIF

              CALL QUADRI_COEFFS(X,Y,COEFFS(1,LTLTLTIGO),
     &                                  IXY(1,LTLTLTIGO))
           ENDIF

           DO JT=1,NSIM     ! LOOP NSIM TIMES TO FILL RING
              FI = DFI * JT
              X  = SIN(FI) * YQ
              Y  = COS(FI) * YQ

              CALL QUADRI_COEFFS(X+CNS2, Y+CNR2,
     &                           COEFFS(1,JT+IGO),  IXY(1,JT+IGO))
              CALL QUADRI_COEFFS(Y+CNS2, -X+CNR2,
     &                           COEFFS(1,JT+LTIGO),IXY(1,JT+LTIGO))
 
              IF (MODE == 'F')  THEN
C                FILL OTHER HALF OF RING
                 CALL QUADRI_COEFFS(-X+CNS2, -Y+CNR2,
     &                       COEFFS(1,JT+LTLTIGO), IXY(1,JT+LTLTIGO))
                 CALL QUADRI_COEFFS(-Y+CNS2, X+CNR2,
     &                       COEFFS(1,JT+LTLTLTIGO),IXY(1,JT+LTLTLTIGO))

             ENDIF
	   ENDDO

 	ENDDO

        END

C     ------------------------- QUADRI_COEFFS ------------------------

      SUBROUTINE QUADRI_COEFFS(X, Y, COEFFS, IXY) 

C     PURPOSE: GIVEN A VALUE FOR: X,Y THIS SUBROUTINE RETURNS
C              COEFF FOR 'QUADRATIC INTERPOLATION': C0,C1,C2,C3,C4,CC
C              ALSO SETS: IXY CONTAINING VALUES OF: X & Y

      IMPLICIT NONE

      REAL, INTENT(IN)     :: X,Y
      REAL, INTENT(OUT)    :: COEFFS(0:5)
      INTEGER, INTENT(OUT) :: IXY(2)

      INTEGER              :: IX,IY
      REAL                 :: DX0,DY0,DXB,DYB

      IX        = IFIX(X)       ! TRUNCATES
      IY        = IFIX(Y)

      DX0       = X - IX        ! DIFFERENCE X (NON INTEGER PART OF X)
      DY0       = Y - IY        ! DIFFERENCE Y (NON INTEGER PART OF Y)

      DXB       = DX0 - 1       ! <= 0
      DYB       = DY0 - 1       ! <= 0

      COEFFS(0) = 1 -DX0 -DY0 -DX0*DXB -DY0*DYB +DX0*DY0
      COEFFS(5) = DX0*DY0
      COEFFS(1) = DX0*(1 +DXB*0.5 -DY0)
      COEFFS(2) = DX0*0.5*DXB
      COEFFS(3) = DY0*0.5*DYB
      COEFFS(4) = DY0*(1 +DYB*0.5 -DX0)

      IXY(1)    = IX
      IXY(2)    = IY

#ifdef NEVER
      if (ix < 1 ) then
         write(6,*)'ix,x:',ix,x
         call errt(102,'Bad ix in coeffs',ix)
         stop
      endif
      if (iy < 1 ) then
         write(6,*)'iy,y:',iy,y
         call errt(102,'Bad iy in coeffs',iy)
         stop
      endif

      !WRITE(6,90) IXY,COEFFS
   90 FORMAT('  (',I5,',',I5,') :',7F8.3)
#endif

      END

