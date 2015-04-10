
C++*********************************************************************
C
C CROSRNG_COM.F  MERGED CROSRMG_DS.F & CROSRNG_DS    AUG 04 ARDEAN LEITH
C                REWRITE FOR FFTW3                   MAR 08 ARDEAN LEITH
C                REMOVED TT FOR CROSRNG_EP           JUN 10 ARDEAN LEITH
C                REMOVED CROSRNG_E*                  JUN 10 ARDEAN LEITH
C                CROSRNG_COM_N PARAMETERS            JUL 10 ARDEAN LEITH
C                CROSRNG_COM_R ADDED                 OCT 10 ARDEAN LEITH
C                CROSRNG_COM_RR   IRAY1,IRA          NOV 11 ARDEAN LEITH
C                CROSRNG_COM_NEW  BAD MAXLTEST       AUG 12 ARDEAN LEITH
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012 Health Research Inc.,                         *
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
C CROSRNG_COM_N(QS,LCIRC,MAXRIN,QMAX,POS_MAX,FFTW3PLAN)
C
C PURPOSE: CROSS CORRELATION OF RADIAL RINGS FOR USE IN ROTATIONAL
C          ALIGNMENT. COMMON CODE FOR CROSRNG_2 CROSRNG_M*
C
C PARAMETERS:
C    QS      - RING FOR FFT                                (SENT)
C    LCIRC   - SIZE OF CIRCS RING ARRAYS                   (SENT)
C    MAXRIN  - LONGEST RING                                (SENT)
C    QMAX    - CC MAX                                      (RETURNED)
C    POS_MAX - POSITION OF CC MAX                          (RETURNED)
C    FFTW3PLAN PLAN FOR FFTW3 REVERSE FOURIER OF MAXRIN    (SENT)
C
C
C    Typical radii:  32, 64, 128, 256, 512  (only)
c               ip: -9    lcirc: 23168
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

cpgi$g opt=3

        SUBROUTINE CROSRNG_COM_N(QS,NRAYS,FFTW3PLAN, QMAX,POS_MAX)

C       COMMON CODE FOR REVERSE FFT AND MAX LOCATION DETERMINATION
C       USES FFT3 OR NATIVE SPIDER FFT
C       NO SGI FFT TT REMAINS 
C       RETURNS: QMAX,POS_MAX
C       CALLERS: crosrng_2.f, crosrng_trans.f


        IMPLICIT NONE

        REAL,             INTENT(INOUT) :: QS(NRAYS+2) ! 2=FFT CMPLX PAD
        INTEGER,          INTENT(IN)    :: NRAYS
        INTEGER*8,        INTENT(IN)    :: FFTW3PLAN  ! STRUCTURE POINTER
        DOUBLE PRECISION, INTENT(OUT)   :: QMAX
        REAL,             INTENT(OUT)   :: POS_MAX

C       AUTOMATIC ARRAYS
        DOUBLE PRECISION                :: T7(-3:3)
        INTEGER                         :: MAXL_ARRAY(1)

        INTEGER                         :: IRTFLG,MAXL,K,J
        REAL                            :: POS

        LOGICAL, PARAMETER              :: SPIDER_SIGN  = .FALSE.
        LOGICAL, PARAMETER              :: SPIDER_SCALE = .FALSE.
        INTEGER, PARAMETER              :: INV = -1   ! REVERSE TRANSFORM

        INTEGER, SAVE :: IDONE = 1

C       REVERSE FOURIER TRANSFORM ON QS

        CALL FMRS(QS,NRAYS,1,1,FFTW3PLAN,
     &            SPIDER_SIGN,SPIDER_SCALE,INV,IRTFLG)
        !call chkreal('real qs',qs,1026, 20,1, 1025)

C       FIND MAXIMUM AND ITS LOCATION INSIDE QS
        MAXL_ARRAY = MAXLOC(QS(1:NRAYS)) ! RETURNS ARRAY OF LENGTH: 1
        MAXL       = MAXL_ARRAY(1)       ! LOCATION OF MAXIMUM
        QMAX       = QS(MAXL)            ! MAXIMUM VALUE

        QMAX       = 1.00048 * QMAX / (NRAYS) ! HACK TO = PRE FFTW3

       
C       INTERPOLATION OVER 3+3 NEIGHBORHOOD AROUND MAXL
        DO K=-3,3
           J     = MOD(MAXL+K+NRAYS-1,NRAYS) + 1
           T7(K) = QS(J)
        ENDDO

        CALL PRB1D(T7,7,POS)
        POS_MAX = FLOAT(MAXL) + POS    ! SUB-PIXEL LOCATION
 
#ifdef DEBUG
        !write(6,*)'nrays:',nrays
        if (idone .le. 2) write(6,901)maxl,qmax
 901    format(' crosrng_com_n; maxl: ',i6,' qmax: ',f16.3)
        idone=idone+1
#endif

        END





C--************************  CROSRNG_COM ******************************

        SUBROUTINE CROSRNG_COM(QS,LCIRC,MAXRIN,QMAX,POS_MAX,TT)

C       COMMON CODE FOR REVERSE TRANSFORM AND MAX LOCATION DETERMINATION
C       FOR USE WITH SPIDER FFT ONLY  CALLED FROM: oracfmskm.f
C       RETURNS: QMAX,POS_MAX,MAXL

        DOUBLE PRECISION              :: QS(MAXRIN + 2)
        INTEGER, INTENT(IN)           :: LCIRC,MAXRIN
        DOUBLE PRECISION, INTENT(OUT) :: QMAX
        REAL, INTENT(OUT)             :: POS_MAX
        DOUBLE PRECISION, INTENT(IN)  :: TT(*)     ! UNUSED

        DOUBLE PRECISION              :: T7(-3:3)
        INTEGER                       :: MAXL_ARRRAY(1)
        
C       REVERSE FOURIER TRANSFORM ON QS   
        IP = -LOG2(MAXRIN)
        CALL FFTR_D(QS,IP)            ! SPIDER FFT

C       FIND MAXIMUM AND ITS LOCATION INSIDE QS
        MAXL_ARRRAY = MAXLOC(QS(1:MAXRIN)) ! RETURNS ARRAY OF LENGTH: 1
        MAXL        = MAXL_ARRRAY(1)       ! LOCATION OF MAXIMUM
        QMAX        = QS(MAXL)             ! MAXIMUM VALUE

        QMAX = QMAX/MAXRIN

C       INTERPOLATION OVER 3+3 NEIGHBORHOOD AROUND MAXL
        DO K=-3,3
           J     = MOD(MAXL+K+MAXRIN-1,MAXRIN) + 1
           T7(K) = QS(J)
        ENDDO

        CALL PRB1D(T7,7,POS)
        POS_MAX = FLOAT(MAXL) + POS    ! SUB-PIXEL LOCATION

        END

C--************************  CROSRNG_COM_RR ******************************

        SUBROUTINE CROSRNG_COM_RR(QS,NRAYS,IRAY1,IRAY2,
     &                           FFTW3PLAN,QMAX,POS_MAX,MAXL)

C       COMMON CODE FOR REVERSE FFT AND MAX LOCATION DETERMINATION
C       USES FFT3 
C       CAN RESTRICT RAY SEARCH RANGE 
C       RETURNS: QMAX,POS_MAX,MAXL

        IMPLICIT NONE

        REAL,        INTENT(INOUT) :: QS(NRAYS+2) ! 2=FFT CMPLX PAD
        INTEGER,     INTENT(IN)    :: NRAYS,IRAY1,IRAY2
        INTEGER*8,   INTENT(IN)    :: FFTW3PLAN  ! STRUCTURE POINTER
        REAL,        INTENT(OUT)   :: QMAX
        REAL,        INTENT(OUT)   :: POS_MAX
        INTEGER,     INTENT(OUT)   :: MAXL

C       AUTOMATIC ARRAYS
        REAL                       :: T7(-3:3)
        INTEGER                    :: MAXL_ARRAY(1)

        INTEGER                    :: IRTFLG,MAXL2,K,J
        REAL                       :: POS

        LOGICAL, PARAMETER         :: SPIDER_SIGN  = .FALSE.
        LOGICAL, PARAMETER         :: SPIDER_SCALE = .FALSE.
        INTEGER, PARAMETER         :: INV = -1   ! REVERSE TRANSFORM

C       REVERSE FOURIER TRANSFORM ON QS
        CALL FMRS(QS,NRAYS,1,1,FFTW3PLAN,
     &            SPIDER_SIGN,SPIDER_SCALE,INV,IRTFLG)
        !call chkreal('real qs',qs,1026, 20,1, 1025)

C       FIND MAXIMUM AND ITS LOCATION INSIDE QS
        IF (IRAY1 > IRAY2) THEN
C          RANGE GOES ACROSS ZERO, BREAK IT IN TWO PARTS
           MAXL_ARRAY = MAXLOC(QS(IRAY1:NRAYS))  ! RETURNS ARRAY OF LENGTH: 1
           MAXL       = IRAY1 + MAXL_ARRAY(1) -1 ! LOCATION OF MAXIMUM

           MAXL_ARRAY = MAXLOC(QS(1:IRAY2))      ! RETURNS ARRAY OF LENGTH: 1
           MAXL2      = MAXL_ARRAY(1)            ! LOCATION OF MAXIMUM

           IF (QS(MAXL2) > QS(MAXL)) MAXL = MAXL2
        ELSE   
           MAXL_ARRAY = MAXLOC(QS(IRAY1:IRAY2)) ! RETURNS ARRAY OF LENGTH: 1
           MAXL       = MAXL_ARRAY(1)           ! LOCATION OF MAXIMUM
        ENDIF

        QMAX       = QS(MAXL)                 ! MAXIMUM VALUE
        QMAX       = 1.00048 * QMAX / (NRAYS) ! HACK TO = PRE FFTW3

C       CAN NOT! USE PRB1 FOR INTERPOLATION OVER 3+3 
C       NEIGHBORHOOD AROUND MAXL SINCE IT IS NOT HIGHEST POINT!!!
        POS_MAX = FLOAT(MAXL)  
 
#ifdef DEBUGNEVER
        !if (maxl .eq. 1) write(6,*) 'bad j:',it
        write(6,*)maxl,pos
 901    format(' crosrng_com_rr; maxl,pos: ',i6,f8.2)
#endif
        END



C--************************  CROSRNG_COM_R ******************************

        SUBROUTINE CROSRNG_COM_R(QS,NRAYS,FFTW3PLAN,QMAX,POS_MAX,MAXL)

C       COMMON CODE FOR REVERSE FFT AND MAX LOCATION DETERMINATION
C       RETURNS: QMAX, POS_MAX, MAXL


        IMPLICIT NONE

        REAL,        INTENT(INOUT) :: QS(NRAYS+2) ! 2=FFT CMPLX PAD
        INTEGER,     INTENT(IN)    :: NRAYS
        INTEGER*8,   INTENT(IN)    :: FFTW3PLAN  ! STRUCTURE POINTER
        REAL,        INTENT(OUT)   :: QMAX
        REAL,        INTENT(OUT)   :: POS_MAX
        INTEGER                    :: MAXL

C       AUTOMATIC ARRAYS
        REAL                       :: T7(-3:3)
        INTEGER                    :: MAXL_ARRAY(1)

        INTEGER                    :: IRTFLG,K,J
        REAL                       :: POS

        LOGICAL, PARAMETER         :: SPIDER_SIGN  = .FALSE.
        LOGICAL, PARAMETER         :: SPIDER_SCALE = .FALSE.
        INTEGER, PARAMETER         :: INV = -1   ! REVERSE TRANSFORM

C       REVERSE FOURIER TRANSFORM ON QS
        CALL FMRS(QS,NRAYS,1,1,FFTW3PLAN,
     &            SPIDER_SIGN,SPIDER_SCALE,INV,IRTFLG)
        !call chkreal('real qs',qs,1026, 20,1, 1025)

C       FIND MAXIMUM AND ITS LOCATION INSIDE QS
        MAXL_ARRAY = MAXLOC(QS(1:NRAYS)) ! RETURNS ARRAY OF LENGTH: 1
        MAXL       = MAXL_ARRAY(1)       ! LOCATION OF MAXIMUM
        QMAX       = QS(MAXL)            ! MAXIMUM VALUE

        QMAX       = 1.00048 * QMAX / (NRAYS) ! HACK TO = PRE FFTW3

C       INTERPOLATION OVER 3+3 NEIGHBORHOOD AROUND MAXL
         DO K=-3,3
           J     = MOD(MAXL+K+NRAYS-1,NRAYS) + 1
           T7(K) = QS(J)
        ENDDO

        CALL PRB1(T7,7,POS)
        POS_MAX = FLOAT(MAXL) + POS    ! SUB-PIXEL LOCATION

        END

























