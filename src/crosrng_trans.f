
C *********************************************************************
C
C CROSRNG_TRANS.F  FROM CROSRNG_E                  APR 10 ARDEAN LEITH
C 
C **********************************************************************
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
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
C   CROSRNG_TRANS(CIRCR,CIRCE,LCIRCD2, NLOCS,NRAYSC,
C                 USE_OMP,FFTW3PLAN, USE_UN,USE_MIR,                         
C                 ISMIRRORED,QMAX,POS_MAX)
C
C PURPOSE: CROSS CORRELATION OF RADIAL RINGS FOR USE IN ROTATIONAL
C          ALIGNMENT. CAN CHECKS ONLY UN-MIRRORED & MIR. POSITION
C          USES TRANSFORMED RAY ORDER RING DATA!
C
C PARAMETERS:
C    CIRCR      - FT OF RINGS MULTIPLIED BY WEIGHTS           (SENT)
C    CIRCE      - FT OF RINGS MULTIPLIED BY WEIGHTS           (SENT)
C    LCIRCD2    - SIZE OF CIRCS ARRAYS                        (SENT)
C    NLOCS      - RAY DIMENSIONS                              (SENT)
C    NRAYSC     - NUMBER OF RAYS                              (SENT)
C    USE_OMP    - USE || OMP                                  (SENT)
C    FFTW3PLAN  - PLAN FOR REVERSE FFT OF RING                (SENT)
C    USE_UN     - USE UN MIRRORED                             (SENT)
C    USE_MIR    - USE MIRRORED                                (SENT)
C    ISMIRRORED -                                             (RETURNED)
C    QMAX       - CC MAX                                      (RETURNED)
C    POS_MAX    - POSITION OF CC MAX                          (RETURNED)
C
C    U = 11*21+ 12*22  (e.g.)
C
C    ABNORMAL LOCATION   <MAXRIN       NVAL=MAXRIN    (U = USUAL)
C     E.g. NVAL:            256             512
C        1:                 11x21          11x21         (same)
C        2:                   0            12x22 
C        NVAL-1               U              U           (same)
C        NVAL                 U              U (maxrin)  (same)
C        NVAL+1             12x22            0 (unused)
C        MAXRIN               0 (unused)     U               
C        MAXRIN + 1:          0              0           (same) (unused) 
C        MAXRIN + 2:          0              0           (same) (unused)    
C        summation          0..257         0..514
C
C  NOTES:     COMPLEX CONJUGATE OF COMPLEX NUMBER = a - bi
C             USUAL COMPLEX MULTIPLICATION: (ac-bd),(ad+bc)
C             OUR   COMPLEX MULTIPLICATION: (ac+bd),(-ad+bc) 
C             (If loaded by fmrs can use usual method)
C             SPIDER_SIGN flag when data originally transformed))
C
C    Typical radii:  32, 64, 128, 256, 512  (only)
c  ip:            -9    lcirc:          23168
C
C  NOTE: JUN 2010 NO IMPROVEMENT USING REAL VARIABLES INSTEAD OF
C        COMPLEX IN CALCULATION
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE CROSRNG_TRANS(CIRCR,CIRCE,LCIRCD2, NLOCS,NRAYSC,
     &                           USE_OMP,FFTW3PLAN,
     &                           USE_UN,USE_MIR,                         
     &                           ISMIRRORED,QMAX,POS_MAX)

        IMPLICIT NONE

        COMPLEX,          INTENT(IN)  :: CIRCE(LCIRCD2), CIRCR(LCIRCD2)
        INTEGER,          INTENT(IN)  :: LCIRCD2
        INTEGER,          INTENT(IN)  :: NLOCS(2,NRAYSC+1)
        INTEGER,          INTENT(IN)  :: NRAYSC
        LOGICAL,          INTENT(IN)  :: USE_OMP
        INTEGER*8,        INTENT(IN)  :: FFTW3PLAN  ! STRUCTURE POINTER
        LOGICAL,          INTENT(IN)  :: USE_UN,USE_MIR
        LOGICAL,          INTENT(OUT) :: ISMIRRORED
        DOUBLE PRECISION, INTENT(OUT) :: QMAX
        REAL,             INTENT(OUT) :: POS_MAX

C       AUTOMATIC ARRAYS
        COMPLEX                       :: QU(NRAYSC),QM(NRAYSC)
        COMPLEX                       :: QUT,QMT

        DOUBLE PRECISION              :: QMAXU,QMAXM
        REAL                          :: POS_MAXU,POS_MAXM
        INTEGER                       :: IRAY,NVAL,JGO,J

        DO IRAY=1,NRAYSC               ! LOOP OVER ALL FFT'D RAYS 

           NVAL  = NLOCS(1,IRAY)       ! # PTS (#RINGS) ON RAY WITH FFT PAD
           JGO   = NLOCS(2,IRAY)       ! INDEX OF RAY START

           !write(6,*) ' iray,nval,jgo:', iray,nval,jgo,nraysc

           QUT   = 0.0              
           QMT   = 0.0              

           IF (USE_UN .AND. USE_MIR) THEN
              DO J= JGO,JGO+NVAL-1
C                NON MIRRORED RING SET 
                 QUT = QUT +       CIRCR(J) * CONJG(CIRCE(J))

C                MIRRORED RING SET 
                 QMT = QMT + CONJG(CIRCR(J)) * CONJG(CIRCE(J))
              ENDDO

              QU(IRAY) = QUT
              QM(IRAY) = QMT

           ELSEIF (USE_UN) THEN
              DO J= JGO,JGO+NVAL-1
C                NON MIRRORED RING SET 
                 QUT = QUT +       CIRCR(J) * CONJG(CIRCE(J))
              ENDDO

              QU(IRAY) = QUT

           ELSEIF (USE_MIR) THEN
              DO J= JGO,JGO+NVAL-1
C                MIRRORED RING SET 
                 QMT = QMT + CONJG(CIRCR(J)) * CONJG(CIRCE(J))
              ENDDO

              QM(IRAY) = QMT

           ENDIF
        ENDDO

        QMAXU = 0.0
        QMAXM = 0.0

        IF (USE_UN) THEN
C          FOR UN-MIRRORED
           CALL CROSRNG_COM_N(QU,NRAYSC*2-2,FFTW3PLAN,QMAXU,POS_MAXU)
        ENDIF

        IF (USE_MIR) THEN
C          FOR MIRRORED
           CALL CROSRNG_COM_N(QM,NRAYSC*2-2,FFTW3PLAN,QMAXM,POS_MAXM)
        ENDIF

        IF (QMAXM .GT. QMAXU) THEN
            ISMIRRORED = .TRUE.
            QMAX       = QMAXM
            POS_MAX    = POS_MAXM
        ELSE
            ISMIRRORED = .FALSE.
            QMAX       = QMAXU
            POS_MAX    = POS_MAXU
        ENDIF
       
        END



C *************************** CROSRNG_TRANS_NOC  **********************

        SUBROUTINE CROSRNG_TRANS_NOC(CIRCR,CIRCE,LCIRCD2, NLOCS,NRAYSC,
     &                           USE_OMP,FFTW3PLAN,
     &                           USE_UN,USE_MIR,                         
     &                           ISMIRRORED,QMAX,POS_MAX)

        IMPLICIT NONE

        REAL,             INTENT(IN)  :: CIRCR(LCIRCD2*2)
        REAL,             INTENT(IN)  :: CIRCE(LCIRCD2*2)

        INTEGER,          INTENT(IN)  :: LCIRCD2
        INTEGER,          INTENT(IN)  :: NLOCS(2,NRAYSC+1)
        INTEGER,          INTENT(IN)  :: NRAYSC
        LOGICAL,          INTENT(IN)  :: USE_OMP
        INTEGER*8,        INTENT(IN)  :: FFTW3PLAN  ! STRUCTURE POINTER
        LOGICAL,          INTENT(IN)  :: USE_UN,USE_MIR
        LOGICAL,          INTENT(OUT) :: ISMIRRORED
        DOUBLE PRECISION, INTENT(OUT) :: QMAX
        REAL,             INTENT(OUT) :: POS_MAX

C       AUTOMATIC ARRAYS
        REAL                          :: QU(NRAYSC*2),QM(NRAYSC*2)

        REAL                          :: FUT1,FUT2,FMT1,FMT2

        DOUBLE PRECISION              :: QMAXU,QMAXM
        REAL                          :: POS_MAXU,POS_MAXM
        INTEGER                       :: IRAY,NVAL,JGO,J

        REAL                          :: C1,C2,D1,D2

        DO IRAY=1,NRAYSC              ! LOOP OVER ALL FFT'D RAYS 

           NVAL  = NLOCS(1,IRAY)*2    ! # PTS (#RINGS) ON RAY INCL FFT PAD
           JGO   = NLOCS(2,IRAY)*2-1  ! INDEX OF RAY START

           !write(6,*) ' iray,nval,jgo:', iray,nval,jgo,nraysc

           FUT1   = 0.0              
           FUT2   = 0.0              
           FMT1   = 0.0              
           FMT2   = 0.0              

           IF (USE_UN .AND. USE_MIR) THEN
              DO J= JGO,JGO+NVAL-1,2
C                NON MIRRORED RING SET 
    	         C1      = CIRCR(J)
 	         C2      = CIRCR(J+1)
                 D1      = CIRCE(J)
                 D2      = CIRCE(J+1)

 	         FUT1    = FUT1 + C1 * D1 + C2 * D2
	         FUT2    = FUT2 - C1 * D2 + C2 * D1

 	         FMT1    = FMT1 + C1 * D1 - C2 * D2
	         FMT2    = FMT2 - C1 * D2 - C2 * D1
              ENDDO

              QU(IRAY*2-1) = FUT1
              QU(IRAY*2)   = FUT2
              QM(IRAY*2-1) = FMT1
              QM(IRAY*2)   = FMT2

           ELSEIF (USE_UN) THEN
              DO J= JGO,JGO+NVAL-1,2
C                NON MIRRORED RING SET 
    	         C1      = CIRCR(J)
 	         C2      = CIRCR(J+1)
                 D1      = CIRCE(J)
                 D2      = CIRCE(J+1)

 	         FUT1    = FUT1 + C1 * D1 + C2 * D2
	         FUT2    = FUT2 - C1 * D2 + C2 * D1
              ENDDO

              QU(IRAY*2-1) = FUT1
              QU(IRAY*2)   = FUT2

           ELSEIF (USE_MIR) THEN
              DO J= JGO,JGO+NVAL-1,2
C                MIRRORED RING SET 
    	         C1      = CIRCR(J)
 	         C2      = CIRCR(J+1)
                 D1      = CIRCE(J)
                 D2      = CIRCE(J+1)

 	         FMT1    = FMT1 + C1 * D1 - C2 * D2
	         FMT2    = FMT2 - C1 * D2 - C2 * D1
              ENDDO

              QM(IRAY*2-1) = FMT1
              QM(IRAY*2)   = FMT2

           ENDIF
        ENDDO

        QMAXU = 0.0
        QMAXM = 0.0

        IF (USE_UN) THEN
C          FOR UN-MIRRORED
           CALL CROSRNG_COM_N(QU,NRAYSC*2-2,FFTW3PLAN,QMAXU,POS_MAXU)
        ENDIF

        IF (USE_MIR) THEN
C          FOR MIRRORED
           CALL CROSRNG_COM_N(QM,NRAYSC*2-2,FFTW3PLAN,QMAXM,POS_MAXM)
        ENDIF

        IF (QMAXM .GT. QMAXU) THEN
            ISMIRRORED = .TRUE.
            QMAX       = QMAXM
            POS_MAX    = POS_MAXM
        ELSE
            ISMIRRORED = .FALSE.
            QMAX       = QMAXU
            POS_MAX    = POS_MAXU
        ENDIF
           

        END
















