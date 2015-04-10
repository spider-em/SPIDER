

C++*********************************************************************
C
C  CROSRNG_2.F   
C              TRAP FOR COMPILER ERROR ON ALTIX   FEB 2005 ARDEAN LEITH
C              REWRITE USING CROSRNG_COM          FEB 2008 ARDEAN LEITH
C              REWRITE USING FFTW3 RINGS          FEB 2008 ARDEAN LEITH
C              REMOVED TT                         JUN 2010 ARDEAN LEITH
C              CONSOLIDATES CROSRNG_E & _M        JUN 2010 ARDEAN LEITH
C              CROSRNG_2C                         OCT 2010 ARDEAN LEITH
C              CROSRNG_2C  SIGN BUG  CONJG        APR 2011 ARDEAN LEITH
C              CROSRNG_2CQ SIGN BUG  CONJG        APR 2011 ARDEAN LEITH
C              QMAXU = -HUGE(QMAXU)               JUL 2011 ARDEAN LEITH
C              IRAY1,IRAY2                        NOV 2011 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
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
C=* merchantability or fitness for a particular purpose.  See the GNU  *
C=* General Public License for more details.                           *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program. If not, see <http://www.gnu.org/licenses> *
C=*                                                                    *
C **********************************************************************
C
C CROSRNG_2(CIRCR,CIRCE,LCIRC, NRING,MAXRIN,NUMR,
C           USE_OMP,FFTW3PLAN,   USE_UN,USE_MIR,                         
C           ISMIRRORED,QMAX,POS_MAX)
C
C PURPOSE: CROSS CORRELATION OF RADIAL RINGS FOR USE IN ROTATIONAL
C          ALIGNMENT.  CHECKS BOTH STRAIGHT & MIRRORED POSITIONS
C          USES NUMR TABLE FOR MAPPING INTO Q ARRAY 
C          USES SIMPLIFIED LOGIC FOR BOUNDARY VALUES, FLOATING PT. ARITH.
C
C PARAMETERS:
C    CIRCR      - FT OF RINGS MULTIPLIED BY WEIGHTS           (SENT)
C    CIRCE      - FT OF RINGS MULTIPLIED BY WEIGHTS           (SENT)
C    LCIRC      - SIZE OF CIRCS ARRAYS                        (SENT)
C    NRING      - NUMBER OF RINGS                             (SENT)
C    MAXRIN     - LONGEST RING                                (SENT)
C    NUMR       - RING LOCATION POINTERS                      (SENT)
C    IRAY1,2    - USE THIS RAY RANGE                          (SENT)
C    USE_OMP    - USE || OMP                                  (SENT)
C    FFTW3PLAN  - PLAN FOR REVERSE FFT OF RING                (SENT)
C    USE_UN     - USE UN MIRRORED                             (SENT)
C    USE_MIR    - USE MIRRORED                                (SENT)
C    ISMIRRORED -                                             (RETURNED)
C    QMAX       - CC MAX                                      (RETURNED)
C    POS_MAX    - POSITION OF CC MAX                          (RETURNED)
C
C  NOTES: AUG 04 ATTEMPTED SPEEDUP USING 
C       PREMULTIPLY  ARRAYS ie( CIRC12 = CIRC1 * CIRC2) much slower
C       VARIOUS  OTHER ATTEMPTS  FAILED TO YIELD IMPROVEMENT
C       THIS IS A VERY IMPORTANT COMPUTE DEMAND IN ALIGNMENT & REFINE.
C       OPTIONAL LIMIT ON ANGULAR SEARCH SHOULD BE ADDED.
C       COMPLEX ARRAY ARE USUALLY SLOWER in 2010
C
C       THE UNEVEN POLAR SAMPLING IN THE FOURIER DOMAIN, IS DUE TO THE 
C       INNER RINGS HAVING HAD DENSER SAMPLING, THE VALUES AT THE 
C       HIGHER FREQUENCIES WOULD HAVE BEEN ZERO ANYWAYS, SO THERE IS 
C       NO NEED TO "INCLUDE" THEM. 

C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

cpgi$g opt=3

C       *************  CROSRNG_2CQR *************************************

C       USES COMPLEX VARIABLES INTERNALLY, USUALLY FASTER
C       USES TRANSFERRED: QU & QM (NON-AUTOMATIC = NOT ON STACK)

        SUBROUTINE CROSRNG_2CQR(CIRCR,CIRCE,LCIRCD2, 
     &                        NRING,NRAYS,NUMR, 
     &                        USE_OMP,FFTW3PLAN,
     &                        USE_UN,USE_MIR,
     &                        QU,QM,
     &                        ISMIRRORED,QMAX,POS_MAX,MAXL)

        IMPLICIT NONE

        COMPLEX,  INTENT(IN)  :: CIRCE(LCIRCD2), CIRCR(LCIRCD2)
        INTEGER,  INTENT(IN)  :: LCIRCD2
        INTEGER,  INTENT(IN)  :: NUMR(3,NRING)
        INTEGER,  INTENT(IN)  :: NRING,NRAYS
        LOGICAL,  INTENT(IN)  :: USE_OMP
        INTEGER*8,INTENT(IN)  :: FFTW3PLAN  ! STRUCTURE POINTER
        LOGICAL,  INTENT(IN)  :: USE_UN,USE_MIR
        COMPLEX,  INTENT(OUT) :: QU(NRAYS/2+1),QM(NRAYS/2+1)
        LOGICAL,  INTENT(OUT) :: ISMIRRORED
        REAL,     INTENT(OUT) :: QMAX
        REAL,     INTENT(OUT) :: POS_MAX
        INTEGER,  INTENT(OUT) :: MAXL

        REAL                  :: QMAXU,QMAXM
        REAL                  :: POS_MAXU,POS_MAXM
        INTEGER               :: I,IGOM1,NVAL,J,JC
        INTEGER               :: MAXLU,MAXLM
        INTEGER               :: NUMTH

        IF (USE_UN) THEN
C          ZERO WHOLE QU ARRAY, STRAIGHT  = CIRCR * CONJG(CIRCE)
	   QU = 0.0
        ENDIF

        IF (USE_MIR) THEN
C          ZERO QM ARRAY,  QM - MIRRORED  = CONJG(CIRCR) * CONJG(CIRCE)
	   QM = 0.0
        ENDIF

        DO I=1,NRING
	   IGOM1 = NUMR(2,I) / 2  
           NVAL  = NUMR(3,I) / 2     

	   IF (USE_UN .AND. USE_MIR)  THEN
	      DO J=1,NVAL
                 JC    = J + IGOM1

C                NON MIRRORED RING SET 
                 QU(J) = QU(J) + CONJG(CIRCR(JC)) *      (CIRCE(JC))

C                MIRRORED RING SET 
                 QM(J) = QM(J) + CONJG(CIRCR(JC)) * CONJG(CIRCE(JC))
 	      ENDDO

           ELSEIF (USE_UN) THEN
	      DO J=1,NVAL
                 JC    = J + IGOM1
                 QU(J) = QU(J) + CONJG(CIRCR(JC)) *      (CIRCE(JC))
  	      ENDDO

           ELSEIF (USE_MIR) THEN
	      DO J=1,NVAL
                 JC    = J + IGOM1
                 QM(J) = QM(J) + CONJG(CIRCR(JC)) * CONJG(CIRCE(JC))
	      ENDDO
           ENDIF
	ENDDO

        QMAXU = -HUGE(QMAXU)
        QMAXM = -HUGE(QMAXM)

        IF (USE_UN) THEN
C          FOR UN-MIRRORED
              CALL CROSRNG_COM_R(QU,NRAYS,FFTW3PLAN,
     &                           QMAXU,POS_MAXU,MAXLU)
        ENDIF

        IF (USE_MIR) THEN
C          FOR MIRRORED
              CALL CROSRNG_COM_R(QM,NRAYS,FFTW3PLAN,
     &                           QMAXM,POS_MAXM,MAXLM)
        ENDIF

        IF (QMAXM .GT. QMAXU) THEN
            ISMIRRORED = .TRUE.
            QMAX       = QMAXM
            POS_MAX    = POS_MAXM
            MAXL       = MAXLM
        ELSE
            ISMIRRORED = .FALSE.
            QMAX       = QMAXU
            POS_MAX    = POS_MAXU
            MAXL       = MAXLU
        ENDIF
           

        END


C       *************  CROSRNG_2CQRR*************************************

C       USES COMPLEX VARIABLES INTERNALLY, USUALLY FASTER
C       USES TRANSFERRED: QU & QM (NON-AUTOMATIC = NOT ON STACK)

        SUBROUTINE CROSRNG_2CQRR(CIRCR,CIRCE,LCIRCD2, 
     &                        NRING,NRAYS,NUMR, IRAY1,IRAY2,
     &                        USE_OMP,FFTW3PLAN,
     &                        USE_UN,USE_MIR,
     &                        QU,QM,
     &                        ISMIRRORED,QMAX,POS_MAX,MAXL)

        IMPLICIT NONE

        COMPLEX,  INTENT(IN)  :: CIRCE(LCIRCD2), CIRCR(LCIRCD2)
        INTEGER,  INTENT(IN)  :: LCIRCD2
        INTEGER,  INTENT(IN)  :: NUMR(3,NRING)
        INTEGER,  INTENT(IN)  :: IRAY1,IRAY2
        INTEGER,  INTENT(IN)  :: NRING,NRAYS
        LOGICAL,  INTENT(IN)  :: USE_OMP
        INTEGER*8,INTENT(IN)  :: FFTW3PLAN  ! STRUCTURE POINTER
        LOGICAL,  INTENT(IN)  :: USE_UN,USE_MIR
        COMPLEX,  INTENT(OUT) :: QU(NRAYS/2+1),QM(NRAYS/2+1)
        LOGICAL,  INTENT(OUT) :: ISMIRRORED
        REAL,     INTENT(OUT) :: QMAX
        REAL,     INTENT(OUT) :: POS_MAX
        INTEGER,  INTENT(OUT) :: MAXL

        REAL                  :: QMAXU,QMAXM
        REAL                  :: POS_MAXU,POS_MAXM
        INTEGER               :: I,IGOM1,NVAL,J,JC
        INTEGER               :: MAXLU,MAXLM
        INTEGER               :: NUMTH

        IF (USE_UN) THEN
C          ZERO WHOLE QU ARRAY, STRAIGHT  = CIRCR * CONJG(CIRCE)
	   QU = 0.0
        ENDIF

        IF (USE_MIR) THEN
C          ZERO QM ARRAY,  QM - MIRRORED  = CONJG(CIRCR) * CONJG(CIRCE)
	   QM = 0.0
        ENDIF

        DO I=1,NRING
	   IGOM1 = NUMR(2,I) / 2  
           NVAL  = NUMR(3,I) / 2     

	   IF (USE_UN .AND. USE_MIR)  THEN
	      DO J=1,NVAL
                 JC    = J + IGOM1

C                NON MIRRORED RING SET 
                 QU(J) = QU(J) + CONJG(CIRCR(JC)) *      (CIRCE(JC))

C                MIRRORED RING SET 
                 QM(J) = QM(J) + CONJG(CIRCR(JC)) * CONJG(CIRCE(JC))
 	      ENDDO

           ELSEIF (USE_UN) THEN
	      DO J=1,NVAL
                 JC    = J + IGOM1
                 QU(J) = QU(J) + CONJG(CIRCR(JC)) *      (CIRCE(JC))
  	      ENDDO

           ELSEIF (USE_MIR) THEN
	      DO J=1,NVAL
                 JC    = J + IGOM1
                 QM(J) = QM(J) + CONJG(CIRCR(JC)) * CONJG(CIRCE(JC))
	      ENDDO
           ENDIF
	ENDDO

        QMAXU = -HUGE(QMAXU)
        QMAXM = -HUGE(QMAXM)

        IF (USE_UN) THEN
C          FOR UN-MIRRORED
           IF (IRAY1 .EQ. 1 .AND. IRAY2 .EQ. NRAYS) THEN
              CALL CROSRNG_COM_R(QU,NRAYS,FFTW3PLAN,
     &                           QMAXU,POS_MAXU,MAXLU)
           ELSE
              CALL CROSRNG_COM_RR(QU,NRAYS,IRAY1,IRAY2,FFTW3PLAN,
     &                            QMAXU,POS_MAXU,MAXLU)
           ENDIF
        ENDIF

        IF (USE_MIR) THEN
C          FOR MIRRORED
           IF (IRAY1 .EQ. 1 .AND. IRAY2 .EQ. NRAYS) THEN
              CALL CROSRNG_COM_R(QM,NRAYS,FFTW3PLAN,
     &                           QMAXM,POS_MAXM,MAXLM)
           ELSE
              CALL CROSRNG_COM_RR(QM,NRAYS,IRAY1,IRAY2,FFTW3PLAN,
     &                            QMAXM,POS_MAXM,MAXLM)
           ENDIF
        ENDIF

        IF (QMAXM .GT. QMAXU) THEN
            ISMIRRORED = .TRUE.
            QMAX       = QMAXM
            POS_MAX    = POS_MAXM
            MAXL       = MAXLM
        ELSE
            ISMIRRORED = .FALSE.
            QMAX       = QMAXU
            POS_MAX    = POS_MAXU
            MAXL       = MAXLU
        ENDIF
           

        END


C********************  CROSRNG_2 *************************************

	SUBROUTINE CROSRNG_2(CIRCR,CIRCE,LCIRC, 
     &                       NRING,MAXRIN,NUMR,
     &                       USE_OMP,FFTW3PLAN,
     &                       USE_UN,USE_MIR,      
     &                       ISMIRRORED,QMAX,POS_MAX)

C       USES NUMR TABLE FOR MAPPING INTO Q ARRAY 
C       USES SIMPLIFIED LOGIC FOR BOUNDARY VALUES, FLOATING PT. ARITH.
C       USES FFT3 OR NATIVE SPIDER, 
C       CAN DO BOTH/EITHER MIR. & UN-MIR.  NO SGI FFT TT REMAINS

        IMPLICIT NONE

        REAL,             INTENT(IN)  :: CIRCR(LCIRC), CIRCE(LCIRC)
        INTEGER,          INTENT(IN)  :: LCIRC
        INTEGER,          INTENT(IN)  :: NUMR(3,NRING)
        INTEGER,          INTENT(IN)  :: NRING,MAXRIN
        LOGICAL,          INTENT(IN)  :: USE_OMP
        INTEGER*8,        INTENT(IN)  :: FFTW3PLAN  ! STRUCTURE POINTER
        LOGICAL,          INTENT(IN)  :: USE_UN,USE_MIR
        LOGICAL,          INTENT(OUT) :: ISMIRRORED
        DOUBLE PRECISION, INTENT(OUT) :: QMAX
        REAL,             INTENT(OUT) :: POS_MAX
 
        DOUBLE PRECISION              :: QMAXU, QMAXM
        INTEGER                       :: I,IGOM1,NVAL,J,JC
        REAL                          :: C1,C2,D1,D2,T1,T3,T2,T4
        REAL                          :: POS_MAXU,POS_MAXM

C       AUTOMATIC ARRAYS
        REAL                          :: QU(MAXRIN+2),QM(MAXRIN+2)

        IF (USE_UN) THEN
C          ZERO WHOLE QU ARRAY, STRAIGHT  = CIRCR * CONJG(CIRCE)
	   QU = 0.0D0
        ENDIF

        IF (USE_MIR) THEN
C          ZERO QM ARRAY,  QM - MIRRORED  = CONJG(CIRCR) * CONJG(CIRCE)
	   QM = 0.0D0
        ENDIF

        DO I=1,NRING
           IGOM1 = NUMR(2,I) - 1
           NVAL  = NUMR(3,I)     

	   IF (USE_UN .AND. USE_MIR)  THEN
	      DO J=1,NVAL,2
	         JC      = J + IGOM1

 	         C1      = CIRCR(JC)
 	         C2      = CIRCR(JC+1)
                 D1      = CIRCE(JC)
                 D2      = CIRCE(JC+1)

  	         T1      = C1 * D1
 	         T3      = C1 * D2
 	         T2      = C2 * D2
 	         T4      = C2 * D1

	         QU(J)   = QU(J)   + T1 + T2
	         QU(J+1) = QU(J+1) - T3 + T4
	         QM(J)   = QM(J)   + T1 - T2
	         QM(J+1) = QM(J+1) - T3 - T4
	      ENDDO

           ELSEIF (USE_UN) THEN

	      DO J=1,NVAL,2
	         JC      = J + IGOM1

 	         C1      = CIRCR(JC)
 	         C2      = CIRCR(JC+1)
                 D1      = CIRCE(JC)
                 D2      = CIRCE(JC+1)

                 T1      = C1 * D1   ! NO FASTER IF T1.. NOT USED
                 T3      = C1 * D2
                 T2      = C2 * D2
                 T4      = C2 * D1

                 QU(J)   = QU(J)   + T1 + T2
                 QU(J+1) = QU(J+1) - T3 + T4
 
	      ENDDO

           ELSEIF (USE_MIR) THEN

	      DO J=1,NVAL,2
	         JC     = J + IGOM1

 	         C1     = CIRCR(JC)
 	         C2     = CIRCR(JC+1)
                 D1     = CIRCE(JC)
                 D2     = CIRCE(JC+1)

                 T1     = C1 * D1
                 T3     = C1 * D2
                 T2     = C2 * D2
                 T4     = C2 * D1

                 QM(J)   = QM(J)   + T1 - T2
                 QM(J+1) = QM(J+1) - T3 - T4

	      ENDDO
           ENDIF
	ENDDO

        QMAXU = -HUGE(QMAXU)
        QMAXM = -HUGE(QMAXM)

        IF (USE_UN) THEN
C          FOR UN-MIRRORED
           CALL CROSRNG_COM_N(QU,MAXRIN,FFTW3PLAN,QMAXU,POS_MAXU)
        ENDIF

        IF (USE_MIR) THEN
C          FOR MIRRORED
           CALL CROSRNG_COM_N(QM,MAXRIN,FFTW3PLAN,QMAXM,POS_MAXM)

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


        !write(6,*) 'qmaxs:',qmaxu,qmaxm
      
	END

C********************  CROSRNG_2C *************************************

C       USES COMPLEX VARIABLES INTERNALLY, USUALLY FASTER
C       USES AUTOMATIC: QU & QM 

        SUBROUTINE CROSRNG_2C(CIRCR,CIRCE,LCIRCD2, 
     &                        NRING,MAXRIN,NUMR,
     &                        USE_OMP,FFTW3PLAN,
     &                        USE_UN,USE_MIR,                         
     &                        ISMIRRORED,QMAX,POS_MAX,MAXL)

        IMPLICIT NONE

        COMPLEX,       INTENT(IN)  :: CIRCE(LCIRCD2), CIRCR(LCIRCD2)
        INTEGER,       INTENT(IN)  :: LCIRCD2
        INTEGER,       INTENT(IN)  :: NUMR(3,NRING)
        INTEGER,       INTENT(IN)  :: NRING,MAXRIN
        LOGICAL,       INTENT(IN)  :: USE_OMP
        INTEGER*8,     INTENT(IN)  :: FFTW3PLAN  ! STRUCTURE POINTER
        LOGICAL,       INTENT(IN)  :: USE_UN,USE_MIR
        LOGICAL,       INTENT(OUT) :: ISMIRRORED
        REAL,          INTENT(OUT) :: QMAX
        REAL,          INTENT(OUT) :: POS_MAX
        INTEGER,       INTENT(OUT) :: MAXL

C       AUTOMATIC ARRAYS
        COMPLEX                    :: QU(MAXRIN/2+1),QM(MAXRIN/2+1)

        REAL                       :: QMAXU,QMAXM
        REAL                       :: POS_MAXU,POS_MAXM
        INTEGER                    :: I,IGOM1,NVAL,J,JC
        INTEGER                    :: MAXLU,MAXLM
        INTEGER                    :: NUMTH

        IF (USE_UN) THEN
C          ZERO WHOLE QU ARRAY, STRAIGHT  = CIRCR * CONJG(CIRCE)
	   QU = 0.0
        ENDIF

        IF (USE_MIR) THEN
C          ZERO QM ARRAY,  QM - MIRRORED  = CONJG(CIRCR) * CONJG(CIRCE)
	   QM = 0.0
        ENDIF

        DO I=1,NRING
	   IGOM1 = NUMR(2,I) / 2  
           NVAL  = NUMR(3,I) / 2     

	   IF (USE_UN .AND. USE_MIR)  THEN
	      DO J=1,NVAL
	         JC    = J + IGOM1

C                NON MIRRORED RING SET 
                 QU(J) = QU(J) + CONJG(CIRCR(JC)) *      (CIRCE(JC))


C                MIRRORED RING SET 
                 QM(J) = QM(J) + CONJG(CIRCR(JC)) * CONJG(CIRCE(JC))
 	      ENDDO

           ELSEIF (USE_UN) THEN
	      DO J=1,NVAL
	         JC    = J + IGOM1
                 QU(J) = QU(J) + CONJG(CIRCR(JC)) *      (CIRCE(JC))
	      ENDDO

           ELSEIF (USE_MIR) THEN
	      DO J=1,NVAL
	         JC    = J + IGOM1
                 QM(J) = QM(J) + CONJG(CIRCR(JC)) * CONJG(CIRCE(JC))
	      ENDDO
           ENDIF
	ENDDO

        QMAXU = -HUGE(QMAXU)
        QMAXM = -HUGE(QMAXM)

        IF (USE_UN) THEN
C          FOR UN-MIRRORED
           CALL CROSRNG_COM_R(QU,MAXRIN,FFTW3PLAN,QMAXU,POS_MAXU,MAXLU)
        ENDIF

        IF (USE_MIR) THEN
C          FOR MIRRORED
           CALL CROSRNG_COM_R(QM,MAXRIN,FFTW3PLAN,QMAXM,POS_MAXM,MAXLM)
        ENDIF

        IF (QMAXM .GT. QMAXU) THEN
            ISMIRRORED = .TRUE.
            QMAX       = QMAXM
            POS_MAX    = POS_MAXM
            MAXL       = MAXLM
        ELSE
            ISMIRRORED = .FALSE.
            QMAX       = QMAXU
            POS_MAX    = POS_MAXU
            MAXL       = MAXLU
        ENDIF
           

        END


C******************************** CROSRNG_SATU  ***********************
C
C  CROSRNG_SATU(CIRCEXP,CIRCREF, NUMR,NRING,LCIRC,FFTW_PLANS,  
C               NLOCS,NRAYSC, TRANS,CPLX,  USE_UN,USE_MIR,
C               ISMIRRORED, CCR, ROTMP, IRTFLG)
 
C  
C  PURPOSE: CARRY OUT  CROSSCORRELATION OF RINGS.

C  SOME PARAMETERS:
C       CIRCEXP             IMAGE IN POLAR RINGS              (INPUT)
C       CIRCREF             IMAGE IN POLAR RINGS              (INPUT)
C       IRTFLG              ERROR FLAG                        (OUTPUT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE CROSRNG_SATU(CIRCEXP,CIRCREF, NUMR,NRING,LCIRC,
     &                          FFTW_PLANS,  NLOCS,NRAYSC,
     &                          TRANS,CPLX,  USE_UN,USE_MIR,
     &                          ISMIRRORED, CCR, ROTMP, IRTFLG)

C       NOTE: OFTEN RUNS WITHIN OMP PARALLEL SECTION OF CODE!

        IMPLICIT NONE

        INTEGER                :: NRING,LCIRC
        INTEGER                :: NUMR(3,NRING) 
        INTEGER *8             :: FFTW_PLANS(*)  ! STRUCTURE POINTERS
        REAL                   :: CIRCEXP(LCIRC)
        INTEGER                :: NLOCS(2,*)
        INTEGER                :: NRAYSC

        CHARACTER (LEN=9)      :: SETRINGS
        LOGICAL                :: TRANS          ! TRANSFORMED RAYS
        LOGICAL                :: CPLX           ! COMPLEX CROSRNG
        REAL                   :: CIRCREF(LCIRC)
        LOGICAL                :: USE_UN,USE_MIR

 	REAL,   INTENT(OUT)    :: ROTMP
        LOGICAL,INTENT(OUT)    :: ISMIRRORED
        INTEGER,INTENT(OUT)    :: IRTFLG

	DOUBLE PRECISION       :: CCD
	REAL,INTENT(OUT)       :: CCR

        LOGICAL, PARAMETER     :: USE_OMP = .FALSE.

        INTEGER                :: MAXRIN,MAXL

        MAXRIN = NUMR(3,NRING) - 2 ! ACTUAL LENGTH OF LONGEST RING

C       CROSS CORRELATION OF POLAR RINGS -----------------------

        IF (TRANS) THEN
C          USING TRANSFORMED RINGS
           IF (CPLX) THEN
C             USING COMPLEX VARIABLES ('AP SHCTC')
              CALL CROSRNG_TRANS(CIRCREF,CIRCEXP, LCIRC/2,
     &                           NLOCS,NRAYSC,
     &                           USE_OMP,FFTW_PLANS(1),
     &                           USE_UN,USE_MIR,
     &                           ISMIRRORED,CCD,ROTMP)
              CCR = CCD    ! CROSRNG_* DOES CALCS IN REAL NOW, SO OK
        ELSE
C             NOT USING COMPLEX VARIABLES ('AP SHCTU')
              CALL CROSRNG_TRANS_NOC(CIRCREF,CIRCEXP, LCIRC/2,
     &                           NLOCS,NRAYSC,
     &                           USE_OMP,FFTW_PLANS(1),
     &                           USE_UN,USE_MIR,
     &                           ISMIRRORED,CCD,ROTMP)
              CCR = CCD    ! CROSRNG_* DOES CALCS IN REAL NOW, SO OK
           ENDIF

        ELSE
C          USING NON-TRANSFORMED RINGS
           IF (CPLX) THEN
C             USING COMPLEX VARIABLES ('AP SHC', 'AP SHCUC')
              CALL CROSRNG_2C(CIRCREF,CIRCEXP,LCIRC/2,
     &                        NRING, MAXRIN,NUMR, 
     &                        USE_OMP,FFTW_PLANS(1),
     &                        USE_UN,USE_MIR,
     &                        ISMIRRORED,CCR,ROTMP,MAXL)
           ELSE
C             NOT USING COMPLEX VARIABLES ('AP SHCUU')
              CALL CROSRNG_2(CIRCREF,CIRCEXP,LCIRC,
     &                       NRING, MAXRIN,NUMR, 
     &                       USE_OMP,FFTW_PLANS(1),
     &                       USE_UN,USE_MIR,
     &                       ISMIRRORED,CCD,ROTMP)
              CCR = CCD    ! CROSRNG_* DOES CALCS IN REAL NOW, SO OK
           ENDIF
        ENDIF

        IRTFLG = 0

        !write(6,'a,2f8.2') ' 0 ccd:',ccd,ccr 
        END





