
C ++********************************************************************
C                                                                      *
C  ANISOF                                                              *
C                           STACKS SUPPORT OCT 2002 ArDean Leith       *
C                           DUMMIES DIMENSION BUG NOV 02 ArDean Leith  *
C
C  ANISOF(LUN1,LUN2,NSAM,NROW,NSLICE,IRTFLG)
C
C  PARAMETERS: LUN1,LUN2   IO UNITS                             (INPUT)
C              NSAM        X DIMENSIONS                         (INPUT)
C              NROW        Y DIMENSIONS                         (INPUT)
C              NSLICE      Z DIMENSIONS                         (INPUT)
C
C  PURPOSE: ALTER CONTRAST IN AN IMAGE/ VOLUME USING ANISOTROPIC       *
C           DIFFUSION                                                  *
C                                                                      *
C  COPYRIGHT (C) 2001 BY ACHILLES FRANGAKIS.                           *
C                                                                      *
C  PORTED FROM C CODE by A. Leith                                      *
C **********************************************************************

        SUBROUTINE ANISOF(LUN1,LUN2,NSAM,NROW,NSLICE,ITER,HT,
     &              FLAMBDA,SIGMA,IRTFLG)

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL, ALLOCATABLE, DIMENSION(:,:) :: BUF

        ALLOCATE(BUF(NSAM+2,NROW+ 2), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'ANISO, OUT...',IER)
           RETURN
        ENDIF

C       LOAD IMAGE INTO CENTRAL PORTION OF BUF 
        J = 2
        DO IREC=1,NROW
           CALL REDLIN(LUN1,BUF(2,J),NSAM,IREC)
           J = J + 1
        ENDDO

C       REPORT MINIMUM, MAXIMUM, MEAN, VARIANCE
        CALL ANALYSE(BUF,NSAM,NROW,FMINT,FMAXT,FAVT,FSIGT)
        K = 0
        IF (VERBOSE) WRITE(NOUT,91)K,FMINT,FMAXT,FAVT,FSIGT
91      FORMAT(' ITER: ',I5,'  MIN: ',G13.7,'  MAX: ',G13.7,
     &         '  AV: ',G13.7,'  SIG: ',G13.7)

        DO K=1,ITER
C          NON-LINEAR ANISOTROPIC DIFFUSION
           CALL EED(HT,NSAM,NROW,FLAMBDA,SIGMA,BUF)

C          REPORT MINIMUM, MAXIMUM, MEAN, VARIANCE
           CALL ANALYSE(BUF,NSAM,NROW,FMINT,FMAXT,FAVT,FSIGT)
           IF (VERBOSE) WRITE(NOUT,91)K,FMINT,FMAXT,FAVT,FSIGT
        ENDDO

        J = 2
        DO IREC=1,NROW
           CALL WRTLIN(LUN2,BUF(2,J),NSAM,IREC)
           J = J + 1
        ENDDO
  
9999    IF (ALLOCATED(BUF)) DEALLOCATE(BUF)
        RETURN
        END

C     -------------------------ANALYSE ----------------------

      SUBROUTINE ANALYSE(BUF,NX,NY,FMINT,FMAXT,FAVT,FSIGT)

      REAL, DIMENSION(NX+2,NY+2) :: BUF

      DOUBLE PRECISION DAV,DAV2,DTOP,FNALL,DTEMP

      NPIX   = NX * NY

      DAV    = 0.0
      DAV2   = 0.0

      FMAXT  = BUF(2,2)
      FMINT  = BUF(2,2)

      DO I = 2,NX + 1
         DO J = 2,NY+1
            B     = BUF(I,J)
            FMAXT = MAX(B,FMAXT)
            FMINT = MIN(B,FMINT)
            DAV   = DAV  + DBLE(B)
            DAV2  = DAV2 + DBLE(B) * DBLE(B)
         ENDDO
      ENDDO

      DTOP  = DAV2 - DAV * DAV / DBLE(NPIX)

      FAVT  = DAV / DBLE(NPIX)
      FSIGT = DSQRT( DTOP / DBLE(NPIX-1.0))

      END


C     ----------------------------- EED -------------------------------*/
C
C  PURPOSE:
C     EDGE-ENHANCING ANISOTROPIC DIFFUSION. 
C     EXPLICIT DISCRETIZATION.
C
C   PARAMETERS:
C     HT            TIME STEP SIZE, 0 < HT <= 0.25
C     NX            IMAGE DIMENSION IN X DIRECTION 
C     NY            IMAGE DIMENSION IN Y DIRECTION 
C     FLAMBDA       CONTRAST PARAMETER
C     SIGMA         NOISE SCALE
C     UBUF          INPUT: ORIGINAL IMAGE  OUTPUT: SMOOTHED
C   VARIABLES:
C     HX, HY              PIXEL SIZE IN X & Y DIRECTION
C     I, J                LOOP VARIABLES
C     RXX, RXY, RYY       TIME SAVERS
C     WN, WNE, WE, WSE    WEIGHTS
C     WS, WSW, WW, WNW    WEIGHTS
C     FBUF                WORK COPY OF U
C     DXX,  DXY, DYY      ENTRIES OF DIFFUSION TENSOR
       
      SUBROUTINE EED(HT, NX, NY, FLAMBDA, SIGMA, UBUF)

      REAL, DIMENSION(NX+2,NY+2) :: UBUF

C     AUTOMATIC ARRAYS
      REAL, DIMENSION(NX+2,NY+2) :: FBUF,DXX,DXY,DYY

      HX = 1.0
      HY = 1.0

C     COPY UBUF INTO FBUF AND ASSIGN DUMMY BOUNDARIES
      FBUF = UBUF
      CALL DUMMIES(FBUF, NX, NY)

C     REGULARIZE FBUF 
      IF (SIGMA .GT. 0.0) THEN
         CALL GAUSS_CONV(SIGMA,NX,NY,HX,HY,3.0,FBUF)
      ENDIF

C     CALCULATE ENTRIES OF DIFFUSION TENSOR 
      CALL DIFF_TENSOR(FBUF,FLAMBDA,NX,NY,HX,HY,DXX,DXY,DYY)

C     CALCULATE EXPLICIT NONLINEAR DIFFUSION OF UBUF 

      RXX  = HT / (2.0 * HX * HX)
      RYY  = HT / (2.0 * HY * HY)
      RXY  = HT / (4.0 * HX * HY)

      CALL  DUMMIES(DXX, NX, NY)
      CALL  DUMMIES(DXY, NX, NY)
      CALL  DUMMIES(DYY, NX, NY)

C     COPY UBUF INTO FBUF AND ASSIGN DUMMY BOUNDARIES
      FBUF = UBUF
      CALL DUMMIES(FBUF, NX, NY)

C     DIFFUSE
      DO I=2, NX+1
         DO J=2, NY+1
C           WEIGHTS

            TMP = ABS(DXY(I,J))

            WE  =   RXX * (DXX(I,J+1) + DXX(I,J) - 
     &              ABS(DXY(I,J+1))   - TMP)

            WW  =   RXX * (DXX(I,J-1) + DXX(I,J) - 
     &              ABS(DXY(I,J-1))   -  TMP)

            WS  =   RYY * (DYY(I+1,J) + DYY(I,J) - 
     &              ABS(DXY(I+1,J))   - TMP)

            WN  =   RYY * (DYY(I-1,J) + DYY(I,J) - 
     &              ABS(DXY(I-1,J))   -  TMP)

            WSE =   RXY * (ABS(DXY(I+1,J+1)) + 
     &              DXY(I+1,J+1) + TMP + DXY(I,J))

            WNW =   RXY * (ABS(DXY(I-1,J-1)) + 
     &              DXY(I-1,J-1) + TMP + DXY(I,J))

            WNE =   RXY * (ABS(DXY(I-1,J+1)) - 
     &              DXY(I-1,J+1) + TMP - DXY(I,J))

            WSW =   RXY * (ABS(DXY(I+1,J-1)) - 
     &              DXY(I+1,J-1) + TMP - DXY(I,J))

C           MODIFY WEIGHTS TO PREVENT FLUX ACROSS BOUNDARIES
            IF (I .EQ. 2) THEN
C              SET WESTERN WEIGHTS ZERO
               WSW = 0.0
               WW  = 0.0
               WNW = 0.0

            ELSEIF (I.EQ.(NX+1)) THEN 
C              SET EASTERN WEIGHTS ZERO
               WNE = 0.0
               WE  = 0.0
               WSE = 0.0
            ENDIF

            IF (J.EQ.2) THEN 
C              SET NORTHERN WEIGHTS ZERO
               WNW = 0.0
               WN  = 0.0
               WNE = 0.0

            ELSEIF (J.EQ.(NY+1)) THEN 
C              SET SOUTHERN WEIGHTS ZERO
               WSE = 0.0
               WS  = 0.0 
               WSW = 0.0
            ENDIF

C           EVOLUTION
            UBUF(I,J) = FBUF(I,J) + 
     &           WE  * (FBUF(I+1,J)   - FBUF(I,J)) + 
     &           WW  * (FBUF(I-1,J)   - FBUF(I,J)) + 
     &           WS  * (FBUF(I,J+1)   - FBUF(I,J)) + 
     &           WN  * (FBUF(I,J-1)   - FBUF(I,J)) +
     &           WSE * (FBUF(I+1,J+1) - FBUF(I,J)) + 
     &           WNW * (FBUF(I-1,J-1) - FBUF(I,J)) + 
     &           WSW * (FBUF(I-1,J+1) - FBUF(I,J)) + 
     &           WNE * (FBUF(I+1,J-1) - FBUF(I,J))
         ENDDO
      ENDDO

      RETURN
      END


C ----------------------------DUMMIES ----------------------

      SUBROUTINE DUMMIES(V, NX, NY)        

C     CREATES DUMMY BOUNDARIES BY MIRRORING

      REAL, DIMENSION(NX+2,NY+2) :: V

      DO I=2,NX+1
        V(I,1)    = V(I,2)
        V(I,NY+2) = V(I,NY+1)
      ENDDO

      DO J=1, NY+2
        V(1,J)    = V(2,J)
        V(NX+2,J) = V(NX-1,J)
      ENDDO

      RETURN
      END


C -------------------------- PA_BACKTRANS --------------------------

      SUBROUTINE PA_BACKTRANS( C, S, EIG1, EIG2, A11, A12, A22)

C      C        1. COMP. OF 1. EIGENVECTOR 
C      S        2. COMP. OF 1. EIGENVECTOR 
C      EIG1     1. EIGENVALUE 
C      EIG2     2. EIGENVALUE 
C      A11      COEFF. OF (2*2)-MATRIX, OUTPUT 
C      A12      COEFF. OF (2*2)-MATRIX, OUTPUT 
C      A22      COEFF. OF (2*2)-MATRIX, OUTPUT 

C      PRINCIPAL AXIS BACKTRANSFORMATION OF A SYMMETRIC (2*2)-MATRIX. 
C      A   = U * DIAG(EIG1, EIG2) * U_TRANSPOSE WITH U = (V1 | V2)     
C      V1  = (C, S) IS FIRST EIGENVECTOR

       A11 = C * C * EIG1 + S * S * EIG2
       A22 = EIG1  + EIG2 -  A11           
       A12 = C * S * (EIG1 - EIG2)

       RETURN
       END


C   ---------------------------- DIFF_TENSOR ----------------------
C
C
C   PURPOSE:  CALCULATES DIFFUSION TENSOR OF EED.
C
C   PARAMETERS:
C      V        IN: REGULARIZED IMAGE, UNCHANGED
C      FLAMBDA  CONTRAST PARAMETER
C      NX       IMAGE DIMENSION IN X DIRECTION
C      NY       IMAGE DIMENSION IN Y DIRECTION
C      HX       PIXEL WIDTH     IN X DIRECTION
C      HY       PIXEL WIDTH     IN Y DIRECTION
C      DXX      OUT: DIFFUSION TENSOR ELEMENT
C      DXY      OUT: DIFFUSION TENSOR ELEMENT 
C      DYY      OUT: DIFFUSION TENSOR ELEMENT 
C
C   VARIABLES:
C      DV_DX, DV_DY        DERIVATIVES OF V
C      TWO_HX, TWO_HY      TIME SAVERS
C      C, S                SPECIFY FIRST EIGENVECTOR
C      GRAD                |GRAD(V)|
C      EIG1, EIG2          EIGENVALUES OF DIFFUSION TENSOR


       SUBROUTINE DIFF_TENSOR(V, FLAMBDA, NX, NY,  HX,  HY, 
     &                        DXX, DXY, DYY)

       REAL, DIMENSION(NX+2,NY+2) :: V,DXX,DYY,DXY

       EIG2   = 1.0
       TWO_HX = 2.0 * HX 
       TWO_HY = 2.0 * HY 
       CALL DUMMIES(V, NX, NY)

       DO I=2, NX+1
          DO J=2, NY+1

C            CALCULATE GRAD(V)
             DV_DX = (V(I+1,J) - V(I-1,J)) / TWO_HX
             DV_DY = (V(I,J+1) - V(I,J-1)) / TWO_HY
             GRAD  = SQRT(DV_DX * DV_DX + DV_DY * DV_DY)

C            CALCULATE FIRST EIGENVECTOR OF DIFFUSION TENSOR (NORMALIZED)
C            CALCULATE EIGENVALUES
             IF (GRAD .GT.  0.0) THEN
                C    = DV_DX / GRAD
                S    = DV_DY / GRAD
                EIG1 = 1.0 - EXP(-3.31488 / ((GRAD/FLAMBDA)**4.0)) 
             ELSE
                C    = 1.0
                S    = 0.0
                EIG1 = 1.0
             ENDIF

C            PRINCIPAL AXIS BACKTRANSFORMATION
             CALL PA_BACKTRANS(C, S, EIG1, EIG2, 
     &                         DXX(I,J), DXY(I,J), DYY(I,J))
          ENDDO
      ENDDO

      RETURN
      END

C   -------------------------- GAUSS_CONV --------------------------
C
C   PURPOSE:
C      GAUSSIAN COONVOLUTION. 
C
C   PARAMETERS:
C      SIGMA              STANDARD DEVIATION OF GAUSSIAN 
C      NX                 IMAGE DIMENSION IN X DIRECTION  
C      NY                 IMAGE DIMENSION IN Y DIRECTION  
C      HX                 PIXEL SIZE IN X DIRECTION 
C      HY                 PIXEL SIZE IN Y DIRECTION 
C      PRECISION          CUTOFF AT PRECISION * SIGMA 
C      FBUF               INPUT: ORIGINAL IMAGE   OUTPUT: SMOOTHED 
C                         0=DIRICHLET, 1=REFLECING, 2=PERIODIC 
C  VARIABLES
C      SUM                FOR SUMMING UP
C      CONV               CONVOLUTION VECTOR 
C      HELP               ROW OR COLUMN WITH DUMMY BOUNDARIES
C
C  BOUNDARY TYPE SET TO: 1 
C

      SUBROUTINE GAUSS_CONV (SIGMA,NX,NY,HX,HY,PRECISION,FBUF)

      INTEGER  ::                        P
      REAL, ALLOCATABLE, DIMENSION(:) :: CONV, HELP
      REAL, DIMENSION(NX+2,NY+2) ::      FBUF
      
C     DIFFUSION IN X DIRECTION ------------

C     CALCULATE LENGTH OF CONVOLUTION VECTOR
      LENGTH = (PRECISION * SIGMA / HX) + 1
      IF (LENGTH .GT. NX) THEN
         CALL ERRT(101,'GAUSS_CONV: SIGMA TOO LARGE',NE) 
         RETURN
      ENDIF

C     ALLOCATE STORAGE FOR CONVOLUTION VECTOR & ROW
      ALLOCATE(CONV(LENGTH+1), HELP(NX+LENGTH+LENGTH), STAT=IRTFLG)
      IF (IRTFLG .NE. 0)  THEN
         CALL  ERRT(101,'UNABLE TO ALLOCATE CONV & ROW',NE)
         RETURN
      ENDIF

C     CALCULATE ENTRIES OF CONVOLUTION VECTOR
      SUM = 0.0
      DO I=1, LENGTH+1
         CONV(I) = 1 / (SIGMA * SQRT(2.0 * 3.1415926)) *
     &       EXP (- ((I-1) * (I-1) * HX * HX) / (2.0 * SIGMA * SIGMA))
         SUM = SUM + 2.0 * CONV(I)
      ENDDO

C     DIVISION OVER ALL OF CONV VECTOR
      CONV = CONV / SUM

      DO  J=2, NY+1

C        COPY IN ROW VECTOR
         DO I=2,NX+1
            HELP(I+LENGTH-1) = FBUF(I,J)
         ENDDO

C        PERIODIC B.C.
         DO P=0, LENGTH-1
C           LEFT BOUNDARY
            HELP(LENGTH-P)      = HELP(NX+LENGTH-P)
C           RIGHT BOUNDARY
            HELP(NX+LENGTH+1+P) = HELP(LENGTH+P+1)
         ENDDO

C        CONVOLUTION STEP
         DO I=LENGTH+1, NX+LENGTH

C           CALCULATE CONVOLUTION
            SUM = 0.0

            DO P=1, LENGTH+1
               SUM = SUM + CONV(P) * (HELP(I+P) + HELP(I-P))
            ENDDO

C           WRITE BACK
            FBUF(I-LENGTH+1,J) = SUM
         ENDDO
       ENDDO

       IF (ALLOCATED(HELP)) DEALLOCATE(HELP)
       IF (ALLOCATED(CONV)) DEALLOCATE(CONV)


C      --------------------- DIFFUSION IN Y DIRECTION ----------------

C     CALCULATE LENGTH OF CONVOLUTION VECTOR
      LENGTH = (PRECISION * SIGMA / HY) + 1
      IF (LENGTH .GT. NY) THEN
         CALL ERRT(101,'GAUSS_CONV: SIGMA TOO LARGE',NE) 
         RETURN
      ENDIF

C     ALLOCATE STORAGE FOR CONVOLUTION VECTOR & ROW
      ALLOCATE(CONV(LENGTH+1),HELP(NY+LENGTH+LENGTH), STAT=IRTFLG)
      IF (IRTFLG .NE. 0)  THEN
         CALL  ERRT(101,'UNABLE TO ALLOCATE CONV & ROW',NE)
         RETURN
      ENDIF

C     CALCULATE ENTRIES OF CONVOLUTION VECTOR
      SUM = 0.0
      DO I=1, LENGTH+1
         CONV(I) = 1 / (SIGMA * SQRT(2.0 * 3.1415926)) *
     &           EXP (- ((I-1) * (I-1) * HY * HY) / 
     &           (2.0 * SIGMA * SIGMA))
         SUM = SUM + 2.0 * CONV(I)
      ENDDO

      CONV = CONV / SUM

      DO  I=2, NX+1

C        COPY IN COL VECTOR
         DO J=1+1,NY+1
            HELP(J+LENGTH-1) = FBUF(I,J)
         ENDDO

C        ASSIGN BOUNDARY CONDITIONS, PERIODIC B.C.
         DO P=0, LENGTH-1
C           LEFT BOUNDARY
            HELP(LENGTH-P)      = HELP(NY+LENGTH-P)
C           RIGHT BOUNDARY
            HELP(NY+LENGTH+1+P) = HELP(LENGTH+P+1)
         ENDDO

C        CONVOLUTION STEP
         DO J=LENGTH+1, NY+LENGTH

C           CALCULATE CONVOLUTION
            SUM = 0.0

            DO P=1, LENGTH+1
               SUM = SUM + CONV(P) * (HELP(J+P) + HELP(J-P))
            ENDDO

C           WRITE BACK
            FBUF(I,J-LENGTH+1) = SUM
         ENDDO
      ENDDO

      IF (ALLOCATED(HELP)) DEALLOCATE(HELP)
      IF (ALLOCATED(CONV)) DEALLOCATE(CONV)
 
      RETURN
      END

 
