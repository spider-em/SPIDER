
C ++********************************************************************
C                                                                      *
C  ANISOF3   
C                                                                      *
C
C  ANISOF3(LUN1,LUN2,NSAM,NROW,NSLICE,IRTFLG)
C
C  PARAMETERS: LUN1,LUN2   IO UNITS                             (INPUT)
C              NSAM        X DIMENSIONS                         (INPUT)
C              NROW        Y DIMENSIONS                         (INPUT)
C              NSLICE      Z DIMENSIONS                         (INPUT)
C
C  PURPOSE: ALTER CONTRAST IN AN IMAGE OR VOLUME USING ANISOTROPIC
C           DIFFUSION
C                                                                      *
C
C  COPYRIGHT (C) 2001 BY ACHILLES FRANGAKIS.
C
C  PORTED FROM C CODE BY A. LEITH
C
C **********************************************************************

	SUBROUTINE ANISOF3(LUN1,LUN2,NSAM,NROW,NSLICE,ITER,HT,
     &                     FLAMBDA,SIGMA,IRTFLG)

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: BUF

        ALLOCATE(BUF(NSAM+2,NROW+2,NSLICE+2), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'ANISO, OUT...',IER)
           RETURN
        ENDIF

        HX = 1
        HY = 1
        Hz = 1

        DO ISLICE = 1,NSLICE
C          LOAD THIS SLICE INTO CENTRAL PORTION OF BUF 
           J = 2
           DO IREC=(ISLICE-1)*NROW+1, ISLICE*NROW
              CALL REDLIN(LUN1,BUF(2,J,ISLICE+1),NSAM,IREC)
              J = J + 1
           ENDDO
        ENDDO

C       REPORT MINIMUM, MAXIMUM, MEAN, VARIANCE
        CALL ANALYSE3(BUF,NSAM,NROW,NSLICE,FMINT,FMAXT,FAVT,FSIGT)
        K = 0
        WRITE(NOUT,91) K,FMINT,FMAXT,FAVT,FSIGT
91      FORMAT(' ITER: ',I5,'  MIN: ',G13.7,'  MAX: ',G13.7,
     &            '  AV: ',G13.7,'  SIG: ',G13.7)

        DO K=1,ITER

C          NON-LINEAR ANISOTROPIC DIFFUSION
           CALL EED3(HT,NSAM,NROW,NSLICE,HX,HY,HZ,SIGMA,FLAMBDA,BUF)
 
C          REPORT MINIMUM, MAXIMUM, MEAN, VARIANCE
           CALL ANALYSE3(BUF,NSAM,NROW,NSLICE,FMINT,FMAXT,FAVT,FSIGT)
           WRITE(NOUT,91) K,FMINT,FMAXT,FAVT,FSIGT

        ENDDO

        DO ISLICE = 1,NSLICE
C          LOAD THIS SLICE INTO CENTRAL PORTION OF BUF 
           J = 2
           DO IREC=(ISLICE-1)*NROW+1, ISLICE*NROW
              CALL WRTLIN(LUN2,BUF(2,J,ISLICE+1),NSAM,IREC)
              J = J + 1
           ENDDO
        ENDDO
  
        IF (ALLOCATED(BUF)) DEALLOCATE(BUF)
        RETURN
        END

C     -------------------------ANALYSE3 ----------------------

      SUBROUTINE ANALYSE3(BUF,NX,NY,NZ,FMINT,FMAXT,FAVT,FSIGT)

      REAL, DIMENSION(NX+2,NY+2,NZ+2) :: BUF

      DOUBLE PRECISION DAV,DAV2,DTOP,FNALL,DTEMP

      NPIX   = NX * NY * NZ

      DAV    = 0.0
      DAV2   = 0.0

      FMAXT  = BUF(2,2,2)
      FMINT  = BUF(2,2,2)
      
      DO I = 2,NZ+1
         DO J = 2,NY+1
            DO K = 2,NX+1
               B     = BUF(K,J,I)
               FMAXT = MAX(B,FMAXT)
               FMINT = MIN(B,FMINT)
               DAV   = DAV  + DBLE(B)
               DAV2  = DAV2 + DBLE(B) * DBLE(B)
            ENDDO
         ENDDO
      ENDDO

      DTOP  = DAV2 - DAV * DAV / DBLE(NPIX)

      FAVT  = DAV / DBLE(NPIX)
      FSIGT = DSQRT( DTOP / DBLE(NPIX-1.0))

      END
 
 

C ----------------------- EED3 -----------------------------------
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

C     EED3
C
C     Thanks to: Jose-Jesus Fernandez, MRC-LMB, Cambridge for fixing
C                some bugs in the FORTRAN port. al.
C           
C     PURPOSE:
C        COHERENCE-ENHANCING ANISOTROPIC DIFFUSION.
C        EXPLICIT DISCRETIZATION.
C
C     PARAMETERS:
C        HT       TIME STEP SIZE, 0 < HT <= 0.25
C        NX       IMAGE DIMENSION IN X DIRECTION
C        NY       IMAGE DIMENSION IN Y DIRECTION
C        NZ       IMAGE DIMENSION IN Z DIRECTION
C        HX       PIXEL SIZE IN X DIRECTION
C        HY       PIXEL SIZE IN Y DIRECTION
C        HZ       PIXEL SIZE IN Z DIRECTION
C        SIGMA    NOISE SCALE
C        FLAMBDA  LAMDA PARAMETER FOR EED3
C        U        INPUT: ORIGINAL IMAGE,  OUTPUT: SMOOTHED
C
C     VARIABLES:
C        RXX, RXY, RXZ, RYY, RYZ, RZZ                  TIME SAVERS
C        WN, WNE, WE, WSE                              WEIGHTS
C        WS, WSW, WW, WNW                              WEIGHTS
C        WB, WSB, WNB, WEB, WWB                        WEIGHTS
C        WF, WSF, WNF, WEF, WWF                        WEIGHTS
C

       SUBROUTINE EED3(HT,NX,NY,NZ,HX,HY,HZ,SIGMA,FLAMBDA,U)

       REAL, DIMENSION(NX+2,NY+2,NZ+2) :: U

       REAL, ALLOCATABLE,DIMENSION(:,:,:) :: F,GRD,DXX,DXY,DXZ,DYY,
     &                                       DYZ,DZZ

C      AUTOMATIC ARRAYS

       ALLOCATE(F(NX+2,NY+2,NZ+2),  GRD(NX+2,NY+2,NZ+2),
     &          DXX(NX+2,NY+2,NZ+2),DXY(NX+2,NY+2,NZ+2),
     &          DXZ(NX+2,NY+2,NZ+2),DYY(NX+2,NY+2,NZ+2),
     &          DYZ(NX+2,NY+2,NZ+2),DZZ(NX+2,NY+2,NZ+2),
     &            STAT=IRTFLG)
       IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'ANISO, V',NE)
           RETURN
       ENDIF

C     COPY U INTO F
      DO I=1+1,NX+1
         DO J=1+1,NY+1
            DO K=1+1,NZ+1
               F(I,J,K) = U(I,J,K)
            ENDDO
         ENDDO
      ENDDO

C     CALCULATE ENTRIES OF STRUCTURE TENSOR FOR EED3
      CALL STRUCT_TENSOR_EED(F, NX, NY, NZ, HX, HY, HZ, SIGMA,
     &                      DXX, DXY, DXZ, DYY, DYZ, DZZ, GRD)

C     CALCULATE ENTRIES OF DIFFUSION TENSOR
      CALL DIFF_TENSOR3(FLAMBDA, NX, NY, NZ, HX, HY, HZ,
     &                SIGMA, GRD, DXX, DXY, DXZ, DYY, DYZ, DZZ)

C     CALCULATE EXPLICIT NONLINEAR DIFFUSION OF U

      RXX  = HT / (2.0 * HX * HX)
      RYY  = HT / (2.0 * HY * HY)
      RZZ  = HT / (2.0 * HZ * HZ)
      RXY  = HT / (4.0 * HX * HY)
      RXZ  = HT / (4.0 * HX * HZ)
      RYZ  = HT / (4.0 * HY * HZ)

      CALL DUMMIES3(DXX, NX, NY, NZ)
      CALL DUMMIES3(DXY, NX, NY, NZ)
      CALL DUMMIES3(DXZ, NX, NY, NZ)
      CALL DUMMIES3(DYY, NX, NY, NZ)
      CALL DUMMIES3(DYZ, NX, NY, NZ)
      CALL DUMMIES3(DZZ, NX, NY, NZ)

C     COPY U INTO F AND ASSIGN DUMMY BOUNDARIES
      DO I=1+1,NX+1
         DO J=1+1,NY+1
            DO K=1+1,NZ+1
               F(I,J,K) = U(I,J,K)
            ENDDO
         ENDDO
       ENDDO

      CALL DUMMIES3(F, NX, NY, NZ)

C     DIFFUSE
      DO I=1+1,NX+1
         DO J=1+1,NY+1
            DO K=1+1,NZ+1

C              WEIGHTS
               WE  = RXX * (DXX(I+1,J,K) + DXX(I,J,K)) - RXY * 
     &               (ABS(DXY(I+1,J,K)) + 
     &                ABS(DXY(I,J,K))) - RXZ * 
     &               (ABS(DXZ(I+1,J,K)) + 
     &                ABS(DXZ(I,J,K)))

               WW  =   RXX * (DXX(I-1,J,K) + DXX(I,J,K)) - RXY * 
     &               (ABS(DXY(I-1,J,K)) + 
     &                ABS(DXY(I,J,K))) - RXZ * 
     &               (ABS(DXZ(I-1,J,K)) + 
     &                ABS(DXZ(I,J,K)))

               WS  =   RYY * (DYY(I,J+1,K) + DYY(I,J,K)) - RXY * 
     &               (ABS(DXY(I,J+1,K)) + 
     &                ABS(DXY(I,J,K))) - RYZ * 
     &               (ABS(DYZ(I,J+1,K)) + 
     &                ABS(DYZ(I,J,K)))

               WN  =   RYY * (DYY(I,J-1,K) + DYY(I,J,K)) - RXY * 
     &               (ABS(DXY(I,J-1,K)) + 
     &                ABS(DXY(I,J,K))) - RYZ * 
     &               (ABS(DYZ(I,J-1,K)) + 
     &                ABS(DYZ(I,J,K)))

               WB  =   RZZ * (DZZ(I,J,K-1) + DZZ(I,J,K)) - RYZ * 
     &               (ABS(DYZ(I,J,K-1)) + 
     &                ABS(DYZ(I,J,K))) - RXZ * 
     &               (ABS(DXZ(I,J,K-1)) + 
     &               ABS(DXZ(I,J,K)))

               WF  =   RZZ * (DZZ(I,J,K+1) + DZZ(I,J,K)) - RYZ * 
     &               (ABS(DYZ(I,J,K+1)) + 
     &                ABS(DYZ(I,J,K))) - RXZ * 
     &               (ABS(DXZ(I,J,K+1)) + 
     &                ABS(DXZ(I,J,K)))

               WSE =  RXY * (   DXY(I+1,J+1,K) + DXY(I,J,K) + 
     &                ABS(DXY(I+1,J+1,K)) + 
     &                ABS(DXY(I,J,K)))

               WNW =   RXY * (   DXY(I-1,J-1,K) + DXY(I,J,K) + 
     &               ABS(DXY(I-1,J-1,K)) + 
     &               ABS(DXY(I,J,K)))

               WNE =   RXY * ( - DXY(I+1,J-1,K) - DXY(I,J,K) + 
     &               ABS(DXY(I+1,J-1,K)) + 
     &               ABS(DXY(I,J,K)))

               WSW =   RXY * ( - DXY(I-1,J+1,K) - DXY(I,J,K) + 
     &               ABS(DXY(I-1,J+1,K)) + 
     &               ABS(DXY(I,J,K)))

               WSF =   RYZ * (   DYZ(I,J+1,K+1) + DYZ(I,J,K) + 
     &               ABS(DYZ(I,J+1,K+1)) + 
     &               ABS(DYZ(I,J,K)))

               WNF =   RYZ * ( - DYZ(I,J-1,K+1) - DYZ(I,J,K) + 
     &               ABS(DYZ(I,J-1,K+1)) + 
     &               ABS(DYZ(I,J,K)))

               WEF =   RXZ * (   DXZ(I+1,J,K+1) + DXZ(I,J,K) + 
     &               ABS(DXZ(I+1,J,K+1)) + 
     &               ABS(DXZ(I,J,K)))

               WWF =   RXZ * ( - DXZ(I-1,J,K+1) - DXZ(I,J,K) + 
     &               ABS(DXZ(I-1,J,K+1)) + 
     &               ABS(DXZ(I,J,K)))

               WSB =   RYZ * ( - DYZ(I,J+1,K-1) - DYZ(I,J,K) + 
     &               ABS(DYZ(I,J+1,K-1)) + 
     &               ABS(DYZ(I,J,K)))

               WNB =   RYZ * (   DYZ(I,J-1,K-1) + DYZ(I,J,K) + 
     &               ABS(DYZ(I,J-1,K-1)) + 
     &               ABS(DYZ(I,J,K)))

               WEB =   RXZ * ( - DXZ(I+1,J,K-1) - DXZ(I,J,K) + 
     &               ABS(DXZ(I+1,J,K-1)) + 
     &               ABS(DXZ(I,J,K)))

               WWB =   RXZ * (   DXZ(I-1,J,K-1) + DXZ(I,J,K) + 
     &               ABS(DXZ(I-1,J,K-1)) + 
     &               ABS(DXZ(I,J,K)))


C              MODIFY WEIGHTS TO PREVENT FLUX ACROSS BOUNDARIES
               IF (I .EQ. 1+1) THEN
C                 SET WESTERN WEIGHTS ZERO
                  WSW = 0.0
                  WW  = 0.0
                  WNW = 0.0
                  WWB = 0.0
                  WWF = 0.0
               ENDIF

               IF (I .EQ. NX+1) THEN
C                 SET EASTERN WEIGHTS ZERO
                  WNE = 0.0
                  WE  = 0.0
                  WSE = 0.0
                  WEB = 0.0
                  WEF = 0.0
               ENDIF

               IF (J .EQ. 1+1) THEN
C                 SET NORTHERN WEIGHTS ZERO
                  WNW = 0.0
                  WN  = 0.0
                  WNE = 0.0
                  WNB = 0.0
                  WNF = 0.0
               ENDIF

               IF (J .EQ. NY+1) THEN
C                 SET SOUTHERN WEIGHTS ZERO
                  WSE = 0.0
                  WS  = 0.0
                  WSW = 0.0
                  WSB = 0.0
                  WSF = 0.0
               ENDIF

               IF (K .EQ. NZ+1) THEN
C                 SET FORWARD WEIGHTS ZERO
                  WEF = 0.0
                  WF  = 0.0
                  WWF = 0.0
                  WNF = 0.0
                  WSF = 0.0
               ENDIF

               IF (K .EQ. 1+1) THEN
C                 SET BACKWARD WEIGHTS ZERO
                  WSB = 0.0
                  WB  = 0.0
                  WNB = 0.0
                  WEB = 0.0
                  WWB = 0.0
               ENDIF

C              EVOLUTION
               U(I,J,K) = F(I,J,K)
     &            + WE  * (F(I+1,J,K)   - F(I,J,K))
     &            + WW  * (F(I-1,J,K)   - F(I,J,K))
     &            + WS  * (F(I,J+1,K)   - F(I,J,K))
     &            + WN  * (F(I,J-1,K)   - F(I,J,K))
     &            + WB  * (F(I,J,K-1)   - F(I,J,K))
     &            + WF  * (F(I,J,K+1)   - F(I,J,K))
     &            + WSE * (F(I+1,J+1,K) - F(I,J,K))
     &            + WNW * (F(I-1,J-1,K) - F(I,J,K))
     &            + WSW * (F(I-1,J+1,K) - F(I,J,K))
     &            + WNE * (F(I+1,J-1,K) - F(I,J,K))
     &            + WNB * (F(I,J-1,K-1) - F(I,J,K))
     &            + WNF * (F(I,J-1,K+1) - F(I,J,K))
     &            + WEB * (F(I+1,J,K-1) - F(I,J,K))
     &            + WEF * (F(I+1,J,K+1) - F(I,J,K))
     &            + WWB * (F(I-1,J,K-1) - F(I,J,K))

     &            + WWF * (F(I-1,J,K+1) - F(I,J,K))
     &            + WSB * (F(I,J+1,K-1) - F(I,J,K))
     &            + WSF * (F(I,J+1,K+1) - F(I,J,K))
            ENDDO
         ENDDO
      ENDDO


C     DEALLOCATE STORAGE
      DEALLOCATE(F,GRD,DXX,DXY,DXZ,DYY,DYZ,DZZ)

      RETURN
      END



C ----------------------------DUMMIES3 ----------------------
C
C     Thanks to: Jose-Jesus Fernandez, MRC-LMB, Cambridge for fixing
C                some bugs in the FORTRAN port. al.
C           

      SUBROUTINE DUMMIES3(V, NX, NY, NZ)

C     CREATES DUMMY BOUNDARIES BY PERIODIC CONTINUATION

      REAL, DIMENSION(NX+2,NY+2,NZ+2) :: V

      DO I=2,NX+1
         DO J=2,NY+1
            V(I,J,1)    = V(I,J,NZ+1)
            V(I,J,NZ+2) = V(I,J,2)
        ENDDO
      ENDDO

      DO J=2,NY+1
         DO K=2,NZ+1
           V(1,J,K)    = V(NX+1,J,K)
           V(NX+2,J,K) = V(2,J,K)
         ENDDO
      ENDDO

      DO K=2,NZ+1
         DO I=2,NX+1
           V(I,1,K)    = V(I,NY+1,K)
           V(I,NY+2,K) = V(I,2,K)
         ENDDO
      ENDDO

      RETURN
      END

C -------------------------- PA_BACKTRANS3 --------------------------

      SUBROUTINE PA_BACKTRANS3(EV11,EV12,EV13,EV21,EV22,EV23,EV31,
     &                        EV32, EV33, EIG1, EIG2, EIG3,
     &                        A11,A12,A13,A22,A23,A33) 

C     EV11   1. COMP. OF 1. EIGENVECTOR 
C     EV12   2. COMP. OF 1. EIGENVECTOR 
C     EV13   3. COMP. OF 1. EIGENVECTOR
C     EV21   1. COMP. OF 2. EIGENVECTOR 
C     EV22   2. COMP. OF 2. EIGENVECTOR 
C     EV23   3. COMP. OF 2. EIGENVECTOR
C     EV31   1. COMP. OF 3. EIGENVECTOR 
C     EV32   2. COMP. OF 3. EIGENVECTOR 
C     EV33   3. COMP. OF 3. EIGENVECTOR
C     EIG1   1. EIGENVALUE 
C     EIG2   2. EIGENVALUE 
C     EIG3   3. EIGENVALUE 
C     A11   COEFF. OF (3*3)-MATRIX, OUTPUT 
C     A12   COEFF. OF (3*3)-MATRIX, OUTPUT 
C     A13   COEFF. OF (3*3)-MATRIX, OUTPUT 
C     A22   COEFF. OF (3*3)-MATRIX, OUTPUT 
C     A23   COEFF. OF (3*3)-MATRIX, OUTPUT 
C     A33   COEFF. OF (3*3)-MATRIX, OUTPUT 


C      PRINCIPAL AXIS BACKTRANSFORMATION OF A SYMMETRIC (3*3)-MATRIX. 
C      A   = U * DIAG(EIG1, EIG2, EIG3) * U_TRANSPOSE WITH 
C      U = (V1 | V2 | V3)     
C      V1  = (EV11, 3V12, EV13) IS FIRST EIGENVECTOR

       A11 = EIG1*EV11*EV11 + EIG2*EV21*EV21 + EIG3*EV31*EV31
       A12 = EIG1*EV11*EV12 + EIG2*EV21*EV22 + EIG3*EV31*EV32
       A13 = EIG1*EV11*EV13 + EIG2*EV21*EV23 + EIG3*EV31*EV33
       A22 = EIG1*EV12*EV12 + EIG2*EV22*EV22 + EIG3*EV32*EV32
       A23 = EIG1*EV12*EV13 + EIG2*EV22*EV23 + EIG3*EV32*EV33
       A33 = EIG1*EV13*EV13 + EIG2*EV23*EV23 + EIG3*EV33*EV33

       RETURN
       END
C   ---------------------------- DIFF_TENSOR3 ----------------------
C
C
C   PURPOSE:  CALCULATES DIFFUSION TENSOR OF EED.
C
C   PARAMETERS:
C     FLAMBDA  EDGE ENHANCING DIFFUSION COEFFICIENT
C     NX       IMAGE DIMENSION IN X DIRECTION
C     NY       IMAGE DIMENSION IN Y DIRECTION
C     NZ       IMAGE DIMENSION IN Z DIRECTION
C     HX       PIXELSIZE IN X DIRECTION
C     HY       PIXELSIZE IN Y DIRECTION
C     HZ       PIXELSIZE IN Z DIRECTION
C     SIGMA    NOISE SCALE
C     GRD      GRADIENT OF THE IMAGE
C     DXX      IN: STRUCTURE TENSOR EL., OUT: DIFF. TENSOR EL.
C     DXY      IN: STRUCTURE TENSOR EL., OUT: DIFF. TENSOR EL. 
C     DXZ      IN: STRUCTURE TENSOR EL., OUT: DIFF. TENSOR EL. 
C     DYY      IN: STRUCTURE TENSOR EL., OUT: DIFF. TENSOR EL.
C     DYZ      IN: STRUCTURE TENSOR EL., OUT: DIFF. TENSOR EL. 
C     DZZ      IN: STRUCTURE TENSOR EL., OUT: DIFF. TENSOR EL. 
C
C   VARIABLES:
C      DV_DX, DV_DY        DERIVATIVES OF V
C      TWO_HX, TWO_HY      TIME SAVERS
C      C, S                SPECIFY FIRST EIGENVECTOR
C      GRAD                |GRAD(V)|
C      EIG1, EIG2          EIGENVALUES OF DIFFUSION TENSOR
C      EV11, EV12, EV13    SPECIFY FIRST EIGENVECTOR
C      EV21, EV22, EV23    SPECIFY SECOND EIGENVECTOR
C      EV31, EV32, EV33    SPECIFY THIRD EIGENVECTOR
C      EMU1, EMU2, EMU3    EIGENVALUES OF STRUCTURE TENSOR
C      EIG1, EIG2, EIG3    EIGENVALUES OF DIFFUSION TENSOR
C      V                   IMAGE
C      AA                  REAL TENSOR MATRIX
C      D                   MATRIX WITH UNORDERED EIGENVALUES
C      V_NORM              MATRIX WITH UNORDERED EIGENVECTORS
C      FNROT               NUMBER OF ITERATIONS THAT WERE REQUIRED


       SUBROUTINE DIFF_TENSOR3(FLAMBDA,NX,NY,NZ,HX,HY,HZ,SIGMA,GRD,
     &                        DXX,DXY,DXZ,DYY,DYZ,DZZ)

       REAL, DIMENSION(NX+2,NY+2,NZ+2) :: DXX,DXY,DXZ,DYY,DYZ,DZZ
       REAL, DIMENSION(NX+2,NY+2,NZ+2) :: GRD

C      ALLOCATABLE
       REAL, ALLOCATABLE,DIMENSION(:,:,:) :: V

C      AUTOMATIC ARRAYS
       REAL, DIMENSION(3,3) :: AA,V_NORM
       REAL, DIMENSION(3)   :: D

       ALLOCATE(V(NX+2,NY+2,NZ+2), STAT=IRTFLG)
       IF (IRTFLG .NE. 0) THEN 
          CALL ERRT(46,'ANISO, V',NE)
          RETURN
       ENDIF

       DO I=2, NX+1
          DO J=2, NY+1
             DO K=2, NZ+1
                CALL READ_STRUCT_TENS(DXX(I,J,K), DXY(I,J,K), 
     &                                DXZ(I,J,K), DYY(I,J,K), 
     &                                DYZ(I,J,K), DZZ(I,J,K), AA)

                CALL JACOBI(AA, 3, D, V_NORM, FNROT)
                CALL EIGSRT(D, V_NORM, 3)
                CALL READ_EI(D, V_NORM, 
     &                       EV11, EV12, EV13, EV21, EV22, EV23, 
     &                       EV31, EV32, EV33, EMU1, EMU2, EMU3)
         
                IF (GRD(I,J,K) .GT. 0.0) THEN
                   EIG1 = 1.0-EXP(-3.31488/(GRD(I,J,K)/FLAMBDA**16.0)) 
                   EIG2 = 1.0-EXP(-3.31488/(GRD(I,J,K)/FLAMBDA**16.0))
                ELSE
                   EIG1 = 1
                   EIG2 = 1
                ENDIF

                EIG3 = 1
                CALL PA_BACKTRANS3(EV11,EV12,EV13,EV21,EV22,EV23,
     &                             EV31,EV32,EV33,EIG1,EIG2,EIG3, 
     &                             DXX(I,J,K), DXY(I,J,K), 
     &                             DXZ(I,J,K), DYY(I,J,K), 
     &                             DYZ(I,J,K), DZZ(I,J,K)) 

             ENDDO
          ENDDO
      ENDDO

      IF (ALLOCATED(V)) DEALLOCATE(V)

      RETURN
      END

C ----------------------- EIGSRT -------------------------------

      SUBROUTINE EIGSRT(D,V,N)

      REAL, DIMENSION(3,3) :: V 
      REAL, DIMENSION(3)   :: D

      DO I=1,N
         P = D(I)
         K = I
         DO J=I+1,N
            IF (D(J) .GE. P) THEN
               P = D(J)
               K = J
            ENDIF
         ENDDO

	 IF (K .NE. I) THEN
            D(K) = D(I)
            D(I) = P
            DO J=1,N
               P      = V(J,I)
               V(J,I) = V(J,K)
               V(J,K) = P
            ENDDO
         ENDIF
      ENDDO
      END
C ----------------------- STRUCT_TENSOR_EED -------------------------------

C   PURPOSE:  DEFINES STRUCTURAL TENSOR 
C
C   PARAMETERS:
C     V        IMAGE
C     NX       IMAGE DIMENSION IN X DIRECTION
C     NY       IMAGE DIMENSION IN Y DIRECTION
C     NZ       IMAGE DIMENSION IN Z DIRECTION
C     HX       PIXELSIZE IN X DIRECTION
C     HY       PIXELSIZE IN Y DIRECTION
C     HZ       PIXELSIZE IN Z DIRECTION
C     SIGMA    NOISE SCALE
C     DXX      IN: STRUCTURE TENSOR EL., OUT: DIFF. TENSOR EL.
C     DXY      IN: STRUCTURE TENSOR EL., OUT: DIFF. TENSOR EL. 
C     DXZ      IN: STRUCTURE TENSOR EL., OUT: DIFF. TENSOR EL. 
C     DYY      IN: STRUCTURE TENSOR EL., OUT: DIFF. TENSOR EL.
C     DYZ      IN: STRUCTURE TENSOR EL., OUT: DIFF. TENSOR EL. 
C     DZZ      IN: STRUCTURE TENSOR EL., OUT: DIFF. TENSOR EL. 
C     GRD      GRADIENT OF THE IMAGE

C      BUILDING TENSOR PRODUCT

      SUBROUTINE STRUCT_TENSOR_EED(V,NX,NY,NZ,HX,HY,HZ,SIGMA,
     &                        DXX,DXY,DXZ,DYY,DYZ,DZZ,GRD)

      REAL, DIMENSION(NX+2,NY+2,NZ+2) :: DXX,DXY,DXZ,DYY,DYZ,DZZ
      REAL, DIMENSION(NX+2,NY+2,NZ+2) :: GRD,V

      TWO_HX = 2.0 * HX  
      TWO_HY = 2.0 * HY    
      TWO_HZ = 2.0 * HZ    
      CALL DUMMIES3(V, NX, NY, NZ)

      DO I=2,NX+1
        DO J=2,NY+1
          DO K=2,NZ+1

             DV_DX      = (V(I+1,J,K) - V(I-1,J,K)) / TWO_HX
             DV_DY      = (V(I,J+1,K) - V(I,J-1,K)) / TWO_HY
             DV_DZ      = (V(I,J,K+1) - V(I,J,K-1)) / TWO_HZ
             DXX(I,J,K) = DV_DX * DV_DX
             DXY(I,J,K) = DV_DX * DV_DY 
             DXZ(I,J,K) = DV_DZ * DV_DX
             DYY(I,J,K) = DV_DY * DV_DY 
             DYZ(I,J,K) = DV_DZ * DV_DY
             DZZ(I,J,K) = DV_DZ * DV_DZ
             GRD(I,J,K) = SQRT(DV_DX*DV_DX + DV_DY*DV_DY + DV_DZ*DV_DZ)
           ENDDO
        ENDDO
      ENDDO

C     SMOOTH THE GRADIENT
      IF (SIGMA .GT.  0.0) THEN 
         CALL GAUSS_CONV3(SIGMA, NX, NY, NZ, HX, HY, HZ, 2, DXX)
         CALL GAUSS_CONV3(SIGMA, NX, NY, NZ, HX, HY, HZ, 2, DXY)
         CALL GAUSS_CONV3(SIGMA, NX, NY, NZ, HX, HY, HZ, 2, DXZ)
         CALL GAUSS_CONV3(SIGMA, NX, NY, NZ, HX, HY, HZ, 2, DYY)
         CALL GAUSS_CONV3(SIGMA, NX, NY, NZ, HX, HY, HZ, 2, DYZ)
         CALL GAUSS_CONV3(SIGMA, NX, NY, NZ, HX, HY, HZ, 2, DZZ)
      ENDIF

      RETURN
      END
C ---------------------------- READ_STRUCT_TENS -----------------------

C     CALCULATING EIGENVECTORS AND EIGENVALUES 
C     DEFINE THE STRUCT_TENS MATRIX

      SUBROUTINE READ_STRUCT_TENS(V_DXX,V_DXY,V_DXZ,V_DYY,V_DYZ,
     &                            V_DZZ,A) 

      REAL, DIMENSION(3,3) :: A

      A(1,1) = V_DXX
      A(1,2) = V_DXY
      A(1,3) = V_DXZ
      A(2,1) = V_DXY
      A(2,2) = V_DYY
      A(2,3) = V_DYZ
      A(3,1) = V_DXZ
      A(3,2) = V_DYZ
      A(3,3) = V_DZZ

      RETURN
      END

C ----------------------- ROTATE -----------------------------------

      SUBROUTINE ROTATE(A,I,J,K,L,TAU,S)
 
      REAL, DIMENSION(3,3) :: A

      G      = A(I,J)
      H      = A(K,L)
      A(I,J) = G - S * (H+G*TAU)
      A(K,L) = H + S * (G-H*TAU)

      END

C ----------------------- JACOBI -----------------------------------
 
      SUBROUTINE JACOBI(A, N, D, V, FNROT)

      REAL, DIMENSION(3,3) :: A,V
      REAL, DIMENSION(3)   :: D

C     AUTOMATIC ARRAYS
      REAL, DIMENSION(3)   :: B,Z

C     INITIALIZE WHOLE V,Z ARRAY
      V = 0.0
      Z = 0.0

      DO IP=1,N
	 V(IP,IP) = 1.0
      ENDDO

      DO IP=1,N
         B(IP) = A(IP,IP)
         D(IP) = A(IP,IP)
      ENDDO

      FNROT = 0
      DO I=1,50

         SM = 0.0
         DO IP = 1,N-1
            DO IQ = IP+1,N
               SM = SM + ABS(A(IP,IQ))
            ENDDO
         ENDDO

C        END ITERATIONS IF SM IS ZERO
         IF (SM .EQ. 0.0) RETURN

         IF (I .LT. 4) THEN
            TRESH = 0.2*SM / (N*N)
         ELSE
            TRESH = 0.0
         ENDIF

         DO IP=1,N-1
            DO IQ=IP+1,N
               G = 100.0 * ABS(A(IP,IQ))
               IF (I .GT. 4 .AND. 
     &            (ABS(D(IP))+G) .EQ. (ABS(D(IP))) .AND.
     &            (ABS(D(IQ))+G) .EQ. (ABS(D(IQ)))) THEN
                  A(IP,IQ) = 0.0
                  
               ELSEIF (ABS(A(IP,IQ)) .GT.  TRESH) THEN
                  H = D(IQ) - D(IP)
                  IF ((ABS(H)+G) .EQ. ABS(H)) THEN
                     T = (A(IP,IQ))/H
                  ELSE 
                     THETA = 0.5 * H / (A(IP,IQ))
                     T     = 1.0/(ABS(THETA) + SQRT(1.0+THETA*THETA))
                     IF (THETA .LT.  0.0) T = -T
                  ENDIF
 
                  C        = 1.0 / SQRT(1+T*T)
                  S        = T * C
                  TAU      = S / (1.0+C)
                  H        = T * A(IP,IQ)
                  Z(IP)    = Z(IP) - H
                  Z(IQ)    = Z(IQ) + H
                  D(IP)    = D(IP) - H
                  D(IQ)    = D(IQ) + H
                  A(IP,IQ) = 0.0

                  DO J=1,IP-1
                     CALL ROTATE(A,J,IP,J,IQ,TAU,S)
                  ENDDO
                  DO J=IP+1,IQ-1
                     CALL ROTATE(A,IP,J,J,IQ,TAU,S)
                  ENDDO
                  DO J=IQ+1,N
                     CALL ROTATE(A,IP,J,IQ,J,TAU,S)
                  ENDDO
                  DO J=1,N
                     CALL ROTATE(V,J,IP,J,IQ,TAU,S)
                  ENDDO
                  FNROT = FNROT + 1
               ENDIF
            ENDDO
         ENDDO

         DO IP=1,N
            B(IP) = B(IP) + Z(IP)
            D(IP) = B(IP)
            Z(IP) = 0.0
         ENDDO
      ENDDO

      CALL ERRT(101,"TOO MANY ITERATIONS IN  JACOBI",NE)

      RETURN
      END

C ---------------- READ_EI ------------------------------------------
       SUBROUTINE READ_EI(D,E_VEC,EV11,EV12,EV13,EV21,EV22,EV23,
     &                           EV31,EV32,EV33,EIG1,EIG2,EIG3)

C      PURPOSE:  READ THE EIGENVALUE AND EIGENVECTORS
C
C      PARAMETERS:
C         D,      INPUT, EIGENVALUES IN DESCENDING ORDER
C         E_VEC   INPUT, CORRESPONDING EIGENVECTORS
C         EV11    1. COMP. OF 1. EIGENVECTOR 
C         EV12    2. COMP. OF 1. EIGENVECTOR 
C         EV13    3. COMP. OF 1. EIGENVECTOR
C         EV21    1. COMP. OF 2. EIGENVECTOR 
C         EV22    2. COMP. OF 2. EIGENVECTOR 
C         EV23    3. COMP. OF 2. EIGENVECTOR
C         EV31    1. COMP. OF 3. EIGENVECTOR 
C         EV32    2. COMP. OF 3. EIGENVECTOR 
C         EV33    3. COMP. OF 3. EIGENVECTOR
C         EIG1    1. EIGENVALUE 
C         EIG2    2. EIGENVALUE 
C         EIG3    3. EIGENVALUE 

      REAL, DIMENSION(3) :: D
      REAL, DIMENSION(3,3) :: E_VEC

      EV11 = E_VEC(1,1)
      EV12 = E_VEC(2,1)
      EV13 = E_VEC(3,1)
      EV21 = E_VEC(1,2)
      EV22 = E_VEC(2,2)
      EV23 = E_VEC(3,2)
      EV31 = E_VEC(1,3)
      EV32 = E_VEC(2,3)
      EV33 = E_VEC(3,3)

      EIG1 = D(1)
      EIG2 = D(2)
      EIG3 = D(3)

      RETURN
      END

C     ------------------------ GAUSS_CONV3 ---------------------------
C
C     Thanks to: Jose-Jesus Fernandez, MRC-LMB, Cambridge for fixing
C                some bugs in the FORTRAN port. al.
C           
      SUBROUTINE GAUSS_CONV3(SIGMA, NX, NY, NZ, HX, HY, HZ,
     &           PRECISION, FBUF)

C     SIGMA       STANDARD DEVIATION OF GAUSSIAN
C     NX          IMAGE DIMENSION IN X DIRECTION
C     NY          IMAGE DIMENSION IN Y DIRECTION
C     NZ          IMAGE DIMENSION IN Z DIRECTION
C     HX          PIXEL SIZE IN X DIRECTION
C     HY          PIXEL SIZE IN Y DIRECTION
C     HZ          PIXEL SIZE IN Z DIRECTION
C     PRECISION   CUTOFF AT PRECISION * SIGMA
C     F           INPUT: ORIGINAL IMAGE   OUTPUT: SMOOTHED


C     PURPOSE: GAUSSIAN COONVOLUTION.


C     VARIABLES:
C         LENGTH     CONVOLUTION VECTOR: 0..LENGTH
C         SUM        FOR SUMMING UP
C         CONV       CONVOLUTION VECTOR
C         HELP       ROW OR COLUMN WITH DUMMY BOUNDARIES


      INTEGER  ::                        P,PRECISION
      REAL, ALLOCATABLE, DIMENSION(:) :: HELP
      REAL, DIMENSION(NX+2,NY+2,NZ+2) :: FBUF

C     AUTOMATIC ARRAYS
      REAL, DIMENSION(PRECISION + 1) :: CONV

C     -------------------- DIFFUSION IN X DIRECTION --------------
C234567
C     CALCULATE LENGTH OF CONVOLUTION VECTOR
      LENGTH = PRECISION + 1

      IF (LENGTH > NX .OR. LENGTH > NY .OR. LENGTH > NZ) THEN
         CALL ERRT(101,'GAUSS_CONV: SIGMA TOO LARGE',NE)
         RETURN
      ENDIF

C     CALCULATE ENTRIES OF CONVOLUTION VECTOR
      SUM = 0.0
      DO I=1, LENGTH
         CONV(I) = 1 / (SIGMA * SQRT(2.0 * 3.1415927)) *
     &             EXP (- ((I-1) * (I-1) * HX * HX) /
     &             (2.0 * SIGMA * SIGMA))
         SUM = SUM + 2.0 * CONV(I)
      ENDDO

C     DIVISION OVER ALL OF CONV VECTOR
      CONV = CONV / SUM

C     ALLOCATE STORAGE FOR ROW
      ALLOCATE(HELP(NX+LENGTH+LENGTH), STAT=IRTFLG)
      IF (IRTFLG .NE. 0)  THEN
         CALL  ERRT(101,'UNABLE TO ALLOCATE HELP',NE)
         RETURN
      ENDIF

      DO J=1+1, NY+1
         DO K=1+1, NZ+1

C           COPY IN ROW VECTOR
            DO I=1+1, NX+1
               HELP(I+LENGTH-1) = FBUF(I,J,K)
            ENDDO

C           ASSIGN BOUNDARY CONDITIONS
            DO P=1, LENGTH
               HELP(LENGTH+1-P)  = HELP(NX+LENGTH+1-P)
               HELP(NX+LENGTH+P) = HELP(LENGTH+P)
            ENDDO

C           CONVOLUTION STEP
            DO I=LENGTH+1, NX+LENGTH

C                CALCULATE CONVOLUTION
                 SUM = CONV(1) * HELP(I)

                 DO P=1, LENGTH-1
                    SUM = SUM + CONV(1+P) * (HELP(I+P) + HELP(I-P))
                 ENDDO

C                WRITE BACK
                 FBUF(I-LENGTH+1,J,K) = SUM
            ENDDO
          ENDDO
       ENDDO

       IF (ALLOCATED(HELP)) DEALLOCATE(HELP)

C     ------------------------ DIFFUSION IN Y DIRECTION ------------

C     CALCULATE ENTRIES OF CONVOLUTION VECTOR
      SUM = 0.0
      DO J=1, LENGTH
         CONV(J) = 1 / (SIGMA * SQRT(2.0 * 3.1415927)) *
     &             EXP (- ((J-1) * (J-1) * HY * HY) /
     &             (2.0 * SIGMA * SIGMA))
         SUM = SUM + 2.0 * CONV(J)
      ENDDO

C     DIVISION OVER ALL OF CONV VECTOR
      CONV = CONV / SUM

C     ALLOCATE STORAGE FOR ROW
      ALLOCATE(HELP(NY+LENGTH+LENGTH), STAT=IRTFLG)
      IF (IRTFLG .NE. 0)  THEN
         CALL  ERRT(101,'UNABLE TO ALLOCATE HELP',NE)
         RETURN
      ENDIF

      DO I=1+1, NX+1
         DO K=1+1, NZ+1

C           COPY IN COLUMN VECTOR
            DO J=1+1, NY+1
               HELP(J+LENGTH-1) = FBUF(I,J,K)
            ENDDO

C           ASSIGN BOUNDARY CONDITIONS
            DO P=1, LENGTH
               HELP(LENGTH+1-P)  = HELP(NY+LENGTH+1-P)
               HELP(NY+LENGTH+P) = HELP(LENGTH+P)
            ENDDO

C           CONVOLUTION STEP
            DO J=LENGTH+1, NY+LENGTH

C              CALCULATE CONVOLUTION
               SUM = CONV(1) * HELP(J)

               DO P=1, LENGTH-1
                  SUM = SUM + CONV(1+P) * (HELP(J+P) + HELP(J-P))
               ENDDO

C              WRITE BACK
               FBUF(I,J-LENGTH+1,K) = SUM
            ENDDO
         ENDDO
      ENDDO

      IF (ALLOCATED(HELP)) DEALLOCATE(HELP)


C     -------------- DIFFUSION IN Z DIRECTION ----------------


C     CALCULATE ENTRIES OF CONVOLUTION VECTOR
      SUM = 0.0
      DO K=1, LENGTH
         CONV(K) = 1 / (SIGMA * SQRT(2.0 * 3.1415927)) *
     &             EXP (- ((K-1) * (K-1) * HZ * HZ) /
     &             (2.0 * SIGMA * SIGMA))
         SUM = SUM + 2.0 * CONV(K)
      ENDDO

C     DIVISION OVER ALL OF CONV VECTOR
      CONV = CONV / SUM

C     ALLOCATE STORAGE FOR ROW
      ALLOCATE(HELP(NZ+LENGTH+LENGTH), STAT=IRTFLG)
      IF (IRTFLG .NE. 0)  THEN
         CALL  ERRT(101,'UNABLE TO ALLOCATE HELP',NE)
         RETURN
      ENDIF

      DO I=1+1, NX+1
         DO J=1+1, NY+1

C           COPY IN COLUMN VECTOR
            DO K=1+1, NZ+1
               HELP(K+LENGTH-1) = FBUF(I,J,K)
            ENDDO

C           ASSIGN BOUNDARY CONDITIONS
            DO P=1, LENGTH
               HELP(LENGTH+1-P)  = HELP(NZ+LENGTH+1-P)
               HELP(NZ+LENGTH+P) = HELP(LENGTH+P)
            ENDDO

C           CONVOLUTION STEP
            DO K=LENGTH+1, NZ+LENGTH

C              CALCULATE CONVOLUTION
               SUM = CONV(1) * HELP(K)

               DO P=1, LENGTH-1
                  SUM = SUM + CONV(1+P) * (HELP(K+P) + HELP(K-P))
               ENDDO

C              WRITE BACK
               FBUF(I,J,K-LENGTH+1) = SUM
            ENDDO
         ENDDO
      ENDDO

      IF (ALLOCATED(HELP)) DEALLOCATE(HELP)
      RETURN
      END





