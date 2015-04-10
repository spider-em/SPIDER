C++*********************************************************************
C                                                                      *
C ALRQ_MS_FBS  FBS USED                          JUL 2011 ARDEAN LEITH *
C              FBS2                              DEC 2011 ARDEAN LEITH *
C              FBS2 NXP                          DEC 2011 ARDEAN LEITH *
C              OUT OF BOUNDS TRAP ENLARGED       FEB 2013 ARDEAN LEITH
C                                                                      *
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
C                                                                      *
C ALRQ_MS_FBS(XIM, NX,NY, CNS2,CNR2, NUMR,CIRC,LCIRC,                  *
C             NRING,MODE,NEWFFT,USE_OMP_PARALLEL, AVO,VRIN             *
C                                                                      *
C PURPOSE: INTERPOLATES: XIM IMAGE INTO POLAR COORDINATES AND STACKS   *
C          THE RINGS INTO: CIRC OUTPUT ARRAY.                          * 
C          FFTW3 WILL BE USED ON RINGS.                                *
C                                                                      *
C PARAMETERS:                                                          *
C             XIM              IMAGE ARRAY                     (INPUT) *
C             NX,NY            IMAGE DIMENSIONS                (INPUT) *
C             CNS2,CNR2        PREDEFINED CENTERS              (INPUT) *
C             NUMR             RING CONTROL ARRAY              (INPUT) *
C             CIRC             POLAR ARRAY                     (OUTPUT *)
C             LCIRC            LENGTH OF CIRC ARRAY            (INPUT) *
C             NRING            NUMBER OF RINGS                 (INPUT) *
C             MODE             HALF OR FULL CIRCLE             (INPUT) *
C             NEWFFT                                           (INPUT) *
C             USE_OMP_PARALLEL                                 (INPUT) * 
C             AVO                                              (INPUT) *
C             VRIN                                             (INPUT) *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE ALRQ_MS_FBS(XIM, NX,NY, CNS2,CNR2, 
     &                         NUMR,CIRC,LCIRC,
     &                         NRING,MODE,NEWFFT,USE_OMP_PARALLEL,
     &                         AVO,VRIN)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'

        REAL, INTENT(IN)             :: XIM(NX,NY)
        INTEGER                      :: NX,NY
        REAL                         :: CNS2,CNR2
        REAL, INTENT(INOUT)          :: CIRC(LCIRC)
        INTEGER                      :: LCIRC,NRING
        INTEGER, INTENT(IN)          :: NUMR(3,NRING)
        CHARACTER(LEN=1)             :: MODE
        LOGICAL, INTENT(IN)          :: NEWFFT,USE_OMP_PARALLEL
        DOUBLE PRECISION, INTENT(IN) :: AVO,VRIN

        REAL, ALLOCATABLE            :: F0(:,:),X1(:,:),Y1(:,:),XY2(:,:)

        DOUBLE PRECISION             :: PI,DFI
        REAL                         :: CIRCT

        INTEGER                      :: NXLD,IRTFLG,MWANT,IT,INR,IGO
        INTEGER                      :: NVAL,JT,NE
        INTEGER                      :: LT,LTIGO,LTLTIGO,LTLTLTIGO,NSIM
        REAL                         :: X,Y,YQ,XOLD,YOLD,FI,XT,YT

        REAL                         :: fbs2

C       CNS2 AND CNR2 ARE PREDEFINED CENTERS

        NXLD   = NX + 2 - MOD(NX,2)

        ALLOCATE (F0 (NXLD, NY),
     &            X1 (NXLD, NY),
     &            Y1 (NXLD, NY),
     &            XY2(NXLD, NY),
     &            STAT=IRTFLG)
        IF (IRTFLG /= 0) THEN 
            MWANT = 4*NXLD*NY 
            CALL ERRT(46,'ALRQ_MS_FBS; F0...',MWANT)
            GOTO 9999
        ENDIF  

C       PAD XIM INTO F0 
        F0(1:NX,1:NY) = XIM(1:NX,1:NY)

C       CALCULATE F0 DERIVATIVES: X1,Y1,XY2
        CALL FBS2_PREP(F0, X1,Y1,XY2, NXLD,NX,NY, IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       CNS2 AND CNR2 ARE PREDEFINED CENTERS

C       USING NORMALIZATION JUST FROM RINGS GIVES DIFFERENT ANSWERS
C       FROM NORMASC, PROBABLY OK BUT I WILL NOT USE IT.

        PI     = 2 * DATAN(1.0D0)  ! NOT REALLY PI!!

C       FILL ALL THE RINGS

        IF (USE_OMP_PARALLEL) THEN  ! --------------------------------

c$omp   parallel do private(it,inr,yq,igo,nval,lt,
c$omp&  ltigo,ltltigo,ltltltigo,nsim,dfi,xold,yold,
c$omp&  circt,jt,fi,x,y)
        DO  IT=1,NRING 

           INR  = NUMR(1,IT)        ! RADIUS OF THE CURRENT RING
           YQ   = INR               ! FLOATING POINT RADIUS
           IGO  = NUMR(2,IT)        ! STARTING LOCATION FOR RING 
           NVAL = NUMR(3,IT)        ! LENGTH OF THIS RING

           IF (NEWFFT) THEN
C             THE ACTUAL, POWER-OF-TWO LENGTH IS NUMR(3,I)-2, ADDITIONAL
C             TWO LOCATIONS ARE ONLY FOR THE NEW FFT.
              CIRC(IGO+NVAL-1) = 0.0
              CIRC(IGO+NVAL-2) = 0.0
              NVAL             = NVAL - 2
           ENDIF

           IF (MODE == 'H')  THEN
              LT = NVAL / 2
           ELSEIF (MODE.EQ.'F') THEN
              LT = NVAL / 4
           ENDIF

           LTIGO        = LT + IGO
           LTLTIGO      = LT + LT + IGO
           LTLTLTIGO    = LT + LT + LT + IGO

           NSIM         = LT - 1
           DFI          = PI / (NSIM+1)

C          TO AVOID SLOW BOUNDARY TESTS IN FBS2, PUT THEM HERE
           X            = CNS2
           Y            = INR + CNR2
           IF (X < 2.0 .OR. X  >  (FLOAT(NX)-1.0) .OR. 
     &         Y < 2.0 .OR. Y  >  (FLOAT(NY)-1.0) ) THEN
                  WRITE(NOUT,*) 'For image size1: ',NX,NY,INR
               WRITE(NOUT,90) X,Y
90             FORMAT('  FOR LOCATION: (',F7.1,',',F7.1,')')
               CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
C              RETURN NOT POSSIBLE WITH INTEL PARALLEL COMPILER
               STOP
           ENDIF

           CIRCT     = FBS2(X,Y, NXLD,NX,NY, 
     &                      XIM,NX, X1,Y1,XY2, .FALSE.)
           CIRC(IGO) = (CIRCT - AVO) * VRIN

C          TO AVOID SLOW BOUNDARY TESTS IN FBS2, PUT THEM HERE
           X            = INR + CNS2
           Y            =     + CNR2
           IF (X  <  2.0 .OR. X  >  (FLOAT(NX)-1.0) .OR. 
     &         Y  <  2.0 .OR. Y  >  (FLOAT(NY)-1.0) ) THEN
                  WRITE(NOUT,*) 'For image size2: ',NX,NY,INR
               WRITE(NOUT,90) X,Y
               CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
C              RETURN NOT POSSIBLE WITH INTEL PARALLEL COMPILER
               STOP
           ENDIF

           CIRCT       = FBS2(X,Y, NXLD,NX,NY, 
     &                        XIM,NX, X1,Y1,XY2, .FALSE.)
           CIRC(LTIGO) = (CIRCT - AVO) * VRIN

           IF (MODE .EQ. 'F')  THEN
C             FILL OTHER HALF OF CIRCLE

C             TO AVOID SLOW BOUNDARY TESTS IN FBS2, PUT THEM HERE
              X     = 0.0  + CNS2
              Y     = -INR + CNR2
              IF (X  <  2.0 .OR. X  >  (FLOAT(NX)-1.0) .OR. 
     &            Y  <  2.0 .OR. Y  >  (FLOAT(NY)-1.0) ) THEN
                  WRITE(NOUT,*) 'For image size3: ',NX,NY,INR
                  WRITE(NOUT,90) X,Y
                  CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
C                 RETURN NOT POSSIBLE WITH INTEL PARALLEL COMPILER
                  STOP
              ENDIF

              CIRCT         =  FBS2(X,Y, NXLD,NX,NY,
     &                              XIM,NX, X1,Y1,XY2, .FALSE.)
              CIRC(LTLTIGO) = (CIRCT - AVO) * VRIN

C             TO AVOID SLOW BOUNDARY TESTS IN FBS2, PUT THEM HERE
              X     = -INR + CNS2
              Y     =  0.0 + CNR2
              IF (X  <  2.0 .OR. X  >  (FLOAT(NX)-1.0) .OR. 
     &            Y  <  2.0 .OR. Y  >  (FLOAT(NY)-1.0) ) THEN
                  WRITE(NOUT,*) 'For image size4: ',NX,NY,INR
                  WRITE(NOUT,90) X,Y
                  CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
C                 RETURN NOT POSSIBLE WITH INTEL PARALLEL COMPILER
                  STOP
              ENDIF
              CIRCT           = FBS2(X,Y, NXLD,NX,NY,
     &                             XIM,NX, X1,Y1,XY2, .FALSE.)
              CIRC(LTLTLTIGO) = (CIRCT - AVO) * VRIN
           ENDIF

           DO JT=1,NSIM     ! LOOP NSIM TIMES TO FILL RING
              FI    = DFI * JT
              X     = SIN(FI) * YQ
              Y     = COS(FI) * YQ

              CIRCT        = FBS2(X+CNS2,Y+CNR2,  NXLD,NX,NY,
     &                            XIM,NX, X1,Y1,XY2, .FALSE.)
              CIRC(JT+IGO) = (CIRCT - AVO) * VRIN

              CIRCT        = FBS2(Y+CNS2,-X+CNR2, NXLD,NX,NY,
     &                            XIM,NX, X1,Y1,XY2, .FALSE.)
              CIRC(JT+LTIGO) = (CIRCT - AVO) * VRIN

              IF (MODE .EQ. 'F')  THEN
C                FILL OTHER HALF OF CIRCLE
                 CIRCT = FBS2(-X+CNS2,-Y+CNR2, NXLD,NX,NY,
     &                        XIM,NX, X1,Y1,XY2, .FALSE.)
                 CIRC(JT+LTLTIGO) = (CIRCT - AVO) * VRIN

                 CIRCT = FBS2(-Y+CNS2,X+CNR2, NXLD,NX,NY,
     &                        XIM,NX, X1,Y1,XY2, .FALSE.)
                 CIRC(JT+LTLTLTIGO) = (CIRCT - AVO) * VRIN
              ENDIF
	   ENDDO
	ENDDO


        ELSE ! -------------------------------- NOT IN OMP --------

C       FILL ALL THE RINGS 
        DO  IT=1,NRING 

           INR  = NUMR(1,IT)        ! RADIUS OF THE CURRENT RING
           YQ   = INR               ! FLOATING POINT RADIUS
           IGO  = NUMR(2,IT)        ! STARTING LOCATION FOR RING 
           NVAL = NUMR(3,IT)        ! LENGTH OF THIS RING

           IF (NEWFFT) THEN
C             THE ACTUAL, POWER-OF-TWO LENGTH IS NUMR(3,I)-2, ADDITIONAL
C             TWO LOCATIONS ARE ONLY FOR THE NEW FFT.
              CIRC(IGO+NVAL-1) = 0.0
              CIRC(IGO+NVAL-2) = 0.0
              NVAL             = NVAL - 2
           ENDIF

           IF (MODE .EQ. 'H')  THEN
              LT = NVAL / 2
           ELSEIF (MODE.EQ.'F') THEN
              LT = NVAL / 4
           ENDIF

           LTIGO        = LT + IGO
           LTLTIGO      = LT + LT + IGO
           LTLTLTIGO    = LT + LT + LT + IGO

           NSIM         = LT - 1
           DFI          = PI / (NSIM+1)

C          TO AVOID SLOW BOUNDARY TESTS IN FBS2, PUT THEM HERE
           X            = CNS2
           Y            = INR + CNR2

           IF (X  <  2.0 .OR. X  >  (FLOAT(NX)-1.0) .OR. 
     &         Y  <  2.0 .OR. Y  >  (FLOAT(NY)-1.0) ) THEN
               !WRITE(NOUT,*) 'For image size1: ',NX,NY,INR
               WRITE(NOUT,90) X,Y
               CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
               RETURN
           ENDIF

           CIRCT        = FBS2(X,Y, NXLD,NX,NY, 
     &                         XIM,NX, X1,Y1,XY2,.FALSE.)
           CIRC(IGO)    = (CIRCT - AVO) * VRIN

C          TO AVOID SLOW BOUNDARY TESTS IN FBS2, PUT THEM HERE
           X            = INR + CNS2
           Y            =     + CNR2
           IF (X  <  2.0 .OR. X  >  (FLOAT(NX)-1.0) .OR. 
     &         Y  <  2.0 .OR. Y  >  (FLOAT(NY)-1.0) ) THEN
               !WRITE(NOUT,*) 'For image size2: ',NX,NY,INR
               WRITE(NOUT,90) X,Y
               CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
               STOP
           ENDIF

           CIRCT        = FBS2(X,Y, NXLD,NX,NY,
     &                         XIM,NX, X1,Y1,XY2, .FALSE.)
           CIRC(LTIGO)  = (CIRCT - AVO) * VRIN

           IF (MODE .EQ. 'F')  THEN
C             FILL OTHER HALF OF CIRCLE

C             TO AVOID SLOW BOUNDARY TESTS IN FBS2, PUT THEM HERE
              X     = 0.0  + CNS2
              Y     = -INR + CNR2
              IF (X  <  2.0 .OR. X  >  (FLOAT(NX)-1.0) .OR. 
     &            Y  <  2.0 .OR. Y  >  (FLOAT(NY)-1.0) ) THEN
                  !WRITE(NOUT,*) 'For image size3: ',NX,NY,INR
                  WRITE(NOUT,90) X,Y
                  CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
                  STOP
              ENDIF

              CIRCT =  FBS2(X,Y, NXLD,NX,NY, 
     &                      XIM,NX, X1,Y1,XY2,.FALSE.)

              CIRC(LTLTIGO) = (CIRCT - AVO) * VRIN

C             TO AVOID SLOW BOUNDARY TESTS IN FBS2, PUT THEM HERE
              X     = -INR + CNS2
              Y     =  0.0 + CNR2
              IF (X  <  2.0 .OR. X  >  (FLOAT(NX)-1.0) .OR. 
     &            Y  <  2.0 .OR. Y  >  (FLOAT(NY)-1.0) ) THEN
                  !WRITE(NOUT,*) 'For image size4: ',NX,NY,INR
                  WRITE(NOUT,90) X,Y
                  CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
                  STOP
              ENDIF
              CIRCT = FBS2(X,Y, NXLD,NX,NY, 
     &                     XIM,NX, X1,Y1,XY2, .FALSE.)
              CIRC(LTLTLTIGO) = (CIRCT - AVO) * VRIN
           ENDIF
          
           DO JT=1,NSIM     ! LOOP NSIM TIMES TO FILL RING
              FI = DFI * JT
              X  = SIN(FI) * YQ
              Y  = COS(FI) * YQ

              XT    = X+CNS2
              YT    = Y+CNR2
              CIRCT = FBS2(XT,YT, NXLD,NX,NY, 
     &                     XIM,NX, X1,Y1,XY2,.FALSE.)
              CIRC(JT+IGO) = (CIRCT - AVO) * VRIN
 
              XT =  Y+CNS2
              YT = -X+CNR2

              CIRCT = FBS2(XT,YT, NXLD,NX,NY, 
     &                     XIM,NX, X1,Y1,XY2,.FALSE.)
              CIRC(JT+LTIGO) = (CIRCT - AVO) * VRIN

              IF (MODE .EQ. 'F')  THEN
C                FILL OTHER HALF OF CIRCLE
                 XT    = -X+CNS2
                 YT    = -Y+CNR2
                 CIRCT = FBS2(XT,YT, NXLD,NX,NY,
     &                        XIM,NX, X1,Y1,XY2, .FALSE.)
                 CIRC(JT+LTLTIGO) = (CIRCT - AVO) * VRIN

                 XT    = -Y+CNS2
                 YT    =  X+CNR2
                 CIRCT = FBS2(XT,YT, NXLD,NX,NY,
     &                        XIM,NX, X1,Y1,XY2, .FALSE.)
                 CIRC(JT+LTLTLTIGO) = (CIRCT - AVO) * VRIN
              ENDIF
	   ENDDO
	ENDDO

        ENDIF

9999    IF (ALLOCATED(F0))  DEALLOCATE (F0)
        IF (ALLOCATED(X1))  DEALLOCATE (X1)
        IF (ALLOCATED(Y1))  DEALLOCATE (Y1)
        IF (ALLOCATED(XY2)) DEALLOCATE (XY2)

        END
