
C++*********************************************************************
C
C APRINGS_TRANS     ADAPTED FOR TRANS             JUN 10 ARDEAN LEITH
C **********************************************************************
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* Authors: J. Frank & A. Leith                                        *
C=* Copyright 1985-2010  Health Research Inc.                          *
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
C  APRINGS_TRANS_REFORM
C
C  PURPOSE: CONVERTS RING ORDER 'RINGS' ARRAY TO A RAY ORDER ARRAY
C
C  PARAMETERS: 
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

cpgi$g opt=O3

      SUBROUTINE APRINGS_TRANS_REFORM(NUMR,NRING,  NLOCS,NRAYSC,
     &                                COLD,LCIRCOD2, CNEW,LCIRCND2, WR)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'

      INTEGER, INTENT(IN)      :: NUMR(3,NRING)
      INTEGER, INTENT(IN)      :: NRING
      INTEGER, INTENT(IN)      :: NLOCS(2,NRAYSC+1)
      INTEGER, INTENT(IN)      :: NRAYSC  
      INTEGER, INTENT(IN)      :: LCIRCOD2,LCIRCND2
      COMPLEX, INTENT(IN)      :: COLD(LCIRCOD2)
      COMPLEX, INTENT(OUT)     :: CNEW(LCIRCND2)
      REAL,    INTENT(IN)      :: WR(NRING)

      INTEGER                  :: MAXRIN,J,IRING,NR,IRAY,IT,ILOC,I
      LOGICAL                  :: WEIGHT
      REAL                     :: WRT

      WEIGHT = (WR(1) .NE. 0)            ! WEIGHTING FLAG
      MAXRIN = NUMR(3,NRING)/2           ! # CMPLX RAYS ON LONGEST RING
      J      = 0                         ! ZERO COLD INDEX
      DO IRING = 1,NRING                 ! LOOP OVER ALL RINGS

         NR  = NUMR(3,IRING) / 2         ! # CMPLX RAYS ON THIS RING
         WRT = WR(IRING)                 ! WEIGHT FOR THIS RING

         DO IRAY = 1,NR                  ! LOOP OVER RAYS ON RING
            J   = J + 1                  ! COLD INDEX

            IT  = IRING - (NRING - NLOCS(1,IRAY))

C           NLOCS IS INDEXED FOR COMPLEX NUMBER PAIRS
            ILOC = NLOCS(2,IRAY)+IT-1        ! CNEW INDEX

            IF (.NOT. WEIGHT) THEN
               CNEW(ILOC) = COLD(J)          ! MOVE
            ELSE
               CNEW(ILOC) = COLD(J)  * WRT   ! MOVE & WEIGHT
            ENDIF

            !if(iray == 1 )write(6,90)nr,it, j,iloc, cold(j) 
            !if(iring == 1 .and. iray == 1)write(6,90)nr,it,j,iloc,cold(j)
90          format(i5,i6,i6,' -->',i6,' (',f10.2,',',f10.2,')')

         ENDDO   ! END OF: DO IRAY  = 1,NR

         IF (WEIGHT .AND. NR .NE. MAXRIN) THEN
C           IF RING LENGTH IS LESS THAN THE MAX. RING LENGTH THEN
C           WEIGHT OF FINAL FFT COEF FOR THIS RING IS HALF USUAL WEIGHT.
C           WHY? al
            !write(6,*) ' wr:',nr,maxrin,j,iloc
            CNEW(ILOC) = CNEW(ILOC) * 0.5 
         ENDIF

       ENDDO     ! END OF: DO IRING = 1, NRING

       END

C     --------------------- APRINGS_TRANS_LOCS ------------------------

      SUBROUTINE APRINGS_TRANS_LOCS(NUMR,NRING, NLOCS,NRAYSC)

      IMPLICIT NONE

      INTEGER, INTENT(IN )     :: NRING,NRAYSC
      INTEGER, INTENT(IN)      :: NUMR(3,NRING)
      INTEGER, INTENT(OUT)     :: NLOCS(2,NRAYSC+1)

      INTEGER                  :: ILOC,IRING,NVT,IRAY

      NVT   = 0                        ! NUMBER OF VALUES TO CC
      NLOCS = 0                        ! ZERO WHOLE NLOC ARRAY

C     FIND NUMBER OF RINGS ON EACH RAY (1...NRAYSC)
      DO ILOC = 1,NRAYSC               ! LOOP OVER RAYS
         DO IRING = 1,NRING            ! LOOP OVER RINGS
            IF (ILOC .LE. (NUMR(3,IRING) / 2))
     &         NLOCS(1,ILOC) = NLOCS(1,ILOC)+1
         ENDDO
         !write(6,'i4,": ",i5,i6') iloc,nloc(iloc,1)
         NVT = NVT + NLOCS(1,ILOC)
      ENDDO

      !WRITE(6,97)'NUMBER OF RINGS:',MAXVAL(NLOCS(1,:))
      !WRITE(6,97)'NUMBER OF RING POINTS:',NVT
 97   FORMAT('  ',A,T30,I10)

C     FIND STARTING LOCATION OF EACH RAY (1...MAXRAY) AND FINAL RAY LEN.
      ILOC = 1 

      DO IRAY = 1,NRAYSC
         !IM = MOD(ILOC-1,2)
         !IF (IM == 0) THEN
         !NLOCS(2,IRAY) = ILOC
         !ELSE
         !IPAD = 2 - IM
         !NLOCS(2,IRAY) = ILOC + IPAD
         !ENDIF

         NLOCS(2,IRAY) = ILOC
         ILOC          = ILOC + NLOCS(1,IRAY) ! NEXT RAY START INDEX
      ENDDO
      NLOCS(2,NRAYSC+1  )= ILOC               ! TO FIND LENGTH OF LAST RA

      !write (6,*)'  nlocs(1,:) -- (# of rings on ray) --nraysc:',nraysc
      !write(6,*) nlocs(1,:)
      !write (6,97) '-------------------------'
      !write (6,*) '  nlocs(2,:) -- (starting index) ---------'
      !write(6,*) nlocs(2,:)
      !write (6,97) '-------------------------'

      END


C     --------------------- APRINGS_TRANS_ONE ------------------------

C       PURPOSE: USED FOR TRANSFORMED DATA WITHOUT COEF
C                USED WHEN ONLY ONE USAGE OF COEF WOULD BE INEFFICIENT

        SUBROUTINE APRINGS_TRANS_ONE(XIM,  NSAM,NROW, CNS2,CNR2, 
     &                       NUMR,NRING,   NLOCS,NRAYSC,
     &                       MODE,USE_OMP, WR, FFTW_PLANS,
     &                       CIRC,LCIRC,   CIRCN,LCIRCC)

        IMPLICIT  NONE

        REAL,        INTENT(IN)      :: XIM(NSAM,NROW)
        INTEGER,     INTENT(IN)      :: NSAM,NROW
        REAL,        INTENT(IN)      :: CNS2,CNR2
        INTEGER,     INTENT(IN)      :: NUMR(3,NRING),NRING
        INTEGER,     INTENT(IN)      :: NLOCS(2,NRAYSC+1),NRAYSC
        CHARACTER*1, INTENT(IN)      :: MODE
        LOGICAL,     INTENT(IN)      :: USE_OMP
        REAL,        INTENT(IN)      :: WR(*)
        INTEGER*8,   INTENT(IN)      :: FFTW_PLANS(*) ! STRUCTURE POINTERS
        REAL,        INTENT(OUT)     :: CIRC(LCIRC)   ! WORK ARRAY
        INTEGER,     INTENT(IN)      :: LCIRC
        COMPLEX,     INTENT(OUT)     :: CIRCN(LCIRCC)
        INTEGER,     INTENT(IN)      :: LCIRCC


        DOUBLE PRECISION             :: AVO,VRINV
        DOUBLE PRECISION             :: PI,DFI
        REAL                         :: CIRCT
        INTEGER                      :: NSB,NSE,NRB,NRE,IT,INR,IGO,NVAL
        INTEGER                      :: LT,LTIGO,LTLTIGO,LTLTLTIGO,NSIM
        INTEGER                      :: NE,JT,INDX,IRTFLG
        REAL                         :: YQ,X,Y,FI,XT,YT

        REAL                         :: QUADRI_FAST
        INTEGER                      :: LOG2

        LOGICAL, PARAMETER           :: SPIDER_SIGN = .FALSE.

        INCLUDE 'CMBLOCK.INC'

C       FIND PARAMETERS TO NORMALIZE UNDER THE MASK,  
C       TRIED DOING THIS COMPLETELY ON THE
C       POLAR RINGS BUT IT GIVES SOME DIFFERENT REF. CHOICES. al

C       CNS2 AND CNR2 ARE PREDEFINED CENTERS
C       CALCULATE DIMENSIONS FOR NORMALIZING MASK
	NSB  = -CNS2
	NSE  =  NSB + NSAM - 1
	NRB  = -CNR2
	NRE  =  NRB + NROW - 1

C       GET PARAMETERS AVO & VRINV TO NORMALIZE UNDER CIRCULAR MASK
        CALL NORMASC(XIM, NSB,NSE,NRB,NRE, NUMR,NUMR(1,NRING),
     &               AVO,VRINV,USE_OMP)

C       INTERPOLATE INTO POLAR COORDINATES & APPLY NORMALIZATION
C       CREATING CIRC (RADIAL IMAGE CIRCLES) FOR THIS IMAGE POSITION
C       FOURIER TRANSFORM CIRC RINGS
C       OPTIONAL: WEIGHT TRANSFORMED CIRC RINGS USING  FACTORS FROM WR

        PI = 2 * DATAN(1.0D0)

C       FILL ALL THE RINGS

        IF (USE_OMP) THEN

c$omp      parallel do private(it,inr,yq,igo,nval,lt,ltigo,ltltigo,
c$omp&                     ltltltigo,nsim,dfi,x,y,circt,jt,fi,indx)
           DO  IT=1,NRING 

              INR  = NUMR(1,IT)        ! RADIUS OF THE CURRENT RING
              YQ   = INR               ! FLOATING POINT RADIUS
              IGO  = NUMR(2,IT)        ! STARTING LOCATION FOR RING 
              NVAL = NUMR(3,IT)        ! LENGTH OF THIS RING

C             ACTUAL, POWER-OF-TWO LENGTH IS NUMR(3,I)-2, ADDITIONAL
C             TWO LOCATIONS ARE ONLY FOR THE NEW FFT.
              CIRC(IGO+NVAL-1) = 0.0
              CIRC(IGO+NVAL-2) = 0.0
              NVAL             = NVAL - 2

              IF (MODE .EQ. 'H')    THEN
                 LT = NVAL / 2
              ELSEIF (MODE .EQ.'F') THEN
                 LT = NVAL / 4
              ENDIF

              LTIGO     = LT + IGO
              LTLTIGO   = LT + LT + IGO
              LTLTLTIGO = LT + LT + LT + IGO

              NSIM      = LT - 1
              DFI       = PI / (NSIM+1)

C             AVOID SLOW BOUNDARY TESTS IN QUADRI_FAST, PUT THEM HERE
              X         = CNS2
              Y         = INR + CNR2
              IF (X .LE. 2.0 .OR. X .GE. (FLOAT(NSAM)-1.0) .OR. 
     &            Y .LE. 2.0 .OR. Y .GE. (FLOAT(NROW)-1.0) ) THEN
                 WRITE(NOUT,*) 'For image size1: ',NSAM,NROW,INR
                 WRITE(NOUT,90) X,Y
90               FORMAT('  FOR LOCATION: (',F7.1,',',F7.1,')')
                 CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
C                RETURN NOT POSSIBLE WITH INTEL PARALLEL COMPILER
                 STOP
              ENDIF

              CIRCT     = QUADRI_FAST(X,Y,NSAM,NROW,XIM)
              CIRC(IGO) = (CIRCT - AVO) * VRINV

C             AVOID SLOW BOUNDARY TESTS IN QUADRI_FAST, PUT THEM HERE
              X  = INR + CNS2
              Y  =     + CNR2
              IF (X .LE. 2.0 .OR. X .GE. (FLOAT(NSAM)-1.0) .OR. 
     &            Y .LE. 2.0 .OR. Y .GE. (FLOAT(NROW)-1.0) ) THEN
                 WRITE(NOUT,*) 'For image size2: ',NSAM,NROW,INR
                 WRITE(NOUT,90) X,Y
                 CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
C                RETURN NOT POSSIBLE WITH INTEL PARALLEL COMPILER
                 STOP
              ENDIF

              CIRCT       = QUADRI_FAST(X,Y,NSAM,NROW,XIM)
              CIRC(LTIGO) = (CIRCT - AVO) * VRINV

              IF (MODE .EQ. 'F')  THEN
C                FILL OTHER HALF OF CIRCLE

C                TO AVOID SLOW BOUNDARY TESTS IN QUADRI_FAST, PUT THEM HERE
                 X  = 0.0  + CNS2
                 Y  = -INR + CNR2
                 IF (X .LE. 2.0 .OR. X .GE. (FLOAT(NSAM)-1.0) .OR. 
     &               Y .LE. 2.0 .OR. Y .GE. (FLOAT(NROW)-1.0) ) THEN
                    WRITE(NOUT,*) 'For image size3: ',NSAM,NROW,INR
                    WRITE(NOUT,90) X,Y
                    CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
C                   RETURN NOT POSSIBLE WITH INTEL PARALLEL COMPILER
                    STOP
                 ENDIF

                 CIRCT         =  QUADRI_FAST(X,Y,NSAM,NROW,XIM)
                 CIRC(LTLTIGO) = (CIRCT - AVO) * VRINV

C                AVOID SLOW BOUNDARY TESTS IN QUADRI_FAST, PUT THEM HERE
                 X = -INR + CNS2
                 Y =  0.0 + CNR2
                 IF (X .LE. 2.0 .OR. X .GE. (FLOAT(NSAM)-1.0) .OR. 
     &               Y .LE. 2.0 .OR. Y .GE. (FLOAT(NROW)-1.0) ) THEN
                    WRITE(NOUT,*) 'For image size4: ',NSAM,NROW,INR
                    WRITE(NOUT,90) X,Y
                    CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
C                   RETURN NOT POSSIBLE WITH INTEL PARALLEL COMPILER
                    STOP
                 ENDIF
                 CIRCT           = QUADRI_FAST(X,Y,NSAM,NROW,XIM)
                 CIRC(LTLTLTIGO) = (CIRCT - AVO) * VRINV
              ENDIF

              DO JT=1,NSIM     ! LOOP NSIM TIMES TO FILL RING
                 FI           = DFI * JT
                 X            = SIN(FI) * YQ
                 Y            = COS(FI) * YQ

                 CIRCT        = QUADRI_FAST(X+CNS2,Y+CNR2,NSAM,NROW,XIM)
                 CIRC(JT+IGO) = (CIRCT - AVO) * VRINV

                 CIRCT        = 
     &                         QUADRI_FAST(Y+CNS2,-X+CNR2,NSAM,NROW,XIM)
                 CIRC(JT+LTIGO) = (CIRCT - AVO) * VRINV

                 IF (MODE .EQ. 'F')  THEN
C                   FILL OTHER HALF OF CIRCLE
                    CIRCT = QUADRI_FAST(-X+CNS2,-Y+CNR2,NSAM,NROW,XIM)
                    CIRC(JT+LTLTIGO) = (CIRCT - AVO) * VRINV

                    CIRCT = QUADRI_FAST(-Y+CNS2,X+CNR2,NSAM,NROW,XIM)
                    CIRC(JT+LTLTLTIGO) = (CIRCT - AVO) * VRINV
                 ENDIF
	      ENDDO

C             INPLACE FORWARD FOURIER TRANSFORM ON THIS CIRC RING
              INDX = LOG2(NVAL) - 1       ! INDEX FOR PLAN
              
              CALL FMRS(CIRC(IGO), NVAL,1,1, FFTW_PLANS(INDX),
     &                  SPIDER_SIGN,.FALSE.,+1,IRTFLG)
              !write(6,*) ' fftd circ(1):',circ(1)
	   ENDDO
       ELSE

C          FILL ALL THE RINGS
           DO IT=1,NRING 

              INR  = NUMR(1,IT)        ! RADIUS OF THE CURRENT RING
              YQ   = INR               ! FLOATING POINT RADIUS
              IGO  = NUMR(2,IT)        ! STARTING LOCATION FOR RING 
              NVAL = NUMR(3,IT)        ! LENGTH OF THIS RING

C             THE ACTUAL, POWER-OF-TWO LENGTH IS NUMR(3,I)-2, ADDED
C             TWO LOCATIONS ARE ONLY FOR THE NEW FFT.
              CIRC(IGO+NVAL-1) = 0.0
              CIRC(IGO+NVAL-2) = 0.0
              NVAL             = NVAL - 2

              IF (MODE .EQ. 'H')  THEN
                 LT = NVAL / 2
              ELSEIF (MODE.EQ.'F') THEN
                 LT = NVAL / 4
              ENDIF

              LTIGO      = LT + IGO
              LTLTIGO    = LT + LT + IGO
              LTLTLTIGO  = LT + LT + LT + IGO

              NSIM       = LT - 1
              DFI        = PI / (NSIM+1)

C             AVOID SLOW BOUNDARY TESTS IN QUADRI_FAST, PUT THEM HERE
              X  = CNS2
              Y  = INR + CNR2
              IF (X .LE. 2.0 .OR. X .GE. (FLOAT(NSAM)-1.0) .OR. 
     &            Y .LE. 2.0 .OR. Y .GE. (FLOAT(NROW)-1.0) ) THEN
                 !WRITE(NOUT,*) 'For image size1: ',NSAM,NROW,INR
                 WRITE(NOUT,90) X,Y
                 CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
                 RETURN
              ENDIF

              CIRCT     = QUADRI_FAST(X,Y,NSAM,NROW,XIM)
              CIRC(IGO) = (CIRCT - AVO) * VRINV
              !write(6,*)it,':1',igo,circ(igo)

C             AVOID SLOW BOUNDARY TESTS IN QUADRI_FAST, PUT THEM HERE
              X  = INR + CNS2
              Y  =     + CNR2
              IF (X .LE. 2.0 .OR. X .GE. (FLOAT(NSAM)-1.0) .OR. 
     &            Y .LE. 2.0 .OR. Y .GE. (FLOAT(NROW)-1.0) ) THEN
                 !WRITE(NOUT,*) 'For image size2: ',NSAM,NROW,INR
                 WRITE(NOUT,90) X,Y
                 CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
                 STOP
              ENDIF

              CIRCT       = QUADRI_FAST(X,Y,NSAM,NROW,XIM)
              CIRC(LTIGO) = (CIRCT - AVO) * VRINV
              !write(6,*)it,':',ltigo,circ(ltigo)

              IF (MODE .EQ. 'F')  THEN
C                FILL OTHER HALF OF CIRCLE

C                AVOID SLOW BOUNDARY TESTS IN QUADRI_FAST, PUT THEM HERE
                 X  = 0.0  + CNS2
                 Y  = -INR + CNR2
                 IF (X .LE. 2.0 .OR. X .GE. (FLOAT(NSAM)-1.0) .OR. 
     &               Y .LE. 2.0 .OR. Y .GE. (FLOAT(NROW)-1.0) ) THEN
                    !WRITE(NOUT,*) 'For image size3: ',NSAM,NROW,INR
                    WRITE(NOUT,90) X,Y
                    CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
                    STOP
                 ENDIF
  
                 CIRCT         =  QUADRI_FAST(X,Y,NSAM,NROW,XIM)
                 CIRC(LTLTIGO) = (CIRCT - AVO) * VRINV
                 !write(6,*)it,':',ltltigo,circ(ltltigo)

C                VOID SLOW BOUNDARY TESTS IN QUADRI_FAST, PUT THEM HERE
                 X  = -INR + CNS2
                 Y  =  0.0 + CNR2
                 IF (X .LE. 2.0 .OR. X .GE. (FLOAT(NSAM)-1.0) .OR. 
     &               Y .LE. 2.0 .OR. Y .GE. (FLOAT(NROW)-1.0) ) THEN
                    !WRITE(NOUT,*) 'For image size4: ',NSAM,NROW,INR
                    WRITE(NOUT,90) X,Y
                    CALL ERRT(101,'RING GOES OUTSIDE IMAGE',NE)
                    STOP
                 ENDIF
                 CIRCT           = QUADRI_FAST(X,Y,NSAM,NROW,XIM)
                 CIRC(LTLTLTIGO) = (CIRCT - AVO) * VRINV
                 !write(6,*)it,':',ltltltigo,circ(ltltltigo)

              ENDIF
          
              DO JT=1,NSIM     ! LOOP NSIM TIMES TO FILL RING
                 FI             = DFI * JT
                 X              = SIN(FI) * YQ
                 Y              = COS(FI) * YQ

                 XT             = X+CNS2
                 YT             = Y+CNR2
                 CIRCT          = QUADRI_FAST(XT,YT,NSAM,NROW,XIM)
                 CIRC(JT+IGO)   = (CIRCT - AVO) * VRINV
 
                 XT             =  Y+CNS2
                 YT             = -X+CNR2

                 CIRCT          = QUADRI_FAST(XT,YT,NSAM,NROW,XIM)
                 CIRC(JT+LTIGO) = (CIRCT - AVO) * VRINV

                 IF (MODE .EQ. 'F')  THEN
C                   FILL OTHER HALF OF CIRCLE
                    XT                 = -X+CNS2
                    YT                 = -Y+CNR2
                    CIRCT              = QUADRI_FAST
     &                                   (XT,YT,NSAM,NROW,XIM)
                    CIRC(JT+LTLTIGO)   = (CIRCT - AVO) * VRINV

                    XT                 = -Y+CNS2
                    YT                 =  X+CNR2

                    CIRCT              = QUADRI_FAST
     &                                   (XT,YT,NSAM,NROW,XIM)
                    CIRC(JT+LTLTLTIGO) = (CIRCT - AVO) * VRINV
                 ENDIF
	      ENDDO

C             INPLACE FORWARD FOURIER TRANSFORM ON THIS CIRC RING
              INDX = LOG2(NVAL) - 1       ! INDEX FOR PLAN
              
              CALL FMRS(CIRC(IGO), NVAL,1,1, FFTW_PLANS(INDX),
     &                  SPIDER_SIGN,.FALSE.,+1,IRTFLG)
              !write(6,*) ' fftd circ(1):',circ(1)

 	   ENDDO
        ENDIF

C       REFORM RINGS BY RAY ORDER AND APPLY OPTIONAL WEIGHTING
        CALL APRINGS_TRANS_REFORM(NUMR,NRING, NLOCS,NRAYSC,  
     &                            CIRC,LCIRC/2, CIRCN,LCIRCC, WR)

        END

