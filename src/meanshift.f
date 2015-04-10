
C ++********************************************************************
C                                                                      *
C MEANSHIFT                                                            *
C                                                                      *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NROW 12204.      *
C=* Email: spider@wadsworth.org                                        *
C=*                                                                    *
C=* SPIDER is free software; you can redistribute it and/or            *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) aNROW later version.                    *
C=*                                                                    *
C=* SPIDER is distributed in the hope that it will be useful,          *
C=* but WITHOUT ANROW WARRANTY; without even the implied warranty of     *
C=* merchantability or fitness for a particular purpose.  See the GNU  *
C=* General Public License for more details.                           *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program. If not, see <http://www.gnu.org/licenses> *
C=*                                                                    *
C **********************************************************************
C   PURPOSE: MEAN SHIFT DENOISING AND SMOOTHING FILTER (2D/3D FILES)   *
C                                                                      *
C      THE MEAN SHIFT FILTER INTRODUCED FIRST BY COMANICIU AND MEER    *
C      (2002) IS A DYNAMIC NONLINEAR FILTER, THAT ACHIEVES A HIGH      *
C      QUALITY EDGE-PRESERVING FILTERING                               *
C                                                                      *
C         USER:                                                        *
C     1. SET THE RADIUS OF CIRCULAR OR SPHERICAL KERNEL IN PIXELS.     *
C     2. SET THE VALUE OF DENSITY DISTANCE.                            *
C         ~ THE KERNEL RADIUS IN THE RANGE 3-7 PIXELS AND THE DENSITY  *
C        DENSITY DISTANCE IN THE RANGE 1/2-2 STANDARD DEVIATION ARE    *
C        RECOMMENDED                                                   *
C                                                                      *
C         ALGORITHM:                                                   *
C     1. START AT A POINT 'A'.                                         *
C     2. SELECT THE PIXELS (VOXELS) FALLING WITHIN A CERTAIN           *
C         SPATIAL DISTANCE (KERNEL) AND WITHIN CERTAIN DENSITY         *
C         DISTANCE.                                                    *
C     3. CALCULATE THE CENTER OF MASS OF THE SET OF SELECTED PIXELS    *
C        (VOXELS), AND DISPLACE THE CENTER OF KERNEL TO THE CENTER     *
C        OF MASS.                                                      *
C     4. REPEAT ITERATIVELY UNTIL THE SPATIAL VARIATION IS LOWER       *
C        THAN 2/3 OF AKERNEL RADIUS                                    *                                           *
C     5. APPLY THE MEAN DENSITY OF SELECTED PIXELS AROUND THE FINAL    *
C        POINT 'B' (THE MODE) TO STARTING POINT 'A'.                   *
C        ~ CALCULATION OF MEAN DENSITY INSIDE KERNEL INVOLVES ONLY     *
C        PIXELS/VOXELS WITHIN AN APPOINTED DENSITY DISTANCE            *
C        ~ IN SOME OTHER SOFTWARE THE DENSITY OF FINAL POINT 'B'       *
C                      ITSELF IS USED.                                 *
C     6. REPEAT FOR ALL PIXELS (VOXELS) IN THE IMAGE.                  *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

       SUBROUTINE MEANSHIFT(LUN1,LUN2,NSAM,NROW,NSLICE,SIG1,IRTFLG)

       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC'

       INTEGER               :: LUN1,LUN2
       INTEGER               :: NSAM,NROW,NSLICE
       REAL                  :: SIG1
       INTEGER               :: IRTFLG

       INTEGER               :: KERNEL
       REAL, ALLOCATABLE     :: BUF1(:,:,:),BUF2(:,:,:)
       REAL                  :: DFACTOR,GRAD
       LOGICAL               :: erri2


       KERNEL = 3
       CALL RDPRI1S(KERNEL,NOT_USED,'RADIUS',IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 9999
       IF (ERRI2(KERNEL,IDUM,1,1,30,0,0)) GOTO 9999
       
       DFACTOR = 1.0
       CALL RDPRM1S(DFACTOR,NOT_USED,
     &       'DENSITY DISTANCE FACTOR (OR <CR> FOR 1)',IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 9999
       IF (DFACTOR .LE. 0) THEN
          CALL ERRT(101,'FACTOR MUST BE > 0',NE)
          GOTO 9999
       ENDIF

       GRAD = SIG1 * DFACTOR

       ALLOCATE (BUF1(NSAM,NROW,NSLICE), 
     &           BUF2(NSAM,NROW,NSLICE),STAT=IRTFLG)
       IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'MEANSHIFT; BUF1,BUF2',2*NSAM*NROW*NSLICE)
           GOTO 9999
       ENDIF

       CALL REDVOL(LUN1,NSAM,NROW,1,NSLICE,BUF1,IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 9999

       IF (NSLICE < 2) THEN
          CALL MEANSHIFT2(BUF1,BUF2,NSAM,NROW,
     &                    SIG1,KERNEL,GRAD,IRTFLG)
       ELSE
          CALL MEANSHIFT3(BUF1,BUF2,NSAM,NROW,NSLICE,
     &                    SIG1,KERNEL,GRAD,IRTFLG)
       ENDIF
       
       CALL WRTVOL(LUN2,NSAM,NROW,1,NSLICE,BUF2,IRTFLG)

9999   IF (ALLOCATED(BUF1))   DEALLOCATE(BUF1)                          
       IF (ALLOCATED(BUF2))   DEALLOCATE(BUF2)                          
       CLOSE(LUN1)
       CLOSE(LUN2)
            
       END

  
C      ******************************* MEANSHIFT2 *******************

       SUBROUTINE MEANSHIFT2(BUF1,BUF2,NX,NY,
     &                        SIG1,KERNEL,GRAD,IRTFLG)

       INTEGER :: J, I, J2,I2
       INTEGER    NX,NY,NJ,NI
       INTEGER :: JTMP,ITMP
       INTEGER :: KERNEL
       INTEGER :: JCNTR, ICNTR
       INTEGER :: NN,  TEMP
       INTEGER :: WEIGHT
       INTEGER :: FILT(-30:30,-30:30)
       REAL    :: BUF1(NX,NY)
       REAL    :: BUF2(NX,NY)
       REAL    :: GRAD
       REAL    :: RJCNTR,  RICNTR
       REAL    :: SIG1
       REAL    :: DENS, INCNTR
       REAL    :: RAD

       !WRITE(6,*) 'SIG1 =',SIG1
       !WRITE(6,*) 'NY,NX =',NY,NX
       !WRITE(6,*) 'KERNEL =',KERNEL

C      CALCULATE THE CIRCULAR AREA OF KERNEL
C      INVOLVING A ROUNDING OF BOUNDARIES

       DO I=-KERNEL,KERNEL
          DO J=-KERNEL,KERNEL
             RAD = I**2 + J**2
             RAD = SQRT(RAD)
             RAD = INT(RAD+0.5)
             IF (RAD <= KERNEL) THEN
	        FILT(J,I) = 1
             ELSE
	        FILT(J,I) = 0		   
             ENDIF
          ENDDO
       ENDDO

! SIMPLE OMP IS SLOWER!!
!c$omp  parallel do private(i,j,ni,nj,nn,icntr,jcntr,weight,dens,itmp,
!c$omp&                     jtmp,incntr,i2,j2,temp)  
       DO I=1,NY
          DO J=1,NX
               NI     = I
               NJ     = J
	       NN     = 0

14             ICNTR  = 0
               JCNTR  = 0
               WEIGHT = 0
 	       DENS   = 0
               ITMP   = 1 + MODULO(NI-1,NY)
               JTMP   = 1 + MODULO(NJ-1,NX)
               INCNTR = BUF1(JTMP,ITMP)

               DO I2=NI-KERNEL, NI+KERNEL
                  ITMP = 1 + MODULO(I2-1,NY)

                  DO J2=NJ-KERNEL, NJ+KERNEL
                     JTMP = 1 + MODULO(J2-1,NX)

                    IF (ABS(BUF1(JTMP,ITMP)-INCNTR) < GRAD) THEN
                       ICNTR  = ICNTR  + I2*FILT(J2-NJ,I2-NI)
                       JCNTR  = JCNTR  + J2*FILT(J2-NJ,I2-NI)
                       WEIGHT = WEIGHT +    FILT(J2-NJ,I2-NI)
                       DENS   = DENS   + BUF1(JTMP,ITMP) *
     &                                   FILT(J2-NJ,I2-NI)
                    ENDIF
                ENDDO
             ENDDO
             NN    = NN+1
             ICNTR = INT(REAL(ICNTR) / REAL(WEIGHT))
             JCNTR = INT(REAL(JCNTR) / REAL(WEIGHT))
             TEMP  = ABS(JCNTR-NJ) + ABS(ICNTR-NI)

             IF (3*TEMP > 2*KERNEL) THEN
                
                NI = ICNTR
                NJ = JCNTR
                IF (NN < 10) GOTO 14
             ENDIF

             BUF2(J,I) = DENS / REAL(WEIGHT)
          ENDDO
       ENDDO
       
       END

C      ******************************* MEANSHIFT3 *******************

       SUBROUTINE MEANSHIFT3(BUF1,BUF2,NX,NY,NZ,
     &                        SIG1,KERNEL,GRAD,IRTFLG)

       INTEGER :: S, J, I, S2,J2,I2
       INTEGER :: NZ,NX,NY,NS,NJ,NI
       INTEGER :: STMP,JTMP,ITMP
       INTEGER :: KERNEL
       INTEGER :: SCNTR, JCNTR, ICNTR
       INTEGER :: NN,  TEMP
       INTEGER :: WEIGHT
       INTEGER :: FILT(-30:30,-30:30,-30:30)
       REAL    :: BUF1(NX,NY,NZ)
       REAL    :: BUF2(NX,NY,NZ)
       REAL    :: GRAD
       REAL    :: RSCNTR,  RJCNTR,  RICNTR
       REAL    :: SIG1
       REAL    :: DENS, INCNTR
       REAL    :: RAD
	 
       !WRITE(6,*) 'SIG1 =',SIG1
       !WRITE(6,*) 'NX,NY,NZ =',NX,NY,NZ
       !WRITE(6,*) 'KERNEL =',KERNEL

C      CALCULATE THE SPHERICAL AREA OF KERNEL
C       INVOLVING A ROUNDING OF BOUNDARIES

       DO S=-KERNEL,KERNEL
           DO I=-KERNEL,KERNEL
              DO J=-KERNEL,KERNEL
                 RAD = S**2 + I**2 + J**2
                 RAD = INT(RAD + 0.5)
                 IF (RAD <= KERNEL) THEN
	            FILT(J,I,S) = 1
                 ELSE
	            FILT(J,I,S) = 0		   
                 ENDIF
              ENDDO
           ENDDO
       ENDDO

! SIMPLE OMP IS SLOWER!!
!c$omp  parallel do private(s,i,j, ns,ni,nj, nn, rscntr,ricntr,rjcntr
!c$omp&                     weight,dens,stmp,itmp,jtmp,incntr,s2,i2,j2,
!c$omp&                     temp)  

       DO S=1,NZ
          DO I=1,NY

             DO J=1,NX
               NS     = S
               NI     = I
               NJ     = J
	       NN     = 0

14             RSCNTR = 0
               RICNTR = 0
               RJCNTR = 0
               WEIGHT = 0
 	       DENS   = 0
               STMP   = 1 + MODULO(NS-1,NZ)
               ITMP   = 1 + MODULO(NI-1,NY)
               JTMP   = 1 + MODULO(NJ-1,NX)
               INCNTR = BUF1(JTMP,ITMP,STMP)

               DO S2=NS-KERNEL, NS+KERNEL
                  STMP = 1 + MODULO(S2-1,NZ)

                  DO I2=NI-KERNEL, NI+KERNEL
                     ITMP = 1 + MODULO(I2-1,NY)

                     DO J2=NJ-KERNEL, NJ+KERNEL
                        JTMP = 1 + MODULO(J2-1,NX)

                        IF (ABS(BUF1(JTMP,ITMP,STMP)-INCNTR) < GRAD)
     &                  THEN
                           RSCNTR = RSCNTR + S2*FILT(J2-NJ,I2-NI,S2-NS)
                           RICNTR = RICNTR + I2*FILT(J2-NJ,I2-NI,S2-NS)
                           RJCNTR = RJCNTR + J2*FILT(J2-NJ,I2-NI,S2-NS)
                           WEIGHT = WEIGHT +    FILT(J2-NJ,I2-NI,S2-NS)
                           DENS   = DENS   + BUF1(JTMP,ITMP,STMP) * 
     &                                       FILT(J2-NJ,I2-NI,S2-NS)
                        ENDIF
                    ENDDO
                ENDDO
             ENDDO
             NN    = NN+1
             SCNTR = INT(RSCNTR/WEIGHT)
             ICNTR = INT(RICNTR/WEIGHT)
             JCNTR = INT(RJCNTR/WEIGHT)
             TEMP  = ABS(JCNTR-NJ) + ABS(ICNTR-NI) + ABS(SCNTR-NS)

             IF (3*TEMP > 2*KERNEL) THEN
                NS = SCNTR
                NJ = JCNTR
                NI = ICNTR
                IF (NN < 40) GOTO 14
             ENDIF

16           BUF2(J,I,S)  =  DENS / REAL(WEIGHT)
          ENDDO
        ENDDO
      ENDDO
       
      END
