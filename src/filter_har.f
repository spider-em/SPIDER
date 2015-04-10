C ++********************************************************************
C                                                                      *
C FILTER_HAR.F            CREATED MAY 01 ARDEAN LEITH                  * 
C                                                                      *
C **********************************************************************
C=*                                                                    *
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
C  FILTER_HAR(LUN1,LUN2,NSAM,NROW,NSLICE,)
C
C  PARAMETERS: LUN1,LUN2   IO UNITS                             (INPUT)
C              NSAM        X DIMENSIONS                         (INPUT)
C              NROW        Y DIMENSIONS                         (INPUT)
C              NSLICE      Z DIMENSIONS                         (INPUT)
C              FMINT       INPUT IMAGE MINIMUM                  (INPUT)
C              FMAXT       INPUT IMAGE MAXIMUM                  (INPUT)

C  PURPOSE: ALTER CONTRAST IN AN IMAGE OR VOLUME USING HARALICK TEXTURE
C           CONVOLUTION METHODS
C                                                                      *
C **********************************************************************

	SUBROUTINE FILTER_HAR(LUN1,LUN2,NSAM,NROW,NSLICE,FMINT,FMAXT)

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

	INTEGER, ALLOCATABLE, DIMENSION(:,:)  :: ICOOC

C       AUTOMATIC ARRAY
        REAL, DIMENSION(NSAM)                 :: VIN,VOUT
 	INTEGER, DIMENSION(NSAM,NROW)         :: IVIN

        IF (NSLICE .GT. 1) THEN
           WRITE(NOUT,*) 'THIS FILTER IS NOT IMPLEMENTED IN 3D!' 
           WRITE(NOUT,*) 'VOLUME WILL BE PROCESSED SLICE-BY-SLICE' 
        ENDIF

C       SET DEFAULT NEIGHBORS
        LX    = 7
        LY    = 7
10      CALL RDPRIS(LX,LY,NOT_USED,'NEIGHBORHOOD X & Y',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (LX .LT. 3 .OR. MOD(LX,2) .EQ. 0 .OR. 
     &      LY .LT. 3 .OR. MOD(LY,2) .EQ. 0) THEN
           CALL ERRT(101,'DIMENSIONS MUST BE ODD  AND > 2',IDUM) 
           GOTO 10
        ENDIF
        LXD2  = LX / 2
        LYD2  = LY / 2
        NEIGH = LX * LY

11      CALL RDPRI1S(INTEN,NOT_USED,'INTENSITY LEVELS',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 10

        IF (INTEN .LT. 3) THEN
           CALL ERRT(101,'NUMBER OF INTENSITY LEVELS MUST BE > 2',IDUM) 
           GOTO 11
        ENDIF

12      IDY = NROW + 1
        CALL RDPRIS(IDX,IDY,NOT_USED,'OFFSET IN X & Y',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 11
        IF (IDY .GT. NROW) IDY = IDX

        IF (IDX .EQ. 0  .AND. IDY .EQ. 0) THEN
           CALL ERRT(101,'AT LEAST ONE OFFSET MUST NON ZERO',IDUM) 
           GOTO 12
        ENDIF

13      CALL RDPRI1S(IMODE,NOT_USED,'MODE NUMBER (1...6)',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 12

        IF (IMODE .LE. 0  .OR. IMODE .GT. 6) THEN
           CALL ERRT(101,'MODE MUST BE 1...6',IDUM) 
           GOTO 13
        ENDIF

C       ALLOCATE SPACE FOR CO-OCURANCE MATRIX & IMAGE
 	ALLOCATE(ICOOC(INTEN,INTEN),STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'CE, ICOOC',IER)
           GOTO 9999
        ENDIF

C       FIND SCALING FACTOR TO CONVERT INTENSITY TO: 1....INTEN
        CALL SCALE32TO8(FMINT,FMAXT,1,INTEN,SCAL,OFFSET,IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           CALL ERRT(101,'BLANK IMAGE',IER)
           GOTO 9999
        ENDIF

        WRITE(NOUT,90) NEIGH
90      FORMAT(' NUMBER OF NEIGHBORS: ',I3)
        FNEIGH  = 1.0 / NEIGH
        FNEIGH2 = 1.0 / (NEIGH**2)

C       PROCESS SLICE BY SLICE (ADDED LOOP)
        DO ISLICE = 1,NSLICE
        IRECGO = (ISLICE-1) * NROW

C       LOAD INPUT SLICE
        DO K = 1,NROW
           CALL REDLIN(LUN1,VIN,NSAM,IRECGO+K)
           DO I = 1,NSAM
              IVIN(I,K) = VIN(I) * SCAL + OFFSET
           ENDDO
        ENDDO

        DO IY=1,NROW
           DO IX = 1,NSAM
C             CREATE CO-OCCURANCY MATRIX
	      CALL FILTER_HAR2(NSAM,NROW,IVIN,LXD2,LYD2,NEIGH,INTEN,
     &                        IX,IY,IDX,IDY,ICOOC)

C             APPLY CO-OCCURANCY MATRIX STATISTICS
              IF (IMODE .EQ. 1) THEN
C                HOMOGENEITY USING SECOND MOMENT
C                CASTLEMAN REFERS TO THIS AS "ENERGY"
                 F1 = 0.0
                 DO I = 1,INTEN
                    DO J = 1,INTEN
                       F1 = F1 + (ICOOC(I,J))**2
                    ENDDO
                 ENDDO
                 VOUT(IX) = F1 * FNEIGH2
 
              ELSEIF (IMODE .EQ.2) THEN
C                CONTRAST USING DIFFERENCE MOMENT
                 F2 = 0.0
                 DO N = 1,INTEN-1
                    F2T = 0.0
                    DO I = 1,INTEN-N
                       J   = I + N
                       F2T = F2T + ICOOC(I,J) + ICOOC(J,I)
                    ENDDO
                    F2 = F2 + F2T * N**2
                 ENDDO
                 VOUT(IX) = F2 * FNEIGH
 
              ELSEIF (IMODE .EQ.3) THEN
C                WEIGHTED AVERAGE ABSOLUTE DISTANCE FROM DIAGONAL
C                SAID TO MEASURE "LACK OF SMOOTHNESS" IN GOSE,
C                JOHNSONBAUGH & JOST, "PATTERN RECOGNITION AND IMAGE
C                ANALYSIS"
C                CASTLEMAN REFERS TO THIS AS "INERTIA"

                 D1 = 0.0
                 D2 = 0.0 
                 DO I = 1,INTEN
                    DO J = 1,INTEN
                       D1 = D1 + ABS(I-J) * ICOOC(I,J)
                       D2 = D2 + ICOOC(I,J)
                    ENDDO
                 ENDDO
                 VOUT(IX) =  D1/D2

              ELSEIF (IMODE .EQ.4) THEN
C                CASTLEMAN REFERS TO THIS AS "ENTROPY"

                 E1 = 0.0
                 DO I = 1,INTEN
                    DO J = 1,INTEN
                       IF (ICOOC(I,J) .NE. 0) THEN
                          E1 = E1 + LOG(FLOAT(ICOOC(I,J))) * ICOOC(I,J)
                       ENDIF
                    ENDDO
                 ENDDO
                 VOUT(IX) =  E1

              ELSEIF (IMODE .EQ.5) THEN
C                INTENSITY OF MAXIMUM COCCURANCE PROBABILITY

                 MAX1 = 0
                 DO I = 1,INTEN
                    DO J = 1,INTEN
                       IF (ICOOC(I,J) .GT. MAX1) THEN
                          MAX1 = ICOOC(I,J)
                          MAXJ = J
                       ENDIF
                    ENDDO
                 ENDDO
                 VOUT(IX) =  MAXJ

              ELSEIF (IMODE .EQ.6) THEN
C                RUSS REFERS TO THIS AS "LINEAR DEPENDENCE OF BRIGHTNESS"
C                BUT I HAVE NEGLECTED CONSTANT TERMS FROM SUMMATION
                 F3 = 0.0
                 DO I = 1,INTEN
                    DO J = 1,INTEN
                       F3 = F3 + I * J * ICOOC(I,J)
                    ENDDO
                 ENDDO
                 VOUT(IX) =  F3
              ENDIF
           ENDDO
C          END LOOP:  DO IX=1,NSAM

C          OUTPUT IMAGE
           CALL WRTLIN(LUN2,VOUT,NSAM,IRECGO+IY)

        ENDDO
        ENDDO

9999    IF (ALLOCATED(ICOOC)) DEALLOCATE(ICOOC)

        END



C       ------------------------- FILTER_HAR2----------------------------

	SUBROUTINE FILTER_HAR2(NSAM,NROW,IVIN,LXD2,LYD2,NEIGH,INTEN,
     &                        IX,IY,IDX,IDY,ICOOC)

C       SURE, THIS IS SLOW, BUT IT IS NOT IMPORTANT!

	INTEGER, DIMENSION(NSAM,NROW) ::   IVIN
	INTEGER, DIMENSION(INTEN,INTEN) :: ICOOC

C       COMPUTE CO-OCCURANCE MATRIX FOR INPUT IMAGE (IVIN) HAVING
C       A SET NUMBER OF INTENSITES (INTEN).  CO-OCCURANCE IS FOR
C       PIXELS WITH OFFSET (IDX,IDY)  WITH NEIGHBORHOOD 
C       (-LXD2...LXD2,-LYD2...LYD2) AROUND A CENTRAL PIXEL (IX,IY)

C       ZERO THE WHOLE CO-OCCURANCE MATRIX 
        ICOOC = 0

        DO MY=-LYD2,LYD2
           IYT  = MOD(IY+MY+NROW-1,NROW)+1
           IYTD = MOD(IY+IDY+MY+NROW-1,NROW)+1
           DO MX=-LXD2,LXD2
              IXT  = MOD(IX+MX+NSAM-1,NSAM)+1
              IXTD = MOD(IX+IDX+MX+NSAM-1,NSAM)+1

              IV1  = IVIN(IXT,IYT)
              IV1D = IVIN(IXTD,IYTD)
              ICOOC(IV1,IV1D) = ICOOC(IV1,IV1D) + 1
           ENDDO
        ENDDO
          
        END


C       ------------------------- SCALE32TO8 ---------------------------

      SUBROUTINE SCALE32TO8(FMIN,FMAX,NMINT,NMAXT,SCAL,OFFSET,IRTFLG)

C     USAGE:  IVAL = FVAL * SCAL + OFFSET

      EPS = EPSILON(EPS)

      IF (ABS(FMAX - FMIN) < EPS) THEN
         IRTFLG = 1
         RETURN
      ENDIF

C     CONVERSION FACTOR FROM FLOATING POINT TO INTEGER RANGE 
      SCAL   = FLOAT(NMAXT - NMINT) / (FMAX - FMIN)
      OFFSET = -FMIN * SCAL + NMINT + 0.5

      IRTFLG = 0
      END
