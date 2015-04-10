
C++*********************************************************************
C
C  ENHANC.F    FIXED UNDEFINED BOTTOM BUG            APR 02 ARDEAN LEITH
C              REMOVED 'CE L" (BUGGY FOR 16 YRS)            ARDEAN LEITH
C              SETPRMB PARAMETERS                    MAY 09 ARDEAN LEITH
C              NBINS, AUTO ? BUG, FORMAT             MAR 11 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2011  Health Research Inc.,                         *
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
C    ENHANC(FILNAM,LUNI,LUNO,NSAM,NROW)
C
C    PURPOSE: IMAGE ENHANCEMENT ROUTINE.  CALLS OTHER ROUTINES FOR
C             HISTOGRAM BASED ENHANCEMENTS AND DOES THRESHOLDING.
C
C    PARAMETERS:
C         FILNAM     NAME OF FILE
C         LUNI       LOGICAL UNIT NUMBER OF INPUT FILE
C         LUNO       LOGICAL UNIT NUMBER OF OUTPUT FILE
C         NSAM,NROW  DIMENSIONS OF FILE
C
C--*******************************************************************

      SUBROUTINE ENHANC(FILNAM,LUNI,LUNO,NSAM,NROW,NSLICE)

      COMMON BUF(1)

      INCLUDE 'CMBLOCK.INC'

      REAL             :: VAL(4)
      CHARACTER(LEN=*) :: FILNAM
      CHARACTER(LEN=1) :: ANS
      CHARACTER(LEN=1) :: NULL 
      CHARACTER(LEN=1) :: IRESH = 'N'

      EQUIVALENCE(B1,VAL(1)),(T1,VAL(2)),(B2,VAL(3)),(T2,VAL(4))

      INTEGER          :: NBINS = 128

      NULL = CHAR(0)

      MAPA = NSAM + 128

      CALL RDPRMC(ANS,NC,.TRUE.,
     &  '(S)INGLE, (A)UTOMATIC OR (D)OUBLE MAPPING? (S/A/D)',NULL,IRT)
      IF (IRT .EQ. -1) RETURN

      IF (ANS .EQ. 'A') THEN
         CALL RDPRM(PERC,NOT_USED,'INTEGRAL THRESHOLD PERCENT')

         CALL RDPRMC(IRESH,NC,.TRUE.,'PLOT RESULT HISTOGRAM? (Y/N)',
     &               NULL,IRT)
      ENDIF

      IF (ANS .EQ. 'S') THEN
         BOTTOM = 0        ! ARBITRARY USELESS
         CALL RDPRM1S(BOTTOM,NOT_USED,'BOTTOM DENSITY VALUE',IRTFLG)

         TOP    = 255      ! ARBITRARY USELESS
         CALL RDPRM1S(TOP,   NOT_USED,'TOP DENSITY VALUE',IRTFLG)

         IF (BOTTOM .GE. TOP) GOTO 650
      ENDIF

C     DETERMINE HISTOGRAM FROM IMAGE 
      CALL HIST(LUNI,0,0,NSAM,NROW,NSLICE,HMIN,HMAX,HSIG,HMODE)

      HINCO = (HMAX - HMIN) / FLOAT(NBINS - 1)

      WRITE(NOUT,1111) HMIN,HMAX
1111  FORMAT('  HISTOGRAM MINIMUM =',1PG13.5,
     &       '  HISTOGRAM MAXIMUM =',1PG13.5)

      IF (ANS .EQ. 'D') THEN  ! ----------------------------------- D
C        NON-UNIQUE MAPPING WITH DOUBLE LINEAR MAPPING FUNCTION
100      CALL RDPRM(B1,NOT_USED,'BOTTOM1')
         CALL RDPRM(T1,NOT_USED,'TOP1   ')
         CALL RDPRM(B2,NOT_USED,'BOTTOM2')
         CALL RDPRM(T2,NOT_USED,'TOP2   ')

         IF (B1 .LT. FMIN) B1 = FMIN
         IF (T2 .GT. FMAX) T2 = FMAX

C        CHECK IF B1,T1,B2,T2 ARE ORDERED ACCORDING TO INCREASING VALUE
         A = VAL(1)
         DO I = 1,4
           A = MAX(A,VAL(I))
           IF (A .NE. VAL(I)) GOTO 650
	 ENDDO

         IF (B1 .NE. T1) THEN
            IF (B2 .NE. T2) GOTO 130
            BOTTOM = B1
            TOP    = T1
         ELSE
            BOTTOM = B2
            TOP    = T2
         ENDIF
         GOTO 30


130      NB1 = (B1-HMIN) / HINCO + 1.5
         NT1 = (T1-HMIN) / HINCO + 1.5

C        NEW DENSITY INCREMENT ON LEFT HAND SIDE
         HINCN1 = 2.0 / FLOAT(NT1-NB1)
         NB2    = (B2-FMIN) / HINCO+1.5
         NT2    = (T2-FMIN) / HINCO+1.5
C        NEW DENSITY INCREMENT ON RIGHT HAND SIDE
         HINCN2 = 2.0 / FLOAT(NT2-NB2)

         IF (NB1 .GE. 2) THEN
           DO  K = 1,NB1-1
              BUF(MAPA+K) = 0.
	   ENDDO
         ENDIF
         DO  K = NB1,NT1
            BUF(MAPA+K) = FLOAT(K-NB1)* HINCN1
	 ENDDO
         IF (NT1 .NE. NB2) THEN
             DO  K = NT1+1,NB2-1
                BUF(MAPA+K) = 2.
	     ENDDO
         ENDIF
         DO  K = NB2,NT2
            BUF(MAPA+K) = FLOAT(K-NB2)*HINCN2      
	 ENDDO
         IF (NT2 .NE. NBINS) THEN
            DO K = NT2+1,NBINS
               BUF(MAPA+K) = 2.
	    ENDDO
         ENDIF
         GOTO 500

      ELSEIF (ANS .EQ. 'A') THEN ! ----------------------------- A
        P     = PERC*FLOAT(NSAM)*FLOAT(NROW)*FLOAT(NSLICE)/100.
        ADD   = 0.0

        DO   I = 1,NBINS
           IM  = I
           ADD = ADD + BUF(NSAM+I)
           WRITE(NDAT,2222)I,P,ADD,BUF(NSAM+I)
2222       FORMAT(1X,I3,3F10.5)
           IF (ADD .GT. P) EXIT
	ENDDO

        BOTTOM = FLOAT(IM) * HINCO + HMIN
        ADD    = 0.
        DO  I = NBINS,IM,-1
           IC  = I
           ADD = ADD+BUF(NSAM+I)
           IF (ADD .GT. P) EXIT
	ENDDO

        TOP = FLOAT(IC) * HINCO + HMIN
      ENDIF


      ! ------------------------------------------------------- S & A

30    IF (BOTTOM .LT. FMIN) BOTTOM = FMIN
      IF (TOP    .GT. FMAX) TOP    = FMAX
      IF (BOTTOM .EQ. FMIN .AND. TOP .EQ. FMAX) RETURN

      NB = (BOTTOM-HMIN) / HINCO+1.5
      NT = (TOP   -HMIN) / HINCO+1.5
      WRITE(NOUT,31)BOTTOM,TOP
31    FORMAT('  BOTTOM DENSITY ',1PG10.2,' , TOP DENSITY ',1PG10.2)
      HINCN = 2./FLOAT(NT-NB)

      IF (NB .GE. 2) THEN
         DO  K = 1,NB-1
           BUF(MAPA+K) = 0.
	 ENDDO
      ENDIF

70    DO  K = NB,NT
         BUF(MAPA+K) = FLOAT(K-NB) * HINCN
      ENDDO
      IF (NT .NE. NBINS) THEN
          DO  K = NT+1,NBINS
             BUF(MAPA+K) = 2.
	  ENDDO
      ENDIF

C     APPLY MAPPING FUNCTION TO DATA.
C     RESULT IS NORMALIZED BETWEEN 0. AND 2.

500   CALL GRAPHS(NDAT,BUF(MAPA+1),NBINS,1,0,1.0,IRTFLG)

      AV = 0.0
      DO  I = 1,NROW*NSLICE
         CALL REDLIN(LUNI,BUF,NSAM,I)
         DO  K     = 1,NSAM
            MAP    = (BUF(K)-HMIN)/HINCO+1.5
            T      = BUF(MAPA+MAP)
            AV     = AV + T
            BUF(K) = T
	 ENDDO
         CALL WRTLIN(LUNO,BUF,NSAM,I)
      ENDDO

      FMAX  = 0.0
      FMIN  = 0.0
      CALL SETPRMB(LUNO, 0.0,0.0, 0.0,0.0)

      IF (IRESH .EQ. 'Y') 
     &    CALL HIST(LUNO,0,0,NSAM,NROW,NSLICE,HMIN,HMAX,HSIG,HMODE)

      RETURN


650   CALL ERRT(14,'ENHANC',NE)

      END
