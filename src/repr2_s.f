
C ++********************************************************************
C                                                                      *
C  REPR2_S                                                             *
C           REFORMATTED, LESS OUTPUT          OCT  8 2009 ArDean Leith *
C           ENDLESS LOOP WITH BAD DATA FIXED  FEB  8 2014 ArDean Leith *
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
C                                                                      *
C  REPR2_S(BCKE,BCN,NSAM,LCYL,NSLICE,                                  *
C          NROWL,NROWH,NANG,IPCUBE,NN,PROJ,PRN,IRI,LTB,LTBN,ABA,INPIC) *
C                                                                      *
C  PURPOSE:  ITERATIVE BACK PROJECTION                                 *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE REPR2_S(BCKE,BCN,NSAM,LCYL,NSLICE,
     &	   NROWL,NROWH,NANG,IPCUBE,NN,PROJ,PRN,IRI,LTB,LTBN,ABA,INPIC)

        INCLUDE 'CMBLOCK.INC' 

	REAL             :: BCKE(NSAM*NSLICE*3),BCN(NSAM*NSLICE)
	REAL             :: PROJ(NSAM,NANG),PRN(NSAM,NANG)
	INTEGER          :: IPCUBE(5,NN)
	DOUBLE PRECISION :: ABA
	LOGICAL          :: ACTIVE_MIN,ACTIVE_MAX,ENDIT,DODO,DOD

	DATA IOFF/6/    ! IO UNIT SET IN bps2.f

	CALL RDPRM2(ALA,AIM,NOT_USED,'LAMBDA, CORRECTION LIMIT')

	CALL RDPRMI(MAXIT,INTER,NOT_USED,
     &	            'ITERATION LIMIT, NUMBER OF INTERNAL STEPS')

	CALL RDPRMI(MODE,IDUM,NOT_USED,'MODE')

	CALL RDPRM2(TMIN,TMAX,NOT_USED,'MINIMUM, MAXIMUM')

	TMIN = TMIN - ABA
	TMAX = TMAX - ABA

	WRITE(NOUT,2059) TMIN,TMAX
2059	FORMAT(/,'  MIN. & MAX. AFTER AVERAGE SUBTRACTION:',
     &          2(5X,1PE10.3))

	CALL RDPRM(T,NOT_USED,'SMOOTHING CONSTANT')

	T = MIN(MAX(T,0.0), 1.0)

	IF (.NOT.
     &     (MODE == 1 .OR. MODE == 3 .OR.
     &      MODE == 6 .OR. MODE == 8)) THEN
	   IF (MAXIT > 1)  THEN
	      INTER  = INTER * MAXIT
	      MAXIT  = 1
	      DODO   = .FALSE.
	      WRITE(NOUT,6066) INTER
6066	      FORMAT(/ ' NO SMOOTHING REQUESTED, ITERATION LIMIT ',
     &           'SET TO 1 AND NO. OF INTERNAL STEPS IS SET TO: ',I0,/)
	   ELSE
	      DODO = .TRUE.
	   ENDIF
	ENDIF

	NMAT       = NSAM * NSLICE
        ACTIVE_MIN = .FALSE.
        ACTIVE_MAX = .FALSE.
        SQOLD      = 1.E23
        ENDIT      = .FALSE.
        ITER       = 0

	DO WHILE (.NOT. ENDIT)
            ITER = ITER + 1
            WRITE(NOUT,971) ITER
971         FORMAT('  ITERATION: ',I6)

            SQ = 0.0
            DO LS = NROWL,NROWH
               IF (ITER > 1)  THEN
                  DO IRC=1,NSLICE
                     IRB = 1+(IRC-1)*NSAM
                     CALL REDLIN(INPIC,BCKE(IRB),NSAM,
     &                            (IRC-1)*LCYL+LS-NROWL+1)
                  ENDDO
               ELSE
                  DO I=1,NMAT
                     BCKE(I) = 0.0
                  ENDDO
               ENDIF

C              GET THE PROJECTION DATA, READS FROM MULTIPLE OPEN FILES
               DO J=1,NANG
                  CALL REDLIN(IOFF+J,PROJ(1,J),NSAM,LS)
               ENDDO

C              REMOVE AVERAGE
               DO J=1,NANG
                  DO I=1,NSAM
                     PROJ(I,J) = PROJ(I,J) - ABA
                  ENDDO
               ENDDO

               SQINT = 1.0E23
               DO ITR=1,INTER
                  IF (.NOT. ENDIT)  THEN
                     IF (ITER == 1 .AND. ITR == 1)  THEN
c$omp                   parallel do private(i,j)
                        DO J=1,NANG
                           DO I=1,NSAM
                              PRN(I,J) = PROJ(I,J)
                           ENDDO
                        ENDDO
                     ELSE
C                       BCKE -> PRN
                        CALL  PRJS2
                        IF (DODO) THEN
                           DOD = (ITR  ==  INTER)
                        ELSE
                           DOD = .TRUE.
                        ENDIF

                        IF (DOD) THEN
			   CALL FMIN_S(PRN,LTBN,GMIN)
			   CALL FMAX_S(PRN,LTBN,GMAX)

cal09                      WRITE(NOUT,2062) LS,GMIN,GMAX
2062	                   FORMAT('  MIN &  MAX IN SLICE: ',I0,'= ',
     &                            1PE10.3,3X,1PE10.3)
                        ENDIF
                        IF ((MODE == 2.OR.MODE == 3 .OR.
     &                       MODE == 7.OR.MODE == 8)
     &		           .AND. .NOT. ACTIVE_MIN .AND. DOD)  THEN
                           IF (GMIN < TMIN)  THEN
                              CALL BMIN_S(BCKE,NMAT,IPCUBE,NN,BMIN)

cal09                         WRITE(NOUT,2051) BMIN
2051	                      FORMAT(
     &                       '  Min constraint active, value in 3D:',
     &                        1PE10.3)

			      ACTIVE_MIN = .TRUE.
			   ENDIF
                        ENDIF

	                IF ((MODE == 5.OR.MODE == 6.OR.MODE == 7 .OR.
     &                       MODE == 8) 
     &		           .AND. .NOT.ACTIVE_MAX .AND. DOD)  THEN
                           IF (GMAX > TMAX)  THEN
                              CALL BMAX_S(BCKE,NMAT,IPCUBE,NN,BMAX)
cal09                         WRITE(NOUT,2052) BMAX
2052                          FORMAT(
     &                        '  Max constraint active, value in 3D:',
     &                        1PE10.3)
			      ACTIVE_MAX = .TRUE.
			   ENDIF
                        ENDIF
                         
c$omp                   parallel do private(i,j)
                        DO J=1,NANG
                           DO I=1,NSAM
                             PRN(I,J) = PROJ(I,J) - PRN(I,J)
                           ENDDO
		        ENDDO
                     ENDIF
                     SQT = 0.0

c$omp              parallel do private(i,j),reduction(+:sqt)
                   DO J=1,NANG
		      DO I=1,NSAM
		         SQT = SQT + PRN(I,J) * PRN(I,J)
                      ENDDO
		   ENDDO

                   IF (SQT <= SQINT)  THEN
                      SQINT = SQT
		      IF (DOD) SQ = SQ + SQT

C                     PRN  -> BCN
                      CALL BPRJ2
		      CALL DOCORR3_S(BCKE,BCN,NMAT,IPCUBE,NN,ALA,SBQ)
cal09 	              WRITE(NOUT,*)' Squared correction of structure: ',
cal09     &                     LS,ITR,SBQ
                      IF (MODE .NE. 0)  THEN
		         IF (ACTIVE_MIN) 
     &                     CALL DOMIN3_S(BCKE,NMAT,IPCUBE,NN,BMIN)
                         IF (ACTIVE_MAX) 
     &                     CALL DOMAX3_S(BCKE,NMAT,IPCUBE,NN,BMAX)
                      ENDIF
                   ELSE
C                     BREAK IT UP, THE ERROR INCREASED.
                      ENDIT = .TRUE.
                   ENDIF
                ENDIF
             ENDDO    ! END OF DO-LOOP OVER ITR

             DO IRC=1,NSLICE
                IRB = 1 + (IRC-1)*NSAM
                CALL WRTLIN(INPIC,BCKE(IRB),NSAM,
     &                        (IRC-1)*LCYL+LS-NROWL+1)
             ENDDO
          ENDDO

cal09	  WRITE(NOUT,*)' Squared discrepancies between projections: ',SQ
          IF (MODE == 1 .OR. MODE == 3 .OR. 
     &        MODE == 6 .OR. MODE == 8) THEN

             DO LS=NROWL+1,NROWH-1
                IF (LS  ==  NROWL+1)  THEN
                    DO IRC=1,NSLICE
                       IRB = 1+ (IRC-1)*NSAM
	               CALL REDLIN(INPIC,BCKE(IRB),NSAM,
     &                           (IRC-1)*LCYL+LS-1-NROWL+1)
                    ENDDO

                    DO IRC=1,NSLICE
                       IRB = 1+ (IRC-1) *NSAM+NSAM*NSLICE
                       CALL REDLIN(INPIC,BCKE(IRB),NSAM,
     &                           (IRC-1)*LCYL+LS-NROWL+1)
                    ENDDO
                 ENDIF
		 DO IRC=1,NSLICE
		    IRB = 1 + (IRC-1)*NSAM + 2*NSAM*NSLICE
	            CALL REDLIN(INPIC,BCKE(IRB),NSAM,
     &                          (IRC-1)*LCYL+LS+1-NROWL+1)
		 ENDDO

                 CALL SMT3(T,BCKE,BCN,NSAM,NSLICE,IPCUBE,NN)

                 DO IRC=1,NSLICE
		    IRB = 1 + (IRC-1)*NSAM
	            CALL WRTLIN(INPIC,BCN(IRB),NSAM,
     &                         (IRC-1)*LCYL+LS-NROWL+1)
                 ENDDO
	         CALL CPPB(BCKE(NSAM*NSLICE+1),BCKE,2*NSAM*NSLICE)
	      ENDDO
	   ENDIF

           IF (.NOT. ENDIT) THEN
              IF (SQ > AIM .AND. ITER < MAXIT)  THEN
                 IF (ITER  ==  1)  THEN
C                   CONTINUE ITERATIONS
                 ELSEIF (SQ < SQOLD)  THEN
                    SQOLD = SQ      ! CONTINUE ITERATIONS
	         ELSE
	            ENDIT = .TRUE.  ! STOP ITERATIONS
	         ENDIF
              ELSE
	         ENDIT = .TRUE.     ! STOP ITERATIONS
	      ENDIF
	   ENDIF
	ENDDO

        END         
