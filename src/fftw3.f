C++*********************************************************************
C
C FFTW3.F                 CREATED FROM FMRS3.F    DEC 2007 ARDEAN LEITH
C                         NUMTHWANTR              JAN 2008 ARDEAN LEITH
C                         INSIDE OMP THREAD OK    JAN 2008 ARDEAN LEITH
C                         REWRITE                 JAN 2008 ARDEAN LEITH
C                         MAKEPLANB               MAR 2008 ARDEAN LEITH
C                         ECHO FORMATTING         DEC 2014 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C  FFTW3_MAKEPLANB(BUF,LDA,NX,NY,NZ,NUMTHINOUT, PLAN,INV,IRTFLG)
C
C  PURPOSE:   USES PRE_EXISTING BUFFER, CREATES PLAN
C
C
C  FFTW3_MAKEPLAN (NX,NY,NZ,NUMTHINOUT, PLAN,INV,IRTFLG)
C
C  PURPOSE:   ALLOCATES AND FREES BUFFER, CREATES PLAN
C
C
C  FFTW3_USEPLAN(BUF,NX,NY,NZ, SPIDER_SIGN, SPIDER_SCALE,PLAN,INV,IRTFLG)
C
C  PURPOSE:   USES A FFTW3 PLAN TO DO A FFT ON BUF
C
C
C  FFTW3_KILLPLAN(PLAN,IRTFLG)
C
C  PURPOSE:   DESTROYS A FFTW3 PLAN.
C
C
C  PARAMETERS:  BUF     ARRAY (LDA*NY*NZ)                SENT/RET.
C               NX..    DIMENSONS                             SENT
C               INV     1=REG. FILE, -1= FOURIER FILE         SENT
C               IRTFLG  ERROR FLAG (O=NORMAL)                 RET.
C   
C  NOTES:    REAL MIXED RADIX FFT.
C            INPUT:  A(N) - REAL ARRAY
C            OUTPUT: N EVEN  A(N+2)
C            ORDER OF ELEMENTS:
C            R(0),0.0, R(1), I(1), R(2), I(2), ....., 
C                  R(N/2-1), I(N/2-1), R(N/2),0.0
C
C            N ODD  X(N+1)
C                 R(0),0.0, R(1), I(1), R(2), I(2), ....., 
C                 R(N/2-1), I(N/2-1), R(N/2),I(N/2)
C            (INTEGER DIVISION ROUNDS DOWN, E.G. 5/2 =2)
C
C--*********************************************************************


       MODULE TYPE_KINDS

       INTEGER, PARAMETER, PUBLIC  :: I_4  = SELECTED_INT_KIND(8)  !I*4
       INTEGER, PARAMETER, PRIVATE :: I_8T = SELECTED_INT_KIND(16) !I*8
       INTEGER, PARAMETER, PUBLIC  :: I_8  = (((1 + SIGN(1,I_8T)) / 2)*
     &                     I_8T) + (((1 - SIGN(1,I_8T))/2) * I_4 )! I*8

C       USAGE:  USE TYPE_KINDS
C       USAGE:  INTEGER(KIND=I_8) :: INTPOINTER

        END MODULE TYPE_KINDS


C       ---------------  FFTW3_MAKEPLAN -----------------------------

	SUBROUTINE FFTW3_MAKEPLAN(NX,NY,NZ,NUMTHINOUT,
     &                            PLAN,INV,IRTFLG)

        USE TYPE_KINDS

        INCLUDE 'CMBLOCK.INC'

        INTEGER, INTENT(INOUT) :: NUMTHINOUT
        INTEGER, INTENT(IN)    :: INV
        INTEGER, INTENT(OUT)   :: IRTFLG
	
C       PLAN IS A POINTER TO A STRUCTURE 
        INTEGER(KIND=I_8), INTENT(INOUT) :: PLAN

        INTEGER           :: OMP_GET_NUM_THREADS
        LOGICAL, SAVE     :: MUST_INIT = .TRUE.

	REAL, ALLOCATABLE, DIMENSION(:,:,:) :: AF

C       USING FFTW3 LIBRARY CALLS FOR FFT

#include "FFTW3.INC"

        IRTFLG = 1

#ifdef SP_MP
C       COMPILED FOR OMP USE

        IF (NUMTHINOUT .EQ. 0) THEN
           NUMTHWANT = NUMFFTWTH
        ELSE 
C          INPUT HAS SPECIFIED NUMTHINOUT     
           NUMTHWANT = NUMTHINOUT
        ENDIF

c       write(6,*) ' numthwant:',numthwant

        NUMTH = OMP_GET_NUM_THREADS()

        IF (NUMTH .GT. 1) THEN
C          CURRENTLY INSIDE OMP PARALLEL SECTION --> FFTW3 ERRORS

           WRITE(NOUT,98)NX,NY,NZ,NUMTHWANT,INV
 98        FORMAT('  Trying to create FFTW3 plan for: (',
     &               I0,',', I0,',', I0,') Threads: ',I0,' INV: ',I0)              

           CALL ERRT(102,
     &       'PGM ERROR, CAN NOT CALL FFTW3_MAKEPLAN WITHIN OMP ||:',
     &       NUMTH)
           RETURN
        ENDIF

        IF (NUMTHWANT .GT. 1 .AND. MUST_INIT) THEN        
C          MUST INITIALIZE THREADS ONCE

           WRITE(NOUT,*) ' Initializing FFTW3 for threaded use' 

           IRTFLGL = -4
           CALL SFFTW_INIT_THREADS(IRTFLGL)

	   IF (IRTFLGL .EQ. 0) THEN
              CALL ERRT(101,'FFTW3 MULTIPLE THREADS INIT. FAILED',IER)
              RETURN
           ENDIF
           MUST_INIT = .FALSE.
        ENDIF

#else
        IF (NUMFFTWTH .EQ. 0) WRITE(NOUT,*) ' Single threaded FFTW3 ' 
        NUMFFTWTH = 1
        NUMTHWANT = 1
#endif

        LDA = NX + 2 - MOD(NX,2)
        JH  = LDA/2

C       ALLOCATE TEST ARRAY: AF
        ALLOCATE(AF(LDA,NY,NZ), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'FMRS_FFTW3, AF',LDA*NY*NZ)
           RETURN
        ENDIF

        IF (INV .GT. 0) THEN
C          CREATE PLAN FOR FORWARD TRANSFORM

          !write(6,*) 'fftw3Creating -> plan:',NX,NY,NUMTHWANT,lda

#ifdef SP_MP
          !write(6,*) 'sfftw_plan_with_nthreads:',numthwant
          CALL SFFTW_PLAN_WITH_NTHREADS(NUMTHWANT)
#endif

C         DESTROY EXISTING OLD PLAN
          !write(6,*) 'calling : sfftw_destroy_plan: ',PLAN
          IF (PLAN .GT. 0) CALL SFFTW_DESTROY_PLAN(PLAN)

          IF (NZ .GT. 1) THEN
             CALL SFFTW_PLAN_DFT_R2C_3D(PLAN,NX,NY,NZ,
     &                                  AF,AF,FFTW_MEASURE)
          ELSEIF (NY .GT. 1) THEN
             !write(6,*) 'calling : SFFTW_PLAN_DFT_R2C_2D '
             CALL SFFTW_PLAN_DFT_R2C_2D(PLAN,NX,NY,
     &                                  AF,AF,FFTW_MEASURE)
          ELSE
C            BIMAL SAID USE ESTIMATE, I DO NOT KNOW WHY?
	     CALL SFFTW_PLAN_DFT_R2C_1D(PLAN,NX,
     &                                  AF,AF,FFTW_ESTIMATE)
          ENDIF

          !write(6,90)PLAN,NX,NY,NZ,NUMTHWANT
          IF (VERBOSE) WRITE(NOUT,90)PLAN,NX,NY,NZ,NUMTHWANT
 90       FORMAT('  New Forward FFTW3 Plan: ',I0,' ( ',
     &              I0,', ', I0,', ', I0,') Threads: ',I0)              

       ELSE
C         REVERSE TRANSFORM

          IF (PLAN .GT. 0) CALL SFFTW_DESTROY_PLAN(PLAN)
#ifdef SP_MP
          CALL SFFTW_PLAN_WITH_NTHREADS(NUMTHWANT)
#endif
                  
          IF (NZ .GT. 1) THEN
             CALL SFFTW_PLAN_DFT_C2R_3D(PLAN,NX,NY,NZ,
     &                                  AF,AF,FFTW_MEASURE)
          ELSEIF (NY .GT. 1) THEN
             CALL SFFTW_PLAN_DFT_C2R_2D(PLAN,NX,NY,
     &                                   AF,AF,FFTW_MEASURE)
          ELSE
	     CALL SFFTW_PLAN_DFT_C2R_1D(PLAN,NX, AF,AF,FFTW_MEASURE)
          ENDIF

          !write(6,91)PLAN,NX,NY,NZ,NUMTHWANT
          IF (VERBOSE) WRITE(NOUT,91)PLAN,NX,NY,NZ,NUMTHWANT
 91       FORMAT('  New Reverse FFTW3 Plan: ',I0,' ( ',
     &              I0,', ', I0,', ', I0,') Threads: ',I0)              
        ENDIF

C       FREE MEMORY	   
        IF (ALLOCATED(AF)) DEALLOCATE(AF)

        NUMTHINOUT = NUMTHWANT
        IRTFLG  = 0

	END


C       ---------------  FFTW3_MAKEPLANB -----------------------------

	SUBROUTINE FFTW3_MAKEPLANB(AF,LDA,NX,NY,NZ,
     &                             NUMTHINOUT, PLAN,INV,IRTFLG)

        USE TYPE_KINDS

	REAL,    INTENT(INOUT)           :: AF(LDA,NY,NZ)
        INTEGER, INTENT(INOUT)           :: NUMTHINOUT
        INTEGER, INTENT(IN)              :: INV
        INTEGER, INTENT(OUT)             :: IRTFLG
	
C       PLAN IS A POINTER TO A STRUCTURE 
        INTEGER(KIND=I_8), INTENT(INOUT) :: PLAN

        INTEGER                          :: OMP_GET_NUM_THREADS
        LOGICAL, SAVE                    :: MUST_INIT = .TRUE.

        INCLUDE 'CMBLOCK.INC'


C       USING FFTW3 LIBRARY CALLS FOR FFT

#include "FFTW3.INC"

        IRTFLG = 1

#ifdef SP_MP
C       COMPILED FOR OMP USE

        IF (NUMTHINOUT .EQ. 0) THEN
           NUMTHWANT = NUMFFTWTH
        ELSE 
C          INPUT HAS SPECIFIED NUMTHINOUT     
           NUMTHWANT = NUMTHINOUT
        ENDIF

c       !write(6,*) ' NUMTHWANT:',NUMTHWANT
        NUMTH = OMP_GET_NUM_THREADS()

        !write(6,*) ' Current omp threads in fftw3 call:',NUMTH
        IF (NUMTH .GT. 1) THEN
C          CURRENTLY INSIDE OMP PARALLEL SECTION --> FFTW3 ERRORS

           WRITE(NOUT,98)NX,NY,NZ,NUMTHINOUT,INV
 98        FORMAT('  Trying to create FFTW3 plan for: ( ',
     &               I0,', ',I0,', ',I0,') Threads: ',I0,' INV: ',I0)              

           CALL ERRT(101,
     &       'PGM ERROR, CAN NOT CALL FFTW3_MAKEPLAN WITHIN OMP ||:',
     &        NUMTH)
           RETURN
        ENDIF

        IF (NUMTHWANT .GT. 1 .AND. MUST_INIT) THEN        
C          MUST INITIALIZE THREADS ONCE

           WRITE(NOUT,*) ' Initializing FFTW3 for threaded use' 

           IRTFLGL = -4
           CALL SFFTW_INIT_THREADS(IRTFLGL)

	   IF (IRTFLGL .EQ. 0) THEN
              CALL ERRT(101,
     &           'FFTW3 MULTIPLE THREADS INITIALIZATION FAILED',IER)
              RETURN
           ENDIF
           MUST_INIT = .FALSE.
        ENDIF

#else
        IF (NUMFFTWTH .EQ. 0) WRITE(NOUT,*) ' Single threaded FFTW3 ' 
        NUMFFTWTH = 1
        NUMTHWANT = 1
#endif

        LDA = NX + 2 - MOD(NX,2)
        JH  = LDA/2

        IF (INV .GT. 0) THEN
C          CREATE PLAN FOR FORWARD TRANSFORM
           !write(6,*) 'Creating forward plan:',NX,NY,NUMTHWANT

#ifdef SP_MP
          CALL SFFTW_PLAN_WITH_NTHREADS(NUMTHWANT)
#endif

C         DESTROY EXISTING OLD PLAN
          IF (PLAN .GT. 0) CALL SFFTW_DESTROY_PLAN(PLAN)

          IF (NZ .GT. 1) THEN
             CALL SFFTW_PLAN_DFT_R2C_3D(PLAN,NX,NY,NZ,
     &                                  AF,AF,FFTW_MEASURE)
          ELSEIF (NY .GT. 1) THEN
             CALL SFFTW_PLAN_DFT_R2C_2D(PLAN,NX,NY,
     &                                  AF,AF,FFTW_MEASURE)
          ELSE
C            BIMAL SAID USE ESTIMATE, I DO NOT KNOW WHY?
	     CALL SFFTW_PLAN_DFT_R2C_1D(PLAN,NX,
     &                                  AF,AF,FFTW_ESTIMATE)
          ENDIF

          !write(6,90)PLAN,NX,NY,NZ,NUMTHWANT
          IF (VERBOSE) WRITE(NOUT,90)PLAN,NX,NY,NZ,NUMTHWANT
 90       FORMAT('  New Forward FFTW3 Plan: ',I0,' ( ',
     &              I0,', ', I0,', ', I0,') Threads: ',I0)              

       ELSE
C         REVERSE TRANSFORM

          IF (PLAN .GT. 0) CALL SFFTW_DESTROY_PLAN(PLAN)
#ifdef SP_MP
          CALL SFFTW_PLAN_WITH_NTHREADS(NUMTHWANT)
#endif
                  
          IF (NZ .GT. 1) THEN
             CALL SFFTW_PLAN_DFT_C2R_3D(PLAN,NX,NY,NZ,
     &                                  AF,AF,FFTW_MEASURE)
          ELSEIF (NY .GT. 1) THEN
             CALL SFFTW_PLAN_DFT_C2R_2D(PLAN,NX,NY,
     &                                   AF,AF,FFTW_MEASURE)
          ELSE
	     CALL SFFTW_PLAN_DFT_C2R_1D(PLAN,NX, AF,AF,FFTW_MEASURE)
          ENDIF

c         write(6,91)PLAN,NX,NY,NZ,NUMTHWANT
          IF (VERBOSE) WRITE(NOUT,91)PLAN,NX,NY,NZ,NUMTHWANT
 91       FORMAT('  New Reverse FFTW3 Plan: ',I0,' ( ',
     &              I0,', ', I0,', ', I0,') Threads: ',I0)              
        ENDIF

        NUMTHINOUT = NUMTHWANT
        IRTFLG  = 0

	END




C     ----------------- FFTW3_USEPLAN -------------------------------

	SUBROUTINE FFTW3_USEPLAN(BUF, NX,NY,NZ, 
     &                           SPIDER_SIGN, SPIDER_SCALE, 
     &                           PLAN,INV,IRTFLG)

        USE TYPE_KINDS

	REAL ,INTENT(INOUT)  :: BUF(*)
        LOGICAL, INTENT(IN)  :: SPIDER_SIGN, SPIDER_SCALE
        INTEGER              :: INV
        INTEGER, INTENT(OUT) :: IRTFLG

C       PLAN IS A POINTER TO A STRUCTURE 
        INTEGER(KIND=I_8), INTENT(IN) :: PLAN

C       USING FFTW3 LIBRARY CALLS FOR FFT

        IRTFLG = 1
        IF (PLAN .EQ. 0) THEN
           CALL ERRT(101,'NO EXISTING FFTW3 PLAN',NE)
           RETURN
        ENDIF

        LDA = NX + 2 - MOD(NX,2)
        JH  = LDA / 2

        IF (INV .GT. 0) THEN
C          FORWARD TRANSFORM, USING FFTW3 GURU INTERFACE

c           write(6,90)PLAN
c           IF (VERBOSE) WRITE(NOUT,90)PLAN,NX,NY,NZ 
 90        FORMAT( '  Forward FFTW3, Plan:   ',I0 ' ( ',
     &              I0,',' , I0,', ', I0,')')

           CALL SFFTW_EXECUTE_DFT_R2C(PLAN,BUF,BUF) 

           IF (SPIDER_SIGN) THEN
C             CHANGE FFTW FORMAT TO SPIDER FFT FORMAT 
C             SPIDER IMAGINARY PARTS HAVE OPPOSITE SIGNS FROM FFTW 

c$omp         parallel do private(i)
	      DO I = 2,2*JH*NY*NZ,2	
	         BUF(I) = -BUF(I)           
   	      ENDDO
           ENDIF

        ELSE

C          REVERSE TRANSFORM, USING FFTW3 GURU INTERFACE

c          write(6,91)PLAN
c           IF (VERBOSE) WRITE(NOUT,91)PLAN,NX,NY,NZ
 91        FORMAT( '  Reverse FFTW3, Plan:   ', I0 ' ( ',
     &              I0,', ', I0,', ', I0,') ')
              
           IF (SPIDER_SIGN) THEN
C             CHANGE FFTW FORMAT TO SPIDER FFT FORMAT 
C             SPIDER IMAGINARY PARTS HAVE OPPOSITE SIGNS FROM FFTW 

c$omp         parallel do private(i)
	      DO I = 2,2*JH*NY*NZ,2	
	         BUF(I) = -BUF(I)           
 	      ENDDO
           ENDIF

C          USING FFTW3 GURU INTERFACE FOR BACK TRANSFORM
           CALL SFFTW_EXECUTE_DFT_C2R(PLAN,BUF,BUF)

           IF (SPIDER_SCALE) THEN
C             SCALING NEEDED
              PIX = 1.0 / FLOAT (NX * NY * NZ)

c$omp         parallel do private(i)
              DO I=1,LDA*NY*NZ
                 BUF(I) = BUF(I) * PIX
              ENDDO
           ENDIF

        ENDIF
        IRTFLG = 0

	END

C     ----------------- FFTW3_KILLPLAN -------------------------------

	SUBROUTINE FFTW3_KILLPLAN(PLAN,IRTFLG)

        USE TYPE_KINDS

C       PLAN IS A POINTER TO A STRUCTURE 
        INTEGER(KIND=I_8), INTENT(INOUT) :: PLAN

        INTEGER, INTENT(OUT)             :: IRTFLG

C       USING FFTW3 LIBRARY CALLS FOR FFT, DESTROY EXISTING OLD PLAN
        IF (PLAN .GT. 0) THEN
           CALL SFFTW_DESTROY_PLAN(PLAN)
           PLAN = 0
        ENDIF

        IRTFLG = 0
        END

C     ----------------- FFTW3_KILLPLANS-------------------------------

	SUBROUTINE FFTW3_KILLPLANS(PLANS,NPLANS,IRTFLG)

        USE TYPE_KINDS

C       PLAN IS A POINTER TO A STRUCTURE 
        INTEGER(KIND=I_8), INTENT(INOUT) :: PLANS(NPLANS)

        INTEGER, INTENT(OUT)             :: IRTFLG

C       USING FFTW3 LIBRARY CALLS FOR FFT, DESTROY EXISTING OLD PLAN
        DO N = 1,NPLANS
           IF (PLANS(N) .GT. 0) THEN
              CALL SFFTW_DESTROY_PLAN(PLANS(N))
              PLANS(N) = 0
           ENDIF
        ENDDO

        IRTFLG = 0
        END




C     ----------------- FFTW3_INIT_THREADS ---------------------------

	SUBROUTINE FFTW3_INIT_THREADS(IRTFLG)

        INCLUDE 'CMBLOCK.INC'

#if defined (SP_MP) 

        IF (NUMFFTWTH .EQ. 0) THEN
C          FIRST TIME FFTW3 CALLED, SET NUMFFTWTH TO ALL POSSIBLE
           CALL GETTHREADS(NUMFFTWTH)

C          IF (VERBOSE) WRITE(NOUT,92) NUMFFTWTH
 92        FORMAT( '  Number of OMP threads possible:', I0)
        ENDIF
#else
        NUMFFTWTH = 1
#endif

        IRTFLG = 0
        END

