
C++*********************************************************************
C
C FMRS.F      CREATED FROM FMRS_1 & _2 & _3.F    JAN 2008 ARDEAN LEITH
C             ALTERED                            FEB 2008 ARDEAN LEITH
C             FMRS_PLANB                         APR 2008 ARDEAN LEITH
C             FMRS_PLANB 16TH SLOT BUG           DEC 2008 ARDEAN LEITH
C             IPD(6  FOR -PLAN BUG               MAR 2009 ARDEAN LEITH
C             REMOVED OLD SGI CODE               OCT 2010 ARDEAN LEITH
C             REMOVED SPIDER FFT CALLS           DEC 2010 ARDEAN LEITH
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C  FMRS(BUF, NX,NY,NZ, FFTW3PLAN, 
C       SPIDER_SIGN, SPIDER_SCALE, INV,IRTFLG)
C
C  PARAMETERS: BUF            REAL ARRAY (LDA*NY*NZ)     SENT/RET.
C                                LDA = NX + 2 - MOD(NX,2)
C                                EVEN NX --> EVEN LDA = NX +2
C                                ODD  NX --> EVEN LDA = NX +1
C              NX..           NON-FOURIER DIMENSIONS               SENT
C              FFTW3PLAN      EXISTING FFTW3 PLAN OR ZERO          SENT
C              SPIDER_SIGN    CHANGE SIGN OF IMAGINARY TO SPIDER   SENT
C              SPIDER_SCALE   SCALE THE OUTPUT FLAG                SENT
C              INV            1=REG. DATA, -1= FOURIER DATA        SENT
C              IRTFLG         ERROR FLAG (O=NORMAL)                RET.
C
C  PURPOSE:  REAL MIXED RADIX FFT.
C
C            OUTPUT: 
C            N=EVEN  BUF(N+2)
C               ORDER OF ELEMENTS IN BUF:
C               R(0),0.0, R(1), I(1), R(2), I(2), ....., 
C                     R(N/2-1), I(N/2-1), R(N/2),0.0
C
C            N=ODD  BUF(N+1)
C               ORDER OF ELEMENTS IN BUF:
C                    R(0),0.0, R(1), I(1), R(2), I(2), ....., 
C                    R(N/2-1), I(N/2-1), R(N/2),I(N/2)
C
C            FOLLOWING CONVENTION THAT INTEGER DIVISION 
C            IS ROUNDED DOWN, E.G. 5/2 =2)
C
C--*********************************************************************

      MODULE FMRS_INFO
        USE TYPE_KINDS      
        INTEGER, PARAMETER                   :: NPLANS = 16
        INTEGER(KIND=I_8), DIMENSION(NPLANS) :: PLANS  = 0
        INTEGER, DIMENSION(NPLANS,6),SAVE    :: IPD
      END MODULE FMRS_INFO



C       -------------------- FMRS ------------------------------------

	SUBROUTINE FMRS(BUF, NX,NY,NZ, FFTW3PLAN, 
     &                  SPIDER_SIGN, SPIDER_SCALE, INV,IRTFLG)

        USE TYPE_KINDS      
        INTEGER(KIND=I_8), INTENT(IN) :: FFTW3PLAN  ! STRUCTURE POINTER 
        INTEGER(KIND=I_8)             :: PLAN       ! STRUCTURE POINTER 

	REAL,    INTENT(INOUT) :: BUF(*)
        INTEGER, INTENT(IN)    :: NX,NY,NZ
        LOGICAL, INTENT(IN)    :: SPIDER_SIGN, SPIDER_SCALE 
        INTEGER, INTENT(INOUT) :: INV
        INTEGER, INTENT(OUT)   :: IRTFLG

        INTEGER                :: OMP_GET_NUM_THREADS
        INTEGER                :: OMP_GET_NUM_PROCS
        LOGICAL                :: USEBUF

        INCLUDE 'CMBLOCK.INC'

        LDA = NX + 2 - MOD(NX,2)
        IF (FFTW3PLAN .GT. 0) THEN 
C          USE PLAN SENT FROM CALLER
           PLAN = FFTW3PLAN
        ELSE

           IF (NUMFFTWTH .LE. 0) THEN
C             FIRST TIME FFTW3 USED IN THIS RUN
	      CALL FFTW3_INIT_THREADS(IRTFLG)
           ENDIF

           NUMTHWANT = NUMFFTWTH
#ifdef SP_MP
C          COMPILED FOR OMP USE

C          NUMTHP = OMP_GET_NUM_PROCS()
           NUMTH = OMP_GET_NUM_THREADS()
           IF (NUMTH .GT. 1) THEN
C             INSIDE OMP PARALLEL SECTION, CAN NOT MAKE A NEW PLAN
              NUMTHWANT = 1
           ENDIF
           !write(6,*) 'in fmrs, threads,numthwant: ',numth,numthwant
#endif
                   
           !write(6,908) NX,NY,NZ,NUMTHWANT,INV
 908       format( ' Calling fmrs_planb for: (',
     &              I5,',', I5,',', I5,') ',2i5)

C          USE CACHED PLAN OR CREATE A NEW ONE IF CACHE NOT USEFULL
           USEBUF = .FALSE.
	   CALL FMRS_PLANB(.FALSE.,FDUM,NX,NY,NZ,NUMTHWANT,
     &                     INV,PLAN,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ENDIF

        CALL FFTW3_USEPLAN(BUF, NX,NY,NZ, 
     &                     SPIDER_SIGN, SPIDER_SCALE, 
     &                     PLAN,INV,IRTFLG)

	END





C       -------------------- FMRS_PLAN ------------------------------------

        SUBROUTINE FMRS_PLAN(USEBUF,BUF, NX,NY,NZ,
     &                       NUMTHWANT,INV,IRTFLG)

        LOGICAL, INTENT(IN)            :: USEBUF
	REAL,    INTENT(IN)            :: BUF(*)
        INTEGER, INTENT(IN)            :: NX,NY,NZ
        INTEGER, INTENT(IN)            :: NUMTHWANT
        INTEGER, INTENT(IN)            :: INV
        INTEGER, INTENT(OUT)           :: IRTFLG

        INTEGER *8                     :: PLAN  !STRUCTURE POINTER

C       USING FFTW3 LIBRARY CALLS FOR FFT 

C       CREATE A NEW PLAN IF CACHE NOT USEFULL
C       DOES NOT RETURN PLAN TO CALLER

	CALL FMRS_PLANB(USEBUF,BUF,NX,NY,NZ,
     &                     NUMTHWANT,INV,PLAN,IRTFLG)

        RETURN
 	END




C       -------------------- FMRS_PLANB -------------------------------

       SUBROUTINE FMRS_PLANB(USEBUF,BUF, NX,NY,NZ,
     &                       NTH,INV,PLAN,IRTFLG)

C       USING FFTW3 LIBRARY CALLS FOR FFT 

        USE TYPE_KINDS      
        USE FMRS_INFO      

        LOGICAL, INTENT(IN)            :: USEBUF
	REAL,    INTENT(IN)            :: BUF(*)
        INTEGER, INTENT(IN)            :: NX,NY,NZ
        INTEGER, INTENT(IN)            :: NTH
        INTEGER, INTENT(IN)            :: INV
        INTEGER(KIND=I_8), INTENT(OUT) :: PLAN   !STRUCTURE POINTER 
        INTEGER, INTENT(OUT)           :: IRTFLG

        INCLUDE 'CMBLOCK.INC'

        IRTFLG = 1

        MT         = -1              ! EMPTY PLAN SLOT
        IREP       = -1              ! REPLACEABLE PLAN SLOT
        NUMTHINOUT = NTH

C       SEE IF THERE IS ANY SUITABLE CACHED PLAN
        DO IPLAN = 1,NPLANS
           PLAN = PLANS(IPLAN)
           IF (PLAN == 0 .AND. MT < 0) THEN
C             AN MT SLOT FOR A PLAN
              MT = IPLAN
              !write(6,934) mt
934           format('  Found MT slot: ',i3)

           ELSEIF (PLAN .NE. 0) THEN
C             THIS IS A VALID CACHED PLAN

              IF (NX     == IPD(IPLAN,1) .AND. 
     &            NY     == IPD(IPLAN,2) .AND.
     &            NZ     == IPD(IPLAN,3) .AND. 
     &            NTH    == IPD(IPLAN,4) .AND.
     &            INV    == IPD(IPLAN,5)) THEN
C                THIS IS A SUITABLE CACHED PLAN 

                 !write(6,903) iplan, (ipd(iplan,i),i=1,6), plan
 903             format('  Found suitable plan: ',i2,': (',i5,','
     &                  ,i5,',',i5,'):', i2,' ',i2,' ',i2,' Plan:',i14)

C                PRESERVE CACHED PLAN FOR THIS OPERATIONS USE
                 IF (IPD(IPLAN,6) < 0)  IPD(IPLAN,6) = 1

                 IRTFLG = 0
                 RETURN            ! CAN USE THIS PLAN

              ELSEIF (IPD(IPLAN,6) < 0 .AND. IREP < 0) THEN
C                NOT SUITABLE, BUT PLAN CAN BE REPLACED FOR THIS OPERATION
                 IREP = IPLAN

                 !write(6,901) iplan, (ipd(iplan,i),i=1,6), plan
 901             format('  Found replaceble plan: ',i2,': (',i5,','
     &                  ,i5,',',i5,'):', i2,' ',i2,' ',i2,' Plan:',i14)
              ENDIF
           ENDIF
        ENDDO

C       CREATE A NEW CACHED PLAN

        IF (MT > 0) THEN
C          HAVE AN EMPTY PLAN SLOT
           IPLAN = MT
           PLAN  = PLANS(IPLAN)

        ELSEIF (IREP > 0) THEN
C          HAVE AN REPLACEABLE PLAN SLOT
           IPLAN = IREP
           PLAN  = PLANS(IPLAN)

        ELSE
C          NO REMAINING PLAN SLOT FOR ANOTHER CACHED PLAN
        
           WRITE(NOUT,*) '  PROGRAMMING ERROR, TELL THE PROGRAMMER'
           DO IPLAN = 1,NPLANS
              WRITE(6,900) IPLAN,(IPD(IPLAN,I),I=1,6), PLANS(IPLAN)
 900          FORMAT('  CACHED PLAN: ',I2,': (',I5,',',I5,',',I5,'):',
     &                  I2,' ',I2,' ',I2,' PLAN:',I14)
           ENDDO

           CALL ERRT(102,'AVAILABLE FFTW3 PLANS LIMITED TO',NPLANS)
           IRTFLG = 1
           RETURN
        ENDIF

C       THIS IS A EMPTY OR REPLACABLE PLAN SLOT, FILL IT

        LDA = NX + 2 - MOD(NX,2)
        IF (USEBUF) THEN
c          write(6,899)NX,NY,NZ,numthinout,inv
c899       format('  Calling makeplanb for: ',3i6,i4,i3)
           CALL FFTW3_MAKEPLANB(BUF,LDA,NX,NY,NZ,
     &                          NUMTHINOUT,PLAN,INV,IRTFLG)
        ELSE
c          write(6,890)NX,NY,NZ,numthinout,inv
c890       format('  Calling makeplan for: ',3i6,i4,i3)
           CALL FFTW3_MAKEPLAN(NX,NY,NZ,
     &                         NUMTHINOUT,PLAN,INV,IRTFLG)
        ENDIF
        IF (IRTFLG .NE. 0) RETURN

        !write(6,902),iplan,plan,mt,irep
902     format('  Slot:',i4,' Plan:',i20,'  mt,irep:',2i5)

        PLANS(IPLAN) = PLAN
        IPD(IPLAN,1) = NX
        IPD(IPLAN,2) = NY
        IPD(IPLAN,3) = NZ
        IPD(IPLAN,4) = NUMTHINOUT
        IPD(IPLAN,5) = INV
        IPD(IPLAN,6) = 1     
                 
        END


C       -------------------- FMRS_DEPLAN -------------------------------

        SUBROUTINE FMRS_DEPLAN(IRTFLG)

C       USING FFTW3 LIBRARY CALLS FOR FFT 

        USE FMRS_INFO      

        INTEGER, INTENT(OUT)   :: IRTFLG

C       DE-PLAN ALL CACHED PLANS
        DO IPLAN = 1,NPLANS
           IF (PLANS(IPLAN) .NE. 0) IPD(IPLAN,6) = -1 
        ENDDO
        IRTFLG = 0
 
        END






#ifdef NEVER
C-------------------------- old unused -----------------------

C       -------------------- FMRS_PLANB -------------------------------

        SUBROUTINE FMRS_PLANB(USEBUF,BUF,NX,NY,NZ,
     &                       NTH,INV,PLAN,IRTFLG)

C       USING FFTW3 LIBRARY CALLS FOR FFT 

        USE TYPE_KINDS      

        LOGICAL, INTENT(IN)            :: USEBUF
	REAL,    INTENT(IN)            :: BUF(*)
        INTEGER, INTENT(IN)            :: NX,NY,NZ
        INTEGER, INTENT(IN)            :: NTH
        INTEGER, INTENT(IN)            :: INV
        INTEGER(KIND=I_8), INTENT(OUT) :: PLAN   !STRUCTURE POINTER 
        INTEGER, INTENT(OUT)           :: IRTFLG

        INTEGER(KIND=I_8), SAVE :: PLAN1O=0, PLAN1OR=0
        INTEGER(KIND=I_8), SAVE :: PLAN2O=0, PLAN2OR=0
        INTEGER(KIND=I_8), SAVE :: PLAN3O=0, PLAN3OR=0

        INTEGER, SAVE :: NX3O=0,  NY3O=0, NZ3O=0
        INTEGER, SAVE :: NX2O=0,  NY2O=0
        INTEGER, SAVE :: NX1O=0

        INTEGER, SAVE :: NX3OR=0, NY3OR=0, NZ3OR=0
        INTEGER, SAVE :: NX2OR=0, NY2OR=0
        INTEGER, SAVE :: NX1OR=0

        INTEGER, SAVE :: NUMFFTWTHO1=-1, NUMFFTWTHOR1=-1
        INTEGER, SAVE :: NUMFFTWTHO2=-1, NUMFFTWTHOR2=-1
        INTEGER, SAVE :: NUMFFTWTHO3=-1, NUMFFTWTHOR3=-1

        INCLUDE 'CMBLOCK.INC'

        IRTFLG    = 1
        LDA       = NX + 2 - MOD(NX,2)
        NUMTHWANT = NTH

        IF (INV .GT. 0) THEN
C          FORWARD TRANSFORM

           IF (NZ .GT. 1) THEN
C             3D TRANSFORM
              IF (NX     .NE. NX3O .OR. 
     &           NY      .NE. NY3O .OR. 
     &           NZ    .NE. NZ3O .OR.
     &           NUMTHWANT .NE. NUMFFTWTHO3) THEN
C               SIZE OR THREADING CHANGED, CREATE NEW FORWARD PLAN

                !write(6,899)NX,NY,NZ,numthwant,inv
 899            format(' Calling makeplan for: ',3i6,i4,i2)
 
                IF (USEBUF) THEN
                   CALL FFTW3_MAKEPLANB(BUF,LDA,NX,NY,NZ,
     &                               NUMTHWANT,PLAN3O,INV,IRTFLG)
                ELSE
                   CALL FFTW3_MAKEPLAN(NX,NY,NZ,
     &                               NUMTHWANT,PLAN3O,INV,IRTFLG)
                ENDIF
                IF (IRTFLG .NE. 0) RETURN

                NX3O      = NX
                NY3O      = NY
                NZ3O    = NZ
                NUMFFTWTHO3 = NUMTHWANT
             ENDIF
             PLAN = PLAN3O

          ELSEIF  (NY .GT. 1) THEN
C             2D TRANSFORM
              IF (NX     .NE. NX2O .OR. 
     &           NY      .NE. NY2O .OR. 
     &           NUMTHWANT .NE. NUMFFTWTHO2) THEN
C               SIZE OR THREADING CHANGED, CREATE NEW FORWARD PLAN

                !write(6,899)NX,NY,NZ,numthwant,inv

                IF (USEBUF) THEN
                   CALL FFTW3_MAKEPLANB(BUF,LDA,NX,NY,NZ,
     &                               NUMTHWANT,PLAN2O,INV,IRTFLG)
                ELSE
                   CALL FFTW3_MAKEPLAN(NX,NY,NZ,
     &                               NUMTHWANT,PLAN2O,INV,IRTFLG)
                ENDIF
                IF (IRTFLG .NE. 0) RETURN

                NX2O      = NX
                NY2O      = NY
                NUMFFTWTHO2 = NUMTHWANT
             ENDIF

             PLAN = PLAN2O

          ELSE 
C             1D TRANSFORM
              IF (NX     .NE. NX1O .OR. 
     &           NUMTHWANT .NE. NUMFFTWTHO1) THEN
C                SIZE OR THREADING CHANGED, CREATE NEW FORWARD PLAN

                !write(6,899)NX,NY,NZ,numthwant,inv

                IF (USEBUF) THEN
                   CALL FFTW3_MAKEPLANB(BUF,LDA,NX,NY,NZ,
     &                               NUMTHWANT,PLAN1O,INV,IRTFLG)
                ELSE
                   CALL FFTW3_MAKEPLAN(NX,NY,NZ,
     &                               NUMTHWANT,PLAN1O,INV,IRTFLG)
                ENDIF
                IF (IRTFLG .NE. 0) RETURN

                NX1O      = NX
                NUMFFTWTHO1 = NUMTHWANT
             ENDIF

             PLAN = PLAN1O
          ENDIF

       ELSE

C          REVERSE TRANSFORM
           IF (NZ .GT. 1) THEN
C             3D TRANSFORM
              IF (NX      .NE. NX3OR .OR. 
     &            NY      .NE. NY3OR .OR. 
     &            NZ    .NE. NZ3OR .OR.
     &            NUMTHWANT .NE. NUMFFTWTHOR3) THEN
C                SIZE OR THREADING CHANGED, CREATE NEW FORWARD PLAN

                 !write(6,899)NX,NY,NZ,numthwant,inv

                 IF (USEBUF) THEN
                   CALL FFTW3_MAKEPLANB(BUF,LDA,NX,NY,NZ,
     &                               NUMTHWANT,PLAN3OR,INV,IRTFLG)
                 ELSE
                   CALL FFTW3_MAKEPLAN(NX,NY,NZ,
     &                               NUMTHWANT,PLAN3OR,INV,IRTFLG)
                 ENDIF
                 IF (IRTFLG .NE. 0) RETURN

                 NX3OR        = NX
                 NY3OR        = NY
                 NZ3OR        = NZ
                 NUMFFTWTHOR3 = NUMTHWANT
             ENDIF

             PLAN = PLAN3OR

          ELSEIF  (NY .GT. 1) THEN
C             2D TRANSFORM
              IF (NX       .NE. NX2OR .OR. 
     &            NY       .NE. NY2OR .OR. 
     &           NUMTHWANT .NE. NUMFFTWTHOR2) THEN
C               SIZE OR THREADING CHANGED, CREATE NEW FORWARD PLAN

                !write(6,899)NX,NY,NZ,numthwant,inv

                IF (USEBUF) THEN
                   CALL FFTW3_MAKEPLANB(BUF,LDA,NX,NY,NZ,
     &                               NUMTHWANT,PLAN2OR,INV,IRTFLG)
                ELSE
                   CALL FFTW3_MAKEPLAN(NX,NY,NZ,
     &                               NUMTHWANT,PLAN2OR,INV,IRTFLG)
                ENDIF
               IF (IRTFLG .NE. 0) RETURN

                NX2OR      = NX
                NY2OR      = NY
                NUMFFTWTHOR2 = NUMTHWANT
             ENDIF

             PLAN = PLAN2OR

           ELSE 
C             1D TRANSFORM
              IF (NX       .NE. NX1OR .OR. 
     &           NUMTHWANT .NE. NUMFFTWTHOR1) THEN
C               SIZE OR THREADING CHANGED, CREATE NEW FORWARD PLAN

                !write(6,898) NX, NX1or, numthwant ,numfftwthor1
 898            format(' NX,NX1or, numthwant,numfftwthor1:',4i7)
                !write(6,899)NX,NY,NZ,numthwant,inv

                IF (USEBUF) THEN
                   CALL FFTW3_MAKEPLANB(BUF,LDA,NX,NY,NZ,
     &                               NUMTHWANT,PLAN1OR,INV,IRTFLG)
                ELSE
                   CALL FFTW3_MAKEPLAN(NX,NY,NZ,
     &                               NUMTHWANT,PLAN1OR,INV,IRTFLG)
                ENDIF
                IF (IRTFLG .NE. 0) RETURN

                NX1OR      = NX
                NUMFFTWTHOR1 = NUMTHWANT
             ENDIF

             PLAN = PLAN1OR
          ENDIF

        ENDIF
        IRTFLG = 0

        END

#endif
