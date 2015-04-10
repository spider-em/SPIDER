
C++*********************************************************************
C
C    LOCAL.F                                        J.FRANK JULY 1977
C
C    PROGRAM SEEMS TO HAVE BEEN CHANGED IN SEPT. 86
C    PROGRAM CORRECTED TO BE COMPILABLE (BUT NOT TESTED) 7/30/87 M.R.
C    REVISED SEPT 89 AL-- DIDN'T WORK
C
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
C LOCAL(LUN1,LUN2,NX,NY)  
C
C PURPOSE:  CONTRAST ENHANCEMENT BASED ON LOCAL HISTOGRAM INFORMATION
C
C PARAMETERS:
C             LUN1       LOGICAL UNIT NUMBER OF FILE
C             LUN2       LOGICAL UNIT NUMBER OF FILE
C             NX,NY  DIMENSIONS OF FILE
C
C VARIABLES:  NLOCAL     LOCAL EQUALIZATION AREA
C             KCTR1,2    STARTING & ENDING COLUMNS
C             ICTR1,2    STARTING & ENDING ROWS
C             A1         NO. PIXELS IN LOCAL AREA
C             IM         POINTER TO BUFFER POSITION OF ROW START
C             MAP        POISTION OF CURRENT PIXEL IN HISTOGRAM
C             NPTR       ARRAY OF POINTERS TO ROW POSITION IN BUFF
C          
C--*******************************************************************

      SUBROUTINE LOCAL(LUN1,LUN2,NX,NY,NZ)

      
      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      COMMON BUFF(1)   ! ACTUAL SIZE NLOCAL * NX???

      REAL             :: H(256) 
      INTEGER          :: NPTR(64) 
      CHARACTER(LEN=1) :: NULL = CHAR(0)
      CHARACTER(LEN=1) :: ANS
      INTEGER          :: A1

      PARAMETER    (NLMAX = 64)

      CALL RDPRMC(ANS,NC,.TRUE.,
     &     'GENERALIZED HISTOGRAM OR THRESHOLD (G/T)',NULL,IRT)
      IF (IRT .NE. 0) RETURN

      CALL RDPRMI(NLOCAL,NDUM,NOT_USED,'LOCAL AREA SIZE')

      NLM = NLMAX
      IF (NLOCAL > NLM) THEN
         WRITE(NOUT,11) NLM,NLM
11       FORMAT('  MAXIMUM AREA SIZE: ',I2,' * ',I2,'  ASSUMED')
         NLOCAL = NLM
      ENDIF

      NLH = NLOCAL/2
      IF (IMAMI .NE. 1) CALL NORM3(LUN1,NX,NY,NZ,FMAX,FMIN,AV)

      NL1   = NLH + 1
      ICTR1 = NL1
      KCTR1 = NL1
      ICTR2 = NY*NZ-NL1
      ICTR3 = ICTR2+1
      KCTR2 = NX-NL1
      A1    = NLOCAL*NX
      FNALL = NLOCAL**2
      SCAL  = 2./FNALL
      HINCR = (FMAX-FMIN)/127.

C     TEMPORARY CLEARING OF UPPER AND LOWER MARGIN  BECAUSE OF 
C     DISCREPANCY OF NORMALIZATION

      DO  K = 1,A1 + NX
         BUFF(K) = 0.
      ENDDO

      DO  I = 1,NLH
         CALL WRTLIN(LUN2,BUFF,NX,I)
      ENDDO

C     INITIALIZE LHIST
      MODE = 0
      CALL LHIST(BUFF,NX,NLOCAL,KCTR,NPTR,H,MODE)
      IF (MODE .NE. 0) RETURN

C     INITIALIZE MREAD
      CALL MREAD(-1,BUFF,NX,NLOCAL,NPTR)

      IF (ANS .NE. 'T') THEN
C       LOCAL GENERALIZED HISTOGRAM OPERATION

        DO  ICTR = ICTR1,ICTR2
C         SWITCH LHIST TO NON-INCREMENTAL OPERATION WITH INTEGRATION
          MODE = 3

C         READ NLOCAL LINES INTO BUFF
          CALL MREAD(LUN1,BUFF,NX,NLOCAL,NPTR)

          IM = (NPTR(NL1)-1)*NX

          DO  KCTR = KCTR1,KCTR2
            CALL LHIST(BUFF,NX,NLOCAL,KCTR,NPTR,H,MODE)

C           SWITCH LHIST TO INCREMENTAL OPERATION WITH INTEGRATION
            MODE = 4

C           NOW APPLY MAPPING TO EACH PIXEL FROM KCTR1...KCTR2
            MAP           = (BUFF(IM+KCTR)-FMIN)/HINCR+1.5
            BUFF(A1+KCTR) = H(128+MAP)*SCAL
	  ENDDO

C         WRITE THIS ROW TO OUTPUT FILE
          CALL WRTLIN(LUN2,BUFF(A1+1),NX,ICTR)
	ENDDO

      ELSE

        CALL RDPRM(PERC,NOT_USED,'HISTOGRAM THRESHOLD PERCENTAGE')

        THRESH = FNALL * PERC/100.
        DO  ICTR = ICTR1,ICTR2
          MODE = 1

C         SWITCH LHIST TO NON-INCREMENTAL OPERATION
          CALL MREAD(LUN1,BUFF,NX,NLOCAL,NPTR)

          IM = (NPTR(NL1)-1)*NX

          DO  KCTR = KCTR1,KCTR2
            CALL LHIST(BUFF,NX,NLOCAL,KCTR,NPTR,H,MODE)

C           SWITCH LHIST TO INCREMENTAL OPERATION
            MODE = 2

            DO  L = 1,128
               IF (H(L) > THRESH) GOTO 86
	    ENDDO

86          IF (FLOAT(L)*HINCR < BUFF(IM+KCTR)-FMIN) BUFF(A1+KCTR) = 2.
	  ENDDO

          CALL WRTLIN(LUN2,BUFF(A1+1),NX,ICTR)
	ENDDO
      ENDIF


C     BORDER CLEARING. 
      DO  K = 1,NX
         BUFF(K) = 0.
      ENDDO

      DO  I = ICTR3,NY*NZ
         CALL WRTLIN(LUN2,BUFF,NX,I)
      ENDDO

      IMAMI = 0
      SIG   = -1.
      FMAX  = 2.
      FMIN  = 0.
      CALL SETPRM(LUN2,NX,NY,FMAX,FMIN,AV,'U')

      END

