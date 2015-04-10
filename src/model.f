
C++*********************************************************************
C
C MODEL.F   REWRITTEN                            APR 97 ARDEAN LEITH
C           RDPRAF REMOVED                       DEC 05 ARDEAN LEITH 
C           G2..                                 JAN 12 ARDEAN LEITH 
C           C RADIUS BUG                         JAN 12 ARDEAN LEITH 
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
C   MODEL(LUN1,NX,NY)
C
C   PURPOSE:  PREPARES TEST PICTURES
C
C   PARAMETERS:
C        LUN1        OUTPUT UNIT NUMBER OF FILE
C        NX,NY       DIMENSIONS OF FILE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE MODEL(LUN1,NX,NY)

	INCLUDE 'CMBLOCK.INC'

        INTEGER            :: LUN1,NX,NY,NZ

	INTEGER, PARAMETER :: MAXSIN = 50
	INTEGER            :: NCX(MAXSIN),NCY(MAXSIN)
	REAL               :: CX(MAXSIN)
	REAL               :: PH(MAXSIN),FWA(4)
	REAL               :: A0(NX)
	REAL               :: SCF, FI

	INTEGER            :: IRTFLG,NVAL
	CHARACTER(LEN=1)   :: GA
	CHARACTER(LEN=2)   :: ANS
	CHARACTER(LEN=1)   :: NULL = CHAR(0)

	REAL, PARAMETER    :: QUADPI = 3.1415926535897932
	REAL, PARAMETER    :: TWOPI  = 2 * QUADPI

	IF (COPT == 'I' .OR. VERBOSE) WRITE(NOUT,100)
 100    FORMAT(
     &      ' .MENU: B   -- BLANK, CONSTANT DENSITY IMAGE'/
     &      '        C   -- FILLED CIRCLE (FOR MASKING) '/
     &      '        G   -- GAUSSIAN CIRCLE/ELLIPSE (NOT FOR MASKING)'/
     &      '        G1  -- 1ST ORDER GAUSSIAN CIRCLE/ELLIPSE (0...1)'/
     &      '        G2  -- 2ND ORDER GAUSSIAN CIRCLE/ELLIPSE (0...1)'/
     &      '        G3  -- 3RD ORDER GAUSSIAN CIRCLE/ELLIPSE (0...1)'/
     &      '        R   -- RANDOM DENSITY PATTERN'/
     &      '        S   -- SINE WAVES'/
     &      '        T   -- TWO SINE WAVES'/
     &      '        W   -- DENSITY WEDGE'/)

1010  CALL RDPRMC(ANS,NC,.TRUE.,
     &          'OPTION (B/C/CM/G/G1/G2/G3/R/S/T/W)',NULL,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF (ANS == 'B') THEN
C       BLANK IMAGE --------------------------------------------- BLANK

        BACK = 1.0
        CALL RDPRM1S(BACK,NOT_USED,'BACKGROUND CONSTANT',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CALL BLANK(LUN1,NX,NY,1,BACK)


      ELSEIF (ANS(1:1)  ==  'C' ) THEN
C       CIRCLE ------------------------------------------------- CIRCLE 

        RAD = MIN(NX,NY) / 2.0  - 3
        CALL RDPRM1S(RAD,NOT_USED,
     &               'RADIUS (FLOATING POINT)',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        RAD2  = RAD**2
        SCF   = 2.0 / FLOAT(NY+NX)
        IXCEN = NX  / 2+1
        IYCEN = NY  / 2+1

        VAL   = 1.0

	DO I=1,NY
           FI2 = (FLOAT(I-IYCEN))**2

           DO J=1,NX
             A0(J) = 0.0
             IF (FI2 + FLOAT(J-IXCEN)**2 <= RAD2) A0(J) = VAL
          ENDDO 

          CALL WRTLIN(LUN1,A0,NX,I)
        ENDDO


      ELSEIF (ANS(1:1) == 'W') THEN
C       WEDGE --------------------------------------------------- WEDGE

        SCF = 2.0 / FLOAT(NY+NX)

	DO I=1,NY
            FI = FLOAT(I) * SCF
            DO J = 1,NX
               A0(J) = FI + SCF * FLOAT(J)
            ENDDO

           CALL WRTLIN(LUN1,A0,NX,I)
        ENDDO


      ELSEIF (ANS(1:1)  ==  'G') THEN
C       GAUSSIAN   ------------------------------------------- GAUSSIAN

        DX = (NX/2) + 1
        DY = (NY/2) + 1
	CALL RDPRM2S(DX,DY,NOT_USED,
     &    'CENTER COORDINATES X,Y (or <CR> FOR IMAGE CENTER)',
     &    IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        SX = 1
        SY = 1
	CALL RDPRM2S(SX,SY,NOT_USED,
     &               'CHARACTARISTIC RADII IN X & Y',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

	IF (SX <= 0.0 .OR. SY <= 0.0)  THEN
	   CALL ERRT(101,'RADII MUST BE > 0',NE)
	   RETURN
	ENDIF

        WRITE(NOUT,90) DX,DY, SX,SY
90      FORMAT('  CENTER: (',G8.2,',',G8.2,')    STD:',G8.2,',',G8.2) 

C       SET THE ORDER OF SUPERGAUSSIAN 
        NORDER = 1
        IF (ANS(2:2)  ==  '2') NORDER = 2
        IF (ANS(2:2)  ==  '3') NORDER = 3

	GNM = 1.0 / SX / SY / 2.0 / QUADPI
        IF (NC > 1) GNM = 1.0  ! G? not G

	TNM   = ALOG(1.0 / TINY(GNM))
	SXSQ  = SX * SX
	SYSQ  = SY * SY

        DO I = 1,NY
           DO K = 1,NX
	      EEE = 0.5 * ((K-DX)**2 / SXSQ +
     &                     (I-DY)**2 / SYSQ)
	      IF (EEE >= TNM)  THEN
	         A0(K) = 0.0
	      ELSE
	         EEE   = 0.5 * (2*EEE)**NORDER
                 A0(K) = GNM * EXP(-EEE)
	      ENDIF
           ENDDO
           CALL WRTLIN(LUN1,A0,NX,I)
        ENDDO
C
      ELSEIF (ANS(1:1)  == 'R') THEN
C       PUT RANDOM NUMBERS IN THE IMAGE ---------------------- RANDOM

        CALL RDPRMC(GA,NC,.TRUE.,'GAUSSIAN DISTRIBUTED? (Y/N)',
     &             NULL,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (GA == 'Y') THEN
           XM = 1
           SD = 1
           CALL RDPRM2S(XM,SD,NOT_USED,
     &        'MEAN AND STANDARD DEVIATION OF GAUSSIAN DIST.',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ENDIF

        DO I = 1,NY
           IF (GA == 'Y') THEN
              DO K = 1,NX
                A0(K) = RANN(XM,SD)
              ENDDO
           ELSE
	      CALL RANDOM_NUMBER(HARVEST=A0)
           ENDIF
           CALL WRTLIN(LUN1,A0,NX,I)
        ENDDO

      ELSEIF (ANS == 'S' .OR. ANS  ==  'T') THEN

C       PUT SINE WAVES OF INTENSITY IN THE IMAGE ----------------- SINE
        IF (ANS(1:1) == 'T') THEN
C          ONLY ONE WAVE IS WANTED
           NS     = 1
           CX(1)  = 1.
           NCX(1) = 2
           NCY(1) = 2
           PH(1)  = 0.0

        ELSEIF (ANS(1:1) == 'S') THEN

C         FIND NUMBER OF SINE WAVES WANTED IN THE IMAGE 
          NS = 3
2         CALL RDPRI1S(NS,NOT_USED,'NUMBER OF SINE WAVES',IRTFLG)

          IF (NS > MAXSIN) THEN
             WRITE(NOUT,251) MAXSIN
251          FORMAT(' *** RESTRICTED TO',I3,' SINE WAVES')
             NS = MAXSIN
          ENDIF

          DO I = 1,NS
9006        IRTFLG = 0
            CALL RDPRA('AMPLITUDE, PHASE, SPATIAL FREQ. IN (X,Y)',
     &                 4,0,.FALSE.,FWA,NVAL,IRTFLG)
            IF (IRTFLG  .NE.  0) RETURN

            CX(I) = FWA(1)
            PH(I) = FWA(2)
            ANX   = FWA(3)
            ANY   = FWA(4)
	    
                           NCX(I) = ANX+0.5
            IF (ANX < 0.0) NCX(I) = ANX-0.5
                           NCY(I) = ANY+0.5
            IF (ANY < 0.0) NCY(I) = ANY-0.5

            WRITE(NOUT,9008)I,CX(I),PH(I),NCX(I),NCY(I)
9008        FORMAT(1X,I5,2F8.2,2I6)

            PH(I) = PH(I) * TWOPI / 360.0
          ENDDO
        ENDIF

C       PLACE SINE WAVE(S) OF INTENSITY IN THE IMAGE
        DO I=1,NY

           PHASE = FLOAT(I-1) * TWOPI / FLOAT(NY)

           DO J=1,NX
              A0(J) = 0.0
              DO K  = 1,NS
                 A0(J)= CX(K) * SIN(FLOAT(J-1) * TWOPI * 
     &                  FLOAT(NCX(K)) / FLOAT(NX) +
     &                  PHASE * FLOAT(NCY(K)) + PH(K)) + A0(J)
              ENDDO
           ENDDO
           CALL WRTLIN(LUN1,A0,NX,I)
         ENDDO
      ENDIF

      END
