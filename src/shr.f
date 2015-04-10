
C++************************************************************* 1/22/81
C
C SHR.F 
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
C SHR(LUN1,LUN2,NSAM,NROW)
C
C PURPOSE:  THRESHOLD IMAGE FROM OUTSIDE
C
C **********************************************************************

	SUBROUTINE SHR(LUN1,LUN2,NSAM,NROW)

	COMMON ADUM(80),BUF(1)
	COMMON/FUNCTION/FCHAR

        CHARACTER *80 FCHAR
        CHARACTER     NULL,FOPT

        NULL = CHAR(0)

	CALL RDPRMC(FOPT,NCHAR,.TRUE.,
     &    'RETAIN PIXELS (A)BOVE OR (B)ELOW THRESHOLD? (A/B)',
     &    NULL,IRTFLG)

	CALL RDPRM(THR,NOT_USED,'THRESHOLD')
	FACT = 1.0
	IF (FCHAR(4:4) .EQ. 'M') FACT = 0.
	IF (FOPT .EQ. 'B') GOTO 200

	DO 100 I=1,NROW
          CALL REDLIN(LUN1,BUF,NSAM,I)
          IF (FCHAR(4:4) .EQ. 'M') THEN
            DO  K=1,NSAM
              BUF(NSAM+K)=1.
	    ENDDO
          ENDIF
16        DO  K=1,NSAM
            KK=K
            B = BUF(K)
            IF (B .LT. THR) BUF(K+NSAM)=THR*FACT
            IF (B .GE. THR) GOTO 11
	  ENDDO
          GOTO 12

11        IF (KK .EQ. NSAM) GOTO 100
12	  DO  K=NSAM,1,-1
            KK=K
            B=BUF(K)
            IF(B.LT.THR)BUF(K+NSAM)=THR*FACT
            IF(B.GE.THR)GOTO 100
	  ENDDO
100	CALL WRTLIN(LUN2,BUF(NSAM+1),NSAM,I)
	RETURN

200	DO 300 I=1,NROW
          CALL REDLIN(LUN1,BUF,NSAM,I)
          IF (FCHAR(4:4) .EQ. 'M') THEN
            DO  K=1,NSAM
              BUF(NSAM+K)=1.
	    ENDDO
          ENDIF
220	  DO  K=1,NSAM
            KK=K
            B=BUF(K)
            IF (B.GT.THR)BUF(K+NSAM)=THR*FACT
            IF (B.LE.THR)GOTO 321
	  ENDDO
          GOTO 322

321	  IF (KK.EQ.NSAM) GOTO 300
322       DO  K=NSAM,1,-1
            KK=K
            B=BUF(K)
            IF (B.GT.THR)BUF(K+NSAM)=THR*FACT
            IF (B.LE.THR)GOTO 300
	  ENDDO

300	CALL WRTLIN(LUN2,BUF(NSAM+1),NSAM,I)

	END

