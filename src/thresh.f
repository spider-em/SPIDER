C++*********************************************************************
C
C    THRESH.F   
C                     ADDED NREPL                 OCT 2007 ARDEAN LEITH  
C                     NREPL BUG                   OCT 2009 ARDEAN LEITH
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
C    THRESH(LUN1,LUN2,NSAM,NROW,NSLICE)
C
C    PURPOSE: THRESHOLD AN IMAGE FILE
C
C--********************************************************************

	SUBROUTINE THRESH(LUN1,LUN2,NSAM,NROW,NSLICE)

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        COMMON /IOBUF/ BUF(NBUFSIZ)

        CHARACTER(LEN=1) :: ANS,NULL
	REAL             :: B(2,2)
        INTEGER * 8      :: NREPL

	DATA B/0.0,1.0,1.0,0.0/

        NREPL = 0
        FREPL = 0.0

        IF (FCHAR(4:4) .EQ. 'M') THEN
           CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &     'BLANK OUT  (A)BOVE OR (B)ELOW THRESHOLD? (A/B)',NULL,IRT)
           IF (IRT .NE. 0) RETURN

           THR = 1.0
	   CALL RDPRM1S(THR,NOT_USED,'THRESHOLD',IRT)
           IF (IRT .NE. 0) RETURN

 	   ISW = 1
	   IF (ANS .EQ. 'A') ISW = 2

           DO  J=1, NSLICE     
	      DO  I=1,NROW
	         CALL REDLIN(LUN1,BUF,NSAM,(J-1)*NROW+I)

	         DO  K=1,NSAM
	            IF (BUF(K) .LT. THR) THEN 
	               BUF(K) = B(1,ISW)
 	            ELSE 
	               BUF(K) = B(2,ISW)
	            ENDIF
                    IF (BUF(K) > 0.0) NREPL = NREPL + 1
	         ENDDO

                 CALL WRTLIN(LUN2,BUF,NSAM,(J-1)*NROW+I)
	      ENDDO
	   ENDDO
           FREPL = NREPL
           CALL REG_SET_NSEL(1,1, FREPL,0.0, 0.0, 0.0, 0.0,IRTFLG)
           WRITE(NOUT,*) ' NON-ZERO PIXELS IN MASK: ',NREPL
 	   RETURN
        ENDIF



        CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &   'ALTER (A)BOVE THRESHOLD, (B)ELOW, OR (C) BOTH SIDES (A/B/C)',
     &     NULL,IRT)
        IF (IRT .NE. 0) RETURN

	IF (ANS .EQ. 'C') THEN
           CALL RDPRM2S(TH1,TH2,NOT_USED,'UPPER, LOWER THRESHOLD',IRT)
           IF (IRT .NE. 0) RETURN

        ELSE
           CALL RDPRM1S(THR,NOT_USED,'THRESHOLD',IRT)
           IF (IRT .NE. 0) RETURN

	   IF (ANS .EQ. 'B')  THEN
	      TH2 = THR
	      TH1 = HUGE(TH1)
	   ELSE
	      TH1 = THR
	      TH2 = -HUGE(TH1)
	   ENDIF
        ENDIF

 	FIX1 = TH1
	FIX2 = TH2
	IF (FCHAR(4:4) .EQ. 'F')  THEN
	   CALL RDPRM1S(FIX1,NOT_USED,'FIXUP DENSITY',IRT)
           IF (IRT .NE. 0) RETURN
	   FIX2 = FIX1
	ENDIF

        NREPL = 0

	DO  J=1,NSLICE
	   DO  I=1,NROW
	      CALL REDLIN(LUN1,BUF,NSAM,(J-1)*NROW+I)
	      DO  K=1,NSAM
	         IF (BUF(K) .GE. TH1)  THEN
	            BUF(K) = FIX1
                    NREPL  = NREPL + 1

	         ELSEIF (BUF(K) .LE. TH2)  THEN
	            BUF(K) = FIX2
                    NREPL  = NREPL + 1
	         ENDIF
	      ENDDO
	      CALL WRTLIN(LUN2,BUF,NSAM,(J-1)*NROW+I)
	   ENDDO
	ENDDO

        FREPL = FLOAT(NREPL)
        CALL REG_SET_NSEL(1,1, FREPL,0.0, 0.0, 0.0, 0.0,IRTFLG)

        WRITE(NOUT,*) ' PIXELS REPLACED: ',NREPL
 
	END
