C++*********************************************************************
C
C MASKMU.FOR
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
C  JMC, 10/30/86 CHANGED TO ACCEPT 3-D FILES
C
C  MASKMU:
C
C  PARAMETERS:     LUN1 LOGICAL UNIT NUMBER OF MASK
C                  LUN2 LOGICAL UNIT NUMBER OF IMAGE TO BE MASKED
C
C-**************************************************************************

	SUBROUTINE MASKMU(LUN1,LUN2,NSAM,NROW,NSLICE)

	COMMON   BUF(1)

	INCLUDE 'CMBLOCK.INC' 

        NTOTAL=NROW*NSLICE

	IF (FCHAR(4:4) .EQ. 'C') GOTO 150

	CALL RDPRM(BACK,NOT_USED,'BACKGROUND')

	DO  I =1,NTOTAL
      	   CALL REDLIN(LUN1,BUF,NSAM,I)
	   CALL REDLIN(LUN2,BUF(NSAM+1),NSAM,I)
	   DO  K=1,NSAM
	      IF (BUF(K) .LT. 0.5) BUF(NSAM+K)=BACK
           ENDDO
	   CALL WRTLIN(LUN2,BUF(NSAM+1),NSAM,I)
        ENDDO
	RETURN

C       OPTION FOR CONTINUOUS-VALUED MASKS ('MM C')

150	IF (IMAMI.NE.1) THEN
          CALL NORM3 (LUN2,NSAM,NROW,NSLICE,FMAX,FMIN,AV)      
        ENDIF

	DO I=1,NTOTAL
           CALL REDLIN(LUN1,BUF,NSAM,I)
	   CALL REDLIN(LUN2,BUF(NSAM+1),NSAM,I)
	   DO    K=1,NSAM
  	      BUF(NSAM+K)=(BUF(NSAM+K)-AV)*BUF(K) + AV
           ENDDO
           CALL WRTLIN(LUN2,BUF(NSAM+1),NSAM,I)
        ENDDO

	RETURN
	END
