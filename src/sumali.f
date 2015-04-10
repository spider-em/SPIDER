
C++***************************************************************1/5/82
C
C  SUMALI.F
C                           REMOVED INPUT REGISTERS AUG 01 ARDEAN LEITH
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
C  SUMALI.F    
C
C  PURPOSE:  TO 'SUM' ANGLES,SHIFTS FROM TWO ALIGNMENT CYCLES
C
C--*********************************************************************

	SUBROUTINE SUMALI(P)

	INCLUDE 'CMBLOCK.INC'

	PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
	PARAMETER (DGR_TO_RAD = (QUADPI/180))
	
	LOGICAL       P

        CALL REG_GET_USED(NSEL_USED)

        IF (NSEL_USED .GT. 3) THEN
C         DEPRECATED REGISTER INPUT USED (AUG 01)
          CALL REG_GET_NSEL(1,ANGOLD,XSHOLD,YSHOLD,ANGNEW,XSHNEW,IRTFLG)
          CALL REG_GET_NSEL(6,YSHNEW,FDUM,FDUM,FDUM,FDUM,IRTFLG)
          IREGO = 7
        ELSE
C         PREFERRED INPUT METHOD
          CALL RDPRM3S(ANGOLD,XSHOLD,YSHOLD,NOT_USED,
     &                 'INITIAL ROTATION ANGLE, X & Y SHIFTS',IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

          CALL RDPRM3S(ANGNEW,XSHNEW,YSHNEW,NOT_USED,
     &                 'REFINED ROTATION ANGLE, X & Y SHIFTS',IRTFLG)
          IF (IRTFLG .NE. 0) RETURN
          IREGO = 1
        ENDIF

	C      = COS(ANGNEW*DGR_TO_RAD)
	S      = SIN(ANGNEW*DGR_TO_RAD)
	IF (P)  S = -S
	XSHSUM = XSHNEW+XSHOLD*C-YSHOLD*S
	YSHSUM = YSHNEW+XSHOLD*S+YSHOLD*C
	ANGSUM = ANGNEW+ANGOLD
	DO WHILE(ANGSUM.LT.0.0)
	   ANGSUM=ANGSUM+360.0
	ENDDO
	ANGSUM=AMOD(ANGSUM,360.0)

        CALL REG_SET_NSEL(IREGO,3,ANGSUM,XSHSUM,YSHSUM,0.0 ,0.0 ,IRTFLG)

	END
