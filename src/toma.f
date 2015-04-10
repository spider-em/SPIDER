
C ++********************************************************************
C                                                                      *
C                                                                      *
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
C                                                                      *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE TOMA(NSAM,NROW,NNSAM,NNROW,Y,VART,MVAR)

 	DIMENSION Y(NSAM,NROW),VART(MVAR)

C FIRST WE FIND THE MINIMUM AND MAXIMUM FOR EACH IMAGE AND BEFORE EVALUATION
C OF THE ENTROPY PARAMETER EACH PIXEL VALUE IS CHANGED BY SUBTRACTING THE 
C MIN FROM IT AND DIVIDING IT BY THE DIFF OF MAX AND MIN.

	DO IK=1,MVAR
	   VART(IK)=0.0
	END DO

	AMIN=0.1E+30
	AMAX=0.1E-30
	DO  I=1,NROW
	DO  J=1,NSAM
	IF(Y(I,J).LT.AMIN)AMIN=Y(I,J)
	IF(Y(I,J).GT.AMAX)AMAX=Y(I,J)
	ENDDO
	ENDDO
C
	AAM=ABS(AMAX-AMIN)
	SUM=0.
	AD2=0.
	AD3=0.
	AD4=0.
	AD5=0.
	SUMTT=0.0
	DO  I=1,NROW
	   DO 400 J=1,NSAM
	      T=Y(J,I)
	      SUM=SUM+T
	      U=T*T	
	      AD2=AD2+U
 	      P=U*T
	      AD3=AD3+P
	      QT=T*P
	      AD4=AD4+QT
	      TT=(T-AMIN)/AAM
	      SUMTT=SUMTT+TT
	     IF(TT.LE.0.)GO TO 400
	      AD5=AD5+TT*ALOG10(TT)
400	   CONTINUE
	ENDDO
	X=FLOAT(NSAM)*FLOAT(NROW)
	AVG=SUM/X
	SQ=AVG*AVG
	CB=AVG*SQ
	QD=AVG*CB
C	VARI=AD2-(X*SQ)
	VART(1)=AD2-(X*SQ)
C	SKEW=AD3-3.*AVG*AD2+2.*X*CB
	VART(2)=AD3-3.*AVG*AD2+2.*X*CB
C	AKURT=AD4-4.*AVG*AD3-3.*QD*X+6.*SQ*AD2-3.
	VART(3)=AD4-4.*AVG*AD3-3.*QD*X+6.*SQ*AD2-3.
	IF(SUMTT.LE.0.) THEN
C	ENTP=0.0
	VART(4)=0.0
	GO TO 900
	ELSE
C	ENTP=ALOG10(SUMTT)-AD5/SUMTT
	VART(4)=ALOG10(SUMTT)-AD5/SUMTT
	ENDIF
900	CONTINUE
	CALL TIMA(NSAM,Y,AVG,NNSAM,NNROW,AVAV,AVVR,SDAV,SDVR)
C
	VART(5)=AVAV
	VART(6)=AVVR
	VART(7)=SDAV
	VART(8)=SDVR
        END
