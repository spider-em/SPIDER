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
C  PARAMETERS:
C                                                                      *
C IMAGE_PROCESSING_ROUTINE  
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************


	SUBROUTINE FINT(X,Y,NSAM,NROW,NSAM2,NROW2,LSD,LSD2)

	
	DIMENSION   X(LSD,NROW), Y(LSD2,NROW2)

C       TO KEEP THE EXACT VALUES ON THE PREVIOUS GRID YOU SHOULD USE
C        SQ2     = 2.0. HOWEVER, TOTAL ENERGY WILL NOT BE CONSERVED.

	SQ2    = SQRT(2.0)
	ANORM  = FLOAT(NSAM2*NROW2)/FLOAT(NSAM*NROW)
	NM     = MOD(NSAM,2)

	INV = +1
	CALL  FMRS_2(X,NSAM,NROW,INV)

C       NORMALIZATION REQUIRED
        X = X*ANORM

	IF(NSAM<=NSAM2)  THEN
        INSAM  = NSAM2 - NSAM
        INROW  = NROW2 - NROW
C       ZERO PADDING IS DONE

	Y = 0.0

	DO J = 1,NROW/2+1
	   DO I =1,LSD
	      Y(I,J) = X(I,J) 
	   ENDDO
	ENDDO

	DO J = NROW/2+2+INROW, NROW2
	   DO I =1,LSD
	      Y(I,J) = X(I,J-INROW) 
	   ENDDO
	ENDDO

C       WEIGHTING FACTOR USED FOR EVEN NSAM. REQUIRED SINCE ADDING ZERO FOR
C       INTERPOLATION WILL INTRODUCE A COMPLEX CONJUGATE FOR NSAM/2 TH
C       ELEMENT.

	IF (NM .EQ. 0 .AND. INSAM .NE. 0) THEN
	   DO J = 1,NROW2
	      Y(NSAM+1,J) = (1/SQ2)*Y(NSAM+1,J)
	      Y(NSAM+2,J) = (1/SQ2)*Y(NSAM+2,J)
	   ENDDO
       	   DO I =1,LSD
	      Y(I,NROW/2+1+INROW) = Y(I,NROW/2+1)/SQ2
	      Y(I,NROW/2+1)       = Y(I,NROW/2+1)/SQ2
	   ENDDO
	ENDIF
	ELSE
	Y(:,1:NROW2/2+1)=X(1:LSD2,1:NROW2/2+1)
	DO  J=NROW2,NROW2/2+2,-1
	Y(:,J)=X(1:LSD2,J+NROW-NROW2)
	ENDDO
	ENDIF
	INV = -1
	CALL  FMRS_2(Y,NSAM2,NROW2,INV)


	END
