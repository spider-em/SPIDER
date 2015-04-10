
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
C  MRROTATE(ROT, PHI, THETA, PSI)                                      *
C                                                                      *
C  PURPOSE:                                                            *

C   TAKES IN ROTATION MATRIX AND RETURNS THE ANGLES USED TO
C    UNDO THAT MATRIX, USING EULER ANGLES.  ROTATION DONE
C    ABOUT Z(PHI), THEN NEW Y(THETA), THEN NEW Z(PSI).
C    OUTPUTS ANGLES TO CHANGE, IN ORDER, TO MATCH VPT TO RPT.
C
C INPUT:
C     ROT(3,3) = ROTATION MATRIX
C OUTPUT
C     PHI, THETA, PSI = ANGLES IN RADIANS
C
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE MRROTATE(ROT, PHI, THETA, PSI)

        INCLUDE 'CMBLOCK.INC'

	DIMENSION         ROT(3,3)
	DOUBLE PRECISION  R3(3,3)
	PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
	PARAMETER (DEPS = 1.0E-7)

	DO    I=1,3
	   DO   J=1,3
	      R3(I,J)=ROT(I,J)
	   ENDDO
	ENDDO

C       LIMIT PRECISION
	DO    J=1,3
	   DO    I=1,3
	      IF(DABS(R3(I,J)).LT.DEPS)  R3(I,J)=0.0D0
	      IF(R3(I,J)-1.0D0.GT.-DEPS)  R3(I,J)=1.0D0
	      IF(R3(I,J)+1.0D0.LT.DEPS)  R3(I,J)=-1.0D0
	   ENDDO	
	ENDDO

	IF(R3(3,3).EQ.1.0)  THEN
	   THETA=0.0
	   PSI=0.0
	   IF(R3(1,1).EQ.0.0)  THEN
	      PHI=DASIN(R3(1,2))
	   ELSE
	      PHI=DATAN2(R3(1,2),R3(1,1))
	   ENDIF
	ELSEIF(R3(3,3).EQ.-1.0)  THEN
	   THETA=QUADPI
	   PSI=0.0
	   IF(R3(1,1).EQ.0.0)  THEN
	      PHI=DASIN(-R3(1,2))
	   ELSE
	      PHI=DATAN2(-R3(1,2),-R3(1,1))
	   ENDIF
	ELSE
	   THETA=DACOS(R3(3,3))
	   ST=DSIGN(1.0D0,DBLE(THETA))
           IF(R3(3,1).EQ.0.0)  THEN
               IF(ST.NE.DSIGN(1.0D0,R3(3,2)))  THEN
                  PHI=1.5*QUADPI
               ELSE
                  PHI=0.5*QUADPI
               ENDIF
            ELSE
               PHI=DATAN2(R3(3,2)*ST,R3(3,1)*ST)
            ENDIF
            IF(R3(1,3).EQ.0.0)  THEN
               IF(ST.NE.DSIGN(1.0D0,R3(2,3)))  THEN
                  PSI=1.5*QUADPI
               ELSE
                  PSI=0.5*QUADPI
               ENDIF
            ELSE
               PSI=DATAN2(R3(2,3)*ST,-R3(1,3)*ST)
            ENDIF
         ENDIF
         IF(PSI.LT.0.0)  PSI=PSI+2*QUADPI
         IF(THETA.LT.0.0)  THETA=THETA+2*QUADPI
         IF(PHI.LT.0.0)  PHI=PHI+2*QUADPI
      END
