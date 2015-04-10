C++*********************************************************************
C
C $$ CALD.FOR
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
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************
C
C $$ CALD.FOR
C
         SUBROUTINE  CALD(FI1,FI2,FIO)
         DIMENSION  FI1(3),FI2(3),FIO(3)
         DOUBLE PRECISION  R1(3,3),R2(3,3),R3(3,3)
         DOUBLE PRECISION  PSI,THETA,PHI,DEPS
	DOUBLE PRECISION  QUADPI,DGR_TO_RAD,RAD_TO_DGR
	PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
	PARAMETER (DGR_TO_RAD = (QUADPI/180))
	PARAMETER (RAD_TO_DGR = (180.0/QUADPI))
         DATA  DEPS/1.0D-7/
C
         CALL  BLDR(R2,FI2(1),FI2(2),FI2(3))
         CALL  BLDR(R1,-FI1(3),-FI1(2),-FI1(1))
         DO    I=1,3
         DO    J=1,3
         R3(I,J)=0.0
         DO    K=1,3
         R3(I,J)=R3(I,J)+R2(I,K)*R1(K,J)
	 ENDDO
 	 ENDDO
	 ENDDO
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
C  Limit precision
         DO    J=1,3
         DO    I=1,3
         IF(DABS(R3(I,J)).LT.DEPS)  R3(I,J)=0.0D0
         IF(R3(I,J)-1.0D0.GT.-DEPS)  R3(I,J)=1.0D0
         IF(R3(I,J)+1.0D0.LT.DEPS)  R3(I,J)=-1.0D0
	 ENDDO
	 ENDDO
C
         IF(R3(3,3).EQ.1.0)  THEN
         THETA=0.0
         PSI=0.0
         IF(R3(1,1).EQ.0.0)  THEN
         PHI=RAD_TO_DGR*DASIN(R3(1,2))
         ELSE
         PHI=RAD_TO_DGR*DATAN2(R3(1,2),R3(1,1))
         ENDIF
         ELSEIF(R3(3,3).EQ.-1.0)  THEN
         THETA=180.0
         PSI=0.0
         IF(R3(1,1).EQ.0.0)  THEN
         PHI=RAD_TO_DGR*DASIN(-R3(1,2))
         ELSE
         PHI=RAD_TO_DGR*DATAN2(-R3(1,2),-R3(1,1))
         ENDIF
         ELSE
         THETA=RAD_TO_DGR*DACOS(R3(3,3))
         ST=DSIGN(1.0D0,THETA)
         IF(R3(3,1).EQ.0.0)  THEN
         IF(ST.NE.DSIGN(1.0D0,R3(3,2)))  THEN
         PHI=270.0
         ELSE
         PHI=90.0
         ENDIF
         ELSE
         PHI=RAD_TO_DGR*DATAN2(R3(3,2)*ST,R3(3,1)*ST)
         ENDIF
         IF(R3(1,3).EQ.0.0)  THEN
         IF(ST.NE.DSIGN(1.0D0,R3(2,3)))  THEN
         PSI=270.0
         ELSE
         PSI=90.0
         ENDIF
         ELSE
         PSI=RAD_TO_DGR*DATAN2(R3(2,3)*ST,-R3(1,3)*ST)
         ENDIF
         ENDIF
         IF(PSI.LT.0.0)  PSI=PSI+360.0
         IF(THETA.LT.0.0)  THETA=THETA+360.0
         IF(PHI.LT.0.0)  PHI=PHI+360.0
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         FIO(1)=PSI
         FIO(2)=THETA
         FIO(3)=PHI
         END
