
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

	SUBROUTINE POJ(NW,QRT,PJX,PJY,XX,RJX,RJY,SIG2,N2,VV)

	DIMENSION QRT(NW,NW),PJX(NW),PJY(NW),XX(NW),RJX(N2)
	DIMENSION  RJY(N2),SIG2(NW)

C 	INITIALIZE

	DO  II=1,NW
	   XX(II)=II
	   PJX(II)=0.0
	   PJY(II)=0.0
	ENDDO

	DO  II=1,N2
	   RJX(II)=0.0
	   RJY(II)=0.0
	ENDDO

	DO  II=1,NW
           DO  JJ=1,NW
              T=QRT(JJ,II)
              PJX(II)=PJX(II)+T
              PJY(JJ)=PJY(JJ)+T
           ENDDO	
	ENDDO
	MWT=0
	CALL FITT(XX,PJX,NW,SIG2,MWT,A,B,SIGA,SIGB)
	DO  I=1,NW
           PJX(I)=PJX(I)-B*XX(I)-A
	ENDDO
	CALL FITT(XX,PJY,NW,SIG2,MWT,A,B,SIGA,SIGB)
	DO  I=1,NW
           PJY(I)=PJY(I)-B*XX(I)-A
	ENDDO

	K=0
	DO  I=1,NW,4
           K=K+1
           RJX(K)=PJX(I)+PJX(I+1)+PJX(I+2)+PJX(I+3)
           RJY(K)=PJY(I)+PJY(I+1)+PJY(I+2)+PJY(I+3)
	ENDDO

	CALL STEP(RJX,N2,INIJ)
	CALL STEP(RJY,N2,INIK)

	VV=INIJ*INIK
 	END
