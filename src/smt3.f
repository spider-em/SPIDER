
C ++********************************************************************
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
C
C  SMT3(T,X,Y,NSAM,NSLICE,IPCUBE,NN)
C
C  PURPOSE:  
C                                                                  *
C***********************************************************************

	SUBROUTINE  SMT3(T,X,Y,NSAM,NSLICE,IPCUBE,NN)

	DIMENSION   X(NSAM,NSLICE,3),Y(NSAM,NSLICE)
	DIMENSION  IPCUBE(5,NN)

	DATA  LH/1/,K/2/
	D=T/26.0
	TM=1.0-T-D

cc$omp parallel do  private(kn,mk,mj,mi,l,i,j)
        DO KN=1,NN
           I=IPCUBE(3,KN)-1
           J=IPCUBE(5,KN)
           DO L=IPCUBE(1,KN),IPCUBE(2,KN)
              I=I+1
              IF(J.EQ.1 .OR. J.EQ.NSAM)  THEN
                 Y(I,J)=X(I,J,K)
              ELSE
                 QT=0.0
                 DO MK=-LH,LH
                    DO MJ=-LH,LH
                       DO MI=-LH,LH
                          QT=QT+X(I+MI,J+MJ,K+MK)
                       ENDDO
                    ENDDO
                 ENDDO
	         Y(I,J)=TM*X(I,J,K)+D*QT
              ENDIF
	   ENDDO
	ENDDO
	END	
