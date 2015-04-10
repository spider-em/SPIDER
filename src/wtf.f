
C ++********************************************************************
C                                                                      *
C  WTF.F              ADDED ISELECT PARAMETER              6/03/03 al                                                               *
C                     IF (ABS(Y) .LT.                      6/21/06 pp                                            *
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
C  WTF(ISELECT,PROJ,W,NNNN,NSAM,NROW,SS,NANG,SNR,K)                            *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:   ISELECT     SELECTED PROJECTIONS LIST           (SENT)
C  PARAMETERS:   USESELECT   USE SELECTION LIST                  (SENT)
C                PROJ                                       (SENT/RET.)
C                W                                               (WORK)
C                NNNN                                            (SENT)
C                NSAM,NROW                                       (SENT)
C                SS                                              (SENT)
C                NANG                                            (SENT)
C                SNR                                             (SENT)
C                K         CURRENT PROJ/ANGLE NUMBER             (SENT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE WTF(ISELECT,USESELECT,PROJ,W,NNNN,NSAM,NROW,SS,
     &                 NANG,SNR,K)

        REAL    :: SS(6,NANG),PROJ(NNNN,NROW),W(NNNN/2,NROW)
        INTEGER :: ISELECT(NANG)
        LOGICAL :: USESELECT

c$omp   parallel do private(i,j)
	DO J=1,NROW
	   DO I=1,NNNN/2
	      W(I,J) = 0.0
	   ENDDO
	ENDDO

	NR2 = NROW/2

	DO LT=1,NANG

           L = LT
           IF (USESELECT) L = ISELECT(LT)

           OX = SS(6,K)*SS(4,L)*(-SS(1,L)*SS(2,K)+ SS(1,K)*SS(2,L)) + 
     &      SS(5,K)*(-SS(3,L)*SS(4,K)+SS(3,K)*SS(4,L)*(SS(1,K)*SS(1,L) 
     &      + SS(2,K)*SS(2,L)))

           OY = SS(5,K)*SS(4,L)*(-SS(1,L)*SS(2,K)+ SS(1,K)*SS(2,L)) - 
     &       SS(6,K)*(-SS(3,L)*SS(4,K)+SS(3,K)*SS(4,L)*(SS(1,K)*SS(1,L) 
     &      + SS(2,K)*SS(2,L)))

	   IF (OX .NE. 0.0 .OR. OY.NE.0.0) THEN
C             THIS TEST MORE OFTEN TRUE SO PUT IT FIRST FOR SPEED
c$omp         parallel do private(i,j,jy,y,qt)
              DO J=1,NROW
                 JY = (J-1)
                 IF (JY .GT. NR2) JY=JY-NROW
                 DO I=1,NNNN/2
                    Y = OX * (I-1) + OY * JY

C                   CAN NEGLECT SMALL QT's 
c                   IF (Y .LT. 2.3) THEN EXP(-4*Y*Y) = 6.8e-10
c                   IF (Y .LT. 1.8) THEN EXP(-4*Y*Y) = 2.4e-6
c                   IF (Y .LT. 1.7) THEN EXP(-4*Y*Y) = 9.5e-6
c                   IF (Y .LT. 1.6) THEN EXP(-4*Y*Y) = 3.7e-5
c                   IF (Y .LT. 1.5) THEN EXP(-4*Y*Y) = 1.0e-4

                    IF (ABS(Y) .LT. 1.6) THEN
C                      SYSTEM FAULTS ON VERY SMALL QT's ON ALTIX
                       QT     = EXP(-4*Y*Y)
                       W(I,J) = W(I,J) + QT
                    ENDIF
                 ENDDO
	      ENDDO
	   ELSE
c$omp         parallel do private(i,j)
	      DO J=1,NROW
	         DO I=1,NNNN/2
                    W(I,J) = W(I,J) + 1.0
                 ENDDO
              ENDDO
	   ENDIF
	ENDDO

  	INV = +1
	CALL FMRS_2(PROJ,NSAM,NROW,INV)
        IF (INV .EQ. 0)THEN
           CALL ERRT(38,'WTF',NE)
           RETURN
        ENDIF

	WNRM = W(1,1)
c$omp   parallel do private(i,j,kx,ww)
	DO  J=1,NROW
	   DO  I=1,NNNN,2
	      KX          = (I+1)/2
	      WW          = W(KX,J)/WNRM/((W(KX,J)/WNRM)**2+SNR)
	      PROJ(I,J)   = PROJ(I,J)*WW
	      PROJ(I+1,J) = PROJ(I+1,J)*WW
	   ENDDO
	ENDDO

	INV = -1
	CALL FMRS_2(PROJ,NSAM,NROW,INV)

	END
