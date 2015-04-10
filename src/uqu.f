
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
C  UQU                                                                     *
C                                                                      *
C  PURPOSE:  sloppily written!                                                          *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE UQU(PIT)

        REAL     :: PIT(4)

        PARAMETER  (N=3)
        REAL     ::  P(N+1,N),Y(N+1),PR(N),PRR(N),PBAR(N),AK(N)

        COMMON     /ITERU/ ITER

        EXTERNAL   FCNQ

        DATA       PI/3.1415926/

        ITER = 0

        DO I=1,3
           P(2,I) = PIT(I) / 180.0 * PI
        ENDDO

C       BRACKET BY 20 DEGS
        BRACKET = 0.34906585
        P(1,1)  = P(2,1) - BRACKET
        P(1,2)  = P(2,2) + BRACKET
        P(1,3)  = P(2,3) + BRACKET
        P(3,1)  = P(2,1) + BRACKET
        P(3,2)  = P(2,2) + BRACKET
        P(3,3)  = P(2,3) - BRACKET
        P(4,1)  = P(2,1) - BRACKET
        P(4,2)  = P(2,2) + BRACKET
        P(4,3)  = P(2,3) - BRACKET

C       FIND INITIAL VALUES
        DO J=1,4
           DO    I=1,3
              AK(I) = P(J,I)
           ENDDO
           Y(J) = FCNQ(AK)
        ENDDO

        EPS    = 0.0001
        ITRMAX = 300      
        CALL AMOEBA(P,Y,N,EPS,FCNQ,ITER,ITRMAX,PR,PRR,PBAR)

        DO I=1,3
           PIT(I) = P(4,I) * 180.0/PI
        ENDDO

        PIT(4) = 1.0 - Y(4)

        END

