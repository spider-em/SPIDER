C++*********************************************************************
C
C $$ CENT_Q.FOR
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
C IMAGE_PROCESSING_ROUTINE
C
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************
C
C $$ CENT_Q.FOR
C

        SUBROUTINE  CENT_Q(X,LSD,NSAM,NROW,SNS,SNR)
        DIMENSION  X(LSD,NROW)

        C=0.0
        S=0.0
        P=8*DATAN(1.0D0)/NSAM
c$omp parallel do private(i,t,j),reduction(+:c,s)
        DO    I=1,NSAM
           T=0.0
           DO    J=1,NROW
              T=T+X(I,J)
           ENDDO
           C=C+COS(P*(I-1))*T
           S=S+SIN(P*(I-1))*T
        ENDDO

        FI=ATAN2(S,C)
        IF(FI.LT.0.0)  FI=FI+8*DATAN(1.0D0)
        SNS=FI/P+1.0

        C=0.0
        S=0.0
        P=8*DATAN(1.0D0)/NROW
c$omp parallel do private(i,t,j),reduction(+:c,s)
        DO    J=1,NROW
           T=0.0
           DO    I=1,NSAM
              T=T+X(I,J)
           ENDDO
           C=C+COS(P*(J-1))*T
           S=S+SIN(P*(J-1))*T
        ENDDO

        FI=ATAN2(S,C)
        IF(FI.LT.0.0)  FI=FI+8*DATAN(1.0D0)
        SNR=FI/P+1.0
        END
