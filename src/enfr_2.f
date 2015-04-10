C++*********************************************************************
C
C ENFR_2.F
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
C ENFR_2.F
C
C IMAGE_PROCESSING_ROUTINE                                                                     
C        0         2         3         4         5         6         7 
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        DOUBLE PRECISION FUNCTION ENFR_2(A,LSD,NSAM,NROW)

C       ENERGY IN FOURIER SPACE

        DIMENSION  A(LSD,NROW)
        DOUBLE PRECISION  ENFR

        ISBD=2*MOD(NSAM+1,2)
        IRBD=MOD(NROW,2)
        ENFR=0.0D0
        DO    J=1,NROW
           DO    I=3,LSD-ISBD
              ENFR=ENFR+A(I,J)*DBLE(A(I,J))
           ENDDO
        ENDDO
        DO    J=2,NROW/2+IRBD
           ENFR=ENFR+A(1,J)*DBLE(A(1,J))+A(2,J)*DBLE(A(2,J))
        ENDDO
        IF(MOD(NSAM,2).EQ.0) THEN
           DO    J=2,NROW/2+IRBD
        ENFR=ENFR+A(LSD-1,J)*DBLE(A(LSD-1,J))+A(LSD,J)*DBLE(A(LSD,J))
           ENDDO
        ENDIF
        ENFR=2*ENFR
        ENFR=ENFR+A(1,1)*DBLE(A(1,1))
        IF(MOD(NROW,2).EQ.0) ENFR=ENFR+A(1,NROW/2+1)*DBLE(A(1,NROW/2+1))
        IF(MOD(NSAM,2).EQ.0) THEN
           ENFR=ENFR+A(LSD-1,1)*DBLE(A(LSD-1,1))
            IF(MOD(NROW,2).EQ.0) ENFR=ENFR+
     &       A(LSD-1,NROW/2+1)*DBLE(A(LSD-1,NROW/2+1))
        ENDIF
        ENFR_2=ENFR/FLOAT(NSAM*NROW)
        END
