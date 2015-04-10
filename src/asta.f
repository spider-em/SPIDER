
C++*********************************************************************
C
C  ASTA.F
C              ASTA_D                              JAN 13 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2013  Health Research Inc.,                         *
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
C  ASTA  (X,N,RI,ABA,KLP)
C  ASTA_D(X,N,RI,ABA,DKLP_8)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE ASTA(X,N,RI,ABA,KLP)

        REAL              :: X(N,N)
        DOUBLE PRECISION  :: ABA

C       ESTIMATE AVERAGE OUTSIDE THE CIRCLE
        R  = RI * RI
        NC = N / 2 + 1

        DO   J=1,N
           T  = J-NC
           XX = T*T

           DO   I=1,N
              T = I - NC
              IF (XX+T*T .GT. R)    THEN
                 ABA = ABA + DBLE(X(I,J))
                 KLP = KLP + 1
              ENDIF
           ENDDO
        ENDDO

        END



        SUBROUTINE ASTA_D(X,N,RI,ABA,DKLP)

        IMPLICIT NONE

        REAL              :: X(N,N)
        INTEGER           :: N
        REAL              :: RI
        DOUBLE PRECISION  :: ABA, DKLP

        REAL              :: R,T,XX
        INTEGER           :: NC,J,I

C       ESTIMATE AVERAGE OUTSIDE THE CIRCLE

        R  = RI * RI
        NC = N / 2 + 1

        DO   J=1,N
           T  = J - NC
           XX = T * T
           DO   I=1,N
              T = I - NC
              IF (XX+T*T > R) THEN
                 ABA    = ABA + DBLE(X(I,J))
                 DKLP = DKLP + 1
              ENDIF
           ENDDO
        ENDDO

        END


