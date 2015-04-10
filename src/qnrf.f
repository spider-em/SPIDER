C ++********************************************************************
C                                                                      *
C  QNRF                                                                    *
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
C QNRF(Q1,Q2,KLX,KNX,KLY,KNY,KLZ,KNZ,R,AA,AB)                                *
C                                                                      *
C PURPOSE:                                                             *
C                                                                      *
C PARAMETERS: AA,AB                                         (RETURNED) *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE QNRF(Q1,Q2,KLX,KNX,KLY,KNY,KLZ,KNZ,R,AA,AB)

        REAL              :: Q1(KLX:KNX,KLY:KNY,KLZ:KNZ),
     &                       Q2(KLX:KNX,KLY:KNY,KLZ:KNZ)
        DOUBLE PRECISION  :: AA,AB,SX,SY,SXQ,SYQ

        RR  = R*R
        SX  = 0.0
        SY  = 0.0
        SXQ = 0.0
        SYQ = 0.0
        M   = 0

        DO IZ=KLZ,KNZ

           RZ = IZ * IZ

           DO IY=KLY,KNY

              RY = IY * IY + RZ

              DO IX=KLX,KNX
                 RT = IX * IX + RY

                 IF (RT .LE. RR) THEN
                    SXQ = SXQ + Q1(IX,IY,IZ) * DBLE(Q1(IX,IY,IZ))
                    SYQ = SYQ + Q2(IX,IY,IZ) * DBLE(Q2(IX,IY,IZ))
                    SX  = SX  + Q1(IX,IY,IZ)
                    SY  = SY  + Q2(IX,IY,IZ)
                    M   = M   + 1
                 ENDIF
              ENDDO
           ENDDO
        ENDDO

        AA = DSQRT((SXQ - SX * SX / M) * (SYQ - SY * SY / M))

        AB = SX * SY / M

        END

