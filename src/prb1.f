C ++********************************************************************
C
C   PRB1.F     FROM PRBD1                          OCT 10 ARDEAN LEITH
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
C PURPOSE:  INTERPOLATION OVER NPOINT NEIGHBORHOOD, RETURNS OFFSET
C           OF MAX FROM MID-POINT
C 
C PARAMETERS:  B        NEIGHBORHOOD VALUES                   SENT
C              NPOINT   LENGTH OF NEIGHBORHOOD (ODD)          SENT
C              POS      OFFSET OF INTERPOLATED MAX            RET.    
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE PRB1(B,NPOINT,POS)

        INTEGER :: NPOINT
        REAL    :: B(NPOINT)
        REAL    :: POS

        REAL    :: C2,C3
        INTEGER :: NHALF

        NHALF  = NPOINT/2 + 1
        POS    = 0.0

        IF (NPOINT .EQ. 7) THEN
           C2 = 49.*B(1) + 6.*B(2) - 21.*B(3) - 32.*B(4) - 27.*B(5)
     1       - 6.*B(6) + 31.*B(7)
           C3 = 5.*B(1) - 3.*B(3) - 4.*B(4) - 3.*B(5) + 5.*B(7)
        ELSEIF (NPOINT .EQ. 5) THEN
           C2 = (74.*B(1) - 23.*B(2) - 60.*B(3) - 37.*B(4)
     1        + 46.*B(5) ) / (-70.)
           C3 = (2.*B(1) - B(2) - 2.*B(3) - B(4) + 2.*B(5) ) / 14.
        ELSEIF (NPOINT .EQ. 3) THEN
           C2 = (5.*B(1) - 8.*B(2) + 3.*B(3) ) / (-2.)
           C3 = (B(1) - 2.*B(2) + B(3) ) / 2.
        ELSEIF (NPOINT .EQ. 9) THEN
           C2 = (1708.*B(1) + 581.*B(2) - 246.*B(3) - 773.*B(4)
     1         - 1000.*B(5) - 927.*B(6) - 554.*B(7) + 119.*B(8)
     1         + 1092.*B(9) ) / (-4620.)
           C3 = (28.*B(1) + 7.*B(2) - 8.*B(3) - 17.*B(4) - 20.*B(5)
     1        - 17.*B(6) - 8.*B(7) + 7.*B(8) + 28.*B(9) ) / 924.
        ENDIF

        IF (C3 .NE. 0.0)  POS = C2 / (2.*C3) - NHALF

        END




