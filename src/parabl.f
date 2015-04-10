
C++****************************************************** VAX 9/25/81
C
C PARABL.FOR
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
C PARABL  9/25/81 : PARABOLIC FIT TO 3 BY 3 PEAK NEIGHBORHOOD
C
C FORMULA FOR PARABOLOID TO BE FITTED INTO THE NINE POINTS:
C
C	F = C1 + C2*Y + C3*Y**2 + C4*X + C5*XY + C6*X**2
C
C THE VALUES OF THE COEFFICIENTS C1 - C6 ON THE BASIS OF THE
C NINE POINTS AROUND THE PEAK, AS EVALUATED BY ALTRAN:
C
C--*********************************************************************

	SUBROUTINE PARABL(Z,XSH,YSH,PEAKV)

	REAL Z(3,3)


	C1 = (26.*Z(1,1) - Z(1,2) + 2*Z(1,3) - Z(2,1) - 19.*Z(2,2)
     1        -7.*Z(2,3) + 2.*Z(3,1) - 7.*Z(3,2) + 14.*Z(3,3))/9.
C
	C2 = (8.* Z(1,1) - 8.*Z(1,2) + 5.*Z(2,1) - 8.*Z(2,2) + 3.*Z(2,3)
     1        +2.*Z(3,1) - 8.*Z(3,2) + 6.*Z(3,3))/(-6.)
C
	C3 = (Z(1,1) - 2.*Z(1,2) + Z(1,3) + Z(2,1) -2.*Z(2,2)
     1      + Z(2,3) + Z(3,1) - 2.*Z(3,2) + Z(3,3))/6.
C
	C4 = (8.*Z(1,1) + 5.*Z(1,2) + 2.*Z(1,3) -8.*Z(2,1) -8.*Z(2,2)
     1      - 8.*Z(2,3) + 3.*Z(3,2) + 6.*Z(3,3))/(-6.)
C
	C5 = (Z(1,1) - Z(1,3) - Z(3,1) + Z(3,3))/4.
C
	C6 = (Z(1,1) + Z(1,2) + Z(1,3) - 2.*Z(2,1) - 2.*Z(2,2)
     1    -2.*Z(2,3) + Z(3,1) + Z(3,2) + Z(3,3))/6.

C       THE PEAK COORDINATES OF THE PARABOLOID CAN NOW BE EVALUATED AS:

	YSH   = 0.0
	XSH   = 0.0
	DENOM = 4. * C3 * C6 - C5 * C5
	IF (DENOM .EQ. 0.) RETURN

	YSH   = (C4*C5 - 2.*C2*C6) / DENOM - 2.
	XSH   = (C2*C5 - 2.*C4*C3) / DENOM - 2.

	PEAKV = 4.*C1*C3*C6 - C1*C5*C5 - C2*C2*C6 + C2*C4*C5 - C4*C4*C3
	PEAKV = PEAKV / DENOM

C       LIMIT INTERPLATION TO +/- 1. RANGE
	IF (YSH .LT. -1.) YSH = -1.
	IF (YSH .GT. +1.) YSH = +1.
	IF (XSH .LT. -1.) XSH = -1.
	IF (XSH .GT. +1.) XSH = +1.

	END
