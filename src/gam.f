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

      FUNCTION GAM(X1)

      IMPLICIT  REAL*8  (A-H,O-Z)

      DATA A0/.9189385D0/,A1/.0005952D0/,A2/.0007937D0/
      DATA A3/.0027778D0/,A4/.0833333D0/

      X=X1
      IF(X.LT.7.0)GO TO 1
      F=0.0
      GO TO 2

    1 F=1.0
      Z=X-1.0

    3 Z=Z+1.0
      X=Z
      F=F*Z
      IF(Z.LT.7.0)GO TO 3
      X=X+1.0
      F=-DLOG(F)

    2 Z=1.0/X**2
      GAM=F+(X-0.5)*DLOG(X)-X+A0+(((-A1*Z+A2)*Z-A3)*Z+A4)/X
      RETURN
      END
