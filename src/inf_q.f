
C ++********************************************************************
C                                                                      *
C INF_Q                                                               *
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
C  INF_Q                                                               *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE  INF_Q(X,LRCL,NSAM,NROW,Y,N)

      DOUBLE PRECISION  Y(29,29)
      DIMENSION  X(LRCL,NROW)

      DO    I=1,NSAM
      DO    J=1,NROW
        X(I,J)=0.0
      ENDDO
      ENDDO
      DO    I=1,N
      DO    J=1,N
        X(I,J)=Y(I,J)
      ENDDO
      ENDDO
      DO    I=2,N
	II=NSAM-I+2
	JJ=NROW-I+2
      X(1,JJ)=Y(1,I)
      X(II,1)=Y(I,1)
      ENDDO
      DO    I=2,N
      II=NSAM-I+2
      DO    J=2,N
      Z=Y(I,J)
      JJ=NROW-J+2
      X(II,JJ)=Z
      X(II,J)=Z
      X(I,JJ)=Z
      ENDDO
      ENDDO
	
      END
