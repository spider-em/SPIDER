
C ++********************************************************************
C                                                                      *
C  MRALIGN               DUAL PARAMETER BUG     DEC 2008 ARDEAN LEITH  *                               *
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
C  MRALIGN(XYPTS,PRJ,ANGLE,SHIFT,SCALE,PTACTIVE,NTPT)                                                                   *
C                                                                      *
C  PURPOSE:                                                            *
C  TAKES IMAGES XYPTS AND TRANSFORMS THEM BY ROTATING ABOUT THE
C  GEOMETRICAL CENTER (NOT MASS CENTER) AND TRANSLATING IT BY SHIFT.   *
C
C PARAMETERS:                                                          *
C INPUT:
C     XYPTS(2,NTPT,LV)= COORDS OF POINTS TO BE TRANSFORMED
C     ANGLE(3,LV)= EULER ANGLES. ANGLE(2) IS TILT AND ANGLE(1) IS
C                  FINAL ROTATION CORRECTION (APPLIED HERE)
C     SHIFT(2,LV)= SHIFT (X,Y) THAT THE IMAGE MUST UNDERGO
C     SCALE(LV)  = MULTIPLICATIVE SCALING FACTOR
C     NTPT       = TOTAL NUMBER OF DISTINCT MARKERS USED
C
C OUTPUT:
C     XYPTS(2,NTPT,LV) = COORDS OF CORRECTED POINTS
C        
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE MRALIGN(XYPTS,ANGLE,SHIFT,SCALE,PTACTIVE,NTPT)

      LOGICAL     PTACTIVE(NTPT)
      DIMENSION   XYPTS(2,NTPT), PRJ(2,NTPT)
      DIMENSION   SHIFT(2)

      CA = COS(ANGLE)
      SA = SIN(ANGLE)

C     IMAGE CENTER (CIR) WAS ALREADY SUBTRACTED IN MRGETINFO.
C     ROTATION IS MEANT TO BE AROUND THE CENTER OF THE IMAGE
C     NOT AROUND THE CORNER OR CENTER OF MASS

C     X' = SCALE * R(PSI) * X + SHIFT

      DO  IPT=1,NTPT
        IF (PTACTIVE(IPT)) THEN

	  QT           = XYPTS(1,IPT)
          XYPTS(1,IPT) = SCALE * ( QT * CA + XYPTS(2,IPT)*SA) + SHIFT(1)
          XYPTS(2,IPT) = SCALE * (-QT * SA + XYPTS(2,IPT)*CA) + SHIFT(2)
	ENDIF
      ENDDO

      END
