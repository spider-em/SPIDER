
C ++********************************************************************
C                                                                      *
C  MR2TO3D                                                              *
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
C   MR2TO3D(P3D, PRJ, ANGLE,PTACTIVE,NUMPTS,NTVW,NTPT)                                                                    *
C                                                                      *
C  PURPOSE:                                                            *
C       CALCULATES, USING LINEAR REGRESSION AND ASSUMING A SINGLE
C       TILT AXIS THROUGH THE Y AXIS, THE X, Y, AND Z COORDS OF THE
C       MARKERS.
C
C  PARAMETERS:                                                         *
C  INPUT:
C     PRJ(2,LS,LV)   = COORDS OF MARKER IN PROJECTIONS
C     ANGLE(3,LV)    = ANGLES TILTED ABOUT Z,Y,Z
C     NUMPTS(LV)     = MAX INDEX OF POINT USED IN EACH VIEW
C     PTACTIVE(LS,LV)= BOOLEAN ARRAY HOLDING WHETHER POINT IS IN VIEW.
C     NTVW           = TOTAL NUMBER OF VIEWS
C     NTPT           = TOTAL NUMBER OF MARKERS
C  OUTPUT:
C     P3D(3,LS)      = COORDS OF POINTS IN 3-D SPACE                   
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE MR2TO3D(P3D, PRJ, ANGLE, PTACTIVE,NUMPTS,NTVW,NTPT)

      PARAMETER (LV=300)
      PARAMETER (LS=256)

      LOGICAL     PTACTIVE(LS,LV)
      INTEGER     NUMPTS(LV)

      DIMENSION   P3D(3,LS)
      DIMENSION   PRJ(2,LS,LV), ANGLE(3,LV)
      DIMENSION   XPR(LV), YPR(LV)

C     PREPARE ARRAYS FOR REGRESSION. 
      DO I=1,NTPT
        YSUM  = 0.0
	A2    = 0.0
	B2    = 0.0
	AB    = 0.0
	AP    = 0.0
	BP    = 0.0
        NYSUM = 0
        NDATA = 0

        DO  J=1,NTVW
          IF (PTACTIVE(I,J)) THEN
             YSUM  = YSUM + PRJ(2,I,J)
             NYSUM = NYSUM + 1
             A = COS(ANGLE(2,J))
             B = SIN(ANGLE(2,J))
             P = PRJ(1,I,J)

C            SOLVE SYSTEM OF EQUATIONS:
C            X*COS(THETA(I)) - Z*SIN(THETA(I)) = P(I)
C            THAT IS:
C            [ X*COS(THETA(I)) - Z*SIN(THETA(I)) - P(I)]^2 ->MIN

             A2 = A2 + A * A
             B2 = B2 + B * B
             AB = AB + A * B
             AP = AP + A * P
             BP = BP + B * P

             NDATA = NDATA + 1

           ENDIF
	ENDDO

	IF (NDATA.NE.0)  THEN
           P3D(1,I) = (AB * BP - B2 * AP) / (AB*AB-A2*B2)
           P3D(3,I) = (A2 * BP - AB * AP) / (AB*AB-A2*B2)
           P3D(2,I) = YSUM / NYSUM
	ELSE
	   DO J=1,3
	      P3D(J,I) = 0.0
	   ENDDO
	ENDIF
      ENDDO

      END
