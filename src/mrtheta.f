
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
C
C  PURPOSE:                                                            *
C                                                                      *
C           CALCULATES THETA ANGLE DIFFERENCE BETWEEN
C           TWO VIEWS OF APPROXIMATELY THE SAME THETA
C           returns new theta, not the difference !
C
C PARAMETERS:
C     XYPTS(2,LS,LV) = ARRAY OF POINTS TO BE CHECKED
C     IVIEW = INDEX OF IMAGE TO BE CHECKED
C     P3D(3,LS) = 3-D MODEL
C     SCALE = DIFFERENCE IN SCALE BETWEEN IMAGES
C     PSI = IN-PLANE CORRECTION
C     THETA = TILT ANGLE OF IVIEW
C OUTPUT:
C     THETA = TILT DIFFERENCE
C
C***********************************************************************

      SUBROUTINE MRTHETA(RPT,VPT,IVIEW,P3D,THETA,PTACTIVE,NUMPTS,NTPT)

      PARAMETER (LV=300)
      PARAMETER (LS=256)

      LOGICAL     PTACTIVE(LS,LV)
      INTEGER     NUMPTS(LV)
      DIMENSION   RPT(2,LS),VPT(2,LS),P3D(3,LS)
      DIMENSION   CPJ(2,LS),CXY(2,LS)

      SY = SIN(THETA)
      CY = COS(THETA)

C     COMPUTE PROJECTION WITH X-Z POINTS AND VIEW WITH X
C     SO CAN GET COMMON CENTER OF MASS X.
      DO  I=1,NUMPTS(IVIEW)
	  IF (PTACTIVE(I,IVIEW))  THEN
             CPJ(1,I) = P3D(1,I) * SY + P3D(3,I) * CY
             CPJ(2,I) = RPT(1,I)
             CXY(2,I) = VPT(1,I)
	  ENDIF
      ENDDO

      DO  I=1,NUMPTS(IVIEW)
C         COMPUTE VIEW WITH X-Z POINTS
	  IF (PTACTIVE(I,IVIEW))  THEN
            TOTD2 = CPJ(1,I)**2 + CPJ(2,I)**2
            TMPZ2 = TOTD2 - (CXY(2,I)**2)
            TMPZ  = 0.0
            IF (TMPZ2 .GT. 0.0) TMPZ = SQRT(TMPZ2)
            CXY(1,I) = SIGN(TMPZ,CPJ(1,I))
	  ENDIF
      ENDDO

      CALL MRANG2(CPJ,CXY,IVIEW,THDIF,PTACTIVE,NTPT)
      THETA = THETA - THDIF

      END
