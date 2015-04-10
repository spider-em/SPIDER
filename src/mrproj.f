
C ++********************************************************************
C                                                                      *
C  MRPROJ                                                                    *
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
C  MRPROJ(P3D, PRJ, ANGLE,NTVW,NTPT)                                                                    *
C                                                                      *
C  PURPOSE:  MAKES PSUEDO-PROJ FROM 3-D DATA TO COMPARE WITH REAL
C            PROJECTIONS.
C                                                                      *
C  PARAMETERS:                                                         *
C
C INPUT:
C     P3D(3,LS)   = POINTS IN 3 COORDS
C     ANGLE(3,LV) = EULER ANGLES OF IMAGES (2 BEING TILT)
C     NTPT        = NUMBER TOTAL MARKERS
C     NTVW        = NUMBER OF VIEWS
C OUTPUT:
C     PRJ(2,LS,LV)= COORDS OF MARKERS IN PROJECTION. Z DATA INCLUDED
C        FOR COMPLETENESS
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE MRPROJ(P3D, PRJ, ANGLE,NTVW,NTPT)

      PARAMETER (LV=300)
      PARAMETER (LS=256)

      DIMENSION   P3D(3,LS), PRJ(2,LS,LV), ANGLE(3,LV)

C     ROTATE ABOUT Y, Z2 (THETA, PSI)

      DO  IVIEW=1,NTVW
         CY = COS(ANGLE(2,IVIEW))
         SY = SIN(ANGLE(2,IVIEW))

         DO  I=1,NTPT
            PRJ(1,I,IVIEW) = P3D(1,I) * CY -  P3D(3,I) * SY
            PRJ(2,I,IVIEW) = P3D(2,I)
	 ENDDO
      ENDDO

      END
