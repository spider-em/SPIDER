
C ++********************************************************************
C                                                                      *
C  MRCG3D                                                                    *
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
C  MRCG3D(P3D,NTPT)                                                                     *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C  INPUT:
C     NTPT      = TOTAL NUMBER OF MARKERS
C     P3D(3,LS) = COORDS OF POINTS IN 3-D SPACE
C  OUTPUT:
C     P3D(3,LS) = CENTERED COORDS OF POINTS IN 3-D SPACE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE MRCG3D(P3D,NTPT)

        PARAMETER (LV=300)
        PARAMETER (LS=256)

        DIMENSION   P3D(3,LS)

C       GO OVER X-Y-Z
	DO J=1,3
	   CG = 0.0
           DO I=1,NTPT
              CG = CG + P3D(J,I)
           ENDDO
	   CG = CG / NTPT

C          SUBTRACT CENTER OF GRAVITY
           DO I=1,NTPT
	     	P3D(J,I) = P3D(J,I) - CG
           ENDDO
	ENDDO
      END
