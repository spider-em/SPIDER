
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
C     FIND THE SHIFT BETWEEN POINTS AND SHIFT XYPTS
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE MRSHIFT(PRJ,XYPTS,SHIFT,PTACTIVE,NTVW,NTPT)

      PARAMETER (LV=300)
      PARAMETER (LS=256)

      LOGICAL    PTACTIVE(LS,LV)
      DIMENSION  XYPTS(2,LS,LV),PRJ(2,LS,LV)
      DIMENSION  SHIFT(2,LS)


      DO IVIEW=1,NTVW
	DO J=1,2
	   SHIFT(J,IVIEW) = 0.0
	   NT = 0
	   DO I=1,NTPT
	      IF (PTACTIVE(I,IVIEW))  THEN
                 SHIFT(J,IVIEW) = SHIFT(J,IVIEW) - XYPTS(J,I,IVIEW)+
     &                           PRJ(J,I,IVIEW)
	         NT=NT+1
	      ENDIF
	   ENDDO

	   SHIFT(J,IVIEW) = SHIFT(J,IVIEW)/NT

           DO  I=1,NTPT
              IF (PTACTIVE(I,IVIEW))
     &		XYPTS(J,I,IVIEW) = XYPTS(J,I,IVIEW)+SHIFT(J,IVIEW)
	   ENDDO
         ENDDO
      ENDDO
      END
