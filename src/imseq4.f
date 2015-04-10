
C ++********************************************************************
C                                                                      *
C                                                                      *
C       REMOVED FROM JPMSK2.FOR   DEC 90 AL
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
C  PURPOSE:
C                                                            *
C THIS ROUTINE CHECKS IF NUMBER IMNO IS CONTAINED IN ARRAY IDI.
C IF YES, IT RETURNS THE ADDRESS IN IDI WHERE IMNO HAS BEEN LOCATED.
C IF NO, IT RETURNS 0.
C
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	INTEGER   FUNCTION IMSEQ4(NUMIM,IMNO,IDI)

	INTEGER   IMNO, NUMIM, IDI(NUMIM)

	DO  I=1,NUMIM
          IF  (IMNO .EQ. IDI(I)) THEN
	  IMSEQ4 = I
	  RETURN
	  ENDIF
	ENDDO
	IMSEQ4 = 0
	END

