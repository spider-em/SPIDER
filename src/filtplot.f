C++*********************************************************************
C
C FILTPLOT.F                         
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
C  FILTPLOT(MAXMEM)
C
C  PARAMETERS:     MAXMEM      MAXIMUM MEMORY IN UNLABELED COMMON
C
C--*******************************************************************


	SUBROUTINE FILTPLOT(MAXMEM)

        INCLUDE 'CMBLOCK.INC'
        CHARACTER*1 ANS,NULL

        NULL=CHAR(0)

        CALL RDPRMC(ANS,NC,.TRUE.,
     &  '(F)ERMI, (G)AUSS, (B)UTER, (R)EMEZ FILTER? (F/G/B/R)',NULL,IRT)

        IF (ANS .EQ. 'F') THEN
           CALL FERMP

        ELSEIF (ANS .EQ. 'G') THEN
           CALL GAUSSP

	ELSEIF (ANS.EQ.'B')THEN
	   CALL BUTERP

        ELSEIF (ANS .EQ. 'R')  THEN
           CALL  REMEZP(MAXMEM)
        ENDIF

        END
