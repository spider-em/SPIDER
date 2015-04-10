C++*********************************************************************
C
C   NEGATE.F 
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
C   NEGATE(LUN1,LUN2,NSAM,NROW,NSLICE,FMAXT)
C
C   PURPOSE:  CARRIES OUT NEGATE OPERATION ON IMAGE PIXEL BY PIXEL, THEN
C             ADDS ORIGINAL FMAX TO EACH PIXEL.
C
C   PARAMETERS:
C        LUN1         LOGICAL UNIT NUMBER OF FILE 1
C        LUN2         LOGICAL UNIT NUMBER OF FILE 2
C        NSAM,NROW    X & Y DIMENSIONS OF FILES
C        NSLICE       Z DIMENSION OF FILES
C        FMAX         ORIGINAL IMAGE FMAX
C
C--*******************************************************************

      SUBROUTINE NEGATE(LUN1,LUN2,NSAM,NROW,NSLICE,FMAX)
      DIMENSION  BUF(NSAM)
      DO  IROW = 1,NROW*NSLICE
         CALL REDLIN(LUN1,BUF,NSAM,IROW)
             BUF = FMAX - BUF
         CALL WRTLIN(LUN2,BUF,NSAM,IROW)
      ENDDO

      END
