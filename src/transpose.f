
C
C **********************************************************************
C *  AUTHOR :                                                              *
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
C
C **********************************************************************
C

        SUBROUTINE TRANSPOSE(T,NPIX,NACT,NUMIM,LUV,LSAV,LTRANSP)

C        INTEGER*4 NUMIM, NPIX, NACT
        DIMENSION T(NPIX, NACT), LUV(NUMIM)
C, DAV, VAR, WEIGHT
C        INTEGER LUV(NUMIM)
C        INTEGER LTMP, LSAV, LAC, LTRANSP
        DATA  LTMP /78/, LAC /1/  
C        INTEGER JA, I, J

c  Warning !  luv changes meaning ...

        CALL  REW(LSAV, 1)
        JA = 0
        DO  I = 1, NUMIM
          IF(LUV(I) .EQ. 1)  THEN
            JA = JA + 1
            READ(LSAV)  (T(J, JA), J = 1, NPIX), dav, var, luv(i)
          ELSE
            READ(LSAV)
          END IF
        END DO
        
        WRITE(LTRANSP)  NPIX, NACT
        DO   J = 1, NPIX
          VAR = 0.0
          WEIGHT = 0.0
          DO  I = 1, NACT
            WEIGHT = WEIGHT + T(J, I)
            VAR = VAR + T(J, I)**2
          END DO
          DAV = WEIGHT/NACT
          VAR = (VAR - (DAV**2) * NACT )/(NACT - 1)
          WRITE(LTRANSP)  (T(J, I), I = 1, NACT), WEIGHT, VAR, LAC
        END DO
        RETURN
        END
