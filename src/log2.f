C++*******************************************************************
C
C $$ LOG2.FOR
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
C
C $$ LOG2: RETURNS THE SMALLEST POWER BY WHICH 2 HAS TO BE
C          RAISED TO OBTAIN AN INTEGER .GE.M
C 	   ONLY DEFINED FOR M.GE.1 !
C
C      FUNCTION LOG2(M)
C
C--*******************************************************************
C
C
        FUNCTION LOG2(M)
        LOG2 = 0
        I=1
        DO  K=1,100000
           IF(I.GE.M) RETURN
           I = I*2
           LOG2 = LOG2+1
        ENDDO
        END
