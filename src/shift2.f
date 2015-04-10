
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
C       THIS SUBROUTINE WILL SHIFT A TWO-DIMENSIONAL ARRAY OR PICTURE
C       STORED ON LUN1 CIRCULARLY BY NSAMS,NROWS;  OUTPUT TO BE FOUND 
C       ON LUN2.
C
C       CALL SHIFT2(LUN1,LUN2,BUF,NSAM,NROW,NSAMS,NROWS)
C         LUN1           LOGICAL UNIT NUMBER OF INPUT IMAGE
C         LUN2           LOGICAL UNIT NUMBER OF OUTPUT IMAGE
C         BUF            REAL BUFFER OF SIZE NSAM
C         NSAM,NROW      DIMENSIONS OF IMAGE
C         NSAMS,NROWS    NUMBER OF SAMPLES AND ROWS TO BE SHIFTED
C
C--*******************************************************************

      SUBROUTINE SHIFT2(LUN1,LUN2,NSAM,NROW,NSAMS,NROWS)

      COMMON BUF(1)

      NS = MOD(NSAMS,NSAM)
      NR = MOD(NROWS,NROW)
      IF(NR.LT.0) NR = NR+NROW
	if(ns.ne.0)  goto  300
      DO  I = 1,NROW
      CALL REDLIN(LUN1,BUF,NSAM,I)
      I1 = I+NR
      IF(I1.GT.NROW) I1=I1-NROW
      CALL WRTLIN(LUN2,BUF,NSAM,I1)
      ENDDO
      RETURN
300   DO  I = 1,NROW
      CALL REDLIN(LUN1,BUF,NSAM,I)
      CALL SHIFT1(BUF,buf(1+nsam),NSAM,NSAMS)
      I1 = I+NR
      IF(I1.GT.NROW) I1=I1-NROW
      CALL WRTLIN(LUN2,BUF(1+nsam),NSAM,I1)
      ENDDO
      END
