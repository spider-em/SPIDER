C++*******************************************************************
C
C $$ SHIFTT.FOR
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
C $$ SHIFTT:    SHIFT A TWO-DIMENSIONAL ARRAY (PICTURE)
C		THAT IS PART OF A 3-D FILE
C
C       THIS SUBROUTINE WILL SHIFT A TWO-DIMENSIONAL ARRAY OR PICTURE
C       STORED ON LUN1 CIRCULARLY BY NSAMS,NROWS;  OUTPUT TO BE FOUND 
C       ON LUN2. RECORD NUMBER OFFSETS ARE USED TO GET INPUT FROM 
C	APPROPRIATE SLICE AND PUT IT INTO POSITION WITH Z-OFFSET
C
C       CALL SHIFTT(LUN1,LUN2,NSAM,NROW,NSAMS,NROWS,JROFF,JWOFF)
C         LUN1           LOGICAL UNIT NUMBER OF INPUT IMAGE
C         LUN2           LOGICAL UNIT NUMBER OF OUTPUT IMAGE
C         NSAM,NROW      DIMENSIONS OF IMAGE
C         NSAMS,NROWS    NUMBER OF SAMPLES AND ROWS TO BE SHIFTED
C	  JROFF,JWOFF	 RECORD NUMBER OFFSETS ON LUN1,LUN2 FOR READ
C	  AND WRITE, RESPECTIVELY
C
C--*******************************************************************
C
C
      SUBROUTINE SHIFTT(LUN1,LUN2,NSAM,NROW,NSAMS,NROWS,JROFF,JWOFF)
      COMMON BUF(1)

      NS = MOD(NSAMS,NSAM)
      NR = MOD(NROWS,NROW)
      NA = NSAM+1
      NA2 = NSAM+NSAM
      IF(NR.LT.0) NR = NR+NROW
      IF(NR.EQ.0) GOTO 300
      DO 200 I = 1,NROW
      CALL REDLIN(LUN1,BUF(NA),NSAM,I+JROFF)
      I1 = I+NR
      IF(I1.GT.NROW) I1=I1-NROW
      IF(NS)100,200,130
100   NS1 = -NS
      DO  K = 1,NS1
        BUF(NA2+K) = BUF(NA+K-1)
      ENDDO
      GOTO 200
130   NS1 = NS
      DO  K = 1,NS1
        BUF(NA-K) = BUF(NA2-K+1)
      ENDDO
200   CALL WRTLIN(LUN2,BUF(NA-NS),NSAM,I1+JWOFF)
      RETURN
300   DO 400 I = 1,NROW
      CALL REDLIN(LUN1,BUF(NA),NSAM,I+JROFF)
      IF(NS) 310,400,350
310   NS1 = -NS
      DO  K = 1,NS1
        BUF(NA2+K) = BUF(NA+K-1)
      ENDDO
      GOTO 400
350   NS1 = NS
      DO  K = 1,NS1
        BUF(NA-K) = BUF(NA2-K+1)
      ENDDO
400   CALL WRTLIN(LUN2,BUF(NA-NS),NSAM,I+JWOFF)
      END
