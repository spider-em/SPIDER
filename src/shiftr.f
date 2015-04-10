
C++*******************************************************************
C SHIFTR.FOR
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
C   SHIFTR:    SHIFT A TWO-DIMENSIONAL  PICTURE  9/26/81 NON-INTEGER SHIFT
C
C       THIS SUBROUTINE WILL SHIFT A TWO-DIMENSIONAL ARRAY OR PICTURE
C       STORED ON LUN1 CIRCULARLY BY SAMS,ROWS;  OUTPUT TO BE FOUND 
C       ON LUN2.
C
C       CALL SHIFTR(LUN1,LUN2,BUF,NSAM,NROW,NNROWS,NNROWE,NNROWK,SAMS,ROWS)
C         LUN1           LOGICAL UNIT NUMBER OF INPUT IMAGE
C         LUN2           LOGICAL UNIT NUMBER OF OUTPUT IMAGE
C         BUF            REAL BUFFER OF SIZE NSAM
C         NSAM,NROW      DIMENSIONS OF IMAGE
C         NNROWS,NNROWE
C         NNROWK         STARTING ROW, ENDING, AND SKIPING FACTOR
C         SAMS,ROWS      SHIFT VECTOR COMPONENTS IN SAMPLE AND ROW DIR.
C
C  JMC, 1/1/87, THIS SUBROUTINE HAS BEEN CHANGED IN ORDER TO ACCEPT
C               A GENERAL VALUE FOR THE STARTING AND ENDING ROWS, AS
C               WELL AS A SKIPPING FACTOR. THIS CHANGES ARE DONE
C               FOR ITS USE IN THE 3-D CASE   
C
C--*******************************************************************

      SUBROUTINE SHIFTR(LUN1,LUN2,NSAM,NROW,NNROWS,NNROWE,NNROWK,
     1                  SAMS,ROWS)
      COMMON BUF(1)

	NSAMS=SAMS
	NROWS=ROWS
	DX=SAMS-NSAMS
	IF(SAMS.LT.0.)DX= 1+DX
	DY=ROWS-NROWS
	IF(ROWS.LT.0.)DY= 1+DY

	C1 = (1-DX) * (1-DY)
	C2 = DX * (1-DY)
	C3 = DY * (1-DX)
	C4 = DX * DY
C	WRITE(5,999)DX,DY,C1,C2,C3,C4
C999	FORMAT(6F10.4)
 
        NS = MOD(NSAMS,NSAM)+1
        NR = MOD(NROWS,NROW)+1
	IF(SAMS.LT.0.)NS=NS-1
	IF(ROWS.LT.0.)NR=NR-1

C ADDRESSES USED
        NA = NSAM+1
C NA = INPUT BUFFER I OFFSET +1
	NA3 = 3*NSAM+1
C NA3= INPUT BUFFER II OFFSET +1
	NA5 = 5*NSAM
C NA5=OUTPUT BUFFER OFFSET
	NA1= 2*NSAM+1
        NA2 = 3*NSAM
        IF(NR.LT.0) NR = NR+NROW
C NA6=OUTPUT BUFFER FOR A THREE-DIMENSIONAL SHIFT
        NA6=6*NSAM 

C INITIALIZE FIRST BUFFER LINE
        CALL REDLIN(LUN1,BUF(NA),NSAM,NNROWE)
	CALL REDLIN(LUN1,BUF(1),NSAM,NNROWE)
	CALL REDLIN(LUN1,BUF(NA+NSAM),NSAM,NNROWE)

	I=0
        NNNROW=NNROWS-NNROWK
80	I=I+1
        NNNROW=NNNROW+NNROWK
	NOFF1=NA
	NOFF2=NA3
	IF(MOD(I,2).EQ.0) NOFF2=NA
	IF(MOD(I,2).EQ.0) NOFF1=NA3
	NO1 = NOFF1-NS
	NO2 = NOFF2-NS

C NOFF1  POSITION OF FIRST ELEMENT OF OLD LINE
C NOFF2  POSITION OF FIRST ELEMENT OF NEW LINE
	CALL REDLIN(LUN1,BUF(NOFF2),NSAM,NNNROW)
C FOR NS=0, TWO NEIGHBORS ON EITHER END OF BUFFER ARE NEEDED
	BUF(NOFF2+NSAM)=BUF(NOFF2)
	BUF(NOFF2-1)=BUF(NOFF2+NSAM-1)

        I1 = I+NR-1
        IF(I1.GT.NROW) I1=I1-NROW
        IF(I1.LT.1)I1=I1+NROW
        NNROI1=NNROWS+(I1-1)*NNROWK
        IF(NS)100,145,130
100     NS1 = -NS
        DO  K = 1,NS1+1
	BUF(NOFF2+NSAM+K-1) = BUF(NOFF2+K-1)
	ENDDO
        GOTO 145

130     NS1 = NS
        DO  K = 1,NS1+1
        BUF(NOFF2-K) = BUF(NOFF2+NSAM-K)
	ENDDO
145	DO  K=1,NSAM
        BUF(NA5+K) = BUF(NO2+K-1)*C2 + BUF(NO2+K)*C1 + BUF(NO1+K-1)*C4
     1               +  BUF(NO1+K)*C3
        ENDDO
        CALL WRTLIN(LUN2,BUF(NA5+1),NSAM,NNROI1)
C	IF(I.LT.5) WRITE(3,201)(BUF(K),K=1,6*NSAM)
C201	FORMAT(16F7.3)
	IF(I.LT.NROW)GOTO 80
	END
