C ++********************************************************************
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
C PURPOSE:  CREATE A TRIANGLE
C
C  MTRIAN(LUN,NSAM,NROW,RP,IDIM)
C
C **********************************************************************

	SUBROUTINE MTRIAN(LUN,NSAM,NROW,RP,IDIM)

C       I DO NOT KNOW IF SAVE IS NEEDED FEB 99 al
        SAVE

	DIMENSION IX(3),IY(3),K(3)

	COMMON ADUM(80),BUF(1)

	CALL RDPRMI(IX(1),IY(1),NOT_USED,
     &              'COORDINATES OF FIRST  POINT')
	CALL RDPRMI(IX(2),IY(2),NOT_USED,
     &              'COORDINATES OF SECOND POINT')
	CALL RDPRMI(IX(3),IY(3),NOT_USED,
     &              'COORDINATES OF THIRD  POINT')

	DO  I=1,3
	   BUF(I)=1.0
	ENDDO

	DO  I=1,3
	DO 5 L=1,3
	IF(BUF(L).EQ.0.) GOTO 5
	L1=IY(L)
	K1=L
	GOTO 6
5	CONTINUE
6	CONTINUE
	DO 3 J=1,3
	IF(BUF(J).EQ.0.) GOTO 3
	IF(IY(J).GE.L1) GOTO 3
	L1=IY(J)
	K1=J
3	CONTINUE
	K(I)=K1
	BUF(K1)=0.
	ENDDO

	IF(IY(1).EQ.IY(2).AND.IY(2).EQ.IY(3)) GOTO 9000

	IYSTRT=MAX0(1,IY(K(1)))
	IYEND=MIN0(NROW,IY(K(3)))
	IF(IYSTRT.GT.NROW.OR.IYEND.LE.0) GOTO 9000

	S1=1.
	IF(IX(K(1)).GT.IX(K(2))) S1=-1.
	S2=1.
	IF(IX(K(1)).GT.IX(K(3))) S2=-1.

	XA=FLOAT(IX(K(1)))
	YA=FLOAT(IY(K(1)))
	XB=FLOAT(IX(K(2)))
	YB=FLOAT(IY(K(2)))
	XC=FLOAT(IX(K(3)))
	YC=FLOAT(IY(K(3)))
	AM2=(XC-XA)/(YC-YA)

	IF(IY(K(1)).EQ.IY(K(2))) GOTO 9
	IF(IY(K(2)).LE.0) GOTO 8

	AM1=(XA-XB)/(YA-YB)
	GOTO 9
8	IF(IY(K(2)).EQ.IY(K(3))) GOTO 9000
	AM1=(XC-XB)/(YC-YB)
9	CONTINUE
	DO 10 I=IYSTRT,IYEND
	IF(IY(K(2)).LE.0.OR.I.NE.IY(K(2))+1) GOTO 11
	S1=1.
	IF(IX(K(2)).GT.IX(K(3))) S1=-1.
	AM1=(XC-XB)/(YC-YB)

11	CONTINUE
	A1=XB+AM1*(FLOAT(I)-YB)
	B1=XA+AM2*(FLOAT(I)-YA)
	IA1=IFIX(A1+S1*0.5)
	IB1=IFIX(B1+S2*0.5)
	IA2=MIN0(IA1,IB1)
	IB2=MAX0(IA1,IB1)
	IA3=MAX0(1,IA2)
	IB3=MIN0(NSAM,IB2)
	IF(IA3.GT.NSAM.OR.IB3.LE.0) GOTO 10
	CALL REDLIN(LUN,BUF,NSAM,I)
	IF(IDIM.EQ.1) GOTO 13
	DO  J=IA3,IB3
	BUF(J)=RP
	ENDDO
	GOTO 14
13	IF(IA2.EQ.IA3) BUF(IA2)=RP
	IF(IB2.EQ.IB3) BUF(IB3)=RP
14	CALL WRTLIN(LUN,BUF,NSAM,I)
10	CONTINUE
	RETURN

9000	CALL ERRT(14,'MTRIAN',NF)
	RETURN
	END
