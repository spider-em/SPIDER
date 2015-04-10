
C++*********************************************************************
C
C $$ PLNEDG.FOR
c
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
C $$ PLNEDG(X,Y,B)
C
C    CALLED BY:    SPOTWT
C
C        0         2         3         4         5         6         7     
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C***********************************************************************

	SUBROUTINE PLNEDG(X,Y,B)

 

	COMMON/EDGSUM/PN(3,4)
	COMMON/UNITS/LUN,MIN,NOUT

C	WRITE(NOUT,522)(PN(1,I),I=1,4)
C522	FORMAT( ' PN(1,I)BEFORE,I=1,4 :  ',4(F8.3,1X))
C513	FORMAT( ' PN(1,I)AFTER,I=1,4 : ',4(F8.3,1X))

	PN(1,4)=B*X+PN(1,4)
	PN(2,4)=B*Y+PN(2,4)
	PN(3,4)=B+PN(3,4)
	PN(1,1)=X*X+PN(1,1)
	PN(1,2)=X*Y+PN(1,2)
	PN(1,3)=X+PN(1,3)
	PN(2,1)=PN(1,2)
	PN(2,2)=Y*Y+PN(2,2)
	PN(2,3)=Y+PN(2,3)
	PN(3,1)=PN(1,3)
	PN(3,2)=PN(2,3)
	PN(3,3)=1.+PN(3,3)

C	WRITE(NOUT,513)(PN(1,I),I=1,4)
	RETURN
	END
