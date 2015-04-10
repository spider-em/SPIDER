C++*********************************************************************
C
C    WFTCIRC     BILL TIVOL
C                MAXNAM                            JUL 14 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C    WFTCIRC(XD,YD,RR)    ancient!!
C
C    CALL TREE:   'SP' -->  DIFF1O --> LATTICE
C                                  --> WFTCIRC --> LMDIF1  --> LMDIF
C
C--*******************************************************************

	SUBROUTINE WFTCIRC(XD,YD,RR)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

	COMMON XX,YY,PLIST
        COMMON /SIZE/ISIZ
        DIMENSION XX(20),YY(20),FF(20),X(3),PLIST(7),IW(3),WA(100)

	INTEGER                  M
	CHARACTER                NULL

	CHARACTER(LEN=MAXNAM) :: DOCNAM
	REAL*8                   FF,XX,YY,X,WA
	DOUBLE PRECISION         T
	EXTERNAL                 X1YR0
	
	NULL   = CHAR(0)
	LUNDOC = 13
        LW     = 100
        N      = 3
        T      = 10E-10

	CALL RDPRMI(M,ID,NOT_USED,'NUMBER OF SPOTS (3<M<21)')
        IF (M .LT. 3) RETURN
        IF (M .GT. 20) THEN
           WRITE(NOUT,*) ' ONLY THE FIRST 20 SPOTS WILL BE USED.'
           M = 20
        ENDIF

1	CONTINUE
        NOPEN = 0
	CALL FILERD(DOCNAM,NLET,NULL,'SPOT DOCUMENT',IRTFLG)
        IF (IRTFLG.EQ.-1.AND.DOCNAM.EQ.NULL) THEN
           WRITE(NOUT,*) 'NO DOCUMENT FILE'
           RETURN
        ENDIF
	IF (IRTFLG.EQ.-1) GO TO 1
	DO  IKEY = 1,M
	   CALL UNSAV(DOCNAM,NOPEN,LUNDOC,IKEY,PLIST,2,LERR,1)
	   NOPEN = 1
	   IF(LERR.NE.0) GO TO 1
	   XX(IKEY) = PLIST(1)*ISIZ
	   YY(IKEY) = PLIST(2)*ISIZ
           WRITE(NOUT,1000) IKEY,XX(IKEY),YY(IKEY)
           WRITE(NDAT,1000) IKEY,XX(IKEY),YY(IKEY)
1000       FORMAT (' FOR THE ',I2,' TH SPOT, THE CO-ORDS ARE ',
     $               2F10.2)
	ENDDO
	XAB = XX(1)-XX(2)
	XBC = XX(2)-XX(3)
	XCA = XX(3)-XX(1)
	YAB = YY(1)-YY(2)
	YBC = YY(2)-YY(3)
	YCA = YY(3)-YY(1)
	XAYBC = XX(1)*YBC
	XBYCA = XX(2)*YCA
	XCYAB = XX(3)*YAB
	YAXBC = YY(1)*XBC
	YBXCA = YY(2)*XCA
	YCXAB = YY(3)*XAB
	X(1) = (XX(1)*XAYBC+XX(2)*XBYCA+XX(3)*XCYAB
     1         -YAB*YBC*YCA)/(2*(XAYBC+XBYCA+XCYAB))
	X(2) = (YY(1)*YAXBC+YY(2)*YBXCA+YY(3)*YCXAB
     1         -XAB*XBC*XCA)/(2*(YAXBC+YBXCA+YCXAB))
	X(3) = SQRT((XX(1)-X(1))**2 + (YY(1)-X(2))**2)
        WRITE(NOUT,1002) X(1),X(2),X(3)
        WRITE(NDAT,1002) X(1),X(2),X(3)
1002    FORMAT(' INITIAL GUESS CALCULATED FROM THE FIRST THREE SPOTS:',/,
     1  ' X COORD OF CENTER = ',E12.5,', Y COORD OF CENTER = ',E12.5,/,
     2  ' RADIUS = ',E12.5)
        IF(M.EQ.3) THEN
           WRITE(NOUT,*) ' 3 SPOTS USED; NO FITTING NECESSARY.'
        ELSE
           CALL LMDIF1(X1YR0,M,N,X,FF,T,IN,IW,WA,LW)
           WRITE(NOUT,1003) IN
           WRITE(NDAT,1003) IN
1003       FORMAT(' THE CONVERGENCE CRITERION IS ',I2)
        ENDIF
	XD = X(1)
	YD = X(2)
	RR = X(3)
        DO  I=1,M
           FF(I) = SQRT((XX(I)-XD)**2+(YY(I)-YD)**2)-RR
           WRITE(NOUT,1001) I,FF(I)
           WRITE(NDAT,1001) I,FF(I)
1001    FORMAT (' FOR THE ',I2,' TH SPOT, THE RADIAL ERROR IS ',F15.10)
	ENDDO
        CLOSE(LUNDOC)
	RETURN
	END

	SUBROUTINE X1YR0(M,N,X,F,IFLAG)
	REAL*8 F(20),XX(20),YY(20),X(3)
	COMMON XX,YY
	DO I=1,M
	   F(I) = SQRT((XX(I)-X(1))**2+(YY(I)-X(2))**2)-X(3)
        END DO
	RETURN
	END 
