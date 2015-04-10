C++*********************************************************************
C
C $$ MALFI3.FOR
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
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************
C
C $$ MALFI3.FOR
C
         SUBROUTINE MALFI3(A1,NSAM,NROW,NSLICE)
C A1*CONJG(A1) = |A1|^2
         DIMENSION  A1(NSAM,NROW,NSLICE)
        LOGICAL IFNS,IFNR,IFND
         COMPLEX  QT,FQ

C
	PI=4.0*ATAN(1.0)
        IFNR=MOD(NROW,2).EQ.0
        IFND=MOD(NSAM,2).EQ.0
        IFNS=MOD(NSLICE,2).EQ.0
        IF(IFND)  THEN
        LBD=2
        ELSE
        LBD=1
        ENDIF
        IF(IFNR)  THEN
        LBR=2
        ELSE
        LBR=1
        ENDIF
        IF(IFNS)  THEN
        LBS=2
        ELSE
        LBS=1
        ENDIF
        ND2=NSAM/2
        NR2=NROW/2
        NS2=NSLICE/2
C
         DO    I=1,LBS
	IZ=(I-1)*NS2
	ZZ=2.0*PI*FLOAT(IZ)*FLOAT(NS2)/FLOAT(NSLICE)
         DO    J=1,LBR
	IY=(J-1)*NR2
	YY=2.0*PI*FLOAT(IY)*FLOAT(NR2)/FLOAT(NROW)
         DO    K=1,LBD
	IX=(K-1)*ND2
	XX=2.0*PI*FLOAT(IX)*FLOAT(ND2)/FLOAT(NSAM)
	FQ=CMPLX(COS(XX+YY+ZZ),SIN(XX+YY+ZZ))
C	FQ=CMPLX(1.0,0.0)
	QT=CMPLX(A1(K,J,I),0.0)*CMPLX(A1(K,J,I),0.0)*FQ
	A1(K,J,I)=REAL(QT)
	ENDDO
	ENDDO
	ENDDO
C
         DO    I=1,NSLICE
	IZ=(I-1)
	IF(IZ.GT.NS2)IZ=IZ-NSLICE
	ZZ=2.0*PI*FLOAT(IZ)*FLOAT(NS2)/FLOAT(NSLICE)
         DO    K=1,LBD
	IX=(K-1)*ND2
	XX=2.0*PI*FLOAT(IX)*FLOAT(ND2)/FLOAT(NSAM)
         DO  J=3,NROW-1,2
	IY=(J-1)/2
	YY=2.0*PI*FLOAT(IY)*FLOAT(NR2)/FLOAT(NROW)
	FQ=CMPLX(COS(XX+YY+ZZ),SIN(XX+YY+ZZ))
C	FQ=CMPLX(1.0,0.0)
	QT=CMPLX(A1(K,J,I),A1(K,J+1,I))*
     &  CMPLX(A1(K,J,I),-A1(K,J+1,I))*FQ
         A1(K,J,I)=REAL(QT)
         A1(K,J+1,I)=AIMAG(QT)
	 ENDDO
	IF(.NOT.IFNR)THEN
	IY=NR2
	YY=2.0*PI*FLOAT(IY)*FLOAT(NR2)/FLOAT(NROW)
	FQ=CMPLX(COS(XX+YY+ZZ),SIN(XX+YY+ZZ))
C	FQ=CMPLX(1.0,0.0)
	QT=CMPLX(A1(K,NROW,I),A1(K,2,I))*
     &  CMPLX(A1(K,NROW,I),-A1(K,2,I))*FQ
         A1(K,NROW,I)=REAL(QT)
         A1(K,2,I)=AIMAG(QT)
	ENDIF
	ENDDO
	ENDDO
C
         DO    J=1,LBR
	IY=(J-1)*NR2
	YY=2.0*PI*FLOAT(IY)*FLOAT(NR2)/FLOAT(NROW)
         DO    K=1,LBD
	IX=(K-1)*ND2
	XX=2.0*PI*FLOAT(IX)*FLOAT(ND2)/FLOAT(NSAM)
         DO    I=3,NSLICE-1,2
	IZ=(I-1)/2
	ZZ=2.0*PI*FLOAT(IZ)*FLOAT(NS2)/FLOAT(NSLICE)
	FQ=CMPLX(COS(XX+YY+ZZ),SIN(XX+YY+ZZ))
C	FQ=CMPLX(1.0,0.0)
        QT=CMPLX(A1(K,J,I),A1(K,J,I+1))*
     &  CMPLX(A1(K,J,I),-A1(K,J,I+1))*FQ
	A1(K,J,I)=REAL(QT)
        A1(K,J,I+1)=AIMAG(QT)
	ENDDO
	IF(.NOT.IFNS)THEN
	IZ=NS2
	ZZ=2.0*PI*FLOAT(IZ)*FLOAT(NS2)/FLOAT(NSLICE)
	FQ=CMPLX(COS(XX+YY+ZZ),SIN(XX+YY+ZZ))
C	FQ=CMPLX(1.0,0.0)
        QT=CMPLX(A1(K,J,NSLICE),A1(K,J,2))*
     &  CMPLX(A1(K,J,NSLICE),-A1(K,J,2))*FQ
         A1(K,J,NSLICE)=REAL(QT)
         A1(K,J,2)=AIMAG(QT)
	ENDIF
	ENDDO	
	ENDDO
C
         DO    I=1,NSLICE
	IZ=(I-1)
	IF(IZ.GT.NS2)IZ=IZ-NSLICE
	ZZ=2.0*PI*FLOAT(IZ)*FLOAT(NS2)/FLOAT(NSLICE)
         DO    J=1,NROW
	IY=(J-1)
	IF(IY.GT.NR2)IY=IY-NROW
	YY=2.0*PI*FLOAT(IY)*FLOAT(NR2)/FLOAT(NROW)
         DO    K=3,NSAM-1,2
	IX=(K-1)/2
	XX=2.0*PI*FLOAT(IX)*FLOAT(ND2)/FLOAT(NSAM)
	FQ=CMPLX(COS(XX+YY+ZZ),SIN(XX+YY+ZZ))
C	FQ=CMPLX(1.0,0.0)
	QT=CMPLX(A1(K,J,I),A1(K+1,J,I))*
     &  CMPLX(A1(K,J,I),-A1(K+1,J,I))*FQ
         A1(K,J,I)=REAL(QT)
        A1(K+1,J,I)=AIMAG(QT)
	ENDDO
	IF(.NOT.IFND)THEN
	IX=ND2
	XX=2.0*PI*FLOAT(IX)*FLOAT(ND2)/FLOAT(NSAM)
	FQ=CMPLX(COS(XX+YY+ZZ),SIN(XX+YY+ZZ))
C	FQ=CMPLX(1.0,0.0)
	QT=CMPLX(A1(NSAM,J,I),A1(2,J,I))*
     &  CMPLX(A1(NSAM,J,I),-A1(2,J,I))*FQ
        A1(NSAM,J,I)=REAL(QT)
         A1(2,J,I)=AIMAG(QT)
	ENDIF
	ENDDO
	ENDDO
         END
