C++*********************************************************************
C
C   MATINV
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
C    MATINV
C
c     matrix invertion
C
C--*******************************************************************

	SUBROUTINE MATINV(ARRAY,NTERM,DET)

 
	DIMENSION ARRAY(20,20),IK(20),JK(20)

	DET=1
	DO 100 K=1,NTERM
	AMAX=0
21	DO  I=K,NTERM
	DO  J=K,NTERM
	IF ((ABS(AMAX)-ABS(ARRAY(I,J))) .LE. 0) THEN
	AMAX=ARRAY(I,J)
	IK(K)=I
	JK(K)=J
	ENDIF
	ENDDO	
	ENDDO

	IF(AMAX .EQ. 0) THEN
	DET=0
	GOTO 140
	ENDIF
	I=IK(K)
	IF(I-K) 21,51,43
43	DO  J=1,NTERM
	 SAVE=ARRAY(K,J)
	 ARRAY(K,J)=ARRAY(I,J)
	 ARRAY(I,J)=-SAVE
	ENDDO
51	J=JK(K)
	IF(J-K) 21,61,53
53	DO  I=1,NTERM
	 SAVE=ARRAY(I,K)
	 ARRAY(I,K)=ARRAY(I,J)
	 ARRAY(I,J)=-SAVE
	ENDDO

61	DO 70  I=1,NTERM
	 IF(I-K) 63,70,63
63	 ARRAY(I,K)=-ARRAY(I,K)/AMAX
70	CONTINUE
71	DO 80 I=1,NTERM
	DO 80 J=1,NTERM
	IF(I-K) 74,80,74
74	IF(J-K) 75,80,75
75	ARRAY(I,J)=ARRAY(I,J)+ARRAY(I,K)*ARRAY(K,J)
80	CONTINUE
81	DO 90 J=1,NTERM
	IF(J-K) 83,90,83
83	ARRAY(K,J)=ARRAY(K,J)/AMAX
90	CONTINUE
	ARRAY(K,K)=1./AMAX
100 	DET=DET*AMAX

101	DO 130 L=1,NTERM
	K=NTERM-L+1
	J=IK(K)
	IF(J-K) 111,111,105
105	DO  I=1,NTERM
	SAVE=ARRAY(I,K)
	ARRAY(I,K)=-ARRAY(I,J)
	ARRAY(I,J)=SAVE
	ENDDO
111	I=JK(K)
	IF(I-K) 130,130,113
113	DO  J=1,NTERM
	SAVE=ARRAY(K,J)
	ARRAY(K,J)=-ARRAY(I,J)
	ARRAY(I,J)=SAVE
	ENDDO
130	CONTINUE
140	RETURN
	END
