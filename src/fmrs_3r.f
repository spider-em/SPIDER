C
C++*********************************************************************
C
C $$ FMRS_3R.FOR
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
C  For order of elements see fmr_1.
C
C
C--*********************************************************************
C
C $$ FMRS_3R.FOR
C
	SUBROUTINE FMRS_3R(A,NNNN,NSAM,NROW,NSLICE,INV)	
	DIMENSION A(NNNN,NROW,NSLICE)
c
	NDR=INV*NNNN*NROW
C
	IF(INV.GE.0)  THEN
	 DO   I=1,NSLICE
	  CALL FMRS_2(A(1,1,I),NSAM,NROW,INV)
	  IF(INV.EQ.0)   RETURN
	 ENDDO
	ENDIF
C
c$omp parallel do private(j,i),shared(ndrt)
	DO    J=1,NROW
	DO    I=1,NNNN-1,2
	NDRT=NDR
	CALL FFTMCF(A(I,J,1),A(I+1,J,1),NSLICE,NSLICE,NSLICE,NDRT)
	ENDDO	
	ENDDO
	IF(NDRT.EQ.0)  THEN
	INV=0
	RETURN
	ENDIF
	IF(INV.GT.0)  RETURN
C NORMALIZE FOR INVERSE
	Q=1.0/FLOAT(NSLICE)
c$omp parallel do private(k,j,i)
	DO    K=1,NSLICE
	DO    J=1,NROW
	DO    I=1,NNNN
	A(I,J,K)=A(I,J,K)*Q
	ENDDO
	ENDDO
	ENDDO
	DO   I=1,NSLICE
	CALL FMRS_2(A(1,1,I),NSAM,NROW,INV)
	IF(INV.EQ.0)   RETURN
	ENDDO
	END
