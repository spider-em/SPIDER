C
C++*********************************************************************
C
C FMRS_2R.F                     private(j)          OCT 01 ARDEAN LEITH
 
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
C  Order of elements:
C
C--*********************************************************************

         SUBROUTINE  FMRS_2R(X,NNNN,NSAM,NROW,INV)

         DIMENSION  X(NNNN,NROW)

	INS=INV*NNNN

        IF (INV .GE. 0)  THEN
c$omp      parallel do private(i),shared(invt) 
	   DO I=1,NROW
	      INVT = INV
	      CALL FMRS_1(X(1,I),NSAM,INVT)
	   ENDDO
           IF (INVT .LE. 0)  THEN
	      INV = INVT
	      RETURN
	   ENDIF
	ENDIF

c$omp   parallel do private(i),shared(invt) 
        DO I=1,NNNN,2
	   INVT = INS
	   CALL FFTMCF(X(I,1),X(I+1,1),NROW,NROW,NROW,INVT)
	ENDDO

        IF (INVT.EQ.0)  THEN
	   INV = 0
	   RETURN
	ENDIF
	IF (INV.GT.0)  RETURN

C       NORMALIZE FOR INVERSE
        Q = 1.0/FLOAT(NROW)
c$omp   parallel do private(i,j)
        DO J=1,NROW
           DO I=1,NNNN
              X(I,J)=X(I,J)*Q
	   ENDDO
	ENDDO 

c$omp   parallel do private(i),shared(invt) 
	DO I=1,NROW
	   INVT=INV
	   CALL  FMRS_1(X(1,I),NSAM,INVT)
	ENDDO

        IF(INVT.LE.0)  INV=INVT

        END
