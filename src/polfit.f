C++*********************************************************************
C
C  POLFIT.F           SUMX,SUM DIMENSIONS CORRECTED MAY 99 ARDEAN LEITH
C                     A REMOVED FROM COMMON         MAY 02 ARDEAN LEITH
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
C  POLFIT(X,Y,NORD,N,C)
C
C  PURPOSE:  MAKES LEAST SQUARES FIT OF EXPERIMENTAL DATA IN
C            X(N),Y(N) USING A POLYNOMIAL OF ARBITRARY ORDER>1
C
C  PARAMETERS:
C	  X	A REAL ARRAY DIMENSIONED N CONTAINING THE ABSCISSAE
C	  Y	A REAL ARRAY DIMENSIONED N CONTAINING THE FUNCTION VALUES
C         NORD  ORDER OF POLYNOMIAL. NORD<N.
C	  C	A REAL ARRAY CONTAINING THE NORD+1 COEFFICIENTS OF THE
C	  	POLYNOMIAL IN INCREASING ORDER
C
C--*********************************************************************

      SUBROUTINE POLFIT(X,Y,NORD,N,C)

      DIMENSION X(*),Y(*),C(*)

      DOUBLE PRECISION SUM(22),SUMY(22),A(132)

C     MAXIMUM NSUM: 2*10+1 al
      NSUM  = 2 * NORD + 1

      NEXP  = 2 * NORD
      NORD1 = NORD + 1
      I2    = NORD1 * NORD1
      DO  I=1,NSUM
	 SUMY(I) = 0.0
	 SUM(I)  = 0.0
      ENDDO
      SUM(1) = N
      DO  I=1,N
         DO  IEXP=1,NEXP
	    SUM(IEXP + 1) = SUM(IEXP + 1) + X(I) ** IEXP
	 ENDDO
	 SUMY(1) = SUMY(1) + Y(I)

	 DO  NY=2,NORD1
	    SUMY(NY) = SUMY(NY) + Y(I) * X(I) ** (NY - 1)
	 ENDDO
      ENDDO

C     NOW CONSTRUCT MATRIX
      DO  K=1,NORD1
         DO  I=1,NORD1
C           I1: 0...110 MAXIMUM
            I1      = (I-1)*NORD1
            A(I1+K) = SUM(K+I-1)
         ENDDO
         A(I2+K) = SUMY(K)
      ENDDO

      CALL SOLVE(A,NORD1,NORD1)

      DO  K=1,NORD1
         C(K) = A(I2+K)
C        WRITE(6,*)'C(',K,'): ',C(K)
      ENDDO

      RETURN
      END
