
C++*********************************************************************
C
C  SOLVE.F
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
C  SOLVE(A,NFU,NDM)
C
C  PURPOSE: SOLVES A SET OF SIMULTANEOUS EQUATIONS OF THE
C  FORM: AX=B; WHERE:
C  A = MATRIX OF COEFFICIENTS
C  X = VECTOR OF SOLUTIONS
C  B = VECTOR OF KNOWNS
C
C  PARAMETERS:
C  A IS THE AUGMENTED MATRIX AND NFU IS THE NUMBER OF EQUATIONS.
C  A MUST BE DIMENSIONED (NDM,NDM+1) IN TNE MAIN PROGRAM
C  THE SOLUTION IS RETURNED IN (A(I,NFU+1),I=1,NFU).  THE MATRIX IS
C  DESTROYED AND NO CHECK IS MADE TO SEE IF THE MATRIX IS SINGULAR.
C  A SINGULAR MATRIX WILL PROBABLY RESULT IN A FORTRAN ERROR MESSAGE.
C  THE METHOD OF GAUSSIAN ELIMINATION WITH MAXIMUM PIVOT IS USED.
C
C***********************************************************************

	SUBROUTINE SOLVE(A,NFU,NDM)

	DOUBLE PRECISION A(NDM,*)
	DOUBLE PRECISION F,AF,G

	COMMON/UNITS/LUN,NIN,NOUT

	N   = NFU
	N1  = N + 1
	DO  I=1,N

C         FIND THE LARGEST PIVOT

	  F = A(I,I)
C	  WRITE(NOUT,333)F
C333	  FORMAT( ' F:  ',F8.3)

	  AF = ABS(F)
	  JJ = I
	  DO J=I,N
	    IF (ABS(A(J,I)) .GT. AF) THEN
	       F = A(J,I)
C	       WRITE(NOUT,333)F
	       AF = ABS(F)
	       JJ = J
            ENDIF
         ENDDO

C        IF THE LARGEST PIVOT IS IN A DIFFERENT ROW,
C        THEN SWITCH ALL ELEMENTS OF THESE ROWS

	  IF (JJ .NE. I) THEN
	     DO  J=1,N1
	        G       = A(I,J)
	        A(I,J)  = A(JJ,J)
	        A(JJ,J) = G
	     ENDDO
          ENDIF

C         DIVIDE THE REST OF THE ROW BY THE PIVOT
   	  DO  J=I,N1
C            WRITE(NOUT,334)F
C334         FORMAT( ' F IN DENOMINATOR :  ',F8.3)
	     A(I,J) = A(I,J) / F
	  ENDDO

C         SUBTRACT THIS ROW (MULTIPLIED BY A(K,I)) FROM ALL REMAINING ROWS
	  IF (I .NE. N) THEN
	     I1 = I+1
	     DO  K=I1,N
	        F = A(K,I)
	        DO  J=I,N1
	           A(K,J) = A(K,J) - A(I,J) * F
	        ENDDO
	     ENDDO
          ENDIF
       ENDDO

C      SUBTRACT UPWARDS
       NM1 = N - 1
       DO  KM=1,NM1
          K   = N1-KM
	  KM1 = K-1
	  F   = A(K,N1)
	  DO  I=1,KM1
	     A(I,N1) = A(I,N1) - F * A(I,K)
	  ENDDO
       ENDDO

       RETURN
       END
