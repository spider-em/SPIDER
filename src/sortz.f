C++*******************************************************************
C
C  SORTZ.F 
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
C  SORTZ:    SORT PROGRAM  REAL KEY ARRAY, TWO REAL AND ONE INTEGER
C		ARRAY TO BE SORTED.
C
C  PURPOSE: A SINGLETON SORT PROGRAM TO ORDER B,C AND D USING 'A' 
C           D USING A AS A KEY.  AS OF THE PRESENT TIME (FEB. 1971) 
C           THIS IS THE FASTEST GENERAL PURPOSE SORTING METHOD KNOWN.
C
C  PARAMETERS:   A:   SORTING KEY                           (SENT/RET.)
C                B:   SORTED (REAL)                         (SENT/RET.)
C                C:   SORTED (REAL)                         (SENT/RET.)
C                D:   SORTED (INT.)                         (SENT/RET.)
C                N:   NUMBER OF ELEMENTS TO BE SORTED            (SENT)
C
C--********************************************************************

      SUBROUTINE SORTZ ( A, B, C,D,N)

      REAL A(N),B(N),C(N)
      INTEGER D(N)

      INTEGER IL(16), IU(16),Z
C
      M = 1
      I = 1
      J = N
    5 IF (I .GE. J) GO TO 70

C     ORDER THE TWO ENDS AND THE MIDDLE

   10 K = I
      IJ = (I + J)/2
      T = A(IJ)
      IF (A(I) .LE. T) GO TO 20
      A(IJ) = A(I)
      A(I) = T
      T = A(IJ)
      X = B(I)
      B(I) = B(IJ)
      B(IJ) = X
      Y = C(I)
      C(I) = C(IJ)
      C(IJ) = Y
      Z = D(I)
      D(I) = D(IJ)
      D(IJ) = Z
   20 L = J
      IF (A(J) .GE. T) GO TO 40
      IF (A(J) .LT. A(I)) GO TO 25
      A(IJ) = A(J)
      A(J) = T
      T = A(IJ)
      X = B(IJ)
      B(IJ) = B(J)
      B(J) = X
      Y = C(IJ)
      C(IJ) = C(J)
      C(J) = Y
      Z = D(IJ)
      D(IJ) = D(J)
      D(J) = Z
      GO TO 40
   25 A(IJ) = A(I)
      A(I) = A(J)
      A(J) = T
      T = A(IJ)
      X = B(J)
      B(J) = B(IJ)
      B(IJ) = B(I)
      B(I) = X
      Y = C(J)
      C(J) = C(IJ)
      C(IJ) = C(I)
      C(I) = Y
       Z = D(J)
      D(J ) = D(IJ)
      D(IJ) = D(I)
      D(I) = Z
      GO TO 40

C     SPLIT THE SEQUENCE BETWEEN I AND J INTO TWO SEQUENCES.  THAT
C     SEQUENCE BETWEEN I AND L WILL CONTAIN ONLY ELEMENTS LESS THAN OR
C     EQUAL TO T, WHILE THAT BETWEEN K AND J WILL CONTAIN ONLY ELEMENTS
C     GREATER THAN T.

   30 A(L) = A(K)
      A(K) = TT
      X = B(L)
      B(L) = B(K)
      B(K) = X
      Y = C(L)
      C(L) = C(K)
      C(K) = Y
      Z = D(L)
      D(L) = D(K)
      D(K) = Z
   40 L = L - 1
      IF (A(L) .GT. T) GO TO 40
      TT = A(L)
   50 K = K + 1
      IF (A(K) .LT. T) GO TO 50
      IF (K .LE. L) GO TO 30

C     SAVE THE END POINTS OF THE LONGER SEQUENCE IN IL AND IU, AND SORT
C     THE SHORTER SEQUENCE.

      IF (L - I .LE. J - K) GO TO 60
      IL(M) = I
      IU(M) = L
      I = K
      M = M + 1
      GO TO 80
   60 IL(M) = K
      IU(M) = J
      J = L
      M = M + 1
      GO TO 80

C     RETRIEVE END POINTS PREVIOUSLY SAVED AND SORT BETWEEN THEM.

   70 M = M - 1
      IF (M .EQ. 0) RETURN
      I = IL(M)
      J = IU(M)

C     IF THE SEQUENCE IS LONGER THAN 11 OR IS THE FIRST SEQUENCE, SORT
C     BY SPLITTING RECURSIVELY.

   80 IF (J - I .GE. 11) GO TO 10
      IF (I .EQ. 1) GO TO 5

C     IF THE SEQUENCE IS 11 OR LESS LONG, SORT IT BY A SHELLSORT.

      I = I - 1
   90 I = I + 1
      IF (I .EQ. J) GO TO 70
      T = A(I + 1)
      IF (A(I) .LE. T) GO TO 90
      X = B(I+1)
      Y = C(I+1)
      Z = D(I+1)
      K = I
  100 A(K+1) = A(K)
      B(K+1) = B(K)
      C(K+1) = C(K)
      D(K+1) = D(K)
      K = K - 1
      IF (T .LT. A(K)) GO TO 100
      A(K+1) = T
      B(K+1) = X
      C(K+1) = Y
      D(K+1) = Z
      GO TO 90

      END
