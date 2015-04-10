C++*******************************************************************
C
C $$ FSORT.FOR
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
C $$ SORT:    SINGLETON SORT
C
C     CALL FSORT(A,B,C,N)
C       A    ARRAY TO BE SORTED
C       N    
C
C--*******************************************************************
C
C
      SUBROUTINE FSORT ( A,N)
      REAL A(N)
C
C     SINGLETON SORT PROGRAM TO ORDER B AND C USING A AS A KEY
C     AS OF THE PRESENT TIME (FEB. 1971) THIS IS THE FASTEST GENERAL
C     PURPOSE SORTING METHOD KNOWN.
C MODIFIED VERSION WITH REAL KEY ARRAY. J.FRANK, FEB. 1977
C
      INTEGER IL(16), IU(16)

C
      M = 1
      I = 1
      J = N
    5 IF (I .GE. J) GO TO 70
C
C     ORDER THE TWO ENDS AND THE MIDDLE
C
   10 K = I
      IJ = (I + J)/2
      T = A(IJ)
      IF (A(I) .LE. T) GO TO 20
      A(IJ) = A(I)
      A(I) = T
      T = A(IJ)
   20 L = J
      IF (A(J) .GE. T) GO TO 40
      IF (A(J) .LT. A(I)) GO TO 25
      A(IJ) = A(J)
      A(J) = T
      T = A(IJ)
      GO TO 40
   25 A(IJ) = A(I)
      A(I) = A(J)
      A(J) = T
      T = A(IJ)
      GO TO 40
C
C     SPLIT THE SEQUENCE BETWEEN I AND J INTO TWO SEQUENCES.  THAT
C     SEQUENCE BETWEEN I AND L WILL CONTAIN ONLY ELEMENTS LESS THAN OR
C     EQUAL TO T, WHILE THAT BETWEEN K AND J WILL CONTAIN ONLY ELEMENTS
C     GREATER THAN T.
C
   30 A(L) = A(K)
      A(K) = TT
   40 L = L - 1
      IF (A(L) .GT. T) GO TO 40
      TT = A(L)
   50 K = K + 1
      IF (A(K) .LT. T) GO TO 50
      IF (K .LE. L) GO TO 30
C
C     SAVE THE END POINTS OF THE LONGER SEQUENCE IN IL AND IU, AND SORT
C     THE SHORTER SEQUENCE.
C
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
C
C     RETRIEVE END POINTS PREVIOUSLY SAVED AND SORT BETWEEN THEM.
C
   70 M = M - 1
      IF (M .EQ. 0) RETURN
      I = IL(M)
      J = IU(M)
C
C     IF THE SEQUENCE IS LONGER THAN 11 OR IS THE FIRST SEQUENCE, SORT
C     BY SPLITTING RECURSIVELY.
C
   80 IF (J - I .GE. 11) GO TO 10
      IF (I .EQ. 1) GO TO 5
C
C     IF THE SEQUENCE IS 11 OR LESS LONG, SORT IT BY A SHELLSORT.
C
      I = I - 1
   90 I = I + 1
      IF (I .EQ. J) GO TO 70
      T = A(I + 1)
      IF (A(I) .LE. T) GO TO 90
      K = I
  100 A(K+1) = A(K)
      K = K - 1
      IF (T .LT. A(K)) GO TO 100
      A(K+1) = T
      GO TO 90
C
      END
