
C ++********************************************************************
C                                                                      *
C  DIRNM                                                             *
C                                                                      *
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
C                                                                      *
C DIRNM                                                                *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

       SUBROUTINE DIRNM(A,M,B,X,XL,XT,IQ)
 
       IMPLICIT REAL*8 (A-H,O-Z)
       IMPLICIT INTEGER*2 (I-N)
       DIMENSION A(M,M),B(M,M),X(M,M),XL(M),XT(M),IQ(M)

       IZERO=0
       CALL HDIAG(B,M,M,IZERO,X,NR,XT,IQ)
       DO  I=1,M
        XL(I)=1./SQRT(ABS(B(I,I)))
        DO  J=1,M
         B(J,I)=X(J,I)*XL(I)
        ENDDO
       ENDDO 
       DO  I=1,M
        DO  J=1,M
        X(I,J)=0.0
        DO  K=1,M
         X(I,J)=X(I,J)+B(K,I)*A(K,J)
	ENDDO
	ENDDO
	ENDDO
       DO  I=1,M
        DO  J=1,M
        A(I,J)=0.0
         DO  K=1,M
         A(I,J)=A(I,J)+X(I,K)*B(K,J)
	 ENDDO
	ENDDO
	ENDDO
       CALL HDIAG(A,M,M,IZERO,X,NR,XT,IQ)
       DO  I=1,M
        XL(I)=A(I,I)
         DO  J=1,M
         A(I,J)=0.0
          DO  K=1,M
          A(I,J)=A(I,J)+B(I,K)*X(K,J)
	  ENDDO
	 ENDDO
	ENDDO 
      DO  I=1,M
        SUMV=0.0
        DO  J=1,M
         SUMV=SUMV+A(J,I)*A(J,I)
        ENDDO 
        DEN=1./SQRT(SUMV)
        DO  J=1,M
         X(J,I)=A(J,I)*DEN
	ENDDO 
      ENDDO
       END
