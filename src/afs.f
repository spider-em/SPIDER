C++*********************************************************************
C
C AFS.F          IMPLICIT                        FEB 2020 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C  AFS(X,OUT,NX,NY, A,B,C,D,SHXI,SHYI,LUN2,LB)
C
C--*********************************************************************

        SUBROUTINE  AFS(X,OUT, NX,NY, A,B,C,D, SHXI,SHYI, LUN2,LB)

        IMPLICIT NONE

        REAL    :: X(NX,NY),OUT(NX)
        INTEGER :: NX,NY,LUN2,LB
        REAL    :: A,B,C,D,SHXI,SHYI

        REAL    :: AT,BT,CT,DT,DET,XI,YI,XOLD,YOLD
        INTEGER :: ICENT,KCENT,I,K,IX,IY,ISLICE

        REAL    :: quadri  ! FUNCTION

C       DETERMINANT WAS CHECKED IN THE CALLING PROGRAM
	DET = A*D-B*C

	AT    =  A/DET
	DT    =  D/DET
	BT    = -B/DET
	CT    = -C/DET
        ICENT = NY/2+1
        KCENT = NX/2+1

        DO I=1,NY
           YI = I-ICENT-SHYI

c$omp      parallel do private(k,xi,xold,yold)
           DO K=1,NX
              XI     = K-KCENT-SHXI
              XOLD   = AT*XI+BT*YI+ICENT
              YOLD   = CT*XI+DT*YI+KCENT
              OUT(K) = QUADRI(XOLD, YOLD, NX, NY, X)
	   ENDDO

C          OUTPUT THE ROW
           CALL WRTLIN(LUN2,OUT,NX,LB+I)
	 ENDDO

         END
