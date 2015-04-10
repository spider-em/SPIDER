C++*********************************************************************
C
C AFS.FOR
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
C--*********************************************************************

        SUBROUTINE  AFS(X,OUT,NSAM,NROW,
     &                        A,B,C,D,SHXI,SHYI,LUN2,LB)

        DIMENSION  X(NSAM,NROW),OUT(NSAM)

	DET=A*D-B*C
C       DETERMINANT WAS CHECKED IN THE CALLING PROGRAM

	AT=A/DET
	DT=D/DET
	BT=-B/DET
	CT=-C/DET
        ICENT=NROW/2+1
        KCENT=NSAM/2+1
        DO I=1,NROW
           YI=I-ICENT-SHYI
c$omp parallel do private(k,xi,xold,yold)
           DO K=1,NSAM
              XI=K-KCENT-SHXI
              XOLD=AT*XI+BT*YI+ICENT
              YOLD=CT*XI+DT*YI+KCENT
              OUT(K)=QUADRI(XOLD, YOLD, NSAM, NROW, X)
	   ENDDO
           CALL  WRTLIN(LUN2,OUT,NSAM,LB+I)
	 ENDDO
         END
