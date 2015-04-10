C++*********************************************************************
C
C $$ SHFI_Q.FOR
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
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************
C
C $$ SHFI_Q.FOR
C
         SUBROUTINE  SHFI_Q(X,OUT,LSD,NSAM,NROW,ISHX,ISHY)
         DIMENSION  X(LSD,NROW),OUT(LSD,NROW)
         INTEGER  RN2,SN2,RW2,RS2,XI,YI,SHX,SHY,YOLD,XOLD
C
         SHX=MOD(ISHX,NSAM)
         SHY=MOD(ISHY,NROW)
         ICENT=NROW/2+1
         KCENT=NSAM/2+1
         RN2=-NROW/2
         SN2=-NSAM/2
         RW2=-RN2
         RS2=-SN2
         IF(MOD(NSAM,2).EQ.0)  RW2=RW2-1
         IF(MOD(NROW,2).EQ.0)  RS2=RS2-1
c$omp parallel do private(i,k,yi,yold,xi,xold)
         DO    I=1,NROW
          YI=I-ICENT-SHY
          IF(YI.LT.RN2)  YI=MIN0(RW2+YI-RN2+1,RW2)
          IF(YI.GT.RW2)  YI=MAX0(RN2+YI-RW2-1,RN2)
          YOLD=YI+ICENT
           DO    K=1,NSAM
            XI=K-KCENT-SHX
            IF(XI.LT.SN2)  XI=MIN0(RS2+XI-SN2+1,RS2)
            IF(XI.GT.RS2)  XI=MAX0(SN2+XI-RS2-1,SN2)
            XOLD=XI+KCENT
            OUT(K,I)=X(XOLD, YOLD)
	   ENDDO
	 ENDDO
         END
