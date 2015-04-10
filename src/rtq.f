C++*********************************************************************
C
C RTQ.F
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
C IMAGE_PROCESSING_ROUTINE                                                                     
C        0         2         3         4         5         6         7 
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE  RTQ(X,OUT,NSAM,NROW,THETA)

         DIMENSION  X(NSAM,NROW),OUT(NSAM,NROW)
         PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
         PARAMETER (DGR_TO_RAD = (QUADPI/180))

         IF(THETA.EQ.0.0)  THEN
            CALL  COP(X,OUT,NSAM*NROW)
            RETURN
         ENDIF
         ICENT=NROW/2+1
         KCENT=NSAM/2+1
         RN2=-NROW/2
         SN2=-NSAM/2
         RW2=-RN2
         RS2=-SN2
         IF(MOD(NSAM,2).EQ.0)  RW2=RW2-1.0
         IF(MOD(NROW,2).EQ.0)  RS2=RS2-1.0
         COD=COS(THETA*DGR_TO_RAD)
         SID=SIN(THETA*DGR_TO_RAD)
c$omp parallel do private(i,k,yi,ycod,ysid,xi,xold,yold)
         DO    I=1,NROW
            YI=I-ICENT
            IF(YI.LT.RN2)  YI=AMIN1(RW2+YI-RN2+1.0,RW2)
            IF(YI.GT.RW2)  YI=AMAX1(RN2+YI-RW2-1.0,RN2)
            YCOD=YI*COD+ICENT
            YSID=-YI*SID+KCENT
            DO    K=1,NSAM
               XI=K-KCENT
               IF(XI.LT.SN2)  XI=AMIN1(RS2+XI-SN2+1.0,RS2)
               IF(XI.GT.RS2)  XI=AMAX1(SN2+XI-RS2-1.0,SN2)
               YOLD=XI*SID+YCOD
               XOLD=XI*COD+YSID
               OUT(K,I)=QUADRI(XOLD, YOLD, NSAM, NROW, X)
            ENDDO
         ENDDO
         END

