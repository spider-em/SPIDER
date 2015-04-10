C++*********************************************************************
C
C    PRT
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
C    PRT(LUN51,M1,N,JV,V,A)
C
C--*******************************************************************

       SUBROUTINE PRT(LUN51,M1,N,JV,V,A)

 

       IMPLICIT REAL*8    (A-H,O-Z)
       IMPLICIT INTEGER*2 (I-N)
       INTEGER*4          LUN50,LUN51
       DIMENSION          JV(M1),A(M1,M1)
       CHARACTER*4        V(M1)
       DATA   NG1,NG2/9,10/

       L=N
       L1=1
 1     CONTINUE
       IF (L.LE.0) RETURN

       LL=NG1
       IF (L.LT.NG2) LL=L-1
       L2=L1+LL
       L=L-NG2
       WRITE(LUN51,2)(I,I=L1,L2)
 2     FORMAT(/7X,10(I4,'-FUNCTN'))

       WRITE(LUN51,5 )
 5     FORMAT(1X)

       DO  I=1,M1
         WRITE(LUN51,4) JV(I),V(I),(A(I,J),J=L1,L2)
       ENDDO
       L1=L2+1
       GOTO 1

 4     FORMAT(1X,I4,1X,A4,1X,10(1X,G10.4))
       END
