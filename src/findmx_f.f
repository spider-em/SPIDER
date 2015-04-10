C++*********************************************************************
C
C $$ FINDMX_F.FOR
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
C $$ FINDMX_F.FOR
C
         SUBROUTINE  FINDMX_F(D,NSAM,NROW,NSI,CMX,SX,SY,JACUP)
         DIMENSION  D(NSAM,NROW),Z(-1:1,-1:1)

         JC=NROW/2+1
         IC=NSAM/2+1
         CMX=D(IC,JC)
         DO    JT=-NSI,NSI
            J=JT+JC
            DO    IT=-NSI,NSI
               I=IT+IC
               IF(CMX.LE.D(I,J))  THEN
                  CMX=D(I,J)
                  IX=I
                  IY=J
               ENDIF
            ENDDO         
         ENDDO
         SX=IX-IC
         SY=IY-JC
         IF(IY.LT.2.OR.IY.GT.NROW-1.OR.IX.LT.2.OR.IX.GT.NSAM-1) RETURN
         DO    J=-1,1
            DO    I=-1,1
               Z(I,J)=D(IX+I,IY+J)
            ENDDO
         ENDDO
         CALL  PARABL(Z,XSH,YSH,CMX)
        K=INT(100.0/FLOAT(JACUP+1))
        SX=SX
     &      +REAL(IFIX(XSH))
     &  +FLOAT(K)*FLOAT(INT(XSH*100.0/FLOAT(K)))/100.0
        SY=SY
     &      +REAL(IFIX(YSH))
     &  +FLOAT(K)*FLOAT(INT(YSH*100.0/FLOAT(K)))/100.0
         END
