C++*********************************************************************
C
C $$ ROTS3.FOR
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
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE ROTS3(LUN2,Q1,KLX,KNX,KLY,KNY,KLZ,KNZ,PSI,THETA,PHI)

         DIMENSION  Q1(KLX:KNX,KLY:KNY,KLZ:KNZ),Q2(KLX:KNX)
         DIMENSION  IM(3)
         DOUBLE PRECISION  AV,RM(3,3),QR(3),DX,DY,DZ

C        EQUIVALENCE  (IM(1),IX),(IM(2),IY),(IM(3),IZ)

         LEX=KNX-KLX+1

         IF (THETA.EQ.0.0.AND.PHI.EQ.0.0.AND.PSI.EQ.0.0)  THEN
            IBUF=0
            DO IZ=KLZ,KNZ
               DO IY=KLY,KNY
                  IBUF = IBUF + 1
                  CALL WRTLIN(LUN2,Q1(KLX,IY,IZ),LEX,IBUF)
               ENDDO
            ENDDO
         RETURN
         ENDIF

C     AV=0.0
C     DO  1  IZ=KLZ,KNZ
C     DO  1  IY=KLY,KNY
C     DO  1  IX=KLX,KNX
C1    AV=AV+Q1(IX,IY,IZ)
C     AV=AV/FLOAT(KNX-KLX+1)/FLOAT(KNY-KLY+1)/FLOAT(KNZ-KLZ+1)

        CALL BLDR(RM,PSI,THETA,PHI)

         IBUF=0
         DO    IZ=KLZ,KNZ
         DO    IY=KLY,KNY

        QR(1)=RM(1,1)*KLX+RM(2,1)*IY+RM(3,1)*IZ
        QR(2)=RM(1,2)*KLX+RM(2,2)*IY+RM(3,2)*IZ
        QR(3)=RM(1,3)*KLX+RM(2,3)*IY+RM(3,3)*IZ

         DO     IX=KLX,KNX

C         DO  3  I3=1,3
C         QR(I3)=0.0
C         DO  3  I2=1,3
C3        QR(I3)=QR(I3)+RM(I2,I3)*IM(I2)

         IOX=QR(1)+FLOAT(1-KLX)
         DX=QR(1)+FLOAT(1-KLX)-IOX
         DX=DMAX1(DX,1.0D-5)
         IOX=IOX+KLX-1
         IOY=QR(2)+FLOAT(1-KLY)
         DY=QR(2)+FLOAT(1-KLY)-IOY
         DY=DMAX1(DY,1.0D-5)
         IOY=IOY+KLY-1
         IOZ=QR(3)+FLOAT(1-KLZ)
         DZ=QR(3)+FLOAT(1-KLZ)-IOZ
         DZ=DMAX1(DZ,1.0D-5)
         IOZ=IOZ+KLZ-1

         IF(IOX.GE.KLX.AND.IOX.LT.KNX)  THEN
         IF(IOY.GE.KLY.AND.IOY.LT.KNY)  THEN
         IF(IOZ.GE.KLZ.AND.IOZ.LT.KNZ)  THEN

C     Q2(IX)=
C     &  +(1-DX)*(1-DY)*(1-DZ)*Q1(IOX,IOY,IOZ)
C     &  +   DX *(1-DY)*(1-DZ)*Q1(IOX+1,IOY,IOZ)
C     &  +(1-DX)*   DY *(1-DZ)*Q1(IOX,IOY+1,IOZ)
C     &  +(1-DX)*(1-DY)*   DZ *Q1(IOX,IOY,IOZ+1)
C     &  +   DX *   DY *(1-DZ)*Q1(IOX+1,IOY+1,IOZ)
C     &  +   DX *(1-DY)*   DZ *Q1(IOX+1,IOY,IOZ+1)
C     &  +(1-DX)*   DY *   DZ *Q1(IOX,IOY+1,IOZ+1)
C     &  +   DX *   DY *   DZ *Q1(IOX+1,IOY+1,IOZ+1)
C
C faster version :

         A1 = Q1(IOX,IOY,IOZ)
         A2 = Q1(IOX+1,IOY,IOZ) - A1
         A3 = Q1(IOX,IOY+1,IOZ) - A1
         A4 = Q1(IOX,IOY,IOZ+1) - A1
         A5 = -A2 - Q1(IOX,IOY+1,IOZ) + Q1(IOX+1,IOY+1,IOZ)
         A6 = -A2 - Q1(IOX,IOY,IOZ+1) + Q1(IOX+1,IOY,IOZ+1)
         A7 = -A3 - Q1(IOX,IOY,IOZ+1) + Q1(IOX,IOY+1,IOZ+1)
         A8 = -A5 + Q1(IOX,IOY,IOZ+1) - Q1(IOX+1,IOY,IOZ+1)
     &   - Q1(IOX,IOY+1,IOZ+1) + Q1(IOX+1,IOY+1,IOZ+1)
         Q2(IX)= A1 + DZ*(A4 + A6*DX + (A7 + A8*DX)*DY) + A3*DY
     &   + DX*(A2 + A5*DY)

         GOTO  5
         ENDIF
         ENDIF
         ENDIF

C        Q2(IX)=AV
         Q2(IX)=Q1(IX,IY,IZ)
5        CONTINUE
         QR(1) = QR(1) + RM(1,1)
         QR(2) = QR(2) + RM(1,2)
         QR(3) = QR(3) + RM(1,3)
         ENDDO
         IBUF=IBUF+1
         CALL  WRTLIN(LUN2,Q2,LEX,IBUF)
         ENDDO
         ENDDO

         END



