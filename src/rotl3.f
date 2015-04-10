C++*********************************************************************
C
C ROTL3.F
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
C  ROTL3(LUN2,Q1,NSAM,NROW,NSLICE,X1,X2,ALPHA)
C
C  ROTATE ABOUT LINE WITH TRI-LINEAR INTERPOLATION
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE ROTL3(LUN2,Q1,NSAM,NROW,NSLICE,X1,X2,ALPHA)

         DIMENSION  Q1(NSAM,NROW,NSLICE)
         DIMENSION  Q2(NSAM)
         DIMENSION  X1(3),X2(3)

         DOUBLE PRECISION  AV,R1(3,3),R2(3,3),R3(3,3),QR(3),DX,DY,DZ
	 PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
	 PARAMETER (RAD_TO_DGR = (180.0/QUADPI))

	x=x2(1)-x1(1)
	y=x2(2)-x1(2)
	z=x2(3)-x1(3)
        IF(ALPHA.EQ.0.0.OR.(X.EQ.0.0.AND.Y.EQ.0.0.AND.Z.EQ.0.0)) THEN
           IBUF=0
           DO IZ=1,NSLICE
              DO IY=1,NROW
                 IBUF=IBUF+1
                 CALL WRTLIN(LUN2,Q1(1,IY,IZ),NSAM,IBUF)
	      ENDDO
	   ENDDO
           RETURN
        ENDIF

	PSI   = -RAD_TO_DGR*ATAN2(Y,X)
	THETA = RAD_TO_DGR*ATAN2(Z,SQRT(X*X+Y*Y))
	CALL  BLDR(R1,PSI,THETA,90.0)
	CALL  BLDR(R2,0.0,ALPHA,0.0)

         DO I=1,3
            DO J=1,3
               R3(J,I) = 0.0
               DO K=1,3
                  R3(J,I)=R3(J,I)+R2(K,I)*R1(J,K)
	       ENDDO
	    ENDDO
	 ENDDO
         DO I=1,3
            DO J=1,3
               R2(J,I)=0.0
               DO K=1,3
                  R2(J,I)=R2(J,I)+R1(I,K)*R3(J,K)
	       ENDDO
	    ENDDO
	 ENDDO

	IBUF=0
	DO IZ=1,NSLICE
	   ZZ=IZ-X1(3)
	   DO IY=1,NROW
	      YY=IY-X1(2)
               DO 5  IX=1,NSAM
                 XX = IX-X1(1)

                 QR(1)=R2(1,1)*XX+R2(2,1)*YY+R2(3,1)*ZZ
                 QR(2)=R2(1,2)*XX+R2(2,2)*YY+R2(3,2)*ZZ
                 QR(3)=R2(1,3)*XX+R2(2,3)*YY+R2(3,3)*ZZ
                 QR(1)=QR(1)+X1(1)
                 QR(2)=QR(2)+X1(2)
                 QR(3)=QR(3)+X1(3)

                  IOX=QR(1)
                  DX=QR(1)-IOX
                  DX=DMAX1(DX,1.0D-5)
                  IOY=QR(2)
                  DY=QR(2)-IOY
                  DY=DMAX1(DY,1.0D-5)
                  IOZ=QR(3)
                  DZ=QR(3)-IOZ
                  DZ=DMAX1(DZ,1.0D-5)

                  IF(IOX.GE.1.AND.IOX.LT.NSAM)  THEN
                     IF(IOY.GE.1.AND.IOY.LT.NROW)  THEN
                        IF(IOZ.GE.1.AND.IOZ.LT.NSLICE)  THEN

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
C
                        A1 = Q1(IOX,IOY,IOZ)
                        A2 = Q1(IOX+1,IOY,IOZ) - Q1(IOX,IOY,IOZ)
                        A3 = Q1(IOX,IOY+1,IOZ) - Q1(IOX,IOY,IOZ)
                        A4 = Q1(IOX,IOY,IOZ+1) - Q1(IOX,IOY,IOZ)
         A5 = Q1(IOX,IOY,IOZ) - Q1(IOX+1,IOY,IOZ) - Q1(IOX,IOY+1,IOZ)
     &   + Q1(IOX+1,IOY+1,IOZ)
         A6 = Q1(IOX,IOY,IOZ) - Q1(IOX+1,IOY,IOZ) - Q1(IOX,IOY,IOZ+1)
     &   + Q1(IOX+1,IOY,IOZ+1)
         A7 = Q1(IOX,IOY,IOZ) - Q1(IOX,IOY+1,IOZ) - Q1(IOX,IOY,IOZ+1)
     &   + Q1(IOX,IOY+1,IOZ+1)
         A8 = Q1(IOX+1,IOY,IOZ) + Q1(IOX,IOY+1,IOZ)+ Q1(IOX,IOY,IOZ+1)
     &   - Q1(IOX,IOY,IOZ)- Q1(IOX+1,IOY+1,IOZ) - Q1(IOX+1,IOY,IOZ+1)
     &   - Q1(IOX,IOY+1,IOZ+1) + Q1(IOX+1,IOY+1,IOZ+1)
         Q2(IX)= A1 + DZ*(A4 + A6*DX + (A7 + A8*DX)*DY) + A3*DY
     &   + DX*(A2 + A5*DY)
C**********************************************************
                        GOTO  5
                     ENDIF
                  ENDIF
               ENDIF

               Q2(IX)=Q1(IX,IY,IZ)
5              CONTINUE
	       IBUF = IBUF+1
	       CALL WRTLIN(LUN2,Q2,NSAM,IBUF)
	   ENDDO
	ENDDO
	END
