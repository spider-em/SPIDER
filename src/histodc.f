C++*********************************************************************
C
C HISTODC.F
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
C
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE  HISTODC(XR,N,LENH)

        INCLUDE 'CMBLOCK.INC'

        DIMENSION  XR(N),H1(LENH),H2(LENH)

C       GENERATE GAUSSIAN DISTRIBUTION - THIS WILL BE EXPECTED DISTRIBUTION
C       OF PIXEL VALUES AFTER THE TRANSFORMATION.

        SPI=0.0
        DO    I=1,LENH
           X=6.0*FLOAT(I-1)/FLOAT(LENH-1)-3.0
           H1(I)=EXP(-X*X/2.0)
           SPI=SPI+H1(I)
        ENDDO
        DO    I=1,LENH
           H1(I)=H1(I)/SPI
        ENDDO

        XRMI=XR(1)
        XRMA=XRMI
        DO    I=2,N
           XRMI=AMIN1(XRMI,XR(I))
           XRMA=AMAX1(XRMA,XR(I))
        ENDDO
C   B
        XCEN=(XRMA-XRMI)/2.0+XRMI
        TOL=0.001
        CHIMIN=GOLDEN(XRMI,XCEN,XRMA,TOL,B,N,LENH,XR,H1,H2,A,C)

         WRITE(NOUT,206)  A,B,C,CHIMIN
206      FORMAT(' The transformation is  A*(ALOG10(x+B)+C)',/,
     &   ' Parameters found     A =',1pe12.5,'   B =',1pe12.5,
     &   '   C =',1pe12.5,/,
     &   ' Chi-square     =',1pe12.5)

C        IF(NSEL(1).NE.0)  PARAM(NSEL(1))=A
C        IF(NSEL(2).NE.0)  PARAM(NSEL(2))=B
C        IF(NSEL(3).NE.0)  PARAM(NSEL(3))=C
         CALL REG_SET_NSEL(1,3,A,B,C,0.0,0.0,IRTFLG)

         END



      FUNCTION GOLDEN(AX,BX,CX,TOL,XMIN,N,LENH,XI,H1,H2,A,AVR)

      DIMENSION  XI(N),H1(LENH),H2(LENH)
      REAL GOLDEN,AX,BX,CX,TOL,XMIN,F,R,C
      PARAMETER (R=.61803399,C=1.-R)
      REAL F1,F2,X0,X1,X2,X3

      X0=AX
      X3=CX
      IF(ABS(CX-BX).GT.ABS(BX-AX))THEN
         X1=BX
         X2=BX+C*(CX-BX)
      ELSE
         X2=BX
         X1=BX-C*(BX-AX)
      ENDIF

      F1=FHODT(X1,N,LENH,XI,H1,H2,A,AVR)
      F2=FHODT(X2,N,LENH,XI,H1,H2,A,AVR)
1     IF (ABS(X3-X0).GT.TOL*(ABS(X1)+ABS(X2))) THEN
         IF (F2.LT.F1) THEN
            X0=X1
            X1=X2
            X2=R*X1+C*X3
            F1=F2
            F2=FHODT(X2,N,LENH,XI,H1,H2,A,AVR)
         ELSE
            X3=X2
            X2=X1
            X1=R*X2+C*X0
            F2=F1
            F1=FHODT(X1,N,LENH,XI,H1,H2,A,AVR)
         ENDIF
         GOTO 1
      ENDIF
      IF(F1.LT.F2)THEN
         GOLDEN=F1
         XMIN=X1
      ELSE
         GOLDEN=F2
         XMIN=X2
      ENDIF
      RETURN
      END
