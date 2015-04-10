C++*********************************************************************
C
C $$ FFTC_D.FOR
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
C
C $$ FFTC_D.FOR
C
         SUBROUTINE  FFTC_D(BR,BI,LN,KS)
         DOUBLE PRECISION  BR(1),BI(1)
         DOUBLE PRECISION  RNI,SGN,TR1,TR2,TI1,TI2
         DOUBLE PRECISION  CC,C,SS,S,T,X2,X3,X4,X5
         INTEGER B3,B4,B5,B6,B7,B56
         DOUBLE PRECISION  TAB1(15)
#ifdef SP_MP
        TAB1(1)=9.58737990959775D-05
        TAB1(2)=1.91747597310703D-04
        TAB1(3)=3.83495187571395D-04
        TAB1(4)=7.66990318742704D-04
        TAB1(5)=1.53398018628476D-03
        TAB1(6)=3.06795676296598D-03
        TAB1(7)=6.13588464915449D-03
        TAB1(8)=1.22715382857199D-02
        TAB1(9)=2.45412285229123D-02
        TAB1(10)=4.90676743274181D-02
        TAB1(11)=9.80171403295604D-02
        TAB1(12)=1.95090322016128D-01
        TAB1(13)=3.82683432365090D-01
        TAB1(14)=7.07106781186546D-01
        TAB1(15)=1.00000000000000D+00
#else
       DATA  TAB1/           9.58737990959775D-05, 1.91747597310703D-04,
     1 3.83495187571395D-04, 7.66990318742704D-04, 1.53398018628476D-03,
     2 3.06795676296598D-03, 6.13588464915449D-03, 1.22715382857199D-02,
     3 2.45412285229123D-02, 4.90676743274181D-02, 9.80171403295604D-02,
     4 1.95090322016128D-01, 3.82683432365090D-01, 7.07106781186546D-01,
     5   1.00000000000000D+00/
#endif


         N=2**LN
         K=IABS(KS)
         L=16-LN
         B3=N*K
         B6=B3
         B7=K
         IF(KS.GT.0)  THEN
            SGN=1.0
         ELSE
            SGN=-1.0
            RNI=1.0D0/FLOAT(N)
            J=1
            DO    I=1,N
               BR(J)=BR(J)*RNI
               BI(J)=BI(J)*RNI
               J=J+K
            ENDDO
         ENDIF

12       B6=B6/2
         B5=B6
         B4=2*B6
         B56=B5-B6

14       TR1=BR(B5+1)
         TI1=BI(B5+1)
         TR2=BR(B56+1)
         TI2=BI(B56+1)

         BR(B5+1)=TR2-TR1
         BI(B5+1)=TI2-TI1
         BR(B56+1)=TR1+TR2
         BI(B56+1)=TI1+TI2

         B5=B5+B4
         B56=B5-B6
         IF(B5.LE.B3)  GOTO  14
         IF(B6.EQ.B7)  GOTO  20

         B4=B7
         CC=2.0*TAB1(L)**2
         C=1.0-CC
         L=L+1
         SS=SGN*TAB1(L)
         S=SS

16       B5=B6+B4
         B4=2*B6
         B56=B5-B6

18       TR1=BR(B5+1)
         TI1=BI(B5+1)
         TR2=BR(B56+1)
         TI2=BI(B56+1)
         BR(B5+1)=C*(TR2-TR1)-S*(TI2-TI1)
         BI(B5+1)=S*(TR2-TR1)+C*(TI2-TI1)
         BR(B56+1)=TR1+TR2
         BI(B56+1)=TI1+TI2

         B5=B5+B4
         B56=B5-B6
         IF(B5.LE.B3)  GOTO  18
         B4=B5-B6
         B5=B4-B3
         C=-C
         B4=B6-B5
         IF(B5.LT.B4)  GOTO  16
         B4=B4+B7
         IF(B4.GE.B5)  GOTO  12

         T=C-CC*C-SS*S
         S=S+SS*C-CC*S
         C=T
         GOTO  16

20       IX0=B3/2
         B3=B3-B7
         B4=0
         B5=0
         B6=IX0
         IX1=0
         IF(B6.EQ.B7)  RETURN

22       B4=B3-B4
         B5=B3-B5
         X2=BR(B4+1)
         X3=BR(B5+1)
         X4=BI(B4+1)
         X5=BI(B5+1)
         BR(B4+1)=X3
         BR(B5+1)=X2
         BI(B4+1)=X5
         BI(B5+1)=X4
         IF(B6.LT.B4)  GOTO  22

24       B4=B4+B7
         B5=B6+B5
         X2=BR(B4+1)
         X3=BR(B5+1)
         X4=BI(B4+1)
         X5=BI(B5+1)
         BR(B4+1)=X3
         BR(B5+1)=X2
         BI(B4+1)=X5
         BI(B5+1)=X4
         IX0=B6

26       IX0=IX0/2
         IX1=IX1-IX0
         IF(IX1.GE.0)  GOTO   26

         IX0=2*IX0
         B4=B4+B7
         IX1=IX1+IX0
         B5=IX1
         IF(B5.GE.B4)  GOTO  22
         IF(B4.LT.B6)  GOTO  24

         END
 
