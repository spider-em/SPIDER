C++*********************************************************************
C
C $$ FFTR_D.FOR
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
C IMAGE_PROCESSING_ROUTINE
C
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************
C
C $$ FFTR_D.FOR
C
         SUBROUTINE  FFTR_D(X,NV)
         DOUBLE PRECISION  X(2,1)
         DOUBLE PRECISION  RNI,TR1,TR2,TI1,TI2,TR,TI
         DOUBLE PRECISION  CC,C,SS,S,T
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

C
         NU=IABS(NV)
         INV=NV/NU
         NU1=NU-1
         N=2**NU1
         ISUB=16-NU1
         SS=-TAB1(ISUB)
         CC=-2.0*TAB1(ISUB-1)**2
         C=1.0
         S=0.0
         N2=N/2
         IF(INV.GT.0)  THEN
            CALL  FFTC_D(X(1,1),X(2,1),NU1,2)
            TR=X(1,1)
            TI=X(2,1)
            X(1,1)=TR+TI
            X(2,1)=TR-TI
            DO    I=1,N2
               I1=I+1
               I2=N-I+1
               TR1=X(1,I1)
               TR2=X(1,I2)
               TI1=X(2,I1)
               TI2=X(2,I2)
               T=(CC*C-SS*S)+C
               S=(CC*S+SS*C)+S
               C=T
               X(1,I1)=0.5*((TR1+TR2)+(TI1+TI2)*C-(TR1-TR2)*S)
               X(1,I2)=0.5*((TR1+TR2)-(TI1+TI2)*C+(TR1-TR2)*S)
               X(2,I1)=0.5*((TI1-TI2)-(TI1+TI2)*S-(TR1-TR2)*C)
               X(2,I2)=0.5*(-(TI1-TI2)-(TI1+TI2)*S-(TR1-TR2)*C)
            ENDDO
         ELSE
            TR=X(1,1)
            TI=X(2,1)
            X(1,1)=0.5*(TR+TI)
            X(2,1)=0.5*(TR-TI)
            DO    I=1,N2
               I1=I+1
               I2=N-I+1
               TR1=X(1,I1)
               TR2=X(1,I2)
               TI1=X(2,I1)
               TI2=X(2,I2)
               T=(CC*C-SS*S)+C
               S=(CC*S+SS*C)+S
               C=T
               X(1,I1)=0.5*((TR1+TR2)-(TR1-TR2)*S-(TI1+TI2)*C)
               X(1,I2)=0.5*((TR1+TR2)+(TR1-TR2)*S+(TI1+TI2)*C)
               X(2,I1)=0.5*((TI1-TI2)+(TR1-TR2)*C-(TI1+TI2)*S)
               X(2,I2)=0.5*(-(TI1-TI2)+(TR1-TR2)*C-(TI1+TI2)*S)
            ENDDO
            CALL  FFTC_D(X(1,1),X(2,1),NU1,-2)
         ENDIF
         END
