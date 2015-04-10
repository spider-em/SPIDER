C++*********************************************************************
C
C $$ CROSRNG_Q.FOR
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
C $$ CROSRNG_Q.FOR
C
         SUBROUTINE  CROSRNG_Q
     &   (CIRC1,CIRC2,LCIRC,NRING,T,Q,MAXRIN,JACUP,NUMR,QN,TOT,MODE)
C
C  INPUT - Fourier transforms of rings!!!
C
C  NUMR(3,I) - ACTUAL LENGTH +2 (to use new FFT)
C  CIRC - FT of rings stored as:
C  Re(0), 0, Re(1) Im(1), ..., Re((NUMR(3,I)-2)/2), 0.
        INTEGER  NUMR(3,NRING),NUMR2I,NUMR3I,MAXRIN
        DIMENSION  CIRC1(LCIRC),CIRC2(LCIRC)
        DOUBLE PRECISION  T(MAXRIN+2),Q(MAXRIN+2)
        DOUBLE PRECISION  QN,QT,PI,T7(-3:3)
        CHARACTER*1  MODE

        PI=4.0*DATAN(1.0D0)
        IF(MODE.EQ.'F')  PI=2*PI

        Q=0.0D0

        DO    I=1,NRING
           NUMR3I=NUMR(3,I)-2
           NUMR2I=NUMR(2,I)
        WR=REAL(NUMR(1,I))*PI/REAL(NUMR3I)*REAL(MAXRIN)/REAL(NUMR3I)

           DO   J=1,NUMR3I+2,2
              JC=J+NUMR2I-1
        T(J)=DBLE(CIRC1(JC))*CIRC2(JC)+DBLE(CIRC1(JC+1))*CIRC2(JC+1)
        T(J+1)=-DBLE(CIRC1(JC))*CIRC2(JC+1)+DBLE(CIRC1(JC+1))*CIRC2(JC)
           ENDDO

           IF(NUMR3I.LT.MAXRIN)  T(NUMR3I+1)=T(NUMR3I+1)/2.0
C###
           Q(1:NUMR3I+1)=Q(1:NUMR3I+1)+T(1:NUMR3I+1)*WR
        ENDDO
C
        INV=-1
        IP=MAXRIN
        CALL  FMRS_1D(Q,IP,INV)
C
         QN=-1.0D20
         DO    J=1,MAXRIN
            IF(Q(J).GE.QN)  THEN
               QN=Q(J)
               JTOT=J
            ENDIF
         ENDDO
C
         DO  K=-3,3
            J=MOD(JTOT+K+MAXRIN-1,MAXRIN)+1
            T7(K)=Q(J)
         ENDDO
         CALL  PRB1D(T7,7,POS)
         TOT=REAL(JTOT)+POS
         END
