C++*********************************************************************
C
C $$ CROSRMG.FOR
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
C IMAGE_PROCESSING_ROUTINE
C
C        1         2         3         4         5         6         7
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************
C
C $$ CROSRMG.FOR
C
         SUBROUTINE  CROSRMG
     &   (CIRC1,CIRC2,LCIRC,NRING,T,Q,MAXRIN,JACUP,NUMR,QN,TOT,MODE)
C
C  INPUT - Fourier transforms of rings!!!
C  First set is conjugated (mirrored)
C
         INTEGER  NUMR(3,NRING),MAXRIN
         DIMENSION  CIRC1(LCIRC),CIRC2(LCIRC)
         DOUBLE PRECISION  T(MAXRIN),Q(MAXRIN)
         DOUBLE PRECISION  QN,QT,PI,T7(-3:3)
         CHARACTER*1  MODE

         PI=4.0*DATAN(1.0D0)
         IF(MODE.EQ.'F')  PI=2*PI
         IP=-LOG2(MAXRIN)
C
         Q=0.0D0
C
         DO    I=1,NRING
            WR=REAL(NUMR(1,I))*PI/REAL(NUMR(3,I))
     &         *REAL(MAXRIN)/REAL(NUMR(3,I))
C
            T(1)=DBLE(CIRC1(NUMR(2,I)))*CIRC2(NUMR(2,I))
            IF(NUMR(3,I).EQ.MAXRIN)  THEN
               T(2)=DBLE(CIRC1(NUMR(2,I)+1))*CIRC2(NUMR(2,I)+1)
               DO    J=3,MAXRIN,2
                  JC=J+NUMR(2,I)-1
         T(J)=DBLE(CIRC1(JC))*CIRC2(JC)-DBLE(CIRC1(JC+1))*CIRC2(JC+1)
      T(J+1)=-DBLE(CIRC1(JC))*CIRC2(JC+1)-DBLE(CIRC1(JC+1))*CIRC2(JC)
	       ENDDO
               Q=Q+T*WR
            ELSE
         T(NUMR(3,I)+1)=DBLE(CIRC1(NUMR(2,I)+1))*CIRC2(NUMR(2,I)+1)/2.0
               T(2)=0.0
               DO    J=3,NUMR(3,I),2
                  JC=J+NUMR(2,I)-1
         T(J)=DBLE(CIRC1(JC))*CIRC2(JC)-DBLE(CIRC1(JC+1))*CIRC2(JC+1)
      T(J+1)=-DBLE(CIRC1(JC))*CIRC2(JC+1)-DBLE(CIRC1(JC+1))*CIRC2(JC)
	       ENDDO
C  ###
              Q(1:NUMR(3,I)+1)=Q(1:NUMR(3,I)+1)+T(1:NUMR(3,I)+1)*WR
           ENDIF
	 ENDDO
C
         CALL  FFTR_D(Q,IP)
C
         QN=-1.0D20
         DO    J=1,MAXRIN
            IF(Q(J).GE.QN)  THEN
               QN=Q(J)
               JTOT=J
            ENDIF
	 ENDDO
C
       IF(JACUP.EQ.0)  THEN
          TOT=JTOT
       ELSE
          DO    K=-3,3
              J=MOD(JTOT+K+MAXRIN-1,MAXRIN)+1
              T7(K)=Q(J)
	  ENDDO
          CALL  PRB1D(T7,7,POS)
          TOT=REAL(JTOT)+POS
          K=INT(100.0/REAL(JACUP+1))
          TOT=REAL(JTOT)
     &     +REAL(IFIX(POS))
     &      +REAL(K)*REAL(INT(POS*100.0/REAL(K)))/100.0
         ENDIF
         END
