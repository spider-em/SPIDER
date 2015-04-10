C++*********************************************************************
C
C $$ FHODT.FOR
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
C
C
C--*********************************************************************
C
C $$ FHODT.FOR
C
        FUNCTION  FHODT(B,N,LENH,XI,H1,H2,A,AVR)
        DIMENSION  XI(N),H1(LENH),H2(LENH)

C       AVR is "-average"

        AVR=0.0
        DO    I=1,N
           AVR=AVR-ALOG10(XI(I)+B)
        ENDDO
        AVR=AVR/N
        SVR=0.0
        DO    I=1,N
           Q=ALOG10(XI(I)+B)+AVR
           SVR=SVR+Q*Q
        ENDDO
        A=SQRT(N/SVR)

        DO    I=1,LENH
           H2(I)=0.0
        ENDDO
        DO    I=1,N
           Q=A*(ALOG10(XI(I)+B)+AVR)
           L=MAX0(1,MIN0(LENH,INT((Q+3.0)/6.0*FLOAT(LENH-1))+1))
           H2(L)=H2(L)+1
        ENDDO

        FHT=0.0
        DO    I=1,LENH
           IF(H2(I).NE.0.0)  FHT=FHT+(H1(I)-H2(I)/N)**2
        ENDDO
        FHODT=FHT
        END
