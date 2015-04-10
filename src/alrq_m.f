
C ++********************************************************************
C                                                                      *
C                                                                      *
C                                                                      *
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
C                                                                      *
C                                                                      *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C
C IMAGE_PROCESSING_ROUTINE                                             *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE  ALRQ_M
     &  (XIM,NSAM,NROW,CNS2,CNR2,NUMR,CIRC,LCIRC,NRING,MODE)

        DIMENSION  XIM(NSAM,NROW),CIRC(LCIRC)
        INTEGER  NUMR(3,NRING)
        CHARACTER*1  MODE
        DOUBLE PRECISION  PI,DFI

C
C  INTERPOLATION INTO POLAR COORDINATES
C
C        CNS2 and CNR2 are predefined centers
CC no need to set to zero, all elements are defined
Cc$omp parallel do private(i)
C        DO  10  I=1,LCIRC
C 10     CIRC(I)=0.0

        PI=2*DATAN(1.0D0)
c$omp parallel do private(i,j,inr,yq,l,lt,nsim,dfi,kcirc,
c$omp& xold,yold,fi,x,y)
        DO  I=1,NRING

C  RADIUS OF THE RING
           INR=NUMR(1,I)
           YQ=INR

           L=NUMR(3,I)
           IF(MODE.EQ.'H')  THEN
              LT=L/2
           ENDIF
           IF(MODE.EQ.'F')  THEN
              LT=L/4
           ENDIF
           NSIM=LT-1
           DFI=PI/(NSIM+1)
           KCIRC=NUMR(2,I)
           XOLD=0.0
           YOLD=INR
        CIRC(KCIRC)=QUADRI(XOLD+CNS2,YOLD+CNR2,NSAM,NROW,XIM)
           XOLD=INR
           YOLD=0.0
        CIRC(LT+KCIRC)=QUADRI(XOLD+CNS2,YOLD+CNR2,NSAM,NROW,XIM)
           IF(MODE.EQ.'F')  THEN
              XOLD=0.0
              YOLD=-INR
        CIRC(LT+LT+KCIRC)=QUADRI(XOLD+CNS2,YOLD+CNR2,NSAM,NROW,XIM)
              XOLD=-INR
              YOLD=0.0
        CIRC(LT+LT+LT+KCIRC)=QUADRI(XOLD+CNS2,YOLD+CNR2,NSAM,NROW,XIM)
           ENDIF
           DO   J=1,NSIM
              FI=DFI*J
              X=SIN(FI)*YQ
              Y=COS(FI)*YQ

              XOLD=X
              YOLD=Y
        CIRC(J+KCIRC)=QUADRI(XOLD+CNS2,YOLD+CNR2,NSAM,NROW,XIM)
              XOLD=Y
              YOLD=-X
        CIRC(J+LT+KCIRC)=QUADRI(XOLD+CNS2,YOLD+CNR2,NSAM,NROW,XIM)
              IF(MODE.EQ.'F')  THEN
                 XOLD=-X
                 YOLD=-Y
        CIRC(J+LT+LT+KCIRC)=QUADRI(XOLD+CNS2,YOLD+CNR2,NSAM,NROW,XIM)
                 XOLD=-Y
                 YOLD=X
        CIRC(J+LT+LT+LT+KCIRC)=QUADRI(XOLD+CNS2,YOLD+CNR2,NSAM,NROW,XIM)
              ENDIF
	   ENDDO
	ENDDO
c 20    CONTINUE
C
        END
