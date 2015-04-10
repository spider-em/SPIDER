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
C                                                                      *
C IMAGE_PROCESSING_ROUTINE                                             *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE  OUTRNG
     &  (CIRC1,CIRC2,CIRTMP,LCIRC,NRING,NUMR,TOT,MAXRIN,NIMA)

        DIMENSION  CIRC1(LCIRC),CIRC2(LCIRC),CIRTMP(LCIRC)
        INTEGER    NUMR(3,NRING),MAXRIN
        COMPLEX  C

        PI2=8.0D0*DATAN(1.0D0)
        QI=1./REAL(NIMA-1)
        DO    I=1,NRING
           NSIRT=NUMR(3,I)
           CIRTMP(NUMR(2,I))=
     &          (CIRC1(NUMR(2,I))*NIMA-CIRC2(NUMR(2,I)))*QI
           CIRTMP(NUMR(2,I)+1)=
     &          (CIRC1(NUMR(2,I)+1)*NIMA-CIRC2(NUMR(2,I)+1)*
     &           COS(PI2*(TOT-1.0)/2.0*REAL(NSIRT)/REAL(MAXRIN)))*QI
     
           DO    J=3,NSIRT,2
              ARG=PI2*(TOT-1.0)*REAL(J/2)/REAL(MAXRIN)
              C=CMPLX(CIRC2(NUMR(2,I)+J-1),CIRC2(NUMR(2,I)+J))*
     &                 CMPLX(COS(ARG),SIN(ARG))
              CIRTMP(NUMR(2,I)+J-1)=
     &                (CIRC1(NUMR(2,I)+J-1)*NIMA-REAL(C))*QI
              CIRTMP(NUMR(2,I)+J)=
     &                (CIRC1(NUMR(2,I)+J)*NIMA-AIMAG(C))*QI
           ENDDO
        ENDDO
        END
