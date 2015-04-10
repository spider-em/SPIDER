C ++********************************************************************
C                                                                      *
C                                                                      *
C                           PGI BUG            FEB 10 2006 ArDean Leith
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
C  PURPOSE:                                                            *
C
C  CALCULATES CIRCULAR CROSCORRELATION 
C    BETWEEN 180 DEG ROTATED IMAGE AND 0 DEGS
C  NON-POWER-OF-TWO DIMENSIONS
C  INPUT - X,Y FOURIER TRANSFORMS
C  OUTPUT -  O=F(CONJG(X)*CONJG(Y))
C
C  PARAMETERS:                                                         *
C
C IMAGE_PROCESSING_ROUTINE                                                                      
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE  CR180_2(X,Y,O,NNNN,NSAM,NROW)

        COMPLEX  X(NNNN,NROW),Y(NNNN,NROW),O(NNNN,NROW)
        DOUBLE PRECISION  PI2
        COMPLEX :: CTEMP

        PI2=8.0*DATAN(1.0D0)
        ITMP=NSAM/2
        SX=PI2*FLOAT(ITMP)/FLOAT(NSAM)
        ITMP=NROW/2
        SY=PI2*FLOAT(ITMP)/FLOAT(NROW)

c$omp   parallel do private(i,j,ix,iy,argy,arg,ctemp)
        DO J=1,NROW
           IY = J-1
           IF (IY .GT. NROW/2)  IY=IY-NROW
           ARGY = SY * IY

           DO I=1,NNNN
              IX     = I-1
              ARG    = SX * IX + ARGY
              CTEMP  = CMPLX(COS(ARG),SIN(ARG))
              O(I,J) = CONJG(X(I,J)) * CONJG(Y(I,J)) * CTEMP
           ENDDO
        ENDDO

        INS = -1
        CALL  FMRS_2(O,NSAM,NROW,INS)

        END
