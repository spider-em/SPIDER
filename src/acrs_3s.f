C++*********************************************************************
C
C ACRS_3
C            MOD PGI COMPILER BUG              FEB 19 2008 ArDean Leith
C            X RETURNS O                       APR 24 2009 ArDean Leith
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
C  Calculates circular autocorrelation, non-power-of-two dimensions
C  Input - X Fourier transform
C  Output -  O=F(X*conjg(X))
C
C  PARAMETERS: X    FOURIER TRANSFORM                          (SENT)
C                      X=F(X )                                 (RET.)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE ACRS_3S(X, LS,NSAM,NROW,NSLICE)

        COMPLEX           :: X(LS/2,NROW,NSLICE)
        COMPLEX           :: O(LS/2,NROW,NSLICE)

        NNNN = LS / 2

        PI2  = 8.0 * DATAN(1.0D0)
        ITMP = NSAM / 2
        SX   = PI2 * FLOAT(ITMP) / FLOAT(NSAM)
        ITMP = NROW / 2
        SY   = PI2 * FLOAT(ITMP) / FLOAT(NROW)
        ITMP = NSLICE / 2
        SZ   = PI2 * FLOAT(ITMP) / FLOAT(NSLICE)

c$omp   parallel do private(i,j,k,ix,iy,iz,arg,argy,argz)
        DO K=1,NSLICE
           IZ = K-1
           IF (IZ .GT. NSLICE/2)  IZ = IZ-NSLICE
           ARGZ = SZ * IZ

           DO J=1,NROW
              IY = J-1
              IF (IY .GT. (NROW/2))  IY = IY - NROW
              ARGY = SY * IY + ARGZ
              DO I=1,NNNN
                 IX       = I - 1
                 ARG      = SX * IX + ARGY
                 X(I,J,K) = CABS(X(I,J,K)) * CMPLX(COS(ARG),SIN(ARG))
              ENDDO 
           ENDDO
        ENDDO

        INS = -1
        CALL FMRS_3(X,NSAM,NROW,NSLICE,INS)

        END
