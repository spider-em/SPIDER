C++*********************************************************************
C
C  CNRD_3.F 
C                PGI BUG                       FEB 10 2006 ArDean Leith
C                MOD PGI COMPILER BUG          FEB 19 2008 ArDean Leith
C                X RETURNS O                   APR 24 2009 ArDean Leith
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
C  PURPOSE: CALCULATES CONVOLUTION
C
C  PARAMETERS: X  FOURIER TRANSFORM                           (SENT)
C                    X=F(X*Y)                                 (RET.)
C              BUF  IO  
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE CNRD_3(LUN2,X,BUF, LS,NSAM,NROW,NSLICE)

        COMPLEX          :: X(LS/2,NROW,NSLICE)
        COMPLEX          :: BUF(LS/2)

        DOUBLE PRECISION :: PI2
        COMPLEX          :: CTEMP
   
        NNNN = LS / 2
        LRS  = 2 * NNNN

        PI2  = 8.0 * DATAN(1.0D0)
        ITMP = NSAM / 2
        SX   = PI2 * FLOAT(ITMP) / FLOAT(NSAM)
        ITMP = NROW / 2
        SY   = PI2 * FLOAT(ITMP) / FLOAT(NROW)
        ITMP = NSLICE / 2
        SZ   = PI2 * FLOAT(ITMP) / FLOAT(NSLICE)

        NR = 0

        DO K=1,NSLICE
              IZ=K-1
              IF (IZ .GT. (NSLICE/2))  IZ = IZ - NSLICE
              ARGZ = SZ * IZ
              DO J=1,NROW
                 IY = J-1
                 IF (IY .GT. (NROW/2))  IY = IY - NROW
                 ARGY = SY * IY + ARGZ
                 NR   = NR + 1
                 CALL  REDLIN(LUN2,BUF,LRS,NR)

                 DO I=1,NNNN
                    IX       = I - 1
                    ARG      = SX * IX + ARGY
                    CTEMP    = CMPLX(COS(ARG),SIN(ARG))
                    x(I,J,K) = X(I,J,K) * BUF(I) * CTEMP
                 ENDDO
              ENDDO
        ENDDO

        INS = -1
        CALL FMRS_3(X,NSAM,NROW,NSLICE,INS)

        END
