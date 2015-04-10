        SUBROUTINE RMIWF(A,NSAM,NSAMF,C,PST,PINC,NPRO,IFIL)
      

C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2000  M. Radermacher                                  *
C=*                                                                    *
C=* Email:  spider@wadsworth.org                                       *
C=*                                                                    *
C=* This program is free software; you can redistribute it and/or      *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* This program is distributed in the hope that it will be useful,    *
C=* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C=* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
C=* General Public License for more details.                           *
C=*                                                                    *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program; if not, write to the                      *
C=* Free Software Foundation, Inc.,                                    *
C=* 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.      *
C=*                                                                    *
C=**********************************************************************

C       A        ARRAY FOR WEIGHTING FUNCTION
C       NSAM     REAL SPACE DIMENSION (USED TO CALCULATE FOURIER UNITS)
C       NSAMF    FOURIER DIMENSION, NORMALLY NSAM+2
C       DIAMETER DIAMETER OF RECONSTRUCTED OBJECT
C       PST      START ANGLE (IN RADIANS)
C       PINC     ANGULAR INCEMENT
C       NPRO     NUMBER OF PROJECTIONS

        DIMENSION A(*)
        DOUBLE PRECISION DINC,DARG
      
        UNIT=NSAM
        NSAMFH=NSAMF/2
        IPIC=IFIL
        DINC=1./FLOAT(NSAMFH-1)
C       A(1)=1
        IF(NPRO.NE.0) THEN
           A(1)=1./FLOAT(NPRO)
        ELSE
           A(1)=1.
        ENDIF
        GOTO (11,12,13,14) IPIC
11      DO 1 I=2,NSAMFH
1       A(I)=1
        GOTO 100
12      DO 2 I=2,NSAMFH
           DARG=DFLOAT(I-1)*DINC
2       A(I)=DSQRT(DARG)
        GOTO 100
13      DO 3 I=2,NSAMFH
3       A(I)=DFLOAT(I-1)*DINC
        GOTO 100
14      DO 4 I=1,NSAMFH
4       A(I)=1.      
100     RETURN
        END
