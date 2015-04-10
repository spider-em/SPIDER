C ++********************************************************************
C
C   NORMAS.F                ALTRIX SPECIAL CODE    FEB 2005 ARDEAN LEITH 
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
C  NORMAS(X,NS1,NS2,NR1,NR2,IR1,IR2)                                   *
C                                                                      *
C  PURPOSE:    NORMALIZES RING DATA.  COVERED AREA IS: IR1....IR2      *
C                                                                      *
C  PARAMETERS:                                                         *
C
C  NOTE   :    I THINK THIS IS FOR PARALLEL USE ONLY, BECAUSE NORMASS
C              IS QUICKER FOR NON_PARALLEL USE!! al Sept 01
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE NORMAS(X,NS1,NS2,NR1,NR2, IR1,IR2)

        DIMENSION  X(NS1:NS2,NR1:NR2)
        REAL*8     AV,VR

        I1SQ = IR1 * IR1
        I2SQ = IR2 * IR2

        AV   = 0.0
        VR   = 0.0
        N    = 0

#ifndef __ia64
c$omp   parallel do private(j,j2,i,ir),reduction(+:av,vr,n)
#endif
        DO J=NR1,NR2
           J2 = J*J
           DO I=NS1,NS2
              IR = J2 + I*I
              IF (IR .GE. I1SQ .AND. IR .LE. I2SQ) THEN
                 N  = N + 1
                 AV = AV + X(I,J)
                 VR = VR + X(I,J) * X(I,J)
              ENDIF
           ENDDO
        ENDDO

#ifndef __ia64
c$omp   end parallel do 
#endif

        AV = AV / N

C       MULTIPLICATION IS FASTER
        VR = 1.0 / (DSQRT((VR-N*AV*AV) / (N-1)))

#ifdef __ia64
C        ALTIX SEG FAULTS IF THIS IS A || LOOP (I'VE TRIED EVERYTHING)
C        FEB 2005 al

         X  = (X - AV) * VR
#else
c$omp   parallel do private(i,j)
        DO J=NR1,NR2
           DO I=NS1,NS2
              X(I,J) = (X(I,J) - AV ) * VR
           ENDDO
        ENDDO
c$omp   end parallel do 
#endif

        END

