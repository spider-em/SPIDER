C
C++*********************************************************************
C
C NORMASC.F
C                     SPEED-UP                    OCT 2001 ARDEAN LEITH
C                     SQRT NEGATIVE TRAP          MAR 2008 ARDEAN LEITH
C                     REWRITE                     MAR 2008 ARDEAN LEITH
C **********************************************************************
C=* FROM: SPIDER - MODULAR IMAGE PROCESSING SYSTEM.   AUTHOR: J.FRANK  *
C=* Copyright (C) 1985-2008  Health Research Inc.                      *
C=*                                                                    *
C=* HEALTH RESEARCH INCORPORATED (HRI),                                *   
C=* ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455.                   *
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
C **********************************************************************
C
C NORMASC(X,Y,NS1,NS2,NR1,NR2,IR1,IR2)
C
C PURPOSE: DETERMINE NORMALIZATION PARAMETERS: AVG & VARIANCE
C
C  NOTE : FOR OMP PARALLEL USE NORMAS INSTEAD al Sept 01
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE NORMASC(X,NS1,NS2,NR1,NR2,IR1,IR2, AVT,VRINV,USE_OMP)

        IMPLICIT NONE

        REAL, INTENT(IN)              :: X(NS1:NS2,NR1:NR2)
        INTEGER, INTENT(IN)           :: NS1,NS2,NR1,NR2,IR1,IR2
        DOUBLE PRECISION, INTENT(OUT) :: AVT,VRINV
        LOGICAL, INTENT(IN)           :: USE_OMP

        DOUBLE PRECISION              :: AVG, VR, DTEMP
        INTEGER                       :: N,I1SQ,I2SQ, J,JSQ,I,IR

        I1SQ = IR1 * IR1
        I2SQ = IR2 * IR2
        AVG  = 0.0
        VR   = 0.0
        N    = 0

        IF (USE_OMP) THEN
c$omp     parallel do private(j,jsq,i,ir),reduction(+:n,avg,vr)
          DO J=NR1,NR2
              JSQ = J * J
              DO I=NS1,NS2
                 IR = JSQ + I*I
                 IF (IR .GE. I1SQ .AND. IR .LE. I2SQ)  THEN
                    N   = N + 1
                    AVG = AVG + X(I,J)
                    VR  = VR  + X(I,J) * X(I,J)
                 ENDIF
              ENDDO
           ENDDO
c$omp     end parallel do

        ELSE
           DO J=NR1,NR2
              JSQ = J * J
              DO I=NS1,NS2
                 IR = JSQ + I*I
                 IF (IR .GE. I1SQ .AND. IR .LE. I2SQ)  THEN
                    N   = N + 1
                    AVG = AVG + X(I,J)
                    VR  = VR  + X(I,J) * X(I,J)
                 ENDIF
              ENDDO
           ENDDO
        ENDIF

        AVT   = AVG / N
        DTEMP = (VR - N * AVT * AVT)

C       MULTIPLICATION IS FASTER
        IF (DTEMP .GT. 0) THEN
           VR    = DSQRT(DTEMP / DBLE(N-1))
           VRINV = 1.0 / VR
        ELSE
C          TRAP FOR BLANK IMAGE AREA
           VRINV = 0
        ENDIF

        END
