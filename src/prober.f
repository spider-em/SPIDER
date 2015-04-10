
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
C       FROM MPI-5/20/80
C       CALL PROBER(A,N,LUN2,LUN3,NSAM,NROW)
C
C
C  PROBER EVALUATES THE PROBABILITY OF THE ERROR AT EVERY
C  POINT (X,Y), THE MAGNITUDE OF THE CONFIDENCE INTERVAL
C  BEING GIVEN
C
C
C           A:  THE MAGNITUDE OF THE CONFIDENCE INTERVAL
C           N:  NUMBER OF SAMPLE FILES ADDED
C        LUN2:  LUN OF V (VARIANCE) FILE
C        LUN3:  LUN OF FILE CONTAINING THE PROBABILITY OF THE ERROR
C
C SUPPORT_ROUTINE
C
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C--***************************************************************

        SUBROUTINE PROBER(A,N,LUN2,LUN3,NSAM,NROW,NZERO)

        DIMENSION AIMG(NSAM)

        NZERO=0
        AN1=0.70710678*A*SQRT(FLOAT(N))
        DO  I=1,NROW
           CALL REDLIN(LUN2,AIMG,NSAM,I)
           DO  J=1,NSAM
              IF(AIMG(J).LE.0) THEN
                 AIMG(J)=1.0
                 NZERO=NZERO+1
              ELSE
                 ZALPHA=AN1/SQRT(AIMG(J))
                 AIMG(J)=ORCDF(ZALPHA)
              ENDIF
           ENDDO
           CALL WRTLIN(LUN3,AIMG,NSAM,I)
        ENDDO
        END
