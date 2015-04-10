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
C  PURPOSE:   EVALUATES CONFIDENCE INTERVALS AT EACH IMAGE POINT (X,Y)
C
C N:            NUMBER OF IMAGE FILES ADDED
C ALPHA:        PROBABILITY OF ERROR
C LUN1:         LUN CONTAINING AVERAGE FILE
C LUN2:         LUN CONTAINING VARIANCE FILE
C LUN3:         LUN OF FILE TO RECEIVE THE UPPER CONFIDENCE IMAGE
C LUN4:         LUN OF FILE TO RECEIVE THE LOWER CONFIDENCE IMAGE
C
C WRITTEN BY W.HAENICKE, MPI FUER BIOPHYSIKALISCHE CHEMIE, GOETTINGEN.
C
C  PARAMETERS:                                                         *
C 
C SUPPORT_ROUTINE 
C                                                                    *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE CONF(N,ALPHA,LUN1,LUN2,LUN3,LUN4,NSAM,NROW)

C        COMMON ADUM(80),BUF(1)

        DIMENSION AIMG(NSAM), BIMG(NSAM)
        

        AN1=1./SQRT(FLOAT(N))
        PHIX=1.-ALPHA*0.005
        CALL NORPPF(PHIX,ZALPHA,IER)
        IF(IER.EQ.1) ZALPHA=4.

        DO  I=1,NROW
          CALL REDLIN(LUN1,AIMG,NSAM,I)
           CALL REDLIN(LUN2,BIMG,NSAM,I)
           DO  J=1,NSAM
              AM=AIMG(J)
              BM=ZALPHA*AN1*SQRT(BIMG(J))
              AIMG(J)=AM+BM
              BIMG(J)=AM-BM
           ENDDO
           CALL WRTLIN(LUN3,AIMG,NSAM,I)
           CALL WRTLIN(LUN4,BIMG,NSAM,I)
        ENDDO
        END


