C++*********************************************************************
C
C   ARITHSCA.F 
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
C   ARITH(LUN1,LUN2,NSAM,NROW,NSLICE,FMIN,FMAX)
C
C   PURPOSE:  SCALES IMAGE PIXEL BY PIXEL
C
C   PARAMETERS:
C        LUN1         LOGICAL UNIT NUMBER OF FILE 1
C        LUN2         LOGICAL UNIT NUMBER OF FILE 2
C        NSAM,NROW    X & Y DIMENSIONS OF FILES
C        NSLICE       Z DIMENSION OF FILES
C
C--*******************************************************************

      SUBROUTINE ARITHSCA(LUN1,LUN2,NSAM,NROW,NSLICE,
     &                   FMINT,FMAXT,FLOW,FHI)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      COMMON /IOBUF/ BUF(NBUFSIZ)
  
      IF (FMAXT .EQ. FMINT) THEN
         CALL ERRT(101,'BLANK FILE SKIPPED',NE)
         RETURN
      ENDIF
  
      RANGENEW  = FHI   - FLOW
      RANGEOLD  = FMAXT - FMINT
      CON       = RANGENEW / RANGEOLD
      CON2      = FLOW - CON * FMINT

      DO  ISL=1,NSLICE
        IOFF = (ISL-1) * NROW
        DO  I = 1,NROW
           IROW = IOFF + I
           CALL REDLIN(LUN1,BUF,NSAM,IROW)
           DO  K = 1,NSAM
C             BUF(K) = FLOW + (FHI - FLOW) * (BUF(K) - FMINT) / (FMAXT - FMINT)
C             BUF(K) = FLOW + RANGENEW * (BUF(K) - FMINT) / RANGEOLD
C             BUF(K) = FLOW + CON * (BUF(K) - FMINT)
C             BUF(K) = FLOW + CON * BUF(K) - CON * FMINT
              BUF(K) = CON2 + CON * BUF(K)
	   ENDDO
           CALL WRTLIN(LUN2,BUF,NSAM,IROW)
        ENDDO
      ENDDO
      RETURN
      END
