C++*********************************************************************
C
C   ARITHSCA.F   IRTFLG SUPPORT                  1/30/20  ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2020  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C   ARITHSCA(LUN1,LUN2,NX,NY,NZ,FMIN,FMAX,IRTFLG)
C
C   PURPOSE:  SCALES IMAGE PIXEL BY PIXEL
C
C   PARAMETERS:
C        LUN1         I/O UNIT NUMBER OF FILE 1
C        LUN2         I/O UNIT NUMBER OF FILE 2
C        NX,NY,NZ     DIMENSIONS OF FILES
C        IRTFLG       ERROR FLAG
C
C--*******************************************************************

      SUBROUTINE ARITHSCA(LUN1,LUN2,NX,NY,NZ,
     &                   FMINT,FMAXT,FLOW,FHI,IRTFLG)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      REAL        :: BUF
      COMMON /IOBUF/ BUF(NBUFSIZ)      ! NBUFSIZ FROM CMLIMITS.INC
  
      INTEGER     :: LUN1,LUN2,NX,NY,NZ,IRTFLG
      REAL        :: FMINT,FMAXT,FLOW,FHI

      INTEGER     :: ISL,IOFF,I,IROW,K,NDUM
      REAL        :: CON2,CON,RANGEOLD,RANGENEW

      IF (FMAXT == FMINT) THEN
C        MAYBE YOU SHOULD HANDLE BLANK IMAGE IN CALLER NOW
         CALL ERRT(101,'BLANK FILE SKIPPED',NDUM)
         IRTFLG = -99
         RETURN
      ENDIF
  
      RANGENEW  = FHI   - FLOW
      RANGEOLD  = FMAXT - FMINT
      CON       = RANGENEW / RANGEOLD
      CON2      = FLOW - CON * FMINT

      DO  ISL=1,NZ
        IOFF = (ISL-1) * NY
        DO  I = 1,NY
           IROW = IOFF + I
           CALL REDLIN(LUN1,BUF,NX,IROW)
           DO  K = 1,NX
C             BUF(K) = FLOW + (FHI - FLOW) * (BUF(K) - FMINT) / (FMAXT - FMINT)
C             BUF(K) = FLOW + RANGENEW * (BUF(K) - FMINT) / RANGEOLD
C             BUF(K) = FLOW + CON * (BUF(K) - FMINT)
C             BUF(K) = FLOW + CON * BUF(K) - CON * FMINT
              BUF(K) = CON2 + CON * BUF(K)
	   ENDDO
           CALL WRTLIN(LUN2,BUF,NX,IROW)
        ENDDO
      ENDDO

      IRTFLG = 0
      END
