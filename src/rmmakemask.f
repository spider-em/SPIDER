
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


        SUBROUTINE RMMAKEMASK(AMASK,NSAM,NROW,C,RAD)
        
C       PROGRAM TO CALCULATE RADON TRANSFORM FOURIER FILTRATION AMASK FOR NTN
C       REDUCTION.
        DIMENSION AMASK(NSAM,NROW)
C       CALCULATE SOME PRELIMIMARY VALUES:
        PI=3.14159265
        NSAMP=NSAM-2
        NSAMPH=NSAMP/2
        RADMAX=RAD/NSAMPH
        CC=ABS(C)
        DENOM=PI*RADMAX*CC
C       THE ORIGIN LINE IN THE FT IS ALWAYS THERE
        DO I=1,NSAM
           AMASK(I,1)=1.
        ENDDO
        DO I=2,NROW/2+1
           BIGRAD=FLOAT(I)/NSAMPH
C          89.9 DEGREES CORRESPONDS TO 0.00174533
           IF(CC.GT.0.0001) THEN
              CUTOFF=BIGRAD/DENOM*NSAMP
           ELSE
              CUTOFF=NSAM
           ENDIF
           IY=NROW+2-I
           DO K=1,NSAM,2
              IP=(K-1)/2
              P=FLOAT(IP)
              IF(P.LT.CUTOFF) THEN
                 AMASK(K,I) = 0.
                 AMASK(K+1,I) = 0.
                 AMASK(K,IY) = 0.
                 AMASK(K+1,IY) = 0.
              ELSE
                 AMASK(K,I) = 1.
                 AMASK(K+1,I) = 1.
                 AMASK(K,IY) = 1.
                 AMASK(K+1,IY) = 1.
              ENDIF
           ENDDO
        ENDDO
        RETURN
        END
