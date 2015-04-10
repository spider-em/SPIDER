C++*********************************************************************
C
C    CCRS_3.F 
C      PGI BUG                                  FEB 10 2006 ArDean Leith
C      REAL CALCULATIONS                        FEB 11 2006 Bimal Rath
C      MOD PGI COMPILER BUG                     FEB 19 2008 ArDean Leith
C      CCRS_3_PH                                JUN 30 2011 ArDean Leith
C
C **********************************************************************
C=* FROM: SPIDER - MODULAR IMAGE PROCESSING SYSTEM.   AUTHOR: J.FRANK  *
C=* Copyright (C) 1985-2011  Health Research Inc.                      *
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
C  CCRS_PH_3(X,Y, LS,NSAM,NROW,NSLICE),IRTFLG
C
C  PURPOSE: CALCULATES CIRCULAR CROSCORRELATION, NON-POWER-OF-TWO 
C           DIMENSIONS
C  
C  PARAMETERS:   X    FOURIER TRANSFORMS                      (SENT/RET)
C                Y    FOURIER TRANSFORMS                      (SENT)
C                LS   NSAM+2-MOD(NSAM,2)                      (SENT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE CCRS_PH_3(X,Y, LS,NSAM,NROW,NSLICE,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'

        REAL             :: X(LS,NROW,NSLICE)
        REAL             :: Y(LS,NROW,NSLICE)
        INTEGER          :: LS,NSAM,NROW,NSLICE,IRTFLG

        REAL             :: SX,SY,ARGZ,ARGY,ARG,TMPR,TMPI,FNRM,SZ,BTM
        INTEGER          :: K,IZ,J,IY,I,INS,ITMP,NE

        LOGICAL          :: BADBTM 

	REAL, PARAMETER  :: QUADPI = 3.1415926535897932384
	REAL, PARAMETER  :: PI2    = 2*QUADPI


        ITMP   = NSAM / 2
        SX     = PI2 * FLOAT(ITMP)/FLOAT(NSAM)
        ITMP   = NROW / 2
        SY     = PI2 * FLOAT(ITMP )/ FLOAT(NROW)

        BADBTM = .FALSE.

c$omp  parallel do private(i,j,k,iy,iz,arg,argy,argz,tmpr,tmpi,fnrm,
c$omp&                     badbtm,btm)
        DO K=1,NSLICE
           IZ = K - 1
           IF (IZ .GT. (NSLICE/2))  IZ = IZ-NSLICE
           ARGZ = SZ * IZ

           DO J=1,NROW
              IY = J - 1
              IF (IY .GT. (NROW/2)) IY = IY - NROW
              ARGY = SY * IY

              DO I=1,LS,2
                 ARG        = SX * (I-1) / 2 + ARGY

    	         TMPR       = X(I,  J,K) * Y(I,  J,K) + 
     &                        X(I+1,J,K) * Y(I+1,J,K)

	         TMPI       = X(I+1,J,K) * Y(I,  J,K) - 
     &                        X(I,  J,K) * Y(I+1,J,K)

	         X(I,J,K)   = TMPR * COS(ARG) - TMPI * SIN(ARG)
	         X(I+1,J,K) = TMPI * COS(ARG) + TMPR * SIN(ARG)

                 BTM = SQRT(X(I,J,K)**2 + X(I+1,J,K)**2)
                 IF (BTM .LE. 0.0) THEN
                    BADBTM = .TRUE.
                    EXIT
                 ENDIF
                  
	         FNRM       = 1.0 / BTM

	         X(I,  J,K) = X(I,  J,K) * FNRM
	         X(I+1,J,K) = X(I+1,J,K) * FNRM
              ENDDO
           ENDDO
        ENDDO

        IF (BADBTM) THEN
           WRITE(NOUT,90) X(I,J,K),X(I+1,J,K)
90         FORMAT(' COMPLEX VALUES GIVE DIV. BY ZERO:',1PG12.3)
           CALL ERRT(101,'AVOIDED DIV BY ZERO',NE)
           IRTFLG = 1
           RETURN
        ENDIF
   
        INS    = -1
        CALL FMRS_3(X,NSAM,NROW,NSLICE,INS)

        IRTFLG = 0
        END


#ifdef NEVER
		 OREAL = TMPR * COS(ARG) - TMPI * SIN(ARG)
		 OIMG  = TMPI * COS(ARG) + TMPR * SIN(ARG)
		 
		 O(I,J,K)   = OREAL
		 O(I+1,J,K) = OIMG		 
#endif



 
