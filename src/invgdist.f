
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
C	 RETURNS A NORMALLY DISTRIBUTED DEVIATE WITH ZERO MEAN AND UNIT
C	 VARIANCE USING UNIFORM_DEVIATE AS THE SOURCE OF THE UNIFORM
C	 DEVIATE FOR THE GENERATION OF THE GAUSSIAN
C
C	 UNIFORM_DEVIATE IS ASSUMED TO BE BETWEEN 0.0 AND 1.0
C
C	 BASED ON ABRAMOWITZ AND STEGEN P 953 AND 933
C
C	 SDF 27-SEP-89
C
C........................................................................
 	
        SUBROUTINE INVGDIST(UNIFORM_DEVIATE, GAUSSIAN_DEVIATE)

        INCLUDE 'CMBLOCK.INC'

        DOUBLE PRECISION C0,C1,C2,D1,D2,D3
 	PARAMETER (C0 = 2.515517)
 	PARAMETER (C1 = 0.802853)
 	PARAMETER (C2 = 0.010328)
 	PARAMETER (D1 = 1.432788)
 	PARAMETER (D2 = 0.189269)
 	PARAMETER (D3 = 0.001308)
 
 	REAL*4 UNIFORM_DEVIATE, GAUSSIAN_DEVIATE
 	DOUBLE PRECISION  T, A, B
 
 	IF (UNIFORM_DEVIATE .GE. 0.0
     1      .AND. UNIFORM_DEVIATE .LE. 1.0) THEN	    
 	    T = DSQRT(-2.0*DLOG(UNIFORM_DEVIATE+1.0D-20))
 	    A = C0 + T*(C1 + C2*T)
 	    B = 1.0 + T*(D1 + T*(D2 + D3*T))
 	    GAUSSIAN_DEVIATE = T-(A/B)
 	ELSE
 	   WRITE(NOUT,*)
     1    '('' INVGDIST: ERROR - UNIFORM_DEVIATE OUT OF RANGE  '',G)',
     1	     UNIFORM_DEVIATE 
 	    CALL ERRT(100,'INVGDIST',NE)
 	END IF
 	
 	RETURN
 	END
