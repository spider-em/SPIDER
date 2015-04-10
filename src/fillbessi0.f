 
C++*********************************************************************
C 
C  FILLBESSI0.F        NEW                        MAR  07  ARDEAN LEITH     3/22/07
C                      FILLBESSIL                 DEC  08  ArDean Leith
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
C PURPOSE:   THIS ROUTINE GIVES THE REGULAR MODIFIED CYLINDRICAL 
C            BESSEL FUNCTION OF ZERO'TH ORDER 
C
C  PARAMETERS:
C        NSAM      SIZE                                         (SENT)
C        LTABI     ARRAY SIZE                                   (SENT)
C        LNB       -3                                           (RET.)
C        LNE       3                                            (RET.)
C        FLTB      LTABI/ 1.25  / 3                             (RET.)
C        ALPHA     1.75                                         (RET.)
C        RRR       NSAM/2                                       (RET.)
C        V         3/(2*NSAM)                                   (RET.)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12 
C--********************************************************************* 

        SUBROUTINE FILLBESSI0(NSAM,LTABI,LNB,LNE,FLTB,TABI,ALPHA,RRR,V)
 
        REAL            :: TABI(0:LTABI)

	REAL, PARAMETER  :: QUADPI = 3.1415926535897932384626D0
	REAL, PARAMETER  :: TWOPI = 2*QUADPI 

C       KAISER-BESSEL WINDOW ACCORDING TO SCHOMBERG

        LN    = 6
        LNB   = -INT(LN/2)                ! ALWAYS: -3
	LNE   = INT(LN/2)                 ! ALWAYS   3

        NPAD  = 2 * NSAM
	LTAB  = NINT(REAL(LTABI) / 1.25)
	RRR   = NSAM / 2
	ALPHA = 1.75                      ! CONSTANT
	V     = REAL(LN) / 2.0 / REAL(NPAD)

C       ADJUST V TO MAKE SURE THAT IT IS NOT ZERO AT THE WINDOW BORDER
	VADJUST = 1.0 * V

C       GENERATE TABLE WITH INTERPOLANTS
	IF (LTABI .GT. LTAB)  TABI(LTAB+1:LTABI) = 0.0

	B0   = BESI0(TWOPI * ALPHA * RRR * VADJUST)
        FLTB = REAL(LTAB) / REAL(LNE)

cc$omp  parallel do private(i,s,xx)
        DO I=0,LTAB
	   S = REAL(I) / FLTB / NPAD
	   IF (S .LT. VADJUST)  THEN
	      XX      = SQRT(1.0 - (S / VADJUST)**2)
	      TABI(I) = BESI0(TWOPI * ALPHA * RRR * VADJUST * XX) / B0
	   ELSE
	      TABI(I) = 0.0
	   ENDIF
        ENDDO

        END

C       -------------------- FILLBESSIL -----------------------------

        SUBROUTINE FILLBESSIL(N2,LN2,LTAB,FLTB,TABI,ALPHA,AAAA,NNN)

        !! COMMON  /TABS/ LN2,FLTB,TABI(0:LTAB)
 
        REAL      :: TABI(0:LTAB)

C       GENERALIZED KAISER-BESSEL WINDOW ACCORDING TO LEWITT
        LN    = 5                           ! ALWAYS=5
        LN2   = LN / 2                      ! ALWAYS=2
	V     = REAL(LN-1) / 2.0 / REAL(N2) ! ALWAYS=4*N2
	ALPHA = 6.5                         ! ALWAYS=6.5
	AAAA  = 0.9*V                       ! ALWAYS=.9*4*N2                   
	NNN   = 3                           ! ALWAYS=3

C       GENERATE TABLE WITH INTERPOLANTS
 	B0   = SQRT(ALPHA) * BESI1(ALPHA)

        FLTB = REAL(LTAB) / REAL(LN2+1)

cc$omp  parallel do private(i,s,xt)
        DO I=0,LTAB
	   S = REAL(I) / FLTB / N2
	   IF (S .LE. AAAA)  THEN
	      XT      = SQRT(1.0 - (S/AAAA)**2)
	      TABI(I) = SQRT(ALPHA*XT) * BESI1(ALPHA*XT) / B0 
	   ELSE
	      TABI(I) = 0.0
	   ENDIF
        ENDDO

        END
