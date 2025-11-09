C **********************************************************************
C *  GETDEFOCUS.F
C *  THIS CODE IS A MESS!  SEGFAULTS and COMPILER ERRORS!! 
C *  I ABANDONED TRYING TO FIX IT! al
C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2002, L. Joyeux & P. A. Penczek                       *
C=*                                                                    *
C=* University of Texas - Houston Medical School                       *
C=*                                                                    *
C=* Email:  pawel.a.penczek@uth.tmc.edu                                *
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
C **********************************************************************

C   value2(number) : input power spectrum
C  ps (pixel size), cs (constant aberration)
C  contrast (amplitude contrast ration)
C  env(number) : input, contain envelope function squared
C  dzlow1, dzhigh1 : interval search
C  start : value(start:number) are only used

      FUNCTION GETDEFOCUS(VALUE2, X,NUMBER, PS, CS, 
     &                LAMBDA, CONTRAST, XSCORE,ENV, ISTART,ISTOP,ISWI)

CC==X is X coordinates!	  
      	IMPLICIT NONE
      	INCLUDE 'CMBLOCK.INC'
      	INCLUDE 'F90ALLOC.INC'

      	INTEGER   NUMBER,NMB
      	REAL       VALUE(NUMBER), ENV(NUMBER), VALUE2(NUMBER),X(NUMBER)
      	REAL      PS, CS, LAMBDA, CONTRAST, DZLOW1, DZHIGH1, DZINC, DZ
     	 REAL      DZLOW, DZHIGH,XXX,XSCORE

C        GFORT WILL NOT ACCEPT REAL PARAMETERS FOR LOOPING   Sept 2025 al
      	INTEGER      IDZLOW, IDZHIGH,IDZINC, IDZ
      	REAL      GETDEFOCUS
	COMMON /SEARCH_RANGE/ DZLOW1,DZHIGH1,NMB
      	INTEGER          I, J, L, ISTART,ISTOP,ISWI
      	DOUBLE PRECISION SUMA, SUMB, DIFF, DIFFMIN, DZMIN      

      	DZLOW  = DZLOW1 
     	DZHIGH = DZHIGH1
      	DZINC  = (DZHIGH-DZLOW)/10.0 
      	IF (DZINC>1000.0) DZINC = 1000.0
	IF (ISWI.EQ.0) THEN
		DIFFMIN=HUGE(DIFFMIN)
        ELSE
 		DIFFMIN=-HUGE(DIFFMIN)
        ENDIF

        DO
C           GFORT WILL NOT ACCEPT REAL PARAMETERS FOR LOOPING   Sept 2025 al
            IDZLOW  = DZLOW    ! buggy change of limits ,,, redefined below al
	    IDZHIGH = DZHIGH 
	    IDZINC  = DZINC
	    
            DO IDZ  = IDZLOW, IDZHIGH, IDZINC
	      DZ    = IDZ
	      
              CALL CTF_SIGNAL1(VALUE,X,NUMBER, PS, CS, 
     &                      LAMBDA, DZ, CONTRAST)

              IF (ISWI.EQ.0) THEN
                  DO L = ISTART, ISTOP
		  
CCC                     DIFF=DIFF+(VALUE(L)*DBLE(VALUE(L))*ENV(L)-
CCC=CHECK ENVELOPE EFFECT  was buggy?? sept 2025

 		        DIFF = DIFF+(VALUE(L)*DBLE(VALUE(L))-
     &                      DBLE(VALUE2(L)))**2
                  ENDDO
                  IF (DIFF<DIFFMIN .OR. DZ .EQ. DZLOW) THEN
                     GETDEFOCUS = DZ
                     DIFFMIN = DIFF
                  END IF
               ELSE 
                  VALUE(ISTART:ISTOP) = VALUE(ISTART:ISTOP)**2
     &				*ENV(ISTART:ISTOP)
                  DIFF = 0.0
                  SUMA = 0.0
                  SUMB = 0.0
                  DO L = ISTART, ISTOP
                     DIFF = DIFF + VALUE(L)*DBLE(VALUE2(L))
                     SUMA = SUMA + VALUE(L)*DBLE(VALUE(L))
                     SUMB = SUMB + VALUE2(L)*DBLE(VALUE2(L))
                  END DO
                  CLOSE(79)
                  DIFF = DIFF  / (DSQRT(SUMA*SUMB)*(ISTOP-ISTART+1))
                  IF(DIFF>DIFFMIN .OR. DZ .EQ. DZLOW) THEN
                     GETDEFOCUS = DZ
                     DIFFMIN = DIFF
                     XSCORE=DIFFMIN
                  ENDIF
               ENDIF
            ENDDO
            DZLOW = GETDEFOCUS -DZINC*2
            IF(DZLOW<0.0) DZLOW = 0.0
            DZHIGH = GETDEFOCUS +DZINC*2
            DZINC = DZINC / 10.0
            IF(DZINC<1.0) RETURN
	ENDDO

      	END
