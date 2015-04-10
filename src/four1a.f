C++*********************************************************************
C
C     FOUR1A.F
C                                OPFILEC          FEB 03 ARDEAN LEITH    
C                                COS              JUL 12 G. KISHCHENKO    
C                                FQ Q PARAMS      OCT 12 ARDEAN LEITH    
C                                RAISED SINC      OCT 13 ARDEAN LEITH    
C	
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2013  Health Research Inc.,                         *
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
C  FOUR1A()
c
C  PURPOSE:  APPLIES FOURIER FILTERS TO 2-D OR 3-D REAL PICTURES
C            'FQ' : QUICK FILTERING (IN CORE, 2-D OR 3-D)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

	SUBROUTINE FOUR1A

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC'

        CHARACTER (LEN=MAXNAM)              :: FILNAM

        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: BB
       
        INTEGER                             :: MAXIM,NX,NY,NZ
        INTEGER                             :: IOPT,IRTFLG,NOT_USED
        INTEGER                             :: N2X,N2Y,N2Z,LSD
 
        INTEGER, PARAMETER                  :: LUN1 = 21
        INTEGER, PARAMETER                  :: LUN2 = 22

        MAXIM = 0
   	CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NX,NY,NZ,
     &  	   MAXIM,'INPUT',.FALSE.,IRTFLG)
	IF (IRTFLG .NE. 0) RETURN

        IF (IFORM .NE. 1 .AND. IFORM .NE. 3) THEN
           CALL ERRT(101,'OPERATION WORKS ON REAL INPUT',IDUM)
	   GOTO 999

	ELSEIF (IFORM == 1 .AND. NY == 1) THEN
           CALL ERRT(101,'OPERATION DOES NOT WORK ON 1D INPUT',IDUM)
	   GOTO 999
	ENDIF
	
        MAXIM = 0
	CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',IFORM,NX,NY,NZ,
     &		     MAXIM,'OUTPUT',.FALSE.,IRTFLG)
	IF (IRTFLG. NE. 0) GOTO 999

 
 1000   WRITE(NOUT,1009)
 1009   FORMAT
     &  ('  1 - LOW-PASS,                2 - HIGH-PASS'         ,/,
     &   '  3 - GAUSS.  LOW-PASS,        4 - GAUSS.  HIGH-PASS' ,/,
     &   '  5 - FERMI                    6 - FERMI'             ,/,
     &   '  7 - BUTTER. LOW-PASS,        8 - BUTTER. HIGH-PASS' ,/,
     &   '  9 - RAISED COS. LOW-PASS,   10 - RAISED COS. HIGH-PASS' ,/,
     &   '  9 - RAISED COS. LOW-PASS,   10 - RAISED COS. HIGH-PASS' ,/,
     &   ' 11 - BUTTER. ELLIP LOW-PASS, 12 - BUTTER. ELLIP HIGH-PASS',/,
     &   ' 13 - RAISED SINC WINDOW')

        CALL RDPRI1S(IOPT,NOT_USED,'FILTER TYPE (1-13)',IRTFLG)
	IF (IRTFLG. NE. 0) GOTO 999
	
        IF (IOPT < 1 .OR. IOPT > 13) THEN
           CALL ERRT(102,'ILLEGAL VALUE FOR FILTER TYPE',IOPT)
	   GOTO 1000
	ENDIF

       IF (NZ > 1 .AND. IOPT == 13) THEN
           CALL ERRT(101,'OPERATION NOT IMPLEMENTED FOR VOLUMES',IOPT)
	   GOTO 1000
	ENDIF

	IF (FCHAR(4:5) == 'NP')  THEN
	    N2X = NX
            N2Y = NY
            N2Z = NZ
        ELSE
           N2X = 2*NX
           N2Y = 2*NY
	   IF (NZ > 1)  THEN
	      N2Z = 2*NZ
	   ELSE
	      N2Z = 1
           ENDIF
	ENDIF

	LSD = N2X + 2 - MOD(N2X,2)

        ALLOCATE (BB(LSD,N2Y,N2Z), STAT=IRTFLG)           
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'FOUR1A; BB',LSD*N2Y*N2Z)
           GOTO 999
        ENDIF

        IF (IFORM == 1)  THEN
           WRITE(NOUT,"(A,I6,' x',I6)") 
     &                '  Dimensions used:',N2X,N2Y

           CALL FQ_Q(IOPT,LUN1,LUN2, BB, 
     &               LSD,N2X,N2Y, NX,NY, IRTFLG)
	   IF (IRTFLG .NE. 0) GOTO 999

        ELSE
           WRITE(NOUT,"(A,I6,' x',I6,' x',I6)") 
     &                '  Dimensions used:',n2x,n2y,n2z
          
           CALL FQ3_P(LUN1,LUN2, BB,
     &                LSD,N2X,N2Y,N2Z,NX,NY,NZ,IOPT)

	   IF (IOPT < 0) THEN
	      CALL ERRT(101,'USING FFTW',NE)
              GOTO 999
	   ENDIF

        ENDIF


999     IF (ALLOCATED(BB)) DEALLOCATE(BB)

        CLOSE(LUN1)
        CLOSE(LUN2)
       
        END

