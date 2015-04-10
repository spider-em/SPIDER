
C++*********************************************************************
C
C  COPYPOS          ALTERED FOR SPIDER           FEB 90   ARDEAN LEITH
C	            MODIFIED FOR F90             10/22/97 yl
C                   MODIFIED FOR F90             FEB 99   ARDEAN LEITH     
C                   MAXNAM                       JUL 14   ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C  COPYPOS(FILOLD,LUNSPI,LUNPOS,NX,NY,NZ)
C
C  PARAMETERS:      
C
C  PURPOSE: CONVERTS A SPIDER IMAGE TO A ASCII POSTSCRIPT FILE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE COPYPOS(FILOLD,LUNSPI,LUNPOS,NX,NY,NZ)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

	CHARACTER(LEN=*)      :: FILOLD
        INTEGER               :: LUNSPI,LUNPOS,NX,NY,NZ

        INTEGER, PARAMETER    :: NXAX = 17008 
        REAL                  :: FARRAY(NXAX)
        CHARACTER(LEN=MAXNAM) :: PSFILE

        INTEGER*1             :: BBUF(NXAX)
        CHARACTER(LEN=1)      :: NULL = CHAR(0)
	REAL	              :: LENGTH,WIDTH
        LOGICAL               :: FLIP

        DATA FLTZER/10E-30/

	IF (IMAMI .NE.1) CALL NORM3(LUNSPI,NX,NY,NZ,FMAX,FMIN,AV)

        IF ( (FMAX - FMIN) < FLTZER) THEN
           WRITE(NOUT,*) ' *** BLANK FILE SKIPPED'
           RETURN
        ENDIF

C       FORMATTED, SEQUENTIAL FILE FOR POSTSCRIPT
8       CALL OPAUXFILE(.TRUE.,PSFILE,DATEXC,LUNPOS,0,'N',
     &                       'POSTSCRIPT OUTPUT',.TRUE.,IRTFLG)
        NLET2 = LNBLNKN(PSFILE)

	WRITE(NOUT,1150) FMIN,FMAX
1150	FORMAT(/'  DENSITY RANGE: ',G12.4,'...',G12.4,/)

	AMIN = FMIN
	AMAX = FMAX
	WRITE(NOUT,900)
900     FORMAT('  (TO REVERSE CONTRAST, MAKE MAX < MIN.)')
        CALL RDPRM2S(AMIN,AMAX,NOT_USED,
     &     'MIN AND MAX DENSITIES FOR THRESHOLDING (OR <CR>)',
     &     IRTFLG)
        IF (IRTFLG == -1) GOTO 8

	IF (AMIN == AMAX) THEN
	    AMIN  = FMIN
	    AMAX  = FMAX
            FLIP = .FALSE.

	ELSEIF (AMIN > AMAX) THEN
            FLIP  = .TRUE.
	    FTEMP = AMIN
	    AMIN  = AMAX
	    AMAX  = FTEMP 
	ENDIF	
    
3	SCALE = 255.0 / (AMAX-AMIN)
	WRITE(NOUT,1151) AMIN, AMAX
1151	FORMAT(/,'  DENSITY VALUES SET: ',G12.4,'....',G12.4)

	ISKIP   = 1
	NX0     = NX
	NY0     = NY
        RATIM   = FLOAT(NX0) / FLOAT(ISKIP)
        RATPA   = FLOAT(NY0) * 18.0 / 23.0
	IF (RATIM >= RATPA) THEN
	   WIDTH = 18.0
	ELSE
	   WIDTH = 23.0 * FLOAT(NX0 / ISKIP) / FLOAT(NY0)
	ENDIF
	LENGTH = WIDTH * FLOAT(NY0) / FLOAT(NX0 / ISKIP)

	WRITE(NOUT,1600) WIDTH, WIDTH/2.54, LENGTH, LENGTH/2.54
1600	FORMAT('  DEFAULT IMAGE SIZE IS ',F5.2,' CM (',
     &     F5.2,' IN) BY ', F5.2,' CM (',F5.2,' IN)',/)

        CALL RDPRM(QWIDTH,NOT_USED,
     &      'DESIRED WIDTH IN CM (OR <CR>)')
        IF (QWIDTH <= 0.0) QWIDTH = WIDTH
	QLENGTH = QWIDTH * FLOAT(NY0) / FLOAT(NX0)

	IF (QWIDTH  .NE. 0 .AND. QWIDTH < WIDTH .AND.
     &      QLENGTH <= LENGTH) WIDTH = QWIDTH

        WRITE(NOUT,902) WIDTH
902     FORMAT('  SELECTED WIDTH: ',F5.2,' CM')

        IDBG = -1
        CALL RDPRI1S(IDBG,NOT_USED,
     &    'BACKGROUND VALUE 0 (BLACK) - 255 (WHITE), (-1 - SKIP)',
     &    IRTFLG)

        WRITE(LUNPOS,94)PSFILE(1:NLET2)
 94     FORMAT(
     &     '%!PS-Adobe-2.0 EPSF-2.0',/,
     &     '%%Title:',A,/,
     &     '%%Creator: SPIDER copypos (CP PO)',/,
     &     '%%Pages: 1',/,
     &     '%%DocumentFonts:',/,
     &     '%%EndComments',/,
     &     '%%EndProlog',/)
      
        WRITE(LUNPOS,95)
95      FORMAT(
     &     '% remember original state',/,
     &     ' /origstate save def',/)

        WRITE(LUNPOS,96)
96      FORMAT(' ',
     &     '% build a temporary dictionary',/,
     &     ' 20 dict begin',/)
     
        WRITE(LUNPOS,97)
97      FORMAT(
     &     ' /Helvetica findfont ',/,
     &     ' 12 scalefont setfont ',/)

	DO IZ = 1,NZ
C           LOOP THROUGH EACH SLICE OF THE IMAGE
	    WRITE(LUNPOS,*) ' 52 112 moveto '

            IF (NZ > 1) THEN
               WRITE(LUNPOS,90) IZ
90             FORMAT('  (Slice ',I4,10X,'-- X -->) show ')
            ENDIF

	    WRITE(LUNPOS,*) ' 55 98 moveto '

            NLET1 = LNBLNKN(FILOLD)
            WRITE(LUNPOS,91) FILOLD(1:NLET1),DATEXC(1:3),PSFILE(1:NLET2)
91          FORMAT(' ( Image: ',A,'.',A,'  Postscript: ',A,') show')

	    CALL COPYPOS2(LUNPOS,0.75,1.7,WIDTH,NX0,NY0)

	    DO IY =  1, NY
               CALL REDLIN(LUNSPI,FARRAY,NX,IY)
	       CALL COPYPOS3(LUNPOS,SCALE,AMIN,AMAX,FLIP,NX0,IDBG)
	    END DO

	    IF (NZ > 1) WRITE(NOUT,*) ' SLICE:',IZ
	    WRITE(LUNPOS,*) ' showpage'
	ENDDO

	WRITE(LUNPOS,*) '% stop using temporary dictionary'
	WRITE(LUNPOS,*) ' end'
	WRITE(LUNPOS,*) ' '

	WRITE(LUNPOS,*) '% restore original state'
	WRITE(LUNPOS,*) ' origstate restore'
	WRITE(LUNPOS,*) ' '

	WRITE(LUNPOS,*) '%%Trailer'
	WRITE(LUNPOS,*) ' '
                           
        CLOSE(LUNSPI)
	CLOSE(LUNPOS)

C       SEE IF USER WANTS A PRINT OUT OF PS FILE NOW
        CALL POPRINT(PSFILE(1:NLET2))

	END



C++*********************************************************************
C
C   COPYPOS2(LUNPOS,XPOS,YPOS,WIDTH,NNX,NY)
C
C   PURPOSE:      ROUTINE FOR INITIALIZING POSTSCRIPT OUTPUT 
C	          SDF 7-JULI-88
C                 ALTERED FOR SPIDER FEB 90 AL
C          
C--*********************************************************************

	SUBROUTINE COPYPOS2(LUNPOS,XPOS,YPOS,WIDTH,NNX,NY)

        DATA  IA,IB,ID,IPIXEL/0,0,0,8/
	
C       USE WIDTHT IN CASE CONSTANT IS PASSED TO HERE
	WIDTHT   = WIDTH / 2.54
	YLENGTH  = WIDTHT * FLOAT(NY) / FLOAT(NNX)
       
        IC = -NY
	WRITE(LUNPOS,10) NNX,NY,IPIXEL,
     &       NNX,IA,IB,IC,ID,NY,XPOS,YPOS,WIDTHT,YLENGTH

10	FORMAT('  /picstr 1 string def'/
     &    '  /grey { ',3I7/
     &	  ' [ ',6I7,' ] '/
     &    '  { currentfile picstr readhexstring pop } '/
     &    '  image } def'/
     &    '  /inch { 72 mul } def'/
     &    2(1X,F4.1,' inch '),' translate'/
     &    1X,F6.2,' inch ',F6.2,' inch scale'/
     &    '  grey')

        END



C++*********************************************************************
C
C  COPYPOS3.FOR -- FEB 90 
C
C   COPYPOS3(LUNPOS,SCALE,AMIN,AMAX,FLIP,NX,IDBG)
C
C   PURPOSE:      CONVERTS A LINE FROM SPIDER IMAGE FILE TO A
C                 POSTSCRIPT (8 BIT NORMALIZED) READABLE FORMAT 
C          
C--*********************************************************************

       SUBROUTINE COPYPOS3(LUNPOS,SCALE,AMIN,AMAX,FLIP,NX,IDBG)

       PARAMETER      (NXAX = 32000)
       COMMON         FARRAY(NXAX),BBUF(NXAX)
 
       INTEGER*1      BVALUE,BBUF,BVAL(4),BVAL4
       EQUIVALENCE    (BVALUE,IVAL), (BVAL, IVAL), (BVAL4,BVAL(4))
       LOGICAL        FLIP

C      OUTPUT ONE LINE OF DATA

       DO IX = 1,NX
          FVAL = FARRAY(IX)
          IF     (FVAL > AMAX) THEN
             RVAL = 255.0
          ELSEIF (FVAL < AMIN) THEN
             RVAL = 0.0
          ELSE
             RVAL = (FVAL - AMIN) * SCALE
          ENDIF

C         GREATER OF (0 AND (THE SMALLER OF RVAL AND 255))
          IVAL = MAX(0,MIN(RVAL, 255.))

          IF (FLIP)   THEN
             IVAL = 255 - IVAL
             IF (IDBG >= 0 .AND. IVAL == 255)  IVAL = IDBG
          ELSE
             IF (IDBG >= 0 .AND. IVAL == 0)   IVAL = IDBG
          ENDIF

#if defined (__osf__) || defined (SP_NT) || defined (__linux__)
C         DEC & NT  USE OTHER BYTE ORDERING
          BBUF(IX) = BVALUE
#else
C         E.G. SGI   
          BBUF(IX) = BVAL(4)
#endif

C         BVALUE = IVAL   ! Through equivalence

        ENDDO

        WRITE(LUNPOS,'(1X,64Z2.2)') (BBUF(IX),IX=1,NX)

        END
