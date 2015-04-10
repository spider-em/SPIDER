
C++*********************************************************************
C
C  SGRAF      LONG FILE NAMES                    FEB   89 ARDEAN LEITH
C             SGRAF1 ROUTINE MOVED TO SGRAF1.FOR APRIL 89 ARDEAN LEITH
C             PUT SGRAF1, EPUR4, & CLASYMB BACK  AUG   99 ARDEAN LEITH
C             INCREASED IMAGE NUMBER TO 7 DIGITS AUG   99 ARDEAN LEITH
C             PUT EIGEN IN HERE                  AUG   99 ARDEAN LEITH
C             FIXED TXT FOR NT ERROR             NOV   99 ARDEAN LEITH
C             REMOVED TITLE                      NOV   00 ARDEAN LEITH
C             LUNDOCREDALL PARAMETERS CHANGED    DEC   00 ARDEAN LEITH
C             LKUPDC FIXED                       APR   01 ARDEAN LEITH
C             OPENDOC PARAMETERS                 JUL   03 ARDEAN LEITH
C             REWRITTEN                          OCT 2003 ARDEAN LEITH
C             SEIGEN OPTIONAL                    MAR 2004 ARDEAN LEITH
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
c
C   SGRAF(LUNIP,LTMP,LUNE,LUND)
C
C   PURPOSE: CREATE CORRESPONDENCE ANALYSIS MAPS USING COORDINATE FILES
C
C   NOTES:
C
C      S -- SYMBOL PRINTED IS ONE CHARACTER LONG. ONE IS USED 
C           FOR EACH ACTIVE IMAGE, AND SECOND CHAR. FOR EACH INACTIVE
C           IMAGE
C
C      C -- SYMBOL PRINTED IS ONE CHARACTER LONG, AND IS LOOKED UP IN A
C           CLUSTER FILE
C
C      D -- SYMBOL PRINTED IS ONE CHARACTER LONG, AND IS
C           LOOKED UP IN A DOCUMENT FILE UNDER KEY=ID
C
C      I -- SYMBOL PRINTED IS UP TO 7 CHARACTERS LONG, AND REPRESENTS
C           THE ORIGINAL IMAGE ID
C
C **********************************************************************

	SUBROUTINE SGRAF(LUNIP,LTMP,LUNE,LUND)
	
	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'

	COMMON          BUF(1)
        REAL, ALLOCATABLE, DIMENSION(:)    :: COO,X,Y

	CHARACTER(LEN=MAXNAM)              :: FIL_EIG,FIL_INP,FILPRE
        CHARACTER(LEN=1)                   :: NULL,FOPT,LFLIP,ANS  
        CHARACTER(LEN=1)                   :: PTYPE,TSYMB

        NULL = CHAR(0)

        CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &        '(I)MAGE OR (P)IXEL COORDINATES',NULL,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CALL FILERD(FILPRE,NLET,NULL,
     &       'COORDINATE FILE PREFIX (WITHOUT TRAILING _)',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (ANS .EQ. 'P') THEN 
           FIL_INP = FILPRE(1:NLET) // '_PIX'//NULL
        ELSE
           FIL_INP = FILPRE(1:NLET) // '_IMC'//NULL
        ENDIF
        FIL_EIG = FILPRE(1:NLET) // '_EIG'//NULL

        CALL RDPRI1S(NHOR,NOT_USED,
     &              'NUMBER OF HORIZONTAL PATCHES',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

  	CALL RDPRIS(NFCT1,NFCT2,NOT_USED,
     &           'TWO FACTOR NUMBERS FOR MAP (e.g: 1,5)',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       OPEN IMAGE OR PIXEL COORDINATE FILE AS SEQUENTIAL FORMATED
        CALL OPAUXFILE(.FALSE.,FIL_INP,DATEXC,LUNIP,0,
     &                       'O',' ',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

82      CALL RDPRMC(FOPT,NCHAR,.TRUE.,
     &       '(S)YMBOL, (A)SSIGN SYMBOL, (C)LASS, (D)OC, (I)D',
     &       NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        TSYMB = ' '
        IF (INDEX('SACDI',FOPT) .EQ. 0) THEN
            CALL ERRT(101,'UNKNOWN OPTION',NDUM)
            GOTO 82

        ELSEIF (FOPT .EQ. 'S') THEN
           CALL RDPRMC(TSYMB,NCHAR,.TRUE.,
     &          '1 CHAR. SYMBOL FOR ACTIVE IMAGES/PIXELS',
     &          NULL,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 999
        ENDIF

        IF (FOPT .EQ. 'I') THEN
           WRITE(NOUT,*) 
     &  ' IF IMAGE NUMBERS EXCEED: 999 ONLY POSTSCRIPT OUTPUT AVAILABLE'
        ENDIF

        CALL RDPRMC(PTYPE,NCHAR,.TRUE.,
     &          'PREPARE POSTSCRIPT MAP? (Y/N)',NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        NPAGE = 1
        IF (PTYPE .NE. 'Y') THEN
          CALL RDPRI1S(NPAGE,NOT_USED,
     &               'NUMBER OF PAGES (1,2,3) OR <CR>=1',IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 999
	  NPAGE = MAX(NPAGE,1)

C         NLINE = 0 WILL ACTIVATE DEFAULT SIZE LATER
          NLINE = 0
          CALL RDPRI1S(NLINE,NOT_USED,
     &               'NUMBER OF LINES OR <CR> FOR DEFAULT',IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 999
        ENDIF
      
        PEX = 2.3    
        CALL RDPRM1S(PEX,NOT_USED,'NUMBER OF SD OR <CR>=2.3',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999
        IF (PEX .LE. 0.0) THEN
           CALL ERRT(101,'SD MUST BE > 0.0',NDUM)
           GOTO 999
        ENDIF

        CALL RDPRMC(LFLIP,NCHAR,.TRUE.,
     &    '1=FLIP #1/ 2=FLIP #2/ 3=FLIP 1+2/ <CR>=NO FLIP',NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

C       READ HEADER OF IMAGE/PIXEL COORDINATE FILE
        IF (ANS .EQ. 'P') THEN 
	   READ(LUNIP,*) NITEM,NFAC,NDUM1,NDUM2,NDUM3
        ELSE
	   READ(LUNIP,*) NITEM,NFAC,NDUM1,NDUM2,NDUM3,NDUM4
        ENDIF

	OPEN(LTMP,STATUS='SCRATCH',FORM='UNFORMATTED')

C       OPEN EIGENVECTOR FILE AS SEQUENTIAL FORMATED
        CALL OPAUXFILE(.FALSE.,FIL_EIG,DATEXC,LUNE,0,
     &                       'O',' ',.FALSE.,IRTFLG)
        IF (IRTFLG .EQ. 0) THEN
C          PRINT EIGENVALUES
           CALL SEIGEN(LUNE,NDAT)
        ENDIF
        CLOSE(LUNE)
                                                
C       PRINT OUT PATCH MATRIX IF REQUESTED
	IF (NHOR .NE. 0) THEN
           I1 = 1
           I2 = NHOR
	   WRITE(NDAT,9003)
9003	   FORMAT(/,' PATCH CODE LOOKUP',//)

84	   WRITE(NDAT,9002) ((I),I=I1,I2)
9002	   FORMAT(/,1X,20I4)
           IF (I2 .GE. NITEM) GOTO 85
           I1 = I1+NHOR
           I2 = I2+NHOR
           IF (I2 .GT. NITEM) I2 = NITEM
           GOTO 84
	ENDIF

85	NDIM = NITEM + 2
        ALLOCATE (COO(NFAC),X(NDIM),Y(NDIM),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MWANT = NFAC + NDIM + NDIM 
           CALL ERRT(46,'COO...',MWANT)
           GOTO 999
        ENDIF

	CALL SGRAF1(NDIM,NITEM,NFAC,NFCT1,NFCT2,LUNIP,LTMP,LUND,
     &              FOPT,LFLIP,TSYMB,PTYPE,NLINE,NPAGE,PEX,ANS,
     &              COO,X,Y)


999     CLOSE(LUNIP)
	CLOSE(LTMP)
        CLOSE(LUNE)

C       FREE UP ALLOCATIONS
        IF (ALLOCATED(COO))  DEALLOCATE(COO)
        IF (ALLOCATED(X))  DEALLOCATE(X)
        IF (ALLOCATED(Y))  DEALLOCATE(Y)

        RETURN
	END

C++****************************** SGRAF1 *********************************

        SUBROUTINE SGRAF1(NDIM,ITEM,NFAC,N1,N2,LUNIP,LTMP,LUND,
     &             FOPT,LFLIP,TSYMB,PTYPE,NLINE, NPAGE, PEX, ANS, 
     &             COO, X, Y)

	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'
        INCLUDE 'F90ALLOC.INC'

        INTEGER                           :: XFLIP,YFLIP
	INTEGER, DIMENSION(20)            :: ITO
        CHARACTER(LEN=1)             :: ANS,FOPT,PTYPE,NULL,TSYMB,LFLIP
        CHARACTER(LEN=20)                 :: CSYMB
        CHARACTER(LEN=7)                  :: CBUF
        CHARACTER(LEN=7), DIMENSION(NDIM) :: ID
	REAL                              :: COO(NFAC),X(NDIM),Y(NDIM)
	INTEGER                           :: ICLOOK(NDIM)
        REAL, DIMENSION(:,:), POINTER     :: FLOOK
	CHARACTER(LEN=MAXNAM)             :: DOCNAM

        NULL       = CHAR(0)

        CSYMB(1:1) = TSYMB

C	NO ROTATION OF COORDINATES AROUND ANY AXIS
        XFLIP = 1.0
        YFLIP = 1.0
        IF (LFLIP .EQ. '1') THEN
C	   ROTATE COORDINATES AROUND Y AXIS
           XFLIP = -1.0

        ELSEIF (LFLIP .EQ. '2') THEN
C	   ROTATE COORDINATES AROUND X AXIS
           YFLIP = -1.0

        ELSEIF (LFLIP .EQ. '3') THEN
C	   ROTATE COORDINATES AROUND BOTH AXIS
           XFLIP = -1.0
           YFLIP = -1.0
        ENDIF

	IF (FOPT .EQ. 'S') THEN
           MODE = 1

	ELSEIF (FOPT .EQ. 'C') THEN
           MODE = 1

	ELSEIF (FOPT .EQ. 'I') THEN
           MODE = 4

	ELSEIF (FOPT .EQ. 'D') THEN
C          GET 1 CHAR. SYMBOL FROM DOC FILE
           MODE = 1
           MAXX = 2
           MAXY = 0
           CALL GETDOCDAT('SYMBOL DOCUMENT',.TRUE.,DOCNAM,
     &                    LUND,.TRUE.,MAXX,MAXY,FLOOK,IRTFLG)

        ELSEIF (FOPT .EQ. 'A') THEN
           MODE = 1

C          GET ALL SYMBOLS (MAXIMUM 20)
           WRITE(NOUT,*) 'ENTER UP TO 20 SYMBOLS, ONE BY ONE'

           II    = 1
           CSYMB = '                    '
140        CALL RDPRMC(CSYMB(II:II),NCHAR,.FALSE.,
     &                 'SYMBOL',NULL,IRTFLG)
           IF (NCHAR .LE. 0) CSYMB(II:II) = ' '
           CALL RDPRMI(ITO(II),IDUM,NOT_USED,'USED UNTIL ID#')
           IF (ITO(II) .GT. 0) THEN
C             GOT A SYMBOL ASSIGNMENT
              II = II + 1
              IF (II .LE. 20) GOTO 140
           ENDIF
        ENDIF



C       ITEMA COUNTS THE NUMBER OF ACTIVE IMAGES. ONLY THESE ARE FLAGGED
C       WHEN THE CLUSTER LOOKUP OPTION IS USED
        ITEMA = 0

        II = 1
	DO I = 1, ITEM
          READ(LUNIP,90) (COO(K),K=1,NFAC),RWGT,DIO,FIIII,FACT
90        FORMAT(10(G12.5,1X))

          IIII = FIIII
          X(I) = COO(N1) * XFLIP
          Y(I) = COO(N2) * YFLIP

	  IF (FOPT .EQ. 'S') THEN
C           USE ONE CHAR. SYMBOL FROM CSYMB FOR ACTIVE/INACTIVE IMAGES
            ID(I) = CSYMB(1:1)

	  ELSEIF (FOPT .EQ. 'A') THEN
C           SYMBOL FOR PLOT IS LOOKED UP IN SYMBOL STRING
            IF (IIII .LE. ITO(II)) THEN
               ID(I) = CSYMB(II:II)
            ELSE
               II    = II + 1
               ID(I) = CSYMB(II:II)
            ENDIF

         ELSEIF (FOPT.EQ.'C') THEN
C           SYMBOL FOR PLOT IS LOOKED UP IN CLASS. FILE
	    ITEMA     = ITEMA + 1
	    ID(ITEMA) = CHAR(IIII)
	    X(ITEMA)  = COO(N1) * XFLIP
	    Y(ITEMA)  = COO(N2) * YFLIP

	  ELSEIF (FOPT.EQ.'I') THEN
C           SYMBOL FOR PLOT IS IMAGE NUMBER (6 DIGITS) & A INDICATOR

            IF (PTYPE .NE. 'Y' .AND. IIII .GT. 999) THEN
               WRITE(NOUT,*) 
     &            'IMAGE NUMBERS > 999 CAN NOT USE RESULTS FILE PLOT'
               CALL ERRT(100,'SGRF',NDUM)
               RETURN
            ENDIF

C           RIGHT JUSTIFIED NUMBER
            CBUF(1:7) = '       '
            CALL INTTOCHAR(IIII,CBUF(2:7),NLET,0)
            
	    IF (ANS .EQ. 'I') THEN
C              ANS  = I = IMAGE (NOT PIXEL COORDINATES)
               CBUF(1:1) = 'A'
               IF (FACT .GT. 0) CBUF(1:1) = 'I'
	    ENDIF
            ID(I) = CBUF
            IF (FACT .GT. 0)  ITEMA = ITEMA + 1

	  ELSEIF (FOPT.EQ.'D') THEN
C           GET 1 CHAR. SYMBOL FROM DOC FILE
            ID(I) = CHAR(INT(FLOOK(2,IIII))) // '      ' 
c           write(6,*) i ,' : ',iiii,':',FLOOK(2,IIII),' -->',id(i) 
	  ENDIF
 	ENDDO

	IF (FOPT.EQ.'I' .AND. ANS .EQ. 'I' .AND. ITEMA .EQ. ITEM) THEN
           DO I = 1,ITEM
              ID(I)(1:7) = ID(I)(2:7) // CHAR(0)
           ENDDO
        ENDIF


	IF (FOPT .EQ. 'C') THEN
C          CLASS SYMBOL LOOKUP REQUESTED
           CALL CLASYM(ITEMA,ICLOOK,ID)

C          ID CONTAINS OLD IMAGE ID WHEN CLASYM IS CALLED, AND UPON
C          RETURN FROM CLASYM, ID CONTAINS THE CLASS SYMBOL
           ITEM = ITEMA
	ENDIF	

	IF (PTYPE .EQ. 'Y') THEN
C          CONTINUOUS SCALE POSTSCRIPT PLOT. NO LINE PLOTTER OUTPUT
           CALL HISMAP4(NDIM,ITEM,X,Y,ID,MODE,PEX)

	ELSE
C          PRINTER LINE PLOT ONLY
           CALL HISMAP(NDIM,ITEM,N1,N2,X,Y,ID,
     &                 MODE,NLINE,NPAGE,PEX,NDAT,LTMP)
	ENDIF

9999    IF (ASSOCIATED(FLOOK)) DEALLOCATE(FLOOK)
	RETURN
 	END


C ************************** CLASYM **************************************

        SUBROUTINE CLASYM(NUMIM,IBUF,ID)

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 

        CHARACTER(LEN=MAXNAM) ::  CLUSO
        COMMON /COMMUN/           CLUSO

        CHARACTER * 4  ID(1)
        CHARACTER*1    LC
        INTEGER        IBUF(1)

        NCLU = 15

C       IBUF -- TEMPORARY STORAGE OF CLASS NUMBERS
C       ID   -- INPUT: OLD ID
C               OUTPUT: NEW ID

C       OPEN FILE AS SEQUENTIAL UNFORMATTED
        CALL OPAUXFILE(.TRUE.,CLUSO,DATEXC,-NCLU,0,
     &                       'O','CLUSTER',.TRUE.,IRTFLG)

        REWIND NCLU

        READ(NCLU) ICARD, NFAC, NKLA, KFAC
        READ(NCLU) (IBUF(I), I = 1, ICARD)
        CLOSE(NCLU)

        IF (NUMIM .NE. ICARD) THEN
           WRITE(NOUT,101) NUMIM, ICARD
101        FORMAT(' *** CLUSTER FILE INCONSISTENT WITH IMAGE SET'/,
     &            '     EXPECTED ',I5,'  ENCOUNTERED ',I5)
           RETURN
        ENDIF
  
        DO I = 1, ICARD
           IC = IBUF(I)
           LC = '*'
           IF (IC .LT. 10) LC = CHAR(IC+48)
           IF (IC .GE. 10 .AND. IC .LT. 36) LC = CHAR(IC+64-9)
           ID(I) = LC
        ENDDO

        RETURN
        END


C ************************** SEIGEN **************************************
C
C   PURPOSE: CORRESPONDENCE ANALYSIS EIGEN VALUE OUTPUT
C
C   PARAMETERS:
C	 LUNE = EIGENVECTOR FILE, ASSUMED TO BE OPEN
C	 NDAT = RESULTS FILE  
C
C **********************************************************************

	SUBROUTINE SEIGEN(LUNE,NDAT)

	COMMON D(1)

	READ(LUNE,*) NFAC, SUMW, TRACE

	WRITE(NDAT,91) NFAC
 91	FORMAT('  NUMBER OF FACTORS: ',I8)

	WRITE(NDAT,92) SUMW, TRACE
 92	FORMAT('  TOTAL WEIGHT:', E12.5,/,
     &         '  TRACE: ',E12.5,/)

	WRITE(NDAT,*) ' EIGENVALUES:'
        DO I = 1,NFAC
	   READ(LUNE,*) D(I)
        ENDDO

	CUL = 0.
	DO I = 1, NFAC
           PER = 100.0 * D(I) / TRACE
           CUL = CUL + PER

           WRITE(NDAT,93) I, D(I), PER, CUL
 93        FORMAT(1X,I4,2X,E12.5,2X,E12.5,2X,E12.5)
        ENDDO

	RETURN
	END
