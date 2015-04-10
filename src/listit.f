
C++*********************************************************************
C
C  LISTIT.F   REWRITTEN & REMOVED FROM OPENPIC     AUG 96 ARDEAN LEITH
C             LONG FILENAMES                       DEC 88 ARDEAN LEITH
C             IMAGE  LIST BUG FIXED                SEP 98 ARDEAN LEITH
C             FIXED KEY FOR DOCFILE OF COL OUTPUT  AUG 99 ARDEAN LEITH
C             USED NBUFSIZ                         AUG 02 ARDEAN LEITH
C             LUNRED                               FEB 03 ARDEAN LEITH
C             'Pixel' COL & ROW INTERCHANGE BUG    JAN 07 ARDEAN LEITH
C             INCORE COMMENT WRITE BUG             FEB 07 ARDEAN LEITH
C             'LI T x' ADDED                       AUG 14 ARDEAN LEITH
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
C    LISTIT()
C
C    PURPOSE:   TO LIST CONTENTS OF SELECTED IMAGE ELEMENTS 
C
C    PARAMETERS:
C         FILNAM     FILE NAME (ONLY USED FOR TITLE)
C         LUN        LOGICAL UNIT NUMBER OF FILE
C         NX,NY  IMAGE DIMENSIONS
C         NZ     IMAGE DIMENSIONS SLICES
C         DOCPRNT    PRINT TO DOC FILE
C         TERMPRNT   PRINT TO TERMINAL
C                    OTHERWISE PRINT TO RESULTS FILE
C
C--*********************************************************************

	SUBROUTINE LISTIT(FILNAM,LUN,NX,NY,NZ,DOCPRNT,TERMPRNT)

	INCLUDE 'CMBLOCK.INC' 
	INCLUDE 'CMLIMIT.INC' 

	INTEGER                :: IROW(NBUFSIZ),ICOL(NBUFSIZ)
	INTEGER                :: ISLI(NBUFSIZ)
	INTEGER                :: IBUF(NBUFSIZ),JOUT(NBUFSIZ)
	REAL                   :: AOUT(NBUFSIZ),A0(2*NBUFSIZ+1)

        COMPLEX                :: BCOM(1)
        EQUIVALENCE              (A0,BCOM)

        CHARACTER (LEN=MAXNAM) :: DOCNAM
        CHARACTER (LEN=160)    :: COMMENT

        CHARACTER (LEN=*)      :: FILNAM
        CHARACTER (LEN=1)      :: ANS,OTYPE,FTYPE
        LOGICAL                :: INTPRNT, DOCPRNT, TERMPRNT, PHAMPRNT
        LOGICAL                :: NEWFILE
        REAL                   :: DLIST(10)
        CHARACTER (LEN=1)      :: NULL = CHAR(0)

	INTEGER, PARAMETER     :: NDOC = 76


        NLET = INDEX(FILNAM,NULL) - 1

C       DETERMINE OUTPUT OPTION
	NLIS     =  NDAT
        IF (TERMPRNT) NLIS = NOUT
        INTPRNT = .FALSE.


        IF (DOCPRNT) THEN
C           OPEN THE DOC FILE NOW
            CALL OPENDOC(DOCNAM,.TRUE.,NLETD,NDOC,NDOCOUT,.TRUE.,
     &          'OUTPUT DOC FILE',.FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 999

            COMMENT = ' FROM IMAGE FILE: '// FILNAM(1:NLETD) // CHAR(0)
            NL      = LNBLNKN(COMMENT)
            CALL LUNDOCPUTCOM(NDOCOUT,COMMENT(1:NL),IRTFLG)

        ELSE
            WRITE(NLIS,91)  FILNAM(1:NLET)
 91         FORMAT(' -------- LISTING FROM IMAGE FILE: ',A)
        ENDIF

C       SET FOR ALL COLUMS, ROWS, FIRST SLICE ONLY
        NIROW   = 0
        NICOL   = 0
        NISLI   = 1
        ISLI(1) = 1

        CALL RDPRMC(OTYPE,NCHAR,.TRUE.,
     &    'HEADER, PIXEL, ROW, COLUMN, IMAGE, OR WINDOW (H/P/R/C/I/W)',
     &     NULL,IRTFLG)

        IF (OTYPE == 'H') THEN
C          WANT HEADER ONLY
           NICOL = NBUFSIZ
           CALL RDPRAI(ICOL,NBUFSIZ,NICOL,1,256,
     &              'HEADER POSITION(S)', NULL,IRTFLG)

        ELSEIF (OTYPE == 'P') THEN
C          WANT ONE PIXEL ONLY
           CALL RDPRIS(IROW(1),ICOL(1),NOT_USED,
     &        'ROW & COLUMN',IRTFLG)
           NIROW = 1
           NICOL = 1

        ELSEIF (OTYPE == 'R') THEN
C          WANT ONE ROW ONLY
           NIROW = NBUFSIZ
           CALL RDPRAI(IROW,NBUFSIZ,NIROW,1,NY,'ROW(S)',
     &                 NULL,IRTFLG)
           DO I = 1,NX
              ICOL(I) = I
           ENDDO
           NICOL = NX

        ELSEIF (OTYPE == 'C') THEN
C          WANT ONE COL. ONLY
           NICOL = NBUFSIZ
           CALL RDPRAI(ICOL,NBUFSIZ,NICOL,1,NX,'COLUMN(S)',
     &                 NULL,IRTFLG)
           DO I = 1,NY
              IROW(I) = I
           ENDDO
           NIROW = NY

        ELSEIF (OTYPE == 'I') THEN 
C          WANT ONE WHOLE IMAGE 
 
           IF (NZ > 1) THEN
C             WANT ONE IMAGE ONLY FROM THE VOLUME
              NISLI = NBUFSIZ
              CALL RDPRAI(ISLI,NBUFSIZ,NISLI,1,NZ,'SLICE(S)',
     &            NULL,IRTFLG)
           ENDIF

           DO I = 1,NY
              IROW(I) = I
           ENDDO
           NIROW = NY
           DO I = 1,NX
              ICOL(I) = I
           ENDDO
           NICOL = NX

        ELSEIF (OTYPE == 'W') THEN
C          WANT WINDOW FROM IMAGE
           NICOL = NBUFSIZ
           CALL RDPRAI(ICOL,NBUFSIZ,NICOL,1,NX,'COLUMNS(S)',
     &                 NULL,IRTFLG)
           NIROW = NBUFSIZ
           CALL RDPRAI(IROW,NBUFSIZ,NIROW,1,NY,'ROW(S)',
     &                 NULL,IRTFLG)
           IF (NZ > 1) THEN
C             MAY WANT MORE THAN ONE SLICE
              NISLI = NBUFSIZ
              CALL RDPRAI(ISLI,NBUFSIZ,NISLI,1,NZ,'SLICE(S)',
     &                    NULL,IRTFLG)
           ENDIF
        ENDIF            

        IF (.NOT. DOCPRNT .AND. IFORM > 0 .AND.
     &      OTYPE .NE. 'H') THEN
           CALL RDPRMC(FTYPE,NCHAR,.TRUE.,
     &         'FLOATING POINT OR INTEGER FORMAT OUTPUT (F/I)', 
     &          NULL,IRTFLG)
           IF (FTYPE == 'I') INTPRNT = .TRUE.
        ENDIF

        PHAMPRNT = .FALSE.
        IF (IFORM < 0 .AND. OTYPE .NE. 'H') THEN
C          FOURIER FILE
           CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &          'PHASE & MODULI LISTING? (Y/N)', NULL,IRTFLG)
           IF (ANS .NE. 'N') PHAMPRNT = .TRUE.
        ENDIF

        IF (OTYPE == 'H') THEN
C         PRINT OUT HEADER STUFF

C         READ 256 HEADER POSITIONS INTO A0 BUFFER
          CALL LUNGETVALS(LUN,1,256,A0,IRTFLG)

          IF (DOCPRNT) THEN
C           SAVE IN DOC FILE, DO NOT APPEND
            COMMENT = '  FILE HEADER VALUES ---------------------'
            NL      = LNBLNKN(COMMENT)
            CALL LUNDOCPUTCOM(NDOCOUT,COMMENT(1:NL),IRTFLG)

             DO I = 1, NICOL
                CALL LUNDOCWRTDAT(NDOCOUT,ICOL(I),A0(I),1,IRTFLG)
             ENDDO
          ELSE
C            FLOATING POINT HEADER LISTING
             WRITE(NLIS,*) ' '
             WRITE(NLIS,*) 'FILE HEADER VALUES -----------------'

             WRITE(NLIS,302) (A0(ICOL(J)),J=1,NICOL)
302          FORMAT(5(1X,G12.4))
          ENDIF

          CALL REG_GET_USED(NSEL_USED)
          IF (NSEL_USED > 0) THEN
             NSEL_USED = MIN(NSEL_USED,NICOL)
             CALL REG_SET_NSELA(NSEL_USED,A0,.TRUE.,IRTFLG)
          ENDIF

          GOTO 999
        ENDIF

        IF (PHAMPRNT) THEN
           IF (DOCPRNT) THEN
              COMMENT = ' LISTING OF FOURIER MODULI AND PHASES ------'
              NL      = LNBLNKN(COMMENT)
              CALL LUNDOCPUTCOM(NDOCOUT,COMMENT(1:NL),IRTFLG)
           ELSE
              WRITE(NLIS,*) ' '
              WRITE(NLIS,*)'LISTING OF FOURIER MODULI AND PHASES ------'
           ENDIF
        ELSE
           IF (DOCPRNT) THEN
              COMMENT = ' LISTING OF FILE VALUES ------'
              NL      = LNBLNKN(COMMENT)
              CALL LUNDOCPUTCOM(NDOCOUT,COMMENT(1:NL),IRTFLG)
           ELSE
              WRITE(NLIS,*) ' '
              WRITE(NLIS,*) 'LISTING OF FILE VALUES ------'
           ENDIF
        ENDIF

        IKEY = 0
C       FOR EACH SLICE IN THE REQUESTED LISTING
        DO  KSLICE = 1, NISLI
           ISLICET = ISLI(KSLICE)
      
C          FOR EACH ROW IN THE REQUESTED LISTING
           DO  KROW = 1,NIROW
              IROWT = IROW(KROW)

C             RECOVER THE IMAGE ROW
              IROWNOW = (ISLICET - 1) * NY + IROWT
              CALL REDLIN(LUN,A0,NX,IROWNOW)

              IF (DOCPRNT) THEN
C                OUTPUT TO DOC FILE
                 WRITE(COMMENT,902) ISLICET,IROWT
902              FORMAT('    SLICE: ',I5,'   ROW: ',I5)
                 NL      = LNBLNKN(COMMENT)
                 CALL LUNDOCPUTCOM(NDOCOUT,COMMENT(1:NL),IRTFLG)
              ELSE
                 WRITE(NLIS,903) ISLICET,IROWT
 903             FORMAT(' SLICE: ',I5,'   ROW: ',I5)
              ENDIF

C             GET THE PIXELS WANTED FROM THIS LINE, PUT IN AOUT
              IVALS = 0
              DO KSAM = 1,NICOL
                 IVALS = IVALS + 1
                 JOUT(IVALS) = ICOL(KSAM)
                 AOUT(IVALS) = A0(ICOL(KSAM))
              ENDDO

              IF (PHAMPRNT) THEN
C                WANT PHASE & MODULI FOR FOURIER
                 ANGF  = 180. / 3.14159
                 IVALS = 0
                 DO KSAM = 1,NICOL
                    IF (KSAM .LE. (NX / 2)) THEN
                       L  = ICOL(KSAM)
                       BR = A0(2*(L-1)+1)
                       BI = A0(2*L)
                       AM = SQRT(BI**2+BR**2)
                       PH = 0.0
                       IF (BI .NE. 0. .OR. BR .NE. 0.) 
     &                     PH = ATAN2(BI,BR)*ANGF
                       IVALS       = IVALS + 1
                       AOUT(IVALS) = AM
                       JOUT(IVALS) = ICOL(KSAM)
                       IVALS       = IVALS + 1
                       AOUT(IVALS) = PH
                    ENDIF
                 ENDDO
              ENDIF

              IF (DOCPRNT) THEN
C                PRINT INTO DOC FILE ONLY
                 DO I = 1, IVALS 
C                   REG. 0: KEY
                    IKEY     = IKEY + 1
                    DLIST(1) = IKEY

C                   REG. 1: VALUE
                    DLIST(2) = AOUT(I)

C                   REG. 2: COLUMN
                    DLIST(3) = JOUT(I)

C                   REG. 3: ROW
                    DLIST(4) = IROWT

C                   REG. 4: SLICE
                    DLIST(5) = ISLICET
                    CALL LUNDOCWRTDAT(NDOCOUT,IKEY,DLIST(2),4,IRTFLG)
                 ENDDO

              ELSEIF (INTPRNT) THEN
C                PRINT AS INTEGERS
                 DO  J=1,IVALS
                    IBUF(J) = AOUT(J)
                 ENDDO

                 WRITE(NLIS,602)(IBUF(J),J=1,IVALS)
602              FORMAT(10(1X,I10.6))

              ELSE
C                FLOATING POINT LISTING 
                 WRITE(NLIS,106)(AOUT(J),J=1,IVALS)
 106             FORMAT(5G12.4)
              ENDIF
           ENDDO
        ENDDO
        IF (.NOT. DOCPRNT) WRITE(NLIS,*) ' '

999     CONTINUE
        IF (DOCPRNT) CLOSE(NDOC)

        END


