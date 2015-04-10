C++*********************************************************************
C
C    HELS.F
C                       REWRITE                    JUN 2009 ARDEAN LEITH
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
C  HELS(NDOC,NDOUT)
C
C  PURPOSE:  CREATES SELECTION DOCUMENT FILES FOR A GIVEN THRESHOLD 
C            FROM A DENDROGRAM. 
C            ONLY WORKS IN CONJUNCTION WITH HIERARCHICAL 
C            CLUSTERING OUTPUT FROM 'CL HC' OR 'CL CLA'.  
C
C--*********************************************************************

        SUBROUTINE HELS(NDOC,NDOUT)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        INTEGER, PARAMETER       :: NDLI = 4
        REAL                     :: DLISTIN(NDLI),DLISTOUT(NDLI)

        CHARACTER(LEN=MAXNAM)    :: DOCNAM,OUTDOC,FILPAT
        CHARACTER(LEN=23)        :: COMMENT
        CHARACTER(LEN=1)         :: NULL
        LOGICAL                  :: EMPTY,TRUNCIT,NEEDNEW

        NULL = CHAR(0)
C                  123456789 123456789 1234
        COMMENT = 'THRESHOLD LEVEL:       '
    
        THRESH = 0.0
        CALL RDPRM1S(THRESH,NOT_USED,'THRESHOLD % (0 .. 100)', IRTFLG)
	IF (IRTFLG .NE. 0)  RETURN
        TRUNCIT  = (THRESH .GT. 0)   

C       OPEN INPUT DOC FILE NAME
        CALL OPENDOC(DOCNAM,.TRUE.,NLET,NDOC,NICDOC,.TRUE.,
     &               'DENDROGRAM DOCUMENT',
     &               .TRUE.,.FALSE.,.FALSE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0)  RETURN
  
C       GET TEMPLATE FOR OUTPUT DOC FILE NAME      
        NMAX = 0
        CALL  FILSEQP(FILPAT,NLET,IDUM,NMAX,IDUM2,
     &                'TEMPLATE FOR SELECTION DOC',IRTFLG)
        IF (IRTFLG .NE. 0)  RETURN

        NTOTAL    = 0
        ICOL      = 2
        EMPTY     = .TRUE.
        NCLAS     = 0
        NOBJ      = 0
        NEEDNEW   = .TRUE.
        LASTCLASS = -1

        DO    ! ENDLESS LOOP
           CALL LUNDOCREDNXT(NICDOC,IKEYGOT,DLISTIN,NDLI-1,
     &                       IGO,ICOUNT,IRTFLG)
           !write(6,*) ikeygot,icount,dlistin(1),irtflg
 
           IF (IRTFLG .EQ. 1)  THEN
              CALL ERRT(101,'ERROR READING DOC. FILE',NE)
              GOTO 9999

           ELSEIF (IRTFLG .EQ. 2 .AND. EMPTY) THEN
C             END OF DOC FILE AND NO KEYS FOUND
              CALL ERRT(101,'DID NOT FIND ANY KEYS IN DOC. FILE',NE)
              GOTO 9999

           ELSEIF (IRTFLG .NE. 0) THEN
C             HAVE FINISHED ALL KEYS IN DOC FILE
              EXIT

           ELSEIF (ICOUNT .LT. ICOL) THEN
              CALL ERRT(102,'REGISTER MISSING IN DOC. FILE',ICOL)
              GOTO 9999
           ENDIF
           EMPTY  = .FALSE.

           ICLASS = DLISTIN(1)
           HITE   = DLISTIN(2)
           IMGNUM = DLISTIN(3)
           NTOTAL = NTOTAL + 1
           NOBJ   = NOBJ + 1

           IF (.NOT. TRUNCIT .AND. ICLASS .NE. LASTCLASS) THEN
              NEEDNEW   = .TRUE.
              NOBJ      = 1
              LASTCLASS = ICLASS
           ENDIF
           
           IF (NEEDNEW) THEN
C             CREATE NEW OUTPUT DOC FILE 
              NCLAS = NCLAS + 1
              IF (TRUNCIT) THEN
                 NFILE = NCLAS
                 CALL HELS_NEWDOC(NCLAS,NFILE,FILPAT,NLET,
     &                         NDOUT,NICDOCO,THRESH,IRTLFG)
              ELSE
                 NFILE = ICLASS
                 CALL HELS_NEWDOC(NFILE,NFILE,FILPAT,NLET,
     &                         NDOUT,NICDOCO,THRESH,IRTLFG)
              ENDIF
              IF (IRTFLG .NE. 0)  GOTO 9999
              NEEDNEW = .FALSE.
           ENDIF

           IHITE = HITE
           NC    = ICLASS
           IF (TRUNCIT) NC = NCLAS
           IF (VERBOSE) WRITE(NOUT,90) NC,NFILE,NOBJ,IHITE,IMGNUM
 90        FORMAT('  Class:',  I4,'   File: ',I5,
     &            '   Element:',I5,'   Index:',I4,
     &            '   ID:', I7)

C          WRITE CURRENT DATA TO DOC FILE
           DLISTOUT(1) = IMGNUM
           DLISTOUT(2) = HITE
           CALL LUNDOCWRTDAT(NICDOCO,NOBJ, DLISTOUT,2,IRTFLG)

           IF (HITE .GE. THRESH .AND. TRUNCIT) THEN          
C              START A NEW CLASS FOR FURTHER BRANCHES
               NOBJ    = 0                   ! RESET OBJECT COUNTER

C              CREATE NEW OUTPUT DOC FILE FOR NEXT LINE
               NEEDNEW = .TRUE.
               CLOSE(NDOUT)
           ENDIF 
        ENDDO

        WRITE(NOUT,93)  NTOTAL,NCLAS
93      FORMAT('  Total number of objects=',I7,
     &         '  Number of classes=',I5)

9999    CLOSE(NDOC)
        CLOSE(NDOUT)
        END

C       ------------------- HELS_NEWDOC ------------------------------

        SUBROUTINE HELS_NEWDOC(NCLAS,NFILE,FILPAT,NLET,
     &                         NDOUT,NICDOCO,THRESH,IRTLFG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=*)         :: FILPAT
        CHARACTER(LEN=MAXNAM)    :: OUTDOC
        CHARACTER(LEN=34)        :: COMMENT


C       CREATE NEW OUTPUT DOC FILE NAME USING FILE
        CALL  FILGET(FILPAT,OUTDOC,NLET,NFILE,IRTFLG)
        IF (IRTFLG .NE. 0)  RETURN

C       OPEN OUTPUT DOC FILE 
        IRTFLG = -9    ! QUIET
        CALL OPENDOC(OUTDOC,.TRUE.,NLET,NDOUT,NICDOCO,.FALSE.,
     &                   '',.FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0)  RETURN

C                  123456789 123456789 123456789 123456
        COMMENT = 'THRESHOLD LEVEL:        CLASS:    '

        WRITE(COMMENT(18:23),FMT='(F6.2)') THRESH
        WRITE(COMMENT(31:34),FMT='(I4)') NCLAS
        CALL LUNDOCPUTCOM(NICDOCO,COMMENT,IRTFLG)
        CALL LUNDOCPUTCOM(NICDOCO,'KEY,     ID,        INDEX', IRTFLG)

        IRTFLG = 0
        END
