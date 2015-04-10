C++*********************************************************************
C
C  HDLS.F     
C              REWRITE                             JUN 2009 ARDEAN LEITH
C              IF (HITE .GE. THRESH                NOV 2009 ARDEAN LEITH 
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
C  HDLS(NDOC,NDOUT)
C
C  PURPOSE:  CREATES DOCUMENT FILE LISTING NUMBER OF CLASSES
C            FOR A GIVEN THRESHOLD FROM A DENDROGRAM. 
C            ONLY WORKS IN CONJUNCTION WITH HIERARCHICAL 
C            CLUSTERING OUTPUT FROM 'CL HC' OR 'CL CLA'.
C  
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE HDLS(NDOC,NDOUT)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        INTEGER, PARAMETER       :: NDLI = 3
        REAL                     :: DLISTIN(NDLI),DLISTOUT(NDLI)
C       INTEGER, PARAMETER       :: NVALS = 1000
C       INTEGER                  :: INUM(NVALS)

        CHARACTER(LEN=MAXNAM)    :: DOCNAM,OUTDOC
        LOGICAL                  :: EMPTY,TRUNCIT

        THRESH = 0.0
        CALL RDPRM1S(THRESH,NOT_USED,'THRESHOLD % (0 .. 100)', IRTFLG)
	IF (IRTFLG .NE. 0)  RETURN
        TRUNCIT = (THRESH .NE. 0.0)

C       OPEN INPUT DOC FILE NAME
        CALL OPENDOC(DOCNAM,.TRUE.,NLET,NDOC,NICDOC,.TRUE.,
     &               'DENDROGRAM DOCUMENT',
     &               .TRUE.,.FALSE.,.FALSE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0)  RETURN
  
C       OPEN OUTPUT DOC FILE NAME
        CALL OPENDOC(OUTDOC,.TRUE.,NLET,NDOUT,NICDOCO,.TRUE.,
     &         'OBJECT DOCUMENT',.FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9999

        CALL LUNDOCPUTCOM(NICDOCO,
     &                   'KEY = CLASS,    ELEMENTS ',IRTFLG)

        NTOTAL    = 0
        ICOL      = 2
        EMPTY     = .TRUE.
        NCLAS     = 1
        IF (.NOT. TRUNCIT) NCLAS = 0

        NOBJ      = 0
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

           IF (TRUNCIT) THEN
              NOBJ = NOBJ + 1
              IF (HITE .GE. THRESH) THEN
C                PUT OUT STATISTICS ON LAST CLASS
                 IF (VERBOSE) WRITE(NOUT,90) NCLAS,NOBJ
 90              FORMAT('  Class:',  I4,'   Elements:',I5)

C                WRITE FINISHED CLASS DATA TO DOC FILE
                 DLISTOUT(1) = NOBJ
                 CALL LUNDOCWRTDAT(NICDOCO,NCLAS, DLISTOUT,1,IRTFLG)

                 NCLAS = NCLAS + 1
                 NOBJ  = 0
              ENDIF
           ELSE
              IF (ICLASS .NE. LASTCLASS) THEN
C                HAVE STARTED NEW CLASS
                 IF (NOBJ .GT. 0) THEN
C                   PUT OUT STATISTICS ON FINISHED CLASS
                    WRITE(NOUT,90) LASTCLASS,NOBJ

C                   WRITE FINISHED CLASS DATA TO DOC FILE
                    DLISTOUT(1) = NOBJ
                    CALL LUNDOCWRTDAT(NICDOCO,LASTCLASS, 
     &                                DLISTOUT,1,IRTFLG)
                 ENDIF

                 NCLAS     = NCLAS + 1
                 NOBJ      = 1
                 LASTCLASS = ICLASS
              ELSE
                 NOBJ  = NOBJ + 1
              ENDIF
           ENDIF
        ENDDO         

        IF (NOBJ .GT. 0) THEN
C          WRITE DATA FOR LAST CLASS
           NC = NCLAS
           IF (.NOT. TRUNCIT) NC = LASTCLASS

           WRITE(NOUT,90) NC,NOBJ

C          WRITE LAST CLASS DATA TO DOC FILE
           DLISTOUT(1) = NOBJ
           CALL LUNDOCWRTDAT(NICDOCO,NC, DLISTOUT,1,IRTFLG)
        ENDIF

        WRITE(NOUT,93)  NTOTAL,NCLAS
93      FORMAT('  Total number of objects=',I7,
     &         '  Number of classes=',I5)

C       SET OPERATION LINE REGISTERS TO NCLAS
        CALL REG_SET_NSEL(1,1,FLOAT(NCLAS), 0.0,0.0,0.0,0.0, IRTFLG)

9999    CLOSE(NDOC)
        CLOSE(NDOUT)
        END


