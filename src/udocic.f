

C++*********************************************************************
C
C UDOCIC.F      USED LUNDOC                         JUNE 99 ARDEAN LEITH
C               REMOVED ALLOCIT                     MAY  00 ARDEAN LEITH
C               LUNDOCREDALL PARAMETERS CHANGED     DEC  00 ARDEAN LEITH
C               ADDED 'SD IC' CAPABILITY            JUNE 03 ARDEAN LEITH
C               ADDED 'UD NEXT' CAPABILITY          FEB  07 ARDEAN LEITH
C               ADDED 'UD FIND' CAPABILITY          JUN  08 ARDEAN LEITH
C               NEXTKEY SET TO 1                    MAY  09 ARDEAN LEITH
C               ICOLWANT NLIST BUG                  NOV  09 ARDEAN LEITH
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
C   UDOCIC(NDOCT,ENDIT)
C
C   PURPOSE:  'UD IC' OPERATION DUPLICATES FUNCTION OF 'UD'
C             BUT IT USES IN-CORE STORAGE, WHICH MAKES SUBSEQUENT READ 
C             ACCESSES QUICK.  IT SOLICITS DOC FILE NAME.  THE
C             REGISTERS ARE SPECIFED ON THE COMMAND LINE
C   
C   PARAMETERS:  NDOCT     LUN FOR INPUT DOC FILE                (SENT)
C                ENDIT     FLAG FOR ENDING USE                   (SENT)
C
C   USAGE: 'UD IC 11,X11'       - FIRST ACCESS: READ DOC FILE INTO CORE
C          'UD IC 11,X11'       - SUBSEQUENT ACCESSES: READ FROM CORE
C          'UD ICE'             - TERMINATE CORE ACCESS OF CURRENT DOC
C          'UD NEXT [key],[r1]' - SUBSEQUENT ACCESSES: READ FROM CORE
C          'UD NEXT E'          - TERMINATE CORE ACCESS OF CURRENT DOC
C          'UD FIND [key],[r1]' - SUBSEQUENT ACCESSES: READ FROM CORE
C          'UD FIND E'          - TERMINATE CORE ACCESS OF CURRENT DOC
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE UDOCIC(NDOCT,ENDIT)

        USE DOCIC_INFO

        INCLUDE 'CMBLOCK.INC' 

        CHARACTER(LEN=MAXNAM)       :: DOCNAM
        REAL, DIMENSION(:), POINTER :: IPQ

C       MAXIMUM NUMBER OF REGISTERS RETURNED ON OPERATION LINE
        INTEGER,PARAMETER :: MAXLIST=100
        REAL              :: PLIST(MAXLIST)

        LOGICAL           :: FIRST,ENDIT,UDNEXT,UDFIND
        CHARACTER(LEN=1)  :: NULL

        CALL SET_MPI(ICOMM,MYPID,IRTFLG) ! SETS ICOMM AND MYPID

        NULL   = CHAR(0)
        UDNEXT = (FCHAR(4:6) .EQ. 'NEX') .OR. (FCHAR(4:6) .EQ. 'NXT')
        UDFIND = (FCHAR(4:4) .EQ. 'F')  

        CALL FILERD(DOCNAM,NLET,NULL,
     &              'DOCUMENT FILE OR ARRAY NAME~',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
C       SEE IF THIS FILE IS ALREADY OPEN IN-CORE
        CALL ISDOCINCORE(DOCNAM,NIC,MT,IRTFLG)

	IF (ENDIT) THEN
C          WANT TO CEASE USING A DOC FILE -----------------------------
           IF (NIC .GT. 0) THEN
C              CHANGE OLDNAM SO IT CANNOT MATCH ANY FILE NAME.
	       OLDNAM(NIC)(1:1)  = NULL 
               NLETOLDNAM(NIC)   = 0

C              DEALLOCATE CORRESPONDING RUNTIME ARRAY
               IPQ => LOCDOC(NIC)%IPT
               DEALLOCATE(IPQ)
               NULLIFY(LOCDOC(NIC)%IPT)
           ELSE
              IF (MYPID .LE. 0) WRITE(NOUT,91) DOCNAM
91            FORMAT('  DOCUMENT FILE NOT IN-CORE: ',A)
           ENDIF

           IGOT = 0
           DO I = 1,MAXICDOCS 
              IF (OLDNAM(I)(1:1) .NE. NULL) THEN
                 IGOT = IGOT + 1
                 IF (IGOT .EQ. 1 .AND. MYPID .LE. 0) 
     &              WRITE(NOUT,*)' DOCUMENT FILES LEFT IN-CORE:'
                 ILEN =  NLETOLDNAM(I)
                 IF (MYPID .LE. 0) WRITE(NOUT,92) OLDNAM(I)(1:ILEN)
92               FORMAT(5X,A)
              ENDIF
           ENDDO
           IF (MYPID .LE. 0) THEN
              IF (IGOT .LE. 0)
     &           WRITE(NOUT,*)' NO DOCUMENT FILES LEFT IN-CORE'
              WRITE(NOUT,*) ' '
           ENDIF
           RETURN
	ENDIF      ! END OF: IF (ENDIT)	


C       WANT TO GET REGISTERS FROM A DOC FILE ------------------------
        FIRST = .FALSE.
        IF (NIC .EQ. 0 .AND. MT .LE. 0) THEN
C          ERROR -- NAME NOT FOUND, AND NO SPACE LEFT IN DOC NAME LIST
           IF (MYPID .LE. 0) WRITE(NOUT,93) 
93         FORMAT('  *** NO SPACE AVAILABLE IN DOC. NAME LIST, ',
     &            ' CLOSE ANOTHER INCORE FILE FIRST')
           CALL ERRT(100,'UDOCIC',IDUM)
           RETURN
        
        ELSEIF (NIC .EQ. 0 .AND. MT .GT. 0) THEN
C          NAME NOT FOUND, SO FILE IS NOT YET IN-CORE, SPACE AVAILABLE
           FIRST  = .TRUE.
           NIC    = MT

C          DELAY SETTING OLDNAM UNTIL EVERYTHING IS OK ON RETRIEVAL
           OLDNAM(NIC)(1:1) = NULL
           NLETOLDNAM(NIC)  = 0
        ENDIF

C       DOCUMENT NAME FOUND OR BUFFER SPACE AVAILABLE ----------------

C       REGISTER LINE ALREADY LOADED IN RDPR 
        IF (UDNEXT .OR. UDFIND) THEN
           CALL REG_GET_USED(NLIST)
           NLIST = NLIST - 1
           IKEY  = 1
        ELSE
C          PARSE REGISTER LINE TO GET IKEY & NLIST 
           IX = 6
           CALL REG_DOC_PARSE(FCHAR(IX:),COMOUT,IKEY,NLIST,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ENDIF

C       REGISTER LIST OK 
        IF (IKEY .LT. 0) THEN
C          DESIRE A COMMENTED KEY. 
           CALL ERRT(101,'USE <UD> TO RETRIEVE COMMENT KEYS',NE)
           RETURN

        ELSEIF (NLIST .GT. MAXLIST) THEN
           MAXLISTT = MAXLIST
           CALL ERRT(102,'MAX. NUMBER FOR IN-CORE REGISTERS',MAXLISTT)
           RETURN

        ELSEIF (FIRST) THEN
C          OPEN THE PHYSICAL DOC FILE USING EXTENSION DATEXC
           CALL OPENDOC(DOCNAM,.TRUE.,NLET,NDOCT,NDOC,.FALSE.,' ',
     &               .TRUE.,.FALSE.,.FALSE.,NEWFILE,IRTFLGT)
           IF (IRTFLGT .NE. 0) RETURN

C          FIND MAXY BY READING FROM THE ON-DISK DOC FILE
           CALL LUNDOCINFO(NDOC,MAXY,MAXREGS,KEYSINUSE,.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          SET ARRAY DIMENSIONS, EACH ARRAY LINE INCLUDES A KEY REGISTER 
C          ALLOCATE THIS INCORE ARRAY (ONLY NEED MAXY ROWS)
           MEMWANT = (MAXREGS+1) * MAXY
           ALLOCATE(IPQ(MEMWANT),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'IPQ',NE)
              GOTO 999
           ENDIF

C          KEEP ARRAY POINTER IN LOCDOC
           LOCDOC(NIC)%IPT => IPQ

C          RECOVER DOC FILE CONTENTS AND PUT THEM IN IPQ
           MAXCOLS = MAXREGS + 1
           CALL LUNDOCREDALL(NDOC,IPQ(1),MAXCOLS,MAXY,.TRUE.,
     &                       NGOT,IRTFLG)
           CLOSE(NDOCT)
           IF (IRTFLG .NE. 0) THEN
C              DEALLOCATE THE CORRESPONDING RUNTIME ARRAY
               DEALLOCATE(IPQ)
               NULLIFY(LOCDOC(NIC)%IPT)
               RETURN
           ENDIF

C          SAVE ARRAY SIZE SIZES
           OLDNAM(NIC)     = DOCNAM(1:NLET) 
           NLETOLDNAM(NIC) = NLET
           NUMCOLS(NIC)    = MAXCOLS
           NUMKEYS(NIC)    = MAXY

C          INITIALIZE RETRIEVAL KEY FOR LUNDOCGETNEXT USE
           NEXTKEY(NIC) = 1

        ELSE
C          GET ARRAY SIZE SAVED WHEN FILLED
           MAXCOLS = NUMCOLS(NIC)
           MAXY    = NUMKEYS(NIC)

C          RECOVER DOC FILE INFO FROM INCORE IMAGE 
           IPQ => LOCDOC(NIC)%IPT
        ENDIF

        IF (.NOT. UDFIND .AND. NLIST .GT. (MAXCOLS - 1)) THEN
           IF (MYPID .LE. 0) WRITE(NOUT,1447) NLIST,MAXCOLS-1
1447       FORMAT('  *** NUMBER OF REGISTERS REQUESTED: ',I2,
     &               ' NUMBER AVAILABLE IN-CORE: ',I2,/)
           CALL ERRT(100,'UDOCIC',NE)
           NLIST = MAXCOLS - 1
        ENDIF

        IF (UDFIND) THEN
C          GET INPUT FOR 'UD FIND' OPERATION
           IF (NLIST .GT. (MAXCOLS - 1)) THEN
           IF (MYPID .LE. 0) WRITE(NOUT,1446) NLIST,MAXCOLS-1
1446          FORMAT('  NUMBER OF REGISTERS REQUESTED: ',I2,
     &               ' SET TO NUMBER AVAILABLE IN-CORE: ',I2,/)
              NLIST = MAXCOLS - 1
           ENDIF

           COLWANT  = 1.0
           VALWANT  = 1.0
           CALL RDPRM2S(COLWANT,VALWANT,NOT_USED,
     &                  'REGISTER COLUMN & VALUE WANTED',IRTFLG)
           IF (IRTFLG .NE. 0) THEN
C             END OF FILE REACHED, SET REGISTER CONTENTS
              NLIST    = 1
              PLIST(1) = 0
              IRTFLG   = 0
              RETURN
           ENDIF

           ICOLWANT = COLWANT + 1.5  ! KEY IS IN COL: 1

           IF (ICOLWANT .GT. MAXCOLS) THEN
              IF (MYPID .LE. 0) WRITE(NOUT,1445) ICOLWANT -1
1445          FORMAT('  *** REGISTER COLUMN: ',I2, ' NOT AVAILABLE ',/)
              CALL ERRT(102,'REGISTER NOT AVAILABLE',ICOLWANT)
              NLIST    = 1
              PLIST(1) = 0
              IRTFLG   = 1
              RETURN
           ENDIF
           ! write(6,*) ' icolwant,valwant: ',icolwant,valwant

C          READ NEXT LINE OF DOC FILE UNTIL DESIRED LINE IS FOUND
           DO WHILE (.TRUE.)      ! ENDLESS LOOP

              CALL LUNDOCREDNXT(-NIC,IKEY,PLIST(2),MAXLIST-1,
     &                          IDUM,NLIST,IRTFLG)
C             write(6,*)'k,igo,1,2,3,:',ikey,igo,(plist(i),i=1,3),irtflg
  
              IF (IRTFLG .EQ. 2) THEN
C                END OF FILE REACHED, SET REGISTER CONTENTS
                 NLIST    = 1
                 PLIST(1) = 0
                 IRTFLG   = 0
                 NEXTKEY(-NIC) = 1
                 EXIT
              ELSEIF (PLIST(ICOLWANT) .EQ. VALWANT) THEN
C                DESIRED REGISTER LINE FOUND
                 NLIST    = NLIST + 1
                 PLIST(1) = IKEY
                 EXIT
              ENDIF
c             write(6,*) ' plist: ',(plist(i),i=1,3)
           ENDDO

C          SET RETRIEVAL KEY FOR LUNDOCGETNEXT USE
           NEXTKEY(NIC) = 1

        ELSEIF (UDNEXT) THEN
C          READ NEXT LINE OF DOC FILE
           CALL LUNDOCREDNXT(-NIC,IKEY,PLIST(2),MAXLIST-1,
     &                       IDUM,NLIST,IRTFLG)
c          write(6,*) ' ikey,igo,nlist,irtflg: ',ikey,igo,nlist,irtflg
 
           IF (IRTFLG .EQ. 2) THEN
C             END OF FILE REACHED, SET REGISTER CONTENTS
              NLIST         = 1
              PLIST(1)      = 0
              IRTFLG        = 0
              NEXTKEY(-NIC) = 1
           ELSE
              NLIST    = NLIST + 1
              PLIST(1) = IKEY
           ENDIF
c          write(6,*) ' plist: ',(plist(i),i=1,3)

        ELSE
C          READ IKEY CONTENTS FROM IPQ AND PLACE IN PLIST
           CALL LUNDOCGETKEY(NDOC,IPQ(1),MAXCOLS,MAXY,IKEY,PLIST,
     &                       NLIST,.TRUE.,IRTFLG)
        ENDIF
        IF (IRTFLG .NE. 0) RETURN

C       SET REGISTER CONTENTS FROM PLIST --------------------------
        CALL REG_SET_NSELA(NLIST,PLIST,.TRUE.,IRTFLG)

        NICLAST = NIC

        IF (.NOT. SILENT .AND. MYPID .LE. 0) WRITE(NOUT,*) ' '
	RETURN

999     CLOSE(NDOCT)
#ifdef USE_MPI
        CALL MPI_BARRIER(ICOMM, MPIERR)
#endif
        RETURN

	END

