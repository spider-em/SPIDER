C++*********************************************************************
C
C DENLST.F          REWRITE                      NOV. 08 ARDEAN LEITH
C                   ALLOW VERY SMALL VALUES      MAY  09 ARDEAN LEITH
C                   ADD FILE PROMPT              JUN. 09 ARDEAN LEITH
C                   REWRITE                      JUN. 09 ARDEAN LEITH
C                   FILE FILE TYPO               NOV  09 ARDEAN LEITH
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
C   DENLST(LUNDOC,NBR,HIT,ICLASS,ID,IRTFLG)
C
C--*********************************************************************

         SUBROUTINE DENLST(LUNDOC,NBR,HIT,ICLASS,ID,IRTFLG)

         INCLUDE 'CMLIMIT.INC' 

         REAL                   :: HIT(NBR),DLIST(3)
         INTEGER                :: ICLASS(NBR), ID(NBR)
         CHARACTER (LEN=MAXNAM) :: DOCNAM
         LOGICAL                :: NEWFILE
         CHARACTER (LEN=1)      :: NULL

         NULL = CHAR(0)

C        CAN LIST THE DENDORGRAM SHAPE IN DOCUMENT FILE IF DESIRED
         CALL FILERD(DOCNAM,NLET,NULL,'DENDROGRAM DOC.',IRTFLG)

         IF ((NLET .EQ. 1 .AND. DOCNAM(1:1) .EQ. 'Y') .OR.
     &       (NLET .EQ. 1 .AND. DOCNAM(1:1) .EQ. 'y') .OR.
     &       (NLET .EQ. 3 .AND. DOCNAM(1:3) .EQ. 'YES') .OR.
     &       (NLET .EQ. 3 .AND. DOCNAM(1:3) .EQ. 'yes')  ) THEN
C           MUST ASK FOR FILE NAME AGAIN
            CALL FILERD(DOCNAM,NLET,NULL,'DENDROGRAM DOC.',IRTFLG)

         ELSEIF ((NLET .EQ. 1 .AND. DOCNAM(1:1).EQ. 'N') .OR.
     &           (NLET .EQ. 1 .AND. DOCNAM(1:1).EQ. 'n') .OR.
     &           (NLET .EQ. 2 .AND. DOCNAM(1:2).EQ. 'NO') .OR.
     &           (NLET .EQ. 2 .AND. DOCNAM(1:2).EQ. 'no') ) THEN
C           DO NOT WANT OUTPUT FILE
            IRTFLG = 1
         ENDIF
         IF (IRTFLG .NE. 0 .OR. DOCNAM(1:1) .EQ. '*') RETURN

         CALL OPENDOC(DOCNAM,.TRUE.,NLET,LUNDOC,LUNDOCT,.FALSE.,'',
     &               .FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
         IF (IRTFLG .NE. 0)RETURN

         CALL LUNDOCPUTCOM(LUNDOCT,
     &              'KEY,    CLASS,      HEIGHT        ID',IRTFLG)

         DO  I = 1,NBR
            DLIST(1) = ICLASS(I)
            DLIST(2) = HIT(I) 
            DLIST(3) = ID(I) 

C           STORE THIS LINE IN DOC FILE
            CALL LUNDOCWRTDAT(LUNDOCT,I,DLIST,3,IRTFLG)
	 ENDDO

         CLOSE(LUNDOCT)
         IRTFLG = 0

         END



