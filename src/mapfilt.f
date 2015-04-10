
C++*********************************************************************
C
C MAPFILT.FOR             CREATED                  NOV  90 ARDEAN LEITH
C                         OPFILEC                  FEB  03 ARDEAN LEITH
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
C    MAPFILT(IRTFL)
C
C    PURPOSE:       READS A DOC. FILE CONTAINING CLUSTER NUMBER +1 (KEY)
C                   AND NUMBER OF VOXELS (REG 1) IN THE CLUSTER.  IF 
C                   NUMBER OF VOXELS < MIN OR > NMAX THEN THE VOXELS
C                   BELONGING TO THIS CLUSTER ARE SET TO ZERO.
C
C    PARAMETERS     IRTFLG       ERROR RETURN
C
C    CALLED BY:     UTIL6
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

       SUBROUTINE MAPFILT(IRTFLG)

       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC' 

       PARAMETER  (NSAMAX=16000)
       PARAMETER  (MAXREG=7)
       PARAMETER  (MAXKEY=16000)
       PARAMETER  (MAXLUT=16000)

       COMMON BUF(NSAMAX),TABLE(MAXLUT)

       COMMON /DOC_BUF/ DBUF(MAXREG,MAXKEY)

       REAL                     :: PLIST(7)
       LOGICAL                  :: USEMAX ,DEBUGGING
       CHARACTER(LEN=MAXNAM)    :: DOCNAM,IMFILE,OUTFILE
       CHARACTER *1             :: NULL = CHAR(0)

       DATA LUND/20/,LUNIM/21/,LUNOUT/22/

       IRTFLG     = 1
       DEBUGGING  = .FALSE.

10     CALL FILERD(DOCNAM,NLETD,DATEXC,'CLUSTER DOC',IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

       NLIST = 2
       IKEY  = 1
       ICALL = 0
       CALL UNSDAL(DOCNAM,ICALL,LUND,IKEY,PLIST,NLIST,DBUF,
     &             MAXKEY,MAXREG,NKEY,IERR)

       WRITE(NOUT,*) ' KEYS FOUND:', NKEY
       IF (NKEY > MAXLUT) THEN
          ITEMP = MAXLUT
          CALL ERRT(102,'TOO MANY KEYS FOR TABLE SIZE',ITEMP)
          RETURN
       ELSEIF (NKEY <= 0) THEN
          CALL ERRT(101,'NO ENTRIES IN DOC. FILE',NE)
          RETURN
       ENDIF

20     MAXIM = 0
       CALL OPFILEC(0,.TRUE.,IMFILE,LUNIM,'O',IFORM,NSAM,NROW,NSLICE,
     &     MAXIM,'CLUSTER INPUT',.FALSE.,IRTFLG)
       IF (IRTFLG .EQ. -1) THEN
           CLOSE (LUND)
           GOTO 10
       ENDIF
       IF (IRTFLG .NE. 0) GOTO 999

30     MAXIM = 0
       CALL OPFILEC(0,.TRUE.,OUTFILE,LUNOUT,'U',IFORM,NSAM,NROW,NSLICE,
     &     MAXIM,'CLUSTER OUTPUT',.FALSE.,IRTFLG)
       IF (IRTFLG .EQ. -1) THEN
           CLOSE (LUND)
           CLOSE (LUNIM)
           GOTO 20
       ENDIF
       IF (IRTFLG .NE. 0) GOTO 999

       USEMAX = .FALSE.
       NMAX   = 0
       CALL RDPRIS(NMIN,NMAX,NOT_USED,
     &   'STARTING AND ENDING CLUSTER SIZES RETAINED:',IRTFLG)
       IF (IRTFLG .EQ. -1) THEN
          CLOSE(LUNOUT)
          GOTO 30
       ENDIF
       IF (NMAX > NMIN) USEMAX = .TRUE.

C      SET THE LOOK-UP-TABLE SO THAT ALL VALUES ARE RETAINED (DEFAULT)
       DO ICLUS = 1,MAXLUT
          TABLE(ICLUS) = ICLUS 
       ENDDO

       NEWVAL  = 0
       MAPINGS = 0
       NOMAP   = 0
       NUNDER  = 0
       NOVER   = 0
 
       DO  100 ICLUSP1 = 1,NKEY

C        CHECK TO SEE IF THIS LINE OF DOC FILE IS IN USE
         ILIST   = DBUF(1,ICLUSP1) 
         IF (ILIST <= 0) GOTO 100

         MAPINGS = MAPINGS + 1

C        FIND NUMBER OF VOXELS IN CLUSTER
         INUM   = DBUF(2,ICLUSP1) 
         ICLUS  = ICLUSP1 - 1

         IF (INUM .LT. 0) THEN
C           IMPOSSIBLE NUMBER OF VALUES IN CLUSTER
            WRITE(NOUT,*) ' *** INUM :',INUM,' CAN NOT BE < 0'
            NOMAP = NOMAP + 1
            GOTO 999

         ELSEIF (INUM .LT. NMIN) THEN
C           UNDER MINIMUM VALUES NEEDED PER CLUSTER, DISCARD CLUSTER
            NUNDER = NUNDER + 1
            IF (ICLUS > 0) TABLE(ICLUS) = 0.0

         ELSEIF (USEMAX .AND. INUM > NMAX) THEN
C           OVER MAXIMUM VALUES ALLOWED PER CLUSTER, DISCARD CLUSTER
            NOVER = NOVER + 1
            IF (ICLUS > 0) TABLE(ICLUS) = 0.0
         ENDIF
100    CONTINUE

c****************DEBUGGING
       IF (DEBUGGING) THEN
          DO I = 1,NKEY
            IF( TABLE(I) .NE. 0.0) THEN
              IT = TABLE(I)
              WRITE(NOUT,*) I,' ->',IT
            ENDIF
          END DO
       ENDIF
C******************************

       WRITE(NOUT,96) MAPINGS
96     FORMAT(/,'  MAPS:',I7)

       WRITE(NOUT,97) NMIN,NUNDER
97     FORMAT('  DISCARDED VALUES <',I7,' =',I7)

       IF (USEMAX) THEN
         WRITE(NOUT,98) NMAX,NOVER
98       FORMAT('  DISCARDED VALUES >',I7,' =',I7)
       ENDIF

       NREC2 = NROW * NSLICE
       CALL MAPIM(LUNIM,LUNOUT,NSAM,1,NREC2,TABLE,NDUM,BUF,IRTFLG)
      
       IRTFLG = 0 
999    CLOSE(LUNIM)
       CLOSE(LUNOUT)
       CLOSE(LUND)
         
       END



