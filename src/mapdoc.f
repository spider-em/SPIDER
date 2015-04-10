
C++*********************************************************************
C
C MAPDOC.F -- CREATED NOV 90 
C              OPFILEC                  FEB  03 ARDEAN LEITH
C **********************************************************************
C *  AUTHOR:  ArDean Leith
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
C    MAPDOC(IRTFL)
C
C    PURPOSE:       READS A DOC FILE FOR MAPPING AN IMAGE FILE INTO
C                   A NEW IMAGE FILE. DOC. FILE CONTAINS X,Y & Z 
C                   LOCATIONS OF A VOXEL IN THE IMAGE FILE.  ALL
C                   VOXELS IN THE FILE HAVING THE CLUSTER NUMBER
C                   FOUND AT THIS LOCATION ARE CARRIED INTO A NEW
C                   FILE.  ALL OTHER VOXELS ARE SET TO ZERO. 
C
C    PARAMETERS     IRTFLG       ERROR RETURN
C
C    CALLED BY:     
C
C        0         2         3         4         5         6         7     
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

       SUBROUTINE MAPDOC(IRTFLG)

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
       IF (NKEY .GT. MAXLUT) THEN
          ITEMP = MAXLUT
          CALL ERRT(102,'TOO MANY KEYS FOR TABLE SIZE',ITEMP)
          RETURN
       ELSEIF (NKEY .LE. 0) THEN
          CALL ERRT(101,'NO ENTRIES IN DOC. FILE',NE)
          RETURN
       ENDIF

20      MAXIM = 0
        CALL OPFILEC(0,.TRUE.,IMFILE,LUNIM,'O',IFORM,
     &     NSAM,NROW,NSLICE,MAXIM,'CLUSTER INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .EQ. -1) THEN
           CLOSE (LUND)
           GOTO 10
        ENDIF
        IF (IRTFLG .NE. 0) GOTO 999

30      MAXIM = 0
        CALL OPFILEC(LUNIM,.TRUE.,OUTFILE,LUNOUT,'U',IFORM,
     &     NSAM,NROW,NSLICE,MAXIM,'CLUSTER OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .EQ. -1) THEN
           CLOSE (LUND)
           CLOSE (LUNIM)
           GOTO 20
        ENDIF
        IF (IRTFLG .NE. 0) GOTO 999

C      ZERO THE LOOK-UP-TABLE
       DO N = 1,MAXLUT
          TABLE(N) = 0
       END DO

       IRECL   = 0
       NEWVAL  = 0
       MAPINGS = 0
       NOMAP   = 0
  
       DO N = 1,NKEY

C        CHECK TO SEE IF THIS LINE OF DOC FILE IS IN USE
         IVAL   = DBUF(1,N) 
         IF (IVAL .LE. 0) CYCLE

         MAPINGS = MAPINGS + 1

C        FIND LOCATION OF VOXEL TO BE MAPPED
         IX     = DBUF(2,N) + 1
         IY     = DBUF(3,N) + 1
         ISLICE = DBUF(4,N)

C        IGNORE ANY SLICES NOT IN THE FILE
         IF (ISLICE > NSLICE) CYCLE

C        FIND RECORD NUMBER FOR THIS VOXEL
         IREC = (ISLICE -1) * NROW + IY
         IF (IREC .NE. IRECL) THEN
C            MUST READ A NEW RECORD OF THE IMAGE FILE
             CALL REDLIN(LUNIM,BUF,NSAM,IREC)
             IRECL = IREC
         ENDIF

C        FIND VALUE OF THE VOXEL IN OLD FILE
         IGOT = BUF(IX)

         IF (IGOT > MAXLUT) THEN
C           OVERFLOW OF LUT TABLE POSSIBLE
            WRITE(NOUT,*) ' *** TABLE OVERFLOW IN DOCMAP'
            GOTO 999

         ELSEIF (IGOT .EQ. 0) THEN
C           UNMAPPABLE CLUSTER 
            NOMAP = NOMAP + 1

         ELSEIF (TABLE(IGOT) .EQ. 0.0) THEN
C           THIS CLUSTER NOT YET MAPPED TO ANYTHING

            NEWVAL = NEWVAL + 1
            TABLE(IGOT) = NEWVAL
         
         ENDIF
       ENDDO


c****************DEBUGGING
       DO I = 1,NEWVAL
         IF(TABLE(I) .NE. 0.0) THEN
           IT = TABLE(I)
           WRITE(NOUT,*) I,' ->',IT
         ENDIF
       END DO
C******************************

       WRITE(NOUT,96) MAPINGS,NEWVAL
96     FORMAT('  MAPS:',I7,' REMAPPED VALUES:',I7)

       WRITE(NOUT,97) NOMAP
97     FORMAT('  UNMAPPED VALUES:',I7)

       NREC2 = NROW * NSLICE
       CALL MAPIM(LUNIM,LUNOUT,NSAM,1,NREC2,TABLE,NDUM,BUF,IRTFLG)
      
       IRTFLG = 0 
       
999    CLOSE(LUNIM)
       CLOSE(LUNOUT)
       CLOSE(LUND)
         
       END



