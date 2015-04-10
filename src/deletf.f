
C++*********************************************************************
C
C    DELETF.F       REMOVED & ALTERED FROM UTIL1    DEC 88 ArDean Leith
C                   VERBOSE                         MAR 02 ArDean Leith
C                   SGI LEAK ON INTERNAL FMT        AUG 02 ArDean Leith
C                   INDEXED STACK                   JAN 03 ARDEAN LEITH
C                   OPFILEC                         FEB 03 ARDEAN LEITH
C                   SPIREOUT                        JUL 05 ARDEAN LEITH
C                   OPENINLN KIND                   OCT 10 ARDEAN LEITH
C                   MPI HEADER NEEDED               MAR 11 ARDEAN LEITH
C                   STACK@ ACCEPTED                 MAR 12 ARDEAN LEITH
C                   FORMATTING                      SEP 13 ARDEAN LEITH
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
C    DELETF(FILNAM,LUN)
C
C    PARAMETERS:     FILNAM    CHAR. VARIABLE FOR FILENAME (EMPTY)
C                    LUN       UNIT FOR FILE OPENING
C
C    PURPOSE:        DELETE SPIDER FILE(S)
C
C **********************************************************************

        SUBROUTINE DELETF(FILNAM,LUN)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=*)             :: FILNAM
        CHARACTER(LEN=MAXNAM)        :: FILN
        CHARACTER(LEN=MAXNAM+20)     :: MESG
        CHARACTER(LEN=1)             :: NULL

        LOGICAL                      :: GOT_IMAGE,ISDIGI

        INTEGER, PARAMETER           :: I_8 = SELECTED_INT_KIND(12)
        INTEGER(KIND=I_8), PARAMETER :: ZERO_8 = 0

#ifdef USE_MPI
        include 'mpif.h'
#endif
        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID 

        NULL = CHAR(0)

200     IF (FCHAR(4:4) == 'A') THEN
           IF (MYPID <= 0) WRITE(NOUT,903) 
903        FORMAT(/,' WARNING, YOU ARE DELETING WHOLE FILE SERIES!'/)
           CALL FILERD(FILNAM,NLET,NULL,'FIRST',IRTFLG)
           IF (IRTFLG .NE. 0 .OR. FILNAM(1:1) .EQ. '*') RETURN

C          MULTIPLE FILE DELETION
           CALL GETFILENUM(FILNAM(1:NLET),IFILE,IDIG,.FALSE.,IRTFLG)

C          IDIG IS NUMBER OF CONSECUTIVE DIGITS AT END OF THE FIRST FILE NAME
           LASTFI = 10**IDIG - 1
           IGO    = NLET - IDIG + 1

        ELSE
           CALL FILERD(FILNAM,NLET,NULL,'DELETE',IRTFLG)
           IF (IRTFLG .NE. 0 .OR. FILNAM(1:1) == '*') RETURN
C          SINGLE FILE DELETION WANTED
           IFILE  = 1
           LASTFI = 1
        ENDIF

        IFIRST = IFILE

        FILN   = FILNAM(1:NLET) // '.' // DATEXC(1:3) // NULL
        NLETA  = NLET + 4

	NOFIND = 0

20      LOCAT     = INDEX(FILNAM,'@')      
        GOT_IMAGE = (LOCAT .GT. 0 .AND. 
     &               LOCAT .LT. NLET .AND. 
     &               ISDIGI(FILNAM(LOCAT+1:LOCAT+1)))

        IF (LOCAT > 0 .AND. GOT_IMAGE) THEN
C          DELETE IMAGE IN AN IMAGE STACK
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'O',IFORM,NSAM,NROW,NSLICE,
     &                MAXIM,' ',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0)  RETURN

C          UPDATE IMAGE IN-USE VARIABLE
           CALL LUNGETINUSE(LUN,IMGNUM,IRTFLG)
           CALL LUNSETINUSE(LUN,0,IRTFLG)
           CALL LUNWRTHED(LUN,NSAM,IMGNUM,IRTFLG)

C          FIND MAXIM IN OVERALL HEADER
           CALL LUNGETMAXIM(LUN,MAXIM,IRTFLG)

           IF (IMGNUM .EQ. MAXIM) THEN
C             DELETED IMAGE IS MAXIM, MUST FIND NEW MAXIM
              CALL FINDMAXIM(LUN,NSAM,MAXIM,MAXIMNEW,IRTFLG)

C             UPDATE MAXIM IN OVERALL HEADER
              CALL LUNSAVMAXIM(LUN,NSAM,MAXIMNEW,IRTFLG)
           ENDIF

           IF (VERBOSE .AND. MYPID .LE. 0) 
     &         WRITE(NOUT,*)'  DELETED STACKED IMAGE: ', FILNAM(1:NLET)

           IF (MYPID .LE. 0 .AND. USE_SPIRE) THEN
              MESG = '.DELETE FILE: ' // FILNAM(1:NLET)  
              CALL SPIREOUT(MESG,IRTFLG)
           ENDIF

           RETURN

        ELSEIF (FILNAM(1:1) == '_') THEN
C          INLINED BUFFER DELETE WANTED, GET BUFFER NUMBER
           CALL INLNBUF(FILNAM,NLET,NUMBUF,IRTFLG)

C          FREE UP THE INLINE BUFFER (DEALLOCATES IF NECESSARY)
           CALL OPENINLN(LUN,NUMBUF,.TRUE.,0,ZERO_8,.FALSE.,IRTFLG)

           IF (MYPID <= 0 .AND. VERBOSE .AND. IRTFLG == 0) THEN
              WRITE(NOUT,*) ' DELETED INLINE FILE: ',FILNAM(1:NLET)
           ENDIF
           RETURN

        ELSEIF(LOCAT == NLET .AND. NLET > 1) THEN        
C          WHOLE STACK DELETE WANTED, GET RID OF @
           FILN  = FILNAM(1:NLET-1) // '.' // DATEXC(1:3) // NULL
           NLETA = NLET + 3

           IF (MYPID <= 0) THEN
              OPEN(LUN,FILE=FILN(1:NLETA),STATUS='OLD',IOSTAT=IER)
              IF (IER == 0) CLOSE(LUN,STATUS='DELETE',IOSTAT=IER)
           ENDIF
#ifdef USE_MPI
           CALL MPI_BCAST(IER, 1, MPI_INTEGER, 0, ICOMM, IERR)
#endif

        ELSE
C          DESTROY SIMPLE FILE  
        
           IF (MYPID <= 0) THEN
              OPEN(LUN,FILE=FILN(1:NLETA),STATUS='OLD',IOSTAT=IER)
              IF (IER == 0) CLOSE(LUN,STATUS='DELETE',IOSTAT=IER)
           ENDIF

#ifdef USE_MPI
           CALL MPI_BCAST(IER, 1, MPI_INTEGER, 0, ICOMM, IERR)
#endif
        ENDIF

        IF (IER == 0) THEN
           IF (MYPID <= 0) THEN 
              IF (VERBOSE) WRITE(NOUT,*)' DELETED:      ',FILN(1:NLETA)
              IF (USE_SPIRE) THEN
                 MESG = '.DELETE FILE: ' // FILN(1:NLETA)  
                 CALL SPIREOUT(MESG,IRTFLG)
              ENDIF
           ENDIF
           NOFIND = 0
        
        ELSE
           IF (MYPID <= 0) THEN
              WRITE(NOUT,*) ' NO SUCH FILE: ',FILN(1:NLETA)
           ENDIF

           IF (IFILE .EQ. IFIRST) NOFIND = 10
           NOFIND = NOFIND + 1
        ENDIF

        IF (IFILE < LASTFI .AND. NOFIND < 10) THEN
C          DELETE NEXT FILE IN SERIES

C          CREATE NEXT FILE NAME
           IFILE = IFILE + 1
           CALL INTTOCHAR(IFILE,FILN(IGO:NLET),NNN,IDIG)
           GOTO 20
        ENDIF

999     IF (VERBOSE .AND. MYPID <= 0) WRITE(NOUT,*) ' '

        END




         SUBROUTINE FINDMAXIM(LUN,NSAM,MAXIMOLD,MAXIMNEW,IRTFLG)

C        FIND HIGHEST NUMBER IMAGE STILL IN THIS STACK

C        SET ERROR RETURN
         IRTFLG = 1

C        MUST SET ISBARE FOR GETOLDSTACK TO WORK
         CALL LUNSETISBARE(LUN,.TRUE.,IRTFLG)

C        START WITH MAXIMOLD IMAGE
         IMGNUM = MAXIMOLD - 1

         DO WHILE (IMGNUM .GT. 0)
C           GET NEXT IMAGE FROM STACK
            CALL GETOLDSTACK(LUN,NSAM,IMGNUM,
     &                      .FALSE.,.FALSE.,.FALSE.,IRTFLGT)
            IF (IRTFLGT .GT. 0) GOTO 999

            IF (IRTFLGT .EQ. 0) EXIT  !FOUND NEW MAXIM

C           TRY NEXT IMAGE DOWN
            IMGNUM = IMGNUM - 1
         ENDDO

100      MAXIMNEW = IMGNUM
         IRTFLG   = 0

999      CONTINUE
C        MUST SET ISBARE FOR GETOLDSTACK TO WORK
         CALL LUNSETISBARE(LUN,.FALSE.,IRTFLG)
         RETURN
         END
