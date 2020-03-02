
C++*********************************************************************
C
C  DELETF.F        REMOVED & ALTERED FROM UTIL1    DEC 88 ArDean Leith
C                  VERBOSE                         MAR 02 ArDean Leith
C                  SGI LEAK ON INTERNAL FMT        AUG 02 ArDean Leith
C                  INDEXED STACK                   JAN 03 ArDean Leith
C                  OPFILEC                         FEB 03 ArDean Leith
C                  SPIREOUT                        JUL 05 ArDean Leith
C                  OPENINLN KIND                   OCT 10 ArDean Leith
C                  MPI HEADER NEEDED               MAR 11 ArDean Leith
C                  STACK@ ACCEPTED                 MAR 12 ArDean Leith
C                  FORMATTING                      SEP 13 ArDean Leith
C                  MRC SUPPORT                     SEP 19 ArDean Leith
C                  REMOVED GETOLDSTACK             FEB 20 ArDean Leith
C **********************************************************************
C=*                                                                    *
C=* Author: ArDean Leith                                               *                                                            *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2020  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C    PURPOSE:        DELETE SPIDER OR MRC FILE(S)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C **********************************************************************

        SUBROUTINE DELETF(FILNAM,LUN)

        IMPLICIT  NONE
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=*)             :: FILNAM
        INTEGER                      :: LUN

        CHARACTER(LEN=MAXNAM)        :: FILN
        CHARACTER(LEN=MAXNAM+20)     :: MESG
        CHARACTER(LEN=1)             :: NULL = CHAR(0)

        LOGICAL                      :: GOT_MRC_STK_IMAGE,IS_MRC
        LOGICAL                      :: GOT_SPI_STK_IMAGE,ISDIGI
        INTEGER                      :: ICOMM,MYPID,MPIERR,NLET,IRTFLG
        INTEGER                      :: IFILE,IDIG,IFIRST,NLETA,NOFIND
        INTEGER                      :: LOCAT,MAXIM,NX,NY,NZ,LASTFI,IGO
        INTEGER                      :: IMGNUM,MAXIMNEW,NUMBUF,IER,NNN
        INTEGER                      :: IDUM

        INTEGER, PARAMETER           :: I_8 = SELECTED_INT_KIND(12)
        INTEGER(KIND=I_8), PARAMETER :: ZERO_8 = 0

        LOGICAL                      :: ISMRCFILE    ! FUNCTION


        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID 

        IF (FCHAR(4:4) == 'A') THEN      
C          OPERATION: 'FI A' WORKS ON WHOLE SERIES

           IF (MYPID <= 0) WRITE(NOUT,903) 
903        FORMAT(/,' WARNING, YOU ARE DELETING WHOLE FILE SERIES!'/)

           CALL FILERD(FILNAM,NLET,NULL,'FIRST',IRTFLG)
           IF (IRTFLG .NE. 0 .OR. FILNAM(1:1) == '*') RETURN

           IF (ISMRCFILE(FILNAM)) THEN
             CALL ERRT(101,'OPERATION NOT SUPPORTED FOR MRC FILES',IDUM)
             RETURN
           ENDIF

C          MULTIPLE FILE DELETION
           CALL GETFILENUM(FILNAM(1:NLET),IFILE,IDIG,.FALSE.,IRTFLG)

C          IDIG IS NUMBER OF CONSECUTIVE DIGITS AT END OF THE FIRST FILE NAME
           LASTFI = 10**IDIG - 1
           IGO    = NLET - IDIG + 1

        ELSE
C          SINGLE FILE DELETION WANTED

C          GET FILE NAME FOR DELETION
           CALL FILERD(FILNAM,NLET,NULL,'DELETE',IRTFLG)
           IF (IRTFLG .NE. 0 .OR. FILNAM(1:1) == '*') RETURN

           IFILE  = 1
           LASTFI = 1
        ENDIF

        IS_MRC =  (ISMRCFILE(FILNAM))
        IF (IS_MRC) THEN
C          MRC FILE DELETION
           FILN   = FILNAM (1:NLET)
           NLETA  = NLET
        ELSE

C          ADD EXTENSION TO FILNAM --> FILN
           FILN = FILNAM(1:NLET) // '.' // DATEXC(1:3) 
           NLETA  = NLET + 4
        ENDIF
          
        IFIRST = IFILE

	NOFIND = 0

20      LOCAT  = INDEX(FILNAM,'@')    ! STACK INDICATOR 

C       SEE IF STACKED SPIDER IMAGE     
        GOT_MRC_STK_IMAGE = (IS_MRC .AND. 
     &                       LOCAT > 1 .AND. 
     &                       ISDIGI(FILNAM(LOCAT-1:LOCAT-1)))

C       SEE IF STACKED SPIDER IMAGE     
        GOT_SPI_STK_IMAGE = (.NOT. IS_MRC .AND. LOCAT > 0 .AND. 
     &                       LOCAT < NLET .AND. 
     &                       ISDIGI(FILNAM(LOCAT+1:LOCAT+1)))

        IF (GOT_MRC_STK_IMAGE) THEN
C          MRC STACKED IMAGE
           CALL ERRT(101,'CAN NOT DELETE AN IMAGE FROM MRC STACK',IDUM)
           IRTFLG = 1
           RETURN
        ENDIF

        

        IF (GOT_SPI_STK_IMAGE) THEN
C          DELETE IMAGE FROM IMAGE STACK
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'O',IFORM,NX,NY,NZ,
     &                MAXIM,' ',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0)  RETURN

C          UPDATE IMAGE IN-USE VARIABLE
           CALL LUNGETINUSE(LUN,IMGNUM,IRTFLG)
           CALL LUNSETINUSE(LUN,0,IRTFLG)
           CALL LUNWRTHED(LUN,NX,IMGNUM,IRTFLG)

C          FIND MAXIM IN OVERALL HEADER
           CALL LUNGETMAXIM(LUN,MAXIM,IRTFLG)
           !write(3,*) '  In deletf; imgnum,maxim:',imgnum,maxim


           IF (IMGNUM == MAXIM) THEN
C             DELETED IMAGE IS MAXIM, MUST FIND NEW MAXIM
              CALL FINDMAXIM(LUN,NX,MAXIM,MAXIMNEW,IRTFLG)

C             UPDATE MAXIM IN OVERALL HEADER
              CALL LUNSAVMAXIM(LUN,NX,MAXIMNEW,IRTFLG)
           ENDIF

           IF (VERBOSE .AND. MYPID <= 0) 
     &         WRITE(NOUT,*)'  DELETED STACKED IMAGE: ', FILNAM(1:NLET)

           IF (MYPID <= 0 .AND. USE_SPIRE) THEN
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

        ELSEIF (NLET > 1 .AND. LOCAT == NLET) THEN        
C          WHOLE SPIDER STACK DELETE WANTED, GET RID OF @
           FILN  = FILNAM(1:NLET-1) // '.' // DATEXC(1:3) // NULL
           NLETA = NLET + 3

           IF (MYPID <= 0) THEN
              OPEN(LUN,FILE=FILN(1:NLETA),STATUS='OLD',IOSTAT=IER)
              IF (IER == 0) CLOSE(LUN,STATUS='DELETE', IOSTAT=IER)
           ENDIF

           CALL BCAST_MPI('DELETF','IER',IER,1,'I', ICOMM)

        ELSE
C          DESTROY SIMPLE FILE  
        
           IF (MYPID <= 0) THEN
              OPEN(LUN,FILE=FILN(1:NLETA),STATUS='OLD',IOSTAT=IER)
              IF (IER == 0) CLOSE(LUN,STATUS='DELETE', IOSTAT=IER)
           ENDIF

           CALL BCAST_MPI('DELETF','IER',IER,1,'I', ICOMM)

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

           IF (IFILE == IFIRST) NOFIND = 10
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


C        -------------------- FINDMAXIM ----------------------------

         SUBROUTINE FINDMAXIM(LUN,NX,MAXIMOLD,MAXIMNEW,IRTFLG)

         IMPLICIT NONE

         INTEGER   :: LUN,NX,MAXIMOLD,MAXIMNEW,IRTFLG

         INTEGER   :: NDUM,IMGNUM,IMUSED
         LOGICAL   :: IS_MRC

C        PURPOSE:  FIND HIGHEST NUMBER IMAGE STILL IN THIS STACK

C        SET ERROR RETURN
         IRTFLG = 1

C        DETERMINE IF MRC OR SPIDER 
         CALL LUNGETIS_MRC(LUN,IS_MRC,IRTFLG)

         IF (IS_MRC) THEN
C          MRC STACK
           CALL ERRT(101,'THIS SUBROUTINE DOES NOT ACCEPT MRC',NDUM)
           IRTFLG = 1
           RETURN
         ENDIF

C        START WITH MAXIMOLD IMAGE
         IMGNUM = MAXIMOLD - 1

         DO WHILE (IMGNUM > 0)
C           GET NEXT IMAGE FROM STACK

C           GET SPECIFIED IMAGE HEADER FROM STACK FILE LOCATION
            CALL LUNREDHED(LUN,NX,0,     .FALSE.,IRTFLG)
            CALL LUNREDHED(LUN,NX,IMGNUM,.FALSE.,IRTFLG)
            !write(3,*) '  In findmax 1, imgnum,irtflg:',imgnum,irtflg


            IF (IRTFLG .NE. 0) EXIT  
 
C           NEED IMUSED FROM THIS STACKED IMAGE
            CALL LUNGETINUSE(LUN,IMUSED,IRTFLG)
            !write(3,*) '  In findmax 2, imgnum,irtflg:',imgnum,irtflg

            IF (IRTFLG == 0 .AND. IMUSED > 0 ) EXIT  !FOUND NEW MAXIM

C           TRY NEXT IMAGE DOWN
            IMGNUM = IMGNUM - 1
         ENDDO

         MAXIMNEW = IMGNUM
         IRTFLG   = 0

         END
