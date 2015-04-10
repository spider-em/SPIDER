
C++*********************************************************************
C
C OPENFIL.F       ADAPTED FROM OPEN3.FOR            OCT 88 ArDean Leith
C                 LONG FILENAMES                    OCT 88 ArDean Leith
C                 MERGED WITH OPENF                 AUG 96 ArDean Leith
C                 -1, -3 & -7 FOURIER REMOVED       AUG 96 ArDean Leith
C                 F90 CHANGES                       APR 98 ArDean Leith
C                 USED LUNHDR                       JAN 99 ArDean Leith
C                 INDEXED STACKS                    JAN 03 ArDean Leith
C                 HEADER COPY                       FEB 03 ArDean Leith
C                 ENDEDNESS                         FEB 03 ArDean Leith
C                 REMOVED IRTFLG INPUT              APR 04 ARDEAN LEITH
C                 MPI OPEN BUG & BCAST CHANGES      OCT 08 ARDEAN LEITH
C                 OPENINLN INTEGER*8                OCT 10 ARDEAN LEITH
C                 NZ LIMITED TO 100000 BUG          JAN 15 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2015  Health Research Inc.,                         *
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
C    OPENFIL(LUNT,FILNAM,LUN,NX,NY,NZ,NSTACK,DISP,IRTFLG)
C
C    PURPOSE:       OPEN NEW OR EXISTING DATA FILE FOR RANDOM 
C                   ACCESS READING/WRITING.  CALLED FOR REGULAR FILES
C                   AND OVERALL STACK HEADERS.  THIS ROUTINE IS NOT
C                   CALLED FOR CREATING STACKED IMAGES WITHIN A
C                   STACK (ONLY FOR THE OVERALL STACK).  NEITHER
C                   IS IT USED FOR INLINE STACKS! 
C
C    PARAMETERS:
C        LUNT       UNIT TO COPY HEADER VALUES FROM               (SENT)
C        FILNAM     CHARACTER ARRAY, CONTAINING FILE NAME         (SENT)
C        LUN        LOGICAL UNIT NUMBER FOR FILNAM.               (SENT)
C        NX,NY      DIMENSIONS OF FILE                (SENT OR RETURNED)
C        NZ         NUMBER OF SLICES IN IMAGE         (SENT OR RETURNED)
C        NSTACK     STACK INDICATOR                    (SENT / RETURNED)
C                   ON INPUT: 
C                      FOR NEW FILE  0 : NOT A STACK, 
C                                   >0 : REGULAR STACK
C                                   <0 : INDEXED STACK, MAX TOTAL IMAGES
C                   ON OUTPUT: 
C                      -2  :  NOT A STACK 
C                      >=0 :  MAX IMAGE NUMBER NOW IN STACK.
C 
C        DISP       CHARACTER VARIABLE CONTAINING DISPOSITION     (SENT)
C                   SEE OPFIL.F FOR POSSIBLE VALUES            
C        KEEPEXT    TO KEEP EXTENSION                             (SENT)
C        IRTFLG     ERROR RETURN FLAG.                        (RETURNED)
C                     0 : NORMAL RETURN
C                     1 : ERROR RETURN
C                     5 : NOT A SPIDER INPUT IMAGE
C
C    CODING:   BASED ON IMAGE PARAMETERS NX,NY, & NZ,  A
C              NEW FILE IS OPENED WITH IREC RECORDS, EACH NX * 4 
C              BYTES LONG.  IREC ALLOWS SPACE FOR THE 3-D DISTRIBUTION,
C              PLUS NORMAL IMAGE HEADER.
C
C    VARIABLES: ITYPE (TYPE)  FILE TYPE SPECIFIER. 
C               +1    R     2-D IMAGE
C               +3    R3    3-D VOLUME FILE
C               -11   O2    2-D FOURIER TRANSFORM, MIXED RADIX ODD
C               -12   E2    2-D FOURIER TRANSFORM, MIXED RADIX EVEN
C               -21   O3    3-D FOURIER TRANSFORM, MIXED RADIX ODD
C               -22   E3    3-D FOURIER TRANSFORM, MIXED RADIX EVEN
C
C    CALL TREE:  SEE OPFILEC
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

	SUBROUTINE OPENFIL(LUNT,FILNAM,LUN,NX,NY,NZ,NSTACK,
     &                     ITYPE,DISP,KEEPEXT,IRTFLG)

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'LABLOCK.INC'

        COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100)

        CHARACTER(LEN=*)             :: FILNAM
        CHARACTER(LEN=MAXNAM)        :: FILNM
        CHARACTER(LEN=2*MAXNAM)      :: MSG
        CHARACTER(LEN=1)             :: NULL,DISP
        LOGICAL                      :: EX,CALLERRT,OPENED,KEEPEXT
        INTEGER, PARAMETER           :: I_8 = SELECTED_INT_KIND(12)
	INTEGER(KIND=I_8)            :: NWORDS8
        INTEGER(KIND=I_8), PARAMETER :: ZERO_8 = 0

#ifdef USE_MPI
        include 'mpif.h' 
        LOGICAL            :: ONLYONE_RED,ONLYONE_WRT
        COMMON /COMM_MPI/ONLYONE_RED,ONLYONE_WRT

        ICOMM = MPI_COMM_WORLD
        CALL MPI_COMM_RANK(ICOMM, MYPID, MPIERR)
        IF (.NOT. ONLYONE_RED) MYPID = -1  !KLUDGE FOR 'AP REF'
#else 
        MYPID = -1    
#endif    
        NULL = CHAR(0)

C       SET FLAG FOR ERRONEOUS RETURN
        IRTFLG = 1
       
C	CHECK IF USER WANTS TO USE THE IN-LINE BUFFER.
        INLNED = 0
        IF (FILNAM(1:1) == '_') THEN
           CALL INLNBUF(FILNAM,NLET,INLNED,IRTFLGT)
           IF (IRTFLGT .NE. 0)  RETURN

	   IF (DISP == 'U' .OR. DISP == 'N') THEN
	      EX = .FALSE.
	   ELSE
	      EX = .TRUE.
	   ENDIF
           FILNM = FILNAM(1:NLET)

        ELSE
C          USE REGULAR FILE
           INLNED = 0
	   LUNC   = LUN

C          NULLIFY THE INLINE POINTER (CLOSEINLN IS INSIDE OPENINLN.F)
           CALL CLOSEINLN(LUN,IRTFLGT)

           IF (KEEPEXT) THEN
              FILNM = FILNAM
              NLET  = lnblnk(FILNM)
           ELSE
              CALL FILNAMANDEXT(FILNAM,DATEXC,FILNM,NLET,.TRUE.,IRTFLGT)
	      IF (IRTFLGT .NE.0) RETURN
           ENDIF

C          SEE IF FILE EXISTS NOW
           IF (MYPID <= 0) THEN
              INQUIRE(FILE=FILNM(1:NLET),IOSTAT=IER,EXIST=EX)
           ENDIF

#ifdef USE_MPI
           IF (ONLYONE_RED) THEN
              CALL BCAST_MPI('OPENFIL','EX', EX,1,'L',ICOMM)
              CALL BCAST_MPI('OPENFIL','IER',IER,1,'I',ICOMM)
           ENDIF
#endif
           IF (IER .NE. 0) THEN
              MSG =  'OPENFIL; INQUIRY ERROR: ' // FILNM(1:NLET)
              CALL ERRT(101,MSG,NE)
              RETURN
           ENDIF
        ENDIF

C       SET ACTUAL LUN FOR THIS FILE
        LUNARB(LUN) = LUN

        CALLERRT = (DISP(1:1) .NE. 'Z' .AND. DISP(1:1) .NE. 'E')

10	IF (DISP == 'U' .OR. DISP == 'N') THEN
C         USER WANTS A NEW FILE TO WRITE INTO ---------------------- NEW

          IF (EX) THEN
C           FILE IS IS OPENED WITH 'UNKNOWN' BUT ALREADY EXISTS; IT WILL
C           BE REPLACED.  OLD FILES ARE DELETED FIRST 
            IF (MYPID <= 0) THEN
               OPEN(UNIT=LUN,FILE=FILNM(1:NLET),STATUS='OLD')
               CLOSE(UNIT=LUN,STATUS='DELETE',IOSTAT=IER)
            ENDIF
#ifdef USE_MPI
            IF (ONLYONE_RED)
     &         CALL BCAST_MPI('OPENFIL','IER',IER,1,'I',ICOMM)
#endif

	    IF (IER .NE. 0) THEN
               MSG =  'OPENFIL; DELETING FILE: ' // FILNM(1:NLET)
               CALL ERRT(101,MSG,NE)
               RETURN
            ENDIF
          ENDIF

C         INITIALIZE HEADER OBJECT FOR NEW IMAGE 
          CALL LUNSETHDR(LUNT,LUN,NX,NY,NZ,ITYPE,NSTACK,IRTFLGT)

C         SET FMIN.... TO AVOID FLT. PT. ERROR ON DEC
          CALL LUNSETSTAT(LUN,0,0.0,0.0,0.0,0.0,IRTFLGT)

	  IF (INLNED == 0) THEN
C	     REGULAR FILE, NOT INLINE FILE

             LENREC = LENOPENFILE(NX*4)
             IF (MYPID <= 0) THEN
                OPEN(LUN,FILE=FILNM(1:NLET),STATUS='NEW',
     &               FORM='UNFORMATTED',
     &               ACCESS='DIRECT',IOSTAT=IER,RECL=LENREC)
             ENDIF
#ifdef USE_MPI
             IF (ONLYONE_RED) THEN
                CALL BCAST_MPI('OPENFIL','IER',IER,1,'I',ICOMM)
             ENDIF
#endif
             IF (IER .NE. 0) THEN        
                MSG = 'OPENFIL; OPENING FILE: ' // FILNM(1:NLET)
                CALL ERRT(101,MSG,NE)
                RETURN
             ENDIF

C            FOR USING NEW FILE ON OLDER SPIDER RELEASES
             CALL LUNSET25(LUN,-1,IRTFLGT)

          ELSE
C            SET UP INLINE BUFFER AND TIE IT TO LUN

C            SAVE ISTACK IN NON-VOLATILE PART OF HEADER OBJECT
             CALL LUNGETISTACK(LUN,ISTACK,IRTFLGT)
             CALL LUNSETSTKALL(LUN,ISTACK,IRTFLGT)

             CALL LUNGETLAB(LUN,LABREC,INDXREC,NRECS,LABBYT,
     &                      LENBYT,IRTFLGT)
             IF (IRTFLGT .NE. 0) THEN
                write(6,*) 'nrecs,indxrec,lenbyt,errorflag:'
                write(6,*)  nrecs,indxrec,lenbyt,irtflgt
                RETURN
             ENDIF

             NWORDS8 = (NRECS + INDXREC) * LENBYT / 4

             CALL OPENINLN(LUN,INLNED,.TRUE.,NX,NWORDS8,
     &                     .TRUE.,IRTFLGT)
             IF (IRTFLGT .NE. 0)  RETURN
	  ENDIF

C         PUSH HEADER OBJECT INFO INTO NEW FILE
          CALL LUNWRTHED(LUN,NX,0,IRTFLGT)

          GOTO 2000

C         --------------------------------------------------------- OLD

	ELSEIF (DISP == 'O' .OR. DISP == 'K' .OR.
     &          DISP == 'Z' .OR. 
     &          DISP == 'E' .OR. DISP == 'M') THEN

C         FILE EXISTS, AND IS ACCESSED WITH 'OLD', OPEN THE FILE

	  IF (.NOT. EX) THEN
C            ERROR -- FILE DOES NOT EXIST, BUT BEING OPENED WITH 'OLD
             IF (MYPID <= 0)
     &             WRITE(NOUT,*) '*** FILE NOT FOUND: ',FILNM(1:NLET)

#ifdef USE_MPI
c            write(6,*) ' openfil; ',mypid, ' not found: ',filnm(1:nlet)
#endif

C	     FOR DISP=Z, DO NOT STOP THE BATCH JOB BY CALLING ERRT
             IF (CALLERRT) THEN
                MSG = 'OPENFIL; FILE NOT FOUND: ' // FILNM(1:NLET)
                CALL ERRT(101,MSG,NE)
             ENDIF
             RETURN
          ENDIF

          IF (INLNED .NE. 0) THEN
C            USE EXISTING INLINE BUFFER, TIE IT TO LUN & RETRIEVE NX
             CALL OPENINLN(LUN,INLNED,.FALSE.,NX,ZERO_8,
     &                     CALLERRT,IRTFLGT)
             IF (IRTFLGT .NE. 0)  RETURN
          ELSE
C            REGULAR FILE ACCESS

             LENREC = LENOPENFILE(256*4)
             IF (MYPID <= 0) THEN
                OPEN(LUN,FILE=FILNM(1:NLET),STATUS='OLD',
     &               ACCESS='DIRECT',
     &              FORM='UNFORMATTED',RECL=LENREC,IOSTAT=IER)
             ENDIF
#ifdef USE_MPI
             IF (ONLYONE_RED) THEN
                CALL BCAST_MPI('OPENFIL','IER',IER,1,'I',ICOMM)
             ENDIF
#endif

             IF (IER .NE. 0) THEN        
C               SEE IF FILE ALREADY OPEN ON DIFFERENT LUN
                IF (MYPID <= 0) THEN
                   INQUIRE(FILE=FILNM(1:NLET),IOSTAT=IER,OPENED=OPENED,
     &                     NUMBER=LUNOPENED)
                ENDIF
#ifdef USE_MPI
                IF (ONLYONE_RED) THEN
                  CALL BCAST_MPI('OPENFIL','IER',IER,1,'I',ICOMM)
                  CALL BCAST_MPI('OPENFIL','OPENED',OPENED,1,'L',ICOMM)
                  CALL BCAST_MPI('OPENFIL','LUNOPENED',LUNOPENED,1,
     &                           'I',ICOMM)
                ENDIF
#endif
                IF (OPENED) THEN
C                  MUST REDIRECT OLD LUNS TO CURRENT LUN, 
                   IF (MYPID <= 0)  CLOSE(LUNOPENED)
                   DO I = 1,100
                      IF (LUNARB(I) == LUNOPENED) LUNARB(I) = LUN
                   ENDDO

C                  NOW WE CAN REOPEN THIS SAME FILE ON LUN
                   IF (MYPID <= 0) THEN
                      OPEN(LUN,FILE=FILNM(1:NLET),STATUS='OLD',
     &                     ACCESS='DIRECT',
     &                     FORM='UNFORMATTED',RECL=LENREC,IOSTAT=IER)
                   ENDIF
#ifdef USE_MPI
                   IF (ONLYONE_RED) THEN
                      CALL BCAST_MPI('OPENFIL','IER',IER,1,'I',ICOMM)
                   ENDIF
#endif
                ELSE
C                  UNKNOWN ERROR
                   MSG = 'OPENFIL; OPENING FILE: ' // FILNM(1:NLET)
                   CALL ERRT(101,MSG,NE)
                   RETURN
                ENDIF
             ENDIF
     	  ENDIF

C         READ OVERALL HEADER FROM FILE 
          CALL LUNREDHED(LUN,256,0,.TRUE.,IRTFLGT)
          IF (IRTFLGT .NE. 0 .AND. MYPID <= 0) 
     &       WRITE(NOUT,*) '*** ERROR READING HEADER OF: ',FILNM(:NLET)
          IF (IRTFLGT .NE. 0) RETURN

C         NEED ITYPE
          CALL LUNGETTYPE(LUN,ITYPE,IRTFLGT)

          IF (ITYPE == -1 .OR. ITYPE == -3 .OR. ITYPE == -7) THEN
C           WANT TO READ OBSOLETE FORMAT FOURIER FILE
            WRITE(NOUT,96)
96          FORMAT(' *** CAN NOT READ OBSOLETE FOURIER FORMAT',/,
     &          '*** CONVERT FOURIER FILE TO REAL FORMAT USING ',
     &          'ORIGINAL VERSION OF SPIDER.'/)
             CALL ERRT(100, 'OPENFIL',NE)
             RETURN
          ENDIF

C         NEED NX VALUE 
          CALL LUNGETSIZE(LUN,NX,NY,NZ,IRTFLGT)
          IF (IRTFLGT .NE. 0) RETURN

          !print *, 'READING HEADER OF: ',FILNM(:NLET)
          !print *, 'NX...: ',LUN,NX,NY,NZ,itype

          IF (ITYPE == 0 .OR. 
     &       NX <= 0     .OR. NY <= 0     .OR. NZ <= 0 .OR.
     &       NX > 10000000 .OR. NY > 10000000 .OR. 
     &       NZ > 10000000) THEN
C           NOT A NATIVE SPIDER IMAGE

C           FLIP BYTES IN HEADER OBJECT 
            CALL LUNFLIPHDR(LUN,IRTFLGT)
            CALL LUNSETFLIP(LUN,1,IRTFLG)

C           NEED ITYPE & SIZE
            CALL LUNGETTYPE(LUN,ITYPE,IRTFLGT)
            CALL LUNGETSIZE(LUN,NX,NY,NZ,IRTFLGT)
            IF (IRTFLGT .NE. 0) RETURN

            IF (ITYPE == 0 ) THEN
C              PROBABLY NOT A SPIDER IMAGE
               IF (CALLERRT) THEN
                 CALL ERRT(102,'INVALID HEADER, BAD FILE TYPE',ITYPE)
               ENDIF
               IRTFLG = 5   ! RETURN ERROR FLAG FOR NON-SPIDER IMAGE

            ELSEIF ( NY <= 0 ) THEN
C              PROBABLY NOT A SPIDER IMAGE
               IF (CALLERRT) THEN
                 CALL ERRT(101,'INVALID HEADER, NX',NX)
               ENDIF
               IRTFLG = 5   ! RETURN ERROR FLAG FOR NON-SPIDER IMAGE

            ELSEIF ( NY <= 0 ) THEN
C              PROBABLY NOT A SPIDER IMAGE
               IF (CALLERRT) THEN
                 CALL ERRT(101,'INVALID HEADER, NY',NY)
               ENDIF
               IRTFLG = 5   ! RETURN ERROR FLAG FOR NON-SPIDER IMAGE

            ELSEIF ( NZ <= 0 ) THEN
C              PROBABLY NOT A SPIDER IMAGE
               IF (CALLERRT) THEN
                 CALL ERRT(101,'INVALID HEADER, NZ',NX)
               ENDIF
               IRTFLG = 5   ! RETURN ERROR FLAG FOR NON-SPIDER IMAGE

             ELSEIF (NX > 1000000) THEN
C              PROBABLY NOT A SPIDER IMAGE
               IF (CALLERRT) THEN
                 CALL ERRT(102,
     &              'INVALID HEADER, NX (LIMIT=10000000)',NX)
               ENDIF
               IRTFLG = 5   ! RETURN ERROR FLAG FOR NON-SPIDER IMAGE

             ELSEIF (NY > 1000000) THEN
C              PROBABLY NOT A SPIDER IMAGE
               IF (CALLERRT) THEN
                 CALL ERRT(102,
     &              'INVALID HEADER, NY (LIMIT=10000000)',NY)
               ENDIF
               IRTFLG = 5   !RETURN ERROR FLAG FOR NON-SPIDER IMAGE

             ELSEIF (NZ > 1000000) THEN
C              PROBABLY NOT A SPIDER IMAGE
               IF (CALLERRT) THEN
                 CALL ERRT(102,
     &              'INVALID HEADER, NZ (LIMIT=10000000)',NZ)
               ENDIF
               IRTFLG = 5   ! RETURN ERROR FLAG FOR NON-SPIDER IMAGE
            ENDIF

            IF (IRTFLG == 5) THEN
C              RETURN ERROR FLAG FOR NON-SPIDER IMAGE
               FILNAM = FILNM(1:NLET)
               RETURN
            ENDIF 

            IF (VERBOSE) 
     &         WRITE(NOUT,*) ' NON-NATIVE BYTE ORDERED SPIDER FILE'
          ENDIF


C         ON UNIX, REOPEN THE FILE WITH FINAL RECORD LENGTH
          IF (INLNED == 0) THEN
             CLOSE(LUN)
             LENREC = LENOPENFILE(NX*4)
             IF (MYPID <= 0) THEN
  	        OPEN(LUN,FILE=FILNM(1:NLET),STATUS='OLD',
     &               FORM='UNFORMATTED',
     &               ACCESS='DIRECT',IOSTAT=IER,RECL=LENREC)
             ENDIF
#ifdef USE_MPI
             IF (ONLYONE_RED) THEN
                CALL BCAST_MPI('OPENFIL','IER',IER,1,'I',ICOMM)
             ENDIF
#endif
             IF (IER .NE. 0) THEN        
                MSG = 'OPENFIL; REOPENING FILE: ' // FILNM(1:NLET)
                CALL ERRT(101,MSG,NE)
                RETURN
             ENDIF
          ENDIF

        ELSE
           MSG = 'OPENFIL; PGM. ERROR, UNKNOWN DISP: ' // DISP
           CALL ERRT(101,MSG,NE)
           RETURN
        ENDIF

2000	CONTINUE

C       SAVE ISTACK IN NON-VOLATILE PART OF HEADER OBJECT
        CALL LUNGETISTACK(LUN,ISTACK,IRTFLGT)
        CALL LUNSETSTKALL(LUN,ISTACK,IRTFLGT)

C       SAVE MAXIM IN NON-VOLATILE PART OF HEADER OBJECT
        CALL LUNCOPYMAXIM(LUN,MAXIM,IRTFLGT)

C       SET FINAL LUNARA OFFSET VALUE FOR USE BY REDLIN/WRTLIN.
C       FOR INDEXED STACKS THIS INCLUDES INDEX OFFSET
        IF (ISTACK .NE. 0) THEN
           CALL LUNSETIMGOFF(LUN,0,NX,IRTFLGT)
        ELSE
           IVAL = 1
           CALL LUNSETIMGOFF(LUN,IVAL,NX,IRTFLGT)
        ENDIF

C       PUT COMMON VALUES INTO COMMON AREA (NOT NEEDED IN FUTURE?)
        CALL LUNSETCOMMON(LUN,IRTFLGT)

C       WRITE OUT FILE OPENING INFO
        CALL LUNSAYINFO(LUN,IRTFLGT)
      
C       NEED TO RETURN NSTACK & MAXIM 
        IF (ISTACK .NE. 0) THEN
C          RETURN CURRENT HIGHEST IMAGE NUMBER IN NSTACK
           NSTACK = MAXIM
        ELSE
C          NOT A STACK, RETURN -2
           NSTACK = -2
        ENDIF

C       SET FLAG FOR NORMAL RETURN	
        IRTFLG = 0
        RETURN

	END
