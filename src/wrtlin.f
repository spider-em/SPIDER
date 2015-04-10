
C++*******************************************************************
C   WRTLIN.F                        IOSTAT ADDED             DEC 97 al
C                                   FILE NAME FOR ERROR      SEP 02 al
C                                   ENDEDNESS                FEB 03 al
C                                   ONLYONE_WRT              NOV 08 al
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
C     WRTLIN(LUNT,BUF,NB,NREC)
C
C     PURPOSE:  WRITE A LINE OF FLOATING POINT NUMBERS INTO FILE
C
C     PARAMETERS:
C        LUNT   LOGICAL UNIT NUMBER OF FILE BEING WRITTEN
C        BUF    BUFFER WHERE RECORD IS READ IN
C        NB     NUMBER OF VALUES IN RECORD TO BE READ
C        NREC   RECORD TO BE READ
C
C        IERR   ERROR CODE 1 IS RETURNED IN CASE OF ERROR
C               IERR IS DEFINED IN COMMON /IOERR/IERR
C 
C--*******************************************************************

      MODULE WRTLIN_PIPE_STUFF
         SAVE
         LOGICAL            :: IMGPIPE    = .FALSE.
         INTEGER, PARAMETER :: LUNIMGPIPE = 303
      END MODULE WRTLIN_PIPE_STUFF



      SUBROUTINE WRTLIN(LUNT,BUF,NB,NREC)

C     HAS  IMG PIPING INFO
      USE WRTLIN_PIPE_STUFF

C     USE INLINE BUFFER COMMON AREA
      INCLUDE 'INLN_INFO.INC'

      REAL            BUF(NB)
      COMMON /IOERR/  IERR
      COMMON /UNITS/  LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT
      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

C     AUTOMATIC ARRAY 
      REAL,DIMENSION(NB) :: FLIPBUF(NB)

      CHARACTER(LEN=80)  :: FILOPEND
      LOGICAL            :: ISOPEN

#ifdef USE_MPI
      LOGICAL            :: ONLYONE_RED,ONLYONE_WRT
      COMMON /COMM_MPI/ONLYONE_RED,ONLYONE_WRT
      include 'mpif.h'
      ICOMM   = MPI_COMM_WORLD
      CALL MPI_COMM_RANK(ICOMM, MYPID, MPIERR)
      IF (.NOT. ONLYONE_WRT) MYPID = -1
#else
      MYPID = -1 
#endif

      IF (ISINLINE(LUNT)) THEN
C        USE INLINED BUFFER FOR I/O (SEE OPENINLN.F)
         CALL INLN_WRTLIN(LUNT,BUF,NB,NREC)
         RETURN

      ELSEIF (IMGPIPE) THEN
C        USE PIPE FOR OUTPUT
         CALL WRTLIN_PIPE(BUF,NB,NREC,IRTFLG)
         RETURN
      ENDIF

      LUN = LUNARB(LUNT)

C     ADD LUNARA(LUN) (FOR LABEL REC) AND LUNSTK (FOR STACK OFFSET)
C     TO NREC TO GET THE CORRECT RECORD NUMBER
      I = NREC + LUNARA(LUN) + LUNSTK(LUNT)

      IERR = 0
      IF (MYPID .LE. 0) THEN
         IF (LUNFLIP(LUN) .EQ. 0) THEN
            WRITE(LUN,REC=I,IOSTAT=IERR) BUF
         ELSE
            CALL FLIPBYTES(BUF,FLIPBUF,NB,IRTFLG)
            WRITE(LUN,REC=I,IOSTAT=IERR) FLIPBUF
         ENDIF

         IF (IERR .NE. 0) THEN
            INQUIRE(UNIT=LUN,OPENED=ISOPEN,NAME=FILOPEND)
            LENT = lnblnkn(FILOPEND)
            WRITE(NOUT,99) IERR,I,NB,LUN,FILOPEND(:LENT)
99          FORMAT( '  *** ERROR(',I4,') WRITING RECORD: ',I6,
     &              ' LENGTH: ', I5,' ON UNIT: ',I3,' TO FILE: ',A)
         ENDIF
      ENDIF

#ifdef USE_MPI
      IF (ONLYONE_WRT) THEN
         CALL BCAST_MPI('WRTLIN','IERR',IERR,1,'I', ICOMM,MPIERR)
         CALL MPI_BARRIER(ICOMM,MPIERR)
      ENDIF
#endif

      END





C++*********************************************************************
C
C WRTLIN_OPENPIPE                   NEW           JAN 2006 ARDEAN LEITH
C
C **********************************************************************
C
C    WRTLIN_OPENPIPE(PIPENAME,IRTFLG)
C
C    PURPOSE:     OPENS PIPE FOR IMAGE/VOLUME    
C
C    PARAMETERS:  PIPENAME  PIPE NAME                           (SENT)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C--*******************************************************************

      SUBROUTINE WRTLIN_OPENPIPE(PIPENAME,IRTFLG)

      USE WRTLIN_PIPE_STUFF

      INCLUDE 'CMLIMIT.INC'
      CHARACTER(LEN=MAXNAM)    :: PIPENAME
      CHARACTER(LEN=MAXNAM+24) :: MSG

#ifdef sgi
C     SETS NAME FOR ASSIGN OBJECT FILE
      CALL SETENV('FILENV','jnkASSIGN1',IRTFLG)
C     CLEAR ANY EXISTING ASSIGN OBJECT FILE
      CALL ASNRM(IRTFLG)
C     INITIALIZE THE ASSIGN OBJECT FILE FOR WRITING
      MSG = '-s u -a ' // PIPENAME
      CALL ASNUNIT(LUNREGPIPE,MSG,IRTFLG)
      IF (IRTFLG .NE. 0) THEN
         CALL ERRT(102,'ASNUNIT TO SET PIPE, RETURNS:',IER)
         RETURN
      ENDIF
#endif

      OPEN(UNIT=LUNIMGPIPE, 
     &    FILE=PIPENAME,
     &    FORM='UNFORMATTED',
     &    ACCESS='SEQUENTIAL',
     &    STATUS='OLD',
     &    ACTION='WRITE',
     &    IOSTAT=IRTFLG)

      IF (IRTFLG .NE. 0) THEN 
         MSG = 'FAILED TO OPEN PIPE: ' // PIPENAME
         CALL ERRT(101,MSG,IRTFLG)
         RETURN
      ENDIF

      IMGPIPE = .TRUE.
      IRTFLG  = 0

      RETURN
      END


C++*********************************************************************
C
C WRTLIN_PIPE                         NEW           JAN 2006 ARDEAN LEITH
C
C **********************************************************************
C
C    WRTLIN_PIPE(BUF,LENT,NB,IRTFLG)
C
C    PURPOSE:    SENDS IMAGE DOWN LUNIMGPIPE   
C
C    PARAMETERS:  BUF   DATA BUFFER LINE                        (SENT)
C                 NB    LENGTH OF DATA BUFFER LINE              (SENT)
C                 NREC  NUMBER OF DATA BUFFER LINE              (SENT)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C--*********************************************************************

      SUBROUTINE WRTLIN_PIPE(BUF,NB,NREC,IRTFLG)

      USE WRTLIN_PIPE_STUFF

      REAL,DIMENSION(NB) :: BUF

      IF (IMGPIPE) THEN
C        WRITE IMAGE/VOL. LINE TO NAMED PIPE
         WRITE(LUNIMGPIPE,IOSTAT=IRTFLG) BUF
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(102,'PIPING IMAGE LINE',NREC)
            RETURN
         ENDIF
      ELSE
         CALL ERRT(101,'NO IMAGE PIPE OPEN',NE)
         IRTFLG = 1
      ENDIF
 
      RETURN
      END

C++*********************************************************************
C
C WRTLIN_PIPE_TOG                      NEW        JAN 2006 ARDEAN LEITH
C
C **********************************************************************
C
C    WRTLIN_PIPE_TOG()
C
C    PURPOSE:    SETS IMAGE PIPLINE TOGGLE    
C
C    PARAMETERS:  
C
C--*********************************************************************

      SUBROUTINE WRTLIN_PIPE_TOG()

      USE WRTLIN_PIPE_STUFF
      INCLUDE 'CMLIMIT.INC' 

      CHARACTER (LEN=MAXNAM), SAVE :: PIPENAME = CHAR(0)
      CHARACTER (LEN=1)            :: NULL = CHAR(0)
      LOGICAL                      :: ISOPEN = .FALSE.
      LOGICAL, SAVE                :: NEEDOPEN = .TRUE.

      INCLUDE 'CMBLOCK.INC'

      IF (IMGPIPE) THEN
C        WRITE IMAGE/VOL. LINE TO FILE(S) 
         WRITE(NOUT,*) '  IMAGES SENT TO FILES'
         IMGPIPE = .NOT. IMGPIPE

      ELSE
C        WRITE IMAGE/VOL. LINE TO NAMED PIPE

         IF (PIPENAME .EQ. NULL) THEN
C           NO FILE OPEN
            CALL FILERD(PIPENAME,NLET,NULL,'PIPE NAME~',IRTFLG)
            IF (IRTFLG .NE. 0) RETURN

            CALL WRTLIN_OPENPIPE(PIPENAME(1:NLET),IRTFLG)
            IF (IRTFLG .NE. 0) RETURN


         ELSE
C          SEE IF THIS PIPE EXISTS, (RETURNS ISOPEN, LUNOP)
           NLET = lnblnkn(PIPENAME)
           INQUIRE(FILE=PIPENAME(:NLET),OPENED=ISOPEN,NUMBER=LUNOP,
     &             IOSTAT=IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (.NOT. ISOPEN) THEN
C             MUST RE-OPEN PIPE
              CALL WRTLIN_OPENPIPE(PIPENAME(1:NLET),IRTFLG)
              IF (IRTFLG .NE. 0) RETURN
           ENDIF
           IMGPIPE = .TRUE.

           WRITE(NOUT,*)'  IMAGES PIPED TO: PIPENAME(1:NLET) ',
     &                  ' INSTEAD OF FILES'
         ENDIF
      ENDIF
 
      END
