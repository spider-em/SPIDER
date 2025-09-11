
C++*********************************************************************
C   WRTLIN.F             IOSTAT ADDED              DEC 97 ArDean Leith
C                        FILE NAME FOR ERROR       SEP 02 ArDean Leith
C                        ENDEDNESS                 FEB 03 ArDean Leith
C                        ONLYONE_WRT               NOV 08 ArDean Leith
C                        MRC SUPPORT               JUL 19 ArDean Leith
C                        MORE MRC SUPPORT          DEC 19 ArDean Leith
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2019  Health Research Inc.,                         *
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
C     MRC_MODE  TYPE                           MY_NBYT
C         0       : INTEGER*1 (UNSIGNED BYTES)   -1 
C         1       : INTEGER*2 (SIGNED)            2 
C         2       : REALS                         4
C         6       : INTEGER*2 (UNSIGNED)         -2
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      MODULE WRTLIN_PIPE_STUFF
         SAVE
         LOGICAL            :: IMGPIPE    = .FALSE.
         INTEGER, PARAMETER :: LUNIMGPIPE = 303
      END MODULE WRTLIN_PIPE_STUFF


C     ---------------------- WRTLIN -------------------------------------

      SUBROUTINE WRTLIN(LUNT,BUF,NB,IREC)

C     HAS  IMG PIPING INFO
      USE WRTLIN_PIPE_STUFF

C     USE INLINE BUFFER COMMON AREA
      INCLUDE 'INLN_INFO.INC'
      INCLUDE 'CMLIMIT.INC'

      INTEGER         :: LUNT,NB,IREC
      REAL            :: BUF(NB)

      INTEGER         :: IERR,LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT
      INTEGER         :: LUNARA,LUNSTK,LUNARB,LUNFLIP
      COMMON /IOERR/  IERR
      COMMON /UNITS/  LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT
      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

      INTEGER *8      :: LUNMRCPOS
      INTEGER         :: LUNMRCNBYT
      COMMON /LUNMRC/    LUNMRCPOS(100),LUNMRCNBYT(100)

C     AUTOMATIC ARRAY 
      REAL                   :: FLIPBUF(NB)

      CHARACTER(LEN=MAXNAM)  :: FILOPEND
      LOGICAL                :: ISOPEN

C     INCLUSION FOR OPTIONAL MPI INITIALIZATION.  
      INTEGER            :: MYPID = -1
#include "MPI_INIT.INC"


      IF (ISINLINE(LUNT)) THEN    ! FROM: INLN_INFO.INC
C        USE INLINED BUFFER FOR I/O (SEE OPENINLN.F)
         CALL INLN_WRTLIN(LUNT,BUF,NB,IREC)
         RETURN

      ELSEIF (IMGPIPE) THEN
C        USE PIPE FOR OUTPUT
         CALL WRTLIN_PIPE(BUF,NB,IREC,IRTFLG)
         RETURN
      ENDIF

      LUN = LUNARB(LUNT)

      IF (LUNMRCPOS(LUN) == 0) THEN
C       DIRECT ACCESS SPIDER FILE

C       ADD LUNARA(LUN) (FOR LABEL REC) AND LUNSTK (FOR STACK OFFSET)
C       TO IREC TO GET THE CORRECT RECORD NUMBER
        I = IREC + LUNARA(LUN) + LUNSTK(LUNT)

        IERR = 0
        IF (MYPID <= 0) THEN
           IF (LUNFLIP(LUN) == 0) THEN
              WRITE(LUN,REC=I,IOSTAT=IERR) BUF
           ELSE
              CALL FLIPBYTESIII(BUF,FLIPBUF,NB,IRTFLG)
              WRITE(LUN,REC=I,IOSTAT=IERR) FLIPBUF
           ENDIF

           IF (IERR .NE. 0) THEN
              INQUIRE(UNIT=LUN,OPENED=ISOPEN,NAME=FILOPEND)
              LENT = lnblnkn(FILOPEND)
              WRITE(NOUT,99) IERR,I,NB,LUN,FILOPEND(:LENT)
99            FORMAT( '  *** ERROR(',I4,') WRITING RECORD: ',I6,
     &                ' LENGTH: ', I5,' ON UNIT: ',I3,' TO FILE: ',A)
           ENDIF
        ENDIF

      ELSE

C       STREAM ACCESS MRC FILE
        CALL WRTLIN_MRC(LUN,BUF,NB,IREC, MYPID,IERR) 

      ENDIF

#ifdef USE_MPI
      IF (ONLYONE_WRT) THEN
         CALL BCAST_MPI('WRTLIN','IERR',IERR,1,'I', ICOMM,MPIERR)
         CALL MPI_BARRIER(ICOMM,MPIERR)
      ENDIF
#endif

      END

C----------------------------- WRTLIN_MRC ----------------------------

      SUBROUTINE WRTLIN_MRC(LUN,BUF,NB,IREC, MYPID,IERR)

C     PURPOSE:  WRITE BUF TO STREAM ACCESS MRC IMAGE/VOLUME FILE

      IMPLICIT NONE

      INCLUDE 'CMLIMIT.INC'   ! NEED: NBUFSIZ

      INTEGER         :: LUN,NB,IREC
      REAL            :: BUF(NB)
      INTEGER         :: MYPID,IERR

      INTEGER *8      :: LUNMRCPOS
      INTEGER         :: LUNMRCNBYT
      COMMON /LUNMRC/    LUNMRCPOS(100),LUNMRCNBYT(100)
      INTEGER         :: LUNARA,LUNSTK,LUNARB,LUNFLIP
      COMMON /LUNARA/   LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

      INTEGER *1      :: I1BUF(NBUFSIZ)
      INTEGER *2      :: I2BUF(NBUFSIZ)

      INTEGER         :: I,IX,IRTFLG
      INTEGER         :: NBYT_PER_VAL,NBYT_PER_REC,MRCMODE,NX,NE
      INTEGER *8      :: IPOSMRC

C     LUNARA CARRIES NX (PIXELS IN X)
      NX = ABS(LUNARA(LUN))

C     DIFFERENT MRC MODES DIFFER IN DATA LENGTHS
      MRCMODE      = LUNMRCNBYT(LUN)   
      NBYT_PER_VAL = ABS(MRCMODE)
      NBYT_PER_REC = NBYT_PER_VAL * NX

C     FILE POSITION TO BEGIN WRITING DEPENDS ON IREC AND ORIGIN

C     LUNARA SIGN DENOTES IMAGE ORIGIN
      IF (LUNARA(LUN) > 0) THEN
C        ORIGIN IS LOWER LEFT
         IPOSMRC = LUNMRCPOS(LUN) - IREC * NBYT_PER_REC

      ELSEIF (LUNARA(LUN) < 0) THEN
C        ORIGIN IS UPPER LEFT
         IPOSMRC = LUNMRCPOS(LUN) + IREC * NBYT_PER_REC

      ELSE
         IRTFLG = 1
         RETURN
      ENDIF

      IF (NBYT_PER_VAL == 4) THEN
C        MRC FILE NEEDS 32 BIT REAL VALUES
         IF (MYPID <= 0) WRITE(LUN, POS=IPOSMRC,IOSTAT=IERR) BUF(1:NB)

      ELSEIF (NBYT_PER_VAL == 2) THEN      
C        MRC FILE NEEDS  16 BIT, INTEGER*2 VALUES        

         I2BUF(1:NB) = BUF(1:NB)   
         IF (MYPID <= 0) WRITE(LUN, POS=IPOSMRC,IOSTAT=IERR) I2BUF(1:NB)

      ELSEIF (NBYT_PER_VAL == 1) THEN         
C        MRC FILE NEEDS   8 BIT INTEGER*1 VALUES 

         I1BUF(1:NB) = BUF(1:NB)   
         IF (MYPID <= 0) WRITE(LUN,POS=IPOSMRC,IOSTAT=IERR)I1BUF(1:NB)
 
      ENDIF

      !if (irec==1) write(3,'(A,2X,i9,f8.1)')
      !&            '  In wrtlin, i1buf,buf:',i1buf(1),buf(1)
      ! if (irec==1) write(3,*)' i1buf(1, 1),(1, nx):',i1buf(1),i1buf(nx) 
      ! if (irec==64)write(3,*)' i1buf(64,1),(64,nx):',i1buf(1),i1buf(nx) 
      ! if (irec==63)  write(3,*)'  irec,lunmrcpos,ipos: ',
      ! &                          irec,lunmrcpos(lun),iposmrc
      ! if (irec==63)  write(3,*)'  irec,lunmrcpos,ipos: ',
      ! &                          irec,lunmrcpos(lun),iposmrc
c     CALL ERRT(101, '16 BIT MRC FILE OUTPUT IS NOT SUPPORTED YET',NE)
c     IERR = 1

      END



C++*********************************************************************
C
C WRTLIN_OPENPIPE                   NEW           JAN 2006 ArDean Leith
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

      END


C++*********************************************************************
C
C WRTLIN_PIPE                         NEW        JAN 2006 ArDean Leith
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
 
      END

C++*********************************************************************
C
C WRTLIN_PIPE_TOG                     NEW        JAN 2006 ArDean Leith
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
