
C ++********************************************************************
C
C COPYTOMRC    MODIFIED FROM COPYMRC               FEB 02 ArDean Leith  
C              ISSWAB ADDED                        JUL 02 ArDean Leith
C              FLIP QUESTION                       MAR 03 ArDean Leith
C              BAD IRECMRC4 & FLIP                 SEP 03 ArDean Leith
C              SCALING                             JAN 05 ArDean Leith
C              I*8                                 SEP 08 ArDean Leith
C              NPIX8                               DEC 08 ArDean Leith
C              BOTLEFT OPTION                      MAY 12 ArDean Leith
C              STREAM IO                           FEB 13 ArDean Leith
C              VOL BUG                             JUN 13 ArDean Leith
C              VOL BUG FIXED                       JUL 13 ArDean Leith
C              MODE 6 STACK SUPPORT                SEP 14 ArDean Leith
C              IPOSMRC INTEGER *8                  JAN 15 ArDean Leith
C              BOTLEFT DEFAULT                     JUL 15 ArDean Leith
C              2015 STACK SUPPORT                  JUL 15 ArDean Leith
C              COPYTOMRC_STK ADDED, STACK BUGGY    JUN 16 ArDean Leith
C              BOTLEFT OPTION REINTRODUCED         SEP 16 ArDean Leith
C              CREATES LITTLE ENDED FILES ONLY     JAN 18 ArDean Leith
C              GETLAB PARAMETERS                   NOV 19 ArDean Leith
C              REWRITE                             DEC 19 ArDean Leith
C              FIXED STACKS                        OCT 25 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2025  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email:                                                             *
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
C COPYTOMRC    (LUNSPI,LUNMRC,NX,NY,NZ)
C COPYTOMRC_STK(LUNSPI,LUNMRC,NX,NY,NZ)    !!! NEVER USED
C                                                                    
C PURPOSE: COPY FROM SPIDER TO MRC FILE FORMAT
C
C NOTES: DATA IN MRC FILE
C        MODE   TYPES OF PIXELS IN IMAGE
C               0 : INTEGER*1 (UNSIGNED BYTES) 
C               1 : INTEGER*2 (SIGNED) 
C               2 : REALS
C               6 : INTEGER*2 (UNSIGNED)
C
C  CALL TREE:
C  IF ('CP....')   
C     |
C  COPY1
C     |
C  IF ('CP TO MRC')    SPIDER FILE(S) TO MRC IMAGE(S) OR VOLUME(S)
C     |
C  COPYTOMRC
C     |
C     OPEN FIRST SPIDER INPUT FILE
C     OPFILES(0,LUNSPI,LUNDOC,LUNXM1,  
C     &       .TRUE.,FILNAM1,NLET1, 'E',
C     &       IFORM1,NX1,NY1,NZ1,NSTACK1,
C     &       'SPIDER INPUT FILE NAME OR TEMPLATE (e.g. SPI_STK@*)~~9'
C     &       TRUE., ILIST1,NILMAX, 
C     &       NOT_USED,NGOT1,IMG1, IRTFLG) 
C     |
C     FILERD (FILNAM2,NLET2,NULL,  
C     &       'MRC OUTPUT FILE NAME OR TEMPLATE (e.g. MRC_STK@*.mrc)~'
C     &       IRTFLG)
C     OPEN FIRST MRC OUTPUT FILE 
C     OPFILES(LUNSPI,LUNMRC,LUNDOC,LUNXM2, 
C     &       .FALSE.,FILNAM2,NLET2,DISP,
C     &       IFORM1,NX1,NY1,NZ1,NSTACK2,
C     &       FILNAM2,
C     &       .TRUE., ILIST2,NILMAX, 
C     &       NOT_USED,NGOT2,IMG2, IRTFLG) 
C       |
C     LOOP
C       REDLIN
C        |
C       WRTLIN
C        |
C       OPEN NEXT SPIDER INPUT & MRC OUTPUT FILES
C        |
C       NEXTFILES(NINDX1,NINDX2,  ILIST1,ILIST2,   
C               .FALSE., LUNXM1,LUNXM2,
C               NGOT1,NGOT2,    NSTACK1,NSTACK2,  
C               LUNSPI,LUNSPI,LUNMRC, FILNAM1,FILNAM2,
C               IMG1,IMG2, IRTFLG)
C        |
C     END LOOP
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE COPYTOMRC(LUNSPI,LUNMRC,
     &                     LUNDOC,LUNXM1,LUNXM2,IRTFLG)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      INTEGER                  :: IERR
      COMMON /IOERR/ IERR

      INTEGER                  :: BUF
      COMMON /IOBUF/ BUF(NBUFSIZ)

      INTEGER                  :: LUNSPI,LUNMRC,LUNDOC,LUNXM1,LUNXM2
      INTEGER                  :: IRTFLG

      INTEGER *2               :: I2MAX,I2VAL
      INTEGER *4               :: I4VAL
      INTEGER                  :: NC,MRCMODE,IGO,IEND,IDUM
      INTEGER                  :: NINDX1,NINDX2
      REAL                     :: FN,FNCON,FMINT,FMAXT,FAVT,FSIGT
      REAL                     :: SCALEX,SCALEY,SCALEZ

      INTEGER,ALLOCATABLE      :: ILIST1(:),ILIST2(:)

      CHARACTER (LEN=MAXNAM)   :: PROMPT 
      CHARACTER (LEN=MAXNAM)   :: FILNAM1,FILNAM2,MRCFILE
      CHARACTER (LEN=2*MAXNAM) :: COMMAN
      LOGICAL                  :: VERBOSE_SAVE,IS_MRC,IS_BARE
      CHARACTER (LEN=1)        :: NULL = CHAR(0)
      CHARACTER (LEN=1)        :: DISP
      CHARACTER (LEN=4)        :: CAXIS,CSTR
      INTEGER                  :: ICOMM,MYPID,MPIERR
      INTEGER                  :: IFORM1,NX1,NY1,NZ1,NSTACK1,NSTACK2
      INTEGER                  :: NGOT1,NGOT2,IMG1,IMG2,NILMAX,IVAL
      INTEGER                  :: MAXIM1,MAXIM2,NLET1,NLET2
      INTEGER                  :: LENT,NOT_USED,IRECIN

      INTEGER                  :: lnblnkn      ! FUNCTION
        
      CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

      VERBOSE_SAVE = VERBOSE           ! SAVE CURRENT VERBOSITY

      NILMAX  = NIMAX             ! FROM CMLIMIT
      ALLOCATE(ILIST1(NIMAX),
     &         ILIST2(NIMAX),
     &         STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
         CALL ERRT(46,'COPYTOMRC; ILIST...',2*NIMAX)
         RETURN
      ENDIF

C     OPEN FIRST INPUT FILE, DISP = 'E' DOES NOT STOP ON ERROR
      MAXIM1 = 0
      IMG1   = 0
      PROMPT =
     &     'SPIDER INPUT FILE NAME OR TEMPLATE (e.g. SPI_STK@*)~~9'
      CALL OPFILES(0,LUNSPI,LUNDOC,LUNXM1,  
     &               .TRUE.,FILNAM1,NLET1, 'O',
     &               IFORM1,NX1,NY1,NZ1,NSTACK1,
     &               PROMPT,
     &              .TRUE., ILIST1,NILMAX, 
     &               NOT_USED,NGOT1,IMG1, IRTFLG) 


#if defined(SP_DBUGIO)
      write(3,*)' ------- RETURNED TO: COPYTOMRC --------'
      write(3,*)' In copytomrc; nx1,ny1,nz1: ',nx1,ny1,nz1
      write(3,*)' '

#endif

      CALL LUNGETIS_MRC(LUNSPI,IS_MRC,IRTFLG)      
      CALL LUNGETISBARE(LUNSPI,IS_BARE,IRTFLG)

      IF (IS_MRC) THEN
         CALL ERRT(101,'OPERATION DOES NOT READ MRC FILES',IDUM)
         GOTO 999
      ENDIF

#if defined(SP_DBUGIO)
      !write(3,'(A,4i6)')' In nstack1; ngot1,img1:',nstack1,ngot1,img1
#endif


C     FIND HEADER VALUES FOR THE MRC FILE

      WRITE(NOUT,*) ' MRC DATA MODE 0:  8-BIT SIGNED INT, ',
     &              ' 1: 16-BIT SIGNED INT  '
      WRITE(NOUT,*)  '               2: 32-BIT REAL,       ',
     &              ' 6: 16-BIT UNSIGNED INT'

      MRCMODE = 2
      CALL RDPRI1S(MRCMODE,NOT_USED,'MRC MODE (0/1/2/6)',IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 999

C     LEGACY INPUT TRAP
      IF (MRCMODE ==  8) MRCMODE = 0   !  8 BIT UNSIGNED INT
      IF (MRCMODE == 16) MRCMODE = 6   ! 16 BIT UNSIGNED INT
      IF (MRCMODE == 32) MRCMODE = 2   ! 32 BIT FLOATING POINT

       WRITE(NOUT,*) ' MRC DATA ORIGIN     UL: UPPER LEFT, ',
     &              ' LL: LOWER LEFT  '
      IF (NZ1 > 1) THEN
         WRITE(NOUT,*) ' MRC HANDEDNESS  L: LEFT, ','        R:  RIGHT '
      ENDIF

      CAXIS = 'UL  '        ! DEFAULT

      IF (NZ1 <= 1) THEN
         CALL RDPRMC(CSTR,NC,.TRUE.,'MRC DATA ORIGIN (UL/LL)',
     &               NULL,IRTFLG)
         IF (IRTFLG .NE. 0)  RETURN

         IF (CSTR(1:NC) == 'LL')  CAXIS(1:2) = 'LL'
 
      ELSE
C        VOLUME OUTPUT
         CALL RDPRMC(CSTR,NC,.TRUE.,  
     &          'DATA ORIGIN (UL/LL) & HANDEDNESS (L/R)',NULL,IRTFLG)
         IF (IRTFLG .NE. 0)  RETURN

C        GET FIRST TOKEN (CHAR. STRING DELIMITED BY A ", ( ) ] -")
         CALL GETNEXTTOKEN_N(CSTR,NC,1,IGO,IEND)
         IF (IGO <= 0) THEN
            CALL ERRT(101,'INVALID INPUT',IDUM)
            GOTO 999
         ENDIF

         IF (CSTR(IGO:IGO+1) == 'LL') CAXIS(1:2) = 'LL'

         IF (CAXIS(1:2) .NE. 'UL' .AND. CAXIS(1:2) .NE. 'LL')THEN
            CALL ERRT(101,'INVALID ORIGIN',IDUM)
            RETURN
         ENDIF

C        GET SECOND TOKEN (CHAR. STRING DELIMITED BY A ", ( ) ] -")
         CALL GETNEXTTOKEN_N(CSTR,NC,IEND+1,IGO,IEND)
         IF (IGO <= 0) THEN
            CALL ERRT(101,'INVALID INPUT',IDUM)
            GOTO 999
         ENDIF
         CAXIS(4:4) = 'L'

         IF     (CSTR(IGO:IGO) == '0') THEN
            CAXIS(4:4) = 'R' 
         ENDIF 
      
         IF (CAXIS(4:4) .NE. 'L' .AND. CAXIS(4:4) .NE. 'R')THEN
            CALL ERRT(101,'INVALID HANDEDNESS',IDUM)
            RETURN
         ENDIF
      ENDIF

C     FIND HEADER VALUES FOR THE MRC FILE

      IF (MRCMODE == 2 .AND. MAXIM1 < 0) THEN
C        32 BIT FLOATING POINT  IMAGE OR VOLUME (NOT WHOLE STACK)

C        NEED FMIN & FMAX (SIG IN COMMON: MASTER)
         IF (IMAMI .NE. 1) CALL NORM3(LUNSPI,NX1,NY1,NZ1,FMAX,FMIN,AV)

         FMINT = FMIN
         FMAXT = FMAX
         FAVT  = AV
         FSIGT = SIG
C        UNFININISHED HERE !!!!!!!!!!!!!!!!!!!

      ELSEIF (MRCMODE == 0  .AND. MAXIM1 < 0 ) THEN
C        8 BIT SIGNED INTEGER IMAGE OR VOLUME (NOT WHOLE STACK)

C        NEED FMIN & FMAX (SIG IN COMMON: MASTER)
         IF (IMAMI .NE. 1) CALL NORM3(LUNSPI,NX1,NY1,NZ1,FMAX,FMIN,AV)
         IF (FMIN < -128) THEN
            I4VAL = FMIN
            CALL ERRT(102,'POSSIBLE 8 BIT INTEGER UNDERFLOW',I4VAL)
            GOTO 999
         ELSEIF (FMAX > 127) THEN 
            I4VAL = FMAX
            CALL ERRT(102,'POSSIBLE 8 BIT INTEGER OVERFLOW',I4VAL)
            GOTO 999
         ENDIF

      ELSEIF (MRCMODE == 1  .AND. MAXIM1 < 0 ) THEN
C        16 BIT SIGNED INTEGER IMAGE OR VOLUME (NOT WHOLE STACK)

C        NEED FMIN & FMAX (SIG IN COMMON: MASTER)
         IF (IMAMI .NE. 1) CALL NORM3(LUNSPI,NX1,NY1,NZ1,FMAX,FMIN,AV)

         IF (FMIN < -HUGE(I2VAL)) THEN
            I4VAL = FMIN
            CALL ERRT(102,'POSSIBLE 16 BIT INTEGER UNDERFLOW',I4VAL)
            GOTO 999
         ELSEIF (FMAX > HUGE(I2VAL)) THEN 
            I4VAL = FMAX
            CALL ERRT(102,'POSSIBLE 16 BIT INTEGER OVERFLOW',I4VAL)
            GOTO 999
         ENDIF

      ELSEIF (MRCMODE == 6  .AND. MAXIM1 < 0 ) THEN
C        16 BIT UNSIGNED INTEGER (NOT WHOLE STACK)

C        NEED FMIN & FMAX (SIG IN COMMON: MASTER)
         IF (IMAMI .NE. 1) CALL NORM3(LUNSPI,NX1,NY1,NZ1,FMAX,FMIN,AV)

         IF (FMIN < -HUGE(I2VAL)) THEN
            I4VAL = FMIN
            CALL ERRT(102,'POSSIBLE 16 BIT INTEGER UNDERFLOW',I4VAL)
            GOTO 999

         ELSEIF (FMAX > HUGE(I2VAL)) THEN 
            I4VAL = FMAX
            CALL ERRT(102,'POSSIBLE 16 BIT INTEGER OVERFLOW',I4VAL)
            GOTO 999
         ENDIF
      ENDIF


C     TRY TO GET SPIDER IMAGE SCALE VALUE (USUALLY NOT SET)
      CALL GETLAB(LUNSPI,21,1,SCALEX,IRTFLG)

      IF (SCALEX > 0) THEN
           WRITE(NOUT,*) ' PIXEL SIZE (A) FOR ALL AXES IN ',
     &                   ' SPIDER HEADER(21): ',SCALEX
           SCALEY  = SCALEX
           SCALEZ  = SCALEX
      ELSE
C          NOTHING IN SPIDER HEADER?
           SCALEX = 1.0
           SCALEY = 1.0
           SCALEZ = 1.0
      ENDIF

        
      IF (NZ1 > 1) THEN
         CALL RDPRM3S(SCALEX,SCALEY,SCALEZ,NOT_USED,
     &              'PIXEL SIZE (A) FOR X, Y, & Z AXES',IRTFLG)
      ELSE
         CALL RDPRM2S(SCALEX,SCALEY,NOT_USED,
     &              'PIXEL SIZE (A) FOR X &  Y AXES',IRTFLG)
      ENDIF
      IF (IRTFLG .NE. 0) RETURN

      WRITE(NOUT,*) ' ' 
      NSTACK2 =  1   ! UNUSED, NOT SENT
      DISP    = 'N'  ! NEW OUTPUT FILE        

C     OPEN FIRST OUTPUT FILE 
      IMG2   = IMG1         ! DEFAULT
      PROMPT =
     &    'MRC OUTPUT FILE NAME OR TEMPLATE (e.g. *@file_stk.mrc)~'

      CALL OPFILES(LUNSPI,LUNMRC,LUNDOC,LUNXM2, 
     &              .TRUE.,FILNAM2,NLET2,DISP,
     &              IFORM1,NX1,NY1,NZ1,NSTACK2,
     &              PROMPT,
     &              .FALSE., ILIST2,NILMAX, 
     &              NOT_USED,NGOT2,IMG2, IRTFLG) 

#if defined (SP_DBUGIO)
      write(3,*)' In copytomrc; nx1,ny1,nz1:    ', nx1,ny1,nz1
      write(3,*)' In copytomrc; nstack2,nilmax: ', nstack2,nilmax
      write(3,*)' In copytomrc; ngot2,img2:     ', ngot2,img2
      write(3,*)' '
#endif

      CALL LUNSETMODE_MRC(LUNMRC,MRCMODE,IRTFLG)

      CALL LUNSETHAND_MRC(LUNMRC,CAXIS,IRTFLG)

      CALL LUNSETPIXSIZES_MRC(LUNMRC,SCALEX,SCALEY,SCALEZ,IRTFLG)

      CALL LUNWRTHED_MRC(LUNMRC,IRTFLG)

      CALL LUNSETPOS_MRC(LUNMRC,IMG1,IRTFLG)

C     SETS LUNMRCNBYT(LUN) = NBYT    ! COMMENTED OUT??

C     DO NOT REPORT FILE INFO IF WHOLE STACK
      IF (NSTACK1 > 0 .AND. NSTACK2 >= 0) VERBOSE = .FALSE. 

      NINDX1 = 1
      NINDX2 = 1
      DO                ! LOOP OVER ALL IMAGES or IMAGE STACKS

C        COPY THE DESIRED NUMBER OF DATA RECORDS FROM EACH FILE
         DO IRECIN = 1,NY1*NZ1   ! NUMBER OF RECORDS COPIED

C           READ EACH ROW OF SPIDER INPUT FILE  (READS INPUT DATA)
            CALL REDLIN(LUNSPI,BUF,NX1,IRECIN)

C           WRITE EACH ROW OF SPIDER INPUT FILE 
            IF (MRCMODE == 2) THEN
C              32 BIT FLOATING POINT            (WRITES OUTPUT DATA)
               CALL WRTLIN_MRC(LUNMRC,BUF,NX1,IRECIN,MYPID,IERR)

            ELSEIF (MRCMODE == 0 ) THEN
C              8 BIT SIGNED INTEGER            (WRITES OUTPUT DATA)
               CALL WRTLIN_MRC(LUNMRC,BUF,NX1,IRECIN,MYPID,IERR)
 
            ELSEIF (MRCMODE == 1) THEN
C              16 BIT SIGNED INTEGER            (WRITES OUTPUT DATA)
               CALL WRTLIN_MRC(LUNMRC,BUF,NX1,IRECIN,MYPID,IERR)
 
            ELSEIF (MRCMODE == 6) THEN
C              16 BIT UNSIGNED INTEGER            (WRITES OUTPUT DATA)
               CALL WRTLIN_MRC(LUNMRC,BUF,NX1,IRECIN,MYPID,IERR)
 
            ENDIF

            IF (IERR .NE. 0) THEN
              LENT = lnblnkn(MRCFILE)
              WRITE(NOUT,99) IERR,IRECIN,NX1,LUNMRC,MRCFILE(:LENT)
99            FORMAT( '  *** ERROR(',I4,') WRITING RECORD: ',I12,
     &                ' LENGTH: ', I9,' UNIT: ',I3,' FILE: ',A)
              CALL ERRT(101,'ON MRC FILE OUTPUT',IERR)
              GOTO 999
            ENDIF

         ENDDO     ! END OF: DO IRECIN = 1,NY*NZ


C        OPEN NEXT SET OF I/O FILES, UPDATES NINDX1 & NINDX2

#if defined (SP_DBUGIO)
         write(3,*)' In copytomrc; calling nextfiles:',nindx1,nindx2
         write(3,*) ' ==================  next set ===================='
#endif
                    
         CALL NEXTFILES(NINDX1,NINDX2,  ILIST1,ILIST2, 
     &                  .FALSE., LUNXM1,LUNXM2,
     &                  NGOT1,NGOT2,    NSTACK1,NSTACK2,  
     &                  LUNSPI,LUNSPI,LUNMRC, FILNAM1,FILNAM2,
     &                  IMG1,IMG2, IRTFLG)

#if defined (SP_DBUGIO)
         write(3,*)' After nextfiles; ngot1,nstack1,img1: ',
     &                                ngot1,nstack1,img1
         write(3,*)' After nextfiles; ngot2,nstack2,img2: ',
     &                                ngot2,nstack2,img2
         write(3,*)' After nextfiles; irtflg:             ', irtflg
#endif

         IF (IRTFLG .NE. 0) EXIT      ! ERROR / END OF INPUT STACK
      ENDDO

      IRTFLG = 0
   
999   CLOSE(LUNSPI)
      CLOSE(LUNMRC)

      VERBOSE = VERBOSE_SAVE          ! RESTORE VERBOSITY 
      IF (ALLOCATED(ILIST1)) DEALLOCATE(ILIST1)
      IF (ALLOCATED(ILIST2)) DEALLOCATE(ILIST2)

      END

