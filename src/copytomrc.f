
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
C              REWRITE                             DEC 25 ArDean Leith
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
C COPYTOMRC(LUNSPI,LUNMRC,LUNDOC,LUNXM1,LUNXM2,IRTFLG)
C
C PURPOSE: COPY FROM SPIDER TO MRC FILE FORMAT
C
C NOTES: DATA IN MRC FILE
C        MODE   TYPES OF PIXELS 
C               0 : INTEGER*1 (UNSIGNED BYTES) 
C               1 : INTEGER*2 (SIGNED) 
C               2 : REALS
C               6 : INTEGER*2 (UNSIGNED)
C
C  CALL TREE:
C
C  IF ('CP....')   		
C      |
C  COPY1
C      |
C  IF ('CP TO MRC')    SPIDER FILE(S) TO MRC IMG(S) / VOL(S)
C      |
C  COPYTOMRC
C      |
C      |->OPFILES          OPEN FIRST SPIDER INPUT FILE
C           |-> FILERD
C      |         ` ->  ECHONAME
C      |-> FILELIST
C      |-> LUNGETIS_MRC      
C      |-> LUNGETISBARE
C      |-> FILERD

C      |-> OPENFIL_MRC     OPEN FIRST MRC OUTPUT FILE
C           |-> GET_FILNAM_INFO
C           |-> LUNNEWHED
C           |-> OPENFIL_N_MRC 
C               |-> OPSTREAMFIL     
C               |-> LUNSET_FILE_MRC
C               |-> LUNSET_MODE_MRC 
C               |-> LUNSET_SIZE_MRC 
C               |-> LUNSET_VIN_MRC   
C               |-> LUNSET_XXX_MRC 
C               |-> LUNZERO_STATS_MRC
C           |-> LUNSET_STK_260
C           |-> LUNSET_ISBARE_MRC
C           |-> WHICH_HAND_MRC --> LUNSET_HAND_MRC
C           |-> LUNSET_POS_MRC
C           |-> LUNWRTHED_MRC
C           |-> LUNGET_TYPE_MRC
C           |-> LUNSET_COMMON_MRC
C           |-> LUNSAYINFO_MRC
C      |-> LUNSET_MODE_MRC      NEEDED HERE, SINCE NOT SENT TO OPENFIL_N
C      |-> LUNSET_HAND_MRC      NEEDED HERE, SINCE NOT SENT TO OPENFIL_N
C      |-> LUNSET_PIXSIZES_MRC  NEEDED HERE or somewhere above
C      |-> LUNWRT_HED_MRC       NEEDED HERE  SINCE NOT SENT TO OPENFIL_N
C      |-> LUNSET_POS_MRC       ?? NEEDED HERE?

C THESE CALLS ONLY USED IN:

C  copytomrc.f:      CALL LUNSET_MODE_MRC originally set there
C  openfil_n_mrc.f:  CALL LUNSET_MODE_MRC may change, needed here also

C  copytomrc.f:      CALL LUNSET_HAND_MRC
C  copyfrommrc.f:    CALL LUNSET_HAND_MRC
C  lunsetmrchdr.f:   CALL LUNSET_HAND_MRC  2x

C  copytomrc.f:      CALL LUNSET_PIXSIZES_MRC   NOT IN SPIDER
C  lunsethdr.f:      CALL LUNSET_PIXSIZ_MRC   IN SPIDER!!  DIFFERENT??

C  copytomrc.f:      CALL LUNSET_POS_MRC
C  copyfrommrc.f:    CALL LUNSET_POS_MRC 
C  openfil_mrc.f:    CALL LUNSET_POS_MRC
C  opfiles_mrc.f:    CALL LUNSET_POS_MRC  2x

C     LOOP
C        REDLIN
C        |
C        WRTLIN
C        |
C        NEXTFILES    OPEN NEXT SPI INPUT & MRC OUTPUT FILE
C        |
C     END LOOP
C
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

      INTEGER *2               :: I2VAL
      INTEGER *4               :: I4VAL
      INTEGER                  :: NC,MRCMODE,IGO,IEND,IDUM
      INTEGER                  :: NINDX1,NINDX2
      REAL                     :: FN,FNCON,FMINT,FMAXT,FAVT,FSIGT
      REAL                     :: SCALEX,SCALEY,SCALEZ

      INTEGER,ALLOCATABLE      :: ILIST1(:),ILIST2(:)

      CHARACTER (LEN=MAXNAM)   :: PROMPT 
      CHARACTER (LEN=MAXNAM)   :: FILNAM1,FILNAM2,MRCFILE
      CHARACTER (LEN=2*MAXNAM) :: COMMAN
      LOGICAL                  :: VERBOSE_SAVE,IS_MRC,ASKNAM
      LOGICAL                  :: IS_BARE_2,IS_BARE_1
      CHARACTER (LEN=1)        :: NULL = CHAR(0)
      CHARACTER (LEN=1)        :: DSP

      CHARACTER (LEN=4)        :: CAXIS,CSTR
      INTEGER                  :: ICOMM,MYPID,MPIERR
      INTEGER                  :: IFORM1,NX1,NY1,NZ1,NSTK1,NSTK2
      INTEGER                  :: NLIST1,NLIST2, INUM1,INUM2, NILMAX
      INTEGER                  :: MAXIM1, NLET1,NLET2, NREC
      INTEGER                  :: LENT,NOT_USED,IRECIN, ITYPE
      INTEGER                  :: NSTK,ISTK1,ISTK2, irep

      INTEGER                  :: lnblnkn      ! FUNCTION
        
      CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

      VERBOSE_SAVE = VERBOSE           ! SAVE CURRENT VERBOSITY

      NILMAX  = NIMAX                  ! FROM CMLIMIT
      ALLOCATE(ILIST1(NILMAX),
     &         ILIST2(NILMAX),
     &         STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
         CALL ERRT(46,'COPYTOMRC; ILIST...',2*NIMAX)
         RETURN
      ENDIF
 
C     OPEN FIRST SPIDER FILE, DSP = 'E' DOES NOT STOP ON ERROR
      MAXIM1 = 0
      INUM1  = 0
      PROMPT =
     &     'SPIDER INPUT FILE NAME OR TEMPLATE (e.g. SPI_STK@*)~'
      CALL OPFILES(0,LUNSPI,LUNDOC,LUNXM1,  
     &             .TRUE.,FILNAM1,NLET1, 'O',
     &              IFORM1,NX1,NY1,NZ1,NSTK1,
     &              PROMPT,
     &             .FALSE., ILIST1,NILMAX, 
     &              NOT_USED,NLIST1,INUM1, IRTFLG) 


#if defined(SP_DBUGIO)
      write(3,*)' ------- RETURNED TO: COPYTOMRC --------'
      write(3,*)' In copytomrc; nx1,ny1,nz1: ',nx1,ny1,nz1
      write(3,*)' In copytomrc; nstk1:       ',nstk1
      write(3,*)' In copytomrc; inum1:       ',inum1
      write(3,*)' '
#endif

      CALL LUNGETIS_MRC(LUNSPI,IS_MRC,IRTFLG)      
      CALL LUNGETISBARE(LUNSPI,IS_BARE_1,IRTFLG)

      IF (IS_MRC) THEN
         CALL ERRT(101,'OPERATION DOES NOT READ MRC FILES',IDUM)
         GOTO 999
      ENDIF

#if defined(SP_DBUGIO)
      !write(3,'(A,4i6)')' In nstk1; NLIST1,INUM1:',nstk1,NLIST1,INUM1
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

         if (CSTR(IGO:IGO) == '0') THEN
            CAXIS(4:4) = 'R' 
         ENDIF 
      
         IF (CAXIS(4:4) .NE. 'L' .AND. CAXIS(4:4) .NE. 'R')THEN
            CALL ERRT(101,'INVALID HANDEDNESS',IDUM)
            RETURN
         ENDIF
      ENDIF

C     FIND HEADER VALUES FOR THE NEW MRC FILE

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
      IF (LUNMRC <= 0 .OR. LUNMRC > 100) THEN
         CALL ERRT(102,'PGM. ERROR: LUN MUST BE 1...100',LUNMRC)
         RETURN
      ENDIF

C     OPEN FIRST MRC OUTPUT FILE 
      MAXIM1 = 0
      INUM2  = 0
      DSP    = 'N'      ! BUT MAY ALREADY EXIST IF A STACKED FILE
      ITYPE  =  0       ! UNUSED BY OPENFIL_MRC  ON INPUT
      NSTK2  =  1       ! UNUSED BY OPENFIL_MRC  ON INPUT

      PROMPT =
     &    'MRC OUTPUT FILE NAME OR TEMPLATE (e.g. *@file_stk.mrc)~~9'

      CALL OPFILES(0,LUNMRC,LUNDOC,LUNXM2,  
     &              .TRUE.,FILNAM2,NLET2, 'N',
     &               ITYPE, NX1,NY1,NZ1, NSTK2,
     &               PROMPT,
     &              .TRUE., ILIST2,NILMAX, 
     &               NOT_USED,NLIST2,INUM2, IRTFLG) 

      CALL LUNGETIS_MRC(LUNMRC,IS_MRC,IRTFLG)      
      CALL LUNGETISBARE(LUNMRC,IS_BARE_2,IRTFLG)

      IF (.NOT. IS_MRC) THEN
         CALL ERRT(101,'OPERATION ONLY CREATES MRC FILES',IDUM)
         GOTO 999
      ENDIF
    
#if defined (SP_DBUGIO)
      write(3,*)' In copytomrc; nx1,ny1,nz1:  ', nx1,ny1,nz1
      write(3,*)' in copytomrc; nstk2,inum2:  ', nstk2,inum2
      write(3,*)' in copytomrc; is_bare_2:  ', is_bare_2
      write(3,*)' '
#endif

      CALL LUNSET_PIXSIZES_MRC(LUNMRC,SCALEX,SCALEY,SCALEZ,IRTFLG)
      CALL LUNSET_HAND_MRC(LUNMRC,CAXIS,IRTFLG)
      CALL LUNSET_MODE_MRC(LUNMRC,MRCMODE,IRTFLG)
      CALL LUNWRTHED_MRC(LUNMRC,IRTFLG)
      CALL LUNSET_POS_MRC(LUNMRC,INUM2,IRTFLG)

C     DO NOT REPORT FILE INFO IF WHOLE STACK
      IF (NSTK1 > 0 .AND. NSTK2 >= 0) VERBOSE = .FALSE. 

      NINDX1 = 1
      NINDX2 = 1

      NREC   = NY1*NZ1  ! NUMBER OF RECORDS IN IMG/VOL

      irep = 0

      DO                ! LOOP OVER ALL IMAGES or IMAGE STACKS
         irep = irep + 1

C        COPY THE DESIRED NUMBER OF DATA RECORDS FROM EACH FILE
         DO IRECIN = 1,NREC   ! NUMBER OF RECORDS COPIED

C           READ EACH ROW OF SPIDER INPUT FILE  (READS INPUT DATA)
            CALL REDLIN(LUNSPI,BUF,NX1,IRECIN)

C           WRITE EACH ROW FROM SPIDER INPUT FILE 
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
         write(3,*)' In copytomrc; irep:   ', irep
         write(3,*)' In copytomrc; records:', nrec
         write(3,*)' In copytomrc; calling nextfiles:',nindx1,nindx2
         write(3,*) ' ==================  Next set ===================='
         
         !if (irep > 1) exit
         exit

#endif
 
                   
         CALL NEXTFILES(NINDX1,NINDX2,  ILIST1,ILIST2, 
     &                  .FALSE., LUNXM1,LUNXM2,
     &                  NLIST1,NLIST2,    NSTK1,NSTK2,  
     &                  LUNSPI,LUNSPI,LUNMRC, FILNAM1,FILNAM2,
     &                  INUM1,INUM2, IRTFLG)

#if defined (SP_DBUGIO)
         write(3,*)' After nextfiles; NLIST1,nstk1,inum1: ',
     &                                NLIST1,nstk1,inum1
         write(3,*)' After nextfiles; NLIST2,nstk2,inum2: ',
     &                                NLIST2,nstk2,inum2
         write(3,*)' After nextfiles; irtflg:           ', irtflg
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

