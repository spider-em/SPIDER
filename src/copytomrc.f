
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
C
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
C COPYTOMRC    (LUNSPI,LUNMRC,NX,NY,NZ)
C COPYTOMRC_STK(LUNSPI,LUNMRC,NX,NY,NZ) 
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
      INTEGER                  :: NDUM,NGOT1,NGOT2,IMG1,IMG2,NILMAX,IVAL
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
      PROMPT =
     &     'SPIDER INPUT FILE NAME OR TEMPLATE (e.g. SPI_STK@**)~~9'
      CALL OPFILES(0,LUNSPI,LUNDOC,LUNXM1,  
     &               .TRUE.,FILNAM1,NLET1, 'E',
     &               IFORM1,NX1,NY1,NZ1,NSTACK1,
     &               PROMPT,
     &              .TRUE., ILIST1,NILMAX, 
     &               NDUM,NGOT1,IMG1, IRTFLG) 


      CALL LUNGETIS_MRC(LUNSPI,IS_MRC,IRTFLG)      
      CALL LUNGETISBARE(LUNSPI,IS_BARE,IRTFLG)

      IF (IS_MRC) THEN
         CALL ERRT(101,'OPERATION DOES NOT READ MRC FILES',NDUM)
         GOTO 999
      ENDIF

      !write(3,'(A,4i6)')' In nstack1,ngot1,img1:',nstack1,ngot1,img1

      CALL FILERD(FILNAM2,NLET2,NULL,
     &       'MRC OUTPUT FILE NAME OR TEMPLATE (e.g. **@MRC_STK.mrc)~',
     &       IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 999

      IF (INDEX(FILNAM2,'.mrc') <= 0 .AND.
     &    INDEX(FILNAM2,'.MRC') <= 0) THEN
         CALL ERRT(101,
     &     'OUTPUT FILE NAME MUST HAVE MRC/mrc EXTENSION',NDUM)
         GOTO 999
      ENDIF


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

      WRITE(NOUT,*) ' ' 
      WRITE(NOUT,*) ' MRC DATA ORIGIN     UL: UPPER LEFT, ',
     &              ' LL: LOWER LEFT  '
      IF (NZ1 > 1) THEN
         WRITE(NOUT,*) ' MRC HANDEDNESS  L: LEFT, ',
     &                 '        R:  RIGHT '
      ENDIF

      CAXIS(1:2) = 'UL'        ! DEFAULT

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
C     UNFININISHED HERE !!!!!!!!!!!!!!!!!!!

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

      NSTACK2 =  1   ! UNUSED
      DISP    = 'U'  ! NEW OUTPUT FILE        

C     OPEN FIRST OUTPUT FILE 
      IMG2 = IMG1
      CALL OPFILES(LUNSPI,LUNMRC,LUNDOC,LUNXM2, 
     &              .FALSE.,FILNAM2,NLET2,DISP,
     &              IFORM1,NX1,NY1,NZ1,NSTACK2,
     &              FILNAM2,
     &              .TRUE., ILIST2,NILMAX, 
     &              NDUM,NGOT2,IMG2, IRTFLG) 

      !write(3,'(A,4i6)')' In copytomrc nstack2,ngot2:',nstack2,ngot2
      !write(3,'(A,i6,A)')' In copytomrc img2:',img2
      !write(3,*)         ' In copytomrc caxis:',caxis

      CALL LUNSETMODE_MRC(LUNMRC,MRCMODE,IRTFLG)
      CALL LUNSETHAND_MRC(LUNMRC,CAXIS,IRTFLG)
      CALL LUNSETPIXSIZES_MRC(LUNMRC,SCALEX,SCALEY,SCALEZ,IRTFLG)
      CALL LUNWRTHED_MRC(LUNMRC,IRTFLG)

C     SETS LUNMRCNBYT(LUN) = NBYT
      CALL LUNSETPOS_MRC(LUNMRC,IMG1,IRTFLG)

C     DO NOT REPORT FILE INFO IF WHOLE STACK
      IF (NSTACK1 > 0 .AND. NSTACK2 >= 0) VERBOSE = .FALSE. 

      NINDX1 = 1
      NINDX2 = 1
      DO                ! LOOP OVER ALL IMAGES/STACKS

C        COPY THE DESIRED NUMBER OF DATA RECORDS FROM EACH FILE
         DO IRECIN = 1,NY1*NZ1

C           READ EACH ROW OF SPIDER INPUT FILE 
            CALL REDLIN(LUNSPI,BUF,NX1,IRECIN)

C           WRITE EACH ROW OF SPIDER INPUT FILE 
            IF (MRCMODE == 2) THEN
C              32 BIT FLOATING POINT  
               CALL WRTLIN_MRC(LUNMRC,BUF,NX1,IRECIN,MYPID,IERR)

            ELSEIF (MRCMODE == 0 ) THEN
C              8 BIT SIGNED INTEGER
               CALL WRTLIN_MRC(LUNMRC,BUF,NX1,IRECIN,MYPID,IERR)
 
            ELSEIF (MRCMODE == 1) THEN
C              16 BIT SIGNED INTEGER
               CALL WRTLIN_MRC(LUNMRC,BUF,NX1,IRECIN,MYPID,IERR)
 
            ELSEIF (MRCMODE == 6) THEN
C              16 BIT UNSIGNED INTEGER
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
         !write(3,*)' In copytomrc, calling nextfiles:',nindx1,nindx2 

         CALL NEXTFILES(NINDX1,NINDX2,  ILIST1,ILIST2, 
     &                  .FALSE., LUNXM1,LUNXM2,
     &                  NGOT1,NGOT2,    NSTACK1,NSTACK2,  
     &                  LUNSPI,LUNSPI,LUNMRC, FILNAM1,FILNAM2,
     &                  IMG1,IMG2, IRTFLG)
         IF (IRTFLG .NE. 0) EXIT      ! ERROR / END OF INPUT STACK
      ENDDO

      IRTFLG = 0
   
999   CLOSE(LUNSPI)
      CLOSE(LUNMRC)

      VERBOSE = VERBOSE_SAVE          ! RESTORE VERBOSITY 
      IF (ALLOCATED(ILIST1)) DEALLOCATE(ILIST1)
      IF (ALLOCATED(ILIST2)) DEALLOCATE(ILIST2)

      END





C **********************************************************************
C
C COPYTOMRC_STK  'CP TO MRC' STACKS               JUN 16 ArDean Leith
C
C **********************************************************************
C 
C COPYTOMRC_STK(LUNSPI,LUNMRC,SPIPAT,ILIST,NIMG, MAXIM,NX,NY,NZ)
C                                                                      
C PURPOSE: CONVERTS SPIDER IMAGE STACK TO MRC FORMAT, CRUDELY WRITTEN!!
C          CREATES NATIVE X386 BYTE ORDER MRC FILES 
C
C NOTES: DATA MODE IN MRC FILE : 2  (32 BIT FLOATING POINT)
C        ONLY CREATES 32 BIT REAL STACKS AS NORMALIZATION DIFFICULT FOR
C        LARGE MRC STACKS.
C
C       CREATE A NEW HEADER FOR THE MRC FILE 2014
C       DMAX < DMIN                         MAX & MIN UNDETERMINED
C       DMEAN < (SMALLER OF DMIN and DMAX)  DMEAN     UNDETERMINED
C       RMS < 0.0                           RMS       UNDETERMINED
C
C       TYPE           ISPG       MZ
C       IMAGE            0         1
C       IMAGE STACK      0      >= 1
C       VOLUME           1      = NZ
C       VOLUME STACK    401     = NZ / # OF SINGLE VOLUMES
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE COPYTOMRC_STK(LUNSPI,LUNMRC,SPIPAT,
     &                           ILIST,NIMG,MAXIM,IMGNUMOUT,
     &                           NX,NY,NZ)
 
C       COPY STACKS FROM SPIDER TO MRC FILE FORMAT 

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
 
        INTEGER                 :: IERR
        COMMON /IOERR/  IERR       ! LEGACY REDLIN ERROR FLAG

        REAL                    :: BUFIN
        COMMON /IOBUF/  BUFIN(NBUFSIZ)

        INTEGER LUNARA,LUNSTK,LUNARB,LUNFLIP
        COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

        INTEGER                 :: LUNSPI,LUNMRC, NIMG,MAXIM,NX,NY,NZ
        CHARACTER(LEN=*)        :: SPIPAT
        INTEGER                 :: ILIST(*)
        INTEGER                 :: IMGNUMOUT

        REAL                    :: FIXLENBUF(256)
        REAL                    :: FLIPBUF(NX)

        CHARACTER(LEN=MAXNAM)   :: MRCFILE
        CHARACTER(LEN=4)        :: CSAXIS
        LOGICAL                 :: ISSWABT

        LOGICAL                 :: isswab,BOTLEFT
        INTEGER                 :: lnblnkn

        INTEGER *8              :: IPOSMRC
        INTEGER                 :: ISTARTX,ISTARTY,ISTARTZ
        REAL                    :: ORIGX,ORIGY,ORIGZ
        INTEGER                 :: IRTFLG,MODE,NX4
        INTEGER                 :: NSYMBT,MACHST,NE,IOFFSET 
        INTEGER                 :: NCHAR,LENT,ILEN,IBOTLEFT
        INTEGER                 :: NOT_USED
        INTEGER                 :: IRECIN,ILOCIN,IRECINT
        REAL                    :: RMS,FMINT,FMAXT,FAVT,FSIGT
        REAL                    :: UNUSED,SCALE,SCALEX,SCALEY,SCALEZ
        INTEGER                 :: NLET,IGO,MZ,IVERSION
        INTEGER                 :: NINDX,NGOT,IBITS,NSYMBYT,ISPG

        LOGICAL, PARAMETER      :: FOUROK = .FALSE.
        INTEGER, PARAMETER      :: LUNXM  = 0


C       OPEN NEW MRC STACK FILE FOR STREAM ACCESS, LITTLE ENDED
        CALL OPSTREAMFILE(.TRUE.,MRCFILE,DATEXC,LUNMRC,
     &                    'UNFORMATTED','NL',
     &                    'MRC OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IVERSION = 20140
        IBOTLEFT = 1               
        CALL RDPRI2S(IVERSION,IBOTLEFT,NOT_USED,
     &    'MRC VERSION (0 or 20140), ORIGIN AT BOTTOM LEFT (YES==1)',
     &    IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        BOTLEFT = (IBOTLEFT == 1)   ! USUAL MRC BOTTOM!!
        IF (IVERSION .NE. 20140 ) IVERSION = 0
  
C       MODE IS 32 BIT FLOATING POINT MRC IMAGE
        MODE    = 2                ! FLOATING POINT IMAGE

        ISTARTX = 0
        ISTARTY = 0
        ISTARTZ = 0

C       ORIGIN ON X,Y & Z AXIS  (MAKE THIS ALWAYS ZERO) 
        ORIGX   = 0.0
        ORIGY   = 0.0
        ORIGZ   = 0.0

C       TRY TO GET SPIDER IMAGE SCALE VALUE (MAY NOT BE SET?)
        CALL GETLAB(LUNSPI,21,1,SCALEX,IRTFLG)
        SCALEY  = SCALEX
        SCALEZ  = SCALEX

        IF (SCALEX > 0) THEN
           WRITE(NOUT,*) ' PIXEL SIZE (A) FOR ALL AXES (HEADER ',
     &                   ' LOCATION=21): ',SCALEX
        ELSEIF (NZ == 1) THEN

           SCALEX = 1.0
           SCALEY = 1.0
           SCALEZ = 1.0
           CALL RDPRM2S(SCALEX,SCALEY,NOT_USED,
     &              'PIXEL SIZE (A) FOR  X & Y AXES',IRTFLG)
        ELSE
           SCALEX = 1.0
           SCALEY = 1.0
           SCALEZ = 1.0
           CALL RDPRM3S(SCALEX,SCALEY,SCALEZ,NOT_USED,
     &              'PIXEL SIZE (A) FOR X, Y, & Z AXES',IRTFLG)
        ENDIF

C       SET UNDETERMINED IMAGE STATISTICS (FASTER)
        IF (IVERSION >= 20140) THEN
           FMINT =  0    
           FMAXT = -1
           FAVT  = -2
           FSIGT = -1
        ELSE
           FMINT =  0    ! UNDETERMINED SHUD BE: FMINT =  0 FOR MRC 2014
           FMAXT =  0    ! UNDETERMINED SHUD BE: FMAXT = -1
           FAVT  =  0    ! UNDETERMINED SHUD BE: FAVT  = -2
           FSIGT = -1    ! UNDETERMINED SHUD BE: FSIGT = -1
        ENDIF

 
        NSYMBYT = 0            ! NO ADDED HEADER LOCATIONS
        ISSWABT = .FALSE.

C       CREATE MRC HEADER. 
        CALL SETHEDCCP4_NEW(FIXLENBUF, NX,NY,NZ, NIMG, MODE,IVERSION,
     &                  ISTARTX,ISTARTY,ISTARTZ, 
     &                  ORIGX,ORIGY,ORIGZ, 
     &                  SCALEX,SCALEY,SCALEZ, 
     &                  FMINT,FMAXT,FAVT,FSIGT,
     &                  ISSWABT,NSYMBYT, IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       WRITE HEADER OF 1024 BYTES (256 FLOATS) TO MRC FILE
        IPOSMRC = 1  
        WRITE(LUNMRC, POS=IPOSMRC,IOSTAT=IERR) FIXLENBUF

        IF (IERR .NE. 0) THEN
           LENT = lnblnkn(MRCFILE)
           ILEN = 256
           WRITE(NOUT,99) IERR,IPOSMRC,ILEN,LUNMRC,MRCFILE(:LENT)
99         FORMAT( '  *** ERROR(',I4,') WRITING POSITION: ',I12,
     &              ' LENGTH: ', I9,' UNIT: ',I3,' FILE: ',A)
           GOTO 9999
        ENDIF

        IOFFSET = 1024
        NINDX   = IMGNUMOUT
        NGOT    = NIMG
        NX4     = NX * 4
          
C       LOOP OVER ALL IMAGES IN SPIDER STACK
        DO             ! LOOP OVER ALL IMAGES IN SPIDER STACK

C          SET STARTING POSITION FOR MRC DATA
           IPOSMRC = (IMGNUMOUT-1)  ! KLUDGE FOR INTEGER *8 PRESERVE
           IPOSMRC = IPOSMRC * NX * NY 
           IPOSMRC = IPOSMRC * 4
           IPOSMRC = IOFFSET + IPOSMRC + 1

c          write(6,'(2x,A,3i6,i12)') '  nindx,ilist,imgnumout,iposmrc:',
c     &                          nindx,ilist(nindx),imgnumout,iposmrc

           DO IRECIN = 1,NY*NZ

C             READ EACH ROW OF SPIDER INPUT FILE 


              IF (BOTLEFT) THEN
                 IRECINT = (NY * NZ) - IRECIN + 1    ! BOTTEM LEFT 
              ELSE
                 IRECINT = IRECIN
              ENDIF  
              CALL REDLIN(LUNSPI,BUFIN,NX,IRECINT)

C             PUT ROW OUT TO MRC FILE
              WRITE(LUNMRC, POS=IPOSMRC,IOSTAT=IERR) BUFIN(1:NX)

              IF (IERR .NE. 0) THEN
                LENT = lnblnkn(MRCFILE)
                WRITE(NOUT,99) IERR,IPOSMRC,NX,LUNMRC,MRCFILE(:LENT)
                GOTO 9999
              ENDIF
              IPOSMRC = IPOSMRC + NX4

           ENDDO     ! END OF: DO IRECIN = 1,NY*NZ

C          OPEN NEXT SPIDER INPUT FILE, UPDATES NINDX 
           CALL NEXTFILE(NINDX, ILIST, 
     &                  FOUROK, LUNXM,
     &                  NGOT,   MAXIM,   
     &                  LUNSPI, 0,
     &                  SPIPAT, 'O',
     &                  IMGNUMOUT, IRTFLG) 
         
           IF (NINDX > MAXIM) EXIT    ! FINISHED 
           IF (IRTFLG .NE. 0) EXIT    ! ERROR

        ENDDO
                      
9999    CLOSE(LUNSPI)
        CLOSE(LUNMRC)

        END

