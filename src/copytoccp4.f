
C ++********************************************************************
C
C COPYTOCCP4   MODIFIED FROM COPYMRC               FEB 02 ArDean Leith         
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
C              COPYTOCCP4_STK ADDED, STACK BUGGY   JUN 16 ArDean Leith
C              BOTLEFT OPTION REINTRODUCED         SEP 16 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2016  Health Research Inc.,                         *
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
C COPYTOCCP4(LUNSPI,LUNMRC,NX,NY,NZ)
C COPYTOCCP4_STK(LUNSPI,LUNMRC,NX,NY,NZ) 
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

        SUBROUTINE COPYTOCCP4(LUNSPI,LUNMRC,NX,NY,NZ)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
 
        COMMON /IOERR/  IERR       ! LEGACY REDLIN ERROR FLAG

        REAL                    :: BUFIN
        COMMON /IOBUF/  BUFIN(NBUFSIZ)

        INTEGER                 :: LUNSPI,LUNMRC, NX,NY,NZ

        REAL,       ALLOCATABLE :: STREAMBUF(:)
        INTEGER *1, ALLOCATABLE :: I1STREAMBUF(:)
        INTEGER *2, ALLOCATABLE :: I2STREAMBUF(:)
        INTEGER,    ALLOCATABLE :: ILIST(:)
        REAL                    :: BUF(NBUFSIZ),FIXLENBUF(256)
        INTEGER *1              :: I1BUF(1024)
        COMMON                     BUF,FIXLENBUF,I1BUF

        CHARACTER(LEN=MAXNAM)   :: MRCFILE,FILPAT,FILOUT
        CHARACTER(LEN=8)        :: ANS
        CHARACTER(LEN=80)       :: PROMPT
        LOGICAL                 :: FLIP,ISSWABT,BOTLEFT
        LOGICAL                 :: isswab
        INTEGER                 :: IVAL
        INTEGER *1              :: I1VAL
        INTEGER *2              :: I2VAL

        INTEGER *2              :: I2V
        INTEGER *1              :: I1V(2),I1TMP
        EQUIVALENCE                (I2V,I1V)

        REAL    *4              :: R4VALIN,R4VALOUT
        INTEGER *1              :: I1VALIN(4),I1VALOUT(4)
        EQUIVALENCE                (R4VALIN,I1VALIN),(R4VALOUT,I1VALOUT)

        INTEGER                 :: IERR,LENOPENB,LENOPENF,IRTFLG,MODE
        INTEGER                 :: NSYMBT,MACHST,NE,MAXIM,IOFFSET 
        INTEGER                 :: LENOPEN,NCHAR,IRECMRC
        INTEGER                 :: IBOTLEFT,NOT_USED,ILOCOUT
        INTEGER                 :: IRECIN,ILOCIN,IRECINT,NSYMBYT
        REAL                    :: RMS,FMINT,FMAXT,FAVT,FSIGT,FN,FNCON
        REAL                    :: UNUSED,SCALE,SCALEX,SCALEY,SCALEZ
        INTEGER                 :: IX,IY,IZ,NLET
        INTEGER                 :: I,NSTACKT,ITYPE,NUNUSED,NSTACKOUT
        INTEGER                 :: IMGNUMOUT,NSTACK,IGO
        INTEGER                 :: NIMG,MZ

        LOGICAL                 :: FOUROK

        IERR = 0

C       OPEN NEW MRC FILE FOR DIRECT ACCESS, RECORD LENGTH 1024 BYTES
        LENOPENB = 1024
        LENOPENF = LENOPENB / 4
        CALL OPAUXFILE(.TRUE.,MRCFILE,DATEXC,LUNMRC,LENOPENB,'U',
     &                 'MRC OUTPUT',.TRUE.,IRTFLG)

        IVAL     = 8
        IBOTLEFT = 1
        CALL RDPRI2S(IVAL,IBOTLEFT,NOT_USED,
     &   'MRC DATA LENGTH (8/32 BITS), ORIGIN AT BOTTOM LEFT (YES==1)',
     &   IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        MODE = 2
        IF (IVAL == 8) MODE = 0

        BOTLEFT = (IBOTLEFT == 1)   ! USUAL MRC BOTTOM!!

C       CREATE A NEW HEADER FOR THE CCP4 (MRC) FILE
        FMINT = FMIN
        FMAXT = FMAX
        FAVT  = AV
        FSIGT = SIG

        IF (MODE == 0) THEN
           FN    = (255.0 - 0.0) / (FMAXT - FMINT)
           FNCON = 0.0 - FN * FMINT

           FMINT = 0.0
           FMAXT = 255.0
           I2VAL = FMINT * FN + FNCON
           FAVT  = I2VAL
C          FSIGT IS NOT RIGHT!!!!
           FSIGT = -1.0
        ENDIF

C       TRY TO GET SCALE VALUE (MAY NOT BE USED)
        CALL GETLAB(LUNSPI,NX,UNUSED,21,1,SCALE,IRTFLG)
        SCALEX  = SCALE
        SCALEY  = SCALE
        SCALEZ  = SCALE

C       FIND IF CURRENTLY SWAPPING BYTES DURING FILE OUTPUT
C       THIS MAY BE DONE BY COMPILER, SO HAVE TO ACTUALLY TEST OUTPUT
        ISSWABT = ISSWAB(99)

C       FLIP BYTES DURING MRC FILE OUTPUT
        IF (ISSWABT) THEN
           CALL LUNSETFLIP(LUNMRC,1,IRTFLG)
           ISSWABT = .FALSE.
        ENDIF

C       CREATE HEADER. (NOTE: FMIN, FMAX, AV ARE SAME AS SPIDER IMAGE)
        NIMG    = 1
        MZ      = NZ
        NSYMBYT = 0
        CALL SETHEDCCP4(FIXLENBUF, NX, NY, NZ,
     &            FMINT,FMAXT,FAVT,FSIGT,
     &            SCALEX,SCALEY,SCALEZ,MODE,
     &            ISSWABT,NSYMBYT, NIMG,MZ,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       WRITE HEADER OF 1024 BYTES (256 FLOATS) TO MRC FILE
        CALL WRTLIN(LUNMRC,FIXLENBUF,LENOPENF,1) 

C       SET STARTING RECORD FOR MRC DATA
        IRECMRC = 1   !SKIPS ONE HEADER RECORD
        ILOCOUT = 0

        IF (MODE == 2) THEN
C          FLOATING POINT OUTPUT

           DO  IRECIN = 1,NY * NZ

C             READ EACH ROW OF SPIDER INPUT FILE 
              IF (BOTLEFT) THEN
                 IRECINT = (NY * NZ) - IRECIN + 1  
              ELSE
                 IRECINT = IRECIN  
              ENDIF
              CALL REDLIN(LUNSPI,BUFIN,NX,IRECINT)

C             PUT ROW OUT TO MRC FILE
              DO ILOCIN=1,NX
                ILOCOUT            = ILOCOUT + 1
                FIXLENBUF(ILOCOUT) = BUFIN(ILOCIN)

                IF (ILOCOUT >= LENOPENF) THEN
C                  PUT OUT COMPLETED RECORD

                   IRECMRC = IRECMRC + 1

                   CALL WRTLIN(LUNMRC,FIXLENBUF,LENOPENF,IRECMRC)
                   !write(6,*) '  Rec: ',irecint,'-->',irecmrc

                   ILOCOUT = 0
                ENDIF
              ENDDO
           ENDDO

           IF (ILOCOUT > 0) THEN
C             PUT OUT ANY REMAINING MRC RECORD
              IRECMRC = IRECMRC + 1
              CALL WRTLIN(LUNMRC,FIXLENBUF,ILOCOUT,IRECMRC)
              !write(6,*) '  Rec: ',irecint,'-->',irecmrc
           ENDIF
        
        ELSEIF (MODE == 0) THEN
C          COPY FROM SPIDER TO MRC 8 BIT FILE FORMAT 

           DO IRECIN = 1,NY * NZ
C             READ EACH ROW OF SPIDER INPUT FILE 
              IF (BOTLEFT) THEN
                 IRECINT = (NY * NZ) - IRECIN + 1  
                 CALL REDLIN(LUNSPI,BUFIN,NX,IRECINT)
              ELSE
                 CALL REDLIN(LUNSPI,BUFIN,NX,IRECIN)
              ENDIF

C             PUT ROW OUT TO MRC FILE
              DO ILOCIN=1,NX
                ILOCOUT        = ILOCOUT + 1
                I2VAL          = BUFIN(ILOCIN) * FN + FNCON
                I1BUF(ILOCOUT) = I2VAL

                IF (ILOCOUT >= LENOPENB) THEN
C                  PUT OUT COMPLETED RECORD
                   IRECMRC = IRECMRC + 1

                   CALL WRTLIN8(LUNMRC,I1BUF,LENOPENB,IRECMRC)
                   IF (IERR .NE. 0) THEN
                      CALL ERRT(102,'WRITING RECORD',IRECMRC)
                      GOTO 9999
                   ENDIF
                   ILOCOUT = 0
                ENDIF
              ENDDO
           ENDDO

           IF (ILOCOUT > 0) THEN
C             PUT OUT REMAINING RECORD
              IRECMRC = IRECMRC + 1
              CALL WRTLIN8(LUNMRC,I1BUF,ILOCOUT,IRECMRC)
          ENDIF

        ELSE
           CALL ERRT(102,'CAN NOT CREATE MRC MODE',MODE)
           GOTO 9999
        ENDIF

        IF (IERR .NE. 0) THEN
           CALL ERRT(102,'WRITING RECORD',IRECIN)
           GOTO 9999
        ENDIF

        
9999    CLOSE(LUNSPI)
        CLOSE(LUNMRC)
        IF(ALLOCATED(STREAMBUF)) DEALLOCATE(STREAMBUF)

        END


C **********************************************************************
C
C COPYTOCCP4_STK  'CP TO MRC' STACKS               JUN 16 ArDean Leith
C
C **********************************************************************
C 
C COPYTOCCP4_STK(LUNSPI,LUNMRC,SPIPAT,ILIST,NIMG, MAXIM,NX,NY,NZ)
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

        SUBROUTINE COPYTOCCP4_STK(LUNSPI,LUNMRC,SPIPAT,
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
        LOGICAL                 :: ISSWABT, SWABT, BOTLEFT

        LOGICAL                 :: isswab
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


C       OPEN NEW MRC STACK FILE FOR STREAM ACCESS
        CALL OPSTREAMFILE(.TRUE.,MRCFILE,DATEXC,LUNMRC,
     &                    'UNFORMATTED','N',
     &                    'MRC OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       FIND IF CURRENTLY SWAPPING BYTES DURING FILE OUTPUT
C       THIS MAY BE DONE IN COMPILATION, HAVE TO ACTUALLY TEST OUTPUT
        ISSWABT = ISSWAB(99)
        
c       write(6,*) '  isswabt,lunflip:', isswabt,lunflip(lunmrc)
        swabt = .false.    

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
        CALL GETLAB(LUNSPI,NX,UNUSED,21,1,SCALEX,IRTFLG)
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
              IF (SWABT) THEN
C                FLIP BYTE ORDERING IN BUFIN
                 CALL FLIPBYTES(BUFIN,FLIPBUF,NX,IRTFLG)
                 WRITE(LUNMRC, POS=IPOSMRC,IOSTAT=IERR) FLIPBUF(1:NX)
              ELSE
                 WRITE(LUNMRC, POS=IPOSMRC,IOSTAT=IERR) BUFIN(1:NX)
              ENDIF

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

