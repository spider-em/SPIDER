
C ++********************************************************************
C
C COPYTOCCP4               MODIFIED FROM COPYMRC   FEB 02 ArDean Leith         
C                          ISSWAB ADDED            JUL 02 ArDean Leith
C                          FLIP QUESTION           MAR 03 ArDean Leith
C                          BAD IRECMRC4 & FLIP     SEP 03 ArDean Leith
C                          SCALING                 JAN 05 ArDean Leith
C                          I*8                     SEP 08 ArDean Leith
C                          NPIX8                   DEC 08 ArDean Leith
C                          BOTLEFT OPTION          MAY 12 ArDean Leith
C                          STREAM IO               FEB 13 ArDean Leith
C                          VOL BUG                 JUN 13 ArDean Leith
C                          VOL BUG FIXED           JUL 13 ArDean Leith
C                          MODE 6 STACK SUPPORT    SEP 14 ArDean Leith
C                          IPOSMRC INTEGER *8      JAN 15 ArDean Leith
C                          BOTLEFT DEFAULT         JUL 15 ArDean Leith
C                          2015 STACK SUPPORT      JUL 15 ArDean Leith
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
C COPYTOCCP4(LUNSPI,LUNMRC,NX,NY,NZ)
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
        INTEGER                 :: IBOTLEF,NOT_USED,ILOCOUT
        INTEGER                 :: IRECIN,ILOCIN,IRECINT,NSYMBYT
        REAL                    :: RMS,FMINT,FMAXT,FAVT,FSIGT,FN,FNCON
        REAL                    :: UNUSED,SCALE,SCALEX,SCALEY,SCALEZ
        INTEGER                 :: IX,IY,IZ,NLET
        INTEGER                 :: I,NSTACKT,ITYPE,NUNUSED,NSTACKOUT
        INTEGER                 :: IMGNUMOUT,NSTACK,IGO
        INTEGER                 :: NIMG,MZ

        LOGICAL                 :: ASKNAM,FOUROK,WANTSTACK,FOLD

        IERR = 0

C       OPEN NEW MRC FILE FOR DIRECT ACCESS, RECORD LENGTH 1024 BYTES
        LENOPENB = 1024
        LENOPENF = LENOPENB / 4
        CALL OPAUXFILE(.TRUE.,MRCFILE,DATEXC,LUNMRC,LENOPENB,'U',
     &                 'MRC OUTPUT',.TRUE.,IRTFLG)

        IVAL    = 8
        IBOTLEF = 1
        CALL RDPRI2S(IVAL,IBOTLEF,NOT_USED,
     &               'MRC DATA LENGTH (8/32 BITS)',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        MODE = 2
        IF (IVAL == 8) MODE = 0

        BOTLEFT = (IBOTLEF > 0)   ! USUAL MRC BOTTOM!!

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
                 CALL REDLIN(LUNSPI,BUFIN,NX,IRECINT)
              ELSE
                 CALL REDLIN(LUNSPI,BUFIN,NX,IRECIN)
              ENDIF

C             PUT ROW OUT TO MRC FILE
              DO ILOCIN=1,NX
                ILOCOUT            = ILOCOUT + 1
                FIXLENBUF(ILOCOUT) = BUFIN(ILOCIN)

                IF (ILOCOUT >= LENOPENF) THEN
C                  PUT OUT COMPLETED RECORD

                   IRECMRC = IRECMRC + 1

                   CALL WRTLIN(LUNMRC,FIXLENBUF,LENOPENF,IRECMRC)
                   ILOCOUT = 0
                ENDIF
              ENDDO
           ENDDO

           IF (ILOCOUT > 0) THEN
C             PUT OUT REMAINING RECORD
              IRECMRC = IRECMRC + 1
              CALL WRTLIN(LUNMRC,FIXLENBUF,ILOCOUT,IRECMRC)
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


C ++********************************************************************
C
C COPYTOCCP4_STK  'CP TO MRC' STACKS               JUN 15 ArDean Leith
C
C **********************************************************************
C 
C COPYTOCCP4_STK(LUNSPI,LUNMRC,SPIPAT,ILIST,NIMG, MAXIM,NX,NY,NZ)
C                                                                      
C PURPOSE: CONVERTS SPIDER IMAGE STACK TO MRC FORMAT, CRUDELY WRITTEN!!
C
C NOTES: DATA MODE IN MRC FILE : 2  (32 BIT FLOATING POINT)
C        ONLY CREATES 32 BIT REAL STACKS AS NORMALIZATION DIFFICULT FOR
C        LARGE MRC STACKS.
C
C        DMAX < DMIN                         MAX & MIN UNDETERMINED
C        DMEAN < (SMALLER OF DMIN and DMAX)  DMEAN     UNDETERMINED
C        RMS < 0.0                           RMS       UNDETERMINED
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE COPYTOCCP4_STK(LUNSPI,LUNMRC,SPIPAT,
     &                          ILIST,NIMG,MAXIM,
     &                          NX,NY,NZ)
 
C       COPY STACKS FROM SPIDER TO MRC FILE FORMAT 

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
 
        INTEGER                 :: IERR
        COMMON /IOERR/  IERR       ! LEGACY REDLIN ERROR FLAG

        REAL                    :: BUFIN
        COMMON /IOBUF/  BUFIN(NBUFSIZ)

        INTEGER                 :: LUNSPI,LUNMRC, NIMG,MAXIM,NX,NY,NZ
        CHARACTER(LEN=*)        :: SPIPAT
        INTEGER                 :: ILIST(*)

        REAL                    :: FIXLENBUF(256)

        CHARACTER(LEN=MAXNAM)   :: MRCFILE
        LOGICAL                 :: ISSWABT
        LOGICAL                 :: isswab

        INTEGER                 :: LENOPENB,LENOPENF,IRTFLG,MODE
        INTEGER                 :: NSYMBT,MACHST,NE,IOFFSET 
        INTEGER                 :: LENOPEN,NCHAR,IRECMRC
        INTEGER                 :: NOT_USED,ILOCOUT
        INTEGER                 :: IRECIN,ILOCIN,IRECINT,NZMRC
        REAL                    :: RMS,FMINT,FMAXT,FAVT,FSIGT,FN,FNCON
        REAL                    :: UNUSED,SCALE,SCALEX,SCALEY,SCALEZ
        INTEGER                 :: NLET,IMGNUMOUT,IGO,MZ
        INTEGER                 :: NINDX1,NGOT,IBITS,NSYMBYT,ISPG

        LOGICAL                 :: ASKNAM,WANTSTACK,FOLD
        LOGICAL, PARAMETER      :: FOUROK = .FALSE.
        INTEGER, PARAMETER      :: LUNXM = 0

C       OPEN NEW MRC FILE FOR DIRECT ACCESS, RECORD LENGTH 1024 BYTES
        LENOPENB = 1024
        LENOPENF = LENOPENB / 4
        CALL OPAUXFILE(.TRUE.,MRCFILE,DATEXC,LUNMRC,LENOPENB,'U',
     &                 'MRC OUTPUT',.TRUE.,IRTFLG)

        IGO   = 1
        IBITS = 32               ! FOR FUTURE
        CALL RDPRI2S(IGO,IBITS,NOT_USED,
     &               'INITIAL MRC IMAGE NUMBER', IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        MODE  = 2                ! FLOATING POINT IMAGE

C       FIND IF CURRENTLY SWAPPING BYTES DURING FILE OUTPUT
C       THIS MAY BE DONE IN COMPILATION, HAVE TO ACTUALLY TEST OUTPUT
        ISSWABT = ISSWAB(99)

        IF (ISSWABT) THEN
C          FLIP BYTES DURING MRC FILE OUTPUT
           CALL LUNSETFLIP(LUNMRC,1,IRTFLG)
           ISSWABT = .FALSE.
        ENDIF

C       CREATE A NEW HEADER FOR THE MRC FILE
C       DMAX < DMIN                         MAX & MIN UNDETERMINED
C       DMEAN < (SMALLER OF DMIN and DMAX)  DMEAN     UNDETERMINED
C       RMS < 0.0                           RMS       UNDETERMINED
C       ISPG   == 0                IMAGE OR IMAGE STACK
C       ISPG   == 1                VOLUMES
C       ISPG   == 401              STACK OF EM VOLUMES

        FMINT =  0    ! UNDETERMINED
        FMAXT =  0    ! UNDETERMINED
        FAVT  = -1    ! UNDETERMINED
        FSIGT = -1    ! UNDETERMINED

C       TRY TO GET SPIDER IMAGE SCALE VALUE (MAY NOT BE SET?)
        CALL GETLAB(LUNSPI,NX,UNUSED,21,1,SCALEX,IRTFLG)
        SCALEY  = SCALEX
        SCALEZ  = SCALEX

        NSYMBYT = 0            ! NO ADDED HEADER LOCATIONS

C       CREATE MRC HEADER. 
        CALL SETHEDCCP4(FIXLENBUF, NX,NY,NZ, FMINT,FMAXT,FAVT,FSIGT,
     &                  SCALEX,SCALEY,SCALEZ, MODE,
     &                  ISSWABT,NSYMBYT,NIMG,MZ,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       WRITE HEADER OF 1024 BYTES (256 FLOATS) TO MRC FILE
        CALL WRTLIN(LUNMRC,FIXLENBUF,LENOPENF,1) 

C       SET STARTING RECORD FOR MRC DATA
        IRECMRC = IGO + 1   ! SKIPS ONE HEADER RECORD
        ILOCOUT = 0
        IERR    = 0

        NINDX1  = 1

        DO      ! LOOP OVER ALL IMAGES IN SPIDER STACK

           DO IRECIN = 1,NY * NZ

C             READ EACH ROW OF SPIDER INPUT FILE 
              IRECINT = (NY * NZ) - IRECIN + 1  
              CALL REDLIN(LUNSPI,BUFIN,NX,IRECINT)

C             PUT ROW OUT TO MRC FILE
              DO ILOCIN=1,NX
                 ILOCOUT            = ILOCOUT + 1
                 FIXLENBUF(ILOCOUT) = BUFIN(ILOCIN)

                 IF (ILOCOUT >= LENOPENF) THEN
C                   PUT OUT COMPLETED RECORD

                    IRECMRC = IRECMRC + 1

                    CALL WRTLIN(LUNMRC,FIXLENBUF,LENOPENF,IRECMRC)
                    ILOCOUT = 0
                 ENDIF
              ENDDO
           ENDDO

           IF (ILOCOUT > 0) THEN
C             PUT OUT REMAINING RECORD IN THIS IMAGE
              IRECMRC = IRECMRC + 1
              CALL WRTLIN(LUNMRC,FIXLENBUF,ILOCOUT,IRECMRC)
           ENDIF
        
           IF (IERR .NE. 0) THEN
              CALL ERRT(102,'WRITING RECORD',IRECIN)
              GOTO 9999
           ENDIF

C          OPEN NEXT SPIDER INPUT FILE, UPDATES NINDX1 
           CALL NEXTFILE(NINDX1, ILIST, 
     &                  FOUROK, LUNXM,
     &                  NGOT,   MAXIM,   
     &                  LUNSPI, 0,
     &                  SPIPAT, 'O',
     &                  IMGNUMOUT, IRTFLG) 

           IF (NINDX1 > MAXIM) EXIT    ! FINISHED 
           IF (IRTFLG .NE. 0)  EXIT    ! ERROR

       ENDDO
       
9999   CLOSE(LUNSPI)
       CLOSE(LUNMRC)

       END




