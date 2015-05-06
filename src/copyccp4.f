
C ++********************************************************************
C
C COPYCCP4                 MODIFIED FROM COPYMRC   FEB 02 ArDean Leith         
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
C COPYCCP4(LUNSPI,LUNMRC,NX,NY,NZ)
C                                                                      
C PURPOSE: CONVERTS SPIDER IMAGES TO OR FROM CCP4 FORMAT
C          CRUDELY WRITTEN!!!
C
C NOTES: DATA IN MRC FILE
C	 MODE   TYPES OF PIXEL IN IMAGE
C               0 : INTEGER*1 (UNSIGNED BYTES) 
C               1 : INTEGER*2 (SIGNED) 
C               2 : REALS
C               6 : INTEGER*2 (UNSIGNED)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE COPYCCP4(LUNSPI,LUNMRC,NX,NY,NZ)

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
        LOGICAL                 :: FLIP,ISSWABT,ISSWAB,BOTLEFT
        INTEGER                 :: IVAL
        INTEGER *1              :: I1VAL
        INTEGER *2              :: I2VAL

        INTEGER *2              :: I2V
        INTEGER *1              :: I1V(2),I1TMP
        EQUIVALENCE                (I2V,I1V)

        REAL    *4              :: R4VALIN,R4VALOUT
        INTEGER *1              :: I1VALIN(4),I1VALOUT(4)
        EQUIVALENCE                (R4VALIN,I1VALIN),(R4VALOUT,I1VALOUT)

	CHARACTER(LEN=1)        :: NULL = CHAR(0)
        INTEGER                 :: IERR,LENOPENB,LENOPENF,IRTFLG,MODE
        INTEGER                 :: NSYMBT,MACHST,NE,MAXIM,IOFFSET 
        INTEGER                 :: LENOPEN,NCHAR,IRECMRC,IRECSPI
        INTEGER                 :: IBOTLEF,NOT_USED,IRECINC,ILOCOUT
        INTEGER                 :: IRECIN,ILOCIN,IRECINT
        REAL                    :: RMS,FMINT,FMAXT,FAVT,FSIGT,FN,FNCON
        REAL                    :: UNUSED,SCALE
        INTEGER                 :: IX,IY,IZ,NLET,LOCAT,LOCAST
        INTEGER                 :: I,NSTACKT,ITYPE,NUNUSED,NSTACKOUT
        INTEGER                 :: IMGNUMOUT,NSTACK,IGO,ISTACK,IRECSTK

        INTEGER *8              :: IPOSMRC

        LOGICAL                 :: ASKNAM,FOUROK,WANTSTACK,FOLD
        INTEGER, PARAMETER      :: LUNDOCSEL = 0
        INTEGER, PARAMETER      :: LUNXM     = 0

        integer   :: imax   = -1000000, iymax  = 0, imax2  = 0
        integer   :: ixmax2 = 0, iymax2 = 0
        integer   :: imin   = 100000000, ixmin = 0,iymin = 0, imin2 = 0
        integer   :: ixmin2 = 0, iymin2 = 0

        IERR = 0

C       FIND IF CURRENTLY SWAPPING BYTES DURING FILE OUTPUT
C       THIS MAY BE DONE BY COMPILER, SO HAVE TO ACTUALLY TEST OUTPUT

        ISSWABT = ISSWAB(99)

        IF (FCHAR(4:5) == 'TO')      GOTO 1000


C       COPY FROM MRC TO SPIDER FILE FORMAT --------------- FROM MRC

C       OPEN MRC FILE AS DIRECT ACCESS, UNFORMATTED, RECL=1024 BYTES
        LENOPENB = 1024
        LENOPENF = 1024 / 4
        CALL OPAUXFILE(.TRUE.,MRCFILE,DATEXC,LUNMRC,LENOPENB,'O',
     &                       'MRC (CCP4) INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       READ MRC HEADER 
        CALL REDLIN(LUNMRC,FIXLENBUF,LENOPENF,1)
 
C	PARSE MRC HEADER	
	CALL GETHEDCCP4(FIXLENBUF,NX,NY,NZ,MODE,FMIN,FMAX,
     &                   AV,RMS,NSYMBT,ISSWABT,FLIP,MACHST,IRTFLG)
        IF (IRTFLG == 2) THEN
            CALL ERRT(101,'NOT CURRENT MRC FORMAT, MAY BE PRE 1999 MRC',
     &                NE)
           GOTO 9999
        ENDIF
C       CLOSE MRC FILE
        CLOSE(LUNMRC)

        WANTSTACK = .FALSE. 
        NSTACK    = 0

        IF (NZ > 1) THEN
C          MRC FILE MAY BE VOLUME OR STACK (STUPID FILE FORMAT FOR STACK)
           PROMPT = 
     &      'OUTPUT VOLUME OR TEMPLATE FOR IMAGE STACK (E.G. STK@*)~'

           CALL FILERD(FILPAT,NLET,NULL,PROMPT,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           LOCAT     = INDEX(FILPAT(1:NLET),'@')   
           LOCAST    = INDEX(FILPAT(1:NLET),'*')
           WANTSTACK = (LOCAST > LOCAT ) 
           ASKNAM    = .FALSE.    ! ALREADY ASKED OUTPUT FILENAME

           IF (WANTSTACK) THEN
C             OPEN FIRST SPIDER OUTPUT FILE	

              IGO = 1
              CALL RDPRI1S(IGO,NOT_USED,
     &                    'FIRST IMAGE NUMBER IN STACK',IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999

              NSTACK = NZ
              ALLOCATE(ILIST(NSTACK),STAT=IRTFLG)
              IF (IRTFLG > 0) THEN
                 CALL ERRT(46,'COPYCCP4; ILIST',NSTACK)
                 GOTO 9999
              ENDIF

C             MAKE LIST OF STACKED FILE NUMBERS
              DO I = IGO,NSTACK
                 ILIST(I) = I
              ENDDO

              NSTACKT = -NSTACK    ! USE ILIST FOR STACKED IMG NUMBERS
              FOUROK  = .FALSE.    ! NOT FOURIER
              NZ      = 1          ! NOT A VOLUME
              ITYPE   = 1

C             OPEN SPIDER STACK OUTPUT FILE	
              CALL OPFILES(0,LUNSPI,LUNDOCSEL,LUNXM,
     &            ASKNAM,FILOUT,NLET, 'U',
     &            ITYPE,NX,NY,NZ,MAXIM,
     &            FILPAT,
     &            FOUROK, ILIST,NSTACKT, 
     &            NUNUSED,NSTACKOUT, IMGNUMOUT, IRTFLG)
              !write(6,*) 'nstackout, imgnumout:', nstackout, imgnumout 

           ELSE
C             OPEN SPIDER VOLUME OUTPUT FILE	
              ITYPE  = 3
              MAXIM  = 0
              CALL OPFILEC(0,ASKNAM,FILPAT,LUNSPI,'U',ITYPE,NX,NY,NZ,
     &                 MAXIM,' ',.FALSE.,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999
           ENDIF
        ELSE
C          OPEN SPIDER IMAGE OUTPUT FILE	
           ITYPE  = 1
           MAXIM  = 0
           CALL OPFILEC(0,.TRUE.,FILPAT,LUNSPI,'U',ITYPE,NX,NY,NZ,
     &                 MAXIM,'SPIDER IMAGE OUTPUT',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
        ENDIF


C       EXTRACT DATA FROM MRC FILE AFTER HEADER & PUT IN SPIDER FILE
        IOFFSET = 1024 + NSYMBT

        IF (MODE == 0) THEN
           ALLOCATE(I1STREAMBUF(NX), STREAMBUF(NX),STAT=IRTFLG)

        ELSEIF (MODE == 1 .OR. MODE == 6)  THEN
           ALLOCATE(I2STREAMBUF(NX), STREAMBUF(NX),STAT=IRTFLG)

        ELSEIF (MODE == 2 )  THEN
           ALLOCATE(STREAMBUF(NX),STAT=IRTFLG)

        ELSE 
           CALL ERRT(102,'UNSUPPORTED MRC MODE',MODE)
        ENDIF 

        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'COPYCCP4; **STREAMBUF',NX)
           GOTO 9999
        ENDIF

C       OPEN MRC FILE FOR STREAM ACCESS
        CALL OPSTREAMFILE(.FALSE.,MRCFILE,NULL,LUNMRC,
     &                    'UNFORMATTED','O',' ',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &           'FLIP BYTE ORDERING? (Y/N), INVERT TOP/BOTTOM? (Y/N)',
     &            NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (ANS(1:1) == 'Y') FLIP = .NOT. FLIP

        ! BOTLEFT IS USUAL MRC BOTTOM IN CURRENT FORMAT !!
        BOTLEFT = (NCHAR > 2 .AND. (INDEX(ANS(3:NCHAR),'Y')) > 0) 

        ISTACK  = 1

        DO  ! POSSIBLE LOOP OVER MRC STACK -------------------------
         IRECSPI = 0
         IRECINC = 1
         IRECSTK = (ISTACK-1) * NY

         IF (BOTLEFT) THEN
C           INVERT TOP & BOT OF EACH IMAGE OR EACH IMAGE WITHIN VOLUME
            IRECSPI = NY*NZ + 1    ! NZ HAS BEEN SET = 1 IF STACK
            IRECINC = -1
         ENDIF

         IF (MODE == 0) THEN
C          SIGNED 8 BIT INTEGER CCP4 (MRC) INPUT FILE

           DO IZ = 1,NZ
             DO IY = 1,NY
                IREC    = IRECSTK + (IZ  -1) * NY + IY

C               IPOSMRC = IOFFSET + (IREC-1) * NX + 1

                IPOSMRC = (IREC-1)    ! KLUDGE FOR INTEGER *8 PRESERVE
                IPOSMRC = IPOSMRC * NX
                IPOSMRC = IOFFSET + IPOSMRC + 1

                READ(LUNMRC, POS=IPOSMRC,IOSTAT=IERR) I1STREAMBUF

                DO IX = 1,NX
                   IVAL = I1STREAMBUF(IX)
                   IF (IVAL < 0) IVAL = 256 + IVAL
                   STREAMBUF(IX) = IVAL
                ENDDO
  
                !if (irec == 1) then
                !   write(6,*)' ',iposmrc,i1streambuf(1),streambuf(1)
                !endif

                IRECSPI = IRECSPI + IRECINC
                !write(6,*) ' irecspi:',irecspi
                CALL WRTLIN(LUNSPI,STREAMBUF,NX,IRECSPI)

              ENDDO ! END OF: IY = 1,NY
           ENDDO    ! END OF: IZ = 1,NZ

        ELSEIF (MODE == 1 .OR. MODE == 6) THEN
C         16 BIT INTEGER CCP4 (MRC) INPUT FILE

          FOLD = (MODE == 1)
          !write(6,*) 'Flip & fold:',flip,fold

          DO IZ = 1,NZ
             DO IY = 1,NY
                IREC    = IRECSTK + (IZ  -1) * NY + IY

c               IPOSMRC = IOFFSET + (IREC-1) * NX * 2 + 1

                IPOSMRC = (IREC-1)    ! KLUDGE FOR INTEGER *8 PRESERVE
                IPOSMRC = IPOSMRC * NX
                IPOSMRC = IPOSMRC * 2
                IPOSMRC = IOFFSET + IPOSMRC + 1

                READ(LUNMRC, POS=IPOSMRC,IOSTAT=IERR) I2STREAMBUF

                IF (FLIP .AND. FOLD) THEN
C                  INVERT BYTE ORDER & CONVERT SIGNED INTEGER TO UNSIGNED
                   DO IX = 1,NX
                      I2V           = I2STREAMBUF(IX)
                      I1TMP         = I1V(1)
                      I1V(1)        = I1V(2)
                      I1V(2)        = I1TMP

C                     FOLD CONVERTS SIGNED INTEGER TO UNSIGNED
                      IF (I2V < 0) I2V = 65536 + I2V

                      STREAMBUF(IX) = I2V
                   ENDDO

                ELSEIF (FOLD) THEN
C                  CONVERT SIGNED INTEGER TO UNSIGNED
                   DO IX = 1,NX
                      I2V = I2STREAMBUF(IX)
                      IF (I2V < 0) I2V = 65536 + I2V
                      STREAMBUF(IX) = I2V
                   ENDDO

                ELSEIF (FLIP) THEN
C                  INVERT BYTE ORDER
                   DO IX = 1,NX
                      I2V           = I2STREAMBUF(IX)
                      I1TMP         = I1V(1)
                      I1V(1)        = I1V(2)
                      I1V(2)        = I1TMP
                      STREAMBUF(IX) = I2V
                   ENDDO

                ELSE
C                  NO CONVERSION
                   STREAMBUF = I2STREAMBUF
                ENDIF

#ifdef NEVER
                do ix = 1,nx
                   if (streambuf(ix) < imin) then
                      imin2  = imin
                      ixmin2 = ixmin
                      iymin2 = iymin
                      imin   = streambuf(ix)
                      ixmin  = ix
                      iymin  = iy
                   endif
                   if (streambuf(ix) > imax) then
                      imax2  = imax
                      ixmax2 = ixmax
                      iymax2 = iymax
                      imax   = streambuf(ix)
                      ixmax  = ix
                      iymax  = iy
                   endif
                enddo
#endif

C               PUT OUT COMPLETED RECORD
                IRECSPI = IRECSPI + IRECINC
                CALL WRTLIN(LUNSPI,STREAMBUF,NX,IRECSPI)

              ENDDO
           ENDDO

#ifdef NEVER
           write(6,*) 'maxs: ',ixmax,iymax,imax
           write(6,*) 'maxs2:',ixmax2,iymax2,imax2
           write(6,*) 'mins: ',ixmin,iymin,imin
           write(6,*) 'mins2:',ixmin2,iymin2,imin2
#endif

        ELSEIF (MODE == 2) THEN
C          32 BIT FOATING POINT CCP4 MRC INPUT FILE

           DO IZ = 1,NZ
             DO IY = 1,NY
                IREC    = IRECSTK + (IZ  -1) * NY + IY
c               IPOSMRC = IOFFSET + (IREC-1) * NX * 4 + 1
                      
                IPOSMRC = (IREC-1)    ! KLUDGE FOR INTEGER *8 PRESERVE
                IPOSMRC = IPOSMRC * NX
                IPOSMRC = IPOSMRC * 4
                IPOSMRC = IOFFSET + IPOSMRC + 1

                READ(LUNMRC, POS=IPOSMRC,IOSTAT=IERR) STREAMBUF
                !if (irec ==1) write(6,*) ' Val:',iposmrc, streambuf(1)
                !write(6,*) ' irec,iposmrc:',irec,iposmrc

                IF (FLIP) THEN
C                  INVERT BYTE ORDER
                   DO IX = 1,NX
                      R4VALIN       = STREAMBUF(IX)
                      I1VALOUT(1)   = I1VALIN(4)
                      I1VALOUT(2)   = I1VALIN(3)
                      I1VALOUT(3)   = I1VALIN(2)
                      I1VALOUT(4)   = I1VALIN(1)
                      STREAMBUF(IX) = R4VALOUT
                   ENDDO
                ENDIF

C               PUT OUT COMPLETED RECORD
                IRECSPI = IRECSPI + IRECINC
                CALL WRTLIN(LUNSPI,STREAMBUF,NX,IRECSPI)

              ENDDO
           ENDDO

        ELSE
           CALL ERRT(102,'CAN NOT COPY MRC MODE',MODE)
        ENDIF

        IF (.NOT. WANTSTACK) EXIT    ! FINISHED IF NOT A MRC STACK

C       OPEN NEXT STACKED OUTPUT FILE 
        CALL NEXTFILE(ISTACK, ILIST, 
     &                FOUROK,LUNXM,
     &                NSTACK,MAXIM,   
     &                LUNSPI,0,
     &                FILPAT,'N',
     &                IMGNUMOUT, IRTFLG) 

        IF (ISTACK > NSTACK) EXIT    ! FINISHED 
        IF (IRTFLG == -99) THEN
           CALL ERRT(102,'INSUFFICIENT OUTPUT FILE NAMES',NSTACK)
           EXIT         
        ELSEIF (IRTFLG .NE. 0) THEN
           EXIT                      ! ERROR
        ENDIF

       ENDDO   ! END OF STACK LOOP --------------------------

       GOTO 9999

	



C      COPY FROM SPIDER TO MRC FILE FORMAT ----------------- TO MRC

1000   CONTINUE

C	OPEN NEW MRC FILE FOR DIRECT ACCESS, RECORD LENGTH 1024 BYTES
        LENOPENB = 1024
        LENOPENF = LENOPENB / 4
        CALL OPAUXFILE(.TRUE.,MRCFILE,DATEXC,LUNMRC,LENOPENB,'U',
     &                 'MRC OUTPUT',.TRUE.,IRTFLG)

C       FLIP BYTES DURING MRC FILE OUTPUT
        IF (ISSWABT) THEN
           CALL LUNSETFLIP(LUNMRC,1,IRTFLG)
           ISSWABT = .FALSE.
        ENDIF


        IVAL    = 8
        IBOTLEF = 0
        CALL RDPRI2S(IVAL,IBOTLEF,NOT_USED,
     &  'MRC DATA LENGTH (8/32 BITS), INVERT IMAGE TOP/BOTTOM =1 (0/1)',
     &     IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        MODE = 2
        IF (IVAL == 8) MODE = 0
        BOTLEFT = (IBOTLEF > 0)   ! USUAL MRC BOTTOM!!

C	CREATE A NEW HEADER FOR THE CCP4 FILE
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

C	CREATE HEADER. (NOTE: FMIN, FMAX, AV ARE SAME AS SPIDER IMAGE)
        CALL SETHEDCCP4(FIXLENBUF, NX, NY, NZ,
     &            FMINT,FMAXT,FAVT,FSIGT,SCALE,MODE,ISSWABT,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C	WRITE HEADER OF 1024 BYTES (256 FLOATS) TO MRC FILE
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

C             PUT ROW OUT TO CCP4 FILE
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


C       -------------- ISSWAB ----------------------------------------

        LOGICAL FUNCTION ISSWAB(LUN)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
 
        INTEGER,DIMENSION(3) ::  IVAL
        CHARACTER(LEN=12) ::     CVAL,CVALIN
        EQUIVALENCE(IVAL,CVAL)

        CHARACTER(LEN=MAXNAM) :: FILNAM
        LOGICAL ::               VERBOSE_SAVE

        CHARACTER(LEN=1) ::      NULL = CHAR(0)

C       DO NOT ECHO FILE OPENING
        VERBOSE_SAVE = VERBOSE
        VERBOSE      = .FALSE.

        FILNAM = 'TMP_JNK_SCRATCH'
        CALL OPAUXFILE(.FALSE.,FILNAM,NULL,LUN,12,'U',' ',.TRUE.,IRTFLG)

        CVAL(1:1)   = CHAR(0)
        CVAL(2:2)   = CHAR(0)
        CVAL(3:3)   = CHAR(0)
        CVAL(4:4)   = CHAR(4)
        CVAL(5:5)   = CHAR(48)
        CVAL(6:6)   = CHAR(48)
        CVAL(7:7)   = CHAR(49)
        CVAL(8:8)   = CHAR(50)
        CVAL(9:9)   = CHAR(0)
        CVAL(10:10) = CHAR(0)
        CVAL(11:11) = CHAR(0)
        CVAL(12:12) = CHAR(4)

        CALL WRTLIN(LUN,IVAL,3,1)
        CLOSE(LUN)

        CALL OPAUXFILE(.FALSE.,FILNAM,NULL,LUN,0,'O',' ',.TRUE.,IRTFLG)

        READ(LUN,*) CVALIN

        CLOSE(LUN,STATUS='DELETE')

c       WRITE(NOUT,*) 'CVALIN: ',CVALIN,' == ',CVAL
 
        ISSWAB   = (CVALIN(8:8) .NE. CVAL(8:8))

c        IF (ISSWAB) THEN
c           WRITE(NOUT,*) 'NON-NATIVE BYTE ORDER '
c        ELSE
c           WRITE(NOUT,*) 'NATIVE BYTE ORDER'
c        ENDIF

        VERBOSE = VERBOSE_SAVE

        END

C       --------------  OPSTREAMFILE ------------------------------

         SUBROUTINE OPSTREAMFILE(ASKNAME,FILNAM,EXTENT,LUNT,
     &                           FORMVAR, DISP, 
     &                           PROMPTT,CALLERRT,IRTFLG)


        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'

        
        LOGICAL           :: ASKNAME
        CHARACTER(LEN=*)  :: FILNAM,EXTENT
        INTEGER           :: LUNT
        CHARACTER(LEN=11) :: FORMVAR
        CHARACTER(LEN=*)  :: DISP,PROMPTT
        LOGICAL           :: CALLERRT
        INTEGER           :: IRTFLG

        LOGICAL           :: EX
        CHARACTER(LEN=96) :: PROMPT
        CHARACTER(LEN=80) :: EXTEN
        CHARACTER(LEN=7)  :: STATVAR

        INTEGER           :: ICOMM,MYPID,MPIERR,LENP,NCHAR
        INTEGER           :: LNBLNKN
        INTEGER           :: LENE,IRTFLGT,LUN,IDUM,LENOPEN,LENOPENFILE
        INTEGER           :: LENOPN,LENREC

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID  #ifdef USE_MPI

C       SET DEFAULT ERROR RETURN
        IRTFLG = 1

C       DO NOT WANT TO RETURN EXTEN
        EXTEN = EXTENT

C       INPUT FILE NAME (IF EXTEN EXISTS IT IS ADDED)

        IF (ASKNAME) THEN
C          SET PROMPT TO ALLOW FILE EXTENSION ON INPUT
           LENP   = LEN(PROMPTT)
           LENP   = MIN(LENP,93)
           PROMPT = PROMPTT(1:LENP) // '~9' 

           CALL FILERD(FILNAM,NCHAR,EXTEN,PROMPT(1:LENP+2),IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ELSE
C          MAY WANT TO ADD EXTENT TO FILNAM
           NCHAR = LNBLNKN(FILNAM)
           LENE  = LNBLNKN(EXTENT)
           IF (LENE > 0) THEN
C             ADD THE EXTENSION THAT IS SENT TO FILNAM
              CALL FILNAMANDEXT(FILNAM,EXTEN,FILNAM,NCHAR,
     &                          .TRUE.,IRTFLGT)
           ENDIF
        ENDIF

        LUN = ABS(LUNT)
        IF ((LUN <= 0 .OR. LUN > 100) .AND.
     &     (LUN .NE. 103)) THEN
C          LUN=103 USED IN  SYMPARTEXT 
           CALL ERRT(102,'IN SOURCE CODE, LUN MUST BE 1...100',LUN)
           RETURN
        ENDIF

        IF (LUN > 0 .AND. LUN <= 100) THEN
C          ZERO THE FLAGS USED IN REDLIN/WRTLIN
           CALL LUNSETLUNS(LUN,0,0,LUN,0,IRTFLGT)
 
C          MAKE SURE THIS IS NOT TREATED AS INLINE FILE
           CALL CLOSEINLN(LUN,IRTFLGT)
        ENDIF

C       SET STATUS FOR OPEN
        STATVAR = 'NEW'

        IF (DISP(1:1) == 'N' .OR. DISP(1:1) == 'U') 
     &     STATVAR = 'REPLACE'

        IF (DISP(1:1) == 'S') STATVAR = 'SCRATCH'

        IF (DISP(1:1) == 'O') THEN
C          CHECK FOR FILE EXISTENCE 
           IF (MYPID <= 0) THEN
              INQUIRE (FILE=FILNAM(1:NCHAR),EXIST=EX,IOSTAT=IRTFLGT) 
           ENDIF

#ifdef USE_MPI
           CALL BCAST_MPI('OPSTREAMFILE','EX',           EX,1,'L',ICOMM)
           CALL BCAST_MPI('OPSTREAMFILE','IRTFLGT', IRTFLGT,1,'I',ICOMM)
#endif

           IF (IRTFLGT .NE. 0) THEN
              WRITE(NOUT,*) '*** INQUIRY ERROR'
              IF (CALLERRT)  CALL ERRT(4,'OPSTREAMFILE',IDUM)
              RETURN
        
           ELSEIF (.NOT. EX) THEN
              WRITE(NOUT,*) '*** FILE DOES NOT EXIST: ',FILNAM(1:NCHAR)
              IF (CALLERRT)  CALL ERRT(100,'OPSTREAMFILE',IDUM)
              RETURN

           ENDIF
           STATVAR = 'OLD'
        ENDIF

C       OPEN FILE FOR STREAM ACCESS

C       COMPUTE RECL UNITS (DIFFERS WITH OS &A COMPILER FLAGS)
        LENOPN = LENOPENFILE(LENREC)

        IF (MYPID <= 0) THEN
           IF (STATVAR == 'SCRATCH') THEN
	      OPEN(UNIT=LUN,STATUS=STATVAR,
     &             FORM=FORMVAR, ACCESS='STREAM',
     &             IOSTAT=IRTFLGT)
           ELSE
	      OPEN(UNIT=LUN,FILE=FILNAM(1:NCHAR),STATUS=STATVAR,
     &             FORM=FORMVAR, ACCESS='STREAM', 
     &             IOSTAT=IRTFLGT)
           ENDIF
        ENDIF

#ifdef USE_MPI
        CALL BCAST_MPI('OPSTREAMFILE','IRTFLGT', IRTFLGT,1, 'I',ICOMM)
#endif


        IF (IRTFLGT .NE. 0) THEN
           WRITE(NOUT,90) FORMVAR(1:1), FILNAM(:NCHAR)
 90        FORMAT(' ERROR OPENING (',A1,'): ',A)
           IF (CALLERRT) CALL ERRT(102,'OPSTREAMFILE',IRTFLGT)
           RETURN
        ENDIF

        IF (VERBOSE .AND. MYPID <= 0) THEN
           WRITE(NOUT,91) FORMVAR(1:1), FILNAM(:NCHAR)
 91        FORMAT('  OPENED (',A1,'): ',A)
        ENDIF

        IRTFLG = 0

        END


