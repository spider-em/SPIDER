
C ++********************************************************************
C
C COPYFROMCCP4   MODIFIED FROM COPYMRC             FEB 02 ArDean Leith         
C                ISSWAB ADDED                      JUL 02 ArDean Leith
C                FLIP QUESTION                     MAR 03 ArDean Leith
C                BAD IRECMRC4 & FLIP               SEP 03 ArDean Leith
C                SCALING                           JAN 05 ArDean Leith
C                I*8                               SEP 08 ArDean Leith
C                NPIX8                             DEC 08 ArDean Leith
C                BOTLEFT OPTION                    MAY 12 ArDean Leith
C                STREAM IO                         FEB 13 ArDean Leith
C                VOL BUG                           JUN 13 ArDean Leith
C                VOL BUG FIXED                     JUL 13 ArDean Leith
C                MODE 6 STACK SUPPORT              SEP 14 ArDean Leith
C                IPOSMRC INTEGER *8                JAN 15 ArDean Leith
C                BOTLEFT DEFAULT                   JUL 15 ArDean Leith
C                2015 STACK SUPPORT                JUL 15 ArDean Leith
C                STACK END BUG                     OCT 15 ArDean Leith
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
C COPYFROMCCP4(LUNSPI,LUNMRC,NX,NY,NZ)
C                                                                      
C PURPOSE: CONVERTS SPIDER IMAGES TO OR FROM MRC FORMAT
C          CRUDELY WRITTEN!!!
C
C NOTES: DATA IN MRC FILE
C        MODE   TYPES OF PIXEL IN IMAGE
C               0 : INTEGER*1 (UNSIGNED BYTES) 
C               1 : INTEGER*2 (SIGNED) 
C               2 : REALS
C               6 : INTEGER*2 (UNSIGNED)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE COPYFROMCCP4(LUNSPI,LUNMRC,NX,NY,NZ)

C       COPY FROM MRC TO SPIDER FILE FORMAT --------------- FROM MRC

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
 
        COMMON /IOERR/  IERR       ! LEGACY REDLIN ERROR FLAG

        REAL                    :: BUFIN
        COMMON /IOBUF/  BUFIN(NBUFSIZ)

        INTEGER                 :: LUNSPI,LUNMRC, NXT,NYT,NZT

        INTEGER                 :: NX,NY,NZ
        REAL,       ALLOCATABLE :: STREAMBUF(:)
        INTEGER *1, ALLOCATABLE :: I1STREAMBUF(:)
        INTEGER *2, ALLOCATABLE :: I2STREAMBUF(:)
        INTEGER,    ALLOCATABLE :: ILIST(:)
        REAL                    :: BUF(NBUFSIZ),FIXLENBUF(256)
        INTEGER *1              :: I1BUF(1024)
        COMMON                     BUF,FIXLENBUF,I1BUF

        CHARACTER(LEN=MAXNAM)   :: MRCFILE,FILOUT
        CHARACTER(LEN=8)        :: ANS
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
        INTEGER                 :: LENOPEN,NCHAR,IRECSPI
        INTEGER                 :: IBOTLEF,NOT_USED,IRECINC,ILOCOUT
        INTEGER                 :: IRECIN,ILOCIN,IRECINT
        REAL                    :: RMS,FMINT,FMAXT,FAVT,FSIGT,FN,FNCON
        REAL                    :: UNUSED,SCALE
        INTEGER                 :: IX,IY,IZ,NLET,LOCAT,NIMG
        INTEGER                 :: I,NSTACKT,ITYPE,NUNUSED,NSTACKOUT
        INTEGER                 :: IMGNUMOUT,IGO,NINDX,IRECSTK
        INTEGER                 :: ISPG,MZ

        INTEGER *8              :: IPOSMRC

        LOGICAL                 :: ASKNAM,FOUROK,WANTSTACK,FOLD
        INTEGER, PARAMETER      :: LUNDOCSEL = 0
        INTEGER, PARAMETER      :: LUNXM     = 0

        integer   :: imax   = -1000000, iymax  = 0, imax2  = 0
        integer   :: ixmax2 = 0, iymax2 = 0
        integer   :: imin   = 100000000, ixmin = 0,iymin = 0, imin2 = 0
        integer   :: ixmin2 = 0, iymin2 = 0

        IERR = 0

        NX   = NXT
        NY   = NYT
        NZ   = NZT

C       COPY FROM MRC TO SPIDER FILE FORMAT --------------- FROM MRC

C       OPEN MRC FILE AS DIRECT ACCESS, UNFORMATTED, RECL=1024 BYTES
        LENOPENB = 1024
        LENOPENF = 1024 / 4
        CALL OPAUXFILE(.TRUE.,MRCFILE,DATEXC,LUNMRC,LENOPENB,'O',
     &                       'MRC INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       READ MRC HEADER 
        CALL REDLIN(LUNMRC,FIXLENBUF,LENOPENF,1)
 
C       FIND IF CURRENTLY SWAPPING BYTES DURING FILE OUTPUT
C       THIS MAY BE DONE BY COMPILER, SO HAVE TO ACTUALLY TEST OUTPUT

        ISSWABT = ISSWAB(99)

C       PARSE MRC HEADER        
        CALL GETHEDCCP4(FIXLENBUF,NX,NY,NZ,MODE,FMIN,FMAX,
     &                  AV,RMS,NSYMBT,ISSWABT,FLIP,MACHST,
     &                  ISPG,MZ,IRTFLG)
        IF (IRTFLG == 2) THEN
            CALL ERRT(101,
     &          'NOT CURRENT MRC FORMAT, CAN NOT COPY PRE-2000 MRC', NE)
           GOTO 9999
        ENDIF
C       CLOSE MRC FILE
        CLOSE(LUNMRC)

        WANTSTACK = .FALSE. 

        IF     (ISPG == 0 .AND. MZ == 1 .AND.  NZ > 1) THEN
C          SPIDER IMAGE STACK OUTPUT FILE (SOME PRE 2015 STACKS)
           WANTSTACK = .TRUE. 
           NIMG      = NZ
           NZ        = 1
           ITYPE     = 1
           MAXIM     = 1

        ELSEIF ((ISPG ==   0 .AND. MZ == 1) .OR.
     &          (ISPG ==   1 .AND. MZ == 1)) THEN
C          SPIDER IMAGE OUTPUT         
           NIMG   = 1
           NZ     = MZ

C          OPEN SPIDER IMAGE OUTPUT FILE     
           ITYPE  = 1
           MAXIM  = 0
           CALL OPFILEC(0,.TRUE.,FILOUT,LUNSPI,'U',ITYPE,NX,NY,NZ,
     &                  MAXIM,'SPIDER OUTPUT',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

        ELSEIF (ISPG ==   0 .AND. MZ >  1) THEN
C          SPIDER IMAGE STACK OUTPUT (FOLLOWS MRC-2015 CONVENTION)
           WANTSTACK = .TRUE. 
           NIMG      = MZ
           NZ        = 1
           ITYPE     = 1
           MAXIM     = 1

        ELSEIF (ISPG ==   1 .AND. MZ > 1) THEN
C          VOLUME (MAY BE VOLUME STACK IN PRE-2015)
           NIMG   = 1
           NZ     = MZ

C          OPEN SPIDER VOLUME OUTPUT FILE     
           ITYPE  = 3
           MAXIM  = 0
           CALL OPFILEC(0,.TRUE.,FILOUT,LUNSPI,'U',ITYPE,NX,NY,NZ,
     &                  MAXIM,'SPIDER OUTPUT',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

        ELSEIF (ISPG == 401 .AND. MZ >  1) THEN
C          POST 2015 FORMAT VOLUME STACK 
           WANTSTACK = .TRUE. 
           NIMG      = NZ / MZ
           NZ        = MZ
           ITYPE     = 3
           MAXIM     = 1

        ELSE 
           CALL ERRT(102,'BAD STACK OR VOLUME PARAMETERS',ISPG)
           IRTFLG = 1
           RETURN
        ENDIF


        IF (NIMG > 1) THEN
C          SPIDER IMAGE OR VOLUME STACK OUTPUT

           CALL FILERD(FILOUT,NLET,NULL,'SPIDER OUTPUT STACK~',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          HACK TO USE OPFILES
           LOCAT = INDEX(FILOUT(1:NLET),'@')
           IF (LOCAT == 0) THEN
              NLET              = NLET + 1
              FILOUT(NLET:NLET) = '@'
           ELSEIF (LOCAT .NE. NLET) THEN
              NLET = LOCAT
              FILOUT(NLET+1:) = ' ' 
           ENDIF

           IGO = 1
           CALL RDPRI1S(IGO,NOT_USED,
     &              'FIRST IMAGE NUMBER IN SPIDER STACK',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999


           ALLOCATE(ILIST(NIMG),STAT=IRTFLG)
           IF (IRTFLG > 0) THEN
              CALL ERRT(46,'COPYCCP4; ILIST',NIMG)
              GOTO 9999
           ENDIF

C          MAKE LIST OF STACKED FILE NUMBERS
           DO I = IGO,NIMG
              ILIST(I) = I
           ENDDO

           NSTACKT   = -NIMG    ! USE ILIST FOR STACKED IMG NUMBERS
           FOUROK    = .FALSE.  ! NOT FOURIER
           IMGNUMOUT = IGO
           ASKNAM    = .FALSE.
           !write(6,*) 'filout(1:nlet):', filout(1:nlet)

C          OPEN SPIDER STACK OUTPUT FILE     
           CALL OPFILES(0,LUNSPI,LUNDOCSEL,LUNXM,
     &            ASKNAM,FILOUT,NLET, 'N',
     &            ITYPE,NX,NY,NZ,MAXIM,
     &            FILOUT(1:NLET),
     &            FOUROK, ILIST,NSTACKT, 
     &            NUNUSED,NSTACKOUT, IMGNUMOUT, IRTFLG)
ccccccc           write(6,*) 'nstackout, imgnumout:', nstackout, imgnumout 
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
     &           'FLIP BYTE ORDERING? (Y/N)', NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (ANS(1:1) == 'Y') FLIP = .NOT. FLIP

        ! BOTLEFT IS USUAL MRC BOTTOM IN CURRENT FORMAT !!!!!
        BOTLEFT = .TRUE.
        IF (NCHAR > 2 .AND. (INDEX(ANS(3:NCHAR),'N')) > 0) 
     &     BOTLEFT = .FALSE.
 
        NINDX  = 1

        DO  ! POSSIBLE LOOP OVER MRC STACK -------------------------
         IRECSPI = 0
         IRECINC = 1
         IRECSTK = (NINDX-1) * NY

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
C          32 BIT FOATING POINT CCP4 (MRC) INPUT FILE

           DO IZ = 1,NZ
             DO IY = 1,NY
                IREC    = IRECSTK + (IZ  -1) * NY + IY
c               IPOSMRC = IOFFSET + (IREC-1) * NX * 4 + 1
                      
                IPOSMRC = (IREC-1)    ! KLUDGE FOR INTEGER *8 PRESERVE
                IPOSMRC = IPOSMRC * NX
                IPOSMRC = IPOSMRC * 4
                IPOSMRC = IOFFSET + IPOSMRC + 1

                READ(LUNMRC, POS=IPOSMRC,IOSTAT=IERR) STREAMBUF
                IF (IERR .NE. 0) THEN
                   CALL ERRT(102,'READ ERROR ',IERR)
                ENDIF

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

                !if (irec ==1) write(6,*) ' Val:',iposmrc, streambuf(1)
                !if (iz==1 .and. iy==1)
     &          !        write(6,*) ' Val:',iposmrc, streambuf(1)

C               PUT OUT COMPLETED RECORD
                IRECSPI = IRECSPI + IRECINC
                CALL WRTLIN(LUNSPI,STREAMBUF,NX,IRECSPI)

              ENDDO
           ENDDO

        ELSE
           CALL ERRT(102,'CAN NOT COPY MRC MODE',MODE)
        ENDIF

        IF (.NOT. WANTSTACK) EXIT    ! FINISHED IF NOT A MRC STACK
        IF (NINDX >= NIMG)   EXIT    ! ALREADY FINISHED WITH STACK

C       OPEN NEXT STACKED OUTPUT FILE 
        CALL NEXTFILE(NINDX, ILIST, 
     &                FOUROK,LUNXM,
     &                NIMG,MAXIM,   
     &                LUNSPI,0,
     &                FILOUT,'N',
     &                IMGNUMOUT, IRTFLG) 
        
       IF (IRTFLG == -99) THEN
           CALL ERRT(102,'INSUFFICIENT OUTPUT FILE NAMES',NINDX)
           EXIT         
        ELSEIF (IRTFLG .NE. 0) THEN
           EXIT                      ! ERROR
        ENDIF

       ENDDO   ! END OF STACK LOOP --------------------------

      !iposmrc = iposmrc + nx  
      !write(6,*) ' Last iposmrc:',iposmrc, streambuf(1)
        
9999   CLOSE(LUNSPI)
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

