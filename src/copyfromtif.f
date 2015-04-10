C ++********************************************************************
C          
C COPYFROMTIF      NEW                          MAR 2014 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C COPYFROMTIF(LUNIN,LUNSPI,MAXMEM) 
C 
C PURPOSE: CONVERT TIFF IMAGE TO SPIDER IMAGE
C
C PARAMETER:               
C       COMMON TIFF TAGS 
C	IMAGEWIDTH		     256
C	IMAGELENGTH		     257
C	BITSPERSAMPLE		     258
C	COMPRESSION		     259
C	PHOTOMETRIC INTERPRETATION   262
C	STRIPOFFSETS		     273
C	ROWSPERSTRIP		     278
C	STRIPBYTECOUNTS		     279
C	RESOLUTIONUNIT		     296
C
C WARNING: IS DEPENDENT ON: byteswapio COMPILATION IF USE ANYTHING OTHER
C          THAN READ I1 OR CHAR INPUT!!!
C
C NOTES:   DEBUG    = (FCHAR(12:12) == 'D') 
C          FLIPDATA = (FCHAR(12:12) == 'S') 
C
C          TRY: od -t u1 -N200 RSp614.tif
C
C *********************************************************************

#ifdef osx86

        SUBROUTINE COPYFROMTIF(LUNIN,LUNSPI,IRTFLG)

        IMPLICIT NONE
        INTEGER          :: LUNIN,LUNSPI,IRTFLG
        INTEGER          :: NE

        CALL ERRT(101,'CAN NOT COMPILE ON PRE-FORTRAN 2008 COMPILER',NE)

        END

#else
        SUBROUTINE COPYFROMTIF(LUNIN,LUNSPI,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

 
        INTEGER                       :: LUNIN,LUNSPI,IRTFLG

        CHARACTER(LEN=MAXNAM)         :: FILOLD,FILNEW
        CHARACTER(LEN=1)              :: NULL = CHAR(0)
        CHARACTER(LEN=1)              :: CFLIP
        CHARACTER(LEN=2)              :: C2,CMAGIC

        INTEGER, PARAMETER            :: NLISTMAX = 128000 
 
        INTEGER(KIND=1),  ALLOCATABLE :: IFD1BUF(:)
        INTEGER(KIND=1),  ALLOCATABLE :: I1BUF(:)
        INTEGER(KIND=4),  ALLOCATABLE :: IHEADSTRIPS(:)
        INTEGER(KIND=4),  ALLOCATABLE :: IBYTESTRIPS(:)
        REAL,             ALLOCATABLE :: BUF(:)

        INTEGER(KIND=1),  ALLOCATABLE :: I1VALS(:)
        INTEGER(KIND=4),  ALLOCATABLE :: I4VALS(:)

        CHARACTER(LEN=1), ALLOCATABLE :: CVAL(:)

        CHARACTER(LEN=24)             :: CDATE

        INTEGER(KIND=1)               :: I1VAL
        INTEGER(KIND=2)               :: I2VAL
        INTEGER                       :: I4VAL
        INTEGER(KIND=8)               :: IPOS

        INTEGER(KIND=1)               :: I1(2), I1_TMP

        LOGICAL               :: FLIP,FOLD,DEBUG,FLIPDATA
        LOGICAL               :: GOTXRES,GOTYRES,GOTDATE,GOTRESU

        INTEGER               :: NX,NY,NZ,IHEAD
        INTEGER               :: IX,IY,IVAL,IERR,IRECSPI,IRECINC
        INTEGER               :: MAGIC,IFDHEADER, NTAGS, NLEN, NUNUSED
        INTEGER               :: IDIR,   IOFF, ITAG, ITYPE, ICOUNT, IXT
        INTEGER               :: NBITS,  ICOMP,  NSTRIPS, NYPERSTRIP
        INTEGER               :: IXRES1, IXRES2, IYRES1,  IYRES2
        INTEGER               :: IRESU,  MAXIM,  IBYTES,  NBYTES,ISTRIP
        REAL                  :: XRES,YRES

        INTEGER               :: IFD1_LEN,I1VALS_LEN,NCHAR,ILOC,IBUF_LEN
        INTEGER               :: I4VALS_LEN,I1BUF_LEN,ICVAL_LEN

C       MAX NUMBER OF TAGS ALLOWED HERE IS: 1000!!!
        INTEGER, PARAMETER    :: IFD1_TAGS = 1000

C       OPEN TIF FILE FOR STREAM ACCESS
        CALL OPSTREAMFILE(.TRUE.,FILOLD,NULL,LUNIN,
     &                    'UNFORMATTED','O', 'TIFF INPUT',
     &                    .TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       MAX NUMBER OF TAGS ALLOWED HERE IS: 1000
        IFD1_LEN = IFD1_TAGS * 12
        ALLOCATE(IFD1BUF(IFD1_LEN),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'COPYFROMTIF; IFD1BUF',IFD1_LEN)
            GOTO 9999
        ENDIF

C       READ FIRST 8 BYTES FROM FILE
        READ(LUNIN, POS=1,IOSTAT=IERR) IFD1BUF(1:8)
        IF (IERR .NE. 0) THEN
            CALL ERRT(102,'OVERALL HEADER I/O',IERR) 
            GOTO 9999
        ENDIF

        WRITE(NOUT,*) ' '

        C2 = TRANSFER(IFD1BUF(1:2),C2)
        WRITE(NOUT,90) C2
90      FORMAT('  ENDEDNESS:',T32,A7) 

        FLIP = (C2 .NE. 'II')  ! OTHER IS: 'MM'
        WRITE(NOUT,91) FLIP
91      FORMAT('  FLIP BYTES:',T32,L7) 

        !write(6,*) '  IFD1BUF : ',ifd1buf(1:8)
        IF (FLIP) CALL FLIP_I2(2,IFD1BUF)
        !write(6,*) '  IFD1BUF : ',ifd1buf(1:8)

        MAGIC = TRANSFER(IFD1BUF(3:4),I2VAL)
        !write(6,*) ' magic i1, flip : ',ifd1buf(3:4),flip

        WRITE(NOUT,92) MAGIC
92      FORMAT('  MAGIC NUMBER (SHOULD = 42):',T32,I7) 

        IF (MAGIC .NE. 42) THEN
           CALL ERRT(102,
     &       'INVALID TIFF FILE MAGIC NUMBER:',MAGIC)
           GOTO 9999
        ENDIF

        IF (FLIP) CALL FLIP_I4(1,IFD1BUF(5))
        !write(6,*) '  ifd1buf : ',ifd1buf(1:8)

        IFDHEADER = TRANSFER(IFD1BUF(5:8),I4VAL)

        WRITE(NOUT,94) IFDHEADER
94      FORMAT('  TIFF IFD HEADER: ',T30,I9) 

C       READ FIRST 2 BYTE INTEGER FROM IFD = NTAGS
        READ(LUNIN, POS=(IFDHEADER+1),IOSTAT=IERR) IFD1BUF(1:8)
        IF (IERR .NE. 0) THEN
            CALL ERRT(102,'IFD1BUF(1:8) I/O',IERR) 
            GOTO 9999
        ENDIF

        IF (FLIP) CALL FLIP_I2(1,IFD1BUF)
        NTAGS = TRANSFER(IFD1BUF(1:2),I2VAL)
        WRITE(NOUT,95) NTAGS
95      FORMAT('  NUMBER OF TAGS: ',T32,I7,/) 

C       READ WHOLE FIRST IFD TAG SET, 12 BYTES PER TAG
        NLEN = NTAGS * 12 

        IF (NLEN > IFD1_LEN) THEN
            CALL ERRT(102,'COPYFROMTIF; IFD1_LEN MUST BE >',NLEN)
            GOTO 9999
        ENDIF

        READ(LUNIN, POS=(IFDHEADER+1),IOSTAT=IERR) IFD1BUF(1:NLEN)
        IF (IERR .NE. 0) THEN
            CALL ERRT(102,'IFD1BUF TAGS I/O',IERR) 
            GOTO 9999
        ENDIF

        NUNUSED  = 0
        GOTXRES  = .FALSE.
        GOTYRES  = .FALSE.
        GOTDATE  = .FALSE.
        GOTRESU  = .FALSE.
        DEBUG    = (FCHAR(12:12) == 'D') 
        FLIPDATA = (FCHAR(12:12) == 'S') 

        I1VALS_LEN = 0
        ICVAL_LEN  = 0

        DO IDIR = 1,NTAGS
           IOFF   = (IDIR - 1) * 12 + 3  ! 3 IS OFFSET FOR:NTAGS

C          FIND TAG
           ITAG   = TRANSFER(IFD1BUF(IOFF:IOFF+1), I2VAL)
           IF (FLIP) CALL FLIP_I2(1,ITAG)
           IF (ITAG < 0) ITAG = 65536 + IVAL
           !WRITE(6,'(2X,6(A,1X,I0,2X))') 'tag:',itag

C          FIND TYPE: 1=BYTE, 2=ASCII, 3=SHORT, 4=LONG, 5=2*LONG=NUM/DEN
           ITYPE  = TRANSFER(IFD1BUF(IOFF+2:IOFF+3), I2VAL)
           IF (FLIP) CALL FLIP_I2(1,ITYPE)

C          FIND COUNT (4BYTES)
           ICOUNT = TRANSFER(IFD1BUF(IOFF+4:IOFF+7), I4VAL)
           !write(6,*) ' IFD1BUF(IOFF+4:',IFD1BUF(IOFF+4:IOFF+7),icount
           IF (FLIP) CALL FLIP_I4(1,ICOUNT)
 
C          FIND VALUE (2BYTES OR 4BYTES)
           IF (ITYPE == 4 ) THEN
C             4 BYTES VALUE

              IVAL = TRANSFER(IFD1BUF(IOFF+8:IOFF+11),I4VAL)
              IF (FLIP) CALL FLIP_I4(1,IVAL)

           ELSEIF (ITYPE == 3 ) THEN
C             2 BYTES VALUE

              I2VAL = TRANSFER(IFD1BUF(IOFF+8:IOFF+9),I2VAL)

              IF (FLIP) CALL FLIP_I2(1,I2VAL)
              IVAL = I2VAL

           ENDIF          


           IF ((ITYPE == 4 .AND. ICOUNT > 1) .OR. ITYPE == 5 ) THEN

C             READ ICOUNT INTEGERS FROM POINTER LOCATION
              IF (ITYPE == 5 ) ICOUNT = 2

C             SET LENGTH FOR VALUES BUFFER
              IF ( (ICOUNT * 4) > I1VALS_LEN) THEN
                IF (ALLOCATED(I1VALS)) DEALLOCATE(I1VALS)
                IF (ALLOCATED(I4VALS)) DEALLOCATE(I4VALS)
                I1VALS_LEN = ICOUNT * 4
                I4VALS_LEN = ICOUNT 
                ALLOCATE(I1VALS(I1VALS_LEN),
     &                   I4VALS(I4VALS_LEN),STAT=IRTFLG)
                
                IF (IRTFLG .NE. 0) THEN
                   CALL ERRT(46,'COPYFROMTIF; I1VALS...',I1VALS_LEN)
                   GOTO 9999
                 ENDIF
              ENDIF

              IPOS = IVAL + 1 
              READ(LUNIN, POS=IPOS, IOSTAT=IERR) I1VALS(1:ICOUNT*4)
 
              IF (IERR .NE. 0) THEN
                   WRITE(NOUT,'(2X,6(A,1X,I0,2X))')
     &                 'POS:',IPOS, 'CNT:',ICOUNT, 'ERR:',IERR
c                   CALL ERRT(102,'I/O',IERR) 
c                   GOTO 9999
              ENDIF
              IF (FLIP) CALL FLIP_I4(ICOUNT,I1VALS)
              DO IX = 1,ICOUNT
                 IXT        =  (IX-1) *4 + 1  
                 I4VAL      = TRANSFER(I1VALS(IXT:ixt+3),I4VAL)
                 I4VALS(IX) = I4VAL
              ENDDO  

              IF (DEBUG)
     &            WRITE(NOUT,*)' POINTED VALUES:',I4VALS(1:4), IPOS

           ELSEIF  (ITYPE == 2)  THEN

C             ASCII VALUES, SET LENGTH FOR CVAL BUFFER
              IF ( ICOUNT  > ICVAL_LEN) THEN
                 IF (ALLOCATED(CVAL)) DEALLOCATE(CVAL)
                 ICVAL_LEN = ICOUNT
                 ALLOCATE(CVAL(1)(ICVAL_LEN), STAT=IRTFLG)
                
                 IF (IRTFLG .NE. 0) THEN
                   CALL ERRT(46,'COPYFROMTIF; CVAL...',ICVAL_LEN)
                   GOTO 9999
                 ENDIF
              ENDIF

              IPOS = IVAL + 1 
              READ(LUNIN, POS=IPOS,IOSTAT=IERR) CVAL(1:ICOUNT)
              IF (DEBUG) WRITE(NOUT,*)' Text: ',CVAL(1:ICOUNT)

           ENDIF

           WRITE(NOUT,96)ITAG,ITYPE,ICOUNT,IVAL
96         FORMAT('  TAG: ',I6,
     &            '  TYPE:',I3,
     &            '  LENGTH: ',I6,
     &            '   VALUE:',I10, I8) 

           IF (ITAG == 256) THEN
              NX = IVAL
c             WRITE(NOUT,'(2X,A,1X,I0)') 'NX:',NX

           ELSEIF (ITAG == 257) THEN
              NY = IVAL
c             WRITE(NOUT,'(2X,A,1X,I0)') 'NY:',NY

           ELSEIF (ITAG == 258) THEN
              NBITS = IVAL
c             WRITE(NOUT,'(2X,A,1X,I0)') 'NBITS',NY

           ELSEIF (ITAG == 259) THEN
              ICOMP = IVAL
              IF (ICOMP .NE. 1) THEN
                 CALL ERRT(102,'COMPRESSION NOT SUPPORTED',ICOMP)
                 GOTO 9999
              ENDIF

           ELSEIF (ITAG == 262) THEN   !PHOTOMETRIC INTERPRETATION 
              ITYPE = IVAL 	
              IF (ITYPE .NE. 1) THEN
                 CALL ERRT(102,'NOT GREY-SCALE DATA',ITYPE)
                 GOTO 9999
              ENDIF

           ELSEIF (ITAG == 273) THEN             ! STRIP HEADERS
              IF (IVAL == 0) ICOUNT = 1
              NSTRIPS = ICOUNT

              ALLOCATE(IHEADSTRIPS(NSTRIPS),STAT=IRTFLG)
              IF (IRTFLG .NE. 0) THEN
                 CALL ERRT(46,'COPYFROMTIF; IHEADSTRIPS',NSTRIPS)
                 GOTO 9999
              ENDIF

             IF (IVAL == 0) THEN
                 IHEADSTRIPS(1) = 8
                 IF (DEBUG) WRITE(NOUT,'(2X,A,1X,I0)')
     &                'NSTRIPS:',NSTRIPS
                 IF (DEBUG) WRITE(NOUT,'(2X,A,1X,I0)')
     &                'STRIP HEADER DEFAULT: ',IHEADSTRIPS(1)

              ELSEIF (ICOUNT == 1) THEN
                 IHEADSTRIPS(1) = IVAL
                 IF (DEBUG)WRITE(NOUT,'(2X,A,1X,I0)')'NSTRIPS:',NSTRIPS
                 IF (DEBUG)WRITE(NOUT,'(2X,A,1X,9(I0,1X))') 
     &                'STRIP HEADER:',IHEADSTRIPS(1)
              ELSE
                 IHEADSTRIPS(1:ICOUNT) = I4VALS  ! ARRAY OP
                 IF (DEBUG)WRITE(NOUT,'(2X,A,1X,I0)')'NSTRIPS:',NSTRIPS
                 IF (DEBUG)WRITE(NOUT,'(2X,A,1X,9(I0,1X))') 
     &                'STRIP HEADERS:',IHEADSTRIPS(1:4)
              ENDIF

           ELSEIF (ITAG == 277) THEN
              IF (DEBUG) WRITE(NOUT,'(2X,A,1X,I0)')
     &                   'SAMPLES / PIXEL:',IVAL
              IF (IVAL .NE. 1) THEN
                 CALL ERRT(102,
     &               'IMAGE SAMPLES/PIXEL NOT SUPPORTED',IVAL)
                 GOTO 9999
              ENDIF

           ELSEIF (ITAG == 278) THEN
              NYPERSTRIP = IVAL
              IF (DEBUG) WRITE(NOUT,'(2X,A,1X,I0)') 
     &                   'ROWS PER STRIP:',IVAL

           ELSEIF (ITAG == 279) THEN             ! STRIP BYTES

              ALLOCATE(IBYTESTRIPS(NSTRIPS),STAT=IRTFLG)
              IF (IRTFLG .NE. 0) THEN
                 CALL ERRT(46,'COPYFROMTIF; IBYTESTRIPS',NSTRIPS)
                 GOTO 9999
              ENDIF

              IBYTESTRIPS(1) = NBITS / 8 * NX *NYPERSTRIP

              IF (IVAL == 0) THEN
                 IF (DEBUG) WRITE(NOUT,'(2X,A,1X,I0)')
     &                'STRIP BYTES DEFAULT: ',IBYTESTRIPS(1)

              ELSEIF (ICOUNT == 1) THEN
                 IBYTESTRIPS(1) = IVAL
                 IF (DEBUG) WRITE(NOUT,'(2X,A,1X,9(I0,1X))') 
     &                'STRIP BYTES:',IBYTESTRIPS(1)
              ELSE
                 IBYTESTRIPS(1:ICOUNT) = I4VALS  ! ARRAY OP
                 IF (DEBUG) WRITE(NOUT,'(2X,A,1X,9(I0,1X))') 
     &                'STRIP BYTES:',IBYTESTRIPS(1:4)
              ENDIF

           ELSEIF (ITAG == 282) THEN           ! X RESOLUTION
              IXRES1  = I4VALS(1)
              IXRES2  = I4VALS(2)
              XRES    = FLOAT(IXRES1) / FLOAT(IXRES2)
              GOTXRES = .TRUE.

              IF (DEBUG) WRITE(NOUT,'(2X,A,1X,I0,1X,I0,2x,F12.2)')
     &                   'X RESOLUTIONS:', IXRES1,IXRES2,XRES

           ELSEIF (ITAG == 283) THEN           ! Y RESOLUTION
              IYRES1  = I4VALS(1)
              IYRES2  = I4VALS(2)
              YRES    = FLOAT(IYRES1) / FLOAT(IYRES2)
              GOTYRES = .TRUE.

              IF (DEBUG) WRITE(NOUT,'(2X,A,1X,I0,1X,I0,2x,F12.2)')
     &                   'Y RESOLUTIONS:', IYRES1,IYRES2,YRES

           ELSEIF (ITAG == 284) THEN
c             IF (DEBUG) WRITE(NOUT,'(2X,A,1X,I0)') 
c                             'PLANAR CONFIG.:',IVAL

           ELSEIF (ITAG == 296) THEN
              IF (DEBUG) WRITE(NOUT,'(2X,A,1X,I0)') 
     &                         'RESOLUTION UNIT:',IVAL
              IRESU = IVAL
              GOTRESU = .TRUE.

           ELSEIF (ITAG == 306) THEN
             CDATE = CVAL(1)(1:22)
c            IF (DEBUG) WRITE(NOUT,'(2X,A,1X,200A)')'DATE & TIME:',
c     &                                    CVAL(1)(1:ICOUNT)
             GOTDATE = .TRUE.

           ELSE
c             WRITE(NOUT,'(2X,A,T32,I7)') 'UNUSED TAG:',ITAG
              NUNUSED = NUNUSED + 1
           ENDIF
        ENDDO

        IF (DEBUG) WRITE(NOUT,'(2X,A,1X,I0)') 'UNUSED TAGS:',NUNUSED

        WRITE(NOUT,*) ' '
        WRITE(NOUT,'(2X,A,1X,I0)') 'NUMBER OF STRIPS:',NSTRIPS
        IF (NSTRIPS == 1) WRITE(NOUT,'(2X,A,1X,I0)') 
     &             'BYTES/STRIP:',IBYTESTRIPS(1) 
        IF (GOTXRES)WRITE(NOUT,'(2X,A,1X,I0)') 
     &             'RESOLUTION UNIT:',IRESU
        IF (GOTXRES)WRITE(NOUT,'(2X,A,1X,F12.2)') 'X RESOLUTION:',XRES
        IF (GOTXRES)WRITE(NOUT,'(2X,A,1X,F12.2)') 'Y RESOLUTION:',YRES
        IF (GOTDATE)WRITE(NOUT,'(2X,A,1X,A)')'DATE & TIME:',CDATE(1:19)
        WRITE(NOUT,*)' ' 
        WRITE(NOUT,97) NX,NY,NBITS,IHEADSTRIPS(1)
97      FORMAT('  NX: ',I0,
     &         '  NY: ',I0,
     &         '    BITS: ',I0, 
     &         '    HEADER BITS: ',I0) 
        WRITE(NOUT,*) ' '

C       OPEN NEW OUTPUT SPIDER IMAGE FILE WITH SPECIFIED SIZE
        IFORM  = 1
        NZ     = 1
        MAXIM  = 0 
        CALL OPFILEC(0,.TRUE.,FILNEW,LUNSPI,'U',IFORM,NX,NY,NZ,
     &              MAXIM,'SPIDER OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        FOLD   = .TRUE.

        IF (DEBUG) THEN
           WRITE(NOUT,'(2X,A,T32,I10)') 'NSTRIPS:',NSTRIPS
           WRITE(NOUT,'(2X,A,T32,I10)') 'NYPERSTRIP:',NYPERSTRIP
           WRITE(NOUT,'(2X,A,T32,I10)') 'NBYTES:',NBYTES
           WRITE(NOUT,'(2X,A,T32,L10)') 'FLIP:',FLIP
           WRITE(NOUT,'(2X,A,T32,L10)') 'FLIPDATA:',FLIPDATA
           WRITE(NOUT,'(2X,A,T32,L10)') 'FOLD:',FOLD
        ENDIF

        CALL FLUSHRESULTS ! -------------------------------------

        IRECSPI = 0
        IRECINC = 1

C       FIND MAX NUMBER OF BYTES IN READ BUFFER
        I1BUF_LEN = MAXVAL(IBYTESTRIPS(1:NSTRIPS)) 

        IF (NBITS == 8) THEN
C          8 BIT INTEGER TIFF INPUT FILE

           IBUF_LEN  = I1BUF_LEN   
           ALLOCATE(I1BUF(I1BUF_LEN), BUF(IBUF_LEN),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'COPYFROMTIF; I1BUF..',I1BUF_LEN+IBUF_LEN)
              GOTO 9999
           ENDIF

C          8 BIT INTEGER TIFF INPUT FILE
           DO ISTRIP = 1,NSTRIPS
C             NEW IMAGE STRIP

              IPOS   = IHEADSTRIPS(ISTRIP) + 1
              IBYTES = IBYTESTRIPS(ISTRIP) 

              READ(LUNIN, POS=IPOS,IOSTAT=IERR) I1BUF(1:IBYTES)
              IF (IERR .NE. 0) THEN
                 WRITE(3,'(2X,6(A,I0,2X))')
     &               'POS: ',IPOS,'BYTES: ',IBYTES,' ERR: ',IERR
                 ! CALL ERRT(102,'I/O',IERR) 
                 RETURN
              ENDIF
              !write(6,*) ' Got pos, ibites:',ipos,ibytes

              IF (FOLD) THEN

                DO IX = 1,IBYTES
                   I1VAL = I1BUF(IX)
                   IF (I1VAL < 0) THEN
                      BUF(IX) = 256 + I1VAL
                   ELSE
                      BUF(IX) = I1VAL
                   ENDIF
                   !if (istrip == 4 .and. I1VAL <0 .and. ix>500.and.ix<600 ) 
                   !write(nout,*) 'ix,i1val:',I1VAL,buf(ix)
               ENDDO 
  
             ELSE
                DO IX = 1,IBYTES
                   BUF(IX) = I1BUF(IX)
                   !write(nout,*) 'ix,i1val:',ix,I1BUF(IX:IX),buf(ix)
               ENDDO 
             ENDIF

             ILOC = 1 - NX
             DO IY = 1,NYPERSTRIP
C               NEW ROW WITHIN STRIP

                ILOC    = ILOC + NX
                IRECSPI = IRECSPI + IRECINC
                CALL WRTLIN(LUNSPI,BUF(ILOC),NX,IRECSPI)
                
              ENDDO   ! END OF: DO IY     = 1,NYPERSTRIP
           ENDDO      ! END OF: DO ISTRIP = 1,NSTRIPS

       ELSEIF (NBITS == 16) THEN
C          16 BIT INTEGER TIFF INPUT FILE

           IBUF_LEN  = I1BUF_LEN / 2 
           ALLOCATE(I1BUF(I1BUF_LEN), BUF(IBUF_LEN),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'COPYFROMTIF; I1BUF..',I1BUF_LEN+IBUF_LEN)
              GOTO 9999
           ENDIF

           DO ISTRIP = 1,NSTRIPS
C             NEW IMAGE STRIP

              IPOS   = IHEADSTRIPS(ISTRIP) + 1
              IBYTES = IBYTESTRIPS(ISTRIP) 

              READ(LUNIN, POS=IPOS,IOSTAT=IERR) I1BUF(1:IBYTES)
              IF (IERR .NE. 0) THEN
                 WRITE(3,'(2X,6(A,I0,2X))')
     &               'POS: ',IPOS,'BYTES: ',IBYTES,' ERR: ',IERR
                 ! CALL ERRT(102,'I/O',IERR) 
                 RETURN
              ENDIF
             !write(6,*) 'got pos,ibites:',ipos,ibytes

              IF (FLIPDATA) THEN
C               INVERT BYTE ORDER

                DO IX = 1,IBUF_LEN
                   IXT   =  IX * 2 - 1

                   I1(1:2) = TRANSFER(I1BUF(IXT:IXT+1),I2VAL)
                   I1_TMP  = I1(1)
                   I1(1)   = I1(2)
                   I1(2)   = I1_TMP
                   I4VAL   = TRANSFER(I1(1:2), I2VAL)

                   IF (I4VAL < 0) THEN
                       BUF(IX) = 65536 + I4VAL
                       !write(6,*) 'ix,i4val:',i4val,buf(ix)
                    ELSE
                       BUF(IX) = I4VAL
                   ENDIF
                ENDDO   ! END OF: DO IX = 1,IBUF_LEN

              ELSEIF (FOLD) THEN

                DO IX = 1,IBUF_LEN
                   IXT   =  IX * 2 - 1
                   I2VAL = TRANSFER(I1BUF(IXT:IXT+1),I2VAL)

                   IF (I2VAL < 0) THEN
                       BUF(IX) = 65536 + I2VAL
                       !write(6,*) 'ix,i2val:',i2val,buf(ix)
                    ELSE
                       BUF(IX) = I2VAL
                   ENDIF
                ENDDO   ! END OF: DO IX = 1,IBUF_LEN
  
             ELSE
                DO IX = 1,IBUF_LEN
                   IXT     =  IX * 2 - 1
                   BUF(IX) = TRANSFER(I1BUF(IXT:IXT+1),I2VAL)
                ENDDO  
             ENDIF

             ILOC = 1 - NX
             DO IY = 1,NYPERSTRIP
C               NEW ROW WITHIN STRIP

                ILOC    = ILOC + NX
                IRECSPI = IRECSPI + IRECINC
                CALL WRTLIN(LUNSPI,BUF(ILOC),NX,IRECSPI)
                
              ENDDO   ! END OF: DO IY     = 1,NYPERSTRIP

           ENDDO      ! END OF: DO ISTRIP = 1,NSTRIPS

        ELSE
           CALL ERRT(102,'UNSUPPORTED BIT LENGTH',NBITS)

        ENDIF

9999    CLOSE(LUNIN)
        CLOSE(LUNSPI)

        IF (ALLOCATED(IFD1BUF))     DEALLOCATE(IFD1BUF)
        IF (ALLOCATED(I1BUF))       DEALLOCATE(I1BUF)
        IF (ALLOCATED(IHEADSTRIPS)) DEALLOCATE(IHEADSTRIPS)
        IF (ALLOCATED(IBYTESTRIPS)) DEALLOCATE(IBYTESTRIPS)
        IF (ALLOCATED(I1VALS))      DEALLOCATE(I1VALS)
        IF (ALLOCATED(I4VALS))      DEALLOCATE(I4VALS)
        IF (ALLOCATED(BUF))         DEALLOCATE(BUF)

        END

#endif



C       ------------------------ FLIP_I2 ---------------------

      SUBROUTINE FLIP_I2(N,I2)

      IMPLICIT NONE

      INTEGER(KIND=2) :: I2(N)
      INTEGER(KIND=4) :: N

      INTEGER(KIND=2) :: I2_VAL
      INTEGER(KIND=1) :: I1(2), I1_TMP
      INTEGER(KIND=4) :: I


      DO I = 1, N
         I1 = TRANSFER(I2(I), I1)

         I1_TMP = I1(1)
         I1(1)  = I1(2)
         I1(2)  = I1_TMP
         !write(6,*) ' i,i2(i),i1',i,i2(i), i1

         I2(i) = TRANSFER(I1(1:2), I2_VAL)
         !write(6,*) ' i2(i)', i2(i)

      ENDDO

      END

C     ------------------------ FLIP_I4 ---------------------

      SUBROUTINE FLIP_I4(N,I4)

      IMPLICIT NONE

      INTEGER(KIND=4) :: I4(N)
      INTEGER(KIND=4) :: N

      INTEGER(KIND=4) :: I4_VAL
      INTEGER(KIND=1) :: I1(4), I1_TMP(4)
      INTEGER(KIND=4) :: I,J


      DO I = 1, N
        I1(1:4) = TRANSFER (I4(I), I1)

         
        !write(6,*) ' i,i4(i),i1',  i,i4(i),':', i1

        I1_TMP = I1(1:4)

        DO J = 1, 4
           I1(J) = I1_TMP(5-J)
        END DO

        I4(I) = TRANSFER (I1(1:4), I4_VAL)

        !write(6,*) ' i4(i)', i4(i)

      ENDDO

      END


#ifdef NEVER
              !WRITE(6,'(2X,6(A,1X,I0,2X))') 'tag:',itag,
              !  'type:',itype, 'Cnt:',icount,'Ival:',ival
              !IVAL = TRANSFER(IFD1BUF(IOFF+8:IOFF+9),I4VAL)
              !write(6,*) ' IFD1BUF(IOFF+8:',IFD1BUF(IOFF+8:IOFF+11),ival
              !IF (FLIP) CALL FLIP_I4(2,IFD1BUF(IOFF+8))

              !IVAL = TRANSFER(IFD1BUF(IOFF+8:IOFF+11),I4VAL)
              !write(6,*) ' IFD1BUF(IOFF+8:',IFD1BUF(IOFF+8:IOFF+11),ival
              !write(6,*) ' ival:',IFD1BUF(IOFF+8:IOFF+11),ival
              !IF (FLIP) CALL FLIP_I2(2,IVAL)
              !write(6,*) ' ival:',ival

              WRITE(6,'(2X,6(A,1X,I0,2X))') 'tag:',itag,
     &          'type:',itype, 'Cnt:',icount,'Ival:',ival

#endif                                                                       
