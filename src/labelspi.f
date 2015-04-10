
C++*********************************************************************
C
C  LABELSPI.F     NEW                              OCT 02 ARDEAN LEITH
C                 NROWOUT BUG                      MAY 03 ARDEAN LEITH
C                 ~9                               APR 04 ARDEAN LEITH
C                 \ REMOVED FOR LINUX              JUN 07 ARDEAN LEITH
C                 BLANK IMAGE SUPPORT              OCT 10 ARDEAN LEITH
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
C   LABELSPI(LUN1,LUN2,LUN3,NSAM,NROW,NSLICE,FMIN1,FMAX1,BOTTEM)
C
C   PURPOSE: PATCH, OR INSERT LABEL IN IMAGE
C
C   PARAMETERS:
C        LUN1,LUN2      INPUT IMAGE & FONT IMAGE 
C        LUN3           OUTPUT IMAGE  
C        NSAM,NROW      DIMENSIONS OF INPUT IMAGE
C        NSAMF,NROWF    DIMENSIONS OF FONT IMAGE
C        FMIN2,FMAX2    MIN & MAX FOR INPUT IMAGE
C
C   VARIABLES
C        FMINF,FMAXF    MIN & MAX FOR FONT IMAGE
C        WIDEF          WIDTH OF LETTER IN FONT IMAGE
C        OFFCON         START OF LETTERS IN FONT IMAGE 
C
C--********************************************************************

      SUBROUTINE LABELSPI(LUN1,LUN2,LUN3,NSAM,NROW,NSLICE,FMIN1,FMAX1,
     &                     BOTTEM)

      INCLUDE 'CMBLOCK.INC'    
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=256)    :: LABELSTR
      REAL, ALLOCATABLE     :: BUF(:,:),BUFFONT(:,:)
      CHARACTER(LEN=1)      :: NULL = CHAR(0)
      CHARACTER(LEN=93)     :: LETS
      CHARACTER(LEN=MAXNAM) :: FONTDIR,FONTNAM
      LOGICAL               :: BLACKFONT,BOTTEM

C                    123456789 123456789 123456789 123456789 1234567
      LETS(1:47)  = '! #$ & ()*+,-./0123456789:<=>?@ABCDEFGHIJKLMNOP'
C                    89 123456789 123456789 123456789 123456789 123
C     LETS(48:93) = 'QRSTUVWXYZ[\]^_ abcdefghijklmnopqrstuvwxyz{|}~'
      LETS(48:93) = 'QRSTUVWXYZ[ ]^_ abcdefghijklmnopqrstuvwxyz{|}~'

CC    FONT WIDTH (DEPENDS ON FONT IN USE!!!
CC    CALL RDPRM2S(WIDEF,OFFCON,NOT_USED,
CC     &           'FONT WIDTH & OFFSET',IRTFLG)
CC    IF (IRTFLG .NE. 0) GOTO 9999

CC    FONT WIDTH (DEPENDS ON FONT IN USE!!!
CC    CALL RDPRM2S(AA,AB,NOT_USED,
CC     &          'OFFSET A & A',IRTFLG)
CC     IF (IRTFLG .NE. 0) GOTO 9999

C      KLUDGE TO GET WORKING FAST, FIX IT SOMETIME
       WIDEF  = 12.3
       OFFCON = 18.0
       AA     = 385.5
       AB     = 754.0

       IXADD  = -1
       IYADD  = -25

       BLACKFONT = (FCHAR(4:4) .EQ. 'B' .OR. FCHAR(5:5) .EQ. 'B')  

C      GET DIR FOR FONT INPUT IMAGE 
       CALL MYGETENV('SPPROC_DIR',FONTDIR,NCHART,
     &              'dir-for-proc-files',IER)
       IF (IER .NE. 0) THEN
         CALL ERRT(101,'NO ENVIRONMENT VARIABLE',NE)
         GOTO 9999
       ENDIF

       IF (BLACKFONT) THEN 
         FONTNAM = FONTDIR(:NCHART) // 'black_font.img' // CHAR(0)
       ELSE 
         FONTNAM = FONTDIR(:NCHART) // 'white_font.img' // CHAR(0) 
       ENDIF

C      OPEN FONT IMAGE, KEEP EXTENSION (~9)
       MAXIM2 = 0
       CALL OPFILEC(0,.FALSE.,FONTNAM,LUN2,'O',IFORM,
     &              NSAMF,NROWF,NSLICEF, MAXIM2,'DUM~9',.FALSE.,IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 9999

       IF (IMAMI .NE. 1)
     &    CALL NORM3(LUN2,NSAMF,NROWF,NSLICEF,FMAX,FMIN,AV)
       FMINF = FMIN
       FMAXF = FMAX

C      FIND OUTPUT IMAGE SIZE
       NROWOUT = NROW
       IF (BOTTEM) NROWOUT = NROW + NROWF

C      ALLOCATE IMAGE & FONT BUFFERS
       ALLOCATE (BUF(NSAM,NROWOUT), BUFFONT(NSAMF,NROWF), STAT=IRTFLG)
       IF (IRTFLG .NE. 0)  THEN
          MWANT = NSAM * NROWOUT + NSAMF * NROWF
          CALL  ERRT(46,'BUF...',MWANT)
          GOTO 9999
       ENDIF

       IF (BOTTEM) THEN
C        MUST CLEAR REST OF OUTPUT IMAGE BUFFER 
         FCLEAR = FMIN1
         IF (BLACKFONT) FCLEAR = FMAX1

         DO IROW = NROW+1,NROW+NROWF
            DO ISAM = 1,NSAM
               BUF(ISAM,IROW) = FCLEAR
            ENDDO
         ENDDO

       ELSEIF (NSLICE .GT. 1) THEN
C        FILL OUTPUT VOLUME WITH INPUT VOLUME
         DO ISLICE = 1,NSLICE
            CALL REDVOL(LUN1,NSAM,NROW,ISLICE,ISLICE,BUF,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999
            CALL WRTVOL(LUN3,NSAM,NROWOUT,ISLICE,ISLICE,BUF,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999
         ENDDO
       ENDIF

C      LOAD INPUT IMAGE INTO BUF
       ISLICE = 1
       CALL REDVOL(LUN1,NSAM,NROW,ISLICE,ISLICE,BUF,IRTFLG)
       IF (IRTFLG .NE. 0)  GOTO 9999

C      LOAD FONT IMAGE INTO BUFFONT
       CALL REDVOL(LUN2,NSAMF,NROWF,1,1,BUFFONT,IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 9999

C      SCALE FONT INTENSITY TO IMAGE RANGE
C      BUF = FVAL * SCALE + FMIN1

       IF (FMAX1 .EQ. FMIN1 .AND. BLACKFONT) THEN
C         BLANK IMAGE WITH BLACK FONT

          SCALE  = 1 
          FMIN1  = -1
          !write(6,91) 'Black: ',fmin1,fmax1,scale,fminf,fmaxf
91        format(' ',A,f5.2,'...',f5.2,f9.2,f5.2,'..',f8.2)
 
       ELSEIF (FMAX1 .EQ. FMIN1) THEN
C         BLANK IMAGE WITH WHITE FONT
          SCALE  = 255.0 / (FMAXF - FMINF)
          !write(6,91) 'White: ',fmin1,fmax1,scale,fminf,fmaxf
       ELSE
C         NON-BLANK IMAGE
          SCALE  = (FMAX1 - FMIN1) / (FMAXF - FMINF)
          !write(6,91) 'Input: ',fmin1,fmax1,scale,fminf,fmaxf
       ENDIF

10     CONTINUE
C      DO NOT UPPERCASE THE INPUT
       IRTFLG = -999
       CALL RDPRMC(LABELSTR,NLET,.TRUE.,'LABEL (<CR> TO END)',
     &            NULL,IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 9999
       IF (NLET .LE. 0) GOTO 20

       IF (BOTTEM) THEN
C        CHECK LABEL LENGTH & CENTER ACROSS WIDTH OF IMAGE
         IXGO = ((NSAM - (NLET * WIDEF)) / 2.0) 
         IF (IXGO .LE. 0) THEN
            NLET = NSAM / WIDEF
            WRITE(NOUT,*) ' LABEL TRUNCATED TO:',NLET,'  CHARACTERS'
            IXGO = ((NSAM - (NLET * WIDEF)) / 2.0) 
         ENDIF
         IYGO = NROW + 1

       ELSE
C        LABELING INSIDE IMAGE, NEED TO GET LOCATION
         IF (NSLICE .GT. 1) THEN
C           VOLUME INPUT  
15          CALL RDPRI3S(IX,IY,IZ,NOT_USED,'LOCATION, X, Y & Z',IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999

            IF (IZ .LT. 1 .OR. IZ .GT. NSLICE) THEN
               CALL ERRT(102,'ILLEGAL SLICE:',IZ)
               GOTO 15

            ELSEIF (IZ .NE. ISLICE) THEN
C              MUST LOAD INPUT SLICE INTO BUFFER
               ISLICE = IZ
               CALL REDVOL(LUN1,NSAM,NROW,ISLICE,ISLICE,BUF,IRTFLG)
            ENDIF
         ELSE
C           2D IMAGE
            CALL RDPRIS(IX,IY,NOT_USED,'LOCATION, X & Y',IRTFLG)
            IZ = 1
         ENDIF
         IF (IRTFLG .NE. 0) GOTO 9999

C        FIND Y STARTING LOCATION IN IMAGE 
         IF (IY .GT. NROW) THEN
            IY = NROW
            WRITE(NOUT,*)' LABEL OFF IMAGE, MOVED UP TO:',NROW
         ENDIF

         IF (IY .LT. 1) THEN
            IY   = 1
            WRITE(NOUT,*)' LABEL OFF IMAGE, MOVED DOWN TO:',IY
         ENDIF

C        CHECK X STARTING LOCATION FOR 1'ST CHAR. IN IMAGE
         IXGO = IX + IXADD

         NT =  (NSAM - IX) / WIDEF
         NLETCAN = MIN(NLET,NT)
         IF (NLETCAN .LT. NLET) THEN
            WRITE(NOUT,*)' LABEL TRUNCATED TO:',NLETCAN,'  CHARACTERS'
         ENDIF
         IYGO = IY + IYADD
      ENDIF
      
      DO I = 1, NLET
         ILET = INDEX(LETS,LABELSTR(I:I))

         IF (ILET .LE. 0) THEN
            WRITE(NOUT,*) ' UNKNOWN CHARACTER: ',LABELSTR(I:I),
     &                    'REPLACED WITH A BLANK.'
            ILET = 2
         ENDIF

C        STARTING LOCATION OF ILET IN FONT IMAGE 
         IF (ILET .LE. 30) THEN
            IXF = (ILET - 1) * WIDEF + OFFCON
         ELSEIF (ILET .GT. 30 .AND. ILET .LE. 60) THEN
            IXF = (ILET - 31) * WIDEF + AA
         ELSE
            IXF = (ILET - 61) * WIDEF + AB
         ENDIF

         IYF = 1

C        STARTING LOCATION FOR LETTER IN IMAGE 
         IXI = IXGO + (I - 1) * WIDEF

         DO IY = 0,NROWF-1
            DO IX = 0,WIDEF-1
               FVAL = BUFFONT(IXF+IX,IYF+IY)

               IF (BOTTEM) THEN
                  BUF(IXI+IX,IY+IYGO) = FVAL * SCALE + FMIN1

               ELSE
                  IYI = IY + IYGO
                  IF (IYI .GT. 0 .AND. IYI .LE. NROW) THEN
C                    INSIDE IMAGE
 
                     IF (BLACKFONT .AND. FVAL .LE. 254) THEN
                        BUF(IXI+IX,IY+IYGO) = FVAL * SCALE + FMIN1
 
                     ELSEIF (.NOT. BLACKFONT .AND. FVAL .GT. 0) THEN
                        BUF(IXI+IX,IY+IYGO) = FVAL * SCALE + FMIN1
                     ENDIF 
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      IF (.NOT. BOTTEM) GOTO 10

20    CONTINUE
C     DUMP IMAGE BUFFER TO FILE
      !write(6,*)'  NSAM,NROWOUT,ISLICE:', NSAM,NROWOUT,ISLICE
      CALL WRTVOL(LUN3,NSAM,NROWOUT,ISLICE,ISLICE,BUF,IRTFLG)
      IF (IRTFLG .NE. 0)  GOTO 9999

9999  IF (ALLOCATED(BUF))     DEALLOCATE(BUF)
      IF (ALLOCATED(BUFFONT)) DEALLOCATE(BUFFONT)

      END

