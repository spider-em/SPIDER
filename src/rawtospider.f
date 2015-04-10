
C++*************************************************************************
C
C  RAWTOSPIDER.F   -- CREATED               JUL 95 ARDEAN LEITH
C                     ADDED -32           APRIL 00 ARDEAN LEITH
C                     REWRITTEN           JUNE  01 ARDEAN LEITH
C                     ALTERED             MAR   02 ARDEAN LEITH
C                     RETURNED IOSTAT     JAN   03 ARDEAN LEITH
C                     QUESTION ORDER      FEB   03 ARDEAN LEITH
C                     RDPR PARAMETERS     04/14/05 ARDEAN LEITH
C                     -33 CORRECTED       07/30/07 ARDEAN LEITH
C                     64 ADDED            07/30/07 ARDEAN LEITH
C                     NPIX8               12/18/08 ARDEAN LEITH
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
C  RAWTOSPIDER(LUNOLD,LUNNEW,IRTFLG)
C
C  PURPOSE:  CONVERTS "RAW" INTEGER TO SPIDER FORMAT
C  copyemi,copymrc,copyfromnik
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

       SUBROUTINE RAWTOSPIDER(LUNOLD,LUNNEW,IRTFLG)

       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC'

       COMMON           BUFRAW(NBUFSIZ)

       CHARACTER(LEN=MAXNAM) :: FILOLD,FILNEW,PROMPT

       CHARACTER(LEN=1)      :: CDUM,ANS,NULL
       LOGICAL               :: FOLD,FLIP
       LOGICAL               :: BIGENDARCH,BIGENDED
       INTEGER * 8           :: NPIX8
 
       NULL = CHAR(0)

C      GET FILENAME FOR EXISTING RAW IMAGE FILE
       CALL FILERD(FILOLD,NLET,DATEXC,'EXISTING RAW~9',IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

C      GET INPUT MODE NOW
       MODE = 8
       CALL RDPRI1S(MODE,NOT_USED,
     &       'BITS / PIXEL IN INPUT IMAGE (8, 16 OR 32)',IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

       MODEA = ABS(MODE)
       IF (MODEA .EQ. 33) MODEA = 32
       IF (MODEA .EQ. 65) MODEA = 64

       IF (MODEA .NE. 8  .AND. 
     &     MODEA .NE. 16 .AND. 
     &     MODEA .NE. 32 .AND.
     &     MODEA .NE. 64) THEN
          CALL ERRT(100,'MUST BE (8,16,32,-32,33,64,-64, OR 65)!',NE)
          RETURN
       ENDIF

        NSLICE = 1
        NSAM   = -1
        CALL RDPRI3S(NSAM,NROW,NSLICE,NOT_USED,
     &               'COLUMNS, ROWS, & SLICES',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

22     IOFFSET = 0
       CALL RDPRI1S(IOFFSET,NOT_USED,
     &             'HEADER BYTES TO BE SKIPPED',IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 9999

       IF ( MODEA .EQ. 16 .OR. MODEA .EQ. 64) THEN
23         ISIGB   = 1
           CALL RDPRI1S(ISIGB,NOT_USED,
     &        'MOST SIGNIFICANT BYTE (1 OR 2)',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
           IF (ISIGB .NE. 1 .AND. ISIGB .NE. 2) THEN
              CALL ERRT(16,'RAWTOSPIDER',NE)
              GOTO 23
           ENDIF

C          GET CURRENT ARCHITECTURE ENDED-NESS
           BIGENDARCH = BIGENDED(0)

           FLIP = .TRUE.
           IF ((ISIGB .EQ. 1 .AND. BIGENDARCH) .OR. 
     &         (ISIGB .EQ. 2 .AND. .NOT. BIGENDARCH)) FLIP = .FALSE. 

           CALL RDPRMC(ANS,NCHAR,.TRUE.,'FOLD NEGATIVES? (N/Y)',
     &                 CDUM,IRTFLG)
           IF (IRTFLG .EQ. -1) GOTO 23
           FOLD = (ANS .EQ. 'Y')
       ENDIF

C      CODE IS INVOLVED IF HEADER IS NOT SAME WORD LENGTH AS DATA

C      OPEN NEW OUTPUT SPIDER IMAGE FILE WITH SPECIFIED SIZE
       IFORM  = 1
       IF (NSLICE .GT. 1) IFORM = 3
       MAXIM = 0 
       CALL OPFILEC(0,.TRUE.,FILNEW,LUNNEW,'U',IFORM,NSAM,NROW,NSLICE,
     &              MAXIM,'NEW SPIDER IMAGE',.FALSE.,IRTFLG)
       IF (IRTFLG .EQ. -1) GOTO 22
       IF (IRTFLG .NE. 0) GOTO 9999

       NPIX8 = NSAM * NROW 
       NPIX8 = NPIX8 * NSLICE

       IOFFMOD2 = MOD(IOFFSET,2)
       IOFFMOD4 = MOD(IOFFSET,4)

C      OPEN RAW FILE AS DIRECT ACCESS, UNFORMATTED REC. 
       LENOPEN = NSAM * (MODEA / 8)
       IF (MODEA .EQ. 33) LENOPEN = NSAM * (32 / 8)
       CALL OPAUXFILE(.FALSE.,FILOLD,CHAR(0),LUNOLD,LENOPEN,'O',
     &                   ' ',.TRUE.,IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 9999

       IF (MODE .EQ. 8) THEN
C         FOR 8 BIT RAW INTEGER INPUT 
          CALL RAW8TOSPI(LUNOLD,LUNNEW,NSAM,NPIX8,IOFFSET,.TRUE.,
     &                   LENOPEN,BUFRAW,IRTFLG)
       
       ELSEIF (MODE .EQ. 16 .AND. IOFFMOD2 .EQ. 0) THEN
C         16 BIT RAW INTEGER WITH HEADER LENGTH DIVISABLE BY TWO
          FOLD = .NOT. FOLD
          CALL RAW16TOSPI(LUNOLD,LUNNEW,NSAM,NPIX8,IOFFSET,FLIP,
     &                   FOLD,LENOPEN,BUFRAW,IRTFLG)

       ELSEIF (MODEA .EQ. 32) THEN
C         COPY/AND ALTER BYTE ORDER OF 32 BIT RAW FLOATING POINT IMAGES

          NFLIP = 0
          IF (MODE .EQ. -32) NFLIP = -1  ! INVERTS BYTES
          IF (MODE .EQ. -33) NFLIP =  1  ! FLIPS BYTES WITHIN WORDS
          IF (MODE .EQ.  33) NFLIP =  2  ! FLIPS BYTES & WORDS

          CALL RAW32TOSPI(LUNOLD,LUNNEW,NSAM,NPIX8,IOFFSET,
     &                    NFLIP,LENOPEN,BUFRAW,IRTFLG)

       ELSEIF (MODEA .EQ. 64 .OR. MODEA .EQ. 65) THEN
C         COPY/AND ALTER BYTE ORDER OF 32 BIT RAW INTEGER IMAGES

          NFLIP = 0
          IF (MODE .EQ. -64)  NFLIP = -1  ! INVERTS BYTES
          IF (MODE .EQ. -65)  NFLIP =  1  ! FLIPS BYTES WITHIN WORDS
          IF (MODE .EQ.  65)  NFLIP =  2  ! FLIPS BYTES & WORDS

          FOLD = .NOT. FOLD
          CALL RAW32INTTOSPI(LUNOLD,LUNNEW,NSAM,NPIX8,IOFFSET,
     &                    NFLIP,FOLD,LENOPEN,IRTFLG)

       ELSE
C         OTHER RAW IMAGE FORMATS OR CONVERSIONS (SLOW)

C         OPEN RAW FILE AS DIRECT ACCESS, UNFORMATTED REC. = 1 BYTE
          CLOSE(LUNOLD)
          LENOPEN = 1
          CALL OPAUXFILE(.FALSE.,FILOLD,CHAR(0),LUNOLD,LENOPEN,'O',
     &                   ' ',.TRUE.,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9999

          CALL RAWTOSPI(LUNOLD,LUNNEW,NSAM,NPIX8,IOFFSET,
     &                 MODE,FOLD,FLIP,ISIGB,IRTFLG)
       ENDIF


9999   CLOSE(LUNOLD)
       CLOSE(LUNNEW)

       RETURN
       END
     


C ------------------------- RAW8TOSPI -------------------------------

       SUBROUTINE RAW8TOSPI(LUNRAW,LUNSPI,NSAM,NPIX8,IOFFSET,FOLD,
     &                      LENOPEN,BUFRAW,IRTFLG)

C      FOR 8 BIT INPUT DATA

       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC'

       COMMON /IOBUF/  BUFSPI(NBUFSIZ)

       INTEGER * 8 :: NPIX8,NEED8,IGOT8,NT8
       INTEGER *1  :: BUFRAW(LENOPEN)
       LOGICAL     :: FOLD

       IRECRAW = 0
       IRECSPI = 0
       LOCSPI  = 0
       IGOT8   = 0
       LOCRAW  = LENOPEN 

C      NEED8 IS NUMBER OF BYTE VALUES TO BE READ INCLUDING HEADER
       NEED8 = NPIX8 + IOFFSET

       DO WHILE (IGOT8 .LT. NEED8)
C         PROCESS EACH BYTE VALUE IN HEADER AND IMAGE
          LOCRAW = LOCRAW + 1

          IF (LOCRAW .GT. LENOPEN) THEN
C            NEED TO READ NEW RECORD FROM RAW INPUT FILE
             NVAL    = LENOPEN
             NT8     = NEED8 - IGOT8 
             IF (NT8 .LT. NVAL) NVAL = NT8
             IRECRAW = IRECRAW + 1
             CALL REDLIN8(LUNRAW,BUFRAW,NVAL,IRECRAW,IRTFLG)

             IF (IRTFLG .NE. 0 .AND. IRTFLG .NE. 253) THEN
                WRITE(NOUT,90) IRTFLG,IRECRAW,NVAL
90              FORMAT('*** IO ERROR: ',I5,' RECORD: ',I6,' NVAL: ',I6)
                CALL ERRT(101,'READING FILE',IDUM)
                RETURN
             ENDIF

             LOCRAW = 1
          ENDIF

          IGOT8 = IGOT8 + 1
          IF (IGOT8 .GT. IOFFSET) THEN
C            WANT THIS RAW VALUE FOR AN OUTPUT PIXEL

             LOCSPI = LOCSPI + 1

             IVAL   = BUFRAW(LOCRAW)
             IF (IVAL .LT. 0 .AND. FOLD) THEN
                BUFSPI(LOCSPI) = 256 + IVAL
             ELSE
                BUFSPI(LOCSPI) = IVAL
             ENDIF

             IF (LOCSPI .GE. NSAM) THEN
C               PUT OUT COMPLETED SPIDER RECORD
                IRECSPI = IRECSPI + 1
                CALL WRTLIN(LUNSPI,BUFSPI,NSAM,IRECSPI)
                LOCSPI = 0
             ENDIF
          ENDIF

      ENDDO

      IF (LOCSPI .GT. 0) THEN
C        PUT OUT UNFINISHED SPIDER RECORD
         IRECSPI = IRECSPI + 1
         NVAL    = MIN(NSAM,LOCSPI)
         CALL WRTLIN(LUNSPI,BUF,NVAL,IRECSPI)
      ENDIF

      RETURN
      END

C ------------------------- RAW16TOSPI ----------------------------------

       SUBROUTINE RAW16TOSPI(LUNOLD,LUNNEW,NSAM,NPIX8,IOFFSET,FLIP,
     &                      FOLD,LENOPEN,INBUFB,IRTFLG)

C      FOR 16 BIT INPUT DATA

       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC'

       COMMON /IOBUF/  BUF(NBUFSIZ)

       INTEGER * 8 :: NPIX8,NEED8,IGOT8,NT8
       INTEGER *1  :: INBUFB(LENOPEN)
       LOGICAL     :: FOLD,FLIP

       INTEGER *2  :: I2VAL
       INTEGER *1  :: I1VAL(2)
       EQUIVALENCE     (I2VAL,I1VAL)
 
       IF (MOD(LENOPEN,2) .NE. 0) THEN
          CALL ERRT(102,'PGM ERROR, LENOPEN MUST BE EVEN',LENOPEN)
          RETURN
       ELSEIF (MOD(IOFFSET,2) .NE. 0) THEN
          CALL ERRT(102,'PGM ERROR, OFFSET MUST BE EVEN',IOFFSET)
          RETURN
       ENDIF

       IRECRAW = 0
       IRECSPI = 0
       ILOCSPI = 0
       IRECSPI = 0
       IGOT8   = 0
       ILOCRAW = LENOPEN

C      NEED8 IS NUMBER OF BYTE VALUES TO BE READ INCLUDING HEADER
       NEED8 = NPIX8 * 2  
       NEED8 = NEED8 + IOFFSET 

       DO WHILE (IGOT8 .LT. NEED8)
          ILOCRAW = ILOCRAW + 2
          IF (ILOCRAW .GT. LENOPEN) THEN
C            NEED TO READ NEW RECORD FROM INPUT
             NVAL    = LENOPEN
             NT8     = NEED8 - IGOT8 
             IF (NT8 .LT. NVAL) NVAL = NT8

             IRECRAW = IRECRAW + 1
             CALL REDLIN8(LUNOLD,INBUFB,NVAL,IRECRAW,IRTFLG)

             IF (IRTFLG .NE. 0 .AND. IRTFLG .NE. 253) THEN
                WRITE(NOUT,90) IRTFLG,IRECRAW,NVAL
90              FORMAT('*** IO ERROR: ',I5,' RECORD: ',I6,' NVAL: ',I6)
                CALL ERRT(101,'READING FILE',IDUM)
                RETURN
             ENDIF
             ILOCRAW = 1
          ENDIF

          IF (IGOT8 .GT. IOFFSET) THEN
C            WANT THIS VALUE FOR AN OUTPUT PIXEL

             IF (FLIP) THEN
                I1VAL(2) = INBUFB(ILOCRAW)
                I1VAL(1) = INBUFB(ILOCRAW+1)
             ELSE
                I1VAL(1) = INBUFB(ILOCRAW)
                I1VAL(2) = INBUFB(ILOCRAW+1)
             ENDIF

             ILOCSPI = ILOCSPI + 1
             IF (FOLD .AND. (I2VAL .LT. 0)) THEN
                BUF(ILOCSPI) = 65536 + I2VAL
             ELSE
                BUF(ILOCSPI) = I2VAL
             ENDIF

             IF (ILOCSPI .GE. NSAM) THEN
C               PUT OUT COMPLETED RECORD
                IRECSPI = IRECSPI + 1
                CALL WRTLIN(LUNNEW,BUF,NSAM,IRECSPI)
                ILOCSPI = 0
             ENDIF
          ENDIF
          IGOT8 = IGOT8 + 2

      ENDDO

      IF (ILOCSPI .GT. 0) THEN
C        PUT OUT COMPLETED RECORD
         IRECSPI = IRECSPI + 1
         NVAL    = MIN(NSAM,ILOCSPI)
         CALL WRTLIN(LUNNEW,BUF,NSAM,IRECSPI)
      ENDIF

      RETURN
      END


C------------------------- RAW32TOSPI ----------------------------------

       SUBROUTINE RAW32TOSPI(LUNOLD,LUNNEW,NSAM,NPIX8,IOFFSET,
     &                       NFLIP,LENOPEN,INBUFB,IRTFLG)

C      FOR 32 BIT INPUT DATA

       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC'

       COMMON /IOBUF/  BUF(NBUFSIZ)

       INTEGER * 8  :: NPIX8,NEED8,IGOT8,NT8
        LOGICAL     :: FLIP
       INTEGER *1      INBUFB(LENOPEN)

       REAL *4         R4VAL
       INTEGER *4      I4VAL
       INTEGER *2      I2VAL(4)
       INTEGER *1      I1VAL(4)
       EQUIVALENCE     (R4VAL,I1VAL)
       EQUIVALENCE     (I4VAL,I1VAL)
       EQUIVALENCE     (I2VAL,I1VAL)
 
C      PREVENT NaN ERROR IN DEBUG MODE
       I4VAL     = 0

       IF (MOD(LENOPEN,4) .NE. 0) THEN
          CALL ERRT(102,'PGM ERROR, LENOPENNOT DIVISABLE BY 4!',LENOPEN)
          RETURN
       ENDIF

       IRECRAW = 0
       IRECSPI = 0
       ILOCSPI = 0
       IRECSPI = 0
       IGOT8    = 0
       ILOCRAW = LENOPEN + 1
       LOCB    = 0
       IBYTE   = 0

C      NEED8 IS TOTAL NUMBER OF BYTES TO BE READ INCLUDING HEADER
       NEED8 = NPIX8 * 4 
       NEED8 = NEED8 + IOFFSET

       DO WHILE (IGOT8 .LT. NEED8)
          ILOCRAW = ILOCRAW + 1
          IF (ILOCRAW .GT. LENOPEN) THEN
C            NEED TO READ NEW RECORD FROM INPUT
             NVAL    = LENOPEN
             NT8     = NEED8 - IGOT8 
             IF (NT8 .LT. NVAL) NVAL = NT8

             IRECRAW = IRECRAW + 1
             CALL REDLIN8(LUNOLD,INBUFB,NVAL,IRECRAW,IRTFLG)

             IF (IRTFLG .NE. 0 .AND. IRTFLG .NE. 253) THEN
                WRITE(NOUT,90) IRTFLG,IRECRAW,NVAL
90              FORMAT('*** IO ERROR: ',I5,' RECORD: ',I6,' NVAL: ',I6)
                CALL ERRT(101,'READING FILE',IDUM)
                RETURN
             ENDIF

             ILOCRAW = 1
          ENDIF
C         IGOT8 POINTS TO CURRENT BYTE IN INPUT FILE
          IGOT8 = IGOT8 + 1

          IF (IGOT8 .GT. IOFFSET) THEN
C            WANT THIS VALUE FOR AN OUTPUT PIXEL
             IBYTE = IBYTE + 1

             IF (NFLIP .EQ. 0) THEN
C               NO FLIP
                I1VAL(IBYTE) = INBUFB(ILOCRAW)

             ELSEIF (NFLIP .EQ. -1) THEN
C               INVERT BYTE ORDER
                I1VAL(5-IBYTE) = INBUFB(ILOCRAW)

             ELSEIF (NFLIP .EQ. 1) THEN
C               FLIP BYTES WITHIN WORDS
                IF (IBYTE .EQ. 1) THEN
                   I1VAL(2) = INBUFB(ILOCRAW)
                ELSEIF (IBYTE .EQ. 2) THEN
                   I1VAL(1) = INBUFB(ILOCRAW)
                ELSEIF (IBYTE .EQ. 3) THEN
                   I1VAL(4) = INBUFB(ILOCRAW)
                ELSEIF (IBYTE .EQ. 4) THEN
                   I1VAL(3) = INBUFB(ILOCRAW)
                ENDIF

             ELSEIF (NFLIP .EQ. 2) THEN
C               FLIP BYTES AND WORDS
                IF (IBYTE .EQ. 1) THEN
                   I1VAL(3) = INBUFB(ILOCRAW)
                ELSEIF (IBYTE .EQ. 2) THEN
                   I1VAL(4) = INBUFB(ILOCRAW)
                ELSEIF (IBYTE .EQ. 3) THEN
                   I1VAL(1) = INBUFB(ILOCRAW)
                ELSEIF (IBYTE .EQ. 4) THEN
                   I1VAL(2) = INBUFB(ILOCRAW)
                ENDIF
             ENDIF

             IF (IBYTE .EQ. 4) THEN
C               TIME TO OUTPUT THE BYTES
                ILOCSPI      = ILOCSPI + 1
                BUF(ILOCSPI) = R4VAL
                IBYTE        = 0
#ifdef NEVER
                IF ((ILOCSPI .LE. 4) .AND. IRECSPI .EQ. 1) THEN
                   WRITE(6,*) 'ILOCSPI : ',ILOCSPI,' I1VALS: ',
     &                 I1VAL(1),I1VAL(2),I1VAL(3),I1VAL(4)
                   WRITE(6,*) 'I2VAL: ',I2VAL(1),I2VAL(2)
                   WRITE(6,*) 'R4VAL: ',R4VAL
                   WRITE(6,*) 'I4VAL: ',I4VAL
                   WRITE(6,*) ' '
                ENDIF
#endif
             ENDIF

             IF (ILOCSPI .GE. NSAM) THEN
C               PUT OUT COMPLETED RECORD
                IRECSPI = IRECSPI + 1
                CALL WRTLIN(LUNNEW,BUF,NSAM,IRECSPI)
                ILOCSPI = 0

             ENDIF
          ENDIF
      ENDDO

      IF (ILOCSPI .GT. 0) THEN
C        PUT OUT REMAINING RECORD
         IRECSPI = IRECSPI + 1
         CALL WRTLIN(LUNNEW,BUF,NSAM,IRECSPI)
      ENDIF

      RETURN
      END

C dd if=avg000.dat of=tes032.dat count=16384 bs=1
C tail +1025c avg000.dat > tes032.dat

#ifdef NEVER
             if (icount .le. 1) THEN
               write(6,*) 'inbufb(',ILOCRAW,'): ',inbufb(ILOCRAW)
               it = ilocin + 1
               write(6,*) 'inbufb(',it,'): ',inbufb(it)
               write(6,*) 'i1val(1):  ', i1val(1),' i1val(2): ',i1val(2)
               write(6,*) 'i1val(3):  ', i1val(3),' i1val(4): ',i1val(4)
               write(6,*) 'i4val: ',i4val

c               write(6,*) '  buf(',ilocout,'): ',buf(ilocout)
c               write(6,*) 'i2val(1):  ', i2val(1),' i2val(2): ',i2val(2)
c               write(6,*) 'ilocin: ', ilocin,' ilocout: ',ilocout
c               write(6,*) 'recin:  ', irecin,' irow: ',irow
               write(6,*) ' -------------',icount
               icount = icount + 1
             endif
#endif



C------------------------- RAW32INTTOSPI ----------------------------------

       SUBROUTINE RAW32INTTOSPI(LUNOLD,LUNNEW,NSAM,NPIX8,IOFFSET,
     &                       NFLIP,FOLD,LENOPEN,IRTFLG)

C      FOR 32 BIT INTEGER INPUT DATA (SOME OF THIS IS PROBABLY 
C      REDUNDANT AND NOT NECESSARY!! al

       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC'

       INTEGER * 8  :: NPIX8,NEED8,IGOT8,NT8
       REAL         :: RBUF(NBUFSIZ)
       LOGICAL      :: FOLD                    

       INTEGER *1   :: INBUFB(LENOPEN)

       INTEGER *4   :: I4VAL
       INTEGER *1   :: I1VAL(4)
       EQUIVALENCE (I4VAL,I1VAL)
 
C      PREVENT NaN ERROR IN DEBUG MODE
       I4VAL     = 0

       IF (MOD(LENOPEN,4) .NE. 0) THEN
          CALL ERRT(102,'PGM ERROR, LENOPENNOT DIVISABLE BY 4!',LENOPEN)
          RETURN
       ENDIF

       IRECRAW = 0
       IRECSPI = 0
       ILOCSPI = 0
       IRECSPI = 0
       IGOT8   = 0
       ILOCRAW = LENOPEN + 1
       LOCB    = 0
       IBYTE   = 0

C      NEED8 IS TOTAL NUMBER OF BYTES TO BE READ INCLUDING HEADER
       NEED8 = NPIX8 * 4 
       NEED8 = NEED8 + IOFFSET

       DO WHILE (IGOT8 .LT. NEED8)
          ILOCRAW = ILOCRAW + 1
          IF (ILOCRAW .GT. LENOPEN) THEN
C            NEED TO READ NEW RECORD FROM INPUT
             NVAL    = LENOPEN
             NT8     = NEED8 - IGOT8 
             IF (NT8 .LT. NVAL) NVAL = NT8
             IRECRAW = IRECRAW + 1
             CALL REDLIN8(LUNOLD,INBUFB,NVAL,IRECRAW,IRTFLG)

             IF (IRTFLG .NE. 0 .AND. IRTFLG .NE. 253) THEN
                WRITE(NOUT,90) IRTFLG,IRECRAW,NVAL
90              FORMAT('*** IO ERROR: ',I5,' RECORD: ',I6,' NVAL: ',I6)
                CALL ERRT(101,'READING FILE',IDUM)
                RETURN
             ENDIF

             ILOCRAW = 1
          ENDIF
C         IGOT8 POINTS TO CURRENT BYTE IN INPUT FILE
          IGOT8 = IGOT8 + 1

          IF (IGOT8 .GT. IOFFSET) THEN
C            WANT THIS VALUE FOR AN OUTPUT PIXEL
             IBYTE = IBYTE + 1

             IF (NFLIP .EQ. 0) THEN
C               NO FLIP  (64)
                I1VAL(IBYTE) = INBUFB(ILOCRAW)

             ELSEIF (NFLIP .EQ. -1) THEN
C               INVERT BYTE ORDER   (-64)
                I1VAL(5-IBYTE) = INBUFB(ILOCRAW)

             ELSEIF (NFLIP .EQ. 1) THEN
C               FLIP BYTES WITHIN WORDS   (65)
                IF (IBYTE .EQ. 1) THEN
                   I1VAL(2) = INBUFB(ILOCRAW)
                ELSEIF (IBYTE .EQ. 2) THEN
                   I1VAL(1) = INBUFB(ILOCRAW)
                ELSEIF (IBYTE .EQ. 3) THEN
                   I1VAL(4) = INBUFB(ILOCRAW)
                ELSEIF (IBYTE .EQ. 4) THEN
                   I1VAL(3) = INBUFB(ILOCRAW)
                ENDIF

             ELSEIF (NFLIP .EQ. 2) THEN
C               FLIP BYTES AND WORDS   (-65)
                IF (IBYTE .EQ. 1) THEN
                   I1VAL(3) = INBUFB(ILOCRAW)
                ELSEIF (IBYTE .EQ. 2) THEN
                   I1VAL(4) = INBUFB(ILOCRAW)
                ELSEIF (IBYTE .EQ. 3) THEN
                   I1VAL(1) = INBUFB(ILOCRAW)
                ELSEIF (IBYTE .EQ. 4) THEN
                   I1VAL(2) = INBUFB(ILOCRAW)
                ENDIF
             ENDIF


             IF (IBYTE .EQ. 4) THEN
C               TIME TO OUTPUT THE BYTES
                ILOCSPI       = ILOCSPI + 1
                IF (FOLD .AND. (I4VAL .LT. 0)) THEN
#if defined (SP_GFORTRAN)
                   RBUF(ILOCSPI) = 2.0**32 + I4VAL
#else
                   RBUF(ILOCSPI) = 2**32 + I4VAL
#endif
                ELSE
                   RBUF(ILOCSPI) = I4VAL
                ENDIF
                IBYTE         = 0
             ENDIF


             IF (ILOCSPI .GE. NSAM) THEN
C               PUT OUT COMPLETED RECORD
                IRECSPI = IRECSPI + 1
                CALL WRTLIN(LUNNEW,RBUF,NSAM,IRECSPI)
                ILOCSPI = 0
             ENDIF
          ENDIF
      ENDDO

      IF (ILOCSPI .GT. 0) THEN
C        PUT OUT REMAINING RECORD
         IRECSPI = IRECSPI + 1
         CALL WRTLIN(LUNNEW,RBUF,NSAM,IRECSPI)
      ENDIF

      RETURN
      END

C ------------------------- RAWTOSPI ----------------------------------

       SUBROUTINE RAWTOSPI(LUNOLD,LUNNEW,NSAM,NPIX8,IOFFSET,
     &                     MODE,FOLD,FLIP,ISIGB,IRTFLG)

       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC'

       COMMON /IOBUF/  FBUF(NBUFSIZ)

       INTEGER * 8      :: NPIX8,NBYTES8
       CHARACTER(LEN=1) :: CINBUF
c       EQUIVALENCE     (CINBUF,I1BUF)

       REAL *4         F4BUF
       INTEGER * 1     I1BUF(4)
       EQUIVALENCE     (F4BUF,I1BUF)

       LOGICAL  ::     FOLD,FLIP
       INTEGER * 1     I1VAL(4)
       INTEGER * 4     I4VAL
       EQUIVALENCE     (I4VAL,I1VAL(1))
       EQUIVALENCE     (FVAL,I1VAL(1))
 

C      CLEAN THE BUFFER, THIS ASSUMES THE INTENSITIES IN THE RAW IMAGE 
C      ARE ALL POSITIVE NUMBERS				---- liy
       DO I = 1, 4
	  I1VAL(I) = 0
       ENDDO

C      ILOCSPI IS POINTER TO CURRENT NUMBER OF IMAGE VALUES OUTPUT TO
C      CURRENT OUTPUT RECORD
       ILOCSPI = 0

C      IRECSPI IS NUMBER OF CURRENT OUTPUT RECORD
       IRECSPI   = 0

C      BYTE POINTER FOR 16 & 32 BIT INPUT
       IBITE  = 0

       NBYTES8 = NPIX8 * ABS(MODE)
       NBYTES8 = NBYTES8 / 8

#if  defined (SP_NT) || defined (__osf__) || defined (__linux__) 
C      MUST ADJUST BYTE ORDER IN I4VAL
       LOC1 = 1
       LOC2 = 2
#else
       LOC1 = 3
       LOC2 = 4
#endif

       IF (IOFFSET .GT. 0) THEN
C         SKIP HEADER OFFSET

          DO ILOC = 1,IOFFSET
C            THE (1) IS ESSENTIAL ON THESE READS
             READ(LUNOLD,REC=ILOC,IOSTAT=IRTFLG) I1BUF(1)
             IF (IRTFLG .NE. 0) THEN
                WRITE(NOUT,90) IRTFLG,ILOC
90              FORMAT('*** ERROR: (',I5,
     &                 ') READING HEADER BYTE: ',I4)
                CALL ERRT(100,'RAWTOSPI',NE)
                RETURN
             ENDIF
          ENDDO
       ENDIF

C      DEC 08 possible overflow of iloc for large nbytes8???? al
       DO ILOC =  IOFFSET + 1,IOFFSET + NBYTES8
          READ(LUNOLD,REC=ILOC,IOSTAT=IRTFLG) I1BUF(1)

          IF (IRTFLG .NE. 0) THEN
             WRITE(NOUT,91) IRTFLG, ILOC - IOFFSET
91           FORMAT('*** ERROR: (',I5,') READING IMAGE BYTE: ',I8)
             CALL ERRT(100,'RAWTOSPI',NE)
             RETURN
          ENDIF

          IF (MODE .EQ. 8) THEN
C            INPUT IS 8 BIT ------------------------------- 8 BIT IN
C            CONVERT TO 32 BIT FLOATING POINT 
             ILOCSPI = ILOCSPI + 1
             IVAL    = I1BUF(1)
             IF (IVAL .LT. 0) IVAL = 256 + IVAL
             FBUF(ILOCSPI) = IVAL

          ELSE IF (MODE .EQ. 16) THEN
C            INPUT IS 16 BIT ----------------------------- 16 BIT IN
C            CONVERT TO 32 BIT FLOATING POINT 
             IBITE = IBITE + 1
             IF (ISIGB .EQ. 1 .AND. IBITE .EQ. 1) THEN
                I1VAL(LOC1) = I1BUF(1)
             ELSE IF (ISIGB .EQ. 1 .AND. IBITE .EQ. 2) THEN
                I1VAL(LOC2) = I1BUF(1)
             ELSE IF (IBITE .EQ. 1) THEN
                I1VAL(LOC2) = I1BUF(1)
             ELSE IF (IBITE .EQ. 2) THEN
                 I1VAL(LOC1) = I1BUF(1)
             ENDIF                 
             IF (IBITE .EQ. 2) THEN
C               PUT OUT COMPLETED VALUE
                IF (FOLD) THEN
                   IF (I4VAL .LE. 32767) THEN
                      I4VAL = I4VAL + 32768
                   ELSE
                      I4VAL = I4VAL - 32768
                   ENDIF
                ENDIF
                ILOCSPI       = ILOCSPI + 1
                FBUF(ILOCSPI) = I4VAL
                IBITE         = 0
             ENDIF

          ELSE IF (MODE .EQ. 32) THEN
C            INPUT IS 32 BIT ---------------------------- 32 BIT IN
C            CONVERT TO 32 BIT FLOATING POINT 
             IBITE        = IBITE + 1
             IF (FLIP) THEN
                I1VAL(5 - IBITE) = I1BUF(1)
             ELSE
                I1VAL(IBITE) = I1BUF(1)
             ENDIF

             IF (IBITE .GE. 4) THEN
C               PUT OUT COMPLETED VALUE
                ILOCSPI       = ILOCSPI + 1
                FBUF(ILOCSPI) = FVAL
                IBITE         = 0
             ENDIF

         ENDIF

         IF (ILOCSPI .EQ. NSAM) THEN
C           PUT OUT COMPLETED RECORD
            IRECSPI = IRECSPI + 1
            CALL WRTLIN(LUNNEW,FBUF,NSAM,IRECSPI)
            ILOCSPI = 0
         ENDIF
      ENDDO

      IF (ILOCSPI .GT. 0) THEN
C        PUT OUT COMPLETED RECORD
         IRECSPI = IRECSPI + 1
         CALL WRTLIN(LUNNEW,FBUF,NSAM,IRECSPI)
      ENDIF

C     SET NO ERROR FLAG
      IRTFLG = 0

      RETURN
      END

