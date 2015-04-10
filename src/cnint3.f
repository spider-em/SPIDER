
C **********************************************************************
C
C   CNINT3.F    -- CREATED MAY 87
C                  CHANGED TO POSTSCRIPT MAR 99 ARDEAN LEITH
C                  MAXNAM                JUL 14 ARDEAN LEITH
C
C **********************************************************************
C * AUTHOR: ArDean Leith                                               *
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
C      CNINT3(LUNIM,LUNPOS,MAXDIM)
C
C      PURPOSE:     READS  MULTIPLE SPIDER PICTURE FILES, EXTRACTS CONTOURS, 
C                   AND TRANSFERS THEM TO A POSTSCRIPT FILE
C
C      PARAMETERS:  LUNIM   = SPIDER IMAGEE INPUT UNIT           (SENT)
C                   LUNPOS  = POSTSCRIPT OUTPUT UNIT             (SENT)
C                   MAXDIM  = COMMON BUFFER SIZE                 (SENT)
C
C      CALLED BY:   PLOT1
C
C      CALL TREE:   PLOT1..CNINT3..CNTUR..CNSCAN..CNTRCE..CNSTUFF
C                                               ..CNCALC
C                      
C23456789012345678901234567890123456789012345678901234567890123456789012
C--********************************************************************

       SUBROUTINE CNINT3(LUNIM,LUNPOS,NSAM,NROW,NSLICE,
     &                   FMINT,FMAXT,MAXDIM)
  
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

C-------- START OF EM-PLOTT-COMMON-------------------------------------
C     INTEGERS
      COMMON /CONT1/ ICALL, IDIDIT, IDONE, IDX, IDY, ILINE, INTT,
     &               IRCD, ISS, ISTART, ISUM1, ISUM2, ISUM3, IT, IV, 
     &               IXX1, IXX2, IXX3, IX, IY, JSUM1, JSUM2, JSUM3, JT,
     &               LEVEL, LW, M, MF, MI, MT, N, NDIV, NF, NI, NT, NW

C     FLOATING POINT
      COMMON /CONT2/ APDIV, APDIVX, CV, DL, PY, RA, RC, RS, SKALE, THE,
     &               SX, SY, DENSL

C     ARRAYS
      COMMON /CONT3/ INCX(3), IORGX(3), INX(8),
     &               INY(8),  IPT(3,3), IMAP(12), NG(3), NP(3)

      COMMON /CONT4/ CTRI(6),FCTR(6),CTRDEL(6),ICNDSH(6),ICNCEL

C--------END OF EM-PLOTT-COMMON----------------------------------------

        PARAMETER        (MAXIRR  = 80000)
        PARAMETER        (NSIZE   = 2000)

        COMMON           BUF(MAXIRR),X(NSIZE),Y(NSIZE),AM(1)
        COMMON /SPI_BUF/ DATA(3,NSIZE),WORK(2,NSIZE)

        COMMON /POLY/    MINPTS,ISLICE,ICNT
      
        CHARACTER(LEN=MAXNAM) :: FILPOS,FILPAT

        CHARACTER *1     LABEL1
        CHARACTER *10    LABEL2
        CHARACTER *28    LABEL3
        CHARACTER *28    LABEL4
        CHARACTER *12    PROMPT

C*      NMAXI IS DIMENSION OF X AND Y ARRAYS
        NMAXI = NSIZE

C       MINIMUM NUMBER OF POINTS SAVED IN CONTOUR FILE / CONTOUR
        MINPTS = 4

C       MAX NUMBER OF POINTS DESIRED ON CONTOUR FILE CONTOUR
        MAXPTS  = 800

        MAXPIX = MAXDIM - 3 * NSIZE 
        IF ((NSAM * NROW) .GT. MAXPIX) THEN
            WRITE(NOUT,9945) MAXPIX
 9945       FORMAT(' SORRY, PGM LIMITED TO: ',I8,' PIXELS'/)
            CALL ERRT(100,'CNINT3',NE)
            GOTO 9999
        ENDIF

C       DISPLAY MAX AND MIN VALUE OF PICTURE , ASK FOR THE CONTOUR LEVELS
        WRITE(NOUT,106) FMINT,FMAXT
  106   FORMAT(' IMAGE RANGE: ',1PG11.3,'....',1PG11.3)

   11   CALL RDPRM2S(BLEV,ELEV,NOT_USED,
     &             'STARTING AND ENDING CONTOUR LEVELS',IRTFLG)
        IF (IRTFLG .EQ. -1) RETURN

   12   CALL RDPRM1S(RINC,NOT_USED,
     &    'CONTOUR LEVEL INCREMENT (USE INCR. > END FOR ONE LEVEL)',
     &     IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 11

        IF (FCHAR(4:4) .EQ. 'S') THEN
C          GET NAME OF OUTPUT FILE AND OPEN AS SEQUENTIAL FORMATTED
           CALL OPAUXFILE(.TRUE.,FILPOS,'ssr',LUNPOS,0,'N',
     &                    'OUTPUT',.TRUE.,IRTFLG)
           IF (IRTFLG .EQ. -1) GOTO 12
           IF (IRTFLG .NE. 0) GOTO 9999
           ICNT = 0
          
        ELSE 
           IF (NSLICE .GT. 1) THEN
C             GET NAME OF POSTSCRIPT FILE TEMPLATE
              CALL FILSEQP(FILPAT,NLET,ILIST,0,NUM,
     &                  'OUTPUT FILE TEMPLATE (E.G. SLI***)',IRTFLG)

C             FIND FILE NAME FOR FIRST SLICE
              CALL FILGET(FILPAT,FILPOS,NLET,ISLICE,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999

C             OPEN THE FIRST SLICE OUTPUT FILE (SEQ. FORMATTED)
              CALL OPAUXFILE(.FALSE.,FILPOS,'ps',LUNPOS,0,'N',
     &                    ' ',.TRUE.,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999

           ELSE
C             OPEN POSTSCRIPT FILE
              CALL OPAUXFILE(.TRUE.,FILPOS,'ps',LUNPOS,0,'N',
     &                    'OUTPUT',.TRUE.,IRTFLG)
           ENDIF
           IF (IRTFLG .EQ. -1) GOTO 12
           IF (IRTFLG .NE. 0) GOTO 9999

           CALL POSTRT(-LUNPOS)

C          INITIALIZE & SET SCALING PARAMETERS FOR POSTSCRIPT
           XLL = 0  
           YLL = -30
           XUR = NSAM
           YUR = NROW

           CALL POSCALE(LUNPOS, 1.0,1.0,  XLL,YLL,  XUR,YUR)
           CALL POGETSCALE(LUNPOS,SCALE,SCALED)

        ENDIF


C       SET X AND Y EXTENT FOR CONTOUR FILE
        XMIN  = 0.0
        YMIN  = 0.0
        XMAX  = NSAM
        YMAX  = NROW

C       PUT IN FRAME AROUND PLOT
        DATA(1,1) = 1.0
        DATA(2,1) = 1.0

        DATA(1,2) = NROW
        DATA(2,2) = 1.0

        DATA(1,3) = NROW
        DATA(2,3) = NSAM

        DATA(1,4) = 1.0
        DATA(2,4) = NSAM

        DATA(1,5) = 1.0
        DATA(2,5) = 1.0

        NDATA     = 5

        IF (FCHAR(4:4) .EQ. 'S') THEN
C          PUSH FRAME DATA INTO FILE
           CALL SSPUSH(LUNPOS,DATA,NDATA,0,IRTFLG)
        ELSE
C          SET TEXT CHARACTARISTICS FOR LABEL
           ITSIZE = 12
           ITANGL = 0
           JUST   = 0

C          SCALING DATA LABELS
           YPOS1  = -6.0 * SCALED
           XPOS1  = 1.0
           LABEL1 = '1'

           XPOS2  = NSAM
           YPOS2  = -6.0* SCALED
           CALL INTTOCHAR(NSAM,LABEL2,NLET2,1)

C          SLICE LABEL
           XPOS3       = NSAM / 2.0
           YPOS3       = -60.0 * SCALED
           LABEL3(1:7) = 'SLICE: '

C          THRESHOLD LABEL
           XPOS4        = XPOS3
           YPOS4        = -30.0 * SCALED
           LABEL4(1:16) = 'STARTING LEVEL: '
           WRITE(LABEL4(17:),8000,IOSTAT=IERR) BLEV
8000       FORMAT(1PG10.3)
           NLET4 = 27
        ENDIF

        DO ISLICE = 1,NSLICE

          IF (FCHAR(4:4) .NE. 'S') THEN
C            POSTSCRIPT OUTPUT

             IF (ISLICE .GT. 1) THEN
C               FIND FILE NAME FOR THIS SLICE
                CALL FILGET(FILPAT,FILPOS,NLET,ISLICE,IRTFLG)
                IF (IRTFLG .NE. 0) GOTO 9999

C               OPEN THE SLICE OUTPUT FILE (SEQ., FORMATTED)
                CALL OPAUXFILE(.FALSE.,FILPOS,'ps',LUNPOS,0,'N',
     &                    ' ',.TRUE.,IRTFLG)
                IF (IRTFLG .NE. 0) GOTO 9999

C               INITIALIZE POSTSCRIPT FILE
                CALL POSTRT(-LUNPOS)
                CALL POSCALE(LUNPOS, 1.0,1.0,  XLL,YLL,  XUR,YUR)
             ENDIF

C            PUSH FRAME DATA INTO FILE
             CALL POARAYF(LUNPOS,DATA,NDATA,.FALSE.,.FALSE.)

C            LABEL SCALING
             CALL POTEX(LUNPOS,LABEL1,1,XPOS1,YPOS1,ITSIZE,ITANGL,JUST)
             CALL POTEX(LUNPOS,LABEL2,NLET2,XPOS2,YPOS2,
     &                  ITSIZE,ITANGL,JUST)

             IF (NSLICE .GT. 1) THEN
C               LABEL SLICE
                CALL INTTOCHAR(ISLICE,LABEL3(8:),NLET3,1)
                CALL POTEX(LUNPOS,LABEL3,NLET3+7,XPOS3,YPOS3,
     &                     ITSIZE,ITANGL,1)
             ENDIF

C            LABEL THRESHOLD
             CALL POTEX(LUNPOS,LABEL4,NLET4,XPOS4,YPOS4,
     &                  ITSIZE,ITANGL,1)
             
          ENDIF

C         READ THE SPIDER FILE INTO AM ARRAY
          IREC1 = (ISLICE-1)*NROW
          ILOC  = 1
C         INVERT FOR CORRECT MIRRORING AS A SPIDER IMAGE
          DO  IREC = IREC1+NROW,IREC1,-1
             CALL REDLIN (LUNIM,AM(ILOC),NSAM,IREC)
             ILOC = ILOC + NSAM
          ENDDO

C         CNTUR EXTRACTS THE CONTOUR FROM THE SPIDER FILE AND PLACES
C         IT IN THE OUTPUT FILE
          CALL CNTUR(AM,NSAM,NROW,BLEV,ELEV,RINC,BUF,X,Y,
     &                NMAXI,LUNPOS,.FALSE.,MAXPTS,MAXIRR)

C         ECHO OUTPUT FILE NAME TO TERMINAL
          NLETP = LNBLNKN(FILPOS)
          WRITE(NOUT,*) ' PLOT PLACED IN: ',FILPOS(1:NLETP)

          IF (FCHAR(4:4) .NE. 'S') THEN
C            CLOSE THE POSTSCRIPT-FILE
             CALL POEND(LUNPOS)
             CLOSE(LUNPOS)
C            GIVE IMMEDIATE POSTSCRIPT PLOT (IF DESIRED)
             CALL POPRINT(FILPOS(1:NLETP))
          ENDIF
        ENDDO


9998    CLOSE(LUNPOS)
9999    CLOSE(LUNIM)
        RETURN
        END
