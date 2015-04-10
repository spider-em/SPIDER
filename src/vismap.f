C ++********************************************************************
C
C VISMAP      FILENAME OVERWRITE BUG FIXED       NOV 2013 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2013  Health Research Inc.,                         *
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
C                                                                      *
C  VISMAP(LUN,LUNO,NDOC1,NDOC2)                                                              *
C                                                                      *
C  PURPOSE:  CREATE VISUAL MAP IN SPIDER IMAGE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE VISMAP(LUN,LUNO,NDOC1,NDOC2)

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC'

        PARAMETER (MAXREG   = 25)
        PARAMETER (MAXKEY   = 9992)
        PARAMETER (MAXIMDIM = 128)

C       WORKSTATION DISPLAY SIZE:
        PARAMETER (MAXXDP = 1024)
        PARAMETER (MAXYDP = 800)

        CHARACTER (LEN=MAXNAM) :: DOCF1,DOCF2,FILPAT,FILN1,FIL2 	
        CHARACTER (LEN=1)      :: NULL = CHAR(0)
        LOGICAL                :: WINDX,WINDY,PADX,PADY
        DOUBLE PRECISION       :: DAVX,DAVY,DAVX2,DAVY2,DAVD,DAVD2

 	REAL    :: BUF(1024), DBUF(MAXREG,MAXKEY),
     &             BUFX(MAXKEY),BUFY(MAXKEY),PLIST(MAXREG),
     &             PORTION(MAXXDP,132)

 	INTEGER :: IBINK(MAXKEY),IBIN(MAXKEY) ,NKEY(MAXKEY)

        MAXXD = MAXXDP
        MAXYD = MAXYDP

        WRITE(NOUT,401) MAXXD,MAXYD
401     FORMAT('  DEFAULT OUTPUT IMAGE SIZE:',I0,' BY ',I0)

        CALL RDPRMI(NEWMAXXD,NEWMAXYD,NOT_USED,
     &             'NEW IMAGE SIZE: X,Y (<RET>=DEF)')

        IF (NEWMAXXD.NE.0) MAXXD = NEWMAXXD
        IF (NEWMAXYD.NE.0) MAXYD = NEWMAXYD
        NSAMO = MAXXD
        CALL FILERD(DOCF1,NLET,DATEXC,'INPUT DOCUMENT',IRTFLG)       
        CALL FILERD(DOCF2,NLET,DATEXC,'OUTPUT DOCUMENT',IRTFLG)
       
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILPAT,LUN,'O',IFORM,NSAM,NROW,NSLICE,
     &                   MAXIM,'INPUT',.FALSE.,IRTFLG)
        CLOSE(LUN)

        WRITE(NOUT,103)
103     FORMAT('  RESULT WILL BE AUTOMATICALLY WINDOWED/PADDED')

        CALL RDPRMI(IX,IY,NOT_USED,
     &         'NUMBER OF DIVISIONS IN X,Y')
 
        CALL RDPRM(STNDRT,NOT_USED,
     &         'STANDARD DEVIATION (2.3=DEF)')

        CALL RDPRM2(RANGEU,RANGEL,NOT_USED,
     &         'UPPER/LOWER IMAGE THRESHOLD IN UNITS OF SIGMA') 

        IF (STNDRT.EQ.0) STNDRT=2.3

        CALL FILERD(FIL2,NLET,NULL,'OUTPUT',IRTFLG)
        NSLIC = 1

        CALL RDPRMI(KEY1,KEY2,NOT_USED,
     &             'FIRST & LAST KEY=IMAGE NUMBERS')
        IF (KEY2.GT.MAXKEY) THEN 
           WRITE(NOUT,102)
102        FORMAT('  ENDING KEY TOO LARGE, YOU MIGHT CONSIDER',/
     &            ' TO ASK YOUR PROGRAMMER, IF HE MIGHT BE SO'/
     &            ' KIND TO CHANGE THE PARAMETERS IN VISMAP.F')
           RETURN
        ENDIF

	CALL RDPRMI(ICOLX,ICOLY,NOT_USED,
     &        'COLUMN #S IN DOC. FILE USED FOR X,Y COORD.')
        NREG=MAX(ICOLX,ICOLY)

C       PREPARE WINDOWING AND PADDING
        WINDX = .FALSE.
        WINDY = .FALSE.
        PADX  = .FALSE.
        PADY  = .FALSE.

C       CALCULATE AVAILABLE SPACE FOR IMAGES (BORDERS 2 PIXELS):
        IXSEC   = MAXXD/IX
        IYSEC   = MAXYD/IY
        IYSEC   = MIN(IYSEC,130)
        IXDIM   = MAXXD/IX-2
        IYDIM   = MAXYD/IY-2
        NROWOUT = IYSEC*IY

        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FIL2,LUNO,'U',IFORM,MAXXD,NROWOUT,NSLIC,
     &                   MAXIM,' ',.FALSE.,IRTFLG)

C       FIND OUT IF PADDING/WINDOWING NECESSARY.
        IF (IXDIM .LE. NSAM) WINDX = .TRUE.
        IF (IYDIM .LE. NROW) WINDY = .TRUE.
        IF (.NOT. WINDX) PADX =.TRUE.
        IF (.NOT. WINDY) PADY =.TRUE.

C        DETERMINE STARTING/ENDING INDICES 
        IXA = 1
        IXE = NSAM

        IF (WINDX) THEN
           IXA = (NSAM-IXDIM)/2+1
           IXE = IXA+IXDIM
        ENDIF
        IYA = 1
        IYE = NROW
        IF (WINDY) THEN 
           IYA = (NROW-IYDIM)/2+1
           IYE = IYA+IYDIM
        ENDIF
        IXSTRT = 2
        IXEND = IXDIM+1
        IF (PADX) THEN
           IXSTRT = (IXDIM-NSAM)/2+1
           IXEND  = IXSTRT+NSAM
        ENDIF
        IYSTRT = 2
        IYEND  = IYDIM+1
        IF (PADY) THEN
           IYSTRT = (IYDIM-NROW)/2+1
           IYEND  = IYSTRT+NROW
        ENDIF

        !WRITE(NDAT,221) 
!     &     IXA,IXE,IXSTRT,IXEND,IXDIM,IYA,IYE,IYSTRT,IYEND,IYDIM
221     FORMAT(' ',
     &    'WINDOW/PAD X: INPUT: IXA,IXE, OUTPUT IXSTRT,IXEND,IXDIM'/10X
     &    ,6(1X,I6))
 
C       GET ALL THE COORDINATES:

        ILINE = 0
	NOPEN = 0

C	ADD AN EXTENSION TO THE FILENAME

        DO I=KEY1,KEY2
           CALL UNSDAL(DOCF1,NOPEN,NDOC1,I,PLIST,NREG,
     &         DBUF,MAXKEY,MAXREG,LKEY,LERR) 

           IF (DBUF(1,I) .NE. 0) THEN
              ILINE       = ILINE + 1
              NKEY(ILINE) = I
              BUFX(ILINE) = PLIST(ICOLX)
              BUFY(ILINE) = PLIST(ICOLY)
              !WRITE(NOUT,101) BUFX(ILINE),BUFY(ILINE),ILINE,I
101           FORMAT(' ',2(2X,F12.6),2(2X,I5))
           ENDIF

        END DO

        NOPEN = 0
        CLOSE(NDOC1)

C       FIRST SORT INTO BINS
        XMIN = BUFX(1)
        XMAX = BUFX(1)
        YMIN = BUFY(1)
        YMAX = BUFY(1)

C       DETERMINE MAXIMUM/MINIMUM X AND Y COORDINATES
        DAVX  = 0
        DAVX2 = 0
        DAVY  = 0
        DAVY2 = 0

        DO K=2,ILINE
           BX    = BUFX(K)
           BY    = BUFY(K)
           DAVX  = DAVX+BX
           DAVX2 = DAVX2+BX*DBLE(BX)
           DAVY  = DAVY+BY
           DAVY2 = DAVY2+BY*DBLE(BY)
           IF (BX .LT. XMIN) XMIN = BX
           IF (BX .GT. XMAX) XMAX = BX
           IF (BY .LT. YMIN) YMIN = BY
           IF (BY .GT. YMAX) YMAX = BY
        END DO

        YR = (YMAX-YMIN)*1.01
        XR = (XMAX-XMIN)*1.01                                 

C       CALCULATE SIGMA AND AVERAGE
        FLINE = FLOAT(ILINE)
        DAVX  = DAVX/FLINE
        DAVY  = DAVY/FLINE

        SIGX  = DSQRT((DAVX2-DAVX*DAVX/FLINE)/DBLE(FLINE-1))
        SIGY  = DSQRT((DAVY2-DAVY*DAVY/FLINE)/DBLE(FLINE-1))

        XRA   = 2.*STNDRT*SIGX
        YRA   = 2.*STNDRT*SIGY

        IF (XRA.GT.XR) THEN 
           XRANGE = XR
           TLX    = SIGN(1.,XMIN)*(ABS(XMIN)*1.01)
           TUX    = SIGN(1.,XMAX)*(ABS(XMAX)*1.01)
        ELSE
           XRANGE  = XRA
           XCENTER = DAVX
           XDEV    = XRANGE/2
           TLX     = XCENTER-XDEV
           TUX     = XCENTER+XDEV
        ENDIF

        IF (YRA.GT.YR) THEN 
           YRANGE  = YR
           TLY     = SIGN(1.,YMIN)*(ABS(YMIN)*1.01)
           TUY     = SIGN(1.,YMAX)*(ABS(YMAX)*1.01)
        ELSE
           YRANGE  = YRA
           YCENTER = DAVY
           YDEV    = YRANGE/2
           TLY     = YCENTER-YDEV
           TUY     = YCENTER+YDEV
        ENDIF

        XINTER = XRANGE/FLOAT(IX)
        YINTER = YRANGE/FLOAT(IY) 

        WRITE(NDAT,*) ' '

C       WRITE(NOUT,345) DAVX,DAVY,XDEV,YDEV,XRANGE,YRANGE,TLX,TUX,
C     &                 TLY,TUY,XINTER,YINTER
        WRITE(NDAT,345) XMAX,XMIN,XR,YMAX,YMIN,YR, 
     &                  DAVX,DAVY,XDEV,YDEV,XRANGE,YRANGE,TLX,TUX
     &                  ,TLY,TUY,XINTER,YINTER,IX,IY
345     FORMAT(
     &  '  X-DRECTION MAX: ',E10.4,' MIN: ',E10.4,' RANGE ',E10.4,/
     &  '  Y-DRECTION MAX: ',E10.4,' MIN: ',E10.4,' RANGE ',E10.4, /
     &  '  AVERAGES              X: ',D10.4,' Y: ',D10.4,/
     &  '  DEVIATION             X: ',E10.4,' Y: ',E10.4, /
     &  '  RANGES                X: ',E10.4,' Y: ',E10.4, /
     &  '  LOWER/UPPER BOUDARIES X:',2(1X,E10.4),' Y: ',2(1X,E10.4)/
     &  '  INTERVALS             X: ',E10.4,' Y: ',E10.4,
     &  '  PARITIONS             X: ', I5  ,' Y: ', I5 )
        WRITE(NDAT,*) ' '
C
C       PUT INTO BINS
        DO K=1,ILINE
           INTX     = ((BUFX(K)-TLX)/XRANGE)*IX+1
           INTY     = ((BUFY(K)-TLY)/YRANGE)*IY+1
           IBINK(K) = NKEY(K)
           IBIN(K)  = (INTY-1)*IX+INTX
C          WRITE(NOUT,346) BUFX(K),BUFX(K)-TLX,
C     &                    INTX,BUFY(K),BUFY(K)-TLY,INTY 
   
c           WRITE(NDAT,346) BUFX(K),BUFX(K)-TLX,INTX,BUFY(K)
c     &        ,BUFY(K)-TLY,INTY,NKEY(K),IBIN(K)    
346        FORMAT(' ',
     &  'BUFX(K)   ,BUFX(K)-TLX,INTX,BUFY(K),BUFY(K)-TLY,INTY,NKEY,IBIN'    
     &     ,/2(1X,E10.4),I5,2(1X,E10.4),3I5)
           IF (INTX .LT. 1 .OR. INTX .GT. IX) IBIN(K)=10000 
           IF (INTY .LT. 1 .OR. INTY .GT. IY) IBIN(K)=10000 
        END DO

        CALL SORTINT(IBIN,IBINK,ILINE) 
        NRUN = 0
        DO K=1,ILINE
           PLIST(1) = K
           PLIST(2) = IBINK(K)
           PLIST(3) = IBIN(K)
           CALL SAVDN1(NDOC2,DOCF2,PLIST,3,NRUN,0)
           NRUN = 1
        END DO

        CLOSE(NDOC2)

C       NOW CREATE THE IMAGE:
        DO K=IYEND+1,IYSEC
           DO I=1,MAXXD
              PORTION(I,K)=0.
           END DO
        END DO
        
        DO K=1,IYSTRT
           DO I=1,MAXXD
              PORTION(I,K)=0.
           END DO
        END DO

C       START LOOP OVER Y STRIPS

        IBCOUNT =1
        DO 201 K=1,IY
C          FIRST LINE # FOR WRITING THE SECTION
           IYOFF=(K-1)*IYSEC
C          CLEAR ARRAY:
           DO KK=IYSTRT+1,IYSTRT+IYE
              DO I=1,MAXXD
                 PORTION(I,KK)=0 
              END DO
           END DO       
           IF (IBCOUNT.GT.ILINE) GOTO 217 
C          END CLEAR ARRAY

C          LOOP OVER X-DIRECTION
           ICOUNT = 0
           NLETI  = lnblnkn(FILPAT)

           DO 202 L=1,IX
203          CONTINUE

             ICURR = IBIN(IBCOUNT)
             IQ    = MOD(ICURR,IX)

C            FOR IQ=0 THIS IS THE LAST SECTION IN X
             IF  (IQ .EQ. 0) IQ=IX
             IQY = (ICURR-IQ) / IY+1

C            IF IQ .NE. L THEN WE ARE IN THE WRONG QUADRANT, 
C            GOTO NORMALIZE LAST SECTION

             !WRITE(NDAT,301) IBCOUNT,ICURR,IQ,IQY,IBINK(IBCOUNT),K,L
301          FORMAT(
     &          ' IBCOUNT,ICURR  ,IQ     ,IQY,IBINK   ,K,L'/7(3X,I5))

             IF (IQ  .NE. L) GOTO 207
             IF (IQY .NE. K) GOTO 207

C            WE HAVE THE CORRECT QUADRANT, ADD IMAGE GET FILE
             NUMBER = IBINK(IBCOUNT)    
             CALL FILGET(FILPAT,FILN1,NLETI,NUMBER,IRTFLG)

             !write(6,*) ' FILPAT:',FILPAT(1:30)
             !write(6,*) ' nleti,number:',nleti,number

             MAXIM = 0
             CALL OPFILEC(0,.FALSE.,FILN1,LUN,'O',IFORM,NSAM,NROW,NSLIC,
     &                    MAXIM,' ',.FALSE.,IRTFLG)

             IF (IMAMI .NE. 1) 
     &          CALL NORM3(LUN,NSAM,NROW,NSLIC,FMAX,FMIN,AV)

             IBCOUNT = IBCOUNT+1
             ICOUNT  = ICOUNT+1
             IYCOUNT = IYSTRT
             DO I=IYA,IYE
                IYCOUNT = IYCOUNT+1
                CALL REDLIN(LUN,BUF,NSAM,I)           
                ILC   = 0
                IXIND = (L-1)*IXSEC+IXSTRT
                DO KK=IXA,IXE
                  ILC = ILC+1
                  ILP = ILC+IXIND
                  PORTION(ILP,IYCOUNT)=PORTION(ILP,IYCOUNT)+(BUF(KK)-AV)
C                IF((I.EQ.IYA.AND.KK.EQ.IXA).OR.(I.EQ.IYE.AND.KK.EQ.IXE)) 
C     &          WRITE(NDAT,321) ILP,IYCOUNT,ICURR,IQ,KK,I,ILC,IXIND,IXSTRT
C321             FORMAT(' ',
C     &          'UPPER LEFT CORNER,ICURR,IQ,IYCOUNT,KK,I,ILC,IXIND,IXSTRT'/
C     &          9(1X,I5))
               END DO
            END DO

            CLOSE(LUN)
            IF (IBCOUNT.GT.ILINE) GOTO 207 
            GOTO 203

C          NORMALIZE
207        CONTINUE 

           !WRITE(NDAT,302) ICOUNT 
302        FORMAT('  SCALING WITH: 1/',I4)
 
           IF (ICOUNT .EQ. 0) GOTO 202
           IYCOUNT = IYSTRT
           PMAX    = -.999E20
           PMIN    = -PMAX
           DAVD    = 0
           DAVD2   = 0
           PNALL   =  (IYE-IYA+1)*(IXE-IXA+1)

           DO I=IYA,IYE
              IYCOUNT = IYCOUNT+1
              ILC     = 0
              IXIND   = (L-1)*IXSEC+IXSTRT

              DO KK=IXA,IXE
                 ILC = ILC + 1
                 ILP = ILC + IXIND

C                IF ((I.EQ.IYA.OR.I.EQ.IYE).AND.(KK.EQ.IXE.OR.KK.EQ.IXA)) 
C     &             WRITE(NDAT,322) ILP,IYCOUNT,IXIND,IXSEC,IXSTRT
C     &             ,ILC,IXA,IXE                      
322              FORMAT(X,'ILP,IYCOUNT,IXIND,IXSEC,IXSTRT,ILC,IXA,IXE'
     &              /8(1X,I5))  
                    
                 PORTION(ILP,IYCOUNT) = PORTION(ILP,IYCOUNT) / 
     &                                    FLOAT(ICOUNT)
                 B     = PORTION(ILP,IYCOUNT)
                 PMAX  = AMAX1(B,PMAX)
                 PMIN  = AMIN1(B,PMIN)
                 DAVD  = DAVD+B
                 DAVD2 = DAVD2+DBLE(B)*DBLE(B)
               END DO
           END DO

           PSIG = DSQRT((DAVD2-DAVD*DAVD/PNALL)/DBLE(PNALL-1.0))

           THUP   = RANGEU*PSIG
           THDOWN = -RANGEL*PSIG
           IF (PMAX .GT. THUP) PMAX=THUP
           IF (PMIN .LT. THDOWN) PMIN=THDOWN
           PNORMAL = PMAX-PMIN
           IYCOUNT = IYSTRT

           DO I=IYA,IYE
              IYCOUNT = IYCOUNT+1
              ILC     = 0
              IXIND   = (L-1)*IXSEC+IXSTRT

              DO KK=IXA,IXE
                 ILC = ILC+1
                 ILP = ILC+IXIND
                 B   = PORTION(ILP,IYCOUNT) 
                 IF (B.GT.THUP)   PORTION(ILP,IYCOUNT) = THUP  
                 IF (B.LT.THDOWN) PORTION(ILP,IYCOUNT) = THDOWN
                 PORTION(ILP,IYCOUNT) = PORTION(ILP,IYCOUNT)/PNORMAL 
              END DO
           END DO

C          -NORMALIZE END

179        ICOUNT = 0
           IF (IBCOUNT .GT. ILINE) GOTO 217 
202     CONTINUE

C       X-QUADRANTS FINISHED, WRITE ARRAY
217     CONTINUE
        !WRITE(NDAT,303) IYOFF+1,IYOFF+IYSEC
303     FORMAT('  WRITING FROM LINE: ',I0,' TO LINE ',I0) 

        DO I=1,IYSEC
           II=I+IYOFF
           CALL WRTLIN(LUNO,PORTION(1,I),NSAMO,II)
        ENDDO

201     CONTINUE

        CLOSE(LUNO)

        END



