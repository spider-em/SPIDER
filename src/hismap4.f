
C++*********************************************************************
C
C HISMAP4  -- CREATED FEB 88  BY ARDEAN LEITH
C             ADAPTED FROM HISMAP.FOR       FEB 25 88 ARDEAN LEITH
C             CHANGED OUTPUT TO POSTSCRIPT  MAR    99 ARDEAN LEITH
C             USED LNBLNKN                  AUG    99 ARDEAN LEITH
C             NO ASK                        OCT    03 ARDEAN LEITH
C             KLIC OVERFLOW BUG             JAN    14 ARDEAN LEITH
C **********************************************************************
C * AUTHOR: A. LEITH                                                   *                                                *
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
C    HISMAP4(IDIM,NPTS,X,Y,ID,MOD,PEX)
C
C    PURPOSE:  PREPARES POSTSSCRIPT PLOT  
C
C    PARAMETERS:  IDIM     DIMENSION FOR ARRAYS                   (SENT)
C                 NPTS     NO. OF POINTS ON MAP                   (SENT)
C                 X,Y      COORDINATES OF POINTS                  (SENT)
C                 ID       ID OF POINTS                           (SENT)
C                 MOD      MODE FOR SYMBOLS OR LABELS             (SENT)
C                 PEX      STANDARD DEVIATIONS                    (SENT)
C
C    NOTE:
C      COORDINATES X(*) FOR HORIZONTAL AXIS JX, Y(*) FOR VERTICAL AXIS JY 
C      LABELS ARE IN ID(*), FORMAT A1 IF MOD=1, FORMAT A4 IF MOD=4  
C      POINTS AT MORE THAN PEX STANDARD DEVIATIONS ARE POSITIONED ON THE
C      EDGES OF THE GRAPH (SUBROUTINE EPUR4).      
C      WARNING:  X(*), Y(*), ID(NPTS+1) ARE DESTROYED UPON RETURN
C      GRAPH IS ABORTED IF MORE THAN 264 POINTS ARE ON THE EDGES
C
C   CALLED BY:  SGRAF  
C
C **********************************************************************

      SUBROUTINE HISMAP4(IDIM,NPTS,X,Y,ID,MOD,PEX)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=MAXNAM) :: FILPOS
      LOGICAL               :: ASK
      CHARACTER(LEN=1)      :: NULL = CHAR(0)
      CHARACTER(LEN=1)      :: CCHAR

      REAL                  :: X(IDIM), Y(IDIM)
      CHARACTER(LEN=7)      :: ID(IDIM)      
      CHARACTER(LEN=7)      :: CDATA,CID
      CHARACTER(LEN=2)      :: AXTYPE
      INTEGER               :: IDUM
      INTEGER, PARAMETER    :: LUNPOS  = 80

C     CHECK TO SEE THAT NUMBER OF POINTS IS NOT EXCESSIVE
      IF (NPTS > IDIM) THEN
         CALL ERRT(102,'HISMAP4; NPTS EXCEEDS ARRAY DIMENSION',IDIM)
         RETURN
      ENDIF         

C     GET NAME OF POSTSCRIPT FILE AND OPEN AS SEQUENTIAL FORMATTED
10    CALL OPAUXFILE(.TRUE.,FILPOS,'ps',LUNPOS,0,'N',
     &               'POSTSCRIPT OUTPUT',.TRUE.,IRTFLGT)
      IF (IRTFLGT .NE. 0) RETURN
      NLETP = LNBLNKN(FILPOS)

C     GET TEXT SIZE
11    ITSIZA = 10
      ITSIZD = 9
      CALL RDPRIS(ITSIZA,ITSIZD,NOT_USED,
     &   'TEXT SIZE FOR AXIS AND DATA (USE <CR> FOR DEFAULT = 10,9)',
     &   IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

      ASK    = (ITSIZA > 0)
      ITSIZA = ABS(ITSIZA)
      ITSIZD = ABS(ITSIZD)

C     FIND POINTS ON BOUNDARY OF MAP
      CALL EPUR4(IDIM,NPTS, X,Y,ID,MOD,PEX,KP,IDUM,IRTFLG,NDAT)
      IF (IRTFLG == 1)  GOTO 9999     

C     FIND MIN/MAX
      XMINT = MINVAL(X(1:NPTS))
      XMAXT = MAXVAL(X(1:NPTS))
      YMINT = MINVAL(Y(1:NPTS))
      YMAXT = MAXVAL(Y(1:NPTS))

      !CALL BORNS(NPTS,X,XMINT,XMAXT)  ! altered 2014 al
      !CALL BORNS(NPTS,Y,YMINT,YMAXT)  ! altered 2014 al             
                  
      CALL POSTRT(-LUNPOS)
      CALL POSCALE(LUNPOS,1.0,1.0,  -12.0,-7.0,  125.0,107.0)

C     ADD AXIS TO PLOT
      XORG   = 0.0
      YORG   = 0.0

      XEND   = 120.0
      YEND   = 100.0

      AXTYPE = 'XO'
21    IF (.NOT. ASK) IRTFLG = -9
      CALL POSAXIS(AXTYPE, XMINT,XMAXT, XORG,YORG, XEND,YEND,XFACTR,
     &            LUNPOS,IRTFLG)
      IF (IRTFLG == -1) GOTO 11

      AXTYPE = 'YO'
      IF (.NOT. ASK) IRTFLG = -9
      CALL POSAXIS(AXTYPE, YMINT,YMAXT, XORG,YORG, XEND,YEND,YFACTR,
     &            LUNPOS,IRTFLG)
      IF (IRTFLG == -1) GOTO 21
 
C     PLOT X = 0 ORIGIN LINE
      X1 = XORG
      Y1 = (0.0 - YMINT) * YFACTR
      X2 = XEND
      Y2 = Y1
      CALL POSEG(LUNPOS,X1,Y1,X2,Y2)

C     PLOT Y ORIGIN LINE
      X1 = (0.0 - XMINT) * XFACTR
      Y1 = YORG
      X2 = X1
      Y2 = YEND
      CALL POSEG(LUNPOS,X1,Y1,X2,Y2)

C     PUT FILENAME AT TOP
      XPOS   = 0.0
      YPOS   = 107.0

      ITANGL = 0
      ITSIZE = ITSIZA
      JUST   = 0
      CALL POTEX(LUNPOS,FILPOS,NLETP,XPOS,YPOS, ITSIZE,ITANGL,JUST)

C     SET PARAMETERS FOR TEXT CONTOURS
      ITSIZE = ITSIZD
      ITANGL = 0
      JUST   = 1

      DO  IPT = 1,NPTS 

        CDATA = ID(IPT)
        NCHAR = LNBLNKN(CDATA)

C       LOCATION OF ID ON MAP
        XPOS = (X(IPT) - XMINT) * XFACTR
        YPOS = (Y(IPT) - YMINT) * YFACTR
        CALL POTEX(LUNPOS,CDATA,NCHAR,XPOS,YPOS, ITSIZE,ITANGL,JUST)
      ENDDO

C     CLOSE THE POSTSCRIPT-FILE 
      CALL POEND(LUNPOS)

      WRITE(NOUT,*) ' GRAPH PLACED IN: ',FILPOS(1:NLETP)

9999  CLOSE(LUNPOS)

      END                                                                     
