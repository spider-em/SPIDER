
C ++********************************************************************
C                                                                      *
C  WINAVE.F    INTERP. FOR NONINTEGER AVERAGING. 5/2/83 MR
C              LONG FILE NAMES                   JAN 89 ARDEAN LEITH
C              OPFILEC                           FEB 03 ARDEAN LEITH
C              SETPRMB PARAMETERS                MAY 09 ARDEAN LEITH
C                                                                      *
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
C WINAVE(LUN1,LUN2,NDOC,NSAM1,NROW1,NSAM2,NROW2)
C
C PARAMETERS:
C	  LUN1		LOGICAL UNIT NUMBER OF LARGE FILE
C	  LUN2		LOGICAL UNIT NUMBER OF SMALL FILE
C	  NDOC		LOGICAL UNIT NUMBER FOR DOCUMENT FILE
C	  NSAM1,NROW1	DIMENSIONS OF IMAGE FILE TO BE AVERAGED
C	  NSAM2,NROW2	DIMENSIONS OF SMALL FILE
C
C   INTERPOLATION FOR NONINTEGER AVERAGING INTRODUCED. 5/2/83 MR
C
C--*******************************************************************

	SUBROUTINE WINAVE(LUN1,LUN2,NDOC,NSAM1,NROW1,NSAM2,NROW2)

	INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 

        CHARACTER(LEN=MAXNAM)     ::  DOCFIL,FILNAM,FILPAT

	COMMON DUM(80),NPEAK(4000),PLIST(10),BUF(1)

        CHARACTER *1   NULL,FIX,Y,CONT,CDUM
 	LOGICAL        CFLAG,C2FLAG
	INTEGER        XBIG,YBIG,ILOW,IHI,IRTFLG

	DATA   Y/'Y'/,NBUF/4000/

        NULL = CHAR(0)
	LUN3 = 20

C       12/12/85 SEQUENTIAL SEARCH OPTION
	NSS = 0
	IF (FCHAR(4:4) .EQ. 'S') THEN
           NSS = 1
           WRITE(NOUT,1)
1	   FORMAT(' ** WARNING: SEQUENTIAL DOCUMENT SEARCH SPECIFIED')
        ENDIF

	NOPEN=0
	NSW=1
	NOFF1=31
	NOF11=NOFF1+NSAM1
	NOFF2=NOFF1+2*NSAM1
	NOFF3=NOFF2+NSAM2+1
	IB1=NOFF1
	IB2=NOF11

C       MAKE OUTPUT FILE BLANK
	DO  IS=1,NSAM2
	   BUF(NOFF1+IS-1) = 0.
	ENDDO
	DO  IR=1,NROW2
          CALL WRTLIN(LUN2,BUF(NOFF1),NSAM2,IR)
	ENDDO

	NSAMC = NSAM1 - NSAM2
	NROWC = NROW1 - NROW2
	CALL SETPRMB(LUN2 ,0.,0., 0.,0.)  ! UNNEEDED??

	CALL FILERD(DOCFIL,NLETD,NULL,'DOCUMENT',IRTFLG)
	IF (DOCFIL(1:1) .EQ. '*') RETURN

10 	CONTINUE
        IDONE = 1
        NUM   = NBUF
	ILOW  = 0
	IHI   = NBUF 
        CALL RDPRAI(NPEAK,NBUF,NUM,ILOW,IHI,'PEAK NUMBERS',CDUM,IRTFLG)

	CALL RDPRMC(FIX,NCHAR,.TRUE.,
     &          'NUMBER OF PEAKS FIXED? (Y/N)',NULL,IRTFLG)
	IFIX=0
	IF (FIX.EQ.'Y') IFIX=1

	CALL RDPRMI(ISTEP,IDUM,NOT_USED,'PEAK NUMBER INCREMENT')
	CALL RDPRMC(CONT,NCHAR,.TRUE.,'CONTROL WINDOWS? (Y/N)',
     &              NULL,IRTFLG)
	CFLAG = .FALSE.
	ICONT = 1

	IF (CONT .EQ. 'Y') THEN
	  CFLAG =.TRUE.
	  CALL RDPRMI(ICONT,IDUM,NOT_USED,'CONTROL INTERVAL')

C         GET FILE NAME PATTERN
          CALL FILSEQP(FILPAT,NLET,IDUM,0,IDUM,
     &                 'PREFIX OF CONTROL WINDOWS',IRTFLG)
          IF (IRTFLG .NE. 0) RETURN
        ENDIF

	IF (ISTEP.LE.0) ISTEP=1
	ICONT  = ICONT*ISTEP
        ICOUNT = 0
        ILOWX  = (NSAM1-NSAMC)/2+1
        ILOWY  = (NROW1-NROWC)/2+1
        IHIGHX = NSAM1/2+NSAMC/2
        IHIGHY = NROW1/2+NROWC/2
	WRITE (NOUT,1111) ILOWX,ILOWY,IHIGHX,IHIGHY
1111    FORMAT(' ILOWX,-Y,IHIGHX,-Y',4I6)

C	HANDMADE DOLOOP:
	I=1-ISTEP
120 	I=I+ISTEP
 
C       COORDS FROM DOCUMENT FILE
 
	ICK=I-1
	IDO=MOD(ICK,ICONT)
	C2FLAG = .FALSE.
	IF (CFLAG .AND. IDO.EQ.0 .OR. IDO.EQ.ISTEP) C2FLAG=.TRUE.
C       12/12/85: NSS SUBSTITUTED FOR 0 IN LAST ARGUMENT TO ENABLE 
C       SEQUENTIAL SEARCH OPTION IF SPECIFIED ('WV S')
        IKEY = NPEAK(I)

C       *****************88DEBUG
        WRITE(NOUT,*) 'DEBUGGING: I,IKEY,NPEAK(2):',I,IKEY,NPEAK(2)
C        *****************

	CALL UNSAV(DOCFIL,NOPEN,NDOC,IKEY,PLIST,2,LERR,NSS)
	NOPEN = 1
	IF (LERR .NE. 0) THEN
           WRITE (NDAT,9301) NPEAK(I)
9301       FORMAT(' ***',I4,' OUTSIDE PEAK FILE RANGE')
           GOTO 9500
        ENDIF

        X1=PLIST(1)
        Y1=PLIST(2)

        IF (X1.LT.ILOWX .OR. X1.GT.IHIGHX .OR. 
     &      Y1.LT.ILOWY .OR. Y1.GT.IHIGHY)THEN
           IF (IFIX.NE.0) THEN
             DO  III=NUM+1,NUM+ISTEP
               NPEAK(III) = NPEAK(III-1)+1
	     ENDDO
             NUM=NUM+ISTEP
           ENDIF

           IF (I.LT.NUM) GOTO 120
           IF (IDONE.EQ.1) GOTO 9500
           NSW = 2
           GOTO 10
        ENDIF

	IF (C2FLAG) THEN
          IWNUM = NPEAK(I)
          CALL FILGET(FILPAT,FILNAM,NLET,IWNUM,IRTFLG)

          MAXIM  = 0
          NSLICE = 1
          CALL OPFILEC(0,.FALSE.,FILNAM,LUN3,'U',IFORM,
     &             NSAM2,NROW2,NSLICE,MAXIM,' ',.TRUE.,IRTFLG)
        ENDIF

        IX1=X1
        IY1=Y1
        DX=1-(X1-IX1)
        DY=1-(Y1-IY1)
        C1=(1-DX)*(1-DY)
        C2=DX*(1-DY)
        C3=DY*(1-DX)
        C4=DX*DY

C	THE FOLLOWING LINES TAKE INTO ACCOUNT THE INTERPOLATION
	IXBIG=NSAM1+2-IX1-NSAM2/2-1
	IYBIG=NROW1+2-IY1-NROW2/2-1
        ICOUNT=ICOUNT+1
	WRITE (NDAT,11) NPEAK(I),X1,Y1,IXBIG,DX,IYBIG,DY
11	FORMAT(' PEAK #:',I4,' X:',F7.2,' Y:',F7.2,
     &     'UPPER LEFT CORNER: X=',I5,'+',F5.2,' Y=',I5,'+',F5.2)

	IF (IYBIG+NROW2-1.GT.NROW1 .OR. IXBIG+NSAM2-1.GT.NSAM1)THEN
           WRITE(NOUT,91)NPEAK(I)
91         FORMAT(' *** PEAK NO.',I3,' OUT OF LIMITS')

110        IF (IFIX.NE.0) THEN
             DO  III=NUM+1,NUM+ISTEP
                NPEAK(III)=NPEAK(III-1)+1
 	     ENDDO 
             NUM=NUM+ISTEP
           ENDIF

           IF (I.LT.NUM) GOTO 120
           IF (IDONE.EQ.1) GOTO 9500
           NSW=2
           GOTO 10
        ENDIF

	CALL REDLIN(LUN1,BUF(IB1),NSAM1,IYBIG)
	DO  J=1,NROW2
          LINE=IYBIG+J
          IF(LINE.GT.NROW1)LINE=NROW1
          CALL REDLIN(LUN1,BUF(IB2),NSAM1,LINE)
          CALL REDLIN(LUN2,BUF(NOFF2),NSAM2,J)
          DO  K=1,NSAM2
            IX2=IXBIG+K
            IF(IX2.GT.NSAM1)IX2=NSAM1
            I1=IB1-1+IXBIG+K-1
            I2=IB1-1+IX2
            I3=IB2-1+IXBIG+K-1
            I4=IB2-1+IX2
            VAL=BUF(I1)*C1+BUF(I2)*C2+BUF(I3)*C3+BUF(I4)*C4
            IF (C2FLAG) BUF(NOFF3+K-1)=VAL
            BUF(NOFF2+K-1)=BUF(NOFF2+K-1)+VAL
	  ENDDO
          IZ=IB1
          IB1=IB2
          IB2=IZ
          IF (C2FLAG) CALL WRTLIN(LUN3,BUF(NOFF3),NSAM2,J)
          CALL WRTLIN(LUN2,BUF(NOFF2),NSAM2,J)
	ENDDO

	IF (C2FLAG) CLOSE(LUN3)
        IF (I.LT.NUM)  GOTO 120
        IF (IDONE.EQ.1)  GOTO 9500
        NSW = 2
        GOTO 10

9500	CLOSE(NDOC)
	WRITE(NDAT,12) NSAMC,NROWC
12      FORMAT(' DIMENSIONS OF CCF WINDOW USED:',2I6)
        WRITE(NOUT,14) ICOUNT
14      FORMAT(' ',I5,' PEAKS USED FOR AVERAGING')

	RETURN
	END

