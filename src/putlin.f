
C++*********************************************************************
C
C PUTLIN.F                          NEW          OCT 94 ArDean Leith 
C                                   FILNAMANDEXT MAR 03 ArDean Leith
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
C  PUTLIN(LUN,LUNDOC,MAXDIM)
C
C  PURPOSE: PUTS LINES IN SPIDER IMAGE FROM DOC FILE. THIS SHOULD USE
C           BRESENHAMS ALGORITM FOR SPEED BUT I AM IN A HURRY SO I 
C           JUST BORROWED SOME EXISTING CODE!!!!! 
C            
C	  LUN 		LOGICAL UNIT NUMBER OF INPUT FILE
C	  LUNDOC	LOGICAL UNIT NUMBER OF DOCUMENT FILE
C	  MAXDIM        MAX. BUFFER SIZE
C
C--*********************************************************************
 
	SUBROUTINE PUTLIN(LUN,LUNDOC,MAXDIM)

	INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=MAXNAM)   ::   DOCNAM,FILNAM,FILTST
	COMMON /COMMUN/DOCNAM,FILNAM,FILTST

C       ACTUAL BUF LENGTH IS MAXDIM
	COMMON BUF(1)

        PARAMETER (MAXKEY=9999) 
        PARAMETER (MAXREG=7)    
        COMMON /DOC_BUF/ DBUF(MAXREG,MAXKEY*2)

	DIMENSION     PLIST(MAXREG+1)
        CHARACTER     NULL,DISP
        LOGICAL       EX,NEWCNT
        INTEGER       CNTNUM

        NULL = CHAR(0)

C       OPEN DOC FILE THAT CONTAINS COORDINATES, ~9 ALLOWS EXTENSION
12      CALL FILERD(DOCNAM,NLET,DATEXC,'COORDINATE DOCUMENT~9',IRTFLG)

        ICALL = 0
        CALL UNSDAL(DOCNAM,ICALL,LUNDOC,1,PLIST,1,DBUF,
     &              MAXKEY,MAXREG,NKEY,LER)
        IF (LER .GT. 0) GOTO 999

C       GET NAME FOR EXISTING OR NEW IMAGE FILE
10      CALL FILERD(FILNAM,NLET,NULL,'OUTPUT IMAGE',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

C       APPEND EXTENSION
        CALL FILNAMANDEXT(FILNAM,DATEXC,FILTST,NLET,.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

C       SEE IF THE IMAGE ALREADY EXISTS
        INQUIRE(FILE=FILTST,ERR=999,EXIST=EX)
        DISP = 'U'
        IF (EX) DISP = 'O'

C       OPEN IMAGE FILE
        IFORM  = 1
        NSAM   = 0
        NROW   = 0
        NSLICE = 0
        MAXIM  = 0
        CALL OPFILEC(0,.FALSE.,FILNAM,LUN,DISP,IFORM,NSAM,NROW,NSLICE,
     &             MAXIM,'IMAGE',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 10
        IF (NSAM * NROW .GT. MAXDIM) THEN
           CALL ERRT(6,'PUTLIN',NE)
           GOTO 999
        ENDIF


14      ICOLX = 1
        ICOLY = 2 
        CALL RDPRIS(ICOLX,ICOLY,NOT_USED,'X-COL., Y-COL.',IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 12
        ICOLX = ICOLX + 1
        ICOLY = ICOLY + 1
        IF ((ICOLX .LT. 0      .OR. ICOLY .LT. 0) .OR.
     &      (ICOLX .GT. MAXREG .OR. ICOLY .GT. MAXREG)) THEN
           CALL ERRT(101,'COLUMN OUT OF REGISTER RANGE',NE)
           GOTO 10
        ENDIF


C       ICOLI IS COLUMN OF DOC FILE CONTAINING INTENSITY
16      ICOLI = -55
        IF (IRTFLG .EQ. -1) GOTO 12
        CALL RDPRI1S(ICOLI,NOT_USED,
     &     'LINE INTENSITY COL. (< 0 ASKS FOR INTENSITY INPUT)',IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 14

17      IF (ICOLI .LT. 0) THEN
          CALL RDPRM1S(FOREGR,NOT_USED,'LINE INTENSITY',IRTFLG)
           IF (IRTFLG .EQ. -1) GOTO 16

           FOREGR2 = FOREGR
           ICOLI   = -1
        ELSE
           ICOLI = ICOLI + 1
           IF (ICOLI .LT. 0 .OR. ICOLI .GT. MAXREG) THEN
              CALL ERRT(101,'COLUMN OUT OF MAX. REGISTER RANGE',NE)
              GOTO 16
           ENDIF
        ENDIF

18      IF (.NOT. EX) THEN
           BACKGR = 0.0
           CALL RDPRM2S(BACKGR,NDUM,NOT_USED,'BACKGROUND INTENSITY',
     &                  IRTFLG)
          IF (IRTFLG .EQ. -1 .AND. ICOLI .EQ. -1) GOTO 17
          IF (IRTFLG .EQ. -1 ) GOTO 16
        ENDIF

C       ICOLCNT IS COLUMN OF DOC FILE CONTAINING CONTOUR NUMBER
181     ICOLCNT = 0
        CALL RDPRIS(ICOLCNT,IDUM,NOT_USED,
     &       'CONTOUR NO. COL. (= 0 FOR NONE)',IRTFLG)
        IF (IRTFLG .EQ. -1 .AND. .NOT. EX)      GOTO 18
        IF (IRTFLG .EQ. -1 .AND. ICOLI .EQ. -1) GOTO 17
        IF (IRTFLG .EQ. -1) GOTO 16
        ICOLCNT = ICOLCNT + 1
      
19      XFACT = 1.0
        YFACT = 1.0
        CALL RDPRM2S(XFACT,YFACT,NOT_USED,'X-FACTOR, Y-FACTOR',IRT)
        IF (IRT .EQ. -1) GOTO 181

20      XOFF = 0.0
        YOFF = 0.0
	CALL RDPRM2S(XOFF,YOFF,NOT_USED,'X-OFFSET, Y-OFFSET',IRT)
        IF (IRT .EQ. -1) GOTO 19

        IF (EX) THEN
C          FILL BUFFER WITH EXISTING IMAGE
           DO IREC = 1, NROW
              J = (IREC -1) * NSAM
              CALL REDLIN(LUN,BUF(J),NSAM,IREC)
           ENDDO
        ELSE
C          FILL BUFFER WITH BACKGROUND VALUE
           DO IVOX = 1, NSAM*NROW
              BUF(IVOX) = BACKGR
           ENDDO
        ENDIF

        LASTCNT = -1
        NEWCNT  = .TRUE.
     
	DO I=1,NKEY
C          GET COORDS FROM DOCUMENT FILE

           IX2 = DBUF(ICOLX,I) * XFACT + XOFF
           IY2 = DBUF(ICOLY,I) * YFACT + YOFF
           IF (IY2   .LT. 0) IY2 = NROW + IY2 + 1
           IF (ICOLI .GT. 0) FOREGR2 = DBUF(ICOLI,I)

           IF (ICOLCNT .GT. 1) THEN
C             READ CONTOUR NUMBER FROM DOC FILE
              CNTNUM = DBUF(ICOLCNT,I)
c*************************
c          write(6,*) 'cntnum,lastcnt:',cntnum,lastcnt,icolcnt,dbuf(4,i),
c     &                DBUF(5,I)
c**********************
              IF (CNTNUM .NE. LASTCNT) THEN
                 LASTCNT = CNTNUM
                 NEWCNT  = .TRUE.
              ENDIF
           ENDIF

           IF (NEWCNT) THEN
C             START A NEW CONTOUR

              IX1 = IX2
              IY1 = IY2
              IF (IY1 .LT. 0) IY1 = NROW - IY1 + 1
              IF (ICOLI .GT. 0) FOREGR = DBUF(ICOLI,11)

              IF ((IX1 .GT. NSAM .OR. IX1 .LE. 0) .OR.
     &            (IY1 .GT. NROW .OR. IY1 .LE. 0)) THEN
                  IONE = 1
                  WRITE(NOUT,91) IONE,IX1,IY1
              ENDIF
              NEWCNT = .FALSE.

           ELSEIF ((IX2 .GT. NSAM .OR. IX2 .LE. 0) .OR.
     &         (IY2 .GT. NROW .OR. IY2 .LE. 0)) THEN
              WRITE(NOUT,91) I,IX2,IY2
91            FORMAT('*** POINT NO.',I4,':(',I4,',',I4,
     &               ') OUT OF IMAGE LIMITS')

           ELSE
C             POINT IS WITHIN IMAGE
              IF (IY1 .EQ. IY2) THEN
C                HORIZONTAL LINE WOULD CAUSE DIVISION BY ZERO
                 ICON = (IY1 -1) * NSAM
                 IGO  = MIN(IX1,IX2)
                 IEND = MAX(IX1,IX2)
                 DO IX = IGO,IEND
                    BUF(ICON + IX) = FOREGR
                 ENDDO

              ELSE
                 FACT =  FLOAT(IX2-IX1) /  FLOAT(IY2-IY1)
                 FCON =  - FACT * IY1 + 0.5
                 IF (IX1 .GT. IX2) FCON = - FACT * IY1 - 0.5
                 IGO  = MIN(IY1,IY2)
                 IEND = MAX(IY1,IY2)

                 IXL = IX1
                 IF (IY2 .LT. IY1) IXL = IX2

                 DO IY = IGO,IEND
C                  FIND X VALUE FOR THIS Y COORDINATE
                   IX = IX1 + IFIX(FACT * FLOAT(IY) + FCON)

C                  SET BUFFER AT THIS LOCATION TO FOREGROUND
                   BUF((IY -1) * NSAM + IX) = FOREGR

                   IF (IABS(IX - IXL) .GT. 1) THEN
C                    MUST ADD IN INTERPOLATED X VALUES 
                     IGOX  = MIN(IX,IXL)
                     IENDX = MAX(IX,IXL)
                     IHALF = IGOX + (IENDX - IGOX) / 2

                     IYT   = IYL
                     IYEND = IY
                     IF (IX .LT. IXL) THEN
                         IYT   = IY
                         IYEND = IYL
                     ENDIF

                     DO IXT = IGOX,IENDX
                         BUF((IYT -1) * NSAM + IXT) = FOREGR
                        IF (IXT .EQ. IHALF) IYT = IYEND                       
                     ENDDO
                   ENDIF
                   IXL = IX
                   IYL = IY
                 ENDDO 
              ENDIF
           ENDIF
C********************************
c         WRITE(6,977) IX1,IY1,ix2,iy2
c977      FORMAT(' (',I4,',',I4,')-->(',I4,',',I4,')')
C*********************************

C          SET NEW STARTING POINT
           IX1    = IX2
           IY1    = IY2
           FOREGR = FOREGR2
        ENDDO

C       PLACE BUFFER BACK IN IMAGE FILE
        DO IREC = 1, NROW
           J = (IREC - 1) * NSAM
           CALL WRTLIN(LUN,BUF(J),NSAM,IREC)
        ENDDO

        CALL SETPRM(LUN,NSAM,NROW,0.,0.,0.,'U') 
          
999	CLOSE(LUN)
	RETURN
	END
