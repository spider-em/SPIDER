
C++*********************************************************************
C
C PUTLIN.F            NEW                        OCT 1994 ArDean Leith 
C                     FILNAMANDEXT               MAR 2003 ArDean Leith
C                     NO MAXDIM                  JUN 2018 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2018  Health Research Inc.,                         *
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
C  PUTLIN(LUN,LUNDOC)
C
C  PURPOSE: PUTS LINES IN SPIDER IMAGE FROM DOC FILE. THIS SHOULD USE
C           BRESENHAMS ALGORITM FOR SPEED BUT I AM IN A HURRY SO I 
C           JUST BORROWED SOME EXISTING CODE!!!!! 
C 
C  PARAMETERS:           
C	    LUN         LOGICAL UNIT NUMBER OF INPUT FILE
C	    LUNDOC      LOGICAL UNIT NUMBER OF DOCUMENT FILE
C
C--*********************************************************************
 
	SUBROUTINE PUTLIN(LUN,LUNDOC)

	IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC'

        INTEGER               :: LUN,LUNDOC

        CHARACTER(LEN=MAXNAM) :: DOCNAM,FILNAM,FILTST
        CHARACTER(LEN=1)      :: NULL = CHAR(0)
        REAL, ALLOCATABLE     :: BUF(:)

        INTEGER, PARAMETER    :: MAXKEY=9999 
        INTEGER, PARAMETER    :: MAXREG=7
        REAL                  :: DBUF    
        COMMON /DOC_BUF/ DBUF(MAXREG,MAXKEY*2)

	REAL                  :: PLIST(MAXREG+1)

        CHARACTER(LEN=1)      :: DISP
        LOGICAL               :: EX,NEWCNT
        INTEGER               :: CNTNUM
        INTEGER               :: ICALL,LER,IRTFLG,NKEY,MAXIM
        INTEGER               :: nlet,nx,ny,nz,ne,mwant,icolx,icoly
        integer               :: not_used,icoli,ndum,icolcnt,idum,iy
        real                  :: foregr,foregr2,backgr,xfact,yfact
        real                  :: xoff,yoff,fact,fcon
        integer               :: igox,iendx,ihalf,ixt,ixl,iyt,iyl,iy1
        integer               :: iyend,ixj,ivox,lastcnt,i,ix2,iy2,ix1
        integer               :: ione,icon,igo,j,iend,ix

C       OPEN DOC FILE THAT CONTAINS COORDINATES, ~9 ALLOWS EXTENSION
        CALL FILERD(DOCNAM,NLET,DATEXC,'COORDINATE DOCUMENT~9',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        ICALL = 0
        CALL UNSDAL(DOCNAM,ICALL,LUNDOC,1,PLIST,1,DBUF,
     &              MAXKEY,MAXREG,NKEY,LER)
        IF (LER > 0) GOTO 999

C       GET NAME FOR EXISTING OR NEW IMAGE FILE
        CALL FILERD(FILNAM,NLET,NULL,'OUTPUT IMAGE',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

C       APPEND EXTENSION
        CALL FILNAMANDEXT(FILNAM,DATEXC,FILTST,NLET,.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

C       SEE IF THE IMAGE ALREADY EXISTS
        INQUIRE(FILE=FILTST,ERR=999,EXIST=EX)
        DISP = 'U'
        IF (EX) DISP = 'O'

C       OPEN IMAGE FILE
        IFORM = 1
        NX    = 0
        NY    = 0
        NZ    = 0
        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FILNAM,LUN,DISP,IFORM,NX,NY,NZ,
     &              MAXIM,'IMAGE',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        ALLOCATE(BUF(NX*NY), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MWANT = NX*NY
           CALL ERRT(46,'PUTLIN; BUF',MWANT)
           GOTO 999
        ENDIF

        ICOLX = 1
        ICOLY = 2 
        CALL RDPRIS(ICOLX,ICOLY,NOT_USED,'X-COL., Y-COL.',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        ICOLX = ICOLX + 1
        ICOLY = ICOLY + 1
        IF ((ICOLX < 0      .OR. ICOLY < 0) .OR.
     &      (ICOLX > MAXREG .OR. ICOLY > MAXREG)) THEN
           CALL ERRT(101,'COLUMN OUT OF REGISTER RANGE',NE)
           GOTO 999
        ENDIF


C       ICOLI IS COLUMN OF DOC FILE CONTAINING INTENSITY
        ICOLI = -55
        CALL RDPRI1S(ICOLI,NOT_USED,
     &     'LINE INTENSITY COL. (< 0 ASKS FOR INTENSITY INPUT)',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        IF (ICOLI < 0) THEN
          CALL RDPRM1S(FOREGR,NOT_USED,'LINE INTENSITY',IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 999

           FOREGR2 = FOREGR
           ICOLI   = -1
        ELSE
           ICOLI = ICOLI + 1
           IF (ICOLI < 0 .OR. ICOLI > MAXREG) THEN
              CALL ERRT(101,'COLUMN OUT OF MAX. REGISTER RANGE',NE)
              GOTO 999
           ENDIF
        ENDIF

       IF (.NOT. EX) THEN
           BACKGR = 0.0
           CALL RDPRM2S(BACKGR,NDUM,NOT_USED,'BACKGROUND INTENSITY',
     &                  IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 999 
        ENDIF

C       ICOLCNT IS COLUMN OF DOC FILE CONTAINING CONTOUR NUMBER
        ICOLCNT = 0
        CALL RDPRIS(ICOLCNT,IDUM,NOT_USED,
     &       'CONTOUR NO. COL. (= 0 FOR NONE)',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        ICOLCNT = ICOLCNT + 1
      
        XFACT = 1.0
        YFACT = 1.0
        CALL RDPRM2S(XFACT,YFACT,NOT_USED,'X-FACTOR, Y-FACTOR',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        XOFF = 0.0
        YOFF = 0.0
	CALL RDPRM2S(XOFF,YOFF,NOT_USED,'X-OFFSET, Y-OFFSET',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        IF (EX) THEN
C          FILL BUFFER WITH EXISTING IMAGE
           DO IREC = 1, NY
              J = (IREC -1) * NX
              CALL REDLIN(LUN,BUF(J),NX,IREC)
           ENDDO
        ELSE
C          FILL BUFFER WITH BACKGROUND VALUE
           DO IVOX = 1, NX*NY
              BUF(IVOX) = BACKGR
           ENDDO
           buf = backgr

        ENDIF

        LASTCNT = -1
        NEWCNT  = .TRUE.
     
	DO I=1,NKEY
C          GET COORDS FROM DOCUMENT FILE

           IX2 = DBUF(ICOLX,I) * XFACT + XOFF
           IY2 = DBUF(ICOLY,I) * YFACT + YOFF
           IF (IY2   < 0) IY2     = NY + IY2 + 1
           IF (ICOLI > 0) FOREGR2 = DBUF(ICOLI,I)

           !write(6,*) ' foregr:',foregr,foregr2

           IF (ICOLCNT > 1) THEN
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
              IF (IY1 < 0)   IY1    = NY - IY1 + 1
              IF (ICOLI > 0) FOREGR = DBUF(ICOLI,I)

              IF ((IX1 > NX .OR. IX1 <= 0) .OR.
     &            (IY1 > NY .OR. IY1 <= 0)) THEN
                  IONE = 1
                  WRITE(NOUT,91) IONE,IX1,IY1
              ENDIF
              NEWCNT = .FALSE.

           ELSEIF ((IX2 > NX .OR. IX2 <= 0) .OR.
     &             (IY2 > NY .OR. IY2 <= 0)) THEN
              WRITE(NOUT,91) I,IX2,IY2
91            FORMAT('*** POINT NO.',I4,':(',I4,',',I4,
     &               ') OUT OF IMAGE LIMITS')

           ELSE
C             POINT IS WITHIN IMAGE
              
              IF (IY1 == IY2) THEN
C                HORIZONTAL LINE WOULD CAUSE DIVISION BY ZERO
                 ICON = (IY1 -1) * NX
                 IGO  = MIN(IX1,IX2)
                 IEND = MAX(IX1,IX2)
                 DO IX = IGO,IEND
                    BUF(ICON + IX) = FOREGR
                    !write(6,*) ' flat:',ix,iy1,foregr,backgr
                 ENDDO

              ELSE
                 FACT =  FLOAT(IX2-IX1) /  FLOAT(IY2-IY1)
                 FCON =  - FACT * IY1 + 0.5
                 IF (IX1 > IX2) FCON = - FACT * IY1 - 0.5
                 IGO  = MIN(IY1,IY2)
                 IEND = MAX(IY1,IY2)

                 IXL = IX1
                 IF (IY2 < IY1) IXL = IX2

                 DO IY = IGO,IEND
C                  FIND X VALUE FOR THIS Y COORDINATE
                   IX = IX1 + IFIX(FACT * FLOAT(IY) + FCON)

C                  SET BUFFER AT THIS LOCATION TO FOREGROUND
                   BUF((IY -1) * NX + IX) = FOREGR
                   !write(6,*) ' buf:',ix,iy,foregr,backgr

                   IF (IABS(IX - IXL) > 1) THEN
C                    MUST ADD IN INTERPOLATED X VALUES 
                     IGOX  = MIN(IX,IXL)
                     IENDX = MAX(IX,IXL)
                     IHALF = IGOX + (IENDX - IGOX) / 2

                     IYT   = IYL
                     IYEND = IY
                     IF (IX < IXL) THEN
                         IYT   = IY
                         IYEND = IYL
                     ENDIF

                     DO IXT = IGOX,IENDX
                         BUF((IYT -1) * NX + IXT) = FOREGR
                         !write(6,*) ' buf:',ixt,iyt,foregr,backgr
                         IF (IXT == IHALF) IYT = IYEND                       
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

        !call chkmaxloc('bufval  ', buf,  nx*ny)


C       PLACE BUFFER BACK IN IMAGE FILE
       CALL WRTVOL(LUN,NX,NY,1,1,BUF,IRTFLG)

        CALL SETPRM(LUN,NX,NY,0.,0.,0.,'U') 
          
999	CLOSE(LUN)
	
	END
