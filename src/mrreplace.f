
C ++********************************************************************
C                                                                      *
C    MRREPLACE                                                                  *
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
C  MRREPLACE(LUN1,LUN2)
C                                                                   *
C  PURPOSE:                                                            *
C    SUBROUTINE TO REPLACE THE 3D FOURIER TRANSFORM IN THE RANGE OF 
C    MEASURED DATA (+- ALPHA FOR SINGLE AXIS TILTING, OUTSIDE THE
C    MISSING CONE FOR CONICAL DATA).
C    PART OF MODULAR POCS IN SPIDER. 
C    THE PROGRAM USES THE SIMPLE 3D FOURIER FORMAT (-9)
C    M.R. DEC.91
C
C    CHANGED TO 'NEW' 3D FOURIER FORMAT, MAI 97, M.R.
C **********************************************************************

      SUBROUTINE MRREPLACE(LUN1,LUN2)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      PARAMETER       (LINDIM=2050)
      COMMON          FBUFW(LINDIM),FBUFD(LINDIM)

      CHARACTER(LEN=MAXNAM) :: FIL1,FIL2

      CHARACTER*1              NULL,SC
      LOGICAL                  SFLAG,CFLAG,NFLAG

      NULL = CHAR(0)

      PI=3.141593
      DTOR=1/180.*PI

      WRITE(NOUT,*) 'THIS FILE WILL BE OVERWRITTEN:'
      MAXIM = 0
      CALL OPFILEC(0,.TRUE.,FIL1,LUN1,'O',IFORM,NSAMSW,NROWW,NSLICEW,
     &               MAXIM,'WORK',.TRUE.,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF (NSAMSW.GT.LINDIM) THEN 
         CALL ERRT(102,'X-DIMENSION TOO LARGE, MAXIMUM IS',LINDIM)
         CLOSE(LUN1)
         RETURN
      ENDIF

      MAXIM = 0
      CALL OPFILEC(0,.TRUE.,FIL2,LUN2,'O',IFORM,NSAMD,NROWD,NSLICED,
     &               MAXIM,'DATA',.TRUE.,IRTFLG)
      IF (IRTFLG .NE. 0) THEN
         CLOSE(LUN1)
         RETURN
      ENDIF

      IF (NSAMD.NE.NSAMSW.OR.NROWD.NE.NROWW.OR.NSLICED.NE.NSLICEW) THEN
         CLOSE(LUN1)
         CLOSE(LUN2)
         CALL ERRT(1,'MRREPLACE',NDUM)
      RETURN

      ENDIF
      NA=1
110   CALL RDPRMC(SC,NA,.TRUE.,
     &    '(S)SINGLE AXIS OR (C)ONICAL, (N)TN REDUCTION',NULL,IRT)
      SFLAG=.FALSE.
      CFLAG=.FALSE.
      NFLAG=.FALSE.
      IF (SC.EQ.'S') SFLAG=.TRUE.
      IF (SC.EQ.'C'.OR.SC.EQ.'N') CFLAG=.TRUE.
      IF (SC.EQ.'N') NFLAG=.TRUE.
      IF(.NOT.SFLAG.AND..NOT.CFLAG) THEN
          WRITE(NOUT,*) 'TRY ANSWERING THIS QUESTION AGAIN, THANKS'
         GOTO 110
      ENDIF

      IF (NFLAG) THEN
          IF ((NSAMD-2).NE.NROWD) THEN
              CALL ERRT(101,'IMPLEMENTED ONLY FOR SQUARE PROJECTIONS',
     &                   NDUM)
              CLOSE(LUN1)
              CLOSE(LUN2)
              RETURN
          ENDIF
          CALL RDPRM(FACT,NOT_USED,
     &        'FACTOR FOR ZERO RANGE CALCULATION (DEF=0.875)')
          IF (FACT.LE.0) FACT=0.875
          FACTL=(FACT/PI)**2
      ENDIF
      IF (SFLAG) THEN 
C         CORRECTION: DALPHU IS LOWER TILT ANGLE, DALPHL IS HIGHER TILT
C             ANGLE. SIGN IS CHAGED TO MAKE IT CONSISTENT WITH 
C             WHAT COMES FROM THE MICROSCOPE.
         CALL RDPRM2(DALPHU,DALPHL,NOT_USED,
     &              'LOWER,UPPER TILT ANGLE')
         DALPHU=-DALPHU
         DALPHL=-DALPHL
      ENDIF

      IF (CFLAG) THEN 
         CALL RDPRM2(DALPHL,rmax,NOT_USED, 'TILT ANGLE')
         DALPHU=DALPHL 
      ENDIF

      ALPHL=DALPHL*DTOR
      ALPHU=DALPHU*DTOR
      COSALP2=COS(ALPHL)**2
      IREADO=0
      NS2=NSLICEW/2
      NR2=NROWW/2

C     DO THE REPLACE FOR SINGLE AXIS:
      BORDERU=1024
      BORDERL=1024
      SU=SIN(ALPHU)
      SL=SIN(ALPHL)
      IF(SU.NE.0) BORDERU=COS(ALPHU)/SU
      IF(SL.NE.0) BORDERL=COS(ALPHL)/SL
      IXA=1
      DO  L=1,NSLICEW
      LZ=L-NSLICEW/2-1
      IF (LZ.GE.0) XBORDER=ABS(LZ*BORDERU)
      IF (LZ.LT.0) XBORDER=ABS(LZ*BORDERL)
      RAD2=XBORDER*XBORDER
      IXB=XBORDER+1+.5
      IXB=2*IXB-1
C     LINEZ=(L-1)*NROWW
      LLZ=LZ
      IF (LZ.EQ.-NS2) LLZ=-LZ
      IF (NFLAG) THEN 
          RAD2=FACTL*LLZ*LLZ 
       ENDIF
      DO   K=1,NROWW
      KY=K-NROWW/2-1
      RY2=KY*KY
c     LINE=LINEZ+K
      KKY=KY
      IF (KKY.EQ.-NR2) KKY=-KKY

      LINE=LINE3DF(KKY,LLZ,NROWW,NSLICEW)
      CALL REDLIN(LUN1,FBUFW,NSAMSW,LINE)
      CALL REDLIN(LUN2,FBUFD,NSAMD,LINE)
      
      IF (SFLAG) IXA=IXB
      DO 13 I=IXA,NSAMSW,2
C      FOR CONICAL TILTS MAKE A CIRCLE:
       IF (CFLAG) THEN
         RX=INT(I/2.)
         RX2=RX*RX
         IF(NFLAG)RX2=RX2*COSALP2
         RAD=RY2+RX2
         IF(NSAMSW.LE.10) WRITE(NOUT,401) RX,RX2,RY2,RAD,RAD2
401      FORMAT(' RX,RX2,RY2,RAD,RAD2:',5F12.4)

         IF(RAD.LT.RAD2) GOTO 13
       ENDIF
      FBUFW(I)=FBUFD(I)
      FBUFW(I+1)=FBUFD(I+1)
13    CONTINUE

      IF (NSAMSW.LE.10) WRITE(NOUT,400) 
     &     KY,LZ,LINE,IXB,(FBUFW(LL),LL=1,NSAMSW)
400   FORMAT(1X,'KY:',I3,' LZ:',I3,' LI:',I3,' IXB:',I3,
     &     /' B:',10(1X,E10.4))
      CALL WRTLIN(LUN1,FBUFW,NSAMSW,LINE)
      ENDDO
      ENDDO
      CLOSE(LUN1)
      CLOSE(LUN2)

      RETURN
      END

C     -----------------------------------------------------------------


      INTEGER FUNCTION LINE3DF(KY,KZ,NROW,NSLICE)
      
C     FUNCTION TO FIND A LINE IN THE CURRENT (5/97) 3D FOURIER TRANSFORM
C     IN SPIDER
C     AUTHOR: M.RADERMACHER, 5/97 
      
      LINZ    = (KZ+NSLICE)
      LINZ2   = MOD(LINZ,NSLICE)
      LINZ3   = LINZ2*NROW
     
      LINE    = (KY + NROW) 
      LINE2   = MOD(LINE,NROW)
      LINE3   = LINZ3+LINE2+1
      LINE3DF = LINE3

      RETURN
      END

