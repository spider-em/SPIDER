C++*********************************************************************
C
C MRRSURF.F        FILENAMES LENGTHENED            JAN 89 al
C                  OPFILEC                         FEB 03 ARDEAN LEITH
C                  MAXNAM                          JUL 14 ARDEAN LEITH
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
C  MRRSURF
C
C  PURPOSE:
C      CREATE SURFACE REPRESENTATION OF THREE-DIMENSIONAL VOLUME
C
C  AUTHOR: M.RADERMACHER, DEZ.1982
C
C--*******************************************************************

      SUBROUTINE MRRSURF

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      COMMON ADUM(512),BUFA(512),BUF(512),BC(512,2),BOX(512),A(262144)

      CHARACTER(LEN=MAXNAM) :: FLN1,FLN3D

      CHARACTER *1  NULL,YZ
      INTEGER       LUN1,LUN2

      NULL = CHAR(0)
      LUN1 = 10
      LUN2 = 11

C----- READ INPUT -----------------------------------------------

      MAXIM = 0
      CALL OPFILEC(0,.TRUE.,FLN3D,LUN1,'O',IFORM,NSAM,NROW,NSLICE,
     &                   MAXIM,'3-D',.FALSE.,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      CALL RDPRMC(YZ,NCHAR,.TRUE.,
     &     'ROTATION AXIS (Y)(=DEF.)  OR (Z)',NULL,IRTFLG)
      IF (YZ.EQ.'Z') THEN 
         NROWP=NSLICE
         IYY=NROW
      ELSE
         NROWP=NROW
         IYY=NSLICE
      ENDIF

      CALL FILERD(FLN1,NLET,NULL,'OUTPUT',IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      CALL RDPRMI(NSAMP,IDUM,NOT_USED,
     &    'HORIZONTAL DIMENSION OF SURFACE DEPTH BUFFER')

      MAXIM = 0
      IFORM = 1
      CALL OPFILEC(LUN1,.FALSE.,FLN1,LUN2,'U',IFORM,NSAMP,NROWP,1,
     &                   MAXIM,'3-D',.FALSE.,IRTFLG)
      IF (IRTFLG .NE. 0) THEN
         CLOSE(LUN1)
         RETURN
      ENDIF

      CALL RDPRMI(NS1,NS2,NOT_USED,'DEPTH FROM, TO')
      NDEEP=NS2-NS1
      CALL RDPRM(PHI,NOT_USED,'VIEWING ANGLE')
      PHI=PHI/180.*3.14159265
      CALL RDPRM(SW,NOT_USED,'THRESHOLD')
      WRITE(NOUT,501) SW
501   FORMAT(' PRECISE THRESHOLD:',G10.4)
      FS=SIGN(1.,SW)
      SW=ABS(SW)
      CALL RDPRMI(IBACK,IDUM,NOT_USED,
     & 'BACKGROUND (0)ZERO,(1)MIN,(2)MAX,(3)LOCAL MIN,(4)LOCAL MAX')
      IF (IBACK.EQ.3.OR.IBACK.EQ.4) THEN 
          CALL RDPRMI(LBOX,IDUM,NOT_USED,
     &         'BOX LENGTH FOR BACKGROUND FILTER')
          CALL RDPRM(BOFF,NOT_USED,'BACKGROUND OFFSET')
      END IF
70    CONTINUE

C----- DETERMINE START VALUES ------------------------------------------
      DO  K=1,NSAMP
         BUF(K)=0.
      ENDDO
      BMIN=1000.
      BMINT=1000.
      BMAX=0.
      BMAXT=0.
      
C     CLEAR BUFFER FOR BACKGROUND CORRECTION:
      DO  L=1,NROWP
        BC(L,1)=0.
62      BC(L,2)=0.
      END DO
      
C     READ SLICE:
      IF (YZ.EQ.'Z') THEN 
        NWHAT=NSLICE
      ELSE
        NWHAT=NROW
      ENDIF
      DO  LAUF=1,NWHAT
         IF (YZ.NE.'Z') THEN
           DO  I=1,NSLICE
             N0=(I-1)*NROW+LAUF
             CALL REDLIN(LUN1,BUFA,NSAM,N0)
             IND0=(I-1)*NSAM
             DO  K=1,NSAM
               IND=IND0+K
               A(IND)=BUFA(K)
	     ENDDO
	   ENDDO
         ENDIF
         IF (YZ.EQ.'Z') THEN
            DO  I=1,NROW
               N0=(LAUF-1)*NROW+I
               IND=(I-1)*NSAM+1
               CALL REDLIN(LUN1,ADUM,NSAM,N0)
               DO  K=1,NSAM
                  A(IND-1+K)=ADUM(NSAM-K+1)
	       ENDDO
	    ENDDO
         ENDIF
         CPHI=COS(PHI)
         SPHI=SIN(PHI)
         DO 5 K=1,NSAMP
C        KSP=NSAMP+1-K
C       CHANGE SIDEDNESS OF IMAGE:
        KSP=K
        IF1=1
        IF2=1
        IF3=1
        IF4=1
        XP=K-1-NSAMP/2
        DO 6 L=1,NDEEP
          YP=L-1.+NS1
          X0=XP*CPHI-YP*SPHI+NSAM/2+1.
          Y0=XP*SPHI+YP*CPHI+IYY/2+1
          IF(X0.GE.(NSAM-1).OR.X0.LE.2) GOTO 6
          IF(Y0.GE.(IYY-1).OR.Y0.LE.2) GOTO 6
          IX0=INT(X0)
          IY0=INT(Y0)
          IND1=IX0+(IY0-1)*NSAM
          IND2=IND1+1
          IND3=IND1+NSAM
          IND4=IND3+1
          IF(A(IND1)*FS.GE.SW) IF1=0
          IF(A(IND2)*FS.GE.SW) IF2=0
          IF(A(IND3)*FS.GE.SW) IF3=0
          IF(A(IND4)*FS.GE.SW) IF4=0
          IG=IF1*IF2*IF3*IF4
          IF(IG.EQ.1) GOTO 6
C         WRITE(NOUT,100)IND1,IND2,IND3,IND4,A(IND1),A(IND2)
C    $        ,A(IND3),A(IND4)
100       FORMAT(' VALUE FOUND',4I4,4F12.4)
C-----------------------------------------------------------------------
C         INTERPOLATION
C-----------------------------------------------------------------------
          AI1=A(IND1)+(A(IND2)-A(IND1))*(X0-IX0)
          AI2=A(IND3)+(A(IND4)-A(IND3))*(X0-IX0)
          AI3=(AI1+(AI2-AI1)*(Y0-IY0))*FS
          IF(AI3.LT.SW) GOTO 6
          X2=XP*CPHI-(YP-1)*SPHI+NSAM/2+1
          Y2=XP*SPHI+(YP-1)*CPHI+IYY/2+1
          IX2=INT(X2)
          IY2=INT(Y2)
          IND1=IX2+(IY2-1)*NSAM
          IND2=IND1+1
          IND3=IND1+NSAM
          IND4=IND3+1
C-----------------------------------------------------------------------
C         CALCULATE VALUES OF NEXT POINTS
C-----------------------------------------------------------------------
          AJ1=A(IND1)+(A(IND2)-A(IND1))*(X2-IX2)
          AJ2=A(IND3)+(A(IND4)-A(IND3))*(X2-IX2)
          AJ3=(AJ1+(AJ2-AJ1)*(Y2-IY2))*FS
          IF(AI3.LE.AJ3) GOTO 7
          DIF=(SW-AJ3)/(AI3-AJ3)
          IF(DIF.LT.0) DIF=ABS((SW-AJ3)/(AI3+AJ3))
          BUF(KSP)=L+DIF-1
C     WRITE(NDAT,101) K,BUF(KSP),DIF,XP,YP,X0,Y0
101   FORMAT(' K=',I4,'BUF(K)=',F12.4,'DIF=',E12.4,
     &'XP,YP ',2F14.4,' X0,Y0 ',2F14.4)
          GOTO 5
7         BUF(KSP)=L
          GOTO 5
6       CONTINUE
5       CONTINUE
        CALL WRTLIN(LUN2,BUF,NSAMP,LAUF)
        DO  K=1,NSAMP
        IF(BMAX.LT.BUF(K)) BMAX=BUF(K)
        IF(BMIN.GT.BUF(K).AND.BUF(K).GT.0.) BMIN=BUF(K)
        BUF(K)=0.
	ENDDO
        BC(LAUF,1)=BMIN-BOFF
        BC(LAUF,2)=BMAX+BOFF
        IF(BMAXT.LT.BMAX) BMAXT=BMAX
        IF(BMINT.GT.BMIN) BMINT=BMIN
        BMAX=0.
        BMIN=1000
      ENDDO
C----- CALCULATE BACKGROUND --------------------------------------------

      IF(IBACK.EQ.0) GOTO 17



C----- CALCULATE BACKGROUND --------------------------------------------

      IF(IBACK.NE.1) GOTO 24

C----- MINIMUM TOTAL ---------------------------------------------------

      BMIN=BMINT-1
24    CONTINUE
      IF(IBACK.NE.2) GOTO 25

C----- MAXIMUM TOTAL ---------------------------------------------------

      BMAX=BMAXT+1.
25    CONTINUE
      IF (IBACK.NE.3.AND.IBACK.NE.4) GOTO 26

C----- LOCAL MINIMUM AND MAXIMUM ---------------------------------------

      IF(IBACK.EQ.3)J=1
      IF(IBACK.EQ.4)J=2
      DO  K=1,NROW
      BZ=BC(K,J)
      IZ=K
      IF(BZ.GT.0.) GOTO 32
      ENDDO
32    BC(1,J)=BZ
      DO  K=2,NROW
      IF(BC(K,J).LE.0.) BC(K,J)=BC(K-1,J)
      ENDDO
      WRITE (NOUT,999) J,(BC(K,J),K=1,NROW)
999   FORMAT(' BC',I3,1X,/,25(/1X,10F10.2))

C----- BOX CONVOLUTION -------------------------------------------------

      LEN=(LBOX-0.5)/2
      DO  K=1,NROW
      BOX(K)=BC(K,J)
      DO  L=1,LEN
      IND1=K-L
      IND2=K+L
      IF(IND1.LT.1)IND1=1
      IF(IND2.GT.NROW)IND2=NROW
      BOX(K)=BOX(K)+BC(IND1,J)+BC(IND2,J)
      ENDDO
      ENDDO
      DO  K=1,NROW
       BC(K,J)=BOX(K)/(2*LEN+1)
      ENDDO
      WRITE(NOUT,999) J,(BC(K,J),K=1,NROW)

C-----------------------------------------------------------------------

26    CONTINUE
      DO  K=1,NROW
        CALL REDLIN(LUN2,BUF,NSAMP,K)
        IF(IBACK.EQ.1) THEN

C----- MINIMUM TOTAL ----------------------------------------------
           DO  L=1,NSAMP
             IF(BUF(L).LE.0.0) BUF(L)=BMINT
	   ENDDO
        END IF
        IF(IBACK.EQ.2) THEN

C----- MAXIMUM TOTAL ----------------------------------------------
          DO  L=1,NSAMP
            IF(BUF(L).LE.0.0) BUF(L)=BMAXT
	  ENDDO
	END IF
        IF(IBACK.EQ.3.OR.IBACK.EQ.4) THEN

C----- LOCAL MINIMUM AND MAXIMUM ---------------------------------
          DO  L=1,NSAMP
             IF(BUF(L).GT.0.) BUF(L)=BUF(L)-BC(K,J)
	  ENDDO
        END IF
        CALL WRTLIN(LUN2,BUF,NSAMP,K)
      ENDDO
17    CONTINUE

C-----------------------------------------------------------------------

C     WRITE(6,200) NL1,NL2,NS1,NS2,SW,NSAM,NROW,NSLICE,CS,CQ
200   FORMAT(' NL1,NL2,NS1,NS2:',4I5,' SW:',F10.2,
     &       'NSAM,NROW,NSLICE:',3I5,'CS,CQ:',F10.2)
      CLOSE(LUN1)
      CLOSE(LUN2)
      RETURN
      END
