C ++********************************************************************
C                                                                      *
C  FGR                                                                 *
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
C                                                                      *
C  FGR                                                                 *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

       SUBROUTINE FGR(LUN51,IP,M,D,NG,NMAX,TMEAN,JG,N
     &               ,MD,XX,AR,JV,VV,MXM,E,IHISTI,XT,MDT,LEST)

       IMPLICIT REAL*8 (A-H,O-Z)
       IMPLICIT INTEGER*2 (I-N)
       INTEGER*4  LUN51,LERC
       DIMENSION  XX(M,M),AR(NMAX,MD),N(NMAX),LIN(90)
       DIMENSION  IHISTI(NMAX,NMAX),XT(M)
       DIMENSION  E(NMAX)
       DIMENSION  D(M,2),MAP(30,42),VV(MXM)
       DIMENSION  TMEAN(M),JG(NMAX),JV(M)
       COMMON  /HFGR/  MAP,LINE
       CHARACTER*4  NG(NMAX)
       CHARACTER*1  LINE(90),IG(12),NAM(12),IX,JX,KX,LX,IKR,MX,NX
     &              ,IGWZ
       CHARACTER*10 IBI
       LOGICAL*1    IFR

       DATA IG/'1','2','3','4','5','6','7','8','9','0','A','B'/
       DATA NAM/'A','B','C','D','E','F','G','H','I','J','K','L'/
       DATA IX,JX,KX,LX /'+','I','-',' '/
       DATA IGWZ /'*'/
       DATA IKR /'.'/

       M1   = M
       NSUM = 0

       WRITE(51)M1,MD,NMAX
       WRITE(51)(TMEAN(I),I=1,M1)
       DO  J=1,M1
          WRITE(51)(D(J,K),K=1,MD)
       ENDDO
       DO  I=1,NMAX
          WRITE(51)(AR(I,K),K=1,MD)
       ENDDO
       WRITE(51)(E(I),I=1,NMAX)
       WRITE(51)(JV(I),I=1,M1)

       DO  I=1,30
          DO  J=1,42
             MAP(I,J)=0
          ENDDO
       ENDDO
       IF (LEST .LT.100) LERC=0

C      CALL  WRTXT( 'CLASSIFICATION, RECORD NUMBER:',35,17,15,1)

       REWIND   4
       NMR = 0
  1    READ(4,5687,END=1000)  VV
5687   FORMAT(2X,F3.1,4(1X,1PE14.7),/,4(1X,1PE14.7),/,1(1X,1PE14.7))
       NMR = NMR+1
       IF (LEST.LT.100)  THEN
           LERC = LERC + 1
           READ(10,REC=LERC) IFR
           IF (IFR)  GOTO  1
       ENDIF
       KG = VV(IP)
       DO I=1,NMAX
           IF(KG .EQ. JG(I)) GOTO 6
       ENDDO
       GOTO 1

 6     KG=I

       WRITE(IBI,7023) NSUM+1
 7023  FORMAT(I8)
C      CALL  WRTXT( IBI,8,52,15,3)

       DO  I=1,M1
          J     = JV(I)
          XT(I) = VV(J)
       ENDDO
       NSUM = NSUM+1
       X=0.0
       Y=0.0
       DO I=1,M1
          Z=XT(I)-TMEAN(I)
          X=X+Z*D(I,1)
          Y=Y+Z*D(I,2)
       ENDDO
       CALL RYS1(X,Y,KG,MAP)
       CALL DIST4(M1,NMAX,MDT,KG,XX,XT,TMEAN,
     &           AR,N,E,IHISTI,VV)
       GOTO 1

1000  CONTINUE
       CALL DIST1(NMAX,MD,AR,N,LIN,E)
       DO  I=1,42
          DO  J=1,90
            LINE(J)=LX
          ENDDO
          CALL DIST2(I,NMAX,MD,AR,N,LIN,LINE,E)
          NX=JX
          IF(I.NE.2.AND.I.NE.41) GOTO 7
          MX=KX
          NX=IX
          DO  J=1,90
             LINE(J)=MX
          ENDDO

 7        CONTINUE
          LINE(3)=NX
          LINE(88)=NX
          CALL RYS3(I,LINE,MAP,NAM)
          DO J=1,NMAX
             IF(N(J).LT.1) GOTO 600
             IY=21.-7.*AR(J,2)
             IF(IY.LT.1) IY=1
             IF(IY.GT.42) IY=42
             IF(IY.NE.I) GOTO 600
             IZ=45.+15.*AR(J,1)
             IF(IZ.GT.90) IZ=90
             IF(IZ.LT.1) IZ=1
             IF(IZ.EQ.90) IZ=89
             IF(LINE(IZ+1).EQ.LX) LINE(IZ+1)=LINE(IZ)
             IF(LINE(IZ+1) .EQ. LINE(IZ)) GOTO 601
             IF(LINE(IZ+1) .NE. IKR) LINE(IZ+1)=IGWZ

 601         CONTINUE
             LINE(IZ)=IG(J)

 600         CONTINUE
          ENDDO
          WRITE(LUN51,101)   LINE
 101      FORMAT(1X,90A1)
       ENDDO

       WRITE(LUN51,102)
       DO  I=1,NMAX
          WRITE(LUN51,102) IG(I),NAM(I),JG(I),NG(I)
       ENDDO
 102   FORMAT(3X,A1,' , ',A1,' = ',I4,2X,A4)

       RETURN
       END
