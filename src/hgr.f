C ++********************************************************************
C                                                                      *
C   HGR                                                               *
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
C   HGR                                                              *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

       SUBROUTINE HGR(LUN51,IP,M,D,NG,NMAX,TMEAN,JG,
     &    AR,JV,VV,MXM,E,IHISTI,XT,LEST)

        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*2 (I-N)
	INTEGER*4 LUN51
       DIMENSION AR(NMAX),IHISTI(NMAX,NMAX),XT(M),E(NMAX)
       DIMENSION D(M),VV(MXM),IH(90,12),LBAR(12)
       DIMENSION TMEAN(M),JG(NMAX),JV(M)
       COMMON /HFGR/IH,LINE
       CHARACTER*4  NG(NMAX)
       CHARACTER*1  LINE(90),IG(12),NAM(12),LX,IKR

       CHARACTER*10 IBI
       LOGICAL*1    IFR
       INTEGER*4    LERC

       DATA IG/'1','2','3','4','5','6','7','8','9','0','A','B'/
       DATA NAM/'P','N','J','D','E','F','G','H','I','J','K','L'/
       DATA LX /' '/
       DATA IKR /'.'/

       MIM(X)=MAX0(MIN0(90,IDINT(X*15.+45.)),1)

       DO    J=1,NMAX
       DO    I=1,90
        IH(I,J)=0.0
       ENDDO
       ENDDO
       M1=M
       NSUM=0
         IF(LEST.LT.100) LERC=0
C                       CALL  WRTXT(
C     &'CLASSIFICATION     , RECORD NUMBER:',35,17,15,1)
         REWIND   4
 1      READ(4,9787,END=1000)  VV
9787	FORMAT(2X,F3.1,4(1X,1PE14.7),/,4(1X,1PE14.7),/,1(1X,1PE14.7))
         IF(LEST.LT.100)  THEN
         LERC=LERC+1
         READ(10,REC=LERC) IFR
         IF(IFR)  GOTO  1
         ENDIF
        KG=VV(IP)
       DO  I=1,NMAX
        IF(KG.EQ.JG(I)) GOTO 6
       ENDDO
       GOTO 1
 6     KG=I
       WRITE(IBI,7023) NSUM+1
 7023  FORMAT(I8)
C                       CALL  WRTXT(
C     & IBI,8,52,15,3)
       DO  I=1,M1
       J=JV(I)
       XT(I)=VV(J)
       ENDDO
       NSUM=NSUM+1
       X=0.0
       DO  I=1,M1
       Z=XT(I)-TMEAN(I)
       X=X+Z*D(I)
       ENDDO
       DS=1.D30
       DO  77  I=1,NMAX
       Y=(X-AR(I))**2-E(I)
       IF(Y.GE.DS)  GOTO  77
       DS=Y
       L=I
 77    CONTINUE
       IHISTI(L,KG)=IHISTI(L,KG)+1
       J=MIM(X)
       IH(J,L)=IH(J,L)+1
       GOTO 1
 1000  CONTINUE
       DS=1.D30
       X=-3.0
       DO  71  J=1,NMAX
       Y=(X-AR(J))**2-E(J)
       IF(Y.GE.DS)  GOTO 71
       DS=Y
       LA=J
 71    CONTINUE
       KL=0
       DO  72  I=2,90
       DS=1.D30
       X=FLOAT(I-45)/15.
       DO  73  J=1,NMAX
       Y=(X-AR(J))**2-E(J)
       IF(Y.GE.DS)  GOTO 73
       DS=Y
       LB=J
 73    CONTINUE
       IF(LA.EQ.LB)  GOTO  72
       LA=LB
       KL=KL+1
       LBAR(KL)=I
 72    CONTINUE
       KI=0
       DO    I=1,90
        DO    J=1,NMAX
	 ILH=IH(I,J)
         KI=MAX(ILH,KATC)
	 KATC=KI
	ENDDO
       ENDDO
       Z=KI
       WRITE(LUN51,102)
       WRITE(LUN51,801)
 801   FORMAT(' +',44('-'),'+',45('-'),'+')
       DO    KI=1,40
       IK=41-KI
       DO    J=1,90
        LINE(J)=LX
       ENDDO
       DO    J=1,KL
       JK=LBAR(J)
       LINE(JK)=IKR
       ENDDO
       DO    J=1,90
       IT=41
       LM=0
       DO    I=1,NMAX
       K=DBLE(IH(J,I))/Z*40.0
       IF(K.GE.IK .AND. K.LT.IT) THEN
       IT=K
       LM=I
       ENDIF
       ENDDO
       IF(LM.NE.0)  LINE(J)=NAM(LM)
       ENDDO
       WRITE(LUN51,802)  LINE
 802   FORMAT(' I',90A1,'I')
       ENDDO
       WRITE(LUN51,801)
       DO    I=1,90
        LINE(I)=LX
       ENDDO	
       DO    I=1,NMAX
       J=45.+AR(I)*15.
       LINE(J)=IG(I)
       ENDDO
       WRITE(LUN51,803)  LINE
 803   FORMAT('  ',90A1)
       WRITE(LUN51,102)
       DO    I=1,NMAX
        WRITE(LUN51,102)  IG(I),NAM(I),JG(I),NG(I)
       ENDDO
 102     FORMAT(3X,A1,' , ',A1,' = ',I4,2X,A4)
       RETURN
       END
C
