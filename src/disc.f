C ++********************************************************************
C                                                                      *
C DISC                                                                 *
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

       SUBROUTINE DISC(LUN50,LUN51)

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       IMPLICIT INTEGER*2 (I-N)
       PARAMETER  (L40=3,L50=10)
       INTEGER*4  LUN50,LUN51,N1,N2,N10,N20
       DIMENSION ALL(L50*L50*L40),T(L50*L50),W(L50*L50),A(L50*L50)
       DIMENSION  AA(L50,L50),BB(L50,L50),TMEAN(L50)
       DIMENSION JG(3),N(L40)
       DIMENSION JV(L50),XMIN(L50*L40),XMAX(L50*L40),XMEAN(L50*L40)
       DIMENSION VV(L50),XQ(L50)
       DIMENSION EE(L40),IHISTI(L40,L40),XT(L50),IQ(L50)
       CHARACTER*4  NG(3),K1,NV(10),NTMP,NVT
C      CHARACTER*1  NWZ(4),NVW(4)
       CHARACTER*10 IBI
c      CHARACTER *80 IFOR
C      EQUIVALENCE (NTMP,NWZ),(NVT,NVW)
       CHARACTER*1  NULL


        DATA UNIV /32000.0D0/
C	DATA  NWZ/'G','R','0','0'/,NVW/' ','V','0','0'/
c       DATA NSP/'    '/
	DATA NV(1)/'GRP '/,NV(2)/'VAR '/,NV(3)/'SKEW'/,NV(4)/'KURT'/
     &  ,NV(5)/'ENTP'/,NV(6)/'AVAV'/,NV(7)/'AVVR'/,NV(8)/'SDAV'/
     &  ,NV(9)/'SDVR'/,NV(10)/'VV  '/
	DATA NG(1)/'PRT '/,NG(2)/'NSE '/,NG(3)/'JNK '/
	DATA JG/1,2,3/

       NULL=CHAR(0)
       L1000=L50*L50
       L200=L50*L40
       IGR=L1000*L40
       DO  I=1,IGR
          ALL(I)=0.0
       ENDDO
       DO I=1,L1000
          T(I)=0.0
          W(I)=0.0
       ENDDO
       DO  I=1,L200
          XMIN(I)=1.E30
          XMAX(I)=-1.E30
          XMEAN(I)=0.0
       ENDDO 
       DO  I=1,L50
          JV(I)=I
          TMEAN(I)=0.0
       ENDDO
       DO  I=1,L40
          N(I)=0
       ENDDO

       WRITE(LUN51,1)
1      FORMAT(////' MANOVA  &  DISCRIMINANT  ANALYSIS       ' /)
       LEST=100
       MGR=3
       M=10
       M1=M-1

       KG=1
       IP=KG
       K0=JV(KG)
       JV(KG)=JV(M)
       JV(M)=K0
       K1=NV(KG)
       NV(KG)=NV(M)
       NV(M)=K1
	
       MDMAX=L40
       MD=M1
       IF (MD.GE.MGR) MD=MGR-1
       IF (MD.GT.MDMAX) MD=MDMAX
 9202  DO   K=1,MGR
          EE(K)=0.0
       ENDDO

       CALL MYREAD(M,MGR,KG,ALL,T,XMEAN,TMEAN,XMIN,XMAX,N,JG,LEST)
       MGR=0
       NSUM=0
       NMAX=0
       FA1S=0.0
       GA1S=0.0
       DO 11 K=1,L40
          I=N(K)
          IF(I.EQ.0) GOTO 11
          NMAX=K
          MGR=MGR+1
          NSUM=NSUM+I
          IF(I.EQ.1) GOTO 11
          X=I-1
          FA1S=FA1S+1.0/X
          GA1S=GA1S+1.0/(X*X)
 11    CONTINUE
       MTEMP = MGR - 1
       MDT=MIN(MTEMP,M1)

       WRITE(LUN51,97) M1,MGR
 97    FORMAT(//,1X,'ANALYSIS FOR',I4,' VARIABLES AND',I4,' GROUPS' )
 98    FORMAT(//,1X,26('-'),'MANOVA.ANALYSIS OF VARIANCE' ,26('-'))
 100   FORMAT(//,1X,'VARIABLE-CRITERION' ,5X,I4,1X,A4)
 99    FORMAT(//,1X,40('- '))
       WRITE(LUN51,98)
C.......................................................................
       DETL=0.0
C      CALL  WRTXT(
C     &   'CALCULATION OF COVARIANCE MATRIX FOR GROUP:',431,17,6,1)
       DO  K=1,NMAX
          WRITE(IBI,7023) K
 7023     FORMAT(I4)
C         CALL  WRTXT(IBI,4,58,6,3)
          J=N(K)
          WRITE(LUN51,100) K0,K1
          IF(J.EQ.0) GOTO 130
          WRITE(LUN51,96) JG(K),NG(K),J
 96       FORMAT(' GROUP: ',I4,1X,A4,10X,I5,' SUBJECTS' /)
          CALL MYWR0(LUN51,K,M,J,DET,JV,NV,XMEAN,ALL,XMIN,XMAX,W)
          IF(DET.LE.0.0D0)  THEN
             WRITE(LUN51,3071)
 3071       FORMAT(' ',/,'  DETERMINANT NEGATIVE - ANALYSIS ABORTED !!',
     *  //,'   TOO FEW DATA OR VARIABLES ARE STRONGLY CORRELATED')
             GOTO  3075
          ENDIF
          DETL=DETL+ ((J-1)*DLOG(DET))
          GOTO 132
 130      WRITE(LUN51,131) JG(K),NG(K)
 131      FORMAT(' GROUP: ',I4,1X,A4,13X,'NO SUBJECTS' /)
 132      CONTINUE
C         WRITE(LUN51,99)
       ENDDO
       WRITE(LUN51,100) K0,K1
       WRITE(LUN51,95)
 95    FORMAT(/' MEANS FOR TOTAL SAMPLE  ',
     &         'AND  POOLED-SAMPLES STANDARD DEVIATION  ' )
       CALL MYWR1(LUN51,M,MGR,NMAX,JV,NV,TMEAN,ALL,W,NSUM)
C      WRITE(LUN51,99)
C      WRITE(LUN51,94)
C 94   FORMAT(//19H CORRELATION MATRIX,/)
       CALL  COREL(T,AA,BB,M,M1,NSUM)
       CALL MTPR(M,M1)
       WRITE(LUN51,99)
C      WRITE(LUN51,93)
C 93   FORMAT(//37H SIGNIFICATIONS OF CORRELATION COEFF. /)
       CALL MTPR(M,M1)
       JV(M)=K0
C      WRITE(LUN51,99)
C      WRITE(LUN51,206)
C 206  FORMAT(//' MATRICES OF MAHALANOBIS  SQUARE-DISTANCES',
C     & ' AND THEIR SIGNIFICANCES'/)
C                       CALL  WRTXT(
C     &'CALCULATION OF MAHALANOBIS DISTANCES',36,23,9,1)
       MGR1=NMAX+1
       CALL MAHAL(M,MGR1,W,A,XMAX,XMEAN,N,XMIN)
       CALL MTPR(MGR1,NMAX)
       CALL MTPR(MGR1,NMAX)
C      WRITE(LUN51,99)
 207   CONTINUE
C      WRITE(LUN51,92)
C 92   FORMAT(//' TESTS FOR THE GIVEN STRUCTURE OF MEAN VALUE'/,
C     & '    NULL HYPOTHESIS H0: MEAN(GR,VAR)=MEAN(VAR) '/)
       CALL MYWR2(MGR,M,NSUM,W,T,A)
C      WRITE(LUN51,99)
       N01=1
       N02=M1
       L=1
       DO   I=1,NMAX
         DO    K=N01,N02
            XMEAN(L)=XMEAN(K)
            L=L+1
	 ENDDO
         N01=N02+2
         N02=N02+M
       ENDDO
       CALL  PRP(AA,BB,W,T,M,M1)
       CALL MTNV(W,M,DETW)
       CALL MTNV(T,M,DETT)
       X=NSUM-MGR
       XMM=X*(-M1*DLOG(X)+DLOG(DETW))
       XMM=XMM-DETL
       EK=MGR
       EN=NSUM
       EM=M1
       EK1=MGR-1
       EM1=M1-1
       F1=0.5*EK1*EM*(EM+1.)
       A1A=(FA1S-1.0/X)*(2.*(EM*EM)+3.0*EM-1.0)
       A1=A1A/(6.0*EK1*(EM+1.0))
       A2=(GA1S-1.0/(X*X))*EM1*(EM+2.0)/(6.0*EK1)
       DIF=A2-A1*A1
       IF(DIF.GT.0.0) GOTO 25
       F2=-(F1+2.0)/DIF
       B1=F2/(1.0-A1+(2.0/F2))
       F=(F2*XMM)/(F1*(B1-XMM))
       GOTO 45
 25    F2=(F1+2.0)/DIF
       B1=F1/(1.0-A1-(F1/F2))
       F=XMM/B1
 45    CONTINUE
       IF(F1.GT.UNIV) F1=UNIV
       IF(F2.GT.UNIV) F2=UNIV
       N1=F1
       N2=F2
       ALPH=ALPHAINT(F,N1,N2)
C      WRITE(LUN51,46)
C 46   FORMAT(//' TEST FOR THE HYPOTHESIS OF THE EQUALITY OF GROUP',
C     & ' DISPERSION MATRICES'/)
C      WRITE(LUN51,47) XMM,N1,N2,F,ALPH
C 47    FORMAT(//32H TESTED BY BARTLETT-BOX TEST M = ,G14.4 /
C     & 30H F-TEST APPROXIMATION WITH DF1,I6,8H AND DF2,I8/4H F =,G14.4,
C     & 5X,6HSIGN =,F7.3)
       CALL MULT(M,TRACE,A,W)
       DIV=EN-EK1*EM-2.0
       IF(DIV.LE.0.0) GOTO 200
       F1=(EK1*EM*(X-EM))/DIV
       F2=X-EM1
       GOTO 201
 200   F1=UNIV
       F2=X-EM1
       F2=F2-(F2-2.)*(F2-4.)*DIV/((X-1.)*(EN-EM1-1.))
 201   F0=(X-EM-1.0)*F2*TRACE/(EM*EK1*(F2-2.0))
       N10=F1
       N20=F2
 202  FORMAT(//' TESTED BY HOTELLING"S TRACE T =',G14.4/
     &' F-TEST APPROXIMATION WITH DF1',I6,' AND DF2',I8/' F =',G14.4,
     &5X,'SIGN =',F7.3)
C       WRITE(LUN51,99)
C       WRITE(LUN51,210)
C 210   FORMAT(//' TESTS FOR THE GIVEN STRUCTURE OF DISPERSION VALUES'/,
C     & '      NULL HYPOTHESIS H0: DISPERSION(GR,VAR)=DISPERSION(VAR)'/)
       CALL DISTEST(M,NSUM,MGR,NMAX,FA1S,ALL,N)
C       WRITE(LUN51,99)
       XL=DETW/DETT
C       WRITE(LUN51,49) XL
C 49     FORMAT(//' TEST FOR THE HYPOTHESIS OF THE EQUALITY OF GROUP ',
C     & 'MEAN VECTORS'///' TESTED BY WILKS  LAMBDA TEST W =',G14.4)
       IF(M1.GT.2.AND.MGR.GT.3) GOTO 50
       IF(M1.EQ.1.OR.MGR.EQ.2) GOTO 5001
       YL=SQRT(XL)
       F1=2*M1
       IF(M1.EQ.2) F1=2*MGR-2
       F2=2.0*EN-F1-4.0
       GOTO 51
 5001  YL=XL
       F1=M1
       IF(M1.EQ.1) F1=MGR-1
       F2=EN-F1-1.0
        GOTO 51
 50    SL=SQRT(((EM*EM)*(EK1*EK1)-4.0)/((EM*EM)+(EK1*EK1)-5.0))
       YL=XL**(1.0/SL)
       PL=(EN-1.)-((EM+EK)/2.0)
       QL=-((EM*EK1)-2.0)/2.0
       F1=EM*EK1
       F2=(PL*SL)+QL
  51    IF(F1.GT.UNIV) F1=UNIV
        IF(F2.GT.UNIV) F2=UNIV
        N1=F1
       N2=F2
       F=((1.0-YL)/YL)*(F2/F1)
       ALPH=ALPHAINT(F,N1,N2)
C      WRITE(LUN51,52) N1,N2,F,ALPH
C 52   FORMAT(30H F-TEST APPROXIMATION WITH DF1,I6,8H AND DF2,I8/4H F =,
C     &G14.4,5X,6HSIGN =,F7.3)
       ALPH=ALPHAINT(F0,N10,N20)
C      WRITE(LUN51,202) TRACE,N10,N20,F0,ALPH
C      WRITE(LUN51,98)
C
       CALL DISCRIM(LUN51,M,NMAX,NSUM,MGR,IP,N,JG,NG,NV,M1,AA,BB,
     1                XMEAN,TMEAN,JV,MDT,W,W,T,ALL,
     2                XQ,XMIN,IHISTI,EE,IQ,XT,VV,MD,LEST)
 3075 CLOSE(2)
      IF (LEST .LT. 100)  CLOSE(UNIT=10,STATUS='DELETE')
      CLOSE(51)
      CLOSE(LUN50)
      CLOSE(LUN51)

      RETURN
      END
