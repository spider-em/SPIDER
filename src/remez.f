

C ++********************************************************************
C                                                                      *
C REMEZ                                                                *
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
C   REMEZ(EDGE,NBANDS,NOUT)
C                                                                   *
C***********************************************************************

          SUBROUTINE  REMEZ(EDGE,NBANDS,NOUT)

          COMMON TDM(2000),
     &	  PI2,AD,DEVI,X,Y,GRID,DES,WT,ALPHA,IEXT,NFCNS,NGRID 

          DIMENSION  IEXT(66),AD(66),ALPHA(66),X(66),Y(66)
          DIMENSION  EDGE(20) 
          DIMENSION  DES(1056),GRID(1056),WT(1056)
          DIMENSION  A(66),P(65),Q(65)
          DOUBLE PRECISION  PI2,DNUM,DDEN,DTEMP,A,P,Q 
          DOUBLE PRECISION  AD,DEVI,X,Y
          DOUBLE PRECISION  DEE,GEE
          EXTERNAL  DEE,GEE

          ARCOS(QT)=1.57079632-ATAN(QT/SQRT(1.0-QT*QT)) 

          ITRMAX=25 
          DEVL=-1.0 
          NZ=NFCNS+1
          NZZ=NFCNS+2 
          NITER=0 

 100      CONTINUE
          IEXT(NZZ)=NGRID+1 
          NITER=NITER+1 
          IF (NITER.GT.ITRMAX)  GO TO 400
          DO J=1,NZ 
             JJ=IEXT(J)
             DTEMP=GRID(JJ)
             DTEMP=DCOS(DTEMP*PI2) 
             X(J)=DTEMP
	  ENDDO
          JET=(NFCNS-1)/15+1
          DO J=1,NZ 
             AD(J)=DEE(J,NZ,JET) 
	  ENDDO
          DNUM=0.0
          DDEN=0.0
          K=1 
          DO    J=1,NZ 
            L=IEXT(J) 
            DTEMP=AD(J)*DES(L)
            DNUM=DNUM+DTEMP 
            DTEMP=FLOAT(K)*AD(J)/WT(L)
            DDEN=DDEN+DTEMP 
            K=-K
	  ENDDO
          DEVI=DNUM/DDEN 
          NU=1
          IF (DEVI.GT.0.0)  NU=-1 
          DEVI=-FLOAT(NU)*DEVI
          K=NU
          DO    J=1,NZ 
            L=IEXT(J) 
            DTEMP=FLOAT(K)*DEVI/WT(L)
            Y(J)=DES(L)+DTEMP 
            K=-K
	  ENDDO
          IF (DEVI .LT. DEVL)  THEN
             WRITE(NOUT,1111)
1111         FORMAT(' ********** FAILURE TO CONVERGE **********'/
     &              '0PROBABLE CAUSE IS MACHINE ROUNDING ERROR'/ 
     &              '0THE IMPULSE RESPONSE MAY BE CORRECT'/
     &              '0CHECK WITH A FREQUENCY RESPONSE')
              GO TO 400
          ENDIF 

          DEVL=DEVI
          JCHNGE=0
          K1=IEXT(1)
          KNZ=IEXT(NZ)
          KLOW=0
          NUT=-NU 
          J=1 


 200          IF (J.EQ.NZZ)   YNZ=COMP 
          IF (J.GE.NZZ)  GO TO 300 
          KUP=IEXT(J+1) 
          L=IEXT(J)+1 
          NUT=-NUT
          IF (J.EQ.2)  Y1=COMP 
          COMP=DEVI
          IF (L.GE.KUP)  GO TO 220 
          ERR=GEE(L,NZ) 
          ERR=(ERR-DES(L))*WT(L)
          DTEMP=FLOAT(NUT)*ERR-COMP 
          IF (DTEMP.LE.0.0)  GO TO 220 
          COMP=FLOAT(NUT)*ERR 
 210      L=L+1 
          IF (L.GE.KUP)  GO TO 215 
          ERR=GEE(L,NZ) 
          ERR=(ERR-DES(L))*WT(L)
          DTEMP=FLOAT(NUT)*ERR-COMP
          IF (DTEMP.LE.0.0)   GO TO 215
          COMP=FLOAT(NUT)*ERR
          GO TO 210

 215      IEXT(J)=L-1
          J=J+1
          KLOW=L-1
          JCHNGE=JCHNGE+1
          GO TO 200

 220      L=L-1 
 225      L=L-1 
          IF (L.LE.KLOW) GO TO 250
          ERR=GEE(L,NZ)
          ERR=(ERR-DES(L))*WT(L)
          DTEMP=FLOAT(NUT)*ERR-COMP
          IF (DTEMP.GT.0.0) GO TO 230
          IF (JCHNGE.LE.0) GO TO 225
          GO TO 260
 
 230      COMP=FLOAT(NUT)*ERR 
 235      L=L-1 
          IF (L.LE.KLOW)  GO TO 240
          ERR=GEE(L,NZ) 
          ERR=(ERR-DES(L))*WT(L)
          DTEMP=FLOAT(NUT)*ERR-COMP 
          IF (DTEMP.LE.0.0) GO TO 240 
          COMP=FLOAT(NUT)*ERR 
          GO TO 235 

 240      KLOW=IEXT(J)
          IEXT(J)=L+1 
          J=J+1 
          JCHNGE=JCHNGE+1 
          GO TO 200 

 250      L=IEXT(J)+1 
          IF (JCHNGE.GT.0)  GO TO 215
 255      L=L+1 
          IF (L.GE.KUP)  GO TO 260 
          ERR=GEE(L,NZ) 
          ERR=(ERR-DES(L))*WT(L)
          DTEMP=FLOAT(NUT)*ERR-COMP 
          IF (DTEMP.LE.0.0)  GO TO 255 
          COMP=FLOAT(NUT)*ERR 
          GO TO 210
 
 260      KLOW=IEXT(J)
          J=J+1 
          GO TO 200 

 300      IF (J.GT.NZZ)  GO TO 320 
          IF (K1.GT.IEXT(1))  K1=IEXT(1) 
          IF (KNZ.LT.IEXT(NZ))   KNZ=IEXT(NZ)
          NUT1=NUT
          NUT=-NU 
          L=0 
          KUP=K1
          COMP=YNZ*1.00001
          LUCK=1
 310      L=L+1 
          IF (L.GE.KUP)  GO TO 315 
          ERR=GEE(L,NZ) 
          ERR=(ERR-DES(L))*WT(L)
          DTEMP=FLOAT(NUT)*ERR-COMP 
          IF (DTEMP.LE.0.0)  GO TO 310 
          COMP=FLOAT(NUT)*ERR 
          J=NZZ 
          GO TO 210
 
 315      LUCK=6
          GO TO 325 

 320      IF (LUCK.GT.9) GO TO 350
          IF (COMP.GT.Y1) Y1=COMP 
          K1=IEXT(NZZ)
 325      L=NGRID+1 
          KLOW=KNZ
          NUT=-NUT1
          COMP=Y1*1.00001 
 330      L=L-1 
          IF (L.LE.KLOW) GO TO 340
          ERR=GEE(L,NZ) 
          ERR=(ERR-DES(L))*WT(L)
          DTEMP=FLOAT(NUT)*ERR-COMP 
          IF (DTEMP.LE.0.0) GO TO 330 
          J=NZZ 
          COMP=FLOAT(NUT)*ERR 
          LUCK=LUCK+10
          GO TO 235 

 340      IF (LUCK.EQ.6)  GO TO 370
          DO J=1,NFCNS
            NZJ=NZ-J
            NZZJ=NZZ-J
            IEXT(NZZJ)=IEXT(NZJ)
	  ENDDO
          IEXT(1)=K1
          GO TO 100
 
 350          KN=IEXT(NZZ)
          DO    J=1,NFCNS
           IEXT(J)=IEXT(J+1) 
	  ENDDO
          IEXT(NZ)=KN 
          GO TO 100
 
 370          IF (JCHNGE.GT.0)  GO TO 100
 400          CONTINUE
          NM1=NFCNS-1 
          FSH=1.0E-06 
          GTEMP=GRID(1) 
          X(NZZ)=-2.0 
          CN=2*NFCNS-1
          DELF=1.0/CN 
          L=1 
          KKK=0 
          IF (EDGE(1).EQ.0.0 .AND. EDGE(2*NBANDS).EQ.0.5)  KKK=1 
          IF (NFCNS.LE.3)  KKK=1 
          IF (KKK.EQ.1)  GO TO 405 
          DTEMP=DCOS(PI2*GRID(1)) 
          DNUM=DCOS(PI2*GRID(NGRID))
          AA=2.0/(DTEMP-DNUM) 
          BB=-(DTEMP+DNUM)/(DTEMP-DNUM)
 405      CONTINUE
          DO    J=1,NFCNS
          FT=FLOAT(J-1)*DELF
          XT=DCOS(PI2*FT) 
          IF (KKK.EQ.1)  GO TO 410 
          XT=(XT-BB)/AA 
          FT=ARCOS(XT)/PI2
 410      XE=X(L) 
          IF (XT.GT.XE)  GO TO 420 
          IF ((XE-XT).LT.FSH)  GO TO 415 
          L=L+1 
          GO TO 410
 
 415      A(J)=Y(L) 
          GO TO 425 

 420      IF ((XT-XE).LT.FSH) GO TO 415 
          GRID(1)=FT
          A(J)=GEE(1,NZ)
 425          CONTINUE
          IF (L.GT.1) L=L-1 
	  ENDDO
          GRID(1)=GTEMP 
          DDEN=PI2/CN 
          DO J=1,NFCNS
            DTEMP=0.0 
            DNUM=FLOAT(J-1)*DDEN
            IF (NM1.LT.1) GO TO 505
              DO K=1,NM1
                 DTEMP=DTEMP+A(K+1)*DCOS(DNUM*FLOAT(K))
	      ENDDO
 505        DTEMP=2.0*DTEMP+A(1)
            ALPHA(J)=DTEMP
	  ENDDO
          DO J=2,NFCNS
             ALPHA(J)=2.0*ALPHA(J)/CN
	  ENDDO
          ALPHA(1)=ALPHA(1)/CN
          IF (KKK.EQ.1)  GO TO 545 
          P(1)=2.0*ALPHA(NFCNS)*BB+ALPHA(NM1) 
          P(2)=2.0*AA*ALPHA(NFCNS)
          Q(1)=ALPHA(NFCNS-2)-ALPHA(NFCNS)
          DO  540  J=2,NM1
          IF (J.LT.NM1)  GO TO 515 
          AA=0.5*AA 
          BB=0.5*BB
 
 515      CONTINUE
          P(J+1)=0.0
          DO K=1,J
             A(K)=P(K) 
             P(K)=2.0*BB*A(K)
	  ENDDO
          P(2)=P(2)+A(1)*2.0*AA 
          JM1=J-1 
          DO K=1,JM1
             P(K)=P(K)+Q(K)+AA*A(K+1)
	  ENDDO
          JP1=J+1 
          DO    K=3,JP1
             P(K)=P(K)+AA*A(K-1) 
	  ENDDO
          IF (J.EQ.NM1)  GO TO 540 
          DO K=1,J
             Q(K)=-A(K)
	  ENDDO
          NF1J=NFCNS-1-J
          Q(1)=Q(1)+ALPHA(NF1J) 
 540          CONTINUE
          DO J=1,NFCNS
             ALPHA(J)=P(J) 
	  ENDDO
 545      CONTINUE

          IF (NFCNS.GT.3)  RETURN

          ALPHA(NFCNS+1)=0.0
          ALPHA(NFCNS+2)=0.0

          END 
