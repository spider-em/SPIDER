
C ++********************************************************************
C                                                                      *
C                                                                      *
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
C                                                                      *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE  ROT3_P(POLTMP,NDIM,ROT)

	PARAMETER  (MAXROT=29)
	IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
	REAL POLTMP(65)
	DIMENSION  ROT(MAXROT,MAXROT,MAXROT)
        COMMON  DUMMY(26000),POL(MAXROT),AMAT(MAXROT,MAXROT,MAXROT)
     &  ,FAT(MAXROT)
     & ,AM(MAXROT),BM(MAXROT),CM(MAXROT),DM(MAXROT),EM(MAXROT)
     & ,FM(MAXROT),GM(MAXROT),HM(MAXROT)
     & ,BB(MAXROT),W(132)
     & ,RPOT(MAXROT),RK(MAXROT),RP(MAXROT),RZ(MAXROT)
	INTEGER  I,K,L,M,P,Q,R,S
C      EQUIVALENCE  (AMAT,POL),(RPOT,FAT),(RK,BB),(RP,DM)
C      EQUIVALENCE  (QQ,AMAT),(POLTMP,POL)
       REAL  IP,IQ,IB
       DATA  A,B,C,D,E,F,G,H/.25,.25,.25,.25,.25,.25,.25,-.75/

	ONE=1.0
	ZERO=0.0
      AM(1)=ONE
      BM(1)=ONE
      CM(1)=ONE
      DM(1)=ONE
      EM(1)=ONE
      FM(1)=ONE
      GM(1)=ONE
      HM(1)=ONE
      FAT(1)=ONE
C	OPEN(20,FILE='FILTER.OUT',FORM='UNFORMATTED',STATUS='UNKNOWN')
C       READ(20)  NDIM
C       READ(20)   (POLTMP(J),J=NDIM,1,-1)
C	CLOSE(20)
	DO    J=NDIM,1,-1
 	POL(NDIM-J+1)=POLTMP(J)
	ENDDO
       NDM1=NDIM-1
       DO  2  J=1,NDM1
       AM(J+1)=AM(J)*A
       BM(J+1)=BM(J)*B
       CM(J+1)=CM(J)*C
       DM(J+1)=DM(J)*D
       EM(J+1)=EM(J)*E
       FM(J+1)=FM(J)*F
       GM(J+1)=GM(J)*G
       HM(J+1)=HM(J)*H
       FATT=ONE
       IF(J.LE.1)  GOTO  2
       DO    K=1,J
         FATT=FATT*FLOAT(K)
       ENDDO
 2     FAT(J+1)=FATT
       DO    J=2,NDIM
         POL(J)=2.0*POL(J)
       ENDDO
C
       CALL  TCNP(ONE,ZERO,BB,NDIM,POL,W)

C
       DO    I=1,NDIM
       DO    J=1,NDIM
	DO    K=1,NDIM
       AMAT(I,J,K)=ZERO
       ROT(I,J,K)=ZERO
       ENDDO
	ENDDO
	ENDDO

       DO    I=1,NDIM
       DO    K=1,I
       COEF1=FAT(I)/FAT(K)
       K1=I-K+1
       DO    L=1,K1              
       COEF2=COEF1/FAT(L)
       K2=I-K-L+2
       DO    M=1,K2              
       COEF3=COEF2/FAT(M)
       K3=I-K-L-M+3
       DO    P=1,K3              
       COEF4=COEF3/FAT(P)
       K4=I-K-L-M-P+4
       DO    Q=1,K4              
       COEF5=COEF4/FAT(Q)
       K5=I-K-L-M-P-Q+5
       DO    R=1,K5              
       COEF6=COEF5/FAT(R)
       K6=I-K-L-M-P-Q-R+6
       DO    S=1,K6              
       N=I-K-L-M-P-Q-R-S+7
       COEF=COEF6/FAT(S)/FAT(N)
       COEFF=AM(K)*BM(L)*CM(M)*DM(P)*EM(Q)*FM(R)*GM(S)*HM(N)
       KM1=K+P+Q+S-3
       LM1=L+P+R+S-3
       MM1=M+Q+R+S-3
       AMAT(KM1,LM1,MM1)=AMAT(KM1,LM1,MM1)+COEF*COEFF*BB(I)
       ENDDO
       ENDDO
       ENDDO
       ENDDO
       ENDDO
       ENDDO
       ENDDO
       ENDDO
C
       DO    I=1,NDIM
       IM1=I-1
       RPOT(1)=ONE
       IF(IM1.EQ.0)  GOTO  8
       RPOT(2)=ONE
       IF(IM1.EQ.1)  GOTO  8
       IP=1
       DO    K=2,IM1
       DO    L=2,K
       IQ=RPOT(L)
       IB=IP+IQ
       IP=IQ
       RPOT(L)=IB
       ENDDO
       RPOT(K+1)=ONE
       ENDDO
 8     JK=I/2
       NODD=I-JK*2
       IF(IM1)   9,9,10
 9     RN=ONE
       GO TO 11
 10    RN=2.0**(IM1-1)
 11    IF(NODD)  12,12,14
 12    DO   J=1,JK
         RK(J)=RPOT(J)/RN
       ENDDO
       JJK=JK
       GO TO 18
 14    JKK=JK+1
       DO    J=1,JKK
         RK(J)=RPOT(J)/RN
       ENDDO
       IF(IM1)  17,17,16
 16    RK(JKK)=RK(JKK)/2.0
 17    JJK=JKK
 18    DO    L=1,NDIM
       LM1=L-1
       RPOT(1)=ONE
       IF(LM1.EQ.0)  GO TO 21
       RPOT(2)=ONE
       IF(LM1.EQ.1)  GO TO 21
       IP=1
       DO    K=2,LM1
       DO    M=2,K
       IQ=RPOT(M)
       IB=IP+IQ
       IP=IQ
       RPOT(M)=IB
       ENDDO
       RPOT(K+1)=ONE
       ENDDO
 21    JP=L/2
       NODD=L-JP*2
       IF(LM1)  22,22,23
 22    RN=ONE
       GO TO 24
 23    RN=2.0**(LM1-1)
 24    IF(NODD)   25,25,27
 25    DO    J=1,JP
       RP(J)=RPOT(J)/RN
       ENDDO
       JJP=JP
       GO TO 31
 27    JPP=JP+1
       DO    J=1,JPP
         RP(J)=RPOT(J)/RN
       ENDDO
       IF(LM1)   30,30,29
 29    RP(JPP)=RP(JPP)/2.0
 30    JJP=JPP
C
 31    DO    KZ=1,NDIM
       KM1=KZ-1
       RPOT(1)=ONE
       IF(KM1.EQ.0)  GO TO 41
       RPOT(2)=ONE
       IF(KM1.EQ.1)  GO TO 41
       IP=1
       DO    K=2,KM1
       DO    M=2,K
       IQ=RPOT(M)
       IB=IP+IQ
       IP=IQ
       RPOT(M)=IB
       ENDDO
       RPOT(K+1)=ONE
       ENDDO
 41    JZ=KZ/2
       NODD=KZ-JZ*2
       IF(KM1)  42,42,43
 42    RN=ONE
       GO TO 44
 43    RN=2.0**(KM1-1)
 44    IF(NODD)   45,45,47
 45    DO    J=1,JZ
         RZ(J)=RPOT(J)/RN
       ENDDO
       JJZ=JZ
       GO TO 51
 47    JZZ=JZ+1
       DO    J=1,JZZ
         RZ(J)=RPOT(J)/RN
       ENDDO
       IF(KM1)   50,50,49
 49    RZ(JZZ)=RZ(JZZ)/2.0
 50    JJZ=JZZ
 51    DO    J=1,JJK
       M=IM1-2*J+3
       DO    K=1,JJP
       N=LM1-2*K+3
       DO    JZZZ=1,JJZ
       MZZZ=KM1-2*JZZZ+3
       ROT(M,N,MZZZ)=ROT(M,N,MZZZ)+RK(J)*RP(K)*RZ(JZZZ)*AMAT(I,L,KZ)
       ENDDO
	ENDDO
	ENDDO

       ENDDO
       ENDDO
       ENDDO
	DO    I=2,NDIM
	DO    J=2,NDIM
       ROT(I,J,1)=ROT(I,J,1)/4.0
       ROT(1,I,J)=ROT(1,I,J)/4.0
       ROT(I,1,J)=ROT(I,1,J)/4.0
	ENDDO
	ENDDO
       DO    I=2,NDIM
	ROT(I,1,1)=ROT(I,1,1)/2.0
	ROT(1,I,1)=ROT(1,I,1)/2.0
	ROT(1,1,I)=ROT(1,1,I)/2.0
       DO   J=2,NDIM
       DO  K=2,NDIM
       ROT(I,J,K)=ROT(I,J,K)/8.0
       ENDDO
	ENDDO
	ENDDO
       END
