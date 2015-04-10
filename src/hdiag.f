
C ++********************************************************************
C                                                                      *
C  HDIAG                                                               *
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
C  HDIAG                                                               *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

       SUBROUTINE HDIAG(H,NVAR,NDIM,IEGEN,U,NR,X,IQ)

       IMPLICIT REAL*8 (A-H,O-Z)
       IMPLICIT INTEGER*2 (I-N)
       DIMENSION H(NDIM,NDIM),U(NDIM,NDIM),X(NDIM),IQ(NDIM)

       N=NVAR
       IF(IEGEN) 15,10,15
  10   DO 14 I=1,N
       DO 14 J=1,N
       IF(I-J) 12,11,12
  11   U(I,J)=1.0
       GOTO 14
  12   U(I,J)=0.0
  14   CONTINUE
  15   NR=0
       IF(N-1) 1000,1000,17
  17   NMI1=N-1
       DO 30 I=1,NMI1
       X(I)=0.0
       IPL1=I+1
       DO 30 J=IPL1,N
       IF(X(I)-ABS(H(I,J))) 20,20,30
  20   X(I)=ABS(H(I,J))
       IQ(I)=J
  30   CONTINUE
       RAP=7.450580596D-9
       HDTEST=1.0D37
  40   DO 70 I=1,NMI1
       IF(I-1) 60,60,45
  45   IF(XMAX-X(I)) 60,70,70
  60   XMAX=X(I)
       IPIV=I
       JPIV=IQ(I)
  70   CONTINUE
       IF(XMAX) 1000,1000,80
  80   IF(HDTEST) 90,90,85
  85   IF(XMAX-HDTEST) 90,90,148
  90   HDIMIN=ABS(H(1,1))
       DO 110 I=2,N
       IF(HDIMIN-ABS(H(I,I))) 110,110,100
  100  HDIMIN=ABS(H(I,I))
  110  CONTINUE
       HDTEST=HDIMIN*RAP
       IF(HDTEST-XMAX) 148,1000,1000
  148  NR=NR+1
  150   XLAM=H(IPIV,JPIV)
       XMI=0.5*(H(IPIV,IPIV)-H(JPIV,JPIV))
       XNI=SQRT(XLAM**2+XMI**2)
       COSINE=SQRT((ABS(XMI)+XNI)/XNI*0.5)
       IF(ABS(COSINE) .LT. 1.0E-15)  COSINE=SIGN(1.0D-15 , COSINE)
       ZNAK=1.
       IF(XMI .LT.0.0) ZNAK=-1.0
       SINE=ZNAK*XLAM/(2.0*XNI*COSINE)
       TANG=SINE / COSINE
       HII=H(IPIV,IPIV)
       H(IPIV,IPIV)=COSINE**2*(HII+TANG*(2.*H(IPIV,JPIV)+TANG*H(JPIV,
     &  JPIV)))
       H(JPIV,JPIV)=COSINE**2*(H(JPIV,JPIV)-TANG*(2.*H(  IPIV,JPIV)-
     &TANG*HII))
       H(IPIV,JPIV)=0.0
       H(JPIV,IPIV)=0.0
       IF(H(IPIV,IPIV)-H(JPIV,JPIV)) 152,153,153
  152  HTEMP=H(IPIV,IPIV)
       H(IPIV,IPIV)=H(JPIV,JPIV)
       H(JPIV,JPIV)=HTEMP
       IF(SINE) 94,95,95
   94   HTEMP=COSINE
       GOTO 96
  95   HTEMP=-COSINE
  96   COSINE=ABS(SINE)
       SINE=HTEMP
  153  CONTINUE
       DO 350 I=1,NMI1
       IF(I-IPIV) 210,350,200
  200  IF(I-JPIV) 210,350,210
  210  IF(IQ(I)-IPIV) 230,240,230
  230  IF(IQ(I)-JPIV) 350,240,350
  240  K=IQ(I)
  250  HTEMP=H(I,K)
       H(I,K)=0.0
       IPL1=I+1
       X(I)=0.0
       DO 320 J=IPL1,N
       IF(X(I)-ABS(H(I,J))) 300,300,320
  300  X(I)=ABS(H(I,J))
       IQ(I)=J
  320  CONTINUE
       H(I,K)=HTEMP
  350  CONTINUE
       X(IPIV)=0.0
       X(JPIV)=0.0
       DO 530 I=1,N
       IF(I-IPIV) 370,530,420
  370  HTEMP=H(I,IPIV)
       H(I,IPIV)=COSINE*HTEMP+SINE*H(I,JPIV)
       IF(X(I)-ABS(H(I,IPIV))) 380,390,390
  380  X(I)=ABS(H(I,IPIV))
       IQ(I)=IPIV
  390  H(I,JPIV)=-SINE*HTEMP+COSINE*H(I,JPIV)
       IF(X(I)-ABS(H(I,JPIV))) 400,530,530
  400  X(I)=ABS(H(I,JPIV))
       IQ(I)=JPIV
       GOTO 530
  420  IF(I-JPIV) 430,530,480
  430  HTEMP = H(IPIV,I)
       H(IPIV,I)=COSINE*HTEMP+SINE*H(I,JPIV)
       IF(X(IPIV)-ABS(H(IPIV,I))) 440,450,450
  440  X(IPIV)=ABS(H(IPIV,I))
       IQ(IPIV)=I
  450  H(I,JPIV)=-SINE*HTEMP+COSINE*H(I,JPIV)
       IF(X(I)-ABS(H(I,JPIV))) 400,530,530
  480  HTEMP=H(IPIV,I)
       H(IPIV,I)=COSINE*HTEMP+SINE*H(JPIV,I)
       IF(X(IPIV)-ABS(H(IPIV,I))) 490,500,500
  490  X(IPIV)=ABS(H(IPIV,I))
       IQ(IPIV)=I
  500  H(JPIV,I)=-SINE*HTEMP+COSINE*H(JPIV,I)
       IF(X(JPIV)-ABS(H(JPIV,I))) 510,530,530
  510  X(JPIV)=ABS(H(JPIV,I))
       IQ(JPIV)=I
  530  CONTINUE
       IF(IEGEN) 40,540,40
  540  DO  I=1,N
       HTEMP=U(I,IPIV)
       U(I,IPIV)=COSINE*HTEMP+SINE*U(I,JPIV)
       U(I,JPIV)=-SINE*HTEMP+COSINE*U(I,JPIV)
       ENDDO
       GOTO 40
 1000  RETURN
       END
