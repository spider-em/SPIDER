
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


       SUBROUTINE MYREAD(M,MGR,KG,ALL,T,XMEAN,TMEAN,XMIN,XMAX,N,
     &                    JG,LEST)

       IMPLICIT REAL*8 (A-H,O-Z)
       IMPLICIT INTEGER*2 (I-N)
       DIMENSION ALL(M,M,1),T(M,M)
       DIMENSION XMEAN(M,1),XMIN(M,1),XMAX(M,1)
       DIMENSION N(1),JG(1),TMEAN(1)
       DIMENSION V(50)
       INTEGER * 4 LERC

       CHARACTER*10 IBI
       LOGICAL*1 IFR

      IF (LEST.LT.100)  THEN
         QF=REAL(LEST)/100.0
         LERC=0
         OPEN(UNIT=10,FILE='RANF.DAT',STATUS='UNKNOWN',ACCESS='DIRECT',
     &       RECL=1)
         ISEED=21717
      ENDIF

      M1=M-1
      NS=0
 1    READ(4,5676,END=10) (V(K),K=1,M)
5676  FORMAT(2X,F3.1,4(1X,1PE14.7),/,4(1X,1PE14.7),/,1(1X,1PE14.7))
      IF (LEST.LT.100)  THEN
         IF(QF.GE.RAND_P(ISEED)) THEN
            IFR=.TRUE.
         ELSE
            IFR=.FALSE.
      ENDIF
      LERC=LERC+1
      WRITE(10,REC=LERC)  IFR
         IF(.NOT.IFR) GOTO  1
      ENDIF
       I=V(KG)
       DO  K=1,MGR
          IF(I.EQ.JG(K)) GOTO 4
       ENDDO
       GOTO 1

 4     V(KG)=V(M)
       WRITE(IBI,7023) NS+1
 7023  FORMAT(I8)
C                       CALL  WRTXT(
C     & IBI,8,51,3,3)
       DO  I=1,M1
          XMEAN(I,K)=XMEAN(I,K)+V(I)
          TMEAN(I)=TMEAN(I)+V(I)
          IF(V(I).GT.XMAX(I,K)) XMAX(I,K)=V(I)
          IF(V(I).LT.XMIN(I,K)) XMIN(I,K)=V(I)
          DO  J=I,M1
             X=V(I)*V(J)
             ALL(I,J,K)=ALL(I,J,K)+X
             T(I,J)=T(I,J)+X
          ENDDO
       ENDDO
       N(K)=N(K)+1
       NS=NS+1
       GOTO 1

 10    CONTINUE
       DO 6 I=1,MGR
          EK=N(I)
          IF(N(I).EQ.0) GOTO 6
          DO  J=1,M1
             DO  L=J,M1
                 ALL(J,L,I)=ALL(J,L,I)-XMEAN(J,I)*XMEAN(L,I)/EK
             ENDDO
             XMEAN(J,I)=XMEAN(J,I)/EK
          ENDDO
 6     CONTINUE
       EN=NS
       DO  J=1,M1
          DO  L=J,M1
             T(J,L)=T(J,L)-TMEAN(J)*TMEAN(L)/EN
             T(L,J)=T(J,L)
          ENDDO
          TMEAN(J)=TMEAN(J)/EN
       ENDDO
       END
