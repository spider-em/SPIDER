 
C ++********************************************************************
C                                                                     
C  MRSURF                                                            
C                  LONG FILENAMES                  JAN 89 al
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
C                                                                      *
C  MRSURF                                                             *
C                                                                      *
C  PURPOSE:  Probably limited to images of <= 512 !!                                                           *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE MRSURF

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      COMMON BUF1(512),BUF2(512),BUF3(512),BUF4(512),BUF(512)
     &      ,BUFZ(512),BUFZB(512),BC(512,4),BOX(512)

      CHARACTER(LEN=MAXNAM):: FLN1,FLN2,FLN3,FLN4,FLN3D,FILDUM,FILPAT

      CHARACTER *1  NULL

      NULL = CHAR(0)

C----- READ INPUT -----------------------------------------------

      CALL FILERD(FLN3D,NLET,NULL,'3-D',IRTFLG)

      CALL FILERD(FLN1,NLETP,NULL,'FIRST',IRTFLG)

      CALL FILCAD(FLN1,FILPAT,K1,IRTFLG)

      FLN2   = FLN1
      FLN3   = FLN1
      FLN4   = FLN1
      FILDUM = FLN1

      K2     = K1+1
      K3     = K2+1
      K4     = K3+1

      CALL FILGET(FILDUM,FLN2,NLET,K2,IRTFLG)
      CALL FILGET(FILDUM,FLN3,NLET,K3,IRTFLG)
      CALL FILGET(FILDUM,FLN4,NLET,K4,IRTFLG)

      MAXIM  = 0
      CALL OPFILEC(0,.FALSE.,FLN3D,10,'O',IFORM,NSAM,NROW,NSLICE,
     &                   MAXIM,' ',.FALSE.,IRTFLG)

      IFORM  = 1

      MAXIM  = 0
      CALL OPFILEC(0,.FALSE.,FLN1,11,'U',IFORM,NSAM,NROW,1,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
      MAXIM  = 0
      CALL OPFILEC(0,.FALSE.,FLN2,12,'U',IFORM,NSAM,NROW,1,
     &                   MAXIM,' ',.FALSE.,IRTFLG)

      MAXIM  = 0
      CALL OPFILEC(0,.FALSE.,FLN3,13,'U',IFORM,NSAM,NROW,1,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
      MAXIM  = 0
      CALL OPFILEC(0,.FALSE.,FLN4,14,'U',IFORM,NSAM,NROW,1,
     &                   MAXIM,' ',.FALSE.,IRTFLG)

      CALL RDPRMI(NS1,NS2,NOT_USED,'FIRST,LAST SLICE')
      CALL RDPRMI(NL1,NL2,NOT_USED,'FIRST,LAST X-LINE')
      CALL RDPRM(SW,NOT_USED,'THRESHOLD')
      CALL RDPRMI(IBACK,IDUM,NOT_USED,
     &    'BACKGROUND (0)ZERO,(1)MIN,(2)MAX,(3)LOCAL MIN,(4)LOCAL MAX')
      IF (IBACK.NE.3.AND.IBACK.NE.4) GOTO 70 
      CALL RDPRMI(LBOX,IDUM,NOT_USED,
     &    'BOX LENGTH FOR BACKGROUND FILTER')
      CALL RDPRM(BOFF,NOT_USED,'BACKGROUND OFFSET')
70    CONTINUE
      IF (NS1.EQ.1) NS1=2
      IF (NS2.EQ.NSLICE) NS2=NSLICE-1
      IF (NL1.EQ.1) NL1=2
      IF (NL2.EQ.NSAM) NL2=NSAM-1

C----- DETERMINE START VALUES ------------------------------------------

      CS=NS2-NS1+1
      CQ=NL2-NL1+1
      NCLEAR=MAX(NSAM,NSLICE)
      DO  K=1,NCLEAR
      BUF1(K)=0.
      BUF2(K)=0.
      BUF3(K)=0.
      BUF4(K)=0.
      ENDDO
      B1MIN=1000.
      B2MIN=1000.
      B3MIN=1000.
      B4MIN=1000.
      B1MAX=0.
      B2MAX=0.
      B3MAX=0.
      B4MAX=0.
      DO  K=1,4
      DO  L=1,NROW
         BC(L,K)=0.
      ENDDO
      ENDDO
      NOFF=NROW

C----- START OF THE CALCULATIONS ---------------------------------------

      DO  J=1,NROW
      N0=(NS1-2)*NOFF+J
      CALL REDLIN(10,BUFZ,NSAM,N0)
      N0=(NS1-1)*NOFF+J
      CALL REDLIN(10,BUF,NSAM,N0)
      KK=0
      DO  K=NS1,NS2
      K4=NS2+1-K
      KK=KK+1
      IND=K*NOFF
      N0=IND+J
      CALL REDLIN(10,BUFZB,NSAM,N0)
      LL=0
      DO 5 L=NL1,NL2
      L2=NL2+1-L
      LL=LL+1
      IF(BUF(L).LT.SW.AND.SW.GE.0.) GOTO 5
      IF(BUF(L).GT.SW.AND.SW.LT.0.) GOTO 5
      IF(BUF1(L).GT.0.0)  GOTO 51
      BUF1(L)=FLOAT(KK)

C----- INTERPOLATION ---------------------------------------------------

      IF (BUFZ(L).GE.BUF(L)) GOTO 51
      DIF=(SW-BUFZ(L))/(BUF(L)-BUFZ(L))
       IF(DIF.LT.0.) DIF=ABS((SW-BUFZ(L))/(BUF(L)+BUFZ(L)))
       BUF1(L)=KK-1.+DIF
C      CALL KONTR(L,1,L,DIF,SW,NCLEAR)
C-----------------------------------------------------------------------
51    IF(BUF3(K).GT.0.0) GOTO 52
       BUF3(K)=FLOAT(LL)

C----- INTERPOLATION ---------------------------------------------------

       IF(BUF(L-1).GE.BUF(L)) GOTO 52
       DIF=(SW-BUF(L-1))/(BUF(L)-BUF(L-1))
       IF(DIF.LT.0.) DIF=ABS((SW-BUF(L-1))/(BUF(L)+BUF(L+1)))
       BUF3(K)=LL-1.+DIF
C      CALL KONTR(K,3,L,DIF,SW,NCLEAR)
52    CONTINUE
C-----------------------------------------------------------------------
      BUF2(L2)=CS-FLOAT(KK)
C
C---- INTERPOLATION ----------------------------------------------------
C
      IF(BUFZB(L).GE.BUF(L)) GOTO 5
      DIF=(SW-BUFZB(L))/(BUF(L)-BUFZB(L))
      IF(DIF.LT.0.) DIF=ABS((SW-BUFZB(L))/(BUF(L)+BUFZB(L)))
      BUF2(L2)=CS-KK-1.+DIF
C     CALL KONTR(L2,2,L,DIF,SW,NCLEAR)
C-----------------------------------------------------------------------
5     CONTINUE
      LL=0
      DO 18 JJ=NL2,NL1,-1
      LL=LL+1
      IF(BUF(JJ).LT.SW.AND.SW.GE.0.) GOTO 18
      IF(BUF(JJ).GT.SW.AND.SW.LT.0.) GOTO 18
      IF(BUF4(K4).GT.0.0) GOTO 18
      BUF4(K4)=LL
C
C---- INTERPOLATION ----------------------------------------------------
      IF(BUF(JJ+1).GE.BUF(JJ)) GOTO 18
      DIF=(SW-BUF(JJ+1))/(BUF(JJ)-BUF(JJ+1))
      IF(DIF.LT.0.) DIF=ABS((SW-BUF(JJ+1))/(BUF(JJ)+BUF(JJ+1)))
      BUF4(K4)=LL-1+DIF
C     CALL KONTR(K4,4,JJ,DIF,SW,NCLEAR)
C-----------------------------------------------------------------------
18    CONTINUE
C
C---- SAVE THE LAST TWO LINES FOR THE NEXT SLICE -----------------------
C
      DO  JJ=1,NSAM
      BUFZ(JJ)=BUF(JJ)
      BUF(JJ)=BUFZB(JJ)
      ENDDO
C-----------------------------------------------------------------------
      ENDDO

C
C----- CALCULATE BACKGROUND --------------------------------------------
C
      IF(IBACK.NE.1.AND.IBACK.NE.3) GOTO 16
C
C----- MINIMUM BACKGROUND LOCAL AND TOTAL ------------------------------
C
      DO  JK=1,4
        BC(J,JK)=0.
      ENDDO
      BC1MIN=1000.
      BC2MIN=1000.
      BC3MIN=1000.
      BC4MIN=1000.
      DO  JJ=1,NSAM
      B2=BUF2(JJ)
      B1=BUF1(JJ)
      IF(BC2MIN.GT.B2.AND.B2.GT.0.0) BC2MIN=B2
      IF(BC1MIN.GT.B1.AND.B1.GT.0.0) BC1MIN=B1
      ENDDO 
      IF(BC1MIN.LT.1000) BC(J,1)=BC1MIN-BOFF
      IF(BC2MIN.LT.1000) BC(J,2)=BC2MIN-BOFF
      DO  JJ=1,NSLICE
      B4=BUF4(JJ)
      B3=BUF3(JJ)
      IF(BC4MIN.GT.B4.AND.B4.GT.0.0) BC4MIN=B4
      IF(BC3MIN.GT.B3.AND.B3.GT.0.0) BC3MIN=B3
      ENDDO
      IF(BC3MIN.LT.1000) BC(J,3)=BC3MIN-BOFF
      IF(BC4MIN.LT.1000) BC(J,4)=BC4MIN-BOFF
16    CONTINUE
      IF(B1MIN.GT.BC1MIN) B1MIN=BC1MIN
      IF(B2MIN.GT.BC2MIN) B2MIN=BC2MIN
      IF(B3MIN.GT.BC3MIN) B3MIN=BC3MIN
      IF(B4MIN.GT.BC4MIN) B4MIN=BC4MIN
C
      IF(IBACK.NE.2.AND.IBACK.NE.4) GOTO 23
C
C----- MAXIMUM BACKGROUND LOCAL AND TOTAL ------------------------------
C
      BC1MAX=0.
      BC2MAX=0.
      BC3MAX=0.
      BC4MAX=0.
      DO  JJ=1,NSAM
      B2=BUF2(JJ)
      B1=BUF1(JJ)
      IF(BC2MAX.LT.B2.AND.B2.GT.0.0) BC2MAX=B2
      IF(BC1MAX.LT.B1.AND.B1.GT.0.0) BC1MAX=B1
      ENDDO
      DO  JJ=1,NSLICE
      B4=BUF4(JJ)
      B3=BUF3(JJ)
      IF(BC4MAX.LT.B4.AND.B4.GT.0.0) BC4MAX=B4
      IF(BC3MAX.LT.B3.AND.B3.GT.0.0) BC3MAX=B3
      ENDDO
      IF(BC1MAX.GT.0.) BC(J,1)=BC1MAX+BOFF
      IF(BC2MAX.GT.0.) BC(J,2)=BC2MAX+BOFF
      IF(BC3MAX.GT.0.) BC(J,3)=BC3MAX+BOFF
      IF(BC4MAX.GT.0.) BC(J,4)=BC4MAX+BOFF
      IF(B1MAX.LT.BC1MAX) B1MAX=BC1MAX
      IF(B2MAX.LT.BC2MAX) B2MAX=BC2MAX
      IF(B3MAX.LT.BC3MAX) B3MAX=BC3MAX
      IF(B4MAX.LT.BC4MAX) B4MAX=BC4MAX
C
C-----------------------------------------------------------------------
C
23    CONTINUE
C
      CALL WRTLIN(11,BUF1,NSAM,J)
      CALL WRTLIN(12,BUF2,NSAM,J)
      CALL WRTLIN(13,BUF3,NSLICE,J)
      CALL WRTLIN(14,BUF4,NSLICE,J)
      DO  K=1,NCLEAR
      BUF1(K)=0.
      BUF2(K)=0.
      BUF3(K)=0.
      BUF4(K)=0.
      ENDDO
      ENDDO
      IF(IBACK.EQ.0) GOTO 17


C
C----- CALCULATE BACKGROUND --------------------------------------------
C
      IF(IBACK.NE.1) GOTO 24
C
C----- MINIMUM TOTAL ---------------------------------------------------
C
      B1MIN=B1MIN-1
      B2MIN=B2MIN-1
      B3MIN=B3MIN-1
      B4MIN=B4MIN-1
24    CONTINUE
      IF(IBACK.NE.2) GOTO 25
C
C----- MAXIMUM TOTAL ---------------------------------------------------
C
      B1MAX=B1MAX+1.
      B2MAX=B2MAX+1.
      B3MAX=B3MAX+1.
      B4MAX=B4MAX+1.
25    CONTINUE
      IF(IBACK.NE.3.AND.IBACK.NE.4) GOTO 26
C
C----- LOCAL MINIMUM AND MAXIMUM ---------------------------------------
C
      DO  J=1,4
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
      ENDDO
C
C----- BOX CONVOLUTION -------------------------------------------------
C
      LEN=(LBOX-0.5)/2
      DO  J=1,4
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
      ENDDO
C
C-----------------------------------------------------------------------
C
26    CONTINUE
      DO  K=1,NROW
      CALL REDLIN(11,BUF1,NSAM,K)
      CALL REDLIN(12,BUF2,NSAM,K)
      CALL REDLIN(13,BUF3,NSLICE,K)
      CALL REDLIN(14,BUF4,NSLICE,K)
      IF(IBACK.NE.1) GOTO 40
C
C----- MINIMUM TOTAL ---------------------------------------------------
C
      DO  L=1,NSAM
      IF(BUF1(L).LE.0.0) BUF1(L)=B1MIN
      IF(BUF2(L).LE.0.0) BUF2(L)=B2MIN
      ENDDO
      DO  J=1,NSLICE
      IF(BUF3(J).LE.0.0) BUF3(J)=B3MIN
      IF(BUF4(J).LE.0.0) BUF4(J)=B4MIN
      ENDDO
40    IF(IBACK.NE.2) GOTO 41
C
C----- MAXIMUM TOTAL ---------------------------------------------------
C
      DO  L=1,NSAM
      IF(BUF1(L).LE.0.0) BUF1(L)=B1MAX
      IF(BUF2(L).LE.0.0) BUF2(L)=B2MAX
      ENDDO
      DO  L=1,NSLICE
      IF(BUF3(L).LE.0.0) BUF3(L)=B3MAX
      IF(BUF4(L).LE.0.0) BUF4(L)=B4MAX
      ENDDO
41    CONTINUE
      IF(IBACK.NE.3.AND.IBACK.NE.4) GOTO 42
C
C----- LOCAL MINIMUM AND MAXIMUM ---------------------------------------
C
      DO  L=1,NSAM
      IF(BUF1(L).GT.0.) BUF1(L)=BUF1(L)-BC(K,1)
      IF(BUF2(L).GT.0.) BUF2(L)=BUF2(L)-BC(K,2)
      ENDDO
      DO  L=1,NSLICE
      IF(BUF3(L).GT.0.) BUF3(L)=BUF3(L)-BC(K,3)
      IF(BUF4(L).GT.0.) BUF4(L)=BUF4(L)-BC(K,4)
      ENDDO
42    CONTINUE
      CALL WRTLIN(11,BUF1,NSAM,K)
      CALL WRTLIN(12,BUF2,NSAM,K)
      CALL WRTLIN(13,BUF3,NSLICE,K)
      CALL WRTLIN(14,BUF4,NSLICE,K)
      ENDDO
17    CONTINUE

C-----------------------------------------------------------------------

C     WRITE(6,200) NL1,NL2,NS1,NS2,SW,NSAM,NROW,NSLICE,CS,CQ
200   FORMAT(' NL1,NL2,NS1,NS2:',4I5,' SW:',F10.2,
     &       'NSAM,NROW,NSLICE:',3I5,'CS,CQ:',F10.2)
      CLOSE (10)
      CLOSE (11)
      CLOSE (12)
      CLOSE (13)
      CLOSE (14)
      RETURN
      END
 
