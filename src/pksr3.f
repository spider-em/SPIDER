C++*********************************************************************
C
C PKSR3.F
C
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
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE  PKSR3(LUN,Q,NSAM,NROW,NSLICE,
     &      NSA1,NSA2,NRO1,NRO2,NSL1,NSL2,SGN,PEAK,NPC,RPC,ML,NMAX)

         DIMENSION  Q(NSAM,NROW,-1:1),NPC(3,ML),RPC(3,ML),PEAK(ML)
         DIMENSION  LX(-1:1),LY(-1:1),LZ(-1:1),LS(-1:1)
         LOGICAL    T,NEGV

         NMAX=0
         DO  1  L=NSL1,NSL2
            IF(L.EQ.NSL1)  THEN
               DO   J=-1,1
                  LZ(J)=J
	       ENDDO

               DO    J=-1,1
                  LSJ=MOD(L+J-1+NSLICE,NSLICE)+1
                  CALL RDSL_P(LUN,Q(1,1,J),NSAM,NROW,LSJ)
	       ENDDO
               DO  MZ=-1,1
                  DO  MY=1,NROW
                     DO  MX=1,NSAM
                        Q(MX,MY,MZ)=Q(MX,MY,MZ)*SGN
                     ENDDO
                  ENDDO
               ENDDO
           ELSE
C  l>nsl1
               J=MOD(L-1-NSL1,3)-1
               LSJ=MOD(L+NSLICE,NSLICE)+1
               CALL RDSL_P(LUN,Q(1,1,J),NSAM,NROW,LSJ)
               DO  MY=1,NROW
                  DO  MX=1,NSAM
                     Q(MX,MY,J)=Q(MX,MY,J)*SGN
                  ENDDO
               ENDDO
C             rotate numbering of slices
               DO  J=-1,1
                  LZ(J)=LZ(J)+1
                 IF(LZ(J).EQ.+2)  LZ(J)=-1
               ENDDO
            ENDIF
C
            DO  1  J=NRO1,NRO2
C               IYM=MOD(J-2+NROW,NROW)+1
C               IYP=MOD(J+NROW,NROW)+1
               LY(-1)=MOD(J-2+NROW,NROW)+1
               LY(1)=MOD(J+NROW,NROW)+1
               LY(0)=J
               DO  1  I=NSA1,NSA2
C                  IXM=MOD(I-2+NSAM,NSAM)+1
C                 IXP=MOD(I+NSAM,NSAM)+1
                  LX(-1)=MOD(I-2+NSAM,NSAM)+1
                  LX(1)=MOD(I+NSAM,NSAM)+1
                  LX(0)=I
                  T=.TRUE.
                  QC=Q(I,J,LZ(0))
                  DO    JZ=-1,1
                     DO    JY=-1,1
                        DO    JX=-1,1
                           IF(QC.LT.Q(LX(JX),LY(JY),LZ(JZ)))  T=.FALSE.
	                ENDDO	 
	             ENDDO	 
	          ENDDO
                  IF (T)  THEN
                     IF(NMAX.EQ.0)  THEN
                        NMAX=1
                        PEAK(1)=QC
                        NPC(1,1)=I
                        NPC(2,1)=J
                        NPC(3,1)=L
C                       get floating point shifts
                        NN = 1
                        ASSIGN  501 TO LABA
                        GOTO  5505

501                     CONTINUE
                     ELSE
                     DO    N=NMAX,1,-1
                        IF(QC.LT.PEAK(N)) THEN
                           NN=N+1
                           GOTO  22
                        ENDIF
	             ENDDO
                     NN=1
22                   CONTINUE
C                    no place for more peaks
                     IF(NN.GT.ML)  GOTO 1
                     DO    M=MIN0(NMAX+1,ML),NN+1,-1
                        PEAK(M)=PEAK(M-1)
                        DO  MM=1,3
                           NPC(MM,M)=NPC(MM,M-1)
                           RPC(MM,M)=RPC(MM,M-1)
                        ENDDO
                     ENDDO
                     PEAK(NN)=QC
                     NPC(1,NN)=I
                     NPC(2,NN)=J
                     NPC(3,NN)=L
C                    get floating point shifts
                     ASSIGN  502 TO LABA

                     GOTO  5505
502                  CONTINUE
                     NMAX=MIN0(NMAX+1,ML)
                  ENDIF
               ENDIF
1        CONTINUE
         RETURN

C        Procedure to get the floating point X, Y, Z shifts.
C        Computes 1D interpolation along 3 axes, using the 2 elements
C        on either side of the central voxel of the 3x3x3 cube

5505     NEGV=.FALSE.

         RX=0.0
         RY=0.0
         RZ=0.0

C        the central element
         QMAX = SGN * Q(LX(0),LY(0),LZ(0))

C    X shift
         QA = SGN * Q(LX(-1),LY(0),LZ(0))
         QC = SGN * Q(LX(1),LY(0),LZ(0))
         CALL  PKSR3_SUB(QMAX,QA,QC,SGN,RX)

C    Y shift
         QA = SGN * Q(LX(0),LY(-1),LZ(0))
         QC = SGN * Q(LX(0),LY(1),LZ(0))
         CALL  PKSR3_SUB(QMAX,QA,QC,SGN,RY)

C    Z shift
         QA = SGN * Q(LX(0),LY(0),LZ(-1))
         QC = SGN * Q(LX(0),LY(0),LZ(1))
         CALL  PKSR3_SUB(QMAX,QA,QC,SGN,RZ)

         RPC(1,NN)= RX + I
         RPC(2,NN)= RY + J
         RPC(3,NN)= RZ + L

         GOTO  LABA

         END

C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************
C Returns R, the subpixel amount of shift for the 3 values: QA,QMAX,QC

        SUBROUTINE  PKSR3_SUB(QMAX,QA,QC,SGN,R)

         IF (SGN.GT.0) THEN
            IF(QA.LT.QC) THEN
               QMIN = QA
               QMID = QC
               QF = 0.5
            ELSE
               QMIN = QC
               QMID = QA
               QF = -0.5
            ENDIF
         ENDIF
         IF (SGN.LT.0) THEN 
            IF(QA.GT.QC) THEN
               QMIN = QA
               QMID = QC
               QF = 0.5
            ELSE
               QMIN = QC
               QMID = QA
               QF = -0.5
            ENDIF
         ENDIF

         IF ((QMAX - QMIN).NE.0) THEN
            QM = QF / (QMAX - QMIN)
            R = QM * (QMID - QMIN)
         ENDIF

         RETURN

         END
