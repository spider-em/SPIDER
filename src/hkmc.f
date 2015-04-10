
C ++********************************************************************
C                                                                      *
C  HKMC.F                                                              *
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
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE  HKMC
     &    (CIRSEED,CIRC,CIRNEW,NK,LCIRC,IP,IQ,ES,EC,E,DIST,ROT,NRING,
     &     TEMP,MAXRIN,JACUP,NUMR,MAXIT,NIMA,MODE,TOTAL,NDUM,MIRROR)

        INCLUDE 'CMBLOCK.INC'

        INTEGER           :: NUMR(3,NRING),MAXRIN
        INTEGER*2         :: IQ(NK),IP(NIMA)
        REAL              :: CIRC(LCIRC,NIMA),CIRSEED(LCIRC,NK)
        REAL              :: CIRNEW(LCIRC)
        REAL              :: DIST(NIMA),ROT(NIMA)
        DOUBLE PRECISION  :: TEMP(MAXRIN,2),
     &                       ENER,PEAK,PEAM,SNEW,EC(NIMA),ES(NK),E(NK)
        DOUBLE PRECISION  :: F,DS,WSS,WSSMAX,TOTAL
        LOGICAL           :: CH_ANG,ZEROGR
        CHARACTER*1       :: MODE
        LOGICAL           :: MIRROR

C       ES  - ENERGY OF A SEED
C       E   - SUM OF SQUARES IN A CLUSTER

        DO J=1,NK
           ES(J) = ENER(CIRSEED(1,J),LCIRC,NRING,NUMR,MODE)
        ENDDO

        WSSMAX = 1.0D30

        DO  ITER=1,MAXIT

C          FIND ASSIGNMENT
           IQ = 0

           DO I=1,NIMA
              F = 1.0D30
              DO  J=1,NK
                 CALL CROSRNG(CIRSEED(1,J),CIRC(1,I),
     &                        LCIRC,NRING,TEMP,TEMP(1,2),
     &                        MAXRIN,JACUP,NUMR,PEAK,TOT,MODE)

	         IF (MIRROR)  THEN
                    CALL CROSRMG(CIRSEED(1,J),CIRC(1,I),
     &                        LCIRC,NRING,TEMP,TEMP(1,2),
     &                        MAXRIN,JACUP,NUMR,PEAM,TMT,MODE)
                    IF (PEAM .GT. PEAK)  THEN
                       DS = ES(J)+EC(I)-2.0*PEAM
C
                       IF (DS .LT. F)  THEN
                          ROT(I) = -TMT
                          F      = DS
                          IR     = J
                       ENDIF
                       GOTO 151
                    ENDIF
                 ENDIF

                 DS=ES(J)+EC(I)-2.0*PEAK
C                WRITE(NOUT,*) 'Object #',I,'  Group #',J,'  Dist.=',DS

                 IF (DS .LT. F)  THEN
                    ROT(I) = TOT
                    F      = DS
                    IR     = J
                 ENDIF
151              CONTINUE
              ENDDO
              IP(I)  = IR
              IQ(IR) = IQ(IR)+1
C             WRITE(NOUT,*)  'Set to group #',IP(I)
           ENDDO

          
           WRITE(NOUT,2020) ITER
2020       FORMAT(/,'  Iteration:',I3,
     &              '  Number of objects in clusters:')
           LK  = NK / 10
           MK  = MOD(NK,10)
           IP1 = 1
           IP2 = 0
           DO J=1,LK+MIN0(MK,1)
              IP2 = MIN0(IP1+9,NK)
              WRITE(NOUT,2092) (I,I=IP1,IP2)
              WRITE(NOUT,2092) (IQ(I),I=IP1,IP2)
              !WRITE(NOUT,2092)
              IP1 = IP1+10
           ENDDO
2092       FORMAT(1X,10I5)

C          WRITE(NOUT,  *,' Assignment'
C          WRITE(NOUT,  2002,IP

C          BUILD NEW CENTERS WITH CORRECTIONS
           DO  201 J=1,NK
              IF (IQ(J) .EQ. 1)  THEN
                 DO I=1,NIMA
                    II = I
                    IF (IP(I) .EQ. J) GOTO 468
                 ENDDO
468              CONTINUE
                 CIRSEED(:,J) = CIRC(:,II)
                 ROT(II)      = 0.0
                 DIST(II)     = 0.0
                 ES(J)        = ENER(CIRSEED(1,J),LCIRC,NRING,NUMR,MODE)
                 E(J)         = 0.0
                 GOTO  201
              ENDIF
              CIRSEED(:,J) = 0.0
C
C             BUILD FIRST AVERAGE
C
              IQ(J) = 0
              DO  203  I=1,NIMA
                 IF (IP(I).NE.J)  GOTO  203
                 IQ(J)=IQ(J)+1
                 TOT  =ROT(I)
                 IMI=IQ(J)
                 IF (TOT.GT.0.0)  THEN
                    CALL UPDTC(CIRSEED(1,J),CIRC(1,I),
     &                   LCIRC,NRING,NUMR,TOT,MAXRIN,IMI)
                 ELSE
                    CALL UPDTM(CIRSEED(1,J),CIRC(1,I),
     &                   LCIRC,NRING,NUMR,-TOT,MAXRIN,IMI)
                 ENDIF
203           CONTINUE

              ES(J) = ENER(CIRSEED(1,J),LCIRC,NRING,NUMR,MODE)
              E(J)  = 0.0

C2001         FORMAT(10(1X,F6.1))

C             ITERATIONS TO GET BETTER APPROXIMATION

              ITR = 0
901           CONTINUE
              ITR    = ITR+1
              CH_ANG = .FALSE.
              SNEW   = 0.0
              IQT    = 0
C 
              IF (MOD(ITR,2) .EQ. 1)  THEN
                 K1 = NIMA
                 K2 = 1
                 K3 = -1
              ELSE
                 K1 = 1
                 K2 = NIMA
                 K3 = 1
              ENDIF
              DO  812  IMI=K1,K2,K3
                 IF (MOD(ITR,2) .EQ. 1) THEN
                    KTN = K1-IMI+1
                 ELSE
                    KTN = IMI
                 ENDIF
                 IF (IP(IMI).NE.J) GOTO  812
                 IQT = IQT+1
                 CALL  CROSRNG
     &             (CIRSEED(1,J),CIRC(1,IMI),LCIRC,NRING,TEMP,TEMP(1,2),
     &                  MAXRIN,JACUP,NUMR,PEAK,TOT,MODE)
	         IF (MIRROR)  THEN
                   CALL  CROSRMG
     &             (CIRSEED(1,J),CIRC(1,IMI),LCIRC,NRING,TEMP,TEMP(1,2),
     &                  MAXRIN,JACUP,NUMR,PEAM,TMT,MODE)
                   IF (PEAM .GT. PEAK)  THEN
                      IF (ROT(IMI) .NE. -TMT)  THEN
                         CH_ANG=.TRUE.
                         ROT(IMI)=-TMT
                      ENDIF
                      F = ES(J)+EC(IMI)-2.0*PEAM
                      CALL UPDTM(CIRNEW,CIRC(1,IMI),LCIRC,NRING,NUMR,
     &                           TMT,MAXRIN,IQT)
                      GOTO 152
                   ENDIF
                 ENDIF

                 IF (ROT(IMI) .NE. TOT)  THEN
                    CH_ANG   = .TRUE.
                    ROT(IMI) = TOT
                 ENDIF

                 F=ES(J)+EC(IMI)-2.0*PEAK
                 CALL UPDTC(CIRNEW,CIRC(1,IMI),LCIRC,NRING,NUMR,
     &                      TOT,MAXRIN,IQT)

152              SNEW      = SNEW+F
                 DIST(IMI) = F
812           CONTINUE

c             WRITE(NOUT,5515)  j,itr,snew,CH_ANG
c5515         format(' Group #',i2,'  Iteration #',i4,' SNEW=',1pe13.6,2x,l1)

              IF ((E(J).EQ.0.0.OR.(SNEW.LE.E(J))).AND.CH_ANG)  THEN
                 CIRSEED(:,J) = CIRNEW
                 ES(J) = ENER(CIRNEW,LCIRC,NRING,NUMR,MODE)
                 E(J)  = SNEW
                 GOTO  901
              ENDIF

              E(J)=SNEW
201        CONTINUE

           ZEROGR = .FALSE.
           DO  I=1,NK
              IF (IQ(I) .EQ. 0)  THEN
                 ZEROGR = .TRUE.
                 DMAX   = -1.E30
                 DO  J=1,NIMA
                    IF (DIST(J) .GT. DMAX)  THEN
                       DMAX    = DIST(J)
                       MAX     = J
                       DIST(J) = -1.E30
                    ENDIF
                 ENDDO
                 CIRSEED(:,I)=CIRC(:,MAX)
                 WRITE(NOUT,5992)  I
5992             FORMAT('  WARNING! EMPTY CLUSTER DETECTED #',I3)

                 WRITE(NOUT,5993) MAX
5993             FORMAT('        NEW SEED CREATED FROM OBJECT #',I5)
              ENDIF
           ENDDO
           IF (.NOT. ZEROGR)  THEN

C             CALCULATE  WITHIN SUM OF SQUARES
              WSS = SUM(E)
              WNK = SUM(E/IQ)

              WRITE(NOUT,2060)WSS
2060          FORMAT('  Within sum of squares:',1PE10.3)
              WRITE(NOUT,2008) E
2008          FORMAT('  Within-group sum of squares:',/,5(2X,1PD10.3))
C             WRITE(NOUT,  2001,(ANG(ROT(J),MODE),J=1,NIMA)
C
C             HALTING RULE
C
              IF (ABS(WSS-WSSMAX) / WSSMAX .LT. 0.001)  EXIT
              WSSMAX = WSS
           ENDIF
        ENDDO     !END OF LOOP OVER ITERATIONS

        IF (ITER .GE. MAXIT)
     &     WRITE(NOUT,*) '  MAXIMUM NUMBER OF ITERATIONS REACHED'

        IF (VERBOSE) WRITE(NOUT,*) '  '
        WRITE(NOUT,*) ' NUMBER OF OBJECTS IN CLUSTERS:'
        LK  = NK / 10
        MK  = MOD(NK,10)
        IP1 = 1
        IP2 = 0
        DO J=1,LK+MIN0(MK,1)
           IP2 = MIN0(IP1+9,NK)
           WRITE(NOUT,2092) (I,I=IP1,IP2)
           WRITE(NOUT,2092) (IQ(I),I=IP1,IP2)
           WRITE(NOUT,2092)
           IP1 = IP1+10
        ENDDO 
        WRITE(NOUT,*)' ASSIGNMENT'
        WRITE(NOUT,2002) IP
2002    FORMAT(1X,20I3)

        WRITE(NOUT,  2003)WSS, TOTAL, TOTAL-WSS
2003    FORMAT('  SUM OF SQUARES:',
     &         '  WITHIN:', 1PD10.3,
     &         '  TOTAL:',  1PD10.3,
     &         '  BETWEEN:',1PD10.3)

        WRITE(NOUT,2004)WSS*(TOTAL-WSS),NK*WSS,WNK
2004    FORMAT('  CRITERIA:  W*B,           K*W,    SUM(WK/NK)',/,
     &         4X,3(5X,1PE10.3))
        IF (VERBOSE) WRITE(NOUT,*) '  '

        END

