C++*********************************************************************
C
C DEFO001.F 
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
C   DEFO001(NUM,NP,KP,NA,NSAM,SPMAX)
C
C   USING LEAST SQUARE METHOD TO DETERMINE DEFOCUS, AMPLITUDE CONTRAST
C   X(K,A)=PI*(0.5*CS*LAMBDA**3*K**4-DZ*LAMBDA*K**2)-OFFSET
C   X(K,A)=PI*(0.5*CS*LAMBDA**3*K**4-A1*LAMBDA*K**2)-A2
C
C       NUM: NUMBER OF IMAGES
C       NP(I): NUMBER OF MINIMUM IN EACH IMAGES
C       KP(I,J): ARRAY OF SP. FREQ. POINTS OF MINIMUM
C       NA(I,J): ARRAY OF ABBERATION
C       NSAM: IMAGE DIMENSION
C       SPMAX: MAX OF SP. FREQ.
C       NA: NUMBER OF ABBERATION IN UNIT OF PI
C       
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/

        SUBROUTINE DEFO001(NUM,NP,KP,NA,NSAM,SPMAX)

C       I DO NOT KNOW IF SAVE IS NEEDED FEB 99 al
        SAVE

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 

        CHARACTER (LEN=MAXNAM) :: OUTNAME
        
        DIMENSION     NP(*)
        REAL          NA(20,20),KP(20,20)

        COMMON        C(20,20),B(20),Y1(10,10),Y2(10,10),A0(20),A1(20)
        DIMENSION     AM0(30),AM2(30),NM(30),XM(30)
        DIMENSION     F2(512),WEIGHT(20,10)

        REAL          KM,KS,KF,LAMBDA
        CHARACTER *1  CHO1,NULL

        DATA PI/3.141592654/

        AA0 = -100.0

        LUN2 = 10

C       INPUT EM PARAMETERS 
        WRITE(NOUT,*)' INPUT PARAMETERS OF IMAGES'
  
        A1(1) = 0

        DO I=2,NUM
           WRITE(NOUT,*) '#',I, '   IMAGE'
           CALL RDPRM(A1(I),NOT_USED,
     &            'DEFOCUS INTERVAL TO FIRST ONE [A]')
        ENDDO

        CALL RDPRM1S(LAMBDA,NOT_USED,'WAVELENGTH LAMBDA [A]',IRTFLG)

        CALL RDPRM1S(CS,NOT_USED,'SPHERICAL ABERRATION CS [MM]',IRTFLG)

        IF (CS < 0.0001) CS = 0.0001
        CS = CS * 1.0E07
        KM = SPMAX
        KS = KM / FLOAT(NSAM)

C       GET VALUE OF Y1   
        DO  I=1,NUM
           DO  J=1,NP(I)
              Y1(I,J)=PI*NA(I,J)
           ENDDO
        ENDDO

        DO  I=1,20
           A0(I) = 0
           B(I)  = 0
           DO  J=1,20
              C(I,J) = 0
           ENDDO
        ENDDO

        DO I=1,NUM
           DO J=1,NP(I)
              KF     = KP(I,J) * KS
              C(1,1) = C(1,1) - PI * (LAMBDA * KF**2)**2
              C(2,1) = C(2,1) - PI * LAMBDA * KF**2
              C(1,2) = C(1,2) - LAMBDA * KF**2
              C(2,2) = C(2,2) - 1.0
              B(1) = B(1) - (PI*(0.5 * CS * LAMBDA**3 * KF**4 -
     &                       A1(I) * LAMBDA * KF**2)-
     &                       Y1(I,J)) * (LAMBDA*KF**2)
              B(2) = B(2)-(PI*(0.5*CS*LAMBDA**3*KF**4 -
     &                       A1(I)*LAMBDA*KF**2) - Y1(I,J))
           ENDDO
        ENDDO

        CALL MATINV(C,2,DET)

        DO  I=1,2
           DO  J=1,2
             A0(I) = A0(I)+C(I,J)*B(J)
           ENDDO
        ENDDO

        A2 = A0(2)

        IF (A0(2) .GT. 0.3 .OR. A0(2) .LT. 0) THEN
C          USING THE DEEPEST GRADIENT PROGRAM 
C          CALCULATE THE WEIGHT OF EACH POINTS
           DO  I=1,NUM
              DO  J=1,NP(I)
                 KF = KP(I,J) * KS
                 WEIGHT(I,J) = PI*(2.*CS*LAMBDA**3*KF**3-2.*
     &                        (A0(1)+A1(I))*LAMBDA*KF)*2.*KS
              ENDDO
           ENDDO
     
           DO  K=1,30
              A0(1) = AA0
              A2    = FLOAT(K)*0.01

C             SET ITERATION STEP
              NSTEP = 0

C             CALCULATE VALUE OF Y2
              DO I=1,NUM
                 DO J=1,NP(I)
                    KF      = KP(I,J)*KS
                    Y2(I,J) = PI*(0.5*CS*LAMBDA**3*KF**4-(A0(1)+A1(I))*
     &                       LAMBDA*KF**2)-A2
                 ENDDO
              ENDDO

C             SET INITIAL VALUE X**2(X,A)=SUM((Y1(I)-Y2(I))**2) 
              X1 = 0
              DO  I=1,NUM
                 DO  J=1,NP(I)
                    X1=X1+(Y1(I,J)-Y2(I,J))**2/WEIGHT(I,J)**2
                 ENDDO
              ENDDO

C             CALCULATE THE VALUE OF Y2(I) 
999           DA0 = 0.001*A0(1)
              DA2 = 0.001*A2

              DO  I=1,NUM
                 X = A0(1)+0.1*DA0+A1(I)
                 DO  J=1,NP(I)
                    KF = KP(I,J)*KS
                    Y2(I,J) = PI*(0.5*CS*LAMBDA**3*KF**4 -
     &                        X*LAMBDA*KF**2)-A2
                 ENDDO
              ENDDO

C             CALCULATE DX**2/DA0   
              X2 = 0
              DO  I=1,NUM
                 DO  J=1,NP(I)
                    X2=X2+(Y1(I,J)-Y2(I,J))**2/WEIGHT(I,J)**2
                 ENDDO
              ENDDO

              DXA0 = (X2-X1)/(0.1*DA0)

C             CALCULATE DX**2/DA2 
              DO  I=1,NUM
                 X = A0(1) + A1(I)
                 DO  J=1,NP(I)
                    KF = KP(I,J)*KS
                 Y2(I,J) = PI*(0.5*CS*LAMBDA**3*KF**4-X*LAMBDA*KF**2)-
     &                   (A2+0.1*DA2)
                 ENDDO
              ENDDO

              X2 = 0
              DO  I=1,NUM
                 DO  J=1,NP(I)
                    X2 = X2+(Y1(I,J)-Y2(I,J))**2/WEIGHT(I,J)**2
                 ENDDO
              ENDDO

              DXA2 = (X2-X1)/(0.1*DA2)

C             CALCULATE THE SUM=SUM((DX**2/DAI*DAI)**2) */
              SUM   = SQRT((DXA0*DA0)**2+(DXA2*DA2)**2)
              A0(1) = A0(1)-DXA0*DA0**2/SUM
              A2    =  A2-DXA2*DA2**2/SUM

C             CRITERI FOR ITERATION 
              DO  I=1,NUM
                 X = A0(1)+A1(I)
                 DO  J=1,NP(I)
                    KF = KP(I,J)*KS
                 Y2(I,J)=PI*(0.5*CS*LAMBDA**3*KF**4-X*LAMBDA*KF**2)-A2
                 ENDDO
              ENDDO

              X2 = 0
              DO  I=1,NUM
                 DO  J=1,NP(I)
                    X2 = X2 + (Y1(I,J) - Y2(I,J))**2 / WEIGHT(I,J)**2
                 ENDDO
              ENDDO
              IF ((X1-X2) .LT. 0) THEN
                 A0(1) = A0(1) + 0.5 * DXA0 * DA0**2 / SUM
                 A2    = A2 + 0.5*DXA2*DA2**2/SUM

c                WRITE(NOUT,*) 'INITIAL A2=',FLOAT(K)*0.01,'  STEP',NSTEP
c                WRITE(NOUT,*) 'A0=',A0,'OFFSET(RAD)=',A2,'X**2=',X1
                 Q      = SIN(A2)/COS(A2)*100.
                 XM(K)  = X1
                 AM0(K) = A0(1)
                 AM2(K) = A2
                 NN     = INT(A2/0.01+0.5)
                 NM(NN) = NM(NN)+1
              ELSE
C                SET PARAMETERS FOR NEXT STEP */
                 X1    = X2
                 NSTEP = NSTEP+1
                 GOTO 999
              ENDIF
           ENDDO

C          FIND CONVERGE POINTS
           NN = 0
           DO I=1,30
              IF(NM(I) .GT. 5) NN = NN + 1
           ENDDO

           IF (NN .GT. 1) THEN
C             THERE ARE TWO CONVERGE POINTS
              X0 = 999999999.0
              DO I=1,30
                 IF(NM(I) .GT. 5 .AND. XM(I) .LT. X0) THEN
                    NN=I
                 X0 = XM(I)       
                 ENDIF
              ENDDO
           ELSE
C             THERE IS ONLY ONE CONVERGE POINT
              NN = 0
              DO I=1,30
                 IF (NM(I) .GT. NN .OR. (NM(I) .EQ. NN .AND. 
     &               XM(I) .LT. XM(NN))) NN=I
              ENDDO
           ENDIF

           A0(1) = AM0(NN)
           A2    = AM2(NN)
        ENDIF

        WRITE(NOUT,*)' DEFOCUS ARE [A]',(A0(1)+A1(I), I=1,NUM)

        WRITE(NOUT,140) A2
140     FORMAT('  AMPLITUDE CONTRAST:',F10.6)

C       GENERATE A FILTER FILE FROM CTF W/O ENVELOPE FUNCTION
C       GATE VALUE IS 0.08

        GATE = 0.08
        CALL RDPRMC(CHO1,NUMC,.TRUE.,
     &      'DO YOU WANT TO GENERATE A FILTER? (Y/N)',NULL,IRT)

        IF ( CHO1 .EQ. 'Y') THEN

            DO I=1,NUM
               X = A0(1)+A1(I)

               DO  J=1,NSAM
                  KF = FLOAT(J)*KS
                  X1 = PI*(0.5*CS*LAMBDA**3*KF**4-X*LAMBDA*KF**2)-A2
                  F1 = SIN(X1)
                  IF (F1 .GT. GATE) THEN
                     F2(J) = 1
                  ELSE
                     IF (F1 .LT. -1.* GATE) THEN
                        F2(J) = -1
                     ELSE
                        F2(J) = 0
                     ENDIF
                  ENDIF
               ENDDO

               IFORM = 1
               MAXIM = 0
               CALL OPFILEC(0,.TRUE.,OUTNAME,LUN2,'U',IFORM,NSAM,1,1,
     &               MAXIM,'OUTPUT',.FALSE.,IRTFLG)

               CALL WRTLIN(LUN2,F2,NSAM,1)

               CLOSE(LUN2)
            ENDDO
        ENDIF
                
        END
