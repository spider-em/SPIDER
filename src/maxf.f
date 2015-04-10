
C ++********************************************************************
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
C
C   MAXF(NOUT,LUN1,LUN2,NSAM1,NROW1,NSAM2,NROW2,IXB,IYB)
C
C   PARAMETERS:
C         LUN1          LOGICAL UNIT NUMBER
C         LUN2          LOGICAL UNIT NUMBER 
C         NSAM1,NROW1   DIMENSIONS OF IMAGE
C         NSAM2,NROW2   DIMENSIONS OF IMAGE
C         IXB,IYB
C
C--*********************************************************************

      SUBROUTINE MAXF(LUN1,LUN2,NSAM1,NROW1,NSAM2,NROW2,
     &                IXB,IYB,KX,KY,C,FINT,IDUMP,R,CRIT)

C     DANGEROUS TO USE CMBLOCK DUE TO USE OF AV HERE!
      COMMON /UNITS/LUN,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT

      INCLUDE 'CMLIMIT.INC'
      COMPLEX        BUF(1)
      COMMON /IOBUF/ BUFF(NBUFSIZ)
      EQUIVALENCE    (BUF,BUFF)

C     PEAK ARRAY SIZE
      PARAMETER (NREFLM = 200)
      COMPLEX   CDUM
      INTEGER   HINDDUM,KINDDUM
      REAL      FINTDUM,KXDUM,KYDUM
      COMMON    CDUM(NREFLM),HINDDUM(NREFLM),KINDDUM(NREFLM),
     &          FINTDUM(NREFLM),KXDUM(NREFLM),KYDUM(NREFLM),AREA(1)
C     WARNING:  TVWN3 PASSES STUFF IN UNLABELED COMMON, HERE ONLY
C               AREA IS USED DIRECTLY FROM THE COMMON

      COMPLEX  C
      REAL     IMAX,IMAX1

      IYA   = KX+NROW1/2 - IYB/2
      IXA   = KY+NSAM1/2-  IXB/2+1
      IMAX  = 0
      FINT  = 0.0
      NA    = 0
      R2    = R**2
      NWIN2 = IYB*IXB

C     DO LOOP TO FIND MAX OF WINDOW VALUES LOGICALLY SHOULD BE:
C     IYA-1+K AND IXA-1+I.  BUT IT WORKS RIGHT THIS WAY!

      DO  K=1,IYB
         CALL REDLIN(LUN1,BUFF,NSAM1,IYA+K)
         IF (IDUMP.EQ.1) WRITE(NDAT,99)(BUFF(IXA+I),I=1,IXB)
99       FORMAT(1X,12G10.3)
         DO  I=1,IXB
           B = BUFF(IXA+I)
           FINT = FINT+B
           NA = NA+1
           AREA(NA) = B
           IF(B.GT.IMAX) THEN
             IMAX=B
             KX=K
             KY=I
           ENDIF
	ENDDO
      ENDDO

      IF(R.GT.FLOAT(IYB)) GOTO 174

C     COMPUTE AVERAGE OF AREA OUTSIDE REFLECTION
      NAV = 0
      AV = 0.
      DO  K = 1,IYB
        NOFF = (K-1)*IXB
        K2 = (K-KX)**2
        DO  I = 1,IXB
           IF (FLOAT(K2)+FLOAT((I-KY)**2).GE.R2) THEN
              AV = AV +AREA(NOFF+I)
              NAV = NAV+1
           ELSE
              AREA(NOFF+I) = -1000.
           ENDIF
	ENDDO
      ENDDO
      AV = AV / FLOAT(NAV)

C     NOW COMPUTE STANDARD DEVIATION OF OUTSIDE AREA
C     AND APPLY THRESHOLD CRITERION

      SIG = 0.
      DO  I = 1,NWIN2
         IF (AREA(I).GE.-999.) SIG = SIG+(AREA(I)-AV)**2
      ENDDO

      SIG = SQRT(SIG/FLOAT(NAV))
      IF (IMAX-AV.LT.CRIT*SIG) THEN
         FINT = -1.
         RETURN
      ENDIF

C     WRITE(NDAT,700) KX,KY
C174  KX=KX-NROW1/2+IYA
C     KY=KY-NSAM1/2-1+IXA
174   KY=KY-NROW1/2+IYA
      KX=KX-NSAM1/2-1+IXA
      KXF = KX
      KYF = KY
      IF (KX .LT. 0) THEN
         KYF = -KY
         KXF = -KX
      ENDIF

C     FIND LOCATION IN FOURIER FILE
175   CONTINUE
      IF (KYF .GE. 0)  THEN
         ITREC=KYF
      ELSE
	 ITREC=KYF+NROW1+1
      ENDIF
      CALL REDLIN(LUN2,BUF,NSAM2,ITREC)
C     WRITE(NOUT,700) KXF,KYF,ITREC,NADDR

      C = BUF(KX+1)
C     IF (PHF.NE.1.) C=CONJG(C)
      IF (KXF.EQ.0) C=CONJG(C)
      IMAX1 = IMAX-AV
C     WRITE(NDAT,181)C,IMAX,IMAX1,AV,SIG
181   FORMAT(1X,6G12.3)

      RETURN
      END
