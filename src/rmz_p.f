C ++********************************************************************
C                                                                      *
C  RMZ_P                                                               *
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
C       REMEZ ALGORITHM
C
C **********************************************************************

	SUBROUTINE  RMZ_P(H,NFILT)

	INCLUDE 'CMBLOCK.INC'

        COMMON TDM(2000),
     &	   PI2,AD,DEVI,X,Y,RR(2112),WT,ALPHA,IEXT,NFCNS,NGRID 
     &	   ,LB,EDGE,FX,WTX
	CHARACTER*40      BL
        DIMENSION         IEXT(66),AD(66),ALPHA(66),X(66),Y(66)
        DIMENSION         H(66)
        DIMENSION         DES(1056),GRID(1056),WT(1056)
        DIMENSION         EDGE(20),FX(10),WTX(10),DEVIAT(10) 
	CHARACTER*4       IBT
        EQUIVALENCE       (RR,GRID),(RR(1057),DES)
        EQUIVALENCE       (DEVIAT,WTX)
        DOUBLE PRECISION  AD,DEVI,X,Y
        DOUBLE PRECISION  PI2,PI

	NOUTB=NOUT

        PI=3.1415926539793D0
        PI2=6.283185307179586D0
C       DIMENSION OF ROT MATRIX IN COMMON DEPENDS ON NFMAX (=NFMAX/2+1=29)
        NFMAX=57 
	JTYPE=1
        IF (NFILT.GT.NFMAX.OR.NFILT.LT.3)  THEN
	   CALL ERRT(101,'ERROR IN INPUT DATA',IDUM)
	   RETURN
	ENDIF                         
C       IF (LGRID.LE.0 .OR. LGRID.GT.32)  LGRID=32
	LGRID=32
	CALL  RDPRMI(NBANDS,IDUM,NOT_USED,'NUMBER OF BANDS')
        IF (NBANDS.LE.0)  NBANDS=1 
        JB=2*NBANDS 
	DO J=1,JB,2
	  WRITE(IBT,2202)  J/2+1
 2202	  FORMAT(I2)
     	  BL='BAND #'//IBT(1:2)//' - LOWER AND UPPER EDGES$'
 	  CALL  RDPRM2(EDGE(J),EDGE(J+1),NOT_USED,BL)
	ENDDO
	DO    J=1,NBANDS
	  WRITE(IBT,2202)  J
     	  BL='BAND #'//IBT(1:2)//' - DESIRED VALUE$'
 	  CALL  RDPRM2(FX(J),FX(J+1),NOT_USED,BL)
	ENDDO
	DO    J=1,NBANDS
	  WRITE(IBT,2202)  J
    	  BL='BAND #'//IBT(1:2)//' - WEIGHTING$'
 	  CALL  RDPRM2(WTX(J),WTX(J+1),NOT_USED,BL)
	ENDDO
          NEG=1 
          IF(JTYPE.EQ.1)  NEG=0 
          NODD=NFILT/2
          NODD=NFILT-2*NODD 
          NFCNS=NFILT/2 
          IF(NODD.EQ.1.AND.NEG.EQ.0)  NFCNS=NFCNS+1 
          GRID(1)=EDGE(1) 
          DELF=LGRID*NFCNS
          DELF=0.5/DELF 
          IF(NEG.EQ.0)  GO TO 135 
          IF(EDGE(1).LT.DELF)  GRID(1)=DELF 
 135          CONTINUE
          J=1 
          L=1 
          LBAND=1 
 140      FUP=EDGE(L+1) 
 145      TEMP=GRID(J)
          DES(J)=EFF(TEMP,FX,WTX,LBAND,JTYPE) 
          WT(J)=WATE(TEMP,FX,WTX,LBAND,JTYPE) 
          J=J+1 
          GRID(J)=TEMP+DELF 
          IF(GRID(J).GT.FUP)  GO TO  150
          GO TO 145 
 150          GRID(J-1)=FUP 
          DES(J-1)=EFF(FUP,FX,WTX,LBAND,JTYPE)
          WT(J-1)=WATE(FUP,FX,WTX,LBAND,JTYPE)
          LBAND=LBAND+1 
          L=L+2 
          IF(LBAND.GT.NBANDS)  GO TO  160 
          GRID(J)=EDGE(L) 
          GO TO 140 
 160          NGRID=J-1 
          IF(NEG.NE.NODD)  GO TO 165
          IF(GRID(NGRID).GT.(0.5-DELF))  NGRID=NGRID-1
 165          CONTINUE
          IF(NEG)   170,170,180 
 170          IF(NODD.EQ.1)  GO TO 200
          DO    J=1,NGRID
             CHANGE=DCOS(PI*GRID(J)) 
             DES(J)=DES(J)/CHANGE
             WT(J)=WT(J)*CHANGE
	  ENDDO
          GO TO 200 
 180          IF(NODD.EQ.1)  GO TO 190
          DO    J=1,NGRID
             CHANGE=DSIN(PI*GRID(J)) 
             DES(J)=DES(J)/CHANGE
             WT(J)=WT(J)*CHANGE
	  ENDDO
          GO TO 200 
 190          DO    J=1,NGRID
          CHANGE=DSIN(PI2*GRID(J))
          DES(J)=DES(J)/CHANGE
          WT(J)=WT(J)*CHANGE
	      ENDDO
 200          TEMP=FLOAT(NGRID-1)/FLOAT(NFCNS)
          DO    J=1,NFCNS
             IEXT(J)=IFIX(FLOAT(J-1)*TEMP+1.0) 
	  ENDDO
          IEXT(NFCNS+1)=NGRID 
          NM1=NFCNS-1 
          NZ=NFCNS+1
          CALL  REMEZ(EDGE,NBANDS,NOUT)
          IF(NEG)  300,300,320
 300          IF(NODD.EQ.0)  GO TO 310
          DO    J=1,NM1
              NZMJ=NZ-J 
              H(J)=0.5*ALPHA(NZMJ)
	  ENDDO
          H(NFCNS)=ALPHA(1) 
          GO TO 350 
 310          H(1)=0.25*ALPHA(NFCNS)
          DO    J=2,NM1
             NF2J=NFCNS+2-J
             NZMJ=NZ-J 
             H(J)=0.25*(ALPHA(NZMJ)+ALPHA(NF2J)) 
	  ENDDO
          H(NFCNS)=0.5*ALPHA(1)+0.25*ALPHA(2) 
          GO TO 350 
 320          IF(NODD.EQ.0)  GO TO 330
          H(1)=0.25*ALPHA(NFCNS)
          H(2)=0.25*ALPHA(NM1)
          DO    J=3,NM1
             NZMJ=NZ-J 
             NF3J=NFCNS+3-J
             H(J)=0.25*(ALPHA(NZMJ)-ALPHA(NF3J)) 
	  ENDDO
          H(NFCNS)=0.5*ALPHA(1)-0.25*ALPHA(3) 
          H(NZ)=0.0 
          GO TO 350 
 330          H(1)=0.25*ALPHA(NFCNS)
          DO    J=2,NM1
          NZMJ=NZ-J 
          NF2J=NFCNS+2-J
           H(J)=0.25*(ALPHA(NZMJ)-ALPHA(NF2J)) 
	  ENDDO
          H(NFCNS)=0.5*ALPHA(1)-0.25*ALPHA(2) 
 350          WRITE(NOUT,360)
 360          FORMAT(1X,132('*')//25X,'FINITE IMPULSE RESPONSE (FIR)'/
     125X,'LINEAR PHASE DIGITAL FILTER DESIGN'/ 
     225X,'REMEZ EXCHANGE ALGORITHM'/)
          IF(JTYPE.EQ.1)   WRITE(NOUT,365) 
 365          FORMAT(25X,'BANDPASS FILTER'/)
          IF(JTYPE.EQ.2)  WRITE(NOUT,370)
 370          FORMAT(25X,'DIFFERENTIATOR'/) 
          IF(JTYPE.EQ.3)   WRITE(NOUT,375) 
 375          FORMAT(25X,'HILBERT TRANSFORMER'/)
          WRITE(NOUT,378)   NFILT
 378          FORMAT(//15X,'FILTER LENGTH = ',I3)
          WRITE(NOUT,380)
 380          FORMAT(15X,'***** IMPULSE RESPONSE *****')
          DO    J=1,NFCNS
          K=NFILT+1-J 
          IF(NEG.EQ.0)   WRITE(NOUT,382)  J,H(J),K 
          IF(NEG.EQ.1)   WRITE(NOUT,383)  J,H(J),K 
	  ENDDO
 382          FORMAT(20X,'H(',I3,') = ',E15.8,' = H(',I4,')') 
 383          FORMAT(20X,'H(',I3,') = ',E15.8,' =-H(',I4,')') 
          IF(NEG.EQ.1.AND.NODD.EQ.1)   WRITE(NOUT,384)   NZ
 384          FORMAT(20X,'H(',I3,') = 0.0') 
          DO  450  K=1,NBANDS,4 
          KUP=K+3 
          IF(KUP.GT.NBANDS)  KUP=NBANDS 
          WRITE(NOUT,385)  (J,J=K,KUP) 
 385          FORMAT(/24X,4('BAND',I3,8X))
          WRITE(NOUT,390)  (EDGE(2*J-1),J=K,KUP) 
 390          FORMAT(2X,'LOWER BAND EDGE',5F15.9) 
          WRITE(NOUT,395)   (EDGE(2*J),J=K,KUP)
 395          FORMAT(2X,'UPPER BAND EDGE',5F15.9) 
          IF(JTYPE.NE.2)   WRITE(NOUT,400)  (FX(J),J=K,KUP)
 400          FORMAT(2X,'DESIRED VALUE',2X,5F15.9)
          IF(JTYPE.EQ.2)   WRITE(NOUT,405)   (FX(J),J=K,KUP) 
 405          FORMAT(2X,'DESIRED SLOPE',2X,5F15.9)
          WRITE(NOUT,410)   (WTX(J),J=K,KUP) 
 410          FORMAT(2X,'WEIGHTING',6X,5F15.9)
          DO    J=K,KUP
            DEVIAT(J)=DEVI/WTX(J)
	  ENDDO
          WRITE(NOUT,425)  (DEVIAT(J),J=K,KUP) 
 425          FORMAT(2X,'DEVIATION',6X,5F15.9) 
          IF(JTYPE.NE.1)  GO TO 450 
          DO    J=K,KUP
             DEVIAT(J)=20.0*ALOG10(DEVIAT(J))
	  ENDDO
          WRITE(NOUT,435)   (DEVIAT(J),J=K,KUP)
 435          FORMAT(2X,'DEVIATION IN DB',5F15.9)
 450          CONTINUE
C          DO  456  J=1,NZ
C          JJ=IEXT(J)
C 456      WT(J)=GRID(JJ)
C          WRITE(NOUT,455)  (WT(J),J=1,NZ)
C 455          FORMAT(/2X,'EXTREMAL FREQUENCIES'/,(2X,5F12.7)) 
 458      CONTINUE
          END 

