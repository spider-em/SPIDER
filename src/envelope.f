C **********************************************************************
C
C ENVELOPE.F
C                   OPFILEC                       FEB 03 ARDEAN LEITH
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
C  ENVELOPE(IRTFLG) : CALCULATE THE ENVELOPE OF CTF 
C  USING LEAST SQUARE METHOD(DEEPEST GRADIENT) TO FIT EXPERIMENT PROFILES
C  (DELETED NOISE BACKGROUND)
C  FITTING THE FOURIER COEFFICIENTS WITH EQUATION 
C  f(A1-A4)=A1*SIN(X(K,DZ,CS,A2)*E1*E2*E3
C  X(K,DZ,CS,Q)=PI*(0.5*CS*LAMBDA**3*K**4-DZ*LAMBDA*KF**2)-Q
C  E1=EXP(-2.*PI**2*A2**2*(KF**3*CS*LAMBDA**3-DZ*KF*LAMBDA)**2)
C  A2: SOURCESIZE
C  E2=EXP(-PI**2*A3**2*K**4*LAMBDA**2/(16*LN2))
C  A3: DEFOCUS SPREAD
C  E3=1/[1+(KF/KFILM)**2]
C  KFILM: CHARACTER OF FILM
C  A4: GAUSSIAN ENVOLOPE HALFWIDTH
C  E4=EXP(-(KF/A4)**2)
C
C  NUM: NUMBER OF IMAGE
C  NSAM: IMAGE DIMESION
C  SPMAX: MAX. OF SPATIAL FREQUENCE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	SUBROUTINE ENVELOPE(IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

        COMMON       NF(10,2),Y1(10,512),Y2(512),Y3(512),
     &               Y4(512,10),Y5(512),
     &               A1(10),DA1(10,10),DZ(10),DXA1(10)

	CHARACTER*1  CHOICE,C2,NULL
	REAL         KM,KS,KF,LAMBDA,KFILM
	LOGICAL      FLAG
        CHARACTER(LEN=MAXNAM)              ::   FILNAME,OUTNAME

	DATA PI/3.1415926535/

	LUN1   = 8
	LUN2   = 10
	IRTFLG = 0
	NULL   = CHAR(0)
5	WRITE(NOUT,*)
     &      ' FITTING THE PROFILE INTO f(A1,A2,A3,A4)=A1*SIN(X(KF))*E1',
     &      '(A2)*E2(A3)*E3(KFILM)*E4(A4)'
 	WRITE(NOUT,*)
     &    ' SIN(X(K))=SIN(PI*(0.5*CS*LAMBDA**3*KF**4-DZ*LAMBDA*KF**2-Q)'
	WRITE(NOUT,*)
     &    ' E1=EXP(-2*PI**2*A2**2*(K**3*CS*LAMBDA**3-DZ*K*LAMBDA)**2 )'
	WRITE(NOUT,*)' E2=EXP(-PI**2*A3**2*K**4*LAMBDA**2/16LN2)'
	WRITE(NOUT,*)' E3=1/[1+(KF/KFILM)**2]'
	WRITE(NOUT,*)' E4=EXP(-(KF/A4)**2]'
	CALL RDPRMI(NUM,NDUM,NOT_USED,
     &              'HOW MANY IMAGES IN THE SERIES?')
	WRITE(NOUT,*)' INPUT EM PARAMETERS'
	CALL RDPRM(LAMBDA,NOT_USED,'WAVELENGTH[A]=')
	CALL RDPRM(CS,NOT_USED,'SPHERICAL ABERRATION CS[MM]')
           IF (CS < 0.0001)    CS = 0.0001
	CS=CS*1.E7
	CALL RDPRM(SPMAX,NOT_USED,'MAX. SP. FREQ.[A-1]=')
	CALL RDPRM(AC,NOT_USED,'AMPLITUDE CONTRAST [RAD]=')
	CALL RDPRM(A2,NOT_USED,'SOURCE SIZE [A-1]=')
	CALL RDPRM(A3,NOT_USED,'DEFOCUS SPREAD [A]=')
	CALL RDPRM(KFILM,NOT_USED,'CHARACTER OF THE FILM Kf [A-1]=')
	CALL RDPRM(A4,NOT_USED,'GAUSSIAN ENVELOPE CHARACTER [A-1]=')
	DO 20 I=1,NUM
	   WRITE(NOUT,*)'# ',I,'   IMAGE'
           MAXIM = 0
	   CALL OPFILEC(0,.TRUE.,FILNAME,LUN1,'O',IFORM,
     &          NSAM,NROW,NSLICE,MAXIM,'IMAGE', .FALSE.,IRTFLG)
	   IF (IRTFLG .NE. 0) RETURN

	   WRITE(NOUT,24) NSAM,NROW
24	   FORMAT(' FILE"S DIMENSION:',I4,' X',I4)
	   CALL RDPRM(DZ(I),NOT_USED,'DEFOCUS [A]=')
C	   CALL RDPRM(A1(I),NOT_USED,'AMPLITUDE OF PHASE CONTRAST=')
C	   CALL RDPRM(A4(I),NOT_USED,'GAUSSIAN ENVELOPE CHARACTER =')
	   CALL REDLIN(LUN1,Y3,NSAM,1)
	   DO K=1,NSAM
	      Y1(I,K)=Y3(K)
	   ENDDO
	   CALL RDPRMI(NF(I,1),NF(I,2),NOT_USED,
     &                 'FITTING REGION (N1-N2):$')
	   CLOSE(LUN1)
C          GET AMPLITUDE OF PHASE CONTRAST
	   DO  J=10,NSAM
	      FLAG=.TRUE.
	      IF(Y3(J) .GT. Y3(J+1) .AND. Y3(J) .GT. Y3(J-1)) THEN
	         DO K=J-1,J-10,-1
	            IF(Y3(J) .LT. Y3(K)) FLAG=.FALSE.
	         ENDDO
	         DO K=J+1,J+10
	            IF(Y3(J) .LT. Y3(K)) FLAG=.FALSE.
	         ENDDO
	         IF (FLAG) THEN
	            A1(I)=Y3(J)*(1+5000/DZ(I))
	            GOTO 20
	         ENDIF
	      ENDIF
	   ENDDO

20	   CONTINUE
	   KS=SPMAX/FLOAT(NSAM)
	   DO  I=1,NUM
	      DO  J=1,NUM
	         DA1(I,J)=0
	      ENDDO
	   ENDDO
C.......   GET Y4 VALUE
	   DO I=1,NUM
	      DO J=1,NSAM
	         Y4(J,I)=1
	      ENDDO
	   ENDDO

	   FLAG=.TRUE.
	   NEIB=3
	   DO  I=1,NUM
	      DO  J=1+NEIB,NSAM-NEIB
	         BX=Y1(I,J)
	         IF (BX .GT. Y1(I,J+1)) THEN
	            DO K=J,J-NEIB,-1
	              IF(BX .LT. Y1(I,K)) THEN
	                 FLAG=.FALSE.
	              ENDIF
                 ENDDO
	         DO K=J,J+NEIB
                   IF(BX .LT. Y1(I,K)) THEN
	              FLAG=.FALSE.
	           ENDIF
	         ENDDO
	        IF(FLAG) THEN
	           Y4(J-3,I)=10
	           Y4(J-2,I)=10
	           Y4(J-1,I)=10
	           Y4(J,I)=10
	           Y4(J+1,I)=10
	           Y4(J+2,I)=10
	           Y4(J+3,I)=10
	        ELSE
	           FLAG=.TRUE.
	        ENDIF	
	      ENDIF
	   ENDDO
	ENDDO

C......GET VALUE OF Y2
	X=0.0
	DO I=1,NUM
	CALL Y2VALUE(Y2,NSAM,KS,CS,LAMBDA,DZ(I),AC,KFILM,A1(I),A2,A3,
     &  A4)
C......GET VALUE OF X0=SUM((Y1-Y2)**2)
	DO K=1,NSAM
	Y3(K)=Y1(I,K)
	y5(k)=y4(k,i)
	ENDDO
	CALL XSUM(X0,Y3,Y2,y5,NF(I,1),NF(I,2))
	X=X+X0
	ENDDO
	X0=X
C......SET INITIATION STEP
	NSTEP=0
999	DO I=1,NUM
	DA1(I,I)=0.001*A1(I)
	ENDDO
	DA2=0.001*A2
	DA3=0.001*A3
	DKFILM=0.001*KFILM	
	DA4=0.001*A4
	DAC=0.001*AC
C......GET VALUE OF dX**2/dA1(I)
	 DO  I=1,NUM
	 X=0
	 DO  J=1,NUM
	  CALL Y2VALUE(Y2,NSAM,KS,CS,LAMBDA,DZ(J),AC,KFILM,
     &    A1(J)+0.1*DA1(I,J),A2,A3,A4)
	   DO K=1,NSAM
	    Y3(K)=Y1(J,K)
	    y5(k)=y4(k,j)
	   ENDDO
	  CALL XSUM(X1,Y3,Y2,y5,NF(J,1),NF(J,2))
          X=X+X1
	 ENDDO
	 X1=X
	DXA1(I)=-(X1-X0)/(0.1*DA1(I,I))
	ENDDO
C......GET VALUE OF dX**2/dA2
	X=0
	DO  I=1,NUM
	CALL Y2VALUE(Y2,NSAM,KS,CS,LAMBDA,DZ(I),AC,KFILM,A1(I),
     &               A2+0.1*DA2, A3,A4)
	 DO K=1,NSAM
	  Y3(K)=Y1(I,K)
	  y5(k)=y4(k,i)
	 ENDDO
	CALL XSUM(X1,Y3,Y2,y5,NF(I,1),NF(I,2))
	X=X+X1
	ENDDO
	X1=X
	DXA2=-(X1-X0)/(0.1*DA2)
C......GET VALUE OF dX**2/dA3
	X=0
	DO  I=1,NUM
	CALL Y2VALUE(Y2,NSAM,KS,CS,LAMBDA,DZ(I),AC,KFILM,A1(I),A2,
     &  A3+0.1*DA3,A4)
	  DO K=1,NSAM
	  Y3(K)=Y1(I,K)
	  y5(k)=y4(k,i)
	  ENDDO
	CALL XSUM(X1,Y3,Y2,y5,NF(I,1),NF(I,2))
	X=X+X1
	ENDDO
	X1=X
C	DXA3=-(X1-X0)/(0.1*DA3)
	DXA3=0

C......GET VALUE dX**2/dKFILM
	X=0
	DO  I=1,NUM
	CALL Y2VALUE(Y2,NSAM,KS,CS,LAMBDA,DZ(I),AC,
     &               KFILM+0.1*DKFILM,A1(I),A2, A3,A4)
	  DO K=1,NSAM
	  Y3(K)=Y1(I,K)
	  y5(k)=y4(k,i)
	  ENDDO
	CALL XSUM(X1,Y3,Y2,y5,NF(I,1),NF(I,2))
	X=X+X1
	ENDDO
	X1=X
	DXKFILM=-(X1-X0)/(0.1*DKFILM)
	

C......GET VALUE dX**2/dA4

	X=0
	DO  J=1,NUM
	CALL Y2VALUE(Y2,NSAM,KS,CS,LAMBDA,DZ(J),AC,KFILM,A1(J),A2,A3,
     &  A4+0.1*DA4)
	  DO K=1,NSAM
	   Y3(K)=Y1(J,K)
	   y5(k)=y4(k,j)
	  ENDDO
	CALL XSUM(X1,Y3,Y2,y5,NF(J,1),NF(J,2))
	X=X+X1
	ENDDO
	X1=X
140	DXA4=-(X1-X0)/(0.1*DA4)

C......GET VALUE dX**2/dAC

	X=0
	DO  J=1,NUM
	CALL Y2VALUE(Y2,NSAM,KS,CS,LAMBDA,DZ(J),AC+0.1*DAC,KFILM,A1(J),
     &  A2,A3,A4)
	  DO K=1,NSAM
	   Y3(K)=Y1(J,K)
 	   y5(k)=y4(k,j)
	  ENDDO
	CALL XSUM(X1,Y3,Y2,y5,NF(J,1),NF(J,2))
	X=X+X1
	ENDDO
	X1=X
	DXAC=-(X1-X0)/(0.1*DAC)

C......GET NEW VALUE OF A1,A2,A3,KFILM,A4,AC
	SUM=0
	DO I=1,NUM
	SUM=SUM+(DXA1(I)*DA1(I,I))**2
	ENDDO
	
	SUM=SUM+(DXA4*DA4)**2
	SUM=SUM+(DXA2*DA2)**2+(DXA3*DA3)**2+(DXKFILM*DKFILM)**2
	SUM=SUM+(DXAC*DAC)**2
	SUM=SQRT(SUM)
	DO I=1,NUM
	A1(I)=A1(I)+DXA1(I)*DA1(I,I)**2/SUM
	ENDDO
	A2=A2+DXA2*DA2**2/SUM
	A3=A3+DXA3*DA3**2/SUM
	KFILM=KFILM+DXKFILM*DKFILM**2/SUM
	A4=A4+DXA4*DA4**2/SUM
	AC=AC+DXAC*DAC**2/SUM
	
C......GET NEW VALUE OF Y2
	X=0
	DO I=1,NUM
	CALL Y2VALUE(Y2,NSAM,KS,CS,LAMBDA,DZ(I),AC,KFILM,
     &  A1(I),A2,A3,A4)
C......GET VALUE OF X1=SUM((Y1-Y2)**2)
	  DO K=1,NSAM
	  Y3(K)=Y1(I,K)
	  y5(k)=y4(k,i)
	  ENDDO
	CALL XSUM(X1,Y3,Y2,y5,NF(I,1),NF(I,2))
	X=X+X1
	ENDDO
	X1=X               
C......CRITERIA FOR ITERATION
	D=X0-X1
	IF(D .LT. 0.) THEN
	DO I=1,NUM
	A1(I)=A1(I)-0.5*DXA1(I)*DA1(I,I)**2/SUM
	ENDDO
	A2=A2-0.5*DXA2*DA2**2/SUM
	A3=A3-0.5*DXA3*DA3**2/SUM
	KFILM=KFILM-0.5*DXKFILM*DKFILM**2/SUM
	A4=A4-0.5*DXA4*DA4**2/SUM
	AC=AC-0.5*DXAC*DAC**2/SUM	
C	WRITE(NOUT,*)'STEP=',NSTEP
	WRITE(NOUT,210)
210	FORMAT(' A1=')
	DO I=1,NUM
	   WRITE(NOUT,200)A1(I)
	ENDDO
200	FORMAT(4F16.10)
	WRITE(NOUT,*)'SOURCE SIZE=',A2,'  KFILM=',KFILM

	WRITE(NOUT,220)
220	FORMAT(' GAUSSIAN ENVELOPE HALFWIDTH=')
	WRITE(NOUT,200)A4

	WRITE(NOUT,*)'AMPLITUDE CONTRAST=',AC	
	WRITE(NOUT,*)'X**2=',X1

C......CREATE A DIFFRACTION  PATTERN
1200	WRITE(NOUT,250)
250	FORMAT(' CREATE A POWER SPECTRUM')
	DO I=1,NUM
	  WRITE(NOUT,*)'FILE #', I

          MAXIM = 0
          CALL OPFILEC(0,.TRUE.,OUTNAME,LUN2,'U',IFORM,NSAM,NROW,1,
     &             MAXIM,'OUTPUT',.FALSE.,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

	  DO  K=1,NSAM
	    KF=FLOAT(K)*KS
	    O1=PI*(0.5*CS*LAMBDA**3*KF**4-DZ(I)*LAMBDA*KF**2)-AC
	    O2=-2.0*PI**2*A2**2*(KF**3*CS*LAMBDA**3-DZ(I)*KF*LAMBDA)**2
	    O3=-1.*PI**2*A3**2*KF**4*LAMBDA**2/(16.*ALOG(2.))
	    O4=1./(1.+(KF/KFILM)**2)
	    O5=EXP(-(KF/A4)**2)
	    Y3(K)=ABS(A1(I)*SIN(O1)*EXP(O2+O3)*O4*O5)
	  ENDDO
          CALL WRTLIN(LUN2,Y3,NSAM,1)
	  CLOSE(LUN2)
	ENDDO
	ELSE
	NSTEP=NSTEP+1
	X0=X1
C	WRITE(NOUT,*)'SPET',NSTEP
C	WRITE(NOUT,210)
C	DO I=1,NUM
C	WRITE(NOUT,200)A1(I)
C	ENDDO
C	WRITE(NOUT,*)'A2=',A2,'  A3=',A3
C	WRITE(NOUT,220)
C	DO I=1,NUM
C	WRITE(NOUT,200)A4(I)
C	ENDDO
	GOTO 999
	ENDIF
	RETURN
	END

        SUBROUTINE Y2VALUE(Y2,NSAM,KS,CS,LAMBDA,DZ,Q,KFILM,
     &	A1,A2,A3,A4)
	DIMENSION Y2(*)
	REAL KF,KS,LAMBDA,KFILM
	DATA PI/3.1415926535/
	DO  I=1,NSAM
        KF=FLOAT(I)*KS
	O1=PI*(0.5*CS*LAMBDA**3*KF**4-DZ*LAMBDA*KF**2)-Q
	O2=-2.0*PI**2*A2**2*(KF**3*CS*LAMBDA**3-DZ*KF*LAMBDA)**2
	O3=-1.*PI**2*A3**2*KF**4*LAMBDA**2/(16.*ALOG(2.))
	O4=1./(1.+(KF/KFILM)**2)
	O5=EXP(-(KF/A4)**2)
	Y2(I)=ABS(A1*SIN(O1)*EXP(O2+O3)*O4*O5)
	ENDDO
        RETURN
	END


	SUBROUTINE XSUM(X,Y1,Y2,y3,NM,NSAM)
	DIMENSION Y1(*),Y2(*),y3(512)
	X=0
	DO  I=NM,NSAM
	X=X+(y3(I)*(Y1(I)-Y2(I)))**2
	ENDDO
	RETURN
	END
