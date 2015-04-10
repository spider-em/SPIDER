C	
C ++********************************************************************
C
C PKD
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
C                                                                      *
C PKD                                                                  *
C
C SUPPORT_ROUTINE
C
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012                                                                      *
C***********************************************************************

        SUBROUTINE PKD(LUN1,NSAM,NROW,NQ,XYZ,itmp,THRSH,L,NNS,NNR)

        DIMENSION Q(NSAM,-NQ:NQ),G(-NQ:NQ,-NQ:NQ)
        DIMENSION XYZ(3,itmp) 
        LOGICAL T
C
        KM(K)=MOD(K+3*NQ,2*NQ+1)-NQ
        KQ(K)=MOD(K+NROW,NROW)+1
	
        EPSN=1.E-4
        DD=NQ*NQ
        L=0
        AG=0.0
        BG=0.0
        DO    J=-NQ,NQ-1
           CALL  REDLIN(LUN1,Q(1,J),NSAM,KQ(J))
        ENDDO
C       GAUSSIAN OF HALF-WIDTH=NQ
        XN=FLOAT((2*NQ+1)**2)
        XX=FLOAT(NQ*NQ)
        DO  I=-NQ,NQ
           DO  J=-NQ,NQ
              G(I,J)=EXP(-FLOAT(I*I+J*J)/XX)
              AG=AG+G(I,J)
              BG=BG+G(I,J)*G(I,J)
           ENDDO   
        ENDDO
        AVGG=AG/XN
        STDG=SQRT(BG-AG*AG/XN)
        DO  I=-NQ,NQ
           DO  J=-NQ,NQ
              G(I,J)=(G(I,J)-AVGG)/STDG
           ENDDO
        ENDDO
C

        KNSAM=NNS/2+1
        KNROW=NNR/2+1
	
C
        DO    J=1,NROW
           CALL  REDLIN(LUN1,Q(1,KM(J+NQ)),NSAM,KQ(J+NQ-1))
           DO    I=NQ+1,NSAM-NQ
              T=.TRUE.
              AT=0.0
              BT=0.0
              CT=0.0
              Z=Q(I,KM(J))    
              DO  JT=-NQ,NQ
                 JTM=KM(JT+J)
                 DO  IT=-NQ,NQ
                    IF(IT.EQ.0.AND.JT.EQ.0)GO TO 8
                    IF(Z.LT.Q(I+IT,JTM))GOTO 7
8                   CONTINUE
                 ENDDO
              ENDDO         

	      
              DO  JJ=-NQ,NQ
                 JJM=KM(JJ+J)
                 DO  II=-NQ,NQ
                    AT=AT+Q(I+II,JJM)
                    BT=BT+Q(I+II,JJM)*Q(I+II,JJM)
                    CT=CT+G(II,JJ)*Q(I+II,JJM)
                 ENDDO
              ENDDO
              RT=(BT-AT*AT/XN)
              IF(RT.LE.EPSN) THEN
                 COEF=0.0
                 GOTO 7
              ELSE
                 COEF=CT/SQRT(RT)
                 IF(COEF.GE.THRSH)THEN
                    T=.TRUE.
                 ELSE
                    GOTO 7
                 ENDIF
              ENDIF
	      
              IF (T)THEN
                 IF((I-KNSAM).LE.0.OR.(I+KNSAM).GT.NSAM) GOTO 7
                 IF((J-KNROW).LE.0.OR.(J+KNROW).GT.NROW) GOTO 7
                 XX=FLOAT(I)
                 YY=FLOAT(J)
 
                 IF(L.GT.0)  THEN
C     Check whether there were any other peaks in the vicinity
                    DO  LT=1,L
                       IF(((XYZ(1,LT)-XX)**2+
     &                    (XYZ(2,LT)-YY)**2).LE.DD) THEN
                          IF(XYZ(3,LT).LT.Z)  THEN
                             XYZ(1,LT)=XX
                             XYZ(2,LT)=YY
                             XYZ(3,LT)=Z
                          ENDIF
                          GOTO 7
                       ENDIF
                    ENDDO
                 ENDIF
                 L=L+1
                 XYZ(1,L)=XX
                 XYZ(2,L)=YY
                 XYZ(3,L)=Z

              ENDIF
7             CONTINUE
	
           ENDDO

        ENDDO

        END
			
