C
C++*********************************************************************
C
C $$ PW3SDR.FOR
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
C IMAGE_PROCESSING_ROUTINE
C
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************
C
C $$ PW3SDR.FOR
C
        SUBROUTINE  PW3SDR(LUN1,LUN2,NSAM,NROW,NSLICE,MODE)
        DIMENSION   Q(NSAM+2-MOD(NSAM,2)),T(NSAM)
        CHARACTER*1  MODE
        LOGICAL*1 IFND,IFNL
        
        
        NNNN=NSAM+2-MOD(NSAM,2)
        IFND=MOD(NROW,2).EQ.0
        IFNL=MOD(NSLICE,2).EQ.0
        NSC=NSAM/2+1

        SCL=2.0/FLOAT(NSAM)/FLOAT(NROW)/FLOAT(NSLICE)
        DO    K=1,NSLICE
           DO    J=1,NROW
              NR=J+(K-1)*NROW
              CALL  REDLIN(LUN1,Q,NNNN,NR)
C  Reaplce F(0,0,0) for visualization
	IF(K.EQ.1 .AND. J.EQ.1)  Q(1:2)=Q(3:4)
        IF(MODE.EQ.'2')  THEN
              DO    I=NNNN-1,1,-2
                 Q(NNNN-1-(NNNN-1-I)/2)=
     &          SCL*SCL*(Q(I+1)*Q(I+1)+Q(I)*Q(I))
              ENDDO
        ELSE
              DO    I=NNNN-1,1,-2
                 Q(NNNN-1-(NNNN-1-I)/2)=
     &          SCL*SQRT(Q(I+1)*Q(I+1)+Q(I)*Q(I))
              ENDDO
        ENDIF
              IF(MODE.EQ.'L')  CALL  AL10(Q(NSAM/2),NSAM/2)
              Q(1)=Q(NNNN-1)
              CALL  WRTLIN(LUN2,Q,NSAM,NR)
           ENDDO
        ENDDO

        DO    K=1,NSLICE
           DO    J=1,NROW/2
              JJ=J+NROW/2+MOD(NROW,2)
              NR=J+(K-1)*NROW
              NRJ=JJ+(K-1)*NROW
              CALL  REDLIN(LUN2,Q,NSAM,NR)
              CALL  REDLIN(LUN2,T,NSAM,NRJ)
              CALL  WRTLIN(LUN2,Q,NSAM,NRJ)
              CALL  WRTLIN(LUN2,T,NSAM,NR)
           ENDDO
        ENDDO

        IF(.NOT.IFND)  THEN
           DO    K=1,NSLICE
              NR=NROW/2+1+(K-1)*NROW
              CALL  REDLIN(LUN2,T,NSAM,NR)
              DO    J=NROW/2+1,NROW-1
                 NR=J+(K-1)*NROW
                 NRJ=J+1+(K-1)*NROW
                 CALL  REDLIN(LUN2,Q,NSAM,NRJ)
                 CALL  WRTLIN(LUN2,Q,NSAM,NR)
              ENDDO
              NR=NROW+(K-1)*NROW
              CALL  WRTLIN(LUN2,T,NSAM,NR)
           ENDDO
        ENDIF   

        DO    K=1,NSLICE/2
           KK=K+NSLICE/2+MOD(NSLICE,2)
           DO   J=1,NROW
              NR=J+(K-1)*NROW
              NRK=J+(KK-1)*NROW
              CALL  REDLIN(LUN2,Q,NSAM,NR)
              CALL  REDLIN(LUN2,T,NSAM,NRK)
              CALL  WRTLIN(LUN2,Q,NSAM,NRK)
              CALL  WRTLIN(LUN2,T,NSAM,NR)
           ENDDO
        ENDDO

        IF(.NOT.IFNL)  THEN
           DO    J=1,NROW
              NR=J+(NSLICE/2)*NROW
              CALL  REDLIN(LUN2,T,NSAM,NR)
              DO    K=NSLICE/2+1,NSLICE-1
                 NR=K+(K-1)*NROW
                 NRK=K+K*NROW
                 CALL  REDLIN(LUN2,Q,NSAM,NRK)
                 CALL  WRTLIN(LUN2,Q,NSAM,NR)
              ENDDO
              NR=J+(NSLICE-1)*NROW
              CALL  WRTLIN(LUN2,T,NSAM,NR)
           ENDDO
        ENDIF   

        NSC=NNNN/2-1
        NSL=NSLICE/2
        IF(IFND)  THEN
           IF(IFNL)  THEN
              KB=2
              CALL  REDLIN(LUN2,Q,NSAM,1)
              DO    I=2,NSC
                 II=NNNN-I
                 Q(I)=Q(II)
              ENDDO
              CALL  WRTLIN(LUN2,Q,NSAM,1)
           ELSE
              KB=1
           ENDIF
           JB=2
           DO    K=KB,NSLICE
              KK=NSLICE-K+KB
              NR=1+(K-1)*NROW
              NRK=1+(KK-1)*NROW
              CALL  REDLIN(LUN2,Q,NSAM,NR)
              CALL  REDLIN(LUN2,T,NSAM,NRK)
              Q(1)=T(1)
              DO    I=2,NSC
                 II=NNNN-I
                 Q(I)=T(II)
              ENDDO
              CALL  WRTLIN(LUN2,Q,NSAM,NR)
           ENDDO
        ELSE
           JB=1
        ENDIF
        IF(IFNL)  THEN
           KB=2
           DO    J=JB,NROW
              JJ=NROW-J+JB
              CALL  REDLIN(LUN2,Q,NSAM,J)
              CALL  REDLIN(LUN2,T,NSAM,JJ)
              Q(1)=T(1)
              DO    I=2,NSC
                 II=NNNN-I
                 Q(I)=Q(II)
              ENDDO
              CALL  WRTLIN(LUN2,Q,NSAM,J)     
           ENDDO
        ELSE
           KB=1
        ENDIF
        DO    K=KB,NSLICE
           KK=NSLICE-K+KB
           DO    J=JB,NROW
              JJ=NROW-J+JB
              NR=J+(K-1)*NROW
              NRJ=JJ+(KK-1)*NROW
              CALL  REDLIN(LUN2,Q,NSAM,NR)
              CALL  REDLIN(LUN2,T,NSAM,NRJ)
              Q(1)=T(1)
              DO    I=2,NSC
                 II=NNNN-I
                 Q(I)=Q(II)
              ENDDO
              CALL  WRTLIN(LUN2,Q,NSAM,NR)    
           ENDDO
        ENDDO
C
        END
