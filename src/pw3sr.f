C++*********************************************************************
C
C PW3SR.F                       HALF BUG FIXED FEB 02 ArDean Leith
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
C PW3SR(Q,NSAM,NROW,NSLICE,MODE)
C
C PURPOSE: POWER SPECTRUM OF VOLUME
C
C IMAGE_PROCESSING_ROUTINE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE PW3SR(Q,NSAM,NROW,NSLICE,MODE)

        DIMENSION    Q(NSAM+2-MOD(NSAM,2),NROW,NSLICE)
        CHARACTER*1  MODE
        LOGICAL*1    IFND,IFNL

        NNNN  = NSAM+2-MOD(NSAM,2)
        IFND  = MOD(NROW,2).EQ.0
        IFNL  = MOD(NSLICE,2).EQ.0
        NSC   = NSAM/2+1

        SCL = 2.0/FLOAT(NSAM)/FLOAT(NROW)/FLOAT(NSLICE)

        DO K=1,NSLICE
	   IF (MODE .EQ. '2')  THEN
              DO J=1,NROW
                 DO I=NNNN-1,1,-2
                    Q(NNNN-1-(NNNN-1-I)/2,J,K) =
     &               SCL*SCL*(Q(I+1,J,K)*Q(I+1,J,K)+Q(I,J,K)*Q(I,J,K))
                 ENDDO
              ENDDO

	   ELSE
              DO J=1,NROW
                 DO I=NNNN-1,1,-2
                    Q(NNNN-1-(NNNN-1-I)/2,J,K) =
     &               SCL*SQRT(Q(I+1,J,K)*Q(I+1,J,K)+Q(I,J,K)*Q(I,J,K))
                 ENDDO
              ENDDO
	   ENDIF
        ENDDO

        DO K=1,NSLICE
           DO J=1,NROW/2
              JJ = J+NROW/2+MOD(NROW,2)
              DO I=NSC,NNNN-1
                 TEMP      = Q(I,J,K)
                 Q(I,J,K)  = Q(I,JJ,K)
                 Q(I,JJ,K) = TEMP
              ENDDO
           ENDDO
        ENDDO

        IF (.NOT. IFND)  THEN
           DO K=1,NSLICE
              DO I=NSC,NNNN-1
                 TEMP = Q(I,NROW/2+1,K)
                 DO J=NROW/2+1,NROW-1
                    Q(I,J,K) = Q(I,J+1,K)
                 ENDDO
                 Q(I,NROW,K) = TEMP
              ENDDO
           ENDDO
        ENDIF   

        DO K=1,NSLICE/2
           KK=K+NSLICE/2+MOD(NSLICE,2)
           DO J=1,NROW
              DO I=NSC,NNNN-1
                 TEMP      = Q(I,J,K)
                 Q(I,J,K)  = Q(I,J,KK)
                 Q(I,J,KK) = TEMP
              ENDDO
           ENDDO
        ENDDO

        IF (.NOT. IFNL)  THEN
           DO J=1,NROW
              DO I=NSC,NNNN-1
                 TEMP=Q(I,J,NSLICE/2+1)
                 DO K=NSLICE/2+1,NSLICE-1
                    Q(I,J,K)   = Q(I,J,K+1)
                 ENDDO
                 Q(I,J,NSLICE) = TEMP
              ENDDO
           ENDDO
        ENDIF   

        NSC = NNNN/2-1
        NSL = NSLICE/2
        IF (IFND)  THEN
           IF (IFNL)  THEN
              KB=2
              DO I=1,NSC
                 II=NNNN-I
                 Q(I,1,1) = Q(II,1,1)
              ENDDO
           ELSE
              KB=1
           ENDIF

           JB = 2
           DO K=KB,NSLICE
              KK=NSLICE-K+KB
              DO I=1,NSC
                 II       = NNNN-I
                 Q(I,1,K) = Q(II,1,KK)
              ENDDO
           ENDDO
        ELSE
           JB=1
        ENDIF

        IF (IFNL)  THEN
           KB=2
           DO J=JB,NROW
              JJ=NROW-J+JB
              DO I=1,NSC
                 II       = NNNN-I
                 Q(I,J,1) = Q(II,JJ,1)     
              ENDDO
           ENDDO
        ELSE
           KB=1
        ENDIF
        DO K=KB,NSLICE
           KK = NSLICE-K+KB
           DO J=JB,NROW
              JJ=NROW-J+JB
              DO I=1,NSC
                 II       = NNNN-I
                 Q(I,J,K) = Q(II,JJ,KK)
              ENDDO
           ENDDO
        ENDDO

        IF (MODE .EQ. 'L') CALL AL10(Q,NNNN*NROW*NSLICE)

	Q(NSAM/2+1,NROW/2+1,NSLICE/2+1) = Q(NSAM/2,NROW/2,NSLICE/2)

        END
