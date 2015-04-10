
C++*********************************************************************
C
C PW2SR.F
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2011  Health Research Inc.,                         *
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
C PW2SR(Q,NSAM,NROW,MODE)
C
C PURPOSE: POWER SPECTRUM OF IMAGE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE   PW2SR(Q,NSAM,NROW,MODE)

        DIMENSION    :: Q(NSAM+2-MOD(NSAM,2),NROW)
        CHARACTER*1  :: MODE
        LOGICAL      :: IFND

        NNNN = NSAM+2-MOD(NSAM,2)
        IFND = MOD(NROW,2).EQ.0
        NSC  = NSAM/2+1

        SCL = 2.0/FLOAT(NSAM)/FLOAT(NROW)

        DO J=1,NROW
           IF (MODE .EQ. '2')  THEN  ! THIS IS NEVER USED JUN 2011!
              DO I=NNNN-1,1,-2
                 Q(NNNN-1-(NNNN-1-I)/2,J)=
     &            SCL*SCL*(Q(I+1,J)*Q(I+1,J)+Q(I,J)*Q(I,J))
              ENDDO
           ELSE
              DO I=NNNN-1,1,-2
                 Q(NNNN-1-(NNNN-1-I)/2,J)=
     &            SCL*SQRT(Q(I+1,J)*Q(I+1,J)+Q(I,J)*Q(I,J))
              ENDDO
           ENDIF

           IF (MODE .EQ. 'L') CALL AL10(Q(NSAM/2,J),NSAM/2)
        ENDDO

        DO J=1,NROW/2
           JJ = J+NROW/2+MOD(NROW,2)
           DO I=NSC,NNNN-1
              TEMP    = Q(I,J)
              Q(I,J)  = Q(I,JJ)
              Q(I,JJ) = TEMP
           ENDDO
        ENDDO

        IF (.NOT.IFND)  THEN   !ODD ROW LENGTH
           DO I=NSC,NNNN-1
              TEMP = Q(I,NROW/2+1)
              DO J=NROW/2+1,NROW-1
                 Q(I,J) = Q(I,J+1)
              ENDDO
              Q(I,NROW) = TEMP
          ENDDO
        ENDIF   

        NSC = NNNN/2-1
        IF (IFND) THEN    ! EVEN ROW LENGTH
           JB = 2
           DO I=1,NSC
              II     = NNNN-I
              Q(I,1) = Q(II,1)
           ENDDO
        ELSE
           JB = 1
        ENDIF

        DO J=JB,NROW
           JJ = NROW-J+JB
           DO I=1,NSC
              II     = NNNN-I
              Q(I,J) = Q(II,JJ)
           ENDDO
        ENDDO

         Q(NSAM/2+1,NROW/2+1) = Q(NSAM/2,NROW/2)

         END
