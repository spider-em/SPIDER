C ++********************************************************************
C                                                                      
C                                                                      
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
C  PURPOSE:                                                            
C
C  Calculates circular croscorrelation, non-power-of-two dimensions
C  Input - X,Y Fourier transforms
C  Output -  O=F(X*conjg(Y))
C
C  PARAMETERS:                                                         
C
C IMAGE_PROCESSING_ROUTINE
C
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012 
C***********************************************************************

        SUBROUTINE  CRSM_2(X,Y,O,NSAM,NROW,WRK)
C  WRK(MAX0(NSAM,NROW))

        DIMENSION  X(NSAM,NROW),Y(NSAM,NROW),O(NSAM,NROW),WRK(*)
        LOGICAL  IFNS,IFNR

        IFNS=MOD(NSAM,2).EQ.0
        IFNR=MOD(NROW,2).EQ.0

C       INS=1
C       CALL  Fmr_2(X,NSAM,NROW,WRK,INS)
C       IF(INS.EQ.0)  RETURN
C       CALL  Fmr_2(Y,NSAM,NROW,WRK,INS)

        I=NSAM/2-1
c$omp parallel do private(j)
        DO   J=1,NROW
           CALL  MLC(X(3,J),Y(3,J),O(3,J),I)
        ENDDO

        IF(IFNS)  THEN
c$omp parallel do private(i,j,kj)
           DO    I=1,2
              DO    J=2,NROW/2
                 KJ=2*J-1
                 O(I,KJ+1)=-X(I,KJ)*Y(I,KJ+1)+X(I,KJ+1)*Y(I,KJ)
                 O(I,KJ)=X(I,KJ)*Y(I,KJ)+X(I,KJ+1)*Y(I,KJ+1)
              ENDDO
              IF(.NOT.IFNR)  THEN
                 O(I,2)=-X(I,NROW)*Y(I,2)+X(I,2)*Y(I,NROW)
                 O(I,NROW)=X(I,NROW)*Y(I,NROW)+X(I,2)*Y(I,2)
              ENDIF
           ENDDO
           O(2,1)=X(2,1)*Y(2,1)
           IF(IFNR)  THEN
              O(1,2)=X(1,2)*Y(1,2)
              O(2,2)=X(2,2)*Y(2,2)
           ENDIF
        ELSE    

c$omp parallel do private(j,kj)
           DO    J=2,NROW/2
              KJ=2*J-1
              O(1,KJ+1)=-X(1,KJ)*Y(1,KJ+1)+X(1,KJ+1)*Y(1,KJ)
              O(1,KJ)=X(1,KJ)*Y(1,KJ)+X(1,KJ+1)*Y(1,KJ+1)
           ENDDO
           IF(IFNR)  THEN
              O(1,2)=X(1,2)*Y(1,2)
           ELSE
              O(1,2)=-X(1,NROW)*Y(1,2)+X(1,2)*Y(1,NROW)
              O(1,NROW)=X(1,NROW)*Y(1,NROW)+X(1,2)*Y(1,2)
           ENDIF
c$omp parallel do private(j)
           DO    J=1,NROW
              O(2,J)=-X(NSAM,J)*Y(2,J)+X(2,J)*Y(NSAM,J)
              O(NSAM,J)=X(NSAM,J)*Y(NSAM,J)+X(2,J)*Y(2,J)
           ENDDO
        ENDIF
        O(1,1)=X(1,1)*Y(1,1)

        INS=-1
        CALL  Fmr_2(O,NSAM,NROW,WRK,INS)

        NR2=NROW/2
        NS2=NSAM/2
        IF(IFNS.AND.IFNR)  THEN
           DO    J=1,NR2
              JJ=J+NR2
              DO    I=1,NS2
                 II=I+NS2
                 Q=O(I,J)
                 O(I,J)=O(II,JJ)
                 O(II,JJ)=Q
                 Q=O(I,JJ)
                 O(I,JJ)=O(II,J)
                 O(II,J)=Q
              ENDDO
           ENDDO
        ELSEIF(.NOT.IFNS .AND. .NOT.IFNR)  THEN
           K=0
           DO    I=1,NS2+1
              K=K+1
              WRK(K)=O(I,NR2+1)
           ENDDO
           DO    J=1,NR2
              K=K+1
              WRK(K)=O(NS2+1,J)
           ENDDO
           DO J=1,NR2+1
              DO    I=1,NS2
                 O(I,J+NR2)=O(I+NS2+1,J)
              ENDDO
              IF(J.EQ.NR2+1)  GOTO  23
              DO    I=1,NS2+1
                 O(I+NS2,J)=O(I,J+NR2+1)
              ENDDO
23            CONTINUE
           ENDDO
           DO    J=1,NR2
              DO    I=1,NS2
                 O(I+NS2,J+NR2)=O(I,J)
                 O(I,J)=O(I+NS2+1,J+NR2+1)
              ENDDO
           ENDDO
           K=0
           DO    I=1,NS2+1
              K=K+1
              O(I+NS2,NROW)=WRK(K)
           ENDDO
           DO    J=1,NR2
              K=K+1
              O(NSAM,J+NR2)=WRK(K)
           ENDDO

        ELSEIF(.NOT.IFNS .AND.IFNR)  THEN
           DO    J=1,NR2
              WRK(J)=O(NS2+1,J)
           ENDDO

           DO I=1,NS2+1
              DO    J=1,NR2
                 O(I+NS2,J)=O(I,J+NR2)
              ENDDO
              IF(I.EQ.NS2+1)  GOTO  52
              DO    J=1,NR2
                 O(I,J+NR2)=O(I+NS2+1,J)
              ENDDO
52            CONTINUE
           ENDDO
           DO    J=1,NR2
              DO    I=1,NS2
                 O(I+NS2,J+NR2)=O(I,J)
                 O(I,J)=O(I+NS2+1,J+NR2)
              ENDDO
           ENDDO
           DO    J=1,NR2
              O(NSAM,J+NR2)=WRK(J)
           ENDDO

        ELSEIF(IFNS .AND. .NOT.IFNR)  THEN
           DO    I=1,NS2
              WRK(I)=O(I,NR2+1)
           ENDDO

           DO J=1,NR2+1
              DO    I=1,NS2
                 O(I,J+NR2)=O(I+NS2,J)
              ENDDO
              IF(J.EQ.NR2+1)  GOTO  62
              DO    I=1,NS2
                 O(I+NS2,J)=O(I,J+NR2+1)
              ENDDO
62            CONTINUE
           ENDDO
           DO    I=1,NS2
              DO    J=1,NR2
                 O(I+NS2,J+NR2)=O(I,J)
                 O(I,J)=O(I+NS2,J+NR2+1)
              ENDDO
           ENDDO
           DO    I=1,NS2
              O(I+NS2,NROW)=WRK(I)
           ENDDO
        ENDIF
        END
