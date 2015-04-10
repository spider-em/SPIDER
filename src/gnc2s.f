C ++********************************************************************
C                                                                      *
C  GNC2S                                                               *
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
C  GNC2S                                                              *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C
C IMAGE_PROCESSING_ROUTINE
C
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C **********************************************************************

        SUBROUTINE  GNC2S(X1,U1,N,M,QL,H0,EPS,LUNO)

        DIMENSION  X1(N,M),U1(N,M)
        COMMON  /GP_D/  C,QL2,R,Q
        REAL, ALLOCATABLE, DIMENSION(:,:) :: U2
        DO    I=1,N
        DO    J=1,M
        U1(I,J)=X1(I,J)
        ENDDO
        ENDDO
C
        CS=0.25
        ALPHA=H0*H0*QL/2.0
        W=2.0*(1.0-1.0/SQRT(2.0)/QL)
        QL2=QL*QL
C
        WRITE(LUNO,71)  QL,H0,EPS,ALPHA
 71     FORMAT(//' GNC ALGORITHM 2-D'
     *  ,/,' PARAMETERS: LAMBDA    H0       EPS       ALPHA'
     *  ,/,12X,F5.0,3X,3G10.3,/)
C
        P=1.0

 
        ALLOCATE (U2(N,M), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'GNC2S, U2',IER)
           RETURN
        ENDIF


 100    CONTINUE
C
        C=CS/P
        R=SQRT(ALPHA*(2.0/C+1.0/QL2))
        Q=ALPHA/QL2/R
	
	

C
C  Start from (1,1)
        U2(1,1)=U1(1,1)-W*(2.0*(U1(1,1)-X1(1,1))+
     1          GP(U1(1,1)-U1(2,1))+
     2          GP(U1(1,1)-U1(1,2)))/
     3          (2.0+4.0*QL2)
C
        WT=W/(2.0+6.0*QL2)
        DO    J=2,M-1
        U2(1,J)=U1(1,J)-WT*(2.0*(U1(1,J)-X1(1,J))+
     1          GP(U1(1,J)-U2(1,J-1))+
     2          GP(U1(1,J)-U1(2,J))+
     3          GP(U1(1,J)-U1(1,J+1)))
        ENDDO
C
        U2(1,M)=U1(1,M)-W*(2.0*(U1(1,M)-X1(1,M))+
     1          GP(U1(1,M)-U2(1,M-1))+
     2          GP(U1(1,M)-U1(2,M)))/
     3          (2.0+4.0*QL2)
C
        WT=W/(2.0+6.0*QL2)
        DO    I=2,N-1
        U2(I,1)=U1(I,1)-WT*(2.0*(U1(I,1)-X1(I,1))+
     1          GP(U1(I,1)-U2(I-1,1))+
     2          GP(U1(I,1)-U1(I,2))+
     3          GP(U1(I,1)-U1(I+1,1)))
        ENDDO
C
        U2(N,1)=U1(N,1)-W*(2.0*(U1(N,1)-X1(N,1))+
     1          GP(U1(N,1)-U2(N-1,1))+
     2          GP(U1(N,1)-U1(N,2)))/
     3          (2.0+4.0*QL2)
C
        WT=W/(2.0+8.0*QL2)
        DO    I=2,N-1
        DO    J=2,M-1
        U2(I,J)=U1(I,J)-WT*(2.0*(U1(I,J)-X1(I,J))+
     1          GP(U1(I,J)-U2(I-1,J))+
     2          GP(U1(I,J)-U2(I,J-1))+
     3          GP(U1(I,J)-U1(I+1,J))+
     4          GP(U1(I,J)-U1(I,J+1)))
        ENDDO
        ENDDO
C
        WT=W/(2.0+6.0*QL2)
        DO    I=2,N-1
        U2(I,M)=U1(I,M)-WT*(2.0*(U1(I,M)-X1(I,M))+
     1          GP(U1(I,M)-U2(I-1,M))+
     2          GP(U1(I,M)-U2(I,M-1))+
     3          GP(U1(I,M)-U1(I+1,M)))
        ENDDO

C
        WT=W/(2.0+6.0*QL2)
        DO    J=2,M-1
        U2(N,J)=U1(N,J)-WT*(2.0*(U1(N,J)-X1(N,J))+
     1          GP(U1(N,J)-U2(N,J-1))+
     2          GP(U1(N,J)-U2(N-1,J))+
     3          GP(U1(N,J)-U1(N,J+1)))
        ENDDO
C
        U2(N,M)=U1(N,M)-W*(2.0*(U1(N,M)-X1(N,M))+
     1          GP(U1(N,M)-U2(N-1,M))+
     2          GP(U1(N,M)-U2(N,M-1)))/
     3          (2.0+4.0*QL2)
C
        ER=ERC(U1,U2,N*M)
        IF(ER.LT.EPS)  GOTO  1000
C
C  Start from (1,M)
        U2(1,M)=U1(1,M)-W*(2.0*(U1(1,M)-X1(1,M))+
     1          GP(U1(1,M)-U1(1,M-1))+
     2          GP(U1(1,M)-U1(2,M)))/
     3          (2.0+4.0*QL2)
C
        WT=W/(2.0+6.0*QL2)
        DO    J=M-1,2,-1
        U2(1,J)=U1(1,J)-WT*(2.0*(U1(1,J)-X1(1,J))+
     1          GP(U1(1,J)-U1(1,J-1))+
     2          GP(U1(1,J)-U1(2,J))+
     3          GP(U1(1,J)-U2(1,J+1)))
        ENDDO
C
        U2(1,1)=U1(1,1)-W*(2.0*(U1(1,1)-X1(1,1))+
     1          GP(U1(1,1)-U1(2,1))+
     2          GP(U1(1,1)-U2(1,2)))/
     3          (2.0+4.0*QL2)
C
        WT=W/(2.0+6.0*QL2)
        DO    I=2,N-1
        U2(I,M)=U1(I,M)-WT*(2.0*(U1(I,M)-X1(I,M))+
     1          GP(U1(I,M)-U2(I-1,M))+
     2          GP(U1(I,M)-U1(I,M-1))+
     3          GP(U1(I,M)-U1(I+1,M)))
        ENDDO
C
        U2(N,M)=U1(N,M)-W*(2.0*(U1(N,M)-X1(N,M))+
     1          GP(U1(N,M)-U1(N-1,M))+
     2          GP(U1(N,M)-U2(N,M-1)))/
     3          (2.0+4.0*QL2)
C
        WT=W/(2.0+8.0*QL2)
        DO    I=2,N-1
        DO    J=M-1,2,-1
        U2(I,J)=U1(I,J)-WT*(2.0*(U1(I,J)-X1(I,J))+
     1          GP(U1(I,J)-U2(I-1,J))+
     2          GP(U1(I,J)-U1(I,J-1))+
     3          GP(U1(I,J)-U1(I+1,J))+
     4          GP(U1(I,J)-U2(I,J+1)))
        ENDDO
        ENDDO
C
        WT=W/(2.0+6.0*QL2)
        DO    J=M-1,2,-1
        U2(N,J)=U1(N,J)-WT*(2.0*(U1(N,J)-X1(N,J))+
     1          GP(U1(N,J)-U1(N,J-1))+
     2          GP(U1(N,J)-U2(N-1,J))+
     3          GP(U1(N,J)-U2(N,J+1)))
        ENDDO
C
        WT=W/(2.0+6.0*QL2)
        DO    I=2,N-1
        U2(I,1)=U1(I,1)-WT*(2.0*(U1(I,1)-X1(I,1))+
     1          GP(U1(I,1)-U2(I-1,1))+
     2          GP(U1(I,1)-U2(I,2))+
     3          GP(U1(I,1)-U1(I+1,1)))
        ENDDO
C
        U2(N,1)=U1(N,1)-W*(2.0*(U1(N,1)-X1(N,1))+
     1          GP(U1(N,1)-U2(N-1,1))+
     2          GP(U1(N,1)-U2(N,2)))/
     3          (2.0+4.0*QL2)
C
        ER=ERC(U1,U2,N*M)
        IF(ER.LT.EPS)  GOTO  1000
C
C  Start from (N,M)
        U2(N,M)=U1(N,M)-W*(2.0*(U1(N,M)-X1(N,M))+
     1          GP(U1(N,M)-U1(N-1,M))+
     2          GP(U1(N,M)-U1(N,M-1)))/
     3          (2.0+4.0*QL2)
C
        WT=W/(2.0+6.0*QL2)
        DO    J=M-1,2,-1
        U2(N,J)=U1(N,J)-WT*(2.0*(U1(N,J)-X1(N,J))+
     1          GP(U1(N,J)-U1(N,J-1))+
     2          GP(U1(N,J)-U1(N-1,J))+
     3          GP(U1(N,J)-U2(N,J+1)))
        ENDDO
C
        WT=W/(2.0+6.0*QL2)
        DO    I=N-1,2,-1
        U2(I,M)=U1(I,M)-WT*(2.0*(U1(I,M)-X1(I,M))+
     1          GP(U1(I,M)-U1(I-1,M))+
     2          GP(U1(I,M)-U1(I,M-1))+
     3          GP(U1(I,M)-U2(I+1,M)))
        ENDDO
C
        U2(1,M)=U1(1,M)-W*(2.0*(U1(1,M)-X1(1,M))+
     1          GP(U1(1,M)-U1(1,M-1))+
     2          GP(U1(1,M)-U2(2,M)))/
     3          (2.0+4.0*QL2)
C
        U2(N,1)=U1(N,1)-W*(2.0*(U1(N,1)-X1(N,1))+
     1          GP(U1(N,1)-U1(N-1,1))+
     2          GP(U1(N,1)-U2(N,2)))/
     3          (2.0+4.0*QL2)
C
        WT=W/(2.0+8.0*QL2)
        DO    I=N-1,2,-1
        DO    J=M-1,2,-1
        U2(I,J)=U1(I,J)-WT*(2.0*(U1(I,J)-X1(I,J))+
     1          GP(U1(I,J)-U1(I-1,J))+
     2          GP(U1(I,J)-U1(I,J-1))+
     3          GP(U1(I,J)-U2(I+1,J))+
     4          GP(U1(I,J)-U2(I,J+1)))
        ENDDO
        ENDDO
C
        WT=W/(2.0+6.0*QL2)
        DO    I=N-1,2,-1
        U2(I,1)=U1(I,1)-WT*(2.0*(U1(I,1)-X1(I,1))+
     1          GP(U1(I,1)-U1(I-1,1))+
     2          GP(U1(I,1)-U2(I,2))+
     3          GP(U1(I,1)-U2(I+1,1)))
        ENDDO
C
        WT=W/(2.0+6.0*QL2)
        DO    J=M-1,2,-1
        U2(1,J)=U1(1,J)-WT*(2.0*(U1(1,J)-X1(1,J))+
     1          GP(U1(1,J)-U1(1,J-1))+
     2          GP(U1(1,J)-U2(2,J))+
     3          GP(U1(1,J)-U2(1,J+1)))
        ENDDO
C
        U2(1,1)=U1(1,1)-W*(2.0*(U1(1,1)-X1(1,1))+
     1          GP(U1(1,1)-U2(2,1))+
     2          GP(U1(1,1)-U2(1,2)))/
     3          (2.0+4.0*QL2)
C
        ER=ERC(U1,U2,N*M)
        IF(ER.LT.EPS)  GOTO  1000
C
C  Start from (N,1)
        U2(N,1)=U1(N,1)-W*(2.0*(U1(N,1)-X1(N,1))+
     1          GP(U1(N,1)-U1(N-1,1))+
     2          GP(U1(N,1)-U1(N,2)))/
     3          (2.0+4.0*QL2)
C
        WT=W/(2.0+6.0*QL2)
        DO    I=N-1,2,-1
        U2(I,1)=U1(I,1)-WT*(2.0*(U1(I,1)-X1(I,1))+
     1          GP(U1(I,1)-U1(I-1,1))+
     2          GP(U1(I,1)-U1(I,2))+
     3          GP(U1(I,1)-U2(I+1,1)))
        ENDDO
C
        WT=W/(2.0+6.0*QL2)
        DO    J=2,M-1
        U2(N,J)=U1(N,J)-WT*(2.0*(U1(N,J)-X1(N,J))+
     1          GP(U1(N,J)-U2(N,J-1))+
     2          GP(U1(N,J)-U1(N-1,J))+
     3          GP(U1(N,J)-U1(N,J+1)))
        ENDDO
C
        U2(N,M)=U1(N,M)-W*(2.0*(U1(N,M)-X1(N,M))+
     1          GP(U1(N,M)-U1(N-1,M))+
     2          GP(U1(N,M)-U2(N,M-1)))/
     3          (2.0+4.0*QL2)
C
        U2(1,1)=U1(1,1)-W*(2.0*(U1(1,1)-X1(1,1))+
     1          GP(U1(1,1)-U2(2,1))+
     2          GP(U1(1,1)-U1(1,2)))/
     3          (2.0+4.0*QL2)
C
        WT=W/(2.0+8.0*QL2)
        DO    I=N-1,2,-1
        DO    J=2,M-1
        U2(I,J)=U1(I,J)-WT*(2.0*(U1(I,J)-X1(I,J))+
     1          GP(U1(I,J)-U1(I-1,J))+
     2          GP(U1(I,J)-U2(I,J-1))+
     3          GP(U1(I,J)-U2(I+1,J))+
     4          GP(U1(I,J)-U1(I,J+1)))
        ENDDO
        ENDDO
C
        WT=W/(2.0+6.0*QL2)
        DO    J=2,M-1
        U2(1,J)=U1(1,J)-WT*(2.0*(U1(1,J)-X1(1,J))+
     1          GP(U1(1,J)-U2(1,J-1))+
     2          GP(U1(1,J)-U2(2,J))+
     3          GP(U1(1,J)-U1(1,J+1)))
        ENDDO
C
        WT=W/(2.0+6.0*QL2)
        DO    I=N-1,2,-1
        U2(I,M)=U1(I,M)-WT*(2.0*(U1(I,M)-X1(I,M))+
     1          GP(U1(I,M)-U1(I-1,M))+
     2          GP(U1(I,M)-U2(I,M-1))+
     3          GP(U1(I,M)-U2(I+1,M)))
        ENDDO
C
        U2(1,M)=U1(1,M)-W*(2.0*(U1(1,M)-X1(1,M))+
     1          GP(U1(1,M)-U2(1,M-1))+
     2          GP(U1(1,M)-U2(2,M)))/
     3          (2.0+4.0*QL2)
C
        ER=ERC(U1,U2,N*M)
        IF(ER.LT.EPS)  GOTO  1000
        GOTO  100
 1000   WRITE(LUNO,24)  P,ER
C       PRINT *,P,ER
        P=P/2.0
        IF(P.GE.1.0/2.0/QL)  GOTO  100

        DEALLOCATE(U2)
        RETURN
 24     FORMAT(10X,'STEP=',F12.6,'  ERROR=',G10.3)
        
        END
