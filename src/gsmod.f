C
C **********************************************************************
C
C
C **********************************************************************
C *  AUTHOR :                                                              *
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
C     ORTHONORMALIZATION OF THE  JCARD (FIRST) COLUMNS  OF                      
C     X(ICARD,*) BY THE METHOD GRAM - SCHMIDT (MODIFIED) .
C     INPUT         1/ IDIM   RESERVED DIMENSION FOR X(IDIM,*)
C                   2/ ICARD  ACTUAL NUMBER OF ROWS FOR X(ICARD,*)
C                   3/ JCARD  NBR. OF COLUMNS TO PROCESS X(ICARD, JCARD)
C                   4/ P(*)   WEIGHT VECTOR , DIMENSION P(IDIM) 
C     INPUT - OUTPUT  5/ X(*,*) INPUT  THE MATRIX TO BE PROCESSED
C                             OUTPUT ORTHONORMALIZED MATRIX (P METRIC)          
C     OUTPUT        6/ KRANG  RANK OF MATRIX X(ICARD, JCARD)  
C                   7, 8/ T(*), V(*) WORKING ARRAYS T(IDIM), V(IDIM)
C  IF THERE IS COLINEARITY,  THE CORRESPONDING COLUMN IS SET TO ZERO 
C
C
C **********************************************************************
C
        
        SUBROUTINE  GSMOD  (IDIM, ICARD, JCARD, P, X, KRANG, V)

        DIMENSION  X(IDIM, JCARD) ,  P(IDIM) ,  V(IDIM)
        DATA  EPS / 1.0 E - 10 /

        KRANG   =  JCARD 

C     INITIAL NORMS. ORTHONORMALIZATION OF X(* , 1)                           
        DO   J  =  1, JCARD 
          V(J)    =  0.0                          
          DO  I  =  1, ICARD 
             V(J)    =  V(J)  +  P(I) * X(I, J) * X(I, J)
          END DO
          V(J) = AMAX1(V(J), EPS)
        END DO

        C  =  1.0 / SQRT(V(1))                                                 
        DO  I  =  1, ICARD                                                      
          X(I, 1)  =  C * X(I, 1)                     
        END DO
        IF (JCARD  .EQ.  1)  RETURN

C     ORTHOGONALISATION OF X(*, J1). MODIFICATION OF THE NEXT.               
        DO  J  =  1, JCARD  -  1
          J1  =  J  +  1                      
          DO   JJ  =  J1, JCARD                
            TJJ   =  0.0                        
            DO   I  =  1, ICARD                
              TJJ   =  TJJ  +  P(I) * X(I, JJ) * X(I, J)
            END DO
            DO   I  =  1, ICARD             
              X(I, JJ) =  X(I, JJ)  -  TJJ * X(I, J)
            END DO
          END DO

C      COLINEARITY TEST. NORMALIZATION OF X(*, J1).                           
          C  =  0.0 
          DO  I  =  1, ICARD             
            C  =  C  +  P(I) * X(I, J1) * X(I, J1)
          END DO
          IF (C/V(J1) .LE. EPS)  THEN
            KRANG   =  KRANG  -  1
            DO   I  =  1, ICARD
              X(I, J1) = 0.0
            END DO
          ELSE
            C  =  1.0 / SQRT(C)
            DO  I  =  1, ICARD
              X(I, J1) =  C * X(I, J1)
            END DO
          END IF
        END DO

c
c	these following lines  are here for debugging the program
c	and are discarded when done.
c
        if(jcard.eq.jcard)  return
        print  *,'  GSMOD'
        do  j=1,jcard
        qt=0.0
        do  i=1,idim
        qt=qt+p(i)*x(i,j)*x(i,j)
        enddo
        print *,' Norm of column #',j,' =',qt
        enddo
        do  j=1,jcard-1
        do  k=j+1,jcard
        qt=0.0
        do  i=1,idim
        qt=qt+p(i)*x(i,j)*x(i,k)
        enddo
        print *,' Cosine between columns #',j,k,' =',qt
        enddo
        enddo
                            
        RETURN
        END
