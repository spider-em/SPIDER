
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
C     ORTHONORMALIZATION OF THE CHOSEN COLUMNS  OF X(IDIM,*)
C     INPUT         1/ IDIM   RESERVED DIMENSION FOR X(IDIM,*)
C                   3/ JCARD  NBR. OF COLUMNS TO PROCESS X(IDIM, JCARD)         
C                   4/ P(*)   WEIGHT VECTOR , DIMENSION P(IDIM)
C     INPUT - OUTPUT  5/ X(*,*) INPUT  THE MATRIX TO BE PROCESSED
C                             OUTPUT ORTHONORMALIZED MATRIX (P METRIC)
C      IF INB(J) = .TRUE.  THEN J'TH COLUMN HAS TO BE PROCESSED
C      IF INB(J) = .FALSE.  THEN J'TH COLUMN IS ALREADY ORTHONORMALIZED
C
C
C **********************************************************************
C

        SUBROUTINE  GSMODl(IDIM, JCARD, P, X, KRANG, v, inb)

        DIMENSION  X(IDIM, JCARD) ,  P(IDIM)
        logical  inb(jcard), v(jcard)
        DATA  EPS / 1.0 E - 10 /

        v(1) = .true.
        do  j = 1, jcard
          if(.not.inb(j))  v(1) = .false.
        enddo
        if(v(1))  then
          CALL GSMOD(idim, idim,  JBASE,  PJ,  x,  KRANG, v)
          return
        endif
        KRANG   =  JCARD
        do  j = 1, jcard
          v(j) = .not.inb(j)
        enddo
        DO  J  =  1, JCARD
          if(inb(j))  then
            DO JJ  =  1, JCARD
              if(v(jj))  then
                TJJ   =  0.0
                DO   I  =  1, IDIM
                  TJJ   =  TJJ  +  P(I) * X(I, JJ) * X(I, J)
                END DO
                DO   I  =  1, IDIM
                  X(I, J) =  X(I, J)  -  TJJ * X(I, Jj)
                END DO
                Vv = 0.0
                DO  I  =  1, IDIM
                  Vv =  Vv  +  P(I) * X(I, J) * X(I, J)
                END DO
                Vv  =  AMAX1(Vv, EPS)
                C  =  1.0 / SQRT(Vv)
                DO  I  =  1, IDIM
                  X(I, j)  =  C * X(I, j)
                END DO
              endif
            END DO
            v(j) = .true.
          endif
        END DO
c
c	for debugging, discard after that.
c
#ifdef NEVER
        if(jcard.eq.jcard)  return
        print  *,'  GSMODL'
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
#endif
        END
