C ++********************************************************************
C                                                                      *
C  MD3                                                                 *
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
C  MD3(X,NSAM,NROW,NSLICE,L,K,MODE,LUNO)                                                                  *
C                                                                      *
C  PURPOSE:  MEDIAN FILTER ON VOLUME                                                          *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C  IMAGE_PROCESSING_ROUTINE                                             *                                                                     *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE  MD3(X,NSAM,NROW,NSLICE,L,K,MODE,LUNO)

        REAL,    DIMENSION(NSAM,NROW,NSLICE) :: X
C       AUTOMATIC ARRAYS
        REAL,    DIMENSION(NSAM)             :: Y
        REAL,    DIMENSION(K)                :: A
        CHARACTER(LEN=1)                     :: MODE

        LH  = L / 2
        K21 = K / 2 + 1

        DO N=1,NSLICE                    
           DO J=1,NROW                      
              IF (MODE .EQ. 'C')  THEN
                 DO I=1,NSAM
                    LB = 0
                    DO M=-LH,LH
                      LB = LB + 1
                      IF (M .NE. 0)  THEN
                          A(LB) = X(I,MOD(J+M+NROW-1,NROW)+1,N)

                          LB    = LB + 1
                          A(LB) = X(MOD(I+M+NSAM-1,NSAM)+1,J,N)

                          LB    = LB + 1
                          A(LB) = X(I,J,MOD(N+M+NSLICE-1,NSLICE)+1)
                       ELSE
                          A(LB) = X(I,J,N)
                       ENDIF
                    ENDDO

                    CALL FSORT(A,K)
                    Y(I) = A(K21)
                 ENDDO
              ELSE
                 DO I=1,NSAM
                    LB = 0
                    DO MN=-LH,LH
                       MNM = MOD(N+MN+NSLICE-1,NSLICE)+1
                       DO MJ=-LH,LH
                          MJM = MOD(J+MJ+NROW-1,NROW)+1
                          DO MI=-LH,LH
                             LB    = LB+1
                             A(LB) = X(MOD(I+MI+NSAM-1,NSAM)+1,MJM,MNM)
                          ENDDO
                       ENDDO
                    ENDDO
                    CALL FSORT(A,K)
                    Y(I) = A(K21)
                 ENDDO
              ENDIF
              CALL WRTLIN(LUNO,Y,NSAM,J+(N-1)*NROW)
           ENDDO
        ENDDO
        END  


