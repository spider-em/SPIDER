C ++********************************************************************
C                                                                      *
C  MD2                                                                                       *
C      MD2_NOLUN                                 AUG 2012 ARDEAN LEITH                                                                                *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C  MD2(X,NX,NY,L,K,MODE,LUNO)                                                                 *
C                                                                      
C  PURPOSE:  MEDIAN FILTER ON IMAGE                                                            *
C                                                                      
C  PARAMETERS: L = BOX SIZE
C              K = NUMBER OF PIXELS TO BE TESTED (BOX = L**2)                                                        
C                                                                      
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE  MD2(X, NX,NY,L,K,MODE, LUNO)

        INTEGER           :: NX,NY,L,K
        REAL              :: X(NX,NY)
        CHARACTER(LEN=1)  :: MODE
        INTEGER           :: LUNO

C       AUTOMATIC ARRAYS
        REAL              :: Y(NX)
        REAL              :: A(K)

        LH  = L/2
        K21 = K/2+1

        DO J=1,NY

           IF (MODE == 'C')  THEN
              DO I=1,NX
                 LB = 0

                 DO M=-LH,LH
                    IF (M .NE. 0)  THEN
                       LB    = LB + 1
                       A(LB) = X(I,MOD(J+M+NY-1,NY)+1)
                       LB    = LB + 1
                       A(LB) = X(MOD(I+M+NX-1,NX)+1,J)
                    ELSE
                       LB    = LB + 1
                       A(LB) = X(I,J)
                    ENDIF
                 ENDDO

                 CALL FSORT(A,K)

                 Y(I) = A(K21)
              ENDDO

           ELSE
              DO I=1,NX
                 LB = 0

                 DO MJ=-LH,LH
                    MJM = MOD(J+MJ+NY-1,NY)+1
                    DO MI=-LH,LH
                       LB    = LB + 1
                       A(LB) = X(MOD(I+MI+NX-1,NX)+1,MJM)
                    ENDDO
                 ENDDO

                 CALL FSORT(A,K)

                 Y(I) = A(K21)
              ENDDO
           ENDIF

           CALL  WRTLIN(LUNO,Y,NX,J)
        ENDDO

        END 



C ****************************** MD2_NOLUN *****************************

        SUBROUTINE  MD2_NOLUN(XIN,XOUT, NX,NY, L,K, MODE)

        IMPLICIT NONE
        INTEGER           :: NX,NY,L,K
        REAL              :: XIN (NX,NY)
        REAL              :: XOUT(NX,NY)
        CHARACTER(LEN=1)  :: MODE

C       AUTOMATIC ARRAYS
        REAL              :: Y(NX)
        REAL              :: A(K)

        INTEGER           :: LH,K21,J,I,LB,M,MJ,MJM,MI

        LH  = L / 2
        K21 = K / 2 + 1

        DO J=1,NY

           IF (MODE == 'C')  THEN
C             CROSS SHAPED BOX
              DO I=1,NX
                 LB = 0

                 DO M=-LH,LH
                    IF (M .NE. 0)  THEN
                       LB    = LB + 1
                       A(LB) = XIN(I, MOD(J+M+NY-1,NY)+1)

                       LB    = LB + 1
                       A(LB) = XIN(MOD(I+M+NX-1,NX)+1, J)
                    ELSE
                       LB    = LB + 1
                       A(LB) = XIN(I,J)
                    ENDIF
                 ENDDO

                 CALL FSORT(A,K)

                 Y(I) = A(K21)
              ENDDO
           ELSE
C             SQUARE SHAPED BOX
              DO I=1,NX
                 LB = 0

                 DO MJ=-LH,LH
                    MJM = MOD(J+MJ+NY-1,NY)+1

                    DO MI=-LH,LH
                       LB    = LB + 1
                       A(LB) = XIN(MOD(I+MI+NX-1,NX)+1, MJM)
                    ENDDO
                 ENDDO

                 CALL FSORT(A,K)

                 Y(I) = A(K21)
              ENDDO
           ENDIF

           XOUT(1:NX,J) = Y(1:NX)

        ENDDO

        END 
