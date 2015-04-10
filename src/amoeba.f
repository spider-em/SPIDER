C++*********************************************************************
C
C AMOEBA.F
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
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE AMOEBA(P,Y,NDIM,FTOL,FUNK,ITER,ITMAX,PR,PRR,PBAR)

         INCLUDE 'CMBLOCK.INC'

         PARAMETER (ALPHA=1.0,BETA=0.5,GAMMA=2.0)
         DIMENSION P(NDIM+1,NDIM),Y(NDIM+1),PR(NDIM)
         DIMENSION PRR(NDIM),PBAR(NDIM)

         EXTERNAL  FUNK

         MPTS = NDIM+1
         ITER = 0

1        ILO = 1
         IF (Y(1) .GT. Y(2)) THEN
            IHI  = 1
            INHI = 2
         ELSE
            IHI  = 2
            INHI = 1
         ENDIF

         DO I=1,MPTS
            IF (Y(I) .LT. Y(ILO)) ILO=I
            IF (Y(I) .GT. Y(IHI)) THEN
               INHI=IHI
               IHI=I
            ELSE IF (Y(I) .GT. Y(INHI)) THEN
               IF (I .NE. IHI) INHI=I
            ENDIF
         ENDDO

         RTOL = 2. * ABS(Y(IHI)-Y(ILO))/(ABS(Y(IHI))+ABS(Y(ILO)))
         IF (RTOL .LT. FTOL) RETURN

         IF (ITER == ITMAX) THEN
            WRITE(NOUT,*) ' Amoeba exceeding maximum iterations.'
            RETURN
         ENDIF

         ITER = ITER + 1
         DO J=1,NDIM
            PBAR(J) = 0.
         ENDDO

         DO I=1,MPTS
            IF (I .NE. IHI) THEN
               DO  J=1,NDIM
                  PBAR(J) = PBAR(J) + P(I,J)
               ENDDO
            ENDIF
         ENDDO

         DO J=1,NDIM
            PBAR(J) = PBAR(J)/NDIM
            PR(J)   = (1.+ALPHA)*PBAR(J)-ALPHA*P(IHI,J)
         ENDDO

         YPR = FUNK(PR)
         IF (YPR .LE. Y(ILO)) THEN
            DO  J=1,NDIM
               PRR(J) = GAMMA * PR(J) + (1.-GAMMA) * PBAR(J)
            ENDDO

            YPRR = FUNK(PRR)
            IF (YPRR .LT. Y(ILO)) THEN
               DO J=1,NDIM
                  P(IHI,J) = PRR(J)
               ENDDO

               Y(IHI) = YPRR
            ELSE
               DO  J=1,NDIM
                  P(IHI,J) = PR(J)
               ENDDO

               Y(IHI) = YPR
            ENDIF
         ELSE IF (YPR .GE. Y(INHI)) THEN
            IF (YPR .LT. Y(IHI)) THEN
               DO  J=1,NDIM
                  P(IHI,J) = PR(J)
               ENDDO

               Y(IHI)=YPR
            ENDIF
            DO J=1,NDIM
               PRR(J) = BETA * P(IHI,J) + (1.-BETA) * PBAR(J)
            ENDDO

            YPRR = FUNK(PRR)
            IF (YPRR .LT. Y(IHI)) THEN
               DO J=1,NDIM
                  P(IHI,J) = PRR(J)
               ENDDO

               Y(IHI) = YPRR
            ELSE
               DO  I=1,MPTS
                  IF (I .NE. ILO) THEN
                     DO  J=1,NDIM
                        PR(J)  = 0.5 * (P(I,J) + P(ILO,J))
                        P(I,J) = PR(J)
                     ENDDO

                     Y(I) = FUNK(PR)
                  ENDIF
               ENDDO
            ENDIF
         ELSE
            DO  J=1,NDIM
               P(IHI,J) = PR(J)
            ENDDO

            Y(IHI) = YPR
         ENDIF

         GO TO 1

         END

