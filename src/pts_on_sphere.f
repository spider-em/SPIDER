C ++********************************************************************
C                                                                      *
C   PTS_ON_SPHERE                                                      *
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
C
C   SUBROUTINE PTS_ON_SPHERE(NT,ITER,THETA,PHI,HALF,NGOT)
C
C   PURPOSE: CREATE N POINTS ON A SPHERE APROXIMATELY EQUI-DISTANT 
C            FROM EACH OTHER.  N POINTS ARE RANDOMLY PLACED ON SPHERE 
C            AND THEN MOVED AROUND UNTIL THE MINIMAL DISTANCE 
C            BETWEEN THE CLOSEST TWO POINTS IS MINIMIZED.
C            ADAPTED FROM ALGORITHM PROPOSED BY: PAUL BOURKE
C
C   VARIABLES: 
C               NT      NUMBER OF POINTS                       (SENT)
C               ITER    NUMBER OF ITERATIONS (SHOULD BE > 100) (SENT)
C               THETA   EULER ANGLE FOR VECTOR                 (RET.)
C               PHI     EULER ANGLE FOR VECTOR                 (RET.)
C               HALF    REPORT +Z ONLY                         (SENT)
C               NGOT    NUMBER OF POINTS RETURNED              (RET.)
C                       (MAY BE LESS THAN NT IF HALF)
C
C **********************************************************************

      SUBROUTINE PTS_ON_SPHERE(N,ITER,X,Y,Z,
     &                        HALF,SPIRAL,MOVE,NGOT,IRTFLG)

      INCLUDE 'CMBLOCK.INC'

      DOUBLE PRECISION :: D, DMIND, DMAXD, VALL, BOTT
      DOUBLE PRECISION :: X(N),Y(N),Z(N)
      DOUBLE PRECISION :: XO(N),YO(N),ZO(N)
      REAL             :: RRAND
      LOGICAL          :: HALF,SPIRAL,MOVE,DEBUGGING

      DEBUGGING = .FALSE.
C     DEBUGGING = .TRUE.

      ITERT = ITER
      IF (ITER .LT. 1000) ITERT = 1000

      IRTFLG = 1

      IF (SPIRAL) THEN
         SVAL =  3.6 / SQRT(FLOAT(N))
         PHIT =  0.0

         X(1) =  0.0
         Y(1) =  0.0
         Z(1) = -1.0

         DO K = 2 , N-1
            ZVAL = -1 + 2 * FLOAT(K) / FLOAT(N)
            RVAL = SQRT(1 - ZVAL * ZVAL)
            PHIT = PHIT + SVAL / RVAL

            X(K) = COS(PHIT) * RVAL
            Y(K) = SIN(PHIT) * RVAL
            Z(K) = ZVAL
            IF (DEBUGGING)write(6,90)X(K),Y(K),Z(K),ZVAL,RVAL,SVAL,PHIT
         ENDDO
         X(N) = 0.0
         Y(N) = 0.0
         Z(N) = 1.0

      ELSE

C        CREATE AN INITIAL 'RANDOM' CLOUD 
         DO I=1,N 
            CALL RANDOM_NUMBER(RRAND)
            X(I) =  RRAND * 1000 - 500

            CALL RANDOM_NUMBER(RRAND)
            Y(I) = RRAND * 1000  - 500

            CALL RANDOM_NUMBER(RRAND)
            Z(I) = RRAND * 1000  - 500

            VALL = 1 / SQRT( X(I)*X(I) + Y(I)*Y(I) + Z(I)*Z(I) )
            X(I) = X(I) * VALL
            Y(I) = Y(I) * VALL
            Z(I) = Z(I) * VALL
cc          write(6,90) x(i),y(i),Z(i)
90          format(3(f8.2,' '),' ',3(f8.2,' ')' ',3(f8.2,' '))
         ENDDO
       ENDIF

       IF (MOVE) THEN
         XO     = X
         YO     = Y
         ZO     = Z

         ITERNOW = 1
         DO WHILE (ITERNOW .LT. ITERT)
C           FIND THE CLOSEST TWO POINTS
            MINP1  = 1
            MINP2  = 2
            DMIND  = SQRT( (X(MINP1)-X(MINP2))**2 +
     &                     (Y(MINP1)-Y(MINP2))**2 +
     &                     (Z(MINP1)-Z(MINP2))**2)
            DMAXD  = DMIND

            DO I=1,N
               DO J=I+1,N
                  D = SQRT( (X(I)-X(J))**2 +
     &                      (Y(I)-Y(J))**2 +
     &                      (Z(I)-Z(J))**2)
                   IF (D .LT.  DMIND) THEN
                      DMIND = D 
                      MINP1 = I 
                      MINP2 = J 
                  ENDIF
                   IF (D .GT. DMAXD) THEN
                      DMAXD = D        
                   ENDIF
               ENDDO
            ENDDO

C           MOVE TWO MINIMAL POINTS APART BY 1%
            X(MINP2)  = X(MINP1) + 1.01 * (X(MINP2) - X(MINP1))
            Y(MINP2)  = Y(MINP1) + 1.01 * (Y(MINP2) - Y(MINP1))
            Z(MINP2)  = Z(MINP1) + 1.01 * (Z(MINP2) - Z(MINP1))

            X(MINP1)  = X(MINP1) - 0.01 * (X(MINP2) - X(MINP1))
            Y(MINP1)  = Y(MINP1) - 0.01 * (Y(MINP2) - Y(MINP1))
            Z(MINP1)  = Z(MINP1) - 0.01 * (Z(MINP2) - Z(MINP1))

            BOTT     = SQRT(X(MINP1)*X(MINP1) + Y(MINP1)*Y(MINP1) +
     &                      Z(MINP1)*Z(MINP1))
            IF (BOTT .EQ. 0) THEN
                WRITE(NOUT,*) 'BAD BOTT:',BOTT
                GOTO 9999
            ENDIF
            VALL     = 1/ BOTT

            X(MINP1) = X(MINP1) * VALL
            Y(MINP1) = Y(MINP1) * VALL
            Z(MINP1) = Z(MINP1) * VALL

            BOTT     = SQRT(X(MINP2)*X(MINP2) + Y(MINP2)*Y(MINP2) +
     &                      Z(MINP2)*Z(MINP2))
            IF (BOTT .EQ. 0) THEN
                WRITE(NOUT,*) 'BAD BOTT:',BOTT
                GOTO 9999
            ENDIF

            VALL     = 1 / BOTT

            X(MINP2) = X(MINP2) * VALL
            Y(MINP2) = Y(MINP2) * VALL
            Z(MINP2) = Z(MINP2) * VALL

            ITERNOW  = ITERNOW + 1
          ENDDO
       ENDIF

       NGOT = N

       IF (HALF) THEN
C         REMOVE BOTTEM HALF OF SPHERE
          NGOT = 0
          DO I = 1,N
             IF ( Z(I) .GE. 0.0) THEN

                IF (MOVE .AND. DEBUGGING) THEN
                   D = SQRT( (XO(I)-X(I))**2 +
     &                       (YO(I)-Y(I))**2 +
     &                       (ZO(I)-Z(I))**2)
                   WRITE(NOUT,92) I,D
92                 FORMAT(' DISTANCE MOVED(',I8,'): ',G10.2)
                ENDIF

                NGOT    = NGOT + 1
                X(NGOT) = X(I)
                Y(NGOT) = Y(I)
                Z(NGOT) = Z(I)
             ENDIF
          ENDDO
       ELSEIF (DEBUGGING .AND. MOVE) THEN
          DO I = 1,N
             D = SQRT( (XO(I)-X(I))**2 +
     &                 (YO(I)-Y(I))**2 +
     &                 (ZO(I)-Z(I))**2)
             WRITE(NOUT,92) I,D
          ENDDO
       ENDIF

       IRTFLG = 0

9999   RETURN
       END




