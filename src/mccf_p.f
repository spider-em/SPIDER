C++*********************************************************************
C
C MCCF_P.F
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
C PURPOSE:  CROSS CORRELATION - MASKED AND NORMALIZED 
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE MCCF_P(NSAM,NROW,X,Y,LSD,IRA,MODE)

        INCLUDE 'CMBLOCK.INC'

        DIMENSION         X(LSD,2*NROW),Y(LSD,2*NROW)
        DOUBLE PRECISION  AVE
        CHARACTER*1         MODE

        R   = IRA
        NS2 = NSAM/2+1
        NR2 = NROW/2+1

c$omp    parallel do private(j,i)
         DO J=1,2*NROW
          DO I=NSAM+1,2*NSAM
           X(I,J)=0.0
           Y(I,J)=0.0
          ENDDO  
         ENDDO
        
c$omp    parallel do private(j,i)
         DO J=NROW+1,2*NROW
          DO I=1,NSAM
           X(I,J)=0.0
           Y(I,J)=0.0
          ENDDO
         ENDDO
         AVE=0.0
         ILE=0
         AVEY=0.0
         ILEY=0
c$omp    parallel do private(j,i,a,tr) reduction(+:AVE,ILE,AVEY,ILEY)
         DO J=1,NROW
          A=FLOAT(J-NR2)**2
           DO    I=1,NSAM
            TR=SQRT(FLOAT(I-NS2)**2+A)
            IF(TR.GT.R)  THEN
             X(I,J)=0.0
             Y(I,J)=0.0
            ELSE
             AVE=AVE+X(I,J)
             ILE=ILE+1
             AVEY=AVEY+X(I,J)
             ILEY=ILEY+1
            ENDIF
           ENDDO
         ENDDO
         AVE=AVE/ILE
         AVEY=AVEY/ILEY
        
c$omp    parallel do private(j,i,a,tr)
         DO J=1,NROW
          A=FLOAT(J-NR2)**2
          DO I=1,NSAM
           TR = SQRT(FLOAT(I-NS2)**2+A)
           IF (TR.LE.R)  THEN
                X(I,J)=X(I,J)-AVE
                Y(I,J)=Y(I,J)-AVEY
           ENDIF
          ENDDO
         ENDDO 

         INS=1
         CALL FMRS_2(X,2*NSAM,2*NROW,INS)
         IF (INS .EQ. 0)  THEN
            CALL ERRT(38,'CC MS ',NE)
            RETURN
         ENDIF

         INS=1
         CALL FMRS_2(Y,2*NSAM,2*NROW,INS)
         IF (INS .EQ. 0)  THEN
            CALL ERRT(38,'CC MS ',NE)
            RETURN
         ENDIF

         LSC = 2*NSAM+2-MOD(2*NSAM,2)
         CALL CCRS_2I(X,Y, LSC,2*NSAM,2*NROW)

         IF (MODE .EQ. 'F')  THEN
           NRL=1
           NRU=2*NROW
           NSL=1
           NSU=2*NSAM
         ELSE
           NRL=NR2
           NRU=NROW+NR2-1
           NSL=NS2
           NSU=NSAM+NS2-1
         ENDIF
        
c$omp    parallel do private(j,i,qt,a,t,m)
         DO J=NRL,NRU
            QT = FLOAT(J-(NROW+1))**2
            DO I=NSL,NSU
               A = SQRT(FLOAT(I-(NSAM+1))**2+QT)/2.0
               IF (A .EQ. 0.0)  THEN
                  X(I,J) = X(I,J)/NINT(3.1415926*R*R)*ILE
               ELSE
                  IF (R.GT.A)  THEN
                     T = 2.0*ATAN(SQRT((R/A)**2-1.0))
                     M = NINT(R*R*(T-SIN(T)))

C                    NORMALIZATION IS APPLIED TO THESE AC COEFF. WHICH 
C                    WERE ESTIMATED FROM AT LEAST  5 PIXELS
C                    OTHERWISE AC COEFFS. ARE SET TO ZERO.

                    IF (M.GE.5)  THEN
                       X(I,J) = X(I,J)/FLOAT(M)*ILE
                    ELSE
                       X(I,J) = 0.0
                    ENDIF
                 ELSE
                    X(I,J) = 0.0
                 ENDIF
              ENDIF
           ENDDO        
        ENDDO
        END
