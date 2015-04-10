
C++*********************************************************************
C
C MACF_P.F
C                        ACRS_ CALLS ADDED LS      FEB 2008 ArDean Leith
C                        ACRS PARAMETERS           APR 2009 ArDean Leith
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
C PURPOSE:         CROSS CORRELATION - MASKED AND NORMALIZED 

C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE MACF_P(NSAM,NROW,X, LSD,IRA,MODE,SMD)

        INCLUDE 'CMBLOCK.INC'

        DIMENSION          X(LSD,2*NROW)
        DOUBLE PRECISION   AVE
        CHARACTER*1        MODE,SMD

        R=IRA
        NS2=NSAM/2+1
        NR2=NROW/2+1

c$omp    parallel do private(j,i)
         DO J=1,2*NROW
          DO I=NSAM+1,2*NSAM
           X(I,J)=0.0
	  ENDDO	 
	 ENDDO
	
c$omp    parallel do private(j,i)
         DO J=NROW+1,2*NROW
          DO I=1,NSAM
           X(I,J)=0.0
	  ENDDO
	 ENDDO
         AVE=0.0
         ILE=0

c$omp    parallel do private(j,i,a,tr) reduction(+:AVE,ILE)
         DO J=1,NROW
          A=FLOAT(J-NR2)**2
           DO I=1,NSAM
            TR=SQRT(FLOAT(I-NS2)**2+A)
            IF (TR.GT.R)  THEN
             X(I,J)=0.0
            ELSE
             AVE=AVE+X(I,J)
             ILE=ILE+1
            ENDIF
	   ENDDO
	 ENDDO
         AVE=AVE/ILE
	
c$omp    parallel do private(j,i,a,tr)
         DO J=1,NROW
          A=FLOAT(J-NR2)**2
          DO I=1,NSAM
           TR=SQRT(FLOAT(I-NS2)**2+A)
           IF (TR.LE.R)  X(I,J)=X(I,J)-AVE
	  ENDDO
	 ENDDO 

         INS=1
         CALL FMRS_2(X,2*NSAM,2*NROW,INS)
         IF (INS.EQ.0)  THEN
            CALL ERRT(38,'AC MS ',NE)
            RETURN
         ENDIF

         LS = (2*NSAM)+2 
 	 IF (SMD .EQ. 'S')  THEN
            CALL ACRS_2S(X, LS,2*NSAM,2*NROW)
	 ELSE
            CALL ACRS_2(X, LS,2*NSAM,2*NROW)
	 ENDIF

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
	 D1=1./REAL(NINT(3.1415926*R*R)*ILE)
c$omp    parallel do private(j,i,qt,a,t,m),shared(d1)
         DO J=NRL,NRU
            QT=FLOAT(J-(NROW+1))**2
            DO I=NSL,NSU
               A=SQRT(FLOAT(I-(NSAM+1))**2+QT)/2.0
               IF(A.EQ.0.0)  THEN
                  X(I,J)=X(I,J)*D1
               ELSE
                  IF(R.GT.A)  THEN
                     T=2.0*ATAN(SQRT((R/A)**2-1.0))
C                    Should be NINT without +0.5, but omp won't take it...
                     M=INT(R*R*(T-SIN(T))+0.5)

C                   NORMALIZATION IS APPLIED TO THESE AC COEFF. WHICH WERE
C                   ESTIMATED FROM AT LEAST  5 PIXELS
C                   OTHERWISE AC COEFFS. ARE SET TO ZERO.

                    IF (M.GE.5)  THEN
                       X(I,J)=X(I,J)/FLOAT(M)*ILE
                    ELSE
                       X(I,J)=0.0
                    ENDIF
                 ELSE
                    X(I,J)=0.0
                 ENDIF
              ENDIF
	   ENDDO	
	ENDDO
        END

