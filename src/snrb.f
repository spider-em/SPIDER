
C **********************************************************************
C *  SNRB.F
C=**********************************************************************
C=* Copyright (C)2002, L. Joyeux & P. A. Penczek                       *
C=* University of Texas - Houston Medical School                       *
C=* Email:  pawel.a.penczek@uth.tmc.edu                                *
C=*                                                                    *
C=* This file is part of:                                              * 
C=* SPIDER - Modular Image Processing System.   Author: J. FRANK       *
C=* Copyright 1985-2009  Health Research Inc.,                         *
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
C=* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
C=* General Public License for more details.                           *
C=*                                                                    *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program. If not, see <http://www.gnu.org/licenses> *
C=*                                                                    *
C=**********************************************************************

      SUBROUTINE SNRB

      IMPLICIT   NONE

      INCLUDE 'CMBLOCK.INC' 
      INCLUDE 'CMLIMIT.INC' 
      INCLUDE 'F90ALLOC.INC'

      REAL, DIMENSION(:,:), POINTER :: PANG, FREQ
      CHARACTER (LEN=1)             :: ANS,NULL 
      CHARACTER (LEN=MAXNAM)        :: FILE

      INTEGER          MAXXT, MAXYT, IRTFLG, NA
      REAL             FP, FS
      LOGICAL          FITLOW
      
      MAXXT = 0
      MAXYT = 0
      CALL GETDOCDAT('FSC DOC',.TRUE.,FILE,
     &      77,.TRUE.,MAXXT,MAXYT,PANG,IRTFLG)

      ALLOCATE(FREQ(2,MAXYT))
      FREQ(1,:) = PANG(2,:)
      FREQ(2,:) = PANG(4,:)

      DEALLOCATE(PANG)
      
      NULL = CHAR(0)
      CALL RDPRMC(ANS,NA,.TRUE.,
     &       'LOW PASS / HIGH PASS (L/H)',NULL,IRTFLG)    
     
      FP = 0.1
      FS = 0.9
      FITLOW = (ANS .EQ. 'L') .OR. (ANS .EQ. 'L')
      CALL FITBUFFCT(FREQ, MAXYT, FP, FS, FITLOW)

      CALL REG_SET_NSEL(1,2,FP,FS,0.0, 0.0, 0.0, IRTFLG)

      DEALLOCATE(FREQ)
      RETURN
      END
      
CCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE FITBUFFCT(FREQ, NB, FP, FS, FITLOW)
      IMPLICIT   NONE
      
      INTEGER    NB
      REAL      FREQ(2,NB), FP, FS
      REAL      EPS,AA
      REAL      DELTA, MIN, BUT
      LOGICAL   FITLOW
      INTEGER  MINJ, MINI, I, J, K
      REAL      FPN, FSN
      REAL      BUTDIFF
      
      DELTA = 0.05

      MIN = BUTDIFF(FREQ, NB, FP, FS, FITLOW)

      DO K=1,20
         DO
            MINI = 0
            MINJ = 0
            DO I=-1,1
               DO J=-1,1
                  FPN = FP+DELTA*I
                  FSN = FS+DELTA*J
                  IF(FPN.LE.0 .OR. FSN.GE.1.0) GOTO 100
                  BUT = BUTDIFF(FREQ, NB, FPN, FSN , FITLOW)
                  IF(BUT<MIN) THEN
                     MINI = I
                     MINJ = J
                     MIN = BUT
                  END IF
100               CONTINUE
               END DO
            END DO
C            WRITE(NOUT,*) __FILE__, __LINE__, " MINI MINJ ", MINI, MINJ
            IF(MINI.EQ.0 .AND. MINJ.EQ.0) GOTO 200
            FP = FP + DELTA*MINI
            FS = FS + DELTA*MINJ
         END DO
200      DELTA = DELTA / SQRT(2.0)
C         WRITE(NOUT,*) __FILE__, __LINE__, " FP FS ", FP, FS
      END DO
      
      RETURN 
      END
      
CCCCCCCCCCCCCCCCCCCCCC
      
      FUNCTION BUTDIFF(FREQ, NB, FP, FS, FITLOW)
      IMPLICIT   NONE

      REAL BUTDIFF
      INTEGER    NB
      REAL      FREQ(2,NB), FP, FS, DIFF
      LOGICAL   FITLOW

       REAL   EPS, AA
      PARAMETER(EPS=0.882)
      PARAMETER(AA=10.624)      
      INTEGER I
      REAL   ORD, RAD
      REAL   A, B

      ORD = 2*ALOG10(EPS/SQRT(AA**2-1))/ALOG10(FP/FS)
      RAD = FP/EPS**(2/ORD)
      
      IF(FITLOW) THEN
         A = 0
         B = 1
      ELSE
         A = 1
         B = -1
      END IF
      
      BUTDIFF = 0 
      DO I=1, NB
         DIFF = FREQ(2,I) - (A+B/SQRT(1+(FREQ(1,I)/RAD)**ORD))
         IF(DIFF>0) DIFF = DIFF * 5
         BUTDIFF = BUTDIFF + DIFF**2
      END DO
      
      RETURN
      END
#if 0
BUF(X)=SQRT(1/(1+(X/RAD(FP,FS))**ORD(FP,FS)))
RAD(FP,FS)=FP/EPS**(2/ORD(FP,FS))
ORD(FP,FS)=2*LOG10(EPS/SQRT(AA**2-1))/LOG10(FP/FS)
H(X)=(SGN(X)+1)/2
FSC2SNR(X)=X/(1-X)
S0(X)=H(X)*X
#endif
      
