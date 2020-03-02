C++*********************************************************************
C
C  AF.F                               
C                 GETNEWSTACK PARAM.             FEB 03   ARDEAN LEITH
C                 SETPRMB PARAMETERS             MAY 09   ARDEAN LEITH 
C                 GETNEWSTACK PARAM.             OCT 10   ARDEAN LEITH
C                 REMOVED NEWSTACK               FEB 20   ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2020  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C   AF(MAXDIM, LUN1,LUN2, NX,NY,NZ, IRTFLG)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE AF(MAXDIM, LUN1,LUN2, NX,NY,NZ, IRTFLG)

      IMPLICIT NONE
      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      INTEGER    :: MAXDIM,LUN1,LUN2, NX,NY,NZ, IRTFLG 

      REAL       :: BUF
      COMMON        BUF(1)

      INTEGER    :: L,J,ISLICE,MAXX,MAXY,NXNY2,MEMTOT,NE
      INTEGER    :: NOT_USED,NUMB,LB,K1
      REAL       :: AAA,BBB,CCC,DDD,SHXI,SHYI,DET


         IRTFLG = 1

C        FIND MEMORY NEEED 
         K1     = 1 + NX * NY
         MEMTOT = K1 + NX

         IF (MEMTOT  > MAXDIM)  THEN
C           OVERFLOW WILL OCCUR
            CALL  ERRT(6,'AF',NE)
            RETURN
         ENDIF

C        SINGLE IMAGE OPERATION
         AAA = 1.0
         BBB = 0.0
         CALL RDPRM2S(AAA,BBB,NOT_USED,
     &             'TRANSFORMATION PARAMETERS A, B',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         CCC = 0.0
         DDD = 1.0
         CALL RDPRM2S(CCC,DDD,NOT_USED,
     &             'TRANSFORMATION PARAMETERS C, D',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

	 SHXI = 0.0
	 SHYI = 0.0
         CALL  RDPRM2S(SHXI,SHYI,NOT_USED,
     &                'SHIFTS IN X AND Y',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

C        CHECK THAT DETERMINANT OF TRANSFORMATION IS NOT TOO SMALL
	 DET = AAA * DDD - BBB * CCC
	 IF (ABS(DET)  <  1.0E-4)  THEN
	    CALL  ERRT(101,'DETERMINANT OVERFLOW',NE)
	    RETURN
	 ENDIF

C        LOAD AND TRANSFORM SLICE BY SLICE

         DO L=1,NZ
           LB = (L-1)*NY
           DO J=1,NY
              CALL  REDLIN(LUN1,BUF(1+(J-1)*NX),NX,LB+J)
           ENDDO

C          TRANSFORM THIS SLICE
           CALL AFS(BUF,BUF(K1), NX,NY,
     &		    AAA,BBB,CCC,DDD, SHXI,SHYI, LUN2,LB)
         ENDDO

C        RESET FILE HEADER FOR ALTERATIONS IN STATISTICS
         CALL SETPRMB(LUN2, 0.0,0.0, 0.0,0.0)

         IRTFLG = 0

         END

