
C++*********************************************************************
C
C LHIST.FOR
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
C  PURPOSE: LOCAL HISTOGRAM GENERATION
C
C  LHIST(BO,NSAM,NLOCAL,KCTR,NPTR,H,MODE)
C         BO
C         NSAM       NUMBER OF SAMPLES
C         NLOCAL     LOCAL AREA DIMENSION (X IS = TO Y)
C         KCTR       COLUMN POINTER
C         NPTR       ROW    POINTER
C         H          HISTOGRAM BUFFER
C         MODE       1 NON-INCREMENTAL OPERATION, NO INTEGRATION
C                    2 INCREMENTAL OPERATION, NO INTEGRATION
C                    3 NON-INCREMENTAL OPERATION WITH INTEGRATION
C                    4 INCREMENTAL OPERATION WITH INTEGRATION
C*****************************************************

      SUBROUTINE LHIST(B0,NSAM,NLOCAL,KCTR,NPTR,H,MODE)

 
      INCLUDE 'CMBLOCK.INC'

      REAL      H(1024),B0(*)
      DIMENSION NPTR(1)

      SAVE FVAL

      IF (MODE .EQ. 0) THEN
        IF (FMAX .EQ. FMIN) THEN
          CALL ERRT(5,'LHIST ',NE)
          MODE = -1
          RETURN
        ENDIF

        FVAL = 127./ (FMAX-FMIN)
        FN   = 1.  / FLOAT(NLOCAL)**2
        RETURN
      ENDIF

      KCTR1 = KCTR-NLOCAL / 2
      GOTO (10,200,8,200),MODE
     
C     INITIALIZE HISTOGRAM
8     DO  K = 129,256
        H(K) = 0.
      ENDDO

10    DO  K = 1,128
        H(K) = 0.
      ENDDO

      DO  I = 1,NLOCAL
         II = (NSAM * (NPTR(I) - 1) + KCTR1) - 1
         DO  K = 1,NLOCAL
            J = INT((B0(II+K)-FMIN)*FVAL) + 1.5
            IF (J .LE. 1) THEN
              H(1) = H(1) + 1.0
            ELSEIF (J. GE. 128) THEN
              H(128) =H(128) + 1.0
            ELSE
              H(J) = H(J) + 1.0
            ENDIF
	 ENDDO
      ENDDO

      IF (MODE .EQ. 3) THEN
         H(129) = H(1)
         DO  K = 2,128
           H(128+K) = H(127+K) + H(K)
	 ENDDO
      ENDIF
      RETURN



C     INCREMENTAL UPDATING OF HISTOGRAM
200   DO  I = 1,NLOCAL
         NRUN =  1.0
         HADD = -1.0
C        FIRST COLUMN SUBTRACTED IN FIRST PASS, LAST COL. ADDED IN SECOND PASS
         II = NSAM*(NPTR(I)-1)+KCTR1

201      J = INT((B0(II)-FMIN)*FVAL)+1.5
         IF (J .LE. 1) THEN
           H(1) = H(1) + HADD
         ELSEIF (J .GE. 128) THEN
           H(128) = H(128) + HADD
         ELSE
           H(J) = H(J) + HADD
         ENDIF

         IF (NRUN .EQ. 1) THEN
            II = II + NLOCAL -1
            NRUN = 2
            HADD = +1.
            GOTO 201
         ENDIF
      ENDDO


      IF (MODE .EQ. 4) THEN
         H(129) = H(1)
         DO  K = 2,128
           H(128+K) = H(127+K)+H(K)
	 ENDDO
      ENDIF

      RETURN
      END
