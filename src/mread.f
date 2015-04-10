
C++*******************************************************************
C
C    MREAD.FOR
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
C    MREAD(LUN,BUF,NSAM,NREC,NPTR)
C 
C       THIS SUBROUTINE IS USED FOR SCROLLING THROUGH AN IMAGE
C	WHERE ONLY NREC IMAGE ROWS ARE KEPT IN THE BUFFER BUF
C	AT A TIME. EACH CALL TO MREAD CAUSES ONE ROW TO BE READ IN
C	TO REPLACE THE 'OLDEST' ROW CURRENTLY IN THE BUFFER. THE
C	ADDRESSES OF THE ROWS IN THE BUFFER ARE BEING ROTATED
C 	ACCORDING TO THE PERMUTATIONS OF THE NUMBERS 1...NREC
C	AN INDEX ARRAY NPTR MAKES THE CURRENT ADDRESSES AVAILABLE.
C
C  PARAMETERS:
C  LUN		LOGICAL UNIT NUMBER OF IMAGE FILE
C  BUF		BUFFER ARRAY OF LENGTH AT LEAST NSAM*NREC
C  NSAM		NUMBER OF ELEMENTS PER ROW
C  NREC		NUMBER OF ROWS IN IMAGE SEGMENT TO BE READ IN
C  NPTR		ARRAY OF LENGTH NREC. NPTR(I) POINTS TO THE LOCATION
C		OF THE I-TH ROW IN THE BUFFER: E.G. IMAGE ROW #8 OF 
C		CURRENT SEGMENT STARTS IN BUFFER LOCATION (NPTR(I)-1)*NSAM+1
C
C	MREAD HAS TO BE INITIALIZED BY 
C			CALL MREAD(-1,BUF,NSAM,NREC,NPTR)
C	THE FIRST CALL WITH NON-NEGATIVE LUN WILL READ IN ALL NREC
C	RECORDS, AND INITIALIZE THE INDEX ARRAY.
C
C 	EACH SUBSEQUENT CALL TO MREAD WILL CAUSE ONE ADDITIONAL ROW
C	TO BE READ IN, WHICH WILL TAKE THE PLACE OF THE "OLDEST" ROW.
C	E.G., THE SECOND CALL WILL READ IN ROW #(NREC+1) AND STORE IT
C	IN THE BUFFER AT THE LOCATION PREVIOUSLY OCCUPIED BY THE FIRST
C	IMAGE ROW, AND SO ON.
C
C--*******************************************************************

      SUBROUTINE MREAD(LUN,BUF,NSAM,NREC,NPTR)

 

C     I DO NOT KNOW IF SAVE IS NEEDED FEB 99 al
      SAVE

      DIMENSION NPTR(*),BUF(*)

      IF (LUN .LE. 0) THEN
         JR    = 0
         NCALL = 0

      ELSEIF (NCALL .LE. 0) THEN
         NCALL = 1
         N     = 1

         DO  I = 1,NREC
            NPTR(I) = I
            JR      = JR + 1
            CALL REDLIN(LUN,BUF(N),NSAM,JR)
            N = N + NSAM
	 ENDDO

      ELSE
         N = (NPTR(1)-1)*NSAM+1
         DO  I = 1,NREC
             NPTR(I) = MOD(NPTR(I)+1,NREC)
             IF (NPTR(I) .EQ. 0) NPTR(I) = NREC
	 ENDDO
         JR = JR + 1
         CALL REDLIN(LUN,BUF(N),NSAM,JR)
      ENDIF

      RETURN
      END

