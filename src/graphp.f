
C++*****************************************************************1/17/81
C
C GRAPHP.F
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
C  GRAPHP(LUN1,NSAM,NROW)
C
C  PURPOSE:    TO PRINT/TYPE SINGLE IMAGE ROW
C
C **********************************************************************

	SUBROUTINE GRAPHP(LUN1,NSAM,NROW)

	INCLUDE 'CMBLOCK.INC' 

	COMMON BUF(1)

        NSAV = NDAT
	IF (FCHAR(4:4) .EQ. 'T') THEN
	   NDEV = NOUT
C          (IDEV = 0) IS TERMINAL &  (IDEV = 1) IS LINE PRINTER FORMAT.  
  	   IDEV = 0
        ELSE
           NDEV = NDAT
           IDEV = 1
        ENDIF

C       IMAG INDICATES THE YSIZE SCALING FACTOR TO BE USED.
5	CALL RDPRMI(NLINE,IMAG,NOT_USED,
     &           'ROW TO BE DISPLAYED, SCALE FACTOR')

	IF (NLINE .LE. 0 .OR. NLINE .GT. NROW) THEN
C          COMPLETION, RESTORE THE LUN OF THE PRINTER OUTPUT.
 	   IF (FCHAR(4:4) .EQ. 'T') NDAT = NSAV
	   RETURN
        ENDIF

        IF (IMAG .LE. 0) IMAG = 1

C       READ THE LINE
	CALL REDLIN(LUN1,BUF,NSAM,NLINE)

C       PLOT THE LINE
        NLINES = 1
        FMAG   = IMAG
        CALL GRAPHS(NDEV,BUF,NSAM,NLINES,IDEV,FMAG,IRTFLG)
	GOTO 5

	END
