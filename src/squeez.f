
C++*********************************************************************
C
C SQUEEZ.F              LONG FILENAMES JAN89 al
C                       MIRRORING BUG FIXED JUL 02 ArDean Leith
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
C SQUEEZ(LUNI,LUNO,NSAM,NROW,IRTFLG)
C
C **********************************************************************

	SUBROUTINE SQUEEZ(LUNI,LUNO,NSAM,NROW,NSLICE,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        COMMON /IOBUF/ BUF(NBUFSIZ)

	COMMON BUFOUT(1)

	CHARACTER (LEN=MAXNAM) ::  FILNAM
        COMMON /COMMUN/FILNAM

	PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
	PARAMETER (DGR_TO_RAD = (QUADPI/180))

        IF (IMAMI .NE. 1) CALL NORM3(LUNI,NSAM,NROW,NSLICE,FMAX,FMIN,AV)
        AVSAV = AV

 	CALL RDPRM(ALPHA,NOT_USED,
     &             'ANGLE BETWEEN UNIT VECTORS IN DEGREES')
	IF (ALPHA .GE. 90.0 .OR. ALPHA .LE. 0.0)  THEN
          CALL ERRT(101,'ANGLE MUST BE 0...90 DEGREES',NE)
          RETURN
	ENDIF

C       FIND OUTPUT IMAGE WIDTH (PADDED FOR SHEAR ANGLE)
	CTG   = COS(ALPHA*DGR_TO_RAD)/SIN(ALPHA*DGR_TO_RAD)
	NSAMO = NSAM + INT(ABS(CTG*NROW))

        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNO,'U',IFORM,NSAMO,NROW,NSLICE,
     &              MAXIM,'OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        DO ISLICE = 1,NSLICE
C          FIND I/O RECORD
           IREC1 = (ISLICE - 1) * NROW

	   DO  I = 1,NROW
              CALL REDLIN(LUNI,BUF,NSAM,I+IREC1)

C             FIND SHEAR AMOUNT IN PIXELS
              SH      = (NROW - I) * CTG
              KSH     = SH
              EPS     = SH  - KSH
              ONEMEPS = 1.0 - EPS
    
C             FILL OUTPUT BUFFER WITH AVERAGE
              DO K = 1,NSAMO
                 BUFOUT(K) = AVSAV
              ENDDO

C             FIRST PIXEL INTERPOLATED WITH AVSAV VALUE
              BUFOUT(1+KSH) = AVSAV * EPS + BUF(1) * ONEMEPS

C             INTERPOLATE INTERNAL PIXELS
              DO  K = 2,NSAM
                 BUFOUT(K+KSH) = BUF(K-1) * EPS + BUF(K) * ONEMEPS
              ENDDO

C             LAST PIXEL INTERPOLATED WITH AVSAV VALUE
              IF ((NSAM + KSH) .LE. NSAMO) 
     &           BUFOUT(NSAM+KSH) = BUF(NSAM) * EPS + AVSAV * ONEMEPS

C             WRITE OUTPUT BUFFER TO FILE
              CALL WRTLIN(LUNO,BUFOUT,NSAMO,I+IREC1)
	   ENDDO
        ENDDO

        IRTFLG = 0
	END
