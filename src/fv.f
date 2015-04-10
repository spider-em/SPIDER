C++*********************************************************************
C
C FV.F                            AUTHOR: PAWEL PENCZEK
C                                 USED REG_SET AUG 00 AL                        
C                                 ENDLESS BUG SEPT 01 PP
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
C  FV
C
C  PURPOSE: FIND A THRESHOLD IN THE VOLUME THE WILL RESULT
C           IN A SPECIFIED NUMBER OF VOXELS THAT HAVE DENSITIES
C           HIGHER THAN THIS THRESHOLD
C  
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE  FV

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

	CHARACTER(LEN=MAXNAM) :: FILNAM
	REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  VOLIN

	DATA	LUN1/76/

	MAXIM = 0
	CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NSAM,NROW,NSLICE,
     &		   MAXIM,'INPUT',.FALSE.,IRTFLG)
	IF (IRTFLG.NE.0)  RETURN

        ALLOCATE(VOLIN(NSAM,NROW,NSLICE),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL  ERRT(46,'FV, VOLIN',IER)
           RETURN
        ENDIF

	CALL  RDPRMI(ILE,IDUM,NOT_USED,'NUMBER OF VOXELS')

	ILE = MIN(NSAM*NROW*NSLICE,(MAX(1,ILE)))

	DO L=1,NSLICE
	   DO J=1,NROW
              CALL  REDLIN(LUN1,VOLIN(1,J,L),NSAM,J+(L-1)*NROW)
	   ENDDO
	ENDDO
	CLOSE(LUN1)

	THR1 = MAXVAL(VOLIN)
	THR3 = MINVAL(VOLIN)
	THR2 = (THR1-THR3)/2 + THR3
	CALL FTH(VOLIN,NSAM*NROW*NSLICE,THR1,THR2,THR3,ILE,THR)

	DEALLOCATE(VOLIN)
	WRITE(NOUT,101)  THR,ILE
101	FORMAT(' Threshold =',G12.3,'  Number of voxels=',I12)

        CALL REG_SET_NSEL(1,2,THR,FLOAT(ILE),0.0, 0.0,0.0,IRTFLG)

        RETURN
	END



	SUBROUTINE  FTH(V,N,AX,BX,CX,ILE,THR)

	DIMENSION  V(N)
	PARAMETER (R=.61803399,C=1.-R)

	X0 = AX
	X3 = CX
	IF (ABS(CX-BX) .GT. ABS(BX-AX))THEN
           X1 = BX
           X2 = BX+C*(CX-BX)
	ELSE
           X2 = BX
           X1 = BX-C*(BX-AX)
	ENDIF
	LF1 = COUNT(V .GE. X1) - ILE
        F1  = LF1 * REAL(LF1)
	LF2 = COUNT(V .GE. X2) - ILE
        F2  = LF2 * REAL(LF2)

        DO WHILE((.NOT.(LF1.EQ.0.OR.LF2.EQ.0).AND.IABS(LF1-LF2).GE.1) 
     &          .AND.
     &         (ABS(X1-X2).GT.1.0E-5 .AND. ABS(X1-X3).GT.1.0E-5 .AND.
     &         ABS(X2-X3).GT.1.0E-5))

           IF (F2 .LT. F1)THEN
              X0  = X1
              X1  = X2
              X2  = R * X1 + C * X3
              F1  = F2
	      LF2 = COUNT(V .GE. X2) - ILE
	      F2  = LF2 * REAL(LF2)
           ELSE
              X3  = X2
              X2  = X1
              X1  = R * X2 + C * X0
              F2  = F1
	      LF1 = COUNT(V .GE. X1) - ILE
              F1  = LF1 * REAL(LF1)
           ENDIF
	ENDDO

	IF (F1 .LT. F2) THEN
	   ILE = LF1+ILE
           THR = X1
	ELSE
	   ILE = LF2+ILE
           THR = X2
	ENDIF

	END
