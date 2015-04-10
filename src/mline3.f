C ++********************************************************************
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
C SUBROUTINE MLINE TO CREATE A LINE in 3D DEFINED BY TWO END POINTS
C
C	MLINE3(LUN,NSAM,NROW,NSLICE,RP)
C
C	LUN	         : LOGICAL UNIT NUMBER
C	NSAM,NROW,NSLICE : FILE DIMENSIONS
C	RP	         : VALUE TO BE CORRECTED
C
C **********************************************************************

	SUBROUTINE MLINE3(LUN,NSAM,NROW,NSLICE,RP)

	INCLUDE 'CMBLOCK.INC' 

	COMMON BUF(1)
	INTEGER  X1,X2,Y1,Y2,Z1,Z2

	CALL RDPRMI(X1,Y1,NOT_USED,
     &             'X,Y COORDINATES OF FIRST POINT')
	CALL RDPRMI(Z1,IDUM,NOT_USED,
     &             'Z COORDINATE OF FIRST POINT')
	CALL RDPRMI(X2,Y2,NOT_USED,
     &             'X,Y COORDINATES OF SECOND POINT')
	CALL RDPRMI(Z2,IDUM,NOT_USED,
     &             'Z COORDINATE OF SECOND POINT')


	IF(X2.EQ.X1)  THEN
	  IF(Y2.EQ.Y1)  THEN
	    IF(Z2.GE.Z1)  THEN
            IC=1
            ELSE
            IC=-1
            ENDIF
            DO  IZ=Z1,Z2,IC
            IREC=(IZ-1)*NROW+Y1
            CALL REDLIN(LUN,BUF,NSAM,IREC)
	    BUF(X1)=RP
            CALL WRTLIN(LUN,BUF,NSAM,IREC)
            ENDDO
          ELSE
            IF(Y2.GE.Y1)  THEN
            IC=1
            ELSE
            IC=-1
            ENDIF
            DO  IY=Y1,Y2,IC
            IZ=REAL( (Z2-Z1)*(IY-Y1) )/REAL(Y2-Y1) + Z1
             IF(IZ.GT.0.AND.IZ.LE.NSLICE)  THEN
	      IREC=(IZ-1)*NROW + IY
              CALL REDLIN(LUN,BUF,NSAM,IREC)
	      BUF(X1)=RP
              CALL WRTLIN(LUN,BUF,NSAM,IREC)
             ENDIF
            ENDDO
          ENDIF
        ELSE
          IF(X2.GE.X1)  THEN
          IC=1
          ELSE
          IC=-1
          ENDIF
	  DO  IX=X1,X2,IC
           IY=REAL( (Y2-Y1)*(IX-X1) )/REAL(X2-X1) + Y1
           IF(IY.GT.0.AND.IY.LE.NROW)  THEN
             IZ=REAL( (Z2-Z1)*(IX-X1) )/REAL(X2-X1) + Z1
             IF(IZ.GT.0.AND.IZ.LE.NSLICE)  THEN
	      IREC=(IZ-1)*NROW + IY
              CALL REDLIN(LUN,BUF,NSAM,IREC)
	      BUF(IX)=RP
              CALL WRTLIN(LUN,BUF,NSAM,IREC)
             ENDIF
           ENDIF
          ENDDO
        ENDIF
	END
