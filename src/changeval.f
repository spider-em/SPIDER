C++*********************************************************************
C
C    CHANGEVAL.F  
C                     ADDED NREPL                 OCT 2007 ARDEAN LEITH  
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
C    CHANGEVAL(LUN1,LUN2,NSAM,NROW,NSLICE)
C
C    PURPOSE:  CHANGE A SPECIFIC VALUE IN AN IMAGE FILE
C
C    CALLED BY:   UTIL2
C
C--********************************************************************

	SUBROUTINE CHANGEVAL(LUN1,LUN2,NSAM,NROW,NSLICE)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        COMMON /IOBUF/ BUF(NBUFSIZ)

	REAL        NEWVAL

	CALL RDPRM2S(OLDVAL,NEWVAL,NOT_USED,'OLD AND NEW VALUES',IRTFLG)
	IF (IRTFLG .LT. 0)  RETURN

        NREPL = 0
        IRECT = 0

        DO J=1, NSLICE    
	   DO I=1,NROW
              IRECT = IRECT + 1
	      CALL REDLIN(LUN1,BUF,NSAM,IRECT)

	      DO K=1,NSAM
	         IF (BUF(K) .EQ. OLDVAL)  THEN
                    BUF(K) = NEWVAL
                    NREPL  = NREPL + 1
                 ENDIF
               ENDDO

              CALL WRTLIN(LUN2,BUF,NSAM,IRECT)
           ENDDO
        ENDDO

        CALL SETPRM(LUN2,NSAM,NROW,0.,0.,0.,'U')

        FREPL = FLOAT(NREPL)
        CALL REG_SET_NSEL(1,1, FREPL,0.0, 0.0, 0.0, 0.0,IRTFLG)

        WRITE(NOUT,*) ' PIXELS REPLACED: ',NREPL
 	END
