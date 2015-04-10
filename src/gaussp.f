
C ++********************************************************************
C                                                                      *
C                                                                      *
C                                                                      *
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
C                                                                      *
C                                                                      *
C                                                                      *
C  PURPOSE: PLOT FERMI DISTRIBUTIONS AND SUMS OF IT                                                           *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE GAUSSP

	INCLUDE 'CMBLOCK.INC'

        DIMENSION  BUF(256)
        CHARACTER*1    NULL,WHAT,MULAD

        NULL=CHAR(0)

        CALL RDPRMI(IDIM,IDUM,NOT_USED,'PLOT X-DIMENSION')
	IDIM=MIN0(256,IDIM)
        CALL RDPRMC(WHAT,NCHAR,.TRUE.,
     $    '(L)OWPASS, (H)IGHPASS, (B)ANDPASS (L/H/B)',NULL,IRTFLG)
        IF (WHAT.EQ.'H' .OR. WHAT.EQ.'L') 
     $    CALL RDPRM2(RAD,TEMP,NOT_USED,'RADIUS')
        IF (WHAT.EQ.'B') CALL RDPRM2(RAD,RADH,NOT_USED,
     $  'RADIUS FOR LOWPASS,RADIUS FOR HIGHPASS')

        CALL RDPRMC(MULAD,NCHAR,.TRUE.,
     $     '(M)ULTIPLICATIVE,(A)DDITIVE',NULL,IRTFLG)
	RAD=2.*RAD**2
	IF(WHAT.EQ.'B')  RADH=2.*RADH**2
        DO  I=1,IDIM
          X=(I-1)/FLOAT(2*IDIM)
          IF(WHAT.EQ.'L') THEN
	    BUF(I)=EXP(-X**2/RAD)*50.
          ELSEIF(WHAT.EQ.'H') THEN
            BUF(I)=(1-EXP(-X**2/RAD))*50.
	  ELSE
C	IF(WHAT.EQ.'B') THEN
          FLOW= EXP(-X**2/RAD)
          IF(MULAD.EQ.'A') BUF(I)=(FLOW+(1-EXP(-X**2/RADH)))*50.    
          IF(MULAD.EQ.'M') BUF(I)=(FLOW*(1-EXP(-X**2/RADH)))*50.    
         ENDIF
	ENDDO

        CALL MRKUR3(BUF,IDIM,0.,0,60)

        END
