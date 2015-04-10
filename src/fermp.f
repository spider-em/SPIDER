
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
C  PURPOSE:  PLOT FERMI DISTRIBUTIONS AND SUMS OF IT
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C***********************************************************************

        SUBROUTINE FERMP

        INCLUDE 'CMBLOCK.INC'
        COMMON A(80),BUF(1024)
        CHARACTER *1 WHAT,MULAD,NULL

        NULL=CHAR(0)

        CALL RDPRMI(IDIM,IDUM,NOT_USED,'PLOT X-DIMENSION')
6666    CALL RDPRMC(WHAT,NCHAR,.TRUE.,
     $     '(L)OWPASS, (H)IGHPASS, OR (B)ANDPASS? (L/H/B)',NULL,IRT)

        IF (WHAT .EQ. 'H' .OR. WHAT .EQ. 'L') THEN
            CALL RDPRM2(RAD,TEMP,NOT_USED,
     $     'FERMI CUTOFF RADIUS, TEMP. FACTOR')

        ELSEIF (WHAT .EQ. 'B') THEN
           CALL RDPRM2(RAD,TEMP,NOT_USED,
     $     'FERMI CUTOFF RADIUS, TEMP. FACTOR FOR LOWPASS')
           CALL RDPRM2(RADH,TEMPH,NOT_USED,
     $     'FERMI CUTOFF RADIUS, TEMP. FACTOR FOR HIGHPASS')
           CALL RDPRMC(MULAD,NCHAR,.TRUE.,
     $     '(M)ULTIPLICATIVE, OR (A)DDITIVE? (M/A)',NULL,IRT)
        ELSE
           GOTO 6666
        ENDIF


        DO  I=1,IDIM
          X=(I-1)/FLOAT(2*IDIM)
          IF (WHAT.EQ.'L') BUF(I)=1./(1.+EXP((X-RAD)/TEMP))*50.
          IF (WHAT.EQ.'H') BUF(I)=1./(1.+EXP(-(X-RAD)/TEMP))*50.
          IF (WHAT.EQ.'B') THEN
             FLOW=1./(1.+EXP((X-RAD)/TEMP))
             IF (MULAD .EQ. 'A') 
     $          BUF(I) = (FLOW+1./(1.+EXP(-(X-RADH)/TEMPH)))*50.
             IF (MULAD .EQ. 'M') 
     $          BUF(I) = (FLOW*1./(1.+EXP(-(X-RADH)/TEMPH)))*50.
          ENDIF
	ENDDO

        CALL MRKUR3(BUF,IDIM,0.,0,60)

        END
