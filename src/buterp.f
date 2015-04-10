C ++********************************************************************
C                                                                      *
C BUTERP                                                               *
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
C
C***********************************************************************

        SUBROUTINE BUTERP

        INCLUDE 'CMBLOCK.INC'

        COMMON       A(80),BUF(1024)
        CHARACTER    NULL,WHAT,ANS

        EPS = 0.882
        AA  = 10.624

        NULL = CHAR(0)

        CALL RDPRMI(IDIM,IDUM,NOT_USED,'PLOT X-DIMENSION')
        CALL RDPRMC(WHAT,NCHAR,.TRUE.,
     &    '(L)OWPASS, (H)IGHPASS, (L/H)',NULL,IRTFLG)

        IF (WHAT.EQ.'H' .OR. WHAT.EQ.'L') 
     &     CALL RDPRM2(RAD1,RAD2,NOT_USED,
     &	              'LOWER & UPPER LIMITING FREQUENCIES')

        ORD=2.*ALOG10(EPS/SQRT(AA**2-1.0))
        ORD=ORD/ALOG10(RAD1/RAD2)
        RAD=RAD1/(EPS)**(2./ORD)

        DO  I=1,IDIM
          XX=0.5*(I-1)/FLOAT(IDIM-1)
          IF(WHAT.EQ.'L') BUF(I)=SQRT(1./(1.+(XX/RAD)**ORD))*50. 
          IF(WHAT.EQ.'H') BUF(I)=(1-SQRT(1.0/(1.0+(XX/RAD)**ORD)))*50.    
	ENDDO

        CALL MRKUR3(BUF,IDIM,0.,0,60)

        CALL RDPRMC(ANS,NC,.TRUE.,
     &    'LIKE AN EXAMPLE WITH STEP FUNCTION?(Y/N)',NULL,IRT)

	IF (ANS .EQ. 'Y') THEN

	   NDIM=LOG2(IDIM)
	   MDIM=2**NDIM
	   IF (MDIM.NE.IDIM)THEN
	      WRITE(NOUT,*)
     &          'WORKS FOR DIMENSIONS EQUAL TO POWERS OF TWO ONLY'
              CALL ERRT(100,'BUTERP',NET)
	      RETURN
	   ENDIF

           DO I=1,IDIM
	      BUF(I)=1.0*50.0
	      IF (I.GT.IDIM/2)BUF(I)=0.0
	   END DO
 
	   WRITE(NOUT,*) 'STEP FUNCTION'
 
           CALL MRKUR3(BUF,IDIM,0.,0,60)
 
           CALL FFTR_Q(BUF,NDIM)
	   X1=IDIM**2

           DO  I=1,IDIM
              IF(I.LE.2)IX=(I-1)*64
              IX=(I-1)/2
              ARG=SQRT(FLOAT(IX*IX)/X1)

              IF (WHAT.EQ.'L')
     &           BUF(I)=BUF(I)*SQRT(1./(1.+(ARG/RAD)**ORD))
              IF (WHAT.EQ.'H') BUF(I)=
     &           BUF(I)*(1.-SQRT(1./(1.+(ARG/RAD)**ORD)))
	  ENDDO

          CALL FFTR_Q(BUF,-NDIM)

	  WRITE(NOUT,*) 'FILTERED FUNCTION'

	  CALL MRKUR3(BUF,IDIM,0.,0,60)
	
	ENDIF
	END
