
C++*********************************************************************
C
C COPYR.FOR                  LONG FILENAMES       JAN 89 ARDEAN LEITH
C                            OPFILEC              FEB 03 ARDEAN LEITH
C                            REWRITE             JAN 14 ARDEAN LEITH
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C COPYR(LUN2
C
C--*********************************************************************

	SUBROUTINE COPYR(LUN2)

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC' 

        INTEGER                :: LUN2

        CHARACTER (LEN=MAXNAM) :: FILNAM
	REAL                   :: BUF(NBUFSIZ)  ! FROM CMLIMIT
        INTEGER                :: MAXIM,NX,NY,NZ,ITYPE,IRTFLG,I,NGOT,NE
        INTEGER                :: J,NC
C                                            123456789 123456789 123456789 
        CHARACTER (LEN=MAXNAM) :: PROMPTS = 'SLICE('
        CHARACTER (LEN=MAXNAM) :: PROMPTL = 'ROW('

        MAXIM  = 0
        NX     = 0
        NY     = 0
        NZ     = 0
        ITYPE  = 3
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'U',IFORM,NX,NY,NZ,
     &             MAXIM,'OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        WRITE(NOUT,*) ' '

	DO J=1,NZ              ! LOOP OVER SLICES IN IMAGE/VOL.

          IF (NZ > 1) THEN
             CALL INTTOCHAR(J,PROMPTS(7:),NC,1)
             WRITE(NOUT,90) PROMPTS(1:6+NC)
 90          FORMAT(2X,A,') -----')
          ENDIF

	  DO I=1,NY           ! LOOP OVER LINES IN IMAGE/VOL.
            CALL INTTOCHAR(I,PROMPTL(5:),NC,1)
            PROMPTL(5+NC:5+NC) = ')'

            CALL RDPRA(PROMPTL(1:5+NC), NX,0,.FALSE.,BUF,NGOT,IRTFLG)

            IF (IRTFLG .NE. 0 ) THEN
              CALL ERRT(101,'INPUT CAN NOT BE UNDERSTOOD',NE)
              EXIT
            ELSEIF (NGOT < NX) THEN
              CALL ERRT(102,'INSUFFICIENT VALUES FOR THIS LINE',NGOT)
              EXIT
            ENDIF

            CALL WRTLIN(LUN2,BUF,NX,(J-1)*NY+I)

C           WRITE(NOUT,9000)(BUF(K),K=1,NX)
9000        FORMAT(16F8.3)
          ENDDO
        ENDDO

        WRITE(NOUT,*) ' '
        CLOSE(LUN2)

	END

                                    
