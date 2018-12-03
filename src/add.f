


C++*********************************************************************
C
C  ADD.F                    CHANGED              7/21/86 MN
C                           OUTPUT FILES CHANGED AUG 96  ARDEAN LEITH
C                           REWRITTEN            MAR 99  ARDEAN LEITH
C                           DIV BUG              SEP 03  ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C  ADD(VOLBUF,LUNIN,ITYPE,NX,NY,NZ,SIGT)
C
C  PURPOSE:  ADD, SUBTRACT, OR MULTIPLY IMAGES
C            OR MULTIPLY FOURIER TRANSFORMS
C            BOTH IMAGES HAVE TO HAVE THE SAME SIZE.
C
C  PARAMETERS:     
C        VOLBUF       INPUT VOLUME (#1)
C        LUNIN        I/O UNIT NUMBER OF INPUT FILE #2
C        ITYPE       IFORM OF INPUT VOLUME
C        NX,NY    X & Y DIMENSIONS OF IMAGES
C        NZ       Z DIMENSION OF IMAGES
C        SIGT        +1    1 IS ADDED TO 2
C                    -1    2 IS SUBTRACTED FROM 1 
C                    +2    1 IS MULTIPLIED WITH 2 
C                    -2    2 IS DIVIDED BY 1,
C                          OR COMPLEX FOURIER MULTIPLICATION WITH CONJUGATE
C                    -3    COMPLEX 2 IS DIVIDED BY COMPLEX 1
C                    +5    ARITHMETIC OR OF 1 WITH 2 
C                          
C    VARIABLES: IFORM (TYPE)  FILE TYPE SPECIFIER. 
C	         +1    R    2-D IMAGE
C                +3   R3    3-D IMAGE
C               -11   O2    2-D FOURIER TRANSFORM, MIXED RADIX ODD
C               -12   E2    2-D FOURIER TRANSFORM, MIXED RADIX EVEN
C               -21   O3    3-D FOURIER TRANSFORM, MIXED RADIX ODD
C               -22   E3    3-D FOURIER TRANSFORM, MIXED RADIX EVEN
C
CC23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE ADD(VOLBUF,LUNIN,ITYPE,NX,NY,NZ,SIGN)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL       :: VOLBUF(1)
        INTEGER    :: LUNIN,ITYPE,NX,NY,NZ 
        REAL       :: SIGN

        REAL       :: SIGT

        COMMON /IOBUF/ BUF(NBUFSIZ)

C       DOES NOT WORK ON SOME ODD FILE FORMATS NO LONGER IN USE
        IF (ITYPE == 0  .OR. ITYPE == 8  .OR. ITYPE == 11 .OR.
     &	    ITYPE == 12 .OR. ITYPE == 16 .OR. ITYPE == -9)  THEN
            CALL ERRT(102,'ADD; UNSUPPORTED FORMAT',ITYPE)
            RETURN
        ENDIF

C       SIGT MAY BE A CONSTANT
        SIGT = SIGN

        IF ( INDEX(FCHAR(4:) ,'D') > 0 .OR. 
     &             FCHAR(1:2) == '12')   THEN
C           DIVIDE
            IF (ITYPE < 0) THEN
C              FOURIER FILES
               SIGT = -3.0
               WRITE(NOUT,*) 'COMPLEX DIVISION'
            ELSE
               SIGT = -2.0
               WRITE(NOUT,*) 'DIVISION'
            ENDIF

        ELSEIF ( INDEX(FCHAR(4:) ,'M') > 0 ) THEN
           SIGT = -2.0
           WRITE(NOUT,*) 'COMPLEX MULTIPLICATION -- (X) * CONJG(Y)'

        ELSEIF ( INDEX(FCHAR(4:) ,'O') > 0 ) THEN
           SIGT = +5.0
           WRITE(NOUT,*) 'ARITHMETIC OR'
        ENDIF


	NREC  = NY * NZ
        ILOC  = 0

	IF (SIGT == +1.0)  THEN
C           SIGT        +1    2 IS ADDED TO 1 --------------------- ADD
            DO IREC=1,NREC
		CALL REDLIN(LUNIN,BUF,NX,IREC)

            	DO  IX=1,NX
                   ILOC         = ILOC + 1
  		   VOLBUF(ILOC) = VOLBUF(ILOC) + BUF(IX)
                ENDDO
            ENDDO
           	
	ELSEIF (SIGT == -1.0)  THEN
C          SIGT        -1    2 IS SUBTRACTED FROM 1 --------------- SUB
           DO IREC=1,NREC
	       CALL REDLIN(LUNIN,BUF,NX,IREC)
          	DO  IX=1,NX
                  ILOC         = ILOC + 1
                  VOLBUF(ILOC) = VOLBUF(ILOC) - BUF(IX)
               ENDDO
           ENDDO

	ELSEIF (SIGT == +2.0)  THEN
C          SIGT        +2    2 IS MULTIPLIED BY 1 ---------------- MULT

	   IF (IFORM > 0)  THEN
C              REAL FILES
               DO IREC=1,NREC
		  CALL REDLIN(LUNIN,BUF,NX,IREC)
            	  DO  IX=1,NX
                     ILOC         = ILOC + 1
  		     VOLBUF(ILOC) = VOLBUF(ILOC) * BUF(IX)
                  ENDDO
              ENDDO

	   ELSEIF (IFORM < 0)  THEN
C             FOURIER FILES
              CALL CADD(VOLBUF,LUNIN,ITYPE,NX,NY,NZ,SIGT)	
           ENDIF
	
	ELSEIF (SIGT  ==  -2.0)  THEN
C          SIGT        -2    1 IS DIVIDED BY 2 ----------------- DIVIDE
C          OR COMPLEX FOURIER MULTIPLICATION WITH CONJUGATE
           WRITE(NOUT,*) 'DIVISION'

	   IF (IFORM > 0)  THEN
	       IZCOUN = 0
               DO IREC=1,NREC
		  CALL REDLIN(LUNIN,BUF,NX,IREC)
            	  DO  IX=1,NX
                     ILOC         = ILOC + 1
		     IF (BUF(IX) .NE. 0) THEN
		        VOLBUF(ILOC) = VOLBUF(ILOC) / BUF(IX)
		     ELSE
		        VOLBUF(ILOC) = 0.0
		        IZCOUN = IZCOUN + 1
		     ENDIF
                  ENDDO
              ENDDO

	      IF (IZCOUN .GT. 0) WRITE(NOUT,160) IZCOUN
160	      FORMAT(' --- WARNING: ',I7,' OUTPUT PIXELS SET TO 0.0 ',
     &               'WHEN DIVISION BY 0 ATTEMPTED')

	   ELSEIF (IFORM .LT. 0)  THEN
C             COMPLEX FOURIER MULTIPLICATION WITH CONJUGATE

              CALL CADD(VOLBUF,LUNIN,ITYPE,NX,NY,NZ,SIGT)	
	   ENDIF


   	ELSEIF (SIGT == -3.0)  THEN
C          SIGT        -3    COMPLEX 2 IS DIVIDED BY COMPLEX 1 ---- DIV
           WRITE(NOUT,*) 'COMPLEX DIVISION'

           CALL CADD(VOLBUF,LUNIN,ITYPE,NX,NY,NZ,SIGT)	


	ELSEIF (SIGT == +5.0)  THEN
C          SIGT        +5   ARITHMETIC OR OF 1 WITH 2 -------------- OR

	   IF  (IFORM .GT. 0)  THEN
               DO IREC=1,NREC
		  CALL REDLIN(LUNIN,BUF,NX,IREC)
            	  DO  IX=1,NX
                     ILOC = ILOC + 1
                     B    = BUF(IX)
                     IF (B  ==  0.0) B = VOLBUF(ILOC)
                     VOLBUF(ILOC) = B
                  ENDDO
              ENDDO

	   ELSE
C              NOT IMPLEMENTED FOR FOURIER
	       CALL ERRT(101,'ADD, OR NOT IMPLEMENTED FOR FOURIER',NE)
	   ENDIF
	ELSE
C          UNKNOWN SIGT
	   CALL ERRT(23,'ADD',NE)
	ENDIF

        RETURN
        END
 
C++*********************************************************************
C
C  CADD.F                                   WRITTEN MAR 99 ARDEAN LEITH
C
C **********************************************************************
C
C       CADD IS KLUDGE TO AVOID EQUIVALENCING VOLBUF TO  CVOLBUFF
C
C **********************************************************************

        SUBROUTINE CADD(CVOLBUF,LUNIN,ITYPE,NX,NY,NZ,SIGT)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        COMPLEX         :: CVOLBUF(1)

        COMMON /IOBUF/ BUF(NBUFSIZ)
        COMPLEX        :: CBUF(1)
        EQUIVALENCE       (BUF(1),CBUF(1))

	NREC  = NY * NZ
        ILOC  = 0

	IF (SIGT == +2.0)  THEN
C          SIGT        +2    2 IS MULTIPLIED BY 1 --------------- MULT

C          FOURIER FILES
	   NSH = NX / 2

           DO IREC=1,NREC
	      CALL REDLIN(LUNIN,BUF,NX,IREC)
              DO  IX=1,NSH
                 ILOC          = ILOC + 1
                 CVOLBUF(ILOC) = CVOLBUF(ILOC) * CBUF(IX)
              ENDDO
           ENDDO

	ELSEIF (SIGT == -2.0)  THEN

C          COMPLEX 1 FOURIER MULTIPLICATION WITH 2 CONJUGATE ---- CMULT
           NSH = NX / 2

           DO IREC=1,NREC
              CALL REDLIN(LUNIN,BUF,NX,IREC)
              DO  IX=1,NSH
                 ILOC          = ILOC + 1
                 CVOLBUF(ILOC) = CVOLBUF(ILOC) * CONJG(CBUF(IX))
              ENDDO
           ENDDO

	ELSEIF (SIGT == +3.0)  THEN
C          SIGT        +3    1 IS SQUARED -------------------------- SQ

	   NSH = NX / 2

           DO IX=1,NREC*NSH
  	      CVOLBUF(IX) = CVOLBUF(IX) *CONJG(CVOLBUF(IX))
           ENDDO

	ELSEIF (SIGT == -3.0)  THEN
C          SIGT        -3    COMPLEX 1 IS DIVIDED BY COMPLEX 2 --- CDIV

	   NSH = NX / 2

           DO IREC=1,NREC
	      CALL REDLIN(LUNIN,BUF,NX,IREC)
              DO  IX=1,NSH
                 ILOC          = ILOC + 1
  	         CVOLBUF(ILOC) =  CVOLBUF(ILOC) / CBUF(IX)
              ENDDO
           ENDDO

        ENDIF

        END

