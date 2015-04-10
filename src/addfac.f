
C++*********************************************************************
C
C  ADDFAC.F                 NEW                    MAR 03  ARDEAN LEITH
C                           NAN TRAP               OCT 09  ARDEAN LEITH
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
C  ADDFAC(VOLBUF,LUNIN,IFORMT,NSAM,NROW,NSLICE,SIGN,FACT1,FACT2)
C
C  PURPOSE:  ADD, SUBTRACT, OR MULTIPLY IMAGES WITH FACTORS
C            BOTH IMAGES HAVE TO HAVE THE SAME SIZE.
C
C  PARAMETERS:     
C        VOLBUF       INPUT VOLUME (#1)
C        LUNIN        I/O UNIT NUMBER OF INPUT FILE #2
C        IFORMT       IFORM OF INPUT VOLUME
C        NSAM,NROW    X & Y DIMENSIONS OF IMAGES
C        NSLICE       Z DIMENSION OF IMAGES
C        SIGN        +1000  1 IS ADDED/SUBTRACTED TO 2
C                    +2000  1 IS MULTIPLIED WITH 2 with ratio
C                          
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE ADDFAC(VOLBUF,LUNIN,IFORMT,NSAM,NROW,NSLICE,SIGN,
     &                 FACT1,FACT2)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        COMMON /IOBUF/ BUF(NBUFSIZ)

        REAL,DIMENSION(*) :: VOLBUF
        INTEGER *8        :: ILOC

	NREC  = NROW * NSLICE
        ILOC  = 0

	IF (SIGN .EQ. 1000.0)  THEN
C          a TIMES FACT1 IS ADDED TO b * FACT2

           WRITE(NOUT,*) ' DENOTING FILES AS a, b, & c'
           WRITE(NOUT,90) FACT1,FACT2
90         FORMAT('  c = ',1PG12.4,' * a + ', 1PG12.4,' * b',/)

           DO IREC=1,NREC
              CALL REDLIN(LUNIN,BUF,NSAM,IREC)
              DO  ISAM=1,NSAM
                  ILOC         = ILOC + 1
                  VOLBUF(ILOC) = VOLBUF(ILOC) * FACT1 + 
     &                           BUF(ISAM)    * FACT2
                ENDDO
           ENDDO
           	
	ELSE
C          RATIO OF a TO b TIMES FACT2 IS ADDED TO a * FACT1

           WRITE(NOUT,*) ' DENOTING FILES AS a, b, & c'
           WRITE(NOUT,91) FACT1,FACT2
91         FORMAT('  c = ',1PG12.4,' * a + ', 1PG12.4,' * (a/b)',/)

           DO IREC=1,NREC
              CALL REDLIN(LUNIN,BUF,NSAM,IREC)
              DO  ISAM=1,NSAM
                ILOC         = ILOC + 1
                IF (BUF(ISAM) .EQ. 0.0) THEN
                    CALL ERRT(102,
     &                      'CAN NOT DIVIDE BY ZERO IN RECORD',IREC)
                    RETURN
                ELSE 
                   VOLBUF(ILOC) = VOLBUF(ILOC) * FACT1 + 
     &                           (VOLBUF(ILOC) / BUF(ISAM)) * FACT2
                ENDIF
              ENDDO
            ENDDO
        ENDIF

        RETURN
        END
 
