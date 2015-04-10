
C++*********************************************************************
C
C  IMSQ.F                         SPLIT FROM ADD  MAR 99 ARDEAN LEITH
C                                 USED REDVOL DEC 2000   ARDEAN LEITH
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
C  IMSQ(ROOT,LUNIN,LUNOUT,ITYPE,NSAM,NROW,NSLICE,IRTFLG)
C
C  PURPOSE:  SQUARE OR SQRT OF AN IMAGE
C
C  PARAMETERS:     
C        ROOT         FLAG FOR SQRT                             (SENT)
C        LUNIN        I/O UNIT NUMBER OF INPUT                  (SENT)
C        LUNOUT       I/O UNIT NUMBER OF OUTPUT                 (SENT)
C        ITYPET       IFORM OF INPUT VOLUME                     (SENT)
C        NSAM,NROW    X & Y DIMENSIONS OF IMAGES                (SENT)
C        NSLICE       Z DIMENSION OF IMAGES                     (SENT)
C        IRTFLG       ERROR FLAG                                (RET.)
C                          
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE IMSQ(ROOT,LUNIN,LUNOUT,ITYPE,NSAM,NROW,NSLICE,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        LOGICAL :: ROOT

        REAL,    ALLOCATABLE, DIMENSION(:) :: VOLBUF
        COMPLEX, ALLOCATABLE, DIMENSION(:) :: CVOLBUF

        COMMON /IOBUF/ BUF(NBUFSIZ)
        COMPLEX        CBUF(1)
        EQUIVALENCE    (BUF(1),CBUF(1))

        IRTFLG = 1

C       DOES NOT WORK ON SOME ODD FILE FORMATS NO LONGER IN USE
        IF (ITYPE.EQ.0 .OR. ITYPE.EQ.8 .OR. ITYPE.EQ.11 .OR.
     &	     ITYPE.EQ.12 .OR. ITYPE.EQ.16 .OR. ITYPE.EQ.-9)  THEN
            CALL ERRT(39,'IMSQ',NE)
            GOTO 9999
        ENDIF

C       LOAD VOLUME FROM FIRST FILE INTO VOLBUF 
        IF (ITYPE .LT. 0) THEN
           NSH = NSAM / 2
           ALLOCATE(CVOLBUF(NSH*NROW*NSLICE), STAT=IRTFLGT)
           IF (IRTFLGT .NE. 0) THEN
              CALL ERRT(46,'IMSQ,CVOLBUF',NDUM)
              GOTO 9999
           ENDIF
           CALL REDVOL(LUNIN,NSAM,NROW,1,NSLICE,CVOLBUF,IRTFLGT)
        ELSE
           ALLOCATE(VOLBUF(NSAM*NROW*NSLICE), STAT=IRTFLGT)
           IF (IRTFLGT .NE. 0) THEN
              CALL ERRT(46,'IMSQ,CVOLBUF',NDUM)
               GOTO 9999
           ENDIF
           CALL REDVOL(LUNIN,NSAM,NROW,1,NSLICE,VOLBUF,IRTFLGT)
        ENDIF
        IF (IRTFLGT .NE. 0) GOTO 9998

        IF (ROOT) THEN
C          WANT SQRT OF INPUT IMAGE ------------------------------ SQRT

	   IF  (ITYPE .GT. 0)  THEN
              DO  ISAM=1,NSAM * NROW * NSLICE
		 IF (VOLBUF(ISAM) .LT. 0.0) THEN
                    JSLICE = (ISAM / NSLICE) + 1

                    JROW   = (ISAM - (JSLICE - 1) * NSLICE) / NROW + 1

                    JSAM   = (ISAM - (JSLICE - 1) * NSLICE) -
     &                               (JROW   - 1) * NROW
                    WRITE(NOUT,90) JSLICE,JROW,JSAM,VOLBUF(ISAM)
90                  FORMAT(' *** AT: (',I6,',',I6,',',I6,') VALUE: ',
     &                     1PG11.4) 
		    CALL ERRT(101,'ATTEMPTED SQRT OF NEG. NUMBER',NE)
                    RETURN
		 ENDIF
  	         VOLBUF(ISAM) = SQRT(VOLBUF(ISAM)) 
              ENDDO

	   ELSE
C             NOT IMPLEMENTED FOR FOURIER
	      CALL ERRT(101,'NOT IMPLEMENTED FOR FOURIER FILES',NE)
              GOTO 9998
	   ENDIF

        ELSE
C          WANT TO SQUARE IMAGE ----------------------------------- SQ
	   IF (ITYPE .GT. 0)  THEN
C             REAL IMAGE

              DO  ISAM=1,NSAM * NROW * NSLICE
  		 VOLBUF(ISAM) = VOLBUF(ISAM) * VOLBUF(ISAM)
              ENDDO

           ELSEIF (ITYPE .LT. 0)  THEN	
C             FOURIER IMAGE
              DO ISAM=1,NSH * NROW * NSLICE
  	         CVOLBUF(ISAM) = CVOLBUF(ISAM) * CONJG(CVOLBUF(ISAM))
              ENDDO

           ENDIF
        ENDIF

C       WRITE TO OUTPUT FILE ON LUNOUT
        ILOC = 1
        DO IREC=1,NROW*NSLICE
           IF (ITYPE .LT. 0) THEN
              CALL WRTLIN(LUNOUT,CVOLBUF(ILOC),NSAM,IREC)
           ELSE
              CALL WRTLIN(LUNOUT,VOLBUF(ILOC),NSAM,IREC)
           ENDIF
           ILOC = ILOC + NSAM
        ENDDO
        IRTFLG = 0

9998    IF (ITYPE .LT. 0)  THEN
           DEALLOCATE(CVOLBUF)
        ELSE
           DEALLOCATE(VOLBUF)
        ENDIF

9999    CLOSE(LUNOUT)
        CLOSE(LUNIN)

        RETURN
        END

