
C++*******************************************************************
C
C WINDOW.F
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
C   WINDOW(LUN1,LUN2,BUF,NSAM,NROW,NSLICE,NSAMW1,NROWW1,
C              NSLIW1,NSAM2,NROW2,NSLIC2,BACKG)
C
C   PURPOSE:    CUTS OUT A RECTANGULAR IMAGE SECTION
C
C       THIS SUBROUTINE CUTS OUT A RECTANGULAR IMAGE SECTION FROM
C       PICTURE 1 WITH A SPECIFIED SIZE NSAM2,NROW2 AT A SPECIFIED
C       LOCATION NSAMW1,NROWW1 AND WRITES OUT AN IMAGE.
C
C       WINDOW(LUN1,LUN2,BUF,NSAM,NROW,NSLICE,NSAMW1,NROWW1,
C              NSLIW1,NSAM2,NROW2,NSLIC2)
C         LUN1           	LOGICAL UNIT NUMBER OF INPUT IMAGE
C         LUN2           	LOGICAL UNIT NUMBER OF OUTPUT IMAGE
C         BUF            	BUFFER ARRAY OF SIZE NSAM + NSAM2
C         NSAM,NROW,NSLICE 	DIMENSIONS OF INPUT VOLUME/PICTURE
C         NSAMW1,NROWW1,NSLIW1  COORDINATES, WITH RESPECT TO INPUT VOLUME, 
C                        	OF TOP LEFT CORNER OF WINDOW
C         NSAM2,NROW2,NSLIC2    DIMENSIONS OF WINDOW = DIMENSIONS OF OUTPUT 
C                        	VOLUME
C
C NOTE: THIS APPEARS TO BE POORLY DESIGNED FOR EFFICIENCY, BUT I AM
C       RELUCTANT TO TAMPER WITH IT. agl NOV 95
C
C       FINALLY FIXED THIS ROUTINE. IT HAS BEEN BUGGY FOR YEARS
C       AND SLOW ALSO.  agl AUG 96
C
C--*******************************************************************

      SUBROUTINE WINDOW(LUN1,LUN2, NSAM,  NROW,  NSLICE, 
     &                             NSAMW1,NROWW1,NSLIW1, 
     &                             NSAM2, NROW2, NSLIC2, BACKG)

      INCLUDE 'CMBLOCK.INC'


      COMMON   BUF(1)
C     ASSUMES THAT COMMON BUFFER > NSAM + 2 * NSAMW1 

      IF (NSAMW1 .GT. NSAM .OR. NROWW1 .GT. NROW .OR. 
     &    NSLIW1 .GT. NSLICE) THEN
          WRITE(NOUT,*) '*** WINDOW OUTSIDE OF IMAGE.'
          CALL ERRT(100,'WINDOW',NE)
          RETURN
      ENDIF

C     FIND NUMBER OF BLANK SLICES AT TOP OF WINDOW
      ITOPS = 0
      IF (NSLIW1 .LT. 1) ITOPS = -NSLIW1 + 1
C     FIND NUMBER OF BLANK SLICES AT BOTTEM OF WINDOW
      IBOTS = 0
      NSLIW2 = NSLIW1 + NSLIC2 - 1
      IF (NSLIW2 .GT. NSLICE) IBOTS = NSLIW2 - NSLICE 

      NSLS   = 1
      NSLE   = 1
      IF (NSLICE .GT. 1) THEN
C        INPUT/OUTPUT FILE IS A VOLUME
         NSLS = MAX(1,NSLIW1)
         NSLE = MIN(NSLIW2,NSLICE)
      ENDIF

C     FIND NUMBER OF BLANK ROWS AT TOP OF WINDOW
      ITOP = 0
      IF (NROWW1 .LT. 1) ITOP = -NROWW1 + 1
C     FIND NUMBER OF BLANK ROWS AT BOTTEM OF WINDOW
      IBOT = 0
      NROWW2 = NROWW1 + NROW2 -1
      IF (NROWW2 .GT. NROW) IBOT = NROWW2 - NROW 

C     FIND FIRST ROW READ FROM IMAGE
      NS = MAX(NROWW1,1)
C     FIND LAST ROW READ FROM IMAGE
      NE = MIN(NROWW2,NROW)

C     CLEAR BUFFERS ONCE
      IF ((ITOP  .GT. 0 .OR. IBOT  .GT. 0) .OR.
     &    (ITOPS .GT. 0 .OR. IBOTS .GT. 0)) THEN
C        FILL BLANKING BUFFER FOR TOP / BOTTEM OVERFLOW
         DO  I=1,NSAM2
            BUF(I) = BACKG
         ENDDO
      ENDIF

C     SET BUFFER LOCATION FOR READ/WRITE OF IMAGE DATA
      KBUF1 = NSAM2 + 1
      IIN   = KBUF1
      IIOUT = KBUF1 + NSAMW1 - 1

      IF (NSAMW1 .LT. 1) THEN
C        WINDOW STARTS BEFORE SOURCE IMAGE ON LEFT
         IIN   = KBUF1 - NSAMW1 + 1
         IIOUT = KBUF1
      ENDIF
 
      IF (NSAMW1 .LT. 1 .OR. (NSAMW1 + NSAM2 - 1) .GT. NSAM) THEN
C        FILL OUTBUT BUFFER FOR LEFT / RIGHT OVERFLOW
         DO  I=KBUF1,KBUF1+NSAM2+NSAM+NSAM2
            BUF(I) = BACKG
         ENDDO
      ENDIF

      IRECT = 1
      IF (ITOPS .GT. 0) THEN
C        FILL BLANK SLICES BEFORE WINDOW
         DO  L = 1,ITOPS*NROW2
             CALL WRTLIN(LUN2,BUF,NSAM2,IRECT)
             IRECT = IRECT + 1
         ENDDO
      ENDIF                 
  
      DO  L = NSLS,NSLE
        IRECIN  = (L-1)*NROW

        IF (ITOP .GT. 0) THEN
C         WINDOW STARTS BEFORE IMAGE AT TOP, FILL ROWS WITH BACKG
          DO I = 1,ITOP
             CALL WRTLIN(LUN2,BUF,NSAM2,IRECT)
             IRECT = IRECT + 1
          ENDDO
        ENDIF

        DO I = NS,NE
           CALL REDLIN(LUN1,BUF(IIN),  NSAM, I+IRECIN)
           CALL WRTLIN(LUN2,BUF(IIOUT),NSAM2,IRECT)
           IRECT = IRECT + 1 
        ENDDO

        IF (IBOT .GT. 0) THEN
C         WINDOW GOES OFF IMAGE AT BOTTEM, FILL THESE WITH BACKG
          DO I = 1,IBOT
             CALL WRTLIN(LUN2,BUF,NSAM2,IRECT)
             IRECT = IRECT + 1
          ENDDO
        ENDIF
      ENDDO

      IF (IBOTS .GT. 0) THEN
C        FILL BLANK SLICES AFTER WINDOW
         DO  L = 1,IBOTS*NROW2
             CALL WRTLIN(LUN2,BUF,NSAM2,IRECT)
             IRECT = IRECT + 1
         ENDDO
      ENDIF  
               
      RETURN
      END
