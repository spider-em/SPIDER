
C++*********************************************************************
C
C PUTPT1.F       NEW                               SEP 98 ARDEAN LEITH                 
C                RDPRAF REMOVED                    DEC 05 ARDEAN LEITH 
C                PROMPT                            MAY 13 ARDEAN LEITH

C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2013  Health Research Inc.,                         *
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
C PUTPT1(LUN,NX,NY,NZ)
C
C PURPOSE:  SUPERIMPOSE PIXELS ONTO AN IMAGE, 
C	    PIXEL LOCATIONS READ FROM  TERMINAL    
C            
C PARAMETERS: LUN	  UNIT NUMBER OF I/O FILE
C	      NX,NY,NZ    DIMENSIONS OF INPUT FILE
C
C--*********************************************************************
 
	SUBROUTINE PUTPT1(LUN,NX,NY,NZ)

	INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC'

        INTEGER, PARAMETER :: MAXNUM = 400
        INTEGER            :: IXCOOR(MAXNUM),IYCOOR(MAXNUM)
        INTEGER            :: IZCOOR(MAXNUM)
        REAL               :: FHEIGHT(MAXNUM)

        COMMON /IOBUF/BUF(NBUFSIZ)

        CHARACTER          :: NULL   = CHAR(0)
        DIMENSION          :: FLIST(4)
        LOGICAL            :: KEEPGO

        K      = 0
        KEEPGO = .TRUE.

        DO WHILE (KEEPGO)
           IF (IFORM .EQ. 3) THEN
C             VOLUME
              FLIST(4) = 1.0
              CALL RDPRA('X, Y, Z, & INTENSITY (0 TO HALT INPUT)',
     &            4,0,.FALSE.,FLIST,NVAL,IRTFLG)
              NVAL = 4

           ELSE
C             2-D IMAGE
              NVAL     = 3
              FLIST(3) = 1.0
              CALL RDPRM3S(FLIST(1),FLIST(2),FLIST(3),NOT_USED,
     &           'X, Y, & INTENSITY (0 TO HALT INPUT)',IRTFLG)
              IF (IRTFLG .NE. 0) RETURN
           ENDIF
           IF (IRTFLG .NE. 0) RETURN

           K          = K + 1
           IXCOOR(K)  = FLIST(1)
           IYCOOR(K)  = FLIST(2)
           IZCOOR(K)  = FLIST(3)
           FHEIGHT(K) = FLIST(NVAL)
           IF (NZ .EQ. 1) IZCOOR(K) = 1 

C          IS THIS END OF INPUT?
           IF (IXCOOR(K) <= 0 .OR. IYCOOR(K) <= 0) THEN
              KEEPGO = .FALSE.
              KGOT = K - 1

           ELSEIF (K >= MAXNUM) THEN
C             ARRAY OVERFLOW WILL OCCUR NEXT INPUT
              KEEPGO = .FALSE.
              WRITE(NOUT,*) '*** INPUT HALTED TO AVOID BUFFER OVERFLOW'
           ENDIF 
        ENDDO

C       ALL COORDINATES HAVE BEEN INPUT
        NUMSET = 0
	DO  I=1,KGOT
           IXCOR  = IXCOOR(I)
           IYCOR  = IYCOOR(I) 
           IZCOR  = IZCOOR(I)
           HEIGHT = FHEIGHT(I)

           IF ((IXCOR .GT. NX   .OR. IXCOR .LE. 0) .OR.
     &         (IYCOR .GT. NY   .OR. IYCOR .LE. 0) .OR.
     &         (IZCOR .GT. NZ .OR. IZCOR .LE. 0)) THEN
               WRITE(NOUT,721) IXCOR,IYCOR,IZCOR
721            FORMAT('  *** LOCATION: (',I5,',',I5,',',I5,
     &                ') OUTSIDE IMAGE, CONTINUING')
           ELSE
C             THIS COORDINATE IS OK, PUT IT IN FILE
              IREC = (IZCOR -1) * NY + IYCOR
              CALL REDLIN(LUN,BUF,NX,IREC)

              BUF(IXCOR) = HEIGHT
              CALL WRTLIN(LUN,BUF,NX,IREC)
              NUMSET = NUMSET + 1
           ENDIF
        ENDDO

9300    WRITE(NOUT,90) NUMSET
 90     FORMAT('  NUMBER OF LOCATIONS SET: ',I4/)

        IF (NUMSET > 0) THEN
C          SET FMIN, FMAX AS UNDETERMINED
           CALL SETPRM(LUN,NX,NY,0.0,0.0,0.0,'U')
        ENDIF
        
        CLOSE(LUN)

	END
