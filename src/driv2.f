
C++*********************************************************************
C
C    DRIV2                USED OPFILEC              MAR 03 ARDEAN LEITH
C                         ADDED WARP3               MAR 04 ARDEAN LEITH
C                         ADDED TOSPIRE             SEP 05 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2015  Health Research Inc.,                         *
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
C   DRIV2(MAXDIM)
C
C   PURPOSE:   CONTAINS ROUTINES REMOVED FROM DRIVER IN MAR 93
C              ALSO CONTAINS SOME ROUTINES RECENTLY ADDED
C
C   PARAMETERS: MAXDIM     MAX LENGTH OF COMMON BUFFER
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE DRIV2(MAXDIM)

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 

        CHARACTER (LEN=MAXNAM)   ::  FILNAM
        CHARACTER (LEN=2*MAXNAM) ::  CSTRING
        CHARACTER (LEN=1)        ::  CDUM

C       DATA MENU/'WA','WT','SE','TS'/

	DATA  LUNIN,LUNOUT,LUNDOC/21,22,71/

        MAXIM  = 0
        IRTFLG = 1

C       SELECT THE OPERATION
        SELECT CASE(FCHAR(1:2))

        CASE('WA') 
C          WARPING ------------------------------------------------ WA
           CALL WARP3(LUNDOC,LUNIN,LUNOUT)

        CASE('WT') 
C          LATTICE REFLECTION PICKING --------------------------- WT TV
C          MANUAL PAGE MISSING 2015
           CALL TVWN3(MAXDIM)

        CASE('SE') 
C          SEED FILL ---------------------------------------------- SE 
C          OPEN INPUT FILE, NO FOURIER INPUT ALLOWED 
201	   CALL OPFILEC(0,.TRUE.,FILNAM,LUNIN,'O',IFORM,
     &               NSAM1,NROW1,NSLICE1,MAXIM,'INPUT',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 5000

           IF (NSLICE1 .GT. 1) THEN
              CALL ERRT(101,'CAN NOT USE 3-D FILES',NE)
              GOTO 5000

           ELSEIF ((2 * NSAM1 * NROW1 + 4) .GE. MAXDIM) THEN
              CALL ERRT(9,'DRIV2 ',NE)
              GOTO 5000
           ENDIF
           FMIN1 = FMIN
           IF (IMAMI .NE. 1) 
     &      CALL NORM3(LUNIN,NSAM1,NROW1,NSLICE1,FMAX1,FMIN1,AVR1)

C          OPEN THE OUTPUT FILE
202        IFORM = 1
           MAXIM = 0
           CALL OPFILEC(LUNIN,.TRUE.,FILNAM,LUNOUT,'U',IFORM,
     &              NSAM1,NROW1,NSLICE1,MAXIM,'OUTPUT',.FALSE.,IRTFLG)
           IF (IRTFLG .EQ. -1) THEN
              CLOSE(LUNIN)
              GOTO 201
           ENDIF
           IF (IRTFLG .NE. 0) GOTO 5000
 
           CALL SEEDFILL(LUNIN,LUNOUT,NSAM1,NROW1,FMIN1,MAXDIM,IRTFLG)
           CLOSE(LUNOUT)
           IF (IRTFLG .EQ. -1)  GOTO 202
           CLOSE(LUNIN)

        CASE('TS') 
C          SPIRE OUTPUT -------------------------------------------- TS 

           IRTFLG = -999    ! NO UPPERCASING
           IF (FCHAR(4:4) .EQ. 'F') THEN
              CALL RDPRMC(CSTRING,NCHAR,.FALSE.,
     &                    'SPIRE OUTPUT FILE STRING',CDUM,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 5000

C             PDATES ADDS DATE AND CALLS SPIREOUT
              CSTRING = CSTRING(1:NCHAR) // '  WRITTEN  '
              CALL PDATES(CSTRING,-1)

           ELSE
              CALL RDPRMC(CSTRING,NCHAR,.FALSE.,'SPIRE OUTPUT STRING',
     &                    CDUM,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 5000
              CALL SPIREOUT(CSTRING,IRTFLG)
           ENDIF

        END SELECT

5000    RETURN
        END     


