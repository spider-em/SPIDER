
C++*********************************************************************
C
C SETMODE_MRC              NEW FOR 'MD MRC'        NOV 19 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2019  Health Research Inc.,                         *
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
C  SETMODE_MRC
C
C  PURPOSE:  SET VARIOUS OPTIONAL MRC MODES (CURRENTLY ONLY ONE MODE)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE SETMODE_MRC()

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'
 
C     NUMBER OF OPERATIONS IN MRC MODE MENU
      INTEGER, PARAMETER    :: IMOFNC = 3
      CHARACTER(LEN=12)     :: MOMENU(IMOFNC)

      CHARACTER(LEN=12)     :: MODE
      CHARACTER(LEN=2)      :: CORIGIN
      CHARACTER(LEN=1)      :: CHAND
      CHARACTER(LEN=1)      :: NULL = CHAR(0)

      INTEGER               :: NLET,IRTFLG,IDUM

      DATA MOMENU/'ME      ',
     &            'ST      ',
     &            'AX      '/ 



C     MODE SWITCH OPERATION,  READ IN THE MODE
      CALL RDPRMC(MODE,NLET,.TRUE.,'MRC MODE',NULL,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      SELECT CASE(MODE(1:2))

      CASE("ME")
C       MENU ------------------------------------------------------- ME
        WRITE(NOUT,9610)
9610    FORMAT(/
     &  '  ME          ',T10, ' MRC MODE MENU    '/
     &  '  STA         ',T10, ' STATUS OF MODES  '/
     &  '  AX          ',T10, ' SET MRC DATA AXES')
        WRITE(NOUT,*) ' '

      CASE("ST")
C       DETERMINE STATUS ------------------------------------------ STA
        WRITE(NOUT,*) ' DATA ORIGIN:     ',MRC_AXIS(1:2)
        WRITE(NOUT,*) ' DATA HANDEDNESS: ',MRC_AXIS(4:4)

      CASE("AX")
C       SET MRC DATA AXIS ------------------------------------ SET AXIS
        CALL RDPRMC(CORIGIN,NLET,.TRUE.,  
     &              'DATA ORIGIN  (UL/LL)',NULL,IRTFLG)
        IF (IRTFLG .NE. 0)  RETURN
        IF (CORIGIN(1:2) .NE. 'UL' .AND. CORIGIN(1:2) .NE. 'LL')THEN
           CALL ERRT(101,'INVALID ORIGIN',IDUM)
           RETURN
        ENDIF

        CALL RDPRMC(CHAND,NLET,.TRUE.,
     &              'DATA HANDEDNESS (L/R)',NULL,IRTFLG)
        IF (IRTFLG .NE. 0)  RETURN
        IF (CHAND(1:1) .NE. 'L' .AND. CHAND(1:1) .NE. 'R')THEN
           CALL ERRT(101,'INVALID HANDEDNESS',IDUM)
           RETURN
        ENDIF
 
        MRC_AXIS(1:4) = CORIGIN(1:2) // ' ' // CHAND(1:1)

      CASE DEFAULT
        WRITE(NOUT,*) '  *** UNKNOWN MRC MODE'

      END SELECT

      END

