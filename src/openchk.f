
C++*********************************************************************
C
C OPENCHK.F -- NEW NOVEMBER 1996  AUTHOR: ARDEAN LEITH
C              NX, NY, NZ ALL ON ONE LINE         AUG 2002 ARDEAN LEITH
C              'MS' MODIFICATIONS                 AUG 2002 ARDEAN LEITH
C              NX                                 DEC 2011 ARDEAN LEITH
C              NO MORE 'ENTER'                    JUN 2013 ARDEAN LEITH
C           
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
C    OPENCHK(NX,NY,NZ,ITYPE,NXMAX,IRTFLG)
C
C    PURPOSE:       CHECK ON INFO TO OPEN A NEW SPIDER IMAGE.
C
C    PARAMETERS:
C        NX,NY      DIMENSIONS OF FILE                     (SENT)
C        NZ         NUMBER OF PLANES                       (SENT)
C        ITYPE      IFORM                                  (SENT/RET.)                    
C        NXMAX      MAX. RECORD SIZE                       (SENT)                    
C        IRTFLG     ERROR RETURN FLAG.                     (RET.)
C                   IRTFLG = 0    NORMAL RETURN
C                   IRTFLG = 1    ERROR RETURN
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE OPENCHK(NX,NY,NZ,ITYPE,NXMAX,IRTFLG)

      INCLUDE 'CMBLOCK.INC'
	
C     CHECK THAT INPUT PARAMETERS ARE OK
C     IF NX OR NY IS NOT LOGICAL ASK FOR THEM

      IF ((ITYPE == 0  .OR. ITYPE == -1 ) .AND.  NZ <= 0) THEN
C        NEED TO INQUIRE DIMENSIONS (FROM 'MS')
         CALL RDPRI3S(NX,NY,NZ,NOT_USED,
     &              'X, Y, & Z DIMENSIONS',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         IF (ITYPE == -1 .AND. NZ <= 1)  THEN
C          FOURIER IMAGE STACK
           ITYPE  = -12 + MOD(NX,2)
           NZ = 1
           
         ELSEIF (ITYPE == -1)  THEN
C          FOURIER VOLUME STACK
           ITYPE = -22 + MOD(NX,2)

         ELSEIF (NZ <= 1) THEN
C           IMAGE STACK
            NZ = 1
            ITYPE  = 1
         ELSE
C           VOLUME STACK
            ITYPE  = 3
         ENDIF

      ELSEIF (ITYPE == 3 .AND. NZ <= 0) THEN
C        NEED TO INQUIRE DIMENSIONS FROM OTHER OPS THAN 'MS'
         CALL RDPRI3S(NX,NY,NZ,NOT_USED,
     &              'X, Y, & Z DIMENSIONS',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         IF (NZ .EQ. -9999) THEN
C           SET NZ TO 1
            NZ = 1

         ELSEIF (NZ <= 0) THEN
C           STILL NEED TO INQUIRE AS TO NZ
            CALL RDPRI1S(NZ,NOT_USED,
     &                 'Z DIMENSION',IRTFLG)
         ENDIF
         IF (IRTFLG .NE. 0) RETURN

         IF (NZ <= 1) ITYPE = 1

      ELSEIF (NX <= 0 .OR. NY <= 0) THEN
C        NEED TO FIND NX & NY FOR NEW FILE
         NY = 0
         CALL RDPRIS(NX,NY,NOT_USED,
     &               'X & Y DIMENSIONS',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
         IF (ITYPE == 0) ITYPE = 1
      ENDIF
      IRTFLG = 1
      IF (NY .LE. 0) NY = NX

      IF (NX <= 0 .OR. NY <= 0) THEN
C        FOR NEW FILES, NON-ZERO NX, NY NEED TO BE SUPPLIED
         CALL ERRT(101,'NX OR NY NOT GIVEN',NE)
         RETURN

      ELSEIF (NX > NXMAX) THEN
C        EXCESSIVE ROW LENGTH 
         CALL ERRT(102,'ROWLENGTH MUST BE <=:',NXMAX)
         RETURN

      ELSEIF (NZ == 0) THEN
C        NEW FILES NEED A NON ZERO SLICE NUMBER
         NZ = 1

      ELSEIF (NZ < 0) THEN
C        WANT TO CREATE OBSOLETE SHORT HEADER FILE
         WRITE(NOUT,*) 
     &            'WARNING:  SHORT HEADER FILES NO LONGER SUPPORTED.'
         NZ = -NZ
      ENDIF

      IF (ITYPE == -1 .OR. ITYPE == -3 .OR. ITYPE == -7) THEN
C        WANT TO CREATE OBSOLETE FOURIER FILE
         CALL ERRT(101,
     &         'CAN NOT CREATE OBSOLETE FOURIER FORMAT FILE',NE)
         RETURN

      ELSEIF (ITYPE == 0) THEN
C        WANT TO CREATE ???? TYPE
         CALL ERRT(101, 'PGM. ERROR, UNKNOWN FILE FORMAT',NE)
         RETURN

      ENDIF

      IRTFLG = 0

      END

