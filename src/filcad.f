C++*********************************************************************
C
C    FILCAD.FOR -- CREATED JAN 23 1989
C                  REPLACES FILCOD.FOR     USES LONG FILE NAMES
C
C **********************************************************************
C    AUTHOR:  ARDEAN LEITH
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
C     FILCAD(FILNAM,FILPAT,NUMRET,IRTFLG)
C
C     PURPOSE:   CONVERTS ALL DIGITS AT END OF FILE NAME INTO INTEGER -- 
C                E.G., PIC064 WILL RESULT IN:  NUMRET = 64.
C                AND ALSO RETURNS A PATTERN FOR THE FILE NAME.
C                OFTEN USED WITH FILGET(FILPAT,FILNAM,NLET,INUM,IRTFLG)
C                TO CREATE A SERIES OF FILES
C
C     PARAMETERS:
C        FILNAM    CHAR. VARIABLE HOLDING INPUT FILE NAME   (INPUT)
C        FILPAT    CHAR. VARIABLE HOLDING FILE NAME PATTERN (RETURNED)
C        NUMRET    INTEGER VARIABLE TO RETURN NUMBER        (RETURNED)
C        IRTFLG    ERROR FLAG (0 IS NORMAL)                 (RETURNED)
C
C--*******************************************************************

        SUBROUTINE FILCAD(FILNAM,FILPAT,NUMRET,IRTFLG)
 
        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'

        CHARACTER (LEN=*) :: FILNAM,FILPAT
        INTEGER           :: NUMRET,IRTFLG

        CHARACTER (LEN=1) :: CHARI
        CHARACTER (LEN=1) :: NULL = CHAR(0)
        CHARACTER (LEN=9) :: DIGITS
        CHARACTER (LEN=9) :: ASTS = '*********'
        INTEGER           :: NLET,I,IGO,NDIGITS,IERR

        NLET = INDEX(FILNAM,NULL) - 1

C       GO BACKWARDS THRU FILENAME TO LOCATE DIGITS
        DO I=NLET,1,-1
          CHARI = FILNAM(I:I)
          IF (CHARI < '0' .OR. CHARI > '9' ) EXIT
        ENDDO

        IGO     = I + 1
        NDIGITS = NLET - IGO + 1
        IRTFLG  = 1

        IF (NDIGITS > 0 .OR. NDIGITS < 9) THEN

           READ(FILNAM(IGO:NLET),'(I9)',IOSTAT=IERR) NUMRET
           IF (IERR .NE. 0) THEN
              WRITE(NOUT,*) ' *** ILLEGAL FILE NUMBER'
              RETURN
           ENDIF

           IF (LEN(FILPAT) > (IGO+NDIGITS)) THEN
              FILPAT  = FILNAM(1:IGO-1) // ASTS(1:NDIGITS) // ' '
              IRTFLG  = 0
           ELSE
              WRITE(NOUT,*)' *** PGM. ERROR, FILPAT TOO SHORT IN FILCAD'
           ENDIF
        ENDIF           
        
        END
