
C++*********************************************************************
C
C GETFILENUM.F  -- NEW JAN 1999                   AUTHOR: ARDEAN LEITH
C                  EXTRACTED FROM LUNSETHDR        AUG 02 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors: Joachim Frank & ArDean Leith & Tapu Shaikh  *
C=* Copyright 1985-2010  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email:                                                             *
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
C    GETFILENUM(FILNAM,IMGNUM,NDIGITS,CALLERRT,IRTFLG)  
C
C    PURPOSE:    FINDS FILE NUMBER AT END OF FILENAME
C    
C    PARAMETERS: FILNAM    CHAR. VARIABLE FILE NAME             (SENT)
C                IMGNUM    NUMBER IN FILE NAME                   (RET)
C                NDIGITS   NUMBER OF DIGITS                      (RET)
C                CALLERRT  CALL ERRT IF ERROR                   (SENT)
C                IRTFLG    ERROR FLAG                            (RET)
C
C    CALLED BY:  deletf.f   filgen.f    openinstk.f 
C                openstk.f  to_peaks.f 
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************
 
       SUBROUTINE GETFILENUM(FILNAM,IMGNUM,NDIGITS,CALLERRT,IRTFLG)

       IMPLICIT  NONE

       INCLUDE 'CMBLOCK.INC'

       CHARACTER *(*) :: FILNAM
       LOGICAL        :: CALLERRT
       INTEGER        :: IMGNUM,NDIGITS,IRTFLG

       CHARACTER *1   :: CHARI
       INTEGER        :: NLET,IGO,I,NE
       INTEGER        :: lnblnkn


C      FIND NUMBER OF CHAR. IN FILNAM
       NLET   = LNBLNKN(FILNAM)

       IGO    = NLET - 1
       IMGNUM = 0

C      EXTRACT IMGNUM FROM FILENAME
       DO I = NLET,1,-1
          CHARI = FILNAM(I:I)
          IF (CHARI .LT. '0' .OR. CHARI .GT. '9') EXIT
          IGO = I 
       ENDDO

       NDIGITS = NLET - IGO + 1
       IF (NDIGITS .LE. 0 .OR. NDIGITS .GT. 10) THEN
C         NO NUMBER AT END OF FILNAM OR > 10 DIGITS
          IRTFLG = -1
          RETURN
       ENDIF

       READ(FILNAM(IGO:NLET),'(I10)',ERR=999) IMGNUM
       IRTFLG = 0
       RETURN
 
    
999    WRITE(NOUT,*) '*** CAN NOT GET FILE NUMBER FROM: ',FILNAM(:NLET)
       IF (CALLERRT) THEN
          CALL ERRT(100,'GETFILENUM',NE)
       ENDIF
       IRTFLG = 1

       RETURN
       END




#ifdef NEVER
      PURE LOGICAL FUNCTION LITTLE_ENDIAN()
      
      INTEGER(INT8)  :: J(2)
      INTEGER(INT16) :: I
      
      I = 1
      J = TRANSFER(SOURCE=I,MOLD=J,SIZE=2)
      
      IF (J(1) == 1) THEN
         LITTLE_ENDIAN = .TRUE.
      ELSE
         LITTLE_ENDIAN = .FALSE.
      ENDIF
      
      END FUNCTION
      
      
      SUBROUTINE ENDIAN(LITEND)

C     CHECKS IF THIS IS A LITTLE ENDIAN MACHINE
C     RETURNS LITEND =.TRUE. IF IT IS, LITEND =.FALSE. IF NOT

      INTEGER*1 J(2)
      INTEGER*2 I
      EQUIVALENCE (I,J)
      
      LOGICAL LITEND

      I = 1
      IF (J(1) .EQ. 1) THEN
         LITEND = .TRUE.
      ELSE
         LITEND = .FALSE.
      ENDIF

      END
      
      PROGRAM TESTI
      ! Posted by Perseus in comp.lang.fortran on 4 July 2005.
      ! and Paul Van Delst and David Flower on 5 July 2005.

      LOGICAL, PARAMETER :: BIGEND = IACHAR(TRANSFER(1,"A")) == 0

      IF (BIGEND) THEN
        PRINT *, "BIG ENDIAN"
      ELSE
        PRINT *, "LITTLE ENDIAN"
      ENDIF

END PROGRAM TESTI

PROGRAM DETERMINE_ENDIANNESS
    IMPLICIT NONE
    INTEGER :: UNIT, IOSTAT
    INTEGER(4) :: VALUE
    CHARACTER(LEN=4) :: BYTES

    ! OPEN THE FILE IN UNFORMATTED MODE
    OPEN(UNIT=10, FILE="FILE.DAT", FORM="UNFORMATTED", ACCESS="STREAM", STATUS="OLD", IOSTAT=IOSTAT)
    IF (IOSTAT /= 0) THEN
        PRINT *, "ERROR OPENING FILE."
        STOP
    END IF

    ! READ THE FIRST 4 BYTES (ASSUMING THE FILE STARTS WITH A KNOWN INTEGER)
    READ(10) VALUE
    CLOSE(10)

    ! CONVERT THE INTEGER TO BYTES FOR ANALYSIS
    CALL MOVE_ALLOC(TRANSFER(VALUE, BYTES), BYTES)

    ! CHECK THE BYTE ORDER
    IF (BYTES == 'EXPECTED_BIG_ENDIAN_BYTES') THEN
        PRINT *, "FILE IS BIG ENDIAN."
    ELSE IF (BYTES == 'EXPECTED_LITTLE_ENDIAN_BYTES') THEN
        PRINT *, "FILE IS LITTLE ENDIAN."
    ELSE
        PRINT *, "UNABLE TO DETERMINE ENDIANNESS."
    END IF
END PROGRAM DETERMINE_ENDIANNESS



#endif


