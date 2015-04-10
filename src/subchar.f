
C++*********************************************************************
C
C  SUBCHAR.F -- CREATED                            SEP 00 ARDEAN LEITH 
C               ADDED SUBCHARS                     NOV 09 ARDEAN LEITH
C               SUBCHARS REMOVE BUG                MAR 14 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C    SUBCHAR(INSERT,ORIGINAL,LOC1,LOC2,LENUSED,IRTFLG)
C
C    PURPOSE:       STRING SUBSTITUTION
C
C    PARAMETERS:    INSERT     TO BE INSERTED                    (SENT)
C                   ORIGINAL   ORIGINAL STRING              (SENT/RET.)
C                   LOC1       START OF INSERT                   (SENT)
C                   LOC2       END OF INSERT                     (SENT)
C                   LENUSED    TOTAL LEN OF NEW STRING           (RET.)
C                   IRTFLG     RETURN FLAG (0 IS NORMAL)         (RET.)
C   
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

       SUBROUTINE SUBCHAR(INSERT,ORIGINAL,LOC1,LOC2,LENUSED,IRTFLG)

       CHARACTER *(*) INSERT, ORIGINAL

       IRTFLG   = 1

       LENI     = LEN(INSERT)
       LENO     = LEN(ORIGINAL)
       LENAFTER = 0
       IF (LOC2 < LENO) LENAFTER = LNBLNKN(ORIGINAL(LOC2+1:))

       LENUSED  = LOC1 - 1 + LENI + LENAFTER

       IF (LENUSED > LENO) THEN
C         OVERFLOW
          CALL ERRT(102,'STRING HAS TOO MANY CHARACTERS',LENUSED)
          IRTFLG = 1
          RETURN
       ENDIF
 
C      MOVE SUFFIX AFTER INSERT SO THAT STRING HOLDS INSERT 
       ORIGINAL(LOC1+LENI:) = ORIGINAL(LOC2+1:)

C      PLUG IN INSERT AT LOC1 BEFORE SUFFIX
       IF ((LOC1+LENI-1) >= LOC1) 
     &        ORIGINAL(LOC1:LOC1+LENI-1)  = INSERT

       IRTFLG = 0

       END


C --------------------------- SUBCHARS --------------------------------


C ************************** SUBCHARS *********************************
C
C    SUBCHARS(ORIGINAL,REMOVE,INSERT,LENOUT,IRTFLG)
C
C    PURPOSE:       STRING SUBSTITUTION
C
C    PARAMETERS:    ORIGINAL   ORIGINAL STRING              (SENT/RET.)
C                   REMOVE     START OF INSERT                   (SENT)
C                   INSERT     STRING TO BE INSERTED             (SENT)                   (SENT)
C                   LENOUT     TOTAL LEN OF NEW STRING           (RET.)
C                   IRTFLG     RETURN FLAG (0 IS NORMAL)          (RET.)
C   
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************
 
      SUBROUTINE SUBCHARS(ORIGINAL,REMOVE,INSERT,LENOUT,IRTFLG)

      CHARACTER(LEN=*) ::  ORIGINAL,REMOVE,INSERT

      LENOUT    = lnblnkn(ORIGINAL)
      LENR      = LEN(REMOVE)
      LENINSERT = LEN(INSERT)

      DO 
         ILOC = INDEX(ORIGINAL(1:LENOUT),REMOVE)
         IF (ILOC <= 0) EXIT

         IF (LENINSERT <= 0) THEN
C           SIMPLE REMOVAL OF REMOVE STRING
            ORIGINAL(ILOC:) = ORIGINAL(ILOC+LENR:)
            LUNOUT = LENOUT - LENR
         ELSE
C           SUBSTITUTE INSERT FOR REMOVE
            CALL SUBCHAR(INSERT,ORIGINAL,ILOC,ILOC+LENR-1,LENOUT,IRTFLG)
            IF (IRTFLG .NE. 0) RETURN
         ENDIF
      ENDDO

      IRTFLG = 0
      END

