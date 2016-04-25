C++*********************************************************************
C
C    GETNEXTTOKEN.F   ACCEPTS TAB                  NOV 99 ARDEAN LEITH
C                     GETNEXTTOKEN_N               MAY 07 ARDEAN LEITH
C                     GETNEXTTOKEN2                JUN 09 ARDEAN LEITH
C                     INQUOTE1 BUG                 AUG 09 ARDEAN LEITH
C                     , INQUOTE1 BUG 9.TOKEN2)     NOV 09 ARDEAN LEITH
C                     ! COMMENT DELIMITER          DEC 11 ARDEAN LEITH
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2011  Health Research Inc.,                         *
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
C    GETNEXTTOKEN(STRING,NCHAR,IFIRST,IGO,IEND)
C
C    PURPOSE:           GETS NEXT TOKEN DELIMITED BY A COMMA, BLANK, 
C                       SEMICOLON ETC, FROM STRING(IGO:IEND)
C                       IGO...IEND IS RETURNED TOKEN POSITION, IF
C                       NO MORE TOKENS LEFT IN STRING, IGO IS ZERO.
C
C    PARAMETERS:
C	  STRING	INPUT STRING                             (SENT)
C	  IFIRST        FIRST CHAR TO BE SCANNED FOR TOKEN       (SENT)
C	  IGO		START OF TOKEN                           (RET.)
C	  IEND		END OF TOKEN                             (RET.)
C	  NCHAR		LAST CHAR IN STRING                      (SENT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

C     ---------------------  GETNEXTTOKEN_NOT --------------------

      SUBROUTINE GETNEXTTOKEN_NOT(STRING,VALID,IFIRST,IGO,IEND)

C     VALID TOKEN CONTENTS GIVEN IN: VALID

      IMPLICIT NONE
      CHARACTER(LEN=*) :: STRING
      CHARACTER(LEN=*) :: VALID
      INTEGER          :: IFIRST,IGO,IEND

      CHARACTER        :: CTEMP
      INTEGER          :: ILAST,I,ID

C     SET DEFAULT RETURN VALUES
      IGO   = 0
      IEND  = 0

C     FIND LAST CHAR POSITION IN STRING
      ILAST = LEN(STRING)

      DO I = IFIRST,ILAST
         CTEMP = STRING(I:I)

         ID = INDEX(VALID,CTEMP)   ! >0 is within VALID

         IF (ID == 0 .AND. IGO > 0) THEN
C           THIS IS TOKEN END
            EXIT

         ELSEIF (ID > 0 .AND. CTEMP > CHAR(31)) THEN
C           ARE INSIDE A TOKEN, SET IGO IF THIS IS FIRST CHAR INSIDE
            IF (IGO == 0) IGO = I
            IEND = I
            !write(6,*) ' Token start:',id,ctemp

         ENDIF
      ENDDO
      !write(6,*) ' Token range:',igo,iend,string(igo:iend)

C     TOKEN FOUND OR RAN OFF END OF STRING

      END

C     -------------------------  GETNEXTTOKEN_D --------------------

      SUBROUTINE GETNEXTTOKEN_D(STRING,DELIM,IFIRST,IGO,IEND)

C     TOKEN DELIMITER IS GIVEN IN: DELIM

      IMPLICIT NONE
      CHARACTER(LEN=*) :: STRING
      CHARACTER(LEN=*) :: DELIM
      INTEGER          :: IFIRST,IGO,IEND

      CHARACTER        :: CTEMP
      INTEGER          :: ILAST,I,ID

C     SET DEFAULT RETURN VALUES
      IGO   = 0
      IEND  = 0

C     FIND LAST CHAR POSITION IN STRING
      ILAST = LEN(STRING)

      DO I = IFIRST,ILAST
         CTEMP = STRING(I:I)

         ID = INDEX(DELIM,CTEMP)

         IF (ID > 0 .AND. IGO > 0) THEN
C           THIS IS TOKEN END
             EXIT

         ELSEIF (ID == 0 .AND. CTEMP > CHAR(31)) THEN
C           ARE INSIDE A TOKEN, SET IGO IF THIS IS FIRST CHAR INSIDE
            IF (IGO == 0) IGO = I
            IEND = I
            !write(6,*) ' token start:',id,ctemp

         ENDIF
      ENDDO

C     TOKEN FOUND OR RAN OFF END OF STRING

      END

C      ********************* GETNEXTTOKEN_N **************************


      SUBROUTINE GETNEXTTOKEN_N(STRING,NCHAR,IFIRST,IGO,IEND)

C     TOKEN DELIMITER IS A COMMA, BLANK, ;, !, ], OR NEG. SIGN

      CHARACTER *(*) STRING
      CHARACTER      CTEMP

C     SET DEFAULT RETURN VALUES
      IGO   = 0
      IEND  = 0

      INVAR = 0

      DO I = IFIRST,NCHAR
         CTEMP = STRING(I:I)
         IF (CTEMP .EQ. '[') INVAR = INVAR + 1
         IF (CTEMP .EQ. ']') INVAR = INVAR - 1

         IF (CTEMP == ';' .OR. CTEMP == '!') THEN
C           START OF COMMENT, THIS IS TOKEN END
            RETURN

         ELSEIF (CTEMP == ']') THEN
C           TOKEN END ENCOUNTERED
            IF (IGO .EQ. 0) IGO = I
            IEND = I
            RETURN

         ELSEIF (CTEMP == '-' .AND. IGO == 0) THEN
C           TOKEN END ENCOUNTERED,  JUST  A SINGLE INITIAL -
            IGO  = I
            IEND = I
            RETURN

         ELSEIF (CTEMP == '-' .AND. INVAR .LE. 0) THEN
C           TOKEN END ENCOUNTERED,  IS A -
            RETURN

         ELSEIF (CTEMP .NE. ',' .AND. CTEMP .NE. ' ' .AND.
     &           CTEMP .NE. '(' .AND. CTEMP .NE. ')' .AND.
     &           CTEMP .GT. CHAR(31)) THEN
C           ARE INSIDE A TOKEN, SET IGO IF THIS IS FIRST CHAR INSIDE
            IF (IGO .EQ. 0) IGO = I
            IEND = I

         ELSEIF (IGO > 0) THEN
C           TOKEN END ENCOUNTERED
            RETURN
         ENDIF
      ENDDO

C     TOKEN FOUND OR RAN OFF END OF STRING

      RETURN
      END


C     -------------------------  GETNEXTTOKEN --------------------

      SUBROUTINE GETNEXTTOKEN(STRING,IFIRST,IGO,IEND)

C     TOKEN DELIMITER IS A COMMA, BLANK, !, OR SEMICOLON

      CHARACTER *(*) STRING
      CHARACTER      CTEMP

C     SET DEFAULT RETURN VALUES
      IGO   = 0
      IEND  = 0

C     FIND LAST CHAR POSITION IN STRING
      ILAST = LEN(STRING)

      DO I = IFIRST,ILAST
         CTEMP = STRING(I:I)

         IF (CTEMP == ';' .OR. CTEMP == '!') THEN
C           START OF COMMENT, THIS IS TOKEN END
            EXIT

         ELSEIF (CTEMP == ']') THEN
C           TOKEN END ENCOUNTERED
            IF (IGO .EQ. 0) IGO = I
            IEND = I
            RETURN

         ELSEIF (CTEMP .NE. ',' .AND. CTEMP .NE. ' ' .AND.
     &           CTEMP .NE. '(' .AND. CTEMP .NE. ')' .AND.
     &           CTEMP .GT. CHAR(31)) THEN
C           ARE INSIDE A TOKEN, SET IGO IF THIS IS FIRST CHAR INSIDE
            IF (IGO .EQ. 0) IGO = I
            IEND = I

         ELSEIF (IGO > 0) THEN
C           TOKEN END ENCOUNTERED
            RETURN
         ENDIF
      ENDDO

C     TOKEN FOUND OR RAN OFF END OF STRING

      RETURN
      END


C     -------------------------  GETNEXTTOKEN2 --------------------

      SUBROUTINE GETNEXTTOKEN2(STRING,IFIRST,IGO,IEND)

C     TOKEN END DELIMITER IS A COMMA, BLANK, ), ], ", ', !, OR ;

      CHARACTER(LEN=*)    :: STRING
      CHARACTER(LEN=1)    :: CTEMP,NQ1,NQ2
      LOGICAL             :: INQUOTE1,INQUOTE2

      NQ1   = CHAR(39)   ! '
      NQ2   = CHAR(34)   ! "

C     SET DEFAULT RETURN VALUES
      IGO      = 0
      IEND     = 0

      INQUOTE1 = .FALSE.
      INQUOTE2 = .FALSE.

C     FIND LAST CHAR POSITION IN STRING
      ILAST = LEN(STRING)

      DO I = IFIRST,ILAST
         CTEMP = STRING(I:I)
         
         IF (CTEMP == '=' .OR.
     &      (CTEMP == ']' .AND. 
     &      (.NOT. INQUOTE1 .AND. .NOT. INQUOTE2)) .OR.  
     &      (INQUOTE1 .AND. CTEMP == NQ1)   .OR. 
     &      (INQUOTE2 .AND. CTEMP == NQ2)) THEN
C           TOKEN END (" ') ENCOUNTERED, INCLUDE THIS SYMBOL IN TOKEN
C           TOKEN END ENCOUNTERED, INCLUDE THIS SYMBOL IN TOKEN
            IF (IGO .EQ. 0) IGO = I
            IEND = I
            EXIT

         ELSEIF ((CTEMP == ' ' .OR. CTEMP == ',' .OR. 
     &            CTEMP == '(' .OR. CTEMP == ')' ) .AND. 
     &          (INQUOTE1 .OR. INQUOTE2)) THEN
C           BLANK, COMMA,(, OR ) INSIDE QUOTED STRING IS NOT END OF TOKEN
            CYCLE

         ELSEIF (CTEMP .NE. ',' .AND. CTEMP .NE. ' ' .AND.
     &           CTEMP .NE. '(' .AND. CTEMP .NE. ')' .AND.
     &           CTEMP .NE. ';' .AND. CTEMP .NE. '!' .AND. 
     &           CTEMP .GT. CHAR(31)) THEN
c           NOT A: COMMA, BLANK, SEMI, !, OR PARENTHESIS
C           ARE INSIDE A TOKEN, SET IGO IF THIS IS FIRST CHAR INSIDE
            IF (IGO .EQ. 0) IGO = I
            IEND = I
 
C           SET INQUOTE FOR THIS TYPE OF QUOTE DELIMITER
            IF (CTEMP .EQ. NQ1 .AND. 
     &          .NOT. INQUOTE1 .AND. .NOT. INQUOTE2) INQUOTE1 = .TRUE. 
            IF (CTEMP .EQ. NQ2 .AND. 
     &          .NOT. INQUOTE1 .AND. .NOT. INQUOTE2) INQUOTE2 = .TRUE. 
            IF (.NOT. INQUOTE1 .AND. .NOT. INQUOTE2 .AND. 
     &          CTEMP .EQ. NQ1) INQUOTE1 = .TRUE.

         ELSEIF (IGO > 0) THEN
C           TOKEN HAS ENDED, DO NOT INCLUDE CHAR IN TOKEN
            EXIT
         ENDIF
      ENDDO

C     TOKEN FOUND OR RAN OFF END OF STRING
      !if (igo > 0) write(6,*) ' Token:',string(igo:iend),':',IGO,IEND

      END


C     -------------------------  GET_TOKENS --------------------


       SUBROUTINE GET_TOKENS(RECLIN,NLET,MAXLENTOK,NTOKSMAX,
     &                      INSIDE,     DELIM, TOKENS, NGOT,   IRTFLG)

C      FILLS TOKENS ARRAY WITH CONSECUTIVE TOKENS FROM A LINE

       IMPLICIT NONE
       INCLUDE 'CMBLOCK.INC'

       CHARACTER(LEN=NLET)      :: RECLIN          ! TOKEN LINE  (SENT)
       INTEGER                  :: NLET            ! LINE LENGTH (SENT)
       INTEGER                  :: MAXLENTOK,NTOKSMAX  ! (SENT)
       CHARACTER(LEN=MAXLENTOK) :: TOKENS(NTOKSMAX)! TOKEN FIELDS
       LOGICAL                  :: INSIDE          ! DELIM IS INSIDE (SENT)
       CHARACTER(LEN=*)         :: DELIM           ! TOKEN DELIMITERS  (SENT)
       INTEGER                  :: NGOT            ! NO. OF TOKENS (RETURNED
       INTEGER                  :: IRTFLG          ! ERROR FLAG     (RETURNED)

       INTEGER                  :: IGO,IEND,IFIRST,NCHAR         
       INTEGER                  :: I

       NGOT   = 0
       IRTFLG = 1
       IEND   = 0

       DO      ! ----------------- LOOP --------------------------
         IFIRST = IEND + 1

         IF (INSIDE) THEN
C           GET TOKEN (CHAR. STRING CONTAINS SYMBOLS IN: DELIM)
            CALL GETNEXTTOKEN_NOT(RECLIN(1:NLET),
     &                            DELIM,IFIRST,IGO,IEND)
         ELSE

C          GET TOKEN (CHAR. STRING DELIMITED BY SYMBOLS IN: DELIM)
            CALL GETNEXTTOKEN_D(RECLIN(1:NLET),
     &                            DELIM,IFIRST,IGO,IEND)
         ENDIF

C        SEE IF ALL TOKENS FROM RECLIN STRING HAVE BEEN EVALUATED
         IF (IGO <= 0) EXIT

C        TOKEN RETURNED, SET TOKEN IN TOKENS
         NGOT = NGOT + 1
         IF (NGOT > NTOKSMAX) THEN
            CALL ERRT(102,'MAX NO. OF TOKENS',NTOKSMAX)
            RETURN
         ENDIF

C        TOKEN RETURNED
         !write(6,'(A,i0,a,a)')'  Token(',ngot,'):',reclin(igo:iend)
         !write(6,*) '  igo...iend: ',igo,iend

         NCHAR = IEND - IGO + 1
         IF (NCHAR > MAXLENTOK) THEN
            CALL ERRT(102,'MAX LENGTH OF TOKEN',MAXLENTOK)
            RETURN
         ENDIF

C        PUT TOKEN IN LIST FOR RETURN
         TOKENS(NGOT) = RECLIN(IGO:IEND)

         IF (NGOT >= NTOKSMAX) EXIT
      ENDDO

      IRTFLG = 0

      END

      SUBROUTINE GET_TOKENS_INTS(RECLIN,NLET,MAXLENTOK,NTOKSMAX,
     &                      INSIDE,     DELIM, TOKENS, NGOT,   IRTFLG)

C      FILLS TOKENS ARRAY WITH CONSECUTIVE TOKENS FROM A LINE. SINCE
C      >7 DIGIT INTEGERS MAY NOT BE EXACTLY REPRESENTED BY A FLOAT THIS
C      WILL PARTITION SUCH FIELDS INTO MULTIPLE TOKENS. 

       IMPLICIT NONE
       INCLUDE 'CMBLOCK.INC'

       CHARACTER(LEN=NLET)      :: RECLIN          ! TOKEN LINE  (SENT)
       INTEGER                  :: NLET            ! LINE LENGTH (SENT)
       INTEGER                  :: MAXLENTOK,NTOKSMAX  ! (SENT)
       CHARACTER(LEN=MAXLENTOK) :: TOKENS(NTOKSMAX)! TOKEN FIELDS
       LOGICAL                  :: INSIDE          ! DELIM IS INSIDE (SENT)
       CHARACTER(LEN=*)         :: DELIM           ! TOKEN DELIMITERS  (SENT)
       INTEGER                  :: NGOT            ! NO. OF TOKENS (RETURNED
       INTEGER                  :: IRTFLG          ! ERROR FLAG     (RETURNED)

       INTEGER                  :: IGO,IEND,IFIRST,NCHAR         
       INTEGER                  :: I,I1,I2

       NGOT   = 0
       IRTFLG = 1
       IEND   = 0

       !write(6,*) '  reclin: ',RECLIN(1:NLET)

       DO      ! ----------------- LOOP --------------------------
         IFIRST = IEND + 1

         IF (INSIDE) THEN
C           GET TOKEN (CHAR. STRING CONTAINS SYMBOLS IN: DELIM)
            CALL GETNEXTTOKEN_NOT(RECLIN(1:NLET),
     &                            DELIM,IFIRST,IGO,IEND)
         ELSE

C          GET TOKEN (CHAR. STRING DELIMITED BY SYMBOLS IN: DELIM)
            CALL GETNEXTTOKEN_D(RECLIN(1:NLET),
     &                            DELIM,IFIRST,IGO,IEND)
         ENDIF
         !write(6,*) '  ifirst,go,end,nlet: ',ifirst, igo, iend, nlet

C        SEE IF ALL TOKENS FROM RECLIN STRING HAVE BEEN EVALUATED
         IF (IGO <= 0) EXIT

C        TOKEN RETURNED, SET TOKEN IN TOKENS

         NCHAR = IEND - IGO + 1
         I1    = IGO

         IF (NCHAR > 7) THEN
C           INEXACT IN DOC FILE, SPLIT

            !IFIELDS = NCHAR / 7
            !IF ((MOD(NCHAR,7)) > 0) IFIELDS = IFIELDS + 1

             DO 
               I2 = I1 + 6
               IF (I2 > IEND) I2 = IEND
               NGOT = NGOT + 1

               IF (NGOT > NTOKSMAX) THEN
                  CALL ERRT(102,'MAX NO. OF TOKENS',NTOKSMAX)
                  RETURN
               ENDIF

C              PUT TOKEN IN LIST FOR RETURN
               TOKENS(NGOT) = RECLIN(I1:I2)

               !write(6,*) '  nchar,i1...i2: ',nchar,i1,i2
               !write(6,'(A,i0,a,a)')'  Token(',ngot,'):',reclin(i1:i2)

               IF (I2 >= IEND) EXIT
               I1 = I2 + 1

            ENDDO

         ELSE
            NGOT = NGOT + 1
            IF (NGOT > NTOKSMAX) THEN
               CALL ERRT(102,'MAX NO. OF TOKENS',NTOKSMAX)
               RETURN
            ENDIF

            !write(6,'(A,i0,a,a)')'  Token(',ngot,'):',reclin(igo:iend)
            !!write(6,*) '  igo...iend: ',igo,iend

            NCHAR = IEND - IGO + 1
            IF (NCHAR > MAXLENTOK) THEN
               CALL ERRT(102,'MAX LENGTH OF TOKEN',MAXLENTOK)
               RETURN
            ENDIF

C           PUT TOKEN IN LIST FOR RETURN
            TOKENS(NGOT) = RECLIN(IGO:IEND)

            IF (NGOT >= NTOKSMAX) EXIT
            I1 = IEND + 1
            !write(6,*) ' ---------- got, i1:', i1
         ENDIF

      ENDDO

      IRTFLG = 0

      END

#ifdef NEVER
#endif
