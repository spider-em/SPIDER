
C++*********************************************************************
C
C   EXPRESS3Q.F             FILENAMES LENGTHENED   JAN 89 ArDean Leith
C                           MERGED WITH EXPRESSQ   FEB 99 ArDean Leith
C                           NGOT ADDED             MAY 01 ArDean Leith
C                           FROM EXPRESS3Q         NOV 05 ARDEAN LEITH
C                           REWRITTEN              NOV 05 ARDEAN LEITH
C                           ALLOW BLANKS           JAN 06 ARDEAN LEITH
C                           ! COMMENT DELIMITER    DEC 11 ARDEAN LEITH
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
C     EXPRESS3Q(STRING,ISTACK,NVAL,VALUES,NGOT,INPARLOOP,IRTFLG)
C
C     PURPOSE:    EVALUATES EXPRESSIONS IN STRING.
C
C     PARAMETERS: STRING        INPUT STRING                     (SENT)
C                 ISTACK        REGISTER STACK LEVEL             (SENT)
C                 NVAL          MAX NUMBER OF RETURNED VALUES    (SENT)
C                 VALUES        RETURNED VALUES                  (RET.)
C                 NGOT          NUMBER OF RETURNED VALUES        (RET.)
C                 INPARLOOP     INSIDE () FOR LOOP               (RET.)
C                 IRTFLG        ERROR FLAG                       (RET.)
C
C--*******************************************************************
       
	SUBROUTINE EXPRESS3Q(STRING,ISTACK,NVAL,
     &                     VALUES,NGOT,INPARLOOP,IRTFLG)

	INCLUDE 'CMBLOCK.INC' 

        CHARACTER(LEN=*)      :: STRING
        REAL                  :: VALUES(NVAL)
	LOGICAL               :: INPARLOOP

        INTEGER,PARAMETER     :: LENEXP  = 80  ! EXPR LENGTH LIMIT
        CHARACTER(LEN=LENEXP) :: EXPR
	LOGICAL               :: HASCOMMA

        NGOT  = 0

C       ASSUME COMMENT (AFTER SEMICOLON) ALREADY REMOVED
        NCHAR = LEN(STRING)

C       AN EXPRESSION CAN BE DELIMINATED BY BLANK, ',', OR ';'
        LOCCOMMA  = INDEX(STRING,',')
        HASCOMMA  = (LOCCOMMA .GT. 0)

        INPARLOOP = .FALSE.
        I         = 1
        DO WHILE (I .LE. NCHAR)

C         GET NEXT EXPRESSION IN STRING
          CALL GETNEXTEXP(STRING(1:NCHAR),I,HASCOMMA,
     &                    IGORET,EXPR,NLET,IRTFLG)
          IF (IRTFLG .NE. 0) THEN
             WRITE(NOUT,*) ' *** INVALID EXPRESSION IN INPUT: ',
     &                     STRING(1:NCHAR)
             RETURN
          ENDIF

C         EXIT LOOP IF NO MORE EXPRESSION(S) RETURNED
          IF (NLET .LE. 0) EXIT

C         CHECK FOR INTIAL OR FINAL UNBALANCED ()
          NLEFPAR = 0
          NRITPAR = 0
          DO J=1,NLET
            IF (EXPR(J:J) .EQ. '(') NLEFPAR = NLEFPAR + 1
            IF (EXPR(J:J) .EQ. ')') NRITPAR = NRITPAR + 1
          ENDDO

          IF (NLEFPAR .GT. NRITPAR .AND. NGOT .EQ. 0 .AND.
     &        EXPR(1:1) .EQ. '(') THEN
              INPARLOOP = .TRUE.
C             STRIP EXTRA INITIAL PARENTHESIS
              EXPR = EXPR(2:NLET)
              NLET = NLET - 1

          ELSEIF (NLEFPAR .LT. NRITPAR .AND. EXPR(NLET:NLET).EQ.')')THEN
C             STRIP EXTRA FINAL PARENTHESIS
              EXPR = EXPR(1:NLET-1)
              NLET = NLET - 1

          ENDIF

C         EVALUATE THE EXPRESSION IN THIS SUB-STRING
          CALL EVALNQ(ISTACK,EXPR(1:NLET),NLET,RETVALUE,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

C         EXPRESSION RETURNS A VALID VALUE
          NGOT         = NGOT + 1
          VALUES(NGOT) = RETVALUE
          I            = IGORET
          IF (NGOT .GE. NVAL) RETURN
       ENDDO

       IF (STRING(1:1) .EQ. '(' .AND. STRING(NCHAR:NCHAR) .EQ. ')'.AND.
     &      NGOT .EQ. 1) INPARLOOP = .TRUE.

       END
       

C++************************ GETNEXTEXP ***************************
C
C    GETNEXTEXP.F  GETNEXTEXP ACCEPTS TAB NOV 99 ARDEAN LEITH
C
C **********************************************************************
C
C    GETNEXTEXP(STRING,IFIRST,HASCOMMA,
C               INEXT,STRINGOUT,NCHAR,IRTFLG)
C
C    PURPOSE:           GETS NEXT EXPRESSION DELIMITED BY A COMMA, 
C                       SEMICOLON, OR OPTIONALLY BY A BLANK, 
C                       FROM STRING, RETURNS EXPRESSION STRING
C                       NO MORE EXP. LEFT IN STRING, NCHAR IS ZERO.
C
C    PARAMETERS:
C	  STRING	INPUT STRING                             (SENT)
C	  IFIRST        FIRST CHAR TO BE SCANNED FOR EXP.        (SENT)
C	  INEXT		START OF NEXT EXPRESSION                 (RET.)
C	  HASCOMMA	HAS , BETWEEN EXPRESSIONS                (SENT)
C	  STRINGOUT	RETURNED EXPRESSION                      (RET.)
C	  NCHAR		LENGTH OF RETURNED EXPRESSION            (RET.)
C	  IRTFLG	ERROR FLAG                               (RET.)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE GETNEXTEXP(STRING,IFIRST,HASCOMMA,
     &                      INEXT,STRINGOUT,NCHAR,IRTFLG)

C     EXPRESSION DELIMITER IS A COMMA, SEMICOLON, & OPTIONALLY A BLANK

      CHARACTER(LEN=*) :: STRING,STRINGOUT
      CHARACTER(LEN=1) :: CTEMP
      LOGICAL          :: HASCOMMA

C     SET DEFAULT RETURN VALUES
      INEXT   = 0
      NCHAR   = 0

C     FIND LAST CHAR POSITION IN STRING
      ILEN = LEN(STRING)

      DO I =IFIRST,ILEN 
         CTEMP = STRING(I:I)
         INEXT = I + 1

         IF (CTEMP .EQ. ';' .OR. CTEMP .EQ. '!') THEN
C           START OF COMMENT, THIS IS END
            EXIT

         ELSEIF (CTEMP .EQ. ',') THEN
C           EXPRESSION END ENCOUNTERED
            EXIT

         ELSEIF (HASCOMMA .AND.
     &           NCHAR .GT. 0 .AND.
     &           CTEMP .LE. CHAR(31) ) THEN
C           END AT A NON-PRINTING CHAR.
            EXIT

         ELSEIF (.NOT. HASCOMMA .AND.
     &          (NCHAR .GT. 0 .AND.
     &          (CTEMP .EQ. ' ' .OR. CTEMP .LE. CHAR(31)))) THEN
C           END AT THE BLANK, OR NON-PRINTING CHAR.
            EXIT

         ELSEIF (CTEMP .NE. ' ' .AND. CTEMP .GT. CHAR(31)) THEN
C           INSIDE AN EXPRESSION (IGNORES BLANKS)
            NCHAR                  = NCHAR + 1
            STRINGOUT(NCHAR:NCHAR) = CTEMP 
         ENDIF
      ENDDO

C     EXPRESSION FOUND OR RAN OFF END OF STRING
      IRTFLG = 0

      RETURN
      END


        
     
