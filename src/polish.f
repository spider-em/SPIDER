C++********************************************************************* 
C
C POLISH.F                    
C              CHANGED SOME ARRAYS FROM CHAR. TO INT FOR IBM JAN 2000 AL
C              READ(EXPR,FMTR LEAK                           AUG 2002 AL
C              SIMPLIFIED NUMBER INTERPRETATION              AUG 2002 AL
C              STACK LEVEL & LOWERCASE                       DEC 2005 AL
C              REMOVED: 'NO REGSITER VAR.......              NOV 2009 AL
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
C  POLISH(ILEVEL,EXPR,NLET,IPOLSH,NPOL,VAL,NVAL,IRTFLG)
C
C  PURPOSE: PARSE INFIX EXPESSION INTO POSTFIX EXPESSION
C
C  PARAMETERS:
C	ILEVEL     STACK LEVEL                                   (SENT)
C       EXPR       CHARACTER STRING CONTAINING EXPRESSION        (SENT)
C	NLET       LENGTH OF EXPR                                (SENT)
C	IPOLSH     INT ARRAY RETURNS POSTFIX EXPRESSION      (RETURNED)
C       NPOL       NO. OF ELEMENTS IN IPOLSH ARRAY           (RETURNED)
C       VAL        ARRAY STORES VALUES WHICH INDEX BY        (RETURNED)
C                        POLISH'S ELEMENTS
C       NVAL       NO. OF ELEMENTS USED IN VAL               (RETURNED)
C       IRTFLG     ERROR FLAG                                (RETURNED)
C
C FOR VAL INDEX
C ASCII(N) - ASCII(0) FOR N FROM 1 TO Z ON THE ASCII TABLE
C FOR MATH FUNCTION ARGUMENTS
C       PAD -> a       
C       SIN -> b
C       EXP -> c
C       LOG -> d
C       COS -> e
C       SQR -> f
C       LON -> g (NATURAL LOGARITHM)
C       INT -> h
C       ABS -> i
C       ATA -> j
C       ASI -> k
C       ACO -> l 
C       TAN -> m
C       RAN -> n   [0,1] uniform distribution
C       RNN -> o   (0,1) normal distribution
C       CHANGE THE SIGN -> p
C
C THE SUBROUTINE TRANSFORMS THE EXPRESSION FROM A TOKEN TO A SINGLE
C FOR EXAMPLE EXPRESSION 43+COS(6) -> 1+e(2)
C CONVERT FROM INFIX TO POSTFIX NOTATION
C  12E+	
C
C  NOTE:  TO ADD A NEW MATH FUNCTION:
C     1. CHOOSE A SEQUENTIAL LETTER FOR SUBSTITUTING A MATH FUNCTION
C
C **********************************************************************

	SUBROUTINE POLISH(ILEVEL,EXPR,NLET,IPOLSH,NPOL,
     &                        VAL,NVAL,IRTFLG)

        COMMON /UNITS/  LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT

        CHARACTER(LEN=*)  :: EXPR

        INTEGER,PARAMETER :: IVALEN  = 40  ! RPN LENGTH LIMIT
        INTEGER,PARAMETER :: IRPNLEN = 80  ! RPN LENGTH LIMIT
        INTEGER,PARAMETER :: NFUNC  = 15   ! 

	INTEGER           :: IPOLSH(IRPNLEN),IEXPR1(IRPNLEN)
	INTEGER           :: ISTACK(IRPNLEN),ISTP(10)
        REAL              :: VAL(IVALEN)
	LOGICAL           :: FP,BGS,UNITARY,POWER

        CHARACTER *80   MSG
        CHARACTER *7    FRMT
        CHARACTER *3    TCHR
        CHARACTER *1    CTEMP,CNEXT,MINUS
	CHARACTER *3    FUNT(NFUNC)

	DATA  FUNT/
     &  'PAD','SIN','EXP','LOG','COS','SQR','LON','INT',
C         a     b     c     d     e     f     g     h
     &  'ABS','ATA','ASI','ACO','TAN','RAN','RNN'/
C         i     j     k     l     m     n     o
C       'p'  RESERVED FOR UNITARY OPERATIONS!
	DATA  MINUS/'p'/

C       SET ERROR RETURN
        IRTFLG = 1

        IEXPR1  = 0
        ISTACK  = 0
        IPOLSH  = 0

	K       = 1
	J       = 1
	I       = 0
	BGS     = .TRUE.
	UNITARY = .FALSE.

C       POSITION ON A STACK (ISTP) OF UNITARY OPERATIONS
	LISTP   = 0

C       PARANTHESIS NESTING LEVEL IS ON THE STACK
        POWER = .FALSE.
        LPOW  = 0

C       -------------- LOOP ----------------------------------------
C       LEXICAL ANALYSIS	
40	I = I + 1
	IF (I .GT. NLET) GOTO 41

        CTEMP = EXPR(I:I)

C       ALL THE SPACES ARE ASSUMED TO BE REMOVED

	IF (BGS) THEN
           IF (CTEMP.EQ. '+')  THEN
C             IGNORE UNITARY +
	      GOTO 41

           ELSEIF (CTEMP.EQ.'-')  THEN
C             UNITARY OPERATION   MINUS SIGN
              IEXPR1(K)   = ICHAR(MINUS)
              IEXPR1(K+1) = ICHAR('(')
              K = K + 2
C             HERE PUSH ON STACK
              LISTP = LISTP+1
              ISTP(LISTP) = 0
              UNITARY     = .TRUE.
              I           = I + 1
	      IF (I .GT. NLET)  THEN
                 CALL ERRT(101,' EXPRESSION CANNOT END WITH -',NE)
                 RETURN
              ENDIF
              CTEMP = EXPR(I:I)

           ELSE
              BGS = .FALSE.           
           ENDIF
	ENDIF

	IF (CTEMP .EQ. '(') THEN
C          LEFT PARENTHESIS
           IF (POWER) LPOW = LPOW+1
           IF (UNITARY) THEN
C             BEGIN NEW EXPRESSION
              ISTP(LISTP) = ISTP(LISTP)+1
              BGS = .TRUE.
           ELSE
              BGS = .TRUE.
           ENDIF
           IEXPR1(K) = ICHAR(CTEMP)
           K          = K + 1
           GOTO 41

	ELSEIF (CTEMP.EQ.')') THEN
C          RIGHT PARENTHESIS
           IEXPR1(K) = ICHAR(CTEMP)
           K          = K + 1
           IF (UNITARY)  THEN
              BGS = .FALSE.
C             POP FROM THE STACK, IF END OF STACK UNITARY=.FALSE.
              ISTP(LISTP) = ISTP(LISTP)-1
              IF (ISTP(LISTP) .EQ. 0)  THEN
                 LISTP = MAX(LISTP-1,0)
                 IF (LISTP .EQ. 0)  UNITARY = .FALSE.
C                SEE IF THE NEXT OPERATION IS '**', IF YES DO NOTHING
                 IF (I+2 .LE. NLET) THEN
                    IF (EXPR(I+1:I+2) .EQ. '**')  THEN
                       IF (POWER)  THEN
                          CALL ERRT(101,
     &                       'TOO MANY NESTED POWER OPERATORS',NE)
                          RETURN
                       ENDIF
                       POWER = .TRUE.
                       LPOW  = 0
                       GOTO 41
                    ENDIF
                 ENDIF
                 IEXPR1(K) = ICHAR(')')
                 K         = K + 1
              ENDIF
           ENDIF
           IF (POWER) THEN
              LPOW = LPOW-1
              IF (LPOW .EQ. 0)  THEN
                 POWER = .FALSE.
                 IEXPR1(K) = ICHAR(')')
                 K         = K + 1
              ENDIF
           ENDIF
           GOTO 41

        ELSEIF (I+1 .LE. NLET)  THEN
           IF (EXPR(I:I+1) .EQ. '**') THEN
              IEXPR1(K) = ICHAR('^')
              K         = K+1
              I         = I+1
              GOTO 41
           ENDIF
        ENDIF

        IF (CTEMP  .EQ. '[') THEN
C          [] IS RESERVED FOR REGISTERS (SYMBOLS ALREADY SUBSTTUTED OUT)

           CALL REG_GET_VAR(ILEVEL,EXPR(I:),.FALSE.,VALDUM,
     &                      IREG,IENDVAR,IER)
	   IF (IER .NE. 0) RETURN

           I         = I + IENDVAR - 1 

          VAL(J)     = IREG
          IEXPR1(K)  = (127+J)
          J          = J+1
          K          = K+1

          IF (UNITARY) THEN
             IF (BGS) THEN
                BGS = .FALSE.
C               POP FROM THE STACK, IF END OF STACK UNITARY= .FALSE.
                LISTP = MAX(LISTP-1,0)
                IF (LISTP .EQ. 0)  UNITARY= .FALSE.
C               SEE IF THE NEXT OPERATION IS '**', IF YES DO NOTHING
                IF (I+2 .LE. NLET) THEN
                   IF (EXPR(I+1:I+2) .EQ. '**')  THEN
                      IF (POWER)  THEN
                         CALL ERRT(101,
     &                   'TOO MANY NESTED POWER OPERATORS',NE)
                         RETURN
                      ENDIF
                      POWER = .TRUE.
                      LPOW  = 0
                      GOTO 41
                   ENDIF
                ENDIF
                IEXPR1(K) = ICHAR(')')
                K         = K + 1
             ENDIF
          ENDIF
          IF (POWER .AND. LPOW .EQ.0) THEN
             POWER = .FALSE.
             IEXPR1(K) = ICHAR(')')
             K         = K + 1
          ENDIF
	
        ELSEIF (CTEMP .GE. 'A') THEN

C          MATH FUNCTIONS PARSING BEGINS HERE
C          SUBSTITUTE LETTERS FOR MATH FUNCTIONS
	
           IF (I+2 .LE. NLET)  THEN
              TCHR = EXPR(I:I+2)
              CALL SSUPCAS(TCHR)  !COULD BE LOWERCASE
              DO  L=1,NFUNC
                 IF (TCHR .EQ. FUNT(L)) GOTO  51
              ENDDO
              GOTO 52

C             CODE FUNCTION BY A SMALL LETTER
51            IEXPR1(K) = (L - 1 + ICHAR('a'))
	      K = K+1
	      I = I+2
              GOTO 41

           ENDIF
52         IF (I+1 .LE. NLET)  THEN
              IF (EXPR(I:I+1) .EQ. 'P1' .OR. EXPR(I:I+1) .EQ. 'p1') THEN
C                PIXEL OPERATIONS FOR 'AR' USE
                 IF (J .GT. IVALEN)  THEN
                    CALL ERRT(101,'EXPRESSION TOO LONG',NE)
                    RETURN
                 ENDIF

C                RESERVE PLACE FOR PIXEL IN VAL ARRAY
                 VAL(J)    = 1.0
                 IEXPR1(K) = (200+J)
                 K = K+1
                 I = I+1
                 J = J+1
                 IF (UNITARY)  THEN
                    IF (BGS) THEN
                       BGS = .FALSE.
C                      POP FROM THE STACK, IF END OF STACK UNITARY= .FALSE.
                       LISTP = MAX(LISTP-1,0)
                       IF (LISTP .EQ. 0)  UNITARY = .FALSE.
C                      SEE IF THE NEXT OPERATION IS '**', 
C                      IF YES DO NOTHING
                       IF (I+2 .LE. NLET) THEN
                          IF (EXPR(I+1:I+2) .EQ. '**')  THEN
                             IF (POWER)  THEN
                                CALL ERRT(101,
     &                          'TOO MANY NESTED POWER OPERATORS',NE)
                                RETURN
                             ENDIF
                             POWER = .TRUE.
                             LPOW  = 0
                             GOTO 41
                          ENDIF
                       ENDIF
                       IEXPR1(K) = ICHAR(')')
                       K         = K + 1
                    ENDIF
                 ENDIF
                 IF (POWER .AND. LPOW .EQ. 0) THEN
                    POWER     = .FALSE.
                    IEXPR1(K) = ICHAR(')')
                    K         = K + 1
                 ENDIF
                 GOTO 41
              ENDIF
           ENDIF

       ELSEIF (CTEMP .EQ. '+' .OR. CTEMP .EQ. '-' .OR.
     &         CTEMP .EQ. '*' .OR. CTEMP .EQ. '/' .OR.
     &         CTEMP .EQ. '^') THEN
C         ARITHMETIC OPERATION +-*/
          IEXPR1(K) = ICHAR(EXPR(I:I))
          K         = K + 1

       ELSEIF (CTEMP .EQ.'.' .OR. (CTEMP.GE.'0'.AND.CTEMP.LE.'9')) THEN
C         A NUMBER IN EXPRESSION

          IGO  = I
          INOT = VERIFY(EXPR(IGO:NLET),'.Ee0123456789')
          I    = NLET
          IF (INOT .GT. 0) THEN
             I = IGO + INOT - 2
             IF (EXPR(I:I) .EQ. 'E' .OR. EXPR(I:I) .EQ. 'e') THEN
C               CAN HAVE INCLUDED '+" or '-'
                INOT  = VERIFY(EXPR(I+2:NLET),'.Ee0123456789')
                IF (INOT .GT. 0) THEN
                   I = I + INOT
                ELSE
                   I = NLET
                ENDIF
             ENDIF
          ENDIF

c          WRITE(NOUT,9009) INOT,IGO,I,EXPR(IGO:I)
c9009      FORMAT('inot: ',i2,'  EXPR(',i2,':',i2,') :'A)

C         EVALUATE THE NUMBER
          IF (J .GT. IVALEN)  THEN
              CALL ERRT(101,'EXPRESSION TOO LONG',NE)
              RETURN
          ENDIF
          IEXPR1(K) = (J+48)
          READ(EXPR(IGO:I),'(F20.0)',IOSTAT=IER) VAL(J)

	   IF (IER .NE. 0)  THEN
              WRITE(NOUT,*) 'IN EXPRESSION: ',EXPR(IGO:I)
              CALL ERRT(101,'READING NUMBER',NE)
              RETURN
           ENDIF

	   J = J + 1
	   K = K + 1

 30	   FORMAT(I1)
        
           IF (UNITARY)  THEN
              IF (BGS) THEN
                 BGS = .FALSE.
C                POP FROM STACK, IF END OF STACK UNITARY= .FALSE.
                 LISTP = MAX(LISTP-1,0)
                 IF (LISTP .EQ.0)  UNITARY= .FALSE.
C                SEE IF THE NEXT OPERATION IS '**', IF YES DO NOTHING
                 IF (I+2 .LE. NLET) THEN
                    IF (EXPR(I+1:I+2) .EQ. '**')  THEN
                       IF (POWER)  THEN
                          CALL ERRT(101,
     &                       'TOO MANY NESTED POWER OPERATORS',NE)
                          RETURN
                       ENDIF
                       POWER= .TRUE.
                       LPOW=0
                       GOTO 41
                    ENDIF
                 ENDIF
                 IEXPR1(K) = ICHAR(')')
                 K          = K + 1
              ENDIF
           ENDIF
           IF (POWER .AND. LPOW .EQ.0) THEN
              POWER= .FALSE.
              IEXPR1(K) = ICHAR(')')
              K          = K + 1
           ENDIF
        ELSE
C          UNEXPECTED CHARACTER IN THE EXPRESSION
           MSG = 'UNEXPECTED CHARACTER IN EXPRESSION: '//CTEMP//CHAR(0)
           NCHARE = lnblnkn(MSG)
           CALL ERRT(101,MSG(:NCHARE),NE)
           RETURN
        ENDIF

C       END OF NUMBER PROCESSING LOOP

41	IF (I .LT. NLET) GOTO 40

C       --------------------- END LOOP ------------------------

C       CONVERT FROM INFIX TO POSTFIX
C       SYNTAX ANALYSIS
        NVAL       = J-1
        IEXPR1(K)  = ICHAR(')')
	NCHAR2     = K+1
	IPNT       = 1
	ITOP       = 1
	ISTACK(1)  = ICHAR('(')
	IRANK      = 0
	I          = 0

4	NEXT = IEXPR1(IPNT)
	IPNT = IPNT+1
	IF (IPNT .GT. NCHAR2) THEN
C          ALL DONE, CAN RETURN NOW
	   IF (ITOP .NE.0) THEN
	      CALL ERRT(43,'POLISH',NE)
              RETURN
	   ENDIF
           J      = 1
           IRTFLG = 0
           NPOL   = I
           RETURN
        ENDIF

	IF (ITOP .LT. 1)  THEN
           WRITE(NOUT,*) ' *** INVALID EXPRESSION: ',EXPR
	   CALL ERRT(101,'INVALID ARITHMETIC EXPRESSION',NE)
           RETURN
	ENDIF

7       CONTINUE
        IF (NEXT .EQ. ICHAR('+') .OR. NEXT.EQ. ICHAR('-')) THEN
          NFQT = 1
        ELSEIF (NEXT .EQ. ICHAR('*') .OR. 
     &          NEXT .EQ. ICHAR('/')) THEN
          NFQT = 3
        ELSEIF (NEXT .EQ. ICHAR('^')) THEN
          NFQT = 6
        ELSEIF (NEXT .EQ. ICHAR('(')) THEN
          NFQT = 9
        ELSEIF (NEXT .EQ. ICHAR(')')) THEN
          NFQT = 0
        ELSE
          NFQT = 7
        ENDIF

        ITEMP = ISTACK(ITOP)
        IF (ITEMP .EQ. ICHAR('+') .OR. ITEMP .EQ. ICHAR('-')) THEN
           NGQT = 2
        ELSEIF (ITEMP .EQ. ICHAR('*') .OR. ITEMP .EQ. ICHAR('/')) THEN
            NGQT = 4
        ELSEIF (ITEMP .EQ. ICHAR('^')) THEN
            NGQT = 5
        ELSEIF (ITEMP .EQ. ICHAR('(')) THEN
            NGQT = 0
        ELSE
            NGQT = 8
        ENDIF

 	IF (NFQT .LE. NGQT) THEN
           IF (NFQT .NE. NGQT) THEN
              I          = I + 1
              IPOLSH(I)  = ITEMP
              NRQT       = 1
              IF (ITEMP .EQ. ICHAR('+') .OR. ITEMP .EQ. ICHAR('-') .OR.
     &            ITEMP .EQ. ICHAR('*') .OR. ITEMP .EQ. ICHAR('/') .OR.
     &            ITEMP .EQ. ICHAR('^')) NRQT = -1

              IRANK = IRANK + NRQT
              IF (IRANK .LE. 0)   THEN
		 CALL ERRT(43,'POLISH',NE)
                 RETURN	
	      ENDIF
              ITOP = ITOP - 1
              GOTO 7
           ENDIF
           ITOP = ITOP - 1
        ELSE
           ITOP         = ITOP+1
           ISTACK(ITOP) = NEXT
        ENDIF
	GOTO 4

	END

