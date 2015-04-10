
C++*********************************************************************
C
C  INQUIREARCH.F   -- CREATED March 4 2002    ARDEAN LEITH
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
C  INQUIREARCH(LUNOLD,FLIP,FOLD,IRTFLG)
C
C  PURPOSE:  PRINTS OUT DIFFERNT WORD REPREPESENTATIONS FROM A RAW
C            FILE
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

       SUBROUTINE INQUIREARCH(LUNOLD,FLIP,FOLD,IRTFLG)

       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC'

       CHARACTER(LEN=MAXNAM) :: FILNAM
       CHARACTER (LEN=1) ::     NULL
       LOGICAL ::               FLIP,FOLD,BIGENDED,BIGEND

       INTEGER *1      I1BUF(4)

       CHARACTER (LEN=1) ::     CCHAR(4)
       REAL *4         R4VAL
       INTEGER *4      I4VAL
       INTEGER *2      I2VAL(2)
       INTEGER *1      I1VAL(4)
       EQUIVALENCE     (R4VAL,I1VAL)
       EQUIVALENCE     (I4VAL,I1VAL)
       EQUIVALENCE     (I2VAL,I1VAL)
       EQUIVALENCE     (CCHAR,I1VAL)

       INTEGER *1      I1REV(4)
       REAL *4         R4REV
       INTEGER *4      I4REV
       EQUIVALENCE     (R4REV,I4REV)
       EQUIVALENCE     (I1REV,I4REV)

       INTEGER *1      I1FLIP(4)
       INTEGER *2      I2FLIP(2)
       REAL *4         R4FLIP
       INTEGER *4      I4FLIP
       EQUIVALENCE     (R4FLIP,I4FLIP)
       EQUIVALENCE     (I2FLIP,I4FLIP)
       EQUIVALENCE     (I1FLIP,I4FLIP)

       INTEGER *2      I1FOLD(4)
       INTEGER *4      I2FOLD(2)

       INTEGER *1      I1FLIPB(4)
       REAL *4         R4FLIPB
       INTEGER *4      I4FLIPB
       EQUIVALENCE     (R4FLIPB,I4FLIPB)
       EQUIVALENCE     (I1FLIPB,I4FLIPB)

       INTEGER *4      I2FF(2)

        NULL = CHAR(0)
        IERR = 0
        NMAX = NIMAX

        BIGEND = BIGENDED(0)
        WRITE(NOUT,*) '  '
        IF (BIGEND) THEN
           WRITE(NOUT,*) ' CURRENT ARCHITECTURE: BIG-ENDED'
        ELSE
           WRITE(NOUT,*) ' CURRENT ARCHITECTURE: LITTLE-ENDED'
        ENDIF

C       OPEN  FILE AS DIRECT ACCESS, UNFORMATTED, RECL= 1 BYTES
        LENREC = 1
        CALL OPAUXFILE(.TRUE.,FILNAM,DATEXC,LUNOLD,LENREC,'O',
     &                       'INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
 
        NVAL = NMAX
        CALL RDPRAI(INUMBR,NMAX,NVAL,1,NMAX,
     &              'ENTER STARTING-ENDING BYTE',NULL,IRTFLG)
        
        DO I = 1,NVAL
           DO IBYTE = 1,4
              I1FLIP(2) = 0
              I1VAL(I)  = 0
              I1FOLD(I) = 0
           ENDDO

           DO IBYTE = 1,4
              ILOC = INUMBR(I) + IBYTE - 1
              READ(LUNOLD,REC=ILOC,IOSTAT=IERR) I1BUF(1)

              IF (IERR .NE. 0) THEN
                 WRITE(NOUT,90) IERR,ILOC
90               FORMAT('*** ERROR: (',I4,') READING LOCATION: ',I8 )
                 CALL ERRT(100,'INQUIREARCH',NE)
                 GOTO 9999
              ENDIF

C             NO FLIP
              I1VAL(IBYTE) = I1BUF(1)

C             INVERT BYTE ORDER
              I1REV(5-IBYTE) = I1BUF(1)

C             FLIP BYTES WITHIN WORDS
              IF (IBYTE .EQ. 1) THEN
                 I1FLIP(2) = I1BUF(1)
              ELSEIF (IBYTE .EQ. 2) THEN
                 I1FLIP(1) = I1BUF(1)
              ELSEIF (IBYTE .EQ. 3) THEN
                 I1FLIP(4) = I1BUF(1)
              ELSEIF (IBYTE .EQ. 4) THEN
                 I1FLIP(3) = I1BUF(1)
              ENDIF

C             FLIP BYTES AND WORDS
              IF (IBYTE .EQ. 1) THEN
                 I1FLIPB(3) = I1BUF(1)
              ELSEIF (IBYTE .EQ. 2) THEN
                 I1FLIPB(4) = I1BUF(1)
              ELSEIF (IBYTE .EQ. 3) THEN
                 I1FLIPB(1) = I1BUF(1)
              ELSEIF (IBYTE .EQ. 4) THEN
                 I1FLIPB(2) = I1BUF(1)
              ENDIF
           ENDDO

C          FOLD WITHIN BYTES 
           DO IBYTE = 1,4 
              I1FOLD(IBYTE) = I1VAL(IBYTE)
              IF (I1FOLD(IBYTE) .LT. 0)
     &           I1FOLD(IBYTE) = 256 + I1FOLD(IBYTE)
           ENDDO

C          FOLD WITHIN 2 BYTE WORDS
           DO IWORD = 1,2 
             I2FOLD(IWORD) = I2VAL(IWORD)
             IF (I2FOLD(IWORD) .LT. 0)
     &           I2FOLD(IWORD) = 65536 + I2FOLD(IWORD)

C            FLIP BYTES WITHIN WORDS THEN FOLD
             I2FF(IWORD) = I2FLIP(IWORD)
             IF (I2FF(IWORD) .LT. 0)
     &           I2FF(IWORD) = 65536 + I2FF(IWORD)
           ENDDO


           WRITE(NOUT,*) '  '
           WRITE(NOUT,*) ' --- STARTING AT LOCATION: ',INUMBR(I),' ----'

           WRITE(NOUT,93) CCHAR(1),CCHAR(2),CCHAR(3),CCHAR(4)
93         FORMAT('  (CHARS):                    (',4A1,')')

           WRITE(NOUT,91) I1VAL(1),I1VAL(2), I1VAL(3),I1VAL(4)
91         FORMAT('  INT*1:                   ',4(I5,'  '))

           WRITE(NOUT,92) I1FOLD(1),I1FOLD(2), I1FOLD(3),I1FOLD(4)
92         FORMAT('  FOLD INT*1:              ',4(I5,'  '))
           WRITE(NOUT,*) ' '

           WRITE(NOUT,*) ' INT*2:                     ',I2VAL(1),
     &                                                  I2VAL(2)
           WRITE(NOUT,*) ' FLIP INT*2:                ',I2FLIP(1),
     &                                                  I2FLIP(2)
           WRITE(NOUT,*) ' FOLD INT*2 (BY WORDS):     ',I2FOLD(1),
     &                                                  I2FOLD(2)
           WRITE(NOUT,*) ' FLIP & FOLD INT*2:         ',I2FF(1),
     &                                                  I2FF(2)
           WRITE(NOUT,*) ' '

           WRITE(NOUT,*) ' INT*4:                     ',I4VAL
           WRITE(NOUT,*) ' FLIP BYTES INT*4:          ',I4FLIP
           WRITE(NOUT,*) ' REVERSE BYTES INT*4:       ',I4REV
CC           WRITE(NOUT,*) ' FOLD INT*4 (BY BYTES):     ',I4FOLD
           WRITE(NOUT,*) ' FLIP BYTES & WORDS INT*4:  ',I4FLIPB
           WRITE(NOUT,*) ' '

           WRITE(NOUT,*) ' REAL*4   :                 ',R4VAL
           WRITE(NOUT,*) ' FLIP BYTES REAL*4:         ',R4FLIP
           WRITE(NOUT,*) ' REVERSE BYTES REAL*4:      ',R4REV
           WRITE(NOUT,*) ' FLIP BYTES & WORDS REAL*4: ',R4FLIPB

           WRITE(NOUT,*) ' '

        ENDDO
       
9999    CLOSE(LUNOLD)
        END


 
