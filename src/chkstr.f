
C++*******************************************************************
C
C    CHKSTR.FOR       ADAPTED FROM CHKSTR        AUG 89  ArDean Leith
C                     IGNORE COMMENT NOW         AUG 99  ArDean Leith
C                     GETNEXTTOKEN ACCEPTS TAB   NOV 99 ARDEAN LEITH
C                     SPLIT GETNEXTTOKEN TO FILE JAN 00 ARDEAN LEITH
C                     'RR' X? BINDING TIME BUG   NOV 01 ARDEAN LEITH
C                     [] REG SUPPORT             NOV 05 ARDEAN LEITH
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
C    CHKSTR(STRING,NCHAR,TYPE,NUML,FNUML,NMAX,INUM,IRTFLG)
C
C    PURPOSE:           EXTRACTS INTERGER OR REAL VALUES FROM STRING.  
C                       CAN SUBSTITUTE FOR REGISTER CONTENT OR RETURN
C                       REGISTER NUMBER LISTING         
C    PARAMETERS:
C	  STRING	INPUT STRING                             (SENT)
C	  NCHAR		# OF CHARACTERS IN INPUT STRING          (SENT)
C	  TYPE		TYPE OF NUMBERS RETURNED                 (SENT)
C                       I  = INTEGER LIST OF REGISTER CONTENTS
C                       R  = FLOATING LIST
C                       IT = REGISTER NUMBERS IN NUML
C                       RE = FLOATING LIST OF REGISTER CONTENTS
C                            AND INTEGER LIST OF REGISTER NUMBERS
C
C	  NUML		INTEGER ARRAY, IF NEEDED                 (RET.)
C	  FNUML		REAL ARRAY, IF NEEDED                    (RET.)
C	  NMAX          SIZE OR NUML AND FNUML                   (SENT)
C	  INUM		# OF NUMBERS RETURNED IN NUML OR FNUML   (RET.)
C	  IRTFLG	ERROR FLAG, =0 IF ALL O.K.               (RET.)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE CHKSTR(STRING,NCHAR,TYPE,NUML,FNUML,NMAX,INUM,IRTFLG)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=*)  :: STRING,TYPE
      INTEGER           :: NUML(NMAX)
      REAL              :: FNUML(NMAX)
      LOGICAL           :: RFLAG

      INUM   = 0
      IRTFLG = 1
      IFIRST = 1

      IF (TYPE .EQ. 'IRRE' .OR. TYPE .EQ. 'IR') THEN
C         OBSOLETE
          WRITE(NOUT,*) ' *** PGM ERROR OBSOLETE TYPE: ',TYPE,
     &                  ' IN CHKSTR'
          CALL ERRT(101,'PGM ERROR OBSOLETE TYPE IN CHKSTR',NE)
          RETURN     
      ENDIF

      DO
c       GET TOKEN (CHAR. STRING DELIMITED BY A ", ()") FROM STRING
        CALL GETNEXTTOKEN(STRING(1:NCHAR),IFIRST,IGO,IEND)
        IF (IGO .LE. 0) THEN
C          ALL TOKENS FROM STRING HAVE BEEN EVALUATED
           IRTFLG = 0
           RETURN
        ENDIF

C       TOKEN RETURNED
          
        RFLAG = .FALSE.
        IF (STRING(IGO:IGO) .EQ. '[')THEN
C          TOKEN STARTS WITH '[', SO IS A REGISTER REFERENCE
           IF (IGO .GE. NCHAR) RETURN
           RFLAG = .TRUE.
        ENDIF

        INUM = INUM + 1
        IF (INUM .GT. NMAX) THEN
C          PREVENT NUML OVERFLOW!!
           WRITE(NOUT,91) NMAX
91         FORMAT('*** ONLY: ',I6,
     &            '  VALUES RETURNED BY CHKSTR TO PREVENT OVERFLOW'/)
           RETURN
        ENDIF

        IF (RFLAG) THEN
C          GET REGISTER VALUE & NUMBER FROM TOKEN
           CALL REG_GET_VAR(0,STRING(IGO:IEND),.FALSE.,FVAL,IREG,
     &                      IEND,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
           IEND = IEND + IGO -1

           IF (TYPE .EQ. 'R') THEN
C             RETURN FLOATING POINT IN FNUML
              FNUML(INUM) = FVAL

           ELSEIF (TYPE .EQ. 'RE') THEN
C             RETURN FLOATING POINT IN FNUML & IREG IN INUML
              FNUML(INUM) = FVAL
              NUML(INUM)  = IREG

           ELSEIF (TYPE .EQ. 'I') THEN
C             EXTRACT FLOATING POINT AND RETURN AS INTEGER 
              NUML(INUM) = FVAL

           ELSEIF (TYPE .EQ. 'IT') THEN
C             RETURN REGISTER NUMBERS IN NUML, (NOT REGISTER CONTENTS)
C             ACCEPTS [] INPUT AND RETURNS REGISTER NUMBER             
              NUML(INUM) = IREG 
           ENDIF

        ELSEIF (TYPE .EQ. 'I') THEN
C          EXTRACT INTEGER      
           READ(STRING(IGO:IEND),*,IOSTAT=IERR) NUML(INUM)
           IF (IERR .NE. 0) THEN
              WRITE(NOUT,*) 'BAD INPUT STRING: ',STRING(IGO:IEND)
              RETURN
           ENDIF
            IF (IERR .NE. 0) RETURN

        ELSEIF (TYPE .EQ. 'R' .OR. TYPE .EQ. 'RE') THEN
C          EXTRACT REAL NUMBER    
           READ(STRING(IGO:IEND),*,IOSTAT=IERR) FNUML(INUM)
           IF (IERR .NE. 0) THEN
              WRITE(NOUT,*) 'BAD INPUT STRING: ',STRING(IGO:IEND)
              RETURN
c             WRITE(NOUT,*) 'NO. CHARACTERS IN STRING: ',NCHAR
c             WRITE(NOUT,*) 'IFIRST: ',IFIRST,'  INUM: ',INUM,
c     &                     '  IGO: ',IGO,' IEND: ',IEND
           ENDIF
           IF (TYPE .EQ. 'RE') NUML(INUM) = 0
        ENDIF
        IFIRST = IEND + 1

        IF (IFIRST .GT. NCHAR) EXIT
      ENDDO

79    IRTFLG = 0

      RETURN
      END



        
             
