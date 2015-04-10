C++************************************************************ 10/11/79
C
C READRQ.F    FCHAR                              AUG 1989 ARDEAN LEITH
C             USED ERRT                          SEP 1997 ARDEAN LEITH
C             HANDLES 20 REGISTERS               JUL 1999 ARDEAN LEITH
C             USED REG_                          AUG 2000 ARDEAN LEITH
C             ADDED PROMPT                       SEP 2000 ARDEAN LEITH
C             SIMPLIFED WITH SETSYMPAR IN RDPR   APR 2001 ARDEAN LEITH
C             NATIVE [] VARIABLES                NOV 2005 ARDEAN LEITH
C             GLOBAL [] VARIABLES                DEC 2005 ARDEAN LEITH
C             MAXNSEL TO 24                      JAN 2006 ARDEAN LEITH
C             'RR S'                             APR 2006 ARDEAN LEITH
C             MAXNVAL TO 150                     DEC 2011 ARDEAN LEITH
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
C READRQ
C
C PURPOSE:   PLACE SOLICITED VALUES INTO REGISTERS
C
C PARAMETERS:   
C
C--*********************************************************************

      SUBROUTINE READRQ()

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC' 

      INTEGER             :: NSEL_USED,IDUM,ILEVEL,IVAL,IRTFLG
      INTEGER             :: NGOT,I,IV,NOT_USED
      REAL                :: FDUM

C     DANGER!! MAXNVAL MUST BE GREATER THAN MAXNSEL, SET IN reg_set.f & spider.f(NPARG)!!
      INTEGER, PARAMETER  :: MAXNVAL = 150    
      REAL                :: VALUES(MAXNVAL) 

C     SETTING A  REGISTER VARIABLE(S) FROM OPERATION LINE

C     COPY VALUES TO LOCAL REGISTER(S)

C     COUNT NUMBER OF REGISTERS USED IN NSEL ARRAY
      CALL REG_GET_USED(NSEL_USED)
      IF (NSEL_USED .LE. 0) THEN
         CALL ERRT(101,'NO REGISTERS ON OPERATION LINE',IDUM)
         RETURN 
      ENDIF

C     CAN FILL OR READ GLOBAL ALSO
      ILEVEL = 0
      IF (FCHAR(1:4) .EQ. 'RR C') ILEVEL = 1

      IF (FCHAR(1:4) .NE. 'RR S') THEN
C        COPY A LIST OF VALUES TO OP. LINE REGISTER(S)

C        INPUT LIST OF VALUES FOR THE SELECTED REGISTERS
         CALL RDPRA('VALUE(S)',NSEL_USED,ILEVEL,.FALSE.,
     &               VALUES,NGOT,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

C        SET DESIRED REGISTERS TO: VALUES(...)
         CALL REG_SET_NSELA(NGOT,VALUES,.FALSE.,IRTFLG)

      ELSE
C        COPY ONLY DESIRED VALUE(S) TO OP. LINE REGISTER(S)

         CALL RDPRA('VALUE(S)',MAXNVAL,ILEVEL,.FALSE.,
     &              VALUES,NGOT,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         IVAL = 1
         CALL RDPRI1S(IVAL,NOT_USED,'POSITION WITHIN LIST',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

C        COPY VALUES TO OP. LINE REGISTER(S)
         DO I =IVAL,IVAL+NSEL_USED-1

           IF (I .LE. NGOT) THEN
              IV = I-IVAL+1
              CALL REG_SET_NSEL(IV,1,VALUES(I),
     &                          FDUM,FDUM,FDUM,FDUM,IRTFLG)
              WRITE(NOUT,*) ' REGISTER SET TO: ',VALUES(I)
           ELSE
             CALL ERRT(102,' POSITION NOT AVAILABLE IN LIST',I)
             RETURN
           ENDIF
         ENDDO
      ENDIF

      END

c     write(6,*) ' nsel_used:',nsel_used
c     write(6,*) ' igot:',igot
c     write(6,*) ' values(1):',values(1)


