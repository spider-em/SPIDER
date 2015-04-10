
C++*********************************************************************
C
C  SEARCHQ.F -- CHANGED TO STRING VARIABLES        AUG  89 ArDean Leith
C               REMOVED READCH                     JAN  97 ArDean Leith
C               USED GETLBNO                       OCT  99 ArDean Leith
C               SPEEDED UP                         JULY 00 ArDean Leith
C               INCORE PROCS                       JAN  01 ArDean Leith
C               DO ...;COMMENT BUG                 NOV  01 ArDean Leith
C               NO LABEL DO                        NOV  06 ArDean Leith
C               READ(NIN,81,...) FCHAR PGI 10.+    OCT  10 ArDean Leith
C               ! COMMENT DELIMITER                DEC  11 ArDean Leith
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
C    SEARCHQ(LBNOWANT,IRTFLG)
C
C    PURPOSE:         READS INPUT FILE UNTIL IT FIND A LINE HAVING 
C                     'DO LB??' WHERE '??' IS EQUAL TO LBNOWANT
C                     OR 'DO ' WITHOUT A 'LB' AND -LBNOWANT IS THE
C                     LINE NUMBER OF UNLABEL 'DO 's IN FILE
C        
c    PARAMETERS:      LBNOWANT  LB?LINE NUMBER WANTED             (SENT)
C                     IRTFLG    ERROR NUMBER                  (RETURNED)
C
C--*******************************************************************

	SUBROUTINE SEARCHQ(LBNOWANT,IRTFLG)

C       FCHAR IS PASSED BACK TO CALLER IN COMMON!
        INCLUDE 'CMBLOCK.INC'
        COMMON /LUNDOECHO/ LUNDONOW,NDOLINE,NINSAVEOF

        LOGICAL  :: DONOLB

        IF (COPT .EQ. 'I') THEN
C          READ FROM FILE OPENED ON NIN (FOR INTERACTIVE LOOP)
           REWIND(NIN) 
           NDOLINE = 0
        ENDIF

        DONOLB = (LBNOWANT .LT. 0)

c       write(6,*) ' Searching for lb: ',lbnowant

C       READ INPUT STRING FROM PROCEDURE ----------------------------
10      CONTINUE

C       READ ANSWER STRING
        IF (COPT .EQ. 'B') THEN
C          INCREMENT PROCEDURE READ POINTER
	   IBCNT = IBCNT + 1       ! LOCATION IN PROCEDURE FILE
           ILINE = IBCNT 

C          READ FROM CURRENT STORED PROCEDURE LINE IBCNT
           CALL PROC_GETPLINE(IBCNT,0,FCHAR,NCHAR1,IRTFLG)

        ELSE
C          READ FROM FILE OPENED ON NIN (FOR INTERACTIVE LOOP)
           READ(NIN,81,IOSTAT=IRTFLG) FCHAR  ! PGI 10+ NEEDS 81
81         FORMAT(A)

           NCHAR1  = lnblnk(FCHAR)
	   NDOLINE = NDOLINE + 1     ! LOCATION IN DOLOOP FILE
           ILINE   = NDOLINE
           !write(6,*)' Read: ',FCHAR(1:NCHAR1),':',NDOLINE,lbnowant
           !write(6,*) ' read, irtflg:',irtflg,':',fchar(1:nchar1)
        ENDIF
        IF (IRTFLG .NE. 0) RETURN

        IF (DONOLB .AND. ILINE .EQ. -LBNOWANT) THEN
C          FOUND DESIRED DO-LOOP, CAN RETURN NOW
           CALL SSUPCAS(FCHAR(1:NCHAR1))

           IRTFLG = 0
           RETURN

        ELSEIF ( .NOT. DONOLB) THEN
C          SEE IF INPUT LINE CONTAINS 'DO LB' 

C          FIND LENGTH OF LINE BEFOR END OR SEMICOLON (COMMENT)
           IGOSEMI = SCAN(FCHAR,';!')
           IEND    = NCHAR1
           IF (IGOSEMI .GT. 0 .AND. IGOSEMI .LT. IEND) IEND = IGOSEMI
           IF (IEND .LE. 0) GOTO 10   ! JUST A COMMENT

           CALL SSUPCAS(FCHAR(1:IEND))

           IGODO  = INDEX(FCHAR(1:IEND),'DO LB')

           IF (IGODO .GT. 0) THEN
C             HAS 'DO LB', GET LABEL NUMBER FROM 'LB##'
	      CALL GETLBNO(FCHAR(1:NCHAR1),LBNO,IRTFLG)
          
              IF (LBNO .EQ. LBNOWANT) THEN
C                FOUND DESIRED DO-LOOP, CAN RETURN NOW
                 IRTFLG = 0
                 RETURN
              ENDIF
           ENDIF
        ENDIF
        GOTO 10

	END


