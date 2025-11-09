
C++*********************************************************************
C
C INTTOCHAR.F   -- NEW 7 Jan 99 ArDean Leith
C                  AVOIDS SGI'S MEMORY LEAK aug 2002  ArDean Leith
C                  DEBUG OUTPUT             sep 2025  ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2025  Health Research Inc.,                         *
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
C    INTTOCHAR(NUMBER,STRING,NLET,MINLEN)
C
C    PURPOSE:       CONVERT INTEGER TO A CHARACTER STRING
C
C    PARAMETERS:
C        NUMBER     INTEGER NUMBER                                (SENT)
C        STRING     OUTPUT STRING                                 (RET.) 
C        NLET       NUMBER OF CHARACTERS IN OUTPUT STRING         (RET.)
C                   (<0 INDICATES ERROR)
C        MINLEN     MINIMUM LENGTH OF OUTPUT STRING               (SENT)
C
C     NOTES: MINLEN AND NLET MUST BE NO MORE THAN 10 DIGITS OR
C            NEGATIVE NLET WILL BE RETURNED 
C            THIS IMPLEMENTATION AVOIDS SGI'S MEMORY LEAK aug 02
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE INTTOCHAR(NUMBER,STRING,NLET,MINLEN)

        INTEGER                :: NUMBER,NLET,MINLEN
        CHARACTER (LEN=*)      :: STRING

        INTEGER, PARAMETER     :: MAXSTR = 10
        CHARACTER (LEN=MAXSTR) :: CTEMP,CZEROS

        INTEGER                :: NZ,IGO


        DATA CZEROS/'0000000000'/

C       FIND NUMBER OF DIGITS TO BE WRITTEN INTO STRING
        NLET = NUMDIG(NUMBER,0)

C       FIND NUMBER OF LEADING ZEROS (IF ANY)
        NZ   = MAX(NLET,MINLEN) - NLET

        IGO  = (MAXSTR-NLET+1)
        NLET = NLET + NZ
        
c        write(6,*) ' Inttochar  number,minlen,nlet,nz: ', 
c     &                          number,minlen,nlet,nz

C       CHECK FOR OVERFLOW OF STRING
        IF (NLET .GT. LEN(STRING) .OR. NLET .GT. MAXSTR) GOTO 999

C       I10 MUST BE SAME AS MAXSTR!!

c       write(3,*)' inttochar maxstr: ',maxstr

        WRITE(CTEMP,'(I10)',ERR=999) NUMBER

        IF (NZ .GT. 0) THEN
C          MUST PREFIX ZERO'S
           STRING = CZEROS(1:NZ) // CTEMP(IGO:MAXSTR) 
        ELSE
           STRING = CTEMP(IGO:MAXSTR) 
        ENDIF

c       write(6,*) ' inttochar number,nlet,string,nlet:',
c       &                      number,nlet,string,nlet

        RETURN

C       ERROR RETURN
999     NLET   = -1
        STRING = CHAR(0)
        RETURN

        END
