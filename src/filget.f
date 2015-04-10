C++*********************************************************************
C
C    FILGET.FOR  - CREATED JAN 89
C                  SUBSTITUTES ANYWHERE IN STRING JULY 99 ARDEAN LEITH  
C                  NO LONGER ALTERS FILPAT FOR {***?} JAN 2001 A. LEITH
C                  FILNAMSUB PARAMETERS CHANGED       DEC 2005 A. LEITH
C                  XMIPP SELFILE                      DEC 2010 A. LEITH
C                  NAME@* SUB                         AUG 2011 A. LEITH
C **********************************************************************
C=* AUTHOR: A. LEITH
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
C    FILGET(FILPAT,FILNAM,NLET,INUM,IRTFLG)
C
C    PURPOSE:        USED AFTER FILSEQ TO CREATE NEXT FILENAME FROM 
C                    FILE-NAME TEMPLATE AND ARRAY OF FILE NUMBERS.
C                    ANALGOUS TO FILNUM.FOR
C
C    PARAMETERS:     FILPAT    CHAR. FILE NAME PATTERN            (SENT)
C                    FILNAM    CHAR. FILE NAME                    (RET.)
C                    NLETT     NUMBER OF LETTERS IN FILE NAME     (SENT)
C                    INUM      FILE NUMBER                        (SENT)
C                                 <0 IS Xmipp selfile             (SENT)
C                    IRTFLG    ERROR FLAG (0 IS NORMAL)           (RET.)
C
C--*********************************************************************

	SUBROUTINE FILGET(FILPAT,FILNAM,NLETT,INUM,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER (LEN=*)      :: FILPAT,FILNAM
        INTEGER                :: NLETT,INUM,IRTFLG

        CHARACTER (LEN=8)      :: FMT = '(I  .  )'
        LOGICAL                :: ISDIGI, HASAST
        CHARACTER (LEN=MAXNAM) :: RECLIN

        INTEGER                :: lnblnk,numdig
        INTEGER                :: NLET,LUNXM,LOCAST,LENAST,IGO,IEND
        INTEGER                :: I,NUMI,NGOT
        INTEGER                :: LOCAT,LASTPOS,LENMAX

        !DATA FMT  /'(I  .  )'/

        NLET   = NLETT
        FILNAM = FILPAT
        IF (NLET .LE. 0) NLET = lnblnk(FILNAM)

        IF (INUM < 0) THEN
C          USING XMIPP SELFILE, -INUM IS: INPUT UNIT = LUNXM
           LUNXM = ABS(INUM)
           CALL GETNEXT_XMSEL(LUNXM,.TRUE.,FILNAM,NLET,IRTFLG)
           write(6,*) ' Filget, xmipp filnam: ',filnam(:nlet)
           RETURN
        ENDIF

C       SUBSTITUTE FOR STRINGS WITHIN {...} PORTION(S) OF THE FILE NAME
        CALL FILNAMSUB(FILNAM,NLET,0,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       FIND LOCATION OF FIRST * IN FILNAM
        LOCAST = INDEX(FILNAM,'*')
        LOCAT  = INDEX(FILNAM,'@')

        IF (LOCAT .EQ. NLET .OR. 
     &      (LOCAT > 0 .AND. LOCAST > LOCAT)) THEN
C          TERMINAL @, APPEND THE STACKED IMAGE NUMBER
           NUMI    = NUMDIG(INUM,0)   ! NUMBER OF DIGITS IN: INUM
           LOCAST  = LOCAT + 1
          
C          CHECK FOR FILE NAME LENGTH OVERFLOW
           LASTPOS = LOCAT + NUMI
           LENMAX  = LEN(FILNAM)
           IF (LASTPOS > LENMAX) THEN
              WRITE(NOUT,92) INUM,FILPAT(1:NLET),LENMAX
92            FORMAT(' *** APPENDING FILE NUMBER:',I6,
     &               ' TO: ',A,' OVERFLOWS NAME LENGTH:',I4)
              IRTFLG = 1
              RETURN
           ENDIF
           !write(6,*) '  filgetat, filnam: ',filnam(1:nlet)

           CALL INTTOCHAR(INUM,FILNAM(LOCAST:LENMAX),NGOT,NUMI)
           IRTFLG = 0
           IF (NGOT .LE. 0) THEN
              WRITE(NOUT,*) ' *** ERROR IN FILENAME'
              IRTFLG = 1
           ENDIF
           !write(6,*) ' filgetat, filnam: ',filnam(:nlet)
           RETURN 


        ELSEIF (LOCAST .LE. 0) THEN
C          NO *, ACCEPT IT ANYWAY IF FILE HAS TERMINAL DIGITS
C          JUST REPLACE THE TERMINAL DIGITS
           LENAST = 0
           I      = NLET
           DO WHILE (I .GT. 0 .AND. ISDIGI(FILNAM(I:I))) 
              LENAST = LENAST + 1
              LOCAST = I
              I      = I - 1
           ENDDO 

C          FIND NUMBER OF DIGITS IN FILE NUMBER. 
C          DO NOT MAKE IT EQUAL LENAST
           NUMI = NUMDIG(INUM,1)

        ELSE
C          FIND NUMBER OF * IN FILNAM
           LENAST = 1
           I      = LOCAST + 1
           DO WHILE (I .LE. NLET .AND. FILNAM(I:I) .EQ. '*') 
              I      = I + 1
              LENAST = LENAST + 1
           ENDDO 

C          FIND NUMBER OF DIGITS IN FILE NUMBER. IF LESS THAN LENAST
C          MAKE IT EQUAL LENAST

           NUMI = NUMDIG(INUM,LENAST)
        ENDIF

C       FIND 1ST CHAR IN FILNAM THAT WILL BE ALTERED
        IGO = LOCAST + LENAST - NUMI

        IF (IGO .LT. 1) THEN
           WRITE(NOUT,90) INUM,FILPAT(1:NLET)
90         FORMAT(' *** FILE NUMBER:',I6,' DESTROYS FILE PATTERN: ',A)
           IRTFLG = 1
           RETURN

        ELSEIF (IGO .LT. LOCAST) THEN
          WRITE(NOUT,91) INUM

91        FORMAT(' WARNING, FILE NUMBER:',I8,' MAY ALTER FILE PATTERN.')
        ENDIF

        IEND = IGO + MAX(LENAST,NUMI) - 1

        CALL INTTOCHAR(INUM,FILNAM(IGO:IEND),NGOT,NUMI)
        IRTFLG = 0
        IF (NGOT .LE. 0) THEN
           WRITE(NOUT,*) ' *** ERROR IN FILENAME'
           IRTFLG = 1
        ENDIF
        !write(6,*) ' filget, filnam: ',filnam(:nlet)

        END


