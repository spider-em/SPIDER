
C++*********************************************************************
C
C  EVALSYMPAR.F -- CREATED 6/8/02 ARDEAN LEITH 
C                 [] DEFAULT FOR VARIABLES         OCT 2005 ARDEAN LEITH
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
C EVALSYMPAR(SYMPARIN,SYMPAROUT,NCHARV,IRTFLG)
C
C PURPOSE: SUBSTITUTE FOR <>,X,I IN SYMPARVAL FROM HIGHER LEVELS
C             
C PARAMETERS:     SYMPARIN     SYMBOLIC PARAMETER VALUE        SENT
C                 SYMPAROUT    SYMBOLIC PARAMETER VALUE        RETURNED
C                 NCHARV       LENGTH OF SYMPAROUT             RETURNED
C                 IRTFLG       ERROR FLAG (0 IS NORMAL)        RETURNED
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE EVALSYMPAR(SYMPARIN,SYMPAROUT,NCHARV,IRTFLG)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=*) ::        SYMPARIN,SYMPAROUT
      CHARACTER(LEN=MAXNAM) ::   SYMPARINT

C     FOR LOCAL SYMBOLIC PARAMETER HANDLING 
      INTEGER, DIMENSION(MAXPRC) :: IPSTACK,IPNUMSTACK,IPARNUM
      COMMON /QSTR_STUFF1/ ISTOP,ITI,ITIN,IWHERE,IPSTACK,
     &                     IPNUMSTACK,IPARNUM

      SYMPAROUT = SYMPARIN
      NCHARV    = LEN(SYMPARIN)

C     SYMPARIN MAY CONTAIN PARAMETER SUBSTITUTION FROM HIGHER PROCS.
C     IF SO, SUBSTITUTE THE VALUES FROM THE HIGHER LEVEL PROCS. NOW
      IF (ISTOP .GT. 1) THEN
         SYMPARINT = SYMPARIN

         DO ILEVEL = ISTOP - 1,1,-1
           SYMPARINT = SYMPAROUT

C          SEE IF ANY VARIABLE STRING ARRIVED IN SYMPARIN 
           IF (INDEX(SYMPARINT(1:NCHARV),'<') .GT. 0) THEN
C             CONVERT OLD STYLE VARIABLE DELIMITERS TO []
              DO I = 1,NCHARV
                 IF (SYMPARINT(I:I) .EQ. '<') SYMPARINT(I:I) = '['
                 IF (SYMPARINT(I:I) .EQ. '>') SYMPARINT(I:I) = ']'
              ENDDO
           ENDIF
C          MAY NEED VARIABLE SUBSTITUTION
           CALL SUBSYMPAR(SYMPARINT,SYMPAROUT,NCHARV,ILEVEL,IRTFLG)

C          SEE IF NEED TO CONVERT OLD x11 REGISTER FORMAT 
           IX = SCAN(SYMPAROUT(1:NCHARV),'xX')
           IF (IX .GT. 0) THEN
C             CONVERT OLD x11 REGISTER FORMAT TO TO NEW: [name] FORMAT
              CALL DEXREG(SYMPAROUT,NCHARV)
           ENDIF

C          MAY WANT TO SUBSTITUTE. FOR REGS? FROM SYMPARINT?
           ISUB = SCAN(SYMPAROUT(1:NCHARV), '{[*')
           IF (ISUB .GT. 0) THEN
C             SUBSTITUTE FOR: {***[]}  {---[]}  ***[]  ${ENV}  .1[] 
C             PASS THE CURRENT REGISTER SET TO FILNAMSUB
C             WILL STOP IN ERRT IN FILNAMSUB IF THERE IS AN ERROR

              CALL FILNAMSUB(SYMPAROUT,NCHARV,-1,IRTFLG)
           ENDIF
           IF (INDEX(SYMPARINT,'<') .LE. 0) EXIT
         ENDDO
      ENDIF
    
      IRTFLG = 0

      RETURN
      END

