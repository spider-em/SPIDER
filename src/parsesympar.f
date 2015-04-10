
C++*********************************************************************
C
C    PARSESYMPAR             
C                 REWRITTEN FROM RDPR             JUNE 2002 ARDEAN LEITH
C                 [] DEFAULT FOR VARIABLES        OCT  2005 ARDEAN LEITH
C                 MPI                             JUN  2009 ARDEAN LEITH
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
C   PARSESYMPAR(PROMPTNID,RESPONSE,PROMPT,NCHARP,SYMPARID,NCHARI,
C               SYMPARVAL,NCHARV,CALLERRT,IRTFLG)
C
C   PURPOSE:   FIND SYMPARID AND SYMPARVAL FROM INPUT STRING(S)
C
C   PARAMETERS: 
C
C   CALL TREE:  DRIV1 
C                v          
C               SYMPAR       
C                v 
C       'FR F'   |-> FILESYMPAR -> ----PARSESYMPAR
C                |                     SETSYMPAR
C                |                                     
C       'FR L'   |-> LOCALSYMPAR   
C                |     v            
C          'FR'  |-> RDPRMC ----->  RDPR -> FRSYMPAR  -> PARSESYMPAR                    ^
C                     ^                                  EVALSYMPAR
C                     ^                                  SETSYMPAR
C      ?..? <ID> -----'                              
C                                                   
C
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

       SUBROUTINE PARSESYMPAR(PROMPTNID,RESPONSE,PROMPT,NCHARP,
     &                       SYMPARID,NCHARI,
     &                       SYMPARVAL,NCHARV,CALLERRT,IRTFLG)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER (LEN=*) ::       PROMPTNID,RESPONSE,PROMPT,SYMPARID
      CHARACTER (LEN=*) ::       SYMPARVAL
      CHARACTER(LEN=1) ::        NULL
      LOGICAL       ::           FROMBACK,CALLERRT

      CALL SET_MPI(ICOMM,MYPID,MPIERR)

      NULL      = CHAR(0)
      IRTFLG    = 1

      SYMPARID  = NULL
      SYMPARVAL = NULL
      PROMPT    = NULL

      IF (PROMPTNID .NE. NULL) THEN
C        SEE IF THERE IS A PROMPT STRING IN PROMPTNID STRING
         CALL CHARINSIDE(PROMPTNID,'?','?',.FALSE.,.FALSE.,
     &                   IPQ1,IPQ2,NCHARP)

         IF (CALLERRT .AND. NCHARP .LT. 3) THEN
            IF (MYPID .LE. 0) THEN
               WRITE(NDAT,*)' *** NO PROMPT (?PROMPT?) IN: ',PROMPTNID
            ENDIF
            CALL ERRT(101,'PARSESYMPAR',NE)
            RETURN
         ENDIF

         IF (NCHARP .GT. 0) PROMPT = PROMPTNID(IPQ1:IPQ2) // NULL

C        SEE IF THERE IS A VARIABLE NAME STRING IN PROMPTNID STRING
C        IF (NCHARP .GT. 0) HAVE ?...?, FIND LAST VARIABLE <LABEL> STRING 
         FROMBACK = (NCHARP .GT. 0)

         CALL CHARINSIDE(PROMPTNID,'[',']',.FALSE.,FROMBACK,
     &                   IP1,IP2,NCHARI)
         IF (CALLERRT .AND. NCHARI .LT. 3) THEN
ccc         IF (MYPID .LE. 0) THEN
ccc            WRITE(NDAT,*)' *** NO VARIABLE NAME ([NAME]) IN: ',
ccc  &                   PROMPTNID
ccc         ENDIF
ccc            CALL ERRT(101,'PARSESYMPAR',NE)
            RETURN
         ENDIF
         IF (NCHARI .GT. 0) THEN
            SYMPARID                = PROMPTNID(IP1:IP2) // NULL
            SYMPARID(1:1)           = '<'
            SYMPARID(NCHARI:NCHARI) = '>'
         ENDIF
      ENDIF


      IF (RESPONSE .NE. NULL) THEN
C        SEE IF THERE IS A VARIABLE ID STRING IN RESPONSE STRING
         CALL CHARINSIDE(RESPONSE,'?','?',.FALSE.,.FALSE.,IPQ,IPQ,NCQ)

C        IF (NCQ .GT. 0) HAVE ?...?, FIND LAST VARIABLE <NAME> STRING 
         FROMBACK = (NCQ .GT. 0)

         CALL CHARINSIDE(RESPONSE,'[',']',.FALSE.,FROMBACK,
     &                   IP1,IP2,NCHARI)
         IF (CALLERRT .AND. NCHARI .LT. 3) THEN
            IF (MYPID .LE. 0) THEN
             WRITE(NDAT,*)' *** NO VARIABLE NAME ([NAME]) IN: ',RESPONSE
            ENDIF
            CALL ERRT(101,'PARSESYMPAR',NE)
            RETURN
         ENDIF
         IF (NCHARI .GT. 0) THEN
            SYMPARID = '<' // RESPONSE(IP1+1:IP2-1) // '>' 
         ENDIF

C        GET VARIABLE VALUE IN RESPONSE STRING
         NCHARV = 0
         IF (IP2 .LT. LEN(RESPONSE)) THEN
            SYMPARVAL = RESPONSE(IP2+1:) // NULL
            NCHARV    = LEN(RESPONSE) - IP2
         ENDIF
      ENDIF

      IRTFLG = 0

      END



