C++*********************************************************************
C
C  INQUIREREG.F                      NEW ROUTINE  JAN 1998 ARDEAN LEITH
C                  FOR VARIABLES                  NOV 2005 ARDEAN LEITH
C                  SAYANS = .TRUE.                JUL 2006 ARDEAN LEITH
C                  [i] BUG                        AUG 2009 ARDEAN LEITH
C                  SHOULD NOT ACCEPT SYMVAR       NOV 2009 ARDEAN LEITH
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
C   INQUIREREG(SAYIT,SHUDSTOP,IRTFLG)    
C
C   PURPOSE: DETERMINES IF A REGISTER/REAL VARIABLE'S VALUE IS CORRECT
C
C   PARAMETERS:  SAYIT    LOGICAL FLAG TO LIST ERROR
C                SHUDSTOP LOGICAL FLAG TO STOP IF WRONG IN BATCH
C                IRTFLG   ERROR FLAG
C                         -1 IS ABORT INPUT
C                          0  IS NORMAL
C                          1 INQUIRY WAS NOT AS EXPECTED
C
C--*********************************************************************

      SUBROUTINE INQUIREREG(SAYIT,SHUDSTOP,IRTFLG)

      INCLUDE 'CMBLOCK.INC' 
 
      LOGICAL            :: SAYIT,SHUDSTOP
      CHARACTER(LEN=80)  :: VARNAME,RESPONSE,REGNAME
      CHARACTER(LEN=160) :: MESG
      CHARACTER(LEN=1)   :: CDUM
      REAL, DIMENSION(2) :: FLIST

      LOGICAL            :: GETANS,UPPER,WANTSUB
      LOGICAL            :: SAYPRMT,SAYANS,ENDATSEMI,STRIP

      GETANS    = .TRUE.
      UPPER     = .FALSE.
      WANTSUB   = .FALSE.
      SAYPRMT   = .TRUE.
      SAYANS    = .TRUE.
      ENDATSEMI = .TRUE.
      STRIP     = .TRUE.

      CALL RDPR('REGISTER VARIABLE, & ITS CORRECT VALUE',NCHAR,RESPONSE,.TRUE.,.FALSE.,.FALSE.,
     &       GETANS,UPPER,WANTSUB,SAYPRMT,SAYANS,ENDATSEMI,STRIP,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SEE IF THERE IS A VARIABLE ID STRING IN RESPONSE STRING
      CALL CHARINSIDE(RESPONSE,'[',']',.TRUE.,.FALSE.,
     &                 IP1,IP2,NCHARI)

      IF (NCHARI .GE. 1) THEN
C         MAKE SURE VARIABLE IS A REGISTER VARIABLE NOT A STRING VAR.  
          CALL ISSYMPAR(RESPONSE(IP1-1:IP2+1),-1,ICVAR,IRTFLG) 
          IF (ICVAR .GT. 0) THEN
C            THIS IS A STRING VARIABLE, NOT A REGISTER VARIABLE!
             MESG = RESPONSE(IP1-1:IP2+1) //'IS NOT A REGISTER VARIABLE'
             IF (SAYIT) WRITE(NOUT,92) MESG 
92           FORMAT('  *** ERROR: ',A)
             IF (SHUDSTOP) CALL ERRT(101,MESG,NE)
             IRTFLG = 1
             RETURN
          ENDIF

C         EVALUATE RESPONSE 
          RESPONSE = '~' // RESPONSE(1:)
          CALL RDPRINC(RESPONSE(1:NCHAR+1),2,.FALSE.,NOT_USED,
     &                 REGVALUE,CORRECTVAL,FDUM,IRTFLG)        
          IF (IRTFLG .NE. 0) RETURN

          IF (REGVALUE .NE. CORRECTVAL) THEN
             VARNAME = RESPONSE(2:IP2+2) 

             IF (SAYIT) WRITE(NOUT,90) VARNAME(1:IP2+1),
     &                  REGVALUE,CORRECTVAL 
90           FORMAT(' *** ',A,' =',G11.3, '  SHOULD BE: ',G11.3)

             IF (SHUDSTOP) CALL ERRT(101,'REGISTER VALUE INCORRECT',NE)
             IRTFLG = 1
             RETURN
          ENDIF

       ELSE
C         NOT A NAMED REGISTER VARIABLE, MAYBE OBSOLETE NUMBER ONLY
          RESPONSE(1:NCHAR+1) = '~' // RESPONSE(1:NCHAR)
          CALL RDPRINC(RESPONSE(1:NCHAR+1),2,.FALSE.,NOT_USED,
     &                 REGVALUE,CORRECTVAL,FDUM,IRTFLG)        
          IF (IRTFLG .NE. 0) RETURN

          IREG         = REGVALUE
          REGNAME(1:2) = '[_'
          CALL INTTOCHAR(IREG,REGNAME(3:),NLET,1)
          REGNAME(NLET+3:NLET+3) = ']'
          CALL REG_GET_VAR(0,REGNAME,.FALSE.,REGVALUE,IREG,IEND,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

          IF (REGVALUE .NE. CORRECTVAL) THEN
             LENT = lnblnkn(REGNAME)
             IF (SAYIT) WRITE(NOUT,91) REGNAME(:LENT),
     &                                 REGVALUE,CORRECTVAL 
91           FORMAT(' ***  ',A,' = ',G11.3,'   SHOULD BE: ',G11.3)

             IF (SHUDSTOP) CALL ERRT(101,'REGISTER VALUE INCORRECT',NE)
             IRTFLG = 1
             RETURN
           ENDIF
        ENDIF

        IRTFLG = 0
        WRITE(NOUT,*) ' '

        END






