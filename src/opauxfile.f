C **********************************************************************
C
C  OPAUXFILE -- NEW  (MERGED SOME OLD FILES)     FEB 1999 ArDean Leith
C               ADDED SCRATCH                    APR 2001 ArDean Leith
C               FIXED INLINE BUG                 SEP 2001 ArDean Leith
C               LUNSETFLIP                       FEB 2003 ArDean Leith
C               LUNSETLUNS                       FEB 2003 ArDean Leith
C               REMOVED IRTFLG INPUT             APR 2004 ArDean Leith
C               SUPPORT FOR LUN=101              NOV 2006 ArDean Leith
C               IRTFLG = 1 IF NOT EXIST          SEP 2014 ArDean Leith
C               CONVERT LITTLE ENDED OPTION      JAN 2018 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2018  Health Research Inc.,                         *
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
C  OPAUXFILE(ASKNAME,FILNAM,EXTENT,LUNT,LENREC,DISP,PROMPTT,
C            CALLERRT,IRTFLG)
C
C  PURPOSE:   OPENS A NON-SPIDER FILE (CAN HAVE EXTENSION OTHER
C             THAN DATEXC)
C
C PARAMETERS: ASKNAME    LOGICAL FLAG TO REQUEST NAME            (SENT)
C             FILNAM     FILE NAME                           (SENT/RET)
C             EXTENT     FILE EXTENSION (OPTIONAL)               (SENT)
C             LUNT       IO UNIT                                 (SENT)
C                        IF < 0 : FLAG FOR UNFORMATTED, SEQUENTIAL  
C             LENREC     RECORD LENGTH FOR OPEN (BYTES)          (SENT)
C                        >0 : LENGTH FOR UNFORMATTED, DIRECT ACCESS
C                        <0 : LENGTH FOR FORMATTED,   DIRECT ACCESS
C                         0 : FORMATTED, SEQUENTIAL ACCESS
C                         0 & LUNT < 0 : UNFORMATTED, SEQUENTIAL ACCESS 
C             DISP       CHAR FLAG THAT FILE IS OLD, ETC         (SENT)
C                        'O'  OLD (MUST EXIST)
C                        'N'  NEW (WILL BE REPLACED IF EXISTS)
C                        'S'  TEMPORARY SCRATCH FILE
C             PROMPTT    PROMPT FOR FILE NAME (USED IF ASKNAME)  (SENT)
C             CALLERRT   LOGICAL FLAG TO CALL ERRT               (SENT)
C             IRTFLG     ERROR FLAG                              (RET)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C **********************************************************************

        SUBROUTINE OPAUXFILE(ASKNAME,FILNAM,EXTENT,LUNT,LENREC,
     &                       DISP, PROMPTT,CALLERRT,IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'

        INTEGER           :: LUNT,LENREC,IRTFLG
        CHARACTER(LEN=*)  :: FILNAM,EXTENT,PROMPTT,DISP
        LOGICAL           :: CALLERRT,EX,ASKNAME,CONVERT_TO_LITTLE
        CHARACTER(LEN=96) :: PROMPT
        CHARACTER(LEN=80) :: EXTEN
        CHARACTER(LEN=11) :: FORMVAR
        CHARACTER(LEN=10) :: ACCVAR
        CHARACTER(LEN=7)  :: STATVAR

        INTEGER           :: LENP,NCHAR,LNBLNKN,LENE,IRTFLGT,LUN
        INTEGER           :: IDUM,LENOPN,LENOPENFILE
        INTEGER           :: ICOMM, MYPID, MPIERR

        CALL SET_MPI(ICOMM, MYPID, MPIERR) ! SETS ICOMM AND MYPID

C       SET DEFAULT ERROR RETURN
        IRTFLG = 1

C       DO NOT WANT TO RETURN EXTEN
        EXTEN = EXTENT

C       INPUT FILE NAME (IF EXTEN EXISTS IT IS ADDED)

        IF (ASKNAME) THEN
C          SET PROMPT TO ALLOW FILE EXTENSION ON INPUT
           LENP   = LEN(PROMPTT)
           LENP   = MIN(LENP,93)
           PROMPT = PROMPTT(1:LENP) // '~9' 

           CALL FILERD(FILNAM,NCHAR,EXTEN,PROMPT(1:LENP+2),IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ELSE
C          MAY WANT TO ADD EXTENT TO FILNAM
           NCHAR = LNBLNKN(FILNAM)
           LENE  = LNBLNKN(EXTENT)
           IF (LENE > 0) THEN
C             ADD THE EXTENSION THAT IS SENT TO FILNAM
              CALL FILNAMANDEXT(FILNAM,EXTEN,FILNAM,NCHAR,
     &                          .TRUE.,IRTFLGT)
           ENDIF
        ENDIF

        LUN = ABS(LUNT)
        IF ((LUN <= 0 .OR. LUN > 100) .AND. (LUN .NE. 103)) THEN
C          LUN=103 USED IN  SYMPARTEXT 
           CALL ERRT(102,'IN SOURCE CODE, LUN MUST BE 1...100',LUN)
           RETURN
        ENDIF

        IF (LUN > 0 .AND. LUN <= 100) THEN
C          ZERO THE FLAGS USED IN REDLIN/WRTLIN
           CALL LUNSETLUNS(LUN,0,0,LUN,0,IRTFLGT)
 
C          MAKE SURE THIS IS NOT TREATED AS INLINE FILE
           CALL CLOSEINLN(LUN,IRTFLGT)
        ENDIF

C       SET STATUS FOR OPEN
        CONVERT_TO_LITTLE = (DISP(2:2) == 'L')

        STATVAR = 'NEW'
        IF (DISP(1:1) == 'N' .OR. DISP(1:1) == 'U') 
     &     STATVAR = 'REPLACE'

        IF (DISP(1:1) == 'S') STATVAR = 'SCRATCH'

        IF (DISP(1:1) == 'O') THEN
C          CHECK FOR FILE EXISTENCE 
           IF (MYPID <= 0) THEN
              INQUIRE (FILE=FILNAM(1:NCHAR),EXIST=EX,IOSTAT=IRTFLGT) 
           ENDIF

#ifdef USE_MPI
           CALL BCAST_MPI('OPAUXFILE','EX',           EX,1, 'L',ICOMM)
           CALL BCAST_MPI('OPAUXFILE','IRTFLGT', IRTFLGT,1, 'I',ICOMM)
#endif

           IF (IRTFLGT .NE. 0) THEN
              WRITE(NOUT,*) '*** INQUIRY ERROR'
              IF (CALLERRT)  CALL ERRT(4,'OPAUXFILE',IDUM)
              IRTFLG = 1
              RETURN
        
           ELSEIF (.NOT. EX) THEN
              WRITE(NOUT,*) '*** FILE DOES NOT EXIST: ',FILNAM(1:NCHAR)
              IF (CALLERRT)  CALL ERRT(100,'OPAUXFILE',IDUM)
              IRTFLG = 1
              RETURN

           ENDIF
           STATVAR = 'OLD'
        ENDIF

        ACCVAR  = 'DIRECT'
        IF (LENREC == 0) ACCVAR = 'SEQUENTIAL'

        FORMVAR = 'UNFORMATTED'
        IF (LENREC <= 0) FORMVAR = 'FORMATTED'
        IF (LUNT   <  0) FORMVAR = 'UNFORMATTED'

        IF (ACCVAR == 'DIRECT') THEN
C          OPEN FILE FOR DIRECT ACCESS

C          COMPUTE RECL UNITS (DIFFERS WITH OS &A COMPILER FLAGS)
           LENOPN = LENOPENFILE(LENREC)

           IF (MYPID <= 0) THEN
              IF (STATVAR == 'SCRATCH') THEN
	         OPEN(UNIT=LUN,STATUS=STATVAR,
     &               FORM=FORMVAR, ACCESS=ACCVAR, RECL=LENOPN,
     &               IOSTAT=IRTFLGT)

              ELSEIF (CONVERT_TO_LITTLE) THEN
C                FORCE OUTPUT TO LITTLE_ENDIAN

	         OPEN(UNIT=LUN,FILE=FILNAM(1:NCHAR),STATUS=STATVAR,
     &               CONVERT='LITTLE_ENDIAN',
     &               FORM=FORMVAR, ACCESS=ACCVAR, RECL=LENOPN,
     &               IOSTAT=IRTFLGT)

              ELSE
	         OPEN(UNIT=LUN,FILE=FILNAM(1:NCHAR),STATUS=STATVAR,
     &               FORM=FORMVAR, ACCESS=ACCVAR, RECL=LENOPN,
     &               IOSTAT=IRTFLGT)
              ENDIF
           ENDIF
         ELSE
C          OPEN FILE FOR  SEQUENTIAL ACCESS
           IF (MYPID <= 0) THEN
              IF (STATVAR == 'SCRATCH') THEN
	         OPEN(UNIT=LUN,STATUS=STATVAR,
     &               FORM=FORMVAR, ACCESS=ACCVAR, 
     &               IOSTAT=IRTFLGT)
              ELSE
                 OPEN(UNIT=LUN,FILE=FILNAM(1:NCHAR),STATUS=STATVAR,
     &               FORM=FORMVAR, ACCESS=ACCVAR, 
     &               IOSTAT=IRTFLGT)
              ENDIF
           ENDIF
        ENDIF

#ifdef USE_MPI
        CALL BCAST_MPI('OPAUXFILE','IRTFLGT', IRTFLGT,1, 'I',ICOMM)
#endif

        IF (IRTFLGT .NE. 0) THEN
           WRITE(NOUT,90) ACCVAR(1:1),FORMVAR(1:1), FILNAM(:NCHAR)
 90        FORMAT(' ERROR OPENING (',A1,A1,'): ',A)
           IF (CALLERRT) CALL ERRT(102,'OPAUXFILE',IRTFLGT)
           RETURN
        ENDIF

        IF (VERBOSE .AND. MYPID <= 0) THEN
           WRITE(NOUT,91) ACCVAR(1:1),FORMVAR(1:1), FILNAM(:NCHAR)
 91        FORMAT('  OPENED (',A1,A1,'): ',A)
        ENDIF

        IRTFLG = 0

        END
