
C **********************************************************************
C
C OPSTREAMFILE  STREAM IO                        FEB 2013 ArDean Leith
C               CONVERT LITTLE ENDED OPTION      JAN 2018 ArDean Leith
C               ADDED DEBUG OUTPUT               OCT 2025 ArDean Leith
C               ADDED DEBUG OUTPUT               DEC 2025 ArDean Leith
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
C OPSTREAMFILE(ASKNAME,FILNAM,EXTENT,LUNT, FORMVAR, DISP, 
C              PROMPTT,CALLERRT,IRTFLG)
C                                                                      
C PURPOSE: OPENS NON-SPIDER STREAM FILE (CAN HAVE EXTENSION OTHER
C          THAN DATEXC)
C
C NOTES:   IF IRTFLG == 999 ON INPUT DO NOT ECHO FILE OPENING INFO
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

         SUBROUTINE OPSTREAMFILE(ASKNAME,FILNAM,EXTENT,LUNT,
     &                           FORMVAR, DISP, 
     &                           PROMPTT,CALLERRT,IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'  
     
        LOGICAL           :: ASKNAME
        CHARACTER(LEN=*)  :: FILNAM,EXTENT
        INTEGER           :: LUNT
        CHARACTER(LEN=11) :: FORMVAR
        CHARACTER(LEN=*)  :: DISP,PROMPTT
        LOGICAL           :: CALLERRT
        INTEGER           :: IRTFLG

        CHARACTER(LEN=MAXNAM)  :: exten4
        CHARACTER(LEN=MAXNAM)  :: trfilnam
        CHARACTER(LEN=1)       :: c22,c23

        LOGICAL           :: EX, CONVERT_TO_LITTLE,CONVERT_TO_BIG
        LOGICAL           :: SAYIT
        CHARACTER(LEN=96) :: PROMPT
        CHARACTER(LEN=80) :: EXTEN
        CHARACTER(LEN=7)  :: STATVAR

        INTEGER           :: ICOMM,MYPID,MPIERR,LENP,NCHAR
        INTEGER           :: LENE,IRTFLGT,LUN,IDUM,LENOPEN,LENOPENFILE
        INTEGER           :: LENOPN,LENREC

        INTEGER           :: lnblnkn  ! FUNCTIONS  

        INTEGER           :: MY_IOSTAT
        CHARACTER(256)    :: MY_IOMSG


        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

        SAYIT = (IRTFLG .NE. 999)

C       SET DEFAULT ERROR RETURN
        IRTFLG = 1

C       DO NOT WANT TO RETURN EXTEN
        EXTEN = EXTENT

C       INPUT FILE NAME (IF EXTEN EXISTS IT IS ADDED)

#if defined(SP_DBUGIO)
        write(3,*)' In opstreamfile 0; askname,lunt: ', askname,lunt
        write(3,*)' In opstreamfile 0; filnam:  |', trim(filnam),'|'
#endif

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

#if defined(SP_DBUGIO)
        write(3,*)' In opstreamfile 1; nchar:  ', nchar
        write(3,*)' In opstreamfile 1; lene:   ', lene
        write(3,*)' In opstreamfile 1; exten:  |', trim(exten), '|'
        write(3,*)' In opstreamfile 1; filnam: |', trim(filnam),'|'
#endif

           IF (LENE > 0) THEN
C             ADD THE EXTENSION THAT IS SENT TO FILNAM
              CALL FILNAMANDEXT(FILNAM,EXTEN,FILNAM,NCHAR,
     &                          .TRUE.,IRTFLGT)
              NCHAR = LNBLNKN(FILNAM)
           ENDIF
        ENDIF


        LUN = ABS(LUNT)
        IF ((LUN <= 0 .OR. LUN > 100) .AND.
     &      (LUN .NE. 103)) THEN
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

        CONVERT_TO_BIG    = (DISP(2:2) == 'B')
        CONVERT_TO_LITTLE = (DISP(2:2) == 'L')
        STATVAR           = 'NEW'

        IF (DISP(1:1) == 'N' .OR. DISP(1:1) == 'U') 
     &     STATVAR = 'REPLACE'

        IF (DISP(1:1) == 'S') STATVAR = 'SCRATCH'
        
#if defined(SP_DBUGIO)
        write(3,*)' In opstreamfile 3; nchar:   ', nchar
        write(3,*)' In opstreamfile 3; tr filnam:|', trim(filnam),'|' 
        write(3,*)' In opstreamfile 3; exten:    |', trim(exten),'|'
        write(3,*)' In opstreamfile 3; disp,mypid :', disp,mypid
        write(3,*)' In opstreamfile 3; to_big,to_little:',
     &                      convert_to_big,',',convert_to_little
        write(3,*)' '
#endif


        IF (DISP(1:1) == 'O') THEN

C          CHECK FOR FILE EXISTENCE 
           IF (MYPID <= 0) THEN
              INQUIRE (FILE=FILNAM(1:NCHAR),EXIST=EX,IOSTAT=IRTFLGT) 
           ENDIF

#ifdef USE_MPI
           CALL BCAST_MPI('OPSTREAMFILE','EX',          EX,1,'L',ICOMM)
           CALL BCAST_MPI('OPSTREAMFILE','IRTFLGT',IRTFLGT,1,'I',ICOMM)
#endif

           IF (IRTFLGT .NE. 0) THEN
              WRITE(NOUT,*) '*** INQUIRY ERROR'
              WRITE (NOUT,*)' ERROR IOSTAT: ', IRTFLGT, 
     &                      ' IOMSG: '//TRIM(MY_IOMSG)
              IF (CALLERRT) CALL ERRT(4,'OPSTREAMFILE',IDUM)
              RETURN
        
           ELSEIF (.NOT. EX) THEN

              WRITE(NOUT,*)'*** FILE DOES NOT EXIST: ',TRIM(FILNAM)
              !WRITE (NOUT,*)' ERROR IOMSG: '// TRIM(MY_IOMSG)
              IF (CALLERRT)  CALL ERRT(100,'OPSTREAMFILE',IDUM)
              RETURN

           ENDIF
           STATVAR = 'OLD'  
        ENDIF

C       OPEN FILE FOR STREAM ACCESS


C       COMPUTE RECL UNITS (DIFFERS WITH OS & COMPILER FLAGS)
        LENOPN = LENOPENFILE(LENREC)

        IF (MYPID <= 0) THEN

           IF (STATVAR == 'SCRATCH') THEN   !----- SCRATCH FILE ----

              IF (CONVERT_TO_LITTLE) THEN
C                FORCE OUTPUT TO LITTLE_ENDIAN

                 OPEN(UNIT=LUN,STATUS=STATVAR,
     &             CONVERT='LITTLE_ENDIAN',
     &             FORM=FORMVAR, ACCESS='STREAM',
     &             IOSTAT=IRTFLGT,IOMSG=MY_IOMSG)

              ELSEIF (CONVERT_TO_BIG) THEN
C                FORCE OUTPUT TO BIG_ENDIAN

                 OPEN(UNIT=LUN,STATUS=STATVAR,
     &             CONVERT='BIG_ENDIAN',
     &             FORM=FORMVAR, ACCESS='STREAM',
     &             IOSTAT=IRTFLGT,IOMSG=MY_IOMSG)

              ELSE

                 OPEN(UNIT=LUN,STATUS=STATVAR,
     &             FORM=FORMVAR, ACCESS='STREAM',
     &             IOSTAT=IRTFLGT,IOMSG=MY_IOMSG)
              ENDIF

           ELSE    ! ----------- NOT A SCRATCH FILE -------------

              IF (CONVERT_TO_LITTLE) THEN
C               FORCE OUTPUT TO LITTLE_ENDIAN

                OPEN(UNIT=LUN,FILE=FILNAM(1:NCHAR),STATUS=STATVAR,
     &             CONVERT='LITTLE_ENDIAN',
     &             FORM=FORMVAR, ACCESS='STREAM', 
     &             IOSTAT=IRTFLGT,IOMSG=MY_IOMSG)

              ELSEIF (CONVERT_TO_BIG) THEN
C                FORCE OUTPUT TO BIG_ENDIAN

                 OPEN(UNIT=LUN,FILE=FILNAM(1:NCHAR),STATUS=STATVAR,
     &             CONVERT='BIG_ENDIAN',
     &             FORM=FORMVAR, ACCESS='STREAM',
     &             IOSTAT=IRTFLGT,IOMSG=MY_IOMSG)

              ELSE
                 OPEN(UNIT=LUN,FILE=FILNAM(1:NCHAR),STATUS=STATVAR,
     &             FORM=FORMVAR, ACCESS='STREAM', 
     &             IOSTAT=IRTFLGT,IOMSG=MY_IOMSG)

             ENDIF
           ENDIF
        ENDIF

#ifdef USE_MPI
        CALL BCAST_MPI('OPSTREAMFILE','IRTFLGT', IRTFLGT,1, 'I',ICOMM)
#endif

        IF (IRTFLGT .NE. 0) THEN

#if defined(SP_DBUGIO)
           write (3,*) ' Open failed iostat: ', irtflgt, 
     &                 ' iomsg = '//trim(my_iomsg)
#endif
 
           WRITE(NOUT,90) FORMVAR(1:1), TRIM(FILNAM)
 90        FORMAT(' ERROR OPENING (',A1,'): ',A)
           WRITE (NOUT,*) ' OPEN FAILED IOMSG: '//TRIM(MY_IOMSG)
           IF (CALLERRT) CALL ERRT(102,'OPSTREAMFILE',IRTFLGT)
           RETURN
        ENDIF

        IF (SAYIT .AND. VERBOSE .AND. MYPID <= 0) THEN
           WRITE(NOUT,91) FORMVAR(1:1), FILNAM(:NCHAR)
 91        FORMAT('  OPENED  (',A1,'): ',A)
        ENDIF
 
        IRTFLG = 0

#if defined(SP_DBUGIO)
        write(3,*)' Leaving opstreamfile; irtflg: ',irtflg
        write(3,*)
#endif

        END


