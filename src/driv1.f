
C++*********************************************************************
C
C    DRIV1
C            CALLS ROUTINES REMOVED FROM DRIVER   MAR 93
C            CHANGED READRQ PARAMETERS PASSED     AUG 99   ARDEAN LEITH
C            'FR' MOVED TO SPIDER                 SEP 2000 ARDEAN LEITH
C            'PO' FOR POLAR CONVERSION            SEP 2000
C            'FR L' ADDED                         MAR 2001 ARDEAN LEITH
C            SIMPLIFED WITH SETSYMPAR IN RDPR     APR 2001 ARDEAN LEITH
C            SYMPAR                               JUN 2002 ARDEAN LEITH
C            MULTILINE 'VMS'                      SEP 2003 ARDEAN LEITH
C            'MD' TO SPIDER.f                     DEC 2003 ARDEAN LEITH
C            'VO IA' & 'VO EPT'                   AUG 2004 ARDEAN LEITH
C            'PI IMG'                             JAN 2006 ARDEAN LEITH
C            'VM CD'                              JAN 2009 ARDEAN LEITH
C            'TM' REGISTER                        SEP 2010 ARDEAN LEITH
C            'TM' NOT RESET TO ZERO BUG           JUN 2011 ARDEAN LEITH
C            'VO TA' REMOVED                      OCT 2012 ARDEAN LEITH
C            'PO R'                               OCT 2013 ARDEAN LEITH
C            'NC' SHOULD NOT CLOSE NDAT           JAN 2014 ARDEAN LEITH
C            FINDRIDGES                           APR 2014 ARDEAN LEITH
C            FINDRIDGES(RIDGESONLY)               MAR 2016 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2016  Health Research Inc.,                         *
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
C   DRIV1(MAXDIM)
C
C   PURPOSE:   CALLS ROUTINES REMOVED FROM DRIVER 
C
C   PARAMETERS: MAXDIM     MAX LENGTH OF COMMON BUFFER
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE DRIV1(MAXDIM)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER (LEN=1)       :: SWITCH
      CHARACTER (LEN=4)       :: CXNUM
      COMMON /DRIV1_COM/ T1,LOOPREG,CXNUM

      CHARACTER(LEN=2*MAXNAM) :: RESPONSE,PROMPT
      CHARACTER(LEN=7)        :: EXTEN
      INTEGER                 :: HRS,MIN,SEC
      LOGICAL                 :: MULTILINE,RIDGESONLY
      INTEGER                 :: ICOMM,MYPID,MPIERR

      CHARACTER(LEN=1)        :: NULL = CHAR(0)

      INTEGER, PARAMETER      :: LUNDOC = 77

c     MENU/'NC','VM','ME','CK','TM',
c     &    'SR','RR','FR','PO','SA',
c     &    'VO','EV','PI','PB','RI'/   -- IS FOR INITILIZING EXTEN

      CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYMPID


      SELECT CASE(FCHAR(1:2))

      CASE('NC','--') !  NEW FILE EXTENSION ----------------------- NC
         !!!!CLOSE(NDAT) jan 2014 al

C        GET THE PROJECT AND DATA EXTENSION
100      IRTFLG = -999

         CALL RDPRMC(EXTEN,NC,.TRUE.,'PROJECT/DATA EXTENSION',
     &             NULL,IRTFLG)

C        MAKE SURE PROJECT EXTENSION IS VALID
         IF (NC < 3 .OR. IRTFLG .NE. 0) THEN
            IF (MYPID .LE. 0)
     &         WRITE(NOUT,*) ' *** EXTENSION MUST BE 3 CHARACTERS'
            GOTO 100
         ENDIF

         PRJEXC(1:3) = EXTEN(1:3)

         IF (EXTEN(4:4) .NE. '/') THEN
            DATEXC(1:3) = PRJEXC(1:3)
         ELSE
            DATEXC(1:3) = EXTEN(5:7)
         ENDIF


       CASE('VM')   !  VMS SPAWN COMMANDS SYSTEM CALL -------------- VM
          IF (FCHAR(4:4) == 'C') THEN
             CALL VMS_CD()
          ELSE
             MULTILINE =  (FCHAR(4:4) == 'M') 
             CALL VMS(MULTILINE)
          ENDIF

       CASE('ME')   !  MENU LISTING -------------------------------- ME
          WRITE(NOUT,*)
     &         ' OPERATION NO LONGER SUPPORTED, USE MANUAL INSTEAD'
          CALL ERRT(100,'DRIV1',NE)


       CASE('CK')   !  CHECKPOINT OPERATION ------------------------ CK
C        WRITE OUT CURRENT OPERATION AND CURRENT ITERATION NUMBER.
         IF (NTRACE.NE.0) THEN
            WRITE(NOUT,9020) FCHAR(1:NALPH)
9020        FORMAT(' ',A)
	    IF (IABSLP .NE. 0 .AND. LOOPREG .NE. 0) THEN
               IF (LOOPREG < 103) THEN 
                  WRITE(NOUT,9040) LOOPREG,IABSLP
9040              FORMAT('  LOOP COUNTER (',I3,') = ',I5)
               ELSE
                  WRITE(NOUT,9041) CHAR(LOOPREG-103),IABSLP
9041              FORMAT('  LOOP INDEX (',A,') = ',I5)
               ENDIF
             ENDIF
          ENDIF


       CASE('TM')   !  TIME OPERATION ------------------------------ TM
C        GET NUMBER OF SEC. SINCE LAST TM, AND COMPUTE HOURS,MIN, & SEC
         TIM  = SECNDS(T1)
	 HRS  = INT(TIM/3600.)
	 HOUR = HRS*3600.
	 MIN  = INT((TIM-HOUR)/60.)
	 SEC  = TIM-HOUR-MIN*60.

	 WRITE(NDAT,9240)HRS,MIN,SEC
9240	 FORMAT(' TIME: ',I2,' HOURS  ',I2,' MINUTES  ',I2,' SECONDS')
       
         CALL REG_SET_NSEL(1,1,TIM,0.0,0.0,0.0,0.0,IRTFLG)
         T1 = SECNDS(0.0)   ! RESET TO ZERO


       CASE('SR')   !  SAVE REGISTERS ------------------------------ SR
          CALL RDPRMC(SWITCH,NLET,.TRUE.,'(S)AVE OR (U)NSAVE',NULL,IRT)
	  IF (IRT .NE. 0) GOTO 9999
           
          CALL ERRT(101,'OBSOLETE OPERATION NO LONGER SUPPORTED',NE)
          GOTO 9999


       CASE('RR')   !  READ REGISTER ------------------------------- RR
          CALL READRQ()
          GOTO 9999


       CASE('FR')   !  FILE READ ----------------------------------- FR
          CALL SYMPAR(LUNDOC)
          GOTO 9999
 
 
       CASE('RI')   !   RIDGE LOCATIONS ---------------------------- RI

          IF (FCHAR(4:5) == 'RV') THEN 
C            FIND RIDGES & VALLEYS 
             RIDGESONLY = .FALSE.
             CALL FINDRIDGES(RIDGESONLY)

          ELSE  

C            FIND VERTICAL RIDGES 
             RIDGESONLY = .TRUE.
             CALL FINDRIDGES(RIDGESONLY)

          ENDIF


       CASE('PO')   !   POLAR CONVERSION --------------------------- PO

          IF     (FCHAR(4:4) == 'R') THEN
C            CONVERT TO POLAR REPRESENTATION - RAYS ALONG X DIMENSION 
             CALL TO_RAYS()

          ELSEIF (FCHAR(4:4) == 'P') THEN   
             CALL TO_PEAKS()

          ELSE
C            CONVERT TO POLAR REPRESENTATION - RAYS ALONG Y DIMENSION 
             CALL TO_POLAR()
          ENDIF


       CASE('SA')   !    SUM ALIGNMENTS ---------------------------- SA
          IF (FCHAR(4:4) == 'P')    THEN
            CALL SUMALI(.TRUE.)

          ELSEIF(FCHAR(4:4) == '3') THEN
            CALL SUMALI3

          ELSEIF(FCHAR(4:4) == 'E') THEN
            CALL SUMEULER

          ELSE
            CALL SUMALI(.FALSE.)

          ENDIF


       CASE('VO')   !    ------------------------------------------- VO

         IF (FCHAR(4:5) == 'EA')      THEN
           CALL VOEA(MAXDIM)

         ELSEIF (FCHAR(4:6) == 'NEA') THEN
           CALL VONEA()

         ELSEIF (FCHAR(4:6) == 'RAS') THEN
           CALL VORA()

         ELSEIF (FCHAR(4:5) == 'RA')  THEN
           CALL VORA()

         ELSEIF (FCHAR(4:5) == 'MD')  THEN
           CALL VOMD

         ELSEIF (FCHAR(4:5) == 'MQ')  THEN
           CALL VOMQ

         ELSEIF (FCHAR(4:5) == 'EP')  THEN
           CALL VOEPT

         ELSEIF (FCHAR(4:5) == 'IA')  THEN
           CALL VOIA

         ELSEIF (FCHAR(4:5) == 'DA')  THEN
           WRITE(NOUT,*) ' UNDOCUMENTED OPERATION'
           CALL VODA

         ELSEIF (FCHAR(4:5) == 'TA')  THEN
           !CALL VOTA(MAXDIM) oct 2012 al
           CALL ERRT(101,
     &          'UNDOCUMENTED OPERATION, NO LONGER SUPPORTED',NE)
         ENDIF


       CASE('EV')   !    ---SET ENVIRONMENTAL VARIABLE ------------ EV
         IRTFLG = -999
C        IRTFLG = -999 DOES NOT CONVERT INPUT TO UPPERCASE

         CALL RDPRMC(PROMPT,NLET1,.TRUE.,'ENVIRONMENT VARIABLE',
     &               NULL,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999
         IRTFLG = -999
         CALL RDPRMC(RESPONSE,NLET2,.TRUE.,'VALUE',NULL,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999

         CALL SETENV(PROMPT(1:NLET1),RESPONSE(1:NLET2),IRTFLG)


       CASE('PI')   ! --- SEND REGISTER DOWN PIPE ------------------ PI

         IF (FCHAR(4:4) == 'R') THEN
C           SEND REGISTER DOWN PIPE  --------------------------- PI REG
            IRTFLG = -999
            CALL RDPRMC(RESPONSE,NCHAR,.TRUE.,'REGISTER VARIABLE',
     &              CDUM,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999

            CALL REG_PIPE(RESPONSE,IRTFLG)
         ELSE
C          SEND IMAGE DOWN PIPE  ------------------------------- PI IMG
           CALL WRTLIN_PIPE_TOG()
         ENDIF


       CASE('PB')   ! ---  MANIPULATE PDB FILES -------------------  PB
         CALL  PDB

       END SELECT

9999   RETURN
       END




