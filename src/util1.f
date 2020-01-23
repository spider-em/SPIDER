C
C++*********************************************************************
C
C    UTIL1.F   TITLE PROCESSING MODIFIED          NOV  87 ArDean Leith
C              LONG FILE NAMES ADDED              DEC  88 ArDean Leith
C	       ALTERED                            4/21/93 Mahieddine Ladjadj
C	       ALTERED                            8/25/93 Jing Su
C	       LI COMMAND REWRITTEN               8/30/96 ArDean Leith
C	       TT COMMAND REWRITTEN               6/28/97 ArDean Leith
C              ST PARAMETERS ALTERED              SEP  98 ArDean Leith
C              'TT COP' ADDED                     JUN  99 ArDean Leith
C              'CA SM' ALTERED                    AUG  99 ArDean Leith
C              'TF CTS' ALTERED                   NOV  00 Haixiao Gao
C              'TF CTF' REMOVED                   JAN  01 ArDean Leith
C              'TF CRF' ADDED                     JAN  11 Paul Penczek
C              'TF ECTF' ADDED                    JUL  31 Paul Penczek
C              'TF ECTF' --> 'TF ED'              JUN  02 Bill Baxter
C              'HI E' ADDED                       FEB  03 ArDean Leith
C              OPFILEC                            FEB  03 ArDean Leith
C              'HI J' ADDED                       MAR  03 ArDean Leith
C              'PK DR' ADDED                      MAR  03 Bimal Rath
C              'CA' REWRITE                       SEP  03 ArDean Leith
C              'TF EA' REMOVED                    NOV  03 Paul Penczek
C              'TF ED' REPLACED                   NOV  03 Paul Penczek
C              'CA SMI' ADDED                     JAN  04 ArDean Leith
C              'HI J' VMIN, VMAX                  FEB  04 ArDean Leith
C               TRAFC & TRAFCT MERGED             MAR  04 ArDean Leith
C               ~7 REPLACES IRTFLG                APR  04 ArDean Leith
C               'PK 3R'                           NOV  04 ArDean Leith
C               'TF SIM' ADDED                    NOV  07 Bimal  Rath 
C               'TF LM4'ADDED                     MAR  06 Zhong  Huang 
C               'HI DOC' DUPLICATES HI D          MAR  06 ArDean Leith 
C               'FI H' ON STACK HEADER            OCT  10 ArDean Leith 
C               'TF COR'                          NOV  10 ArDean Leith 
C               'FI H' INQUIREHEAD ARGS           NOV  10 ArDean Leith 
C               'ST H'                            NOV  10 ArDean Leith 
C	        'FI H' NO FILE BUG                JAN  11 ArDean Leith
C	        CASE                              JAN  11 ArDean Leith
C	        'CTF FIND'                        MAY  12 ArDean Leith
C	        'CTF ED'                          JUN  12 ArDean Leith
C	        'CA EIGPCT'                       JAN  13 ArDean Leith
C	        'FI' PROMPT SHORTENED             JUN  13 ArDean Leith
C               'TF RCTF' REMOVED                 JAN  14 ArDean Leith
C               'LI T x' ADDED                    AUG  14 ArDean Leith
C               'ST H' KLUDGE FOR BAD EM2EM STACK AUG  14 ArDean Leith
C               CALL NOISE PARAMETERS             NOV  14 ArDean Leith
C	        'FI CEN' ADDED                    DEC  14 ArDean Leith
C	        'TF C?' MERGERS                   NOV  15 ArDean Leith
C	        'FI H' CAN OPEN STACK WITHOUT @   FEB  16 ArDean Leith
C	        'TT' NOT IN MRC FILE              JUL  19 ArDean Leith
C	         MRC FILE SUPPORT                 OCT  19 ArDean Leith
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2019  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C   UTIL1(MAXDIM,IRTRET)
C
C   PURPOSE:    ORIGINALLY, DRIVER FOR ROUTINES REQUIRING A SINGLE FILE
C
C   PARAMETERS: MAXDIM     MAX LENGTH OF COMMON BUFFER
C               IRTRET     ERROR RETURN FLAG
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE UTIL1(MAXDIM,IRTRET)

	INCLUDE 'CMBLOCK.INC' 
	INCLUDE 'CMLIMIT.INC' 

        INTEGER                    :: MAXDIM,IRTRET

        LOGICAL                    :: IS_MRC
        INTEGER                    :: MAXIM,MAXIM2,IDUM

        INTEGER                    :: ICOMM,MYPID,MPIERR 
	INTEGER, PARAMETER         :: NFUNC = 19
        CHARACTER(LEN=2)           :: FUNC(NFUNC)

        LOGICAL                    :: DOCPRNT,TERMPRNT
        CHARACTER(LEN=MAXNAM)      :: FILNAM,CLINE
        CHARACTER(LEN=1)           :: NULL = CHAR(0)
        CHARACTER(LEN=1)           :: DISP

        INTEGER                    :: NX,NY,NZ,NSEL_USED
        INTEGER                    :: IXC,IYC,IZC

 	INTEGER, PARAMETER         :: LUN1   = 8
	INTEGER, PARAMETER         :: LUN2   = 12
	INTEGER, PARAMETER         :: LUN3   = 7
	INTEGER, PARAMETER         :: LUN4   = 9
	INTEGER, PARAMETER         :: LUN5   = 13
	INTEGER, PARAMETER         :: LUNDOC = 80
	INTEGER, PARAMETER         :: LUNXM  = 81
       
C       DATA FUNC/'DE', 'DU', 'FI', 'HI', 'LI', 
C     &           'MO', 'PK', 'RA', 'RN', 'TT', 
C     &           'ST', 'TF', 'FS', 'CA', 'GR', 
C     &           'CG', 'CV', 'CL', 'HD'/

      IRTRET  = 0
      IRTFLG  = 0

      MAXIM   = 0
      MAXIM2  = 0

      CALL SET_MPI(ICOMM,MYPID,MPIERR)  ! SET MYPID

C     NUMBER OF COMMAND LINE REGISTERS
      CALL REG_GET_USED(NSEL_USED)


      SELECT CASE(FCHAR(1:2))

      CASE ('20') ! ---------------------------------------------- 'HIS'
        CALL HISTOG()
        RETURN


      CASE ('DE') ! ----------------------------------------------- 'DE'
        CALL DELETF(FILNAM,LUN1)
        RETURN
C       DO NOT USE GOTO 9000 HERE AS IT CAUSES DOUBLE CLOSING ERROR

      CASE ('DU') ! ----------------------------------------------  'DU'

C       OPEN INPUT FILE, NO FOURIER INPUT ALLOWED 
 	CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &             NX,NY,NZ,
     &             MAXIM,'OVERWRITES INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

C       DETERMINE HSIG & HMODE
        CALL HIST(LUN1,0,0,NX,NY,NZ,HMIN,HMAX,HSIG,HMODE)

C       REMOVE OUT-LIERS
        CALL DUST(LUN1,NX,NY,NZ,HSIG,HMODE,IRTFLG)

	GOTO 9000

      CASE ('FI') ! ----------- FILE INFO ------------------------ 'FI'
                                       
        IF (FCHAR(4:4) == 'A') THEN

C	   FILE INFO FOR MULTIPLE FILES
           CALL FILERD(FILNAM,NLETI,NULL,'INPUT',IRTFLG)
           IF (IRTFLG == -1) GOTO 9000
           CALL FILGEN(FILNAM,NLETI,LUN1)

        ELSEIF (FCHAR(4:4) == 'C') THEN ! -------------------- 'FI CEN'

C	   FIND CENTER OF FILE 
C          OPEN INPUT FILE, NO FOURIER INPUT ALLOWED 
 	   CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &                  NX,NY,NZ,
     &                  MAXIM,'INPUT',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

           CALL FINDFILCEN(ITYPE,NX,NY,NZ, IXC,IYC,IZC,IRTFLG)

           IF (ITYPE == 1) THEN
              WRITE(NOUT,'(A,I0,A,I0,A)'),
     &                   '  IMAGE CENTER: (',IXC,',',IYC,')'
           ELSEIF (ITYPE == 3) THEN
              WRITE(NOUT,'(A,I0,A,I0,A,I0,A)'),
     &                   ' VOLUME CENTER: (',IXC,',',IYC,',',IZC,')'
           ENDIF

           IF (NSEL_USED > 0) THEN
C              OUTPUT TO SPIDER'S REGISTERS
               CALL REG_SET_NSEL(1,3, 
     &            FLOAT(IXC),FLOAT(IYC),FLOAT(IZC),0,0,IRTFLG)
           ENDIF

        ELSEIF (FCHAR(4:4) == 'H')  THEN ! --------------------- 'FI H'
C          RETRIEVE HEADER VARIABLE CONTENTS FROM SINGLE FILE 
           MAXIM = 2   ! ALLOWS QUERY ON STACK HEADER
C          USE ~7 TO ALLOW STACK HEADER ACCESS
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &         NX,NY,NZ,
     &         MAXIM,'INPUT~7',.TRUE.,IRTFLG)

          IF (IRTFLG .NE. 0) THEN
C             FILE NOT FOUND
              CALL ERRT(101,'OPENING FILE',IDUM)
              GOTO 9000 
          ENDIF
          CALL INQUIREHEAD(LUN1,NX,NY,NZ,IRTFLG)

        ELSEIF (FCHAR(4:4) == 'N') THEN !  --------------------- 'FI N'

C	   FILE INFO FOR SINGLE FILE, NO ERROR STOP IF NOT FOUND
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'Z',ITYPE,
     &             NX,NY,NZ,
     &             MAXIM,'INPUT',.TRUE.,IRTFLG)

           IF (IRTFLG .NE. 0) THEN
C             FILE NOT FOUND, ZERO REG. FOR NX, NY ...
              CALL REG_SET(1,0.0,NULL,IRTFLG)
              CALL REG_SET(2,0.0,NULL,IRTFLG)
              CALL REG_SET(3,0.0,NULL,IRTFLG)
              CALL REG_SET(4,0.0,NULL,IRTFLG)
              CALL REG_SET(7,0.0,NULL,IRTFLG)
           ENDIF

           CALL LUNGETIS_MRC(LUN1,IS_MRC,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000
           IF (IS_MRC) THEN
             CALL FILDAT_MRC(LUN1,NX)
           ELSE
             CALL FILDAT(LUN1,NX)
           ENDIF
 
        ELSEIF (FCHAR(4:4) == 'T') THEN !  --------------------- 'FI T'
C          TEST OF FILENAME SUBSTITUTION MECHANISM

           CALL FILERD(FILNAM,NLET1,NULL,'TEST FILE NAME',IRTFLG)
           IF (IRTFLG == -1) RETURN

           IRTFLG = -999
           CALL RDPRMC(CLINE,NLET2,.TRUE.,'CORRECT NAME',
     &                 NULL,IRTFLG)
           IF (IRTFLG == 0 .AND. 
     &        FILNAM(1:NLET1) .NE. CLINE(1:NLET2)) THEN
              WRITE(NOUT,*) '*** ERROR IN FILENAME FORMATION!!!'
              WRITE(NOUT,9065) FILNAM(1:NLET1),CLINE(1:NLET2)
9065          FORMAT( '*** GOT: ',A,' -- SHOULD BE: ',A/)
              CALL ERRT(100,'UTIL1',NE)
           ENDIF

        ELSE !  -------------------------------------------------- 'FI'

C	   FILE INFO FOR SINGLE FILE, ERROR IF NOT FOUND
           CALL FILERD(FILNAM,NLETI,NULL,'INPUT',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000
           IF (FILNAM(1:1) == '?') THEN
C             OLD FASHIONED "FR" SETTING
              WRITE(NOUT,*) 
     &          '*** OBSOLETE: PLEASE USE OPERATION <FR> NOW'
              IRTRET = 1
              BACKSPACE(NIN)
              IBCNT = IBCNT - 1
              RETURN
           ENDIF
           ILOCAT = INDEX(FILNAM,'@') 
           IF (FCHAR(4:4) == '[' .AND. ILOCAT == NLETI) THEN
C             NEED HEADER LOCATION
              MAXIM = 2
           ENDIF
           IF (FCHAR(4:4) == 'X' .AND. ILOCAT == NLETI) THEN
C             NEED HEADER LOCATION
              MAXIM = 2
           ENDIF
           CALL OPFILEC(0,.FALSE.,FILNAM,LUN1,'O',ITYPE,
     &                 NX,NY,NZ,MAXIM,' ',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

           CALL LUNGETIS_MRC(LUN1,IS_MRC,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000
           IF (IS_MRC) THEN
             CALL FILDAT_MRC(LUN1,NX)
           ELSE
             CALL FILDAT(LUN1,NX)
           ENDIF
        ENDIF
        GOTO 9000

      CASE ('HI') ! -------------- HISTOGRAM --------------------- 'HI'

        IF (FCHAR(4:6) == 'DOC')  THEN
          CALL ERRT(101,"USE OPERATION: 'HD'",IDUM)
          GOTO  9000
        ENDIF


C       OPEN INPUT FILE, FOURIER NOT ALLOWED 
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE1,
     &             NX1,NY1,NZ,
     &             MAXIM,'INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        IF (FCHAR(4:4) == 'E') THEN
C          ENTROPY HISTOGRAM
           CALL ENTROP(LUN1,NX1,NY1,NZ,ENTROPY,IRTFLG)

        ELSEIF (FCHAR(4:4) == 'J') THEN
C          JOINT HISTOGRAM FOR MUTUAL SHARED INFORMATION

C          MAKE SURE STATISTICS ARE CURRENT
           FMIN1 = FMIN
           FMAX1 = FMAX
           IF (IMAMI .NE. 1 .AND. ITYPE1 >= 0) 
     &        CALL NORM3(LUN1,NX1,NY1,NZ,FMAX1,FMIN1,AV)

C          OPEN SECOND INPUT FILE, FOURIER ALLOWED 
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'O',ITYPE2,
     &           NX2,NY2,NZ2, 
     &           MAXIM2,'SECOND INPUT',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 8999

           CALL SIZCHK(UNUSED,NX1,NY1,NZ, 0,
     &                        NX2,NY2,NZ2,0,IRTFLG) 
           IF (IRTFLG .NE. 0) GOTO 8999 

C          MAKE SURE STATISTICS ARE CURRENT
           FMIN2 = FMIN
           FMAX2 = FMAX
           IF (IMAMI .NE. 1.AND. ITYPE1 >= 0) 
     &        CALL NORM3(LUN2,NX2,NY2,NZ2,FMAX2,FMIN2,AV)

           IF (ITYPE1 >= 0) THEN
C             NOT FOURIER
              CALL RDPRI1S(NBINS,NOT_USED,
     &           'NUMBER OF BINS IN HISTOGRAM',IRTFLG)
              IF (NBINS < 1) THEN
                  CALL ERRT(1,'UTIL1',IDUM)
                  GOTO 8999 
              ENDIF

              WRITE(NOUT,*)' FIRST IMAGE RANGE:  ',FMIN1,'.....',FMAX1
              WRITE(NOUT,*)' SECOND IMAGE RANGE: ',FMIN2,'.....',FMAX2
              VMIN = MIN(FMIN1,FMIN2)
              VMAX = MAX(FMAX1,FMAX2)
              CALL RDPRM2S(VMIN,VMAX,NOT_USED,
     &                     'HISTOGRAM RANGE',IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

              CALL JOHIST(LUN1,LUN2,NX1,NY1,NZ,NBINS,
     &                   FMIN1,FMAX1,FMIN2,FMAX2,VMIN,VMAX,IRTFLG)
           ELSE
C             FOURIER
              NBINSA = 128
              NBINSP = 360
              CALL RDPRIS(NBINSA,NBINSP,NOT_USED,
     &           'NUMBER OF AMPLITUDE & PHASE BINS IN HISTOGRAM',
     &            IRTFLG)
              IF (NBINSA < 1 .OR. NBINSP < 1) THEN
                  CALL ERRT(1,'UTIL1',IDUM)
                  GOTO 8999 
              ENDIF

              CALL JOHISTF(LUN1,LUN2,NX1,NY1,NZ,
     &                     NBINSA,NBINSP,IRTFLG)
           ENDIF
           GOTO 8999

        ELSEIF (FCHAR(4:4) .NE. 'M') THEN
C          NORMAL HISTOGRAM
           CALL HIST(LUN1,0,LUN2,NX1,NY1,NZ,HMIN,HMAX,
     &               HSIG,HMODE)
        ELSE
C          HISTOGRAM UNDER MASK --------------------------------- 'HI M'
           IFORM1 = IFORM
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'O',ITYPE,
     &                  NX2,NY2,NZ2,MAXIM2,
     &                 'MASK',.FALSE.,IRTFLG)
           IFORM2 = IFORM
           IF (IRTFLG .NE. 0) GOTO 9000

C          IMAGES NUST HAVE SAME DIMENSIONS
           CALL SIZCHK(UNUSED,NX1,NY1,NZ, IFORM1,
     &                        NX2,NY2,NZ2,IFORM2,IRTFLG) 
           IF (IRTFLG .NE. 0) GOTO 9000 

           CALL HIST(LUN1,LUN2,LUN3,NX1,NY1,NZ,HMIN,HMAX,
     &               HSIG,HMODE)
           CLOSE(LUN2)
        ENDIF
	GOTO 9000
    
      CASE ('HD') ! --- HISTOGRAM OF A DOCUMENT FILE COLUMN ------ 'HD'

20      CALL HISD(LUN3)
        CLOSE(LUN3)
	GOTO 9000
    
      CASE ('LI') ! ---------------------------------------------  'LI'

C       CHECK FOR 'LI R' WHICH DOES NOT USE INPUT FILE JUST REGISTERS
5       IF (FCHAR(4:4) == 'R') THEN
C          FOR 'LI R', 'LI RT',  'LI RD' (listregs.f)
           CALL ERRT(101,'OBSOLETE OPERATION',NE)
           GOTO 9000

        ELSEIF (FCHAR (4:5) == '2D') THEN
C          UNDOCUMENTED COMMAND 'LI 2D' CALLED LISTIM
           CALL ERRT(101,'OBSOLETE OPERATION',NE)
           GOTO 9000
         ENDIF

C       NEED INPUT FILE, USE ~7 TO ALLOW STACK HEADER ACCESS
        DISP = 'O'
        IF (FCHAR(4:4) == 'T') DISP = 'Z'
 
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,DISP,ITYPE,
     &              NX,NY,NZ,
     &              MAXIM,'INPUT~7',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        DOCPRNT  = (FCHAR(4:4) == 'D') 
        TERMPRNT = (FCHAR(4:4) == 'T') 
 
        IF (NSEL_USED > 0 .AND. .NOT. TERMPRNT .AND. .NOT. DOCPRNT) THEN
C          SINGLE NUMBER (REGISTER) OPTION:
           CALL LISTITR(FILNAM,LUN1,NX,NY,NZ)
        ELSE
           CALL LISTIT(FILNAM,LUN1,NX,NY,NZ,DOCPRNT,TERMPRNT)
        ENDIF

        GOTO 9000

      CASE ('MO') ! ---------------- MODEL ----------------------- 'MO'
6	NX2  = 0
	NY2  = 0
	NZ = 1
	IFORM  = 1
	IF (FCHAR(4:4) == '3') THEN
	   NZ = 0
	   IFORM  = 3
	ENDIF
	CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'U',IFORM,
     &             NX2,NY2,NZ,
     &             MAXIM,'OUTPUT',.TRUE.,IRTFLG)

	IF (IRTFLG .NE. 0 .AND. FILNAM(1:1) .NE. '*') GOTO 9000

	IF (FCHAR(4:4) == '3') THEN
           IF (FCHAR(5:5) == 'H') THEN 
               CALL ERRT(101,'OBSOLETE OPERATION',IDUM)
               GOTO 9000 
           ELSE
C             FOR 'MO 3'
              CALL MODEL3(LUN1,LUN2,FILNAM,NX2,NY2,NZ)
              CLOSE(LUN2)
           ENDIF
        ELSE
C          FOR 'MO'
           !!write(3,*) ' in util1,   irtflg!!!!!!:',irtflg
           CALL MODEL(LUN1,NX2,NY2)
	ENDIF
        GOTO 9000
           
      CASE ('PK') ! ---------------- PEAK SEARCH ----------------- 'PK'

 	CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NX1,NY1,NZ1,
     &              MAXIM,'INPUT',.TRUE.,IRTFLG)
	IF (IRTFLG .NE. 0) GOTO 9000

        FMAX1 = FMAX
	IF (FCHAR(4:4) == 'M') THEN
           IF (IMAMI.NE.1)
     &     CALL NORM3(LUN1,NX1,NY1,NZ1,FMAX1,FMIN1,AVR1)

           CALL SPEAKM(LUN1,NX1,NY1,NZ1,FMAX1)

        ELSEIF (FCHAR(4:4) == '3') THEN
           CALL SPEAK3(LUN1,NX1,NY1,NZ1,FCHAR(5:5),LUNDOC)
	
        ELSE
            CALL RDPRMI(ML,NOR,NOT_USED,
     &        'NUMBER OF PEAKS, CENTER ORIGIN OVERRIDE (0/1)')
           IF (ML < 1)   ML = 1

           IF (FCHAR(4:4) == 'C' .OR. FCHAR(5:5) == 'C') THEN
             CALL SPEAKC(FILNAM,LUN1,NX1,NY1,MAXDIM,FCHAR(4:4),
     &                   LUNDOC,ML,NOR)
     
           ELSEIF ( FCHAR(5:5) == 'R') THEN
             CALL SPEAKR(FILNAM,LUN1,NX1,NY1,MAXDIM,FCHAR(4:4),
     &                   LUNDOC,ML,NOR)  
             
           ELSE
             CALL SPEAK(FILNAM,LUN1,NX1,NY1,MAXDIM,FCHAR(4:4),
     &                   LUNDOC,ML,NOR)
           ENDIF
        ENDIF
        GOTO 9000
     

      CASE ('RA') ! ------------- RAMP --------------------------- 'RA'

        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &              NX,NY,NZ,
     &              MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        IF (ITYPE  .NE. 1) THEN
   	   CALL ERRT(101,'OPERATION ONLY WORKS ON IMAGES',NE)
           GOTO 9001
        ENDIF 

        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &              NX,NY,NZ,
     &              MAXIM2,'OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        CALL  RAMP_P(LUN1,LUN2,NX,NY,NOUT)
	GOTO 8999
            
      CASE ('RN') ! ---------------- RENAME ---------------------- 'RN'

        CALL ERRT(101,'OBSOLETE OPERATION',NE)
	GOTO 9001


      CASE ('TT') !   -----------CHANGE TITLE--------------------- 'TT'


        IF (FCHAR(4:4) == 'C') THEN
C          OPEN THE FILE THAT CONTAINS DESIRED TITLE
	   CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &                  NX1,NY1,NZ1,
     &                  MAXIM,'TITLE SOURCE',.TRUE.,IRTFLG)
	   IF (IRTFLG .NE. 0) GOTO 9001
        ENDIF

C       OPEN THE FILE THAT RECEIVES TITLE 
	CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'O',ITYPE,
     &               NX1,NY1,NZ1,
     &               MAXIM,'OUTPUT',.TRUE.,IRTFLG)
	IF (IRTFLG .NE. 0) GOTO 9000

        CALL LUNGETIS_MRC(LUN2,IS_MRC,IRTFLG)
        IF (IS_MRC) THEN
          CALL ERRT(101,'NO TITLE IN MRC FILES',IDUM)
          IRTRET = 1
          RETURN
        ENDIF

C       GET NEW TITLE HERE 
        IF (FCHAR(4:4) == 'C') THEN
           CALL LUNGETTITLE(LUN1,CTIT,LENTIT,IRTFLG)
        ELSE
           CALL RDPRMC(CTIT,LENTIT,.FALSE.,'NEW TITLE',NULL,IRTFLG)
        ENDIF
	IF (IRTFLG .NE. 0) GOTO 9000

C	TITLE ALTERATION CAN PROCEED NOW
        CALL TITLE(LUN2,CTIT,LENTIT,.TRUE.,IRTFLG)
        CLOSE(LUN2)
	GOTO 9000

      CASE ('ST') !  ------ SET BUFFER LOCATIONS ----------------- 'ST'

C       DISP OF "Z" ALLOWS  CORRECTING STACK ERROR

        IF (FCHAR(4:4) == 'H')  THEN  
C          SET HEADER VARIABLES IN FILE, SELECT BY NAME

           MAXIM = 2   ! ALLOWS OPENING OVERALL STACK HEADER
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'Z',ITYPE,
     &             NX,NY,NZ,
     &             MAXIM,'SET HEADER VALUES IN',.TRUE.,IRTFLG)


C          SET HEADER VARIABLES IN FILE, IRTFLG==5 IS TO ALLOW IMGNUM CHANGE 

           IF (IRTFLG == 0 .OR. IRTFLG == 5) THEN
              CALL SETHEAD(LUN1,NX,NY,NZ,IRTFLG)    ! MRC OK
           ENDIF

        ELSEIF (FCHAR(4:4) == 'E')  THEN  
C          CORRECT BAD EM2EM HEADER VARIABLES IN FILE
           MAXIM = 2   ! ALLOWS OPENING OVERALL STACK HEADER
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'Z',ITYPE,
     &             NX,NY,NZ,
     &             MAXIM,'FIX STACKED IMAGE HEADER VALUES IN',
     &            .TRUE.,IRTFLG)

           CALL LUNGETIS_MRC(LUN1,IS_MRC,IRTFLG)
           IF (IS_MRC) THEN
             CALL ERRT(101,'NOT RELEVENT FOR MRC FILES',IDUM)
             IRTRET = 1
             RETURN
           ENDIF

C          CORRECT BAD EM2EM HEADER VARIABLES IN FILE
           CALL SETHEADEM2(LUN1,NY,IRTFLG)    ! MRC WAS ERRT

        ELSE
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'Z',ITYPE,
     &               NX,NY,NZ, 
     &               MAXIM,'INPUT~7',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C  	   SET LABEL VALUES TO SOLICITED INPUT
           CALL SETVAL(LUN1,NX,NY,NZ)    ! MRC OK
        ENDIF
        GOTO 9000

      CASE ('18') ! ------------DETERMINE DEFOCUS --------------- 'CTF'

        LENF     = LNBLNK(FCHAR(4:))
        LENCOMMA = INDEX(FCHAR(4:4+LENF-1),',')
        IF (LENCOMMA > 1) LENF = LENCOMMA-1
        LENBLANK = INDEX(FCHAR(4:4+LENF-1),' ')
        IF (LENBLANK > 1) LENF = LENBLANK-1

        SELECT CASE(FCHAR(4:4+LENF-1))

          CASE ('F','FI','FIN','FIND')
             CALL CTFFIND3()
	     	     
          CASE ('E','ED')
             CALL CTFED()

          CASE DEFAULT
             CALL ERRT(101,'UNIDENTIFIED OPERATION',IDUM)
        END SELECT
        GOTO 9000


      CASE ('TF') ! ------------TRANSFER FUNCTION ---------------- 'TF'

        LENF     = LNBLNK(FCHAR(4:))
        LENCOMMA = INDEX(FCHAR(4:4+LENF-1),',')
        IF (LENCOMMA > 1) LENF = LENCOMMA-1
        LENBLANK = INDEX(FCHAR(4:4+LENF-1),' ')
        IF (LENBLANK > 1) LENF = LENBLANK-1

        SELECT CASE(FCHAR(4:4+LENF-1))

          CASE ('D')
             CALL TRAFD(LUN1)

          CASE ('C3')
             !write(6,*) ' calling: TRAFC3  C  '
             CALL TRAFC3(LUN1,.FALSE.)

          CASE ('COR','CO')
             CALL RCTFONE(LUN1) 

          CASE ('CRF','CR')
             CALL TFCRF

          CASE ('CT3')
             !write(6,*) ' Calling: TRAFC3  CT '
             CALL TRAFC3(LUN1,.TRUE.)

          CASE ('CTS')
             CALL RCTFSS(LUN1,LUN2)  
	  	   
          CASE ('CT')
             !write(6,*) ' Calling: TRAFC  CT'
             CALL TRAFC(LUN1,.TRUE.)

          CASE ('C')
             !write(6,*) ' Calling: TRAFC  C'
             CALL TRAFC(LUN1,.FALSE.)

          CASE ('DDF','DD')
             CALL DEFOCUS(IRTFLG)

          CASE ('DEV','DE')
             CALL ENVELOPE(IRTFLG)

          CASE ('DNS','DN')
             CALL NOISE(LUN1,LUN2)

          CASE ('ED')
             CALL TFED

          CASE ('F','FI','FIN','FIND')
             CALL CTFFIND3()
	     	     
          CASE ('L','LIS')
             CALL TRAFL
          
         CASE ('LM4','LM')
             CALL TFLM4
	        
          CASE ('RCTF','R','RC','RCT')
             CALL ERRT(101,'UNDOCUMENTED OPERATION REMOVED 2014',NE)

          CASE ('SIM','SI')
	     CALL TRAFSIM(LUN1)   
   
          CASE ('SNR','SN')
            CALL TFSNR

          CASE DEFAULT
             CALL TRAF(LUN1)
        END SELECT
        GOTO 9000

      CASE ('FS') ! ------------- FILE STATISTICS ---------------- 'FS'


        SELECT CASE(FCHAR(4:4))

          CASE ('V')
            CALL FINDMINSORMAXS(LUN1,LUNDOC)

          CASE ('L')
            CALL QSTATLOC(LUN1,LUNDOC,LUNXM)

          CASE DEFAULT
            !write(6,*) ' In util1 - lun1,lun2:',lun1,lun2
            CALL QSTAT(LUN1,LUN2,LUNDOC,LUNXM)
            CLOSE(LUN2)

        END SELECT
	GOTO 9000

      CASE ('CA') !  ------------ CLUSTER ANALYSIS   ------------- 'CA'

        IENDOP2 = INDEX(FCHAR(4:),' ') - 1 + 4

        SELECT CASE(FCHAR(4:IENDOP2))

          CASE ('S')
C            FACTOR MAP CALCULATION
 	     CALL JPMSK1()

          CASE ('SM','SME')
C            FACTOR MAP PLOT
             CALL SGRAF(LUN1,LUN2,LUN3,LUN4)
             CLOSE(LUN2)

          CASE ('SMI')
C            INACTIVE FACTOR MAP PLOT
             CALL JPMSK3(LUN1,LUN2,LUN3,LUN4,LUN5)

          CASE ('SR', 'SRD', 'SRI', 'SRA', 'SRE')
C            IMAGE RECONSTITUTION
 	     CALL JPMSK2(LUN1,LUN2,LUN3,LUN4,LUN5)
             CLOSE(LUN3)
             CLOSE(LUN4)
             CLOSE(LUN5)

          CASE ('EIGPCT','NOISE','NOIS','EIG','NOI','NO','N')
C            COPY EIG % TO DOC FILE
 	     CALL EIGPERCENT()

          CASE ('VIS')
C            VISUAL MAP CREATION
             CALL VISMAP(LUN1,LUN2,LUN3,LUN4)

          CASE DEFAULT
C            'CA E', 'CA ES', removed PAP 10/05/99
             CALL ERRT(101,'UNIDENTIFIED OPERATION',IDUM)

        END SELECT
        GOTO 8999

      CASE ('GR') ! --------- GRAPH A ROW IN RESULTS FILE  -------- 'GR'

16      CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &         NX1,NY1,NZ, MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000 

	CALL GRAPHP(LUN1,NX1,NY1)
	GOTO 9000

      CASE ('CG') ! ------------ CENTER OF GRAVITY ---------------- 'CG'

C       3-D CENTER OF GRAVITY AND RADIUS OF GYRATION

17	IF (FCHAR(4:4) == 'S' .OR. FCHAR(4:5) == 'PH') THEN
          CALL FINDCENT()
	  GOTO 9001
        ENDIF

	CALL CENGR3(LUN1)
	GOTO 9000      

      CASE ('CV') ! ---------------------------------------------- 'CV'

C       POCS PROGRAMS (06/05/90), MODULAR POCS PROGRAMS: (12/5/91) M.R.

        IF (FCHAR(4:8) == 'REPL2') THEN 
C          tdfrepl uses old Fourier format and was disabled
           CALL ERRT(101,'OBSOLETE OPERATION',NE)
C	   CALL TDFREPL(LUN1,LUN2,LUN3)
        ELSE
           CALL MRREPLACE(LUN1,LUN2)
        ENDIF
        GOTO 9000


      CASE ('CL') ! ---------------------------------------------- 'CL'

        SELECT CASE (FCHAR(4:5))

           CASE('KM')      ! KMEANS CLUSTERING 
           CALL SUBKMNS(LUN1,LUN2)

           CASE('HC')      ! HIERARCHICAL CLUSTERING 
           CALL HCLS(LUN1,LUN2,LUN3)

           CASE('HD')      ! HIERARCHICAL CLUSTERING, CALCULATE CLASSES 
   	   CALL HDLS(LUN1,LUN2)

           CASE('HE')      ! HIERARCHICAL CLUSTERING, CREATE DOC FILE
           CALL HELS(LUN1,LUN2)

           CASE('CL')      ! HIERARCHICAL CLUSTERING 
           CALL SCLASSI(LUN1,LUN2,LUN3)

           CLOSE(LUN3)

        END SELECT
        RETURN

      END SELECT

C       -------------------------------------------------------- END
8999  CLOSE(LUN2)
9000  CLOSE(LUN1)
9001  CONTINUE

      END
