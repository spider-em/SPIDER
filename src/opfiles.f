 
C++*********************************************************************
C                                                                      *
C  OPFILES.F   NEW                            12/15/06  ARDEAN LEITH    
C              BAD NUMBRT() TRAP              05/21/09  ARDEAN LEITH  
C              ASKNAM, PROMPTEX               12/06/10  ARDEAN LEITH 
C              NX...                          03/26/12  ARDEAN LEITH 
C              COPY NON SPIDER INPUT          05/26/14  ARDEAN LEITH 
C              COPY NON SPIDER INPUT          05/26/14  ARDEAN LEITH 
C              LOCAST, ASKLIST FOR ILIST      10/02/14  ARDEAN LEITH
C
C ********************************************************************** 
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C  CONTAINS:  OPFILES, GETOLDIMG, GETNEWIMG, NEXTFILE, NEXTFILES 
C
C  OPFILES(LUNCP,LUNIMG,LUNDOC,LUNXM,  ASKNAM,FILPAT,NLET, DISP,
C          ITYPE,NX,NY,NZ,MAXIM, PROMPT,FOUROK,
C          ILIST,NIMAXT, NTOT,IMGNUM, IRTFLG)
C 
C  PURPOSE: SOLICITS FILE NAME(S) AND OPENS FILE(S)
C           SUPPORT ROUTINE FOR CONVERTING OPERATIONS TO 
C           WORK ON WHOLE STACK OR WITH SELECTION DOC FILE.
C
C  PARAMETERS:
C        LUNCP      UNIT TO COPY HEADER VALUES FROM               (SENT)
C        LUNIMG     UNIT TO OPEN FILE ON                          (SENT)
C        LUNDOC     UNIT TO OPEN LIST DOC FILES ON                (SENT)
C        LUNXM      UNIT TO OPEN XMIPP SELFILE ON                 (SENT)
C        ASKNAM     FLAG TO ASK FOR FILE NAME                     (SENT)
C        FILPAT     FILE NAME PATTERN                             (RET)
C        NLET       CHARS IN FILE NAME PATTERN                    (RET)
C        DISP       CHARACTER CONTAINING ONE OF THE               (SENT) 
C                   FOLLOWING DISPOSITION SPECIFICATIONS:
C                   'O'   -  FILE IS ASSUMED TO EXIST.  DIMENSIONS,
C                            ITYPE AND HEADER INFO (IN COMMON) ARE 
C                            RETURNED TO THE CALLING PROGRAM. 
C                   'B'   -  SAME AS OLD BUT NO LIMIT ON BUFFER LENGTH
C                            FOR OPENCHK. 
C                   'Z/E' -  THE FILE IS ASSUMED TO EXIST.
C                            IF FILE DOES NOT EXIST, THEN BATCH DOES
C                            NOT STOP. (ONLY DIFFERENCE FROM 'O'). 
C                   'N'  -   WANT NEW FILE. SEND NX, NY, NZ & ITYPE.
C                   'U'  -   IT IS NOT KNOWN IF THE FILE EXISTS.  
C                            SEND NX, NY, NZ & ITYPE. IF FILE 
C                            ALREADY EXISTS, IT WILL BE REPLACED.
C        ITYPE      IFORM FOR FILE                         (SENT OR RET) 
C        NX,NY,NZ   IMAGE SIZE                             (SENT OR RET)
C        MAXIM      STACK INDICATOR IF DISP == 'I'                (SENT)
C                   STACK INDICATOR                               (RET)
C        PROMPT     PROMPT FOR FILNAME                            (SENT)
C                     IF NOT (ASKNAM) THIS IS FILE NAME           (SENT)
C                     ~ (TILDE) IN LAST CHAR. SAYS SKIP
C                       "FILE" AT END OF PROMPT
C                     ~9 IN NEXT TO LAST OR
C                        NEXT-TO-NEXT-TO LAST
C                        ACCEPTS AN EXTENSION
C                        (OTHERWISE DISCARDED!)
C                     ~6 KEEPS OLD DATE/TIME
C        FOUROK     CAN USE EXISTING FOURIER FILES?               (SENT)
C        ILIST      IMAGE NUMBER LIST                             (RET)
C                     IF NIMAXT < 0 MUST BE SENT
C                     NOT USED IF SINGLE IMAGE/SELFILE 
C        NIMAXT     MAX LENGTH OF IMAGE NUMBER LIST               (SENT)
C                     <0 MEANS DO NOT ASK FOR LIST
C        UNUSED     UNUSED                                          (?)
C        NTOT       # OF IMAGES IN IMAGE NUMBER LIST              (RET)
C                     ZERO FOR SINGLE IMAGE AND NO * 
C        IMGNUM     IMAGE NUMBER THAT IS CURRENTLY OPEN       (SENT/RET)
C                   ON INPUT: IF (BARESTACK) IS # WANTED
C                   ON OUTPUT:   <0 IS SELFILE IN USE 
C        IRTFLG     ERROR FLAG (0 IS NORMAL)                      (RET)
C                      -1 GOTO PREVIOUS QUESTION
C 
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12 
C--********************************************************************* 

       SUBROUTINE OPFILES(LUNCP,LUNIMG,LUNDOC,LUNXM,
     &                    ASKNAM,FILPAT,NLET, DISP,
     &                    ITYPE,NX,NY,NZ,MAXIM,
     &                    PROMPT,
     &                    FOUROK,ILIST,NIMAXT, 
     &                    UNUSED,NTOT,IMGNUM, IRTFLG) 
 
       INCLUDE 'CMBLOCK.INC' 
       INCLUDE 'CMLIMIT.INC' 
 
       INTEGER                   :: LUNCP,LUNIMG,LUNDOC
       LOGICAL                   :: ASKNAM
       CHARACTER(LEN=*)          :: FILPAT 
       INTEGER                   :: NLET
       CHARACTER(LEN=1)          :: DISP 
       INTEGER                   :: ITYPE,NX,NY,NZ,MAXIM
       CHARACTER(LEN=*)          :: PROMPT 
       LOGICAL                   :: FOUROK
       INTEGER                   :: ILIST(*)
       INTEGER                   :: NIMAXT,NDUM,NTOT,IMGNUM, IRTFLG

       INTEGER                   :: NIMAXP
       LOGICAL                   :: SAYIT,ASKLIST
       CHARACTER (LEN=1)         :: CDUM 
       CHARACTER (LEN=MAXNAM)    :: FILNAM 
       CHARACTER (LEN=2*MAXNAM)  :: MESG 
       CHARACTER (LEN=100)       :: PROMPTEX
       CHARACTER (LEN=1)         :: DISPT 
       CHARACTER (LEN=1)         :: NULL = CHAR(0)

       INTEGER                   :: LUNOP,NLETT

       LOGICAL                   :: ISOPEN,GOTFILE
       CHARACTER (LEN=MAXNAM)    :: FILNAMT 

       LOCTILDE = INDEX(PROMPT,'~') ! unfinished

       IF (ASKNAM .AND. PROMPT == NULL) THEN
C         ASK FOR FILE NAME, CAN ACCEPT EXTENSION

          IF (DISP == 'N' .OR. 
     &        DISP == 'I' .OR.
     &        DISP == 'U') THEN
C            NEW FILE, USE DEFAULT PROMPT
             PROMPTEX = 
     &          'OUTPUT FILE NAME OR TEMPLATE (E.G. STK@****)~~9'
          ELSE
C            OLD FILE, USE DEFAULT PROMPT
             PROMPTEX = 
     &          'INPUT FILE NAME OR TEMPLATE (E.G. STK@****)~~9'
          ENDIF
       
          CALL FILERD(FILPAT,NLET,NULL,PROMPTEX,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN 

       ELSEIF (ASKNAM) THEN
C         ASK FOR FILE NAME USING PROMPT
      
          CALL FILERD(FILPAT,NLET,NULL,PROMPT,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN 

       ELSE
C         USE FILENAME SENT IN: PROMPT
          FILPAT = PROMPT
          NLET   = LNBLNKN(FILPAT)
       ENDIF
 
       IMGNUMIN = IMGNUM
       IMGNUM   = 0 
     
       LOCAT    = INDEX(FILPAT(1:NLET),'@')   
       LOCAST   = INDEX(FILPAT(1:NLET),'*')
       
       !write(6,*)' locast,locat,nlet:',locast,locat,nlet,filpat(1:nlet)
       !write(6,*)' ILIST:',ILIST(1)
       !write(6,*)' locast,locat,nstack:',locast,locat,maxim

       ASKLIST = (NIMAXT >= 0)
       NIMAXP  = ABS(NIMAXT)
 
       !!IF (ASKNAM .AND. LOCAST > 0 .AND. ASKLIST) THEN
       IF (LOCAST > 0 .AND. ASKLIST) THEN
C         GET LIST OF IMAGES FROM DOC. FILE OR INPUT LINE
          NTOT   = 0 
          CALL FILELIST(.FALSE.,LUNDOC,CDUM,NDUM,ILIST,NIMAXP,
     &                  NTOT,' ',IRTFLG)
          IF (IRTFLG .NE. 0) RETURN
       ELSEIF (.NOT. ASKLIST) THEN
          NTOT = NIMAXP
       ENDIF

       IF (LOCAT > 0 .AND. LOCAST > LOCAT) THEN
C         TEMPLATED STACKED FILE: STK@**** -------------- _9@* or STK@**

C         OPEN STACK HEADER
          FILNAM = FILPAT(1:LOCAT)
          MAXIM  = 1  
             
          IF (FILNAM(1:1) .NE. '_') THEN
             CALL FILNAMANDEXT(FILNAM(1:LOCAT-1),DATEXC,
     &                         FILNAMT,NLET,.TRUE.,IRTFLG)
             INQUIRE(FILE=FILNAMT,OPENED=ISOPEN,NUMBER=LUNOP)
             MESG = '  FILE: ' // FILNAMT(1:NLET) //
     &              '  ALREADY OPENED ON' // NULL
             LENE = LNBLNKN(MESG)
             IF (ISOPEN .AND. LUNOP .NE. LUNIMG) THEN
                WRITE(NOUT,'(A,I3)') MESG(1:LENE),LUNOP 
                !CALL ERRT(102,MESG(1:LENE),LUNOP) 
                IRTFLG = -2
                GOTO 9000
             ENDIF
          ENDIF

	  CALL OPFILEC(LUNCP,.FALSE.,FILNAM,LUNIMG,DISP,
     &                 ITYPE,NX,NY,NZ, 
     & 		       MAXIM,' ',FOUROK,IRTFLG) 
	  IF (IRTFLG .NE. 0) GOTO 9000 
          CALL LUNSETISBARE(LUNIMG,.FALSE.,IRTFLG)

C         OPEN FIRST FILE IN STACK SERIES
          IMGWANT  = ILIST(1)
          IF (IMGWANT < 0 .OR. IMGWANT > 10000000) THEN
              CALL ERRT(102,'INVALID IMAGE NUMBER',IMGWANT)
              IRTFLG = 1
              GOTO 9000
          ENDIF

          SAYIT    = .TRUE.
          IF (DISP == 'U' .OR. DISP == 'N') THEN
             CALL GETNEWIMG(LUNCP,LUNIMG,LUNDOC,FILPAT,IMGWANT,
     &                      SAYIT,IMGNUM,IRTFLG)
          ELSE
	     CALL GETOLDIMG(LUNIMG,LUNDOC,FILPAT, IMGWANT,SAYIT,
     &                      FOUROK,IMGNUM,IRTFLG)
          ENDIF
          IF (IRTFLG .NE. 0) GOTO 9000

C         RETRIEVE CURRENT MAXIMUM IMAGE NUMBER FROM OVERALL HEADER
          CALL LUNGETMAXIM(LUNIMG,MAXIM,IRTFLG)

          !write(6,'(A,A,A,i6,a,i6,a,i6)')' Opened templated stack: ',
!     &          FILPAT(1:nlet),' At image: ',IMGNUM, 
!     &          '  Maxim: ',maxim

       ELSEIF (LOCAT == NLET) THEN
C         WHOLE BARESTACK:  STK@  --------------------------_9@ or STK@
          DISPT = DISP
          IF (DISP == 'I') THEN
C             OPEN NEW BARE INDEXED STACK 
              DISPT = 'N'
              MAXIM = -MAXIM  ! FLAG FOR INDEXED STACK
          ELSE
              MAXIM = 1
          ENDIF

          IF (FILNAM(1:1) .NE. '_') THEN
             CALL FILNAMANDEXT(FILPAT(1:LOCAT-1),DATEXC,
     &                         FILNAMT,NLET,.TRUE.,IRTFLG)
             INQUIRE(FILE=FILNAMT,OPENED=ISOPEN,NUMBER=LUNOP)
             MESG = '  FILE: ' // FILNAMT(1:NLET) //
     &              '  ALREADY OPENED ON' // NULL
             LENE = LNBLNKN(MESG)
             IF (ISOPEN .AND. LUNOP .NE. LUNIMG) THEN
                WRITE(NOUT,'(A,I3)') MESG(1:LENE),LUNOP 
                !CALL ERRT(102,MESG(1:LENE),LUNOP) 
                IRTFLG = -2
                GOTO 9000
             ENDIF
          ENDIF

          CALL OPFILEC(LUNCP,.FALSE.,FILPAT,LUNIMG,DISPT,
     &                 ITYPE,NX,NY,NZ, 
     & 		       MAXIM,'INPUT',FOUROK,IRTFLG) 
	  IF (IRTFLG .NE. 0) GOTO 9000

C         OPEN FIRST FILE IN STACK, UNLESS SPECIFIED FOR NEW BARE STACK
          IMGWANT = 1
          SAYIT   = .TRUE.

          IF (DISP == 'U' .OR. 
     &        DISP == 'I' .OR.
     &        DISP == 'N') THEN
C            NEW BARE STACK, OPEN REQUESTED FILE IN STACK

             IF (IMGNUMIN > 0) IMGWANT = IMGNUMIN
             CALL GETNEWIMG(LUNCP,LUNIMG,LUNDOC,FILPAT,IMGWANT,
     &                      SAYIT,IMGNUM,IRTFLG)
             IF (IRTFLG .NE. 0) GOTO 9000

          ELSE
C            EXISTING BARE STACK, OPEN FIRST FILE IN STACK
	     CALL GETOLDIMG(LUNIMG,LUNDOC,FILPAT, IMGWANT,
     &                      SAYIT,FOUROK,IMGNUM,IRTFLG)
             IF (IRTFLG .NE. 0) GOTO 9000

C            CREATE IMAGE NUMBER LIST IN: ILIST
             NTOT = 0
             DO I= 1,MAXIM
                NTOT = NTOT + 1
                IF (NTOT > NIMAXP) THEN
                   CALL ERRT(102,'IMAGE # LIST OVERFLOW AT IMAGE',NTOT)
                   GOTO 9000
                ENDIF
                ILIST(NTOT) = I
             ENDDO
          ENDIF

         !write(6,*)' Opened bare stack:',FILPAT(1:NLET),' Img:',IMGNUM
         
       ELSEIF (LOCAST > 0) THEN
C         A SIMPLE FILE TEMPLATE: IMG*** ----------------------- IMG***

C         OPEN FIRST FILE IN THE SERIES
          IMGNUM = ILIST(1)
          IF (IMGNUM < 0 .OR. IMGNUM > 10000000) THEN
              CALL ERRT(102,'INVALID IMAGE NUMBER',IMGNUM)
              IRTFLG = 1
              GOTO 9000
          ENDIF
          CALL  FILGET(FILPAT,FILNAM,NLET,IMGNUM,IRTFLG)
	  IF (IRTFLG .NE. 0) GOTO 9000 

          MAXIM = 0 
	  CALL OPFILEC(LUNCP,.FALSE.,FILNAM,LUNIMG,DISP,
     &                 ITYPE,NX,NY,NZ, 
     & 		       MAXIM,' ',FOUROK,IRTFLG) 
	  IF (IRTFLG .NE. 0) GOTO 9000 

          !write(6,*)' Opened templated file: ',FILPAT(1:NLET),
      !&              '  for:',NTOT,' images.'

       ELSE
C         SINGLE SIMPLE INPUT FILE: IMG001 --------------------- IMG001
C         OR XMIPP SELFILE LISTING FILE: SELX ----------------- SELFILE
C         OR TRYING TO COPY NON-SPIDER FILE   ----------------- NONSPIFILE
               
C         CHECK FOR XMIPP SELFILE LIST
          IF (LUNXM > 0) THEN
             !write(6,*)' Filpat for openxmsel: ',filpat(:nlet) 
             CALL OPENXMSELFILE(FILPAT,LUNXM,FILNAM,NLET,NTOT,IRTFLG)
             !write(6,*)' Filnam from openxmsel: ',ntot,':',filnam(:nlet) 

             INQUIRE(FILE=FILNAM(:NLET),EXIST=GOTFILE)
           
             IF (NTOT > 0 .AND. GOTFILE) THEN
C               OPEN FIRST FILE IN XMIPP SELFILE LIST
                MAXIM = 0  
	        CALL OPFILEC(LUNCP,.FALSE.,FILNAM(:NLET),LUNIMG,DISP,
     &                       ITYPE,NX,NY,NZ, 
     & 		             MAXIM,'dum~9',FOUROK,IRTFLG) 
	        IF (IRTFLG .NE. 0) GOTO 9000 

                IMGNUM     = -1
                !write(6,*)' Opened selfile image: ',filnam(1:nlet) 
                RETURN

             ENDIF
          ENDIF

C         SINGLE SIMPLE INPUT FILE: IMG001 ------------------ IMG001
          MAXIM = 0  
	  CALL OPFILEC(LUNCP,.FALSE.,FILPAT,LUNIMG,DISP,
     &                 ITYPE,NX,NY,NZ, 
     & 		       MAXIM,PROMPT,FOUROK,IRTFLG) 

C         RETURN FILENAME WITH ANY EXTENSION IF NOT SPIDER IMAGE
          IF (IRTFLG == 5) NLET = lnblnkn(FILPAT)
	  IF (IRTFLG .NE. 0) GOTO 9000 

          NTOT       = 0
          IMGNUM     = 1

          !write(6,*)' Opened simple file: ',FILPAT(1:nlet) 
       ENDIF

9000   RETURN

       END 
 

C++*********************************************************************
C
C GETOLDIMG.F   FROM GETNXTSTK                     JAN 02 ARDEAN LEITH
C
C **********************************************************************
C
C    GETOLDIMG(LUN,LUNXM,FILPAT,NWANT, SAYIT,FOUROK,NGOT,IRTFLG)
C
C    PURPOSE:       TO OPEN A SPECIFIED IMAGE WITHIN STACK FOR RANDOM 
C                   ACCESS READING/WRITING.
C
C    PARAMETERS:
C        LUN        LUN NUMBER FOR FILNAM                         (SENT)
C        LUNXM      LUN FOR XM SELFILE                            (SENT) 
C        FILPAT     FILENAME PATTERN                              (SENT)
C        NWANT      IMAGE NUMBER WANTED (<0 IS SELFILE)           (SENT) 
C        SAYIT      SAY FILE OPENING INFO.                        (SENT)
C        FOUROK     FOURIER INPUT OK FLAG                         (SENT) 
C        NGOT       IMAGE NUMBER FOUND                            (RET.) 
C        IRTFLG     ERROR RETURN FLAG.                            (RET.)
C                   IRTFLG = -1    END OF FILE BEFORE NWANT
C                   IRTFLG =  0    NORMAL RETURN, IMAGE IS STACK
C                   IRTFLG =  2    IMAGE NOT IN USE
C
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************
 
	SUBROUTINE GETOLDIMG(LUN,LUNXM,FILPAT,NWANT, SAYIT,
     &                       FOUROK,NGOT,IRTFLG)

        INCLUDE 'CMLIMIT.INC'

        INTEGER                :: LUN     
        INTEGER                :: LUNXM     
        CHARACTER(LEN=*)       :: FILPAT
        INTEGER                :: NWANT     
        LOGICAL                :: SAYIT
        LOGICAL                :: FOUROK
        INTEGER                :: NGOT     
        INTEGER                :: IRTFLG     

        CHARACTER(LEN=MAXNAM)  :: FILNAM
        CHARACTER(LEN=1)       :: NULL = CHAR(0)

        INTEGER                :: NLET,LOCAST,LOCAT

        NLET   = lnblnkn(FILPAT)

        LOCAST = INDEX(FILPAT(1:NLET),'*')
        LOCAT  = INDEX(FILPAT(1:NLET),'@')

        !write(6,*)' locast,locat:',locast,locat,nlet,filpat(1:nlet)
        !write(6,*)' getoldimg, nwant,: ',nwant,':',filpat(1:nlet)
        
        IF (NWANT < 0) THEN
C          XMIPP SELFILE SIMPLE IMAGE ------------------------- SELAAA

C          RECOVER EXISTING IMAGE SIZE & TYPE
           CALL LUNGETSIZE(LUN,NX1,NY1,NZ1,IRTFLG)
           CALL LUNGETTYPE(LUN,ITYPE1,IRTFLG)
           CLOSE(LUN)    ! USUALLY STILL OPEN
 
C          LOAD FILNAM FROM SELFILE
           CALL GETNEXT_XMSEL(LUNXM,.TRUE.,FILNAM,NLET,IRTFLG)
	   IF (IRTFLG .NE. 0) RETURN 

C          OPEN EXISTING FILE: FILNAM  (HAS EXTENSION)
           MAXIM = 0  
	   CALL OPFILEC(0,.FALSE.,FILNAM(:NLET),LUN,'O',ITYPE,
     &                 NX,NY,NZ, 
     & 		       MAXIM,'~9',FOUROK,IRTFLG) 
	   IF (IRTFLG .NE. 0) RETURN 

C          NEW IMAGE SIZE SHOULD BE SAME AS PREVIOUS FILE
           CALL SIZCHK(NULL,NX1,NY1,NZ1,ITYPE1,
     &                      NX ,NY, NZ, ITYPE, IRTFLG)
	   IF (IRTFLG .NE. 0) RETURN 
 
           !write(6,*)' Opened old xmipp selfile file: ',filnam(1:nlet)

           NGOT   = NWANT   
           IRTFLG = 0
           RETURN

        ELSEIF (LOCAST > 0 .AND. LOCAT <= 0) THEN 
C          TEMPLATED SIMPLE IMAGE --------------------------- IMG***

           CALL LUNGETSIZE(LUN,NX1,NY1,NZ1,IRTFLG)
           CALL LUNGETTYPE(LUN,ITYPE1,IRTFLG)
           CLOSE(LUN)    ! USUALLY STILL OPEN
 
           CALL  FILGET(FILPAT,FILNAM,NLET,NWANT,IRTFLG)
	   IF (IRTFLG .NE. 0) RETURN 

           MAXIM = 0  
	   CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'O',ITYPE,
     &                 NX,NY,NZ, 
     & 		       MAXIM,' ',FOUROK,IRTFLG) 
	   IF (IRTFLG .NE. 0) RETURN 

C          IMAGE SIZE SHOULD BE SAME
           CALL SIZCHK(NULL,NX1,NY1,NZ1,ITYPE1,
     &                      NX ,NY, NZ, ITYPE, IRTFLG)
	   IF (IRTFLG .NE. 0) RETURN 

           !write(6,*)' Opened old templated file: ',filnam(1:nlet),maxim

           NGOT   = NWANT   
           IRTFLG = 0
           RETURN

        ELSEIF (LOCAT  > 0 .AND. LOCAST > LOCAT) THEN
C          TEMPLATED STACKED IMAGE ---------------------------- STK@***

C          MUST LOAD OVERALL HEADER FIRST FOR LUNREDHED (MAY BE MT NOW!)
           CALL LUNGETSIZE(LUN,NX,NY,NZ,IRTFLG)
           CALL LUNREDHED(LUN,NX,0,.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              !write(6,*) ' lunredhed,lun,NX,irtflg:',lun,NX,irtflg
              CALL ERRT(102,'REDHED FAILED ON LUN',LUN) 
              IRTFLG = 2
              RETURN
           ENDIF

C          LOAD SPECIFIED IMAGE HEADER
           CALL LUNREDHED(LUN,NX,NWANT,.FALSE.,IRTFLG)
           IF (IRTFLG == 0) THEN
C             NEED IMUSED FROM THIS STACKED IMAGE
              CALL LUNGETINUSE(LUN,IMUSED,IRTFLG)
           ENDIF
           IF (IRTFLG .NE. 0 .OR. IMUSED == 0) THEN
              CALL ERRT(102,'IMAGE NOT IN STACK',NWANT) 
              IRTFLG = 2
              RETURN
           ENDIF
           NGOT = NWANT

           FILNAM = FILPAT
           CALL LUNSETIMNUM(LUN,FILNAM,NWANT,'O',IRTFLG)
           NLET = lnblnkn(FILNAM)

           !write(6,*)' Opened old templated stacked file: ',FILNAM(:NLET)

        ELSEIF (LOCAT == NLET) THEN
C          WHOLE IMAGE STACK ------------------------------ STK@ or _1@
C          GET SPECIFIED IMAGE HEADER FROM STACK FILE LOCATION
C          DO NOT CALL ERRT IF RUNS OFF END OF FILE

           CALL LUNGETSIZE(LUN,NX,NY,NZ,IRTFLG)

           NGOT = NWANT
           DO
C             LOAD OVERALL HEADER FIRST FOR LUNREDHED (MAY BE MT NOW!)
              CALL LUNREDHED(LUN,NX,0,.FALSE.,IRTFLG)

C             LOAD SPECIFIED IMAGE HEADER
              CALL LUNREDHED(LUN,NX,NGOT,.FALSE.,IRTFLG)
              IF (IRTFLG > 0) THEN
C                PROBABLY RAN OFF END OF STACK FILE
                 IRTFLG = -1
                 NGOT   = 0
                 RETURN
              ENDIF

C             NEED IMUSED FROM THIS STACKED IMAGE
              CALL LUNGETINUSE(LUN,IMUSED,IRTFLG)
              IF (IMUSED > 0) EXIT     ! FOUND NEXT IMAGE

C             THIS IMAGE NOT AN EXISTING IMAGE WITHIN STACK!
C             INCREMENT NGOT AND TRY AGAIN
              NGOT  = NGOT + 1

           ENDDO

           FILNAM = FILPAT
           CALL LUNSETIMNUM(LUN,FILNAM,NGOT,'O',IRTFLG)

           NLET = lnblnkn(FILNAM)

           !write(6,*)' Opened old bare stacked file: ',FILNAM(1:NLET)
           !write(6,*)' ngot,NX:',ngot,NX,lun,imused,irtflg
        ENDIF

C       SET OFFSETS FOR REDLIN/WRTLIN ON THIS LUN
        CALL LUNSETIMGOFF(LUN,NGOT,NX,IRTFLG)

C       WRITE OUT FILE OPENING INFO 
        CALL LUNSAYINFO(LUN,IRTFLG)

C       SET COMMON BLOCK VARIABLES
        CALL LUNSETCOMMON(LUN,IRTFLG)


	END


C **********************************************************************
C
C    GETNEWIMG(LUNCP,LUN,LUNXM,FILPAT,NWANTT, SAYIT,NGOT,IRTFLG)
C
C    PURPOSE:       TO OPEN A SPECIFIED IMAGE WITHIN STACK FOR RANDOM 
C                   ACCESS READING/WRITING.
C
C    PARAMETERS:
C        LUNCP      UNIT NUMBER FOR HEADER TXT COPY               (SENT)
C        LUN        UNIT NUMBER FOR FILNAM.                       (SENT)
C        LUNXM      UNIT NUMBER FOR XM SELFILE                    (SENT)
C        FILPAT     FILENAME PATTERN                              (SENT)
C        NWANT      IMAGE NUMBER WANTED                           (SENT) 
C        SAYIT      SAY FILE OPENING INFO                         (SENT)
C        NGOT       IMAGE NUMBER FOUND                            (RET.) 
C        IRTFLG     ERROR RETURN FLAG.                            (RET.)
C                   IRTFLG = -1    END OF FILE BEFORE NWANT
C                   IRTFLG =  0    NORMAL RETURN, IMAGE IS STACK
C                   IRTFLG =  2    IMAGE NOT IN USE
C
C **********************************************************************

	SUBROUTINE GETNEWIMG(LUNCP,LUN,LUNXM,FILPAT,NWANTT, 
     &                       SAYIT,NGOT,IRTFLG)

        INCLUDE 'CMLIMIT.INC'

        INTEGER                :: LUNCP,LUN,LUNXM
        CHARACTER(LEN=*)       :: FILPAT
        INTEGER                :: NWANTT
        LOGICAL                :: SAYIT
        INTEGER                :: NGOT
        INTEGER                :: IRTFLG

        CHARACTER(LEN=MAXNAM)  :: FILNAM
        CHARACTER(LEN=1)       :: DSP
        LOGICAL                :: FOUROK,ISBARE

        NWANT  = NWANTT          ! NWANTT MAY NOT BE WRITABLE
        NLET   = lnblnkn(FILPAT)
        LOCAST = INDEX(FILPAT(1:NLET),'*')
        LOCAT  = INDEX(FILPAT(1:NLET),'@')

        !write(6,*)' locast,locat,nlet:',locast,locat,nlet,filpat(1:nlet)

        IF (NWANT < 0) THEN
C          XMIPP SELFILE SIMPLE IMAGE ----------------------- SELAAA

C          GET PREVIOUS FILE SIZE AND TYPE (SHOULD BE SAME)
           CALL LUNGETSIZE(LUN,NX1,NY1,NZ1,IRTFLG)
           CALL LUNGETTYPE(LUN,ITYPE1,IRTFLG)
           CLOSE(LUN)           ! USUALLY STILL OPEN
 
C          LOAD FILE NAME FROM SELFILE
           CALL GETNEXT_XMSEL(LUNXM,.TRUE.,FILNAM,NLET,IRTFLG)
	   IF (IRTFLG .NE. 0) RETURN 

C          OPEN FILNAM
           MAXIM = 0  
	   CALL OPFILEC(0,.FALSE.,FILNAM(1:NLET),LUN,'U',ITYPE1,
     &                 NX1,NY1,NZ1, 
     & 		       MAXIM,'~9',FOUROK,IRTFLG) 
	   IF (IRTFLG .NE. 0) RETURN 

           !write(6,*)' Opened new Xmipp selfile file: ',filnam(1:nlet)

           NGOT = NWANT   
           RETURN

        ELSEIF (LOCAT <= 0 .AND. LOCAST > 1) THEN
C          TEMPLATED SIMPLE IMAGE --------------------------- IMG***

C          NEW IMAGE, NEEDS TO KNOW: ITYPE,NX,NY,NZ!
C          GET IT FROM OPFILES OR PREVIOUS CALL
           CALL LUNGETSIZE(LUN,NX,NY,NZ,IRTFLG)
           CALL LUNGETTYPE(LUN,ITYPE,IRTFLG)
           CLOSE(LUN)   ! GET RID OF ANY OPEN IMAGE
 
C          CREATE FILE NAME
           CALL FILGET(FILPAT,FILNAM,NLET,NWANT,IRTFLG)
	   IF (IRTFLG .NE. 0) RETURN 

           MAXIM = 0
           CLOSE(LUN)    ! MAY BE STILL OPEN FROM FIRST CALL  
	   CALL OPFILEC(LUNCP,.FALSE.,FILNAM,LUN,'U',ITYPE,
     &                 NX,NY,NZ, 
     & 		       MAXIM,' ',FOUROK,IRTFLG) 
	   IF (IRTFLG .NE. 0) RETURN 

           !write(6,*)' Opened new templated file: ',FILNAM(1:NLET)

           NGOT   = NWANT   
           IRTFLG = 0
           RETURN

        ELSEIF (LOCAT > 0) THEN
C          STACKED IMAGE ------------------------------- STK@*  or STK@

           CALL LUNGETISBARE(LUN,ISBARE,IRTFLG)
!           IF (ISBARE) THEN
!C             WANT NEXT IMAGE IN STACK, GET CURRENT FILE NUMBER 
!              CALL LUNGETINUSE(LUN,NWANT,IRTFLG)
!              NWANT = NWANT + 1
!              IF (NWANT <= 0) NWANT = 1
!           ENDIF

C          LOAD OVERALL HEADER FIRST FOR LUNREDHED (MAY BE MT NOW!)
           CALL LUNGETSIZE(LUN,NX,NY,NZ,IRTFLG)
           CALL LUNREDHED(LUN,NX,0,.FALSE.,IRTFLG)
 
C          RETRIEVE CURRENT MAXIMUM IMAGE NUMBER FROM OVERALL HEADER
           CALL LUNGETMAXIM(LUN,MAXIM,IRTFLG)

           IF (NWANT > MAXIM) THEN
C             UPDATE OVERALL HEADER WITH MAXIMUM IMAGE NUMBER
              CALL LUNSETMAXIM(LUN,NWANT,IRTFLG)
              CALL LUNSETMAXALL(LUN,NWANT,IRTFLG)
           ENDIF

C          NEED ISTACK 
           CALL LUNCOPYSTK(LUN,ISTACK,IRTFLGT)

          IF (ISTACK < 0) THEN
C             MAKING A NEW INDEXED STACKED FILE, UPDATE INDX LOCATION
              CALL LUNWRTINDX(LUN,NWANT,NX,IRTFLGT)
              IF (IRTFLGT .NE. 0) RETURN
           ENDIF

           IF (NWANT > MAXIM .OR. ISTACK > 2) THEN
C             SAVE OVERALL HEADER NOW TO PRESERVE MAXIM & LASTINDX
              CALL LUNWRTHED(LUN,NX,0,IRTFLGT)
              CALL LUNGETMAXIM(LUN,MAXIM,IRTFLG)
           ENDIF

C          GET FILENAM FROM CURRENT HEADER OBJECT
           CALL LUNGETFILE(LUN,FILNAM,NLET,DSP,IRTFLG)

C          SET NEW FILENAME IN HEADER OBJECT AND GET FILENAME
           CALL LUNSETIMNUM(LUN,FILNAM,NWANT,'N',IRTFLG)
           NLET = lnblnkn(FILNAM)    ! MAY HAVE BEEN ALTERED

C          SET IMGNUM IN HEADER OBJECT
           CALL LUNSETINUSE(LUN,NWANT,IRTFLG)

C          PUT IMAGE ISTACK
           CALL LUNSETISTACK(LUN,0,IRTFLG)

C          PUT COMMON VALUES INTO COMMON AREA (NOT NEEDED IN FUTURE?)
           CALL LUNSETCOMMON(LUN,IRTFLG)

C          PUSH HEADER OBJECT INFO INTO NEW STACKED FILE
           CALL LUNWRTHED(LUN,NX,NWANT,IRTFLG)

C          SET OFFSETS FOR REDLIN/WRTLIN ON THIS LUN
           CALL LUNSETIMGOFF(LUN,NWANT,NX,IRTFLGT)

C          WRITE OUT FILE OPENING INFO TO SCREEN
           CALL LUNSAYINFO(LUN,IRTFLG)

           NGOT = NWANT

           !write(6,*)' Opened new stacked file: ',FILNAM(1:NLET)

        ENDIF

	END




C++*********************************************************************
C
C  NEXTFILES.F  NEW                              12/15/06 ARDEAN LEITH
C               OVERUN OUTPUT LIST = -99          1/15/12 ARDEAN LEITH
C **********************************************************************
C
C NEXTFILES(NINDX1, NINDX2, INUMBR1,INUMBR2, 
C           FOUROK,NGOT1,NGOT2,   
C           MAXIM1,MAXIM2,   
C           LUN1,LUNCP,LUN2, FILPAT1,FILPAT2,
C           IMGNUM1,IMGNUM2,IRTFLG) 
C
C PURPOSE:  GETS NEXT FILES FOR A STACK ORIENTED OPERATION
C           STACKS MUST BE OPENED WITH OPFILES!!!
C
C PARAMETERS: NINDX1,NINDX2    LIST INDICES                   (SENT/RET)
C             INUMBR1,INUMR2   IMAGE NUMBER LISTS                (SENT)
C             FOUROK           FOURIER INPUT IS OK               (SENT)
C             LUNXM1,LUNXM2    LUN FOR SELFILE INPUT             (SENT)
C             NGOT1,NGOT2      # OF IMAGES                       (SENT)
C             MAXIM1,MAXIM2    MAXIM VALUES                      (SENT)
C             LUN1             LUN FOR INPUT  (0 = NO FILE IN)   (SENT)
C             LUNCP            LUN FOR OUTPUT HEADER COPY        (SENT)
C             LUN2             LUN FOR OUTPUT (0 = NO FILE OUT)  (SENT)
C             FILPAT,FILPAT2   FILE NAME PATTERNS                (SENT)
C             IMGNUM1,IMGNUM2  IMAGE NUMBERS                  (SENT/RET)
C             IRTFLG           ERROR (0 IS OK, -1 IS END STACK)  (RET.)
C--*********************************************************************
 
      SUBROUTINE NEXTFILES(NINDX1, NINDX2, INUMBR1,INUMBR2, 
     &                     FOUROK,LUNXM1,LUNXM2,
     &                     NGOT1,NGOT2,     MAXIM1,MAXIM2,   
     &                     LUN1,LUNCP,LUN2, FILPAT1,FILPAT2,
     &                     IMGNUM1,IMGNUM2,IRTFLG) 
 
      IMPLICIT NONE

      INTEGER           :: NINDX1,NINDX2
      INTEGER           :: INUMBR1(NGOT1),INUMBR2(NGOT2)
      LOGICAL           :: FOUROK
      INTEGER           :: LUNXM1,LUNXM2
      INTEGER           :: NGOT1,NGOT2
      INTEGER           :: MAXIM1,MAXIM2,LUN1,LUNCP,LUN2
      CHARACTER(LEN=*)  :: FILPAT1,FILPAT2
      INTEGER           :: IMGNUM1,IMGNUM2,IRTFLG

      INTEGER           :: NWANT1,NWANT2, it
      LOGICAL           :: SAYIT = .TRUE.
      LOGICAL           :: GOTAST1,GOTAST2
      

      NINDX1 = NINDX2 + 1
      NINDX2 = NINDX2 + 1

      !write(6,*) 'nextfiles0  : l1,l2:',lun1,lun2
      IF (LUN1 > 0) THEN  
C        OPEN NEXT INPUT FILE 
         GOTAST1 = (INDEX(FILPAT1,'*') > 0)

         !write(6,'(a,8i5)') '  in:  nindx1,ngot1,maxim1,imgnum1: ',
!     &                              nindx1,ngot1,maxim1,imgnum1

         IF (IMGNUM1 == -1 .AND. LUNXM1 > 0 ) THEN
C           XMIPP SELFILE LISTED IMAGE
            NWANT1 = -1

         ELSEIF (MAXIM1 >= 0 .AND. .NOT. GOTAST1 ) THEN
C           BARE STACK INPUT (NO LIST)
             IF (NINDX1 > NGOT1) THEN
C              FINISHED THE WHOLE STACK
               IRTFLG = -1
               RETURN
            ENDIF
            NWANT1 = NINDX1

         ELSEIF (MAXIM1 == -2  .OR.
     &           MAXIM1 == -1  .OR.
     &           MAXIM1  >    0 ) THEN
C           NON STACKED IMAGE WITH/WITHOUT TEMPLATED LIST
C           STACKED IMAGE WITH/WITHOUT LIST
            IF (NINDX1 > NGOT1) THEN
C              OVERUN INPUT LIST
               IRTFLG = -1
               RETURN
            ENDIF

C           OPEN NEXT INPUT FILE 
            NWANT1 = INUMBR1(NINDX1)

         ENDIF

         CALL GETOLDIMG(LUN1,LUNXM1,FILPAT1,NWANT1,SAYIT, 
     &                  FOUROK,IMGNUM1,IRTFLG)

         !write(6,'(a,8i5)')' Gotoldimg, lun1,nwant1,imgnum1: ',
!     &                                   lun1,nwant1,imgnum1

         IF (IRTFLG .LT. 0) RETURN    ! END OF WHOLE-STACK
         IF (IRTFLG .NE. 0) RETURN    ! ERROR

         IF (MAXIM1 > 0 .AND. .NOT. GOTAST1) THEN
C           INPUT FROM A BARE STACK 
            NINDX1 = IMGNUM1
            NINDX2 = IMGNUM1
         ENDIF
      ENDIF
           ! write(6,*) 'nextfiles1: l1,l2,irtflg:',lun1,lun2,irtflg

      IF (LUN2 > 0) THEN  
C        OPEN NEXT OUTPUT FILE 
         GOTAST2 = (INDEX(FILPAT2,'*') > 0)

         !write(6,'(a,8i5)')'  out: nindx2,ngot2,maxim2,imgnum2: ',
!     &                             nindx2,ngot2,maxim2,imgnum2

!            write(6,*) 'nextfiles2: l2,imgnum2,gotast2:',
!     &                              lun2,imgnum2,gotast2
         IF (IMGNUM2 == -1 .AND. LUNXM2 > 0  ) THEN
C           XMIPP SELFILE LISTED IMAGE
            NWANT2 = -1

         ELSEIF (MAXIM2 >= 0 .AND. .NOT. GOTAST2 ) THEN
C           BARE STACK OUTPUT (NO LIST)
            NWANT2 = NINDX2

            IF (LUN1 > 0) THEN
C              BARE STACK, OUTPUT IMAGE HAS SAME # AS INPUT ALWAYS
               NWANT2 = IMGNUM1
            ENDIF

         ELSEIF (MAXIM2 == -2  .OR.
     &           MAXIM2 == -1  .OR.
     &           MAXIM2  >  0 ) THEN

C           NON-STACK IMAGE WITH/WITHOUT TEMPLATE LIST 
            IF (NINDX2 > NGOT2) THEN
C               OVERUN OUTPUT LIST
                !write(6,*) 'NINDX2 > NGOT2',NINDX2,NGOT2
                IRTFLG = -99
                RETURN
            ENDIF
C           OPEN NEXT OUTPUT FILE 
            NWANT2 = INUMBR2(NINDX2)
         ENDIF

         !write(6,*) ' nextfiles, nwant2,nindx2:',nwant2,nindx2
         !write(6,'(a,8i5)') 
!     &   ' Calling getnew,nwant2,imgnum2,maxim2:',nwant2,imgnum2,maxim2

         CALL GETNEWIMG(LUNCP,LUN2,LUNXM2,FILPAT2,NWANT2,
     &                  SAYIT,IMGNUM2,IRTFLG)

          !write(6,'(a,8i5)') ' Getnew,nwant2,imgnum2,maxim2,irtflg:',
!     &                                nwant2,imgnum2,maxim2,irtflg
         IF (IRTFLG .NE. 0) RETURN     ! ERROR

      ENDIF

      END


C++*********************************************************************
C
C NEXTFILE.F    NEW                              12/15/06 ARDEAN LEITH
C               OVERUN OUTPUT LIST = -99          1/15/12 ARDEAN LEITH
C
C **********************************************************************
C
C NEXTFILE  
C
C PURPOSE:  GETS NEXT FILE FOR A STACK ORIENTED OPERATION
C           STACKS MUST BE OPENED WITH OPFILES!!!
C
C PARAMETERS: NINDX1           LIST INDICES                   (SENT/RET)
C             INUMBR1          IMAGE NUMBER LISTS                (SENT)
C             FOUROK           FOURIER INPUT IS OK               (SENT)
C             LUNXM1           LUN FOR SELFILE INPUT             (SENT)
C             NGOT1            # OF IMAGES                       (SENT)
C             MAXIM1           MAXIM VALUES                      (SENT)
C             LUN1             LUN FOR I/0                       (SENT)
C             LUNCP            LUN FOR OUTPUT HEADER COPY        (SENT)
C             FILPAT1          FILE NAME PATTERN                 (SENT)
C             DISP             IMAGE EXISTANCE                   (SENT)
C             IMGNUM1          IMAGE NUMBER                   (SENT/RET)
C             IRTFLG           ERROR (0 IS OK, -1 IS END STACK)  (RET.)
C--*********************************************************************
 
      SUBROUTINE NEXTFILE(NINDX1,   INUMBR1, 
     &                    FOUROK,   LUNXM1,
     &                    NGOT1,    MAXIM1,   
     &                    LUN1,     LUNCP, 
     &                    FILPAT1,  DISP,
     &                    IMGNUM1,  IRTFLG) 
 
      IMPLICIT NONE

      LOGICAL           :: ISIN
      INTEGER           :: NINDX1 
      INTEGER           :: INUMBR1(NGOT1) 
      LOGICAL           :: FOUROK
      INTEGER           :: LUNXM1 
      INTEGER           :: NGOT1 
      INTEGER           :: MAXIM1, LUN1,LUNCP 
      CHARACTER(LEN=*)  :: FILPAT1
      INTEGER           :: IMGNUM1
      CHARACTER(LEN=1)  :: DISP
      INTEGER           :: IRTFLG

      INTEGER           :: NWANT1, it
      LOGICAL           :: SAYIT = .TRUE.
      LOGICAL           :: GOTAST
      

      NINDX1 = NINDX1 + 1
      GOTAST = (INDEX(FILPAT1,'*') > 0)

      IF (DISP == 'O' .OR. DISP == 'B' .OR. 
     &    DISP == 'Z' .OR. DISP == 'E') THEN 
 
C        OPEN INPUT FILE 
         !write(6,'(a,8i5)') '  in:  nindx1,ngot1,maxim1,imgnum1: ',
c     &                              nindx1,ngot1,maxim1,imgnum1

         IF (IMGNUM1 == -1  .AND. LUNXM1 > 0 ) THEN
C           XMIPP SELFILE LISTED IMAGE
            NWANT1 = -1

         ELSEIF (MAXIM1 >= 0 .AND. .NOT. GOTAST ) THEN
C           BARE STACK INPUT (NO LIST)
             IF (NINDX1 > NGOT1) THEN
C              FINISHED THE WHOLE STACK
               IRTFLG = -1
               RETURN
            ENDIF
            NWANT1 = NINDX1

         ELSEIF (MAXIM1 == -2  .OR.
     &           MAXIM1 == -1  .OR.
     &           MAXIM1  >    0 ) THEN
C           NON STACKED IMAGE WITH/WITHOUT TEMPLATED LIST
C           STACKED IMAGE WITH/WITHOUT LIST
            IF (NINDX1 > NGOT1) THEN
C              OVERUN INPUT LIST
               IRTFLG = -1
               RETURN
            ENDIF

C           OPEN NEXT INPUT FILE 
            NWANT1 = INUMBR1(NINDX1)

         ENDIF
         !write(6,*) ' call gotoldimg, nwant1:',nwant1,nindx1
 
         CALL GETOLDIMG(LUN1,LUNXM1,FILPAT1,NWANT1,SAYIT, 
     &                  FOUROK,IMGNUM1,IRTFLG)

         !write(6,'(a,8i5)')' Gotoldimg, ngot1,nwant1,imgnum1: ',
c     &                           ngot1,nwant1,imgnum1,irtflg,maxim1

         IF (IRTFLG .LT. 0) RETURN    ! END OF WHOLE-STACK
         IF (IRTFLG .NE. 0) RETURN    ! ERROR

         IF (MAXIM1 > 0 .AND. .NOT. GOTAST) THEN
C           BARE STACK 
            NINDX1 = IMGNUM1
         ENDIF

      ELSE   
C        OPEN NEXT OUTPUT FILE 
         !write(6,'(a,8i5)')'  out: nindx1,ngot1,maxim1,imgnum1: ',
c     &                            nindx1,ngot1,maxim1,imgnum1

         IF (IMGNUM1 == -1  .AND. LUNXM1 > 0 ) THEN
C           XMIPP SELFILE LISTED IMAGE
            NWANT1 = -1

         ELSEIF (MAXIM1 >= 0 .AND. .NOT. GOTAST ) THEN
C           BARE STACK OUTPUT (NO LIST)
            NWANT1 = NINDX1

         ELSEIF (MAXIM1 == -2  .OR.
     &           MAXIM1 == -1  .OR.
     &           MAXIM1  >  0 ) THEN
C           NON STACKED IMAGE WITH/WITHOUT TEMPLATED LIST
C           STACKED IMAGE WITH/WITHOUT LIST
            IF (NINDX1 > NGOT1) THEN
C               OVERUN OUTPUT LIST
                !write(6,*) 'NINDX1 > NGOT1',NINDX1,NGOT1
                IRTFLG = -99
                RETURN
            ENDIF
C           OPEN NEXT OUTPUT FILE 
            NWANT1 = INUMBR1(NINDX1)
         ENDIF

         !write(6,*) ' Nextfiles, nwant1,nindx1:',nwant1,nindx1
         !write(6,'(a,8i5)') 
c     &   ' calling Getnew,nwant1,imgnum1,maxim1:',nwant1,imgnum1,maxim1

         CALL GETNEWIMG(LUNCP,LUN1,LUNXM1,FILPAT1,NWANT1,
     &                  SAYIT,IMGNUM1,IRTFLG)

         !write(6,'(a,8i5)') ' Getnew,nwant1,imgnum1,maxim1,irtflg:',
!     &                              nwant1,imgnum1,maxim1,irtflg
         IF (IRTFLG .NE. 0) RETURN     ! ERROR

         IF (MAXIM1 > 0 .AND. .NOT. GOTAST) THEN
C           BARE STACK, OUTPUT IMAGE HAS SAME # AS INPUT ALWAYS
            NINDX1 = IMGNUM1
         ENDIF
      ENDIF

      END


