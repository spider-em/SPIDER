 
C++*********************************************************************
C                                                                      
C  OPFILES.F NEW                              12/15/06  ArDean Leith
C            BAD NUMBRT() TRAP                05/21/09  ArDean Leith 
C            ASKNAM, PROMPTEX                 12/06/10  ArDean Leith 
C            NX...                            03/26/12  ArDean Leith 
C            COPY NON SPIDER INPUT            05/26/14  ArDean Leith 
C            COPY NON SPIDER INPUT            05/26/14  ArDean Leith 
C            LOCAST, ASKLIST FOR ILIST        10/02/14  ArDean Leith
C            DO NOT CHECK .MRC FOR XMIPP       7/25/19  ArDean Leith
C            MRC SUPPORT                       8/05/19  ArDean Leith
C            COMMENTS                          2/06/20  ArDean Leith
C            DEBUG OUTPUT ADDED                9/26/25  ArDean Leith
C            LIST_SIZE SET 0                   9/26/25  ArDean Leith 
C            ADDED CALL TREE                  10/16/25  ArDean Leith
C            PUT NEXTFILE(S) IN nextfiles.f   12/13/25  ArDean Leith
C            RENAMED NTOT -> LIST_SIZE        12/14/25  ArDean Leith
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
C  CONTAINS: OPFILES, GETOLDIMG, GETNEWIMG
C            PUT NEXTFILE(S) IN: nextfiles.f  Dec 2025 

C  OPFILES(LUNCP,LUNIMG,LUNDOC,LUNXM,  
C          ASKNAM,FILPAT,NLET, DISP,
C          ITYPE,NX,NY,NZ,MAXIM, PROMPT,
C          FOUROK, ILIST,NIMAXT, UNUSED
C          LIST_SIZE,IMGNUM, IRTFLG)
C 
C  PURPOSE: SOLICITS FILE NAME(S) AND OPENS FILE(S)
C           SUPPORT ROUTINE FOR CONVERTING OPERATIONS TO 
C           WORK ON WHOLE STACK OR WITH A SELECTION DOC FILE.
C           ONLY USED IN A FEW OPERATIONS.  RAN OUT OF YEARS TO
C           COMPLETE THE IMPLEMENTATION IN MORE OPERATIONS
C           ALSO WORKS FOR MRC FILES NOW
C
C  PARAMETERS:
C        LUNCP     UNIT TO COPY HEADER VALUES FROM               (SENT)
C        LUNIMG    UNIT TO OPEN FILE ON                          (SENT)
C        LUNDOC    UNIT TO OPEN LIST DOC FILES ON                (SENT)
C        LUNXM     UNIT TO OPEN XMIPP SELFILE ON                 (SENT)
C        ASKNAM    FLAG TO ASK FOR FILE NAME                     (SENT)
C        FILPAT    FILE NAME PATTERN                              (RET)
C        NLET      CHARS IN FILE NAME PATTERN                     (RET)
C
C        DISP      CHARACTER CONTAINING ONE OF THE               (SENT)
C                  FOLLOWING DISPOSITION SPECIFICATIONS:
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
C
C        ITYPE     IFORM FOR FILE                        (SENT OR RET) 
C
C        NX,NY,NZ  IMG/VOL SIZE                          (SENT OR RET)
C
C        MAXIM     STACK INDICATOR OR NUMBER FOR IMAGES     (SENT/RET)
C                   ON INPUT: STACK INDICATOR IF DISP == 'I'    (SENT)
C
C                   ON OUTPUT:                                   (RET)
C                       -2 A NON-STACK IMAGE                
C                       -1 A SINGLE STACKED IMAGE                  
C                     >= 0 NSTK:  CURRENT MAX. IMAGE NO. IN STACK 
C                   IF NOT STACK MAXIM IS NUMBER IN LIST
C            
C        PROMPT    PROMPT FOR FILNAME                           (SENT)
C                     IF NOT (ASKNAM) THIS IS FILE NAME!         
C                     AT END: ~  SKIPS "FILE" ON PROMPT
C                             ~6 KEEPS OLD DATE/TIME
C                             ~7 CAN OPEN A STACK WITHOUT @
C                             ~8 SETS OPENED FILE AS READ ONLY (R)
C                             ~9 KEEPS INCOMING EXTENSION
C                               (OTHERWISE DISCARDED FROM FILNAM)
C
C        FOUROK    CAN USE EXISTING FOURIER FILES?              (SENT)
C
C        ILIST     IMAGE NUMBER LIST                             (RET)
C                     IF NIMAXT < 0 MUST BE SENT
C                     NOT USED IF SINGLE IMAGE or SELFILE 
C
C        NIMAXT    MAX LENGTH OF IMAGE NUMBER LIST              (SENT)
C                     <0 MEANS DO NOT ASK FOR LIST
C
C        UNUSED    UNUSED PARAMETER                                (?)
C
C        LIST_SIZE: # OF IMAGES IN IMAGE NUMBER LIST             (RET)
C                   ZERO FOR SINGLE IMAGE AND NO *
C 
C        IMGNUM    IMAGE NUMBER THAT IS CURRENTLY OPEN      (SENT/RET)
C                     ON INPUT:  IF (BARESTACK) IS # WANTED
C                     ON OUTPUT: <0 IS SELFILE IN USE 
C
C        IRTFLG    ERROR FLAG (0 IS NORMAL)                      (RET)
C                      -1 GOTO PREVIOUS QUESTION
C                      -5 NOT A SPIDER OR MRC FILE
C
C  CALL TREE:
C
C     OPFILES Uses (E.G. STK@*)~~9  prompt for filename
C        v    '~9' in last 3 char. of prompt accepts ext. on filename 
C        |
C        |
C        ` -> FILERD  
C        |       ` ->  ECHONAME
C 
C        ` -> FILELIST
C        |
C        |    MRC file:  *.MRC or *.MRCS
C        ` -> OPFILES_MRC
C        | 
C        |    Templated stack: stk@* 
C        ` -> FILNAMANDEXT --> OPFILEC --> GETOLDIMG
C        |                             --> GETNEWIMG
C        |    Whole barestack: stk@     
C        ` -> FILNAMANDEXT --> OPFILEC --> GETOLDIMG
C        |                             --> GETNEWIMG
C        |    File template:   img***  
C        ` -> FILGET      -->  OPFILEC  
C        | 
C        |  Simple file:       img001
C        ` -> ---------------> OPFILEC  
C        | 
C        
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--********************************************************************* 
       SUBROUTINE OPFILES(LUNCP,LUNIMG,LUNDOC,LUNXM,
     &                    ASKNAM,FILPAT,NLET, DISP,
     &                    ITYPE,NX,NY,NZ, MAXIM,
     &                    PROMPT,
     &                    FOUROK,ILIST,NIMAXT, 
     &                    UNUSED,LIST_SIZE,IMGNUM, IRTFLG) 
 
       IMPLICIT NONE
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
       INTEGER                   :: NIMAXT, UNUSED,LIST_SIZE
       INTEGER                   :: IMGNUM, IRTFLG

       INTEGER                   :: ISTK,LIST_VAL,NSTK

       INTEGER                   :: LOCAT,LOCAST
       INTEGER                   :: LUNOP,NLETT,LOCTILDE
       INTEGER                   :: IMGWANT,IMGNUMIN,LENE,I,IFLIPIN
       INTEGER                   :: NIMAXP,LUNXM,NDUM
       LOGICAL                   :: ISOPEN,GOTFILE
       LOGICAL                   :: SAYIT,ASKLIST

       CHARACTER (LEN=MAXNAM)    :: FILNAMT,FILNAM 
       CHARACTER (LEN=2*MAXNAM)  :: MESG 
       CHARACTER (LEN=100)       :: PROMPTEX
       CHARACTER (LEN=1)         :: CDUM 
       CHARACTER (LEN=1)         :: DISPT 
       CHARACTER (LEN=1)         :: NULL = CHAR(0)

       LOGICAL                   :: IS_MRC, IS_BARE
       LOGICAL                   :: ISMRCFILE    ! FUNCTION 

       INTEGER                   :: lnblnkn      ! FUNCTION
  
C      LOCTILDE = INDEX(PROMPT,'~') ! unfinished
C      LENP = lnblnkn(PROMPT)
C      PUT ~7 AT END OF PROMPT TO OPEN* TO OPEN A STACK WITHOUT @
C      OPSTKNOAT =  (INDEX (PROMPT(1:LENP),'~7') > 1)
C      PUT ~9 AT END OF PROMPT TO OPEN* TO KEEP INCOMING EXTENSION
C      KEEPEXT = (INDEX(PROMPT(1:LENP),'~9') > 1)


       IF (ASKNAM .AND. PROMPT == NULL) THEN
C         ASK FOR FILE NAME, CAN ACCEPT EXTENSION

          IF (DISP == 'N' .OR. 
     &        DISP == 'I' .OR.
     &        DISP == 'U') THEN
C            NEW FILE, USE DEFAULT PROMPT
             PROMPTEX = 
     &          'OUTPUT FILE NAME OR TEMPLATE (E.G. stk-file@*)~'
          ELSE
C            OLD FILE, USE DEFAULT PROMPT
             PROMPTEX = 
     &          'INPUT FILE NAME OR TEMPLATE (E.G. stk-file@*)~'
          ENDIF
     
          CALL FILERD(FILPAT,NLET,NULL,PROMPTEX,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN 

#if defined (SP_DBUGIO)
          write(3,*)' In opfiles 111; filpat: ',trim(filpat)
#endif

       ELSEIF (ASKNAM) THEN
C         ASK FOR FILE NAME USING RECEIVED PROMPT
      
          CALL FILERD(FILPAT,NLET,NULL,PROMPT,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN 

#if defined (SP_DBUGIO)
        write(3,*)' In opfiles 222; prompt: ',trim(prompt)
        write(3,*)' In opfiles 222; filpat: ',trim(filpat)
#endif

       ELSE
C         USE FILENAME SENT IN: PROMPT

          NLET   = LNBLNKN(PROMPT)
          FILPAT = TRIM(PROMPT(1:NLET))

#if defined (SP_DBUGIO)
          write(3,*)' In opfiles 333; filpat: ',trim(filpat)
#endif
       ENDIF
 
       LOCAT    = INDEX(FILPAT(1:NLET),'@')   
       LOCAST   = INDEX(FILPAT(1:NLET),'*')
     
       ASKLIST = (NIMAXT > 0)   ! ASK FOR IMAGE NUMBER(S)
       NIMAXP  = ABS(NIMAXT)    ! ILIST DIMENSION, IF ILIST IN USE
 
#if defined (SP_DBUGIO)
       !write(3,*)' In opfiles 444; filpat:       ', trim(filpat)
        write(3,*)' In opfiles 444; locast,locat: ', locast,locat
       !write(3,*)' In opfiles 444; nimaxt:       ', nimaxt
        write(3,*)' In opfiles 444, asklist:      ', asklist 
#endif
       
       IMGWANT   = 1   ! DEFAULT IMAGE NUMBER 
       LIST_SIZE = 0   ! DEFAULT NUMBER UNUSED        
       
       IF (LOCAST > 0 .AND. ASKLIST) THEN
C         GET LIST OF IMAGES FROM DOC. FILE OR INPUT LINE
          LIST_SIZE   = 0 
          CALL FILELIST(.FALSE.,LUNDOC,CDUM,NDUM,ILIST,NIMAXP,
     &                  LIST_SIZE,' ',IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

C         START WITH FIRST FILE IN SERIES
          IMGWANT = ILIST(1)

#if defined (SP_DBUGIO)
          write(3,*)' In opfiles aaa, imgwant:   ',imgwant
          write(3,*)' In opfiles aaa, list_size: ',list_size
#endif

       ELSEIF (LOCAST > 0 .AND. NIMAXT < 0) THEN

C         ILIST CONTAINS TEMPLATED SERIES, 
C         START WITH FIRST FILE IN THIS SERIES
          IMGWANT   = ILIST(1)
          LIST_SIZE = NIMAXP
 
       ELSEIF (.NOT. ASKLIST) THEN
          LIST_SIZE = NIMAXP

       ELSEIF (LOCAST == 0 .AND. ASKLIST .AND. 
     &        (IMGWANT < 0 .OR.  IMGWANT > 10000000)) THEN
C         SIMPLE FILE NAME 
          IMGWANT = 0
       ENDIF


       IF (IMGWANT < 0 .OR. IMGWANT > 10000000) THEN
          CALL ERRT(102,'INVALID IMAGE NUMBER',IMGWANT)
          IRTFLG = 1
          GOTO 9000
       ENDIF

C      SEE IF THIS IS A MRC FILE
       IS_MRC = ISMRCFILE(FILPAT)

#if defined (SP_DBUGIO)
       write(3,*)'   '
       write(3,*)' In opfiles 555; filpat: ',trim(filpat)
       write(3,*)' In opfiles 555; is_mrc:         ', is_mrc
       write(3,*)' In opfiles 555; locast:         ', locast
       write(3,*)' In opfiles 555; asklist,nimaxp: ', asklist,nimaxp 
       write(3,*)' In opfiles 555; list_size,ilist(:2): ',
     &                             list_size,ilist(1:2)
       write(3,*)' In opfiles 555; imgnum,imgwant: ', imgnum,imgwant
#endif


       IF (IS_MRC) THEN
C         OPEN MRC FILE ----------------------------

#if defined (SP_DBUGIO)
          write(3,*)' In opfiles 666 ; Calling opfiles_mrc -----'
#endif

          LIST_VAL  = IMGWANT
          CALL OPFILES_MRC(LUNCP,LUNIMG,LUNDOC,
     &                     FILPAT,NLET, DISP,
     &                     ITYPE,NX,NY,NZ, NSTK, 
     &                     ILIST, LIST_SIZE, LIST_VAL,  
     &                     IRTFLG)
 
          MAXIM  = NSTK       ! RETURNED TO CALLER HERE
          IMGNUM = LIST_VAL   ! RETURNED TO CALLER HERE
C         LIST_SIZE             RETURNED TO CALLER HERE
 
#if defined (SP_DBUGIO)
          write(3,*)' Opfiles 777 returns: nz:              ', nz 
          write(3,*)' Opfiles 777 returns: maxim=nstk:      ', maxim 
          write(3,*)' Opfiles 777 returns: imgnum=list_val: ', imgnum 
          write(3,*)' Opfiles 777 returns: list_size:       ',list_size 
          write(3,*)' Opfiles 777 returns: irtflg:          ', irtflg 
          write(3,*)'   '
                                            
#endif

          RETURN

       ENDIF  ! -------------- End of mrc --------------------


       IMGNUMIN = IMGNUM
       IMGNUM   = 0 

       IF (LOCAT > 0 .AND. LOCAST > LOCAT) THEN
C         TEMPLATED STACKED FILE: STK@* -------------- _9@* or STK@*

          FILNAM = FILPAT(1:LOCAT)
            
          IF (FILNAM(1:1) .NE. '_') THEN
C            CONCATENATE EXTENSION ONTO FILNAM
             CALL FILNAMANDEXT(FILNAM(1:LOCAT-1),DATEXC,
     &                         FILNAMT,NLET,.TRUE.,IRTFLG)
             INQUIRE(FILE=FILNAMT,OPENED=ISOPEN,NUMBER=LUNOP)
             MESG = '  FILE: ' // FILNAMT(1:NLET) //
     &              '  ALREADY OPENED ON' // NULL
             LENE = LNBLNKN(MESG)
             IF (ISOPEN .AND. LUNOP .NE. LUNIMG) THEN
                WRITE(NOUT,'(A,I3)') MESG(1:LENE),LUNOP 
                IRTFLG = -2
                GOTO 9000
             ENDIF
          ENDIF

C         OPEN THE STACK FILE HEADER 
          MAXIM  = 1  
          CALL OPFILEC(LUNCP,.FALSE.,FILNAM,LUNIMG,DISP,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,' ',FOUROK,IRTFLG) 
          IF (IRTFLG .NE. 0) GOTO 9000 

          !write(3,*) ' In opfiles - locast: ',locat,locast,filnam
          !write(3,*) ' In opfiles - imgwant,imgnum2: ',imgwant,imgnum

C         OPEN FIRST FILE IN STACK SERIES
          IMGWANT  = ILIST(1)
          IF (IMGWANT < 0 .OR. IMGWANT > 10000000) THEN
              CALL ERRT(102,'INVALID IMAGE NUMBER',IMGWANT)
              IRTFLG = 1
              GOTO 9000
          ENDIF

C         THIS IS NOT A BARE STACK FILE
          CALL LUNSETISBARE(LUNIMG,.FALSE.,IRTFLG)

          SAYIT    = .TRUE.
          IF (DISP == 'U' .OR. DISP == 'N') THEN
C            WANT A NEW SPIDER STACK FILE FOR OUTPUT

             !write(3,*)' In opfiles  imgwant,imgnum: ',imgwant,imgnum

             CALL GETNEWIMG(LUNCP,LUNIMG,LUNDOC,FILPAT,IMGWANT,
     &                      SAYIT,IMGNUM,IRTFLG)

             !write(3,*) ' In opfiles mgwant,imgnum3:',imgwant,imgnum

          ELSE
C            WANT TO USE EXISTING SPIDER STACK FILE
             CALL GETOLDIMG(LUNIMG,LUNDOC,FILPAT, IMGWANT,SAYIT,
     &                         FOUROK,IMGNUM,IRTFLG)
          ENDIF
          IF (IRTFLG .NE. 0) GOTO 9000

C         RETRIEVE CURRENT MAXIMUM IMAGE NUMBER FROM OVERALL HEADER
          CALL LUNGETMAXIM(LUNIMG,MAXIM,IRTFLG)

          !write(3,'(a,a,a,i6,a,i6,a,i6)')' opened templated stack: ',
          !&      filpat(1:nlet),' at image: ',imgnum,'  maxim: ',maxim


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
     &                 MAXIM,'INPUT',FOUROK,IRTFLG) 
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
             LIST_SIZE = 0
             DO I= 1,MAXIM
                LIST_SIZE = LIST_SIZE + 1
                IF (LIST_SIZE > NIMAXP) THEN
                   CALL ERRT(102,'IMAGE # LIST OVERFLOW AT IMAGE',
     &                            LIST_SIZE)
                   GOTO 9000
                ENDIF
                ILIST(LIST_SIZE) = I
             ENDDO
          ENDIF

          !write(3,*)' Opened bare stack:',filpat(1:nlet),' img:',imgnum
          !write(3,*)' Opened imgnum, stack:',imgnum,filpat(1:nlet)
          !write(3,*)' Opened imgwant,list_size:',  
          !                   imgwant,imgnum,list_size


       ELSEIF (LOCAST > 0) THEN
C         SIMPLE TEMPLATED FILE: IMG*** ---------------------- IMG***
C         NOT A TEMPLATED STACK 

C         SUBSTITUTE IMGWANT INTO FILPAT 
          CALL FILGET(FILPAT,FILNAM,NLET,IMGWANT,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9000 

C         OPEN FIRST FILE IN THE SERIES
          MAXIM = 0 
          CALL OPFILEC(LUNCP,.FALSE.,FILNAM,LUNIMG,DISP,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,' ',FOUROK,IRTFLG) 
          IF (IRTFLG .NE. 0) GOTO 9000 

#if defined(SP_DBUGIO)
          write(3,*)' '
          write(3,*)' In opfiles bbb; imgwant,list_size: ',  
     &                                imgnum,list_size
          write(3,*)' In opfiles bbb; filpat: ', trim(filpat)
          write(3,*)' In opfiles bbb; filnam: ', trim(filnam)
          write(3,*)' '
#endif


       ELSE
C         SINGLE INPUT FILE: IMG001 ------------------------ IMG001
C         SINGLE STACKED FILE: IMGSTK ---------------------- IMGSTK@1
C         OR XMIPP SELFILE LISTING FILE: SELX -------------- SELFILE
C         OR TRY TO COPY NON-SPIDER FILE   ----------------- NONSPIFILE

#if defined(SP_DBUGIO)
          !write(3,*) ' In opfiles ; Simple file: ',filpat(:nlet)
#endif
               
C         CHECK FOR XMIPP SELFILE LIST
          IF (LUNXM > 0 .AND. .NOT. ISMRCFILE(FILPAT) ) THEN
             !write(3,*)' Filpat for openxmsel: ',filpat(:nlet) 
             CALL OPENXMSELFILE(FILPAT,LUNXM,FILNAM,NLET,
     &                          LIST_SIZE,IRTFLG)
#if defined(SP_DBUGIO)
             !write(3,*)' Filnam from openxmsel: ',list_size,':',
             ! &  filnam(:nlet) 
#endif
             INQUIRE(FILE=FILNAM(:NLET),EXIST=GOTFILE)
           
             IF (LIST_SIZE > 0 .AND. GOTFILE) THEN
C               OPEN FIRST FILE IN XMIPP SELFILE LIST
                MAXIM = 0  
                CALL OPFILEC(LUNCP,.FALSE.,FILNAM(:NLET),LUNIMG,DISP,
     &                       ITYPE,NX,NY,NZ, 
     &                              MAXIM,'dum~9',FOUROK,IRTFLG) 
                IF (IRTFLG .NE. 0) GOTO 9000 

                IMGNUM     = -1
                !write(3,*)' Opened selfile image: ',filnam(1:nlet) 
                RETURN

             ENDIF
          ENDIF

C         SINGLE SIMPLE  INPUT FILE: IMG001   ----------------- IMG001
C         SINGLE STACKED INPUT FILE: IMGSTK@3 ----------------- IMGSTK@3

#if defined(SP_DBUGIO)
           write(3,*)' In opfiles, filpat,',trim(filpat)
           write(3,*)' In opfiles, calling opfilec,  disp: ',disp
           write(3,*)'  '
#endif

          MAXIM = 0  
          CALL OPFILEC(LUNCP,.FALSE.,FILPAT,LUNIMG,DISP,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,PROMPT,FOUROK,IRTFLG) 

#if defined(SP_DBUGIO)
          write(3,*)'   '
         !write(3,*)' Back in opfiles, opened:       ', trim(filpat)
          write(3,*)' Back in opfiles, maxim:        ', maxim
          write(3,*)' Back in opfiles, itype,irtflg: ', itype,irtflg
#endif
 
C         RETURN FILENAME WITH ANY EXTENSION IF NOT SPIDER IMAGE
          IF (IRTFLG == 5) NLET = lnblnkn(FILPAT)
          IF (IRTFLG .NE. 0) GOTO 9000 

          LIST_SIZE  = 0
          IMGNUM     = 1

#if defined(SP_DBUGIO_NEVER)
          !write(3,*)' In opfiles: opened simple:  ', trim(filpat)
          !write(3,*)' In opfiles: itype,nimaxt:   ', itype,nimaxt 
          !find native-enddedness of input     
          !call lungetflip(luncp,iflipin,irtflg)
          !write(3,*)' In opfiles: luncp,iflipin: ', luncp,iflipin
#endif
 
      ENDIF

#if defined (SP_DBUGIO)
      write(3,*)' Opfiles 999 returns: maxim =nstk:     ', maxim 
      write(3,*)' Opfiles 999 returns: imgnum=list_val: ', imgnum 
      write(3,*)' Opfiles 999 returns: list_size:       ', imgnum 
      write(3,*)' Opfiles 999 returns: irtflg:          ', irtflg 
      write(3,*)' Opfiles 999 for SPIDER  ----------------------------'
      write(3,*)'    ' 
#endif

9000  RETURN

      END 


C++*********************************************************************
C
C GETOLDIMG.F   FROM GETNXTSTK                     JAN 02 ArDean Leith
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
C        LUNXM      LUN FOR XM SELFILE                            (SENT) C
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
C    CALL TREE:
C
C       GETOLDIMG 
C          |
C          |  Templated simple image:  IMG*** 
C          ` ---> FILGET --> OPFILEC  
C          |                         
C          |  Templated stacked image: STK@*     
C          ` ---> LUN***
C          | 
C          |  Whole image stack: STK@  
C          ` ---> LUNS***
C          | 
C          |  Templated stacked MRC: *@MRC or *@.MRCS
C          ` ---> GETOLDIMG_MRC --> OPFILEC  
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
        LOGICAL                :: IS_MRC

        NLET   = lnblnkn(FILPAT)

        LOCAST = INDEX(FILPAT(1:NLET),'*')
        LOCAT  = INDEX(FILPAT(1:NLET),'@')

C       IS THIS A MRC FILE SET
        CALL LUNGETIS_MRC(LUN,IS_MRC,IRTFLG)

        !write(3,*)' locast,locat:',locast,locat,nlet,filpat(1:nlet)
        !write(3,*)' getoldimg, nwant,: ',nwant,':',filpat(1:nlet)
        
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
     &                        MAXIM,'~9',FOUROK,IRTFLG) 
           IF (IRTFLG .NE. 0) RETURN 

C          NEW IMAGE SIZE SHOULD BE SAME AS PREVIOUS FILE
           CALL SIZCHK(NULL,NX1,NY1,NZ1,ITYPE1,
     &                      NX ,NY, NZ, ITYPE, IRTFLG)
           IF (IRTFLG .NE. 0) RETURN 
 
           !write(3,*)' Opened old xmipp selfile file: ',filnam(1:nlet)

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
     &                        MAXIM,' ',FOUROK,IRTFLG) 
           IF (IRTFLG .NE. 0) RETURN 

C          IMAGE SIZE SHOULD BE SAME
           CALL SIZCHK(NULL,NX1,NY1,NZ1,ITYPE1,
     &                      NX ,NY, NZ, ITYPE, IRTFLG)
           IF (IRTFLG .NE. 0) RETURN 

           !write(3,*)' Opened old templated file: ',filnam(1:nlet), maxim

           NGOT   = NWANT   
           IRTFLG = 0
           RETURN

        ELSEIF (IS_MRC .AND. LOCAT > 0 .AND. LOCAST < LOCAT) THEN

C          TEMPLATED STACKED MRC IMAGE ---------------------- **@STK.mrc
C          BARE STACKED MRC IMAGE      ----------------------   @STK.mrc
           
           !write(3,*)' In getoldimg, nwant: ',nwant,locat,locast,filpat
           CALL GETOLDIMG_MRC(LUN,FILPAT,NWANT,SAYIT,
     &                        FILNAM,NGOT,IRTFLG)
           RETURN

        ELSEIF (LOCAT  > 0 .AND. LOCAST > LOCAT) THEN
C          TEMPLATED STACKED SPIDER IMAGE --------------------- STK@*

C          MUST LOAD OVERALL HEADER FIRST FOR LUNREDHED (MAY BE MT NOW!)
           CALL LUNGETSIZE(LUN,NX,NY,NZ,IRTFLG)
           CALL LUNREDHED(LUN,NX,0,.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              !write(3,*)' In opfiles lunredhed; lun,nx,irtflg:',
c    &                                           lun,nx,irtflg
              CALL ERRT(102,'REDHED FAILED ON LUN',LUN) 
              IRTFLG = 2
              RETURN
           ENDIF

C          LOAD SPECIFIED IMAGE HEADER
           CALL LUNREDHED(LUN,NX,NWANT,.FALSE.,IRTFLG)
           IF (IRTFLG == 0) THEN
C             NEED IMUSED FROM THIS STACKED IMAGE
              CALL LUNGETINUSE(LUN,IMUSED,IRTFLG)
              !write(3,*)' In opfiles lungetinuse, lun,imused:',
c     &                                            lun,imused
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

           !write(3,*)' Opened old templated stacked file: ',
           !filnam(:nlet)

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

#if defined (SP_DBUGIO)
           write(3,*)' In getoldimg; CCCC ', filnam(1:nlet)
           write(3,*)' In getoldimg; locat, nlet: ', locat, nlet
#endif
        ENDIF

#if defined (SP_DBUGIO)
           write(3,*)' In getoldimg; Opened filnam(1:nlet): ',
     &                                      filnam(1:nlet)
           write(3,*)' In getoldimg; locat,nlet,ngot: ', 
     *                               locat,nlet,ngot
#endif
 
C       SET OFFSETS FOR REDLIN/WRTLIN ON THIS LUN
        CALL LUNSETIMGOFF(LUN,NGOT,NX,IRTFLG)

C       WRITE OUT FILE OPENING INFO 
        CALL LUNSAYINFO(LUN,IRTFLG)

C       SET COMMON BLOCK VARIABLES
        CALL LUNSETCOMMON(LUN,IRTFLG)

        END




C **********************************************************************
C
C  GETNEWIMG(LUNCP,LUN,LUNXM,FILPAT,NWANTT, SAYIT,NGOT,IRTFLG)
C
C  PURPOSE:       OPEN A SPECIFIED IMAGE OR WITHIN STACK FOR RANDOM 
C                 ACCESS READING/WRITING.
C
C  PARAMETERS:
C       LUNCP      UNIT NUMBER FOR HEADER TXT COPY               (SENT)
C       LUN        UNIT NUMBER FOR FILNAM.                       (SENT)
C       LUNXM      UNIT NUMBER FOR XM SELFILE                    (SENT)
C       FILPAT     FILENAME PATTERN                              (SENT)
C       NWANT      IMAGE NUMBER WANTED                           (SENT) 
C       SAYIT      SAY FILE OPENING INFO                         (SENT)
C       NGOT       IMAGE NUMBER FOUND                            (RET.) 
C       IRTFLG     ERROR RETURN FLAG.                            (RET.)
C                  IRTFLG = -1    END OF FILE BEFORE NWANT
C                  IRTFLG =  0    NORMAL RETURN, IMAGE IS STACK
C                  IRTFLG =  2    IMAGE NOT IN USE
C
C  CALL TREE:
C
C       GETNEWIMG 
C          |
C          |  Templated simple image:  IMG*** 
C          ` ---> FILGET --> OPFILEC  
C          |                         
C          |  Templated stacked image: STK@*     
C          ` ---> LUN***
C          | 
C          |  Whole image stack: STK@  
C          ` ---> LUNS***  
C          | 
C          |  Templated stacked MRC image: file*@MRC or *@file.MRCS
C          ` ---> GETNEWIMG_MRC    --> OPFILEC  
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--********************************************************************* 
        SUBROUTINE GETNEWIMG(LUNCP,LUN,LUNXM,FILPAT,NWANTT, 
     &                       SAYIT,NGOT,IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMLIMIT.INC'

        INTEGER                :: LUNCP,LUN,LUNXM
        CHARACTER(LEN=*)       :: FILPAT
        INTEGER                :: NWANTT
        LOGICAL                :: SAYIT
        INTEGER                :: NGOT
        INTEGER                :: IRTFLG

        CHARACTER(LEN=MAXNAM)  :: FILNAM
        CHARACTER(LEN=1)       :: DSP
        LOGICAL                :: FOUROK,ISBARE,IS_MRC
        INTEGER                :: NWANT,NLET,LOCAST,LOCAT,NX1,NY1,NZ1
        INTEGER                :: ITYPE1,NX,NY,NZ,ITYPE,ISTK,IRTFLGT
        INTEGER                :: MAXIM

        LOGICAL                :: ISMRCFILE   ! FUNCTION

        INTEGER                :: lnblnkn     ! FUNCTION


        NWANT  = NWANTT          ! NWANTT MAY NOT BE WRITABLE
        NLET   = lnblnkn(FILPAT)
        LOCAST = INDEX(FILPAT(1:NLET),'*')
        LOCAT  = INDEX(FILPAT(1:NLET),'@')
        IS_MRC = ISMRCFILE(FILPAT)

#if defined(SP_DBUGIO)
        write(3,*)' In getnewimg - filpat:       ', trim(filpat)
        write(3,*)' In getnewimg - nwant :       ', nwant 
        write(3,*)' In getnewimg - locast,locat: ', locast,locat 
        write(3,*)' ----------- 55555555555555555555555 ---------------'
#endif

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
     &                 MAXIM,'~9',FOUROK,IRTFLG) 
           IF (IRTFLG .NE. 0) RETURN 

#if defined(SP_DBUGIO)
           write(3,*)' Opened new Xmipp selfile file: ',trim(filnam)
#endif

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
     &                 MAXIM,' ',FOUROK,IRTFLG) 
           IF (IRTFLG .NE. 0) RETURN 

#if defined(SP_DBUGIO)
           write(3,*)' Opened new templated file: ',trim(filnam)
#endif

           NGOT   = NWANT   
           IRTFLG = 0
           RETURN

        ELSEIF (IS_MRC .AND. LOCAT > 0 .AND. LOCAST < LOCAT) THEN

C          TEMPLATED STACKED MRC IMAGE ---------------------- *@STK.mrc
C          BARE MRC IMAGE ------------------------------------ @STK.mrc
           
           CALL GETNEWIMG_MRC(LUN,FILPAT,NWANT,SAYIT,
     &                        FILNAM,NGOT,IRTFLG)
           RETURN


        ELSEIF (LOCAT > 0) THEN
C          STACKED IMAGE ------------------------------- STK@ or @STK@

           CALL LUNGETISBARE(LUN,ISBARE,IRTFLG)

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

C          NEED ISTK 
           CALL LUNCOPYSTK(LUN,ISTK,IRTFLGT)

#if defined(SP_DBUGIO)
           write(3,*)' In getnewimg - isbare,filpat: ',isbare, filpat
           write(3,*)' In getnewimg - nwant,maxim,istk: ',  
     &                                nwant,maxim,istk
#endif

          IF (ISTK < 0) THEN
C             MAKING A NEW INDEXED STACKED FILE, UPDATE INDX LOCATION
              CALL LUNWRTINDX(LUN,NWANT,NX,IRTFLGT)
              IF (IRTFLGT .NE. 0) RETURN
           ENDIF

           IF (NWANT > MAXIM .OR. ISTK > 2) THEN
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

#if defined(SP_DBUGIO)
           write(3,*)' Opened new stacked file: ',filnam(1:nlet)
#endif

        ENDIF

        END









