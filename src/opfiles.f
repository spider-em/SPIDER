 
C++*********************************************************************
C                                                                      
C  OPFILES.F   NEW                              12/15/06  ArDean Leith    
C              BAD NUMBRT() TRAP                05/21/09  ArDean Leith  
C              ASKNAM, PROMPTEX                 12/06/10  ArDean Leith 
C              NX...                            03/26/12  ArDean Leith 
C              COPY NON SPIDER INPUT            05/26/14  ArDean Leith 
C              COPY NON SPIDER INPUT            05/26/14  ArDean Leith 
C              LOCAST, ASKLIST FOR ILIST        10/02/14  ArDean Leith
C              DO NOT CHECK .MRC FOR XMIPP       7/25/19  ArDean Leith
C              MRC SUPPORT                       8/05/19  ArDean Leith
C              COMMENTS                          2/06/20  ArDean Leith
C
C ********************************************************************** 
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2020  Health Research Inc.,                         *
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
C  CONTAINS:  OPFILES, GETOLDIMG, GETNEWIMG, NEXTFILE, NEXTFILES 
C
C  OPFILES(LUNCP,LUNIMG,LUNDOC,LUNXM,  ASKNAM,FILPAT,NLET, DISP,
C          ITYPE,NX,NY,NZ,MAXIM, PROMPT,
C          FOUROK, ILIST,NIMAXT, UNUSED
C          NTOT,IMGNUM, IRTFLG)
C 
C  PURPOSE: SOLICITS FILE NAME(S) AND OPENS FILE(S)
C           SUPPORT ROUTINE FOR CONVERTING OPERATIONS TO 
C           WORK ON WHOLE STACK OR WITH A SELECTION DOC FILE.
C           ONLY USED IN A FEW OPERATIONS.  RAN OUT OF YEARS TO
C           COMPLETE THE IMPLEMENTATION IN MORE OPERATIONS
C           ALSO WORKS FOR MRC FILES NOW
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
C        ITYPE      IFORM FOR FILE                        (SENT OR RET) 
C        NX,NY,NZ   IMAGE SIZE                            (SENT OR RET)
C
C        MAXIM      STACK INDICATOR                          (SENT/RET)
C                   ON INPUT: STACK INDICATOR IF DISP == 'I'     (SENT)
C                   ON OUTPUT:                                    (RET)
C                       -2 A NON-STACK IMAGE                
C                       -1 A SINGLE STACKED IMAGE                  
C                     >= 0 CURRENT MAX. IMAGE NO. IN STACK             
C        PROMPT     PROMPT FOR FILNAME                           (SENT)
C                     IF NOT (ASKNAM) THIS IS FILE NAME          (SENT)
C                     ~ (TILDE) IN LAST CHAR. SAYS SKIP
C                       "FILE" AT END OF PROMPT
C                     ~9 IN NEXT TO LAST OR
C                        NEXT-TO-NEXT-TO LAST
C                        ACCEPTS AN EXTENSION
C                        (OTHERWISE DISCARDED!)
C                     ~6 KEEPS OLD DATE/TIME
C        FOUROK     CAN USE EXISTING FOURIER FILES?              (SENT)
C        ILIST      IMAGE NUMBER LIST                             (RET)
C                     IF NIMAXT < 0 MUST BE SENT
C                     NOT USED IF SINGLE IMAGE/SELFILE 
C        NIMAXT     MAX LENGTH OF IMAGE NUMBER LIST              (SENT)
C                     <0 MEANS DO NOT ASK FOR LIST
C        NTOT:      # OF IMAGES IN IMAGE NUMBER LIST              (RET)
C                   ZERO FOR SINGLE IMAGE AND NO * 
C        IMGNUM     IMAGE NUMBER THAT IS CURRENTLY OPEN      (SENT/RET)
C                   ON INPUT:   IF (BARESTACK) IS # WANTED
C                   ON OUTPUT:  <0 IS SELFILE IN USE 
C        ILIST      IMAGE NUMBER LIST                             (RET)
C                     IF NIMAXT < 0 MUST BE SENT
C                     UNUSED IF SINGLE IMAGE/SELFILE 
C        UNUSED     UNUSED PARAMETER                                (?)
C        NTOT       # OF IMAGES IN IMAGE NUMBER LIST              (RET)
C                     ZERO FOR SINGLE IMAGE AND NO * 
C        IRTFLG     ERROR FLAG (0 IS NORMAL)                      (RET)
C                      -1 GOTO PREVIOUS QUESTION
C                      -5 NOT A SPIDER OR MRC FILE
C 
C  CALL TREE:
C
C     OPFILES 
C        v
C     FILERD
C        v 
C     FILELIST
C        |  Templated stack: STK@*** 
C        ` ---> FILNAMANDEXT --> OPFILEC --> GETOLDIMG
C        |                               --> GETNEWIMG
C        |
C        |  Whole barestack: STK@     
C        ` ---> FILNAMANDEXT --> OPFILEC --> GETOLDIMG
C        |                               --> GETNEWIMG
C        | 
C        |  File template:   IMG***  
C        ` ---> FILGET       --> OPFILEC  
C        | 
C        |  Simple file:    IMG001
C        ` --------------------> OPFILEC  
C        | 
C        |  MRC file:       *.MRC or *.MRCS
C        ` -> OPFILES_MRC
C           |
C           |  Templated stack:  **@STK.MRC
C           ` ---> FILGET_AT --> OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC 
C           |                                            --> OPENFIL_N_MRC
C           |
C           |  Whole barestack:  @STK.MRC 
C           ` ---> FILGET_AT --> OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC
C           |                                            --> OPENFIL_N_MRC
C           | 
C           |  File template:    IMG***.MRC
C           ` ---> FILGET_AT --> OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC
C           |                                            --> OPENFIL_N_MRC
C           | 
C           |  Simple file:      IMG001.MRC
C           ` -----------------> OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC  
C                                                        --> OPENFIL_N_MRC
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12 
C--********************************************************************* 

       SUBROUTINE OPFILES(LUNCP,LUNIMG,LUNDOC,LUNXM,
     &                    ASKNAM,FILPAT,NLET, DISP,
     &                    ITYPE,NX,NY,NZ,MAXIM,
     &                    PROMPT,
     &                    FOUROK,ILIST,NIMAXT, 
     &                    UNUSED,NTOT,IMGNUM, IRTFLG) 
 
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
       INTEGER                   :: NIMAXT,UNUSED,NTOT,IMGNUM, IRTFLG

       INTEGER                   :: LUNOP,NLETT,LOCTILDE,LOCAT,LOCAST
       INTEGER                   :: IMGWANT,IMGNUMIN,LENE,I
       INTEGER                   :: NIMAXP,LUNXM,NDUM
       LOGICAL                   :: ISOPEN,GOTFILE
       LOGICAL                   :: SAYIT,ASKLIST

       CHARACTER (LEN=MAXNAM)    :: FILNAMT,FILNAM 
       CHARACTER (LEN=2*MAXNAM)  :: MESG 
       CHARACTER (LEN=100)       :: PROMPTEX
       CHARACTER (LEN=1)         :: CDUM 
       CHARACTER (LEN=1)         :: DISPT 
       CHARACTER (LEN=1)         :: NULL = CHAR(0)

       LOGICAL                   :: IS_MRC,IS_BARE
       LOGICAL                   :: ISMRCFILE    ! FUNCTION 

       INTEGER                   :: LNBLNKN      ! FUNCTION

       LOCTILDE = INDEX(PROMPT,'~') ! unfinished

       !write(3,*)' In opfiles, imgnum: ',imgnum

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
 
       LOCAT    = INDEX(FILPAT(1:NLET),'@')   
       LOCAST   = INDEX(FILPAT(1:NLET),'*')
       
       !write(3,*)' locast,locat:',locast,locat,filpat(1:nlet)
       !write(6,*)' ilist:',ilist(1)
       !write(6,*)' locast,locat,nstack:',locast,locat,maxim
       !write(3,*)' in opfiles, nimaxt:',nimaxt

       ASKLIST = (NIMAXT > 0)   ! ASK FOR IMAGE NUMBER(S)
       NIMAXP  = ABS(NIMAXT)    ! ILIST DIMENSION, IF ILIST IN USE
 
       !write(3,*)' in opfiles, nimaxt: ',nimaxt
       !write(3,*)' in opfiles, asklist: ',asklist
       !write(3,*)' In opfiles, imgnum: ',imgnum

       IMGWANT = MAX(1,IMGNUM)   ! DEFAULT IMAGE NUMBER     

       IF (LOCAST > 0 .AND. ASKLIST) THEN
C         GET LIST OF IMAGES FROM DOC. FILE OR INPUT LINE
          NTOT   = 0 
          CALL FILELIST(.FALSE.,LUNDOC,CDUM,NDUM,ILIST,NIMAXP,
     &                  NTOT,' ',IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

C         START WITH FIRST FILE IN SERIES
          IMGWANT = ILIST(1)
          !write(3,*)' In opfiles, imgwant aa: ',imgwant

       ELSEIF (LOCAST > 0 .AND. NIMAXT < 0) THEN

C         ILIST CONTAINS TEMPLATED SERIES, 
C         START WITH FIRST FILE IN THIS SERIES
          IMGWANT = ILIST(1)
          NTOT    = NIMAXP
          !write(3,*)' In opfiles, imgwant bb: ',imgwant

       ELSEIF (.NOT. ASKLIST) THEN
          NTOT = NIMAXP
       ENDIF

       IF (IMGWANT < 0 .OR. IMGWANT > 10000000) THEN
          CALL ERRT(102,'INVALID IMAGE NUMBER',IMGWANT)
          IRTFLG = 1
          RETURN
       ENDIF

       !write(3,*)' In opfiles, nimaxt: ',nimaxt
       !write(3,*)' In opfiles, ilist: ',ilist(1:5)
       !write(3,*)' In opfiles, asklist,nimaxp: ',asklist,nimaxp
       !write(3,*)' In opfiles, maxim,ntot,imgwant: ',maxim,ntot,imgwant

C      SEE IF THIS IS A MRC FILE
       IS_MRC = ISMRCFILE(FILPAT)

       !write(3,*)' In opfiles, imgwant cc: ',imgwant
       IF (IS_MRC) THEN
C         OPEN MRC FILE
          CALL OPFILES_MRC(LUNCP,LUNIMG,LUNDOC,
     &                     FILPAT,NLET, DISP,
     &                     ITYPE,NX,NY,NZ,MAXIM,
     &                     IMGWANT, IRTFLG) 
          
          IMGNUM = IMGWANT

          RETURN
       ENDIF

       IMGNUMIN = IMGNUM
       IMGNUM   = 0 
       !write(3,*) ' In opfiles - imgnum 0: ',imgnum

       IF (LOCAT > 0 .AND. LOCAST > LOCAT) THEN
C         TEMPLATED STACKED FILE: STK@**** -------------- _9@* or STK@**

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

          !write(6,*) ' In opfiles - imgwant,imgnum1: ',imgwant,imgnum

C         OPEN THE STACK FILE HEADER 
          MAXIM  = 1  
          CALL OPFILEC(LUNCP,.FALSE.,FILNAM,LUNIMG,DISP,
     &                         ITYPE,NX,NY,NZ, 
     &                         MAXIM,' ',FOUROK,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9000 

          !write(6,*) ' In opfiles - locast: ',locat,locast,filnam
          !write(6,*) ' In opfiles - imgwant,imgnum2: ',imgwant,imgnum


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

             !write(6,*) ' In opfiles - imgwant,imgnum: ',imgwant,imgnum

             CALL GETNEWIMG(LUNCP,LUNIMG,LUNDOC,FILPAT,IMGWANT,
     &                      SAYIT,IMGNUM,IRTFLG)

             !write(6,*) ' In opfiles - imgwant,imgnum3: ',imgwant,imgnum

          ELSE
C            WANT TO USE EXISTING SPIDER STACK FILE
             CALL GETOLDIMG(LUNIMG,LUNDOC,FILPAT, IMGWANT,SAYIT,
     &                         FOUROK,IMGNUM,IRTFLG)
          ENDIF
          IF (IRTFLG .NE. 0) GOTO 9000

C         RETRIEVE CURRENT MAXIMUM IMAGE NUMBER FROM OVERALL HEADER
          CALL LUNGETMAXIM(LUNIMG,MAXIM,IRTFLG)

          !write(6,'(a,a,a,i6,a,i6,a,i6)')' opened templated stack: ',
!     &          filpat(1:nlet),' at image: ',imgnum, 
!     &          '  maxim: ',maxim

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

          !write(3,*)' Opened bare stack:',filpat(1:nlet),' img:',imgnum
          !write(6,*)' Opened imgnum, stack:',imgnum,filpat(1:nlet)
          !write(3,*)' Opened imgwant,ntot:',imgwant,imgnum,ntot

       ELSEIF (LOCAST > 0) THEN
C         A SIMPLE FILE TEMPLATE: IMG*** ----------------------- IMG***

C         FIND IMGNUM FOR FIRST FILE IN THE SERIES

          IMGNUM = ILIST(1)
          IF (IMGNUM < 0 .OR. IMGNUM > 10000000) THEN
              CALL ERRT(102,'INVALID IMAGE NUMBER',IMGNUM)
              IRTFLG = 1
              GOTO 9000
          ENDIF

C         SUBSTITUTE IMGNUM INTO FILPAT 
          CALL  FILGET(FILPAT,FILNAM,NLET,IMGNUM,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9000

C         OPEN FIRST FILE IN THE SERIES
          MAXIM = 0 
          CALL OPFILEC(LUNCP,.FALSE.,FILNAM,LUNIMG,DISP,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,' ',FOUROK,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9000

          !write(6,*)' Opened templated file: ',filpat(1:nlet),
          !&              '  for:',ntot,' images.'

       ELSE
C         SINGLE SIMPLE INPUT FILE: IMG001 -------------------- IMG001
C         OR XMIPP SELFILE LISTING FILE: SELX ----------------- SELFILE
C         OR TRYING TO COPY NON-SPIDER FILE   ----------------- NONSPIFILE

C         CHECK FOR XMIPP SELFILE LIST
          IF (LUNXM > 0 .AND. .NOT. ISMRCFILE(FILPAT) ) THEN
             !write(6,*)' Filpat for openxmsel: ',filpat(:nlet) 
             CALL OPENXMSELFILE(FILPAT,LUNXM,FILNAM,NLET,NTOT,IRTFLG)
             !write(6,*)' Filnam from openxmsel: ',ntot,':',filnam(:nlet) 

             INQUIRE(FILE=FILNAM(:NLET),EXIST=GOTFILE)
           
             IF (NTOT > 0 .AND. GOTFILE) THEN
C               OPEN FIRST FILE IN XMIPP SELFILE LIST
                MAXIM = 0  
                CALL OPFILEC(LUNCP,.FALSE.,FILNAM(:NLET),LUNIMG,DISP,
     &                       ITYPE,NX,NY,NZ, 
     &                       MAXIM,'dum~9',FOUROK,IRTFLG)
                IF (IRTFLG .NE. 0) GOTO 9000

                IMGNUM     = -1
                !write(6,*)' Opened selfile image: ',filnam(1:nlet) 
                RETURN

             ENDIF
          ENDIF

C         SINGLE SIMPLE INPUT FILE: IMG001 ------------------ IMG001

          !write(3,*)' In opfiles, simple file: ',filpat(1:nlet),itype,nx 
          MAXIM = 0  
          CALL OPFILEC(LUNCP,.FALSE.,FILPAT,LUNIMG,DISP,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,PROMPT,FOUROK,IRTFLG)

C         RETURN FILENAME WITH ANY EXTENSION IF NOT SPIDER IMAGE
          IF (IRTFLG == 5) NLET = lnblnkn(FILPAT)
          IF (IRTFLG .NE. 0) GOTO 9000

          NTOT       = 0
          IMGNUM     = 1

          !write(3,*)' In opfiles, Opened simple: ',filpat(1:nlet)
       ENDIF

9000   RETURN

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
C  CALL TREE:
C
C       GETOLDIMG 
C          |
C          |  Templated simple image:  IMG*** 
C          ` ---> FILGET --> OPFILEC  
C          |                         
C          |  Templated stacked image: STK@*(**     
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
     &                 MAXIM,'~9',FOUROK,IRTFLG)
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
     &                 MAXIM,' ',FOUROK,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          IMAGE SIZE SHOULD BE SAME
           CALL SIZCHK(NULL,NX1,NY1,NZ1,ITYPE1,
     &                      NX ,NY, NZ, ITYPE, IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           !write(6,*)' Opened old templated file: ',filnam(1:nlet),maxim

           NGOT   = NWANT   
           IRTFLG = 0
           RETURN

        ELSEIF (IS_MRC .AND. LOCAT > 0 .AND. LOCAST < LOCAT) THEN

C          TEMPLATED STACKED MRC IMAGE ---------------------- **@STK.mrc
C          BARE STACKED MRC IMAGE      ----------------------   @STK.mrc
           
          !write(3,*) ' In getoldimg, nwant: ',nwant,locat,locast,filpat
           CALL GETOLDIMG_MRC(LUN,FILPAT,NWANT,SAYIT,
     &                        FILNAM,NGOT,IRTFLG)
           RETURN

        ELSEIF (LOCAT  > 0 .AND. LOCAST > LOCAT) THEN
C          TEMPLATED STACKED SPIDER IMAGE --------------------- STK@***

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
C  GETNEWIMG(LUNCP,LUN,LUNXM,FILPAT,NWANTT, SAYIT,NGOT,IRTFLG)
C
C  PURPOSE:       OPEN A SPECIFIED IMAGE WITHIN STACK FOR RANDOM 
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
C          |  Templated stacked MRC image: *@MRC or *@.MRCS
C          ` ---> GETNEWIMG_MRC    --> OPFILEC  
C
C **********************************************************************

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
        INTEGER                :: ITYPE1,NX,NY,NZ,ITYPE,ISTACK,IRTFLGT
        INTEGER                :: MAXIM

        LOGICAL                :: ISMRCFILE   ! FUNCTION

        INTEGER                :: lnblnkn     ! FUNCTION


        NWANT  = NWANTT          ! NWANTT MAY NOT BE WRITABLE
        NLET   = lnblnkn(FILPAT)
        LOCAST = INDEX(FILPAT(1:NLET),'*')
        LOCAT  = INDEX(FILPAT(1:NLET),'@')
        IS_MRC = ISMRCFILE(FILPAT)

        !write(3,*)' In getnewimg - nwant,locat:',nwant,locat,filpat

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
     &                 MAXIM,' ',FOUROK,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           !write(6,*)' Opened new templated file: ',filnam(1:nlet)

           NGOT   = NWANT   
           IRTFLG = 0
           RETURN

        ELSEIF (IS_MRC .AND. LOCAT > 0 .AND. LOCAST < LOCAT) THEN

C          TEMPLATED STACKED MRC IMAGE ---------------------- **@STK.mrc
C          BARE MRC IMAGE ------------------------------------ @STK.mrc

           CALL GETNEWIMG_MRC(LUN,FILPAT,NWANT,SAYIT,
     &                        FILNAM,NGOT,IRTFLG)
           RETURN


        ELSEIF (LOCAT > 0) THEN
C          STACKED IMAGE ------------------------------- STK@*  or STK@

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

C          NEED ISTACK 
           CALL LUNCOPYSTK(LUN,ISTACK,IRTFLGT)

           !write(6,*)' In getnewimg - isbare,filpat: ',isbare, filpat
           !write(6,*)' In getnewimg -nwant,maxim,istak: ',nwant,maxim,istack


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
C  NEXTFILES.F  NEW                              12/15/06 ArDean Leith
C               OVERUN OUTPUT LIST = -99          1/15/12 ArDean Leith
C **********************************************************************
C
C NEXTFILES(NINDX1, NINDX2, INUMBR1,INUMBR2, 
C           FOUROK,NLIST1,NLIST2,   
C           NSTACK1,NSTACK2,   
C           LUN1,LUNCP,LUN2, FILPAT1,FILPAT2,
C           IMGNUM1,IMGNUM2,IRTFLG) 
C
C PURPOSE:  GETS NEXT INPUT AND OUTPUT FILES FOR A STACK ORIENTED 
C           OPERATION.  STACKS MUST BE OPENED WITH OPFILES!!!
C
C PARAMETERS: NINDX1,NINDX2    LIST INDICES                 (SENT/RET.)
C             INUMBR1,INUMR2   IMAGE NUMBER LISTS                (SENT)
C             FOUROK           FOURIER INPUT IS OK               (SENT)
C             LUNXM1,LUNXM2    LUN FOR SELFILE INPUT             (SENT)
C             NLIST1,NLIST2    NUMBER OF IMAGES                  (SENT)
C             NSTACK1,NSTACK2  MAX IMAGE IN STACK                (SENT)
C             LUN1             LUN FOR INPUT  (0 = NO FILE IN)   (SENT)
C             LUNCP            LUN FOR OUTPUT HEADER COPY        (SENT)
C             LUN2             LUN FOR OUTPUT (0 = NO FILE OUT)  (SENT)
C             FILPAT,FILPAT2   FILE NAME PATTERNS                (SENT)
C             IMGNUM1,IMGNUM2  IMAGE NUMBERS                (SENT/RET.)
C             IRTFLG           ERROR (0 IS OK, -1 IS END STACK)  (RET.)
C
C  CALL TREE:
C     AFTER OPFILES HAVE BEEN CALLED:
C
C       NEXTFILES 
C          |      
C          |      Old image
C          ` ---> GETOLDIMG    
C          |      
C          |      New image                         
C          ` ---> GETNEWIMG
C
C--*********************************************************************
 
      SUBROUTINE NEXTFILES(NINDX1, NINDX2, INUMBR1,INUMBR2, 
     &                     FOUROK,LUNXM1,LUNXM2,
     &                     NLIST1,NLIST2,   NSTACK1,NSTACK2,   
     &                     LUN1,LUNCP,LUN2, FILPAT1,FILPAT2,
     &                     IMGNUM1,IMGNUM2, IRTFLG) 
 
      IMPLICIT NONE

      INTEGER           :: NINDX1,NINDX2
      INTEGER           :: INUMBR1(NLIST1),INUMBR2(NLIST2)
      LOGICAL           :: FOUROK
      INTEGER           :: LUNXM1,LUNXM2
      INTEGER           :: NLIST1,NLIST2
      INTEGER           :: NSTACK1,NSTACK2,LUN1,LUNCP,LUN2
      CHARACTER(LEN=*)  :: FILPAT1,FILPAT2
      INTEGER           :: IMGNUM1,IMGNUM2,IRTFLG

      INTEGER           :: NWANT1,NWANT2, it
      LOGICAL           :: SAYIT = .TRUE.
      LOGICAL           :: GOTAST1, GOTAST2
      LOGICAL           :: IS_BARE1,IS_BARE2    

      NINDX1 = NINDX2 + 1
      NINDX2 = NINDX2 + 1

      IF (LUN1 > 0) THEN
C        OPEN NEXT INPUT FILE 
         GOTAST1 = (INDEX(FILPAT1,'*') > 0)

C        IS THIS A BARE STACK OPERATION?  (OK FOR SPIDER & MRC)
         CALL LUNGETISBARE(LUN1,IS_BARE1,IRTFLG)

         IF (IMGNUM1 == -1 .AND. LUNXM1 > 0 ) THEN
C           XMIPP SELFILE LISTED IMAGE
            NWANT1 = -1

         ELSEIF ( IS_BARE1 ) THEN
C           MRC OR SPIDER BARE STACK INPUT   (NO LIST)

            IF (IMGNUM1 < 0) IMGNUM1 = 1
            NWANT1 = IMGNUM1 + 1
            IF (NWANT1 > NSTACK1) THEN
C              FINISHED THE WHOLE STACK
               IRTFLG = -1
               RETURN
            ENDIF

         ELSEIF (NSTACK1 == -2  .OR.
     &           NSTACK1 == -1  .OR.
     &           NSTACK1  >    0 ) THEN
C           NON STACKED IMAGE WITH/WITHOUT TEMPLATED LIST
C           STACKED     IMAGE WITH/WITHOUT LIST         

            IF (NINDX1 > NLIST1) THEN
C              OVERUN INPUT LIST
               IRTFLG = -1
               RETURN
            ENDIF

C           OPEN NEXT INPUT FILE 
            NWANT1 = INUMBR1(NINDX1)
         ENDIF

         !write(3,*)' In nextfiles, nwant1: ',nwant1
         !write(6,'(a,8i5)')' In nextfiles, nwant1: ',
         !                                  nwant1,lun1,nwant1,imgnum1
         CALL GETOLDIMG(LUN1,LUNXM1,FILPAT1,NWANT1,SAYIT, 
     &                  FOUROK,IMGNUM1,IRTFLG)

         !write(3,*)' gotoldimg, nstacki1,gotast1: ',nstack1,gotast1

         IF (IRTFLG < 0)    RETURN    ! END OF WHOLE-STACK
         IF (IRTFLG .NE. 0) RETURN    ! ERROR

         IF (IS_BARE1) THEN
C           INPUT FROM A BARE STACK 
            NINDX1 = IMGNUM1
            NINDX2 = IMGNUM1
         ENDIF
      ENDIF
      ! write(6,*) 'In nextfiles1: l1,l2,irtflg:',lun1,lun2,irtflg

      IF (LUN2 > 0) THEN
C        OPEN NEXT OUTPUT FILE 
         GOTAST2 = (INDEX(FILPAT2,'*') > 0)

C        IS THIS IS A BARE STACK OPERATION?  (OK FOR SPIDER & MRC)
         CALL LUNGETISBARE(LUN2,IS_BARE2,IRTFLG)


!        !write(6,'(a,8i5)')'  out: nindx2,ngot2,nstack2,imgnum2: ',
!     &                             nindx2,ngot2,nstack2,imgnum2
!         write(6,*) 'nextfiles2: l2,imgnum2,gotast2:',
!     &                         lun2,imgnum2,gotast2
         IF (IMGNUM2 == -1 .AND. LUNXM2 > 0  ) THEN
C           XMIPP SELFILE LISTED IMAGE
            NWANT2 = -1

         ELSEIF (LUN1 > 0 .AND. IS_BARE1 ) THEN
C           NON-STACK IMAGE WITH/WITHOUT TEMPLATE LIST 
C           MRC OR SPIDER BARE STACK INPUT   (NO LIST AVAILABLE)
C           BARE STACK INPUT , USE SAME  OUTPUT IMAGE NUMBER
            NWANT2 = IMGNUM1

            !write(3,*)'In nextfiles, :',isbare,nstack2,nwant2

         ELSEIF (IS_BARE2 ) THEN
C           MRC OR SPIDER BARE STACK OUTPUT   (NO LIST)
            IF (IMGNUM2 < 0) IMGNUM2 = 1
            NWANT2 = IMGNUM2 + 1

            IF (LUN1 > 0 .AND. IS_BARE1) THEN
C              BARE STACK INPUT , USE SAME  OUTPUT IMAGE NUMBER
                NWANT2 = IMGNUM1
            ENDIF

            !write(3,*)'In nextfiles, :',isbare,nstack2,nwant2


         ELSEIF (NSTACK2 == -2  .OR.
     &           NSTACK2 == -1  .OR.
     &           NSTACK2  >  0 ) THEN

C           NON-STACK IMAGE WITH/WITHOUT TEMPLATE LIST 
            IF (NINDX2 > NLIST2) THEN
C               OVERUN OUTPUT LIST
                !write(6,*) 'nindx2 > nlist2',nindx2,nlist2
                IRTFLG = -99
                RETURN
            ENDIF
C           OPEN NEXT OUTPUT FILE 
            NWANT2 = INUMBR2(NINDX2)
         ENDIF

         !write(3,*)' In nextfiles, nwant2,nindx2:',nwant2,nindx2
         !write(6,'(a,8i5)') 
         !& ' Calling getnew,nwant2,imgnum2,nstack2:',
         !                   nwant2,imgnum2,nstack2

         CALL GETNEWIMG(LUNCP,LUN2,LUNXM2,FILPAT2,NWANT2,
     &                  SAYIT,IMGNUM2,IRTFLG)

         !write(3,'(a,8i5)') ' Getnew,nwant2,imgnum2,nstack2,irtflg:',
      !&                              nwant2,imgnum2,nstack2,irtflg

         IF (IRTFLG .NE. 0) RETURN     ! ERROR

      ENDIF
      !write(3,'(a,5i5)')' In nextfiles, nlist1,nstack1,imgnum1:',
      !&                                 nlist1,nstack1,imgnum1
      !write(3,'(a,5i5)')' In nextfiles, ngot2,nstack2,imgnum2:',
      !&                                 ngot2,nstack2,imgnum2

      END


C++*********************************************************************
C
C NEXTFILE.F    NEW                              12/15/06 ArDean Leith
C               OVERUN OUTPUT LIST = -99          1/15/12 ArDean Leith
C               MRC SUPPORT                       8/22/19 ArDean Leith
C
C **********************************************************************
C
C NEXTFILE(NINDX1, INUMBR1, FOUROK, LUNXM1, NLIST1,   NSTACK1,   
C          LUN1,   LUNCP,   FILPAT1,  DISP, IMGNUM1, IRTFLG)  
C
C PURPOSE:  GETS NEXT INPUT OR OUTPUT FILE FOR BARE STACK INPUT,
C           A NON STACKED IMAGE WITH/WITHOUT TEMPLATED LIST, 
C           OR A  STACKED IMAGE WITH/WITHOUT LIST
C           STACKS MUST BE OPENED WITH OPFILES!!!
C
C PARAMETERS: NINDX1         LIST INDEX                    (SENT/RET.)
C             INUMBR1        IMAGE NUMBER LIST                  (SENT)
C             FOUROK         FOURIER INPUT IS OK                (SENT)
C             LUNXM1         LUN FOR SELFILE INPUT              (SENT)
C             NLIST1         NUMBER OF IMAGES IN LIST           (SENT)
C             NSTACK1        HIGHEST IMAGE IN STACK             (SENT)
C             LUN1           LUN FOR I/0                        (SENT)
C             LUNCP          LUN FOR OUTPUT HEADER COPY         (SENT)
C             FILPAT1        FILE NAME PATTERN                  (SENT)
C             DISP           IMAGE EXISTENCE                    (SENT)
C             IMGNUM1        IMAGE NUMBER                  (SENT/RET.)
C             IRTFLG         ERROR (0 IS OK, -1 IS END STACK)   (RET.)
C
C  CALL TREE:
C     AFTER OPFILES HAVE BEEN CALLED:
C
C       NEXTFILE 
C          |      
C          |      Old image
C          ` ---> GETOLDIMG    
C          |        |
C          |        |  Templated simple image:  IMG*** 
C          |        ` ---> FILGET --> OPFILEC  
C          |        |                         
C          |        |  Templated stacked image: STK@*(**     
C          |        ` ---> LUN***
C          |        | 
C          |        |  Whole image stack:       STK@  
C          |        ` ---> LUNS***
C          |        | 
C          |        |  Templated stacked MRC:  *@MRC or *@.MRCS
C          |        ` ---> GETOLDIMG_MRC --> OPFILEC  
C          |      
C          |      
C          |      New image                         
C          ` ---> GETNEWIMG
C                   |
C                   |  Templated simple image:  IMG*** 
C                   ` ---> FILGET --> OPFILEC  
C                   |                         
C                   |  Templated stacked image: STK@*(**     
C                   ` ---> LUN***
C                   | 
C                   |  Whole image stack:       STK@  
C                   ` ---> LUNS***  
C                   | 
C                   |  Templated stacked MRC:  *@MRC or *@.MRCS
C                   ` ---> GETNEWIMG_MRC --> OPFILEC  
C
C--*********************************************************************
 
      SUBROUTINE NEXTFILE(NINDX1,   INUMBR1, 
     &                    FOUROK,   LUNXM1,
     &                    NLIST1,   NSTACK1,   
     &                    LUN1,     LUNCP, 
     &                    FILPAT1,  DISP,
     &                    IMGNUM1,  IRTFLG) 
 
      IMPLICIT NONE

      INTEGER           :: NINDX1 
      INTEGER           :: INUMBR1(NLIST1) 
      LOGICAL           :: FOUROK
      INTEGER           :: LUNXM1 
      INTEGER           :: NLIST1 
      INTEGER           :: NSTACK1, LUN1,LUNCP 
      CHARACTER(LEN=*)  :: FILPAT1
      CHARACTER(LEN=1)  :: DISP
      INTEGER           :: IMGNUM1
      INTEGER           :: IRTFLG

      INTEGER           :: NWANT1, it
      LOGICAL           :: SAYIT = .TRUE.
      LOGICAL           :: GOTAST,IS_BARE1

      LOGICAL           :: ISMRCFILE    ! FUNCTION
      
      NINDX1 = NINDX1 + 1

      GOTAST = (INDEX(FILPAT1,'*') > 0)

C     IS THIS IS A BARE STACK OPERATION?
      CALL LUNGETISBARE(LUN1,IS_BARE1,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C      write(3,*)' In nextfile, gotast, filpat1:  ', gotast,filpat1
C      write(6,*)' In nextfile, nindx1,isbare:    ',nindx1,is_bare1
C      write(3,*)' In nextfile, ----------------------: '
C      write(3,'(A,5i6)')
C     &         '  In nextfile, num1,nindx1,nlist1,nwant1,nstack1: ',
C     &                      imgnum1,nindx1,nlist1,nwant1,nstack1

      IF (IMGNUM1 == -1  .AND. LUNXM1 > 0 ) THEN
C        XMIPP SELFILE LISTED IMAGE
         NWANT1 = -1

      ELSEIF ( IS_BARE1 ) THEN
C        MRC OR SPIDER BARE STACK INPUT   (NO LIST)

         IF (IMGNUM1 < 0) IMGNUM1 = 1
         NWANT1 = IMGNUM1 + 1

C        write(3,'(A,5i6)')'  In nextfile, nwant1,nstack1: ',
C       &                                  nwant1,nstack1

        IF (NWANT1 > NSTACK1) THEN
C           FINISHED THE WHOLE STACK
            IRTFLG = -1
            RETURN
         ENDIF
          
      ELSEIF (NSTACK1 == -2  .OR.
     &        NSTACK1 == -1  .OR.
     &        NSTACK1  >    0 ) THEN
C        NON STACKED IMAGE WITH/WITHOUT TEMPLATED LIST
C        STACKED     IMAGE WITH/WITHOUT LIST
         IF (NINDX1 > NLIST1) THEN
C           OVERUN I/O LIST
            IRTFLG = -1
            RETURN
         ENDIF

C        OPEN NEXT I/O FILE 
         NWANT1 = INUMBR1(NINDX1)

      ENDIF

      !write(3,*)' In nextfile, nindx1,nlist1,nstack1,imgnum1: ',
      !     &                   nindx1,nlist1,nstack1,imgnum1
      !write(3,*)'  In nextfile, calling get, nwant1,imgnum1,nstack1:',
      !&                                      nwant1,imgnum1,nstack1

      IF (DISP == 'O' .OR. DISP == 'B' .OR. 
     &    DISP == 'Z' .OR. DISP == 'E') THEN 
  
C        OPEN NEXT INPUT FILE
         !write(3,*) ' Call gotoldimg, nwant1: ',nwant1,imgnum1
         CALL GETOLDIMG(LUN1,LUNXM1,FILPAT1,NWANT1,SAYIT, 
     &                  FOUROK,IMGNUM1,IRTFLG)

         !write(3,*)' Gotoldimg, nlist1,nwant1,imgnum1: ',
         !&                      nlist1,nwant1,imgnum1,irtflg,nstack1

         IF (IRTFLG    < 0) RETURN    ! END OF WHOLE-STACK
         IF (IRTFLG .NE. 0) RETURN    ! ERROR

      ELSE   

C        OPEN NEXT OUTPUT FILE 
         !write(3,*) ' Call gotnewimg, nwant1: ',nwant1,imgnum1
         CALL GETNEWIMG(LUNCP,LUN1,LUNXM1,FILPAT1,NWANT1,
     &                  SAYIT,IMGNUM1,IRTFLG)

         !write(3,*) ' Getnewimg, nwant1,imgnum1,nstack1,irtflg:',
         !     &                  nwant1,imgnum1,nstack1,irtflg
         IF (IRTFLG .NE. 0) RETURN   ! ERROR

      ENDIF

      IF (IS_BARE1) THEN
C        BARE STACK, OUTPUT IMAGE HAS SAME # AS INPUT ALWAYS
         NINDX1 = IMGNUM1
      ENDIF

      !   write(3,'(a,5i5)')' In nextfile, nlist1,nstack1,imgnum1:',
      !&                                   nlist1,nstack1,imgnum1
      END


