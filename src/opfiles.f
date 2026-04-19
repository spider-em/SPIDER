
C++*********************************************************************
C                                                                      
C  FILE_OPEN_NOW.F NEW                         04/6/26    ArDean Leith
C 
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2026  Health Research Inc.,                         *
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
C  FILE_OPEN_NOW(FILNAM,LUN,IS_OPEN,IRTFLG)
C
C  PARAMETERS:
C       FILNAM    FILE NAME                                     (SENT)
C       LUN       UNIT TO OPEN FILE ON                          (SENT)
C       IRTFLG    CHARS IN FILE NAME PATTERN                     (RET)
C
C 
C      PURPOSE:  SEE IF FILNAM  ALREADY OPENED ON ANOTHER LUN
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************                                                                  
       SUBROUTINE FILE_OPEN_NOW(FILNAM,LUN,IS_OPEN,IRTFLG)                                  
       IMPLICIT NONE
       INCLUDE 'CMBLOCK.INC' 
       INCLUDE 'CMLIMIT.INC' 
       
       CHARACTER(LEN=*)          :: FILNAM 
       INTEGER                   :: LUN
       LOGICAL                   :: IS_OPEN
       INTEGER                   :: IRTFLG
	
       CHARACTER(LEN=MAXNAM+10)  :: FILNAMPE 
       CHARACTER(LEN=2*MAXNAM)   :: MESG 

       INTEGER                   :: LUNOP,NLET

       INTEGER                   :: lnblnkn


       IRTFLG = 0

       IF (FILNAM(1:1) == '_') RETURN    ! INLINE FILE

       FILNAMPE = TRIM(FILNAM)

       IF (INDEX(FILNAMPE,'.')  <= 0) THEN
          FILNAMPE = TRIM(FILNAM) // '.' // TRIM(DATEXC)
       ENDIF
       NLET = lnblnkn(FILNAMPE)

       INQUIRE(FILE=FILNAMPE,OPENED=IS_OPEN,NUMBER=LUNOP)

       IF (IS_OPEN .AND. LUNOP .NE. LUN) THEN
           MESG = '  FILE: ' // FILNAMPE(1:NLET) //
     &             '  ALREADY OPENED ON: ' 
           NLET = lnblnkn(MESG)

           WRITE(NOUT,'(A,I3)') MESG(1:nlet),LUNOP 

           !CALL ERRT(102,MESG(1:nlet),LUNOP) 
           IRTFLG = -2
        ENDIF

        RETURN
        END


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
C            DEBUG OUTPUT ADDED                4/14/26  ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2026  Health Research Inc.,                         *
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
C
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
C                     >= 0 CURRENT MAX IMAGE # IN STACK 
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
C  CALLED BY: 
C      util_1010,  util_1011, util_1110,   util_11,   util2sup,
C      adds        apscc      bp32f,       bp3f,      bpcg,      bprp3,
C      copy1,      copyd,     copyfrommrc, copytomrc, copytovol, falb,
C      gali,       hali,      interps,     qstat,     qstatloc,  tfed,
C      rotqss,     ssnrb,     softmask   
C
C  CALL TREE:
C
C     OPFILES Uses (E.G. STK@*)~~9  prompt for filename
C        v    '~9' in last 3 char of prompt accepts ext. on filename 
C        |
C        ` -> FILERD  
C        |       ` ->  ECHONAME
C 
C        ` -> FILELIST
C        |
C        |    MRC file:        *.mrc or *.mrcs
C        ` -> OPFILES_MRC
C        | 
C        |    Templated stack: stk@* 
C        ` -> FILNAMANDEXT --> OPFILEC --> GETOLDIMG
C        |                             --> GETNEWIMG
C        |
C        |    Whole barestack: stk@     
C        ` -> FILNAMANDEXT --> OPFILEC --> GETOLDIMG
C        |                             --> GETNEWIMG
C        |
C        |    Templated file:  img***  
C        ` -> FILGET      -->  OPFILEC  
C        | 
C        |    Simple file:     img001
C        ` -> OPFILEC  
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

       INTEGER                   :: ISTK,NSTK_FLG

       INTEGER                   :: LOCAT,LOCAST,  ngot
       INTEGER                   :: LUNOP,NLETT,LOCTILDE, nleti,nletn
       INTEGER                   :: IMGWANT,LENE,I,IFLIPIN
       INTEGER                   :: NIMAXP,LUNXM,NDUM,NLETX
       LOGICAL                   :: IS_OPEN,GOT_FILE,GOT_FILE_X
       LOGICAL                   :: SAYIT,ASKLIST

       LOGICAL                   :: WANT_NEW_IMG_OR_FILE
       LOGICAL                   :: WANT_OLD_IMG_OR_FILE
 
       CHARACTER (LEN=MAXNAM)    :: FILNAMPE,FILNAM 
       CHARACTER (LEN=MAXNAM+1)  :: FILPATP1 
       CHARACTER (LEN=2*MAXNAM)  :: MESG 
       CHARACTER (LEN=100)       :: PROMPTEX
       CHARACTER (LEN=1)         :: CDUM 
       CHARACTER (LEN=1)         :: DISPT 
       CHARACTER (LEN=1)         :: NULL = CHAR(0)

       CHARACTER(LEN=MAXNAM)     :: FIL_NOAT,FIL_DIRS,FIL_BASE,FIL_EXT 
       LOGICAL                   :: IS_MRC, IS_MRCS, IS_BARE
       INTEGER                   :: LIST_VAL

       LOGICAL                   :: ismrcfile    ! FUNCTION 
       INTEGER                   :: lnblnkn      ! FUNCTION
  
C      LOCTILDE = INDEX(PROMPT,'~') ! FUTURE??
C      LENP = lnblnkn(PROMPT)
C      PUT ~7 AT END OF PROMPT TO OPEN* TO OPEN A STACK WITHOUT @
C      OPSTKNOAT =  (INDEX (PROMPT(1:LENP),'~7') > 1)
C      PUT ~9 AT END OF PROMPT TO OPEN* TO KEEP INCOMING EXTENSION
C      KEEPEXT = (INDEX(PROMPT(1:LENP),'~9') > 1)


       DISPT = DISP

       WANT_NEW_IMG_OR_FILE = 
     &     (DISPT == 'N' .OR. DISPT == 'I' .OR. DISPT == 'U')

       WANT_OLD_IMG_OR_FILE = 
     &     (DISPT == 'O' .OR. DISPT == 'I' .OR. DISPT == 'E')

       IF (ASKNAM .AND. PROMPT == NULL) THEN
C         ASK FOR FILE NAME, CAN ACCEPT EXTENSION

          IF (WANT_NEW_IMG_OR_FILE) THEN

C            NEW FILE/IMG, USE DEFAULT PROMPT FOR OUTPUT
             PROMPTEX = 
     &          'OUTPUT FILE NAME OR TEMPLATE (E.G. stk-file@*)~'

          ELSE
C            OLD FILE/IMG, USE DEFAULT PROMPT FOR INPUT
             PROMPTEX = 
     &          'INPUT FILE NAME OR TEMPLATE (E.G. stk-file@*)~'
          ENDIF
     
C         USING PROMPTX, READ IN FILPAT
          CALL FILERD(FILPAT,NLET,NULL,PROMPTEX,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN 

C         NOTE: trim(filpat) can leave a null at end!!!
          NLET  = lnblnkn(FILPAT)

cc#define SP_DBUGIO
#if defined (SP_DBUGIO)
          write(3,*)' '
          write(3,*)' In opfiles 00A; promptex: ', trim(promptex)
          write(3,*)' In opfiles 00A; filpat:   ', filpat(:nlet)
          write(3,*)  '  '
#endif

       ELSEIF (ASKNAM) THEN
C         USING PROMPT, READ IN FILPAT      
          CALL FILERD(FILPAT,NLET,NULL,PROMPT,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN 
          NLET = lnblnkn(FILPAT)

#if defined (SP_DBUGIO)
          write(3,*)'  '
          write(3,*)' In opfiles 00B; prompt: ', trim(prompt)
          write(3,*)' In opfiles 00B; filpat: ', filpat(:nlet)
#endif

       ELSE
C         USE INCOMING 'PROMPT' AS FILENAME 

          NLET   = lnblnkn(PROMPT)
          FILPAT = PROMPT(1:NLET)

#if defined (SP_DBUGIO)
          write(3,*)'  '
          write(3,*)' In opfiles 00C; prompt: ', prompt(1:nlet)
          write(3,*)' In opfiles 00C; filpat: ', filpat
#endif
       ENDIF




 
C      HAVE GOT FILENAME OR TEMPLATE ---------------------------------

C      GET: LOCAT,LOCAST, ASKLIST, NIMAXP, LISTVAL, LIST_SIZE, ILIST

       LOCAT     = INDEX(FILPAT(1:NLET),'@')   
       LOCAST    = INDEX(FILPAT(1:NLET),'*')
     
       ASKLIST   = (NIMAXT > 0)   ! ASK FOR IMAGE NUMBER(S)
       NIMAXP    = ABS(NIMAXT)    ! ILIST DIMENSION, IF ILIST IN USE

       LIST_VAL  = 1   ! DEFAULT IMAGE NUMBER 
       LIST_SIZE = 0   ! DEFAULT NUMBER UNUSED        

#if defined (SP_DBUGIO)
       write(3,*)'  '
       write(3,*)' In opfiles 004, locat:   ', locat
       write(3,*)' In opfiles 004, locast:  ', locast
       write(3,*)' In opfiles 004, asklist:           ', asklist
       write(3,*)' In opfiles 004, nimaxp:  ', nimaxp
#endif


       IF (LOCAST > 0 .AND. ASKLIST) THEN

C         GET LIST OF IMAGE #, FROM DOC. FILE OR INPUT LINE ----------
          CALL FILELIST(.FALSE.,LUNDOC,CDUM,NDUM,ILIST,NIMAXP,
     &                  LIST_SIZE,' ',IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

C         START WITH FIRST FILE/IMG IN SERIES
          LIST_VAL = ILIST(1)

#if defined (SP_DBUGIO)
          write(3,*)'  '
          write(3,*)' In opfiles 005, list_size: ',list_size
          write(3,*)' In opfiles 005, list_val:  ',list_val
#endif

       ELSEIF (LOCAST > 0 .AND. NIMAXT < 0) THEN

C         ILIST ALREADY CONTAINS TEMPLATED SERIES, 
C         START WITH FIRST FILE IN THIS SERIES ----------------------

          LIST_VAL  = ILIST(1)
          LIST_SIZE = NIMAXP

#if defined (SP_DBUGIO)
          write(3,*)' In opfiles 006, list_size: ',list_size
          write(3,*)' In opfiles 006, list_val:  ',list_val
#endif


       ELSEIF (LOCAT > 0 .AND. LOCAST == 0 ) THEN

C         BARE STACK, START WITH FIRST IMAGE? -----------------------
          LIST_VAL  = 1
          LIST_SIZE = 0

#if defined (SP_DBUGIO)
          write(3,*)' '
          write(3,*)' In opfiles 007, list_size: ',list_size
          write(3,*)' In opfiles 007, list_val:  ',list_val
#endif


       ELSEIF (.NOT. ASKLIST) THEN
C         NO IMAGE LIST WANTED  -------------------------------------

          LIST_VAL  = 1
          LIST_SIZE = 0

#if defined (SP_DBUGIO)
         write(3,*)' In opfiles 008, list_size: ',list_size
         write(3,*)' In opfiles 008, list_val:  ',list_val
#endif

       ELSEIF (LOCAST  == 0 .AND. ASKLIST .AND. 
     &        (LIST_VAL < 0 .OR.  LIST_VAL > 10000000)) THEN

C         SIMPLE FILE NAME ??? --------------------------------------
          LIST_VAL = 0

#if defined (SP_DBUGIO)
          write(3,*)' In opfiles 009, list_val:  ',list_val
          write(3,*)' In opfiles 009, list_size: ',list_size
          write(3,*)' In opfiles 009, imgwant:   ',imgwant  
#endif
       ENDIF

       IF (LIST_VAL < 0 .OR. LIST_VAL > 10000000) THEN
          CALL ERRT(102,'INVALID IMAGE NUMBER',IMGWANT)
          IRTFLG = 1
          GOTO 9000
       ENDIF

       IMGWANT = LIST_VAL


#if defined (SP_DBUGIO)
       write(3,*)'  '
       write(3,*)' In opfiles 111; dispt:                 ', dispt
       write(3,*)' In opfiles 111; locat:     ', locat
       write(3,*)' In opfiles 111; locast:    ', locast
       write(3,*)' In opfiles 111; nimaxt:    ', nimaxt
       write(3,*)' In opfiles 111; nimaxp:    ', nimaxp
       write(3,*)' in opfiles 111, list_size: ', list_size
       write(3,*)' in opfiles 111, list_val:  ', list_val  
       write(3,*)' in opfiles 111, imgwant:   ', imgwant 
       !write(3,*)'  ' 
#endif

C      GET: FIL_NOAT, FIL_DIRS, FIL_BASE, FIL_EXT
C           IS_MRC, IS_MRCS, IS_BARE, ISTK  
      
       CALL GET_FILNAM_INFO(FILPAT, LIST_VAL, 
     &                     FIL_NOAT,FIL_DIRS,FIL_BASE,FIL_EXT,
     &                     IS_MRC, IS_MRCS, IS_BARE, GOT_FILE,
     &                     ISTK,   IRTFLG)

       FILNAM = FIL_NOAT  ! FULL PATH TO FILE, MINUS SPIDER EXTENSION
       NLETN  = lnblnkn(FILNAM)

#if defined (SP_DBUGIO)
       write(3,*)'   '
       write(3,*)' In opfiles 666; is_mrc:              ', is_mrc
       write(3,*)' In opfiles 666; is_mrcs:             ', is_mrcs
       write(3,*)' In opfiles 666; is_bare:             ', is_bare
       write(3,*)' In opfiles 666; got_file:            ', got_file
       write(3,*)' In opfiles 666; filpat:     ', filpat(1:nlet)
       write(3,*)' In opfiles 666; fil_noat:   ', trim(fil_noat)
       write(3,*)' In opfiles 666; filnam:     ', filnam(:nletn)
      !write(3,*)' In opfiles 666; fil_ext:    ', trim(fil_ext)

       write(3,*)' In opfiles 666; istk:       ', istk
       write(3,*)' In opfiles 666; list_val:   ', list_val
       write(3,*)' In opfiles 666; list_size:  ', list_size
       write(3,*)' In opfiles 666; ilist(1:2): ', ilist(1:2)
#endif


       IF (IS_MRC) THEN
C         OPEN MRC FILE -------------------- MRC ------------------

#if defined (SP_DBUGIO)
          write(3,*)' B4 opfiles MRC; imgwant: ', imgwant
          write(3,*)' B4 opfiles MRC; disp:    ', disp
          write(3,*)' Calling opfiles_mrc -----'
#endif

          CALL OPFILES_MRC(LUNCP,LUNIMG,LUNDOC,
     &                     FILPAT,NLET, DISP,
     &                     ITYPE,NX,NY,NZ, MAXIM, 
     &                     ILIST, LIST_SIZE, LIST_VAL,  
     &                     IRTFLG)
 
C         MAXIM              ! RETURNED TO CALLER HERE
          IMGNUM = LIST_VAL  ! RETURNED TO CALLER HERE
C         LIST_SIZE          ! RETURNED TO CALLER HERE
 
#if defined (SP_DBUGIO)
          write(3,*)' After Opfiles MRC: nz:              ', nz 
          write(3,*)' After Opfiles MRC: maxim:           ', maxim 
          write(3,*)' After Opfiles MRC: imgnum=list_val: ', imgnum 
          write(3,*)' After Opfiles MRC: list_size:       ', list_size 
          write(3,*)' After Opfiles MRC: irtflg:          ', irtflg 
          write(3,*)'   '                                          
#endif

          RETURN
       ENDIF  ! -------------- End of mrc ---------------------



       IMGNUM   = 0 

ccccccccccccccccccccccccc OLD TEMPLATED STACK    cccccccc OTS cccccccc


       IF (LOCAT  > 0      .AND. 
     &     LOCAST > LOCAT  .AND.
     &     WANT_OLD_IMG_OR_FILE) THEN

C         OLD TEMPLATED STACKED FILE: stk@* ------------- _9@* or stk@*
         
C         ERROR IF THIS FILE IS ALREADY OPENED ON ANOTHER LUN
          CALL FILE_OPEN_NOW(FIL_NOAT,LUNIMG,IS_OPEN,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9000

#if defined (SP_DBUGIO)
          write(3,*)'  '
          nleti = lnblnkn(filnam)
          write(3,*)' In opfiles OTS0, filnam:     ', filnam(:nleti)  
          write(3,*)' In opfiles OTS0, dispt:            ', dispt
          write(3,*)' In opfiles OTS0, Calling opfilec'
          write(3,*)'  '
#endif

C         OPEN THE STACK FILE HEADER 
          MAXIM  = 1  
          CALL OPFILEC(LUNCP,.FALSE.,FILPAT,LUNIMG,DISPT,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,' ',FOUROK,IRTFLG) 
          IF (IRTFLG .NE. 0) GOTO 9000 

#if defined (SP_DBUGIO)
          write(3,*)'  '
          write(3,*)' In opfiles OTS1, filnam: ', filpat(1:nlet)  
          write(3,*)' In opfiles OTS1, maxim:  ', maxim
          write(3,*)'  '
#endif

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

C         WANT TO USE EXISTING SPIDER STACK FILE
          CALL GETOLDIMG(LUNIMG,LUNDOC,FILPAT, IMGWANT,SAYIT,
     &                         FOUROK,IMGNUM,IRTFLG)
#if defined (SP_DBUGIO)
          write(3,*)' '
          write(3,*)' In opfiles OTS2, filpat:  ', filpat(1:nlet)  
          write(3,*)' In opfiles OTS2, imgwant: ', imgwant
          write(3,*)' In opfiles OTS2, imgnum:  ', imgnum
#endif

          IF (IRTFLG .NE. 0) GOTO 9000

C         RETRIEVE CURRENT MAXIMUM IMAGE NUMBER FROM OVERALL HEADER
          CALL LUNGETMAXIM(LUNIMG,MAXIM,IRTFLG)

#if defined (SP_DBUGIO)
          write(3,*) ' '
          write(3,*)' In opfiles OTS3, filpat:  ', filpat(1:nlet)  
          write(3,*)' In opfiles OTS3, imgwant: ', imgwant
          write(3,*)' In opfiles OTS3, imgnum:  ', imgnum
          write(3,*)' In opfiles OTS3, maxim:   ', maxim
#endif


cccccccccccccccccccccc NEW TEMPLATED STACK  ccccccccccc NTS cccccccccc

       ELSEIF (LOCAT > 0      .AND.
     &        LOCAST > LOCAT  .AND.
     &        WANT_NEW_IMG_OR_FILE) THEN
 
C         NEW TEMPLATED STACKED IMG: stk@* ----------- _9@* or stk@*
         
C         ERROR IF THIS FILE IS ALREADY OPENED ON ANOTHER LUN
          CALL FILE_OPEN_NOW(FIL_NOAT,LUNIMG,IS_OPEN,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9000

#if defined (SP_DBUGIO)
          write(3,*)'  '
          write(3,*)' In opfiles NTS1, filnam:     ', filnam(:nlet)  
          write(3,*)' In opfiles NTS1, disp:             ', disp
          write(3,*)' In opfiles NTS1, Calling opfilec'
          write(3,*)'  '
#endif

C         OPEN THE STACK FILE HEADER 
          MAXIM  = 1  
          CALL OPFILEC(LUNCP,.FALSE.,FILPAT,LUNIMG,DISPT,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,' ',FOUROK,IRTFLG) 
          IF (IRTFLG .NE. 0) GOTO 9000 

#if defined (SP_DBUGIO)
          write(3,*)'  '
          write(3,*)' In opfiles NTS2, filnam: ', filnam(:nletn)  
          write(3,*)' In opfiles NTS2, maxim:  ', maxim
#endif

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

#if defined (SP_DBUGIO)
          write(3,*)
          write(3,*)' In opfiles NTS3, filpat:  ', filpat(:nlet)  
          write(3,*)' In opfiles NTS3, imgwant: ', imgwant
          write(3,*)' In opfiles NTS3, imgnum:  ', imgnum
#endif

          CALL GETNEWIMG(LUNCP,LUNIMG,LUNDOC,FILPAT,IMGWANT,
     &                      SAYIT,IMGNUM,IRTFLG)

#if defined (SP_DBUGIO)
          write(3,*)
          write(3,*)' In opfiles NTS4, filpat:  ', filpat(:nlet)  
          write(3,*)' In opfiles NTS4, imgwant: ', imgwant
          write(3,*)' In opfiles NTS4, imgnum:  ', imgnum
#endif

          IF (IRTFLG .NE. 0) GOTO 9000

C         RETRIEVE CURRENT MAXIMUM IMAGE NUMBER FROM OVERALL HEADER
          CALL LUNGETMAXIM(LUNIMG,MAXIM,IRTFLG)

#if defined (SP_DBUGIO)
          write(3,*)
          write(3,*)' In opfiles NTS5, filpat:  ', filpat(:nlet)  
          write(3,*)' In opfiles NTS5, imgwant: ', imgwant
          write(3,*)' In opfiles NTS5, imgnum:  ', imgnum
#endif

ccccccccccccccccccccccccc  OLD BARE  ccccccccccccccc obs cccccccccccc

       ELSEIF (IS_BARE   .AND. 
     &         WANT_OLD_IMG_OR_FILE) THEN

C         WANT OLD SPIDER BARE STACK INPUT:  ---------- stk@,@stk,_9@

          FILNAM = TRIM(FIL_NOAT)
          NLETN  = lnblnkn(FILNAM)

C         ERROR IF THIS FILE ALREADY OPENED ON ANOTHER LUN
          CALL FILE_OPEN_NOW(FILNAM,LUNIMG,IS_OPEN,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9000

#if defined (SP_DBUGIO)
          write(3,*)' '
          write(3,*)' In opfiles OBS1 ; filpat:   ', filpat(1:nlet)
          write(3,*)' In opfiles OBS1 ; filnam:   ', filnam(1:nletn)
          write(3,*)' In opfiles OBS1 ; dispt:             ', dispt
          write(3,*)' In opfiles OBS1 ; Calling opfilec'
#endif

          CALL OPFILEC(LUNCP,.FALSE.,FILPAT,LUNIMG,DISPT,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,'INPUT',FOUROK,IRTFLG) 
          IF (IRTFLG .NE. 0) GOTO 9000

          IF (MAXIM > 0) THEN
C            CREATE IMAGE NUMBER LIST IN: ILIST
             LIST_SIZE = 0
             DO I= 1,MAXIM
                LIST_SIZE = LIST_SIZE + 1
                IF (LIST_SIZE > NIMAXP) THEN
                   CALL ERRT(102,'IMAGE # LIST OVERFLOW AT IMAGE:',
     &                            LIST_SIZE)
                   GOTO 9000
                ENDIF
                ILIST(LIST_SIZE) = I
             ENDDO
          ENDIF

C         EXISTING BARE STACK, OPEN FIRST FILE IN STACK -----------
          IMGWANT = 1
          
#if defined (SP_DBUGIO)
          write(3,*)' '
          write(3,*)' In opfiles OBS2 ; dispt:               ', dispt
          write(3,*)' In opfiles OBS2 ; maxim:       ', maxim 
          write(3,*)' In opfiles OBS2 ; imgwant:     ', imgwant
          write(3,*)' In opfiles OBS2 ; list_size:   ', list_size
          write(3,*)' In opfiles OBS2 ; Calling getoldimg'
#endif

          SAYIT   = .TRUE.
          CALL GETOLDIMG(LUNIMG,LUNDOC,FILPAT, IMGWANT,
     &                      SAYIT,FOUROK,IMGNUM,IRTFLG)

#if defined (SP_DBUGIO)
          write(3,*)
          write(3,*)' In opfiles OBS3; maxim:   ', maxim
          write(3,*)' In opfiles OBS3; imgwant: ', imgwant
          write(3,*)' In opfiles OBS3; imgnum:  ', imgnum
          write(3,*)' In opfiles OBS3; filpat:  ',filpat(1:nlet)  
#endif
          IF (IRTFLG .NE. 0) GOTO 9000
 
cccccccccccccccc   NEW BARE STACK FILE cccccccccccc nbs ccccccccccccc

       ELSEIF (IS_BARE         .AND. 
     &        (.NOT. GOT_FILE) .AND. 
     &        WANT_NEW_IMG_OR_FILE) THEN

C         NEW BARE STACK, OPEN REQUESTED FILE IN STACK -----------

          IF (DISPT == 'I') THEN
C             OPEN NEW BARE INDEXED STACK 
              DISPT = 'N'
              MAXIM = -MAXIM  ! FLAG FOR INDEXED STACK
          ENDIF

#if defined (SP_DBUGIO)
          write(3,*)' '
          write(3,*)' In opfiles NBS1 ; filpat:  ',filpat(1:nlet)  
          write(3,*)' In opfiles NBS1 ; got_file:          ', got_file
          write(3,*)' In opfiles NBS1 ; disp:              ', disp
          write(3,*)' In opfiles NBS1 ; dispt:             ', dispt
          write(3,*)' In opfiles NBS1 ; Calling opfilec'
#endif

          CALL OPFILEC(LUNCP,.FALSE.,FILPAT,LUNIMG,DISPT,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,'INPUT',FOUROK,IRTFLG) 
          IF (IRTFLG .NE. 0) GOTO 9000

#if defined (SP_DBUGIO)
          write(3,*)' '
          write(3,*)' In opfiles NBS2 ; maxim:       ', maxim 
          write(3,*)' In opfiles NBS2 ; imgwant:     ', imgwant
          write(3,*)' In opfiles NBS2 ; Calling getnewimg'
#endif

          CALL GETNEWIMG(LUNCP,LUNIMG,LUNDOC,FILPAT,IMGWANT,
     &                   SAYIT,IMGNUM,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9000

#if defined (SP_DBUGIO)
          write(3,*)'   '
          write(3,*)' In opfiles NBS3; imgwant: ', imgwant
          write(3,*)' In opfiles NBS3; imgnum:  ', imgnum
#endif
 

cccccccccccccc  APPEND NEW IMAGE TO OLD BARE STACK cccc nt0 ccccccc

       ELSEIF (IS_BARE   .AND. 
     &         GOT_FILE  .AND. 
     &         WANT_NEW_IMG_OR_FILE) THEN

          NSTK_FLG = 3

          IF (DISPT == 'I') THEN
C             OPEN  BARE INDEXED STACK 
              DISPT    = 'N'
              NSTK_FLG = -NSTK_FLG  ! FLAG FOR INDEXED STACK
          ENDIF

C         KLUDGE TO OPEN THE OVERALL STACK FILE HEADER 
          FILPATP1 = FILPAT(1:NLET) // '1'

#if defined (SP_DBUGIO)
          write(3,*)' '
          write(3,*)' In opfiles NTO1, filpat:   ', filpat(1:nlet)
          write(3,*)' In opfiles NTO1, filpatp1: ', filpatp1(1:nlet+1)
          write(3,*)' In opfiles NTO1, disp:               ', disp
          write(3,*)' In opfiles NTO1, dispt:              ', dispt
  
          write(3,*)' In opfiles NTO1; nstk_flg: ', nstk_flg
          write(3,*)' In opfiles NTO1; Calling opfilec'
#endif


          CALL OPFILEC(LUNCP,.FALSE.,FILPATP1,LUNIMG,DISPT,
     &                 ITYPE,NX,NY,NZ, 
     &                 NSTK_FLG,' ',FOUROK,IRTFLG) 
          IF (IRTFLG .NE. 0) GOTO 9000 

C         RETRIEVE CURRENT MAXIMUM IMAGE NUMBER FROM OVERALL HEADER
          CALL LUNGETMAXIM(LUNIMG,MAXIM,IRTFLG)

          IMGWANT = MAXIM + 1   ! NEW IMG #

C         THIS IS NOT A BARE STACK FILE
          CALL LUNSETISBARE(LUNIMG,.FALSE.,IRTFLG)

#if defined (SP_DBUGIO)
          write(3,*)
          write(3,*)' In opfiles NTO3, filpat:   ', filpat(:nlet)  
          write(3,*)' In opfiles NTO3, nstk_flg: ', nstk_flg
          write(3,*)' In opfiles NTO3, maxim:    ', maxim
          write(3,*)' In opfiles NTO3, imgwant:  ', imgwant
#endif

          SAYIT = .TRUE.
          CALL GETNEWIMG(LUNCP,LUNIMG,LUNDOC,FILPAT,IMGWANT,
     &                      SAYIT,IMGNUM,IRTFLG)

#if defined (SP_DBUGIO)
          write(3,*)
          write(3,*)' In opfiles NTO4, filpat:  ', filpat(:nlet)  
          write(3,*)' In opfiles NTO4, imgwant: ', imgwant
          write(3,*)' In opfiles NTO4, imgnum:  ', imgnum
#endif

          IF (IRTFLG .NE. 0) GOTO 9000




c-------------------------- unused -----------------------
#if defined (NEVER)

#if defined (SP_DBUGIO)
          write(3,*)' '
          write(3,*)' In opfiles NTO1, filpat:   ', filpat(1:nlet)
          lene = lnblnkn(filnam)
          write(3,*)' In opfiles NTO1, filnam:   ', filnam(1:nletn)
          write(3,*)' In opfiles NTO1, disp:               ', disp
          write(3,*)' In opfiles NTO1, dispt:              ', dispt
  
          write(3,*)' In opfiles NTO1; nstk_flg: ', nstk_flg
          write(3,*)' In opfiles NTO1; Calling opfilec'
#endif

C         OPEN THE STACK FILE HEADER 
          CALL OPFILEC(LUNCP,.FALSE.,FILPAT,LUNIMG,DISPT,
     &                 ITYPE,NX,NY,NZ, 
     &                 NSTK_FLG,'INPUT',FOUROK,IRTFLG) 
          IF (IRTFLG .NE. 0) GOTO 9000


C         RETRIEVE CURRENT MAXIMUM IMAGE NUMBER FROM OVERALL HEADER
          CALL LUNGETMAXIM(LUNIMG,MAXIM,IRTFLG)

          IMGWANT = MAXIM + 1   ! NEW STACK LENGTH

#if defined (SP_DBUGIO)
          write(3,*)'   '
          write(3,*)' In opfiles NTO2; nstk_flg:    ', nstk_flg
          write(3,*)' In opfiles NTO2; maxim:       ', maxim 
          write(3,*)' In opfiles NTO2; imgwant:     ', imgwant
          write(3,*)' In opfiles NTO2; Calling getoldimg'
#endif

          SAYIT   = .TRUE.
          CALL GETOLDIMG(LUNIMG,LUNDOC,FILPAT, IMGWANT,
     &                      SAYIT,FOUROK,IMGNUM,IRTFLG)

#if defined (SP_DBUGIO)
          write(3,*)
          write(3,*)' In opfiles NTO4; maxim:   ', maxim
          write(3,*)' In opfiles NTO4; imgwant: ', imgwant
          write(3,*)' In opfiles NTO4; imgnum:  ', imgnum
          write(3,*)' In opfiles NTO4; filpat:  ',filpat(1:nlet)  
#endif
          IF (IRTFLG .NE. 0) GOTO 9000
 
          CALL GETOLDIMG(LUNIMG,LUNDOC,FILPAT, IMGWANT,
     &                   SAYIT,FOUROK,IMGNUM,IRTFLG)

#if defined (SP_DBUGIO)
          write(3,*)
          write(3,*)' In opfiles NTO3; imgwant: ', imgwant
          write(3,*)' In opfiles NTO3; imgnum:  ', imgnum
          write(3,*)' In opfiles NTO3; irtflg:  ', irtflg
#endif
          IF (IRTFLG .NE. 0) GOTO 9000

          IMGWANT = MAXIM + 1  ! APPEND TO EXISTING STACK

#if defined (SP_DBUGIO)
          write(3,*)
          write(3,*)' In opfiles NTO4; imgwant: ', imgwant
          write(3,*)' In opfiles NTO4; imgnum:  ', imgnum
          write(3,*)' In opfiles NTO4; filpat:  ', filpat(1:nlet)  
          write(3,*)'   '
#endif
          CALL ADD_TO_BARE(LUNIMG,FILPAT,IMGWANT, 
     &                     .TRUE.,NGOT,IRTFLG)

#endif
c ------------------------- end of unused --------------------


#if defined (SP_DBUGIO)
          write(3,*)
          write(3,*)' In opfiles NTS5; irtflg:    ', irtflg
          write(3,*)' In opfiles NTS5; list_size: ',list_size
          write(3,*)' In opfiles NTS5; imgwant:   ', imgwant
          write(3,*)' In opfiles NTS5; imgnum:    ', imgnum
          write(3,*)' In opfiles NTS5; filpat:    ', filpat(1:nlet)  
#endif

          IF (IRTFLG .NE. 0) GOTO 9000

       ELSEIF (LOCAST > 0) THEN

cccccccccccccccccc  SIMPLE TEMPLATED  FILE ccccccc  img***  ccccccccc

C         SIMPLE TEMPLATED FILE: IMG*** ---------------------- img***
C         NOT A TEMPLATED STACK 

C         SUBSTITUTE IMGWANT INTO FILPAT 
          CALL FILGET(FILPAT,FILNAM,NLETN,IMGWANT,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9000 

C         OPEN FIRST FILE IN THE SERIES
          MAXIM = 0 
          CALL OPFILEC(LUNCP,.FALSE.,FILNAM,LUNIMG,DISPT,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,' ',FOUROK,IRTFLG) 
          IF (IRTFLG .NE. 0) GOTO 9000 

#if defined(SP_DBUGIO)
          write(3,*)' '
          write(3,*)' In opfiles ST0; imgwant:   ', imgwant
          write(3,*)' In opfiles ST0; list_size: ', list_size
          write(3,*)' In opfiles ST0; filpat:    ', filpat(:nlet)
          write(3,*)' In opfiles ST0; filnam:    ', filnam(:nletn)
          write(3,*)'  '
#endif


       ELSEIF (LOCAT  > 0      .AND. 
     &         WANT_OLD_IMG_OR_FILE) THEN

ccccccccccccc OLD SINGLE STACKED FILE: imgstk@1 ---------- imgstk@7

#if defined(SP_DBUGIO)
          write(3,*)' '
          write(3,*)' In opfiles OSS1; stacked img: ', filpat(:nlet)
          write(3,*)' In opfiles OSS1; filpat:      ', filpat(:nlet)
          write(3,*)' In opfiles OSS1; dispt:          |', dispt,'|'
          write(3,*)' In opfiles OSS1; maxim:       ', maxim
          write(3,*)' In opfiles OSS1; Calling opfilec '
          write(3,*)'  '
#endif
          CALL OPFILEC(LUNCP,.FALSE.,FILPAT,LUNIMG,DISPT,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,PROMPT,FOUROK,IRTFLG) 

#if defined(SP_DBUGIO)
          write(3,*)'   '
          write(3,*)' In opfiles OSS2; Opened:    ', filpat(:nlet)
          write(3,*)' In opfiles OSS2; maxim:     ', maxim
          write(3,*)' In opfiles OSS2; itype:     ', itype
          write(3,*)' In opfiles OSS2; irtflg:    ', irtflg
#endif
 
C         RETURN FILENAME WITH ANY EXTENSION IF NOT SPIDER IMAGE
          IF (IRTFLG == 5) NLET = lnblnkn(FILPAT)
          IF (IRTFLG .NE. 0) GOTO 9000 

          LIST_SIZE  = 0
          IMGNUM     = 1

#if defined(SP_DBUGIO_NEVER)
          write(3,*)'  '
          write(3,*)' In opfiles OSS3: Opened stacked: ', filpat(:nlet)

#endif

       ELSEIF (LOCAT  > 0      .AND. 
     &        WANT_NEW_IMG_OR_FILE) THEN

ccccccccccc NEW SINGLE STACKED FILE: ccccc imgstk@1 ---------- imgstk@1

#if defined(SP_DBUGIO)
          write(3,*)' '
          write(3,*) ' In opfiles NSS1; stacked img: ',filpat(:nlet)
           write(3,*)' In opfiles NSS1; filpat: ',  filpat(:nlet)
           write(3,*)' In opfiles NSS1; dispt:           |', dispt,'|'
           write(3,*)' In opfiles NSS1; maxim:       ', maxim
           write(3,*)' In opfiles NSS1; Calling opfilec '
           write(3,*)'  '
#endif
          CALL OPFILEC(LUNCP,.FALSE.,FILPAT,LUNIMG,DISPT,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,PROMPT,FOUROK,IRTFLG) 

          IF (IMGWANT > MAXIM) MAXIM = IMGWANT


#if defined(SP_DBUGIO)
          write(3,*)'   '
          write(3,*)' In opfiles NSS2; Opened:    ', filpat(:nlet)
          write(3,*)' In opfiles NSS2; maxim:     ', maxim
          write(3,*)' In opfiles NSS2; itype:     ', itype
          write(3,*)' In opfiles NSS2; irtflg:    ', irtflg

          IF (IMGWANT > MAXIM) MAXIM = IMGWANT
          write(3,*)' In opfiles NSS2; new maxim: ', maxim

#endif
 
          IF (IRTFLG .NE. 0) GOTO 9000 

          LIST_SIZE  = 0
          IMGNUM     = 1

#if defined(SP_DBUGIO_NEVER)
          write(3,*)'  '
          write(3,*)' In opfiles NSS3: Opened simple: ', filpat(:nlet)
#endif

       ELSE

ccccccccccccccccccc ONE SIMPLE FILE  -------- sf --------- img001 
ccccccccccccccccccc or XMIPP SELFILE LISTING FILE: ------- selfile
ccccccccccccccccccc or TRY TO COPY NON-SPIDER FILE ------- anyfile.ext

#if defined(SP_DBUGIO)
          write(3,*)' '
          write(3,*) ' In opfiles SF1; Simple file: ',filpat(:nlet)
#endif
               
C         CHECK FOR XMIPP SELFILE LIST
          IF (LUNXM > 0 .AND. .NOT. ISMRCFILE(FILPAT) ) THEN
             !write(3,*)' Filpat for openxmsel: ',filpat(:nlet) 
             CALL OPENXMSELFILE(FILPAT,LUNXM,FILNAM,NLETX,
     &                          LIST_SIZE,IRTFLG)
#if defined(SP_DBUGIO)
             nletn = lnblnkn(filnam)
             write(3,*)' '
             write(3,*)' In opfiles SF2; selected file: ',filnam(:nletx)
             write(3,*)' In opfiles SF2; list_size: ',list_size
#endif
             INQUIRE(FILE=FILNAM(:NLETX),EXIST=GOT_FILE_X)
           
             IF (LIST_SIZE > 0 .AND. GOT_FILE_X) THEN
C               OPEN FIRST FILE IN XMIPP SELFILE LIST
                MAXIM = 0  
                CALL OPFILEC(LUNCP,.FALSE.,FILNAM(:NLETX),LUNIMG,DISPT,
     &                       ITYPE,NX,NY,NZ, 
     &                       MAXIM,'dum~9',FOUROK,IRTFLG) 
                IF (IRTFLG .NE. 0) GOTO 9000 

                IMGNUM     = -1
#if defined(SP_DBUGIO)
                write(3,*)' Opened XMipp selfile img: ',filnam(1:nletx)
#endif
                RETURN
             ENDIF
          ENDIF

          MAXIM = 0 
          IF (LOCAT > 0) MAXIM = 1   ! >= 1 STACKED IMAGE: e.g. imgstk@3
  
#if defined(SP_DBUGIO)
           write(3,*)'  '
           write(3,*)' In opfiles SF3; filpat: ',  filpat(:nlet)
           write(3,*)' In opfiles SF3; dispt:           |', dispt,'|'
           write(3,*)' In opfiles SF3; maxim:       ', maxim
           write(3,*)' In opfiles SF3; Calling opfilec '
           write(3,*)'  '
#endif
          CALL OPFILEC(LUNCP,.FALSE.,FILPAT,LUNIMG,DISPT,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,PROMPT,FOUROK,IRTFLG) 

#if defined(SP_DBUGIO)
          write(3,*)'   '
          write(3,*)' In opfiles SF4; Opened:    ', filpat(:nlet)
          write(3,*)' In opfiles SF4; maxim:     ', maxim
          write(3,*)' In opfiles SF4; itype:     ', itype
          write(3,*)' In opfiles SF4; irtflg:    ', irtflg
#endif
 
C         RETURN FILENAME WITH ANY EXTENSION IF NOT SPIDER IMAGE
          IF (IRTFLG == 5) NLET = lnblnkn(FILPAT)
          IF (IRTFLG .NE. 0) GOTO 9000 

          LIST_SIZE  = 0
          IMGNUM     = 1

#if defined(SP_DBUGIO_NEVER)
          write(3,*)'  '
          write(3,*)' In opfiles SF5: Opened simple: ', filpat(:nlet)
#endif
 
      ENDIF


#if defined (SP_DBUGIO)
      write(3,*)'   '
      write(3,*)' Opfiles 999 returns: maxim:           ', maxim 
      write(3,*)' Opfiles 999 returns: imgnum=list_val: ', imgnum 
      write(3,*)' Opfiles 999 returns: list_size:       ', list_size 
      write(3,*)' Opfiles 999 returns: irtflg:          ', irtflg 
      write(3,*)'    ' 
#endif


9000  RETURN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
C          |  Templated simple image:  img*** 
C          ` ---> FILGET --> OPFILEC  
C          |                         
C          |  Templated stacked image: stk@*     
C          ` ---> LUN***
C          | 
C          |  Whole image stack: stk@  
C          ` ---> LUNS***
C          | 
C          |  Templated stacked MRC: *@mrc or *@.mrcs
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

           !write(3,*)' Opened old templated file: ', filnam(1:nlet)

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
           write(3,*)' In getoldimg; filnam: ', filnam(1:nlet)
           write(3,*)' In getoldimg; nlet:   ', nlet
           write(3,*)' In getoldimg; locat:  ', locat
#endif
        ENDIF

#if defined (SP_DBUGIO)
           write(3,*)' In getoldimg; ngot:   ', ngot
           write(3,*)' In getoldimg; Opened filnam(1:nlet): ',
     &                                      filnam(1:nlet)
           write(3,*)' In getoldimg; Calling lunsayinfo'
           write(3,*) '  '
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
C          |  Templated simple image:  img*** 
C          ` ---> FILGET --> OPFILEC  
C          |                         
C          |  Templated stacked image: stk@*     
C          ` ---> LUN***
C          | 
C          |  Whole image stack: stk@  
C          ` ---> LUNS***  
C          | 
C          |  Templated stacked MRC image: file*@mrc or *@file.mrcs
C          ` ---> GETNEWIMG_MRC --> OPFILEC  
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
        write(3,*)'  '
        write(3,*)' In getnewimg 0 ; filpat:       ', trim(filpat)
        write(3,*)' In getnewimg 0 ; nwant :       ', nwant 
        write(3,*)' In getnewimg 0 ; locast,locat: ', locast,locat 
#endif

        IF (NWANT < 0) THEN
C          XMIPP SELFILE SIMPLE IMAGE ----------------------- selaaa

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
C          TEMPLATED SIMPLE IMAGE --------------------------- img***

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
           write(3,*)'  '
           write(3,*)' In getnewimg 2 ;  Opened templated img: ', 
     &                                   trim(filnam)
#endif

           NGOT   = NWANT   
           IRTFLG = 0
           RETURN

        ELSEIF (IS_MRC .AND. LOCAT > 0 .AND. LOCAST < LOCAT) THEN

C          TEMPLATED STACKED MRC IMAGE ---------------------- *@stk.mrc
C          BARE MRC IMAGE ------------------------------------ @stk.mrc
           
           CALL GETNEWIMG_MRC(LUN,FILPAT,NWANT,SAYIT,
     &                        FILNAM,NGOT,IRTFLG)

#if defined(SP_DBUGIO)
           write(3,*)'  '
           write(3,*)' In getnewimg 3 ;  Opened mrc img: ', 
     &                                   trim(filnam)
#endif
           RETURN


        ELSEIF (LOCAT > 0) THEN
C          BARE STACK ---------------------------------- stk@ or @stk

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
           write(3,*)'  '
           write(3,*)' In getnewimg 4 ;  maxim:   ', maxim
           write(3,*)' In getnewimg 4 ;  nwant:   ', nwant
           write(3,*)' In getnewimg 4 ;  istk:    ', istk
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
           write(3,*)'  '  
           write(3,*)' In getnewimg 5 ; Opened :   ', filnam(:nlet)
#endif

        ENDIF

        END






        SUBROUTINE ADD_TO_BARE(LUN,FILNAM,IMGWANT, 
     &                        SAYIT,NGOT,IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMLIMIT.INC'

        INTEGER                :: LUN
        CHARACTER(LEN=*)       :: FILNAM
        INTEGER                :: IMGWANT
        LOGICAL                :: SAYIT
        INTEGER                :: NGOT
        INTEGER                :: IRTFLG

        CHARACTER(LEN=1)       :: DSP
        LOGICAL                :: ISBARE
        INTEGER                :: NLET,LOCAT
        INTEGER                :: NX,NY,NZ,ITYPE,ISTK,IRTFLGT
        INTEGER                :: MAXIM

        INTEGER                :: lnblnkn     ! FUNCTION


        !IMGWANT  = IMGWANT          ! IMGWANT MAY NOT BE WRITABLE
        NLET   = lnblnkn(FILNAM)
        LOCAT  = INDEX(FILNAM(1:NLET),'@')

#if defined(SP_DBUGIO)
        write(3,*)'  '
        write(3,*)' In add_to_bare 0 ; filnam:       ', filnam(:nlet)
        write(3,*)' In add_to_bare 0 ; imgwant :     ', imgwant 
        write(3,*)' In add_to_bare 0 ; locast,locat: ', locat 
#endif

C       BARE STACK ---------------------------------- stk@ or @stk

        CALL LUNGETISBARE(LUN,ISBARE,IRTFLG)

C       LOAD OVERALL HEADER FIRST FOR LUNREDHED (MAY BE MT NOW!)
        CALL LUNGETSIZE(LUN,NX,NY,NZ,IRTFLG)
        CALL LUNREDHED(LUN,NX,0,.FALSE.,IRTFLG)
 
C       RETRIEVE CURRENT MAXIMUM IMAGE NUMBER FROM OVERALL HEADER
        CALL LUNGETMAXIM(LUN,MAXIM,IRTFLG)

        IF (IMGWANT > MAXIM) THEN
C          UPDATE OVERALL HEADER WITH MAXIMUM IMAGE NUMBER
           CALL LUNSETMAXIM(LUN,IMGWANT,IRTFLG)
           CALL LUNSETMAXALL(LUN,IMGWANT,IRTFLG)
        ENDIF

C       NEED ISTK 
        CALL LUNCOPYSTK(LUN,ISTK,IRTFLGT)

#if defined(SP_DBUGIO)
        write(3,*)'  '
        write(3,*)' In add_to_bare 4 ; maxim:   ', maxim
        write(3,*)' In add_to_bare 4 ; imgwant: ', imgwant
        write(3,*)' In add_to_bare 4 ; istk:    ', istk
#endif

        IF (IMGWANT > MAXIM .OR. ISTK > 2) THEN
C          SAVE OVERALL HEADER NOW TO PRESERVE MAXIM & LASTINDX
           CALL LUNWRTHED(LUN,NX,0,IRTFLGT)
           CALL LUNGETMAXIM(LUN,MAXIM,IRTFLG)
        ENDIF

C       GET FILENAME FROM CURRENT HEADER OBJECT
        CALL LUNGETFILE(LUN,FILNAM,NLET,DSP,IRTFLG)

C       SET NEW FILENAME IN HEADER OBJECT AND GET FILENAME
        CALL LUNSETIMNUM(LUN,FILNAM,IMGWANT,'N',IRTFLG)
        NLET = lnblnkn(FILNAM)    ! MAY HAVE BEEN ALTERED

C       SET IMGNUM IN HEADER OBJECT
        CALL LUNSETINUSE(LUN,IMGWANT,IRTFLG)

C       PUT IMAGE ISTACK
        CALL LUNSETISTACK(LUN,0,IRTFLG)

C       PUT COMMON VALUES INTO COMMON AREA (NOT NEEDED IN FUTURE?)
        CALL LUNSETCOMMON(LUN,IRTFLG)

C       PUSH HEADER OBJECT INFO INTO NEW STACKED FILE
        CALL LUNWRTHED(LUN,NX,IMGWANT,IRTFLG)

C       SET OFFSETS FOR REDLIN/WRTLIN ON THIS LUN
        CALL LUNSETIMGOFF(LUN,IMGWANT,NX,IRTFLGT)

C       WRITE OUT FILE OPENING INFO TO SCREEN
        CALL LUNSAYINFO(LUN,IRTFLG)

        NGOT = IMGWANT

#if defined(SP_DBUGIO)
        write(3,*)'  '  
        write(3,*)' In add_to_bare 5 ; Opened:  ', filnam(:nlet)
#endif

        END
