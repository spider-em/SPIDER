 
C++*********************************************************************
C                                                                      
C  OPFILES_MRC.F    CREATED FROM OPFILES         7/30/19  ArDean Leith
C                   COMMENTS                     2/06/20  ArDean Leith
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
C  CONTAINS: OPFILES_MRC, GETOLDIMG_MRC, GETNEWIMG_MRC 
C
C  OPFILES_MRC(LUNCP,LUNIMG,LUNDOC,LUNXM,  
C              ASKNAM,FILPAT,NLET, DISP,
C              ITYPE,NX,NY,NZ,MAXIM,FOUROK,
C              IMGNUM, IRTFLG)
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
C
C        MAXIM      STACK INDICATOR                          (SENT/RET)
C                   ON INPUT: STACK INDICATOR IF DISP == 'I'     (SENT)
C                   ON OUTPUT:                                    (RET)
C                       -2 NON-STACK IMAGE                
C                       -1 STACKED IMAGE                  
C                     >= 0 IS CURRENT MAX. IMAGE NO. FOR STACK             
C        PROMPT     PROMPT FOR FILNAME                            (SENT)
C                     IF NOT (ASKNAM) THIS IS FILE NAME           (SENT)
C                     ~ (TILDE) IN LAST CHAR. SAYS SKIP
C                       "FILE" AT END OF PROMPT
C                     ~9 IN NEXT TO LAST OR
C                        NEXT-TO-NEXT-TO LAST
C                        ACCEPTS AN EXTENSION
C                        (OTHERWISE DISCARDED!)
C                     ~6 KEEPS OLD DATE/TIME
C        IMGNUM     IMAGE NUMBER                              (SENT/RET)
C                     ON INPUT:   IMAGE NUMBER THAT IS WANTED
C                     ON OUTPUT:  >=0 IMAGE NUMBER CURRENTLY OPEN 
C                                 <0 IS SELFILE IN USE 
C        IRTFLG     ERROR FLAG (0 IS NORMAL)                      (RET)
C                      -1 GO TO PREVIOUS QUESTION
C 
C  CALL TREE:
C
C        OPFILES 
C           v
C        FILERD
C           v 
C        FILELIST
C           v     If: MRC file: *.MRC or *.MRCS
C           v 
C        OPFILES_MRC
C           |
C           |   Templated stack: **@STK.MRC
C           ` ---> FILGET_AT --> OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC 
C           |                            --> OPENFIL_MRC --> OPENFIL_N_MRC
C           |
C           |   Whole barestack:  @STK.MRC 
C           ` ---> FILGET_AT --> OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC
C           |                            --> OPENFIL_MRC --> OPENFIL_N_MRC
C           | 
C           |   File template:   IMG***.MRC
C           ` ---> FILGET_AT --> OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC
C           |                            --> OPENFIL_MRC --> OPENFIL_N_MRC
C           | 
C           |   Simple file:     IMG001.MRC
C           ` -----------------> OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC  
C                                        --> OPENFIL_MRC --> OPENFIL_N_MRC
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12 
C--********************************************************************* 

       SUBROUTINE OPFILES_MRC(LUNCP,LUNIMG,LUNDOC,
     &                        FILPAT,NLET, DISP,
     &                        ITYPE,NX,NY,NZ,MAXIM,
     &                        IMGNUM, IRTFLG) 
 
       IMPLICIT NONE

       INCLUDE 'CMBLOCK.INC' 
       INCLUDE 'CMLIMIT.INC' 
 
       INTEGER                   :: LUNCP,LUNIMG,LUNDOC
       CHARACTER(LEN=*)          :: FILPAT 
       INTEGER                   :: NLET
       CHARACTER(LEN=1)          :: DISP 
       INTEGER                   :: ITYPE,NX,NY,NZ,MAXIM
       LOGICAL                   :: FOUROK = .FALSE.
       INTEGER                   :: IMGNUM,IRTFLG

       LOGICAL                   :: SAYIT
       CHARACTER (LEN=MAXNAM)    :: FILNAM,FILNAMNOAT 
       CHARACTER (LEN=2*MAXNAM)  :: MESG 
       CHARACTER (LEN=1)         :: NULL = CHAR(0)

       INTEGER                   :: LOCAT,LOCAST,LENE
       INTEGER                   :: LUNOP,NLETT,I

       LOGICAL                   :: ISOPEN,IS_BAREL
       CHARACTER (LEN=MAXNAM)    :: FILNAMT 

       INTEGER                   :: lnblnkn       ! FUNCTION
       integer                   :: itypet,nxt,nyt,nzt,maximt


       LOCAT    = INDEX(FILPAT(1:NLET),'@')   
       LOCAST   = INDEX(FILPAT(1:NLET),'*')
       IS_BAREL = (LOCAT > 1 .AND. FILPAT(LOCAT-1:LOCAT-1) == '/')
       
       !write(3,*)' In opfiles_mrc,locast,locat,:',locast,locat,is_barel
       !write(3,*)' In opfiles_mrc, imgnum:',imgnum,filpat

       IF (LOCAST > 0 .AND. LOCAST < LOCAT) THEN
C         TEMPLATED STACKED MRC FILE: ***@STK.MRC --------- **@STK.MRC

C         PRINT *, "opfiles_mrc.f : 160: Templated stacked MRC file"  ! get segfault
          !write(3,*)' In opfiles_mrc, opening: ',filpat(1:nlet)
        
C         SUBSTITUTE STACKED IMGNUM INTO FILE NAME PATTERN -> FILNAM 
          CALL FILGET_AT(FILPAT,IMGNUM,FILNAM,NLET,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

C         OPEN IMGNUM IN EXISTING MRC STACK FILE
          SAYIT  = .TRUE.
          MAXIM  = 0          ! OPEN SPECIFIC IMAGE  

          !write(3,*)' In opfiles_mrc, opening: ',filnam(1:nlet)

C         PRINT *, "opfiles_mrc.f : 174: OPFILES_MRC: Calling OPFILEC"
          CALL OPFILEC(0,.FALSE.,FILNAM,LUNIMG,DISP,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,' ',FOUROK,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

C         THIS IS NOT A BARE STACK REQUEST
          CALL LUNSETISBARE_MRC(LUNIMG,.FALSE.,IRTFLG)

          !write(3,*)' In opfiles_mrc, maxim: ',maxim,imgnum
          RETURN

       ELSEIF (LOCAST > 0) THEN
C         A SIMPLE FILE TEMPLATE: IMG***.mrc ------------- IMG***.mrc
C         PRINT *, "opfiles_mrc.f : 186: Simple templated MRC files"

C         SUBSTITUTE IMGNUM INTO FILE NAME PATTERN -> FILNAM   
          CALL  FILGET(FILPAT,FILNAM,NLET,IMGNUM,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

C         OPEN THE FILE
          MAXIM = 0 
          CALL OPFILEC(LUNCP,.FALSE.,FILNAM,LUNIMG,DISP,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,' ',FOUROK,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

C         THIS IS NOT A BARE STACK REQUEST
          CALL LUNSETISBARE_MRC(LUNIMG,.FALSE.,IRTFLG)

c          write(3,*)' Opened templated file: ',filnam(1:20),
c     &              '  for:',imgnum

       ELSEIF (LOCAT == 1 .OR. IS_BAREL) THEN
C         WHOLE BARESTACK:  @STK.mrc  ---------------------- @STK.mrc
C         PRINT *, "opfiles_mrc.f : 208: Whole barestack MRC file"

C         SUBSTITUTE STACKED IMGNUM INTO FILE NAME PATTERN -> FILNAM   
          
          IF (LOCAT == 1) THEN
              FILNAMNOAT = FILPAT(2:NLET)
          ELSE
              FILNAMNOAT = FILPAT(1:LOCAT-1) // FILPAT(LOCAT+1:NLET)  
          ENDIF

C         OPEN FIRST FILE IN MRC STACK
          IMGNUM = 1
          MAXIM  = 1
         !write(3,*)' In opfiles_mrc, filnamnoat: ',filnamnoat(1:nlet)

          INQUIRE(FILE=FILNAMNOAT,OPENED=ISOPEN,NUMBER=LUNOP)
          MESG = '  FILE: ' // FILNAMNOAT(1:NLET)//'  ALREADY OPENED ON' 
          LENE = LNBLNKN(MESG)
          IF (ISOPEN .AND. LUNOP .NE. LUNIMG) THEN
             WRITE(NOUT,'(A,I3)') MESG(1:LENE),LUNOP 
             IRTFLG = -2
             RETURN
          ENDIF

C         SUBSTITUTE STACKED IMGNUM INTO FILE NAME PATTERN -> FILNAM   
          CALL FILGET_AT(FILPAT,IMGNUM,FILNAM,NLET,IRTFLG)
          !write(3,*)' In opfiles_mrc, filpat: ',filpat(1:nlet)
          !write(3,*)' In opfiles_mrc, filnam: ',filnam(1:nlet)
        
C         OPEN FIRST FILE
          CALL OPFILEC(LUNCP,.FALSE.,FILNAM,LUNIMG,DISP,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,'UNUSED',FOUROK,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

C         THIS WAS A BARE STACK REQUEST
          CALL LUNSETISBARE_MRC(LUNIMG,.TRUE.,IRTFLG)

          !write(3,*)' In opfiles_mrc, opened: ',filnam(1:nlet)
          !write(3,*)' In opfiles_mrc, imgnum: ',imgnum

       ELSE
C         SINGLE SIMPLE INPUT FILE: ------------------------ IMG001.mrc
C         STACKED IMAGE FILE:  ------------------------------ 2@STK.mrc
                
          MAXIM = 0        ! SPECIFIC IMAGE WANTED  
          CALL OPFILEC(LUNCP,.FALSE.,FILPAT,LUNIMG,DISP,
     &                 ITYPE,NX,NY,NZ, 
     &                 MAXIM,'UNUSED',FOUROK,IRTFLG)

C         RETURN FILENAME WITH ANY EXTENSION IF NOT SPIDER IMAGE
          IF (IRTFLG == 5) NLET = lnblnkn(FILPAT)
          IF (IRTFLG .NE. 0) RETURN

          ! write(3,*)' In opfiles_mrc, imami : ',imami
          ! write(3,*)' In opfiles_mrc, openrf simple file: ',filpat

          IMGNUM = 1

       ENDIF

       END 
 

C **********************************************************************
C
C    GETNEWIMG_MRC( LUN, FILPAT,NWANT,SAYIT,FILNAM,NGOT,IRTFLG)
C
C    PURPOSE:       OPEN A SPECIFIED TEMPLATED MRC IMAGE FOR 
C                   RANDOM ACCESS READING/WRITING.
C
C    PARAMETERS:
C        LUN        I/O UNIT NUMBER FOR FILNAM.                   (SENT)
C        FILPAT     FILENAME PATTERN                              (SENT)
C        NWANT      IMAGE NUMBER WANTED                           (SENT) 
C        SAYIT      SAY FILE OPENING INFO                         (SENT)
C        FILNAM     FILENAME OPENED                               (RET.)
C        NGOT       IMAGE NUMBER FOUND                            (RET.) 
C        IRTFLG     ERROR RETURN FLAG.                            (RET.)
C                   IRTFLG = -1    END OF FILE BEFORE NWANT
C                   IRTFLG =  0    NORMAL RETURN, IMAGE IS STACK
C                   IRTFLG =  2    IMAGE NOT IN USE
C
C **********************************************************************

        SUBROUTINE GETNEWIMG_MRC(LUN,FILPAT,NWANT,SAYIT,
     &                           FILNAM,NGOT,IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        INTEGER                :: LUN
        CHARACTER(LEN=*)       :: FILPAT
        INTEGER                :: NWANT
        LOGICAL                :: SAYIT
        CHARACTER(LEN=*)       :: FILNAM
        INTEGER                :: NGOT
        INTEGER                :: IRTFLG

        INTEGER                :: LUNCP = 0    ! UNUSED FOR MRC
        CHARACTER(LEN=1)       :: DSP
        CHARACTER(LEN=4)       :: CAXIS
        LOGICAL                :: FOUROK = .FALSE.
        INTEGER                :: NLET,LOCAST,LOCAT,NX,NY,NZ,ITYPE
        INTEGER                :: MAXIM,MZ,NSTACK,IMGNUM,LENT,IMGNUMOLD
        LOGICAL                :: WANTUL

        INTEGER                :: lnblnkn    ! FUNCTION

        NLET   = lnblnkn(FILPAT)

        LOCAST = INDEX(FILPAT(1:NLET),'*')
        LOCAT  = INDEX(FILPAT(1:NLET),'@')

        !write(3,*)' !!!locast,locat,filpat:',locast,locat,filpat(1:nlet)

        IMGNUM = NWANT

        IF (LOCAST > 1 .or. (LOCAT  > 0 .AND. LOCAST < LOCAT)) THEN
C          SUBSTITUTE IMGNUM INTO FILE NAME PATTERN -> FILNAM   
           CALL FILGET_AT(FILPAT,IMGNUM,FILNAM,NLET,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ENDIF

        IF (LOCAT <= 0 .AND. LOCAST > 1) THEN
C          TEMPLATED SIMPLE MRC IMAGE ---------------------- IMG***.mrc

C          NEW IMAGE, NEEDS TO KNOW: ITYPE,NX,NY,NZ!
C          GET IT FROM OPFILES OR PREVIOUS CALL

           CALL LUNGETSIZE_MRC(LUN,NX,NY,NZ,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
           CALL LUNGETTYPE_MRC(LUN,ITYPE,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           CLOSE(LUN)    ! MAY BE STILL OPEN FROM FIRST CALL  

C          OPEN NEW FILE
           MAXIM = 0
           CALL OPFILEC(LUNCP,.FALSE.,FILNAM,LUN,'U',ITYPE,
     &                  NX,NY,NZ, 
     &                  MAXIM,'UNUSED',FOUROK,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           !write(3,*)' In getnewimg_mrc opened new templated file: ',
           !&                 filnam(1:nlet)

           NGOT   = NWANT   
           IRTFLG = 0
           RETURN

        ELSEIF (LOCAT  > 0 .AND. LOCAST < LOCAT) THEN

C          TEMPLATED STACKED MRC FILE: ***@STK.MRC --------- **@STK.MRC
C          BARE STACKED MRC FILE: @STK.MRC ----------------- **@STK.MRC
           
C          THIS STACK FILE ALREADY OPEN JUST NEED TO UPDATE FOR IMAGE

C          SET UNDETERMINED IMAGE STATISTICS
           CALL LUNSET_STATSIMG_MRC(LUN,0,IRTFLG)

C          GET CURRENT IMGNUM AND NSTACK FROM STATIC HEADER
           CALL LUNGETSTK_MRC(LUN,MZ,NSTACK,IMGNUMOLD,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (IMGNUM > NSTACK) THEN 
C             SET CURRENT IMGNUM AND NSTACK IN STATIC HEADER
              NSTACK = IMGNUM
              
              CALL LUNSETSTK_MRC(LUN,IMGNUM,NSTACK,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

              CALL LUNGETSIZE_MRC(LUN,NX,NY,NZ,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

              CALL LUNSETNSTACK_MRC(LUN, NZ,NSTACK,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

           ENDIF           

C          SET OFFSETS FOR REDLIN/WRTLIN FOR IMAGE OPENED ON THIS LUN
           CALL LUNSETPOS_MRC(LUN,IMGNUM,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
           !write(3,'(a,8i5)') ' In getnu_mrc22, lun,imgnum:',lun,imgnum


C          PUSH HEADER OBJECT INTO FILE TO SAVE STATS, NSTACK, & IMGNUM
           CALL LUNWRTHED_MRC(LUN,IRTFLG)
           IF (IRTFLG .NE. 0) THEN
             LENT = lnblnkn(FILNAM)
             WRITE(NOUT,99) IRTFLG,LUN,FILNAM(:LENT)
99           FORMAT( '  *** ERROR(',I4,') ON UNIT: ',I3,' FILE: ',A)
             RETURN
           ENDIF

C          PUT COMMON VALUES INTO COMMON AREA (NOT NEEDED IN FUTURE?)
           CALL LUNSETCOMMON_MRC(LUN,IRTFLG)

C          WRITE OUT FILE OPENING INFO TO SCREEN
           CALL LUNSAYINFO_MRC(LUN,.TRUE.,IRTFLG)

           NGOT = IMGNUM

           !write(3,*)' In getnewimg_mrc opened new stacked file: ',
       !&                 filnam(1:nlet)

        ENDIF

        END


C++*********************************************************************
C
C GETOLDIMG_MRC.F   FROM GETNXTSTK               AUG 2019 ArDean Leith
C
C **********************************************************************
C
C    GETOLDIMG_MRC(LUN,FILPAT,NWANT,SAYIT,FILNAM,NGOT,IRTFLG)
C
C    PURPOSE:       TO OPEN A SPECIFIED IMAGE WITHIN MRC STACK FOR 
C                   RANDOM ACCESS READING/WRITING.
C
C    PARAMETERS:
C        LUN        LUN NUMBER FOR FILNAM                         (SENT)
C        FILPAT     FILE NAME PATTERN                             (SENT)
C        NWANT      IMAGE NUMBER WANTED (<0 IS SELFILE)           (SENT) 
C        SAYIT      ECHO IMAGE OPENING DETAILS                    (SENT) 
C        FILNAME    FILE NAME OPENED                              (RET)
C        NGOT       IMAGE NUMBER OPENED                           (RET.) 
C        IRTFLG     ERROR RETURN FLAG.                            (RET.)
C                   IRTFLG =  0    NORMAL RETURN, IMAGE IS STACK
C                   IRTFLG =  1    FILE ENDS BEFORE NWANT
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************
 
        SUBROUTINE GETOLDIMG_MRC(LUN,FILPAT,NWANT,SAYIT,
     &                           FILNAM,NGOT,IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMLIMIT.INC'

        INTEGER               :: LUN     
        CHARACTER(LEN=*)      :: FILPAT
        INTEGER               :: NWANT 
        LOGICAL               :: SAYIT    
        CHARACTER(LEN=MAXNAM) :: FILNAM
        INTEGER               :: NGOT     
        INTEGER               :: IRTFLG     

        CHARACTER(LEN=1)      :: NULL = CHAR(0)
        CHARACTER(LEN=4)      :: CAXIS

        INTEGER               :: NLET,LOCAST,LOCAT,MAXIM,ITYPE
        INTEGER               :: NX,NY,NZ,MZ,NSTACK,IMGNUM
        LOGICAL               :: FOUROK = .FALSE. ! NO MRC FOURIER SUPPORT
        LOGICAL               :: WANTUL

        INTEGER               :: lnblnkn      ! FUNCTION

        NLET   = lnblnkn(FILPAT)

        LOCAST = INDEX(FILPAT(1:NLET),'*')
        LOCAT  = INDEX(FILPAT(1:NLET),'@')

        ! write(3,*)' In getoldimg_mrc locast,locat,nlet,filpat:',
        !&            locast,locat,nlet,filpat(1:nlet)
        !write(3,*)' In getoldimg_mrc, nwant,: ',nwant,':',filpat(1:nlet)

        IF (LOCAST > 0 .OR. (LOCAST == 0 .AND. LOCAT > 0)) THEN 
C          SUBSTITUTE STACKED IMGNUM INTO FILE NAME PATTERN -> FILNAM   
           CALL FILGET_AT(FILPAT,NWANT,FILNAM,NLET,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ENDIF

        IF (LOCAST > 0 .AND. LOCAT == 0) THEN 
C          TEMPLATED SIMPLE MRC IMAGE ---------------------- IMG***.mrc

           CLOSE(LUN)    ! USUALLY STILL OPEN
 
C          OPFILEC CALLS OPENFIL_MRC --> OPENFIL_O_MRC -->
C          LUNSETFILE, LUNSETPOS_MRC LUNSET_STATSIMG_MRC,
C          LUNSETISBARE_MRC, LUNSETCOMMON_MRC, LUNSAYINFO_MRC

           MAXIM = 0  
           CALL OPFILEC(0,.FALSE.,FILNAM,LUN,NULL,ITYPE,
     &                  NX,NY,NZ, 
     &                  MAXIM,' ',FOUROK,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           !write(6,*)' In getoldimg_mrc opened old templated file: ',
           !          filpat(1:nlet), maxim

           NGOT   = NWANT   
           IRTFLG = 0
           RETURN

        ELSEIF (LOCAT  > 0 .AND. LOCAST < LOCAT) THEN
C          TEMPLATED OLD STACKED MRC IMAGE ----------------- **@STK.mrc
C                    OLD BARE MRC STACKED IMAGE ------------   @STK.mrc

C          FIND NUMBER OF IMAGES IN STACK
           CALL LUNGETSTK_MRC(LUN,MZ,NSTACK,IMGNUM,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
                
C          CHECK IF IMAGE CAN FIT IN STACK
           IF (NWANT > NSTACK) THEN
             CALL ERRT(102,'GETOLDIMG_MRC, NO SUCH IMAGE',NWANT)
             IRTFLG = 1
             RETURN
           ENDIF

C          SET IMAGE NUMBER IN STATIC FILE HEADER AREA
           CALL LUNSETIMGNUM_MRC(LUN,NWANT,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          SET FILENAME,  IN LUN - FILENAME
           CALL LUNSETFILE_MRC(LUN,FILNAM(1:NLET),'O',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          SET AXIS ORIGIN LOCATION & VOLUME HANDEDNESS IF NEEDED
           CALL WHICH_HAND_MRC(LUN,FILNAM,CAXIS,IRTFLG)

C          SET IMAGE OFFSET FOR READ/WRITE
           CALL LUNSETPOS_MRC(LUN,NWANT,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          PUT COMMON VALUES INTO COMMON AREA (NOT NEEDED IN FUTURE?)
           CALL LUNSETCOMMON_MRC(LUN,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
           !write(3,*)' In getoldimg_mrc after lunsetcommon: ',filnam(1:nlet)

C          WRITE OUT FILE OPENING INFO IF WANTED
           CALL LUNSAYINFO_MRC(LUN,SAYIT,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

cc         write(3,*)' In getoldimg_mrc opened old templated stk file: ',
cc     &               filnam(1:nlet)

C          RETURN IMAGE NUMBER THAT WAS WANTED
           NGOT   = NWANT
   
           IRTFLG = 0
           RETURN

        ENDIF
        IRTFLG = 1

        END

