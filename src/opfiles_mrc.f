
C++*********************************************************************
C                                                                      
C  OPFILES_MRC.F    CREATED FROM OPFILES         7/30/19  ArDean Leith
C                   COMMENTS                     2/06/20  ArDean Leith
C                   DEBUG OUTPUT ADDED           9/26/25  ArDean Leith
C.                  REWRITEN                    11/28/25  ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2025  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=*                                                                    *
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
C  OPFILES_MRC(LUNCP,LUNIMG,LUNDOC,
C              FILPAT,NLET, DSP,
C              ITYPE,NX,NY,NZ, NSTK, ISTK, IRTFLG)
C
C  PURPOSE: SOLICITS FILE NAME(S) AND OPENS FILE(S)
C           SUPPORT ROUTINE FOR CONVERTING OPERATIONS TO 
C           WORK ON WHOLE STACK OR WITH SELECTION DOC FILE.
C
C  PARAMETERS:
C        LUNCP    UNIT TO COPY HEADER VALUES FROM               (SENT)
C        LUNIMG   UNIT TO OPEN FILE ON                          (SENT)
C        LUNDOC   UNIT TO OPEN LIST DOC FILES ON                (SENT)
C
C        FILPAT   FILE NAME PATTERN                              (RET)
C        NLET     CHARS IN FILE NAME PATTERN                     (RET)
C
C        DSP      CHARACTER CONTAINING ONE OF THE               (SENT)
C                 FOLLOWING DISPOSITION SPECIFICATIONS:
C                   'O'   -  FILE IS ASSUMED TO EXIST.  DIMENSIONS,
C                            ITYPE AND HEADER INFO (IN COMMON) ARE 
C                            RETURNED TO THE CALLING PROGRAM. 
C                   'B'   -  SAME AS OLD BUT NO LIMIT ON BUFFER
C                            LENGTH FOR OPENCHK. 
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
C        NX,NY,NZ  IMAGE SIZE                            (SENT OR RET)
C
C        NSTK      STACK SIZE                            (SENT OR RET)
C
C        LIST_SIZE NUMBER OF ENTRIES IN LIST               SENT / RET)
C        ISTK      CURRENT VALUE FROM ILIST                     (SENT) 
C
C        IRTFLG    ERROR FLAG (0 IS NORMAL)                      (RET)
C                      -1 GO TO PREVIOUS QUESTION
C
C GONE   LIST_NUM  CURRENT LOCATION IN ILIST  (UNUSED)            (?) 
C GONE   ISTK      FILE or ISTK NUMBER      (INTERNAL USED ONLY)  (?)
C
C 
C  CALL TREE:
C
C    OPFILES 
C      |
C      ` -> FILERD  
C      |       ` ->  ECHONAME
C      ` -> FILELIST
C      |
C      |  If: MRC file: *.MRC or *.MRCS
C      ` -> OPFILES_MRC
C             |
c           OPENFIL_MRC
C             | If: Old
C             |    `-> OPENFIL_O_MRC -->
C             |
C             | OR: New
C             |    `-> OPENFIL_N_MRC --> 
C             |
C           LUNSET_ISBARE_MRC
C             |
C             | IF (ISTK > 0) 
C             |     ` -> LUNSET_STK_260_MRC
C             |
C             | IF (ISTK > 0) --> LUNSET_STK_260_MRC
C             |     `LUNSET_STK_260_MRC (different input)
C             |
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

       SUBROUTINE OPFILES_MRC(LUNCP,LUNIMG,LUNDOC,
     &                        FILPAT,NLET, DSP,
     &                        ITYPE,NX,NY,NZ, NSTK, 
     &                        ILIST, LIST_SIZE, ISTK, IRTFLG)

C      OPFILES_MRC     NEEDS: ITYPE,   ISTK
C      OPFILE (CALLER) NEEDS: IF_BARE: ILIST,NTOT=LIST_SIZE,
C                                      IMGNUM=ISTK  

       IMPLICIT NONE

       INCLUDE 'CMBLOCK.INC' 
       INCLUDE 'CMLIMIT.INC' 
 
       INTEGER                   :: LUNCP,LUNIMG,LUNDOC
       CHARACTER(LEN=*)          :: FILPAT 
       INTEGER                   :: NLET
       CHARACTER(LEN=1)          :: DSP 
       INTEGER                   :: ITYPE, NX,NY,NZ, NSTK
       INTEGER                   :: ILIST(*)
       INTEGER                   :: LIST_SIZE
       INTEGER                   :: IRTFLG

       CHARACTER (LEN=MAXNAM)    :: FIL_NOAT,FIL_DIRS,FIL_BASE,FIL_EXT 

       LOGICAL                   :: IS_MRC, IS_MRCS, IS_BARE

       INTEGER                   :: N, ISTK
       INTEGER                   :: lnblnkn       ! FUNCTION

       IF (LUNIMG <= 0 .OR. LUNIMG > 100) THEN
          CALL ERRT(102,'PGM. ERROR: LUN MUST BE 1...100',LUNIMG)
          RETURN
       ENDIF

#if defined (SP_DBUGIO)
        write(3,*)' '
        write(3,*)' In opfiles_mrc; dsp:       ', dsp
        write(3,*)' In opfiles_mrc, filpat:    ', trim(filpat)
       !write(3,*)' In opfiles_mrc; nstk:      ', nstk
        write(3,*)' In opfiles_mrc; list_size: ', list_size
        write(3,*)' In opfiles_mrc; istk:      ', istk
#endif


       NLET   = lnblnkn(FILPAT)

       NSTK   =  0   ! UNUSED BY OPENFIL_MRC  ON INPUT

       CALL OPENFIL_MRC(LUNIMG,FILPAT,DSP, 
     &                  NX,NY,NZ, 
     &                  NSTK, ISTK, IS_BARE,
     &                  ITYPE, IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

     
       IF (IS_BARE) THEN
C         FILL ILIST WITH LIST OF ALL POSSIBLE IMAGES IN BARESTACK
          IF (NSTK > NIMAX) THEN
             CALL ERRT(46,'OPFILES_MRC; ILIST...',NIMAX)
             RETURN
          ENDIF

          DO N = 1,NSTK
             ILIST(N) = N
          ENDDO
      
          LIST_SIZE = NSTK
          ISTK      = 1
       ENDIF

C      SET BARE STACK REQUEST
       CALL LUNSET_ISBARE_MRC(LUNIMG,IS_BARE,IRTFLG)

C      SET STATIC NSTK & ISTK
       CALL LUNSET_STK_260_MRC(LUNIMG,ISTK,NSTK,IRTFLG)


#if defined (SP_DBUGIO)
       write(3,*)' In opfiles_mrc; Returning istk,nstk: ', istk,nstk
       write(3,*)' In opfiles_mrc; Returning nstk:      ', nstk
       write(3,*)' In opfiles_mrc; Returning filpat:    ', trim(filpat)
       write(3,*)
#endif
       END 
 




      SUBROUTINE NEXTFILE_MRC(LIST_NINDX, ILIST, 
     &                        FOUROK,     LUNXM1,
     &                        LIST_SIZE,  NSTK,   
     &                        LUN1,       LUNCP, 
     &                        FILPAT,     DSP,
     &                        ISTK,       IRTFLG) 
 
      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC' 
      INCLUDE 'CMLIMIT.INC' 

      INTEGER           :: LIST_NINDX 
      INTEGER           :: ILIST(LIST_SIZE) 
      LOGICAL           :: FOUROK
      INTEGER           :: LUNXM1 
      INTEGER           :: LIST_SIZE 
      INTEGER           :: NSTK, LUN1,LUNCP 
      CHARACTER(LEN=*)  :: FILPAT
      CHARACTER(LEN=1)  :: DSP
      INTEGER           :: ISTK
      INTEGER           :: IRTFLG

      INTEGER           :: LIST_VAL, it
      LOGICAL           :: SAYIT = .TRUE.
      LOGICAL           :: GOTAT,  GOTAST, IS_BARE, IS_MRC

      CHARACTER (LEN=MAXNAM) :: FILNAM_RET
      LOGICAL           :: ISMRCFILE    ! FUNCTION
   
      integer :: mz,istk_b4 

           
#if defined (SP_DBUGIO)
      write(3,*)' '
      write(3,*)' In nextfile_mrc, list_nindx,list_size: ',
     &                             list_nindx,list_size
#endif

      FILNAM_RET = CHAR(0)       

      IF (LIST_NINDX > LIST_SIZE) THEN
C        OVERUN I/O LIST
         IRTFLG = -1
         RETURN
      ENDIF

      IS_MRC  = ISMRCFILE(FILPAT)
      GOTAT   = (INDEX(FILPAT,'@') > 0)
      GOTAST  = (INDEX(FILPAT,'*') > 0)

C     IS THIS IS A BARE STACK OPERATION?
      CALL LUNGETISBARE(LUN1,IS_BARE,IRTFLG)

#if defined (SP_DBUGIO)
      write(3,*)' In nextfile_mrc, gotat,gotast,is_bare: ',
     &                             gotat,gotast,is_bare
      write(3,*)' In nextfile_mrc, istk,irtflg: ',
     &                             istk,irtflg
#endif
      IF (IRTFLG .NE. 0) RETURN


      IF (ISTK == -1  .AND. LUNXM1 > 0 ) THEN
C        XMIPP SELFILE LISTED IMAGE
         LIST_VAL = -1
      ELSE
         LIST_VAL = ILIST(LIST_NINDX)
      ENDIF

#if defined (SP_DBUGIO)
         write(3,*)' In nextfile_mrc, list_val:', list_val
#endif


      IF (DSP == 'O' .OR. DSP == 'B' .OR. 
     &    DSP == 'Z' .OR. DSP == 'E') THEN 
  
C        OPEN NEXT INPUT FILE



#if defined (SP_DBUGIO)
C       find number of images in stack
        call lunget_stk_260_mrc(lun1,mz,nstk,istk_b4,irtflg)
        if (irtflg .ne. 0) return
                
        write(3,*)' in nextfile_mrc, 6666 ; lun1: ', lun1

        write(3,*)' in nextfile_mrc, 6666 ; nstk,istk_b4: ',
     &                                      nstk,istk_b4
#endif




#if defined (SP_DBUGIO)
         write(3,*)' In nextfile_mrc; call getoldimg_mrc'
#endif

         CALL GETOLDIMG_MRC(LUN1,FILPAT,LIST_VAL,SAYIT,
     &                     FILNAM_RET,ISTK,IRTFLG)

#if defined (SP_DBUGIO)
         write(3,*)' In nextfile_mrc; opened: ', trim(filnam_ret)
         write(3,*)' In nextfile_mrc; irtflg: ', irtflg
#endif

         IF (IRTFLG    < 0) RETURN    ! END OF WHOLE-STACK
         IF (IRTFLG .NE. 0) RETURN    ! ERROR

      ELSE   

C        OPEN NEXT OUTPUT FILE 
#if defined (SP_DBUGIO)
         write(3,*)' In nextfile_mrc; call getnewimg_mrc'
#endif

         CALL GETNEWIMG_MRC(LUN1,FILPAT,LIST_VAL,SAYIT,
     &                            FILNAM_RET, ISTK,IRTFLG)

#if defined (SP_DBUGIO)
         write(3,*)' In nextfile_mrc; opened: ', trim(filnam_ret)
         write(3,*)' In nextfile_mrc; irtflg: ',irtflg
#endif

         IF (IRTFLG .NE. 0) RETURN   ! ERROR

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
C        LUN        LUN NUMBER FOR FILNAM                        (SENT)
C        FILPAT     FILE NAME PATTERN                            (SENT)
C        NWANT      IMAGE NUMBER WANTED (<0 IS SELFILE)          (SENT)
C        SAYIT      ECHO IMAGE OPENING DETAILS                   (SENT)
C        FILNAM     FILE NAME OPENED                             (RET)
C        NGOT       IMAGE NUMBER OPENED                          (RET.)
C        IRTFLG     ERROR RETURN FLAG.                           (RET.)
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

        INTEGER               :: NLET,LOCAST,LOCAT,ITYPE
        INTEGER               :: ISTK_B4,ISTK
        INTEGER               :: NX,NY,NZ,MZ, NSTK
        LOGICAL               :: FOUROK = .FALSE. ! NO MRC FOURIER 
        LOGICAL               :: WANTUL

        INTEGER               :: lnblnkn      ! FUNCTION

        NLET   = lnblnkn(FILPAT)

        LOCAST = INDEX(FILPAT(1:NLET),'*')
        LOCAT  = INDEX(FILPAT(1:NLET),'@')

#if defined (SP_DBUGIO)
         write(3,*)'  '
         write(3,*)' In getoldimg_mrc, filpat: ',trim(filpat)
         write(3,*)' In getoldimg_mrc, filnam: ',trim(filnam)
         write(3,*)' In getoldimg_mrc, nwant,locast,locat: ',
     &                                 nwant,locast,locat
#endif


#if defined (FUTURE) 
        IF (LOCAST > 0 .OR. (LOCAST == 0 .AND. LOCAT > 0)) THEN 
C          SUBSTITUTE STACKED NWANT INTO FILE NAME PATTERN -> FILNAM   
           CALL FILGET_AT(FILPAT,NWANT,FILNAM,NLET,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN 
        ENDIF


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
     &                  NX,NY,NZ, 
     &                  MAXIM,'~9',FOUROK,IRTFLG) 
           IF (IRTFLG .NE. 0) RETURN 

C          NEW IMAGE SIZE SHOULD BE SAME AS PREVIOUS FILE
           CALL SIZCHK(NULL,NX1,NY1,NZ1,ITYPE1,
     &                      NX ,NY, NZ, ITYPE, IRTFLG)
           IF (IRTFLG .NE. 0) RETURN 
 
#if defined (SP_DBUGIO)
           write(3,*)' Opened old xmipp selfile file'
           write(3,*)' In getoldimg_mrc, 1111: ',trim(filnam)
#endif

           NGOT   = NWANT   
           IRTFLG = 0
           RETURN
#endif


        IF (LOCAST > 0 .AND. LOCAT == 0) THEN 
C          TEMPLATED SIMPLE MRC IMAGE ---------------------- IMG***.mrc
C          NOT A STACK!

           CLOSE(LUN)    ! USUALLY STILL OPEN
 
C          OPFILEC CALLS OPENFIL_MRC --> OPENFIL_O_MRC -->
C          LUNSET_FILE, LUNSET_POS_MRC LUNSET_STATSIMG_MRC,
C          LUNSET_ISBARE_MRC, LUNSET_COMMON_MRC, LUNSAYINFO_MRC

           istk = 0  
           CALL OPFILEC(0,.FALSE.,FILNAM,LUN,NULL,ITYPE,
     &                  NX,NY,NZ, 
     &                  NWANT,' ',FOUROK,IRTFLG) 
          IF (IRTFLG .NE. 0) RETURN 

#if defined (SP_DBUGIO)
           write(3,*)' In getoldimg_mrc, 2222: ',trim(filnam)
#endif
           NGOT   = NWANT   
           IRTFLG = 0
           RETURN




        ELSEIF (LOCAT  > 0 .AND. LOCAST < LOCAT) THEN
C          TEMPLATED OLD STACKED MRC IMAGE ----------------- *@STK.mrc

           ISTK = NWANT   ! FOR CLARITY

C          FIND NUMBER OF IMAGES IN STACK
           CALL LUNGET_STK_260_MRC(LUN,MZ,NSTK,ISTK_B4,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
                
#if defined (SP_DBUGIO)
           write(3,*)' In getoldimg_mrc, 6666 ; nstk,istk_b4: ',
     &                                          nstk,istk_b4
           write(3,*)' In getoldimg_mrc,      ; istk,nwant: ',
     &                                          istk,nwant
#endif

C          CHECK IF IMAGE CAN FIT IN STACK
           IF ( ISTK > NSTK) THEN
             CALL ERRT(102,'GETOLDIMG_MRC,  NO SUCH IMAGE',ISTK)
             IRTFLG = 1
             RETURN
           ENDIF

C          SET IMAGE NUMBER IN STATIC FILE HEADER AREA
           CALL LUNSET_IMGNUM_MRC(LUN,ISTK,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          SET FILENAME,  IN LUN - FILENAME
           CALL LUNSET_FILE_MRC(LUN,TRIM(FILNAM), 'O',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          SET AXIS ORIGIN LOCATION & VOLUME HANDEDNESS IF NEEDED
           CALL WHICH_HAND_MRC(LUN,FILNAM,CAXIS,IRTFLG)

C          SET IMAGE  OFFSET FOR READ/WRITE
           CALL LUNSET_POS_MRC(LUN,ISTK,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          PUT COMMON VALUES INTO COMMON AREA (NOT NEEDED IN FUTURE?)
           CALL LUNSET_COMMON_MRC(LUN,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

#if defined (SP_DBUGIO)
           write(3,*)' In getoldimg_mrc, 6666 ; istk,irtflg:',
     &                                          istk,irtflg
           write(3,*)' Opened old templated stackd file: ',trim(filnam)
#endif

C          WRITE OUT FILE OPENING INFO IF WANTED
           CALL LUNSAYINFO_MRC(LUN,SAYIT,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

cc         write(6,*)' In getoldimg_mrc opened old stk file: ',
cc     &               trim(filnam)

C          RETURN IMAGE NUMBER THAT WAS WANTED
           NGOT   = ISTK
   
           IRTFLG = 0
           RETURN

        ENDIF
        IRTFLG = 1

        END





C   | Templated stack:  *@stk.MRC
C   |  ` ---> FILGET_AT --> OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC
C   |                                            `  --> OPENFIL_N_MRC
C   |  Whole barestack:  @stk.MRC 
C      ` ---> FILGET_AT --> OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC
C   |                                               --> OPENFIL_N_MRC
C   |  File template:    img***.MRC
C   |   ` ---> FILGET_AT --> OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC
C   |                                                --> OPENFIL_N_MRC
C.  |. Simple file:      img001.MRC
C   | ` -------> OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC
C                                        --> OPENFIL_N_MRC
C
C
C **********************************************************************
C
C    GETNEWIMG_MRC(LUN, FILPAT,NWANT,SAYIT, FILNAM,NGOT,IRTFLG)
C
C    PURPOSE:       OPEN A SPECIFIED TEMPLATED MRC IMAGE FOR 
C                   RANDOM ACCESS READING/WRITING.
C
C    PARAMETERS:
C        LUN        I/O UNIT NUMBER FOR FILNAM.                   (SENT)
C        FILPAT     FILENAME PATTERN                              (SENT)
C        NWANT      IMAGE NUMBER WANTED                           (SENT) C
C        SAYIT      SAY FILE OPENING INFO                         (SENT)
C        FILNAM     FILENAME OPENED                               (RET.)
C        NGOT       IMAGE NUMBER FOUND                            (RET.) C
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
        INTEGER                :: ISTK,ISTK_OLD,MZ,NSTK
        INTEGER                :: LENT

        LOGICAL                :: WANTUL

        INTEGER                :: lnblnkn    ! FUNCTION

        NLET   = lnblnkn(FILPAT)

        LOCAST = INDEX(FILPAT(1:NLET),'*')
        LOCAT  = INDEX(FILPAT(1:NLET),'@')

        IF (LOCAST > 1 .or. (LOCAT  > 0 .AND. LOCAST < LOCAT)) THEN
C          SUBSTITUTE NWANT INTO FILE NAME PATTERN -> FILNAM   
           CALL FILGET_AT(FILPAT,NWANT,FILNAM,NLET,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN 

        ENDIF

C       MAY NEED TO KNOW: ITYPE,NX,NY,NZ!
C       MAY GET IT FROM PREVIOUS STATIC HEADER??

        CALL LUNGET_SIZE_MRC(LUN,NX,NY,NZ,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN 

       IF (LOCAT <= 0 .AND. LOCAST > 1) THEN
C          TEMPLATED SIMPLE MRC IMAGE ---------------------- IMG***.mrc

           CALL LUNGET_TYPE_MRC(LUN,ITYPE,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN 

           CLOSE(LUN)    ! MAY BE STILL OPEN FROM FIRST CALL  

C          OPEN NEW FILE
           ISTK = 0      ! FOR A SPECIFIC IMAGE, NOT STACKED IMG/VOL
           CALL OPFILEC(LUNCP,.FALSE.,FILNAM,LUN,'U',ITYPE,
     &                  NX,NY,NZ, 
     &                  ISTK,'UNUSED',FOUROK,IRTFLG) 
           IF (IRTFLG .NE. 0) RETURN 

           !write(3,*)' In getnewimg_mrc; new templated file: ',
           !&                 trim(filnam)

           NGOT   = NWANT   
           IRTFLG = 0
           RETURN             !------------------------------

        ELSEIF (LOCAT  > 0 .AND. LOCAST < LOCAT) THEN

C          TEMPLATED STACKED MRC FILE: *@STK.MRC ---------- *@STK.MRC
C          BARE STACKED MRC FILE: @STK.MRC ----------------- @STK.MRC
           
C          THIS STACK FILE ALREADY OPEN JUST NEED TO UPDATE FOR IMG/VOL
           ISTK = NWANT

C          SET UNDETERMINED IMAGE STATISTICS
           CALL LUNSET_STATSIMG_MRC(LUN,0,IRTFLG)

C          GET CURRENT ISTK AND NSTK FROM STATIC HEADER
           CALL LUNGET_STK_260_MRC(LUN,MZ,NSTK,ISTK_OLD,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           ISTK = NWANT
           IF (ISTK > NSTK) THEN 
C             SET CURRENT IMG/VOL AND NSTK IN STATIC HEADER
             
              NSTK = ISTK
              CALL LUNSET_STK_260_MRC(LUN,ISTK,NSTK,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

C             ADJUST NZ AND NSTK, AS NECESSARY
              CALL LUNSET_2014_MRC(LUN, NZ,NSTK, IRTFLG)

           ENDIF           

C          SET OFFSETS FOR REDLIN/WRTLIN FOR IMAGE OPENED ON THIS LUN
           CALL LUNSET_POS_MRC(LUN,ISTK,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
           !write(6,'(a,8i5)') ' In getnu_mrc22, lun,istk:',lun,istk

C          PUSH HEADER OBJECT INTO FILE TO SAVE STATS, NSTK, & ISTK
           CALL LUNWRTHED_MRC(LUN,IRTFLG)

           IF (IRTFLG .NE. 0) THEN
              LENT = lnblnkn(FILNAM)
              WRITE(NOUT,99) IRTFLG,LUN,FILNAM(:LENT)
99            FORMAT( '  *** ERROR(',I4,') ON UNIT: ',I3,' FILE: ',A)
              RETURN
           ENDIF

C          PUT COMMON VALUES INTO COMMON AREA (NOT NEEDED IN FUTURE?)
           CALL LUNSET_COMMON_MRC(LUN,IRTFLG)

C          WRITE OUT FILE OPENING INFO TO SCREEN
           CALL LUNSAYINFO_MRC(LUN,.TRUE.,IRTFLG)

           NGOT = NWANT

#if defined(SP_DBUGIO)
           write(3,*) ' '
           write(3,*)' In getnewimg_mrc; opened stacked img/vol: ',
     &                 filnam(1:nlet)
#endif

        ENDIF

        END








