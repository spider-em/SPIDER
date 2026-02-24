
C++*********************************************************************
C
C GET_FILNAM_INFO   NEW OCT 2005                 Oct 2025 ArDean Leith
C                   ALTERED                      Dec 2025 ArDean Leith
          
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors: Joachim Frank & ArDean Leith & Tapu Shaikh  *
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
C ************************* GET_FILNAM_INFO ****************************
C
C   GET_FILNAM_INFO(FILNAM, LIST_VAL,
C                   FIL_NOAT,FIL_DIRS,FIL_BASE,FIL_EXT,
C                   IS_MRC, IS_MRCS, IS_BARE, 
C                   ISTK,   IRTFLG)
C
C   PURPOSE:    FINDS STACKED IMG NUMBER IN FILENAME
C               AlSO SUBSTITUTES FOR *** IN FILENAME
C    
C   PARAMETERS: FILNAM     FILE NAME                           (SENT)
C               LIST_VAL   NUMBER FROM FILELIST                (SENT)
C                          MAY CONTAIN EITHER ISTK OR IMGNUM
C
C               FIL_NOAT   FILE NAME WITHOUT @ STUFF            (RET)
C               FIL_DIRS   FILE DIRECTORY                       (RET)
C               FIL_BASE   BASE FILE NAME                       (RET)
C               FIL_EXT    FILE EXTENSION                       (RET)
C               IS_MRC     IS A MRC/MRCS/EMD/etc FILE           (RET)
C               IS_MRCS    IS A MRCS  FILE                      (RET)
C               IS_BARE    IS A BARESTACK                       (RET)
C 
C               ISTK       STACK IMG/VOL                        (RET)
C                             0 : NOT STACK or BARE STACK
C                            >1 : STACKED IMAGE/VOL NUMBER       
C
C               IRTFLG     ERROR FLAG                           (RET)
C  
C     IF PRESENT, EXTRACT STACK FILENAME AND NUMBER e.g. from
C     056@file.mrcs       or file3321@100.mrcs  or /dir/fileaa@21 or
C     /dir/fileaa@21.mrcs or /dir/056@file         etc 
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE GET_FILNAM_INFO(FILNAM, LIST_VAL,
     &                      FIL_NOAT,FIL_DIRS,FIL_BASE,FIL_EXT,
     &                      IS_MRC, IS_MRCS, IS_BARE, 
     &                      ISTK,   IRTFLG)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=*)      :: FILNAM
      INTEGER               :: LIST_VAL

      CHARACTER(LEN=*)      :: FIL_NOAT,FIL_DIRS,FIL_BASE,FIL_EXT
      LOGICAL               :: IS_MRC,IS_MRCS, IS_BARE
      INTEGER               :: ISTK, IRTFLG

      CHARACTER(LEN=MAXNAM) :: FIL_NEW
      INTEGER               :: LOCAT,LOCAST,LOCAST_GO,LOCAST_END
      LOGICAL               :: GOT_EXT 

      INTEGER               :: LOCSLASH,NLEN,IGO,LOCDOT,IENDB,IASTSB4

      INTEGER               :: LOCEXT,IERR,NE,IDIGS,ILENT,NAST,NLET

      CHARACTER(LEN=10)     :: CINUM

      INTEGER               :: lnblnkn

      character(len=1)      :: c22,c23

      NLEN = lnblnkn(FILNAM)

C     FIND ENDING LOCATION OF  DIRECTORY(S) IN FILNAM_AT
      LOCSLASH = INDEX(FILNAM(1:NLEN),'/',BACK=.TRUE.)

C     STRIP OFF DIRECTORY(S)
      FIL_DIRS = ''
      FIL_EXT  = ''
      FIL_BASE = FILNAM(1:NLEN)

      IF (LOCSLASH > 0) THEN
          FIL_DIRS = FILNAM(1:LOCSLASH)
          FIL_BASE = FILNAM(LOCSLASH+1:NLEN) 
      ENDIF

C     FIND LOCATION OF EXTENSION IN FILNAM
      LOCDOT   = INDEX(FILNAM(1:NLEN),'.',BACK=.TRUE.)

C     EXTRACT  EXTENSION FROM FILNAM
      GOT_EXT  = (LOCDOT > 0)

      IF (GOT_EXT) THEN
C        GET FILE EXTENSION
         FIL_EXT  = TRIM(FILNAM(LOCDOT+1:))

C        STRIP OFF EXTENSION FROM FIL_BASE
         LOCEXT   = INDEX(FIL_BASE,'.',BACK=.TRUE.)
         FIL_BASE = FIL_BASE(1:LOCEXT-1)
      ENDIF
               
      IS_MRC = .FALSE.
      IF  ( ( INDEX(FIL_EXT,'mrc') > 0) .OR. 
     &      ( INDEX(FIL_EXT,'MRC') > 0) .OR. 
     &      ( INDEX(FIL_EXT,'map') > 0) .OR. 
     &      ( INDEX(FIL_EXT,'MAP') > 0) ) THEN
        IS_MRC = .TRUE.
      ENDIF

      IS_MRCS = .FALSE.
      IF  ( ( INDEX(FIL_EXT,'mrcs') > 0) .OR. 
     &      ( INDEX(FIL_EXT,'MRCS') > 0) ) THEN
        IS_MRCS = .TRUE.
      ENDIF

      IENDB   = lnblnkn(FIL_BASE)
      LOCAT   = INDEX(FIL_BASE,'@')
      LOCAST  = INDEX(FIL_BASE,'*')
      IS_BARE = .FALSE.


#if defined (SP_DBUGIO)
      write(3,*)'  '
      write(3,*)' In get_filnam_info  -1-1-1'
      write(3,*)' In get_filnam_info; filnam:   |', trim(filnam),'||'
      write(3,*)' In get_filnam_info; locslash:  ', locslash
      write(3,*)' In get_filnam_info; locdot:    ', locdot
      write(3,*)' In get_filnam_info; locat:     ', locat
      write(3,*)' In get_filnam_info; locast:    ', locast
      write(3,*)' In get_filnam_info; ismrc,ismrcs: ', is_mrc,is_mrcs

      write(3,*)' In get_filnam_info; list_val:  ', list_val
 
      write(3,*)' In get_filnam_info; fil_dirs: |', trim(fil_dirs),'||'
      write(3,*)' In get_filnam_info; fil_base: |', trim(fil_base),'||'
      write(3,*)' In get_filnam_info; fil_ext:  |', trim(fil_ext),'||'
#endif


      IF (LOCAT  == 0 .AND. 
     &    LOCAST == 0 .AND. (.NOT. IS_MRCS)) THEN      !------------

         ISTK      = 0   ! NOT A STACK

#if defined (SP_DBUGIO)
         write(3,*)' In get_filnam_info 0000; istk:  ',istk      
#endif

         GOTO 999          


      ELSEIF (LOCAT == 0 .AND. IS_MRCS ) THEN ! --------------------

         ISTK     = 1      ! MRCS STACK WITHOUT @

#if defined (SP_DBUGIO)
         write(3,*)' In get_filnam_info 1111; istk:   ',istk    
#endif

         GOTO 999            

      ELSEIF (LOCAT == 1 ) THEN         ! -------------------------
C        BARESTACK AT START OF FILE BASENAME

         IS_BARE   = .TRUE.
         ISTK      = 1
         FIL_BASE  = TRIM(FIL_BASE(2:)) 

#if defined (SP_DBUGIO)
      write(3,*)' In get_filnam_info 2222; is_bare,istk: ',is_bare,istk
#endif

         GOTO 999   

      ELSEIF (LOCAT == IENDB ) THEN      ! -----------------------
C        BARESTACK AT END OF FILE BASENAME

         IS_BARE   = .TRUE.
         ISTK      = 1
         FIL_BASE  = TRIM(FIL_BASE(1:IENDB-1)) 

#if defined (SP_DBUGIO)
      write(3,*)' In get_filnam_info 3333; is_bare,istk: ',is_bare,istk      
#endif
         GOTO 999    

      ENDIF



      IASTSB4 = VERIFY(FIL_BASE(1:LOCAT-1),'*',.TRUE.)
 
      IF ((LOCAT > 1) .AND. (IASTSB4 == 0)) THEN   ! --------------
C        ONLY '*'s BEFORE @ AT START OF FILE BASENAME e.g. (*@..)

         ISTK     = LIST_VAL
         FIL_BASE = TRIM(FIL_BASE(LOCAT+1:)) 
  
#if defined (SP_DBUGIO)
         write(3,*)'  '
         write(3,*)' In get_filnam_info 4444; iastsb4: ',iastsb4  
         write(3,*)' In get_filnam_info 4444; istk:    ',istk      
         write(3,*)' In get_filnam_info 4444; fil_base: |',
     &                                      trim(fil_base),'|' 
#endif

         GOTO 999
    
      ENDIF



      IASTSB4 = VERIFY(FIL_BASE(LOCAT+1:IENDB),'*')
 
      IF ((LOCAT > 1) .AND. (IASTSB4 == 0)) THEN     ! --------------
C        ONLY '*'s  AT END OF FILE BASENAME e.g.(..@*)

         ISTK     = LIST_VAL
         FIL_BASE = TRIM(FIL_BASE(1:LOCAT-1))

#if defined (SP_DBUGIO)
         write(3,*)' In get_filnam_info 5555; iastsb4: ',iastsb4  
         write(3,*)' In get_filnam_info 5555; istk:    ',istk      
         write(3,*)' In get_filnam_info 5555; fil_base: |',
     &                                      trim(fil_base),'|' 
#endif

         GOTO 999  
      ENDIF


      IDIGS = VERIFY(FIL_BASE(1:LOCAT-1),'1234567890') 

      IF ((LOCAT > 1) .AND. (IDIGS == 0)) THEN     ! --------------
C        ONLY DIGITS  AT START OF FILE BASENAME BEFORE @ e.g.(123@..)
 
         READ(FIL_BASE(1:LOCAT-1),'(I10)',IOSTAT=IERR) ISTK
         IF (IERR .NE. 0) THEN
             WRITE(NOUT,*)
     &          '*** CAN NOT GET ISTK NUMBER FROM: ',FIL_BASE(1:LOCAT)
             CALL ERRT(100,'GET_FILNAM_NSTK',NE)
             IRTFLG = 1
             RETURN
         ENDIF

         FIL_BASE  = TRIM(FIL_BASE(LOCAT+1:))

#if defined (SP_DBUGIO)
         write(3,*)' In get_filnam_info 6666; idigs:    ',idigs  
         write(3,*)' In get_filnam_info 6666; istk:     ',istk      
         write(3,*)' In get_filnam_info 6666; fil_base: |',
     &                                      trim(fil_base),'|' 
#endif

         GOTO 999    
      ENDIF


      IDIGS = VERIFY(FIL_BASE(LOCAT+1:IENDB),'1234567890') 

      IF ((LOCAT < IENDB) .AND. (IDIGS == 0)) THEN   ! --------------
C        ONLY DIGITS  AT END OF FILE BASENAME AFTER @ e.g.(...@123)
 
        READ(FIL_BASE(LOCAT+1:IENDB),'(I10)',IOSTAT=IERR) ISTK

        IF (IERR .NE. 0) THEN
            WRITE(NOUT,*)
     &        '*** CAN NOT GET ISTK NUMBER FROM: ',FIL_BASE(LOCAT:IENDB)
            CALL ERRT(100,'get_filnam_info',NE)
            IRTFLG = 1
            RETURN
        ENDIF

        FIL_BASE  = TRIM(FIL_BASE(1:LOCAT-1))
    
#if defined (SP_DBUGIO)
        write(3,*)' In get_filnam_info 7777; idigs:    ',idigs  
        write(3,*)' In get_filnam_info 7777; istk:     ',istk      
        write(3,*)' In get_filnam_info 7777; fil_base: |',
     &                                      trim(fil_base),'|' 
#endif

         GOTO 999    
      ENDIF



      LOCAST_GO   = INDEX(FIL_BASE,'*')
      LOCAST_END  = INDEX(FIL_BASE,'*',BACK=.TRUE.)
      NAST        = LOCAST_END - LOCAST_GO

      IF (LOCAST_GO > 0) THEN                        ! --------------
 
C        SUBSTITUTE LIST_VAL (IMGNUM) INTO FILPAT 
         NLET = lnblnkn(FIL_BASE)
         CALL FILGET(FIL_BASE,FIL_NEW,NLET,LIST_VAL,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         ISTK     = 0
         FIL_BASE = TRIM(FIL_NEW)
  
#if defined (SP_DBUGIO)
        write(3,*)' In get_filnam_info 8888; list_val:',list_val
        write(3,*)' In get_filnam_info 8888; fil_base: |',
     &                                       trim(fil_base),'|' 
        write(3,*)' In get_filnam_info 8888; fil_new: |',
     &                                      trim(fil_new),'|'
#endif
         GOTO 999

      ENDIF

C        REPSUB ( IN, LEFT, RIGHT, STRING, OUT ) ?????????
C        WRITE(MYSTRING,'(I10)') K
c        WRITE (ANS, '(A, I3.3)') CHR1, VAL1
C        WRITE (CINUM ,'(A, I3.3)',IOSTAT=IERR)  
C     &                        FIL_BASE(1:LOCAST_GO-1),LIST_VAL
C        WRITE (CINUM ,'(I3.3)',IOSTAT=IERR) LIST_VAL
C        FIL_BASE  = FIL_BASE(1:LOCAST_GO-1) // TRIM(CINUM) //
C     &               FIL_BASE(LOCAST_END+1:) 
C        IF (IERR .NE. 0) THEN
C           WRITE(NOUT,*)
C    &        '*** CAN NOT WRITE NUMBER INTO CINUM STRING : ',LIST_VAL
C           CALL ERRT(100,'GET_FILNAM_INFO',NE)
C           IRTFLG = 1
C           RETURN
C        ENDIF


999   CONTINUE !-----------------------------------------------------


      FIL_NOAT = ''
      IF (LEN_TRIM(FIL_DIRS) > 0) FIL_NOAT = FIL_DIRS

 
      FIL_NOAT  = TRIM(FIL_NOAT) // TRIM(FIL_BASE) 

      IF (GOT_EXT) FIL_NOAT = TRIM(FIL_NOAT) // '.' // TRIM(FIL_EXT)

#if defined (SP_DBUGIO)
      write(3,*)' '
      write(3,*)' In get_filnam_info 9999; fil_dirs: ',trim(fil_dirs)
      write(3,*)' In get_filnam_info 9999; fil_base: ',trim(fil_base) 
      write(3,*)' In get_filnam_info 9999; fil_ext:  ',trim(fil_ext) 
      write(3,*)' In get_filnam_info 9999; fil_noat: ',trim(fil_noat)

      c22 = fil_noat(22:22)
      c23 = fil_noat(23:23)
      write(3,*)' In get_filnam_info 9999; c22: |',c22,'|'
      write(3,*)' In get_filnam_info 9999; c23: |',c23,'|'
#endif



      END


#if defined (SP_DBUGIO_NEVER)
      write(3,*)' In get_filnam_info 0000; locast_go:  ', locast_go  
      write(3,*)' In get_filnam_info 0000; locast_end: ', locast_end  
      write(3,*)' In get_filnam_info 0000; list_val:   ', list_val
 
      write(3,*) ' '
      write(3,*)' In get_filnam_info; tr fil_dirs: ',trim(fil_dirs) 
      ilent = lnblnkn(fil_dirs)
      write(3,*)' In get_filnam_info; lnblnkn   fil_dirs: ',ilent
      ilent = len_trim(fil_dirs)
      write(3,*)' In get_filnam_info; len_trim  fil_dirs: ',ilent

      write(3,*) ' '
      write(3,*)' In get_filnam_info; tr fil_base: ',trim(fil_base) 
      ilent = lnblnkn(fil_base)
      write(3,*)' In get_filnam_info; lnblnkn  fil_base: ',ilent
      ilent = len_trim(fil_base)
      write(3,*)' In get_filnam_info; len_trim fil_base: ',ilent

      write(3,*) ' '
      write(3,*)' In get_filnam_info; tr fil_ext: ',trim(fil_ext) 
      ilent = lnblnkn(fil_ext)
      write(3,*)' In get_filnam_info; lnblnkn  fil_ext:  ',ilent
      ilent = len_trim(fil_ext)
      write(3,*)' In get_filnam_info; lentrim  fil_ext:  ',ilent

      write(3,*) ' '
      write(3,*)' In get_filnam_info; tr fil_noat: ',trim(fil_noat) 
      ilent = lnblnkn(fil_noat)
      write(3,*)' In get_filnam_info; lnblnkn fil_noat: ',ilent
      ilent = len_trim(fil_noat)

      END
#endif












