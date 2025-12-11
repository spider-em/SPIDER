
C++*********************************************************************
C
C GET_FILNAM_INFO   NEW OCT 2005                 OCT 2025 ArDean Leith
C               
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
C   GET_FILNAM_INFO(FILNAM, NLIST, FIL_NOAT, FIL_DIRS,FIL_BASE,FIL_EXT,
C                   IS_MRC, IS_BARE, INUM, ISTK, IRTFLG)
C
C   PURPOSE:    FINDS STACKED IMG   NUMBER IN FILENAME
C    
C   PARAMETERS: FILNAM     FILE NAME                         !!(SENT)
C               LIST_VAL   NUMBER FROM FILELIST                (SENT)
C
C               FIL_NOAT   FILE NAME WITHOUT @ STUFF            (RET)
C               FIL_DIRS   FILE DIRECTORIES                     (RET)
C               FIL_BASE   FILE NAME                            (RET)
C               FIL_EXT    FILE EXTENSION                       (RET)
C               IS_MRC     IS A MRC/MRCS/EMD/etc FILE           (RET)
C               IS_MRCS    IS A MRCS  FILE                      (RET)
C               IS_MRC     IS A MRC/MRCS/EMD/etc FILE           (RET)
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

      INTEGER               :: LOCAT,LOCAST,LOCAST_GO,LOCAST_END
      LOGICAL               :: GOT_EXT 

      INTEGER               :: LOCSLASH,NLEN,IGO,LOCDOT,IENDB,IASTSB4

      INTEGER               :: LOCEXT,IERR,NE,IDIGS,ILENT,NAST
      INTEGER               :: NSTK

      CHARACTER(LEN=10)     :: CINUM

      CHARACTER(LEN=1)      :: NULL = CHAR(0)
 
      INTEGER               :: lnblnkn


      NLEN = lnblnkn(FILNAM)

C     FIND ENDING LOCATION OF  DIRECTORY(S) IN FILNAM_AT
      LOCSLASH = INDEX(FILNAM(1:NLEN),'/',BACK=.TRUE.)

C     STRIP OFF DIRECTORY(S)
      FIL_DIRS = ''
      FIL_BASE = FILNAM(1:NLEN)

      IF (LOCSLASH > 0) THEN
          FIL_DIRS = FILNAM(1:LOCSLASH)
          FIL_BASE = FILNAM(LOCSLASH+1:NLEN) 
      ENDIF

C     FIND LOCATION OF EXTENSION IN FILNAM
      LOCDOT   = INDEX(FILNAM(1:NLEN),'.',BACK=.TRUE.)

C     EXTRACT  EXTENSION FROM FILNAM
      FIL_EXT  = NULL
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
     &      ( INDEX(FIL_EXT,'MRC') > 0) ) THEN
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
      write(3,*)' In get_ filnam_info  -1-1-1-1'
      write(3,*)' In get_filnam_info; filnam:   ', trim(filnam)
      write(3,*)' In get_filnam_info; locslash: ', locslash
      write(3,*)' In get_filnam_info; locdot:   ', locdot
      write(3,*)' In get_filnam_info; locat:    ', locat
      write(3,*)' In get_filnam_info; locast:   ', locast
      write(3,*)' In get_filnam_info; ismrc,ismrcs:', is_mrc,is_mrcs

      write(3,*)' In get_filnam_info; fil_dirs: ', trim(fil_dirs)
      write(3,*)' In get_filnam_info; fil_base: ', trim(fil_base)
      write(3,*)' In get_filnam_info; fil_ext:  ', trim(fil_ext)
      !ilent = lnblnkn(fil_ext)
      !write(3,*)' In get_filnam_info; lnblnkn of fil_ext: ', ilent     
      write(3,*)'  '
#endif


      IF (LOCAT  == 0 .AND. 
     &    LOCAST == 0 .AND. (.NOT. IS_MRCS)) THEN      !------------

         ISTK      = 0   ! NOT A STACK

#if defined (SP_DBUGIO)
         write(3,*)' In get_filnam_info  0000;  locat: ',locat  
         write(3,*)' In get_filnam_info  0000;  istk:  ',istk      
#endif

         GOTO 999          


      ELSEIF (LOCAT == 0 .AND. IS_MRCS ) THEN ! --------------------

         !NZ      = 1
         ISTK     = 1      ! MRCS STACK WITHOUT @

#if defined (SP_DBUGIO)
         write(3,*)' In get_filnam_info  1111; locast: ',locast  
         write(3,*)' In get_filnam_info  1111; istk:   ',istk      
#endif

         GOTO 999            

      ELSEIF (LOCAT == 1 ) THEN         ! -------------------------
C        BARESTACK AT START OF FILE BASENAME

         IS_BARE   = .TRUE.
         ISTK      = 0

         FIL_BASE  = TRIM(FIL_BASE(2:))

#if defined (SP_DBUGIO)
         write(3,*)' In get_filnam_info  2222; locat: ',locat  
         write(3,*)' In get_filnam_info  2222; istk:   ',istk      
#endif

         GOTO 999   

      ELSEIF (LOCAT == IENDB ) THEN      ! -----------------------
C        BARESTACK AT END OF FILE BASENAME

#if defined (SP_DBUGIO)
         write(3,*)' In get_filnam_info    3333 '
#endif

         IS_BARE   = .TRUE.
         ISTK      = 0

         FIL_BASE  = FIL_BASE(1:IENDB-1) // NULL

#if defined (SP_DBUGIO)
         write(3,*)' In get_filnam_info  3333; locat: ',locat  
         write(3,*)' In get_filnam_info  3333; istk:   ',istk      
#endif

         GOTO 999    

      ENDIF



      IASTSB4 = VERIFY(FIL_BASE(1:LOCAT-1),'*',.TRUE.)
 
      IF ((LOCAT > 1) .AND. (IASTSB4 == 0)) THEN   ! --------------
C        ONLY '*'s  AT START OF FILE BASENAME

c         WRITE(CINUM,'(I10)',IOSTAT=IERR) LIST_VAL  
c         IF (IERR .NE. 0) THEN
c             WRITE(NOUT,*)
c     &          '*** PROBLEM WITH LIST_VAL CONVERSION: ',LIST_VAL
c             CALL ERRT(100,'GET_FILNAM_INFO',NE)
c             IRTFLG = 1
c             RETURN
c         ENDIF

         ISTK     = LIST_VAL
         FIL_BASE = FIL_BASE(LOCAT+1:) 
  
#if defined (SP_DBUGIO)
         write(3,*)'  '
         write(3,*)' In get_filnam_info  4444; locat: ',locat  
         write(3,*)' In get_filnam_info  4444; iastsb4: ',iastsb4  
         write(3,*)' In get_filnam_info  4444; istk:   ',istk      
         write(3,*)' In get_filnam_info  4444; fil_base: ',
     &                                    trim(fil_base) 
#endif

         GOTO 999
    
      ENDIF


      IASTSB4 = VERIFY(FIL_BASE(LOCAT+1:IENDB),'*')
 
      IF ((LOCAT > 1) .AND. (IASTSB4 == 0)) THEN     ! --------------
C        ONLY '*'s  AT END OF FILE BASENAME

         ISTK      = LIST_VAL
         FIL_BASE  = TRIM(FIL_BASE(1:LOCAT-1))

#if defined (SP_DBUGIO)
         write(3,*)' In get_filnam_info  5555; locat:   ',locat  
         write(3,*)' In get_filnam_info  5555; iastsb4: ',iastsb4  
         write(3,*)' In get_filnam_info  5555; istk:    ',istk      
         write(3,*)' In get_filnam_info  5555; fil_base: ',
     &                                    trim(fil_base) 
#endif

         GOTO 999  
      ENDIF


      IDIGS = VERIFY(FIL_BASE(1:LOCAT-1),'1234567890') 

      IF ((LOCAT > 1) .AND. (IDIGS == 0)) THEN     ! --------------
C        ONLY DIGITS  AT START OF FILE BASENAME BEFORE @
 
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
         write(3,*)' In get_filnam_info  6666; locat: ',locat  
         write(3,*)' In get_filnam_info  6666; idigs: ',idigs  
         write(3,*)' In get_filnam_info  6666; istk:   ',istk      
         write(3,*)' In get_filnam_info; 6666; fil_base: ',
     &                                    trim(fil_base) 
#endif

         GOTO 999    
      ENDIF


      IDIGS = VERIFY(FIL_BASE(LOCAT+1:IENDB),'1234567890') 

      IF ((LOCAT < IENDB) .AND. (IDIGS == 0)) THEN   ! --------------
C        ONLY DIGITS  AT END OF FILE BASENAME AFTER @
 
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
         write(3,*)' In get_filnam_info  7777; locat: ',locat  
         write(3,*)' In get_filnam_info  7777; idigs: ',idigs  
         write(3,*)' In get_filnam_info  7777; istk:   ',istk      
         write(3,*)' In get_filnam_info; 7777; fil_base: ',
     &                                    trim(fil_base) 
#endif

         GOTO 999    
 
      ENDIF


 999  CONTINUE !-----------------------------------------------------


      LOCAST_GO   = INDEX(FIL_BASE,'*')
      LOCAST_END  = INDEX(FIL_BASE,'*',BACK=.TRUE.)
      NAST        = LOCAST_END - LOCAST_GO

      IF (LOCAST_GO > 0) THEN                        ! --------------
 
C        REPSUB ( IN, LEFT, RIGHT, STRING, OUT ) ?????????
C        WRITE(MYSTRING,'(I10)') K
c        WRITE (ANS, '(A, I3.3)') CHR1, VAL1
C        WRITE (CINUM ,'(A, I3.3)',IOSTAT=IERR)  
C     &                        FIL_BASE(1:LOCAST_GO-1),LIST_VAL

         WRITE (CINUM ,'(I3.3)',IOSTAT=IERR) LIST_VAL

#if defined (SP_DBUGIO)
         write(3,*) ' '
         write(3,*)' In get_filnam_info; list_val: ',list_val
         write(3,*)' In get_filnam_info; cinum: ',cinum
              
#endif

         IF (IERR .NE. 0) THEN
            WRITE(NOUT,*)
     &        '*** CAN NOT WRITE NUMBER INTO CINUM STRING : ',LIST_VAL
            CALL ERRT(100,'GET_FILNAM_INFO',NE)
            IRTFLG = 1
            RETURN
         ENDIF

         FIL_BASE  = FIL_BASE(1:LOCAST_GO-1) // TRIM(CINUM) //
     &               FIL_BASE(LOCAST_END+1:) 

      ENDIF


      FIL_NOAT = ''
      IF (LEN_TRIM(FIL_DIRS) > 0) FIL_NOAT = FIL_DIRS
 
      FIL_NOAT  = TRIM(FIL_NOAT) // TRIM(FIL_BASE) 

      IF (GOT_EXT) FIL_NOAT = TRIM(FIL_NOAT) // '.' // TRIM(FIL_EXT)

#if defined (SP_DBUGIO)
         write(3,*)' In get_filnam_info; fil_dirs: ',trim(fil_dirs)
         write(3,*)' In get_filnam_info; fil_base: ',trim(fil_base) 
         write(3,*)' In get_filnam_info; fil_ext:  ',trim(fil_ext) 
         write(3,*)' In get_filnam_info; fil_noat: ',trim(fil_noat)
         write(3,*) ' '
#endif


#if defined (SP_DBUGIO_NEVER)
         write(3,*)' In get_filnam_info  8888;  locast_go:  ',
     &                                             locast_go  
         write(3,*)' In get_filnam_info  8888;  locast_end: ',
     &                                             locast_end  
         write(3,*)' In get_filnam_info  8888;  list_val: ',list_val
 
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
 
#endif


      END




















C++*********************************************************************
C
C GETFILENUM.F  -- NEW JAN 1999                   AUTHOR: ARDEAN LEITH
C                  EXTRACTED FROM LUNSETHDR AUG 02 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors: Joachim Frank & ArDean Leith & Tapu Shaikh  *
C=* Copyright 1985-2010  Health Research Inc.,                         *
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
C    GETFILENUM(FILNAM,IMGNUM,NDIGITS,CALLERRT,IRTFLG)  
C
C    PURPOSE:    FINDS FILE NUMBER AT END OF FILENAME
C    
C    PARAMETERS: FILNAM    CHAR. VARIABLE FILE NAME             (SENT)
C                IMGNUM    NUMBER IN FILE NAME                   (RET)
C                NDIGITS   NUMBER OF DIGITS                      (RET)
C                CALLERRT  CALL ERRT IF ERROR                   (SENT)
C                IRTFLG    ERROR FLAG                            (RET)
C
C    CALLED BY:  deletf.f   filgen.f    openinstk.f 
C                openstk.f  to_peaks.f 
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************
 
       SUBROUTINE GETFILENUM(FILNAM,IMGNUM,NDIGITS,CALLERRT,IRTFLG)

       IMPLICIT  NONE

       INCLUDE 'CMBLOCK.INC'

       CHARACTER *(*) :: FILNAM
       LOGICAL        :: CALLERRT
       INTEGER        :: IMGNUM,NDIGITS,IRTFLG

       CHARACTER *1   :: CHARI
       INTEGER        :: NLET,IGO,I,NE
       INTEGER        :: lnblnkn


C      FIND NUMBER OF CHAR. IN FILNAM
       NLET   = LNBLNKN(FILNAM)

       IGO    = NLET - 1
       IMGNUM = 0

C      EXTRACT IMGNUM FROM FILENAME
       DO I = NLET,1,-1
          CHARI = FILNAM(I:I)
          IF (CHARI .LT. '0' .OR. CHARI .GT. '9') EXIT
          IGO = I 
       ENDDO

       NDIGITS = NLET - IGO + 1
       IF (NDIGITS .LE. 0 .OR. NDIGITS .GT. 10) THEN
C         NO NUMBER AT END OF FILNAM OR > 10 DIGITS
          IRTFLG = -1
          RETURN
       ENDIF

       READ(FILNAM(IGO:NLET),'(I10)',ERR=999) IMGNUM
       IRTFLG = 0
       RETURN
 
    
999    WRITE(NOUT,*) '*** CAN NOT GET FILE NUMBER FROM: ',FILNAM(:NLET)
       IF (CALLERRT) THEN
          CALL ERRT(100,'GETFILENUM',NE)
       ENDIF
       IRTFLG = 1

       RETURN
       END


#ifdef NEVER
      PURE LOGICAL FUNCTION LITTLE_ENDIAN()
      
      INTEGER(INT8)  :: J(2)
      INTEGER(INT16) :: I
      
      I = 1
      J = TRANSFER(SOURCE=I,MOLD=J,SIZE=2)
      
      IF (J(1) == 1) THEN
         LITTLE_ENDIAN = .TRUE.
      ELSE
         LITTLE_ENDIAN = .FALSE.
      ENDIF
      
      END FUNCTION
      
      
      SUBROUTINE ENDIAN(LITEND)

C     CHECKS IF THIS IS A LITTLE ENDIAN MACHINE
C     RETURNS LITEND =.TRUE. IF IT IS, LITEND =.FALSE. IF NOT

      INTEGER*1 J(2)
      INTEGER*2 I
      EQUIVALENCE (I,J)
      
      LOGICAL LITEND

      I = 1
      IF (J(1) .EQ. 1) THEN
         LITEND = .TRUE.
      ELSE
         LITEND = .FALSE.
      ENDIF

      END
      
      PROGRAM TESTI
      ! Posted by Perseus in comp.lang.fortran on 4 July 2005.
      ! and Paul Van Delst and David Flower on 5 July 2005.

      LOGICAL, PARAMETER :: BIGEND = IACHAR(TRANSFER(1,"A")) == 0

      IF (BIGEND) THEN
        PRINT *, "BIG ENDIAN"
      ELSE
        PRINT *, "LITTLE ENDIAN"
      ENDIF

END PROGRAM TESTI

PROGRAM DETERMINE_ENDIANNESS
    IMPLICIT NONE
    INTEGER :: UNIT, IOSTAT
    INTEGER(4) :: VALUE
    CHARACTER(LEN=4) :: BYTES

    ! OPEN THE FILE IN UNFORMATTED MODE
    OPEN(UNIT=10, FILE="FILE.DAT", FORM="UNFORMATTED", ACCESS="STREAM", STATUS="OLD", IOSTAT=IOSTAT)
    IF (IOSTAT /= 0) THEN
        PRINT *, "ERROR OPENING FILE."
        STOP
    END IF

    ! READ THE FIRST 4 BYTES (ASSUMING THE FILE STARTS WITH A KNOWN INTEGER)
    READ(10) VALUE
    CLOSE(10)

    ! CONVERT THE INTEGER TO BYTES FOR ANALYSIS
    CALL MOVE_ALLOC(TRANSFER(VALUE, BYTES), BYTES)

    ! CHECK THE BYTE ORDER
    IF (BYTES == 'EXPECTED_BIG_ENDIAN_BYTES') THEN
        PRINT *, "FILE IS BIG ENDIAN."
    ELSE IF (BYTES == 'EXPECTED_LITTLE_ENDIAN_BYTES') THEN
        PRINT *, "FILE IS LITTLE ENDIAN."
    ELSE
        PRINT *, "UNABLE TO DETERMINE ENDIANNESS."
    END IF
END PROGRAM DETERMINE_ENDIANNESS



#endif


