
C++*********************************************************************
C
C GETFILENUM.F  -- NEW JAN 1999                   AUTHOR: ARDEAN LEITH
C                  EXTRACTED FROM LUNSETHDR AUG 02 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
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
C    CALLED BY:      deletf.f    filgen.f     openinstk.f 
C                    openstk.f   to_peaks.f 
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


C++*********************************************************************
C
C GET_FILNAM_NSTK   NEW OCT 2005                 OCT 2025 ArDean Leith
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
C ************************* GET_FILNAM_NSTK ****************************
C
C   GET_FILNAM_NSTK(FILNAM,FIL_NOAT, FIL_DIRS,FIL_BASE,FIL_EXT,
C                   IS_MRC,LOCAST, IS_STK,NSTK_FLG, IRTFLG)
C
C   PURPOSE:    FINDS STACKED IMG  NUMBER IN FILENAME
C    
C   PARAMETERS: FILNAM     CHAR. VARIABLE FILE NAME            (SENT)
C               FIL_NOAT   CHAR. VARIABLE FILE NAME             (RET)
C               FIL_DIRS   CHAR. VARIABLE FILE NAME             (RET)
C               FIL_BASE   CHAR. VARIABLE FILE NAME             (RET)
C               FIL_EXT    CHAR. VARIABLE FILE NAME             (RET)
C               IS_MRC     NUMBER IN FILE NAME                  (RET)
C               LOCAST     LOCATION OF * IN FILE NAME           (RET)
C               IS_STK     IS A STACK                   (SENT or RET)
C               NSTK_FLG   STACK IMG NUMBER FROM FILE NAME      (RET)
C                            -2 :  NOT A STACK   (NO @)
C                            -1 :  AST IN STACK   
C                             0 :  BARE STACK
C                            >1 :  STACKED IMAGE NUMBER
C
C               IRTFLG     ERROR FLAG                       (RET.)
C  
C     IF PRESENT, EXTRACT STACK FILENAME AND NUMBER e.g. from
C     056@file.mrcs       or file3321@100.mrcs  or /dir/fileaa@21 or
C     /dir/fileaa@21.mrcs or /dir/056@file         etc 
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE GET_FILNAM_NSTK(FILNAM,
     &                      FIL_NOAT,FIL_DIRS,FIL_BASE,FIL_EXT,
     &                      IS_MRC, LOCAST, IS_STK, NSTK_FLG,IRTFLG)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=*)              :: FILNAM,FIL_DIRS,FIL_BASE,FIL_EXT
      CHARACTER(LEN=*),INTENT(OUT)  :: FIL_NOAT
      LOGICAL                       :: IS_MRC,IS_STK
      INTEGER                       :: LOCAST,NSTK_FLG,IRTFLG

      LOGICAL               :: GOT_EXT, IS_MRCS
      INTEGER               :: LOCAT

      INTEGER               :: LOCSLASH,NLEN,IGO,LOCDOT,IENDB,IASTSB4
      INTEGER               :: LOCEXT,IERR,NE,IDIGS,ILENT
      CHARACTER(LEN=1)      :: NULL = CHAR(0)
 
      INTEGER               :: lnblnkn


      NLEN = lnblnkn(FILNAM)

C     FIND ENDING LOCATION OF  DIRECTORY(S) IN FILNAM_AT
      LOCSLASH = INDEX(FILNAM(1:NLEN),'/',BACK=.TRUE.)

C     STRIP OFF DIRECTORY(S)
      FIL_DIRS = NULL
      FIL_BASE = NULL

      IF (LOCSLASH > 0) THEN
          FIL_DIRS = FILNAM(1:LOCSLASH)
          FIL_BASE = FILNAM(LOCSLASH+1:NLEN) 
      ENDIF

#if defined (SP_DBUGIO)
      write(3,*)' In get_filnam_nstk  0000000'
      write(3,*)' In get_filnam_nstk; locslash:',locslash

      ilent = lnblnkn(filnam)
      write(3,*)' In get_filnam_nstk; filnam: ',   filnam(1:ilent)
      ilent = lnblnkn(fil_dirs)
      write(3,*)' In get_filnam_nstk; fil_dirs: ', fil_dirs(1:ilent)
#endif


C     FIND LOCATION OF EXTENSION IN FILNAM
      LOCDOT   = INDEX(FILNAM(1:NLEN),'.',BACK=.TRUE.)

C     EXTRACT  EXTENSION FROM FILNAM
      GOT_EXT  = (LOCDOT > 0)

      FIL_EXT  = NULL
 
      IF (GOT_EXT) THEN
C        GET FILE EXTENSION
         FIL_EXT  = FILNAM(LOCDOT+1:NLEN)

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
        IS_MRC = .TRUE.
      ENDIF


#if defined (SP_DBUGIO)
      write(3,*)' In get_filnam_nstk; ismrc,ismrcs:', is_mrc,is_mrcs
 
      ilent = lnblnkn(fil_ext)
      write(3,*)' In get_filnam_nstk; fil_ext: ',  fil_ext(1:ilent) 
      ilent = lnblnkn(fil_base)
      write(3,*)' In get_filnam_nstk; fil_base: ', fil_base(1:ilent) 
#endif

      IENDB   = lnblnkn(FIL_BASE)
      LOCAT   = INDEX(FIL_BASE(1:IENDB),'@')

#if defined (SP_DBUGIO)
      write(3,*)' In get_filnam_nstk; iendb,nlen:',iendb,nlen
      write(3,*)' In get_filnam_nstk; locslash,locdot,locat:',
     &                                locslash,locdot,locat 
        
#endif
      

      IRTFLG   = 0

      IF (LOCAT <= 0 .AND. (.NOT. IS_MRCS)) THEN

#if defined (SP_DBUGIO)
         write(3,*)' In get_filnam_nstk  0000000 '      
#endif

         IS_STK    = .FALSE.
         FIL_NOAT  = FILNAM
         NSTK_FLG  = -2   ! NOT A STACK
         GOTO 999         !RETURN

      ELSEIF (LOCAT <= 0 .AND. IS_MRCS ) THEN

#if defined (SP_DBUGIO)
         write(3,*)' In get_filnam_nstk  1111111 ' 
#endif

         IS_STK    = .TRUE.
         FIL_NOAT  = FILNAM
         NSTK_FLG  = 1      ! MRCS STACK WITHOUT @
         GOTO 999           ! RETURN

      ELSEIF (LOCAT == 1 ) THEN
C        BARESTACK AT START OF FILE BASENAME

#if defined (SP_DBUGIO)
         write(3,*)'  In GET_FILNAM_NSTK  2222222 '
#endif
 
         IS_STK    = .TRUE.
         NSTK_FLG  = 0      ! BARESTACK
         FIL_BASE  = FIL_BASE(2:IENDB)

         ILENT     = lnblnkn(FIL_NOAT)
         IF (GOT_EXT) FIL_NOAT = FIL_NOAT(1:ILENT) // '.' // FIL_EXT
         GOTO 999   

      ELSEIF (LOCAT == IENDB ) THEN
C        BARESTACK AT END OF FILE BASENAME

#if defined (SP_DBUGIO)
         write(3,*)' In get_filnam_nstk    3333333 '
#endif

         IS_STK    = .TRUE.
         NSTK_FLG  = 0      ! BARESTACK
         FIL_BASE  = FIL_BASE(1:IENDB-1) // NULL

         ILENT     = lnblnkn(FIL_DIRS)
         FIL_NOAT  = FIL_DIRS(1:ILENT) // FIL_BASE 
         ILENT     = lnblnkn(FIL_NOAT)
         IF (GOT_EXT) FIL_NOAT = FIL_NOAT(1:ILENT) // '.' // FIL_EXT
         GOTO 999    

      ENDIF

#if defined (SP_DBUGIO)
      !write(3,*)' In get_filnam_nstk  3333333  to   4444444 '
#endif
   
      IASTSB4 = VERIFY(FIL_BASE(1:LOCAT-1),'*',.TRUE.) 
      IF ((LOCAT > 1) .AND. (IASTSB4 == 0)) THEN
C        ONLY '*'s  AT START OF FILE BASENAME

#if defined (SP_DBUGIO)
         write(3,*)' In get_filnam_nstk 4444444 '
#endif

         IS_STK    = .TRUE.
         NSTK_FLG  = -1      ! (*)@

         ILENT     = lnblnkn(FIL_DIRS)
         FIL_NOAT  = FIL_DIRS(1:ILENT) // FIL_BASE(LOCAT+1:IENDB) 

         ILENT     = lnblnkn(FIL_NOAT)
         IF (GOT_EXT) FIL_NOAT = FIL_NOAT(1:ILENT) // '.' // FIL_EXT
         GOTO 999    
      ENDIF

#if defined (SP_DBUGIO)
      !write(3,*)' In get_filnam_nstk    4444444 to  5555555 '
#endif

      IASTSB4 = VERIFY(FIL_BASE(LOCAT+1:IENDB),'*') 
      IF ((LOCAT > 1) .AND. (IASTSB4 == 0)) THEN
C        ONLY '*'s  AT END OF FILE BASENAME

#if defined (SP_DBUGIO)
         !write(3,*)' In get_filnam_nstk   5555555 '
         !write(3,*)' In get_filnam_nstk; ilent: ',  ilent
         !write(3,*)' In get_filnam_nstk; iastsb4,locat: ',    iastsb4,locat
         !write(3,*)' In get_filnam_nstk; fil_base ',fil_base(1:locat)
#endif

         IS_STK    = .TRUE.
         NSTK_FLG  = -1       ! @(*)

         !!FIL_BASE  = FIL_BASE(1:LOCAT-1) // NULL

         ILENT     = lnblnkn(FIL_DIRS)
         FIL_NOAT  = FIL_DIRS(1:ILENT) // FIL_BASE 
         ILENT     = lnblnkn(FIL_NOAT)
         IF (GOT_EXT) FIL_NOAT = FIL_NOAT(1:ILENT) // '.' // FIL_EXT

         GOTO 999  
      ENDIF


#if defined (SP_DBUGIO)
      !write(3,*)' In get_filnam_nstk    5555555   TO  6666666 '
#endif

      IDIGS = VERIFY(FIL_BASE(1:LOCAT-1),'1234567890') 

      IF ((LOCAT > 1) .AND. (IDIGS == 0)) THEN
C        ONLY DIGITS  AT START OF FILE BASENAME BEFORE @
 
#if defined (SP_DBUGIO)
        write(3,*)' In get_filnam_nstk   6666666 '
#endif

         READ(FIL_BASE(1:LOCAT-1),'(I10)',IOSTAT=IERR) NSTK_FLG
         IF (IERR .NE. 0) THEN
             WRITE(NOUT,*)
     &          '*** CAN NOT GET STK NUMBER FROM: ',FIL_BASE(1:LOCAT)
             CALL ERRT(100,'GET_FILNAM_NSTK',NE)
             IRTFLG = 1
             RETURN
         ENDIF

         IS_STK    = .TRUE.

         ILENT     = lnblnkn(FIL_DIRS)
         FIL_NOAT  = FIL_DIRS(1:ILENT) // FIL_BASE(LOCAT+1:)
 
         ILENT     = lnblnkn(FIL_NOAT)
         IF (GOT_EXT) FIL_NOAT = FIL_NOAT(1:ILENT) // '.' // FIL_EXT

         GOTO 999    !RETURN
      ENDIF

#if defined (SP_DBUGIO)
      !write(3,*)' In get_filnam_nstk    6666666 to  7777777  '
#endif


      IDIGS = VERIFY(FIL_BASE(LOCAT+1:IENDB),'1234567890') 
      IF ((LOCAT < IENDB) .AND. (IDIGS == 0)) THEN
C        ONLY DIGITS  AT END OF FILE BASENAME AFTER @
 
#if defined (SP_DBUGIO)
        write(3,*)' In get_filnam_nstk    7777777 '
        write(3,*)' In get_filnam_nstk; locat,idigs,iendb:',
     &                                  locat,idigs,iendb
        ilent = lnblnkn(fil_base)
        write(3,*)' In get_filnam_nstk; fil_base ', 
     &                                  fil_base(locat+1:ilent)
c
c        write(3,*)' In get_filnam_nstk; fil_base(:locat-1):',
c     &                                  fil_base(:locat-1)
c
c        write(3,*)' In get_filnam_nstk; fil_base(locat+1:iendb) ',
c     &                                  fil_base(locat+1:iendb)
c
c        ilent = lnblnkn(fil_base)
c        write(3,*)' In get_filnam_nstk; fil_base ', 
c     &                                  fil_base(locat+1:ilent)

#endif

        READ(FIL_BASE(LOCAT+1:IENDB),'(I10)',IOSTAT=IERR) NSTK_FLG

        IF (IERR .NE. 0) THEN
            WRITE(NOUT,*)
     &        '*** CAN NOT GET STK NUMBER FROM: ',FIL_BASE(LOCAT:IENDB)
            CALL ERRT(100,'get_filnam_nstk',NE)
            IRTFLG = 1
            RETURN
        ENDIF

        IS_STK    = .TRUE.

        ILENT     = lnblnkn(FIL_DIRS)
        FIL_NOAT  = FIL_DIRS(1:ILENT) // FIL_BASE(1:LOCAT-1)
 
        ILENT     = lnblnkn(FIL_NOAT)
        IF (GOT_EXT) FIL_NOAT = FIL_NOAT(1:ILENT) // '.' // FIL_EXT
    
        GOTO 999    
 
      ENDIF

      write(3,*) '  Should not get to this line!!!!!! '

 999  CONTINUE

      LOCAST    = INDEX(FIL_NOAT,'*')

#if defined (SP_DBUGIO)

      write(3,*)'  '
      !ilent = lnblnkn(fil_noat)
      !write(3,*)' In get_filnam_nstk; lnblnkn fil_noat: ', ilent
      !ilent = len(fil_noat)
      !write(3,*)' In get_filnam_nstk; ilent fil_noat: ', ilent
      !write(3,*)' End get_filnam_nstk; fil_noat: ',    fil_noat(:ilent)

      ilent = lnblnkn(fil_noat)
      write(3,*)' End get_filnam_nstk; fil_noat: ', fil_noat(1:ilent)

      ilent = lnblnkn(fil_dirs)
      write(3,*)' End get_filnam_nstk; fil_dirs: ', fil_dirs(1:ilent)

      ilent = lnblnkn(fil_base)
      write(3,*)' End get_filnam_nstk; fil_base: ', fil_base(1:ilent)

      ilent = lnblnkn(fil_ext)
      write(3,*)' End get_filnam_nstk; fil_ext:  ', fil_ext(1:ilent)

      write(3,*)' End get_filnam_nstk; locast:   ', locast

      write(3,*)' End get_filnam_nstk; nstk_flg,irtflg: ',
     &                                 nstk_flg,irtflg
   
      write(3,*)' '



        
#endif

      END



