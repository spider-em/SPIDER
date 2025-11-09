C++********************************************************************
C
C LUNSETMRCHDR.F  DERIVED FROM LUNSETHDR.F        JUN 2019 ArDean Leith
C                 NATIVE MRC FILE SUPPORT         JUL 2019 ArDean Leit
C                 NVIDIA COMPILER CHANGES         SEP 2025 Tapu Shaikh
C                 NVIDIA COMPILER CHANGES UNDONE  SEP 2025 ArDean Leith
C                 DEBUG OUTPUT CHANGED,CTIT_LOCAL SEP 2025 ArDean Leith
C                 NBYT, NBYT_PER_VAL CONFUSED     SEP 2025 ArDean Leith
C                 DEBUG OUTPUT ADDED              OCT 2025 ArDean Leith
C                 LUNSETFILE SET (27)             OCT 2025 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2025 Health Research Inc.,                          *
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
C  PURPOSE: HANDLES ALL INTERACTIONS WITH MRC IMAGE FILE HEADERS. 
C           CONTAINS NUMEROUS SUBROUTINES MOST START WITH PREFIX: LUN
C
C  SUBROUTINES:
C     ------------------------- LUNNEWHED_MRC ---------------------- 
C     ------------------------- LUNREDHED_MRC ---------------------- 
C     ------------------------- LUNWRTHED_MRC ---------------------- 
C     ------------------------- LUNGETOBJ_MRC ---------------------- 
C     ------------------------- LUNGETPOSINFO_MRC ------------------ 
C     ------------------------- LUNSETPOS_MRC ---------------------- 
C     ------------------------- LUNGETIS_MRC ----------------------- 
C     ------------------------- LUNGETMAP_MRC ---------------------- 
C     ------------------------- LUNSETHAND_MRC ---------------------
C     ------------------------- LUNGETHAND_MRC ---------------------
C     ------------------------- LUNSETXXX_MRC ----------------------
C     ------------------------- LUNSETMODSIZ_MRC -------------------
C     ------------------------- LUNGETMODSIZ_MRC -------------------
C     ------------------------- LUNGETMOD_MRC -------------------
C     ------------------------- LUNGETVIN_MRC ----------------------
C     ------------------------- LUNSETVIN_MRC ----------------------
C     ------------------------- LUNSETISBARE_MRC -------------------
C     ------------------------- LUNGETISBARE_MRC -------------------
C     ------------------------- LUNSETFILE_MRC ---------------------
C     ------------------------- LUNSETIMGNUM_MRC -------------------
C     ------------------------- LUNSET_STATSIMG_MRC ----------------
C     ------------------------- LUNGET_STATSIMG_MRC ----------------
C     ------------------------- LUNSETSTATS_MRC --------------------
C     ------------------------- LUNGETSTATS_MRC --------------------
C     ------------------------- LUNZEROSTATS_MRC -------------------
C     ------------------------- LUNGETVERSION_MRC ------------------
C     ------------------------- LUNGETPIXSIZ_MRC -------------------
C     ------------------------- LUNSETPIXSIZ_MRC -------------------
C     ------------------------- LUNSETPIXSIZES_MRC -----------------
C     ------------------------- LUNGETSIZE_MRC ---------------------
C     ------------------------- LUNGETVALS_R_MRC -------------------
C     ------------------------- LUNSETVALS_R_MRC -------------------
C     ------------------------- LUNGETTYPE_MRC --------------------- 
C     ------------------------- LUNGETSTK_MRC ---------------------- 
C     ------------------------- LUNSETSTK_MRC ---------------------- 
C     ------------------------- LUNGETNSTACK_MRC ------------------- 
C     ------------------------- LUNSETNSTACK_MRC ------------------- 
C     ------------------------- LUNSETCOMMON_MRC ------------------- 
C     ------------------------- LUNGETLABELS_MRC ------------------- 
C     ------------------------- LUNSAYINFO_MRC --------------------- 
C     ------------------------- LUNGETMODSIZES_MRC -----------------

C  REAL LOCATIONS:      
C
C      CELLAX   = TRANSFER(MRC_HEADER(11),FVAL)  ! PIXEL SIZE * NX
C      CELLAY   = TRANSFER(MRC_HEADER(12),FVAL)  ! PIXEL SIZE * NY
C      CELLAZ   = TRANSFER(MRC_HEADER(13),FVAL)  ! PIXEL SIZE * NZ
C
C      PIXSIZX  = CELLAX / NX                    ! PIXEL SIZE IN X 
C      PIXSIZY  = CELLAY / NY                    ! PIXEL SIZE IN Y 
C      PIXSIZZ  = CELLAZ / NZ                    ! PIXEL SIZE IN Z  
C
C      CELLBX   = TRANSFER(MRC_HEADER(14),FVAL)
C      CELLBY   = TRANSFER(MRC_HEADER(15),FVAL)
C      CELLBZ   = TRANSFER(MRC_HEADER(16),FVAL)
C
C      DMIN     = TRANSFER(MRC_HEADER(20),FVAL) ! MIN  DENSITY VALUE
C      DMAX     = TRANSFER(MRC_HEADER(21),FVAL) ! MAX  DENSITY VALUE  
C      DMEAN    = TRANSFER(MRC_HEADER(22),FVAL) ! MEAN DENSITY VALUE  
C
C      CAXIS    = TRANSFER(MRC_HEADER(26),CSTR(1:4))
C      EXTTYP   = TRANSFER(MRC_HEADER(27),CSTR(1:4))
C
C      ANG1     = TRANSFER(MRC_HEADER(43),FVAL)
C      ANG2     = TRANSFER(MRC_HEADER(44),FVAL)
C      ANG3     = TRANSFER(MRC_HEADER(45),FVAL)
C      ANG4     = TRANSFER(MRC_HEADER(46),FVAL)
C      ANG5     = TRANSFER(MRC_HEADER(47),FVAL)
C      ANG6     = TRANSFER(MRC_HEADER(48),FVAL)
C
C      ORX      = TRANSFER(MRC_HEADER(50),FVAL)
C      ORY      = TRANSFER(MRC_HEADER(51),FVAL)
C      ORZ      = TRANSFER(MRC_HEADER(52),FVAL)
C
C      MAP      = TRANSFER(MRC_HEADER(53),CSTR(1:4))
C
C      RMS      = TRANSFER(MRC_HEADER(55),FVAL) ! DEVIATION FROM MEAN DENSITY   
C
C
C
C  STATIC LOCATION:   257 -- IDSP
C                     258 -- ISBARE
C                     259 -- IMGNUM (CURRENT)  SPIDER: ISTACK
C                     260 -- NSTACK (CURRENT)  SPIDER: MAXIM
C
C                               MRC!=0        CONTENT
C  COMMON /LUNMRC/ LUNMRCPOS   MRC FILE       I/O OFFSET
C                  LUNMRCNBYT  UNSIGNED       BYTES / DATA VALUE
C  COMMON /LUNARA/ LUNARA      UL    LL       NX VALUE   
C                  LUNSTK                     IMGNUM OFFSET
C                  LUNARB                     DUPLICATE LUN POINTER
C                  LUNFLIP                    O = NOFLIP  1=FLIP BYTES
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      MODULE LUNMRCHDR_INFO
         INTEGER, PARAMETER  :: NUMLUNST = 100
         TYPE INTEGER_POINTER
            INTEGER, POINTER :: IPT(:) 
         END TYPE INTEGER_POINTER

         TYPE(INTEGER_POINTER) :: LUNMRCHDRBUF(NUMLUNST)
      END MODULE LUNMRCHDR_INFO

C     ----------- LUNNEWHED_MRC ---------------------------------------

      SUBROUTINE LUNNEWHED_MRC(LUN,IRTFLG)

C     CREATES STORAGE SPACE FOR A MRC_HEADER OBJECT

#include "LUNHDR.INC"

      INTEGER         :: LUN,IRTFLG
      
      INTEGER * 8     :: LUNMRCPOS 
      INTEGER         :: LUNMRCNBYT
      COMMON /LUNMRC/    LUNMRCPOS(100),LUNMRCNBYT(100)

      MPOINTER => LUNMRCHDRBUF(LUN)%IPT 
      IF (.NOT. ASSOCIATED(MPOINTER)) THEN
C        ALLOCATE SPACE FOR THIS HEADER OBJECT
         ALLOCATE(MPOINTER(LENMRCHDR),STAT=IRTFLG)

         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'LUNNEWHED_MRC, FILE HEADER',LENHDRMRC)
            IRTFLG = 1
            RETURN
         ENDIF
         LUNMRCHDRBUF(LUN)%IPT => MPOINTER
      ENDIF

C     TEMPORARY SETTING TO ENSURE LUN POINTS TO MRC BEFORE LENGTH KNOWN
      LUNMRCPOS(LUN) = 1

      !write(3,*) ' In lunnewhed_mrc: ',lun,LUNMRCHDRBUF(LUN)%ipt 
      IRTFLG = 0
      END

C     ----------- LUNREDHED_MRC ---------------------------------------

      SUBROUTINE LUNREDHED_MRC(LUN,CALLERRT,IRTFLG)

C     READS MRC IMAGE HEADER FROM FILE INTO MRC HEADER OBJECT 

#include "LUNHDR.INC"

      INTEGER      :: LUN,IRTFLG 
      LOGICAL      :: CALLERRT

      INTEGER * 8  :: IPOSMRC

      character(len=4) :: caxis

 
C     INCLUSION FOR OPTIONAL MPI INITIALIZATION.  
      INTEGER      :: MYPID = -1
#include "MPI_INIT.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     READ 1024 BYTEs LONG HEADER FROM MRC FILE INTO MRC_HEADER OBJECT

      IPOSMRC = 1          ! HEADER START POSITION

C     USING IOSTAT; IRTFLG IS SET TO ZERO ON SUCCESSFUL READ. 
      IF (MYPID <= 0) 
     &    READ(LUN, POS=IPOSMRC,IOSTAT=IRTFLG) MRC_HEADER(1:256)


#ifdef USE_MPI
      IF (ONLYONE_RED) THEN
C       ALWAYS BROADCASTS IF ONLYONE_RED IS .TRUE. WHEN USING MPI
        CALL BCAST_MPI('LUNREDHED_MRC','MRC_HEADER',MRC_HEADER,256,'I',
     &                  ICOMM,MPIERR)
        CALL BCAST_MPI('LUNREDHED_MRC','IRTFLG',IRTFLG, 1,'I',
     &                  ICOMM,MPIERR)
      ENDIF
#endif



#if defined(SP_DBUG)
       caxis = transfer(mrc_header(26),caxis) 
       write(3,*)' In lunredhed_mrc, mrc_header: ',mrc_header(26)
       write(3,*)' In lunredhed_mrc, caxis:',caxis
#endif

      IF (MYPID <= 0 .AND. IRTFLG .NE. 0 .AND. CALLERRT) THEN
         CALL ERRT(102,'I/O ERROR ON FILE HEADER',IRTFLG)
         IRTFLG = 1
         RETURN
      ENDIF

      IRTFLG = 0
      END

C     ----------- LUNWRTHED_MRC ---------------------------------------

      SUBROUTINE LUNWRTHED_MRC(LUN,IRTFLG)

C     WRITES HEADER OBJECT INTO MRC FILE'S HEADER (ONLY HAS 1 HEADER)

#include "LUNHDR.INC"

      INTEGER    :: LUN,IRTFLG

      INTEGER *8 :: IPOSMRC = 1    ! LONG INTEGER
      real       :: fval

C     INCLUSION FOR OPTIONAL MPI INITIALIZATION.  
      INTEGER    :: MYPID = -1

#include "MPI_INIT.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     WRITE HEADER RECORD FROM MRC_HEADER INTO MRC FILE 

      IF (MYPID <= 0) 
     &    WRITE(LUN, POS=IPOSMRC,IOSTAT=IRTFLG) MRC_HEADER(1:256)

#if defined(SP_DBUGIO)
      write(3,*) ' '
      write(3,*) ' In lunwrthed_mrc; iposmrc,irtflg: ',iposmrc,irtflg
      write(3,*) ' In lunwrthed_mrc; header(25): ',    mrc_header(25)
      write(3,*) ' '
      !CALL BACKTRACE
#endif


#ifdef USE_MPI
      IF (ONLYONE_RED) THEN
C       ALWAYS BROADCASTS IF ONLYONE_RED IS .TRUE. WHEN USING MPI
        CALL BCAST_MPI('LUNWRTHED_MRC','MRC_HEADER',MRC_HEADER,256,'I',
     &                  ICOMM,MPIERR)
        CALL BCAST_MPI('LUNWRTHED_MRC','IRTFLG',IRTFLG, 1,'I',
     &                  ICOMM,MPIERR)
      ENDIF
#endif

      !fval = transfer(mrc_header(20),fval)

      IF (MYPID <= 0 .AND. IRTFLG .NE. 0) THEN

         CALL ERRT(102,'WRITING TO MRC FILE HEADER',IRTFLG)
         IRTFLG = 1
      ENDIF

      IRTFLG = 0
      END

C     ------------------------- LUNGETOBJ_MRC -------------------------

      SUBROUTINE LUNGETOBJ_MRC(LUN,IPOINTER,IRTFLG)

      USE LUNMRCHDR_INFO

      INTEGER                        :: LUN,IRTFLG
      INTEGER, DIMENSION(:), POINTER :: IPOINTER 

      INTEGER, PARAMETER             :: NUMLUNS = 100

C     POINT TO HEADER OBJECT
 
      IRTFLG = 1
      IF (LUN <= 0 .OR. LUN > NUMLUNS) THEN
         CALL ERRT(102,'PGM ERROR, LUN OUT OF RANGE', LUN)
         RETURN
      ENDIF

      IPOINTER => LUNMRCHDRBUF(LUN)%IPT
      IF (.NOT. ASSOCIATED(IPOINTER)) THEN
         CALL ERRT(101,'UNASSOCIATED MRC_HEADER POINTER', IRTFLG)
         RETURN
      ENDIF

      IRTFLG = 0
      END

C     ------------------ LUNGETPOSINFO_MRC ----------------------------

      SUBROUTINE LUNGETPOSINFO_MRC(LUN,NBYT_PER_VAL,IHEDLEN,IRTFLG)

#include "LUNHDR.INC"
     
      INTEGER         :: LUN,NBYT_PER_VAL,IHEDLEN,IRTFLG

      INTEGER         :: NSYMBT,MRCMODE

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET NSYMBT = HEADER EXTRA LENGTH
      NSYMBT  = MRC_HEADER(24)  ! NO. OF ADDITIONAL HEADER BYTES

C     CALCULATE MRC FILE HEADER LENGTH IN BYTES
      IHEDLEN = 1024 + NSYMBT

C     CALCULATE NUMBER OF BYTES IN AN IMAGE VALUE (<0 == UNSIGNED)

      IRTFLG  = 0
      MRCMODE = MRC_HEADER(4)

      IF     (MRCMODE == 0) THEN
         NBYT_PER_VAL =  1         ! SIGNED    8 BIT INTEGER*1
 
      ELSEIF (MRCMODE == 1) THEN
         NBYT_PER_VAL =  2         ! SIGNED   16 BIT INTEGER*2

      ELSEIF (MRCMODE == 6) THEN
         NBYT_PER_VAL = -2         ! UNSIGNED 16 BIT INTEGER*2

      ELSEIF (MRCMODE == 2) THEN
         NBYT_PER_VAL =  4         ! REAL*4
  
      ELSE
         NBYT_PER_VAL = 0          ! NOT A VALID MRC FILE
         IRTFLG = 1
         RETURN
      ENDIF

      IRTFLG = 0
      END


C     ------------------ LUNSET_HEDPOS_MRC ----------------------------

      SUBROUTINE LUNSET_HEDPOS_MRC(LUN,IRTFLG)

      IMPLICIT NONE
     
      INTEGER        :: LUN,IRTFLG

      INTEGER * 8    :: LUNMRCPOS 
      INTEGER        :: LUNMRCNBYT
      COMMON /LUNMRC/   LUNMRCPOS(100),LUNMRCNBYT(100)

      INTEGER        :: LUNARA,LUNSTK,LUNARB,LUNFLIP
      COMMON /LUNARA/   LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

      INTEGER        :: NX,NY,NZ

      CALL LUNGETSIZE_MRC(LUN,NX,NY,NZ,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      LUNARA(LUN)    = NX

C     NO STACK OFFSET NEEDED
      LUNMRCPOS(LUN) = 1

      IRTFLG = 0
      END

C     ------------------ LUNSETPOS_MRC ----------------------------

      SUBROUTINE LUNSETPOS_MRC(LUN,IMGNUM,IRTFLG)

      IMPLICIT NONE
     
      INTEGER        :: LUN,IMGNUM,IRTFLG

      INTEGER * 8    :: LUNMRCPOS 
      INTEGER        :: LUNMRCNBYT
      COMMON /LUNMRC/   LUNMRCPOS(100),LUNMRCNBYT(100)

      INTEGER        :: LUNARA,LUNSTK,LUNARB,LUNFLIP
      COMMON /LUNARA/   LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

      INTEGER         :: NBYT,IHEDLEN,NBYT_PER_VAL,NBYT_PER_REC
      INTEGER         :: NX,NY,NZ
      CHARACTER(LEN=4):: CAXIS

C     GET AXIS ORIGIN IN FILE
C     THIS IS NOT A STANDARD MRC DEFINED HEADER POSITION!!
      CALL LUNGETHAND_MRC(LUN,CAXIS,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET NUMBER OF BYTES IN AN IMAGE VALUE AND HEADER LENGTH
      CALL LUNGETPOSINFO_MRC(LUN,NBYT_PER_VAL,IHEDLEN,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     NZ IS ADJUSTED FOR RELION .mrcs FILES IF NEEDED
      CALL LUNGETSIZE_MRC(LUN,NX,NY,NZ,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      LUNMRCNBYT(LUN) = NBYT_PER_VAL  ! PRESERVES NEGATIVE FOR SIGN

      NBYT_PER_VAL    = ABS(NBYT_PER_VAL)   
      NBYT_PER_REC    = NBYT_PER_VAL * NX

 
      IF (CAXIS(1:1) == 'U') THEN   ! NON-DEFAULT
C        ORIGIN AT UPPER LEFT, POSITION AT ONE RECORD BEFORE IMAGE

         LUNARA(LUN) = -NX

         IF (IMGNUM <= 1) THEN
C          NO STACK OFFSET NEEDED
           LUNMRCPOS(LUN) = IHEDLEN - NBYT_PER_REC + 1

         ELSE
C          STACKED FILE OFFSET NEEDED 
           LUNMRCPOS(LUN) = IHEDLEN + 
     &                      (IMGNUM -1) * NZ * NY * NBYT_PER_REC -
     &                                              NBYT_PER_REC + 1
         ENDIF   ! END OF: IF (IMGNUM <= 1)
 
      ELSE  ! DEFAULTS IF THIS POSITION IS BLANK
C        ORIGIN AT LOWER LEFT, POSITION AT END OF CURRENT IMAGE +1

         LUNARA(LUN) = NX

         IF (IMGNUM <= 1) THEN
C          NO STACK OFFSET NEEDED
           LUNMRCPOS(LUN) = IHEDLEN + 
     &                      (NZ - 1) * NY * NBYT_PER_REC + 
     &                                 NY * NBYT_PER_REC + 1
         ELSE
C          STACKED FILE OFFSET NEEDED 
           LUNMRCPOS(LUN) = IHEDLEN + 
     &                      (IMGNUM -1) * NZ * NY * NBYT_PER_REC +
     &                      (NZ - 1)    *      NY * NBYT_PER_REC + 
     &                                         NY * NBYT_PER_REC + 1
         ENDIF   ! END OF: IF (IMGNUM <= 1)
      ENDIF

      !write(3,*)' In lunsetpos_mrc; imgnum,lunmrcpos(lun): ',
      !&                             imgnum,lunmrcpos(lun)
      !write(3,*)' In lunsetpos_mrc; ihedlen,nx,ny,nz,nbyt,imgnum: ',
      !&                             ihedlen,nx,ny,nz,nbyt,imgnum
      !write(3,*)' In lunsetpos_mrc; nbyt_per_rec: ', nbyt_per_rec
      !write(3,*)' In lunsetpos_mrc; lunmrcpos: ',lunmrcpos(lun)

      IRTFLG = 0
      END

C     ------------------------- LUNFLIPORG_MRC ------------------------

      SUBROUTINE LUNFLIPORG_MRC(LUN,NX,NY,NZ,RESET,IRTFLG)

#include "LUNHDR.INC"

      INTEGER        :: LUN,NX,NY,NZ,IRTFLG
      LOGICAL        :: RESET

      INTEGER * 8    :: LUNMRCPOS 
      INTEGER        :: LUNMRCNBYT
      COMMON /LUNMRC/   LUNMRCPOS(100),LUNMRCNBYT(100)

      INTEGER        :: LUNARA,LUNSTK,LUNARB,LUNFLIP
      COMMON /LUNARA/   LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

      CHARACTER(LEN=4)  :: CAXIS
      INTEGER * 8, SAVE :: LUNMRCPOS_SAV = 0 
      INTEGER           :: NBYT_PER_REC

      IF (LUNMRCPOS(LUN) == 0) THEN
C       DIRECT ACCESS SPIDER FILE
        RETURN

      ELSEIF (RESET) THEN
C       RESET TO ORIGINAL AXIS ORIGIN IN FILE
        LUNMRCPOS(LUN) = LUNMRCPOS_SAV 
        IRTFLG = 0
        RETURN

      ELSE
C       DIRECT ACCESS MRC FILE

C       GET AXIS ORIGIN IN FILE NOT A STANDARD MRC DEFINED POSITION!!
        CALL LUNGETHAND_MRC(LUN,CAXIS,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (CAXIS(1:1) == 'L') THEN    
C         ORIGIN AT LOWER LEFT, BUT POSITION AT ONE RECORD BEFORE IMAGE
          NBYT_PER_REC   = LUNMRCNBYT(LUN)
          LUNMRCPOS_SAV  = LUNMRCPOS(LUN)
          LUNMRCPOS(LUN) = LUNMRCPOS(LUN) - NZ * NY * NBYT_PER_REC

        ENDIF   ! END OF: IF (CAXIS(1:1) == 'L')
      ENDIF
      IRTFLG = 0

      END

C     ------------------------- LUNGETIS_MRC --------------------------

      SUBROUTINE LUNGETIS_MRC(LUN,IS_MRC,IRTFLG)

C     PURPOSE:  QUERIES HEADER LOCATION THAT SPECIFIES MRC FILE

      IMPLICIT NONE

      INTEGER         :: LUN,IRTFLG
      LOGICAL         :: IS_MRC

      INTEGER * 8     :: LUNMRCPOS 
      INTEGER         :: LUNMRCNBYT
      COMMON /LUNMRC/    LUNMRCPOS(100),LUNMRCNBYT(100)

C     Handle invalid index (LUN=0, etc.)  Sept 2025
      IF (LUN >= LBOUND(LUNMRCPOS, 1) .AND. 
     &    LUN <= UBOUND(LUNMRCPOS, 1)) THEN
          IS_MRC = (LUNMRCPOS(LUN) .NE. 0)
      ELSE
          IS_MRC = .FALSE.
      ENDIF

      IRTFLG = 0
      END

C     ------------------------- LUNGETMAP_MRC -------------------------

      SUBROUTINE LUNGETMAP_MRC(LUN,MAPC,MAPR,MAPS,IRTFLG)

#include "LUNHDR.INC"
      INTEGER         :: LUN,MAPC,MAPR,MAPS,IRTFLG

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      MAPC = MRC_HEADER(17)
      MAPR = MRC_HEADER(18)
      MAPS = MRC_HEADER(19)

#if defined(SP_DBUG)
        write(3,*) ' In lungetmap_mrc; mrc_header 17,18,19:  ',
     &        mrc_header(17),' ',mrc_header(18),' ',mrc_header(19)      
        write(3,*) ' In lungetmap_mrc; map c,r,s: ',mapc,mapr,maps
#endif

      IRTFLG = 0
      END

C     ------------------------- WHICH_HAND_MRC ------------------------

      SUBROUTINE WHICH_HAND_MRC(LUN,FILNAM,CAXIS,IRTFLG)

C     PURPOSE:  SET AXIS ORIGIN LOCATION & VOLUME HANDEDNESS 
C               BEFORE LUNSETPOS_MRC

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'

      CHARACTER(LEN=*)    :: FILNAM
      CHARACTER(LEN=4)    :: CAXIS
      INTEGER             :: LUN,IRTFLG

      LOGICAL             :: ISMRCSFILE

C     SET AXIS ORIGIN LOCATION & VOLUME HANDEDNESS BEFORE 
C              CALLING LUNSETPOS_MRC

C     SEE IF RELION COMPATIBLE IMAGE/STACK (ALWAYS 'UL L')
      ISMRCSFILE = (INDEX(FILNAM,'.MRCS') > 1  .OR.
     &              INDEX(FILNAM,'.mrcs') > 1)
        
#if defined(SP_DBUG)
      write(3,*)' In lunwhichhand_mrc; ismrcsfile: ',ismrcsfile
      write(3,*)' In lunwhichhand_mrc; caxis   1111111: ',caxis
#endif

      IF (ISMRCSFILE) THEN
C        RELION COMPATIBLE IMAGE/STACK (ALWAYS 'UL L')
         CAXIS = 'UL L'

         CALL LUNSETHAND_MRC(LUN,CAXIS,IRTFLG)
         !write(3,*)' In which_hand_mrc; called sethand 22222: ',caxis
      
      ELSE
C        SEE IF FILE HEADER ALREADY HAS CAXIS SET
         CALL LUNGETHAND_MRC(LUN,CAXIS,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         IF (INDEX(CAXIS,'L') == 0) THEN
C           CAXIS DOES NOT CONTAIN: 'UL L', 'UL R', 'LL L' or 'LL R'
C                (ALWAYS HAS L)

C           SET ORIGIN & HANDEDNESS TO CURRENT DEFAULT MRC_AXIS
            CAXIS = MRC_AXIS
            CALL LUNSETHAND_MRC(LUN,CAXIS,IRTFLG)

            !write(3,*)' In which_hand_mrc; called sethand 3333: ',caxis
        ENDIF       
      ENDIF 

      END     

C     ------------------------- LUNSETHAND_MRC ------------------------

      SUBROUTINE LUNSETHAND_MRC(LUN,CAXIS,IRTFLG)

C     SETS  AXIS ORIGIN LOCATION & VOLUME HANDEDNESS (CHARACTERS)
C     THIS IS NOT A STANDARD MRC DEFINED HEADER POSITION!!

#include "LUNHDR.INC"

      INTEGER            :: LUN,IRTFLG
      CHARACTER(LEN=4)   :: CAXIS

      INTEGER            :: I4V

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      MRC_HEADER(26) = TRANSFER(CAXIS(1:4),I4V)

#if defined(SP_DBUG)
      write(3,*)' In lunsethand_mrc; in header(26) caxis: ',caxis
#endif

      IRTFLG = 0
      END

C     ------------------------- LUNGETHAND_MRC ------------------------

      SUBROUTINE LUNGETHAND_MRC(LUN,CAXIS,IRTFLG)

C     GETS  AXIS ORIGIN LOCATION & VOLUME HANDEDNESS (CHARACTERS)
C     THIS IS NOT A STANDARD MRC DEFINED HEADER POSITION!!

#include "LUNHDR.INC"

      INTEGER            :: LUN,IRTFLG
      CHARACTER(LEN=4)   :: CAXIS

      CHARACTER(LEN=4)   :: CSTR

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      CAXIS(1:4) = TRANSFER(MRC_HEADER(26),CSTR(1:4))

      !write(3,*)' In lungethand_mrc; lun,caxis: ',lun,caxis

      IRTFLG = 0
      END

C     ------------------------- LUNSETXXX_MRC -------------------------

      SUBROUTINE LUNSETXXX_MRC(LUN,IRTFLG)

C     SETS  N?START,CELLANGS,AXES,NSYMBYT,MAP,MACHST,LABELS,ETC.
#include "LUNHDR.INC"

      INTEGER            :: LUN,IRTFLG

      INTEGER            :: I,INOW
      INTEGER            :: NSYMBYT,NLABEL
      CHARACTER(LEN=1)   :: BLANK = CHAR(32)
      CHARACTER(LEN=4)   :: CVAL
      CHARACTER(LEN=12)  :: CDATT
      CHARACTER(LEN=36)  :: LABEL1
      INTEGER            :: I4V

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET DEFAULT HEADER POSITIONS

      MRC_HEADER(5)  = 1  ! NXSTART  NUMBER OF FIRST COLUMN  IN MAP
      MRC_HEADER(6)  = 1  ! NYSTART  NUMBER OF FIRST ROW     IN MAP
      MRC_HEADER(7)  = 1  ! NZSTART  NUMBER OF FIRST SECTION IN MAP

C     CELL ANGLES (DEGREES)  (CELLBX..) 
      MRC_HEADER(14) = 90.0
      MRC_HEADER(15) = 90.0
      MRC_HEADER(16) = 90.0

C     STORAGE: FAST, MEDIUM, SLOW AXES 
      MRC_HEADER(17) = 1
      MRC_HEADER(18) = 2
      MRC_HEADER(19) = 3

      NSYMBYT        = 0          ! NO EXTRA BYTES IN HEADER
      MRC_HEADER(24) = NSYMBYT

C     ZERO THE EXTRA POSITIONS  25...49, BUT 27=MRC0, 28=IVERSION
      MRC_HEADER(25:26) = 0
      MRC_HEADER(29:49) = 0

C     PUT IN EXTTYP : 'MRCO'
      CVAL = 'MRCO'
      MRC_HEADER(27) = TRANSFER(CVAL(1:4),I4V)

C     PUT IN 'MAP'
      CVAL = 'MAP '
      MRC_HEADER(53) = TRANSFER(CVAL(1:4),I4V)

C     SET MACHINE STAMP
C     MACHST = 286326784    ! FOR BIG ENDIAN DATA
C     MACHST = 4369         ! FOR LITTLE ENDIAN DATA ??
      MACHST = 16708        ! FOR LITTLE ENDIAN DATA ??
      MRC_HEADER(54) = MACHST

C     SET NUMBER OF LABELS
      NLABEL = 1  
      MRC_HEADER(56) = NLABEL

C     ZERO ALL LABELS WITH BLANKS
      CVAL = BLANK // BLANK // BLANK // BLANK
      DO I = 57,256
         MRC_HEADER(I) = TRANSFER(CVAL(1:4),I4V)
      ENDDO

C     ADD ONE LABEL
      CALL DATE_2K(CDATT)            !12234567890123
C               123456789 123456789 123456789 123456789 
      LABEL1 = 'Created using SPIDER: ' // CDATT
      INOW   = 56
      DO I = 1,36,4
         CVAL = LABEL1(I:I+3)
         INOW = INOW + 1
         MRC_HEADER(INOW) = TRANSFER(CVAL(1:4),I4V)
      ENDDO

      IRTFLG = 0
      END

C     ------------------------- LUNSETMODSIZ_MRC -----------------------

      SUBROUTINE LUNSETMODSIZ_MRC(LUN,MRCMODE,NX,NY,NZMRC,
     &                            MX,MY,MZ,IRTFLG)

#include "LUNHDR.INC"

      INTEGER  :: LUN,MRCMODE,NX,NY,NZMRC,MX,MY,MZ,IRTFLG

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      MRC_HEADER(4)  = MRCMODE
      MRC_HEADER(1)  = NX
      MRC_HEADER(2)  = NY
      MRC_HEADER(3)  = NZMRC ! NO. OF IMAGES IN RELION mrcs FILES!
      MRC_HEADER(8)  = MX    ! NO. OF INTERVALS ALONG X
      MRC_HEADER(9)  = MY    ! NO. OF INTERVALS ALONG Y
      MRC_HEADER(10) = MZ    ! NO. OF INTERVALS ALONG Z (>1 = STACK)

      !write(3,*) ' In lunsetmodsiz_mrc: ',mrcmode,nx,ny,nz,mx,my,mz

      IRTFLG = 0
      END

C     ------------------------- LUNSETMODE_MRC ------------------------

      SUBROUTINE LUNSETMODE_MRC(LUN,MRCMODE,IRTFLG)

#include "LUNHDR.INC"

      INTEGER  :: LUN,MRCMODE,IRTFLG

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      MRC_HEADER(4) = MRCMODE

      !write(3,*) ' In lunsetmode_mrc; mrcmode: ',mrcmode

      IRTFLG = 0
      END


C     ------------------------- LUNGETMODSIZES_MRC --------------------

      SUBROUTINE LUNGETMODSIZES_MRC(LUN,NX,NY,NZ,
     &                                  MX,MY,MZ,IRTFLG)

#include "LUNHDR.INC"

      INTEGER  :: LUN,NX,NY,NZ, MX,MY,MZ,IRTFLG

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      NX      = MRC_HEADER(1)
      NY      = MRC_HEADER(2)
      NZ      = MRC_HEADER(3)  ! NO. OF IMAGES IN RELION mrcs FILES! 
      MX      = MRC_HEADER(8)  ! NO. OF INTERVALS ALONG X
      MY      = MRC_HEADER(9)  ! NO. OF INTERVALS ALONG Y
      MZ      = MRC_HEADER(10) ! NO. OF INTERVALS ALONG Z (>1 = STACK)

      !write(3,*) ' In lungetmodsizes_mrc: ',mrcmode,nx,ny,nz,mx,my,mz

      IRTFLG = 0
      END

C     ------------------------- LUNGETMODE_MRC ----------------------

      SUBROUTINE LUNGETMODE_MRC(LUN,MRCMODE,IRTFLG)

#include "LUNHDR.INC"

      INTEGER  :: LUN,MRCMODE,IRTFLG

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      MRCMODE = MRC_HEADER(4)

      !write(3,*) ' In lungetmod_mrc; mrcmode: ',mrcmode

      IRTFLG = 0
      END

C     ------------------------- LUNGETVIN_MRC -------------------------

      SUBROUTINE LUNGETVIN_MRC(LUN,IVERSION,ISPG,NSYMBT,NLABL,IRTFLG)

#include "LUNHDR.INC"
      INTEGER         :: LUN,IVERSION,ISPG,NSYMBT,NLABL,IRTFLG

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
 
C     (23)  ISPG == 0        IMAGE OR IMAGE STACK
C     (23)  ISPG == 1        VOLUMES
C     (23)  ISPG == 401      STACK OF EM VOLUMES

      IVERSION = MRC_HEADER(28)
      ISPG     = MRC_HEADER(23)
      NSYMBT   = MRC_HEADER(24)  ! NO. OF ADDITIONAL HEADER BYTES
      NLABL    = MRC_HEADER(56)  ! NO. OF LABELS

      ! write(3,*) ' In lungetvin_mrc; iversion', iversion

      IRTFLG = 0
      END

C     ------------------------- LUNSETVIN_MRC -------------------------

      SUBROUTINE LUNSETVIN_MRC(LUN,IVERSION,ISPG,NSYMBT,NLABL,IRTFLG)

#include "LUNHDR.INC"
      INTEGER         :: LUN,IVERSION,ISPG,NSYMBT,NLABL,IRTFLG

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
   
C     (23)  ISPG == 0        IMAGE OR IMAGE STACK
C     (23)  ISPG == 1        VOLUMES
C     (23)  ISPG == 401      STACK OF EM VOLUMES

      MRC_HEADER(28) = IVERSION
      MRC_HEADER(23) = ISPG
      MRC_HEADER(24) = NSYMBT ! NUMBER OF ADDITIONAL HEADER BYTES
      MRC_HEADER(56) = NLABL  ! NO. OF LABELS

      !write(3,*) ' In lunsetvin_mrc; iversion', iversion,ispg

      IRTFLG = 0
      END

C     ------------------------- LUNGETISBARE_MRC ----------------------

      SUBROUTINE LUNGETISBARE_MRC(LUN,ISBARE,IRTFLG)

#include "LUNHDR.INC"

      INTEGER   :: LUN,IRTFLG
      LOGICAL   :: ISBARE

 
C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET ISBARE FROM HEADER OBJECT STATIC LOCATION
      ISBARE = (MRC_HEADER(258) == 1)

      IRTFLG = 0
      END

C     --------------------- LUNSETISBARE_MRC -------------------------

      SUBROUTINE LUNSETISBARE_MRC(LUN,ISBARE,IRTFLG)

#include "LUNHDR.INC"

      INTEGER     :: LUN,IRTFLG
      LOGICAL     :: ISBARE

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET ISBARE IN MRC HEADER OBJECT (STATIC AREA)
      IF (ISBARE) THEN
         MRC_HEADER(258) = 1 
      ELSE
         MRC_HEADER(258) = 0
      ENDIF

      IRTFLG = 0
      END

C     ------------------------- LUNSETFILE_MRC -----------------------

      SUBROUTINE LUNSETFILE_MRC(LUN,FILNAM, DSP,IRTFLG)

#include "LUNHDR.INC"

      INTEGER                :: LUN,IRTFLG
      CHARACTER(LEN=*)       :: FILNAM
      CHARACTER(LEN=1)       :: DSP
      CHARACTER(LEN=1)       :: NULL = CHAR(0)

      INTEGER                :: NLET,LOCDOT
      CHARACTER(LEN=4)       :: FIL_EXT

      IF (FILNAM(1:1) .NE. NULL) THEN
C        SET CURRENT FILENAME IN HEADER OBJECT 
         NLET           = LNBLNKN(FILNAM)
         LUNFILNAM(LUN) = FILNAM(1:NLET)

#if defined(NEVER_SP_DBUGIO)
C23456
         write(3,*) '  '
         ilent = lnblnkn(filnam)
         write(3,*)' In lunsetfile; filnam: ', filnam(1:ilent)
         write(3,*)' In lunsetfile; lunfilnam(lun): ', lunfilnam(lun)
#endif
      ENDIF

      IF (DSP .NE. NULL) THEN

C        POINT TO HEADER OBJECT
         CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

C        SET DSP IN STATIC HEADER OBJECT
         MRC_HEADER(257) = 0
         IF (DSP == 'N') MRC_HEADER(257) = 1

         IF (FILNAM(1:1) .NE. NULL) THEN
            NLET     = LNBLNKN(FILNAM)
            LOCDOT   = INDEX(FILNAM(1:NLET),'.',BACK=.TRUE.)

            IF (LOCDOT > 0) THEN
               FIL_EXT = FILNAM(LOCDOT+1:NLET)
               NLET    = lnblnkn(FIL_EXT)

               IF (NLET < 4) FIL_EXT(4:4) = ' '         
               MRC_HEADER(27) = TRANSFER(FIL_EXT(1:4),I4V)
            ENDIF
         ENDIF


      ENDIF



      IRTFLG = 0
      END

C     ------------------------ LUNGETFILE_MRC -----------------------

      SUBROUTINE LUNGETFILE_MRC(LUN,FILNAM,NLET,DSP,IRTFLG)

#include "LUNHDR.INC"

      INTEGER                :: LUN,NLET,IRTFLG
      CHARACTER(LEN=*)       :: FILNAM
      CHARACTER(LEN=1)       :: DSP

C     RETRIEVE CURRENT FILENAME         ! NOT IN HEADER OBJECT
      FILNAM = LUNFILNAM(LUN)
      NLET   = LNBLNKN(FILNAM)

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET DSP FROM STATIC HEADER OBJECT
      DSP = 'O'
      IF (MRC_HEADER(257) == 1) DSP = 'N'

      IRTFLG = 0
      END

C     ------------------------- LUNSETIMGNUM_MRC ----------------------

      SUBROUTINE LUNSETIMGNUM_MRC(LUN,IMGNUM,IRTFLG)

#include "LUNHDR.INC"

      INTEGER  :: LUN,IMGNUM,IRTFLG

C     POINT TO MRC HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
          
C     SET MRC_HEADER VALUE
      MRC_HEADER(259) = IMGNUM       ! LOCATION: 259  CURRENT IMGNUM
         
      IRTFLG = 0
      END

C     ------------------------- LUNSET_STATSIMG_MRC -------------------

      SUBROUTINE LUNSET_STATSIMG_MRC(LUN,IMGSTATS,IRTFLG)

#include "LUNHDR.INC"

      INTEGER  :: LUN,IMGSTATS,IRTFLG

C     IMGSTATS IS THE CURRENT STACKED OR MAIN IMAGE NUMBER IN STATS

C     POINT TO MRC HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      !write(3,*)' In lunset_statsimg_mrc; imgstats: ',imgstats
          
C     SET MRC_HEADER VALUE
      MRC_HEADER(25) = IMGSTATS ! LOCATION: 25  IMGNUM FOR STATISTICS
         
      IRTFLG = 0
      END

C     ------------------------- LUNGET_STATSIMG_MRC -------------------

      SUBROUTINE LUNGET_STATSIMG_MRC(LUN,IMGSTATS,IMGNUM,NSTACK,IRTFLG)

#include "LUNHDR.INC"

      INTEGER  :: LUN,IMGSTATS,IMGNUM,NSTACK,IRTFLG
C     IMGSTATS IS THE CURRENT STACKED OR MAIN IMAGE NUMBER IN STATS

C     POINT TO MRC HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
   
C     GET MRC_HEADER VALUE
      IMGSTATS = MRC_HEADER(25)  ! LOCATION:  25 IMGNUM FOR STATISTICS
      IMGNUM   = MRC_HEADER(259) ! LOCATION: 259 IMGNUM (<0 != STACK)
      NSTACK   = MRC_HEADER(260) ! LOCATION: 260 NSTACK (<0 != STACK)
         
      IRTFLG = 0
      END

C     ------------------------- LUNSETSTATS_MRC -----------------------

      SUBROUTINE LUNSETSTATS_MRC(LUN,IMAMIT,FMINT,FMAXT,AVT,SIGT,IRTFLG)

#include "LUNHDR.INC"

      INTEGER  :: LUN,IMAMIT,IRTFLG
      REAL     :: FMINT,FMAXT,AVT,SIGT,f21

      INTEGER  :: MZ,NSTACK,IMGNUM,IVERSION,IMGSTATS
      INTEGER  :: I4V

C     IMGSTATS IS THE CURRENT STACKED OR MAIN IMAGE NUMBER IN STATS

C     POINT TO MRC HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF (IMAMIT == 0 .OR. (FMAXT < FMINT) ) THEN

C        SET UNDETERMINED FMIN... STATISTICS
         CALL LUNGETVERSION_MRC(LUN,IVERSION,IRTFLG)  ! MUST BE SET FIRST
         IF (IRTFLG .NE. 0) RETURN

         CALL LUNZEROSTATS_MRC(LUN,IVERSION,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         IMGSTATS = 0

      ELSE
C        SET MRC_HEADER VALUES
         MRC_HEADER(20) = TRANSFER(FMINT,I4V) ! DMAX   MIN  DENSITY VALUE
         MRC_HEADER(21) = TRANSFER(FMAXT,I4V) ! DMIN   MAX  DENSITY VALUE  
         MRC_HEADER(22) = TRANSFER(AVT,  I4V) ! DMEAN  MEAN DENSITY VALUE  
         MRC_HEADER(55) = TRANSFER(SIGT, I4V) ! RMS    MAP DEVIATION FROM MEAN DENSITY          

C        NEED IMGNUM FOR IMGSTATS
         CALL LUNGETSTK_MRC(LUN,MZ,NSTACK,IMGNUM,IRTFLG)
        
         IMGSTATS = IMGNUM

      ENDIF

C     SET NUMBER OF IMAGE IN STATS
      CALL LUNSET_STATSIMG_MRC(LUN,IMGSTATS,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
       

      IRTFLG = 0
      END

C     ------------------------- LUNGETSTATS_MRC -----------------------

      SUBROUTINE LUNGETSTATS_MRC(LUN,IMAMIT,FMINT,FMAXT,AVT,SIGT,IRTFLG)

#include "LUNHDR.INC"

      INTEGER  :: LUN,IMAMIT,IRTFLG
      REAL     :: FMINT,FMAXT,AVT,SIGT,f21

      INTEGER  :: IMGSTATS,IMGNUM,NSTACK

C     POINT TO MRC HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET MRC_HEADER VALUES
      FMINT = TRANSFER(MRC_HEADER(20),FVAL) ! DMAX:  MIN  DENSITY VALUE
      FMAXT = TRANSFER(MRC_HEADER(21),FVAL) ! DMIN:  MAX  DENSITY VALUE  
      AVT   = TRANSFER(MRC_HEADER(22),FVAL) ! DMEAN: MEAN DENSITY VALUE  
      SIGT  = TRANSFER(MRC_HEADER(55),FVAL) ! RMS: DEVIATION FROM MEAN DENSITY          

C     MRC FILE ENCODES IMAMI IN RELATIVE VALUES FOR MIN/MAX      
      IMAMIT = 1
      IF (FMAXT < FMINT) IMAMIT = 0
      IF (AVT   < FMINT)  IMAMIT = 0
      IF (SIGT  < 0.0)   IMAMIT = 0

C     IF ALL ZERO STATISTICS ARE PROBABLY NOT SET!! 
      IF (FMINT == 0.0 .AND. FMAXT == 0.0 .AND. SIGT == 0.0) IMAMIT = 0 
  
C     ARE STATS FOR THE CURRENT MRC IMAGE?
      CALL LUNGET_STATSIMG_MRC(LUN,IMGSTATS,IMGNUM,NSTACK,IRTFLG)

      IF (IMGNUM > 0 .AND. (IMGSTATS .NE. IMGNUM)) THEN
C        THIS IS A MRC STACKED IMAGE BUT IMG NOT = IMG IN STATISTICS
         IMAMIT = 0
         !FMAXT  = FMINT - 1

#if defined(SP_DBUG)
!      f21 = TRANSFER(mrc_header(21),f21)

!      write(3,'(A,I5,2X,F8.1,2X,F8.2)')
!     &     '  In lunGetstats_mrc; lun,fmaxt,mrc_header(21):',
!     &        lun,fmaxt,f21

!      write(3,'(A,I5,2X,3(2x,I8))')
!     &     '  In lunGetstats_mrc; nstack,imgnum,imgstats,imamit:',
!     &                            nstack,imgnum,imgstats,imamit
C      write(3,*)'  '
#endif

         CALL LUNSETSTATS_MRC(LUN,IMAMIT,FMINT,FMAXT,AVT,SIGT,IRTFLG)

      ELSEIF (NSTACK > 1 .AND. IMGNUM < 0 ) THEN
C        THIS IS A MRC STACK BUT NO IMAGE SPECIFIED
         IMAMIT = 0
         FMAXT  = FMINT - 1
      ENDIF

#if defined(SP_DBUG)

!      f21 = TRANSFER(mrc_header(21),f21)
!      write(3,'(A,I5,2X,F8.1,2X,F8.2)')
!     &     '  In lunGetstats_mrc; lun,fmaxt,mrc_header(21):',
!     &        lun,fmaxt,f21

!      write(3,'(A,3I7)')
!     &     '  In lunGetstats_mrc;  mrc_header 17,18,19: ',
!     &         mrc_header(17),mrc_header(18),mrc_header(19)

!      write(3,'(A,I5,2X,3(2x,I8))')
!     &     '  In lunGetstats_mrc; nstack,imgnum,imgstats,imamit:',
!     &                            nstack,imgnum,imgstats,imamit

!      write(3,*)' Leaving lunGetstats_mrc ------'
!      write(3,*)'  '
#endif

      IRTFLG = 0
      END

C     ------------------------- LUNZEROSTATS_MRC ----------------------

      SUBROUTINE LUNZEROSTATS_MRC(LUN,IVERSION,IRTFLG)

#include "LUNHDR.INC"

      INTEGER  :: LUN,IVERSION,IRTFLG
 
      REAL     :: FMINT,FMAXT,FAVT,FSIGT
      INTEGER  :: I4V

C     POINT TO MRC HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET UNDETERMINED IMAGE STATISTICS (FASTER)
      IF (IVERSION >= 20140) THEN
         FMINT =  0    ! UNDETERMINED SHUD BE: FMINT =  0 FOR MRC 2014
         FMAXT = -1    ! UNDETERMINED SHUD BE: FMAXT = -1
         FAVT  = -2    ! UNDETERMINED SHUD BE: FAVT  = -2
         FSIGT = -1    ! UNDETERMINED SHUD BE: FSIGT = -1
      ELSE
         FMINT =  0     
         FMAXT =  0     
         FAVT  =  0    
         FSIGT = -1     
      ENDIF
          
C     SET MRC_HEADER VALUES
      MRC_HEADER(20) = TRANSFER(FMINT,I4V) ! DMAX   MIN  DENSITY VALUE
      MRC_HEADER(21) = TRANSFER(FMAXT,I4V) ! DMIN   MAX  DENSITY VALUE  
      MRC_HEADER(22) = TRANSFER(FAVT, I4V) ! DMEAN  MEAN DENSITY VALUE  
      MRC_HEADER(55) = TRANSFER(FSIGT,I4V) ! RMS    MAP DEVIATION FROM MEAN DENSITY          

C     PUT UNKNOWN STAT'S IMAGE NUMBER IN FILE HEADER         
      CALL LUNSET_STATSIMG_MRC(LUN,0,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

#if defined(SP_BACK)
      CALL BACKTRACE 
#endif
 
#if  defined(SP_DBUG)
      f21 = TRANSFER(mrc_header(21),f21)
!      write(3,'(A,I5,2X,F8.1,2X,F8.2)')
!     &     '  In lunZerostats_mrc; lun,fmax,mrc_header(21):',
!     &        lun,fmax,f21

!      write(3,'(A,3I7)')
!     &     '  In lunZerostats_mrc;  mrc_header 17,18,19: ',
!     &         mrc_header(17),mrc_header(18),mrc_header(19)

!      write(3,*)' Leaving lunZerostats_mrc ------'
!      write(3,*)'  '
#endif



      IRTFLG = 0
      END

C     ------------------------- LUNGETVERSION_MRC ---------------------

       SUBROUTINE LUNGETVERSION_MRC(LUN,IVERSION,IRTFLG)

#include "LUNHDR.INC"

      INTEGER   :: LUN,IVERSION,IRTFLG 

C     POINT TO MRC HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET RETURN VALUE (MUST ALREADY BE SET, USUALLY BY: LUNZEROSTATS_MRC)
      IVERSION = MRC_HEADER(28)

      IRTFLG = 0
      END

C     ------------------------- LUNGETPIXSIZ_MRC ----------------------

      SUBROUTINE LUNGETPIXSIZ_MRC(LUN,PIXSIZ,IRTFLG)

C     PURPOSE:  GET PIXEL SIZE, NOTE THAT MRC HAS 3 SIZES BUT NOT SPIDER
C     LOCATIONS 11-13 CELLA       CELL DIMENSIONS IN ANGSTROMS

#include "LUNHDR.INC"

      REAL              :: PIXSIZ,FVAL
      INTEGER           :: LUN,IRTFLG

      INTEGER           :: NX,NY,NZ

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      CALL LUNGETSIZE_MRC(LUN,NX,NY,NZ,IRTFLG)

C     GET RETURN VALUE
      PIXSIZ = TRANSFER(MRC_HEADER(11),FVAL) 
      PIXSIZ = PIXSIZ / FLOAT(NX) 

      IRTFLG = 0
      END

C     ------------------------- LUNSETPIXSIZE_MRC ---------------------

      SUBROUTINE LUNSETPIXSIZ_MRC(LUN,PIXSIZ,IRTFLG)

C     PURPOSE:  SET PIXEL SIZE, NOTE MRC HAS 3 SIZES BUT NOT SPIDER!
C     LOCATIONS 11-13 CELLAX,CELLAY,CELLAZ   CELL DIMENSIONS IN ANG.

#include "LUNHDR.INC"

      REAL              :: PIXSIZ
      INTEGER           :: LUN,IRTFLG

      INTEGER           :: NX,NY,NZ,IV4
      REAL              :: CELLAX,CELLAY,CELLAZ 
   
C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      CALL LUNGETSIZE_MRC(LUN,NX,NY,NZ,IRTFLG)

      CELLAX = PIXSIZ * NX
      CELLAY = PIXSIZ * NY
      CELLAZ = PIXSIZ * NZ

!      write(3,*) ' In lunsetpixsiz_mrc, nx...: ',nx,ny,nz
!      write(3,*) ' In lunsetpixsiz_mrc, sizex...: ',
!     &                                  sizex,sizey,sizez
!      write(3,*) ' In lunsetpixsiz_mrc, cellax...: ',
!     &                                  cellax,cellay,cellaz


C     SET HEADER VALUES FOR PIXEL SIZE * NO. PIXELS IN CELL
      MRC_HEADER(11) = TRANSFER(CELLAX, IV4)
      MRC_HEADER(12) = TRANSFER(CELLAY, IV4)
      MRC_HEADER(13) = TRANSFER(CELLAZ, IV4)

      IRTFLG = 0
      END

C     ------------------------- LUNSETPIXSIZES_MRC --------------------

      SUBROUTINE LUNSETPIXSIZES_MRC(LUN,SIZEX,SIZEY,SIZEZ,IRTFLG)

C     PURPOSE:  SET PIXEL SIZE, NOTE MRC HAS 3 SIZES BUT NOT SPIDER!
C     LOCATIONS 11-13 CELLAX,CELLAY,CELLAZ   CELL DIMENSIONS IN ANG.

#include "LUNHDR.INC"

      REAL              :: SIZEX,SIZEY,SIZEZ
      INTEGER           :: LUN,IRTFLG

      INTEGER           :: NX,NY,NZ,IV4
      REAL              :: CELLAX,CELLAY,CELLAZ 
   
C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      CALL LUNGETSIZE_MRC(LUN,NX,NY,NZ,IRTFLG)

      CELLAX = SIZEX * NX
      CELLAY = SIZEY * NY
      CELLAZ = SIZEZ * NZ

!      write(3,*) ' In lunsetpixsiz_mrc; nx...: ',nx,ny,nz
!      write(3,*) ' In lunsetpixsiz_mrc; sizex...: ',
!     &                                  sizex,sizey,sizez
!      write(3,*) ' In lunsetpixsiz_mrc; cellax...: ',
!     &                                  cellax,cellay,cellaz


!     SET HEADER VALUES FOR PIXEL SIZE * NO. PIXELS IN CELL
      MRC_HEADER(11) = TRANSFER(CELLAX, IV4)
      MRC_HEADER(12) = TRANSFER(CELLAY, IV4)
      MRC_HEADER(13) = TRANSFER(CELLAZ, IV4)

      IRTFLG = 0
      END


!     ------------------------- LUNSETANGS_MRC ------------------------

      SUBROUTINE LUNSETANGS_MRC(LUN,NANG,ANGS,IRTFLG)   !UNUSED??

C     PURPOSE:  SET ANGLES, NOTE THAT MRC HAS 6, SPIDER HAS 9
C     LOCATIONS 43-48 ANG1...ANG6  USING SPIDER ANGLE CONVENTIONS

#include "LUNHDR.INC"

      INTEGER           :: LUN,NANG,IRTFLG
      REAL              :: ANGS(6)

      INTEGER           :: IV4
   
C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET HEADER VALUES FOR ANGLES
      MRC_HEADER(42) = NANG
      MRC_HEADER(43) = TRANSFER(ANGS(1), I4V)
      MRC_HEADER(44) = TRANSFER(ANGS(2), I4V)
      MRC_HEADER(45) = TRANSFER(ANGS(3), I4V)
      MRC_HEADER(46) = TRANSFER(ANGS(4), I4V)
      MRC_HEADER(47) = TRANSFER(ANGS(5), I4V)
      MRC_HEADER(48) = TRANSFER(ANGS(6), I4V)

      IRTFLG = 0
      END

C     ------------------------- LUNGETANGS_MRC ------------------------

      SUBROUTINE LUNGETANGS_MRC(LUN,NANG,ANGS,IRTFLG)  !! UNUSED

C     PURPOSE:  RETRIEVE ANGLES, NOTE THAT MRC HAS 6 SPIDER HAS 9
C     LOCATIONS 43-48 ANG1...ANG6  USING SPIDER ANGLE CONVENTIONS

#include "LUNHDR.INC"

      REAL              :: ANGS(6)
      INTEGER           :: LUN,NANG,IRTFLG

      REAL              :: FVAL

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     RETRIEVE HEADER VALUES FOR ANGLES
      NANG    = MRC_HEADER(42)
      ANGS(1) = TRANSFER(MRC_HEADER(43),FVAL)
      ANGS(2) = TRANSFER(MRC_HEADER(44),FVAL)
      ANGS(3) = TRANSFER(MRC_HEADER(45),FVAL)
      ANGS(4) = TRANSFER(MRC_HEADER(46),FVAL)
      ANGS(5) = TRANSFER(MRC_HEADER(47),FVAL)
      ANGS(6) = TRANSFER(MRC_HEADER(48),FVAL)

      IRTFLG = 0
      END

C     ------------------------- LUNGETSIZE_MRC ------------------------

      SUBROUTINE LUNGETSIZE_MRC(LUN,NX,NY,NZ,IRTFLG)

#include "LUNHDR.INC"

      INCLUDE 'CMLIMIT.INC'

      INTEGER   :: LUN,NX,NY,NZ,IRTFLG 

      CHARACTER (LEN=MAXNAM)  :: FILNAM
      CHARACTER (LEN=1)       :: DSP
      INTEGER                 :: NLET 

C     POINT TO MRC HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET RETURN VALUES
      NX = MRC_HEADER(1)
      NY = MRC_HEADER(2)
      NZ = MRC_HEADER(3)

      IF (NZ < 0) NZ = -NZ

C     HACK TO HANDLE RELION'S MRCS/mrcs STACK FILE WITHOUT MZ VALUE

C     RETRIEVE CURRENT FILENAME (INCLUDES IMGNUM IF STACKED IMAGE)
      CALL LUNGETFILE_MRC(LUN,FILNAM,NLET,DSP,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
       
      IF (INDEX(FILNAM(1:NLET),'.mrcs') > 0 .OR.
     &    INDEX(FILNAM(1:NLET),'.MRCS') > 0) THEN 
         NZ = 1    ! RELION MRC IMAGE STACK
      ENDIF
      !write(3,*) ' In lunsetmrchdr; filnam: ',filnam(1:nlet), nz

      IRTFLG = 0
      END

C     ------------------------- LUNGETVALS_R_MRC ----------------------

      SUBROUTINE LUNGETVALS_R_MRC(LUN,IGO,NVAL,BUFOUT,IRTFLG)

C     GET RETURN VALUES, ONLY WORKS FOR REALS STORED AS I*4

#include "LUNHDR.INC"

      INCLUDE 'CMLIMIT.INC'    ! PROVIDES IAPLOC
      INTEGER,PARAMETER :: IAPLOC_MRC  = 29

      REAL              :: BUFOUT(NVAL)
      INTEGER           :: LUN,IGO,NVAL,IRTFLG

      REAL              :: FVAL
      INTEGER           :: IEND,I,IGOMRC

      INTEGER           :: JENIS(56)  ! INTEGER=1,REAL=2,CHAR=3

C                1,2,3, 4, 5,6,7, 8,9,10, 11,12,13, 14,15,16,
      JENIS = (/ 1,1,1, 1, 1,1,1, 1,1, 1,  2, 2, 2,  2, 2, 2,   

C                17,18,19, 20,21,22, 23, 24 25,26, 27,28,
     &            1, 1, 1,  2, 2, 2,  1,  1, 1, 3,  3, 1, 

C                29,30,31,32,33,34,35,36,37,38,39,40,41, 42,
     &            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,  1,

C                43,44,45,46,47,48,49, 50,51,52, 53,54, 55, 56
     &            2, 2, 2, 2, 2, 2, 2,  2, 2, 2,  3, 1,  2,  1 /)


C     POINT TO MRC HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     HACK FOR TRANSLATING MRC IAPLOC INSTEAD OF SPIDER LOCATIONS
      IGOMRC = IGO
      IF (IGO == (IAPLOC + 1)) THEN
         IGOMRC = IAPLOC_MRC
      ENDIF
 
      IEND = IGOMRC + NVAL - 1
      IF (IGOMRC < 1 .OR. IEND > 256) THEN
         CALL ERRT(102,'HEADER LOCATION MUST BE 1...256',IEND)
         IRTFLG = 1
         RETURN
      ENDIF

C     GET REAL RETURN VALUES
      DO I = IGOMRC,IEND
         IF (JENIS(I) == 2) THEN
            BUFOUT(I-IGOMRC+1) = TRANSFER(MRC_HEADER(I),FVAL)
         ELSEIF (JENIS(I) == 1) THEN
            BUFOUT(I-IGOMRC+1) = MRC_HEADER(I)
         ELSE
           CALL ERRT(102,'CHARACTER HEADER LOCATION',I)
           IRTFLG = 1
           RETURN
         ENDIF
      ENDDO

      IRTFLG = 0
      END

C     ------------------------- LUNSETVALS_R_MRC ----------------------

      SUBROUTINE LUNSETVALS_R_MRC(LUN,IGO,NVAL,BUFVALS,IRTFLG)

C     SET HEADER VALUES ONLY WORKS FOR REALS STORED AS I*4
#include "LUNHDR.INC"

      INCLUDE 'CMLIMIT.INC'    ! PROVIDES IAPLOC
      INTEGER,PARAMETER :: IAPLOC_MRC  = 29

      REAL              :: BUFVALS(NVAL)
      INTEGER           :: LUN,IGO,NVAL,IRTFLG

      INTEGER           :: IEND,I,IGOMRC
      INTEGER           :: I4V

      INTEGER           :: JENIS(56)
      REAL              :: FVAL

C                1,2,3, 4, 5,6,7, 8,9,10, 11,12,13, 14,15,16,
      JENIS = (/ 1,1,1, 1, 1,1,1, 1,1, 1,  2, 2, 2,  2, 2, 2,   

C                17,18,19, 20,21,22, 23, 24 25,26, 27,28,
     &            1, 1, 1,  2, 2, 2,  2,  2, 2, 2,  2, 2, 

C                29,30,31,32,33,34,35,36,37,38,39,40,41, 42,
     &            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,  1,

C                43,44,45,46,47,48,49, 50,51,52, 53,54, 55, 56
     &            2, 2, 2, 2, 2, 2, 2,  2, 2, 2,  3, 1,  2,  1 /)


C     POINT TO MRC HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     HACK FOR TRANSLATING MRC IAPLOC INSTEAD OF SPIDER LOCATIONS
      IGOMRC = IGO
      IF (IGO == (IAPLOC + 1)) THEN
         IGOMRC = IAPLOC_MRC
      ENDIF
 
      IEND = IGOMRC + NVAL - 1
      IF (IGOMRC < 1 .OR. IEND > 56) THEN
         CALL ERRT(102,'HEADER LOCATION MUST BE 1...56',IEND)
         IRTFLG = 1
         RETURN
      ENDIF

C     SET HEADER VALUES WRITTEN AS I*4
      DO I = IGOMRC,IEND
         !write(3,*) ' Bufvals: ',jenis(i),bufvals(i-igomrc+1)

         IF (JENIS(I) == 2) THEN
            MRC_HEADER(I) = TRANSFER(BUFVALS(I-IGOMRC+1),I4V)   

         ELSEIF (JENIS(I) == 1) THEN
C           CONVERT REAL VALUE TO INTEGER
            MRC_HEADER(I) = BUFVALS(I-IGOMRC+1)   

         ELSE
           CALL ERRT(102,'CHARACTER HEADER LOCATION',I)
           IRTFLG = 1
           RETURN
         ENDIF
      ENDDO

      IRTFLG = 0
      END

C     ------------------------- LUNGETTYPE_MRC ------------------------

      SUBROUTINE LUNGETTYPE_MRC(LUN,ITYPE,IRTFLG)

      IMPLICIT NONE

      INTEGER   :: LUN,ITYPE,IRTFLG
 
      INTEGER   :: NX,NY,NZ

C     NZ IS ADJUSTED FOR RELION .mrcs FILES IF NEEDED
      CALL LUNGETSIZE_MRC(LUN,NX,NY,NZ,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      
C     RETURN FILE TYPE USING SPIDER'S CONVENTION 
      ITYPE = 1
      IF (NZ > 1) ITYPE = 3   ! MRC FOURIER NOT SUPPORTED

      IRTFLG = 0
      END

C     ------------------------- LUNGETSTK_MRC -------------------------

       SUBROUTINE LUNGETSTK_MRC(LUN,MZ,NSTK,ISTK,IRTFLG)

#include "LUNHDR.INC"

      INTEGER   :: LUN,MZ,NSTK,IMGNUM,IRTFLG

C     POINT TO MRC HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET RETURN VALUES
      MZ   = MRC_HEADER(10)    ! CAN BE NSTK 
      ISTK = MRC_HEADER(259)   ! CURRENT IMG/VOL IN HEADER
      NSTK = MRC_HEADER(260)   ! NUMBER OF IMAGES/VOLUMES IN STACK

      IRTFLG = 0
      END

C     ------------------------- LUNSETSTK_MRC--------------------------

      SUBROUTINE LUNSETSTK_MRC(LUN,ISTK,NSTK,IRTFLG)

#include "LUNHDR.INC"

      INTEGER  :: LUN,ISTK,NSTK,IRTFLG

C     POINT TO MRC HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET STACK RELATED STATIC HEADER VALUES
      MRC_HEADER(259) = ISTK  ! CURRENT IMG/VOL IN HEADER
      MRC_HEADER(260) = NSTK  ! NUMBER OF IMAGES/VOLUMES IN STACK
 
      ! write(3,*)' In lunsetstk_mrc; istk,nstk: ',istk,nstk

      IRTFLG = 0
      END

C     ------------------------- LUNGETNSTACK_MRC ----------------------

       SUBROUTINE LUNGETNSTACK_MRC(LUN, NX,NY,NZ, NZ3,NSTK_FLG,IRTFLG)

C      PURPOSE: GET SPIDER'S NSTACK FROM MRC HEADER INFO (NOT FILE)
C               MRC FILES LACK STANDARD ENCODING FOR STACK SIZE
C               AND SOME VOLUMES LOOK LIKE STACKS SO THIS IS COMPLEX!!
C
C      RETURNS:   NX,NY,NZ, NZ3,NSTK_FLG 
C
C      VARIABLES FOR MRC 2014
C                 ISPG (23)   NZ (3)     MZ (10)
C       IMAGE           0      1          1
C       IMAGE STACK     0     NSTK       NSTK
C       VOLUME (Z)      1     NZ         NZ
C       VOLUME STACK  401     NZ*NSTK    NZ/NSTK
       
C      IMPLICIT NONE

#include "LUNHDR.INC"
      INCLUDE 'CMLIMIT.INC'

      INTEGER                :: LUN,NX,NY,NZ, MZ,NSTK_FLG,IRTFLG

      CHARACTER (LEN=MAXNAM) :: FILNAM     
      CHARACTER (LEN=1)      :: DSP
      CHARACTER (LEN=4)      :: cloc27
      INTEGER                :: NZ3,NLET,IVERSION,ISPG

      INTEGER                :: lnblnkn   ! FUNCTION

C     POINT TO MRC HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET HEADER VALUES
      NX       = MRC_HEADER(1) 
      NY       = MRC_HEADER(2) 
      NZ3      = MRC_HEADER(3)   ! NO. OF IMAGES IN RELION mrcs FILES!

      MZ       = MRC_HEADER(10)  ! CAN BE EITHER NSTK OR NZ 
      cloc27   = TRANSFER(MRC_HEADER(27),cloc27)   
      IVERSION = MRC_HEADER(28)
      ISPG     = MRC_HEADER(23)

#if defined(SP_DBUGIO)
C23456
      write(3,*)  '  '
      write(3,*)' In lungetnstack; iversion,ispg: ', iversion,ispg
      write(3,*)' In lungetnstack; nz3,mz:        ', nz3,mz
#endif

C     GET CURRENT FILE NAME FROM LUNFILNAM COMMON
      CALL LUNGETFILE_MRC(LUN,FILNAM,NLET,DSP,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

#if defined(SP_DBUGIO)
C23456
      write(3,*)' In lungetnstack; dsp,nlet:        ', dsp,nlet
      write(3,*)' In lungetnstack; filnam:          ', filnam(1:nlet)
#endif


      IF (IVERSION >= 20140 .AND. ISPG == 0) THEN
C        IMAGE OR IMAGES STACK

         IF (NZ3 == 1 .AND. MZ == 1) THEN ! --------------------
C           SINGLE IMAGE FILE - MRC VERSION 20140
            NZ       = 1
            NSTK_FLG = -2

#if defined(SP_DBUGIO)
            write(3,*)' In lungetnstack 11; nz,nstk_flg:', nz,nstk_flg
#endif

         ELSEIF (NZ3 == 1 .AND. MZ > 1) THEN  ! ------------
C           IMAGE STACK FILE IN - MRC VERSION 20140
            NZ       = 1
            NSTK_FLG = MZ

#if defined(SP_DBUGIO)
            write(3,*)' In lungetnstack 22; nz,nstk_flg:', nz,nstk_flg
#endif
         ENDIF


      ELSEIF (IVERSION >= 20140 .AND. ISPG == 1) THEN  ! ------------
C        SINGLE VOLUME FILE - MRC VERSION 20140
         NZ       = NZ3  ! nz3=348 mz=348  ispg=1 ll 
         NSTK_FLG = -2

#if defined(SP_DBUGIO)
         write(3,*)' In lungetnstack 33; nz,nstk_flg:', nz,nstk_flg
#endif


      ELSEIF (IVERSION >= 20140 .AND. ISPG >= 401 ) THEN ! ----------
C        STACK OF VOLUMES IN FILE - MRC VERSION 20140

         NZ       = MZ
         NSTK_FLG = NZ3 / MZ

#if defined(SP_DBUGIO)
         write(3,*)' In lungetnstack 44; nz,nstk_flg:', nz,nstk_flg
#endif



      ELSEIF (INDEX(FILNAM(1:NLET),'.mrcs') > 0 .OR.
     &        INDEX(FILNAM(1:NLET),'.MRCS') > 0) THEN ! ----------

C        MAY NEED HACK FOR RELION mrcs IMAGE STACKS HAVING NZ > 1
C        IN RELION  THIS WAS A STACK OF IMAGES NOT A VOLUME
C          NZ       = 1        
C          MZ       = 1
C          NSTK_FLG = NZ3      ! WAS NZ IN FILE

         NZ       = 1
         NSTK_FLG = NZ3

#if defined(SP_DBUGIO)
         write(3,*)' In lungetnstack 55; nz,nstk_flg:', nz,nstk_flg
#endif 
 

      ELSEIF (INDEX(FILNAM(1:NLET),'.map') > 0 .OR.
     &        INDEX(FILNAM(1:NLET),'.MAP') > 0) THEN  ! ----------
C        HACK FOR EMD VOLUME 

         NZ       = NZ3      ! 
         NSTK_FLG = -2       ! ???????

#if defined(SP_DBUGIO)
         write(3,*)' In lungetnstack 66; nz,nstk_flg:', nz,nstk_flg
#endif

      ELSEIF (INDEX(FILNAM(1:NLET),'@') >0 .AND. 
     &        NZ3 == 1 .AND. MZ == 1) THEN            ! ----------
C        HACK FOR IMAGE STACKS CONTAINING ONLY ONE IMAGE
         NZ       = 1        ! 
         NSTK_FLG = 1        ! 
 
#if defined(SP_DBUGIO)
         write(3,*)' In lungetnstack 77; nz,nstk_flg:', nz,nstk_flg
#endif

      ELSEIF (MZ > 1 .AND. NZ3 == MZ) THEN            ! ----------
C        SINGLE VOLUME                            
         NZ       = NZ3
         NSTK_FLG = -2     ! NOT STACK, RETURN NEGATIVE (-2)

#if defined(SP_DBUGIO)
         write(3,*)' In lungetnstack 88; nz,nstk_flg:', nz,nstk_flg
#endif

      ELSE                                            ! ----------
C        VOLUME/IMAGE.  NOT A STACK
         NZ       = NZ3
         NSTK_FLG = -2     ! NOT STACK, RETURN NEGATIVE (-2)

#if defined(SP_DBUGIO)
         write(3,*)' In lungetnstack 99; nz,nstk_flg:', nz,nstk_flg
#endif
      ENDIF


#if defined(SP_DBUGIO)
C23456  
      write(3,*)' In lungetnstack 99999; dsp,nlet: ', dsp,nlet
      write(3,*)' In lungetnstack 99999; filnam:   ', filnam(1:nlet)
      write(3,*)' In lungetnstack 99999; iversion,ispg:', iversion,ispg
      write(3,*)' In lungetnstack 99999; loc27:        ', loc27

      write(3,*)' In lungetnstack 99999, nz,mz:    ', nz,mz
      write(3,*)' In lungetnstack 99999, nstk_flg; ', nstk_flg
      write(3,*)  '    '
#endif

      IRTFLG = 0
      END





C     -------------------- LUNSETNSTACK_MRC --------------------------

       SUBROUTINE LUNSETNSTACK_MRC(LUN, NZ,NSTACK,IRTFLG)

C      PURPOSE: UPDATE NSTACK RELATED LOCATIONS IN MRC HEADER. 
C               MRC FILES LACK STANDARD ENCODING FOR STACK SIZE
C               SO THIS IS COMPLEX!!

#include "LUNHDR.INC"
      INCLUDE 'CMLIMIT.INC'    ! FOR MAXNAM

      INTEGER                :: LUN,NZ,NSTACK,IRTFLG

      CHARACTER (LEN=MAXNAM) :: FILNAM     
      CHARACTER (LEN=1)      :: DSP
      INTEGER                :: NZMRC,MZ,NLET,IVERSION,ISPG

C     POINT TO MRC HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IVERSION = MRC_HEADER(28)
      ISPG     = MRC_HEADER(23)

C     GET CURRENT FILE NAME FROM LUNFILNAM COMMON
      CALL LUNGETFILE_MRC(LUN,FILNAM,NLET,DSP,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      !write(3,*)' In lunsetnstack, nlet,filnam:', nlet,filnam(1:nlet)
      !write(3,*)' In lunsetnstack, :', iversion,ispg,nz,nstack

      IF (( (INDEX(FILNAM,'MRCS') > 1  .OR.
     &       INDEX(FILNAM,'mrcs') > 1) .OR.
     &      (NSTACK >= 1 .AND. NZ > 1 ))) THEN
C        A RELION COMPATIBLE IMAGE STACK
         NZMRC  = NSTACK         
         MZ     = 1

      ELSEIF (NSTACK >= 1 .AND. NZ == NSTACK .AND.
     &        IVERSION == 20140 .AND. ISPG == 0) THEN
C        A 2014 STANDARD STACK OF IMAGES
         NZMRC  = NZ                     
         MZ     = NSTACK

      ELSEIF (NSTACK >= 1 .AND. NZ == 1 .AND.
     &        IVERSION == 20140 .AND. ISPG == 0) THEN
C        A 2014 STANDARD STACK OF IMAGES
         NZMRC  = NZ                     
         MZ     = NSTACK

      ELSEIF (NSTACK >= 1  .AND.NZ > 1  .AND.
     &        IVERSION == 20140 .AND. ISPG == 401) THEN
C        A 2014 STANDARD STACK OF EM VOLUMES
         NZMRC  = NZ          
         MZ     = NZ * NSTACK

      ELSEIF (NSTACK < 0 .AND. NZ == 1) THEN
C        A SINGLE IMAGE
         NZMRC  = 1                      
         MZ     = 1
      ELSEIF (NSTACK < 0 .AND. NZ > 1 ) THEN
C        A SINGLE VOLUME
         NZMRC  = NZ        
         MZ     = NZ
      ELSE 
         CALL ERRT(102,'BAD STACK OR VOLUME PARAMETERS',NSTACK)
         IRTFLG = 1
         RETURN
      ENDIF

C     UPDATE HEADER VALUES
      MRC_HEADER(3)  = NZMRC  
      MRC_HEADER(10) = MZ     
 
      IRTFLG = 0
      END

C     ------------------- LUNSETCOMMON_MRC ----------------------------

      SUBROUTINE LUNSETCOMMON_MRC(LUN,IRTFLG)

C     PURPOSE: SET SPIDER COMMON BLOCK AND REGISTER VALUES 

      IMPLICIT NONE
      INCLUDE 'CMBLOCK.INC'

      INTEGER           :: LUN,IRTFLG

      CHARACTER(LEN=12) :: CDUM
      CHARACTER(LEN=1)  :: NULL = CHAR(0)
      INTEGER           :: IVERSION,ISPG,NSYMBT,NLABL
      INTEGER           :: NSTACK,NSLICE,IMGNUM,MZ,NANG
      REAL              :: ANGS(9)

C     INCLUDE 'LABLOCK.INC'
C     LABLOCK HAS BEEN INLINED 
C     LABNEW   COMMON BLOCK FOR FILE HEADER VARIABLES
C              USED IN OPENF, FILGEN, FSTACK, FSTADD, ETC.
C     LABLEN             NUMBER OF FLOATING POINT VARIABLES IN LABEL
C     IANGLE             NUMBER OF ANGLES PRESENT
C     PHI,THETA,PSI      ROTATION
C     XOFF,YOFF,ZOFF     X,Y,Z TRANSLATION
C     SCALE,FLAG
C     KANGLE             NUMBER OF SECONDARY ANGLES PRESENT
C     PHI2,THETA2,PSI2   SECONDARY ANGLES
C     PHI1,THETA1,PSI1   SECONDARY ANGLES
C     HDR_VALS           HEADER LOCATIONS 37...101

      INTEGER    :: LABLEN,IANGLE,KANGLE
      REAL       :: PHI,THETA,PSI,XOFF,YOFF,ZOFF,SCALE,FLAG
      REAL       :: PHI2,THETA2,PSI2,PHI1,THETA1,PSI1
      REAL       :: HDR_VALS(64)

      COMMON/LABNEW/LABLEN,IANGLE,PHI,THETA,PSI,XOFF,YOFF,ZOFF,
     &              SCALE,FLAG,KANGLE,PHI2,THETA2,PSI2,PHI1,
     &              THETA1,PSI1,HDR_VALS


C     RARE ANGLE HEADER VALUES NOT AVAILABLE IN MRC FILES
      CALL LUNGETANGS_MRC(LUN,NANG,ANGS,IRTFLG)  ! UNUSED??
      IANGLE   = 0
      KANGLE   = 0

      PHI      = 0.0
      THETA    = 0.0
      PSI      = 0.0

      PHI1     = 0.0
      THETA1   = 0.0
      PSI1     = 0.0


      PHI2      = 0.0
      THETA2    = 0.0
      PSI2      = 0.0

      IF (NANG > 0) THEN
         IANGLE = 1
         PHI    = ANGS(1)
         THETA  = ANGS(2)
         PSI    = ANGS(3)
      ENDIF

      IF (NANG > 3) THEN
         KANGLE = 1
         PHI1   = ANGS(4)
         THETA1 = ANGS(5)
         PSI1   = ANGS(6)
      ENDIF

C     GET HEADER EXTRA LENGTH = NSYMBT 
      CALL LUNGETVIN_MRC(LUN,IVERSION,ISPG,NSYMBT,NLABL,IRTFLG)
      LABLEN = 1024 + NSYMBT

C     RETRIEVE NSTACK & CURRENT IMGNUM
      CALL LUNGETSTK_MRC(LUN,MZ,NSTACK,IMGNUM,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     RETRIEVE VOLUME/IMAGE SIZE  
      CALL LUNGETSIZE_MRC(LUN,NSAMC,NROWC,NSLICE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET VOLUME/IMAGE STATISTICS: FMIN... IN CMBLOCK COMMON 
      CALL LUNGETSTATS_MRC(LUN,IMAMI,FMIN,FMAX,AV,SIG,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     RECOVER FILE DATE, TIME & TITLE FROM HEADER
C     DATE NOT PASSED IN COMMON ANY MORE AS COMMON IS ONLY 10 CHAR.
      CTIM = ' '
      CTIT = ' '
      CDAT = ' '
 
C     SET RESERVED REGISTER VALUES AS NEEDED (SELDOM USED NOW??)
      CALL REG_SET(1,FLOAT(NSAMC),  NULL, IRTFLG)
      CALL REG_SET(2,FLOAT(NROWC),  NULL, IRTFLG)
      CALL REG_SET(3,FMAX,          NULL, IRTFLG)
      CALL REG_SET(4,FMIN,          NULL, IRTFLG)
      CALL REG_SET(5,AV,            NULL, IRTFLG)
      CALL REG_SET(6,SIG,           NULL, IRTFLG)
      CALL REG_SET(7,FLOAT(NSLICE), NULL, IRTFLG) 
      CALL REG_SET(8,FLOAT(NSTACK), NULL, IRTFLG) 

      IRTFLG = 0
      END

C     ------------------------- LUNGETLABELS_MRC -----------------------

      SUBROUTINE LUNGETLABELS_MRC(LUN,NLABL,LABELS,IRTFLG)

#include "LUNHDR.INC"

      INTEGER            :: LUN,NLABL,IRTFLG 
C     CHARACTER(LEN=800) :: LABELS    ; GFORT COMPILER ERROR SEPT 2025
      CHARACTER(LEN=*)   :: LABELS

      INTEGER            :: INOW,I,IGO,IEND,IFLIP
      INTEGER            :: IGOL,IENDL,IGOH,IENDH,IL

      CHARACTER(LEN=4)   :: CSTR

C     POINT TO MRC HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IRTFLG = 1    ! Default error return Sept 2025

      !write(3,*)' In lungetlabels_mrc, nlabl:',nlabl

      IF (NLABL < 0 .OR. NLABL > 80) THEN
          write(3,*) ' In lungetlabels_mrc, Bad nlabl: ',nlabl
          RETURN  ! Stop if we would overflow
      ENDIF

C     GET THE EXISTING NLABL LABELS

      IGOL  = 1
      IENDL = NLABL * 20

      IGOH  = 57
      IENDH = 56 + NLABL * 20

      IL    = -3

      !write(3,*)' In lungetlabels_mrc, len(labels):',len(labels)
      !write(3,*)' In lungetlabels_mrc, nlabl,igoh,iendh,il: ',
      !&                                 nlabl,igoh,iendh,il 
 
      DO IH = IGOH,IENDH      ! 57...57+2*20

        IL = IL + 4

        !write(3,*)' In lungetlabels_mrc, ih,il:',ih,il
 
        LABELS(IL:IL+3) = TRANSFER(MRC_HEADER(IH),CSTR(1:4))

        !write(3,*)' In lungetlabels_mrc, ih,il,labels: ',
        !&                                ih,il,labels(il:il+3) 
      ENDDO
      IGOL  = 1

      !write(3,*)' In lungetlabels_mrc, labels(1:80): ',
      !&                                labels(1:80) 

      CALL LUNGETFLIP(LUN,IFLIP,IRTFLG)

      IF (IFLIP == 1 ) THEN
c        write(3,*) ' In lungetlabels, iflip: ',iflip

         IEND = NLABL * 80
         CALL REVERSEBYTES(LABELS,IEND,IRTFLG)
      ENDIF

C     STRIP LEADING BLANKS FROM LABELS (UNTIDY HACK)
      DO ILAB = 1,NLABL
         IGO  = (ILAB - 1) * 80 + 1
         IEND = IGO + 80 - 1
         IF (IEND > LEN(LABELS)) RETURN   ! Stop if overflow

         DO I = IGO,IEND
            IF (LABELS(I:I) .NE. ' ') THEN
                LABELS(IGO:IEND) = LABELS(I:IEND)
                EXIT
            ENDIF
         ENDDO
      ENDDO

      IRTFLG = 0
      END

C     ------------------------- LUNSAYINFO_MRC ------------------------

      SUBROUTINE LUNSAYINFO_MRC(LUN,SAYIT,IRTFLG)

      INCLUDE 'CMLIMIT.INC'
      INCLUDE 'CMBLOCK.INC'

      INTEGER              :: LUN,IRTFLG
      LOGICAL              :: SAYIT

      INTEGER              :: ITYPE,NSTACK,NX,NY,NZ,IMGNUM
      INTEGER              :: LENLABEL1,HEDBYT
      INTEGER              :: IVERSION,ISPG,NSYMBT,NLABL
      INTEGER              :: MAXIM,LENT,NLET,LT,MZ,MRCMODE
      INTEGER              :: ICOMM,MYPID,MPIERR

      CHARACTER(LEN=1)     :: NULL = CHAR(0)
      CHARACTER(LEN=MAXNAM):: FILNAM
      CHARACTER(LEN=1)     :: DSP
      CHARACTER(LEN=5)     :: DTP
      CHARACTER(LEN=104)   :: CSTRING
      CHARACTER(LEN=11)    :: TYPEF
      CHARACTER(LEN=80)    :: LABEL1
      CHARACTER(LEN=4)     :: CAXIS
      CHARACTER(LEN=1)     :: MODE_STR

      IRTFLG = 0
      IF (.NOT. SAYIT) RETURN

C     SET ICOMM AND MYPID
      CALL SET_MPI(ICOMM,MYPID,MPIERR) 

C     RETRIEVE SPIDER'S ITYPE  
      CALL LUNGETTYPE_MRC(LUN,ITYPE,IRTFLG)

#if defined(SP_DBUGIO)
C23456      
      !write(3,*) 'In lunsayinfo_mrc   0 | type:',itype
#endif

C     RETRIEVE MRC MODE 
      CALL LUNGETMODE_MRC(LUN,MRCMODE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF     (MRCMODE == 0) THEN
         MODE_STR = '0'
      ELSEIF (MRCMODE == 1) THEN
         MODE_STR = '1'
      ELSEIF (MRCMODE == 6) THEN
         MODE_STR = '6'
      ELSEIF (MRCMODE == 2) THEN
         MODE_STR = '2'
      ENDIF


C     RETRIEVE CURRENT HEADER NSTACK & IMGNUM
      CALL LUNGETSTK_MRC(LUN,MZ,NSTACK,IMGNUM,IRTFLG)

#if defined(SP_DBUGIO)
C23456      
      write(3,*)' In lunsayinfo_mrc  1 | mz,nstack,imgnum: ',
     &                                   mz,nstack,imgnum
#endif

C     RETRIEVE SIZE
      CALL LUNGETNSTACK_MRC(LUN, NX,NY,NZ,MZ,NSTACK,IRTFLG)

#if defined(SP_DBUGIO)
C23456      
      write(3,*)' In lunsayinfo_mrc   2 | nz,mz,nstack:     ',
     &                                    nz,mz,nstack 
#endif

      IF     (ITYPE == 1 .AND. NSTACK > 0) THEN
          TYPEF = 'MRC-' // MODE_STR // ' S2'

      ELSEIF (ITYPE == 3 .AND. NSTACK > 0) THEN
          TYPEF = 'MRC-' // MODE_STR // ' S3'

      ELSEIF (ITYPE == 1) THEN
          TYPEF = 'MRC-' // MODE_STR // ' 2'

      ELSEIF (ITYPE == 3) THEN
          TYPEF = 'MRC-' // MODE_STR // ' 3'
      ENDIF
      LT = lnblnkn(TYPEF)

C     GET DATA ORIGIN & HANDEDNES
      CALL LUNGETHAND_MRC(LUN,CAXIS,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     RETRIEVE CURRENT FILENAME (INCLUDES IMGNUM IF STACKED IMAGE)
      CALL LUNGETFILE_MRC(LUN,FILNAM,NLET,DSP,IRTFLG)

#if defined(SP_DBUGIO)
C23456
      write(3,*)' In lunsayinfo_mrc   3 | nlet,filnam: ',
     &                                    nlet,filnam(1:nlet)

      write(3,*) ' '
#endif
      
C     GET HEADER EXTRA LENGTH = NSYMBT 
      CALL LUNGETVIN_MRC(LUN,IVERSION,ISPG,NSYMBT,NLABL,IRTFLG)
      HEDBYT = 1024 + NSYMBT

C     NO FILE DATE, TIME OR TITLE IN MRC FILE HEADER
C     RECOVER LABEL #1 TO USE FOR  TITLE FROM HEADER INSTEAD
      CALL LUNGETLABELS_MRC(LUN,1,LABEL1,IRTFLG)
      LENLABEL1 = lnblnkn(LABEL1)
     
      DTP = '(OLD)'
      IF (DSP == 'N') DTP = '(NEW)'

C     FOR SPIRE OUTPUT
      IF (USE_SPIRE .AND. DSP == 'N' .AND. FILNAM(1:1) .NE. '_') THEN
         CALL SPIREOUT(FILNAM(:NLET),IRTFLG)

         IF (NSTACK .NE. 0 .AND. IMGNUM == 0 .AND. NZ > 1) THEN
C           OVERALL STACKED VOLUME 
            WRITE(CSTRING,890)TYPEF(:LT),NX,NY,NZ,NSTACK, DTP,HEDBYT

         ELSEIF (NSTACK .NE. 0 .AND. IMGNUM == 0) THEN
C           OVERALL STACKED IMAGE 
            WRITE(CSTRING,900)TYPEF(:LT),NX,NY,   NSTACK,DTP,HEDBYT

         ELSEIF (NSTACK .NE. 0 .AND. IMGNUM > 0 .AND. NZ > 1) THEN
C           STACKED VOLUME 
            WRITE(CSTRING,910)TYPEF(:LT),NX,NY,NZ,IMGNUM, DTP,HEDBYT

         ELSE IF (IMGNUM > 0) THEN
C           STACKED IMAGE
            WRITE(CSTRING,920)TYPEF(:LT),NX,NY,   IMGNUM,DTP

        ELSE IF (NZ > 1) THEN
C           SIMPLE VOLUME
            WRITE(CSTRING,930)TYPEF(:LT),NX,NY,NZ,DTP,HEDBYT

         ELSE
C           SIMPLE IMAGE
            WRITE(CSTRING,940)TYPEF(:LT),NX,NY,   DTP,HEDBYT
         ENDIF

         CALL SPIREOUT(CSTRING,IRTFLG)
      ENDIF     ! END OF USE_SPIRE


#if defined(SP_DBUGIO)
      write(3,*) ' '
      write(3,*)' In lunsayinfo_mrc   4 | imgnum,nstack,mz: ',
     &                                    imgnum,nstack,mz
      write(3,*)' In lunsayinfo_mrc   4 | filnam: ',filnam(:nlet)
      write(3,*) ' '
#endif




      IF (VERBOSE .AND. IFOUND .NE. -4) THEN
C         PRINT FILE OPENING INFORMATION

         LENT = LENLABEL1 + NLET   ! LENGTH OF FILENAME AND LABEL1

         IF (NLET > 0 .AND. LENLABEL1 <= 0) THEN
C           HAS FILENAME BUT NO LABEL1
            IF (MYPID <= 0) THEN
               WRITE(NOUT,*) ' ',FILNAM(:NLET)
            ENDIF

         ELSEIF (LENT > 0 .AND. LENT < 70) THEN
C           HAS FILENAME AND LABEL1 THAT FIT ON ONE LINE
            IF (MYPID <= 0) THEN
               WRITE(NOUT,*) ' ',FILNAM(:NLET),'     /',
     &                       LABEL1(1:LENLABEL1)
            ENDIF

         ELSEIF (LENT > 0) THEN
C           HAS FILENAME AND LABEL1 THAT DO NOT FIT ON SINGLE LINE
            IF (MYPID <= 0 .AND. NLET > 0 ) THEN
               WRITE(NOUT,'(2X,A)') FILNAM(:NLET)
               WRITE(NOUT,88) LABEL1(1:LENLABEL1)
88             FORMAT('     ',A)
            ENDIF
         ENDIF


#if defined(SP_DBUGIO)
         write(3,*)' In lunsayinfo_mrc   5 | imgnum,nstack,mz:',
     &                                       imgnum,nstack,mz
         write(3,*) '  '
#endif



         IF (NSTACK >= 0 .AND. IMGNUM == 0 .AND. NZ > 1) THEN
C           OVERALL STACKED VOLUME FILE
            IF (MYPID <= 0) THEN
               WRITE(NOUT,890)TYPEF(:LT),NX,NY,NZ, NSTACK, DTP,
     &                        HEDBYT,CAXIS
890            FORMAT('  (',A,') ',3(I0,1X),' (.. ',I0,')',2X,A,
     &                ' HEADER BYTES: ',I0,'  AXIS:(',A,')')
            ENDIF
            
         ELSEIF (NSTACK >= 0 .AND. IMGNUM <= 0) THEN
C           OVERALL STACKED IMAGE FILE 
            IF (MYPID <= 0) THEN
               WRITE(NOUT,900)TYPEF(:LT),NX,NY ,NSTACK, DTP,
     &                        HEDBYT,CAXIS(:2)
900            FORMAT('  (',A,') ',2(I0,1X),' (.. ',I0,')',2X,A,
     &                ' HEADER BYTES: ',I0,'  AXIS:(',A,')')
            ENDIF

         ELSEIF (IMGNUM > 0 .AND. NZ > 1) THEN
C           STACKED VOLUME
            IF (MYPID <= 0) THEN
               WRITE(NOUT,910) TYPEF(:LT),NX,NY,NZ,IMGNUM, DTP,
     &                         HEDBYT,CAXIS
910            FORMAT('  (',A,') ',3(I0,1X),' (@',I0,')',2X,A,
     &                ' HEADER BYTES: ',I0,'  AXIS:(',A,')')
            ENDIF

         ELSE IF (IMGNUM == 1 .AND. NSTACK < 0 ) THEN
C           SIMPLE IMAGE

#if defined(SP_DBUGIO)
            write(3,*)' In lunsayinfo_mrc   6 | nz,mz:',
     &                                          nz,mz
            write(3,*)' In lunsayinfo_mrc   6 | imgnum,nstack:',
     &                                          imgnum,nstack
            write(3,*) ' '
#endif

            IF (MYPID <= 0) THEN
               WRITE(NOUT,921)TYPEF(:LT),NX,NY,  DTP,
     &                        HEDBYT,CAXIS(:2)
921            FORMAT('  (',A,') ',2(I0,1X), 2X,A,
     &                ' HEADER BYTES: ',I0,'  AXIS:(',A,')')
            ENDIF



         ELSE IF (IMGNUM > 0) THEN
C           STACKED IMAGE

#if defined(SP_DBUGIO)
            write(3,*)' In lunsayinfo_mrc  7 | nz,mz:',
     &                                         nz,mz
            write(3,*)' In lunsayinfo_mrc  7 | imgnum,nstack:',
     &                                         imgnum,nstack
            write(3,*) ' '
#endif

            IF (MYPID <= 0) THEN
               WRITE(NOUT,920)TYPEF(:LT),NX,NY, IMGNUM, DTP,
     &                        HEDBYT,CAXIS(:2)
920            FORMAT('  (',A,') ',2(I0,1X),' (@',I0,')'2X,A,
     &                ' HEADER BYTES: ',I0,'  AXIS:(',A,')')
            ENDIF


         ELSE IF (NZ > 1) THEN
C           SIMPLE VOLUME
            IF (MYPID <= 0) THEN
               WRITE(NOUT,930)TYPEF(:LT),NX,NY,NZ, DTP,HEDBYT,CAXIS
930            FORMAT('  (',A,') ',3(I0,1X),2X,A, 
     &                ' HEADER BYTES: ',I0,'  AXIS:(',A,')')
            ENDIF

         ELSE
C           SIMPLE IMAGE
            IF (MYPID <= 0) THEN
               WRITE(NOUT,940)TYPEF(:LT),NX,NY,   DTP,HEDBYT,CAXIS(:2)
940            FORMAT('  (',A,') ',2(I0,1X),2X,A,
     &                ' HEADER BYTES: ',I0,'  AXIS:(',A,')')
            ENDIF
         ENDIF
      ENDIF

      IRTFLG = 0
      END
