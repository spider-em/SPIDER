
C***********************************************************************
C
C OPENFIL_MRC F ADAPTED FROM OPENFIL            MAY 2019 ArDean Leith
C               BUGS FIXED lnblnkn              SEP 2025 ArDean Leith 
C               DEBUG OUTPUT                    OCT 2025 ArDean Leith
C               REWRITTEN                       NOV 2025 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2025  Health Research InC ,                         *
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
C  OPENFIL_MRC(LUNT,FILPAT, DSP,   NX,NY,NZ, NSTK, 
C              ISTK,IS_BARE, ITYPE, IRTFLG)
C
C  PURPOSE:        OPEN NEW OR OLD MRC DATA FILE FOR RANDOM ACCESS(IO)
C
C  PARAMETERS:
C        LUNT      UNIT NUMBER FOR FILNAM                       (SENT)
C
C        FILPAT    FULL FILE NAME (MAY HAVE @ OR *)                      C                    (STACKS HAVE @ EXCEPT FOR MRCS)            (SENT)
C
C        DSP       CHAR VAR:  O/R/N FOR OLD/READ/NEW            (SENT)
C
C        NX,NY,NZ  DIMENSIONS OF IMAGE/VOL               (SENT OR RET)
C
C        NSTK      STACK SIZE.                           (SENT OR RET)
C
C        ISTK      STACK IMG/VOL NUMBER                  (SENT OR RET)
C
C        IS_BARE   BARE STACK.                           (SENT OR RET)   C
C        ITYPE     SPIDER FILE TYPE SPECIFIER                    (RET)
C
C        IRTFLG    ERROR RETURN FLAG.                            (RET)
C                    0 : NORMAL RETURN
C                    1 : ERROR RETURN
C  VARIABLES:
C
C  OPENFIL_MRC  CALL TREE: 
C         
C    OPFILEC        OPFILES_MRC   COPYTOMRC
C       |             |              |
C    OPENFIL_MRC <----' <-------------
C       |           
C       |-> GET_FILNAM_INFO
C       |
C       |-> LUNNEWHED
C       |                       
C       IF (OLD OR NEW ISTK) 
C              -> OPENFIL_O_MRC 
C                    -> LOADHED_MRC        
C                    -> LUNGET_MAP_MRC  
C                    -> LUNSETFLIP_MRC 
C                    -> LOADHED_MRC        
C                    -> LUNGET_MAP_MRC  
C                    -> LUNSETFLIP_MRC 
C                    -> LUNSET_FILE_MRC    
C                    -> LUNGET_MODE_MRC 
C                    -> LUNGET_VIN_MRC 
C                    -> LUNGET_LABELS_MRC  
C                    -> LUNGET_MODSIZES_MRC 
C                    -> LUNGET_2014_MRC
C                    -> LUNSET_2014_MRC      (NSTK MAY INCREASE)
C                    -> LUNGET_STATS_MRC
C                IF (DSP/=R) -> LUNSET_STATSIMG_MRC 
C
C       IF (NEW)  
C              -> OPENFIL_N_MRC 
C                    -> OPSTREAMFIL     
C                    -> LUNSET_FILE_MRC
C                    -> LUNSET_MODSIZE_MRC 
C                    -> LUNSET_VIN_MRC   
C                    -> LUNSET_XXX_MRC 
C                    -> LUNZERO_STATS_MRC
C
C       |-> LUNSET_STK_260_MRC
C       |-> LUNSET_ISBARE_MRC
C       |-> WHICH_HAND_MRC
C        IF (DSP == R) 
C                -> LUNSET_HEDPOS_MRC
C        ELSE
C                -> LUNSET_POS_MRC
C                -> LUNWRTHED_MRC
C       |-> LUNGET_TYPE_MRC
C       |-> LUNSET_COMMON_MRC
C       |-> LUNSAYINFO_MRC
C
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE OPENFIL_MRC(LUNT,FILPAT,DSP, 
     &                        NX,NY,NZ, 
     &                        NSTK, ISTK, IS_BARE,
     &                        ITYPE, IRTFLG)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=*)        :: FILPAT
      CHARACTER(LEN=1)        :: DSP
      INTEGER                 :: LUNT, NX,NY,NZ, NSTK, LIST_VAL 
      LOGICAL                 :: IS_BARE
      INTEGER                 :: ITYPE,IRTFLG

      LOGICAL                 :: IS_STK
 
      INTEGER                 :: LUN,NE,I,ilent
      CHARACTER(LEN=4)        :: CAXIS

      INTEGER                 :: MRCMODE,ISTK
      REAL                    :: SCALEX,SCALEY,SCALEZ
      LOGICAL                 :: ISBARE,EX,ISOPEN
      INTEGER                 :: LUNOP,LENT,ISTKT,MZT
      LOGICAL                 :: WANTUL

      LOGICAL                 :: ISDIGI   ! FUNCTIONS
      INTEGER                 :: lnblnkn  ! FUNCTIONS  

      CHARACTER(LEN=MAXNAM)   :: FILNAM
      CHARACTER(LEN=MAXNAM)   :: FIL_NOAT,FIL_DIRS,FIL_BASE,FIL_EXT 
      LOGICAL                 :: IS_MRC, IS_MRCS


C     WANT TO OPEN OLD OR NEW MRC FILE FOR STREAM ACCESS
           
      LUN = LUNT

#if defined(SP_DBUGIO)
       write(3,*)' In openfil_mrc, filpat:   ',trim(filpat)
       write(3,*)' In openfil_mrc, list_val: ',list_val
       write(3,*)' '
#endif
      CALL GET_FILNAM_INFO(FILPAT, LIST_VAL, 
     &                     FIL_NOAT,FIL_DIRS,FIL_BASE,FIL_EXT,
     &                     IS_MRC, IS_MRCS, IS_BARE, 
     &                     ISTK,   IRTFLG)

      FILNAM = FIL_NOAT

#if defined(SP_DBUGIO)
       write(3,*)' In openfil_mrc, filnam:   ',trim(filnam)
      !write(3,*)' In openfil_mrc, nx,ny,nz: ',nx,ny,nz
      !write(3,*)' In openfil_mrc, lun:      ',lun
       write(3,*)' In openfil_mrc, istk:     ',istk
       write(3,*)' In openfil_mrc, list_val: ',list_val
       write(3,*)' '
#endif
     
C     CREATE A MRC HEADER OBJECT FOR THIS LUN (CAN REUSE EXISTING)
      CALL LUNNEWHED_MRC(LUN,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
          
      EX = .FALSE.
      IF (DSP == 'N'  .AND. FILNAM(1:1) .NE. '_') THEN
C        MAY BE OPENING A NEW IMAGE IN AN EXISTING STACK?
C        SEE IF THE MRC STACK ALREADY EXISTS
         INQUIRE(FILE=FILNAM,EXIST=EX,OPENED=ISOPEN,NUMBER=LUNOP)

         IF (ISOPEN .AND. LUNOP .NE. LUN) THEN
            LENT = lnblnkn(FILNAM)
            WRITE(NOUT,98) LUN,FILNAM(:LENT),LUNOP
98          FORMAT('  *** ERROR FILE: ',A, 
     &             ' ALREADY OPEN ON UNIT: ',I3)
            IRTFLG = 1
            RETURN
         ENDIF
      ENDIF




      IF (DSP == 'O' .OR. DSP == 'R' .OR. 
     &   (DSP == 'N' .AND. EX)) THEN  ! ---------------------- OLD
C        OPEN AN ALREADY EXISTING MRC FILE, RETURNS: ISTK, NSTK_OLD

#if defined (SP_DBUGIO)
         write(3,*)' In openfil_mrc; B4 openfil_o_mrc -------'
         write(3,*)' In openfil_mrc, filnam:   ', trim(filnam)
         write(3,*)' In openfil_mrc, dsp:      ', dsp
         write(3,*)' In openfil_mrc, istk:     ', istk
         write(3,*)
#endif

         CALL OPENFIL_O_MRC(LUN,FILNAM,DSP,
     &                     NX,NY,NZ,
     &                     NSTK,ISTK, IS_BARE, ITYPE,IRTFLG)
 
#if defined (SP_DBUGIO)
         write(3,*)'  '
         write(3,*)' Back in openfil_mrc; filnam:      ', trim(filnam)
         write(3,*)' Back in openfil_mrc; nz:          ', nz
         write(3,*)' Back in openfil_mrc; nstk,istk:   ', nstk,istk 
         write(3,*)' Back in openfil_mrc; is_bare:     ', is_bare
         write(3,*)'  '
#endif

         NSTK = MAX(ISTK,NSTK)   ! USED BELOW
        
#if defined (SP_DBUGIO)
         write(3,*)' In openfil_mrc; nstk_new: ', nstk 
         write(3,*)' ' 
#endif


      ELSEIF (DSP == 'N') THEN  !----------------------------- NEW
C        OPEN AND INITIALIZE A NEW MRC FILE AND IMAGE

         SCALEX   = 1.0    ! ?????????
         SCALEY   = 1.0
         SCALEZ   = 1.0
         MRCMODE  = 2      ! USE 32 BIT REALS ALWAYS FOR NOW???

         NSTK     = ISTK
         IF (IS_BARE .AND. NSTK < 1) THEN
             NSTK = 1
             ISTK = 1
         ENDIF

         CALL OPENFIL_N_MRC(LUN,FILNAM,DSP,
     &                      MRCMODE,
     &                      NX,NY,NZ, 
     &                      NSTK,ISTK, IS_BARE, ITYPE,
     &                      SCALEX,SCALEY,SCALEZ,IRTFLG)

      ENDIF     !------------------------------------------- OLD/NEW

#if defined (SP_DBUGIO)
      write(3,*) '     '
      write(3,*)' In openfil_mrc; irtflg:    ', irtflg
      write(3,*)' In openfil_mrc, filnam:    ', trim(filnam)
      write(3,*)' In openfil_mrc; is_bare:   ', is_bare 
      write(3,*)' In openfil_mrc; istk,nstk: ', istk,nstk 
#endif

      IF (IRTFLG .NE. 0) RETURN


C     SET ISTK AND NSTK IN STATIC AREA OF FILE HEADER
      CALL LUNSET_STK_260_MRC(LUN,ISTK,NSTK,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      CALL LUNSET_ISBARE_MRC(LUN,ISBARE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN


C     SET AXIS ORIGIN LOCATION & VOLUME HANDEDNESS BEFORE SETPOS
      CALL WHICH_HAND_MRC(LUN,FILNAM,CAXIS,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
     
#if defined (SP_DBUGIO)
      write(3,*)  '  '
      write(3,*)' In openfil_mrc; After lunsetstk_260  ' 
      write(3,*)' In openfil_mrc; After which_hand_mrc '
      write(3,*)' In openfil_mrc, filnam:   ',trim(filnam)
      write(3,*)' In openfil_mrc, caxis:    ',caxis
      write(3,*)' In openfil_mrc, nstk:     ',nstk
      write(3,*)' In openfil_mrc; dsp:                ',dsp 
      write(3,*)  ' '
#endif

      IF (DSP == 'R') THEN    ! NON-WRITEABLE EXISTING FILE --------

C        DO NOT WRITE HEADER INTO FILE, JUST SET READ/WRITE POSITION
         CALL LUNSET_HEDPOS_MRC(LUN,IRTFLG)
      
         IF (IRTFLG .NE. 0) THEN
            LENT = lnblnkn(FILNAM)
            WRITE(NOUT,99) IRTFLG,LUN,FILNAM(:LENT)
99          FORMAT( '  *** ERROR(',I4,') ON UNIT: ',I3,' FILE: ',A)
            RETURN
         ENDIF


      ELSE                 ! WRITEABLE EXISTING OR NEW FILE --------


#if defined (SP_DBUGIO)
      write(3,*)  '  '
      write(3,*)' In openfil_mrc; calling lunset_pos_mrc ' 
      write(3,*)' In openfil_mrc, istk,nstk:   ', istk,nstk
      write(3,*)  ' '
#endif

C        SET READ/WRITE FILE OFFSETS FOR THIS IMAGE IN LUN COMMON 
         CALL LUNSET_POS_MRC(LUN,ISTK,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

#if defined (SP_DBUGIO)
         write(3,*)' In openfil_mrc; calling lunwrthed: '
#endif

C        WRITE HEADER INTO FILE TO PRESERVE ANY ALTERED VALUES
         CALL LUNWRTHED_MRC(LUN,IRTFLG)

         IF (IRTFLG .NE. 0) THEN
            LENT = lnblnkn(FILNAM)
            WRITE(NOUT,99) IRTFLG,LUN,FILNAM(:LENT)
            RETURN
         ENDIF

#if defined (SP_DBUGIO)
         write(3,*)' In openfil_mrc; After lunwrthed_mrc '
         write(3,*)' ' 
#endif

      ENDIF              ! BOTH OLD AND NEW FILES ------------------


C     GET IMAGE TYPE (2D/3D or FOURIER) (RETURNED TO CALLER)
      CALL LUNGET_TYPE_MRC(LUN,ITYPE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     PUT COMMON VALUES INTO COMMON AREA (NOT NEEDED IN FUTURE?)
      CALL LUNSET_COMMON_MRC(LUN,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     WRITE FILE OPENING INFO TO RESULTS/TERMINAL

#if defined (SP_DBUGIO)
      write(3,*)' In openfil_mrc; istk,nstk: ', istk,nstk 
      write(3,*)' In openfil_mrc; Calling lunsayinfo_mrc '
      write(3,*)' '
#endif

      CALL LUNSAYINFO_MRC(LUN,.TRUE.,IRTFLG)
      
C     SET FLAG FOR NORMAL RETURN	
      IRTFLG = 0

#if defined (SP_DBUGIO)
       write(3,*)'    '
       write(3,*)' Leaving openfil_mrc; istk,nstk : ', istk,nstk 
       write(3,*)'    '
#endif
      
      END



C     ------------------------- ISMRCFILE  ----------------------------

 
       LOGICAL FUNCTION ISMRCFILE(FILNAM)

       IMPLICIT NONE

       CHARACTER(LEN=*) :: FILNAM

       INTEGER          :: NLETI
       INTEGER          :: lnblnk

       NLETI = lnblnk(FILNAM)

C      IF IT IHAS .MRC, .MRCS, OR .MAP IT IS LIKELY TO BE MRC FILE 
       ISMRCFILE = (INDEX(FILNAM(1:NLETI),'.MRC') > 0 .OR.
     &              INDEX(FILNAM(1:NLETI),'.mrc') > 0  .OR.
     &              INDEX(FILNAM(1:NLETI),'.map') > 0 .OR.
     &              INDEX(FILNAM(1:NLETI),'.MAP') > 0)

       END


