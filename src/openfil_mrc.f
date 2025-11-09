C***********************************************************************
C
C OPENFIL_MRC.F ADAPTED FROM OPENFIL             MAY 2019 ArDean Leith
C               BUGS FIXED lnblnkn               SEP 2025 ArDean Leith 
C               DEBUG OUTPUT                     OCT 2025 ArDean Leith
C               REMOVED  GETMRCIMGNUM            OCT 2025 ArDean Leith
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
C  OPENFIL_MRC(LUNT,FILNAM,NLET,NX,NY,NZ, NSTK_FLG, ITYPE,DSPT,IRTFLG
C
C  PURPOSE:        OPEN NEW OR OLD MRC DATA FILE FOR RANDOM ACCESS(IO)
C                  ONLY USED IN: OPFILEC!
C
C  PARAMETERS:
C        LUNT      UNIT NUMBER FOR FILNAM                       (SENT)
C        FILNAM    FULL FILE NAME                      
C                    (STACKS HAVE @ EXCEPT FOR MRCS)            (SENT)
C
C        NX,NY,NZ  DIMENSIONS OF IMAGE/VOL               (SENT OR RET)
C
C        NSTK_FLG  STACK INDICATOR  (NOT SENT)                   (RET)
C                  ON OUTPUT: 
C                     -2 :  NOT A STACK   (NO @)
C                     -1 :  AST FOR NSTK  (@* or *@)   
C                     =0 :  BARE STACK (NO IMAGE SPECIFIED)
C                     >0 :  CURRENT STACK SIZE
C 
C        ITYPE     SPIDER FILE TYPE SPECIFIER.                   (RET)
C        DSPT      CHAR VAR. CONTAINING DISPOSITION             (SENT)
C     
C        IRTFLG    ERROR RETURN FLAG.                            (RET)
C                    0 : NORMAL RETURN
C                    1 : ERROR RETURN
C  VARIABLES:
C
C  CALL TREE:   
C     OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC -->
C                             --> OPENFIL_N_MRC -->
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE OPENFIL_MRC(LUNT,FILNAM,NLET,NX,NY,NZ, NSTK_FLG,
     &                       ITYPE, DSPT,IRTFLG)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=*)        :: FILNAM
      INTEGER                 :: LUNT,NLET, NX,NY,NZ,NSTK_FLG,ITYPE
      INTEGER                 :: IRTFLG
      CHARACTER(LEN=1)        :: DSP,DSPT

      CHARACTER(LEN=MAXNAM)   :: FIL_NOAT,FIL_DIRS,FIL_BASE
      CHARACTER(LEN=4)        :: FIL_EXT
      INTEGER                 :: LOCAT, LOCAST, IMGNUM,NSTK_MAX
      LOGICAL                 :: IS_MRC,IS_STK,ISBARESTK
 
      INTEGER                 :: LUN,NE,I,MAXIM,ilent
      CHARACTER(LEN=2)        :: DISP
      CHARACTER(LEN=4)        :: CAXIS

      INTEGER                 :: MRCMODE
      REAL                    :: SCALEX,SCALEY,SCALEZ
      LOGICAL                 :: ISBARE,EX,ISOPEN
      INTEGER                 :: LUNOP,LENT,IMGNUMT,MZT
      LOGICAL                 :: WANTUL

      LOGICAL                 :: ISDIGI   ! FUNCTIONS
      INTEGER                 :: lnblnkn  ! FUNCTIONS  


C     WANT TO OPEN OLD OR NEW MRC FILE FOR STREAM ACCESS
           
      LUN = LUNT
      DSP = DSPT      ! MAY BE ALTERED BELOW


#if defined(SP_DBUGIO)
      !ilent = len(filnam)
      !write(3,*)' In openfil_mrc, ilent for filnam: ',ilent
      !ilent = len(fil_noat)
      !write(3,*)' In openfil_mrc, ilent for fil_noat: ',ilent
      if (dsp == 'N') then
         write(3,*)' In openfil_mrc, nx,ny,nz: ',nx,ny,nz
         write(3,*)' '
      endif
#endif

C     GET_FILNAM_NSTK RETURNS NSTK_FLG AND OTHER STACK INFO
C                  -2 :  NOT A STACK   (NO @)
C                  -1 :  AST FOR NSTK  (@* or *@)   
C                   0 :  BARE STACK    (@ OR MRCS WITHOUT @)
C                  >0 :  STACKED IMAGE NUMBER

      CALL GET_FILNAM_NSTK(FILNAM,FIL_NOAT,FIL_DIRS,FIL_BASE,FIL_EXT,
     &                     IS_MRC,LOCAST, IS_STK ,NSTK_FLG, IRTFLG)

#if defined(SP_DBUGIO)
      write(3,*)' In openfil_mrc, get_filnam_nstk returned ----------'
      write(3,*)' In openfil_mrc; got is_mrc,locast: ',is_mrc,locast
      write(3,*)' In openfil_mrc; got is_stk,nstk_flg: ',is_stk,nstk_flg
      ilent = len(fil_noat)
      write(3,*)' In openfil_mrc; got fil_noat: ',fil_noat(:ilent)
      write(3,*) '  '
#endif
      
C     CREATE A MRC HEADER OBJECT FOR THIS LUN (CAN REUSE EXISTING)
      CALL LUNNEWHED_MRC(LUN,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
          
C     PUT ISBARE IN STATIC (OFF-FILE) AREA OF THE HEADER OBJECT
      ISBARE = (NSTK_FLG == 0)
      CALL LUNSETISBARE_MRC(LUN,ISBARE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      EX = .FALSE.
      IF (DSP == 'N') THEN
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

      IF (DSP == 'O' .OR. (DSP == 'N' .AND. EX)) THEN  !-----------
C        OPEN AN ALREADY EXISTING MRC FILE, RETURNS NSTK_FLG

#if defined (SP_DBUGIO)
         write(3,*)' In openfil_mrc; B4 openfil_o_mrc -------'
#endif

         IMGNUM = 0    ! NOT USED
C        NSTK_FLG      IS NOT NEEDED BY OR SENT TO: OPENFIL_O_MRC ???
         CALL OPENFIL_O_MRC(LUN,FIL_NOAT,NLET,DSP,IMGNUM,
     &                     NX,NY,NZ, NSTK_FLG,IRTFLG)

#if defined (SP_DBUGIO)
        !ilent = lnblnkn(fil_noat)
        !write(3,*)' Back openfil_mrc; fil_noat: ',fil_noat(:ilent)
         write(3,*)' Back in openfil_mrc; imgnum, nz:', 
     &                                    imgnum, nz
         write(3,*)' Back in openfil_mrc; nstk_flg,irtflg:', 
     &                                    nstk_flg,irtflg
         write(3,*)' Back in openfil_mrc; Calling lungetstk_mrc '
         write(3,*)'  '
#endif

         CALL LUNGETSTK_MRC(LUN,MZT,NSTK_MAX,IMGNUMT,IRTFLG)
    
#if defined (SP_DBUGIO)
         write(3,*)' Back in openfil_mrc; mzt,imgnum: ',mzt,imgnum
         write(3,*)' Back in openfil_mrc; nstk_max,irtflg: ',
     &                                    nstk_max,irtflg
#endif

         IF (NSTK_MAX > 0) THEN
            NSTK_MAX = MAX(IMGNUM,NSTK_MAX)   ! USED BELOW
         ENDIF

#if defined (SP_DBUGIO)
         write(3,*)' In openfil_mrc; nstk_max: ', nstk_max
         write(3,*)' '
#endif


      ELSEIF (DSP == 'N') THEN  !---------------------------------
C        OPEN AND INITIALIZE A NEW MRC FILE AND IMAGE

         SCALEX   = 1.0    ! ?????????
         SCALEY   = 1.0
         SCALEZ   = 1.0
         MRCMODE  = 2      ! USE 32 BIT REALS ALWAYS FOR NOW???

         IMGNUM   = 1
         IF (NSTK_FLG > 0) IMGNUM = NSTK_FLG
         NSTK_MAX = IMGNUM


#if defined (SP_DBUGIO)
         ilent = lnblnkn(fil_noat)
         write(3,*)' In openfil_mrc, B4 openfil_n_mrc; ilent: ',ilent 

         write(3,*)' In openfil_mrc, B4 openfil_n_mrc; fil_noat: ',
     &                                           fil_noat(:ilent)

         write(3,*)' In openfil_mrc, B4 openfil_n_mrc; nz,nstk_flg: ',
     &                                                 nz,nstk_flg
#endif

         CALL OPENFIL_N_MRC(LUN,FIL_NOAT,
     &                      MRCMODE,
     &                      NX,NY,NZ, NSTK_FLG,
     &                      SCALEX,SCALEY,SCALEZ,IRTFLG)

         NSTK_MAX = MAX(IMGNUM,NSTK_MAX)

#if defined (SP_DBUGIO)
         write(3,*) '      '
         write(3,*)' In openfil_mrc; nstk_flg:      ',nstk_flg 
         write(3,*)' In openfil_mrc; irtflg,imgnum: ',irtflg,imgnum 
         write(3,*)' In openfil_mrc; nstk_max       ', nstk_max 
#endif

      ENDIF     !--------------------------------------------

      IF (IRTFLG .NE. 0) RETURN





#if defined (SP_DBUGIO)
      write(3,*) '     '
      ilent = lnblnkn(filnam)
      write(3,*)' In openfil_mrc, filnam:           ',filnam(:ilent)
      write(3,*)' In openfil_mrc; nstk_flg:         ',nstk_flg 
      write(3,*)' In openfil_mrc; irtflg,imgnum:    ',imgnum 
      write(3,*)' In openfil_mrc; nstk_max          ', nstk_max 
#endif

C     SET IMGNUM AND NSTK_MAX IN STATIC AREA OF FILE HEADER
      CALL LUNSETSTK_MRC(LUN,IMGNUM,NSTK_MAX,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
     
#if defined (SP_DBUGIO)
      write(3,*)' In openfil_mrc;  After lunsetstk_mrc ' 
#endif

C     IF (FCHAR(1:2) .NE.  '31' .AND. FCHAR(4:5) .NE.  'HE') THEN
C     IF (FCHAR(1:2) .NE.  '31') THEN
C     IF (FCHAR(4:5) .NE.  'HE') THEN
C
      IF ((FCHAR(1:2) == '31' .AND. FCHAR(4:5) == 'HE') .OR.
     &    (FCHAR(4:9) == 'FROM M')) THEN 

C        DO NOT WRITE HEADER INTO FILE, JUST SET READ/WRITE POSITION
         CALL LUNSET_HEDPOS_MRC(LUN,IRTFLG)
      
#if defined (SP_DBUGIO)
         write(3,*)  ' '
         write(3,*)' In openfil_mrc; After lunset_hedpos_mrc  '
         write(3,*) '  '
#endif
         IF (IRTFLG .NE. 0) THEN
            LENT = lnblnkn(FILNAM)
            WRITE(NOUT,99) IRTFLG,LUN,FILNAM(:LENT)
99          FORMAT( '  *** ERROR(',I4,') ON UNIT: ',I3,' FILE: ',A)
            RETURN
         ENDIF

      ELSE

C        SET AXIS ORIGIN LOCATION & VOLUME HANDEDNESS BEFORE SETPOS
         CALL WHICH_HAND_MRC(LUN,FIL_NOAT,CAXIS,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

#if defined (SP_DBUGIO)
         write(3,*)  '  '
         write(3,*)' In openfil_mrc; After which_hand_mrc '
         ilent = lnblnkn(fil_noat)
         write(3,*)' In openfil_mrc, fil_noat: ',fil_noat(:ilent)
         write(3,*)' In openfil_mrc, caxis:    ',caxis
#endif

C        SET READ/WRITE FILE OFFSETS FOR THIS IMAGE IN LUN COMMON 
         CALL LUNSETPOS_MRC(LUN,NSTK_FLG,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

#if defined (SP_DBUGIO)
         write(3,*)  ' '
         write(3,*)' In openfil_mrc; After lunsetpos_mrc  '
         write(3,*) '  '
#endif

C        WRITE HEADER INTO FILE TO PRESERVE ANY ALTERED VALUES
         CALL LUNWRTHED_MRC(LUN,IRTFLG)

         IF (IRTFLG .NE. 0) THEN
            LENT = lnblnkn(FILNAM)
            WRITE(NOUT,99) IRTFLG,LUN,FILNAM(:LENT)
            RETURN
         ENDIF
      ENDIF

  
#if defined (SP_DBUGIO)
      write(3,*)' In openfil_mrc; After lunwrthed_mrc ' 
      write(3,*)' ' 
#endif

C     GET IMAGE TYPE (2D/3D or FOURIER) (RETURNED TO CALLER)
      CALL LUNGETTYPE_MRC(LUN,ITYPE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     PUT COMMON VALUES INTO COMMON AREA (NOT NEEDED IN FUTURE?)
      CALL LUNSETCOMMON_MRC(LUN,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     WRITE FILE OPENING INFO TO RESULTS/TERMINAL

      !!!!!!!!!!!NSTK_FLG = IMGNUM

#if defined (SP_DBUGIO)
      write(3,*)' In openfil_mrc; imgnum,nstk_flg: ',
     &                            imgnum,nstk_flg 
      write(3,*)' In openfil_mrc; Call lunsayinfo_mrc; lun: ',lun
      write(3,*)' '
#endif

      CALL LUNSAYINFO_MRC(LUN,.TRUE.,IRTFLG)
      
C     SET FLAG FOR NORMAL RETURN	
      IRTFLG = 0

#if defined (SP_DBUGIO)
       write(3,*)'    '
       write(3,*)' LEAVING openfil_mrc; imgnum,nstk_flg : ',
     &                                  imgnum,nstk_flg  
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


