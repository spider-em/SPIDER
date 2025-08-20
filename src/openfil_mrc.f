C++*********************************************************************
C
C OPENFIL_MRC.F ADAPTED FROM OPENFIL             MAY 2019 ArDean Leith
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
C  OPENFIL_MRC(LUNT,FILNAM_AT,NLET,NX,NY,NZ,NSTACK, ITYPE,DSP,IRTFLG
C
C  PURPOSE:       OPEN NEW OR OLD MRC DATA FILE FOR RANDOM 
C                 ACCESS READING/WRITING.  
C
C  PARAMETERS:
C        LUNT       UNIT NUMBER FOR FILNAM_AT                   (SENT)
C        FILNAM_AT  FULL FILE NAME (STACKS HAVE @)              (SENT)
C        NX,NY,NZ   DIMENSIONS OF IMAGE                 (SENT OR RET.)
C        NSTACK     STACK INDICATOR                      (SENT / RET.)
C                     ON INPUT: 
C                        FOR NEW FILE  0 : NOT A STACK, 
C                                     >0 : STACK
C                     ON OUTPUT: 
C                        <0 :  NOT A STACK 
C               ????     =0 :  BARE STACK (NO IMAGE SPECIFIED)????
C                        >0 :  CURRENT STACK SIZE
C 
C        ITYPE      SPIDER FILE TYPE SPECIFIER.                 (RET.)
C        DSPT       CHARACTER VAR. CONTAINING DISPOSITION       (SENT)
C     
C        IRTFLG     ERROR RETURN FLAG.                          (RET.)
C                     0 : NORMAL RETURN
C                     1 : ERROR RETURN
C  VARIABLES:
C
C  CALL TREE:   
C     OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC -->
C                             --> OPENFIL_N_MRC -->
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE OPENFIL_MRC(LUNT,FILNAM_AT,NLET,NX,NY,NZ,NSTACK,
     &                       ITYPE, DSPT,IRTFLG)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=*)        :: FILNAM_AT
      INTEGER                 :: LUNT,NLET, NX,NY,NZ,NSTACK,ITYPE
      INTEGER                 :: IRTFLG
      CHARACTER(LEN=1)        :: DSP,DSPT

      CHARACTER (LEN=MAXNAM)  :: FILNAM
      INTEGER                 :: LUN,NE,I,MAXIM
      CHARACTER(LEN=2)        :: DISP
      CHARACTER(LEN=4)        :: CAXIS

      INTEGER                 :: MRCMODE
      REAL                    :: SCALEX,SCALEY,SCALEZ
      LOGICAL                 :: ISBARE,EX,ISOPEN
      INTEGER                 :: IMGNUM,LOCAT,LUNOP,LENT 
      LOGICAL                 :: WANTUL,ISMRCSFILE

      LOGICAL                 :: ISDIGI,lnblnkn  ! FUNCTIONS


C     WANT TO OPEN OLD OR NEW MRC FILE FOR STREAM ACCESS
           
      LUN = LUNT
      DSP = DSPT      ! MAY BE ALTERED BELOW

C     ALL MRC STACKS WILL HAVE AN @ IN NAME
      LOCAT = INDEX(FILNAM_AT,'@')
C     PRINT *, "opfil_mrc.f : 94: OPENFIL_MRC: LOCAT=", LOCAT

      IF (LOCAT <= 0) THEN
C        NOT SPECIFIED AS A STACK
         IMGNUM = -1
         NSTACK = -2
         FILNAM = FILNAM_AT      
      ELSE
C        GET STACKED IMAGE NUMBER (ZERO IF BARESTACK)
C        OUTPUTS IMGNUM    0 : BARE STACK
C                         -1 : NOT A STACK   (NO @)
C                         >1 : STACKED IMAGE NUMBER

C        PRINT *, "opfil_mrc.f : 107: OPENFIL_MRC: Calling GETMRCIMGNUM"
         CALL GETMRCIMGNUM(FILNAM_AT,LOCAT, FILNAM,IMGNUM,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         NSTACK = 1
         IF (IMGNUM > 0) NSTACK = IMGNUM
         !write(3,*)' In openfil_mrc - imgnum,nstack:',imgnum,nstack

      ENDIF
C     PRINT *, "opfil_mrc.f : 116: OPENFIL_MRC: NSTACK=", NSTACK
C     PRINT *, "opfil_mrc.f : 117: OPENFIL_MRC: IMGNUM=", IMGNUM
      !write(3,*)' In openfil_mrc, filnam :',filnam
      
C     CREATE A MRC HEADER OBJECT FOR THIS LUN (CAN REUSE EXISTING)
      CALL LUNNEWHED_MRC(LUN,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     PUT FILENAME AND DSP IN STATIC AREA OF THE HEADER OBJECT
      CALL LUNSETFILE_MRC(LUN,FILNAM,DSP,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

c     PUT ISBARE IN STATIC (OFF-FILE) AREA OF THE HEADER OBJECT
      ISBARE = (IMGNUM == 0)
      CALL LUNSETISBARE_MRC(LUN,ISBARE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      EX = .FALSE.
      IF (DSP == 'N' .AND. LOCAT > 0) THEN
C        MAY BE OPENING A NEW IMAGE IN AN EXISTING STACK?
C        SEE IF THE MRC STACK ALREADY EXISTS
         INQUIRE(FILE=FILNAM,EXIST=EX,OPENED=ISOPEN,NUMBER=LUNOP)
      ENDIF

      IF (DSP == 'O' .OR. (DSP == 'N' .AND. EX)) THEN
C        OPEN AN ALREADY EXISTING MRC FILE
C        PRINT *, "opfil_mrc.f : 142: Calling OPENFIL_O_MRC"
         CALL OPENFIL_O_MRC(LUN,FILNAM,NLET,DSP,IMGNUM,
     &                      NX,NY,NZ,NSTACK,  IRTFLG)
C        PRINT *, "opfil_mrc.f : 145: OPENFIL_MRC: NSTACK",NSTACK

C        HACK TO HANDLE .mrcs files with bad NZ (REMOVED 6/30/20) al
    
         !write(3,*)' In openfil_mrc, after openfil_o_mrc irtflg:',irtflg
         !write(3,*)' In openfil_mrc, after openfil_o_mrc filnam:',filnam
         !write(3,*)' In openfil_mrc, after openfil_o_mrc imgnum:',imgnum
         !write(3,*)' In openfil_mrc, nstack: ',nstack
         !write(3,*)' In openfil_mrc - nz,nstack:', nz,nstack

      ELSEIF (DSP == 'N') THEN
C        OPEN AND INITIALIZE A NEW MRC FILE

         SCALEX  = 1.0    ! ?????????
         SCALEY  = 1.0
         SCALEZ  = 1.0
         MRCMODE = 2      ! USE 32 BIT REALS ALWAYS FOR NOW???
         
C        PRINT *, "opfil_mrc.f : 162: Calling OPENFIL_N_MRC"
         CALL OPENFIL_N_MRC(LUN,FILNAM,
     &                      MRCMODE,
     &                      NX,NY,NZ,NSTACK,
     &                      SCALEX,SCALEY,SCALEZ,IRTFLG)
C        PRINT *, "opfil_mrc.f : 167: OPENFIL_MRC: NSTACK",NSTACK
      ENDIF
      IF (IRTFLG .NE. 0) RETURN

C     SET IMGNUM AND NSTACK IN STATIC AREA OF FILE HEADER
      CALL LUNSETSTK_MRC(LUN,IMGNUM,NSTACK,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
     
C     SET AXIS ORIGIN LOCATION & VOLUME HANDEDNESS BEFORE LUNSETPOS_MRC
      CALL WHICH_HAND_MRC(LUN,FILNAM,CAXIS,IRTFLG)

C     SET READ/WRITE FILE OFFSETS FOR THIS IMAGE IN LUN COMMON 
      CALL LUNSETPOS_MRC(LUN,IMGNUM,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF (FCHAR(1:2) .NE.  '31' .OR. FCHAR(4:5) .ne.  'HE') THEN
C        WRITE HEADER INTO FILE TO PRESERVE ANY ALTERED VALUES
         CALL LUNWRTHED_MRC(LUN,IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            LENT = lnblnkn(FILNAM)
            WRITE(NOUT,99) IRTFLG,LUN,FILNAM(:LENT)
99          FORMAT( '  *** ERROR(',I4,') ON UNIT: ',I3,' FILE: ',A)
            RETURN
         ENDIF
      ENDIF


C     GET IMAGE TYPE (2D/3D or FOURIER) (RETURNED TO CALLER)
      CALL LUNGETTYPE_MRC(LUN,ITYPE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     PUT COMMON VALUES INTO COMMON AREA (NOT NEEDED IN FUTURE?)
      CALL LUNSETCOMMON_MRC(LUN,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     WRITE FILE OPENING INFO TO RESULTS/TERMINAL
      CALL LUNSAYINFO_MRC(LUN,.TRUE.,IRTFLG)
      
C     SET FLAG FOR NORMAL RETURN	
      IRTFLG = 0

      !write(6,*)' In openfil_mrc, ret. itype,iform: ',itype,iform
      !write(3,*)' In openfil_mrc, ret. irtflg: ',irtflg
      
      END


C     ------------------- GETMRCIMGNUM ---------------------------------

      SUBROUTINE GETMRCIMGNUM(FILNAM_AT,LOCAT, FILNAM,IMGNUM,IRTFLG)

C     IMGNUM         0  :  BARE STACK
C                   -1  :  NOT A STACK   (NO @)
C                   >1  :  STACKED IMAGE NUMBER

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=*)    :: FILNAM_AT,FILNAM
      INTEGER             :: NLET,NSTACK,IMGNUM,MAXIM,IRTFLG

      INTEGER             :: LOCAT,LOCSLASH,NLEN,NDIGITS,IGO  

      LOGICAL             :: ISDIGI,lnblnkn


C     IF PRESENT, EXTRACT STACK FILENAME AND NUMBER e.g. 056@file.mrcs

      LOCAT = INDEX(FILNAM_AT,'@')
      NLEN  = lnblnkn(FILNAM_AT)
      
      IF (LOCAT <= 0) THEN
         IMGNUM = -1           ! NOT A STACK
         IRTFLG = 0
         RETURN

      ELSEIF (ISDIGI(FILNAM_AT(LOCAT-1:LOCAT-1) )) THEN
 
C        EXTRACT IMAGE NUMBER WITHIN STACK FILE

         LOCSLASH = INDEX(FILNAM_AT(1:LOCAT),'/',BACK=.TRUE.)
         IGO      = 1 
         IF (LOCSLASH > 0) IGO = LOCSLASH + 1

C        PRINT *, "opfil_mrc.f : 249: GETMRCIMGNUM: Calling GETFILENUM"
         CALL GETFILENUM(FILNAM_AT(IGO:LOCAT-1),IMGNUM,
     &                   NDIGITS,.TRUE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         !write(6,*)' In getmrcimgnum - imgnum,ndigits:',imgnum,ndigits
         !write(6,*)' In getmrcimgnum - stack name:',filnam_at(1:nlen)
         !write(6,*)' In getmrcimgnum - igo,locat:',igo,locat

         IF (IMGNUM <= 0) THEN
            WRITE(NOUT,*) 'STACK NAME:',FILNAM_AT(1:NLEN)
            CALL ERRT(102,'STACKS START WITH IMAGE: 1 NOT',IMGNUM)
            IRTFLG = 1
            RETURN
         ENDIF

         FILNAM = FILNAM_AT(1:IGO-1) // FILNAM_AT(LOCAT+1:NLEN)

      ELSE
C        NO NUMBER BEFORE @ INDICATES BARE STACK
         IMGNUM = 0
         IF (LOCAT > 1)  THEN
           FILNAM = FILNAM(1:LOCAT-1) // FILNAM_AT(LOCAT+1:NLEN)
         ELSE
           FILNAM = FILNAM_AT(2:NLEN)
         ENDIF
      ENDIF
      
      NLET  = lnblnkn(FILNAM)
      !write(3,*)' In getmrcimgnum - filnam_at:',filnam_at(1:25)
      !write(6,*)' In getmrcimgnum - filnam_at: ',filnam_at(:nlen)
      !write(6,*)' In getmrcimgnum - filnam   : ',filnam(:nlet),irtflg

      IRTFLG = 0
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
     &              INDEX(FILNAM(1:NLETI),'.mrc') > 0 .OR.
     &              INDEX(FILNAM(1:NLETI),'.MRCS') > 0 .OR.
     &              INDEX(FILNAM(1:NLETI),'.mrcs') > 0 .OR.
     &              INDEX(FILNAM(1:NLETI),'.map') > 0 .OR.
     &              INDEX(FILNAM(1:NLETI),'.MAP') > 0)

       END


