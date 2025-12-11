C++*********************************************************************
C
C OPENFIL_N_MRC.F ADAPTED FROM OPENFIL           MAY 2019 ArDean Leith
C                  DEBUG OUTPUT                  OCT 2025 ArDean Leith
C                  REWRITTEN                     NOV 2025 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2025   Health Research Inc.,                        *
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
C  OPENFIL_N_MRC(LUN,FILNAM,DSP,MRCMODE,  NX,NY,NZ, 
C                NSTK, ISTK, IS_BARE,
C                ITYPE, SCALEX,SCALEY,SCALEZ,IRTFLG)
C
C  PURPOSE: OPEN NEW MRC FILE FOR STREAM ACCESS, LITTLE ENDED
C
C  PARAMETERS:
C       LUN           LOGICAL UNIT NUMBER FOR FILNAM            (SENT)
C       FILNAM        CHAR  ARRAY, CONTAINING FILE NAME NO @    (SENT)
C       DSP                  (SENT)
C.      MRCMODE       1/2/3/6 FOR DIFFERENT BITS/VALUE          (SENT)
C       NX,NY,NZ      DIMENSIONS OF IMAGE/VOL             (SENT & RET)
C       NSTK                                              (SENT & RET)
C       ISTK          (UNUSED)
C       IS_BARE       (UNUSED)
C       ITYPE         (UNUSED)
C       SCALEX        (UNUSED) 
C       SCALEY        (UNUSED)
C       SCALEZ        (UNUSED)
C
C       IRTFLG        ERROR RETURN FLAG.                         (RET)
C                        0 : NORMAL RETURN
C                        1 : ERROR RETURN
C
C  VARIABLES: 
C        ITYPE  (TYPE)  FILE TYPE SPECIFIER. 
C         +1     R     2-D IMAGE
C         +3     R3    3-D VOLUME FILE
C
C        ISPG == 0     IMAGE OR IMAGE STACK
C        ISPG == 1     VOLUMES
C        ISPG == 401   STACK OF EM VOLUMES  if Version = 20140
C
C        DMAX  < DMIN                        MAX & MIN UNDETERMINED
C        DMEAN < (SMALLER OF DMIN and DMAX)  DMEAN     UNDETERMINED
C        RMS   < 0.0                         RMS       UNDETERMINED
C
C  MRC FILE HEADER:
C     MZ identifies the # of sections along Z in 
C     each of the volumes that are in a volume stack
C
C     2014 MRC STANDARD (IVERSION = 20140)    
C                       ISPG(23)  NZ(3)          MZ(10)
C       Single image       0       1              1
C       Image stack        0       1              >= 1
C       Single volume      1      >1 (#SLICES)    = NZ
C       Volume stack      401     >1*#VOLS        = NZ / #VOLS  
C
C     RELION STANDARD -  MRCS     (IVERSION = 1)
C                        ISPG     NZ             MZ
C       Single image       0       1              1
C       Image stack        0       1              >= 1
C       Single volume      0      >1 (#SLICES)    = NZ??
C       No Volume stack    0      >1*#VOLS        = NZ / #VOLS  
C
C  CALL TREE:   
C         
C    OPFILEC        OPFILES_MRC   COPYTOMRC
C       |             |              |
C    OPENFIL_MRC <----' <-------------
C       |           
C       |-> GET_FILNAM_INFO
C       |
C       |-> LUNNEWHED
C       |                       
C       IF (NEW)  
C          OPENFIL_N_MRC 
C          |
C          -> OPSTREAMFIL     
C          -> LUNSET_FILE_MRC
C          -> LUNSET_MODSIZE_MRC 
C          -> LUNSET_VIN_MRC   
C          -> LUNSET_XXX_MRC 
C          -> LUNZERO_STATS_MRC
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE OPENFIL_N_MRC(LUN,FILNAM,DSP,
     &                         MRCMODE,
     &                         NX,NY,NZ, 
     &                         NSTK,ISTK,IS_BARE,
     &                         ITYPE, SCALEX,SCALEY,SCALEZ, IRTFLG)
 
C     PURPOSE: OPEN NEW MRC FILE FOR STREAM ACCESS, LITTLE ENDED

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'
 
      INTEGER                     :: LUN 
      CHARACTER(LEN=*),INTENT(IN) :: FILNAM
      CHARACTER(LEN=1)            :: DSP
      INTEGER                     :: MRCMODE
      INTEGER                     :: NX,NY,NZ, NSTK, ISTK 
      LOGICAL                     :: IS_BARE
      INTEGER                     :: ITYPE,IRTFLG
      REAL                        :: SCALEX,SCALEY,SCALEZ

      LOGICAL                     :: IS_MRCS 
      INTEGER                     :: NLABL,NSYMBT,MACHST,NE,IOFFSET
      INTEGER                     :: MX,MY,MZ,IVERSION
      INTEGER                     :: NSYMBYT,ISPG

      CHARACTER(LEN=1)            :: NULL = CHAR(0)

      INTEGER                     :: lnblnkn  ! FUNCTION

#if defined (SP_DBUGIO)
      write(3,*)'  '
      write(3,*)' In openfil_n_mrc; nx,ny,nz,nstk : ',nx,ny,nz,nstk 
      write(3,*)' In openfil_n_mrc; filnam:         ',trim(filnam)
#endif

C     OPEN NEW MRC FILE FOR STREAM ACCESS, LITTLE ENDED

      IRTFLG = 999   ! DO NOT ECHO OPENING INFO
      CALL OPSTREAMFILE(.FALSE.,FILNAM,NULL,LUN,
     &                  'UNFORMATTED','NL',
     &                  'NOASK',.TRUE.,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GRID SAMPLING ON X, Y, Z (MX & MY ) 
C     FOR EM, NUMBER OF SECTIONS ALONG Z IN A VOLUME STACK
C     SO: MZ = NZ / NUMBER OF VOLUMES 

C     I AM TRYING TO HANDLE ALL COMMON FOULED UP MRC STACKS
      MX  = NX
      MY  = NY
 
C     SEE IF USER WANTS RELION COMPATIBLE IMAGE/STACK
      IS_MRCS = (INDEX(FILNAM,'MRCS') > 1  .OR.
     &           INDEX(FILNAM,'mrcs') > 1)

#if defined (SP_DBUGIO)
      write(3,*)' In openfil_n_mrc, is_mrcs: ',is_mrcs
#endif

C     PUT FILENAME IN STATIC AREA OF THE HEADER OBJECT
      CALL LUNSET_FILE_MRC(LUN,FILNAM,NULL,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF ( IS_MRCS .OR. (NZ > 1 .AND. NSTK > 1) )  THEN
C        HACK TO CREATE RELION COMPATIBLE STACKS
         IVERSION = 1
         ISPG     = 0
         NZ       = MAX(NSTK,1)    ! NUMBER OF IMAGES NOT SLICES!
         MZ       = NZ

#if defined (SP_DBUGIO)
         write(3,*)' In openfil_n_mrc  222, is_mrcs:    ',is_mrcs
         write(3,*)' In openfil_n_mrc  222; nz,mz,nstk: ',nz,mz,nstk
#endif
 
      ELSEIF (NZ == 1 .AND. NSTK == 0) THEN
C        JUST A SINGLE IMAGE
         IVERSION = 20140        
         ISPG     = 0    
         NZ       = 1            ! IMAGE (NOT A VOLUME OR STACK)        
         MZ       = 1            ! IMAGE (NOT A VOLUME)
         NSTK     = 0            ! NOT A STACK

#if defined (SP_DBUGIO)
         write(3,*)' In openfil_n_mrc  333; nz,mz,nstk: ',nz,mz,nstk
#endif

      ELSEIF (NZ == 1 .AND. NSTK > 0) THEN
C        A 2014 STANDARD IMAGE STACK
         IVERSION = 20140        
         ISPG     = 0    
         NZ       = 1        ! IMAGE (NOT A VOLUME) (SINGLE SLICE) 
         MZ       = NSTK     ! STACK SIZE, NUMBER OF IMAGES IN STACK

#if defined (SP_DBUGIO)
         write(3,*)' In openfil_n_mrc  444; nz,mz,nstk: ',nz,mz,nstk
#endif


      ELSEIF (NZ > 1 .AND. NSTK == 0 ) THEN
C        A 2014 SINGLE VOLUME
         IVERSION = 20140
         ISPG     = 1    
         !NZ      = NZ          ! VOLUME 
         NSTK     = 0           ! NOT A STACK

#if defined (SP_DBUGIO)
         write(3,*)' In openfil_n_mrc  555; nz,mz,nstk: ',nz,mz,nstk
#endif

      ELSEIF (NZ > 1 .AND. NSTK > 0 ) THEN
C        A 2014 STACK OF EM VOLUMES  (WHAT SOFTWARE USES VOLUME STACKS?)
         IVERSION = 20140          ! SO NZ IS NOT FOR STACK
         ISPG     = 401            ! ? IS THIS RIGHT??
             
         MZ       = NZ * NSTK      ! NZ OF EACH VOLUME ??
         !NZ      = NZ * NSTK_FLG  ! NUMBER OF 'SLICES'
         !NSTK    = NSTK           ! STACK SIZE  (MAY BE WRONG :)

#if defined (SP_DBUGIO)
         write(3,*)' In openfil_n_mrc  666; nz,mz,nstk: ',nz,mz,nstk
#endif

      ELSE 
         CALL ERRT(102,'BAD STACK OR VOLUME PARAMETERS',NSTK)
         IRTFLG = 1
         RETURN
      ENDIF
 
C     MODE IS ALWAYS 32 BIT REAL FOR NEW FILES (FOR NOW)
      CALL LUNSET_MODE_MRC(LUN,MRCMODE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET SIZE (NZ NOT NECESSARILY Z! FOR MRCS, ETC)
      CALL LUNSET_SIZE_MRC(LUN, NX,NY,NZ,
     &                          MX,MY,MZ, IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

#if defined (SP_DBUGIO)
      write(3,*)'  '
      write(3,*)' End openfil_n_mrc; nstk:     ', nstk
      write(3,*)' End openfil_n_mrc; nx,ny,nz: ', nx,ny,nz
      write(3,*)' End openfil_n_mrc; mx,my,mz: ', mx,my,mz
#endif

C     SET EXTRA HEADER BYTES, NO. OF LABELS USED, VERSION, ETC
      NSYMBT = 0                      ! RESET BY LUNSET_XXX_MRC
      NLABL  = 1                      ! RESET BY LUNSET_XXX_MRC

      CALL LUNSET_VIN_MRC(LUN,IVERSION,ISPG,NSYMBT,NLABL,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN


#if defined (SP_DBUGIO)
      !write(3,*)' End openfil_n; iversion,ispg,nlabl: ',
      ! &                         iversion,ispg,nlabl
#endif

C     SET:  N?START,CELLANGS,AXES,NSYMBYT,MAP,MACHST,LABELS,...
      CALL LUNSET_XXX_MRC(LUN,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET IMAGE STATISTICS: FMIN.... AS UNDETERMINED IN NEW FILE
C     DOES NOT SET FMIN.... IN COMMON BLOCK!
      CALL LUNZERO_STATS_MRC(LUN,IVERSION,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      
C     SET FLAG FOR NORMAL RETURN	
      IRTFLG = 0

      END











