C++*********************************************************************
C
C OPENFIL_N_MRC.F ADAPTED FROM OPENFIL           MAY 2019 ArDean Leith
C                 DEBUG OUTPUT                   OCT 2025 ArDean Leith
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
C  OPENFIL_N_MRC(LUN,FILNAM,DATEXT, MRCMODE, NX,NY,NZ, NSTK_FLG,
C                SCALEX,SCALEY,SCALEZ,IRTFLG
C
C  PURPOSE: OPEN NEW MRC FILE FOR STREAM ACCESS, LITTLE ENDED
C
C  PARAMETERS:
C       LUN           LOGICAL UNIT NUMBER FOR FILNAM            (SENT)
C       FILNAM        CHAR  ARRAY, CONTAINING FILE NAME NO @    (SENT)
C       LUN           LOGICAL UNIT NUMBER FOR FILNAM            (SENT)
C       NX,NY,NZ      DIMENSIONS OF IMAGE/VOL             (SENT & RET)
C
C       NSTK_FLG      STACK INDICATOR                     (SENT / RET)
C                     ON INPUT: 
C                       -2 :  AST FOR NSTK  (@* or *@)   
C                       -1 :  NOT A STACK   (NO @)
C                        0 :  BARE STACK    (@)
C                       >1 :  STACKED MRC IMAGE/VOLUME? NUMBER
C
C                     ON OUTPUT: 
C                       -1 :  NOT A STACK 
C                      >=0 :  MAX IMAGE NUMBER NOW IN STACK
C 
C       IRTFLG        ERROR RETURN FLAG.                         (RET)
C                       0 : NORMAL RETURN
C                       1 : ERROR RETURN
C
C  VARIABLES: 
C        ITYPE     (TYPE)  FILE TYPE SPECIFIER. 
C             +1     R     2-D IMAGE
C             +3     R3    3-D VOLUME FILE
C
C        ISPG == 0         IMAGE OR IMAGE STACK
C        ISPG == 1         VOLUMES
C        ISPG == 401       STACK OF EM VOLUMES  if Version = 20140
C
C        MZ   ==  1        IMAGE 
C        MZ   >=  1        IMAGE STACK
C        MZ   == 1            IMAGE STACK
C        MZ   == NZ        VOLUME
C        MZ   NZ/NUMVOLS   VOLUME STACK
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
C
C  CALL TREE:   
C     OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC -->
C                             --> OPENFIL_N_MRC -->
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE OPENFIL_N_MRC(LUN,FILNAM,
     &                         MRCMODE,
     &                         NX,NY,NZ, NSTK_FLG,
     &                         SCALEX,SCALEY,SCALEZ,IRTFLG)
 
C     PURPOSE: OPEN NEW MRC FILE FOR STREAM ACCESS, LITTLE ENDED

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'
 
      INTEGER                 :: LUN,MRCMODE

      CHARACTER(LEN=*),INTENT(IN) :: FILNAM

      INTEGER                 :: NX,NY,NZ, NSTK_FLG
      REAL                    :: SCALEX,SCALEY,SCALEZ
      INTEGER                 :: IRTFLG

      LOGICAL                 :: WANTUL,IS_MRCS 
      INTEGER                 :: NLABL,NSYMBT,MACHST,NE,IOFFSET
      INTEGER                 :: MX,MY,MZ,IVERSION
      INTEGER                 :: NSYMBYT,ISPG,LNBLNKN,ILENT

      CHARACTER(LEN=1)        :: NULL = CHAR(0)


      ILENT = LNBLNKN(FILNAM)

#if defined (SP_DBUGIO)
      write(3,*)'  '
      write(3,*)' In openfil_n_mrc; nx,ny,nz,nstk_flg: ', 
     &                              nx,ny,nz,nstk_flg
      !write(3,*)' In openfil_n_mrc; ilent:  ',ilent
       write(3,*)' In openfil_n_mrc; filnam: ',filnam(1:ilent)
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
      CALL LUNSETFILE_MRC(LUN,FILNAM,NULL,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF ( IS_MRCS .OR. (NZ > 1 .AND. NSTK_FLG > 1) )  THEN
C        HACK TO CREATE RELION COMPATIBLE STACKS
         IVERSION = 1
         ISPG     = 0
         NZ       = MAX(NSTK_FLG,1)    ! NUMBER OF IMAGES NOT SLICES!
         MZ       = NZ

#if defined (SP_DBUGIO)
         write(3,*)' In openfil_n_mrc, is_mrcs,nz: ',is_mrcs,nz
#endif
 
      ELSEIF (NZ == 1 .AND. NSTK_FLG < 0) THEN
C        JUST A SINGLE IMAGE
         IVERSION = 20140        
         ISPG     = 0    
         NZ       = 1            ! IMAGE (NOT A VOLUME OR STACK)        
         MZ       = 1            ! IMAGE (NOT A VOLUME)
         NSTK_FLG = -2           ! NOT A STACK


#if defined (SP_DBUGIO)
         write(3,*)' In openfil_n_mrc, nz,nstk_flg: ',nz,nstk_flg
#endif

      ELSEIF (NZ == 1 .AND. NSTK_FLG >= 1) THEN
C        A 2014 STANDARD IMAGE STACK
         IVERSION = 20140        
         ISPG     = 0    
         NZ       = 1            ! IMAGE (NOT A VOLUME) (SINGLE SLICE) 
         MZ       = NSTK_FLG     ! STACK SIZE, NUMBER OF IMAGES IN STACK

      ELSEIF (NZ > 1 .AND. NSTK_FLG < 0 ) THEN
C        A 2014 SINGLE VOLUME
         IVERSION = 20140
         ISPG     = 1    
C        MZ       = MZ           ! VOLUME 
         NSTK_FLG = -1           ! NOT A STACK

      ELSEIF (NZ > 1 .AND. NSTK_FLG > 0 ) THEN
C        A 2014 STACK OF EM VOLUMES  (WHAT SOFTWARE USES VOLUME STACKS?)
         IVERSION = 20140        ! SO NZ IS NOT FOR STACK
         ISPG     = 401    
C        NZ       = NZ           ! NUMBER OF VOLUMES
         MZ       = NZ * NSTK_FLG
         NSTK_FLG = MZ           ! STACK SIZE  (MASY BE WRONG)

      ELSE 
         CALL ERRT(102,'BAD STACK OR VOLUME PARAMETERS',NSTK_FLG)
         IRTFLG = 1
         RETURN
      ENDIF
 
C     MODE IS ALWAYS 32 BIT REAL FOR NEW FILES (FOR NOW)
      CALL LUNSETMODSIZ_MRC(LUN,MRCMODE,NX,NY,NZ,
     &                                  MX,MY,MZ,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

#if defined (SP_DBUGIO)
      write(3,*)' End openfil_n_mrc; nstk_flg: ', nstk_flg
      write(3,*)' end openfil_n_mrc; nx,ny,NZ: ', nx,ny,nz
      write(3,*)' End openfil_n_mrc; mx,my,mz: ', mx,my,mz
#endif


C     SET EXTRA HEADER BYTES, NO. OF LABELS USED, VERSION, ETC
      NSYMBT = 0                      ! RESET BY LUNSETXXX_MRC
      NLABL  = 1                      ! RESET BY LUNSETXXX_MRC

      CALL LUNSETVIN_MRC(LUN,IVERSION,ISPG,NSYMBT,NLABL,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

#if defined (SP_DBUGIO)
      !write(3,*)' End openfil_n; iversion,ispg,nlabl: ',
      ! &                         iversion,ispg,nlabl
#endif

C     SET:  N?START,CELLANGS,AXES,NSYMBYT,MAP,MACHST,LABELS,...
      CALL LUNSETXXX_MRC(LUN,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET IMAGE STATISTICS: FMIN.... AS UNDETERMINED IN NEW FILE
C     DOES NOT SET FMIN.... IN COMMON BLOCK!
      CALL LUNZEROSTATS_MRC(LUN,IVERSION,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      
C     SET FLAG FOR NORMAL RETURN	
      IRTFLG = 0

      END







