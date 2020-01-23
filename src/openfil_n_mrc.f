C++*********************************************************************
C
C OPENFIL_N_MRC.F ADAPTED FROM OPENFIL           MAY 2019 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2020   Health Research Inc.,                         *
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
C  OPENFIL_N_MRC(LUN,FILNAM,DATEXT, MRCMODE, NX,NY,NZ,NSTACK,
C                SCALEX,SCALEY,SCALEZ,IRTFLG
C
C  PURPOSE: OPEN NEW MRC FILE FOR STREAM ACCESS, LITTLE ENDED
C
C  PARAMETERS:
C       LUN       LOGICAL UNIT NUMBER FOR FILNAM.               (SENT)
C       FILNM     CHARACTER ARRAY, CONTAINING FILE NAME WITH @  (SENT)
C       FILNAM    CHARACTER ARRAY, CONTAINING FILE NAME NO @    (SENT)
C       LUN       LOGICAL UNIT NUMBER FOR FILNAM.               (SENT)
C       NX,NY,NZ  DIMENSIONS OF IMAGE                   (SENT OR RET.)
C       NSTACK    STACK INDICATOR                        (SENT / RET.)
C                   ON INPUT: 
C                     -2 : NOT A STACK  
C                     >0 : MRC STACK
C                   ON OUTPUT: 
C                     -1  :  NOT A STACK 
C                     >=0 :  MAX IMAGE NUMBER NOW IN STACK.
C 
C       IRTFLG    ERROR RETURN FLAG.                            (RET.)
C                   0 : NORMAL RETURN
C                   1 : ERROR RETURN
C
C  VARIABLES: ITYPE (TYPE)  FILE TYPE SPECIFIER. 
C               +1    R     2-D IMAGE
C               +3    R3    3-D VOLUME FILE
C
C        ISPG == 0        IMAGE OR IMAGE STACK
C        ISPG == 1        VOLUMES
C        ISPG == 401      STACK OF EM VOLUMES
C
C        MZ   ==  1       IMAGE 
C        MZ   >=  1       IMAGE STACK
C        MZ   == NZ       VOLUME
C        MZ   NZ/NUMVOLS  VOLUME STACK
C
C        DMAX  < DMIN                        MAX & MIN UNDETERMINED
C        DMEAN < (SMALLER OF DMIN and DMAX)  DMEAN     UNDETERMINED
C        RMS   < 0.0                         RMS       UNDETERMINED
C
C  CALL TREE:   
C     OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC -->
C                             --> OPENFIL_N_MRC -->
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE OPENFIL_N_MRC(LUN,FILNAM,
     &                         MRCMODE,
     &                         NX,NY,NZ,NSTACK,
     &                         SCALEX,SCALEY,SCALEZ,IRTFLG)
 
C     PURPOSE: OPEN NEW MRC FILE FOR STREAM ACCESS, LITTLE ENDED

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'
 
      INTEGER                 :: LUN,MRCMODE
      CHARACTER(LEN=*)        :: FILNAM
      INTEGER                 :: NX,NY,NZ,NSTACK
      REAL                    :: SCALEX,SCALEY,SCALEZ
      INTEGER                 :: IRTFLG

      LOGICAL                 :: WANTUL,IS_MRCS 
      INTEGER                 :: NLABL,NSYMBT,MACHST,NE,IOFFSET
      INTEGER                 :: MX,MY,MZ,IVERSION,NZMRC
      INTEGER                 :: NSYMBYT,ISPG
      CHARACTER(LEN=1)        :: NULL = CHAR(0)


      !write(3,*)' In openfil_n_mrc, filnam,nstack: ',filnam,nz,nstack

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

      !write(3,*)' In openfil_n_mrc, filnam: ',filnam,nx,ny,nz,nstack

C     SEE IF USER WANTS RELION COMPATIBLE IMAGE/STACK
      IS_MRCS = (INDEX(FILNAM,'MRCS') > 1  .OR.
     &           INDEX(FILNAM,'mrcs') > 1)

      !write(3,*)' In openfil_n_mrc, filnam: ',is_mrcs, nstack
      IF ( IS_MRCS .OR. (NZ > 1 .AND. NSTACK > 1) )  THEN
C        HACK TO CREATE RELION COMPATIBLE STACKS
         IVERSION = 1
         ISPG     = 0    
         NZMRC    = MAX(NSTACK,1)        ! NUMBER OF IMAGES 
         MZ       = NSTACK
         !write(3,*)' In openfil_n_mrc, ismrcs,nstack,nzmrc: ', is_mrcs, nstack,nzmrc
      ELSEIF (NZ == 1 .AND. NSTACK < 0) THEN
C        JUST A SINGLE IMAGE
         IVERSION = 20140        
         ISPG     = 0    
         NZMRC    = 1           ! IMAGE           
         MZ       = 1
         NSTACK   = -2
      ELSEIF (NZ == 1 .AND. NSTACK >= 1) THEN
C        A 2014 STANDARD IMAGE STACK
         IVERSION = 20140        
         ISPG     = 0    
         NZMRC    = NZ          ! IMAGE           
         MZ       = NSTACK
      ELSEIF (NZ > 1 .AND. NSTACK < 0 ) THEN
C        SINGLE VOLUME
         IVERSION = 20140
         ISPG     = 1    
         NZMRC    = NZ          ! NUMBER OF IMAGES??
         MZ       = NZ
         NSTACK   = -2
         !write(3,*)' In openfil_n_mrc, nstack,nzmrc: ', nstack,nzmrc,nz
      ELSEIF (NZ > 1 .AND. NSTACK > 0 ) THEN
C        STACK OF EM VOLUMES
         IVERSION = 20140       ! SO NZ IS NOT FOR STACK
         ISPG     = 401    
         NZMRC    = NZ          ! NUMBER OF VOLUMES
         MZ       = NZ * NSTACK
         NSTACK   = MZ          ! STACK SIZE
      ELSE 
         CALL ERRT(102,'BAD STACK OR VOLUME PARAMETERS',NSTACK)
         IRTFLG = 1
         RETURN
      ENDIF

      !write(3,*)' In openfil_n_mrc, set:', nx,ny,nz,MX,MY,MZ

C     MODE IS ALWAYS 32 BIT REAL FOR NEW FILES (FOR NOW)
      CALL LUNSETMODSIZ_MRC(LUN,MRCMODE,NX,NY,NZMRC,
     &                                  MX,MY,MZ,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET EXTRA HEADER BYTES, NO. OF LABELS USED, VERSION, ETC
      NSYMBT = 0                      ! RESET BY LUNSETXXX_MRC
      NLABL  = 1                      ! RESET BY LUNSETXXX_MRC
      !write(3,*)' In oprnfil_n, set:', iversion,ispg,nzmrc,mz,nstack
      CALL LUNSETVIN_MRC(LUN,IVERSION,ISPG,NSYMBT,NLABL,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

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







