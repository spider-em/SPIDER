C++*********************************************************************
C
C OPENFIL_O_MRC.F   ADAPTED FROM OPENFIL         MAY 2019 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2019  Health Research Inc.,                         *
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
C  OPENFIL_O_MRC(LUNT,FILNAM,NLET,IMGNUM, NX,NY,NZ,NSTACK, IRTFLG)
C
C  PURPOSE: OPEN EXISTING MRC DATA FILE FOR READ/WRITE  
C
C  PARAMETERS:
C       LUNT     LOGICAL UNIT NUMBER FOR FILNAM.                (SENT)
C       FILNAM   FILE NAME (NO @)                               (SENT)
C       NLET     LENGTH OF FILE NAME                            (RET.)
C       DSP      OLD OR NEW FILE                                (SENT)
C       IMGNUM   IMAGE NUMBER IF STACKED                        (SENT)
C       NX,NY,NZ DIMENSIONS OF IMAGE                            (RET.)
C       NSTACK   STACK INDICATOR                                (RET.)
C                   -2  :  NOT A STACK 
C                   >=0 :  MAX IMAGE NUMBER NOW IN STACK.
C       IRTFLG   ERROR RETURN FLAG.                             (RET.)
C                   0 : NORMAL RETURN
C                   1 : ERROR RETURN
C
C  VARIABLES:
C       (23) ISPG == 0     CONTAINS IMAGE OR IMAGE STACK
C       (23) ISPG == 1     CONTAINS VOLUMES
C       (23) ISPG == 401   CONTAINS STACK OF EM VOLUMES
C       (24) NSYMBT        NO. OF EXTRA SYMMETRY, ETC. BYTES
C       (56) NLABL         NO. OF 80 CHAR TEXT LABELS IN #57..256
C
C  CALL TREE:   
C     OPFILEC --> OPENFIL_MRC --> OPENFIL_O_MRC -->
C                             --> OPENFIL_N_MRC -->
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE OPENFIL_O_MRC(LUNT,FILNAM,NLET,DSP,IMGNUM,
     &                         NX,NY,NZ,NSTACK, IRTFLG)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=*)        :: FILNAM
      CHARACTER(LEN=1)        :: DSP
      INTEGER                 :: LUNT,NLET,IMGNUM, NX,NY,NZ,NSTACK
      INTEGER                 :: IRTFLG

      INTEGER                 :: LUN,NE
      INTEGER                 :: I,NZMRC,MX,MY,MZ
      INTEGER                 :: IVERSION,NSYMBT,IHEDLEN,NLABL,NBYT
      CHARACTER(LEN=1)        :: ENDED
      CHARACTER(LEN=2)        :: DISP
      CHARACTER(LEN=2*MAXNAM) :: MSG
      CHARACTER(LEN=1)        :: NULL = CHAR(0)

      CHARACTER(LEN=800)      :: LABELS
      INTEGER                 :: MRCMODE,IVAL
      REAL                    :: FVAL
      REAL                    :: SCALEX,SCALEY,SCALEZ
      INTEGER                 :: IMAMIT
      REAL                    :: FMINT,FMAXT,AVT,SIGT
      INTEGER                 :: ISPG,INOW,IGO,IEND,ILAB
      LOGICAL                 :: ENDOK
      INTEGER                 :: LOCAT   
      INTEGER                 :: MAPC,MAPR,MAPS
      LOGICAL                 :: WANTUL,ISMRCSFILE 

      LOGICAL                 :: ISDIGI,lnblnkn  ! FUNCTIONS

C     INCLUSION FOR OPTIONAL MPI INITIALIZATION.  
      INTEGER                 :: MYPID = -1
#include "MPI_INIT.INC"

C     PURPOSE:  OPEN EXISTING MRC FILE FOR STREAM ACCESS
           
      LUN = LUNT

C     SET READ TO CONVERT TO LITTLE ENDEDNESS
      DISP(1:1) = 'O'
      ENDED     = 'L'
      DISP(2:2) = ENDED

      !write(3,*)' In openfil_o_mrc, nlet,filnam:',nlet,filnam(:nlet)
      !write(3,*)' In openfil_o_mrc, imgnum:',imgnum

C     OPEN EXISTING FILE, CREATE HEADER OBJECT, AND CHECK ENDEDNESS 
      !write(6,*)' In openfil_o_mrc, call loadhed_mrc, disp:',disp
      CALL LOADHED_MRC(LUN,FILNAM,NLET,DISP,IRTFLG)
      
C     DETERMINE CURRENT FILE ENDED-NESS AS READ IN
      CALL LUNGETMAP_MRC(LUN,MAPC,MAPR,MAPS,IRTFLG)
      !write(3,*)' In openfil_o_mrc, mapc,r,s: ',mapc,mapr,maps

      ENDOK = ((MAPC == 1) .OR. (MAPR == 1) .OR. (MAPS == 1))
      IF (ENDOK)  CALL LUNSETFLIP(LUN,0,IRTFLG)

      IF (.NOT. ENDOK) THEN
C        NOT RIGHT ENDEDNESS, TRY CONVERT FROM BIG ENDED 
        
         ENDED     = 'B'
         DISP(2:2) = ENDED

C        OPEN FILE, CREATE HEADER OBJECT AND CHECK ENDEDNESS 

         !write(6,*)' In openfil_o_mrc, call loadhed_mrc, disp:',disp
         CLOSE (LUN)
         CALL LOADHED_MRC(LUN,FILNAM,NLET,DISP, IRTFLG)

C        DETERMINE CURRENT FILE ENDED-NESS AS READ IN
         CALL LUNGETMAP_MRC(LUN,MAPC,MAPR,MAPS,IRTFLG)
         !write(3,*)' In openfil_o_mrc, mapc,r,s: ',mapc,mapr,maps

         ENDOK = ((MAPC == 1) .OR. (MAPR == 1) .OR. (MAPS == 1))

         IF (ENDOK) THEN
C           MUST FLIP FILE CONTENTS WHEN READ IN
            CALL LUNSETFLIP(LUN,1,IRTFLG)
         ELSE
            IRTFLG = 1
         ENDIF

      ENDIF
      
      IF (MYPID <= 0 .AND. IRTFLG .NE. 0) THEN        
C         NEITHER ENDEDNESS SEEMS OK --> OPENING ERROR 
          MSG = 'OPENFIL_O_MRC; OPENING FILE: ' // FILNAM(1:NLET)
          CALL ERRT(101,MSG,NE)
          IRTFLG =1
          RETURN
      ENDIF

C     GET IMAGE/VOLUME SIZE PARAMETERS AND MODE
      CALL LUNGETMODSIZ_MRC(LUN,MRCMODE,NX,NY,NZMRC, MX,MY,MZ,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN        

      !write(6,*)' In openfil_o_mrc, mrcmode: ',mrcmode
      !write(3,*)' In openfil_o_mrc, size (x,y,z),mz: ',nx,ny,nzmrc,mz

      IF (MYPID <= 0 .AND. (MRCMODE == 3 .OR. MRCMODE == 4)) THEN
         CALL ERRT(102,'UNSUPPORTED MRC TRANSFORM FILE, MODE',MRCMODE)
         IRTFLG = 1
         RETURN
      ENDIF

      !write(3,*)' In openfil_o_mrc, imgnum,mz: ',imgnum,mz

C     READ NLABL LABELS FROM FILE (<= 800 CHAR)
      CALL LUNGETLABELS_MRC(LUN,NLABL,LABELS,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN        

C     GET MRC FILE PARAMETERS: IVERSION,ISPG,NSYMBT,NLABL
      CALL LUNGETVIN_MRC(LUN,IVERSION,ISPG,NSYMBT,NLABL,IRTFLG)
      !write(3,*)' In openfil_o_mrc, ispg,iversion: ', ispg,iversion
      !write(6,*)' In openfil_o_mrc, nsymbt, nlabl: ', nsymbt, nlabl
      IF (IRTFLG .NE. 0) RETURN        
  
C     GET MRC FILE PARAMETERS:  NZ,MZ, & NSTACK.(HANDLES .mrcs FILES OK)
      PRINT *, "opfil_o_mrc.f : 182: Calling LUNGETNSTACK_MRC"
      CALL LUNGETNSTACK_MRC(LUN, NX,NY,NZ,MZ,NSTACK,IRTFLG)
      PRINT *, "opfil_o_mrc.f : 184: OPENFIL_O_MRC: NSTACK",NSTACK
      IF (IRTFLG .NE. 0) RETURN        

      !write(3,*)' In openfil_o_mrc, filnam: ',filnam(1:nlet)
      !write(3,*)' In openfil_o_mrc, nz,mz,nstack: ',nz,mz,nstack !!!!!
      !write(3,*)' In openfil_o_mrc,imgnum,mz,nstack: ',imgnum,mz,nstack

C     BARE STACK REFERENCE WITHOUT '@' ALLOWED FROM 'ST' & 'LI'
      IF (MZ > 1 .AND. IMGNUM < -1 ) THEN 
C        BARE STACK REFERENCE NOT ALLOWED WITHOUT '@' NORMALLY
         CALL ERRT(101,'STACK INDICATOR (@) MISSING',NE)
         IRTFLG = 3
         RETURN
      ENDIF

C     ENSURE STACKED IMAGE EXISTS
      IF (DSP == 'O' .AND. (IMGNUM > 1 .AND. IMGNUM > NSTACK)) THEN
         CALL ERRT(102,'IMAGE EXCEEDS HIGHEST IMAGE IN STACK',MZ)
         IRTFLG = 1
         RETURN        
      ENDIF

      !write(3,*)' In openfil_o_mrc, dsp :'  , dsp 
      !write(3,*)' In openfil_o_mrc, imgnum:', imgnum 
      !write(3,*)' In openfil_o_mrc, nstack:', nstack

C     UPDATE MZ AND NSTACK IN MRC HEADER FOR NEW STACKED IMAGE
      IF (DSP == 'N' .AND.
     &    IMGNUM > 1 .AND. 
     &    (NSTACK > 0 .or. MZ == 1) .AND. 
     &    IMGNUM > NSTACK) THEN
         ! EXTENDING STACK WITH A NEW IMAGE/VOLUME
         !write(3,*)' In openfil_o_mrc2, dsp,imgnum:', dsp,imgnum,nstack
         NSTACK = IMGNUM      
         CALL LUNSETNSTACK_MRC(LUN, NZ,NSTACK,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN        
      ENDIF

C     SET SPIDER IFORM IN COMMON CMBLOCK
      IFORM = 2
      IF (NZ > 1) IFORM = 3 

C     SET IMGNUM AND NSTACK IN STATIC HEADER
      !write(6,*)' In openfil_o_mrc - imgnum,nstack:',imgnum,nstack 
      CALL LUNSETSTK_MRC(LUN,IMGNUM,NSTACK,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN        

      IF (FCHAR(4:5) ==  'HE' .AND.  FCHAR(1:2) ==  '31') THEN
C        SPECIAL KLUDGE FOR 'MRC HED' SO IT DOESN'T ALTER STATS
C        SET FLAG FOR NORMAL RETURN	
         IRTFLG = 0
         RETURN
      ENDIF

C     ADJUST FMIN...  FOR DIFFERENT MRC HEADER POSSIBILITIES
C     DOES NOT PUT FILE VALUES FOR FMIN... INTO COMMON: CMBLOCK   
C     FILE VALUES FOR FMIN... INTO COMMON: CMBLOCK SET LATER IN LUNSETCOMMON   
      CALL LUNGETSTATS_MRC(LUN,IMAMIT,FMINT,FMAXT,AVT,SIGT,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN        
 
      !write(3,*)' In openfil_o_mrc imami,fmin:',imami,fmin 

C     ENSURE THAT STAT'S IMAGE NUMBER IS SET CORRECTLY     
      IF (IMAMI == 1 .AND. NSTACK <= 0) THEN
C        NOT A STACK, SO FMIN/FMAX SHOULD BE OK FOR THIS MRC FILE
         CALL LUNSET_STATSIMG_MRC(LUN,1,IRTFLG)

      ELSEIF (IMAMI == 0) THEN
C        NO IDENTIFIABLE STATS IN MRC IMAGE 
         !write(3,*)'In openfil_o_mrc, fmin,fmax,imami:',fmin,fmax,imami

         CALL LUNSET_STATSIMG_MRC(LUN,0,IRTFLG)
      ENDIF

      !write(6,*)' In openfil_o_mrc - iform,nz: ',iform,nz,nstack 
      !write(6,*)' In openfil_o_mrc - nx,ny,nz: ',nx,ny,nz,nstack 

C     SET FLAG FOR NORMAL RETURN	
      IRTFLG = 0
      
      END


C     ------------------- LOADHED_MRC ---------------------------------

      SUBROUTINE LOADHED_MRC(LUN,FILNAM,NLET,DISP, IRTFLG)

C     COPY MRC FILE HEADER INTO: MRC HEADER

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=*)      :: FILNAM,DISP
      INTEGER               :: LUN,NLET,IRTFLG

      INTEGER               :: NE,I
      CHARACTER(LEN=1)      :: NULL = CHAR(0)

C     INCLUSION FOR OPTIONAL MPI INITIALIZATION.  
      INTEGER               :: MYPID = -1
#include "MPI_INIT.INC"

C     WANT TO OPEN EXISTING MRC FILE FOR STREAM ACCESS
        
c     write(6,*)' Opening old MRC file: ',FILNAM(1:nlet),
c     &          '  Convert:',disp  

      IRTFLG = 999   ! DO NOT ECHO OPENING INFO
      CALL OPSTREAMFILE(.FALSE.,FILNAM(1:NLET),NULL,LUN,
     &                  'UNFORMATTED',DISP,' ',.TRUE.,IRTFLG)

      !write(6,*)' Opened  old MRC file: ',FILNAM(1:nlet),
c     &           '  Convert:',DISP,' Irtflg: ',irtflg  
      IF (IRTFLG .NE. 0) RETURN

#ifdef USE_MPI
      IF (ONLYONE_RED) THEN
         CALL BCAST_MPI('OPENFIL_MRC','IRTFLG',IRTFLG,1,'I',ICOMM)
      ENDIF
#endif
     
      IF (MYPID <= 0 .AND. IRTFLG .NE. 0) THEN
          WRITE(NOUT,*) '*** ERROR OPENING HEADER OF: ',FILNAM(:NLET)
      ENDIF
      IF (IRTFLG .NE. 0) RETURN

C     READ OVERALL HEADER FROM MRC FILE 
      CALL LUNREDHED_MRC(LUN,.TRUE.,IRTFLG)

      IF (IRTFLG .NE. 0 .AND. MYPID <= 0) THEN
         WRITE(NOUT,*) '*** ERROR READING HEADER OF: ',FILNAM(:NLET)
         IF (IRTFLG .NE. 0) RETURN
      ENDIF
     
      END





