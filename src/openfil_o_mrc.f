C++*********************************************************************
C
C OPENFIL_O_MRC.F   ADAPTED FROM OPENFIL         MAY 2019 ArDean Leith
C                   NLABEL BUG                   SEP 2025 ArDean Leith
C                   DEBUG OUTPUT ADDED           SEP 2025 ArDean Leith
C                   IRTFLG AND NLABL CHANGED     SEP 2025 ArDean Leith
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
C  OPENFIL_O_MRC(LUNT,FILNAM,NLET,DSP,IMGNUM,
C                NX,NY,NZ, NSTK_FLG, IRTFLG)
C
C  PURPOSE: OPEN EXISTING MRC DATA FILE FOR READ/WRITE  
C
C  PARAMETERS:
C       LUNT     LOGICAL UNIT NUMBER FOR FILNAM.                (SENT)
C       FILNAM   FILE NAME (NO @)                               (SENT)
C       NLET     LENGTH OF FILE NAME                            (SENT)
C       DSP      OLD OR NEW FILE                                (SENT)
C       IMGNUM   IMAGE NUMBER IF STACKED                        (SENT)
C
C       NX,NY    X & Y DIMENSIONS OF IMAGE/VOL                  (RET)
C       NZ       Z  DIMENSION OF VOL OR MRCS_NSTK               (RET)
C       NSTK_FLG STACK INDICATOR                                (RET)
C                   -2  :  NOT A STACK 
C                   -1  :  AST FOR NSTK  (@* or *@) 
C                    0  :  BARE STACK    (@)
C                   >=0 :  MAX IMAGE NUMBER NOW IN STACK
C
C       GET_IMG_STK  GETS NSTK_FLG: NSTK or NOT or BARE 
C                  -2 :  NOT A STACK   (NO @)
C                  -1 :  AST FOR NSTK  (@* or *@)   
C                   0 :  BARE STACK    (@)
C                  >1 :  STACKED IMAGE NUMBER
C
C       IRTFLG   ERROR RETURN FLAG.                             (RET)
C                   0 : NORMAL RETURN
C                   1 : ERROR RETURN
C
C  VARIABLES:  
C           (FOR MRC 2014)
C                 ISPG (23)   NZ (1)     MZ (10)
C       IMAGE           0      1          1
C       IMAGE STACK     0     NSTK       NSTK
C       VOLUME (Z)            NZ         NZ
C       VOLUME STACK  401     NZ*NSTK    NZ/NSTK
C
C           (FOR MRCS)
C                 ISPG (23)   NZ (1)     MZ (10)
C       IMAGE           0     1          1
C       IMAGE STACK     0     NSTK       1
C       VOLUME (Z)            NZ         NZ
C
C       (23) ISPG == 1     CONTAINS SINGLE VOLUME
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
     &                         NX,NY,NZ, NSTK_FLG, IRTFLG)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=*)        :: FILNAM
      CHARACTER(LEN=1)        :: DSP
      INTEGER                 :: LUNT,NLET,IMGNUM
      INTEGER                 :: NX,NY,NZ,NSTK_FLG
      INTEGER                 :: IRTFLG

      INTEGER                 :: LUN,NE
      INTEGER                 :: I,MX,MY,MZ
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
      INTEGER                 :: ISPG,INOW,IGO,IEND,ILAB,iflip
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

C     SET READ FOR NATIVE LITTLE ENDIANNESS  
      DISP(1:1) = 'O'
      ENDED     = 'L'
      DISP(2:2) = ENDED

#if defined (SP_DBUGIO)
      write(3,*)' In openfil_o_mrc; imgnum:   ', imgnum
      write(3,*)' In openfil_o_mrc; filnam:   ', filnam(:nlet)
      write(3,*)' In openfil_o_mrc; disp:     ', disp
      write(3,*)'   '
      write(3,*)' In openfil_o_mrc; calling loadhed_mrc    '
      write(3,*)'  '
#endif

C     OPEN EXISTING FILE, READ HEADER & LOAD IN HEADER OBJECT  
      CALL LOADHED_MRC(LUN,FILNAM,NLET,DISP,IRTFLG)

#if defined (SP_DBUGIO)
      write(3,*)' In openfil_o_mrc; After calling loadhed_mrc 111  '
#endif

C     DETERMINE CURRENT FILE ENDED-NESS AS READ IN
      CALL LUNGETMAP_MRC(LUN,MAPC,MAPR,MAPS,IRTFLG)

#if defined (SP_DBUGIO)
      write(3,*)' In openfil_o_mrc; mapc,mapr,maps: ',mapc,mapr,maps
#endif

C     CHECK ENDIANNESS 
      ENDOK = ((MAPC == 1) .OR. (MAPR == 1) .OR. (MAPS == 1))
      IF (ENDOK)  THEN

         CALL LUNSETFLIP(LUN,0,IRTFLG)   ! NO FLIP NEEDED

      ELSE
C        NOT RIGHT ENDEDNESS, TRY CONVERT FROM BIG ENDED 
        
         ENDED     = 'B'
         DISP(2:2) = ENDED

 
C        OPEN FILE, CREATE HEADER OBJECT 

#if defined(SP_DBUGIO)
         write(3,*)' In openfil_o_mrc, Call loadhed_mrc, disp:',disp
#endif

C        OPEN EXISTING FILE, READ HEADER & LOAD IN HEADER OBJECT  
         CLOSE (LUN)
         CALL LOADHED_MRC(LUN,FILNAM,NLET,DISP, IRTFLG)

C        DETERMINE CURRENT FILE ENDED-NESS AS READ IN
         CALL LUNGETMAP_MRC(LUN,MAPC,MAPR,MAPS,IRTFLG)

#if defined(SP_DBUGIO)
         write(3,*)' In openfil_o_mrc; mapc,mapr,maps: ',mapc,mapr,maps
#endif

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

                           WRITE(NDAT,904) MAPC,MAPR,MAPS
          IF (NDAT .NE. 6) WRITE(6,904)    MAPC,MAPR,MAPS
904       FORMAT(' *** HEADER MAPC,MAPR,MAPS: ',I8,' ',I8,' ',I8,
     &           ' ARE NONSENSE ')

          MSG = 'OLD FILE HEADER NOT MRC/MAP:  ' // FILNAM(1:NLET)
          CALL ERRT(101,MSG,NE)
          IRTFLG = 1
          RETURN
      ENDIF

C     PUT FILNAM IN COMMON LUN HEADER
      CALL LUNSETFILE_MRC(LUN,FILNAM(1:NLET),NULL,IRTFLG)

C     GET FILE  MODE
      CALL LUNGETMODE_MRC(LUN,MRCMODE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN        

#if defined(SP_DBUGIO)
      write(3,*)' In openfil_o_mrc; mrcmode: ',mrcmode
#endif

      IF (MYPID <= 0 .AND. (MRCMODE == 3 .OR. MRCMODE == 4)) THEN
         CALL ERRT(102,'UNSUPPORTED MRC TRANSFORM FILE, MODE',MRCMODE)
         IRTFLG = 1
         RETURN
      ENDIF

C     GET MRC FILE PARAMETERS: IVERSION,ISPG,NSYMBT,NLABL
      CALL LUNGETVIN_MRC(LUN,IVERSION,ISPG,NSYMBT,NLABL,IRTFLG)

#if defined (SP_DBUGIO)
      write(3,*)' In openfil_o_mrc; iversion,ispg: ', iversion,ispg
      write(3,*)' In openfil_o_mrc; nsymbt, nlabl: ', nsymbt, nlabl
#endif
      IF (IRTFLG .NE. 0) RETURN 
       
C     READ NLABL LABELS FROM FILE (<= 800 CHAR)
      CALL LUNGETLABELS_MRC(LUN,NLABL,LABELS,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

#if defined (SP_DBUGIO)
      write(3,*)' In openfil_o_mrc; nlabl: ',nlabl
      !write(3,*)' In openfil_o_mrc; First label: ' ,labels(1:80)
      call lungetflip(lun,iflip,irtflg)
      write(3,*)' In openfil_o_mrc, iflip: ',iflip
#endif

C     GET MRC FILE PARAMETERS: NZ,MZ, & NSTK_FLG
C               (SHOULD HANDLE .mrcs FILES OK)

C     GET IMAGE/VOLUME SIZE PARAMETERS AND MODE
      CALL LUNGETMODSIZES_MRC(LUN,NX,NY,NZ, 
     &                        MX,MY,MZ,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN        


      CALL LUNGETNSTACK_MRC(LUN, NX,NY,NZ, MZ, NSTK_FLG,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN        

C     NZ FROM FILE IS NOT VOLUME DIMENSION IF MRC VERSION: 2014*
C                 ISPG (23)   NZ (1)     MZ (10)
C       IMAGE           0      1          1
C       IMAGE STACK     0     NSTK       NSTK
C       VOLUME (Z)            NZ         NZ


#if defined (SP_DBUGIO)
      write(3,*)' In openfil_o_mrc  88 ; nx,ny:    ', nx,ny
      write(3,*)' In openfil_o_mrc; 88 ; ispg:     ', ispg 
      write(3,*)' In openfil_o_mrc  88 ; nz,mz:    ', nz,mz
      write(3,*)' In openfil_o_mrc; 88 ; nstk_flg: ', nstk_flg
#endif

      IF (NSTK_FLG > = 0) THEN
         IMGNUM = NSTK_FLG
      ENDIF

#if defined (SP_DBUGIO)
      write(3,*)' In openfil_o_mrc; imgnum,nstk_flg: ',imgnum,nstk_flg
#endif

C      BARE STACK REFERENCE WITHOUT '@' ALLOWED FROM 'ST' & 'LI'  ?????
C      IF (MZ > 1 .AND. IMGNUM == 0 ) THEN 
C        BARE STACK REFERENCE NOT ALLOWED WITHOUT '@' NORMALLY
C         CALL ERRT(101,'STACK INDICATOR (@) MISSING',NE)
C         IRTFLG = 3
C         RETURN
C      ENDIF

C     ENSURE STACKED IMAGE EXISTS
      IF (DSP == 'O' .AND. (IMGNUM > 1 .AND. IMGNUM > NSTK_FLG)) THEN
         CALL ERRT(102,'IMAGE EXCEEDS HIGHEST IMAGE IN STACK',MZ)
         IRTFLG = 1
         RETURN        
      ENDIF

C     UPDATE MZ AND NSTK_FLG IN MRC HEADER FOR NEW STACKED IMAGE

      IF (DSP == 'N' .AND.
     &    IMGNUM  > 1 .AND. 
     &    (NSTK_FLG > 0 .or. MZ == 1) .AND. 
     &    IMGNUM > NSTK_FLG) THEN

         NSTK_FLG = IMGNUM
      
         ! EXTENDING STACK WITH A NEW IMAGE/VOLUME

#if defined (SP_DBUGIO)
         write(3,*)' In openfil_o_mrc 2 Set|   nz,nstk_flg: ',  
     &                                         nz,nstk_flg
         write(3,*)' In openfil_o_mrc 2; Calling lunsetnstack ----'
#endif

         CALL LUNSETNSTACK_MRC(LUN, NZ, NSTK_FLG,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN        
      ENDIF


C     SET SPIDER IFORM IN COMMON CMBLOCK
      IFORM = 2
      IF (NZ > 1) IFORM = 3 

C     SET IMGNUM AND NSTK_FLG IN STATIC HEADER

#if defined (SP_DBUGIO)
      write(3,*)'  '
      write(3,*)' In openfil_o_mrc 3  Set| imgnum,nstk_flg: ',
     &                                     imgnum,nstk_flg
      write(3,*)' In openfil_o_mrc 3  Calling lunsetstk ----'
#endif
 
      CALL LUNSETSTK_MRC(LUN,IMGNUM,NSTK_FLG,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN        

#if defined (SP_DBUGIO)
      write(3,*)' In openfil_o_mrc| fchar(1:2),fchar(4:5): ',
     &                              fchar(1:2),'  ',fchar(4:5)
      write(3,*) '    '
#endif

      IF (FCHAR(1:2) == '31' .AND. FCHAR(4:5) == 'HE') THEN

C     IF (FCHAR(4:5) ==  'HE' .AND.  FCHAR(1:2) ==  '31') THEN
C        SPECIAL KLUDGE FOR 'MRC HED' SO IT DOESN'T ALTER STATS
C        SET FLAG FOR NORMAL RETURN

#if defined (SP_DBUGIO)
         write(3,*)' In openfil_o_mrc; special return for MRC HED '
         write(3,*)'    '
#endif

         IRTFLG = 0
         RETURN
      ENDIF

C     ADJUST FMIN...  FOR DIFFERENT MRC HEADER POSSIBILITIES
C     DOES NOT PUT FILE VALUES FOR FMIN... INTO COMMON: CMBLOCK   
C     FILE VALUES FOR FMIN... IN COMMON: CMBLOCK SET LATER WITH  
C     LUNSETCOMMON

      CALL LUNGETSTATS_MRC(LUN,IMAMIT,FMINT,FMAXT,AVT,SIGT,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN        
 
#if defined (SP_DBUGIO)
      write(3,*)' In openfil_o_mrc 4  Got imamit,fmaxt: ',imamit,fmaxt 
#endif

C     ENSURE THAT STAT'S IMAGE NUMBER IS SET CORRECTLY 
    
      IF (IMAMIT == 1 .AND. NSTK_FLG <= 0) THEN
C        NOT A STACK, SO FMIN/FMAX SHOULD BE OK FOR THIS MRC FILE
         CALL LUNSET_STATSIMG_MRC(LUN,1,IRTFLG)

      ELSEIF (IMAMIT == 0) THEN
C        NO IDENTIFIABLE STATS IN MRC IMAGE
 
#if defined (SP_DBUGIO)
         write(3,*)' In openfil_o_mrc 5; imamit: ',imamit
         write(3,*)' In openfil_o_mrc 5; Calling lunset_statsimg '
         write(3,*)'  '
#endif

         CALL LUNSET_STATSIMG_MRC(LUN,0,IRTFLG)

      ENDIF

#if defined (SP_DBUGIO)
      write(3,*)' From openfil_o_mrc; iform:    ', iform 
      write(3,*)' from openfil_o_mrc; nx,ny,nz: ', nx,ny,nz 
#endif

C     SET FLAG FOR NORMAL RETURN	
      IRTFLG = 0
      
      END


C     ------------------- LOADHED_MRC ---------------------------------

      SUBROUTINE LOADHED_MRC(LUN,FILNAM,NLET,DISP, IRTFLG)

C     COPY MRC FILE HEADER INTO: MRC HEADER

C     DISP  IS INPUT FOR OPSTREAMFILE USE

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
        
      IRTFLG  = 0
      !IRTFLG = 999   ! DO NOT ECHO OPENING INFO  (sept 2025 al)
      CALL OPSTREAMFILE(.FALSE.,FILNAM(1:NLET),NULL,LUN,
     &                  'UNFORMATTED',DISP,' ',.TRUE.,IRTFLG)

#if defined(NEVER_SP_DBUGIO)
C23456  
      write(3,*)' Opened old file; Disp,irtflg: ',disp,irtflg  
      write(3,*)' Opened old file:              ',filnam(1:nlet)
#endif
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





