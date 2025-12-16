

C++*********************************************************************
C
C OPENFIL_O_MRC.F   ADAPTED FROM OPENFIL         MAY 2019 ArDean Leith
C                   NLABEL BUG                   SEP 2025 ArDean Leith
C                   DEBUG OUTPUT ADDED           SEP 2025 ArDean Leith
C                   IRTFLG AND NLABL CHANGED     SEP 2025 ArDean Leith
C                   REWRITE                      NOV 2025 ArDean Leith
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
C  OPENFIL_O_MRC(LUNT,FILNAM,DSP, NX,NY,NZ, 
C.               NSTK, ISTK, IS_BARE,  ITYPE, IRTFLG)
C
C  PURPOSE: OPEN EXISTING MRC DATA FILE FOR READ/WRITE  
C
C  PARAMETERS:
C       LUNT     LOGICAL UNIT NUMBER FOR FILNAM.                (SENT)
C       FILNAM   FILE NAME (NO @)                               (SENT)
C       DSP      OLD OR NEW FILE  (O/R/N)                       (SENT)
C
C       NX,NY    X & Y DIMENSIONS OF IMAGE/VOL                   (RET)
C       NZ       Z  DIMENSION OF VOL OR MRCS_NSTK                (RET)
C
C       NSTK     STACK SIZE.  (IF STACKED)                       (RET)
C       ISTK     IMAGE NUMBER (IF STACKED)                  (SENT/RET)
C
C.      ITYPE    FOURIER, ETC
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
C          -> OPENFIL_O_MRC 
C              -> LOADHED_MRC        
C              -> LUNGET_MAP_MRC  
C              -> LUNSETFLIP_MRC 
C              -> LOADHED_MRC        
C              -> LUNGET_MAP_MRC  
C              -> LUNSETFLIP_MRC 
C              -> LUNSET_FILE_MRC    
C              -> LUNGET_MODE_MRC 
C              -> LUNGET_VIN_MRC 
C              -> LUNGET_LABELS_MRC  
C              -> LUNGET_MODSIZES_MRC 
C              -> LUNGET_2014_MRC
C              -> LUNSET_2014_MRC      (NSTK MAY INCREASE)
C              -> LUNGET_STATS_MRC
C          IF (DSP/=R) -> LUNSET_STATSIMG_MRC 
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE OPENFIL_O_MRC(LUNT,FILNAM,DSP,
     &                          NX,NY,NZ, 
     &                          NSTK,ISTK,IS_BARE,
     &                          ITYPE, IRTFLG)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=*)        :: FILNAM
      CHARACTER(LEN=1)        :: DSP
      INTEGER                 :: LUNT, NX,NY,NZ, NSTK, ISTK 
      LOGICAL                 :: IS_BARE
      INTEGER                 :: ITYPE,IRTFLG

      INTEGER                 :: LUN,NE
      INTEGER                 :: I,MX,MY,MZ,NZ3
      INTEGER                 :: IVERSION,NSYMBT,IHEDLEN,NLABL,NBYT
      CHARACTER(LEN=1)        :: ENDED
      CHARACTER(LEN=2)        :: DISP
      CHARACTER(LEN=2*MAXNAM) :: MSG

      CHARACTER(LEN=800)      :: LABELS
      INTEGER                 :: MRCMODE,IVAL
      REAL                    :: FVAL
      REAL                    :: SCALEX,SCALEY,SCALEZ
      INTEGER                 :: IMAMIT
      REAL                    :: FMINT,FMAXT,AVT,SIGT
      INTEGER                 :: ISPG,INOW,IGO,IEND,ILAB,iflip
      LOGICAL                 :: ENDOK
      INTEGER                 :: MAPC,MAPR,MAPS

      LOGICAL                 :: lnblnkn  ! FUNCTIONS

C     INCLUSION FOR OPTIONAL MPI INITIALIZATION.  
      INTEGER                 :: MYPID = -1
#include "MPI_INIT.INC"

           
      LUN = LUNT

C     SET READ FOR NATIVE LITTLE ENDIANNESS  
      DISP(1:1) = 'O'
      ENDED     = 'L'
      DISP(2:2) = ENDED

#if defined(SP_DBUGIO)
      write(3,*)' In openfil_o_mrc; istk:   ', istk
      write(3,*)' In openfil_o_mrc; filnam: ', trim(filnam)
      write(3,*)' In openfil_o_mrc; disp:   ', disp
      write(3,*)' In openfil_o_mrc; calling loadhed_mrc    '
      write(3,*)'  '
#endif

C     READ EXISTING HEADER FROM FILE & LOAD IN HEADER OBJECT 
       
      CALL LOADHED_MRC(LUN,FILNAM,DISP,IRTFLG)

#if defined(SP_DBUGIO)
      write(3,*)' In openfil_o_mrc; After calling loadhed_mrc 111  '
#endif

C     DETERMINE CURRENT FILE ENDED-NESS AS READ IN
      CALL LUNGET_MAP_MRC(LUN,MAPC,MAPR,MAPS,IRTFLG)

#if defined(SP_DBUGIO)
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

         CALL LOADHED_MRC(LUN,FILNAM,DISP, IRTFLG)

C        DETERMINE CURRENT FILE ENDED-NESS AS READ IN
         CALL LUNGET_MAP_MRC(LUN,MAPC,MAPR,MAPS,IRTFLG)

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

          MSG = 'OLD FILE HEADER NOT MRC/MAP:  ' // TRIM(FILNAM)
          CALL ERRT(101,MSG,NE)
          IRTFLG = 1
          RETURN
      ENDIF



C     PUT FILNAM IN COMMON LUN HEADER ----------------------------
      CALL LUNSET_FILE_MRC(LUN,TRIM(FILNAM),DSP,IRTFLG)

C     GET FILE  MODE
      CALL LUNGET_MODE_MRC(LUN,MRCMODE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN        

#if defined(SP_DBUGIO)
      write(3,*)' In openfil_o_mrc; mrcmode: ',mrcmode
#endif

      IF (MYPID <= 0 .AND. (MRCMODE == 3 .OR. MRCMODE == 4)) THEN
         CALL ERRT(102,'UNSUPPORTED MRC FOURIER FILE, MODE',MRCMODE)
         IRTFLG = 1
         RETURN
      ENDIF

C     GET FILE PARAMETERS:    IVERSION,ISPG,NSYMBT,NLABL
      CALL LUNGET_VIN_MRC(LUN,IVERSION,ISPG,NSYMBT,NLABL, IRTFLG)

#if defined(SP_DBUGIO)
      write(3,*)' In openfil_o_mrc; iversion,ispg: ', iversion,ispg
      write(3,*)' In openfil_o_mrc; nsymbt, nlabl: ', nsymbt,nlabl
#endif
      IF (IRTFLG .NE. 0) RETURN 
       
C     READ NLABL LABELS FROM FILE (<= 800 CHAR)
      CALL LUNGET_LABELS_MRC(LUN,NLABL,LABELS,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN


C     GET IMAGE/VOLUME SIZE PARAMETERS AND MODE
      CALL LUNGET_MODSIZES_MRC(LUN,NX,NY,NZ, 
     &                             MX,MY,MZ,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN        

C     GET MRC 2014 FILE PARAMETERS: NZ & NSTK
C     (SHOULD HANDLE .mrcs FILES OK)

      CALL LUNGET_2014_MRC(LUN, NX,NY,NZ3, MZ, NZ,NSTK ,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN        

#if defined(SP_DBUGIO)
      write(3,*)' In openfil_o_mrc  88 ; nx,ny:  ', nx,ny
      write(3,*)' In openfil_o_mrc; 88 ; ispg:   ', ispg 
      write(3,*)' In openfil_o_mrc  88 ; nz,mz:  ', nz,mz
      write(3,*)' In openfil_o_mrc; 88 ; nstk:   ', nstk
      write(3,*)' In openfil_o_mrc; 88 ; istk:   ', istk
#endif

C     BARE STACK REFERENCE WITHOUT '@' ALLOWED FROM 'ST' & 'LI'  ?????
C     IF (MZ > 1 .AND. ISTK == 0 ) THEN 
C        BARE STACK REFERENCE NOT ALLOWED WITHOUT '@' NORMALLY
C         CALL ERRT(101,'STACK INDICATOR (@) MISSING',NE)
C         IRTFLG = 3
C         RETURN
C     ENDIF

C     UPDATE MZ AND NSTK IN MRC HEADER OBJECT FOR NEW STACKED IMAGE

      IF (IVERSION >= 2014 .AND.
     &    DSP == 'N' .AND.
     &    ISTK  > 1  .AND. 
     &    (NSTK > 0 .OR. MZ == 1) .AND. 
     &    ISTK > NSTK) THEN

         NSTK = ISTK    ! EXTENDING STACK WITH A NEW IMAGE/VOLUME

#if defined(SP_DBUGIO)
         write(3,*)' In openfil_o_mrc 2 Set|   nz,nstk: ', nz,nstk
         write(3,*)' In openfil_o_mrc 2; Calling lunset_2014 ----'
#endif

         CALL LUNSET_2014_MRC(LUN, NZ,NSTK, IRTFLG)
         IF (IRTFLG .NE. 0) RETURN        
      ENDIF


C     ENSURE STACKED IMAGE EXISTS
      IF (DSP == 'O' .AND. (ISTK > 1 .AND. ISTK > NSTK)) THEN
         CALL ERRT(102,'IMAGE EXCEEDS HIGHEST IMAGE IN STACK',NSTK)
         IRTFLG = 1
         RETURN        
      ENDIF

C     SET SPIDER IFORM IN COMMON CMBLOCK
      IFORM = 2
      IF (NZ > 1) IFORM = 3 

#if defined(SP_DBUGIO)
      write(3,*)'  '
      write(3,*)' In openfil_o_mrc set| lun,istk,nstk: ',lun,istk,nstk
      write(3,*)' In openfil_o_mrc Calling lunset_stk_260 ----'
#endif

C     SET ISTK AND NSTK IN STATIC HEADER OBJECT 
      CALL LUNSET_STK_260_MRC(LUN,ISTK,NSTK,IRTFLG)


C     ADJUST FMIN...  FOR DIFFERENT MRC HEADER POSSIBILITIES
C     FILE VALUES FOR FMIN... IN COMMON: CMBLOCK WILL BE SET IN CALLER 
C.    WITH LUNSETCOMMON

      CALL LUNGET_STATS_MRC(LUN,IMAMIT,FMINT,FMAXT,AVT,SIGT,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN        
 
#if defined(SP_DBUGIO)
      write(3,*)' In openfil_o_mrc 4  Got imamit,fmaxt: ',imamit,fmaxt 
#endif


      IF (DSP .NE. 'R') THEN

C        SPECIAL RETURN FOR 'MRC HED' AND UNWRITABLE MRC FILES 

#if defined(SP_DBUGIO)
         write(3,*)' In openfil_o_mrc 5; istk: ',istk
         write(3,*)' In openfil_o_mrc 5; Calling lunset_statsimg '
         write(3,*)'  '
#endif

C        SET STAT'S IMAGE NUMBER IN IMAGE HEADER 
         CALL LUNSET_STATSIMG_MRC(LUN,ISTK,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN        
      ENDIF     

#if defined(SP_DBUGIO)
      write(3,*)' From openfil_o_mrc; nx,ny,nz: ', nx,ny,nz 
      write(3,*)' From openfil_o_mrc; mx,my,mz: ', nx,ny,nz 
      write(3,*)' From openfil_o_mrc; iform:    ', iform 
      write(3,*)' From openfil_o_mrc; istk:     ', istk 
      write(3,*)' From openfil_o_mrc; nstk:     ', nstk 
      write(3,*)' From openfil_o_mrc; dsp:      ', dsp
#endif

C     SET FLAG FOR NORMAL RETURN	
      IRTFLG = 0
      
      END




C     ------------------- LOADHED_MRC ---------------------------------

      SUBROUTINE LOADHED_MRC(LUN,FILNAM,DISP, IRTFLG)

C     COPY MRC FILE HEADER INTO: MRC HEADER

C     DISP  IS INPUT FOR OPSTREAMFILE USE

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=*)      :: FILNAM,DISP
      INTEGER               :: LUN,IRTFLG

      INTEGER               :: NE,I
      CHARACTER(LEN=1)      :: NULL = CHAR(0)

      INTEGER               :: MYPID = -1
C     INCLUSION FOR OPTIONAL MPI INITIALIZATION.  
#include "MPI_INIT.INC"

C     WANT TO OPEN EXISTING MRC FILE FOR STREAM ACCESS
        
#if defined(SP_DBUGIO)
      write(3,*) ' In loadhed_mrc; disp:    ',disp  
      write(3,*) ' In loadhed_mrc; lun:     ',lun  
      write(3,*) ' In loadhed_mrc; opening: ',trim(filnam)
#endif

      IRTFLG  = 0
      !IRTFLG = 999   ! DO NOT ECHO OPENING INFO  (sept 2025 al)
      CALL OPSTREAMFILE(.FALSE.,FILNAM,NULL,LUN,
     &                  'UNFORMATTED',DISP,' ',.TRUE.,IRTFLG)

#if defined(SP_DBUGIO)  
      write(3,*)' In loadhed_mrc; irtflg:  ',irtflg  
      write(3,*)' In loadhed_mrc; opened: ',trim(filnam)
#endif
      IF (IRTFLG .NE. 0) RETURN


#ifdef USE_MPI
      IF (ONLYONE_RED) THEN
         CALL BCAST_MPI('OPENFIL_MRC','IRTFLG',IRTFLG,1,'I',ICOMM)
      ENDIF
#endif
     
      IF (MYPID <= 0 .AND. IRTFLG .NE. 0) THEN
          WRITE(NOUT,*) '*** ERROR OPENING HEADER OF: ',TRIM(FILNAM)
      ENDIF
      IF (IRTFLG .NE. 0) RETURN

C     READ OVERALL HEADER FROM MRC FILE 

      CALL LUNREDHED_MRC(LUN,.TRUE.,IRTFLG)

      IF (IRTFLG .NE. 0 .AND. MYPID <= 0) THEN
         WRITE(NOUT,*) '*** ERROR READING HEADER OF: ',TRIM(FILNAM)
         IF (IRTFLG .NE. 0) RETURN
      ENDIF

      END





