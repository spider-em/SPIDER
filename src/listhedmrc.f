
C **********************************************************************
C 
C  LISTHEDMRC   NEW FOR MRC SUPPORT             OCT 2019  ArDean Leith
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
C  LISTHEDMRC(LUN,IRTFLG)
C                                                                      
C  PURPOSE:    LIST MRC IMAGE HEADER AS USED BY SPIDER     
C                                            
C  PARAMETERS: LUN           INPUT UNIT
C              IRTFLG        ERROR FLAG
C  
C  MAP/IMAGE HEADER FORMAT:
C 
C  LENGTH = 1024 BYTES, ORGANIZED AS 56 LONG WORDS FOLLOWED
C           BY SPACE FOR 10 80 BYTE TEXT LABELS.
C 
C  1    NX          # OF COLUMNS    (FASTEST CHANGING IN MAP)
C  2    NY          # OF ROWS
C  3    NZ          # OF SECTIONS   (SLOWEST CHANGING IN MAP)
C  4    MODE        DATA TYPE
C                   0   IMAGE:  8-BIT BYTES RANGE -128 -->127 (SIGNED)                  
C                   1   IMAGE: 16-BIT INTEGERS HALFWORDS                
C                   2   IMAGE: 32-BIT REALS                    
C                   6   IMAGE: UNSIGNED 8-BIT BYTES RANGE 0 -->255        
C                   3   TRANSFORM: COMPLEX 16-BIT INTEGERS (UNSUPPORTED)   
C                   4   TRANSFORM: COMPLEX 32-BIT REALS    (UNSUPPORTED)     
C  5    NXSTART     # OF FIRST COLUMN  IN MAP 
C  6    NYSTART     # OF FIRST ROW     IN MAP       
C  7    NZSTART     # OF FIRST SECTION IN MAP       
C  8    MX          # OF INTERVALS ALONG X
C  9    MY          # OF INTERVALS ALONG Y
C 10    MZ          # OF INTERVALS ALONG Z    (>1 FOR STACK)
C 11-13 CELLA       CELL DIMENSIONS IN ANGSTROMS
C 14-16 CELLB       CELL ANGLES IN DEGREES  (PHI,THETA,PSI)                        
C 17    MAPC        AXIS CORRESPONDING TO COLUMNS  (1,2,3 FOR X,Y,Z)
C 18    MAPR        AXIS CORRESPONDING TO ROWS     (1,2,3 FOR X,Y,Z)
C 19    MAPS        AXIS CORRESPONDING TO SECTIONS (1,2,3 FOR X,Y,Z)
C 20    DMIN        MINIMUM DENSITY VALUE
C 21    DMAX        MAXIMUM DENSITY VALUE  
C 22    DMEAN       MEAN    DENSITY VALUE  
C 23    ISPG        SPACE GROUP NUMBER  (IMAGES == 0, VOL. == 1 )
C 24    NSYMBT      # OF BYTES USED FOR SYMMETRY DATA (0 OR 80)
C                     PLUS ANY EXTRA HEADER BYTES
C 25-26 EXTRA       EXTRA USER DEFINED STORAGE SPACE.
C 27    ?           CURRENTLY: 'MRCO'
C 28    IVERSION    VERSION NUMBER (CURRENTLY: 20140)
C 29-41 EXTRA       EXTRA, USER DEFINED STORAGE SPACE.
C 42    IANGLE      FLAG FOR ANGLES PRESENT IN LOCATIONS: 43..48 
C 43    ANG1        PHI     (SPIDER & IMOD DEFINED LOCATION)
C 44    ANG2        THETA   (SPIDER & IMOD DEFINED LOCATION)
C 45    ANG3        PSI     (SPIDER & IMOD DEFINED LOCATION)
C 46    ANG4        PHI1    (SPIDER & IMOD DEFINED LOCATION)
C 47    ANG5        THETA1  (SPIDER & IMOD DEFINED LOCATION)
C 48    ANG6        PSI1    (SPIDER & IMOD DEFINED LOCATION)
C 49    EXTRA       EXTRA USER DEFINED STORAGE SPACE
C 50    ORX         X ORIGIN FOR TRANSFORMS
C 51    ORY         Y ORIGIN FOR TRANSFORMS
C 52    ORZ         Z ORIGIN FOR TRANSFORMS
C 53    MAP         CHARACTER STRING 'MAP ' TO IDENTIFY FILE TYPE   
C 54    MACHST      MACHINE STAMP                                   
C 55    RMS         RMS DEVIATION OF MAP FROM MEAN DENSITY          
C 56    NLABL       NUMBER OF LABELS BEING USED                     
C 57-256            LABEL(20,10) 10 80-CHARACTER TEXT LABELS
C 
C SYMMETRY RECORDS IF ANY, FOLLOW - STORED AS TEXT AS IN INTERNATIONAL
C TABLES, OPERATORS SEPARATED BY * AND GROUPED INTO 'LINES' OF 80
C CHARACTERS (IE. SYMMETRY OPERATORS DO NOT CROSS THE ENDS OF THE
C 80-CHARACTER 'LINES' AND THE 'LINES' DO NOT TERMINATE IN A *).
C 
C DATA RECORDS FOLLOW.
C
C NOTES IN VERSION 20140++ :
C        DMAX  < DMIN                        MAX & MIN UNDETERMINED
C        DMEAN < (SMALLER OF DMIN and DMAX)  DMEAN     UNDETERMINED
C        RMS   < 0.0                         RMS       UNDETERMINED
C
C     Bytes 213 and 214 contain 4 `nibbles' (half-bytes) indicating 
C     the representation of float, complex, integer and character 
C     datatypes. Bytes 215 and 216 are unused. The CCP4 library contains 
C     a general representation of datatypes, but in practice it is 
C     safe to use 0x44 0x44 0x00 0x00 for little endian machines, and 
C     0x11 0x11 0x00 0x00 for big endian machines. The CCP4 library 
C     uses this information to automatically byte-swap data if 
C     appropriate, when tranferring data files between machines.  
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE LISTHEDMRC(LUN,IRTFLG)

C     IMPLICIT NONE

#include "LUNMRCHDR.INC"
      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      INTEGER               :: LUN,IRTFLG
      INTEGER               :: NX,NY,NZ,ITYPE,MAXIM,NE
      LOGICAL               :: IS_MRC
      REAL                  :: FVAL
      INTEGER               :: IEND
      LOGICAL               :: BIGENDARCH,BIGENDED
      LOGICAL               :: BIGENDFILE,SAMEENDFILE
      CHARACTER(LEN=MAXNAM) :: FILNAM
      CHARACTER(LEN=800)    :: LABELS
      CHARACTER(LEN=4)      :: CSTR
      CHARACTER(LEN=1)      :: NULL = CHAR(0)
        
      REAL                  :: ANG1,ANG2,ANG3,ANG4,ANG5,ANG6
      REAL                  :: DMIN,DMAX,DMEAN,RMS
      INTEGER               :: NSYMBT,MRCMODE
      INTEGER               :: MACHST,ISPG
      INTEGER               :: MAPC,MAPR,MAPS,NXSTART
      INTEGER               :: NYSTART,NZSTART,MX,MY,MZ
      REAL                  :: CELLAX,CELLAY,CELLAZ
      REAL                  :: CELLBX,CELLBY,CELLBZ
      REAL                  :: PIXSIZX,PIXSIZY,PIXSIZZ
      REAL                  :: ORX,ORY,ORZ
      INTEGER               :: IANGLE,IMGSTATS,NLABL,IVERSION
      CHARACTER(LEN=4)      :: MAP 
      CHARACTER(LEN=4)      :: CAXIS,EXTTYP

      INTEGER               :: LNBLNKN    ! FUNCTION
      LOGICAL               :: ISMRCFILE  ! FUNCTION


C     OPEN INPUT FILE
      CALL OPFILEC(0,.TRUE.,FILNAM,LUN,'O',ITYPE,
     &             NX,NY,NZ,
     &             MAXIM,'MRC',.FALSE.,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF (.NOT. ISMRCFILE(FILNAM)) THEN
         CALL ERRT(101,'OPERATION ONLY WORKS ON MRC FILES',NE)
         IRTFLG = 1
         RETURN
      ENDIF

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ_MRC(LUN,MRC_HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      !write(3,*)' listhedmrc, header1:  ',mrc_header(1)
      !write(3,*)' listhedmrc, header2:  ',mrc_header(2)
      !write(3,*)' listhedmrc, header3:  ',mrc_header(3)
      !write(3,*)' listhedmrc, header4:  ',mrc_header(4)
      !write(3,*)' listhedmrc, header17: ',mrc_header(17)
      !write(3,*)' listhedmrc, header18: ',mrc_header(18)
      !write(3,*)' listhedmrc, header19: ',mrc_header(19)
      !call lungetsize_mrc(lun,nx,ny,nz,irtflg)
      !write(3,*)' listhedmrc3,  nx,ny,nz: ',nx,ny,nz !ok

      NX       = MRC_HEADER( 1)
      NY       = MRC_HEADER( 2)
      NZ       = MRC_HEADER( 3)

      MRCMODE  = MRC_HEADER( 4)

      NXSTART  = MRC_HEADER( 5)
      NYSTART  = MRC_HEADER( 6)
      NZSTART  = MRC_HEADER( 7)

      MX       = MRC_HEADER( 8)
      MY       = MRC_HEADER( 9)
      MZ       = MRC_HEADER(10)

      CELLAX   = TRANSFER(MRC_HEADER(11),FVAL)  ! PIXEL SIZE * NX
      CELLAY   = TRANSFER(MRC_HEADER(12),FVAL)  ! PIXEL SIZE * NY
      CELLAZ   = TRANSFER(MRC_HEADER(13),FVAL)  ! PIXEL SIZE * NZ

      PIXSIZX  = CELLAX / NX                    ! PIXEL SIZE IN X 
      PIXSIZY  = CELLAY / NY                    ! PIXEL SIZE IN Y 
      PIXSIZZ  = CELLAZ / NZ                    ! PIXEL SIZE IN Z  

      CELLBX   = TRANSFER(MRC_HEADER(14),FVAL)
      CELLBY   = TRANSFER(MRC_HEADER(15),FVAL)
      CELLBZ   = TRANSFER(MRC_HEADER(16),FVAL)

      MAPC     = MRC_HEADER(17)
      MAPR     = MRC_HEADER(18)
      MAPS     = MRC_HEADER(19)

      DMIN     = TRANSFER(MRC_HEADER(20),FVAL) ! MIN  DENSITY VALUE
      DMAX     = TRANSFER(MRC_HEADER(21),FVAL) ! MAX  DENSITY VALUE  
      DMEAN    = TRANSFER(MRC_HEADER(22),FVAL) ! MEAN DENSITY VALUE  

      ISPG     = MRC_HEADER(23) 
      NSYMBT   = MRC_HEADER(24) 
      IMGSTATS = MRC_HEADER(25)                     ! SPIDER DEFINED

      CAXIS    = TRANSFER(MRC_HEADER(26),CSTR(1:4)) ! SPIDER DEFINED

C     GET EXTTYP (NOW 'MRCO')
      EXTTYP   = TRANSFER(MRC_HEADER(27),CSTR(1:4))

C     GET VERSION NUMBER 
      IVERSION = MRC_HEADER(28) 

C     GET IANGLE
      IANGLE   = MRC_HEADER(42)

C     GET ANGLES (6 IMOD DEFINED LOCATIONS)
C                (PHI1,THETA1,PSI1),(PHI2,THETA2,PSI2)
      ANG1     = TRANSFER(MRC_HEADER(43),FVAL) 
      ANG2     = TRANSFER(MRC_HEADER(44),FVAL)
      ANG3     = TRANSFER(MRC_HEADER(45),FVAL)

      ANG4     = TRANSFER(MRC_HEADER(46),FVAL)
      ANG5     = TRANSFER(MRC_HEADER(47),FVAL)
      ANG6     = TRANSFER(MRC_HEADER(48),FVAL)

C     GET IMAGE ORIGINS (3  ? DEFINED LOCATIONS)
      ORX      = TRANSFER(MRC_HEADER(50),FVAL)
      ORY      = TRANSFER(MRC_HEADER(51),FVAL)
      ORZ      = TRANSFER(MRC_HEADER(52),FVAL)

C     GET MAP TYPE 
      MAP      = TRANSFER(MRC_HEADER(53),CSTR(1:4))

      MACHST   = MRC_HEADER(54) 
      RMS      = TRANSFER(MRC_HEADER(55),FVAL) ! DEVIATION FROM MEAN       
      NLABL    = MRC_HEADER(56)

      IF (NLABL > 0) THEN
C        GET LABELS
         CALL LUNGETLABELS_MRC(LUN,NLABL,LABELS,IRTFLG)
      ENDIF

C     GET CURRENT ARCHITECTURE ENDED-NESS
      BIGENDARCH = BIGENDED(0)

      SAMEENDFILE = ((MAPC == 1) .OR. (MAPR == 1) .OR. (MAPS == 1))

      BIGENDFILE = ((      SAMEENDFILE .AND.       BIGENDARCH) .OR.
     &              (.NOT. SAMEENDFILE .AND. .NOT. BIGENDARCH))

C     WRITE OUT CONVERSION INFORMATION
      WRITE(NOUT,*) ' '

      IF (BIGENDARCH .AND. BIGENDFILE ) THEN
        WRITE(NOUT,*)' On big ended architecture,',
     &               ' reading big ended file.'
      ELSE IF (BIGENDARCH ) THEN
        WRITE(NOUT,*)' On big ended architecture,',
     &               ' reading little ended file.'
      ELSE IF (BIGENDFILE) THEN
         WRITE(NOUT,*)' On little ended architecture,',
     &                ' reading big ended file.'
      ELSE
         WRITE(NOUT,*)' On little ended architecture,',
     &                ' Reading little ended file.'
      ENDIF


C     WRITE OUT HEADER INFORMATION
      WRITE(NOUT,*) ' *  = Derived from actual header '
      WRITE(NOUT,*) ' ** = From SPIDER specific header location  '

      IF     (MRCMODE == 0) THEN
        WRITE(NOUT,*) ' Data type (Mode) ..........................' //
     &               '    0: 8-bit signed integers, Range -128 --> 127'
      ELSEIF (MRCMODE == 1) THEN
        WRITE(NOUT,*) ' Data type .................................' //
     &                '     1: 16-bit signed integers'
      ELSEIF (MRCMODE == 2) THEN
        WRITE(NOUT,*) ' Data type ................................ ' //
     &                '     2: 32-bit reals '
      ELSEIF (MRCMODE == 3) THEN
        WRITE(NOUT,*) ' Data type  ................................' //
     &                '     3: Complex 16-bit integers'
      ELSEIF (MRCMODE == 4) THEN
        WRITE(NOUT,*) ' Data type  ................................' //
     &                '     4: Complex 32-bit reals'
      ELSEIF (MRCMODE == 6) THEN
        WRITE(NOUT,*) ' Data type  ................................' //
     &                '    6 : 16-bit unsigned integers'
      ENDIF

      WRITE(NOUT,1001) NX,NY,NZ
1001  FORMAT(2X,'Columns, rows, sections ..................  ',3(I7,1X))

      WRITE(NOUT,1003) NXSTART,NYSTART,NZSTART
1003  FORMAT(2X,'Starting columns, rows, sections ......... ',3I7)

      WRITE(NOUT,1004) MX,MY,MZ 
1004  FORMAT(2X,'Intervals (MX,MY,MZ) .....................  ',3I7)

      WRITE(NOUT,1005) PIXSIZX,PIXSIZY,PIXSIZZ 
1005  FORMAT(2X,'Pixel sizes (X,Y,Z) * .................... ',3F10.2)

      WRITE(NOUT,1026) CELLAX,CELLAY,CELLAZ 
1026  FORMAT(2X,'Cell sizes ............................... ',3F10.2)

      WRITE(NOUT,1006) CELLAB,CELLAB,CELLAB 
1006  FORMAT(2X,'Cell angles .............................. ',3F10.2)

      WRITE(NOUT,1007) MAPC,MAPR,MAPS 
1007  FORMAT(2X,'Fast, medium, slow axes .................. ',3I7)

      WRITE(NOUT,1008) DMIN,DMAX  
1008  FORMAT(2X,'Min & max  density ....................... ',2F23.11)

      WRITE(NOUT,1009) DMEAN,RMS
1009  FORMAT(2X,'Mean & RMS density ....................... ',2F23.11)

      WRITE(NOUT,1092) ANG1,ANG2,ANG3,ANG4,ANG5,ANG6
1092  FORMAT(2X,'Angles (Phi,Theta,Psi) ...............     ',6F8.2)

      WRITE(NOUT,1010) ORX,ORY,ORZ
1010  FORMAT(2X,'Origins .................................. ',3F10.2)

      WRITE(NOUT,1011) ISPG,NSYMBT 
1011  FORMAT(2X,'Space group, extra header bytes .......... ',2I7)

      WRITE(NOUT,1091) IMGSTATS 
1091  FORMAT(2X,'Image used for stats ** .................. ',2I7)
      WRITE(NOUT,1012) CAXIS 
1012  FORMAT(2X,'Data origin ** ..........................        ',A)

      WRITE(NOUT,1013) MACHST
1013  FORMAT(2X,'Machine stamp ............................ ',I12)

      WRITE(NOUT,1014) MAP
1014  FORMAT(2X,'Map type..................................       ',A)

      WRITE(NOUT,1015) IVERSION
1015  FORMAT(2X,'MRC version ..............................     ',I7)

      WRITE(NOUT,1016) EXTTYP 
1016  FORMAT(2X,'ExtType .....................................    ',A)

      WRITE(NOUT,1017) NLABL 
1017  FORMAT(2X,'Number of labels ......................... ',I7)

      IF (NLABL > 0) THEN
         IEND = lnblnkn(LABELS)

         WRITE(NOUT,1031)
1031     FORMAT('  Labels:')

         WRITE(NOUT,1032) LABELS(1:IEND)
1032     FORMAT(3X,100(A80))
      ENDIF

      WRITE(NOUT,*) ' '

      CLOSE(LUN)

      END
