
C **********************************************************************
c
C SETHEDCCP4    ISSWAB ADDED                      JUL  02 ARDEAN LEITH
C               ANGSTROMS/PIXEL                   JAN  05 ARDEAN LEITH 
C               STARTING                          FEB  06 ARDEAN LEITH 
C               MRC PROMPTS                       MAY  12 ARDEAN LEITH 
C               NSYMBYT,ISPG                      JUN  15 ARDEAN LEITH 
C               ISSWAB RENAMED ISSWABT            JAN  16 ARDEAN LEITH
C               SETHEDCCP4_NEW ADDED              JUN  16 ARDEAN LEITH
C               SETSTAMP NOW ALWAYS LITTLE ENDED  JAN  18 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2018  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@wadsworth.org                                        *
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
C SETHEDCCP4(HEADBUF, NX, NY, NZ,
C           DMIN,DMAX,DMEAN,DSIG, SCALEX,SCALEY,SCALEZ,
C           IMODE,ISSWABT,NSYMBYT,NIMG,MZ,IRTFLG)
C
C PURPOSE: CREATE NEW CCP4 HEADER. ALL OF THE STANDARD IMAGE 
C          DEFAULTS ARE SET UP GIVEN THE REQUESTED INFORMATION. 
C          HEADER NOT WRITTEN.
C
C          NOTE: THE STARTING POINT FOR COLUMNS,ROWS,SECTIONS
C          ARE SET TO 0 BY DEFAULT.!!!!!!!!!
C
C PARAMETERS:
C       NX...    # OF INTERVALS COLUMNS, ROWS, SECTIONS
C       IMODE    DATA STORAGE MODE (1-4)
C                0 = IMAGE               8               INTEGERS
C                1 = IMAGE               16 BIT SIGNED   INTEGERS 
C                2 = IMAGE               32 BIT REALS
C                3 = FOURIER TRANSFORM   16 BIT INTEGER
C                4 = FOURIER TRANSFORM   32 BIT REALS      
C                6 = IMAGE               16 BIT UNSIGNED INTEGERS 
C
C NOTE:  SEE REDHEDCCP4 FOR SOME MORE DETAILS
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
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE SETHEDCCP4_NEW(HEADBUF, 
     &                        NX, NY, NZ, NIMG, IMODE,IVERSION,
     &                        ISTARTX,ISTARTY,ISTARTZ, 
     &                        ORIGX,ORIGY,ORIGZ, 
     &                        SCALEX,SCALEY,SCALEZ, 
     &                        DMIN,DMAX,DMEAN,DSIG, 
     &                        ISSWABT,NSYMBYT, IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'

        REAL               :: HEADBUF(*)
        INTEGER            :: NX,NY,NZ, IMODE,IVERSION
        INTEGER            :: ISTARTX,ISTARTY,ISTARTZ
        REAL               :: SCALEX,SCALEY,SCALEZ
        REAL               :: ORIGX,ORIGY,ORIGZ
        REAL               :: DMIN,DMAX,DMEAN,DSIG,SCALE
        LOGICAL            :: ISSWABT
        INTEGER            :: NSYMBYT,NIMG,ISPG,MZ,IRTFLG

        INTEGER            :: MX,MY,MACHST,NZMRC 
        INTEGER            :: NLABEL,I,INOW,NOT_USED
        CHARACTER(LEN=44)  :: LABEL1 
        CHARACTER(LEN=4)   :: CVAL
        CHARACTER(LEN=1)   :: BLANK = CHAR(32)
        real               :: cellax,cellay,cellaz

C       NOTE: COPY IMAGE DATA IN REAL*4 FORMAT. 

        WRITE(NOUT,*)' '
        
        NZMRC = NZ
        IF (IVERSION < 20140) THEN
           NZMRC = NIMG
        ENDIF

        CALL CCPMVI(HEADBUF(1),NX,1)
        CALL CCPMVI(HEADBUF(2),NY,1)
        CALL CCPMVI(HEADBUF(3),NZMRC,1)

        CALL CCPMVI(HEADBUF(4),IMODE,1)

C       SET START POINTS ON COLUMNS, ROWS, SECTIONS 
        CALL CCPMVI(HEADBUF(5),ISTARTX,1)
        CALL CCPMVI(HEADBUF(6),ISTARTY,1)
        CALL CCPMVI(HEADBUF(7),ISTARTZ,1)

C       GRID SAMPLING ON X, Y, Z (MX & MY ) 
C       FOR EM, NUMBER OF SECTIONS ALONG Z IN A VOLUME STACK
C       SO: MZ = NZ / NUMBER OF VOLUMES 

        MX = NX
        MY = NY
        IF (NZ == 1 .AND. NIMG == 1) THEN
C          IMAGE
           ISPG  = 0    
           MZ    = 1
           NZMRC = 1            
        ELSEIF (NZ == 1 .AND. NIMG > 1 ) THEN
C          IMAGE STACK
           ISPG  = 0    
           MZ    = NIMG
           NZMRC = NIMG        ! NUMBER OF IMAGES 
        ELSEIF (NZ > 1 .AND. NIMG == 1 ) THEN
C          VOLUME
           ISPG  = 1    
           MZ    = NZ
           NZMRC = NZ          ! NUMBER OF VOLUMES 
        ELSEIF (NZ > 1 .AND. NIMG > 1 ) THEN
C          STACK OF EM VOLUMES
           ISPG  = 401    
           MZ    = NZ * NIMG
           NZMRC = NIMG       ! NUMBER OF VOLUMES
        ELSE 
           CALL ERRT(102,'BAD STACK OR VOLUME PARAMETERS',NIMG)
           IRTFLG = 1
           RETURN
        ENDIF

        CALL CCPMVI(HEADBUF(8), MX,1)
        CALL CCPMVI(HEADBUF(9), MY,1)
        CALL CCPMVI(HEADBUF(10),MZ,1)

C       CELL SIZE (ANGSTROMS) (CELLAX..) 
        CELLAX      = SCALEX * NX
        CELLAY      = SCALEX * NX
        CELLAZ      = SCALEX * NX
        HEADBUF(11) = CELLAX
        HEADBUF(12) = CELLAY
        HEADBUF(13) = CELLAZ  

C       CELL ANGLES (DEGREES)  (CELLBX..) 
        HEADBUF(14) = 90.0
        HEADBUF(15) = 90.0
        HEADBUF(16) = 90.0

C       STORAGE: FAST, MEDIUM, SLOW AXES 
        CALL CCPMVI(HEADBUF(17),1,1)
        CALL CCPMVI(HEADBUF(18),2,1)
        CALL CCPMVI(HEADBUF(19),3,1)
 
C       MIN, MAX. & MEAN DENSITY 
        HEADBUF(20) = DMIN
        HEADBUF(21) = DMAX
        HEADBUF(22) = DMEAN

C       SPACE GROUP, # BYTES SYMMETRY 
        CALL CCPMVI(HEADBUF(23), ISPG,  1)
        CALL CCPMVI(HEADBUF(24), NSYMBYT,1)

C       ZERO THE EXTRA POSITIONS
        DO I = 25,49
           CALL CCPMVI(HEADBUF(I), 0,1)
        ENDDO

C       PUT IN EXTTYP : 'MRCO'
        CVAL = 'MRCO'
        CALL CCPMVC(HEADBUF(27),CVAL,ISSWABT)

C       VERSION NUMBER 
        CALL CCPMVI(HEADBUF(28), IVERSION,1)
                
        HEADBUF(50) = ORIGX
        HEADBUF(51) = ORIGY
        HEADBUF(52) = ORIGZ

C       PUT IN 'MAP'
        CVAL = 'MAP '
        CALL CCPMVC(HEADBUF(53),CVAL,ISSWABT)

C       SET MACHINE STAMP
C       MACHST = 286326784    ! FOR BIG ENDIAN DATA
C       MACHST = 4369         ! FOR LITTLE ENDIAN DATA
        MACHST = 16708        ! FOR LITTLE ENDIAN DATA
        CALL CCPMVI(HEADBUF(54),MACHST,1)

c       write(nout,*) ' map:',cval,headbuf(53)
c       write(nout,*) ' machst:',machst,headbuf(54)

C       SET RMS (ASSUMING THIS DEFINED SAME AS SIG IN SPIDER??)
        HEADBUF(55) = DSIG

C       SET NUMBER OF LABELS
        NLABEL = 1  
        CALL CCPMVI(HEADBUF(56), NLABEL,1)

C       ZERO ALL LABELS WITH BLANKS
        CVAL = BLANK // BLANK // BLANK // BLANK
        DO I = 57,256
           CALL CCPMVC(HEADBUF(I),CVAL,ISSWABT)
        ENDDO

C       ADD ONE LABEL
        LABEL1 = 'Converted from SPIDER to MRC using SPIDER '
        INOW   = 56
        DO I = 1,44,4
           CVAL = LABEL1(I:I+3)
           INOW = INOW + 1
           CALL CCPMVC(HEADBUF(INOW),CVAL,ISSWABT)
        ENDDO
        
        IRTFLG = 0
        END

#ifdef NEVER
C          WRITE OUT HEADER INFORMATION
           WRITE(6,1000) NX,NY,NZ,IMODE,
     &       ISTARTX,ISTARTY,ISTARTZ, MX,MY,MZ,
     &       CELLAX,CELLAY,CELLAZ, 
     &       DMIN,DMAX,DMEAN,DSIG,ORIGX,ORIGY,ORIGZ,ISPG,NSYMBYT,
     &       MACHST,IVERSION, CVAL

1000       FORMAT(
     &     2X,'Columns, rows, sections .................. ',3(I7,1X)/
     &     2X,'Mode ..................................... ',I6/
     &     2X,'Start points on columns, rows, sections .. ',3I7/
     &     2X,'Grid sampling on x, y, z ................. ',3I7/
     &     2X,'Cell axes ................................ ',3F10.2/
     &     2X,'Minimum density .......................... ',F25.12/
     &     2X,'Maximum density .......................... ',F25.12/
     &     2X,'Mean density ............................. ',F25.12/
     &     2X,'RMS deviation ............................ ',F25.12/
     &     2X,'Origins .................................. ',3F10.2/
     &     2X,'Space group, # bytes symmetry ............ ',2I7/
     &     2X,'Machine stamp ............................ ',I12/
     &     2X,'Version ..................................    ',I7/
     &     2X,'ExtType .....................................   ',A)

#endif




C ******************************** SETHEDCCP4 *****************************


        SUBROUTINE SETHEDCCP4(HEADBUF, NX, NY, NZ,
     &                        DMIN,DMAX,DMEAN,DSIG, 
     &                        SCALEX,SCALEY,SCALEZ, IMODE,
     &                        ISSWABT,NSYMBYT,NIMG,MZ,IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'

        REAL               :: HEADBUF(*)
        INTEGER            :: NX,NY,NZ, IMODE
        REAL               :: DMIN,DMAX,DMEAN,DSIG,SCALE
        LOGICAL            :: ISSWABT
        INTEGER            :: NSYMBYT,NIMG,MZ,IRTFLG

        INTEGER            :: ISTARTX,ISTARTY,ISTARTZ,NOT_USED
        INTEGER            :: MX,MY,IVERSION,MACHST,ISPG,NZMRC 
        REAL               :: SCALEX,SCALEY,SCALEZ
        REAL               :: ORIGX,ORIGY,ORIGZ
        INTEGER            :: NLABEL,I,INOW
        CHARACTER(LEN=44)  :: LABEL1 
        CHARACTER(LEN=4)   :: CVAL
        CHARACTER(LEN=1)   :: BLANK = CHAR(32)

c       NOTE: COPY IMAGE DATA IN REAL*4 FORMAT. 
        
        CALL CCPMVI(HEADBUF(1),NX,1)
        CALL CCPMVI(HEADBUF(2),NY,1)
C       NZ = MZ * NUMBER OF VOLUMES
        CALL CCPMVI(HEADBUF(3),NZ,1)

C       COPY IMAGE DATA IN REAL*4 FORMAT.       
        CALL CCPMVI(HEADBUF(4),IMODE,1)

C       SET START POINTS ON COLUMNS, ROWS, SECTIONS

        WRITE(NOUT,*)' '
        IF (NZ > 1) THEN
C          VOLUME
           WRITE(NOUT,*)' STARTING (X,Y,Z) DEFAULT VALUE: (-(NX/2)+ 1,',
     &                 ' -(NY/2)+1), -(NZ/2)+1'
           WRITE(NOUT,*)' (+1 ADDED ONLY IF LENGTH IS ODD) ',
     &               ' USE <CR> FOR DEFAULT VALUE'
                
           ISTARTX = -(NX/2) + (MOD(NX,2))
           ISTARTY = -(NY/2) + (MOD(NY,2))
           ISTARTZ = -(NZ/2) + (MOD(NZ,2))
 
           CALL RDPRI3S(ISTARTX,ISTARTY,ISTARTZ,NOT_USED,
     &                 'STARTING X, Y, & Z FOR MRC DATA',IRTFLG)

        ELSE
C          IMAGE
           ISTARTZ = -(NZ/2) + (MOD(NZ,2))
           WRITE(NOUT,*)' STARTING (X & Y) DEFAULT VALUE: (-(NX/2)+ 1,',
     &                  ' -(NY/2)+1)'
           WRITE(NOUT,*)' (+1 ADDED ONLY IF LENGTH IS ODD) ',
     &               ' USE <CR> FOR DEFAULT VALUE'
                
           ISTARTX = -(NX/2) + (MOD(NX,2))
           ISTARTY = -(NY/2) + (MOD(NY,2))
           ISTARTZ = -(NZ/2) + (MOD(NZ,2))
 
           CALL RDPRI2S(ISTARTX,ISTARTY,NOT_USED,
     &                 'STARTING X & Y FOR MRC DATA',IRTFLG)
        ENDIF

        CALL CCPMVI(HEADBUF(5),ISTARTX,1)
        CALL CCPMVI(HEADBUF(6),ISTARTY,1)
        CALL CCPMVI(HEADBUF(7),ISTARTZ,1)

C       GRID SAMPLING ON X, Y, Z (MX..) SAME AS NX, NY, NZ
        MX = NX
        MY = NY
C       FOR EM, NUMBER OF SECTIONS ALONG Z IN A VOLUME STACK
C       SO: MZ = NZ / NUMBER OF VOLUMES 

        IF (NZ == 1 .AND. NIMG == 1) THEN
C          IMAGE
           ISPG  = 0    
           MZ    = 1
           NZMRC = 1            
        ELSEIF (NZ == 1 .AND. NIMG > 1 ) THEN
C          IMAGE STACK
           ISPG  = 0    
           MZ    = NIMG
           NZMRC = NIMG        ! NUMBER OF IMAGES 
        ELSEIF (NZ > 1 .AND. NIMG == 1 ) THEN
C          VOLUME
           ISPG  = 1    
           MZ    = NZ
           NZMRC = NZ          ! NUMBER OF VOLUMES 
        ELSEIF (NZ > 1 .AND. NIMG > 1 ) THEN
C          STACK OF EM VOLUMES
           ISPG  = 401    
           MZ    = NZ * NIMG
           NZMRC = NZ          ! NUMBER OF VOLUMES
        ELSE 
           CALL ERRT(102,'BAD STACK OR VOLUME PARAMETERS',NIMG)
           IRTFLG = 1
           RETURN
        ENDIF


        CALL CCPMVI(HEADBUF(8), NX,1)
        CALL CCPMVI(HEADBUF(9), NY,1)
        CALL CCPMVI(HEADBUF(10),NZMRC,1)

C       CELL AXES (CELLAX..) 
        IF (SCALEX > 0) THEN
           WRITE(NOUT,*) ' ANGSTROMS/PIXEL FOR ALL AXES (HEADER ',
     &                   ' LOCATION=21): ',SCALEX
        ELSE
           SCALEX = 1.0
           SCALEY = 1.0
           SCALEZ = 1.0
        ENDIF

        IF (NZ > 1) THEN
           CALL RDPRM3S(SCALEX,SCALEY,SCALEZ,NOT_USED,
     &              'ANGSTROMS/PIXEL FOR  X, Y, & Z AXIS',IRTFLG)
        ELSE
           CALL RDPRM2S(SCALEX,SCALEY,NOT_USED,
     &              'ANGSTROMS/PIXEL FOR  X &  Y AXES',IRTFLG)
        ENDIF
        IF (IRTFLG .NE. 0) RETURN

        HEADBUF(11) = SCALEX * NX
        HEADBUF(12) = SCALEY * NY
        HEADBUF(13) = SCALEZ * NZ

C       CELL ANGLES (CELLBX..) 
        HEADBUF(14) = 90.0
        HEADBUF(15) = 90.0
        HEADBUF(16) = 90.0

C       FAST, MEDIUM, SLOW AXES 
        CALL CCPMVI(HEADBUF(17),1,1)
        CALL CCPMVI(HEADBUF(18),2,1)
        CALL CCPMVI(HEADBUF(19),3,1)
 
C       MIN, MAX. & MEAN DENSITY 
        HEADBUF(20) = DMIN
        HEADBUF(21) = DMAX
        HEADBUF(22) = DMEAN

C       SPACE GROUP, # BYTES SYMMETRY 
        CALL CCPMVI(HEADBUF(23), ISPG,  1)
        CALL CCPMVI(HEADBUF(24), NSYMBYT,1)

C       ZERO THE EXTRA POSITIONS
        DO I = 25,49
           CALL CCPMVI(HEADBUF(I), 0,1)
        ENDDO

C       PUT IN 'MRCO'
        CVAL = 'MRCO'
        CALL CCPMVC(HEADBUF(27),CVAL,ISSWABT)

C       VERSION NUMBER 
        IVERSION = 20140
        CALL CCPMVI(HEADBUF(28), IVERSION,1)

C       ORIGIN ON X,Y & Z AXIS 
        ORIGX = (NX/2) + (MOD(NX,2))
        ORIGY = (NY/2) + (MOD(NY,2))
        ORIGZ = (NZ/2) + (MOD(NZ,2))

        WRITE(NOUT,*)' '
        IF (NZ < 2) THEN
          WRITE(NOUT,*)' ORIGIN (X,Y) DEFAULT VALUE:((NX/2)+ 1, ',
     &                 '(NY/2)+1)'
          WRITE(NOUT,*)'   (+1 ADDED ONLY IF AXIS LENGTH IS ODD)',
     &                 ' USE <CR> FOR DEFAULT VALUE'
          CALL RDPRM2S(ORIGX,ORIGY,NOT_USED,
     &               'X, & Y ORIGIN FOR MRC DATA',IRTFLG)
        ELSE 
          WRITE(NOUT,*)' ORIGIN (X,Y,Z) DEFAULT VALUE:((NX/2)+ 1, ',
     &                 '(NY/2)+1), (NZ/2)+1'
          WRITE(NOUT,*)'   (+1 ADDED ONLY IF AXIS LENGTH IS ODD)',
     &                 ' USE <CR> FOR DEFAULT VALUE'
          CALL RDPRM3S(ORIGX,ORIGY,ORIGZ,NOT_USED,
     &               'X, Y, & Z ORIGIN FOR MRC DATA',IRTFLG)
        ENDIF
        IF (IRTFLG .NE. 0) RETURN
                
        HEADBUF(50) = ORIGX
        HEADBUF(51) = ORIGY
        HEADBUF(52) = ORIGZ

C       PUT IN 'MAP'
        CVAL = 'MAP '
        CALL CCPMVC(HEADBUF(53),CVAL,ISSWABT)

C       SET MACHINE STAMP
C       MACHST = 286326784    ! FOR BIG ENDIAN DATA
C       MACHST = 4369         ! FOR LITTLE ENDIAN DATA
        MACHST = 16708        ! FOR LITTLE ENDIAN DATA
        CALL CCPMVI(HEADBUF(54),MACHST,1)
 
c        write(nout,*) ' map:',cval,headbuf(53)
c        write(nout,*) ' machst:',machst,headbuf(54)

C       SET RMS (ASSUMING THIS DEFINE SAME AS SIG IN SPIDER??)
        HEADBUF(55) = DSIG

C       SET NUMBER OF LABELS
        NLABEL = 1  
        CALL CCPMVI(HEADBUF(56), NLABEL,1)

C       ZERO ALL LABELS WITH BLANKS
        CVAL = BLANK // BLANK // BLANK // BLANK
        DO I = 57,256
           CALL CCPMVC(HEADBUF(I),CVAL,ISSWABT)
        ENDDO

C       NOW GO BACK AND ADD IN ONE LABEL
        LABEL1 = 'Converted from SPIDER to MRC using SPIDER '
        INOW   = 56
        DO I = 1,44,4
           CVAL = LABEL1(I:I+3)
           INOW = INOW + 1
           CALL CCPMVC(HEADBUF(INOW),CVAL,ISSWABT)
        ENDDO
        
        IRTFLG = 0

        END





C ------------------------- CCPMVC --------------------------

C     PURPOSE:     PACK 4 CHAR INTO A INTEGER OR FLOAT

C     PARAMETERS:
C     I1ARRAY      ARRAY TO WHICH CHARS ARE TO BE COPIED       RET.
C     CIN          CHAR STRING TO BE COPIED INTO I1ARRAY       SENT
C     REVERSE      FLAG TO INVERT STRING                       SENT

      SUBROUTINE CCPMVC(I1ARRAY,CIN,REVERSE)

      CHARACTER(LEN=4) :: CIN
      INTEGER * 1      :: I1ARRAY(4)
      LOGICAL          :: REVERSE

      IF (REVERSE) THEN
         DO I = 1,4
            I1ARRAY(I) = ICHAR(CIN(5-I:5-I))
         ENDDO
      ELSE
         DO I = 1,4
            I1ARRAY(I) = ICHAR(CIN(I:I))
         ENDDO
      ENDIF

      END

