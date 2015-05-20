
C **********************************************************************
C                                                                      
C  GETHEDCCP4                                                                      *
C                INTEL BYTE_ORDER                  JUL 09 ARDEAN LEITH
C                FORMATTING                        SEP 14 ARDEAN LEITH
C                FORMATTING                        JAN 15 ARDEAN LEITH
C                                                                      
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2015  Health Research Inc.,                         *
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
C  GETHEDCCP4(HEADBUF,NX,NY,NZ,IMODE,DMIN,DMAX,
C             DMEAN,RMS,NSYMBT,ISSWABT,FLIP,MACHST,IRTFLG)
C                                                                      
C  PURPOSE:     DECODE CCP4 (MRC IMAGE 2000) HEADER                                        *
C                                            
C  PARAMETERS: 
C  
C  MAP/IMAGE HEADER FORMAT
C 
C 	LENGTH = 1024 BYTES, ORGANIZED AS 56 LONG WORDS FOLLOWED
C                BY SPACE FOR 10 80 BYTE TEXT LABELS.
C 
C  1	NX		# OF COLUMNS	(FASTEST CHANGING IN MAP)
C  2	NY		# OF ROWS
C  3	NZ		# OF SECTIONS 	(SLOWEST CHANGING IN MAP)
C  4	MODE		DATA TYPE
C			0	IMAGE : SIGNED 8-BIT BYTES RANGE -128 	
C					TO 127				
C			1	IMAGE : 16-BIT HALFWORDS		
C			2	IMAGE : 32-BIT REALS			
C			3	TRANSFORM : COMPLEX 16-BIT INTEGERS	
C			4	TRANSFORM : COMPLEX 32-BIT REALS	
C  5	NXSTART		NUMBER OF FIRST COLUMN  IN MAP 
C  6	NYSTART		NUMBER OF FIRST ROW     IN MAP       
C  7	NZSTART		NUMBER OF FIRST SECTION IN MAP       
C  8	MX		NUMBER OF INTERVALS ALONG X
C  9	MY		NUMBER OF INTERVALS ALONG Y
C 10	MZ		NUMBER OF INTERVALS ALONG Z
C 11-13	CELLA           CELL DIMENSIONS IN ANGSTROMS
C 14-16	CELLB           CELL ANGLES IN DEGREES				
C 17	MAPC		AXIS CORRESPONDING TO COLUMNS  (1,2,3 FOR X,Y,Z)
C 18	MAPR		AXIS CORRESPONDING TO ROWS     (1,2,3 FOR X,Y,Z)
C 19	MAPS		AXIS CORRESPONDING TO SECTIONS (1,2,3 FOR X,Y,Z)
C 20	DMIN		MINIMUM DENSITY VALUE
C 21	DMAX		MAXIMUM DENSITY VALUE
C 22	DMEAN		MEAN    DENSITY VALUE    
C 23	ISPG		SPACE GROUP NUMBER 0 OR 1 (DEFAULT=0)
C 24	NSYMBT		NUMBER OF BYTES USED FOR SYMMETRY DATA (0 OR 80)
C 25-49 EXTRA           EXTRA, USER DEFINED STORAGE SPACE. 29 WORDS MAX.
C 50-52	ORIGIN          ORIGIN IN X,Y,Z USED FOR TRANSFORMS		
C 53	MAP	        CHARACTER STRING 'MAP ' TO IDENTIFY FILE TYPE	
C 54	MACHST          MACHINE STAMP					
C 55	RMS	        RMS DEVIATION OF MAP FROM MEAN DENSITY		
C 56	NLABL	        NUMBER OF LABELS BEING USED			
C 57-256	        LABEL(20,10) 10 80-CHARACTER TEXT LABELS
C 
C SYMMETRY RECORDS FOLLOW - IF ANY - STORED AS TEXT AS IN INTERNATIONAL
C TABLES, OPERATORS SEPARATED BY * AND GROUPED INTO 'LINES' OF 80
C CHARACTERS (IE. SYMMETRY OPERATORS DO NOT CROSS THE ENDS OF THE
C 80-CHARACTER 'LINES' AND THE 'LINES' DO NOT TERMINATE IN A *).
C 
C DATA RECORDS FOLLOW.
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	SUBROUTINE GETHEDCCP4(HEADBUF,NX,NY,NZ,IMODE,DMIN,DMAX,
     &            DMEAN,RMS,NSYMBT,ISSWABT,FLIP,MACHST,IRTFLG)

        INCLUDE 'CMBLOCK.INC'

        REAL                  :: HEADBUF(*)

	CHARACTER(LEN=800)    :: CLABLS
        CHARACTER(LEN=4)      :: MAP
        LOGICAL               :: BIGENDARCH,BIGENDED
        LOGICAL               :: BIGENDFILE,SAMEENDFILE
        LOGICAL               :: FLIP,ISSWABT

        INTEGER               :: LENMAP
	CHARACTER(LEN=1)      :: LXYZ(3)

	DATA LXYZ/'X','Y','Z'/

C       GET CURRENT ARCHITECTURE ENDED-NESS
        BIGENDARCH = BIGENDED(0)

C       DETERMINE CURRENT FILE ENDED-NESS

        CALL MVNFLIP(HEADBUF(17), MAPC, .FALSE.)
        CALL MVNFLIP(HEADBUF(18), MAPR, .FALSE.)
        CALL MVNFLIP(HEADBUF(19), MAPS, .FALSE.)

        !write(6,*) ' mapc,mapr,maps: ',mapc,mapr,maps

	SAMEENDFILE = ((MAPC.EQ.1) .OR. (MAPR.EQ.1) .OR. (MAPS.EQ.1))

        BIGENDFILE = ((      SAMEENDFILE .AND.       BIGENDARCH) .OR.
     &                (.NOT. SAMEENDFILE .AND. .NOT. BIGENDARCH))

C       IF FILE ENDEDNESS DIFFERS FROM THIS MACHINES MUST FLIP BYTES
        FLIP = .NOT. SAMEENDFILE

C       WRITE OUT CONVERSION INFORMATION
        IF (VERBOSE) THEN    	
           IF (BIGENDARCH) THEN
	      WRITE(NOUT,*)' On big ended architecture'
           ELSE
	      WRITE(NOUT,*)' On little ended architecture'
           ENDIF

           IF (BIGENDFILE) THEN
	      WRITE(NOUT,*)' Reading big ended file'
           ELSE
	      WRITE(NOUT,*)' Reading little ended file'
           ENDIF

           IF (.NOT. ISSWABT) THEN
	      WRITE(NOUT,*)' SPIDER I/O Native byte order'
           ELSE
	      WRITE(NOUT,*)' SPIDER I/O Non-native byte order'
           ENDIF

           IF (FLIP) THEN
	      WRITE(NOUT,*)' Flipping byte order'
           ENDIF
        ENDIF

C       GET MAP TYPE FROM INPUT FILE
        CALL MVNREV(HEADBUF(53),MAP,ISSWABT)
          
        LENMAP = lnblnkn(MAP)
        IF (LENMAP <= 0) THEN
C           OLD STYLE, MRC MAP OR UNKNOWN FILE TYPE
            WRITE(NOUT,*) ' *** BLANK MAP STRING IN CCP4 FILE'

        ELSEIF (MAP .NE. 'MAP ') THEN
C           OLD STYLE, MRC MAP OR UNKNOWN FILE TYPE
            WRITE(NOUT,*) ' *** BAD MAP STRING IN CCP4 FILE: ',MAP
        ENDIF

	CALL MVNFLIP(HEADBUF( 1), NX,     FLIP) 
	CALL MVNFLIP(HEADBUF( 2), NY,     FLIP)
	CALL MVNFLIP(HEADBUF( 3), NZ,     FLIP)
	CALL MVNFLIP(HEADBUF( 4), IMODE,  FLIP)

        CALL MVNFLIP(HEADBUF( 5), NXSTART,FLIP)
        CALL MVNFLIP(HEADBUF( 6), NYSTART,FLIP)
        CALL MVNFLIP(HEADBUF( 7), NZSTART,FLIP)
        CALL MVNFLIP(HEADBUF( 8), MX,     FLIP)
        CALL MVNFLIP(HEADBUF( 9), MY,     FLIP)
        CALL MVNFLIP(HEADBUF(10), MZ,     FLIP)

        IF (FLIP) CALL FLIPBYTESI(HEADBUF(11),6,IRTFLG)
        CELLAX  = HEADBUF(11)
        CELLAY  = HEADBUF(12)
        CELLAZ  = HEADBUF(13)
        CELLBX  = HEADBUF(14)
        CELLBY  = HEADBUF(15)
        CELLBZ  = HEADBUF(16)

        CALL MVNFLIP(HEADBUF(17), MAPC,  FLIP)
        CALL MVNFLIP(HEADBUF(18), MAPR,  FLIP)
        CALL MVNFLIP(HEADBUF(19), MAPS,  FLIP)

        IF (FLIP) CALL FLIPBYTESI(HEADBUF(20),3,IRTFLG)
        DMIN    = HEADBUF(20)
        DMAX    = HEADBUF(21)
        DMEAN   = HEADBUF(22)

        CALL MVNFLIP(HEADBUF(23), ISPG,  FLIP)
        CALL MVNFLIP(HEADBUF(24), NSYMBT,FLIP)

        IF (FLIP) CALL FLIPBYTESI(HEADBUF(50),3,IRTFLG)
        ORX     = HEADBUF(50)
        ORY     = HEADBUF(51)
        ORZ     = HEADBUF(52)

        CALL MVNFLIP(HEADBUF(54), MACHST, FLIP)
        IF (FLIP) CALL FLIPBYTESI(HEADBUF(55),1,IRTFLG)
        RMS     = HEADBUF(55)
        CALL MVNFLIP(HEADBUF(56), NLABL, FLIP)

        IF (NLABL > 0) THEN
C          GET LABELS
           INOW = 1
           DO I = 57,56 + NLABL * 20
              CALL MVNREV(HEADBUF(I),CLABLS(INOW:INOW+3),ISSWABT)
              INOW = INOW + 4
	   ENDDO
        ENDIF

        IF (VERBOSE) THEN    	
C          WRITE OUT HEADER INFORMATION
	   WRITE(NOUT,1000) NX,NY,NZ,IMODE,
     &       NXSTART,NYSTART,NZSTART,MX,MY,MZ,
     &	     CELLAX,CELLAY,CELLAZ,CELLBX,CELLBY,CELLBZ,
     &       LXYZ(MAPC),LXYZ(MAPR),LXYZ(MAPS),
     &	     DMIN,DMAX,DMEAN,RMS,ORX,ORY,ORZ,ISPG,NSYMBT,
     &       MACHST,MAP,NLABL

1000	   FORMAT(/
     &     2X,'Number of columns, rows, sections ........ ',3(I7,1X)/
     &     2X,'Pixel mode ............................... ',I6/
     &     2X,'Start points on columns, rows, sections .. ',3I7/
     &     2X,'Grid sampling on x, y, z ................. ',3I7/
     &     2X,'Cell axes ................................ ',3F10.2/
     &     2X,'Cell angles .............................. ',3F10.2/
     &     2X,'Fast, medium, slow axes .................. ',3(4X,A)/
     &     2X,'Minimum density .......................... ',F25.12/
     &     2X,'Maximum density .......................... ',F25.12/
     &     2X,'Mean density ............................. ',F25.12/
     &     2X,'RMS deviation ............................ ',F25.12/
     &     2X,'Origins .................................. ',3F10.2/
     &	   2X,'Space group, # bytes symmetry ............ ',2I7/
     &	   2X,'Machine stamp ............................ ',I12/
     &     2X,'Map ......................................    ',A/
     &     2X,'Number of labels ......................... ',I7)

	   IF (NLABL > 0) THEN
              WRITE(NOUT,1001)
1001          FORMAT('  Labels:')
              WRITE(NOUT,1002) CLABLS(1:NLABL * 80)
1002          FORMAT(3X,100(A80))
           ENDIF
           WRITE(NOUT,*) ' '

        ENDIF

        END

C --------------------------- MVNREV -------------------------------

C     PURPOSE: COPY AN ARRAY OF INTEGER * 1 INTO A CHAR. STRING 
C
C     PARAMETERS:
C     IA        ARRAY HOLDING 4 INTEGER * 1 or 1 FLOAT          SENT
C     COUT      CHAR STRING                                     RET.
C     REVERSE   LOGICAL FLAG TO FLIP ORDER OF ELEMENTS          SENT

      SUBROUTINE MVNREV(I1IN,COUT,REVERSE)
      
      CHARACTER(LEN=4) :: COUT
      INTEGER * 1      :: I1IN(4)
      LOGICAL          :: REVERSE

      IF (REVERSE) THEN
         DO I = 1,4
           COUT(I:I) = CHAR(I1IN(5-I))
         ENDDO
      ELSE
         DO I = 1,4
           COUT(I:I) = CHAR(I1IN(I))
         ENDDO
      ENDIF

      END


C ------------------------ MVNFLIP -------------------------------

      SUBROUTINE MVNFLIP(I1IN,I1OUT,FLIP)

C     ASSIGNS I1IN TO I1OUT AND FLIPS BYTES WITHIN WORDS IF REQUESTED

      LOGICAL     :: FLIP
      INTEGER * 1 :: I1IN(4),I1OUT(4)

      IF (FLIP) THEN
         I1OUT(1) = I1IN(4)
         I1OUT(2) = I1IN(3)
         I1OUT(3) = I1IN(2)
         I1OUT(4) = I1IN(1)
      ELSE
         I1OUT(1) = I1IN(1)
         I1OUT(2) = I1IN(2)
         I1OUT(3) = I1IN(3)
         I1OUT(4) = I1IN(4)
      ENDIF

      END

 
