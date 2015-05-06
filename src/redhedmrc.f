
C ++********************************************************************
C REDHEDMRC
C                              INTEL BYTE_ORDER     JUL 09 ARDEAN LEITH
C                                                                      *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
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
C                                                                      *
C REDHEDMRC(BUF,NX,NY,NZ,MODE,DMIN,DMAX,DMEAN,FLIP)
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C  REDHEDMRC
C
C 		MAP/IMAGE HEADER FORMAT
C 		-----------------------
C 
C 	LENGTH = 1024 BYTES, ORGANIZED AS 56 LONG WORDS FOLLOWED
C                BY SPACE FOR 10 80 BYTE TEXT LABELS.
C 
C  1	NX		# OF COLUMNS	(FASTEST CHANGING IN MAP)
C  2	NY		# OF ROWS
C  3	NZ		# OF SECTIONS 	(SLOWEST CHANGING IN MAP)
C  4	MODE		DATA TYPE
C 			  0 = IMAGE     STORED AS INTEGER*1
C 			  1 = IMAGE	STORED AS INTEGER*2
C 			  2 = IMAGE	STORED AS REALS
C 			  3 = TRANSFORM	STORED AS COMPLEX INTEGER*2
C 			  4 = TRANSFORM	STORED AS COMPLEX REALS
C 
C  5	NXSTART		NUMBER OF FIRST COLUMN  IN MAP (DEFAULT = 0)
C  6	NYSTART		NUMBER OF FIRST ROW     IN MAP      "
C  7	NZSTART		NUMBER OF FIRST SECTION IN MAP      "
C  8	MX		NUMBER OF INTERVALS ALONG X
C  9	MY		NUMBER OF INTERVALS ALONG Y
C 10	MZ		NUMBER OF INTERVALS ALONG Z
C 11	X LENGTH	CELL DIMENSIONS (ANGSTROMS)
C 12	Y LENGTH		     "
C 13	Z LENGTH		     "
C 14	ALPHA		CELL ANGLES     (DEGREES)
C 15	BETA		             "
C 16	GAMMA		             "
C 17	MAPC		WHICH AXIS CORRESPONDS TO COLUMNS  (1,2,3 FOR X,Y,Z)
C 18	MAPR		WHICH AXIS CORRESPONDS TO ROWS     (1,2,3 FOR X,Y,Z)
C 19	MAPS		WHICH AXIS CORRESPONDS TO SECTIONS (1,2,3 FOR X,Y,Z)
C 20	AMIN		MINIMUM DENSITY VALUE
C 21	AMAX		MAXIMUM DENSITY VALUE
C 22	AMEAN		MEAN    DENSITY VALUE    (AVERAGE)
C 23	ISPG		SPACE GROUP NUMBER       (0 FOR IMAGES)
C 24	NSYMBT		NUMBER OF BYTES USED FOR STORING SYMMETRY OPERATORS
C 25	EXTRA		EXTRA, USER DEFINED STORAGE SPACE. 29 WORDS MAX.
C .          "
C .          "
C .          "   (ALL SET TO ZERO BY DEFAULT)
C .          "
C 53	     "
C 54	XORIGIN		X ORIGIN
C 55	YORIGIN		Y ORIGIN
C 56	NLABL		NUMBER OF LABELS BEING USED
C 57-256  LABEL(20,10)	10  80 CHARACTER TEXT LABELS (IE. A4 FORMAT)
C 
C SYMMETRY RECORDS FOLLOW - IF ANY - STORED AS TEXT AS IN INTERNATIONAL
C TABLES, OPERATORS SEPARATED BY * AND GROUPED INTO 'LINES' OF 80
C CHARACTERS (IE. SYMMETRY OPERATORS DO NOT CROSS THE ENDS OF THE
C 80-CHARACTER 'LINES' AND THE 'LINES' DO NOT TERMINATE IN A *).
C 
C 		DATA RECORDS FOLLOW.
C 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       IRDHDR
C
C	Read the header on unit ISTREAM, Print the contents,
C	and return some parameters to the caller:
C
C	INXYZ		  : size of file Columns, Rows, Sections
C	MXYZ		  : # of intervals Columns, Rows, Sections
C	IMODE		  : Data storage mode (1-4)
C				0 = Image		Integer*1
C				1 = Image               Integer*2
C				2 = Image               Reals
C				3 = Fourier Transform   Integer*2
C				4 = Fourier Transform   Reals
C	DMIN,DMAX,DMEAN   : Min, Max, & Mean densities
C
C **********************************************************************

	SUBROUTINE REDHEDMRC(BUF,NX,NY,NZ,MODE,
     &                       DMIN,DMAX,DMEAN,FLIP)

        INCLUDE 'CMBLOCK.INC'


        REAL             :: BUF(*) 
        CHARACTER(LEN=4) :: CLABLS(20,10)

	EQUIVALENCE         (STUFF,ISTUFF)
        INTEGER          :: ISTUFF(31)
        REAL             :: STUFF(31)

        REAL             :: LABELS(20,10)
        INTEGER          :: LABLS(20,10) 

	INTEGER          :: NCRS(3),NCRST(3),NXYZ(3),INXYZ(3),MXYZ(3)
        INTEGER          :: MAPCRS(3)  
        REAL             :: CEL(6),DENMMM(3),ORIGXY(2)  
        CHARACTER(LEN=1) :: LXYZ(3)

        INTEGER          :: MAPC,MAPR,MAPS
        LOGICAL          :: FLIP
	
	DATA LXYZ/'X','Y','Z'/
	
C       DETERMINE CURRENT FILE ENDED-NESS JULY 09 al
        CALL CCPMVI(MAPC,BUF(17),1)
        CALL CCPMVI(MAPR,BUF(18),1)
        CALL CCPMVI(MAPS,BUF(19),1)
	FLIP = ((MAPC.NE.1) .AND. (MAPR.NE.1) .AND. (MAPS.NE.1))
        IF (FLIP) THEN
C          MUST FLIP BUF CONTENTS BEFORE USE
           !write(6,*)' Flipping mrc header bytes now'
           CALL FLIPBYTESI(BUF(1),56,IRTFLG)
           WRITE(NOUT,*)' FLIPPED MRC HEADER BYTES'
        ENDIF

	NB   = 4
	K    = 0
		
        I=1
	CALL  CCPMVI(NCRS, BUF(I), 3)
	I= I + 3 
	
        CALL CCPMVI(MODE, BUF(I), 1)
	I= I + 1
	  
        CALL CCPMVI(NCRST, BUF(I), 3)
	I= I + 3
	  
        CALL CCPMVI(NXYZ, BUF(I), 3)
	I= I + 3
	  
        CALL CCPMVI(CEL, BUF(I), 6)
	I= I + 6

        CALL CCPMVI(MAPCRS, BUF(I), 3)
        !write(6,*) ' mapcrs: ',MAPCRS(1),MAPCRS(2),MAPCRS(3),i

	I= I + 3
        CALL CCPMVI(DENMMM, BUF(I), 3)
	I= I + 3	
	
        CALL CCPMVI(STUFF, BUF(I), 31)
	I= I + 31
 
        CALL CCPMVI(ORIGXY, BUF(I), 2)
	I= I + 2
  
        CALL CCPMVI(NLAB, BUF(I), 1)
	I= I + 1   
 
        CALL CCPMVI(LABLS, BUF(I), 200)
	I= I + 200 -1
 
	CALL CCPMVI(INXYZ,NCRS(1),3) 
	CALL CCPMVI(MXYZ,NXYZ(1),3) 
	DMIN  = DENMMM(1)
	DMAX  = DENMMM(2)
	DMEAN = DENMMM(3)
	ISPG  = ISTUFF(1)
	NBS   = ISTUFF(2)
	IMODE = MODE
	IF (ISPG .GT. 0 .AND. NBS .GT. 0) THEN
	   NBSYM = NBS
	ELSE
	   ISPG  = 0
	   NBSYM = 0
	ENDIF

        CALL CCPIBY(CLABLS,LABLS,800)
	
        NX   = INXYZ(1)
        NY   = INXYZ(2)
        NZ = INXYZ(3)

C       WRITE OUT HEADER INFORMATION
	WRITE(NOUT,1011) (INXYZ(K),K=1,3) 
1011	FORMAT(/
     &  4X,'Number of columns, rows, sections ...........',3I6)

	WRITE(NOUT,1021)  IMODE  
1021	FORMAT(/
     &  4X,'Map mode ....................................',I5)

	WRITE(NOUT,1031) (NCRST(K),K=1,3)
1031	FORMAT(/
     &  4X,'Start points on columns, rows, sections .....',3I6)

	WRITE(NOUT,1041) (MXYZ(K),K=1,3)
1041	FORMAT(/
     &  4X,'Grid sampling in x, y, z ....................',3I6)

	WRITE(NOUT,1002) (CEL(K),K=1,6)
1002	FORMAT(/
     &  4X,'Cell axes ...................................',3G11.4/
     &  4X,'Cell angles .................................',3G11.4)
     
        WRITE(NOUT,1042) (LXYZ(MAPCRS(K)),K=1,3),(ORIGXY(K),K=1,2)
1042	FORMAT(/
     &  4X,'Fast, medium, slow axes .....................',3(4X,A1)/
     &	4X,'Origin in x, y ..............................',2G13.5)

	WRITE(NOUT,1003) DMIN,DMAX,DMEAN 
1003	FORMAT(/
     &  4X,'Minimum density .............................',G13.5/
     &  4X,'Maximum density .............................',G13.5/
     &  4X,'Mean density ................................',G13.5)

	WRITE(NOUT,1033) ISPG,NBSYM
1033	FORMAT(/
     &	4X,'Space group, # bytes symmetry ...............',2I6)

	WRITE(NOUT,1004) NLAB, ((CLABLS(I,K),I=1,20),K=1,NLAB)
1004	FORMAT(/
     &  4X,'Number of titles ............................',I5/
     & ' Titles :'/10(1X,20A4/))
     	
        END


C -----------------------------------------------------------

      SUBROUTINE CCPMVI (ARR1,ARR2,NUM)

C     PURPOSE:  ASSIGNS THE FIRST NUM WORDS OF ARR2 TO ARR1

      INTEGER  :: NUM
      REAL     :: ARR1(*),ARR2(*)
      INTEGER  :: J

      DO J=1,NUM
         ARR1(J) = ARR2(J)
      ENDDO

      END

 

C -----------------------------------------------------------

C  CCPIBY 
 
C COPY AN ARRAY OF INTEGERS INTO AN ARRAY OF UNSIGNED (OR UNSIGNED) BYTES. 
C NOTE: NO OVERFLOW CHECKING IS DONE.
C
C (MUST BE IMPLEMENTED IF CCPBYT FUNCTION RETURNS .TRUE.)
C
C PARAMETERS:
C    IBYT (O)   ARRAY RETURNING BYTE DATA (MAY BE AN INTEGER ARRAY FOR EXAMPLE
C               WITH DATA PACKED INTO ADJACANT BYTES
C      IA (I)   ARRAY HOLDING INTEGER VALUES
C      NB (I)   IF >0, THE NUMBER OF ELEMENTS TO BE COPIED TO UNSIGNED BYTES
C               IF <0, -THE NUMBER OF ELEMENTS TO BE COPIED TO SIGNED BYTES

      SUBROUTINE CCPIBY(IBYT,IA,NB)
      
      INTEGER     :: IA(*)
      INTEGER * 1 :: IBYT(*)
      INTEGER * 1 :: JBYT(4)
      EQUIVALENCE (JA,JBYT(1))
      LOGICAL     :: CALLED, LITEND
      EXTERNAL    :: LITEND
      INTEGER     :: IND
      SAVE        :: CALLED, IND

      DATA    CALLED/.FALSE./

      IF (.NOT. CALLED) THEN
        IF (LITEND(1)) THEN
          IND = 1
        ELSE
          IND = 4
        ENDIF
        CALLED = .TRUE.
      ENDIF
      
      NE = NB
      IF (NE > 0) THEN
         DO  I=1,NE
           JA = IA(I)
           IBYT(I) = JBYT(IND)
	 ENDDO
      ELSE
         NE = -NE
         DO  I=1,NE
           IBYT(I) = IA(I)
	 ENDDO
      ENDIF

      END


C -----------------------------------------------------------

         LOGICAL FUNCTION LITEND(IDUM)

C        CHECK ENDEDNESS, RETURNS TRUE IF LITTLE ENDIAN (VAX, FX2800,
C                                                   ULTRIX, CONVEX)
C                              FALSE IF BIG ENDIAN (IBM,IRIS,ESV)

         INTEGER       :: I, IDUM
         INTEGER * 1   :: B(4)
         EQUIVALENCE   (I,B(1))

C        INITIALISE B

         DO JDO=1,4
            B(JDO) = 0
         ENDDO

         I = 1

         IF (B(1) .NE. 0) THEN
            LITEND = .TRUE.
         ELSE
            LITEND = .FALSE.
         ENDIF

         END
	
