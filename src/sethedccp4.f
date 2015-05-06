
C *********************************************************************
c
C SETHEDCCP4                   ISSWAB ADDED        JUL  02 ARDEAN LEITH
C                              ANGSTROMS/PIXEL     JAN  05 ARDEAN LEITH 
C                              STARTING            FEB  06 ARDEAN LEITH 
C                              MRC PROMPTS         MAY  12 ARDEAN LEITH 
C *********************************************************************
c
C SETHEDCCP4(HEADBUF, NX, NY, NZ,
C            DMIN,DMAX,DMEAN,DSIG, IMODE,ISSWAB,IRTFLG)
C
C PURPOSE: CREATE NEW CCP4 HEADER. ALL OF THE STANDARD IMAGE 
C          DEFAULTS ARE SET UP GIVEN THE REQUESTED INFORMATION. 
C          HEADER NOT WRITTEN.
C
C          NOTE: THE STARTING POINT FOR COLUMNS,ROWS,SECTIONS
C	   ARE SET TO 0 BY DEFAULT.!!!!!!!!!
C
C PARAMETERS:
C	NX...    # OF INTERVALS COLUMNS, ROWS, SECTIONS
C	IMODE    DATA STORAGE MODE (1-4)
C				0 = IMAGE		INTEGER*1
C				1 = IMAGE               INTEGER*2
C				2 = IMAGE               REALS
C				3 = FOURIER TRANSFORM   INTEGER*2
C				4 = FOURIER TRANSFORM   REALS
C
C NOTE:  SEE REDHEDCCP4 FOR DETAILS
C
C *********************************************************************

	SUBROUTINE SETHEDCCP4(HEADBUF, NX, NY, NZ,
     &                        DMIN,DMAX,DMEAN,DSIG, SCALE,IMODE,
     &                        ISSWAB,IRTFLG)

        INCLUDE 'CMBLOCK.INC'

        REAL               :: HEADBUF(*)

        CHARACTER(LEN=44)  :: LABEL1 
        CHARACTER(LEN=4)   :: CVAL
        CHARACTER(LEN=1)   :: BLANK = CHAR(32)
	LOGICAL            :: ISSWAB

c	NOTE: COPY IMAGE DATA IN REAL*4 FORMAT.	
	
        CALL CCPMVI(HEADBUF(1),NX,1)
        CALL CCPMVI(HEADBUF(2),NY,1)
        CALL CCPMVI(HEADBUF(3),NZ,1)

C	COPY IMAGE DATA IN REAL*4 FORMAT.	
        CALL CCPMVI(HEADBUF(4),IMODE,1)

C       SET START POINTS ON COLUMNS, ROWS, SECTIONS
 	WRITE(NOUT,*)' '
 	WRITE(NOUT,*)' STARTING (X,Y,Z) DEFAULT VALUE:(-(NX/2)+ 1, ',
     &               '-(NY/2)+1), -(NZ/2)+1'
 	WRITE(NOUT,*)'   (+1 ADDED ONLY IF LENGTH IS ODD)'
	WRITE(NOUT,*)' USE <CR> FOR DEFAULT VALUE'
	
        ISTARTX = -(NX/2) + (MOD(NX,2))
        ISTARTY = -(NY/2) + (MOD(NY,2))
        ISTARTZ = -(NZ/2) + (MOD(NZ,2))
 
	CALL RDPRI3S(ISTARTX,ISTARTY,ISTARTZ,NOT_USED,
     &              'STARTING X, Y, & Z FOR MRC DATA',IRTFLG)
        CALL CCPMVI(HEADBUF(5),ISTARTX,1)
        CALL CCPMVI(HEADBUF(6),ISTARTY,1)
        CALL CCPMVI(HEADBUF(7),ISTARTZ,1)

C       GRID SAMPLING ON X, Y, Z (MX..) SAME AS NX, NY, NZ
        CALL CCPMVI(HEADBUF(8), NX,1)
        CALL CCPMVI(HEADBUF(9), NY,1)
        CALL CCPMVI(HEADBUF(10),NZ,1)

C       CELL AXES (CELLAX..) 
        IF (SCALE .NE. 0) THEN
           WRITE(NOUT,*) ' ANGSTROMS/PIXEL FOR ALL AXES (HEADER ',
     &                   ' LOCATION=21): ',SCALE
           SCALEX = SCALE
           SCALEY = SCALE
           SCALEZ = SCALE
        ELSE
           SCALEX = 1.0
           SCALEY = 1.0
           SCALEZ = 1.0
        ENDIF

	CALL RDPRM3S(SCALEX,SCALEY,SCALEZ,NOT_USED,
     &              'ANGSTROMS/PIXEL FOR  X, Y, & Z AXIS',IRTFLG)
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
        ISPG   = 1
        NSYMBT = 0
        CALL CCPMVI(HEADBUF(23), ISPG,1)
        CALL CCPMVI(HEADBUF(24), NSYMBT,1)

C       ZERO THE EXTRA POSITIONS
        DO I = 25,49
           CALL CCPMVI(HEADBUF(I), 0,1)
        ENDDO

C	ORIGIN ON X,Y & Z AXIS 
        ORIGX = (NX/2) + (MOD(NX,2))
        ORIGY = (NY/2) + (MOD(NY,2))
        ORIGZ = (NZ/2) + (MOD(NZ,2))
 
 	WRITE(NOUT,*)' '
 	WRITE(NOUT,*)' ORIGIN (X,Y,Z) DEFAULT VALUE:((NX/2)+ 1, ',
     &               '(NY/2)+1), (NZ/2)+1'
 	WRITE(NOUT,*)'   (+1 ADDED ONLY IF AXIS LENGTH IS ODD)'
	WRITE(NOUT,*)' USE <CR> FOR DEFAULT VALUE'
	
	CALL RDPRM3S(ORIGX,ORIGY,ORIGZ,NOT_USED,
     &               'X, Y, & Z ORIGIN FOR MRC DATA',IRTFLG)
	
        HEADBUF(50) = ORIGX
        HEADBUF(51) = ORIGY
        HEADBUF(52) = ORIGZ

C       PUT IN 'MAP'
        CVAL = 'MAP '
        CALL CCPMVC(HEADBUF(53),CVAL,ISSWAB)

C       SET MACHINE STAMP
        CALL SETSTAMP(MACHST,ISSWAB)
        CALL CCPMVI(HEADBUF(54),MACHST,1)

c        write(nout,*) ' map:',cval,headbuf(53)
c        write(nout,*) ' machst:',machst,headbuf(54)

C       SET RMS (ASSUMING THIS SAME AS SIG IN SPIDER???)
        HEADBUF(55) = DSIG

C       SET NUMBER OF LABELS
        NLABEL = 1  
        CALL CCPMVI(HEADBUF(56), NLABEL,1)

C       ZERO ALL LABELS WITH BLANKS
        CVAL = BLANK // BLANK // BLANK // BLANK
        DO I = 57,256
           CALL CCPMVC(HEADBUF(I),CVAL,ISSWAB)
        ENDDO

C       NOW GO BACK AND ADD IN ONE LABEL
        LABEL1 = 'Converted from SPIDER to MRC using SPIDER '
        INOW   = 56
        DO I = 1,44,4
           CVAL = LABEL1(I:I+3)
           INOW = INOW + 1
           CALL CCPMVC(HEADBUF(INOW),CVAL,ISSWAB)
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

C ---------------------- SETSTAMP -----------------------------------
 
      SUBROUTINE SETSTAMP(MACHSTMP,ISSWAB)

C     PURPOSE: SETS MACHINE STAMP FOR THIS ARCHITECTURE

C     NOTE: I HAVE EXTRACTED THIS FROM THE MRC 2000 CODE AND
C           CONVERTED TO FORTRAN. BUT I MAY HAVE BOTCHED IT? al

      INTEGER * 4 :: MACHSTMP
      LOGICAL     :: ISSWAB

      INCLUDE 'CMBLOCK.INC'

c       Little-ended Intel and some MIPS
#if defined(MIPSEL) || defined(i386) || defined(i860)
#  define NATIVEIT 4
#  define NATIVEFT 4
#endif

C       AN attempt at machines using the powerPC chip.             
#if defined (SP_PPC)
#  define NATIVEIT 1
#  define NATIVEFT 1
#endif

C       Compaq alpha's running Unix.             
#ifdef __alpha
#  define NATIVEFT 4
#  define NATIVEIT 4
#endif

C        SGI Altix running GNU/Linux.             
#ifdef __ia64
#  define NATIVEFT 4
#  define NATIVEIT 4
#endif

C       Big-endian ieee includes SGI machines,       
C       HP (68k-based or RISC), RS/6000 and all        
C       Suns except obsolete i386-based ones.       

#if defined(MIPSEB) || defined(__hpux) || defined(_AIX) || defined(m68k) || defined(mc68000) || defined(sparc)
#  define NATIVEIT 1
#  define NATIVEFT 1
#endif

C       Little-endian AMD OPTERON,       
#if defined  __x86_64__ || defined __ORDER_LITTLE_ENDIAN__
#  define NATIVEIT 4
#  define NATIVEFT 4
#endif


#if defined(SP_IBMSP3) 
#  define NATIVEIT 1
#  define NATIVEFT 1
#endif


#ifndef NATIVEFT
#error "Can't determine machine number format"
#endif

      NATIVEFTT = NATIVEFT
      NATIVEITT = NATIVEIT

      IF (ISSWAB) THEN
C        RUNNING WITH NON-NATIVE BYTE-SWAPPING
         IF (NATIVEFTT == 1) THEN
            MACHSTMP = 4369
         ELSE
            MACHSTMP = 286326784
         ENDIF
      ELSE
         IF (NATIVEFTT == 1) THEN
            MACHSTMP = 286326784
         ELSE
            MACHSTMP = 4369
         ENDIF
       ENDIF
c        write(nout,*) ' MACHSTMP:',MACHSTMP


      END


C -----------------------------------------------------------

      SUBROUTINE CCPMVIT(ARR1,ARR2,N)

C     THIS ROUTINE ASSIGNS  ARR2 TO ARR1

      REAL  :: ARR1(N),ARR2(N)

      DO I = 1,N
         ARR1(I) = ARR2(I)
      ENDDO

      END

 
