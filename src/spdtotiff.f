
C++*********************************************************************
C
C  SPDTOTIFF.FOR  --    CREATED, DEC 28, 94 al
C                       OSF REMOVED                JUL 2009 ArDean Leith
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
C
C    SPDTOTIFF(LUNO,LUNN,NX,NY,NZ,IRTFLG)
C
C    PURPOSE:       CONVERT SPIDER IMAGE TO 8BIT TIFF FORMAT
C
C    PARAMETERS:
C        LUNO       LOGICAL UNIT NUMBER TO BE ASSIGNED INPUT.
C        LUNN       LOGICAL UNIT NUMBER TO BE ASSIGNED TO FILNEW
C        IRTFLG     ERROR RETURN FLAG. (0 IS NORMAL)
C
C    VARIABLES:  TIFF TAGS  
C	TIFFTAG_IMAGEWIDTH		256
C	TIFFTAG_IMAGELENGTH		257
C	TIFFTAG_BITSPERSAMPLE		258
C	TIFFTAG_COMPRESSION		259
C	TIFFTAG_PHOTOMETRIC		262
C	TIFFTAG_STRIPOFFSETS		273
C	TIFFTAG_ROWSPERSTRIP		278
C	TIFFTAG_STRIPBYTECOUNTS		279
C	TIFFTAG_RESOLUTIONUNIT		296
C
C    FUTURE      FIND IF CURRENTLY SWAPPING BYTES
C                FLIP = ISSWAB(99)  
C
C--********************************************************************

	SUBROUTINE SPDTOTIFF(LUNO,LUNN,NX,NY,NZ,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        INTEGER                :: IERR
        COMMON /IOERR/ IERR

 	INTEGER                :: LUNO,LUNN,NX,NY,NZ,IRTFLG

        REAL                   :: BUFO(NBUFSIZ) ! FROM CMLIMIT

        INTEGER, PARAMETER     :: LENOPN = 1024 ! ONLY USES: 512
        INTEGER*1              :: LBUF(LENOPN)

	INTEGER*2              :: TIFFINT(70), TIFFDAT(70)
	INTEGER*1              :: TIFFBYTE(140)
	EQUIVALENCE (TIFFINT, TIFFBYTE)

        CHARACTER (LEN=1)      :: NULL = CHAR(0)
        CHARACTER (LEN=MAXNAM) :: FILNEW

        LOGICAL                :: FLIP
	REAL	               :: PVAL
        INTEGER                :: IRTFLGT, ISLICE, NOT_USED
        INTEGER                :: ISTREC, IENDREC, K, IRECOUT, I, NE 
        INTEGER                :: J, NLET 

        INTEGER                :: lnblnkn 

        character (len=2)      :: cvals = 'ii'
        integer * 2            :: i2val
	equivalence (cvals, i2val)
   

C	SET THE TIFFINT INVARIANT HEADER
C       BYTE ORDER, VERSION(=42), OFFSET TO IFD, IFD COUNT,
C       DIRECTORY ENTRIES: TAG #, TYPE(3=SHORT,4=LONG),
C                          LENGTH(4 BYTES), VALUE(4 BYTES)

        !  'II' is: 18761

	DATA TIFFDAT/19789,42,0,8,9, 
     &		       256, 4,0,1,0,0,	257,4,0,1,0,0,
     &		       258, 3,0,1,8,0,	259,3,0,1,1,0,
     &		       262, 3,0,1,1,0,	273,4,0,1,0,1024,
     &		       278, 4,0,1,0,0,	279,4,0,1,0,0,
     &		       296, 3,0,1,3,0,	  0,0,0,0,0,0,
     &                   0, 0,0,0,0/

        IRTFLG = 1

        !WRITE(6,*) ' cvals,i2vals:',cvals,i2val,tiffdat(1)

C	INITIALIZE THE TIFF HEADER WITH CONSTANT VALUES
        TIFFINT = TIFFDAT

        !WRITE(6,*) ' tiffint(1):',tiffint(1)
        !TIFFINT(1:2) = TRANSFER(CVALS(1:2),I2VAL)
        !WRITE(6,*) ' tiffint(1):',tiffint(1)

	    
C       CONVERT SPIDER HEADER TO TIFF HEADER           

C	TIFFTAG_IMAGEWIDTH		256
        TIFFINT(11) = NX

C	TIFFTAG_IMAGELENGTH		257
	TIFFINT(17) = NY

C	TIFFTAG_STRIPOFFSETS		273  ! ONLY ONE STRIP
	TIFFINT(41) = 512

C	TIFFTAG_ROWSPERSTRIP		278
	TIFFINT(47) = NY
	
C	TIFFTAG_STRIPBYTECOUNTS		279  ! ONLY ONE STRIP
	TIFFINT(52) = NX*NY/65536.0
	TIFFINT(53) = NX*NY
	
C       IF NECESSARY FIND STATS FOR INPUT IMAGE 
        IF (IMAMI == 0) CALL NORM3(LUNO,NX,NY,NZ,FMAX,FMIN,AV)

C	COMPUTE THE SCALE VALUE
	PVAL = 255.0 / (FMAX - FMIN)

C       GET OUTPUT FILE NAME
C       OPEN NEW FILE FOR TIFF VERSION, 1024 BYTE RECORD LENGTH
        CALL OPAUXFILE(.TRUE.,FILNEW,'tif',LUNN,LENOPN,'N',
     &                 'TIFF OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C	CHECK THE DATA IS AN IMAGE OR A VOLUME
C	IF IT IS AN IMAGE, CONVERT IT TO A TIFF FILE
C	IF IT IS A VOLUME, ASK THE USER WHICH SLICE SHOULD BE CONVERTED
	IF (ABS(NZ) <= 1) THEN
C	    IT IS AN IMAGE ....
	    ISLICE = 1

	ELSEIF (ABS(NZ) > 1) THEN
C	    IT IS A VOLUME, ASK FOR THE SLICE NO. THAT WILL BE CONVERTED
	    CALL RDPRI1S(ISLICE,NOT_USED, 'SLICE NUMBER',IRTFLG)
            IF (IRTFLG .NE. 0) RETURN
	ENDIF

	ISTREC   = (ISLICE - 1) * NY + 1
	IENDREC  = ISTREC + NY - 1

C	WRITE THE TIFF HEADER 
        LBUF(1:512) = 0

C       FLIP THE HEADER'S BYTES ON INTEL LINUX 
        DO K= 1, 140,2
           LBUF(K+1) = TIFFBYTE(K)
           LBUF(K)   = TIFFBYTE(K+1)
        ENDDO



C       LBUF HAS 512 BYTES IN IT NOW FROM TIFF HEADER
        K       = 512
        IRECOUT = 0
        IERR    = 0
	    
        DO I = ISTREC,IENDREC

C         READ EACH RECORD OF SPIDER FILE
          CALL REDLIN(LUNO,BUFO,NX,I)
          IF (IERR .NE. 0) THEN
             CALL ERRT(102,'READING SPIDER FILE',IERR)
             IRTFLG = 1
             RETURN
          ENDIF

C         CONVERT FLOATING POINT NUMBERS TO BYTE
          DO J=1,NX
             IF (K >= LENOPN) THEN
C               BUFFER IS FULL, PUT IT OUT TO FILE (8bits = no byteswap)
                IRECOUT = IRECOUT + 1
                CALL WRTLIN8(LUNN,LBUF,LENOPN,IRECOUT)
                K = 0
             ENDIF
             K       = K + 1
             LBUF(K) = (BUFO(J) - FMIN) * PVAL
           ENDDO
        ENDDO

        IF (K > 0) THEN
C          BUFFER STILL HAS PIXELS IN IT, PUT THEM OUT TO FILE
           IRECOUT = IRECOUT + 1
           CALL WRTLIN8(LUNN,LBUF,K,IRECOUT)  ! DOES NOT byteswap 
        ENDIF

        NLET = LNBLNKN(FILNEW)
        WRITE(NOUT,*) ' OUTPUT PLACED IN: ',FILNEW(1:NLET)

C       SET FLAG FOR NORMAL RETURN	
        IRTFLG = 0

        END


        INTEGER FUNCTION NEWUNIT(UNIT) RESULT(N)
        ! RETURNS LOWEST I/O UNIT NUMBER NOT IN USE

        INTEGER, INTENT(OUT), OPTIONAL :: UNIT
        LOGICAL                        :: INUSE
        INTEGER, PARAMETER             :: NMIN=200   ! AVOID LOWER NUMBERS 
        INTEGER, PARAMETER :: NMAX=999  ! MAY BE SYSTEM-DEPENDENT

        DO N = NMIN, NMAX
           INQUIRE(UNIT=N, OPENED=INUSE)
           IF (.NOT. INUSE) THEN
              IF (PRESENT(UNIT)) UNIT = N
              RETURN
           ENDIF
        ENDDO

        CALL ERRT(101,'NO newunit AVAILABLE',NE)

        END FUNCTION
