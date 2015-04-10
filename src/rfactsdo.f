C++*********************************************************************
C
C    RFACTSDO.F    ADDED FOURIER INPUT          JUL  2000 ARDEAN LEITH
C                  OPFILEC                      FEB  2003 ARDEAN LEITH
C                  FSCOP                        FEB  2012 ARDEAN LEITH
C                  MASK                         SEP  2012 ARDEAN LEITH
C                  FSCCUT                       SEP  2012 ARDEAN LEITH
C                  WANTSQRTS                    MAY  2014 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C PURPOSE: DIFFERENTIAL CORRELATION COEFFICIENT  AND PHASE RESIDUAL
C          COMPARISON OF 
C          TWO FOURIER TRANSFORMS ON A SERIES OF RINGS.  FOLLOWS
C          PHILOSOPHY OF FRANK ET AL. SCIENCE 214 (1981) 1353-1355. 
C          IN THE CURRENT ROUTINE, ALL RINGS ARE COMPUTED AT ONCE,
C          BUT SCALE SEARCH IS DONE ON EACH RING SEPARATELY.  THUS THE
C          RESULT IS EQUIVALENT TO THE RESULT OF APPLYING A SERIES
C          OF CALLS TO "RF S" WITH SUCCESSIVE RINGS.
C          NOTE THAT THIS APPROACH WILL LEAD TO UNREASONABLE PHASE 
C          RESIDUAL RESULTS IF THE TWO FOURIER TRANSFORMS HAVE 
C          STRONGLY DIFFERENT RADIAL BEHAVIORS.
C
C OPERATIONS:  'RF' and 'FRC'
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C-************************************************************************

        SUBROUTINE RFACTSDO(FSCOP)

        INCLUDE 'CMLIMIT.INC'      
        INCLUDE 'CMBLOCK.INC'      
 
        LOGICAL               :: FSCOP

        INTEGER               :: NX,NY,NZ,NX2,NY2,NZ2
        INTEGER               :: ITYPE1,ITYPE2,MAXIM
        REAL                  :: PIXSIZ,RADMASK,FSCCUT
        LOGICAL               :: WANTSQRTS = .FALSE.

        REAL, ALLOCATABLE     :: AIMG(:,:),BIMG(:,:),CSUM1(:),CSUM(:)
        REAL, ALLOCATABLE     :: PR(:,:),AMP(:,:),AVSUM(:,:)
        REAL, ALLOCATABLE     :: CSUM2(:)
        INTEGER, ALLOCATABLE  :: LR(:)

        CHARACTER (LEN=1)     :: SER
        CHARACTER(LEN=MAXNAM) :: FILNAM1,FILNAM2

        INTEGER, PARAMETER    :: NSCALE = 20
        INTEGER, PARAMETER    :: LUN1   = 21
        INTEGER, PARAMETER    :: LUN2   = 22
        INTEGER, PARAMETER    :: LUNGP  = 23
        INTEGER, PARAMETER    :: LUNDOC = 89

        CALL SET_MPI(ICOMM,MYPID,MPIERR)

C       INPUT FIRST IMAGE
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM1,LUN1,'O',ITYPE1,
     &               NX1,NY1,NZ1,MAXIM,
     &               'FIRST INPUT IMAGE~',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (NZ1 > 1) THEN
           CALL ERRT(101,'OPERATION ONLY FOR IMAGES',NE)
           GOTO 9998
        ELSEIF (ITYPE1 < 0) THEN
C          FOURIER INPUT FILE
           LSD1 = NX1
           NX   = NX1 - MOD(-ITYPE1,10)
        ELSE
C          REAL INPUT FILE
           LSD1 = NX1 + 2 - MOD(NX1,2)
           NX   = NX1
        ENDIF

C       INPUT SECOND IMAGE
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM2,LUN2,'O',ITYPE2,
     &               NX2,NY2,NZ2,MAXIM,
     &               'SECOND INPUT IMAGE~',.TRUE.,IRTFLG)

        IF (IRTFLG .NE. 0) THEN
           GOTO 9998
        ELSEIF (ITYPE2 < 0 .AND. ITYPE1 > 0) THEN
C          FIRST FILE IS NOT FOURIER INPUT FILE
           CALL ERRT(101,'BOTH FILES MUST BE FOURIER FORMAT',NE)
           GOTO 9998
        ELSEIF (ITYPE2 < 0) THEN
C          SECOND FILE IS A FOURIER INPUT FILE
           LSD2  = NX2
        ELSE
C          SECOND FILE IS A REAL INPUT FILE
           LSD2  = NX2 + 2 - MOD(NX2,2)
        ENDIF

        CALL SIZCHK(UNUSED,LSD1,NY1,NZ1,0,
     &                     LSD2,NY2,NZ2,0,IRTFLG)

        LSD        = LSD1
        NY         = NY1 
        NZ         = NZ1 

        PIXSIZ     = 1
        FMAXSPFREQ = 0.5 / PIXSIZ
        RADMASK    = 0
        WI         = 0.5    ! RING WIDTH
        FSCCUT     = 0.5    ! DEFAULT RESOLUTION CUT LINE
 
        IF (FSCOP) THEN 
           SCALE1 = 0.2     ! SCALE FACTOR1
           SCALE2 = 2.0     ! SCALE FACTOR2

           IF (ITYPE1 > 0 .AND. ITYPE2 > 0) THEN
              CALL RDPRM2S(WI, RADMASK, NOT_USED, 
     &          'RING WIDTH (RECIPROCAL SPACE UNITS), MASKING RADIUS',
     &          IRTFLG)
           ELSE
              CALL RDPRM1S(WI, NOT_USED, 
     &          'RING WIDTH (RECIPROCAL SPACE UNITS)',IRTFLG)
           ENDIF
           IF (IRTFLG .NE. 0) GOTO 9999

           WIP    = 1.0 / WI
           CALL RDPRM2S(PIXSIZ, FSCCUT, NOT_USED, 
     &         'PIXEL SIZE (A) & RESOLUTION CUTOFF',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           FMAXSPFREQ = 0.5 / PIXSIZ
        
        ELSE 
           CALL RDPRM1S(WI, NOT_USED, 
     &        'RING WIDTH (RECIPROCAL SPACE UNITS)',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           CALL RDPRM2S(SCALE1,SCALE2,NOT_USED,
     &                'SCALE FACTORS (LOWER,UPPER)',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

        ENDIF

C       MISSING ANGLE STUFF IS ZEROED
        SER    = 'C'
        SSANG  = 90.0

C       FACTOR FOR NOISE COMPARISON
        FACT   = 3.0
        DSCALE = (SCALE2-SCALE1) / FLOAT(NSCALE-1)
        Y1     = FLOAT(MAX(NX,NY))
        INC    = INT(Y1/WI) / 2+1

        ALLOCATE (AIMG(LSD,NY), 
     &            BIMG(LSD,NY), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           MWANT = 2*LSD*NY
           CALL ERRT(46,'RFACTSDO; AIMG & BIMG',MWANT)
           GOTO 9999
        ENDIF

        IF (ITYPE1 > 0) THEN
C          FIRST INPUT FILE IS REAL SPACE, NOT FOURIER
           CALL READV(LUN1,AIMG,LSD,NY,NX1,NY,NZ)
        ELSE
C          FIRST INPUT IS FOURIER ALREADY
           CALL READV(LUN1,AIMG, LSD,NY,NX1,NY,NZ)
        ENDIF
	
        IF (ITYPE2 > 0) THEN
C          SECOND INPUT FILE IS REAL SPACE, NOT FOURIER
           CALL READV(LUN2,BIMG, LSD,NY,NX2,NY,NZ)
        ELSE
C          SECOND INPUT IS FOURIER ALREADY
           CALL READV(LUN2,BIMG, LSD,NY,NX2,NY,NZ)
        ENDIF

        IF (RADMASK > 0) THEN
C          SUPERGAUSSIAN MASKING WANTED

           XCEN       = (NX/2) + 1
           YCEN       = (NY/2) + 1
 	   TNM        = ALOG(1.0 / TINY(TNM))
           RADMASKSQI = 1.0 / (RADMASK**2)

           DO J = 1,NY
              DO I = 1,NX
	        EEE = 0.5 * ((I - XCEN) **2 * RADMASKSQI +
     &                       (J - YCEN) **2 * RADMASKSQI)

	        IF (EEE  >= TNM) THEN
	           AIMG(I,J) = 0.0
	           BIMG(I,J) = 0.0
	        ELSE  
	           EEE       = 0.5 * (2*EEE)**2
                   AIMG(I,J) = EXP(-EEE) * AIMG(I,J)
                   BIMG(I,J) = EXP(-EEE) * BIMG(I,J)
	        ENDIF
             ENDDO
          ENDDO
        ENDIF

        IF (ITYPE1 > 0) THEN
C          FIRST INPUT FILE IS REAL SPACE, NOT FOURIER
           INV = 1
           CALL FMRS_3(AIMG,NX1,NY,NZ,INV)
           IF (INV == 0) THEN
              CALL ERRT(101,'FFT ERROR',NE)
              GOTO 9999
           ENDIF
        ENDIF
	
        IF (ITYPE2 > 0) THEN
C          SECOND INPUT FILE IS REAL SPACE, NOT FOURIER
           INV = 1
           CALL FMRS_3(BIMG,NX2,NY,NZ,INV)
           IF (INV == 0)THEN
              CALL ERRT(101,'FFT ERROR',NE)
              GOTO 9999
           ENDIF 
        ENDIF

        ALLOCATE(PR(NSCALE,INC), AMP(NSCALE,INC),CSUM1(INC),
     &           LR(INC),CSUM(INC),CSUM2(INC), AVSUM(NSCALE,INC),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           MWANT = 3*NSCALE * INC + 4*INC
           CALL ERRT(46,'RFACTSDO; PR... ',MWANT)
           GOTO 9999
        ENDIF

C       CALCULATIONS
        CALL PR3DB(AIMG,BIMG,PR,AMP,CSUM1,LR,CSUM,CSUM2,
     &      AVSUM,LSD,NX,NY,NZ,DSCALE,NSCALE,SCALE1,
     &      SSANG,INC,Y1,WI,SER)

C       WRITE RESULT INTO DOC FILE AND RESULT FILE
        CALL RFACTSD2(PR,AMP,CSUM1,LR,CSUM,CSUM2,AVSUM,
     &                NSCALE,INC,WI,FACT,.TRUE.,
     &                LUNDOC,FSCOP,FMAXSPFREQ,LUNGP,FSCCUT,
     &                WANTSQRTS)

        IF (MYPID <= 0) WRITE(NOUT,*)' '

9999    IF (ALLOCATED(AIMG))  DEALLOCATE (AIMG)
        IF (ALLOCATED(BIMG))  DEALLOCATE (BIMG)
        IF (ALLOCATED(PR))    DEALLOCATE (PR)
        IF (ALLOCATED(AMP))   DEALLOCATE (AMP)
        IF (ALLOCATED(CSUM1)) DEALLOCATE (CSUM1)
        IF (ALLOCATED(LR))    DEALLOCATE (LR)
        IF (ALLOCATED(CSUM))  DEALLOCATE (CSUM)
        IF (ALLOCATED(CSUM2)) DEALLOCATE (CSUM2)
        IF (ALLOCATED(AVSUM)) DEALLOCATE (AVSUM)

9998    CLOSE(LUN1)
        CLOSE(LUN2)
        CLOSE(LUNDOC)
        CLOSE(LUNGP)

        END
