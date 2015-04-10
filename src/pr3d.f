
C++*********************************************************************
C
C PR3D.F           ADDED FOURIER INPUT             JUL 00 ARDEAN LEITH
C                  OPFILEC                         FEB 03 ARDEAN LEITH
C                  OUTPUT FORMAT, ERROR TRAP       NOV 09 ARDEAN LEITH
C                  OUTPUT FORMAT, SIZCHK           FEB 11 ARDEAN LEITH
C                  FSC                             FEB 12 ARDEAN LEITH
C                  MASK                            SEP 12 ARDEAN LEITH
C                  FSCCUT                          SEP 12 ARDEAN LEITH
C                  WANTSQRTS                       MAY 14 ARDEAN LEITH
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
C PURPOSE:  CALCULATE 3-D PHASE RESIDUE OUTSIDE MISSING CONE, 
C           PR OF FOURIER RINGS(RADIUS, DIRECTION RELATIVE TO Z) AND 
C           OF FOURIER SHELLS (RADIUS). 
C
C OPERATIONS:  'RF 3' and 'FSC'
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE PR3D(FSCOP)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        LOGICAL               :: FSCOP

        CHARACTER(LEN=MAXNAM) :: FILNAM1,FILNAM2 

        REAL, ALLOCATABLE     :: AIMG(:,:,:),BIMG(:,:,:)
        REAL, ALLOCATABLE     :: CSUM1(:), CSUM(:)
        REAL, ALLOCATABLE     :: PR(:,:), AMP(:,:), AVSUM(:,:)
        REAL, ALLOCATABLE     :: CSUM2(:)
        INTEGER, ALLOCATABLE  :: LR(:)

        REAL                  :: FSCCUT
        REAL                  :: SSANG, WI, PIXSIZ,FMAXSPFREQ
        REAL                  :: RADMASK
        LOGICAL               :: WANTSQRTS 

        CHARACTER(LEN=1)      :: SER
        CHARACTER(LEN=1)      :: NULL = CHAR(0)

        INTEGER               :: ICOMM,MYPID,MPIERR

        INTEGER               :: MAXIM
        INTEGER, PARAMETER    :: NSCALE = 20
        INTEGER, PARAMETER    :: LUN1   = 21
        INTEGER, PARAMETER    :: LUN2   = 22
        INTEGER, PARAMETER    :: LUNGP  = 23
        INTEGER, PARAMETER    :: LUNDOC = 89

        CALL SET_MPI(ICOMM,MYPID,MPIERR)

        WANTSQRTS = (FSCOP == .TRUE.)

C       INPUT FIRST IMAGE
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM1,LUN1,'O',ITYPE1,
     &               NX1,NY1,NZ1,MAXIM,
     &               'FIRST INPUT VOLUME~',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (NZ1 < 2) THEN
C          ONLY FOR VOLUMES
           CALL ERRT(101,'OPERATION ONLY FOR VOLUMES',NE)
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
     &               'SECOND INPUT VOLUME~',.TRUE.,IRTFLG)

        IF (IRTFLG .NE. 0) THEN
           GOTO 9998
        ELSEIF (ITYPE2 < 0) THEN
C          SECOND FILE IS A FOURIER INPUT FILE
           LSD2 = NX2
        ELSE
           LSD2 = NX2 + 2 - MOD(NX2,2)
        ENDIF

        CALL SIZCHK(UNUSED,LSD1,NY1,NZ1,0,
     &                     LSD2,NY2,NZ2,0,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        LSD        = LSD1
        NY         = NY1  
        NZ         = NZ1
        PIXSIZ     = 1 
        FMAXSPFREQ = 0.5 / PIXSIZ
        WI         = 0.5
        RADMASK    = 0
        FSCCUT     = 0.5 

        IF (FSCOP) THEN 
           SCALE1 = 0.2     ! SCALE FACTOR1
           SCALE2 = 2.0     ! SCALE FACTOR2
           SER    = 'C'     ! CONE
           SSANG  = 90.0    ! MAX. TILT ANGLE
           FACT   = 3       ! FACTOR FOR NOISE COMPARISON

           IF (ITYPE1 > 0 .AND. ITYPE2 > 0) THEN
              CALL RDPRM2S(WI, RADMASK, NOT_USED, 
     &          'SHELL WIDTH (RECIPROCAL SPACE UNITS), MASKING RADIUS',
     &          IRTFLG)
           ELSE
              CALL RDPRM1S(WI, NOT_USED, 
     &          'SHELL WIDTH (RECIPROCAL SPACE UNITS)',IRTFLG)
           ENDIF
           IF (IRTFLG .NE. 0) GOTO 9999

           PIXSIZ = 1 
           CALL RDPRM2S(PIXSIZ, FSCCUT, NOT_USED, 
     &                 'VOXEL SIZE (A) & RESOLUTION CUTOFF',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           FMAXSPFREQ = 0.5 / PIXSIZ

        ELSE 
           CALL RDPRM1S(WI, NOT_USED, 
     &        'SHELL WIDTH (RECIPROCAL SPACE UNITS)',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           CALL RDPRM2S(SCALE1,SCALE2,NOT_USED,
     &                'SCALE FACTOR (LOWER,UPPER)',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           CALL RDPRMC(SER,NUMC,.TRUE.,
     &                 'MISSING CONE/WEDGE ANGLE (C/W)',
     &                 NULL,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           SSANG  = 90.0    ! MAX. TILT ANGLE

           CALL RDPRM1S(SSANG,NOT_USED,
     &       'MAXIMUM TILT ANGLE (90 FOR MOST SINGLE PARTICLE RECONS.)',
     &       IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           IF (SER .NE. 'C' .AND. SER .NE. 'W') SSANG = 90.0

           FACT   = 3       ! DEFAULT FACTOR FOR NOISE COMPARISON
           CALL RDPRM1S(FACT,NOT_USED,
     &               'NOISE COMPARISON FACTOR',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
        ENDIF

        DSCALE = (SCALE2-SCALE1) / FLOAT(NSCALE-1)


        Y1     = FLOAT(MAX(NX,NY,NZ))
        INC    = INT(Y1/WI) / 2+1

        ALLOCATE (AIMG(LSD,NY,NZ),
     &            BIMG(LSD,NY,NZ), 
     &            STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MWANT = 2*LSD*NY*NZ  
           CALL ERRT(46,'PR3D; AIMG & BIMG',MWANT)
           GOTO 9999
        ENDIF

        CALL READV(LUN1,AIMG, LSD,NY,NX1,NY,NZ)
        CALL READV(LUN2,BIMG, LSD,NY,NX2,NY,NZ)

        IF (RADMASK > 0) THEN
C          SUPERGAUSSIAN MASKING WANTED

           XCEN       = (NX/2) + 1
           YCEN       = (NY/2) + 1
           ZCEN       = (NZ/2) + 1
	   TNM        = ALOG(1.0 / TINY(TNM))
           RADMASKSQI = 1.0 / (RADMASK**2)

           DO K = 1,NZ
              DO J = 1,NY
                 DO I = 1,NX
	           EEE = 0.5 * ((I - XCEN) **2 * RADMASKSQI +
     &                          (J - YCEN) **2 * RADMASKSQI +
     &                          (K - ZCEN) **2 * RADMASKSQI)

	           IF (EEE  >= TNM) THEN
	              AIMG(I,J,K) = 0.0
	              BIMG(I,J,K) = 0.0
	           ELSE  
	              EEE         = -0.5 * (2*EEE)**2
                      AIMG(I,J,K) = EXP(EEE) * AIMG(I,J,K)
                      BIMG(I,J,K) = EXP(EEE) * BIMG(I,J,K)
	           ENDIF
                ENDDO
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

        ALLOCATE(PR(NSCALE,INC), AMP(NSCALE,INC), CSUM1(INC),
     &           LR(INC),CSUM(INC),CSUM2(INC), AVSUM(NSCALE,INC),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           MWANT = 3*NSCALE*INC + 4*INC
           CALL ERRT(46,'PR3D; PR...',MWANT)
           GOTO 9999
        ENDIF

C       CALCULATIONS
        CALL PR3DB(AIMG,BIMG,PR,AMP,CSUM1,LR,CSUM,CSUM2,
     &      AVSUM,LSD,NX,NY,NZ,DSCALE,NSCALE,SCALE1,
     &      SSANG,INC,Y1,WI,SER)

#ifdef NEVER
c*************
        INV = -1
        CALL FMRS_3(AIMG,NX,NY,NZ,INV)
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM1,LUN1,'U',ITYPE1,NX,NY,
     &          NZ,MAXIM,'OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CALL WRITEV(LUN1,AIMG,LSD,NY,NX,NY,NZ)
	CLOSE(LUN1)
c*********************
#endif

C       WRITE RESULT INTO DOC FILE AND RESULT FILE
        CALL  RFACTSD2(PR,AMP,CSUM1,LR,CSUM,CSUM2,AVSUM,
     &                 NSCALE,INC,WI,FACT,.FALSE.,
     &                 LUNDOC,FSCOP,FMAXSPFREQ,LUNGP,FSCCUT,
     &                 WANTSQRTS)
       
9999	IF (ALLOCATED(AIMG))  DEALLOCATE (AIMG)
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
