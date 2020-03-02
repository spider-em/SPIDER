C ++********************************************************************
C                                                                      *
C FILTER                  CREATED             MAR 2001 ARDEAN LEITH    * 
C                         ADDED OPERATIONS    MAY 2001 ARDEAN LEITH    *
C                         3D FIXED           JUNE 2001 ARDEAN LEITH    *
C                         TOP HAT NEIGH FIXED APR 2002 ArDean Leith    *
C                         PRE & SOBEL 3D BUG  APR 2002 ArDean Leith    *
C                         LAHE                APR 2002 ArDean Leith    *
C                         VAR BUG FIXED       OCT 2002 ArDean Leith
C                         SOME STACKS SUPPORT OCT 2002 ArDean Leith    *
C                         VS PIXEL REGISTER   OCT 2004 ArDean Leith    *
C                                                                      *
C **********************************************************************
C=*                                                                    *
C=* Author: ArDean Leith                                               *                                                            *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
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
C  FILTER(LUN1,LUN2,NSAM,NROW,NSLICE,MODE)
C
C  PARAMETERS: LUN1,LUN2   IO UNITS                             (INPUT)
C              NSAM        X DIMENSIONS                         (INPUT)
C              NROW        Y DIMENSIONS                         (INPUT)
C              NSLICE      Z DIMENSIONS                         (INPUT)
C              MODE        TYPE SPECIFIER                       (INPUT)
C              FMINT       INPUT IMAGE MINIMUM                  (INPUT)
C                          F:  FREI-CHEN - 3x3 NEIGHBORHOOD
C                          G : GRADIENT  - 3x3 NEIGHBORHOOD
C                          LA: LAPLACIAN - RECTILINEAR NEIGHBORHOOD
C                          LH: LOCAL AREA HISTOGRAM EQUALIZATION
C                          MA: MAXIMUM   - RECTILINEAR NEIGHBORHOOD
C                          MI: MINIMUM   - RECTILINEAR NEIGHBORHOOD
C                          p:  PREWITT   - 3x3 NEIGHBORHOOD
C                          R:  RANGE     - RECTILINEAR NEIGHBORHOOD
C                          RI: RIDGE FOLLOWING
C                          S:  SOBEL     - 3x3 NEIGHBORHOOD
C                          T:  TOPHAT    - CIRCULAR NEIGHBORHOOD
C                          V:  VARIANCE  - RECTILINEAR NEIGHBORHOOD
C                          VS: VARIANCE SMOOTHING  - RECTILINEAR NEIGH.
C
C  PURPOSE: ALTER CONTRAST IN AN IMAGE OR VOLUME USING CONVOLUTION OR 
C           OTHER LOCAL IMAGE SPACE METHODS
C                                                                       
C **********************************************************************

	SUBROUTINE FILTER(LUN1,LUN2,NSAM,NROW,NSLICE,MAXIM,
     &                    MODE,FMINT,FMAXT,SIGT)

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

	REAL, ALLOCATABLE, DIMENSION(:)  ::  VIN,VKERNAL

	REAL, DIMENSION(9)  ::  GRAD2,GRAD2X,GRAD2Y,LAPL2
	REAL, DIMENSION(9)  ::  SOBX,SOBY,PREX,PREY
	REAL, DIMENSION(27) ::  GRAD3,GRAD3X,GRAD3Y,GRAD3Z,LAPL3
	REAL, DIMENSION(27) ::  SOBX3,SOBY3,SOBZ3,PREX3,PREY3,PREZ3

	CHARACTER(LEN=2) ::   MODE 
        LOGICAL          ::   THREED,LOADIT

        REAL, PARAMETER :: sq2 = 0.7071
        REAL, PARAMETER :: sq3 = 0.57735


C       SIMPLE 3X3 LAPLACIAN (DYNAMIC BELOW)
        DATA LAPL2 / -1.0,-1.0,-1.0,
     &               -1.0, 8.0,-1.0,  
     &               -1.0,-1.0,-1.0/

        DATA LAPL3 / -1.0,-1.0,-1.0,  -1.0,-1.0,-1.0,  -1.0,-1.0,-1.0,
     &               -1.0,-1.0,-1.0,  -1.0,26.0,-1.0,  -1.0,-1.0,-1.0,
     &               -1.0,-1.0,-1.0,  -1.0,-1.0,-1.0,  -1.0,-1.0,-1.0/

C       GRADIENT MAGNITUDE
        DATA GRAD2 /0.7071,  1.0, 0.7071,   
     &              1.0000, -1.0, 1.000,  
     &              0.7071,  1.0, 0.7071/

C       GRADIENT X COMPONENT
        DATA GRAD2X/ -.7071, 0.0, 0.7071,  
     &              -1.000,  0.0, 1.000,  
     &               -.7071, 0.0, 0.7071/

C       GRADIENT Y COMPONENT
        DATA GRAD2Y/0.7071, 1.0, 0.7071,  
     &              0.000,  0.0, 0.000,  
     &              -.7071,-1.0, -.7071/

C       GRADIENT MAGNITUDE
        DATA GRAD3/
     &      .57735,.7071,.57735, .7071,1.0,.7071, .57735,.7071,.57735,
     &      .7071,  1.0,  .7071,   1.0, 1.0, 1.0,   .7071, 1.0,  .7071,
     &      .57735,.7071,.57735, .7071,1.0,.7071, .57735,.7071,.57735/

C       GRADIENT X COMPONENT
        DATA GRAD3X/
     &        -.57735,0.0,.57735, -.7071,0.0,.7071, -.57735,0.0,.57735,
     &        -.7071, 0.0,.7071,  -1.0,  0.0, 1.0,  -.7071, 0.0, .7071,
     &        -.57735,0.0,.57735, -.7071,0.0,.7071, -.57735,0.0,.57735/

C       GRADIENT Y COMPONENT
        DATA GRAD3Y/
     &        .57735,.7071,.57735, 0.0,0.0,0.0, -.57735,-.7071,-.57735,
     &        .7071, 1.0,  .7071,  0.0,0.0,0.0, -.7071,  -1.0, -.7071,
     &        .57735,.7071,.57735, 0.0,0.0,0.0, -.57735,-.7071,-.57735/
        
C       GRADIENT Z COMPONENT
        DATA GRAD3Z/
     &          -.57735,-.7071,-.57735,  -.7071,-1.0,-.7071,
     &          -.57735,-.7071,-.57735,
     &          0.0, 0.0, 0.0,   0.0, 0.0,0.0,    0.0, 0.0, 0.0,
     &           .57735, .7071, .57735,   .7071, 1.0,.7071,    
     &           .57735, .7071, .57735/

C       SOBEL X MASK
        DATA SOBX/
     &       -1, 0, 1,
     &       -2, 0, 2,
     &       -1, 0, 1/

C       SOBEL Y MASK
        DATA SOBY/
     &        1, 2, 1,
     &        0, 0, 0,
     &       -1,-2,-1/

C       SOBEL 3D X MASK
        DATA  SOBX3/
     &       -1, 0, 1,
     &       -2, 0, 2,
     &       -1, 0, 1,
     &       -1, 0, 1,
     &       -2, 0, 2,
     &       -1, 0, 1,
     &       -1, 0, 1,
     &       -2, 0, 2,
     &       -1, 0, 1/

C       SOBEL 3D Y MASK
        DATA SOBY3/
     &        1, 2, 1,
     &        0, 0, 0,
     &       -1,-2,-1,

     &        1, 2, 1,
     &        0, 0, 0,
     &       -1,-2,-1,

     &        1, 2, 1,
     &        0, 0, 0,
     &       -1,-2,-1/

C       SOBEL 3D Z MASK
         DATA SOBZ3/
     &        1, 1, 1, 
     &        2, 2, 2, 
     &        1, 1, 1,

     &        0, 0, 0,
     &        0, 0, 0,
     &        0, 0, 0,

     &       -1,-1,-1, 
     &       -2,-2,-2, 
     &       -1,-1,-1/

C       PREWITT X MASK
        DATA PREX/
     &        1, 0,-1,
     &        1, 0,-1,
     &        1, 0,-1/

C       PREWITT Y MASK
        DATA PREY/
     &       -1,-1,-1,
     &        0, 0, 0,
     &        1, 1, 1/

C       PREWITT 3D X MASK
        DATA PREX3/
     &        1, 0,-1,
     &        1, 0,-1,
     &        1, 0,-1,
     &        1, 0,-1,
     &        1, 0,-1,
     &        1, 0,-1,
     &        1, 0,-1,
     &        1, 0,-1,
     &        1, 0,-1/

C       PREWITT 3D Y MASK
        DATA PREY3/
     &       -1,-1,-1,
     &        0, 0, 0,
     &        1, 1, 1,
     &       -1,-1,-1,
     &        0, 0, 0,
     &        1, 1, 1,
     &       -1,-1,-1,
     &        0, 0, 0,
     &        1, 1, 1/
 
C       PREWITT 3D Z MASK
        DATA PREZ3/
     &        1, 1, 1, 
     &        1, 1, 1, 
     &        1, 1, 1,
     &        0, 0, 0,
     &        0, 0, 0,
     &        0, 0, 0,
     &       -1,-1,-1, 
     &       -1,-1,-1, 
     &       -1,-1,-1/

        THREED = (NSLICE  >  1)

C       SET DEFAULT NEIGHBORS
        LX = 3
        LY = 3
        LZ = 1
        IF (THREED) LZ = 3
        NEIGH = LX * LY * LZ

C       CATEGORIES ONLY USED BY HURST
        ICAT = 1

C       SOME FILTERS NEED >1 KERNAL
        NK  = 1
        IF (MODE(1:1)  ==  'S' .OR. MODE(1:1)  ==  'P' .OR.
     &      MODE(1:1)  ==  'H') NK = 2
        IF ((MODE(1:1)  ==  'S' .OR. MODE(1:1)  ==  'P' ) .AND. 
     &       NSLICE  >  1) NK = 3
        IF (MODE(1:1)  ==  'F') NK = 9

        IF (MODE(1:1)  ==  'T') THEN
C          TOP-HAT OVER SPECIFIED NEIGHBORHOOD

11         CALL RDPRIS(IDIN,IDOUT,NOT_USED,
     &                 'INNER & OUTER DIAMETERS',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           IF ((IDIN   <  1 .OR. MOD(IDIN,2)   ==  0) .OR.
     &         (IDOUT  <  1 .OR. MOD(IDOUT,2)  ==  0)) THEN
              CALL ERRT(101,'DIAMETERS MUST BE ODD  AND > 0',IDUM) 
              GOTO 11
           ENDIF

           RADIN  = IDIN  / 2.0 
           RADOUT = IDOUT / 2.0 
           LX     = IDOUT
           LY     = LX
           NEIGH  = LX * LY
           LZ     = 0
           IF (THREED) THEN
              LZ    = LX
              NEIGH = NEIGH * LZ
           ENDIF

        ELSEIF (MODE(1:1)  ==  'H') THEN
C          HURST OVER SPECIFIED DIAMETER

12         CALL RDPRI1S(IDOUT,NOT_USED,'DIAMETER',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           IF (IDOUT  <  1 .OR. MOD(IDOUT,2)  ==  0) THEN
              CALL ERRT(101,'DIAMETER MUST BE ODD  AND > 0',IDUM) 
              GOTO 12
           ENDIF
           LX      = IDOUT 
           LY      = IDOUT
           NEIGH   = LX * LY
           LZ      = 1
           IF (THREED) THEN
              WRITE(NOUT,*) ' THIS FILTER IS NOT IMPLEMENTED IN 3D!' 
              WRITE(NOUT,*) ' VOLUME WILL BE PROCESSED SLICE-BY-SLICE' 
              THREED = .FALSE.
           ENDIF

        ELSEIF (MODE(1:2)  ==  'RA' .OR. MODE(1:1)  ==  'M' .OR.
     &          MODE(1:1)  ==  'L'  .OR. MODE(1:1)  ==  'V' .OR.
     &          MODE(1:1)  ==  '-'  .OR. MODE(1:1)  ==  '-' .OR.
     &          MODE(1:2)  ==  'VS' .OR. MODE(1:2)  ==  'LH') THEN
C          RANGE, MAX, MIN, LAPLACIAN, VARIANCE, LOCAL HISTOGRAM
C          VARIANCE SMOOTHING, FILTER OVER SPECIFIED NEIGHBORHOOD

10         IF (THREED) THEN
              CALL RDPRI3S(LX,LY,LZ,NOT_USED,
     &                    'NEIGHBORHOOD X, Y, & Z',IRTFLG)
           ELSE
              CALL RDPRIS(LX,LY,NOT_USED,'NEIGHBORHOOD X & Y',IRTFLG)
              LZ = 1
           ENDIF
           IF (IRTFLG .NE. 0) GOTO 9999

           IF (LX  <  3 .OR. MOD(LX,2)  ==  0 .OR. 
     &         LY  <  3 .OR. MOD(LY,2)  ==  0 .OR. (THREED .AND.
     &         LZ  <  3 .OR. MOD(LZ,2)  ==  0)) THEN
              CALL ERRT(101,'DIMENSIONS MUST BE ODD  AND > 2',IDUM) 
              GOTO 10
           ENDIF

           NEIGH = LX * LY * LZ

	   IF (MODE(1:2)  ==  'LH') THEN
C             LAHE
              NBINS = 64
              CALL RDPRI1S(NBINS,NOT_USED,'NUMBER OF BINS',IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999
              GOTO 1000

           ENDIF
        ELSEIF (MODE(1:1)  ==  'F') THEN
C          FREI-CHEN FILTER OVER 3x3 NEIGHBORHOOD
           IF (THREED) THEN
              WRITE(NOUT,*) ' THIS FILTER IS NOT IMPLEMENTED IN 3D!' 
              WRITE(NOUT,*) ' VOLUME WILL BE PROCESSED SLICE-BY-SLICE' 
              THREED = .FALSE.
              LZ     = 1
              NEIGH  = LX * LY 
           ENDIF
        ENDIF

C       ALLOCATE SPACE FOR KERNAL(S)
 	ALLOCATE(VKERNAL(NEIGH*NK),STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'CE, VKERNAL',IER)
            RETURN
        ENDIF

        IF (MODE(1:2)  ==  'GX') THEN
C          X GRADIENT OVER A 3X3 NEIGHBORHOOD
           IF (THREED) THEN
              VKERNAL = GRAD3X 
           ELSE
              VKERNAL = GRAD2X
           ENDIF

        ELSEIF (MODE(1:2)  ==  'GY') THEN
C          Y GRADIENT OVER A 3X3 NEIGHBORHOOD
           IF (THREED) THEN
              VKERNAL = GRAD3Y 
           ELSE
              VKERNAL = GRAD2Y
           ENDIF

        ELSEIF (MODE(1:2)  ==  'GZ') THEN
C          Z GRADIENT OVER A 3X3 NEIGHBORHOOD
           VKERNAL = GRAD3Z 

        ELSEIF (MODE(1:1)  ==  'G') THEN
C          GRADIENT OVER A 3X3 NEIGHBORHOOD
           IF (THREED) THEN 
              VKERNAL = GRAD3
           ELSE
              VKERNAL = GRAD2
           ENDIF
           VKERNAL(NEIGH/2 + 1) = - (NEIGH - 1)

        ELSEIF (MODE(1:1)  ==  'V') THEN
C          VARIANCE, VARIANCE SMOOTHING  OVER A RECTILINEAR NEIGHBORHOOD
           CONTINUE

        ELSEIF (MODE(1:1)  ==  'R' .OR.
     &          MODE(1:2)  ==  'MA' .OR. MODE(1:2)  ==  'MI') THEN
C          RANGE, RIDGE, MAX. OR MIN. OVER A RECTILINEAR NEIGHBORHOOD
           DO I = 1, NEIGH
              VKERNAL(I) = 1.0
           ENDDO
           VKERNAL(NEIGH/2 + 1) = 0.0

        ELSEIF (MODE(1:1)  ==  'L') THEN
C          LAPLACIAN  OVER A RECTILINEAR NEIGHBORHOOD
           DO I = 1, NEIGH
              VKERNAL(I) = -1.0
           ENDDO
           VKERNAL(NEIGH/2 + 1)   =  NEIGH - 1

        ELSEIF (MODE(1:1)  ==  'S') THEN
C          SOBEL OVER A 3X3 NEIGHBORHOOD
           IF (.NOT. THREED) THEN
              DO I = 1, NEIGH
                 VKERNAL(I)         = SOBX(I)
                 VKERNAL(I+NEIGH)   = SOBY(I)
              ENDDO
           ELSE
              DO I = 1, NEIGH
                 VKERNAL(I)         = SOBX3(I)
                 VKERNAL(I+NEIGH)   = SOBY3(I)
                 VKERNAL(I+2*NEIGH) = SOBZ3(I)
              ENDDO
           ENDIF

         ELSEIF (MODE(1:1)  ==  'P') THEN
C          PREWITT OVER A 3X3 NEIGHBORHOOD
           IF (.NOT. THREED) THEN
              DO I = 1, NEIGH
                 VKERNAL(I)         = PREX(I)
                 VKERNAL(I+NEIGH)   = PREY(I)
              ENDDO
           ELSE
              DO I = 1, NEIGH
                 VKERNAL(I)         = PREX3(I)
                 VKERNAL(I+NEIGH)   = PREY3(I)
                 VKERNAL(I+2*NEIGH) = PREZ3(I)
              ENDDO
           ENDIF

        ELSEIF (MODE(1:1)  ==  'T') THEN
C          TOP-HAT
           LXD2 = LX / 2
           LZD2 = LZ / 2
	   CALL FILTER_HAT(RADIN,RADOUT,LXD2,LZD2,NEIGH,
     &                     VKERNAL,THREED)

        ELSEIF (MODE(1:1)  ==  'H') THEN
C          HURST
           LXD2 = LX / 2
	   CALL FILTER_HURST(LXD2,THREED,NEIGH,ICAT,
     &             VKERNAL,VKERNAL(NEIGH+1))

        ELSEIF (MODE(1:1)  ==  'F') THEN
C          FREI-CHEN
	   CALL FILTER_FREI(NEIGH,VKERNAL)
        ENDIF

1000    CONTINUE
        ALLOCATE(VIN(NSAM*NROW*NSLICE),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'CE, VIN',IER)
           GOTO 9999
        ENDIF

        IF (MAXIM  <=  0) THEN
C          NOT A WHOLE STACK, LOAD INPUT VOLUME
           CALL REDVOL(LUN1,NSAM,NROW,1,NSLICE,VIN,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
        ENDIF

        LXD2  = LX / 2
        LYD2  = LY / 2
        LZD2  = LZ / 2

        WRITE(NOUT,90) NEIGH
90      FORMAT('  NEIGHBORHOOD: ',I4)

        IF (MODE(1:2)  ==  'LH') THEN
C          LAHE
           CALL FILTER_LAHE(VIN,NSAM,NROW,NSLICE,LXD2,LYD2,LZD2,
     &                     NEIGH,MODE,LUN2,FMINT,FMAXT,NBINS)

        ELSEIF ((NK  >  1) .AND. THREED) THEN
C          USES MULTIPLE KERNALS ON A WHOLE VOLUME
           CALL FILTER3M(VIN,NSAM,NROW,NSLICE,VKERNAL,LXD2,LYD2,LZD2,
     &                  NEIGH,MODE,NK,LUN2,ICAT)

        ELSEIF (NK  >  1) THEN
C          USES MULTIPLE KERNALS ON 2D IMAGE
           DO ISLICE=1,NSLICE
              ILOC = (ISLICE - 1) * NROW * NSAM + 1
              CALL FILTER2M(VIN(ILOC),NSAM,NROW,NSLICE,VKERNAL,
     &                      LXD2,LYD2,1,NEIGH,MODE,NK,LUN2,ICAT,ISLICE)
           ENDDO

        ELSEIF (THREED) THEN
C          OPERATION WORKS ON A WHOLE VOLUME
           CALL FILTER3(VIN,NSAM,NROW,NSLICE,VKERNAL,LXD2,LYD2,LZD2,
     &                  NEIGH,MODE,LUN2,FMINT,SIGT)

        ELSE
C          OPERATION WORKS ON 2D IMAGE OR SLICE-BY-SLICE THRU A VOLUME

           DO ISLICE=1,NSLICE
              ILOC = (ISLICE - 1) * NROW * NSAM + 1
              CALL FILTER2(VIN(ILOC),NSAM,NROW,NSLICE,VKERNAL,
     &                      LXD2,LYD2,NEIGH,MODE,LUN2,FMIN,ISLICE,SIGT)
           ENDDO
        ENDIF

9999    IF (ALLOCATED(VIN))     DEALLOCATE(VIN)
        IF (ALLOCATED(VKERNAL)) DEALLOCATE(VKERNAL)

        END


C       ------------------------- FILTER2M ----------------------------

	SUBROUTINE FILTER2M(VIN,NSAM,NROW,NSLICE,VKERNAL,LXD2,LYD2,LZD2,
     &                     NEIGH,MODE,NKT,LUN2,ICAT,IZ)

	REAL, DIMENSION(NSAM,NROW) :: VIN
	REAL, DIMENSION(-LXD2:LXD2,-LYD2:LYD2,NKT) :: VKERNAL

C       AUTOMATIC ARRAY
	REAL, DIMENSION(NSAM) :: VOUT
	REAL, DIMENSION(ICAT) :: VALSMIN,VALSMAX

	CHARACTER(LEN=2) ::      MODE

        NK = NKT 
        IF (MODE(1:1) == 'H') NK = 1

        DO IY=1,NROW  
                    
           DO IX=1,NSAM
              VALT  = 0.0
              IF (MODE(1:1) == 'R' .OR. MODE(1:1) == 'M' .OR.
     &            MODE(1:1) == 'S' .OR. MODE(1:1) == 'P') THEN
C                RANGE, ETC, USE MIN OR MAX OR BOTH
                 VALMIN = HUGE(VALMIN)
                 VALMAX = -VALMIN

              ELSEIF (MODE(1:1)  ==  'H') THEN
C                HURST,  ZERO DIFFERENCE ARRAYS
                 VALSMIN =  HUGE(VALMIN) 
                 VALSMAX = -VALSMIN 
 
              ELSEIF (MODE(1:1)  ==  'F') THEN
C                FREI CHEN NEEDS TWO SUMS
                 VAL1 = 0.0
                 VAL2 = 0.0
              ENDIF

C             APPLY  MULTIPLE KERNALS
              DO KER = 1,NK
                 IF (MODE  ==  'F') VALT = 0.0
                 DO MY=-LYD2,LYD2
                    IYT = MOD(IY+MY+NROW-1,NROW)+1

                    DO MX=-LXD2,LXD2
C                      VALUE FOR PIXEL UNDER CURRENT KERNAL ELEMENT
                       VOK = VIN(MOD(IX+MX+NSAM-1,NSAM)+1,IYT)

                       IF (MODE(1:1)  ==  'H') THEN
C                         HURST
                          IVOK = VKERNAL(MX,MY,1)
                          IF (IVOK  >  0) THEN
                             VALSMIN(IVOK) = MIN(VALSMIN(IVOK),VOK) 
                             VALSMAX(IVOK) = MAX(VALSMAX(IVOK),VOK) 
                          ENDIF
                       ELSE
                          VALT = VALT + VOK * VKERNAL(MX,MY,KER)
                       ENDIF

                    ENDDO
C                   END LOOP: MX=-LXD2,LXD2 

                 ENDDO
C                END LOOP: DO MY=-LYD2,LYD2 

                 IF (MODE(1:1)  ==  'S' .OR. MODE(1:1)  ==  'P') THEN
C                   SOBEL OR PREWITT
                    VALMAX = MAX(VALMAX,VALT)

                 ELSEIF (MODE(1:1)  ==  'F') THEN
C                   FREI-CHEN
                    IF (KER .GE. 2 .AND. KER  <=  3) THEN
                       VAL1 = VAL1 + VALT**2
                    ELSE
                       VAL2 = VAL2 + VALT**2
                    ENDIF
                 ENDIF

              ENDDO
C             END LOOP:  DO KER = 1,NK

              IF (MODE(1:1)  ==  'S' .OR. MODE(1:1)  ==  'P') THEN
C                SOBEL OR PREWITT
                 VOUT(IX) = VALMAX

              ELSEIF (MODE(1:1)  ==  'H') THEN
C                HURST
	         CALL FILTER_HURST_DO(VKERNAL(-LXD2,-LYD2,2),ICAT,
     &               VALSMIN,VALSMAX,VOUT(IX))

               ELSEIF (MODE(1:1)  ==  'F') THEN
C                FREI-CHEN
                 IF (VAL2 .NE. 0.0) THEN
                    VOUT(IX) = COS(SQRT(VAL1/VAL2))
                 ELSE
C                   I DO NOT KNOW IF THIS IS OK??? al
                    VOUT(IX) = 0.0
                 ENDIF
             ENDIF

           ENDDO
C          END LOOP:  DO IX=1,NSAM
 
C          OUTPUT IMAGE
           IREC = (IZ - 1) * NROW + IY
           CALL WRTLIN(LUN2,VOUT,NSAM,IREC)
         ENDDO

         END	

C       ------------------------- FILTER3M ----------------------------

	SUBROUTINE FILTER3M(VIN,NSAM,NROW,NSLICE,VKERNAL,LXD2,LYD2,LZD2,
     &                     NEIGH,MODE,NKT,LUN2,ICAT)

	REAL, DIMENSION(NSAM,NROW,NSLICE) :: VIN
	REAL, DIMENSION(-LXD2:LXD2,-LYD2:LYD2,-LZD2:LZD2,NKT) :: VKERNAL

C       AUTOMATIC ARRAY
	REAL, DIMENSION(NSAM)      :: VOUT

	CHARACTER(LEN=2) ::           MODE 

        NK = NKT 
        IF (MODE(1:1)  ==  'H') NK = 1

        DO IZ=1,NSLICE

        DO IY=1,NROW  
                    
           DO IX=1,NSAM
              VALT = 0.0
              IF (MODE(1:1)  ==  'R' .OR. MODE(1:1)  ==  'M' .OR.
     &            MODE(1:1)  ==  'S' .OR. MODE(1:1)  ==  'P') THEN
C                RANGE, ETC, USE MIN OR MAX OR BOTH
                 VALMIN = HUGE(VALMIN)
                 VALMAX = -VALMIN
              ENDIF

C             APPLY  MULTIPLE KERNALS
              DO KER = 1,NK
                 DO MZ=-LZD2,LZD2
                    IZT = MOD(IZ+MZ+NSLICE-1,NSLICE)+1 
 
                    DO MY=-LYD2,LYD2
                       IYT = MOD(IY+MY+NROW-1,NROW)+1

                      DO MX=-LXD2,LXD2
C                         VALUE FOR PIXEL UNDER CURRENT KERNAL ELEMENT
                          VOK = VIN(MOD(IX+MX+NSAM-1,NSAM)+1,IYT,IZT)

                          VALT = VALT + VOK * VKERNAL(MX,MY,MZ,KER)
                       ENDDO
C                      END LOOP: MX=-LXD2,LXD2 
                    ENDDO
C                   END LOOP: DO MY=-LYD2,LYD2 
                 ENDDO

                 IF (MODE(1:1)  ==  'S' .OR. MODE(1:1)  ==  'P') THEN
C                   SOBEL OR PREWITT
                    VALMAX = MAX(VALMAX,VALT)
                 ENDIF

              ENDDO
C             END LOOP:  DO KER = 1,NK

              IF (MODE(1:1)  ==  'S' .OR. MODE(1:1)  ==  'P') THEN
C                SOBEL OR PREWITT
                 VOUT(IX) = VALMAX
              ENDIF

           ENDDO
C          END LOOP:  DO IX=1,NSAM
 
C          OUTPUT IMAGE
           IREC = (IZ - 1) * NROW + IY
           CALL WRTLIN(LUN2,VOUT,NSAM,IREC)
         ENDDO
         ENDDO

         END	

C       ------------------------- FILTER3 -----------------------------

	SUBROUTINE FILTER3(VIN,NSAM,NROW,NSLICE,VKERNAL,LXD2,LYD2,LZD2,
     &                     NEIGH,MODE,LUN2,FMINT,SIGT)

	INCLUDE 'CMBLOCK.INC'

	REAL, DIMENSION(NSAM,NROW,NSLICE)                 :: VIN
	REAL, DIMENSION(-LXD2:LXD2,-LYD2:LYD2,-LZD2:LZD2) :: VKERNAL

C       AUTOMATIC ARRAY
        REAL, DIMENSION(NEIGH) :: FMED
	REAL, DIMENSION(NSAM)  :: VOUT

	CHARACTER(LEN=2)       :: MODE 

        IF (MODE(1:1)  ==  'V') THEN
C          VARIANCE OR VARIANCE SMOOTHING
           CONAVG  = 1.0 / FLOAT(NEIGH)
           CONVAR  = 1.0 / FLOAT(NEIGH -1)
           SIGTSQ  = SIGT**2
           KMED    = NEIGH / 2 + 1
           NREPL   = 0
        ENDIF

        DO IZ=1,NSLICE
        DO IY=1,NROW  
                    
           DO IX=1,NSAM
              VALT = 0.0

              IF (MODE(1:1)  ==  'V') THEN
C                VARIANCE
                 VALAVG  = 0
                 VOC     = VIN(IX,IY,IZ)    
                 LMED    = 0

              ELSEIF (MODE(1:1)  ==  'R' .OR. MODE(1:1)  ==  'M' .OR.
     &                MODE(1:1)  ==  'T') THEN
C                RANGE, RIDGE, MIN, MAX, OR TOP-HAT
                 VALMIN  = HUGE(VALMIN)
                 VALMAX  = -VALMIN
                 VALMAX2 = -VALMIN
                 VOC     = VIN(IX,IY,IZ)    
              ENDIF

C             APPLY KERNAL
              DO MZ=-LZD2,LZD2
              IZT = MOD(IZ+MZ+NSLICE-1,NSLICE)+1 
 
              DO MY=-LYD2,LYD2
                 IYT = MOD(IY+MY+NROW-1,NROW)+1

                 DO MX=-LXD2,LXD2
C                   VALUE FOR IMAGE UNDER CURRENT KERNAL ELEMENT
                    VOK = VIN(MOD(IX+MX+NSAM-1,NSAM)+1,IYT,IZT)

                    IF (MODE(1:1)  ==  'V') THEN
C                      USE SQ. FOR VARIANCE
                       VALT   = VALT   + VOK ** 2
                       VALAVG = VALAVG + VOK 

                       IF (MODE(2:2)  ==  'S') THEN
C                          NEED TO GET MEDIAN FOR VARIANCE SMOOTHING
                           LMED       = LMED + 1
                           FMED(LMED) = VOK
                       ENDIF
 
                    ELSEIF (MODE(1:1)  ==  'R' .OR. 
     &                      MODE(1:1)  ==  'M') THEN
C                      RANGE, RIDGE, MIN., OR MAX.
                       VALMIN = MIN(VALMIN,VOK)
                       VALMAX = MAX(VALMAX,VOK)

                    ELSEIF (MODE(1:1)  ==  'T') THEN
C                      TOP-HAT
                       IF (VKERNAL(MX,MY,MZ)  ==  1) THEN
                          VALMAX = MAX(VALMAX,VOK)
                       ELSEIF (VKERNAL(MX,MY,MZ)  ==  2) THEN
                          VALMAX2 = MAX(VALMAX2,VOK)
                       ENDIF

                    ELSE
C                      USE SUM OF KERNAL PRODUCTS
                       VALT = VALT + VOK * VKERNAL(MX,MY,MZ)
                    ENDIF
                 ENDDO
              ENDDO
              ENDDO

              IF (MODE(1:2)  ==  'RA') THEN
C                RANGE
                 VOUT(IX) = VALMAX - VALMIN

              ELSEIF (MODE(1:1)  ==  'R') THEN
C                RIDGE
                 IF (VOC  <  VALMIN .OR. VALMAX  >  VOC) THEN
                     VOUT(IX) = FMINT
                 ELSE
                     VOUT(IX) = VOC
                 ENDIF
              ELSEIF (MODE(1:2)  ==  'MA') THEN
C                MAX
                 VOUT(IX) = VALMAX 

              ELSEIF (MODE(1:2)  ==  'MI') THEN
C                MIN
                 VOUT(IX) = VALMIN 

              ELSEIF (MODE(1:1)  ==  'T') THEN
C                TOP-HAT
                 VOUT(IX) = VALMAX2 - VALMAX 

              ELSEIF (MODE(1:2)  ==  'VS') THEN
C                VARIANCE SMOOTHING
                 SLOC = ABS(CONVAR * 
     &                  (VALT - CONAVG *(VALAVG * VALAVG)))

                 IF (((VOC - VALAVG)**2)  >  SLOC .AND. 
     &                SLOC  >  SIGTSQ) THEN
C                    USE MEDIAN FROM NEIGHBORING VALUES
                     CALL FSORT(FMED,LMED)
                     VOUT(IX) = FMED(KMED)
                     NREPL    = NREPL + 1
                 ELSE
C                    KEEP SAME VALUE
                     VOUT(IX) = VOC
                 ENDIF

              ELSEIF (MODE(1:1)  ==  'V') THEN
C                VARIANCE
                 VOUT(IX) =  CONVAR * (VALT - CONAVG *(VALAVG * VALAVG)) 
              ELSE
                 VOUT(IX) = VALT
              ENDIF
 
          ENDDO

C          OUTPUT VOLUME LINE
           IREC = (IZ - 1) * NROW + IY
           CALL WRTLIN(LUN2,VOUT,NSAM,IREC)
         ENDDO
         ENDDO

         IF (MODE(1:2)  ==  'VS') THEN
C           VARIANCE SMOOTHING

            FREPL = FLOAT(NREPL)
            CALL REG_SET_NSEL(1,1, FNREPL,0.0, 0.0, 0.0, 0.0,IRTFLG)

            WRITE(NOUT,*) ' VOXELS REPLACED: ',NREPL
         ENDIF

         END	


C       ------------------------- FILTER2 -----------------------------

	SUBROUTINE FILTER2(VIN,NSAM,NROW,NSLICE,VKERNAL,LXD2,LYD2,NEIGH,
     &                      MODE,LUN2,FMINT,IZ,SIGT)

	INCLUDE 'CMBLOCK.INC'

	REAL, DIMENSION(NSAM,NROW)             :: VIN
	REAL, DIMENSION(-LXD2:LXD2,-LYD2:LYD2) :: VKERNAL

C       AUTOMATIC ARRAYS
	REAL, DIMENSION(NSAM)                  :: VOUT
        REAL, DIMENSION(NEIGH )                :: FMED

	CHARACTER(LEN=2)                       :: MODE 

        IF (MODE(1:1)  ==  'V') THEN
C          VARIANCE
           CONAVG   = 1.0 / FLOAT(NEIGH)
           CONVAR   = 1.0 / FLOAT(NEIGH - 1)
           SIGTSQ   = SIGT**2
           KMED     = NEIGH / 2 + 1
           NREPL    = 0
        ENDIF

        DO IY=1,NROW  
                    
           DO IX=1,NSAM
              VALT = 0.0

              IF (MODE(1:1)  ==  'V') THEN
C                VARIANCE
                 VALAVG  = 0
                 VOC     = VIN(IX,IY)    
                 LMED    = 0

             ELSEIF (MODE(1:1)  ==  'R' .OR. MODE(1:1)  ==  'M' .OR.
     &                MODE(1:1)  ==  'T') THEN
C                RANGE, RIDGE, MIN, MAX, OR TOP-HAT
                 VALMIN  = HUGE(VALMIN)
                 VALMAX  = -VALMIN
                 VALMAX2 = -VALMIN
                 VOC     = VIN(IX,IY)    
              ENDIF

C             APPLY KERNAL 
              DO MY=-LYD2,LYD2
                 IYT = MOD(IY+MY+NROW-1,NROW)+1

                 DO MX=-LXD2,LXD2
C                   VALUE FOR IMAGE UNDER CURRENT KERNAL ELEMENT
                    VOK = VIN(MOD(IX+MX+NSAM-1,NSAM)+1,IYT)

                    IF (MODE(1:1)  ==  'V') THEN
C                      USE SQ. FOR VARIANCE
                       VALT   = VALT   + VOK ** 2
                       VALAVG = VALAVG + VOK 

                       IF (MODE(2:2)  ==  'S') THEN
C                          NEED TO GET MEDIAN FOR VARIANCE SMOOTHING
                           LMED       = LMED + 1
                           FMED(LMED) = VOK
                       ENDIF
 
                    ELSEIF (MODE(1:1)  ==  'R' .OR. 
     &                      MODE(1:1)  ==  'M') THEN
C                      RANGE, RIDGE, MIN., OR MAX.
                       VALMIN = MIN(VALMIN,VOK)
                       VALMAX = MAX(VALMAX,VOK)

                    ELSEIF (MODE(1:1)  ==  'T') THEN
C                      TOP-HAT
                       IF (VKERNAL(MX,MY)  ==  1) THEN
C                         OUTER RING
                          VALMAX = MAX(VALMAX,VOK)

                       ELSEIF (VKERNAL(MX,MY)  ==  2) THEN
C                         INNER RING
                          VALMAX2 = MAX(VALMAX2,VOK)
                       ENDIF

                    ELSE
C                      USE SUM OF KERNAL * IMAGE PRODUCTS
                       VALT = VALT + VOK * VKERNAL(MX,MY)
                    ENDIF
                 ENDDO
              ENDDO

              IF (MODE(1:2)  ==  'RA') THEN
C                RANGE
                 VOUT(IX) = VALMAX - VALMIN

              ELSEIF (MODE(1:1)  ==  'R') THEN
C                RIDGE
                 IF (VOC  <  VALMIN .OR. VALMAX  >  VOC) THEN
                     VOUT(IX) = FMINT
                 ELSE
                     VOUT(IX) = VOC
                 ENDIF
              ELSEIF (MODE(1:2)  ==  'MA') THEN
C                MAX
                 VOUT(IX) = VALMAX 

              ELSEIF (MODE(1:2)  ==  'MI') THEN
C                MIN
                 VOUT(IX) = VALMIN 

              ELSEIF (MODE(1:1)  ==  'T') THEN
C                TOP-HAT
                 VOUT(IX) = VALMAX2 - VALMAX 

              ELSEIF (MODE(1:2)  ==  'VS') THEN
C                VARIANCE SMOOTHING
C                 VALT   = VALT   - VOC ** 2
C                 VALAVG = VALAVG - VOC 

                 SLOC = ABS(CONVAR * 
     &                  (VALT - CONAVG *(VALAVG * VALAVG)))
                 VT   = VALAVG * CONAVG

                 IF (((VOC - VT)**2)  >  SLOC .AND. 
     &                SLOC  >  SIGTSQ) THEN

C                    USE MEDIAN FROM NEIGHBORING VALUES
                     CALL FSORT(FMED,LMED)
                     VOUT(IX) = FMED(KMED)
C                    VOUT(IX) = (VALAVG - VOC) * CONVAR

CCC                  WRITE(NOUT,*) ' REPLACED(',IX,',',IY,'): ',VOC,
CCC  &                         '  WITH: ',VOUT(IX)

                     NREPL    = NREPL + 1
                 ELSE
C                    KEEP SAME VALUE
                     VOUT(IX) = VOC
                 ENDIF

              ELSEIF (MODE(1:1)  ==  'V') THEN
C                VARIANCE
                 VOUT(IX) =  CONVAR * (VALT - CONAVG*(VALAVG * VALAVG)) 

              ELSE
                 VOUT(IX) = VALT
              ENDIF
           ENDDO

C          OUTPUT IMAGE LINE
           IREC = (IZ - 1) * NROW + IY
           CALL WRTLIN(LUN2,VOUT,NSAM,IREC)
         ENDDO

         IF (MODE(1:2)  ==  'VS') THEN
C           VARIANCE SMOOTHING

            FNREPL = FLOAT(NREPL)
            CALL REG_SET_NSEL(1,1, FNREPL,0.0, 0.0, 0.0, 0.0,IRTFLG)

            WRITE(NOUT,*) ' PIXELS REPLACED: ',NREPL
         ENDIF
         END	



C       ------------------------- FILTER_HAT ----------------------------

	SUBROUTINE FILTER_HAT(RADIN,RADOUT,LXD2,LZD2,NEIGH,
     &                        VKERNAL,THREED)

	REAL, DIMENSION(NEIGH) :: VKERNAL

        LOGICAL :: THREED

        EPS = EPSILON(EPS)

        ILOC = 0

        DO MZ=-LZD2,LZD2
           DO MY=-LXD2,LXD2
              DO MX=-LXD2,LXD2
                 DIST = SQRT(FLOAT(MX)**2 + FLOAT(MY)**2 + FLOAT(MZ)**2)
                 ILOC = ILOC + 1

                 IF (DIST   <  RADIN) THEN
C                   INNER REGION
                    VKERNAL(ILOC) = 2

                 ELSEIF (DIST   <  RADOUT) THEN
C                   OUTER REGION
                    VKERNAL(ILOC) = 1

                 ELSE
C                   BORDER REGION
                    VKERNAL(ILOC) = 0
                 ENDIF
              ENDDO
           ENDDO
        ENDDO

        END


C       ------------------------- FILTER_FREI ----------------------------

	SUBROUTINE FILTER_FREI(NEIGH,VKERNAL)

	REAL, DIMENSION(NEIGH*8) :: VKERNAL

        REAL, DIMENSION(9) :: FC0,FC1,FC2,FC3,FC4,FC5,FC6,FC7,FC8

        DATA FC0/
     &        1, 1, 1,
     &        1, 1, 1,
     &        1, 1, 1/

        DATA FC1/
     &       -1,-1.414,-1,
     &        0,   0,   0,
     &        1, 1.414, 1/

        DATA FC2/
     &      -1,     0, -1,
     &      -1.414, 0, -1.414,
     &      -1,     0,  1/

        DATA FC3/
     &       0,     -1, -1.414,
     &       1,      0, -1 ,
     &      -1.4141, 1,  0/

        DATA FC4/
     &      -1.414, -1,  0,
     &      -1,      0, -1 ,
     &      0 ,      1, -1.414/
        DATA FC5/
     &        0, 1, 0,
     &       -1, 1, 1,
     &        0,-1, 0/

        DATA FC6/
     &       -1, 0, 1,
     &        0, 0, 0,
     &        1, 0,-1/

        DATA FC7/
     &        1,-2, 1,
     &       -2, 4,-2,
     &        1,-2, 1/

        DATA FC8/
     &       -2, 1,-2,
     &        1, 4, 1,
     &       -2, 1,-2/
 
        DO I = 1, NEIGH
           VKERNAL(I+NEIGH*0) = FC0(I)
           VKERNAL(I+NEIGH*1) = FC1(I)
           VKERNAL(I+NEIGH*2) = FC2(I)
           VKERNAL(I+NEIGH*3) = FC3(I)
           VKERNAL(I+NEIGH*4) = FC4(I)
           VKERNAL(I+NEIGH*5) = FC5(I)
           VKERNAL(I+NEIGH*6) = FC6(I)
           VKERNAL(I+NEIGH*7) = FC7(I)
           VKERNAL(I+NEIGH*8) = FC8(I)
        ENDDO

        END

C       ------------------------- FILTER_HURST ----------------------------

	SUBROUTINE FILTER_HURST(LXD2,THREED,NEIGH,ICAT,VKERNAL,VLOGDIST)

C       SURE, THIS IS SLOW, BUT IT IS NOT IMPORTANT!

	REAL, DIMENSION(NEIGH) :: VKERNAL,VLOGDIST

C       AUTOMATIC ARRAYS
        INTEGER, DIMENSION((LXD2+1)**2) :: IUNIQ,IUNIQSORT

C       THREED NOT IMPLEMENTED!!!
        LOGICAL :: THREED,NEWDIST

        LEND  = (2*LXD2+1)**2

C       FIND UNIQUE DISTANCES**2 TO ANY KERNAL ELEMENT
        ICAT = 0
        ILOC = 0
        DO MY=0,LXD2
           DO MX=0,MY
              IVAL     = MX**2 + MY**2
              ILOC     = ILOC + 1 
              NEWDIST  = .TRUE.
              IF (ICAT  >  0) THEN
                 DO I = 1,ICAT
                    IF (IUNIQ(I)  ==  IVAL) THEN
                       NEWDIST = .FALSE.
                       EXIT
                    ENDIF
                 ENDDO
              ENDIF
              IF (NEWDIST .AND. IVAL  >  0) THEN
                 ICAT        = ICAT + 1
                 IUNIQ(ICAT) = IVAL
              ENDIF
           ENDDO
        ENDDO

C       SET KERNAL ELEMENTS TO CATEGORIES
        VKERNAL = 0
        ILOC    = 0
        DO MY = -LXD2,LXD2
           DO MX = -LXD2,LXD2
              ILOC = ILOC + 1
              IVAL = MX**2 + MY**2
              DO I=1,ICAT
                 IF (IVAL  ==  IUNIQ(I)) THEN
                    VKERNAL(ILOC) = I
                    VLOGDIST(I)   = LOG(SQRT(FLOAT(IVAL)))
                    EXIT
                 ENDIF
              ENDDO
           ENDDO
        ENDDO

        END

C       ------------------------- FILTER_HURST_DO -----------------------

	SUBROUTINE FILTER_HURST_DO(VLOGDIST,ICAT,VALSMIN,VALSMAX,VAL)

C       SURE, THIS IS SLOW, BUT IT IS NOT IMPORTANT!

	REAL, DIMENSION(ICAT) :: VALSMIN,VALSMAX,VLOGDIST
	REAL, DIMENSION(2)    :: COEF

        EPS = EPSILON(EPS)

C       FIND RANGE FOR ALL CATEGORIES
        DO I = 1,ICAT
          IF ((VALSMAX(I) - VALSMIN(I))  <=  EPS) THEN
              VAL = 0.0
              RETURN
           ENDIF
           VALSMAX(I) = LOG(VALSMAX(I) - VALSMIN(I))
        ENDDO

C       SUBROUTINE POLFIT(X,Y,NORD,N,C)
C       PURPOSE:  MAKES LEAST SQUARES FIT OF EXPERIMENTAL DATA IN
C                 X(N),Y(N) USING A POLYNOMIAL OF ARBITRARY ORDER>1
C       X     A REAL ARRAY DIMENSIONED N CONTAINING THE ABSCISSAE
C       Y     A REAL ARRAY DIMENSIONED N CONTAINING THE FUNCTION VALUES
C       NORD  ORDER OF POLYNOMIAL. NORD<N.
C       C     A REAL ARRAY CONTAINING THE NORD+1 COEFFICIENTS OF THE
C               POLYNOMIAL IN INCREASING ORDER

        CALL POLFIT(VLOGDIST,VALSMAX,1,ICAT,COEF)

C       RETURN SLOPE
        VAL = COEF(1)
        END
        

C       ------------------------- FILTER_LAHE -------------------------

	SUBROUTINE FILTER_LAHE(VIN,NSAM,NROW,NSLICE,LXD2,LYD2,LZD2,
     &                     NEIGH,MODE,LUN2,FMINT,FMAXT,NBINS)

	REAL, DIMENSION(NSAM,NROW,NSLICE) :: VIN

C       AUTOMATIC ARRAYS
	REAL, DIMENSION(NSAM)     :: VOUT
	INTEGER, DIMENSION(NBINS) :: HIST

	CHARACTER(LEN=2)          :: MODE 

C       HISTOGRAM BIN SIZE
        BINSPVAL = (NBINS - 1) / (FMAXT - FMINT)
        VALPPIX  = (FMAXT - FMINT) / NEIGH

        DO IZ=1,NSLICE
           DO IY=1,NROW  
              DO IX=1,NSAM

C                ZERO THE HIST ARRAY
                 HIST = 0

C                APPLY KERNAL
                 DO MZ=-LZD2,LZD2
                    IZT = MOD(IZ+MZ+NSLICE-1,NSLICE)+1 
 
                    DO MY=-LYD2,LYD2
                       IYT = MOD(IY+MY+NROW-1,NROW)+1

                       DO MX=-LXD2,LXD2
C                         VALUE FOR IMAGE UNDER CURRENT KERNAL ELEMENT
                          VOK = VIN(MOD(IX+MX+NSAM-1,NSAM)+1,IYT,IZT)
                          IBIN       = (VOK - FMINT) * BINSPVAL + 1.5
                          HIST(IBIN) = HIST(IBIN) + 1
                       ENDDO
                    ENDDO
                 ENDDO

C                NUMBER OF BINNED VALUES NEEDED TO GET CURRENT VALUE
C                ITARG=(VIN(IX,IY,IZ)-FMINT)*(NBINS-1)/(FMAXT-FMINT)+1.5
                 ITARG = (VIN(IX,IY,IZ) - FMINT) * BINSPVAL + 1.5
                 IGOT  = 0

                 DO I = 1,ITARG
                    IGOT = IGOT + HIST(I)
                 ENDDO

C                HISTOGRAM EQUALIZED OUTPUT VALUE
C                VOUT(IX) = (IGOT - 1) * ((FMAXT-FMINT)/NEIGH) + FMINT
                 VOUT(IX) = (IGOT - 1) * VALPPIX + FMINT
              ENDDO

C             OUTPUT IMAGE LINE
              IREC = (IZ - 1) * NROW + IY
              CALL WRTLIN(LUN2,VOUT,NSAM,IREC)
            ENDDO
         ENDDO

         END	
