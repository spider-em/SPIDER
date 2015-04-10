
C++*********************************************************************
C  MAKE_PATM    NEW                       MAY 11    GREGORY KISHCHENKO *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C MAKE_PATM(AVG,SD,BUFFFT, LMASK, NX,NY,  N2XLD,N2X,N2Y, 
C           BUFWORK,BUFPATM, IRTFLG)
C
C PURPOSE: RETURN ENLARGED PATTERSON MAP OF BUFIMG.
C          NORMALIZE IMAGE WITHIN MASK.
C          PAD MASKED IMAGE INTO BUFFFT,  INPLACE FORWARD FFT. 
C          CREATE NORMALIZED PATTERSON MAP, SHIFT MAP TO CENTER.
C
C PARAMETERS:
C          AVG             IMAGE AVERAGE                        INPUT
C          SD              IMAGE SD                             INPUT
C          BUFFFT          PADDED IMAGE AVERAGE                 INPUT
C                          FFT OF PADDED IMAGE AVERAGE          OUTPUT
C          LMASK           CIRCULAR LOGICAL MASK                INPUT
C          NX,NY           IMAGE DIMENSIONS                     INPUT        
C          N2XLD,N2X,N2Y   PADDED IMAGE DIMENSIONS              INPUT
C          BUFWORK         WORKING ARRAY FOR CENTERING          OUTPUT
C          BUFPATM         PATTERSON MAP                        OUTPUT
C          IRTFLG          ERROR FLAG                           OUTPUT
c
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************
 
       SUBROUTINE MAKE_PATM(AVG,SD,BUFFFT, CMASK, 
     &                      NX,NY, N2XLD,N2X,N2Y, 
     &                      BUFWORK, BUFPATM, IRTFLG)

       IMPLICIT NONE

       REAL,    INTENT(IN)   :: AVG,SD                    
       INTEGER, INTENT(IN)   :: NX,NY, N2XLD,N2X,N2Y
       REAL,    INTENT(INOUT):: BUFFFT (N2XLD,N2Y)       
       REAL,    INTENT(IN)   :: CMASK  (NX,NY)          
       REAL,    INTENT(OUT)  :: BUFWORK(N2XLD,N2Y)       
       REAL,    INTENT(OUT)  :: BUFPATM(N2X,  N2Y)       
       INTEGER, INTENT(OUT)  :: IRTFLG

       REAL                  :: SDINV                    
       DOUBLE PRECISION      :: DAVPAT,DSIGPAT,DSIGPATINV 
       INTEGER               :: INV,IX,IY,IX2,IY2
       INTEGER               :: IDX,IDY
       INTEGER * 8           :: IPLAN1 = 0     ! STRUCTURE POINTER 
       INTEGER * 8           :: IPLAN2 = 0     ! STRUCTURE POINTER 
       INTEGER               :: INRAD,I,J,  nmax

       LOGICAL               :: SPIDER_SCALE = .FALSE.
       LOGICAL               :: SPIDER_SIGN  = .FALSE.
       LOGICAL               :: USE_OMP      = .TRUE.

      real    :: fmint,fmaxt,flow,fhi,rangenew,rangeold,con,con2 

       !call chkfile('jnkpat00',66,1,n2xld,n2y,1, buffft,irtflg)

C      THE 2X SIZE IS FILLED WITH ZEROS OUTSIDE IMAGE BY GETDATA
C      NORMALIZE IMAGE by SUTRACTING AVG.      (IMPROVES RESULTS) 
C      MASK IMAGE WITH ZEROS                   (IMPROVES RESULTS) 
C      MASK IMAGE WITH AVG                     (POORER   RESULTS) 
C      FOR PWS BEST IF SUBTRACT AVG FROM IMAGE,
C      SET OUTSIDE IMAGE=0, PAD TO 2X WITH AVG,
C      MASK IMAGE WITH ZEROS                   (PWS POOR RESULTS) 
C      BUFFFT(IX,IY) =  (BUFFFT(IX,IY) - AVG) * SDINV (SAME RESULTS)

       !sdinv = 1.0 / sd       ! for speed

c$omp  parallel do private(iy,ix)
       DO IY = 1,NY
          DO IX = 1,NX
              BUFFFT(IX,IY) = (BUFFFT(IX,IY) - AVG ) * CMASK(IX,IY)
          ENDDO
       ENDDO

       !call chkfile('jnkpat01',66,1,n2xld,n2y,1, buffft,irtflg
 
       ! DO NOT USE: FMRS_2() AS IT ALWAYS SCALES!
       INV = +1   ! INPLACE FORWARD FFT 
       CALL FMRS(BUFFFT, N2X,N2Y,1, IPLAN1,
     &            SPIDER_SIGN,SPIDER_SCALE, INV,IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

C      CREATE PATTERSON MAP IN FOURIER SPACE
c$omp  parallel do private(iy,ix)
       DO IY=1,N2Y
          DO IX=0,N2XLD/2-1
             BUFWORK(2*IX+1, IY) = 
     &          ((BUFFFT(2*IX+1, IY)**2 + BUFFFT(2*IX+2, IY)**2)) 
             BUFWORK(2*IX+2, IY) = 0.0

             ! SQRT LINE MAKES RESULTS WORSE!!
             !BUFWORK(2*IX+1, IY) = SQRT(BUFWORK(2*IX+1,IY)) ! POWER->AMPLITUDE
          ENDDO
       ENDDO

       INV = -1   ! INPLACE REVERSE FFT 
       CALL FMRS(BUFWORK, N2X,N2Y,1, IPLAN2, 
     &           SPIDER_SIGN,SPIDER_SCALE, INV,IRTFLG)

       !call chkfile('jnkpat02',66,1,n2xld,n2y,1, bufwork,irtflg)

C      SHIFT NORMALIZED PATTERSON MAP TO SPIDER IMAGE CENTER
C      X-CENTER = N2X/2+1,    Y-CENTER = N2Y/2+1

       IDX    = N2X / 2 - 1
       IDY    = N2Y / 2 - 1

c$omp  parallel do private(iy,iy2,ix,ix2)
       DO IY=1,N2Y
          IY2 = MODULO(IY + IDY, N2Y) + 1

          DO IX=1,N2X
             IX2              = MODULO(IX + IDX, N2X) + 1
             BUFPATM(IX2,IY2) = BUFWORK(IX, IY)
          ENDDO
       ENDDO

C      CENTRAL POINT IS VERY BRIGHT, SO THIS MAKES USEFULL FOR DISPLAY
       BUFPATM(IDX,IDY) = BUFPATM(IDX+1,IDY+1)

       !call chkfile('jnkpat03',66,1,n2x,n2y,1, bufpatm,irtflg)
       !call chkmaxloc2d(' PATM max',bufpatm,n2x,n2y)
       
       END



C      ---------------------- MAKE_CIRC_LMASK -----------------------

       SUBROUTINE MAKE_CIRC_LMASK(LMASK, NX,NY,
     &                      IXCEN,IYCEN,IRAD,  
     &                      LIN,LOUT)

       IMPLICIT NONE
       LOGICAL     :: LMASK(NX,NY)  ! LOGICAL MASK (OUTPUT)
       INTEGER     :: NX,NY, IXCEN,IYCEN,IRAD
       LOGICAL     :: LIN,LOUT      ! INTERNAL & EXTERNAL L VALUES

       INTEGER     :: IY,IYDSQ,IX,IRADSQ

C      CREATE MASK IMAGE WITH LOUT OUTSIDE IRAD FROM (IXCEN,IYCEN)
C  
       IRADSQ = IRAD * IRAD + 1

       DO IY=1,NY

          IYDSQ = (IY - IYCEN)**2

          DO IX=1,NX
             IF (((IX - IXCEN)**2 + IYDSQ) <= IRADSQ) THEN
                LMASK(IX,IY) = LIN
             ELSE
                LMASK(IX,IY) = LOUT
             ENDIF
          ENDDO
       ENDDO

       END


C      ---------------------- MAKE_CIRC_CMASK -----------------------

       SUBROUTINE MAKE_CIRC_CMASK(CMASK, NX,NY, XCEN,YCEN, RAD)

       IMPLICIT NONE
       REAL        :: CMASK(NX,NY)    ! REAL MASK (OUTPUT)
       INTEGER     :: NX,NY
       REAL        :: XCEN,YCEN, RAD

       INTEGER     :: IY,IX
       REAL        :: RADSQ,YRADSQ,RADNOW

C      CREATE CIRCULAR MASK IMAGE WITH CENTER: (XCEN,YCEN)

        RADSQ = RAD**2

        DO IY = 1,NY
           YRADSQ = (IY-YCEN)**2 

           DO IX = 1,NX
	      RADNOW = (IX-XCEN)**2  + YRADSQ

	      IF (RADNOW >= RADSQ)  THEN
	         CMASK(IX,IY) = 0.0
	      ELSE
                 CMASK(IX,IY) = 1.0
	      ENDIF
           ENDDO
       ENDDO

       END


C      ---------------------- MAKE_CIRC_CSMASK -----------------------

       SUBROUTINE MAKE_CIRC_CSMASK(CMASK, NX,NY, XCEN,YCEN, RAD)

       IMPLICIT NONE
       REAL        :: CMASK(NX,NY)    ! REAL MASK (OUTPUT)
       INTEGER     :: NX,NY
       REAL        :: XCEN,YCEN, RAD

       INTEGER     :: NORDER,IY,IX
       REAL        :: TNM,RINV_RADSQ,YRADSQ,RADNOW,VALG

C      CREATE SOFT CIRCULAR MASK IMAGE WITH CENTER: (XCEN,YCEN)

C      SET THE ORDER OF SUPERGAUSSIAN 
       NORDER     = 3
       TNM        = ALOG(1.0 / TINY(TNM))

       RINV_RADSQ = 1.0 / RAD**2

       DO IY = 1,NY
           YRADSQ = (IY-YCEN)**2 * RINV_RADSQ

           DO IX = 1,NX
	      RADNOW = 0.5 * ((IX-XCEN)**2 * RINV_RADSQ + YRADSQ )

	      IF (RADNOW >= TNM)  THEN
	         CMASK(IX,IY) = 0.0
	      ELSE
	         VALG         = 0.5 * (2*RADNOW)**NORDER
                 CMASK(IX,IY) = EXP(-VALG)
	      ENDIF
           ENDDO
       ENDDO

       END



#ifdef NEVER

      SUBROUTINE ARITHSCA_NOLUN(BUFIN,BUFOUT,NX,NY,NZ,
     &                          FMINT,FMAXT,FLOW,FHI)

      IMPLICIT NONE

      REAL    :: BUFIN(NX*NY*NZ)
      REAL    :: BUFOUT(NX*NY*NZ)
      INTEGER :: NX,NY,NZ 
      REAL    :: FMINT,FMAXT,FLOW,FHI 

      REAL    :: RANGENEW,RANGEOLD,CON,CON2

      IF (FMAXT == FMINT) THEN
         BUFOUT = BUFIN
         RETURN
      ENDIF
  
      RANGENEW  = FHI   - FLOW
      RANGEOLD  = FMAXT - FMINT
      CON       = RANGENEW / RANGEOLD
      CON2      = FLOW - CON * FMINT

      BUFOUT    = CON2 + CON * BUFIN

      END
#endif
