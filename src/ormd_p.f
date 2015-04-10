C++*********************************************************************
C
C ORMD_P.F
C             ANGLE BUG                         JAN 2013 ARDEAN LEITH                      
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2013  Health Research Inc.,                         *
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
C  ORMD_P(NX,NY, NRING,LCIRC,MAXRAYS,NUMR,X,NPEAK,MODE,LUNEXP,LUNREF)
C
C  PARAMETERS:   NX,NY                  SIZE                    INPUT
C                LUNEXP,LUNREF,LUNDOC   IO UNITS                INPUT
C
C  NOTE:  THE SUB-PIXEL INTERPOLATION IS DONE DIFFERENTLY FROM
C         'OR SH' THAT IS REASON FOR SMALL VARIANCE IN POSITION
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE ORMD_P(NX,NY, 
     &                    NRING,LCIRC,MAXRAYS,NUMR,X,
     &                    NPEAK,MODE,FFTW_PLANS,
     &                    LUNEXP,LUNREF,LUNDOC)
 
        INCLUDE 'CMBLOCK.INC'

        INTEGER                  :: NX,NY
        INTEGER                  :: NRING,LCIRC,MAXRAYS
        INTEGER                  :: NUMR(3,NRING)
        REAL                     :: X(NX,NY)
        INTEGER                  :: NPEAK
	CHARACTER(LEN=1)         :: MODE 
        INTEGER                  :: LUNEXP,LUNREF,LUNDOC

        REAL, ALLOCATABLE        :: CIRCEXP(:),CIRCREF(:) 
        REAL                     :: Q(MAXRAYS+2)
        REAL                     :: PEAK(2,NPEAK)
        DOUBLE PRECISION         :: T7(-3:3)
        REAL                     :: WRE(NRING)
        REAL                     :: WRR(NRING)
        DOUBLE PRECISION         :: DQMAX
        REAL                     :: QMAX
        REAL                     :: POS_MAX
        INTEGER                  :: MAXL

        INTEGER, PARAMETER       :: NDLI       = 2
        REAL                     :: DLIST(NDLI)
        LOGICAL, PARAMETER       :: USE_OMP    = .FALSE.
        LOGICAL, PARAMETER       :: USE_MIR    = .FALSE.

        ALLOCATE (CIRCEXP(LCIRC), CIRCREF(LCIRC), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'ORMD_P, CIRC',2*LCIRC)
           RETURN
        ENDIF

C       CALCULATE DIMENSIONS FOR NORMALIZING IN APRINGS_ONE
	CNS2 = NX/2+1     ! IMAGE CENTER
	CNR2 = NY/2+1     ! IMAGE CENTER

C       LOAD REF. IMAGE DATA ------------------------------
        CALL REDVOL(LUNREF,NX,NY,1,1,X,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        CLOSE(LUNREF)

C       EXTRACT REF. IMAGE RINGS, NORMALIZE & FFT THEM
        WRR(1) = 0.0     ! ONLY WEIGHT EXP IMAGES
        CALL APRINGS_ONE_NEW(NX,NY, CNS2,CNR2, X, USE_OMP,
     &                       MODE,NUMR,NRING,LCIRC, WRR,FFTW_PLANS,
     &                       CIRCREF,IRTFLG)

C       LOAD EXP. IMAGE DATA ----------------------------
        CALL REDVOL(LUNEXP,NX,NY,1,1,X,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        CLOSE(LUNEXP)

C       FIND  WEIGHTS FOR TRANSFORMED CIRC RINGS 
        CALL RINGWE_NEW(WRE,NUMR,NRING,MAXRAYS)

C       EXTRACT EXP. IMAGE RINGS, NORMALIZE & FFT THEM
        CALL APRINGS_ONE_NEW(NX,NY, CNS2,CNR2, X, USE_OMP,
     &                       MODE,NUMR,NRING,LCIRC, WRE,FFTW_PLANS,
     &                       CIRCEXP,IRTFLG)
    
        CALL CROSRNG_NEW(CIRCREF,CIRCEXP,LCIRC, 
     &                NRING,MAXRAYS,NUMR,
     &                FFTW_PLANS, USE_MIR,
     &                Q,  QMAX,POS_MAX,MAXL)

        IF (NPEAK <= 1)  THEN

C          NORMALIZE PEAK 
           QMAX = 2.0 * Q(MAXL) / MAXRAYS/ MAXRAYS / NRING 

C          CONVERT PEAK LOCATION TO ANGLE
           RANGNEW = ANG_N(POS_MAX,MODE,MAXRAYS)

           WRITE(NOUT,2799)RANGNEW,QMAX    
2799       FORMAT('  Angle: ',F10.4,'    Peak Height: ',G12.5)
          
C          ONE PEAK, FIRST REGISTER IS ANGLE, SECOND IS PEAK VALUE
           CALL REG_SET_NSEL(1,2,RANGNEW,QMAX, 0.0,0.0,0.0,IRTFLG)

        ELSE
C          FIND SPECIFIED NUMBER OF PEAKS
           DO K2=1,NPEAK
              PEAK(1,K2) = -HUGE(PEAK)
              PEAK(2,K2) = -1.0
           ENDDO

           DO J=1,MAXRAYS
              J2 = MOD(J+MAXRAYS-1,MAXRAYS)+1
              J1 = MOD(J+MAXRAYS-2,MAXRAYS)+1
              J3 = MOD(J+MAXRAYS,  MAXRAYS)+1

              IF (Q(J2) > Q(J1) .AND. 
     &            Q(J2) > Q(J3)) THEN
C                LOCAL PEAK HERE

                 DO K2=1,NPEAK
                    IF (Q(J2) > PEAK(1,K2))  THEN
                       IF (NPEAK > 1)  THEN
C                         MOVE THIS PEAK IN LIST
                          DO K3=NPEAK,K2+1,-1
                             PEAK(1,K3) = PEAK(1,K3-1)
                             PEAK(2,K3) = PEAK(2,K3-1)
                          ENDDO
                       ENDIF

                       PEAK(1,K2) = Q(J2)
                       PEAK(2,K2) = J
                       EXIT
                    ENDIF
                 ENDDO
              ENDIF
           ENDDO

           !write(6,*) (peak(2,i),i=1,npeak)

C          CONVERT TO ANGLES AND INTERPOLATE
           DO K2=1,NPEAK
              JTOT = PEAK(2,K2)
              NPT  = K2 
              IF (JTOT == -1) EXIT

              DO K3=-3,3
                 J      = MOD(JTOT+K3+MAXRAYS-1,MAXRAYS) + 1
                 T7(K3) = Q(J)
              ENDDO

C             SUB-PIXEL INTERPOLATION
              CALL PRB1D(T7,7,POS)

              QMAX = 2.0 * PEAK(1,K2) / MAXRAYS/ MAXRAYS / NRING ! HEIGHT
              RANG = PEAK(2,K2) + POS                    ! LOCATION

C             CONVERT PEAK LOCATION TO ANGLE
              RANGNEW = ANG_N(RANG,MODE,MAXRAYS)

              DLIST(1) = RANGNEW     ! ANGLE
              DLIST(2) = QMAX        ! HEIGHT
              CALL LUNDOCWRTDAT(LUNDOC,K2,DLIST,NDLI,IRTFLG)

              !write(6,90) j,pos,rang,peak(2,k2)
  90          format(   ' j,pos,rang,angle: ', i5, f8.2, f8.2, f8.2)

              WRITE(NOUT,2701) PEAK(2,K2),QMAX
2701          FORMAT('  Angle: ',F10.4,'    Peak Height: ',G12.5)
           ENDDO

C          MULTIPLE PEAKS, FIRST REGISTER IS THE NUMBER OF PEAKS SAVED
           CALL REG_SET_NSEL(1,1,FLOAT(NPT), 0.0,0.0,0.0,0.0,IRTFLG)
        ENDIF

9999    IF (ALLOCATED(CIRCEXP)) DEALLOCATE (CIRCEXP)
	IF (ALLOCATED(CIRCREF)) DEALLOCATE (CIRCREF)
       
        END

