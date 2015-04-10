
C **********************************************************************
C *  ORACFMSK.F 
C                  OPFILEC                         FEB 03 ARDEAN LEITH
C                  FMRS_PLAN                       MAY 08 ARDEAN LEITH
C                  CCRS_2I                         APR 09 ARDEAN LEITH
C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2002 & 2009, P. A. Penczek                            *
C=*                                                                    *
C=* University of Texas - Houston Medical School                       *
C=* Email:  pawel.a.penczek@uth.tmc.edu                                *
C=*                                                                    *
C=* This program is free software; you can redistribute it and/or      *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* This program is distributed in the hope that it will be useful,    *
C=* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C=* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
C=* General Public License for more details.                           *
C=*                                                                    *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program; if not, write to the                      *
C=* Free Software Foundation, Inc.,                                    *
C=* 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.      *
C=*                                                                    *
C=**********************************************************************
C
C ORACFMSK
C
C PURPOSE: ORIENTATION SEARCH - 2D
C          DETERMINES ROTATIONAL AND TRANSLATIONAL ORIENTATION BETWEEN 
C          TWO IMAGES USING SELF-CORRELATION FUNCTIONS. OPTIONALLY 
C          WITH ADDITIONAL CHECK OF MIRROR TRANSFORMATION. 
C    
C NOTE:    SLOPPILY WRITTEN.  USES OBSOLETE DUPLICATED ROUTINES.
C          APPEARS TO BE CLONED FROM : ORACFMSKM.F
c
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE ORACFMSK

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=MAXNAM):: FINPIC,FINPAT

        REAL, ALLOCATABLE    :: REF(:,:), ALIGNED(:,:)
        INTEGER, ALLOCATABLE :: NUMR(:,:)
        REAL, ALLOCATABLE    :: REFER_CIRC(:)
        REAL, ALLOCATABLE    :: REFPAD(:,:)

        LOGICAL              :: MIRROR
        CHARACTER(LEN=1)     :: MODE,ASK
        CHARACTER(LEN=3)     :: CMIR

	INTEGER, PARAMETER   :: INPIC   = 77
	INTEGER, PARAMETER   :: INREF   = 76

C       ASK FOR DATA FILES
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FINPIC,INPIC,'O',ITYPE,NX,NY,
     &              NZ,MAXIM,'EXPERMENTAL IMAGE',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FINPIC,INREF,'O',ITYPE,NX,NY,
     &              NZ,MAXIM,'REFERENCE IMAGE',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CLOSE (INPIC)
           RETURN
        ENDIF

        CALL RDPRI2S(MRA,MRR,NOT_USED,
     &    'REAL SPACE MASK RADIUS FOR EXP. AND REF. IMAGES',IRTFLG)
        IF (IRTFLG.NE.0)  GOTO 9999

        IF (MRR <= 0 .OR. 
     &      MRR >= MIN((NX/2),(NY/2))) THEN
           CALL ERRT(102,'INVALID EXP. MASK RADIUS',MRR)
           GOTO 9999
        ENDIF

        IF (MRA <= 0 .OR. 
     &      MRA >= MIN((NX/2),(NY/2))) THEN
           CALL ERRT(102,'INVALID REF MASK RADIUS',MRA)
           GOTO 9999
        ENDIF

        CALL RDPRI2S(NRING,NSHIFT,NOT_USED,
     &    'RADIUS OF THE ACF, MAXIMUM SHIFT',IRTFLG)
        IF (IRTFLG.NE.0)  GOTO 9999

        IF (NRING <= 0. OR.
     &      NRING >= MIN((NX/2),(NY/2)) .OR.
     &      NRING > 2*MIN(MRA,MRR) ) THEN
           CALL ERRT(102,'INVALID ACF RADIUS',MRA)
           GOTO 9999
        ENDIF

        CALL RDPRMC(ASK,NA,.TRUE.,
     &          'CHECK MIRRORED POSITIONS? (Y/N)', NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        MIRROR = (ASK == '1' .OR. ASK == 'Y')

        MODE = 'H'
C       PADDED DIMENSIONS FOR CCF AND SACF ARE ALWAYS EVEN

C       CCF
        I    = MAX(NX,MRA+MRR+NSHIFT)
        NEWS = NEARESTFFTFRIEND(I)
        I    = MAX(NY,MRA+MRR+NSHIFT)
        NEWR = NEARESTFFTFRIEND(I)
        WRITE(NOUT,*)  ' Dimensions for CCF:',NEWS,NEWR

C       ACF
        I     = MAX(4*MAX(MRA,MRR)+2,NX)
        NAS   = NEARESTFFTFRIEND(I)
        I     = MAX(4*MAX(MRA,MRR)+2,NY)
        NAR   = NEARESTFFTFRIEND(I)
        WRITE(NOUT,*) ' Dimensions for ACF:',NAS,NAR

        ALLOCATE(NUMR(3,NRING),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MWANT = 3 * NRING
           CALL ERRT(46,'ORACFMSK;  NUMR',MWANT)
           GOTO 9999
        ENDIF

        DO I=1,NRING
           NUMR(1,I) = I
        ENDDO

C       PUTS CIRCULAR RINGS IN A LINEAR ARRAY, CONCATENATED TOGETHER.
C       RETURNS NUMR & LCIRC.
        CALL ALPRBS(NUMR,NRING,LCIRC,MODE)
        MAXRIN = NUMR(3,NRING)

        ALLOCATE(REF(NX,NY),
     &           ALIGNED(NX,NY),
     &           REFER_CIRC(LCIRC),
     &           REFPAD(NEWS+2,NEWR),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            MWANT = NX*NY*2 + LCIRC + (NEWS+2)*NEWR
            CALL ERRT(46,'ORACFMSK; REFER_CIRC',MWANT)
            GOTO 9999
        ENDIF 

C       INPUT SAMPLE AND REFERENCE IMAGES
        NSL = 1
        CALL READV(INPIC,ALIGNED,NX,NY,NX,NY,NSL)
        CALL READV(INREF,REF,    NX,NY,NX,NY,NSL)

C       DETERMINE NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)
C       WRITE(NOUT,*) ' NUMBER OF OMP THREADS: ',NUMTH

        CALL PREPREF1(REF,NX,NY,NEWS,NEWR,NAS,NAR,MRR,
     &                LCIRC,NUMR,NRING,MAXRIN,NUMTH,
     &                REFER_CIRC,REFPAD,IRTFLG)

        IF (ALLOCATED(REF)) DEALLOCATE (REF)
        IF (IRTFLG .NE. 0) GOTO 9999

        CALL ORACFMSK_1(ALIGNED,NX,NY,
     &          NEWS,NEWR,NAS,NAR,MRA,MRR,NSHIFT,MIRROR,
     &          LCIRC,NUMR,NRING,MAXRIN,NUMTH,REFER_CIRC,REFPAD,
     &          ALPHA,SX,SY,MR,CMA)

        CMIR = 'No '
        IF (MR > 0) CMIR = 'Yes' 
        WRITE(NOUT,*)  ' '
        WRITE(NOUT,90) ALPHA,SX,SY,CMIR,CMA
90      FORMAT('  Rotation:',F7.2,'  Shifts: (',F7.3,',',F7.3,
     &         ')   Mirror: ',A,'  CC:',F8.4)
        IF (VERBOSE) WRITE(NOUT,*)  ' '
 
        CALL  REG_SET_NSEL(1,5,ALPHA,SX,SY,REAL(MR),CMA,IRTFLG)

9999    IF(ALLOCATED(REFER_CIRC)) DEALLOCATE(REFER_CIRC)
        IF(ALLOCATED(REFPAD))     DEALLOCATE(REFPAD)
        IF(ALLOCATED(ALIGNED))    DEALLOCATE(ALIGNED)
        IF(ALLOCATED(NUMR))       DEALLOCATE(NUMR)

        CALL FMRS_DEPLAN(IRTFLG)

        CLOSE(INPIC)
        CLOSE(INREF)

        END

C       ------------  PREPREF1 ----------------------------------

        SUBROUTINE PREPREF1(REF,NX,NY,NEWS,NEWR,NAS,NAR,MRR,
     &                LCIRC,NUMR,NRING,MAXRIN,NUMTH,
     &                REFER_CIRC,REFPAD,IRTFLG)

C       EXTERNAL ARRAYS
        INTEGER :: NUMR(3,NRING)
        REAL    :: REF(NX,NY),REFER_CIRC(LCIRC),REFPAD(NEWS+2,NEWR)

        REAL, ALLOCATABLE :: BUF(:,:)

C       AUTOMATIC ARRAYS
        REAL              :: WR(NRING)

        CHARACTER(LEN=1)  :: MODE = 'H'

        ALLOCATE(BUF(NAS+2,NAR),  STAT=IRTFLG)   !FFTW PADDING
        IF (IRTFLG .NE. 0) THEN
            MWANT = (NAS+2)*NAR
            CALL  ERRT(46,'PREPREF1; BUF',MWANT)
            RETURN
        ENDIF 

        CALL FMRS_PLAN(.TRUE.,BUF, NAS,NAR,1, 0,+1,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        CALL FMRS_PLAN(.TRUE.,BUF, NAS,NAR,1, 0,-1,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        CALL FMRS_PLAN(.TRUE.,BUF, NAS,NAR,1, 1,+1,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        CALL FMRS_PLAN(.TRUE.,BUF, NAS,NAR,1, 1,-1,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       FINDS WEIGHTS FOR RADIAL X-CORRELATION RINGS                                                            *
        CALL RINGWE(WR,NUMR,NRING,MAXRIN)

        BUF(1:NX,1:NY) = REF    ! WHOLE ARRAY COPY, PADDED

C       HAS FORWARD & REVERSE FFT ON: NAS,NAR
        CALL MACF_PL(NX,NY,BUF,NAS,NAR,MRR)

        CALL ALRQ_LSD(BUF,NAS+2,NAS,NAR,NUMR,
     &                REFER_CIRC,LCIRC,NRING,MODE)

C       HAS SPIDER FFT
        CALL FRNGS(REFER_CIRC,LCIRC,NUMR,NRING)

        CALL APPLYWS(REFER_CIRC,LCIRC,NUMR,WR,NRING,MAXRIN)

C       HAS: FORWARD FMRS2 on: NEWS,NEWR 
        CALL PREPFORCCN(REF,REFPAD,NX,NY,NEWS,NEWR,MRR)

9999    IF (ALLOCATED(BUF)) DEALLOCATE(BUF)

        END


C       ************************ ORACFMSK_1 ******************

        SUBROUTINE ORACFMSK_1(ALIGNED,NX,NY,
     &          NEWS,NEWR,NAS,NAR,MRA,MRR,NSHIFT,MIRROR,
     &          LCIRC,NUMR,NRING,MAXRIN,NUMTH,REFER_CIRC,REFPAD,
     &          ALPHA,SX,SY,MR,CMA)

        INTEGER              :: NUMR(3,NRING)
        REAL                 :: ALIGNED(NX,NY)
        REAL                 :: REFER_CIRC(LCIRC)
        REAL                 :: REFPAD(NEWS+2,NEWR)

        REAL, ALLOCATABLE    :: CIRA(:,:),CIRR(:,:),DIVIS(:,:)
        DOUBLE PRECISION     :: TT(1)

        LOGICAL              :: MIRROR

C       PREPARE THE NORMALIZATION FILE FOR THE CCC
        ICPX = INT(NEWS/2)+1
        ICPY = INT(NEWR/2)+1
        MRR2 = MRR**2
        MRA2 = MRA**2

        ALLOCATE(CIRA(NEWS+2,NEWR),
     &           CIRR(NEWS+2,NEWR),
     &           DIVIS(-NSHIFT:NSHIFT,-NSHIFT:NSHIFT),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            MWANT = 2*(NEWS+2)*NEWR + 2*(NSHIFT+1)*2*(NSHIFT+1)
            CALL  ERRT(46,'ORACFMSK_1; CIRA...',MWANT)
            GOTO 9999
        ENDIF

        CALL FMRS_PLAN(.TRUE.,CIRA, NEWS,NEWR,1, 0,+1,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        CALL FMRS_PLAN(.TRUE.,CIRA, NEWS,NEWR,1, 0,-1,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        CALL FMRS_PLAN(.TRUE.,CIRA, NEWS,NEWR,1, 1,+1,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        CALL FMRS_PLAN(.TRUE.,CIRA, NEWS,NEWR,1, 1,-1,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

c$omp   parallel do private (i,j,j2,i2)
        DO J=1,NEWS
          J2 = (J-ICPY)**2
          DO I=1,NEWR
            I2 = J2+(I-ICPX)**2
            IF (I2 > MRA2) THEN
               CIRA(I,J) = 0.0
            ELSE
               CIRA(I,J) = 1.0
            ENDIF

            IF (I2 > MRR2) THEN
                CIRR(I,J) = 0.0
            ELSE
                CIRR(I,J) = 1.0
            ENDIF
          ENDDO
        ENDDO

C       FORWARD FOURIER TRANSFORMS OF: CIRA & CIRR ON: NEWS,NEWR
        INV = +1
        CALL FMRS_2(CIRA,NEWS,NEWR,INV)
        INV = +1
        CALL FMRS_2(CIRR,NEWS,NEWR,INV)

C       CROSS-CORRELATE: CIRA x CIRR, REVERSE FFT: CIRA, RETURNS: CIRA
        LSC = NEWS+2-MOD(NEWS,2)

        CALL CCRS_2(CIRA,CIRR,CIRA, LSC,NEWS,NEWR)
        !!CALL CCRS_2I(CIRA,CIRR, LSC,NEWS,NEWR)

C       divis is only use of cross-correlation step?? al

        TMP   = 1.0 / CIRA(ICPX,ICPY)
        DIVIS = CIRA(ICPX-NSHIFT:ICPX+NSHIFT, 
     &               ICPY-NSHIFT:ICPY+NSHIFT) * TMP

        NIMA = 1
        IDI  = 1

C       OUTPUT PARAMETERS ARE:
C       NUMBER OF THE MOST SIMILAR REFERENCE PROJECTION.
C       CORR COEFF.(D5), ANGLE (D4), SHIFT: XSHSUM, YSHSUM

C       HAS FORWARD & REVERSE FFT ON: NAS,NAR, & ON: NEWS,NEWR
        CALL ORA2D(ALIGNED,REFER_CIRC,REFPAD,TT,NUMR,
     &           NX,NY,NEWS,NEWR, NAS,NAR,MRA,DIVIS,NSHIFT,
     &           LCIRC,NRING,MAXRIN,NIMA,MIRROR,
     &           IDI, CMA,ALPHA,SX,SY)

        IF (IDI > 0) THEN
           MR = 0
        ELSE
           MR = 1
        ENDIF

C       DEALLOCATE LOCAL ARRAYS
9999    IF (ALLOCATED(CIRA)) DEALLOCATE(CIRA)
        IF (ALLOCATED(CIRR)) DEALLOCATE(CIRR)

        END
