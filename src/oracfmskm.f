C **********************************************************************
C *  ORACFMSKM.F 
C                  OPFILEC                         FEB 03 ARDEAN LEITH
C                  FMRS_PLAN                       MAY 08 ARDEAN LEITH
C                  COSMETIC, || BUGGY              MAY 08 ARDEAN LEITH
C                  CCRS_2I                         APR 09 ARDEAN LEITH
C                  REVERT FROM CCRS_2I BUG         JUN 13 ARDEAN LEITH
C                  REWRITE                         JUN 13 ARDEAN LEITH
C
C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2002 & 2013,  P. A. Penczek & ArDean Leith
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
C PURPOSE: 'OR Q' PERFORMS MULTIREFERENCE ALIGNMENT BETWEEN A SERIES OF 
C          IMAGES AND A SET OF REFERENCE IMAGES (TEMPLATES). THE 
C          MIRROR ORIENTATION IS NOT CHECKED. THE OPERATION USES 
C          SELF-CORRELATION FUNCTION. 
C
C NOTE:    SLOPPILY WRITTEN.  USES OBSOLETE & DUPLICATED ROUTINES.
C          LOT OF GARBAGE CODE. al
C          APPEARS TO BE CLONE OF: ORACFMSK.F
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE ORACFMSKM

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 

        LOGICAL               :: MIRROR
        CHARACTER(LEN=1)      :: MODE
        CHARACTER(LEN=1)      :: ASK

        CHARACTER(LEN=MAXNAM) :: FINPAT,FINPIC,FILTOA
        CHARACTER(LEN=MAXNAM) :: FILREF
        INTEGER, ALLOCATABLE  :: NUMR(:,:)
        INTEGER, ALLOCATABLE  :: ILIST(:),IRIST(:)
        REAL,    ALLOCATABLE  :: REFER_CIRC(:,:)
        REAL,    ALLOCATABLE  :: REFPAD(:,:,:)

	INTEGER, PARAMETER    :: LUNREF  = 78
	INTEGER, PARAMETER    :: INPIC   = 77

C       ALLOCATE SPACE FOR REFERENCE IMAGE FILE LIST
        NILMAX = NIMAX
        ALLOCATE(ILIST(NILMAX),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'ORACFMSKM; ILIST',NILMAX)
           GOTO 9999
        ENDIF

C       ASK FOR REFERENCE IMAGE FILE LIST
        CALL FILELIST(.TRUE.,INPIC,FINPAT,NLET,ILIST,NILMAX,NIMA,
     &      'TEMPLATE FOR REFERENCE IMAGE SERIES',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (NIMA < 1)  THEN
           CALL ERRT(101,' No reference images',IDUM)
           GOTO 9999
        ENDIF

        WRITE(NOUT,2001) NIMA
2001    FORMAT('  Number of reference images: ',I6)


C       GET FIRST REFERENCE IMAGE TO DETERMINE DIMS
        CALL FILGET(FINPAT,FINPIC,NLET,ILIST(1),INTFLG)
        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,
     &               NX,NY,NZ,
     &               MAXIM,' ',.FALSE.,IRTFLG)
        IF (IRTFLG.NE.0)  GOTO 9999
        CLOSE(INPIC)

        CALL RDPRI2S(MRA,MRR,NOT_USED,
     &    'REAL SPACE MASK RADIUS FOR EXP. AND REF. IMAGES',IRTFLG)
        IF (IRTFLG.NE.0)  GOTO 9999

        IF (MRR .LE. 0 .OR. 
     &      MRR .GE. MIN((NX/2),(NY/2))) THEN
           CALL ERRT(102,'INVALID EXP. MASK RADIUS',MRR)
           GOTO 9999
        ENDIF

        IF (MRA .LE. 0 .OR. 
     &      MRA .GE. MIN((NX/2),(NY/2))) THEN
           CALL ERRT(102,'INVALID REF MASK RADIUS',MRA)
           GOTO 9999
        ENDIF

        CALL RDPRI2S(NRING,NSHIFT,NOT_USED,
     &    'RADIUS OF THE ACF, MAXIMUM SHIFT',IRTFLG)
        IF (IRTFLG.NE.0)  GOTO 9999

        IF (NRING .LE. 0. OR.
     &      NRING .GE. MIN((NX/2),(NY/2)) .OR.
     &      NRING > 2*MIN(MRA,MRR) ) THEN
           CALL ERRT(102,'INVALID ACF RADIUS',MRA)
           GOTO 9999
        ENDIF

        CALL RDPRMC(ASK,NA,.TRUE.,
     &          'CHECK MIRRORED POSITIONS? (Y/N)', NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        MIRROR = (ASK == '1' .OR. ASK == 'Y')

        MODE = 'H'
C       PADDED DIMENSIONS FOR CCF AND SACF ARE ALWAYS EVEN CCF
        I    = MAX(NX,MRA+MRR+NSHIFT)
        NEWS = NEARESTFFTFRIEND(I)
        I    = MAX(NY,MRA+MRR+NSHIFT)
        NEWR = NEARESTFFTFRIEND(I)
        WRITE(NOUT,*)  ' Dimensions for CCF:',NEWS,NEWR

C       ACF
        I   = MAX(4*MAX(MRA,MRR)+2,NX)
        NAS = NEARESTFFTFRIEND(I)
        I   = MAX(4*MAX(MRA,MRR)+2,NY)
        NAR = NEARESTFFTFRIEND(I)
        WRITE(NOUT,*)  ' Dimensions for ACF:',NAS,NAR

        ALLOCATE(NUMR(3,NRING),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'ORACFMSKM;  NUMR',3*NRING)
           GOTO 9999
        ENDIF

        DO I=1,NRING
           NUMR(1,I) = I
        ENDDO

C       PUTS CIRCULAR RINGS IN A LINEAR ARRAY, CONCATENATED TOGETHER.
C       RETURNS NUMR & LCIRC. PREPARES FOR SPIDER FFT NOT FFTW3
        CALL ALPRBS(NUMR,NRING,LCIRC,MODE)
        MAXRIN = NUMR(3,NRING)

C       ALLOCATE SPACE FOR REFERENCE CIRCLES ARRAY
        ALLOCATE(REFER_CIRC(LCIRC,NIMA),
     &          REFPAD(NEWS+2,NEWR,NIMA),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'ORACFMSKM; REFER_CIRC',IER)
            GOTO 9999
        ENDIF 

C       DETERMINE NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)
C       WRITE(NOUT,*) ' NUMBER OF OMP THREADS: ',NUMTH

        CALL PREPREF(FINPAT,NLET,ILIST,NIMA,
     &               NX,NY,NEWS,NEWR,NAS,NAR,MRR,
     &               LCIRC,NUMR,NRING,MAXRIN,NUMTH,
     &               REFER_CIRC,REFPAD,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
 
C       ALLOCATE SPACE FOR ALIGNED IMAGES FILE LIST
        ALLOCATE(IRIST(NILMAX),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'ORACFMSKM; NUMR',IER)
           GOTO 9999
        ENDIF

C       GET LIST OF SAMPLE IMAGES TO BE ALIGNED
        CALL FILELIST(.TRUE.,INPIC,FILTOA,NLETI,IRIST,NILMAX,NTOTAL,
     &     'TEMPLATE FOR IMAGE SERIES TO BE ALIGNED',IRTFLG)
        IF (IRTFLG.NE.0) GOTO 9999

        IF (NTOTAL < 1)  THEN
           CALL ERRT(101,'  No experimental images!',IDUM)
           GOTO 9999
        ENDIF

        WRITE(NOUT,2002) NTOTAL
2002    FORMAT('  Number of experimental images: ',I6,/)

C       NIMA   IS NUMBER OF REFERENCE IMAGES
C       NTOTAL IS NUMBER OF SAMPLE IMAGES

        CALL ORACFMSK_PS(FILTOA,ILIST,NIMA,IRIST,NTOTAL,NX,NY,
     &          NEWS,NEWR,NAS,NAR,MRA,MRR,NSHIFT,MIRROR,
     &          LCIRC,NUMR,NRING,MAXRIN,NUMTH,REFER_CIRC,REFPAD)

9999    IF (ALLOCATED(IRIST))      DEALLOCATE(IRIST)
        IF (ALLOCATED(ILIST))      DEALLOCATE(ILIST)
        IF (ALLOCATED(NUMR))       DEALLOCATE(NUMR)
        IF (ALLOCATED(REFER_CIRC)) DEALLOCATE(REFER_CIRC)
        IF (ALLOCATED(REFPAD))     DEALLOCATE(REFPAD)

        END


C       --------------- PREPREF -------------------------------------

        SUBROUTINE PREPREF(FINPAT,NLET,ILIST,NIMA,
     &                NX,NY,NEWS,NEWR,NAS,NAR,MRR,
     &                LCIRC,NUMR,NRING,MAXRIN,NUMTH,
     &                REFER_CIRC,REFPAD,IRTFLG)

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

C       EXTERNAL ARRAYS
        INTEGER               :: ILIST(NIMA),NUMR(3,NRING)
        REAL                  :: REFER_CIRC(LCIRC,NIMA)
        REAL                  :: REFPAD(NEWS+2,NEWR,NIMA)
        CHARACTER(LEN=MAXNAM) :: FINPAT,FINPIC

C       AUTOMATIC ARRAYS
        REAL                  :: WR(NRING)

        REAL, ALLOCATABLE     :: A(:,:,:)
        REAL, ALLOCATABLE     :: BUF(:,:,:)

        CHARACTER(LEN=1)      :: MODE

	INTEGER, PARAMETER    :: INPIC   = 58

        IRTFLG = 0
        MODE   = 'H'

        ALLOCATE(A(NX,NY,NUMTH),
     &           BUF(NAS+2,NAR,NUMTH), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            MWANT = NX*NY*NUMTH + (NAS+2)*NAR
            CALL ERRT(46,'PREPREF; A..',MWANT)
            GOTO 9999
        ENDIF 

        CALL RINGWE(WR,NUMR,NRING,MAXRIN)

        CALL FMRS_PLAN(.TRUE.,BUF, NAS,NAR,1, 1,+1,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        CALL FMRS_PLAN(.TRUE.,BUF, NAS,NAR,1, 1,-1,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        CALL FMRS_PLAN(.TRUE.,BUF, NEWS,NEWR,1, 1,+1,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        CALL FMRS_PLAN(.TRUE.,BUF, NEWS,NEWR,1, 1,-1,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       PREPARE CIRCULAR RINGS DATA FOR ALL REFERENCE IMAGES
        DO  IMIT=1,NIMA,NUMTH
           DO  IMI=IMIT,MIN(NIMA,IMIT+NUMTH-1)
              CALL FILGET(FINPAT,FINPIC,NLET,ILIST(IMI),INTFLAG)
              MAXIM = 0
              CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,NXT,NYT,
     &                NZ, MAXIM,' ',.FALSE.,IRTFLG)
              IF (IRTFLG .NE. 0)  GOTO 9999

              CALL READV(INPIC,A(1,1,IMI-IMIT+1),NX,NY,NX,NY,1)
              CLOSE(INPIC)
           ENDDO

c$omp      parallel do private(IMI)
           DO  IMI=IMIT,MIN(NIMA,IMIT+NUMTH-1)
             BUF(1:NX,1:NY,IMI-IMIT+1) = A(:,:,IMI-IMIT+1)

C            HAS: FORWARD FMRS2 on NAS,NAR
             CALL MACF_PL(NX,NY,BUF(1,1,IMI-IMIT+1),NAS,NAR,MRR)

C            INTERPOLATION INTO POLAR COORDINATES
             CALL ALRQ_LSD(BUF(1,1,IMI-IMIT+1),NAS+2,NAS,NAR,NUMR,
     &                  REFER_CIRC(1,IMI),LCIRC,NRING,MODE)

C            HAS SPIDER FFT
             CALL FRNGS(REFER_CIRC(1,IMI),LCIRC,NUMR,NRING)

             CALL APPLYWS(REFER_CIRC(1,IMI),LCIRC,NUMR,WR,NRING,MAXRIN)

C            HAS: FORWARD FMRS2 on: NEWS,NEWR 
             CALL PREPFORCCN(A(1,1,IMI-IMIT+1),
     &                  REFPAD(1,1,IMI),NX,NY,NEWS,NEWR,MRR)
          ENDDO
        ENDDO

C       DEALLOCATE LOCAL ARRAYS
9999    IF (ALLOCATED(A))          DEALLOCATE(A)
        IF (ALLOCATED(BUF))        DEALLOCATE(BUF)

        END

C       ************************ PREPFORCCN *************************

        SUBROUTINE PREPFORCCN(REF,REFPAD,NX,NY,NEWS,NEWR,MRR)

        DIMENSION  REF(NX,NY),REFPAD(NEWS+2,NEWR)
        INTEGER    X31,X32,X33,X34,X35,X36,X37

C       PREPARE REFERENCE IMAGES FOR THE CCN

C       CENTER OF ORIGINAL IMAGES
        X33 = INT(NX/2)+1
        X34 = INT(NY/2)+1


        MRR2 = MRR**2
        AVE1 = 0.0
        ILE1 = 0
        DO J=1,NY
         J2 = (J-X34)**2
          DO I=1,NX
           I2 = J2+(I-X33)**2
           IF (I2 .LE. MRR2)  THEN
              AVE1 = AVE1 + REF(I,J)
              ILE1 = ILE1 + 1
           ENDIF
          ENDDO
        ENDDO
        AVE1 = AVE1 / ILE1

        DO J=1,NY
         J2 = (J-X34)**2
          DO I=1,NX
           I2 = J2+(I-X33)**2
           IF (I2 .LE. MRR2)  THEN
              REF(I,J) = REF(I,J)-AVE1
           ELSE
              REF(I,J) = 0.0
           ENDIF
          ENDDO
        ENDDO

C       FIRST PIXEL TO PAD
        X31    = INT(NEWS/2)-INT(NX/2)+1
        X32    = INT(NEWR/2)-INT(NY/2)+1
        REFPAD = 0.0
        REFPAD(X31:X31+NX-1,X32:X32+NY-1) = REF

C       FORWARD FFT ON: NEWS,NEWR
        INV = +1
        CALL FMRS_2(REFPAD,NEWS,NEWR,INV)

        END
        


C      ******************* ORACFMSK_PS ***********************

        SUBROUTINE ORACFMSK_PS(FILTOA,ILIST,NIMA,IRIST,NTOTAL,NX,NY,
     &          NEWS,NEWR,NAS,NAR,MRA,MRR,NSHIFT,MIRROR,
     &          LCIRC,NUMR,NRING,MAXRIN,NUMTH,REFER_CIRC,REFPAD)

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        CHARACTER(LEN=MAXNAM)         :: FINPIC,FILTOA,DOCNAM
        CHARACTER(LEN=88)             :: COMMENT

        INTEGER                       :: ILIST(NIMA)
        INTEGER                       :: IRIST(NTOTAL)
        INTEGER                       :: NUMR(3,NRING)
        REAL                          :: REFER_CIRC(LCIRC,NIMA)
        REAL                          :: REFPAD(NEWS+2,NEWR,NIMA)

        REAL                          :: DLIST(6,NUMTH)
        INTEGER                       :: NASSIG(NUMTH)
        LOGICAL                       :: MIRROR

        REAL, ALLOCATABLE             :: A(:,:,:)
        REAL, ALLOCATABLE             :: CIRA(:,:),CIRR(:,:),DIVIS(:,:)
        DOUBLE PRECISION              :: TT(1)
        
        LOGICAL                       :: ADDEXT,GETNAME,ISOLD,APPEND
        LOGICAL                       :: NEWFILE,MESSAGE
        INTEGER                       :: NLET,IRTFLG,NDOCO

	INTEGER, PARAMETER            :: NDOC    = 81
	INTEGER, PARAMETER            :: INPIC   = 21

        ADDEXT  = .TRUE.
        GETNAME = .TRUE.
        ISOLD   = .FALSE.
        APPEND  = .FALSE.
        MESSAGE = .TRUE.
        IRTFLG  = -8         ! NO IC USE

        CALL OPENDOC(DOCNAM,ADDEXT,NLET,NDOC,NDOCO,GETNAME,
     &           'ALIGNMENT DOC',ISOLD,APPEND,MESSAGE,
     &            NEWFILE,IRTFLG)

C                  123456789 123456789 123456789 123456789 123456789 123456789 
        COMMENT = '        #-REF        NN-CC           ANGLE       '//
     &            ' X-SHIFT      Y-SHIFT       #-EXP-IMG'
        CALL LUNDOCPUTCOM(NDOCO,COMMENT(1:88),IRTFLG)

        ALLOCATE(A(NX,NY,NUMTH),  STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            MWANT = NX*NY*NUMTH
            CALL ERRT(46,'ORACFMSK_PS; A... ',MWANT)
            RETURN
        ENDIF

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
            CALL  ERRT(46,'ORACFMSK_PS; CIRA...',MWANT)
            GOTO 9999
        ENDIF

        CALL FMRS_PLAN(.TRUE.,CIRA, NEWS,NEWR,1, 0,+1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CALL FMRS_PLAN(.TRUE.,CIRA, NEWS,NEWR,1, 0,-1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CALL FMRS_PLAN(.TRUE.,CIRA, NEWS,NEWR,1, 1,+1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CALL FMRS_PLAN(.TRUE.,CIRA, NEWS,NEWR,1, 1,-1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

c$omp   parallel do private(i,j,j2,i2)
        DO J=1,NEWS
          J2 = (J-ICPY)**2
          DO I=1,NEWR
            I2 = J2+(I-ICPX)**2
            IF (I2 > MRA2)  THEN
                CIRA(I,J) = 0.0
            ELSE
                CIRA(I,J) = 1.0
            ENDIF
            IF (I2 > MRR2)  THEN
                CIRR(I,J) = 0.0
            ELSE
                CIRR(I,J) = 1.0
            ENDIF
          ENDDO
        ENDDO

C       FORWARD FOURIER TRANSFORMS OF: CIRA & CIRR OVER: NEWS,NEWR
        INV = +1
        CALL FMRS_2(CIRA,NEWS,NEWR,INV)
        INV = +1
        CALL FMRS_2(CIRR,NEWS,NEWR,INV)

C       REVERSE FOURIER TRANSFORMS ON: NEWS,NEWR
C       CROSS-CORRELATE: CIRA x CIRR, RETURNS: CIRA
        LSC = NEWS+2-MOD(NEWS,2)

        !CALL CCRS_2I(CIRA,CIRR, LSC,NEWS,NEWR)
        CALL CCRS_2(CIRA,CIRR, CIRA,LSC,NEWS,NEWR)

        TMP   = 1.0 / CIRA(ICPX,ICPY)
        DIVIS = CIRA(ICPX-NSHIFT:ICPX+NSHIFT,ICPY-NSHIFT:ICPY+NSHIFT)*
     &          TMP

C       DEALLOCATE LOCAL ARRAYS
        IF (ALLOCATED(CIRA)) DEALLOCATE(CIRA)
        IF (ALLOCATED(CIRR)) DEALLOCATE(CIRR)

        NLETI = lnblnkn(FILTOA)

C       LOOP OVER IMAGES TO BE ALIGNED
        DO IMIT=1,NTOTAL,NUMTH
           DO IMI=IMIT,MIN(NTOTAL,IMIT+NUMTH-1)

              CALL FILGET(FILTOA,FINPIC,NLETI,IRIST(IMI),INTFLAG)

              MAXIM = 0
              CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,NXT,NYT,
     &             NZ, MAXIM,' ',.FALSE.,IRTFLG)
              IF (IRTFLG.NE.0) GOTO 9999

              IF (NXT .NE. NX .OR. NYT .NE. NY)  THEN
                  CALL ERRT(101,'INCONSISTENT IMAGE DIMENSIONS',NE)
                  CLOSE(INPIC)
                  GOTO 9999
              ENDIF

              CALL READV(INPIC,A(1,1,IMI-IMIT+1),NX,NY,NX,NY,1)
              CLOSE(INPIC)
           ENDDO

C          NUMTH INPUT IMAGES READY TO BE ALIGNED
C          OUTPUT PARAMETERS ARE:
C          NUMBER OF THE MOST SIMILAR REFERENCE PROJECTION.
C          CORR COEFF.(D5), ANGLE (D4), SHIFT: XSHSUM, YSHSUM

c$omp      parallel do private(imi)
           DO  IMI=IMIT,MIN(NTOTAL,IMIT+NUMTH-1)
C             COMPARE EACH IMAGE TO BE ALIGNED WITH ALL REFERENCE IMAGES

C             HAS FORWARD & REVERSE FFT ON: NAS,NAR, & ON: NEWS,NEWR
              CALL ORA2D(A(1,1,IMI-IMIT+1),REFER_CIRC,REFPAD,TT,NUMR,
     &           NX,NY,NEWS,NEWR,NAS,NAR,MRA,DIVIS,NSHIFT,
     &           LCIRC,NRING,MAXRIN,NIMA,MIRROR,
     &           NASSIG(IMI-IMIT+1),
     &           DLIST(2,IMI-IMIT+1),DLIST(3,IMI-IMIT+1),
     &           DLIST(4,IMI-IMIT+1),DLIST(5,IMI-IMIT+1))
           ENDDO

C          OUTPUT
C          1 - NUMBER OF THE MOST SIMILAR REFERENCE PROJECTION.
C          2 - NOT-NORMALIZED CORRELATION COEFFICIENT.
C          3 - PSI ANGLE. (IN=PLANE ROTATION)
C          4 - SX
C          5 - SY
C          6 - INPUT IMAGE NUMBER.

           DO IMI=IMIT,MIN(NTOTAL,IMIT+NUMTH-1)
              !write(6,*) 'imi,imit,nassig(imi-imit+1):', imi,imit,nassig(imi-imit+1)

              IT = IMI-IMIT+1 

              !DLIST(1,IT) = IMI
              DLIST(1,IT) = ISIGN(ILIST(IABS(NASSIG(IT))), NASSIG(IT))
              DLIST(6,IT) = IRIST(IMI)

              CALL LUNDOCWRTDAT(NDOCO,IMI,DLIST(1,IT),6,IRTFLG)

              !CALL SAVD(NDOC,DLIST(1,IT),7,IRTFLG)
           ENDDO

C---------------------------------------------------------------------
        ENDDO

        CLOSE(NDOC)

C       DEALLOCATE LOCAL ARRAYS
9999    IF (ALLOCATED(A)) DEALLOCATE(A)

        END




C      *************** ORA2D.F ************************************

        SUBROUTINE ORA2D(ALIGNED,REFER_CIRC,REFPAD,TT,NUMR,
     &       NX,NY,NEWS,NEWR,NAS,NAR,MRA,DIVIS,NSHIFT,
     &       LCIRC,NRING,MAXRIN,NIMA,MIRROR,
     &       IDI,CCC,PHI,SX,SY)

        REAL                  :: ALIGNED(NX,NY)
        REAL                  :: DIVIS(-NSHIFT:NSHIFT,-NSHIFT:NSHIFT)
        REAL                  :: REFER_CIRC(LCIRC,NIMA)
        INTEGER               :: NUMR(3,NRING)
        REAL                  :: REFPAD(NEWS+2,NEWR,NIMA)
        DOUBLE PRECISION      :: TT(*)
        LOGICAL               :: MIRROR
        CHARACTER*1           :: MODE = 'H'
        DOUBLE PRECISION      :: TOTMIN,TOTMIR
        REAL, ALLOCATABLE     :: A_CIRC(:)
        REAL, ALLOCATABLE     :: ALIGNEDPAD(:,:),ALPADM(:,:),ROTAP(:,:)
        REAL, ALLOCATABLE     :: BUF(:,:)
        INTEGER               :: X31,X32,X33,X34,X35,X36,X37

C       CENTER OF ORIGINAL IMAGES
        X33 = INT(NX/2)+1
        X34 = INT(NY/2)+1

C       FIRST PIXEL TO PAD
        X31 = INT(NEWS/2)-INT(NX/2)+1
        X32 = INT(NEWR/2)-INT(NY/2)+1

        ICPX = INT(NEWS/2)+1
        ICPY = INT(NEWR/2)+1

 !    &           FTMP(NEWS+2,NEWR),
        ALLOCATE(BUF(NAS+2,NAR),
     &           A_CIRC(LCIRC),
     &           ALIGNEDPAD(NEWS+2,NEWR),
     &           ROTAP(NEWS+2,NEWR), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MWANT = (NAS+2)*NAR + LCIRC + 2*(NEWS+2)*NEWR
           CALL ERRT(46,'ORA2D; A_CIRC ...',MWANT)
           RETURN
        ENDIF

        IF (MIRROR)  THEN
           ALLOCATE(ALPADM(NEWS+2,NEWR), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              MWANT = (NEWS+2)*NEWR
              CALL ERRT(46,'ORA2D; ALPADM',MWANT)
              GOTO 9999
           ENDIF
        ENDIF

        BUF(1:NX,1:NY) = ALIGNED   ! BUF = ALIGNED ARRAY COPY

C       HAS FORWARD & REVERSE FFT ON: NAS,NAR
        CALL MACF_PL(NX,NY,BUF,NAS,NAR,MRA)

C       INTERPOLATION INTO POLAR COORDINATES
        CALL ALRQ_LSD(BUF,NAS+2,NAS,NAR,NUMR,A_CIRC,LCIRC,NRING,MODE)

C       HAS SPIDER FFT
        CALL FRNG(A_CIRC,LCIRC,NUMR,NRING)

        IF (ALLOCATED(BUF)) DEALLOCATE(BUF)

C       MASK AND PAD ALIGNED IMAGE FOR THE CCN
        MRA2 = MRA**2
        AVE2 = 0.0
        ILE2 = 0

c$omp   parallel do private (i,j,i2,j2) reduction(+:ave2,ile2)
        DO J=1,NY
          J2 = (J-X34)**2
          DO I=1,NX
             I2 = J2+(I-X33)**2

             IF (I2 > MRA2)  THEN
                ALIGNED(I,J) = 0.0
             ELSE
                AVE2 = AVE2 + ALIGNED(I,J)
                ILE2 = ILE2 + 1
             ENDIF
          ENDDO
        ENDDO

        AVE2 = AVE2 / ILE2
 
c$omp   parallel do private (i,j,i2,j2)
        DO J=1,NY
          J2 = (J-X34)**2
          DO I=1,NX
             I2 = J2+(I-X33)**2
             IF (I2 .LE. MRA2)  ALIGNED(I,J) = ALIGNED(I,J) - AVE2
          ENDDO
        ENDDO

        ALIGNEDPAD = 0.0
        ALIGNEDPAD(X31:X31+NX-1,X32:X32+NY-1) = ALIGNED

        IF (MIRROR) THEN
C          MIRROR AROUND Y
           DO  J=1,NEWR
              ALPADM(1,J) = ALIGNEDPAD(1,J)
              DO  I=2,NEWS
                 ALPADM(NEWS-I+2,J) = ALIGNEDPAD(I,J)
              ENDDO
           ENDDO
        ENDIF

C       COMPARE ALIGNED IMAGE WITH ALL REFERENCE IMAGES
        CCC = -HUGE(AVE2)
        LSC = NEWS+2-MOD(NEWS,2)

        DO  IR=1,NIMA
C          CALCULATE THE ANGLE (IN DEGREES) BETWEEN THE ACFS

           IF (MIRROR) THEN
C             HAS SPIDER FFT
              CALL CROSRNG_MS(A_CIRC,REFER_CIRC(1,IR),LCIRC,NRING,
     &                       MAXRIN,NUMR,TOTMIN,TOT, TOTMIR,TMT, TT)

              ANGLEM = ANGMOR(TMT,MODE,MAXRIN) 
           ELSE
C             HAS SPIDER FFT
              CALL CROSRNG_DS(A_CIRC,REFER_CIRC(1,IR),LCIRC,NRING,
     &                        MAXRIN,NUMR,TOTMIN,TOT,TT)
           ENDIF

           ANGLE = ANGMOR(TOT,MODE,MAXRIN)

C          CHECK ANGLE POSITION
           CALL RTQ_Q(ALIGNEDPAD,ROTAP,NEWS+2,NEWS,NEWR,ANGLE)

C          HAS FORWARD FFT ON: NEWS,NEWR
           INV = +1
           CALL FMRS_2(ROTAP,NEWS,NEWR,INV)

C          HAS REVERSE FFT ON: NEWS,NEWR
          ! CALL CCRS_2I(REFPAD(1,1,IR),ROTAP, LSC,NEWS,NEWR)
           CALL CCRS_2(REFPAD(1,1,IR),ROTAP, ROTAP, LSC,NEWS,NEWR)

           !rotap = tmp

C          FIND THE PEAK
           ROTAP(ICPX-NSHIFT:ICPX+NSHIFT,ICPY-NSHIFT:ICPY+NSHIFT)=
     &     ROTAP(ICPX-NSHIFT:ICPX+NSHIFT,ICPY-NSHIFT:ICPY+NSHIFT)/DIVIS

           CALL FINDMX_Q(ROTAP,NEWS+2,NEWS,NEWR,NSHIFT,CMA,SXA,SYA)

           !write(6,*) 'CMA > CCC:',CMA,CCC,ir
           IF (CMA > CCC)  THEN
               CCC = CMA
               SX  = SXA
               SY  = SYA
               PHI = ANGLE
               IDI = IR
           ENDIF

C          CHECK ANGLE+180 POSITION
           ANGLE = ANGLE + 180.0
           CALL RTQ_Q(ALIGNEDPAD,ROTAP,NEWS+2,NEWS,NEWR,ANGLE)
        
C          HAS FORWARD FFT ON: NEWS,NEWR
           INV = +1
           CALL FMRS_2(ROTAP,NEWS,NEWR,INV)

C          HAS REVERSE FFT ON: NEWS,NEWR
           !CALL CCRS_2I(REFPAD(1,1,IR),ROTAP, LSC,NEWS,NEWR)
           CALL CCRS_2(REFPAD(1,1,IR),ROTAP,ROTAP, LSC,NEWS,NEWR)

C          FIND THE PEAK
           ROTAP(ICPX-NSHIFT:ICPX+NSHIFT,ICPY-NSHIFT:ICPY+NSHIFT) =
     &     ROTAP(ICPX-NSHIFT:ICPX+NSHIFT,ICPY-NSHIFT:ICPY+NSHIFT)/DIVIS

           CALL FINDMX_Q(ROTAP,NEWS+2,NEWS,NEWR,NSHIFT,CMM,SMX,SMY)

           IF (CMM > CCC)  THEN
              CCC = CMM
              SX  = SMX
              SY  = SMY
              PHI = ANGLE
              IDI = IR
           ENDIF
           !write(6,*) 'idi,phi,sx,sy,ccc:', idi,phi,sx,sy,ccc
        
           IF (MIRROR)  THEN

C           CHECK ANGLEM POSITION
            CALL RTQ_Q(ALPADM,ROTAP,NEWS+2,NEWS,NEWR,ANGLEM)

C           HAS FORWARD FFT ON: NEWS,NEWR
            INV = +1
            CALL FMRS_2(ROTAP,NEWS,NEWR,INV)

C           HAS REVERSE FFT ON: NEWS,NEWR
            !CALL CCRS_2I(REFPAD(1,1,IR),ROTAP, LSC,NEWS,NEWR)
            CALL CCRS_2(REFPAD(1,1,IR),ROTAP,ROTAP, LSC,NEWS,NEWR)

C           FIND THE PEAK
            ROTAP(ICPX-NSHIFT:ICPX+NSHIFT,ICPY-NSHIFT:ICPY+NSHIFT) =
     &        ROTAP(ICPX-NSHIFT:ICPX+NSHIFT,ICPY-NSHIFT:ICPY+NSHIFT) /
     &       DIVIS

            CALL FINDMX_Q(ROTAP,NEWS+2,NEWS,NEWR,NSHIFT,CMM,SMX,SMY)
            IF (CMM > CCC)  THEN
                CCC = CMM
                SX  = SMX
                SY  = SMY
                PHI = ANGLEM
                IDI = -IR
            ENDIF

C           CHECK ANGLEM+180 POSITION
            ANGLEM = ANGLEM+180.0
            CALL RTQ_Q(ALPADM,ROTAP,NEWS+2,NEWS,NEWR,ANGLEM)
        
C           HAS FORWARD FFT ON: NEWS,NEWR
            INV = +1
            CALL FMRS_2(ROTAP,NEWS,NEWR,INV)

C           HAS REVERSE FFT ON: NEWS,NEWR
            !CALL CCRS_2I(REFPAD(1,1,IR),ROTAP, LSC,NEWS,NEWR)
            CALL CCRS_2(REFPAD(1,1,IR),ROTAP,ROTAP, LSC,NEWS,NEWR)

C           FIND THE PEAK
            ROTAP(ICPX-NSHIFT:ICPX+NSHIFT,ICPY-NSHIFT:ICPY+NSHIFT)=
     &      ROTAP(ICPX-NSHIFT:ICPX+NSHIFT,ICPY-NSHIFT:ICPY+NSHIFT) /
     &          DIVIS

            CALL FINDMX_Q(ROTAP,NEWS+2,NEWS,NEWR,NSHIFT,CMM,SMX,SMY)
            IF (CMM > CCC)  THEN
               CCC = CMM
               SX  = SMX
               SY  = SMY
               PHI = ANGLEM
               IDI = -IR
            ENDIF
       write(6,*) 'idi,phi,sx,sy,ccc:', idi,phi,sx,sy,ccc

         ENDIF   ! END OF MIRROR CHECK

C       END OF DO-LOOP OVER REFERENCE IMAGES
        ENDDO

9999    IF (ALLOCATED(A_CIRC))     DEALLOCATE(A_CIRC)
        IF (ALLOCATED(ALIGNEDPAD)) DEALLOCATE(ALIGNEDPAD)
        IF (ALLOCATED(ROTAP))      DEALLOCATE(ROTAP)
        IF (ALLOCATED(ALPADM))     DEALLOCATE(ALPADM)

        END


C       -------------------- MACF_PL ------------------------------

        SUBROUTINE MACF_PL(NX,NY,X,NAS,NAR,IRA)

        REAL              :: X(NAS+2,NAR)
        DOUBLE PRECISION  :: AVE

        R   = IRA
        NS2 = NX/2+1
        NR2 = NY/2+1

        IF (NAS > NX)  THEN
c$omp      parallel do private(j,i)
           DO J=1,NAR
             DO I=NX+1,NAS
                X(I,J) = 0.0
             ENDDO       
           ENDDO
        ENDIF

        IF (NAR > NY)  THEN
c$omp      parallel do private(j,i)
           DO J=NY+1,NAR
              DO I=1,NX
                 X(I,J) = 0.0
              ENDDO
           ENDDO
        ENDIF

        AVE = 0.0
        ILE = 0
c$omp   parallel do private(j,i,a,tr) reduction(+:ave,ile)
        DO J=1,NY
           A = FLOAT(J-NR2)**2
           DO I=1,NX
             TR = SQRT(FLOAT(I-NS2)**2+A)
             IF (TR > R)  THEN
                X(I,J) = 0.0
             ELSE
                AVE = AVE+X(I,J)
                ILE = ILE + 1
             ENDIF
           ENDDO
         ENDDO
         AVE = AVE/ILE
        
c$omp    parallel do private(j,i,a,tr)
         DO J=1,NY
            A = FLOAT(J-NR2)**2
            DO I=1,NX
               TR = SQRT(FLOAT(I-NS2)**2+A)
               IF (TR .LE. R)  X(I,J) = X(I,J)-AVE
            ENDDO
         ENDDO 

C        HAS FORWARD FFT ON: NAS,NAR
         INS = +1
         CALL FMRS_2(X,NAS,NAR,INS)

         IF (INS .EQ. 0)  THEN
            CALL ERRT(38,'OR Q',NE)
            RETURN
         ENDIF

C        HAS REVERSE FFT ON: NAS,NAR
         CALL ACRS_2SL(X,X,NAS,NAR)

         NRL = 1
         NRU = NAR
         NSL = 1
         NSU = NAS
         D1  = 1. / REAL(NINT(3.1415926*R*R)*ILE)

c$omp    parallel do private(j,i,qt,a,t,m),shared(d1)
         DO J=NRL,NRU
            QT = FLOAT(J-(NY+1)) ** 2
            DO I=NSL,NSU
               A = SQRT(FLOAT(I-(NX+1))**2+QT)/2.0
               IF (A .EQ. 0.0)  THEN
                  X(I,J) = X(I,J) * D1
               ELSE
                  IF (R > A)  THEN
                     T = 2.0*ATAN(SQRT((R/A)**2-1.0))
C                    SHOULD BE NINT WITHOUT +0.5, BUT OMP WON'T TAKE IT...
                     M = INT(R*R*(T-SIN(T))+0.5)

C                   NORMALIZATION IS APPLIED TO THESE AC COEFF. WHICH WERE
C                   ESTIMATED FROM AT LEAST  5 PIXELS
C                   OTHERWISE AC COEFFS. ARE SET TO ZERO.

                    IF (M .GE. 5)  THEN
                       X(I,J) = X(I,J) / FLOAT(M) * ILE
                    ELSE
                       X(I,J) = 0.0
                    ENDIF
                 ELSE
                    X(I,J) = 0.0
                 ENDIF
              ENDIF
           ENDDO        
        ENDDO
        END

C       *************** ACRS_2SL ***********************************

C       STUPID DUPLICATION OF EXISTING ROUTINE!! al

        SUBROUTINE ACRS_2SL(X,O,NX,NY)

        COMPLEX          :: X((NX+2-MOD(NX,2))/2,NY)
        COMPLEX          :: O((NX+2-MOD(NX,2))/2,NY)
        DOUBLE PRECISION :: PI2

C       CALCULATES CIRCULAR AUTOCORRELATION, NON-POWER-OF-TWO DIMENSIONS
C       INPUT - X FOURIER TRANSFORM
C       OUTPUT -  O=F(X*CONJG(X))

        NNNN = (NX+2-MOD(NX,2))/2

        PI2  = 8.0*DATAN(1.0D0)
        ITMP = NX/2
        SX   = PI2*FLOAT(ITMP)/FLOAT(NX)
        ITMP = NY/2
        SY   = PI2*FLOAT(ITMP)/FLOAT(NY)

c$omp   parallel do private(i,j,ix,iy,argy,arg)
        DO J=1,NY
           IY = J-1
           IF (IY > NY/2)  IY=IY-NY
           ARGY = SY*IY
           DO I=1,NNNN
              IX     = I-1
              ARG    = SX*IX+ARGY
              O(I,J) = CABS(X(I,J)) * CMPLX(COS(ARG),SIN(ARG))
            ENDDO
        ENDDO

C       REVERSE FOURIER TRANSFORM ON: NX,NY
        INS = -1
C       FOR SOME REASON, IN THE PARALLEL MODE THE PROGRAM FAILS WHEN IT HAS TO
C       SWITCH BETWEEN SGI FFTS WITH DIFFERENT DIMENSIONS. SO, USE FORTRAN
C       CODE IN THE PARALLEL MODE.

        CALL FMRS_2(O,NX,NY,INS)

        END     


C       -------------------- ALRQ_LSD ------------------------------
C       STUPID DUPLICATION OF EXISTING ROUTINE!! al

        SUBROUTINE ALRQ_LSD
     &      (XIM,LSD,NX,NY,NUMR,CIRC,LCIRC,NRING,MODE)

        DIMENSION         XIM(LSD,NY),CIRC(LCIRC)
        INTEGER           NUMR(3,NRING)
        CHARACTER*1       MODE
        DOUBLE PRECISION  PI,DFI

C       INTERPOLATION INTO POLAR COORDINATES
C       NO NEED TO SET TO ZERO, ALL ELEMENTS ARE DEFINED

        NS2 = NX/2+1
        NR2 = NY/2+1
        PI  = 2*DATAN(1.0D0)

c$omp   parallel do private(i,j,inr,yq,l,lt,nsim,dfi,kcirc,
c$omp&                      xold,yold,fi,x,y)
        DO  I=1,NRING

C          RADIUS OF THE RING
           INR = NUMR(1,I)
           YQ  = INR

C          THE ACTUAL, POWER-OF-TWO LENGTH IS NUMR(3,I)-2, ADDITIONAL
C          TWO LOCATIONS ARE ONLY FOR THE NEW FFT.
           L = NUMR(3,I)
           IF (MODE .EQ. 'H')  THEN
              LT = L/2
           ENDIF
           IF (MODE .EQ. 'F')  THEN
              LT = L/4
           ENDIF
           NSIM  = LT-1
           DFI   = PI/(NSIM+1)
           KCIRC = NUMR(2,I)
           XOLD  = 0.0
           YOLD  = INR
           CIRC(KCIRC)=QUADRI_Q(XOLD+NS2,YOLD+NR2,LSD,NX,NY,XIM)
           XOLD = INR
           YOLD = 0.0
           CIRC(LT+KCIRC)=QUADRI_Q(XOLD+NS2,YOLD+NR2,LSD,NX,NY,XIM)
           IF (MODE .EQ. 'F')  THEN
              XOLD = 0.0
              YOLD = -INR
        CIRC(LT+LT+KCIRC)=QUADRI_Q(XOLD+NS2,YOLD+NR2,LSD,NX,NY,XIM)
              XOLD  = -INR
              YOLD  = 0.0
              CIRC(LT+LT+LT+KCIRC) =
     &        QUADRI_Q(XOLD+NS2,YOLD+NR2,LSD,NX,NY,XIM)
           ENDIF
           DO J=1,NSIM
              FI = DFI*J
              X = SIN(FI)*YQ
              Y = COS(FI)*YQ

              XOLD=X
              YOLD=Y
        CIRC(J+KCIRC) = QUADRI_Q(XOLD+NS2,YOLD+NR2,LSD,NX,NY,XIM)
              XOLD=Y
              YOLD=-X
        CIRC(J+LT+KCIRC)=QUADRI_Q(XOLD+NS2,YOLD+NR2,LSD,NX,NY,XIM)
              IF (MODE .EQ. 'F') THEN
                 XOLD=-X
                 YOLD=-Y
                 CIRC(J+LT+LT+KCIRC)=
     &           QUADRI_Q(XOLD+NS2,YOLD+NR2,LSD,NX,NY,XIM)
                 XOLD=-Y
                 YOLD=X
                 CIRC(J+LT+LT+LT+KCIRC)=
     &           QUADRI_Q(XOLD+NS2,YOLD+NR2,LSD,NX,NY,XIM)
              ENDIF
           ENDDO
        ENDDO
        END

C       ********************** NEARESTFFTFRIEND *******************************

        INTEGER  FUNCTION NEARESTFFTFRIEND(N)

        LOGICAL   DIVISIBLE
        EXTERNAL  DIVISIBLE

C       MAKE THE ARGUMENT DIVISIBLE BY 2
        L = N + MOD(N,2)

        DO WHILE (.NOT.DIVISIBLE(L))
         L = L + 2
        ENDDO
        NEARESTFFTFRIEND =  L

        END

C       ********************** DIVISIBLE *******************************

        LOGICAL FUNCTION DIVISIBLE(LA)

C       CHECK WHETHER THE ARGUMENT IS DIVISIBLE BY 2,3,5

        PARAMETER  (ND=3)
        DIMENSION  IDV(ND)

        IDV(1) = 5
        IDV(2) = 3
        IDV(3) = 2
        L      = LA

        DO I=1,ND
           DO WHILE(MOD(L,IDV(I)) .EQ. 0)
              L = L / IDV(I)
              IF (L .EQ. 1)  THEN
                 DIVISIBLE = .TRUE.
                 RETURN
              ENDIF
           ENDDO
        ENDDO

        DIVISIBLE = .FALSE.
        END


C       ********************** ANGMOR *******************************

        REAL FUNCTION  ANGMOR(RKK,MODE,MAXRIN)

        IMPLICIT NONE

        REAL         :: RKK
        CHARACTER*1  :: MODE
        INTEGER      :: MAXRIN

        IF (MODE == 'H')  THEN
            ANGMOR = (RKK - 1.0) / MAXRIN * 180.0
            IF (ANGMOR > 0.0)  ANGMOR = 360.0 - ANGMOR
	    ANGMOR = AMOD(ANGMOR+180.0, 180.0)

        ELSEIF (MODE == 'F')  THEN
            ANGMOR = (RKK - 1.0)/ MAXRIN * 360.0
            IF (ANGMOR > 0.0)  ANGMOR = 360.0 - ANGMOR
	    ANGMOR = AMOD(ANGMOR+360.0, 360.0)
        ENDIF

        END
