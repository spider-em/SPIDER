C++*********************************************************************
C                                                                      *
C  WRITPRO_N.F    SPEEDED UP                    FEB 2000 ARDEAN LEITH  *
C                 VERBOSE OUTPUT                NOV 2000 ARDEAN LEITH  *
C                 PUT ANGLES IN HEADER          JUN 2001 ARDEAN LEITH  *
C                 OPFILEC                       FEB 2003 ARDEAN LEITH  *
C                 PJ_RRINGS                     FEB 2005 ARDEAN LEITH  *
C                 APRINGS_NEW                   APR 2008 ARDEAN LEITH  *
C                 APRINGS_INIT_PLANS PARAMS     JUN 2011 ARDEAN LEITH  *
C                 WPRO_FBS                      DEC 2011 G. KISHCHENKO *
C                                                                      *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2011  Health Research Inc.,                         *
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
C WRITPRO_N(PROJ,LUNPROJ,NX,NY,NUMTH,BCKE,NNN,                         *
C         IPCUBE,NN,RI,ISELECT,NANG,MAXKEY,ANGBUF)                     *
C         LUNRINGS,MODE,MR,NR,ISKIP,LDPX,LDPY,LDPZ,IRTFLG              *
C                                                                      *
C PARAMETERS:                                                          *
C             PROJ        ARRAY OF PROJECTIONS                   RET   *
C             NUMTH       NUMBER OF OMP THREADS                  SENT  *
C             BCKE        VOLUME                                 SENT  *
C             RI          RADIUS                                 SENT  *
C                                                                      *
C PURPOSE: COMPUTES A PROJECTION OF A 3D VOLUME ACCORDING TO THE       *
C          THREE EULERIAN ANGLES CAN ALSO SAVE REFERENCE RINGS IN      *
C          POLAR FORM                                                  *
C                                                                      *
C IPCUBE: A RUN LENGTH LIST OF VOXELS ON EACH LINE IN THE              *
C         VOLUME  WHICH ARE WITHIN A SPECIFED RADIUS SPHERE IN VOL.    *
C                1 - BEGINNING VOXEL ON LINE                           *
C                2 - ENDING VOXEL ON LINE                              *
C                3 - IX FOR VOXEL                                      *
C                4 - IY FOR VOXEL                                      *
C                5 - IZ FOR VOXEL                                      *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE WRITPRO_N(PROJ,LUNPROJ, NX,NY,NZ, NUMTH,
     &                BCKE,NNN,IPCUBE,NN,RI,ISELECT,NANG,MAXKEY,ANGBUF,
     &                LUNRINGS,MODE,MR,NR,ISKIP,LDPX,LDPY,LDPZ,
     &                FBS_WANTED,IRTFLG)

        USE TYPE_KINDS                                          

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL                  :: PROJ(NX,NY,NUMTH)
        INTEGER               :: LUNPROJ, NX,NY,NZ, NUMTH
        REAL                  :: BCKE(NX,NY,NZ)
        INTEGER               :: NNN
        REAL                  :: IPCUBE(5,NN)
        REAL                  :: RI
        INTEGER               :: NN,NANG,MAXKEY
        INTEGER               :: ISELECT(MAXKEY)
        REAL                  :: ANGBUF(3,MAXKEY)
        INTEGER               :: LUNRINGS
        CHARACTER(LEN=*)      :: MODE
        INTEGER               :: MR,NR,ISKIP,LDPX,LDPY,LDPZ
        LOGICAL               :: FBS_WANTED
        INTEGER               :: IRTFLG

        REAL                  :: BUFOUT(4)

        CHARACTER(LEN=MAXNAM) :: PROJNAM,PROJPAT,REFNAM
        LOGICAL               :: REFRINGS

        INTEGER, ALLOCATABLE  :: NUMR(:,:)
        REAL, ALLOCATABLE     :: WR(:)

        REAL, ALLOCATABLE     :: XYZ  (:,:,:)   ! BCKE --> XYZ
        REAL, ALLOCATABLE     :: X1   (:,:,:)
        REAL, ALLOCATABLE     :: Y1   (:,:,:)
        REAL, ALLOCATABLE     :: Z1   (:,:,:)
        REAL, ALLOCATABLE     :: XY2  (:,:,:)
        REAL, ALLOCATABLE     :: XZ2  (:,:,:)
        REAL, ALLOCATABLE     :: YZ2  (:,:,:)

        INTEGER, PARAMETER   :: NPLANS = 14
        INTEGER(KIND=I_8)    :: FFTW_PLANS(NPLANS)

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

        REFRINGS = (LUNRINGS > 0) 
        IF (REFRINGS) THEN

C          FIND NUMBER OF REFERENCE-RINGS
           NRING = ((NR - MR) / ISKIP) + 1

C          FILL NUMR WITH RING ADDRESSES
           ALLOCATE(NUMR(3,NRING), WR(NRING),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'NUMR & WR',4*NRING)
              GOTO 9999
           ENDIF

C          NUMBR(1,*) IS RING NUMBER
           NRING = 0
           DO I=MR,NR,ISKIP
              NRING = NRING + 1
              NUMR(1,NRING) = I
	   ENDDO

C          CALCULATES DATA FOR NUMR & LCIRC
           CALL ALPRBS_Q(NUMR,NRING,LCIRC,MODE)

C          RINGWE RETURNS WR
	   CALL RINGWE_NEW(WR,NUMR,NRING,NUMR(3,NRING))
           IF (MODE .EQ. 'H') WR = WR * 0.5

C          INITIALIZE FFTW3 PLANS FOR USE WITHIN OMP || SECTIONS
           CALL APRINGS_INIT_PLANS(NUMR,NRING,
     &                          FFTW_PLANS,NPLANS,0,0,IRTFLG)
         ENDIF

C        GET PROJPAT TEMPLATE ONLY (NOT ILIST)
         NMAX = 0
         CALL FILSEQP(PROJPAT,NLET,ILIST,NMAX,NIMA,
     &        'TEMPLATE FOR 2-D PROJECTION',IRTFLG)
         IF (IRTFLG .NE. 0)  GOTO 9999

         IF (REFRINGS)  THEN
C           OPEN REF RINGS FILE USING SPIDER FORMAT
            MAXIM = 0
            CALL OPFILEC(0,.TRUE.,REFNAM,LUNRINGS,'U',IFORM,
     &            LCIRC,NANG,1,MAXIM,'REFERENCE RINGS',.FALSE.,IRTFLG)
            IF (IRTFLG .NE. 0)  GOTO 9999
         ENDIF

         IF (VERBOSE .AND. MYPID .LE. 0) WRITE(NOUT,*) ' '

         IF (FBS_WANTED) THEN

            NXLD   = NX + 2 - MOD(NX,2)  ! PAD FOR FFTW

C           CREATE SPACE FOR DERIVATIVES
            ALLOCATE (XYZ  (NXLD,NY,NZ),
     &                X1   (NXLD,NY,NZ),
     &                Y1   (NXLD,NY,NZ),
     &                Z1   (NXLD,NY,NZ),
     &                XY2  (NXLD,NY,NZ),
     &                XZ2  (NXLD,NY,NZ),
     &                YZ2  (NXLD,NY,NZ), STAT=IRTFLG)
            IF (IRTFLG .NE. 0) THEN             
               CALL ERRT(46,'WRITPRO_N; XYZ...',7*NXLD*NY*NZ)
               GOTO 9999                          
            ENDIF

C           COPY BCKE INTO XYZ AND PAD EACH ROW
            DO K=1,NZ
               DO J=1,NY
                  DO I=1,NX
                    XYZ(I,J,K) = BCKE(I,J,K)
                  ENDDO
               ENDDO
            ENDDO

C           CALCULATE DERIVATIVES USING 3D FFT
C           RETURNS: XYZ, X1, Y1, Z1, XY2, XZ2, YZ2
            CALL FBS3_PREP(XYZ, NXLD,NX,NY,NZ,
     &                     X1,Y1,Z1,XY2,XZ2,YZ2)

         ENDIF   ! END OF: IF (FBS_WANTED

         IFORM = 1
         NSL   = 1
         IANG  = 1
         DO WHILE (IANG .LE. NANG)

C           BREAK INTO NUMTH CHUNKS, CALCULATE PROJECTIONS
            NEEDED = MIN(IANG+NUMTH-1,NANG)

            IF (FBS_WANTED) THEN

c$omp          parallel do private(i,ifile)
               DO I=IANG,NEEDED
                  IFILE  = ISELECT(I)

                  CALL WPRO_FBS(PROJ(1,1,I-IANG+1),
     &               NXLD,NX,NY,NZ,  BCKE,
     &               X1,Y1,Z1,XY2,XZ2,YZ2,XYZ,
     &               IPCUBE,NN,
     &               ANGBUF(3,IFILE),ANGBUF(2,IFILE),ANGBUF(1,IFILE),
     &               RI, LDPX,LDPY,LDPZ)
               ENDDO
            ELSE
c$omp          parallel do private(i,ifile)
               DO I=IANG,NEEDED
                  IFILE  = ISELECT(I)

                  CALL WPRO_N(PROJ(1,1,I-IANG+1),
     &               NX,NY,NZ,BCKE,
     &               NNN,IPCUBE,NN,
     &               ANGBUF(3,IFILE),ANGBUF(2,IFILE),ANGBUF(1,IFILE),
     &               RI, LDPX,LDPY,LDPZ)
               ENDDO
            ENDIF

C           WRITE PROJECTIONS TO OUTPUT FILES
            DO I=IANG,NEEDED
               IFILE  = ISELECT(I)

C              CREATE OUTPUT FILENAME
               CALL FILGET(PROJPAT,PROJNAM,NLET,IFILE,IRTFLG)

C              OPEN NEXT PROJECTION OUTPUT FILE
               MAXIM = 0
               CALL OPFILEC(0,.FALSE.,PROJNAM,LUNPROJ,'U',IFORM,
     &                      NX,NY,NSL,MAXIM,' ',.FALSE.,IRTFLG)

               CALL WRTVOL(LUNPROJ,NX,NY,1,1,
     &                     PROJ(1,1,I-IANG+1),IRTFLG)

               IF (VERBOSE .AND. MYPID .LE. 0) THEN
                  WRITE(NOUT,90) IFILE,
     &               ANGBUF(1,IFILE),ANGBUF(2,IFILE),ANGBUF(3,IFILE)
90                FORMAT('  PROJECTION:',I6,
     &                   '  PSI:',F6.1,' THETA:',F6.1,' PHI:',F6.1)
               ENDIF

C              PUT ANGLES IN HEADER ALSO
               BUFOUT(1) = ANGBUF(1,IFILE)  !PSI
               BUFOUT(2) = ANGBUF(2,IFILE)
               BUFOUT(3) = ANGBUF(3,IFILE)
               BUFOUT(4) = 1.0
               CALL LUNSETVALS(LUNPROJ,IAPLOC+1,4,BUFOUT,IRTFLG)

               CLOSE(LUNPROJ)
            ENDDO

            IF (REFRINGS)  THEN
C              CALCULATE & WRITE REFERENCE RINGS 
               CALL PJ_RRINGS(NX,NY,
     &                        NRING,LCIRC,NUMR,MODE, 
     &                        NUMTH,LUNRINGS,WR,FFTW_PLANS,
     &                        PROJ,IANG,NEEDED,IRTFLG)
            ENDIF

            IANG = IANG + NUMTH
         ENDDO

9999     IF (REFRINGS) THEN
            IF (ALLOCATED(NUMR))  DEALLOCATE(NUMR)
            IF (ALLOCATED(WR))    DEALLOCATE(WR)
         ENDIF

         IF (ALLOCATED(X1))    DEALLOCATE(X1)
         IF (ALLOCATED(Y1))    DEALLOCATE(Y1)
         IF (ALLOCATED(Z1))    DEALLOCATE(Z1)
         IF (ALLOCATED(XY2))   DEALLOCATE(XY2)
         IF (ALLOCATED(XZ2))   DEALLOCATE(XZ2)
         IF (ALLOCATED(YZ2))   DEALLOCATE(YZ2)
         IF (ALLOCATED(XYZ))   DEALLOCATE(XYZ)

         END



C       ---------------------- PJ_RRINGS  --------------------------

        SUBROUTINE PJ_RRINGS(NX,NY,
     &                       NRING,LCIRC,NUMR, MODE, 
     &                       NUMTH,LUNRINGS,WR, FFTW_PLANS,
     &                       PROJ, IANG1,IANG2, IRTFLG)

        USE TYPE_KINDS                                          

        IMPLICIT NONE

        INTEGER               :: NX,NY,NRING,LCIRC
        INTEGER               :: NUMR(3,NRING)
        CHARACTER(LEN=1)      :: MODE
        INTEGER               :: NUMTH,LUNRINGS
        REAL                  :: WR(NRING) 
        INTEGER(KIND=I_8)     :: FFTW_PLANS(*) ! FFTW_PLANS STRUCTURE
	REAL                  :: PROJ(NX,NY,NUMTH)
        INTEGER               :: IANG1,IANG2,IRTFLG

C       AUTOMATIC ARRAYS
        REAL                  :: CIRCREF(LCIRC,NUMTH)

        INTEGER               :: IMI,IT
        REAL                  :: CNS2,CNR2

        IRTFLG = 0

C       CALCULATE CENTER FOR APRINGS
        CNS2 = NX / 2 + 1 
        CNR2 = NY / 2 + 1

c$omp   parallel do private(imi,it)
	DO IMI=IANG1,IANG2
           IT = IMI - IANG1 + 1

C          NORMALIZE IMAGE VALUES UNDER THE MASK OVER VARIANCE RANGE
C          INTERPOLATE TO POLAR COORDINATES, CREATE FOURIER OF: CIRCREF
C          WEIGHT CIRCREF USING WR
           CALL APRINGS_ONE_NEW(NX,NY, CNS2,CNR2,
     &                          PROJ(1,1,IT),.FALSE.,
     &                          MODE,NUMR,NRING,LCIRC,WR,FFTW_PLANS,
     &                          CIRCREF(1,IT),IRTFLG)
        ENDDO

C       SAVE CIRCREF IN FILE OPENED ON LUNRINGS
        DO IMI=IANG1,IANG2
           IT  = IMI - IANG1 + 1
	   CALL WRTLIN(LUNRINGS,CIRCREF(1,IT),LCIRC,IMI)

c          write(6,*) 'oCIRC(1-100):',CIRCREF(1,1),CIRCREF(100,1)

        ENDDO

	END


#ifdef NEVER
C       I TRIED THIS BUT IT WAS A BIT SLOWER !!!!!!!!!!!--------------
C       -------------- PJ_RRINGSNEW-----------------------------------

        SUBROUTINE PJ_RRINGSNEW(NX,NY,
     &                       NRING,LCIRC,NUMR,MODE, 
     &                       NUMTH,LUNRINGS,WR,
     &                       PROJ,IANG1,IANG2,IRTFLG)

	INCLUDE 'CMBLOCK.INC' 
	INCLUDE 'CMLIMIT.INC' 

        CHARACTER(LEN=1)                    :: MODE

C       EXTERNAL ARRAYS
        INTEGER,DIMENSION(3,NRING)          :: NUMR
	REAL,DIMENSION(NX,NY,NUMTH)     :: PROJ
        REAL,DIMENSION(NRING)               :: WR

C       AUTOMATIC ARRAYS
        REAL,DIMENSION(LCIRC,NUMTH)         :: CIRCREF

        IRTFLG = 0

C       PREPARE CIRCULAR RINGS DATA FOR REFERENCE IMAGES
c$omp   parallel do private(IMI,IT)
	DO IMI=IANG1,IANG2
           IT = IMI - IANG1 + 1     !PROJ & CIRCREF INDEX
C          NORMALIZE, INTERPOLATE INTO POLAR COORDINATES, FFT, & WEIGHT
	   CALL APRINGS_ONE(NX,NY,WR,
     &                      PROJ(1,1,IT),CIRCREF(1,IT),
     &                      MODE,NUMR,NRING,LCIRC, IRTFLG)
        ENDDO

C       SAVE CIRCREF IN FILE OPENED ON LUNRINGS
        DO IMI=IANG1,IANG2
           IT  = IMI - IANG1 + 1
	   CALL WRTLIN(LUNRINGS,CIRCREF(1,IT),LCIRC,IMI)

           if (imi.eq.1) THEN
              write(6,*) 'n(1-100):',CIRCREF(1,1),CIRCREF(100,1),numth
           endif
        ENDDO
	END

#endif


