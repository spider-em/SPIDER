C++*********************************************************************
C                                                                      *
C APRINGS_PATM.F  NEW                              AUG 11 ARDEAN LEITH *
C **********************************************************************
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
C
C  APRINGS_NEW_PATM(ILIST,NUMREF,  NX,NY,
C                  NRING,LCIRC,NUMR,MODE,FFTW_PLANS, 
C                  REFPAT,LUNREF,CIRCREF,CIRCREF_IN_CORE,
C                  LUNRING,SCRFILE, 
C                  AVR,SIGR, IRTFLG)
C 
C  PURPOSE: CREATES SETS OF RADIAL/POLAR RINGS FROM SETS OF IMAGES
C 
C  PARAMETERS: ILIST       LIST OF IMAGE FILE NUMBERS          (INPUT)
C              NUMREF      NO. OF IMAGES                       (INPUT)
C              NX,NY       IMAGE DIMENSIONS                    (INPUT)
C              CMASK       MASK                                (INPUT)
C              NRING       NUMBER OF POLAR RINGS               (INPUT)
C              LCIRC       LENGTH OF CIRCULAR RING ARRAY       (INPUT)
C              NUMR        RING INDICES                        (INPUT)
C              MODE        FULL OR HALF CIRCLE                 (INPUT)
C              FFTW_PLANS  PLAN POINTERS                       (INPUT)
C              REFPAT      IMAGE SERIES FILE TEMPLATE          (INPUT)
C              LUNREF      IMAGE FILE IO UNIT                  (INPUT)
C              CIRCREF     OUTPUT ARRAY                        (OUTPUT)
C              CIRCREF_IN_CORE   NO OUTPUT ARRAY FLAG          (OUTPUT)
C              LUNRING     REF-RINGS FILE IO UNIT              (INPUT)
C              SCRFILE     REF-RINGS FILE                      (INPUT)
C              AVR         AVERAGE                             (OUTPUT)
C              SIGR        SD                                  (OUTPUT)
C              IRTFLG      ERROR FLAG                          (OUTPUT)
C
C NOTE: MOST MEMORY DEMAND DEPENDENT ON LCIRC & NUMREF.  LCIRC IS THE 
C       TOTAL LENGTH OF ARRAY THAT HOLDS THE CIRCULAR RINGS, SO IT IS 
C       DEPENDENT ON NUMBER OF RINGS AND THEIR RADIUS. ARRAY ALLOCATED 
C       IS: CIRCREF(LCIRC,NUMREF) ANOTHER SMALL ALLOCATED ARRAY IS: 
C       A(NX,NY,NUMTH)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

cpgi$g opt=3

        SUBROUTINE APRINGS_NEW_PATM(ILIST,NUMREF,  NX,NY, CMASK,
     &                     NRING,LCIRC,NUMR,MODE,FFTW_PLANS, 
     &                     REFPAT,LUNREF,CIRCREF,CIRCREF_IN_CORE,
     &                     LUNFFT,LUNRING,SCRFILE, 
     &                     FFTPAT, AVR,SIGR, IRTFLG)

        IMPLICIT NONE
	INCLUDE 'CMBLOCK.INC' 
	INCLUDE 'CMLIMIT.INC' 

        INTEGER                  :: ILIST(NUMREF)
        INTEGER                  :: NUMREF,NX,NY,NRING,LCIRC
        REAL                     :: CMASK(NX,NY)
        INTEGER                  :: NUMR(3,NRING) 
	CHARACTER (LEN=1)        :: MODE
        INTEGER*8                :: FFTW_PLANS(*) !POINTERS TO STRUCTURES
        CHARACTER (LEN=*)        :: REFPAT
        INTEGER                  :: LUNREF 
        REAL                     :: CIRCREF(LCIRC,*)
        LOGICAL                  :: CIRCREF_IN_CORE
        INTEGER                  :: LUNFFT,LUNRING
        CHARACTER (LEN=*)        :: SCRFILE
        CHARACTER (LEN=*)        :: FFTPAT
	REAL                     :: AVR(NUMREF),SIGR(NUMREF)
        INTEGER                  :: IRTFLG 

        CHARACTER (LEN=MAXNAM)   :: FILNAM
        LOGICAL                  :: USEREFFILE
        LOGICAL                  :: ISOPEN
        LOGICAL                  :: WEIGHT

        INTEGER                  :: ICOMM,MYPID,MPIERR
        INTEGER                  :: ILOCAT,INULL,NLET,NDUM,LUNOP
        INTEGER                  :: NUMTH,LCIRCT,NUMREFT,NZ,MAXIM
        INTEGER                  :: IMGNUM,NE,IREF
        REAL                     :: DUM
        LOGICAL                  :: INLNED
        INTEGER                  :: lnblnkn,lnblnk
        CHARACTER(LEN=3)         :: TYPET = 'FI'

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID 

C       SEE IF THIS FILE EXISTS, (RETURNS EX, ISOPEN, LUNOP)

C       SEE IF SCRATCH "FILE" EXISTS (MAY BE INCORE FILE)
        ILOCAT = INDEX(SCRFILE,'@')
        INULL  = lnblnkn(SCRFILE)
        IF (INULL > 0 .AND. SCRFILE(1:1) .NE. '_' .AND. 
     &      ILOCAT .EQ. 0) THEN
C          ADD EXTENSION TO PHYSICAL FILENAME
           CALL FILNAMANDEXT(SCRFILE,DATEXC,FILNAM,NLET,.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ELSE
           FILNAM = SCRFILE
           NLET   = lnblnk(FILNAM)
        ENDIF
#ifdef USE_MPI
        IF (MYPID <= 0) THEN
           INQUIRE(FILE=FILNAM,EXIST=USEREFFILE,OPENED=ISOPEN,
     &             NUMBER=LUNOP,IOSTAT=IRTFLG)
        ENDIF
        CALL BCAST_MPI('APRINGS','USEREFFILE',USEREFFILE,1, 'L',ICOMM)
#else
        USEREFFILE = .FALSE. 
        IF (INULL > 0) THEN
           CALL INQUIREIF1(LUNRING,FILNAM,DUM,NDUM,USEREFFILE,ISOPEN,
     &                     LUNOP,INLNED,IMGNUM,IRTFLG)
        ENDIF
#endif

        IF (USEREFFILE .AND. MYPID <= 0) THEN 
           WRITE(NOUT,*) ' Using existing reference rings file: ',
     &                     FILNAM(1:NLET)
        ELSEIF (INULL > 0  .AND. MYPID <= 0) THEN 
           WRITE(NOUT,*) ' No existing reference rings file: ',
     &                     FILNAM(1:NLET)
        ENDIF

        WEIGHT = .TRUE.   ! REFERENCES ARE WEIGHTED

C       NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)
        NUMTH = MIN(NUMTH,NUMREF)

        IF (CIRCREF_IN_CORE  .AND. .NOT. USEREFFILE) THEN
C          CALCULATE REF. RINGS DATA AND FILL CIRCREF ARRAY WITH IT
C
           CALL APRINGS_FILL_PATM(ILIST,NUMREF,  
     &                      NX,NY,NUMTH,CMASK,
     &                      NRING,LCIRC,NUMR,MODE,FFTW_PLANS, 
     &                      REFPAT,LUNREF,LUNFFT,
     &                      CIRCREF,NUMREF,0,WEIGHT,
     &                      FFTPAT,AVR,SIGR,IRTFLG)

          !IF (MYPID .LE. 0) WRITE(NOUT,*) ' Created incore reference rings'

        ELSEIF (CIRCREF_IN_CORE .AND. USEREFFILE) THEN
C          READ EXISTING REF RINGS FILE AND FILL CIRCREF ARRAY WITH ITS DATA

C          OPEN EXISTING REFERENCE RINGS FILE
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,SCRFILE,LUNRING,'O',IFORM,
     &                 LCIRCT,NUMREFT,NZ,MAXIM,' ', .TRUE.,IRTFLG)
	   IF (IRTFLG .NE. 0) RETURN

           IF (LCIRCT .NE. LCIRC .OR. NUMREFT .NE. NUMREF) THEN
               CALL ERRT(101,'REF. RINGS FILE HAS WRONG SIZE',NE)
               IRTFLG = 1
               RETURN
           ENDIF

C          FILL CIRCREF WITH EXISTING RINGS DATA FROM FILE AND RETURN.
           DO IREF=1,NUMREF
              CALL REDLIN(LUNRING,CIRCREF(1,IREF),LCIRC,IREF) 
           ENDDO

           IF (MYPID <= 0) THEN
             WRITE(NOUT,*) ' Loaded reference rings file incore'
           ENDIF

        ELSEIF (.NOT. CIRCREF_IN_CORE .AND. USEREFFILE) THEN
C          OPEN EXISTING REF. RINGS FILE TO READ REF RINGS DATA
 
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,SCRFILE,LUNRING,'O',IFORM,
     &                 LCIRCT,NUMREFT,NZ,MAXIM,' ', .FALSE.,IRTFLG)
	   IF (IRTFLG .NE. 0) RETURN

           IF (LCIRCT .NE. LCIRC .OR. NUMREFT .NE. NUMREF) THEN
               CALL ERRT(101,'REF. RINGS FILE HAS DIFFERENT SIZE',NE)
               IRTFLG = 1
               RETURN
           ENDIF

        ENDIF

        END

C       --------------------- APRINGS_FILL_PATM --------------------------

        SUBROUTINE APRINGS_FILL_PATM(ILIST,NUMREF,
     &                        NX,NY,NUMTH, CMASK,
     &                        NRING,LCIRC,NUMR,MODE,FFTW_PLANS, 
     &                        REFPAT,LUNREF,LUNFFT,
     &                        CIRCREF,ICORE,LUNRING,WEIGHT,
     &                        FFTPAT, AVR,SIGR, IRTFLG)

        IMPLICIT NONE
	INCLUDE 'CMBLOCK.INC' 
	INCLUDE 'CMLIMIT.INC' 

        INTEGER              :: ILIST(NUMREF)
        INTEGER              :: NUMREF,NX,NY,NUMTH,NRING,LCIRC
	REAL                 :: CMASK(NX,NY)
        INTEGER              :: NUMR(3,NRING)
        CHARACTER(LEN=1)     :: MODE
        INTEGER*8            :: FFTW_PLANS(*) ! POINTERS TO STRUCTURES
        CHARACTER (LEN=*)    :: REFPAT 
        INTEGER              :: LUNREF,LUNFFT
        REAL                 :: CIRCREF(LCIRC,ICORE)
        INTEGER              :: ICORE,LUNRING
        LOGICAL              :: WEIGHT                   ! INPUT
        CHARACTER (LEN=*)    :: FFTPAT 
        REAL                 :: AVR(NUMREF),SIGR(NUMREF) ! OUTPUT
        INTEGER              :: IRTFLG,IRT

C       AUTOMATIC ARRAYS
        REAL                 :: WR(NRING)
        REAL                 :: AVRT(NUMTH),SIGRT(NUMTH)  ! TEMP USE

C       ALLOCATABLE ARRAYS
	REAL, ALLOCATABLE    :: BUFFFT (:,:,:)
	REAL, ALLOCATABLE    :: BUFPATM(:,:,:)
	REAL, ALLOCATABLE    :: BUFWORK(:,:)
	REAL, ALLOCATABLE    :: BUFWORK1(:,:)

	!real, allocatable    :: bufwork2(:,:)

        INTEGER              :: MAXRIN,NXLD,MWANT,IGO,IEND
        INTEGER              :: IMI,IT,IPT,IV,IREF,NLET,ITYPE
        INTEGER              :: ICOMM,MYPID,MPIERR 
        INTEGER              :: N2X,N2Y,N2XLD,IDUM 
        REAL                 :: X2CEN,Y2CEN,PADVAL,MAXIM,INV
        INTEGER * 8          :: IPLAN        = 0     ! STRUCTURE POINTER 
        LOGICAL, PARAMETER   :: MPIBCAST     = .TRUE.
        LOGICAL, PARAMETER   :: WANTSTATS    = .TRUE.
        LOGICAL, PARAMETER   :: SPIDER_SCALE = .FALSE.
        LOGICAL, PARAMETER   :: SPIDER_SIGN  = .FALSE.
  
        CHARACTER(LEN=MAXNAM):: FILNAM  

        INTEGER              :: NDUM

        integer :: i,j,ii,jj,kernel  = 4
        real    :: dfactor,sig1,grad
        real    :: avall,sigall
        logical :: savit

        !character(len=10)    :: filPATM = 'jnkPATM0000'
        !integer,save         :: nfil   = 0

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID 
        IRTFLG = 0

        !IF (MYPID <= 0)WRITE(NOUT,'(A,i5)')' Fill ref rings with threads:',NUMTH

        IF (WEIGHT) THEN
           MAXRIN = NUMR(3,NRING) - 2

C          RINGWE RETURNS WR WEIGHTS
	   CALL RINGWE_NEW(WR,NUMR,NRING,MAXRIN)
           IF (MODE == 'H') WR = WR * 0.5
        ELSE
           WR(1) = 0.0
        ENDIF

        N2X   = NX  * 2
        N2Y   = NY  * 2
        NXLD  = NX  + 2 - MOD(NX, 2)
        N2XLD = N2X + 2 - MOD(N2X,2)

        MWANT = N2XLD*N2Y*NUMTH + N2XLD*N2Y + N2X*N2Y*NUMTH +NXLD*NY
        ALLOCATE(BUFFFT (N2XLD,N2Y,NUMTH),
     &           BUFWORK(N2XLD,N2Y),
     &           BUFWORK1(NXLD,NY),
     &           BUFPATM(N2X,  N2Y,NUMTH), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'APRINGS_FILL_PATM; BUFFFT...',MWANT)
           RETURN
        ENDIF 

C       OPEN FFT OUTPUT STACK FILE NOW
        ITYPE = -11
        IF (MOD(NX,2) == 0) ITYPE = -12
        MAXIM = 1            ! FOR STACK
        IREF  = ILIST(1)

        CALL FILGET(FFTPAT,FILNAM,0,IREF,IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9999

        !write(6,*) ' Opening :',NXLD,NROW,MAXIM,ITYPE
        CALL OPFILEC(0,.FALSE.,FILNAM,LUNFFT,'N',ITYPE,
     &               NXLD,NY,1,MAXIM,' ',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

	X2CEN  = N2X / 2 + 1
	Y2CEN  = N2Y / 2 + 1

        PADVAL = 0.0

        !write(6,*) ' numref,numth: ', numref,numth

C       PREPARE CIRCULAR RINGS DATA FOR ALL REFERENCE IMAGES
C       LOOP OVER ALL REF. IMAGES, IN CHUNKS OF NUMTH
        DO IGO=1,NUMREF,NUMTH  
           IEND = MIN(IGO+NUMTH-1,NUMREF)

           !write(6,*) ' numref,numth: ', numref,numth,igo,iend

C          LOAD SET OF REF. IMAGES (NUMTH IMAGES IN SET)
C          PAD INTO BUFFFT FOR FFT, PRESERVE IMAGE AV & SIG.
           CALL AP_GETDATA(ILIST,NUMREF, 
     &                    NX,NY, N2XLD,N2Y,PADVAL,
     &                    NUMTH,REFPAT,LUNREF, IGO,IEND,
     &                    MPIBCAST, BUFFFT, 
     &                    WANTSTATS,AVRT,SIGRT,  IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
           !call chkfile('jnkbuffft',66,1,n2xld,n2y,1,buffft(1,1,iend),irt)

C          CONVERT PADDED REF. IMAGES TO PATTERSON MAPS
C          LOOP OVER NUMTH REF. IMAGES IN THIS CHUNK
	   DO IMI=IGO,IEND                   ! OVERALL INDEX     
              IT        = IMI - IGO + 1      ! NUMTH INDEX
              IREF      = ILIST(IMI)         ! ACTUAL REF IMAGE NUMBER

              AVR (IMI) = AVRT(IT)           ! RETURNED TO CALLER
              SIGR(IMI) = SIGRT(IT)          ! RETURNED TO CALLER
              !write(6,'(a,4i7)')' numref,igo,imi,it:',numref,igo,imi,it
              !write(6,'(a,i4,a,f7.4,f7.4)')'avrt(',imi,'):',avrt(imi),avrt(iref)

C             COPY THE 1X PADDED REF IMAGE INTO TEMP FILE (NOT FFT!)
              BUFWORK1(1:NXLD,1:NY) = BUFFFT(1:NXLD,1:NY,IT)

              !call chkfile('jnkforfft',66,1,nxld,ny,1,bufwork1,irt)
              !if (iref == 1) then
              !call chkfile('jnkrefpato',66,1,n2xld,n2y,1,buffft(1,1,1),irtflg)
              !endif

              INV = +1   ! INPLACE FORWARD FFT ON 1X UNPADDED REF IMAGE 
                         ! DO NOT USE: FMRS_2() AS IT SCALES
              CALL FMRS(BUFWORK1, NX,NY,1, IPLAN,
     &            SPIDER_SIGN,SPIDER_SCALE, INV,IRTFLG)
              IF (IRTFLG .NE. 0)  GOTO 9999
 
C             OPEN REF IMAGE FFT OUTPUT STACKED FILE
              CALL GETNEWIMG(0,LUNFFT,0,FFTPAT,IMI,
     &                       VERBOSE,IDUM,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999

C             PRESERVE THE FFT'D REF IMAGE IN TEMP FILE                          
              CALL WRTVOL(LUNFFT,NXLD,NY,1,1,BUFWORK1,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999

C             DENOISING REF IMAGE WITH MEANSHIFT OR MEDIAN DOES NOT
C             IMPROVE SELECTION!

C             CREATE PATTERSON MAPS FROM PADDED, MASKED REF. IMAGES
C             DOES A FFT ON 2X IMAGE THEN A BACK FFT
              !using: make_patm(avall, sigall, gives same values!

              CALL MAKE_PATM(AVRT(IT), SIGRT(IT),
     &                       BUFFFT(1,1,IT), CMASK,
     &                       NX,NY, N2XLD,N2X,N2Y, 
     &                       BUFWORK, BUFPATM(1,1,IT), IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999

              !write(6,'(a,5i7)')' igo,imi,it,iref:',igo,imi,it,iref
              !call chkfile('jnkrefpat',66,1,n2x,n2y,1,bufpatm(1,1,it),irtflg)
           ENDDO

c$omp      parallel do private(imi,ipt,it)
	   DO IMI=IGO,IEND      
              IT  = IMI - IGO + 1           ! BUFPATM INDEX
              IPT = IMI                     ! CIRCREF INDEX
              IF (LUNRING > 0) IPT = IT     ! CIRCREF INDEX 

C             CONVERT BUFPATM TO POLAR RINGS, FFT, & WEIGHT THE RINGS
 	      CALL APRINGS_ONE_NEW(N2X,N2Y, X2CEN,Y2CEN, 
     &                 BUFPATM(1,1,IT),
     &                .FALSE.,MODE,NUMR,NRING,LCIRC,WR,FFTW_PLANS,
     &                 CIRCREF(1,IPT),IRTFLG)
           ENDDO
c$omp      end parallel do

           !write(6,*) 'circref(1,1):', circref(1,1),bufpatm(30,30,1),bufpatm(60,60,1)

           IF (LUNRING > 0) THEN
C             SAVE CIRCREF IN FILE OPENED ON LUNRING
	      DO IMI=IGO,IEND
                 IV = IMI - IGO + 1 
                 CALL WRTLIN(LUNRING,CIRCREF(1,IV),LCIRC,IMI)
              ENDDO
           ENDIF
        ENDDO

C       DEALLOCATE ARRAYS
9999    IF (ALLOCATED(BUFFFT))   DEALLOCATE(BUFFFT)
        IF (ALLOCATED(BUFWORK))  DEALLOCATE(BUFWORK)
        IF (ALLOCATED(BUFWORK1)) DEALLOCATE(BUFWORK1)
        IF (ALLOCATED(BUFPATM))  DEALLOCATE(BUFPATM)
        CLOSE(LUNFFT)

	END


#ifdef DEBUGNEVER
c-----------------debug
       !write(6,*) ' File:',filnam(1:20),'   sigr:',sigr(iref)
       !call chkavg('fftd ring file1',buffft,nxld*NY)
c---------------------------
#endif


