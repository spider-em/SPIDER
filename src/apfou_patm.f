
C++*********************************************************************
C                                                                      *
C APFOU_PATM.F    ADAPTED FROM: APREF_PM          JUN 11 ARDEAN LEITH  *
C                                                                      *
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
C APFOU_PATM(IREFLIST,NUMREF,IEXPLIST,NUMEXP,
C            NX,NY,NR,RANGE,
C            NRING,LCIRC,NUMR,CIRCREF,CIRCREF_IN_CORE,
C            REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS,
C            REFPAT,EXPPAT,CKMIRROR,CTYPE,ISHRANGE,LUNDOC)
C
C PURPOSE: FIND ROTATIONAL AND SHIFT PARAMETERS TO ALIGN A SERIES OF
C          REFERENCE IMAGES WITH SAMPLE IMAGES
C
C PARAMETERS:
C       IREFLIST            LIST OF REF. IMAGE FILE NUMBERS   (INPUT)
C       NUMREF              NO. OF IMAGES                     (INPUT)
C       IEXPLIST            LIST OF EXP. IMAGE FILE NUMBERS   (INPUT)
C       NUMEXP              NO. OF IMAGES                     (INPUT)
C       NX,NY               IMAGE DIM.                        (INPUT)
C       RANGE               PROJ. ANGLE SEARCH RANGE (DEG)    (INPUT)
C       ANGDIFTHR           ANG. DIFFERENCE THRESHOLD         (INPUT)
C       RADI                MASK RADIUS                       (INPUT)
C       REFANGDOC           REF. ANGLES FILE NAME             (INPUT)
C       EXPANGDOC           EXP. ANGLES FILE NAME             (INPUT)
C       REFPAT              REF. IMAGE SERIES FILE TEMPLATE   (INPUT)
C       EXPPAT              EXP. IMAGE SERIES FILE TEMPLATE   (INPUT)
C       CTYPE               AP OP TYPE == 'FOU'               (INPUT)
C       LUNDOC              OUTPUT UNIT FOR DOC FILE          (INPUT)
C       FBS_WANTED          USE FBS INTERP IF ROTFIRST        (INPUT)
C
C OPERATIONS:  'AP FOU'
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE APFOU_PATM(IREFLIST,NUMREF, IEXPLIST,NUMEXP,
     &          NX,NY,RANGE,ANGDIFTHR,
     &          NRING,LCIRC,NUMR,CIRCREF,CIRCREF_IN_CORE,
     &          REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS,
     &          REFPAT,EXPPAT,CKMIRROR,CTYPE,
     &          ROTFIRST,ISHRANGEX,ISHRANGEY,LUNDOC,FBS_WANTED)

        IMPLICIT NONE

        INCLUDE 'MAKE_CLOSE_LIST.INC'  
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

	INTEGER                      :: IREFLIST(NUMREF) 
	INTEGER                      :: NUMREF 
	INTEGER                      :: IEXPLIST(NUMEXP) 
	INTEGER                      :: NUMEXP 
	INTEGER                      :: NX,NY
        REAL                         :: RANGE,ANGDIFTHR 
	INTEGER                      :: NRING,LCIRC 
        INTEGER                      :: NUMR(3,NRING)
	REAL                         :: CIRCREF(LCIRC,NUMREF)
	LOGICAL                      :: CIRCREF_IN_CORE
        CHARACTER (LEN=*)            :: REFANGDOC
        CHARACTER (LEN=*)            :: EXPANGDOC
        CHARACTER (LEN=*)            :: SCRFILE
        INTEGER *8                   :: FFTW_PLANS(*)
        CHARACTER (LEN=*)            :: REFPAT,EXPPAT 
	LOGICAL                      :: CKMIRROR
        CHARACTER (LEN=*)            :: CTYPE 
        LOGICAL                      :: ROTFIRST
        INTEGER                      :: ISHRANGEX,ISHRANGEY,LUNDOC
        LOGICAL                      :: FBS_WANTED

        CHARACTER (LEN=MAXNAM)       :: FFTPAT 
        CHARACTER (LEN=MAXNAM)       :: FILNAM
        CHARACTER (LEN=74)           :: COMMENT 
        LOGICAL                      :: MIRRORNEW, MIRRORNEWT
	LOGICAL                      :: GOTREFANG,GOTEXPANG
	LOGICAL                      :: ONLYMIRROR
        LOGICAL                      :: LIMITRANGE
        LOGICAL                      :: MIRRORED
        LOGICAL                      :: USE_UN,USE_MIR
        LOGICAL                      :: ANGINHEADER

C       ALLOCATABLE ARRAYS
c        LOGICAL, ALLOCATABLE         :: LMASKPAT(:,:)
        real, allocatable            :: cmaskpat(:,:)
        LOGICAL, ALLOCATABLE         :: LMASKCC(:,:)
        REAL, ALLOCATABLE            :: CIRCEXP(:)
	REAL, ALLOCATABLE            :: BUFEXP(:,:) 
	REAL, ALLOCATABLE            :: BUFPATM(:,:)
	REAL, ALLOCATABLE            :: BUFWORK(:,:)
	REAL, ALLOCATABLE            :: BUFEXP_P(:,:) 
        REAL, ALLOCATABLE            :: BUFREF_P(:,:)
        REAL, ALLOCATABLE            :: F0(:,:),X1(:,:),Y1(:,:),XY2(:,:)

C       AUTOMATIC ARRAYS
 	DOUBLE PRECISION             :: DCCROT(NUMREF)
 	LOGICAL                      :: ISMIRRORED(NUMREF)
 	REAL                         :: ROTPOS(NUMREF) 
	REAL                         :: DLIST(6)
	REAL                         :: ANGEXP(8,NUMEXP)
	REAL                         :: ANGREF(3,NUMREF)
	REAL                         :: REFDIR(3,NUMREF)
	REAL                         :: EXPDIR(3)
	REAL                         :: AVR(NUMREF),SIGR(NUMREF)
	INTEGER, POINTER             :: LCG(:)
        INTEGER, PARAMETER           :: NLISTMAX = 15
        REAL                         :: PARLIST(NLISTMAX) 

	CHARACTER (LEN=1)            :: NULL = CHAR(0) 

        CHARACTER (LEN=1)            :: MODE = 'H'! PS IS ONLY 0--->180
        REAL                         :: DIVAS = 180.0

        REAL                         :: RADPATMASK,RADCCMASK 

        INTEGER                      :: NXLD,MYPID,MAXRIN,NUMTH,MWANT
        INTEGER                      :: IRTFLG,NGOTREF
        REAL                         :: RANGECOS,PEAKV,PADVAL
        REAL                         :: SX,SY,CO,SO 

        INTEGER                      :: NSAID = 0

        INTEGER                      :: ILOC(1)
        LOGICAL                      :: ONLYONEROT
        INTEGER                      :: ITYPET,NXLDT,NGOT,NPROJT
        REAL                         :: CCROTT
        INTEGER                      :: IMIL,IMI,ITHREAD,NGOTPAR,IEXP
        REAL                         :: FDUM,CC,ANGDIF,CCLAS,CCAVG
        INTEGER                      :: IBIGANGDIF,IMPROVCC,IWORSECC
        REAL                         :: ANGDIFAVG,CCIMPROV,CCWORSE
        REAL                         :: F2XCEN,F2YCEN,AVI,SIGI,ROTANGNEW
        REAL                         :: ROTANGNEWT,PEAKVT,CCROT,XSHNEWT
        REAL                         :: YSHNEWT,XSHNEW,YSHNEW,IYDSQ
        INTEGER                      :: NLET,NXT,NYT,NZT
        INTEGER                      :: MAXIM,IMGEXP,IREF,IREFT,NPROJ
        INTEGER                      :: ISHRANGE,IT
        INTEGER                      :: IXCEN,IYCEN,IX,IY
        INTEGER                      :: N2X,N2Y,N2XLD,IRAD,IB
 	REAL                         :: AVRT,SIGRT,ROTPOST 

        LOGICAL, PARAMETER           :: USE_OMP   = .FALSE.
        LOGICAL, PARAMETER           :: WANTSTATS = .TRUE.
        LOGICAL, PARAMETER           :: MPIBCAST  = .TRUE.

        INTEGER                      :: NBORDER = 0     ! UNUSED
        INTEGER                      :: NSUBPIX = 0     ! UNUSED
        INTEGER                      :: lnblnkn,NLEN,IGO

        INTEGER, PARAMETER           :: LUNRING = 50
        INTEGER, PARAMETER           :: LUNFFT  = 51
        INTEGER, PARAMETER           :: LUNPIC  = 77
        INTEGER, PARAMETER           :: LUNANG  = 78

        REAL, PARAMETER              :: QUADPI = 3.1415926535897932384
        REAL, PARAMETER              :: DGR_TO_RAD = (QUADPI/180)         

        real                         :: xcen,ycen,rad

        integer                      :: kernel  = 3
        real                         :: dfactor,sig1,grad
        integer                      :: ndum,i,j,ii,jj,inv
        real, allocatable            :: bufwork1(:,:),bufwork2(:,:)
        real, allocatable            :: bufwork3(:,:)

        MYPID      = -1                 ! NOT USING MPI

        MAXRIN     = NUMR(3,NRING) - 2 ! ACTUAL LENGTH OF LONGEST RING

	RANGECOS   = COS(RANGE*DGR_TO_RAD)
        ONLYONEROT = (CTYPE == 'FOU1')

        N2X   = NX  * 2
        N2Y   = NY  * 2
        NXLD  = NX  + 2 - MOD(NX, 2)  ! FOURIER PAD
        N2XLD = N2X + 2 - MOD(N2X,2)  ! FOURIER PAD

C       SPIDER DOUBLE SIZED IMAGE CENTER
        F2XCEN  = N2X / 2 + 1
        F2YCEN  = N2Y / 2 + 1

C       SPIDER IMAGE CENTER
        IYCEN  = NY / 2 + 1
        IXCEN  = NX / 2 + 1
        !write(6,*) 'NX,NY,FXCEN,FYCEN:',NX,NY,FXCEN,FYCEN

C       FIND NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)
         
        IF (ROTFIRST) THEN
           IF (FBS_WANTED) THEN
            WRITE(NOUT,*)' ALIGNING INPUT IMAGES WITH FBS INTERPOLATION'
           ELSE
           WRITE(NOUT,*)' ALIGNING INPUT IMAGES WITH QUAD INTERPOLATION'
           ENDIF
        ENDIF

!           bufwork1 (nx,  ny), bufwork2 (nx,  ny), 
!     &           LMASKPAT(NX,   NY), 
 	ALLOCATE(BUFEXP  (NX,   NY), 
     &           BUFPATM (N2X,  N2Y), 

     &           bufwork1(N2XLD,N2Y), 

     &           cmaskpat(nx,   ny), 
     &           LMASKCC (NXLD, NY), 
     &           BUFEXP_P(N2XLD,N2Y), 
     &           BUFWORK (N2XLD,N2Y), 
     &           BUFREF_P(NXLD, NY), 
     &           CIRCEXP (LCIRC),  STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           MWANT = 3*NX*NY + N2X*N2Y + 2*N2XLD*N2Y + 2*NXLD*NY + LCIRC
           CALL ERRT(46,'APFOU_PATM; BUFEXP...',MWANT)
           GOTO 9999
        ENDIF

        IF (USE_FBS_INTERP) THEN
           ALLOCATE (F0 (NXLD,  NY),
     &               X1 (NXLD,  NY),
     &               Y1 (NXLD,  NY),
     &               XY2(NXLD,  NY),
     &               STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN 
              MWANT = 4*NXLD*NY 
              CALL ERRT(46,'APFOU_PATM; F0...',MWANT)
              GOTO 9999
           ENDIF  
        ENDIF

C       LOAD REF. PROJ. ANGLES FROM DOC. FILE  OR IMAGE HEAD (IF WANTED)
C       CONVERT REF. ANGLES TO UNITARY DIRECTIONAL VECTORS (REFDIR).
	CALL AP_GETANGAS(IREFLIST,NUMREF,0,REFANGDOC,REFPAT,
     &                   LUNPIC,LUNANG,3,ANGREF,GOTREFANG,NGOTREF,
     &                   .TRUE.,REFDIR,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       CREATE CIRCULAR SOFT  MASK FOR PATTERSON MAP MASKING
C       SOFT MASKING HELPS TINY BIT 
        XCEN = IXCEN
        YCEN = IYCEN
        RAD  = MIN(NX,NY) - MAX(IXCEN,IYCEN) - 2 
        CALL MAKE_CIRC_CSMASK(CMASKPAT, NX,NY,XCEN,YCEN, RAD)
        call chkfile('jnkcmask',66,1,nx,ny,1, cmaskpat,irtflg)


C       CREATE CIRCULAR LOGICAL MASK FOR IMAGE CC
        IRAD = MIN(NX,NY) - MAX(IXCEN,IYCEN) - 3
        CALL MAKE_CIRC_LMASK(LMASKCC, NXLD,NY,IXCEN,IYCEN,
     &                 IRAD, .FALSE.,.TRUE.)

        CALL FILERD(FFTPAT,NLET,NULL,'TEMP FILE PATTERN',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       READ REFERENCE IMAGES PATTERSON MAP INTO REFERENCE RINGS (CIRCREF) ARRAY OR
C       CREATE REFERENCE RINGS FILE FOR LATER READING 

        CALL APRINGS_NEW_PATM(IREFLIST,NUMREF,
     &               NX,NY,CMASKPAT,
     &               NRING,LCIRC,NUMR,MODE,FFTW_PLANS, 
     &               REFPAT,LUNPIC, CIRCREF,CIRCREF_IN_CORE,
     &               LUNFFT,LUNRING,SCRFILE,FFTPAT, 
     &               AVR,SIGR,IRTFLG) 
        IF (IRTFLG .NE. 0) GOTO 9999
        !write(6,*) ' Got ref proj PATM rings'
        !call chkring(6,.false., circref(1,1),lcirc, numr,nring, ndum,ndum)

        IF (NSAID <= 0) THEN
           IF (CIRCREF_IN_CORE) THEN
              WRITE(NOUT,91)NUMTH
91            FORMAT('  Ref. rings in core,  Threads: ',I4)
           ELSE
              WRITE(NOUT,92)NUMTH
92            FORMAT('  Ref. rings not in core,  Threads: ',I4)
           ENDIF
           NSAID = NSAID + 1
        ENDIF

        !write(6,'(a,(3i6,1x,f6.2)/)')' maxrin:',maxrin,nring,lcirc,divas
        !write(6,'(a,/,(3i6))') 'numr:', numr

C       LOAD EXP. ANGLES & ALIGNMENT PARAM. FROM DOC. FILE OR 
C       HEADER (IF WANTED) 
        CALL AP_GETANGAS(IEXPLIST,NUMEXP,0,EXPANGDOC,EXPPAT,
     &                   LUNPIC,LUNANG,8,ANGEXP,GOTEXPANG,NGOTPAR,
     &                   .FALSE.,FDUM,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        !write(6,*) ' Got getangas2'

C       INITIALIZE CC STATISTICS
        CALL AP_STAT_ADD(-1,CC,ANGDIF,ANGDIFTHR,CCLAS,
     &                   CCAVG,IBIGANGDIF,ANGDIFAVG,IMPROVCC,
     &                   CCIMPROV,IWORSECC,CCWORSE)

        PADVAL = 0.0  ! PADS OUTSIDE OF LOADED IMAGE
             
        anginheader = (rotfirst .and. expangdoc .eq. '*') ! UNFINISHED !!!!!!!!!!!

        DO IEXP=1,NUMEXP
C         LOOP OVER ALL EXPERIMENTAL (SAMPLE) IMAGES ------------------
    
C         CONVERT EXP. ANGLE TO UNITARY DIRECTIONAL VECTORS (EXPDIR).
	  CALL AP_GETSATA(ANGEXP(1,IEXP),EXPDIR,8,1,IRTFLG)

C         LOAD CURRENT EXP. IMAGE INTO ARRAY BUFEXP_P
          IF (ROTFIRST) THEN
C            WANT TO ROTATE/SHIFT EXP IMAGE WHEN READING IT

	     CALL AP_GETDATA_RTSQ(IEXPLIST,NUMEXP, 
     &                       NX,NY, N2XLD,N2Y,PADVAL,
     &                       1,EXPPAT,LUNPIC, IEXP,IEXP,
     &                       ANGINHEADER, ANGEXP, 
     &                       MPIBCAST,BUFEXP,BUFEXP_P,
     &                       WANTSTATS,AVI,SIGI,IRTFLG)
          ELSE
C            DO NOT ROTATE/SHIFT EXP IMAGES FIRST
             !write(6,*) ' Calling: AP_GETDATA FOR IMAGE: ',IEXP
  	     CALL AP_GETDATA(IEXPLIST,NUMEXP,
     &                       NX,NY, N2XLD,N2Y, PADVAL,
     &                       1,EXPPAT,LUNPIC, IEXP,IEXP,
     &                       MPIBCAST,BUFEXP_P,
     &                       WANTSTATS,AVI,SIGI,IRTFLG)
          ENDIF
          IF (IRTFLG .NE. 0) GOTO 9999

          !call chkfile('jnkexpimg',66,1,n2xld,n2y,1, bufexp_p,irtflg)

C         SAVE UNPADDED EXP IMAGE FOR ROTATION LATER
          BUFEXP(1:NX,1:NY) = BUFEXP_P(1:NX,1:NY)
          !call chkfile('jnkexpunpad',66,1,nx,ny,1, bufexp,irtflg)

C         MED or MS DENOISING THE EXP IMAGE USED FOR RINGS IS POORER

C         MAKE PATTERSON MAP OF EXP IMAGE (BUFEXP_P), PLACE IN BUFPATM
C         SIMPLE PWS GIVES WORSE RESULTS THAN PATM
 	  CALL MAKE_PATM(AVI,SIGI,BUFEXP_P, CMASKPAT,
     &                   NX,NY, N2XLD,N2X,N2Y,
     &                   BUFWORK,BUFPATM, IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9999
          !call chkfile('jnkexppatm',66,1,n2x,n2y,1, bufpatm,irtflg)
          !call chkring(6,.false., circref(1,1),lcirc, numr,nring, ndum,ndum)
          !write(6,*) ' Calling: APRINGS_ONE FOR IMAGE: ',IEXP

C         EXTRACT POLAR COORD. RINGS, FROM PATTERSON MAP NORMALIZE & FFT THEM
	  CALL APRINGS_ONE_NEW(N2X,N2Y,  F2XCEN,F2YCEN, BUFPATM,.FALSE.,
     &                         MODE,NUMR,NRING,LCIRC, 0.0,FFTW_PLANS,
     &                         CIRCEXP,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9999
          !write(6,*) 'circexp(1):',circexp(1),bufpatm(30,30),bufpatm(60,60)
          !call chkring(6,.false., circexp,lcirc, numr,nring, ndum,ndum)

C         DETERMINE WHICH REF IMAGES ARE TO BE COMPARED
          NPROJ      = NUMREF
          LIMITRANGE = (RANGECOS < 1.0)

C         IF LIMITRANGE, LIST NEARBY REF IMAGES, RETURNS: NPROJ 
          CALL MAKE_CLOSE_LIST(NUMREF,LIMITRANGE,
     &                         REFDIR,EXPDIR,
     &                         RANGECOS, .TRUE., 
     &                         LCG, NPROJ, IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9999

          IF (NPROJ <= 0) THEN
C            NO REF. IMAGE WITHIN COMPARISON ANGLE
             IREF  = 0 
             GOTO 1000
          ENDIF

#ifdef NEVER                                                                  
           limitrange = .true.
           nproj      = 4                                                      
           lcg(1)     = 3339  ! best fou
           lcg(2)     = 5051  ! best ref
           lcg(3)     = 4873  ! best fouless
           lcg(4)     = -4778 ! best fouless 1x
#endif 
                                                                       
C         ALL OR SOME REF. IMAGES FOUND WITHIN COMPARISON RANGE

          IF (CIRCREF_IN_CORE) THEN
C            USE CIRCREF FOR REFERENCE RINGS

c$omp        parallel do private(imil,imi,use_un,use_mir)
             DO IMIL=1,NPROJ        
                IMI = IMIL          ! REF INDEX (IF NOT LIMITED)
             
C               SEARCH BOTH MIRRORED & NON-MIRRORED IF CHKMIR
                USE_UN  = .TRUE.    ! NEEDED FOR OMP ???
                USE_MIR = CKMIRROR  ! NEEDED FOR OMP ???

                IF (LIMITRANGE) THEN
C                  ONLY CHECKING NEARBY PROJECTIONS
	           IMI = ABS(LCG(IMIL)) ! ACTUAL REF INDEX WHEN LIMITED

                   IF (CKMIRROR) THEN
C                     ONLY SEARCH EITHER MIRRORED OR NON-MIRRORED
                      USE_UN  = (LCG(IMIL) >= 0)
                      USE_MIR = (LCG(IMIL) < 0)
                   ENDIF
                ENDIF

C               CHECK MIRRORED/ NON-MIRRORED POSITION 
                CALL CROSRNG_2(CIRCREF(1,IMI),CIRCEXP,
     &                         LCIRC,NRING, MAXRIN,NUMR,
     &                         .FALSE.,FFTW_PLANS(1),
     &                         USE_UN,USE_MIR,
     &                         ISMIRRORED(IMIL),
     &                         DCCROT(IMIL),
     &                         ROTPOS(IMIL))

#ifdef NEVER
                write(6,'(a,i5,a,f8.1,a,f7.2,a,f9.2,a,l)') 
     &            '  Imi:',    imi, 
     &            '  Rotpos:', rotpos(imil),
     &            '  Ang:',    (rotpos(imil)-1)/maxrin*divas,
     &            '  CCrot:',  dccrot(imil),
     &            '  Mir:',    ismirrored(imil)
#endif
	     ENDDO
c$omp        end parallel do
          ELSE
             stop
          ENDIF

C         LOOP OVER ALL RELEVANT REF. IMAGES -------------------------
          PEAKV  = -HUGE(PEAKV)

C         OPEN  REF IMAGE FFT TO START THIS STACK (KEEP OPEN)
          IREF = IREFLIST(1)  ! ACTUAL REF IMAGE #
          CALL FILGET(FFTPAT,FILNAM,0,IREF,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9999
          CALL OPFILEC(0,.FALSE.,FILNAM,LUNPIC,'O',ITYPET,
     &                    NXLDT,NYT,NZT,MAXIM,' ',.TRUE.,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9999
          CALL SIZCHK(NULL,NXLD,NY,1,0, NXLDT,NYT,1,0,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9999

          NPROJT = NPROJ
          IF (ONLYONEROT) THEN
C            ONLY EVALUATE THE REF WITH TOP CCROT VALUE
             NPROJT = 1
             ILOC   = MAXLOC(DCCROT)
             LCG(1) = ILOC(1)
          ENDIF

          DO IMIL=1,NPROJT
             IREFT = IREFLIST(IMIL)                ! ACTUAL REF IMAGE #
	     IF (LIMITRANGE) IREFT = ABS(LCG(IMIL))

             IMI   = IMIL
             IF (ONLYONEROT) THEN
                IREFT = ABS(LCG(IMIL))
                IMI   = IREFT
             ENDIF
             ROTANGNEWT = (ROTPOS(IMI)-1)/MAXRIN*DIVAS ! ROT. ANGLE
             MIRRORNEWT = ISMIRRORED(IMI)              ! MIRROR FLAG
             CCROTT     = DCCROT(IMI)                  ! ROT. CC
             ROTPOST    = ROTPOS(IMI)
             AVRT       = AVR(IREFT)
             SIGRT      = SIGR(IREFT)

C            OPEN REFERENCE IMAGE FFT FILE
             CALL GETOLDIMG(LUNPIC,0,FFTPAT,IREFT,
     &                      VERBOSE,.TRUE.,NGOT,IRTFLG)
             IF (IRTFLG .NE. 0) GOTO 9999
             IF (NGOT .NE. IREFT) THEN
                CALL ERRT(102,'INCORRECT IMAGE NUMBER',NGOT)
                GOTO 9999
             ENDIF

C            READ FFT
             CALL REDVOL(LUNPIC,NXLD,NY,1,1,BUFREF_P,IRTFLG)
             IF (IRTFLG .NE. 0) GOTO 9999

C            ROTATE AND SHIFT THE EXP IMAGE TO MATCH THIS REF.
C            CHECK BOTH ROTANGNEW AND: ROTANGNEW + 180.0
C            RETURNS SUB-PIXEL SHIFTS.

             !call chkfile('jnkexp',66,1,nx,ny,1, bufexp,irtflg)
             !write(6,*)' Calling apshift_fou: ',ireft
             CALL APSHIFT_FOU(LUNPIC, IREFT,
     &                  NX,NY, NXLD,
     &                  BUFEXP,    AVI,SIGI, AVRT,SIGRT, LMASKCC,
     &                  BUFEXP_P, BUFREF_P, F0,X1,Y1,XY2, 
     &                  ISHRANGEX,ISHRANGEY,
     &                  ROTANGNEWT,XSHNEWT,YSHNEWT,
     &                  MIRRORNEWT,PEAKVT,IRTFLG)
             IF (IRTFLG .NE. 0) GOTO 9999

#ifdef NEVER
             !call chkfile('jnkexp_p',66,1,nxld,ny,1, bufexp_p,irtflg)
             !call chkfile('jnkrefgud',66,1,nxld,ny,1, bufref_p,irtflg)
             !write(6,*)' avi,sigr,peak:',avi,sigr(imi),peakvt
             !ref   rotpos  ang     ccrot   shifts           peak  mir
             write(6,
     &     '(a,i4, a,f9.1, a,f7.1, a,f9.2, a,f5.1,a,f5.1,a, af9.4,a,L)')
     &           '  Ref:',       ireft,   
     &           '  Rotpos:',    rotpost,
     &           '  Ang:',       rotangnewt,
     &           '  CCrot:',     ccrott,
     &           '  Shift:(',    xshnewt, 
     &           ',',            yshnewt,')',
     &           '  CC:',        peakvt, 
     &           '  Mir:',       mirrornewt 
#endif
             IF (PEAKVT > PEAKV) THEN
C               GOOD MATCH WITH TOTMIN (MIRRORED OR NOT)  POSITION 
                CCROT     = CCROTT
                CC        = PEAKVT
                PEAKV     = PEAKVT
                ROTANGNEW = ROTANGNEWT
                MIRRORNEW = MIRRORNEWT
                XSHNEW    = XSHNEWT
                YSHNEW    = YSHNEWT
                IREF      = IREFT           ! ACTUAL REF IMAGE #
	     ENDIF
          ENDDO   ! END OF: DO IMIL=1,IEND -------------------------



1000      IMGEXP = IEXPLIST(IEXP)          ! ACTUAL EXP IMAGE #

          IF (IREFT <= 0) THEN
C            NO NEARBY REFERENCE IMAGE
C            IREF IS FOR REFDIR INDEX
             IREF      = 1
             CCROT     = 0.0
             PEAKV     = 0.0
             ROTANGNEW = 0.0
             XSHNEW    = 0.0
             YSHNEW    = 0.0
             MIRRORNEW = .FALSE.
          ENDIF

          !write(6,*) '  iref: ',iref,ireft,iexp
          IF (MOD(NX,2) > 0) XSHNEW = XSHNEW + 1
          IF (MOD(NY,2) > 0) YSHNEW = YSHNEW + 1

C         AP_END WRITES ALIGNMENT PARAMETERS TO DOC FILE 

          ISHRANGE = MAX(ISHRANGEX,ISHRANGEY)
          CALL AP_ENDS(IEXP,IMGEXP,IREF, 
     &         ANGREF(1,IREF),REFDIR(1,IREF),
     &         ANGEXP(1,IEXP), EXPDIR,ISHRANGE,
     &         GOTREFANG, NGOTPAR, CCROT,PEAKV,
     &         ROTANGNEW,XSHNEW,YSHNEW,MIRRORNEW,REFPAT,
     &         NPROJ,CTYPE, LUNDOC,.FALSE.,.FALSE.,PARLIST)

C         WRITE DATA TO IMAGE HEADER 
          CALL AP_END_HEAD(IMGEXP,EXPPAT,LUNPIC,PARLIST,8,IRTFLG)

          CALL  AP_STAT_ADD(NGOTPAR,CC,PARLIST(10),
     &                    ANGDIFTHR,ANGEXP(8,IEXP),
     &                    CCAVG,IBIGANGDIF,ANGDIFAVG,IMPROVCC,
     &                    CCIMPROV,IWORSECC,CCWORSE)

c#ifdef NEVER
            !exp   ref  ang     ccrot   shifts           cc   mir
             write(6,
C dgm     &     '(a,i4, a,i9, a,f7.1,a,f8.2, a,f5.1,a,f5.1,a, af8.4, a,L)')
     &     '(a,i4, a,i9, a,f7.1,a,f8.2, a,f5.1,a,f5.1,a,f8.4, a,L)')
     &           '  Exp:',       iexplist(iexp),
     &           '  Ref:',       iref,   
     &           '     Ang:',    rotangnew,
     &           '  CCrot:',     ccrot,
     &           '  Shift:(',    xshnew, 
     &           ',',            yshnew,')',
     &           '  CC:',        cc, 
     &           '  Mir:',       mirrornew 

c#endif
          
       ENDDO  ! END OF: DO IEXP=1,NUMEXP  LOOP OVER ALL EXP. IMAGES ----


       IF (NUMEXP > 1) THEN
C         SAVE CC & ANGULAR DISPLACEMENT STATISTICS
          CALL AP_STAT(NUMEXP,ANGDIFTHR,IBIGANGDIF,
     &                 ANGDIFAVG, CCAVG,
     &                 IMPROVCC,CCIMPROV,
     &                 IWORSECC,CCWORSE,
     &                 NBORDER,NSUBPIX,LUNDOC)
       ENDIF

9999   CLOSE(LUNANG)
       CLOSE(LUNDOC)
       IF (.NOT. CIRCREF_IN_CORE) CLOSE(LUNRING)

       IF (ALLOCATED(CMASKPAT)) DEALLOCATE(CMASKPAT)
       IF (ALLOCATED(LMASKCC))  DEALLOCATE(LMASKCC)
       IF (ALLOCATED(BUFWORK))  DEALLOCATE(BUFWORK)
       IF (ALLOCATED(BUFEXP_P)) DEALLOCATE(BUFEXP_P)
       IF (ALLOCATED(BUFREF_P)) DEALLOCATE(BUFREF_P)
       IF (ALLOCATED(BUFEXP))   DEALLOCATE(BUFEXP)
       IF (ALLOCATED(BUFPATM))  DEALLOCATE(BUFPATM)
       IF (ALLOCATED(F0))       DEALLOCATE(F0)
       IF (ALLOCATED(X1))       DEALLOCATE(X1)
       IF (ALLOCATED(Y1))       DEALLOCATE(Y1)
       IF (ALLOCATED(XY2))      DEALLOCATE(XY2)
       IF (ALLOCATED(CIRCEXP))  DEALLOCATE(CIRCEXP)
       IF (ASSOCIATED(LCG))     DEALLOCATE(LCG)
       NULLIFY(LCG)

       END
