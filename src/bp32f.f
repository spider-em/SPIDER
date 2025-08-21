
C++*********************************************************************
C
C  BP32F.F                                         JAN 00 P.A. Penczek
C              ANGLES IN HEADER                    JUL 01 ArDean Leith
C              OPFILE 'U'                          FEB 02 ArDean Leith
C              BESSELS                             MAY 02 P.A. Penczek
C              OPFILEC                             FEB 03 ArDean Leith
C              BUILDM PARAMETERS                   JUL 03 ArDean Leith
C              BUILDM PARAMETERS                   SEP 03 ArDean Leith
C              MPI                                 OCT 03 Chao Yang
C              REWRITE                             DEC 06 ArDean Leith
C              MPI CHANGES                         OCT 06 ArDean Leith
C              MPI READV1P BUG                     DEC 09 ArDean Leith
C              OPFILES PARAMETERS                  DEC 10 ArDean Leith
C              PREVIOUSLY NAMED WIW32**            JAN 11 ArDean Leith
C              MPI HAD UNDEFINED NY BUG            MAR 11 ArDean Leith
C              ROT2QS --> RTSQ RENAMED             DEC 11 ArDean Leith
C              NSAM --> NX, RTSQ PARAM.            JAN 12 ArDean Leith
C              ILIST ALLOC                         JAN 13 ArDean Leith
C              COMMENTS & COSMETIC                 FEB 13 ArDean Leith
C
C=**********************************************************************
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright (C)2002,2013 P. A. Penczek & ArDean Leith                *
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
C  BP32F(CALLRTSQ)
C
C  PURPOSE:  FOR: 'BP 32F' IMPROVED WITH LESS MEMORY REQUIRED AND
C            MORE STACK USAGE.
C
C  CONTAINS: WINDKB2A
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE BP32F(CALLRTSQ)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'F90ALLOC.INC'

C       DOC FILE POINTERS
        REAL, POINTER          :: ANGBUF(:,:), ANGSYM(:,:)

        REAL,    ALLOCATABLE   :: DM(:,:),  SM(:,:) 
        REAL,    ALLOCATABLE   :: WE(:,:,:),WO(:,:,:)
        COMPLEX, ALLOCATABLE   :: XE(:,:,:),XO(:,:,:)
        REAL,    ALLOCATABLE   :: TEMP(:)
        REAL,    ALLOCATABLE   :: ROTANG(:),SHX(:),SHY(:)
        INTEGER, ALLOCATABLE   :: ILIST(:)

        LOGICAL                :: ANGINDOC,CALLRTSQ
        CHARACTER(LEN=1)       :: NULL = CHAR(0)
        CHARACTER(LEN=MAXNAM)  :: FILPAT 
        CHARACTER(LEN=MAXNAM)  :: EVENVOL,ODDVOL,BOTHVOL
        CHARACTER(LEN=MAXNAM)  :: ANGDOC

C       COMMON: /TABS/ IS USED IN ONELINE, EXTRACTLINE, PUTLINE3, ETC
        INTEGER, PARAMETER     :: LTAB = 4999
        COMMON  /TABS/ LN2,FLTB,TABI(0:LTAB)

        INTEGER, PARAMETER     :: IOPIC   = 18
        INTEGER, PARAMETER     :: INPROJ  = 19
        INTEGER, PARAMETER     :: LUNDOC  = 80
        INTEGER, PARAMETER     :: LUNXM1  = 0    ! UNSABLE NEED #'s
        !INTEGER, PARAMETER    :: LUNXM2  = 83   ! IN CALLED SUB
        !INTEGER, PARAMETER    :: LUNDOC  = 80   ! IN CALLED SUB
        !INTEGER, PARAMETER    :: LUNROTT = 81   ! IN CALLED SUB 
        INTEGER                :: IMGNUM

        CALL SET_MPI(ICOMM,MYPID,MPIERR)      ! SETS ICOMM AND MYPID

        NILMAX = NIMAXPLUS
        ALLOCATE(ILIST(NILMAX),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'BP32F; ILIST....',NILMAX)
           RETURN
        ENDIF 

C       OPEN INPUT IMAGE FILES 
        IMGNUM = 0    ! Passing uninitialized variable may result in 'INVALID IMAGE NUMBER' error
        CALL OPFILES(0,INPROJ,LUNDOC,LUNXM1,
     &             .TRUE.,FILPAT,NLET, 'O',
     &             ITYPE,NX,NY,NZ,MAXIM1,
     &             'TEMPLATE FOR IMAGE FILES (E.G. STK@****)~~',
     &             .FALSE., ILIST,NILMAX, 
     &             NDUM,NANG,IMGNUM, IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        IF (NX .NE. NY) THEN
           CALL ERRT(101,'ONLY WORKS ON SQUARE IMAGES',NDUM)
           GOTO 999
        ENDIF

        MAXNUM = MAXVAL(ILIST(1:NANG))

C       NANG - TOTAL NUMBER OF IMAGES
        IF (MYPID <= 0) WRITE(NOUT,2001) NANG
2001    FORMAT('  NUMBER OF IMAGES: ',I7)

C       RETRIEVE ARRAY WITH ANGLES DATA IN IT
        ANGINDOC = .TRUE.

C       PSI, THE, PHI, REF#, EXP#, INPLANE, SX, SY  
        MAXXT    = 8 + 1
        MAXYT    = MAXNUM

        CALL GETDOCDAT('ANGLES DOC',.TRUE.,ANGDOC,77,.FALSE.,MAXXT,
     &                 MAXYT,ANGBUF,IRTFLG)
        IF (IRTFLG .NE. 0) ANGINDOC = .FALSE.
        !write(6,*)'gotdocdat, 9999:',angbuf(1,9999),angbuf(1,10000)
                
C       RETRIEVE ARRAY WITH SYMMETRIES DATA IN IT
        MAXXS  = 0
        MAXSYM = 0
        CALL GETDOCDAT('SYMMETRIES DOC',.TRUE.,ANGDOC,77,.TRUE.,MAXXS,
     &                   MAXSYM,ANGSYM,IRTFLG)
        IF (IRTFLG .NE. 0)  MAXSYM = 1
        

        N2   = 2 * NX
        LSD  = N2 + 2 - MOD(N2,2)
        NMAT = LSD * N2 * N2
       
        NDIMANG = 1
        IF (ANGINDOC) NDIMANG = MAXNUM
        NDIMSYM = MAX(1,MAXSYM)

        ALLOCATE(DM(9,NDIMANG), 
     &           SM(9,NDIMSYM), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           MWANT = 9*NDIMANG + 9*NDIMSYM 
           CALL ERRT(46,'BP 32F; DM, SM', MWANT)
           GOTO 999
        ENDIF

        IF (ANGINDOC) THEN
C          GET ANGLES FROM DOCUMENT FILE
           CALL BUILDM1(ILIST,NANG, DM,ANGBUF,9,MAXYT,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 999

C          NO LONGER NEED ANGBUF        
           IF (.NOT. CALLRTSQ) DEALLOCATE(ANGBUF)
        ENDIF

        !write(6,*)'gotdocdat,buildm1:',irtflg,dm(1,5)
        CALL FLUSHRESULTS()

        IF (MAXSYM .GT. 1)  THEN
C          HAVE SYMMETRIES
           CALL BUILDS(SM,MAXSYM,ANGSYM(1,1),IRTFLG)
           DEALLOCATE(ANGSYM)
        ENDIF

        ALLOCATE(XE(0:NX, N2,N2),
     &           WE(0:NX, N2,N2),
     &           XO(0:NX, N2,N2), 
     &           WO(0:NX, N2,N2), STAT=IRTFLG)

        IF (IRTFLG.NE.0) THEN 
C          SOME ARRAYS ARE COMPLEX SO 6 NOT 4
           MWANT = 6 * (NX+1) * N2 * N2
           CALL ERRT(46,'BP 32F; XE, WE, WO, & XO',MWANT)
           GOTO 999
        ENDIF

C       CREATE FFTW3 PLAN FOR REVERSE 3D FFT ON XE USING ALL THREADS
        !write(6,*) ' bp32f;  inv plan on xe: ',n2
        CALL FMRS_PLAN(.TRUE.,XE,N2,N2,N2, 0,-1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
 
C       GENERALIZED KAISER-BESSEL WINDOW ACCORDING TO LEWITT
        LN    = 5
        LN2   = LN / 2                    ! USED IN COMMON /TABS/
        R     = NX / 2
        V     = REAL(LN-1) / 2.0 / REAL(N2)
        ALPHA = 6.5
        AAAA  = 0.9*V
        NNN   = 3

C       GENERATE TABLE WITH INTERPOLANTS
        B0   = SQRT(ALPHA) * BESI1(ALPHA)
        FLTB = REAL(LTAB) / REAL(LN2+1)

cc$omp  parallel do private(i,s,x),shared(mmm)
        DO  I=0,LTAB
           S = REAL(I) / FLTB / N2
           IF (S <= AAAA)  THEN
              X       = SQRT(1.0 - (S/AAAA)**2)
              TABI(I) = SQRT(ALPHA*X) * BESI1(ALPHA*X) / B0
           ELSE
              TABI(I) = 0.0
           ENDIF
        ENDDO

C       BP32FQ(NS,XE,WE,XO,WO,LSD,N,N2, 
C       BP32FQ(NX,XE,WE,XO,WO,LSD,N

        CALL BP32FQ(NX, XE,WE,XO,WO,LSD,N2, CALLRTSQ,FILPAT,
     &               INPROJ,ANGBUF,ILIST,NANG, 
     &               DM,IMGNUM,SM,MAXSYM,ANGINDOC, 
     &               LUNXM1,MAXIM1,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        CALL FILERD(BOTHVOL,NLETI,NULL,'RECONSTRUCTED VOLUME',IRTFLG) 
        IF (IRTFLG .NE. 0) GOTO 999

        CALL FILERD(EVENVOL,NLETI,NULL,'FIRST SAMPLE VOLUME',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999
 
        CALL FILERD(ODDVOL,NLETI,NULL,'SECOND SAMPLE VOLUME',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

C       CBE AND CBO (XE & XO) ARE LSD x N2 x N2
C       CWE AND CWO (WE & WO) ARE (LSD/2) x N2 x N2
        LSDD2 = LSD / 2

        !write(6,*) ' Opening bothvol'
        !call flushresults

C       STORE FOURIER XE TEMPORARILY IN OVERALL VOLUME: BOTHVOL
        MAXIM = 0
        ITYPE = 3
        CALL OPFILEC(0,.FALSE.,BOTHVOL,IOPIC,'U',ITYPE,LSD,N2,N2,
     &             MAXIM,'DUMMY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999
        !write(6,*) ' opened bothvol; irtflg',irtflg

        CALL WRITEV(IOPIC,XE,LSD,N2,LSD,N2,N2)
        IF (MYPID <= 0) CLOSE(IOPIC)

        !write(6,*) ' xe  cached in: ',bothvol
        !call flushresults

C       STORE WEIGHT WE TEMPORARILY IN ODDVOL (ODD VOLUME)
        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,ODDVOL,IOPIC,'U',ITYPE,LSDD2,N2,N2,
     &             MAXIM,'DUMMY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 999

        CALL WRITEV(IOPIC,WE,LSDD2,N2,LSDD2,N2,N2)
        IF (MYPID <= 0) CLOSE(IOPIC)

        !write(6,*) ' we cached in:',oddvol
        !call flushresults

C       PROCESS EVEN VOLUME -------------------------------- EVEN
c        write(6,*) 'vals: even',xe(1,1,1),xe(1,1,1),
c    &                           we(1,1,1),we(2,2,2),
c    &                           xo(2,2,2),xo(2,2,2),
c    &                           wo(2,2,2),wo(2,2,2)
 
C       WEIGHT, FOURIER TRANSFORM, & WINDOW USING: XE AND WE
        !write(6,*) 'BP32F;  nrmw2 xe : ',NX,n2
        CALL NRMW2(XE,WE,NX,N2)
        CALL WINDKB2A(XE,XE,NX,LSD,N2,ALPHA,AAAA,NNN)

C       NOW EVEN VOLUME IS READY, SYMMETRIZE IF NECESSARY
        IF (MAXSYM .GT. 1)  THEN
C          ADDITIONAL SYMMETRIZATION OF VOLUME XE IN REAL SPACE 05/03/02
           ALLOCATE (TEMP(NX*NX*NX), STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN 
              MEMWANT = NX * NX * NX 
              CALL ERRT(46,'BP 32F; TEMP',MEMWANT)
              GOTO 999
           ENDIF
           CALL COP(XE,TEMP,NX*NX*NX)
c$omp      parallel do private(i,j,k)
           DO K=1,N2
              DO J=1,N2
                 DO I=0,NX
                    XE(I,J,K) = CMPLX(0.0,0.0)
                 ENDDO
              ENDDO
           ENDDO
           IF (MOD(NX,2) .EQ. 0)  THEN
              KNX = NX/2-1
           ELSE
              KNX = NX/2
           ENDIF
           KLX = -NX/2
           CALL SYMVOL(TEMP,XE,KLX,KNX,KLX,KNX,KLX,KNX,SM,MAXSYM)
        ENDIF

C       WRITE EVEN VOLUME TO FILE: EVENVOL,  NOTE: NX=NY=NZ 
        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,EVENVOL,IOPIC,'U',ITYPE,NX,NX,NX,
     &             MAXIM,'DUMMY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999
        CALL WRTVOL(IOPIC,NX,NX, 1,NX, XE,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999
        IF (MYPID <= 0) CLOSE(IOPIC)

        !write(6,*) ' even stored in:',evenvol
        !call flushresults

C       RECOVER XE FROM BOTHVOL AND PUT IT IN: XE --------------- BOTH 
        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,BOTHVOL,IOPIC,'O',ITYPE,LSD,N2,N2,
     &               MAXIM,'DUMMY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999
        CALL READV(IOPIC,XE,LSD,N2,LSD,N2,N2)
        CLOSE(IOPIC)

        !write(6,*) ' xe recovered'
        !call flushresults

C       RECOVER WEIGHT WE  TEMP. STORED IN: ODDVOL AND PUT IT IN: WE 
        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,ODDVOL,IOPIC,'O',ITYPE,LSDD2,N2,N2,
     &             MAXIM,'DUMMY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999
        CALL READV(IOPIC,WE,LSDD2,N2,LSDD2,N2,N2)
        CLOSE(IOPIC)

        !write(6,*) ' we recovered'
        !call flushresults

        !write(6,*) 'vals: rec ',xe(1,1,1),xe(1,1,1),
c    &                           we(1,1,1),we(2,2,2),
c    &                           xo(2,2,2),xo(2,2,2),
c    &                           wo(2,2,2),wo(2,2,2)

C       ADD EVEN AND ODD VOLUMES, PLACE IN XE AND WE
        CALL ADDADA(XE,XO,NMAT)
        CALL ADDADA(WE,WO,NMAT/2)

C       WEIGHT, FOURIER TRANSFORM, AND WINDOW BOTHVOL USING: XE AND WE
        !write(6,*) ' BP32F;  nrmw2 xe : ',NX,n2
        CALL NRMW2(XE,WE,NX,N2)
        CALL WINDKB2A(XE,XE,NX,LSD,N2,ALPHA,AAAA,NNN)

C       OVEALL VOLUME IS NOW IN: XE
        IF (MAXSYM .GT. 1)  THEN
C          ADDITIONAL SYMMETRIZATION OF VOLUME TOTAL IN REAL SPACE 05/03/02
           CALL COP(XE, TEMP,NX*NX*NX)
c$omp      parallel do private(i,j,k)
           DO K=1,N2
              DO J=1,N2
                 DO I=0,NX
                    XE(I,J,K) = CMPLX(0.0,0.0)
                 ENDDO
              ENDDO
           ENDDO
           CALL SYMVOL(TEMP,XE,KLX,KNX,KLX,KNX,KLX,KNX,SM,MAXSYM)
        ENDIF

C       WRITE OVERALL VOLUME FROM: XE  TO: BOTHVOL
        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,BOTHVOL,IOPIC,'U',ITYPE,NX,NX,NX,
     &             MAXIM,'DUMMY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 999
        CALL WRTVOL(IOPIC,NX,NX, 1,NX, XE,IRTFLG)
        IF (MYPID <= 0) CLOSE(IOPIC)

        !write(6,*) 'vals: end ',xe(1,1,1),xe(1,1,1),
c    &                          we(1,1,1),we(2,2,2),
c    &                          xo(2,2,2),xo(2,2,2),
c    &                          wo(2,2,2),wo(2,2,2)
        !write(6,*) ' overall stored'
        !call flushresults

C       ODD VOLUME IS STILL UNDAMAGED IN: XO --------------------- ODD

C       WEIGHT, FOURIER TRANSFORM, AND WINDOW ODD VOL. USING: XO AND WO
        !write(6,*) ' BP32F;  nrmw2 xo : ',NX,n2
        CALL NRMW2(XO,WO,NX,N2)
        CALL WINDKB2A(XO,XO,NX,LSD,N2,ALPHA,AAAA,NNN)

        IF (MAXSYM .GT. 1)  THEN
C          ADDITIONAL SYMMETRIZATION OF VOLUME XO IN REAL SPACE 05/03/02
           CALL COP(XO,TEMP, NX*NX*NX)
c$omp      parallel do private(i,j,k)
           DO K=1,N2
              DO J=1,N2
                 DO I=0,NX
                    XO(I,J,K) = CMPLX(0.0,0.0)
                 ENDDO
              ENDDO
           ENDDO
           CALL SYMVOL(TEMP,XO,KLX,KNX,KLX,KNX,KLX,KNX,SM,MAXSYM)
        ENDIF

C       WRITE ODD FILE VOLUME FROM: XO  TO: ODDVOL
        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,ODDVOL,IOPIC,'U',ITYPE,NX,NX,NX,
     &              MAXIM,'DUMMY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999
        CALL WRTVOL(IOPIC,NX,NX, 1,NX, XO,IRTFLG)
        IF (MYPID <= 0) CLOSE(IOPIC)

        !write(6,*) ' odd stored'

999     CONTINUE
        IF (ALLOCATED(ILIST))   DEALLOCATE(ILIST)
        IF (ALLOCATED(DM))      DEALLOCATE(DM)
        IF (ALLOCATED(SM))      DEALLOCATE(SM)
        IF (ALLOCATED(XE))      DEALLOCATE(XE)
        IF (ALLOCATED(WE))      DEALLOCATE(WE)
        IF (ALLOCATED(XO))      DEALLOCATE(XO)
        IF (ALLOCATED(WO))      DEALLOCATE(WO)
        IF (ALLOCATED(TEMP))    DEALLOCATE(TEMP)
        IF (ASSOCIATED(ANGBUF)) DEALLOCATE(ANGBUF)
        IF (ASSOCIATED(ANGSYM)) DEALLOCATE(ANGSYM)

        CLOSE(INPROJ)
        CLOSE(LUNXM1)

        END

C       ---------------- BP32FQ -----------------------------------

        SUBROUTINE BP32FQ(NX,XE,WE,XO,WO,LSD,N, CALLRTSQ,FILPAT,
     &                     INPROJ,ANGBUF,INUMBRT,NANG, 
     &                     DM,IMGNUM,SM,MAXSYM,ANGINDOC,
     &                     LUNXM1,MAXIM1,IRTFLG)

C       NOTE: STUPID TRANSFORM OF N2-->N  !!!!al

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL                  :: DM(3,3,NANG),SM(3,3,MAXSYM),DMS(3,3)

        REAL                  :: ANGBUF(9,NANG)
        LOGICAL               :: ANGINDOC,CALLRTSQ
        CHARACTER(LEN=*)      :: FILPAT
        LOGICAL, ALLOCATABLE  :: RANDLIST(:)
        REAL,    ALLOCATABLE  :: PROJ(:,:), PROJTEMP(:,:)
        COMPLEX, ALLOCATABLE  :: BI(:,:)
        REAL                  :: ANGBUFT(4)
        INTEGER               :: INUMBRT(NANG)
        REAL                  :: WE(0:NX,N,N)
        REAL                  :: WO(0:NX,N,N)
        COMPLEX               :: XE(0:NX,N,N)
        COMPLEX               :: XO(0:NX,N,N)

        CHARACTER(LEN=MAXNAM) :: FILNAM,ODDVOL,FILPATOUT
        INTEGER, ALLOCATABLE  :: INUMBROUT(:)

        DOUBLE PRECISION      :: PI
        LOGICAL               :: ITMP

        INTEGER, PARAMETER    :: LUNDOC  = 80
        INTEGER, PARAMETER    :: LUNROTT = 81 
        INTEGER, PARAMETER    :: LUNXM2  = 81 

#ifdef USE_MPI
c       This MPI version is memory intensive.
c       It requires 6 copies of the 3-D volume,
c       2-D images are read into memory and distributed.
c       Each processor will hold roughly nang/nproc 2-D images.  

        INCLUDE 'mpif.h'
        INTEGER               :: ISTAT(MPI_STATUS_SIZE)
        INTEGER, ALLOCATABLE  :: PSIZE(:), NBASE(:)
                          
        REAL   , ALLOCATABLE  :: PRJLOC(:,:,:)
        REAL   , ALLOCATABLE  :: PRJBUF(:,:,:)
        REAL   , ALLOCATABLE  :: WELOC(:,:,:), WOLOC(:,:,:)
        COMPLEX, ALLOCATABLE  :: XELOC(:,:,:), XOLOC(:,:,:)

        ICOMM = MPI_COMM_WORLD
        CALL MPI_COMM_RANK(ICOMM, MYPID,  MPIERR)
        CALL MPI_COMM_SIZE(ICOMM, NPROCS, MPIERR)
#else
        MYPID = -1
#endif

c$omp   parallel do private(i,j,k)
        DO K=1,N
           DO J=1,N
              DO I=0,NX
                 XE(I,J,K) = CMPLX(0.0,0.0)
                 WE(I,J,K) = 0.0
                 XO(I,J,K) = CMPLX(0.0,0.0)
                 WO(I,J,K) = 0.0
              ENDDO
           ENDDO
        ENDDO

        NSIZE  = 1
        LUNROT = 0    ! IF ZERO --> NO RTSQ OUTPUT WANTED
        NXLD   = NX

        IF (CALLRTSQ) THEN
           IF (USE_FBS_INTERP) NXLD = NX + 2 - MOD(NX,2)

C          READ IMAGE AND ROTATE, SCALE & SHIFT INTO: PROJ
           NSIZE = NX
           ALLOCATE(INUMBROUT(NANG), STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN 
              CALL ERRT(46,'BP32FQ; INUMBROUT', NANG)
              GOTO 9999
           ENDIF

C          OPEN OUTPUT IMAGE(S) FOR RTSQ 
           ITYPE = 1
           CALL OPFILES(INPROJ,LUNROTT,LUNDOC,LUNXM2,
     &            .TRUE.,FILPATOUT,NLETO, 'U',
     &            ITYPE,NX,NX,1, MAXIM2,
     &            'TRANSFORMED OUTPUT IMAGES TEMPLATE (E.G. ROT@****)~',
     &            .FALSE., INUMBROUT,NANG, 
     &            NDUM,NANGT, IMGNUMOUT, IRTFLG) 

C          IF IRTFLG IS NEGATIVE THEN NO RTSQ OUTPUT WANTED
           IF (IRTFLG .EQ. 0) LUNROT = LUNROTT
           IF (IRTFLG .GT. 0) GOTO 9999

           IF (USE_FBS_INTERP) NXLD = NX + 2 - MOD(NX,2)
        ENDIF

        ALLOCATE (PROJ(NX,NX),
     &            PROJTEMP(NX,NSIZE),
     &            BI(0:NX,N),
     &            RANDLIST(NANG), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           MEMWANT = NX*NX + NX*NSIZE + 2*(NX+1)*N + NANG 
           CALL ERRT(46,'BPFQ; PROJ, ...',MEMWANT)
           RETURN
        ENDIF

C       CREATE LIST OF IMAGES FOR EACH RECONSTRUCTION
        RANDLIST(1:NANG/2)      = .TRUE.
        RANDLIST(NANG/2+1:NANG) = .FALSE.

        DO  K=1,NANG
           CALL RANDOM_NUMBER(HARVEST=X)
           IORD           = MIN(NANG,MAX(1,INT(X*NANG+0.5)))
           ITMP           = RANDLIST(IORD)
           RANDLIST(IORD) = RANDLIST(K)
           RANDLIST(K)    = ITMP
        ENDDO


#ifdef USE_MPI

C       --- BEGIN MPI VERSION ----------------------------------- MPI
C       DISTRIBUTE PARTICLES TO PROCESSORS.
C       NANGLOC IS THE NUMBER OF PARTICLES ASSIGNED TO EACH PROCESSOR.

C       CREATE FORWARD FFTW3 PLAN FOR 2D FFT ON BI USING ONE THREAD
        !write(6,*) ' BP32F;  inv plan on bi,n : ',n

        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 1,+1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
   
        ALLOCATE(PSIZE(NPROCS) ,NBASE(NPROCS), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MEMWANT = 2 * NPROCS 
           CALL ERRT(46,'BP 32F; PSIZE & NBASE',MEMWANT)
           RETURN
        ENDIF

C       SETPART RETURNS THE SIZE OF THE LOCAL PIECE AND GLOBAL OFFSET.
C       PARTITIONS NANG DATA ITEMS EVENLY AMONG AVAILABLE PROCESSORS.
C       PSIZE HAS NUMBER OF PARTICLES, NBASE HAS STARTING LOCATION.
        CALL SETPART(NANG, PSIZE, NBASE)

C       NANGLOC IS THE # OF PARTICLES TO BE USED FOR CURRENT PROCESS.
        NANGLOC = PSIZE(MYPID+1)

C       2-D IMAGES ARE DISTRIBUTED AND HELD IN PRJLOC ON EACH PROCESSOR
        ALLOCATE(PRJBUF(NX,NX,PSIZE(1)),
     &           PRJLOC(NX,NX,NANGLOC),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MEMWANT = NX*NX*PSIZE(1) + NX*NX*NANGLOC
           CALL ERRT(46,'BP 32F; PRJBUF, PRJLOC',MEMWANT)
           RETURN
        ENDIF

C       PROCESSOR 0 READS IMAGE FILES AND DISTRIBUTES THEM 
C       SO EACH PROCESSOR GETS DIFFERENT SET OF IMAGE FILES
C       (ASSUMES SUFFICIENT MEMORY TO HOLD NANG/NPROCS IMAGES)
C       FOR IPROC =1,  MYPID=0 READS AND KEEPS IMAGES FOR ITSELF
C       FOR IPROC =1,  MYPID>0 SKIPS READ 
C       FOR IPROC >1,  MYPID=0 READS AND BROADCASTS IMAGES
C       FOR IPROC >1,  MYPID>0 SKIPS READ AND RECEIVES  IMAGES 
           
        !if(mypid <= 0)write(6,'(a,i4)') 'In BP32Fq; nprocs:', nprocs
        !if(mypid <= 0)write(6,'(a,i4)') 'In BP32Fq; psize(1):', psize(1)

        DO IPROC = 1, NPROCS
           NANGLOC = PSIZE(IPROC)

C          READ IMAGES INTO THE PRJBUF BUFFER FIRST, THEN DISTRIBUTE

           DO JLOC = 1, NANGLOC
              JGLB = NBASE(IPROC) + JLOC
            
              NLET = 0     ! TO USE: lnblnkn IN FILGET
              CALL FILGET(FILPAT,FILNAM,NLET,INUMBRT(JGLB),IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999
              MAXIM = 0

              CALL OPFILEC(0,   .FALSE.,  FILNAM, INPROJ, 'O',
     &                     IFORM,   NX,    NY,   NSL,    MAXIM,
     &                     'DUMMY', .FALSE., IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999

              IF (MYPID <= 0) THEN
C                ONLY PROCESS WITH MYPID = 0 NEEDS TO READ
                 DO IROW = 1,NY
                    CALL REDLIN1P(INPROJ,PRJBUF(1,IROW,JLOC),NX,IROW)
                 ENDDO
              ENDIF
              CLOSE(INPROJ)
           ENDDO

           IF (IPROC > 1) THEN
#ifdef MPI_DEBUG
               write(6,'(a,i4)') 'In BP32FQ: iproc: ', iproc
#endif

               IF  (MYPID .EQ. 0) THEN
C                 SEND PROJECTION IMAGES TO OTHER PROCESSOR(S)
#ifdef MPI_DEBUG
                  write(6,*) 'In BP32FQ; Sending to mypid= ', mypid
#endif
                  CALL MPI_SEND(PRJBUF , NX*NX*NANGLOC, MPI_REAL,
     &                          IPROC-1, IPROC-1, ICOMM, MPIERR)
                  IF (MPIERR .NE. 0) THEN
                     WRITE(6,*) ' *** BP32FQ: SEND ERROR!'
                     STOP
                  ENDIF

               ELSEIF (MYPID .EQ. IPROC-1) THEN
C                 RECEIVE PROJECTION IMAGES FROM PROCESSOR 0
#ifdef MPI_DEBUG
                  write(6,'(a,i4)') 'In BP32FQ; Receiving mypid: ', mypid
#endif
                  CALL MPI_RECV(PRJLOC, NX*NX*NANGLOC, MPI_REAL,
     &                          0, MPI_ANY_TAG, ICOMM, ISTAT , MPIERR)
                  IF (MPIERR .NE. 0) THEN
                     WRITE(6,*) ' *** BP32FQ: RECV FAILED'
                     STOP
                  ENDIF
#ifdef MPI_DEBUG
                  write(6,*)' In BP32Fq: mypid: ',mypid,' received' 
#endif
              ENDIF 

           ELSEIF (MYPID .EQ. 0) THEN
C             KEEP PROJECTION IMAGES FOR THIS PROCESS 

#ifdef MPI_DEBUG
              write(6,'(a,i4)') 'In BP32FQ; Keeping mypid: ', mypid
#endif
              DO JLOC = 1, NANGLOC
                 DO ISAM = 1, NX
                     DO JROW = 1, NX
                        PRJLOC(ISAM,JROW,JLOC) = PRJBUF(ISAM,JROW,JLOC)
                     ENDDO
                  ENDDO
               ENDDO
           ENDIF
        ENDDO

        !write(6,'(a,i4)') 'In BP32FQ: deallocating now, pid:',mypid

        IF (ALLOCATED(PRJBUF)) DEALLOCATE(PRJBUF)

        IF (.NOT. ANGINDOC) THEN
C           GET ANGLES FROM HEADER
            ANGBUFT(1) = INUMBRT(K)
            CALL LUNGETVALS(INPROJ,IAPLOC+1,8,ANGBUFT(2),IRTFLG)
            CALL BUILDM(INUMBRT,DM,1,ANGBUFT,.FALSE.,SSDUM,
     &                  .FALSE.,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999
        ENDIF
        CLOSE(INPROJ)

C       PERFORM CALCULATIONS IN PARALLEL NOW
        !if(mypid <= 0)write(6,'(a,i4)') 'In BP32Fq; calculating' 

        ALLOCATE (WELOC(0:NX,N,N), XELOC(0:NX,N,N), 
     &            WOLOC(0:NX,N,N), XOLOC(0:NX,N,N), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MEMWANT = 4*(NX+1)*N*N 
           CALL ERRT(46,'BP 3F; WELOC, XELOC, ...',MEMWANT)
           RETURN
        ENDIF

       !write(6,*) ' BP32Fq: xeloc zero = ', mypid
        DO K=1,N
           DO J=1,N
              DO I=0,NX
                 XELOC(I,J,K) = CMPLX(0.0,0.0)
                 XOLOC(I,J,K) = CMPLX(0.0,0.0)
                 WELOC(I,J,K) = 0.0
                 WOLOC(I,J,K) = 0.0
              ENDDO
           ENDDO
        ENDDO

        INV     = +1      ! FORWARD FOURIER TRANSFORMS
        NANGLOC = PSIZE(MYPID+1)

        DO JLOC = 1, NANGLOC
           JGLB = NBASE(MYPID+1) + JLOC

C          PAD: PRJLOC  TO SIZE: N
           CALL PADD2(PRJLOC(1,1,JLOC),NX,BI,LSD,N)

C          FORWARD FOURIER TRANSFORM OF: BI
           !write(6,*) ' BP32F;  fft on bi,n : ',n
           CALL FMRS_2(BI,N,N,INV)

           DO J=1,N
              DO I=0,NX
                 BI(I,J) = BI(I,J) * (-1)**(I+J+1)
              ENDDO
           ENDDO

           DO  ISYM=1,MAXSYM
              IF (MAXSYM .GT. 1 )  THEN
C                SYMMETRIES, MULTIPLY MATRICES
                 DMS = MATMUL(SM(:,:,ISYM),DM(:,:,JGLB))
              ELSE
                 DMS = DM(:,:,JGLB)
              ENDIF

             IF (RANDLIST(JGLB)) THEN
                 DO J=-NX+1,NX
                    CALL ONELINE(J,N,NX,XELOC,WELOC,BI,DMS)
                 ENDDO
              ELSE
                 DO J=-NX+1,NX
                    CALL ONELINE(J,N,NX,XOLOC,WOLOC,BI,DMS)
                 ENDDO
              ENDIF
           ENDDO
        ENDDO

        !write(6,*) ' before reduce ; n,NX,nangloc: ', n,NX,nangloc
        !write(6,*) ' xe:' , xeloc(0,1,1),xeloc(0,1,2),xeloc(75,75,75)
        !write(6,*) ' xo:' , xoloc(0,1,1),xoloc(0,1,2),xoloc(75,75,75)
        !write(6,*) ' we:' , weloc(0,1,1),weloc(0,1,2),weloc(75,75,75)
        !write(6,*) ' wo:' , woloc(0,1,1),woloc(0,1,2),woloc(75,75,75)

C       SUM UP X AND W FROM THE LOCAL PIECES (XLOC, WLOC)
C       RESIDING ON EACH PROCESSOR
        NVAL = (NX + 1) * N
        DO K3 = 1, N
           CALL ALLREDUCE_MPI('BP32FQ','XELOC',
     &              XELOC(0,1,K3),XE(0,1,K3),NVAL, 'X','S',ICOMM)
           CALL ALLREDUCE_MPI('BP32FQ','XOLOC',
     &              XOLOC(0,1,K3),XO(0,1,K3),NVAL, 'X','S',ICOMM)
           CALL ALLREDUCE_MPI('BP32FQ','WELOC',
     &              WELOC(0,1,K3),WE(0,1,K3),NVAL, 'R','S',ICOMM)
           CALL ALLREDUCE_MPI('BP32FQ','WOLOC',
     &              WOLOC(0,1,K3),WO(0,1,K3),NVAL, 'R','S',ICOMM)
        ENDDO 

        IF (ALLOCATED(PRJLOC))   DEALLOCATE(PRJLOC)
        IF (ALLOCATED(XOLOC))    DEALLOCATE (XOLOC)
        IF (ALLOCATED(XELOC))    DEALLOCATE (XELOC)
        IF (ALLOCATED(WOLOC))    DEALLOCATE (WOLOC)
        IF (ALLOCATED(WELOC))    DEALLOCATE (WELOC)

#ifdef MPI_DEBUG
        write(6,*) 'BP32FQ; Completed global sum, mypid: ', mypid
#endif



C       --------------------  END OF MPI VERSION ----------------------    
#else

C       CREATE FFTW3 PLAN FOR 2D FFT ON BI USING ALL THREADS
        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,+1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
 
#ifdef  SP_MP
        LN1 = LN + 1         ! WHY?? LN ALWAYS EQUALED: 5 (IN CALLER)
        LN1 = 6              
#endif

        NWANT1 = 1
        NWANT2 = 1
        DO 
c          write(6,*)  'Projection #: ',IMGNUM

           IF (CALLRTSQ) THEN
C             READ IMAGE INTO: PROJTEMP

C             Reg. numbers for angle & shift =(6,7,8)
              CALL READV(INPROJ,PROJTEMP, NXLD,NX, NX,NX,1)

              IF (USE_FBS_INTERP) THEN
C                ROTATE & SHIFT FROM: PROJTEMP  INTO: PROJ
	         CALL RTSF(PROJTEMP,PROJ, 
     &                     NXLD,NX,NX,
     &                     ANGBUF(7,IMGNUM), 1.0,
     &                     ANGBUF(8,IMGNUM),ANGBUF(9,IMGNUM),IRTFLG)

              ELSE
C                ROTATE & SHIFT FROM: PROJTEMP  INTO: PROJ
	         CALL RTSQ(PROJTEMP,PROJ, 
     &                     NX,NX, NX,NX,
     &                     ANGBUF(7,IMGNUM), 1.0,
     &                     ANGBUF(8,IMGNUM),ANGBUF(9,IMGNUM),IRTFLG)
              ENDIF

              IF (LUNROT > 0) THEN
                 CALL WRTVOL(LUNROT,NX,NX,1,1, PROJ,IRTFLG)
                 !CLOSE(LUNROT)
              ENDIF

              IF (VERBOSE) WRITE(NOUT,90)IMGNUM,
     &                                   ANGBUF(7,IMGNUM),
     &                                   ANGBUF(8,IMGNUM),
     &                                   ANGBUF(9,IMGNUM)
90            FORMAT('  IMAGE: ',I6,
     &               '  ANGLE: ',G10.3,
     &               '  SHIFT: (',G10.3,',',G10.3,')')
           ELSE
C             READ IMAGE INTO: PROJ
              CALL REDVOL(INPROJ,NX,NX, 1,1, PROJ,IRTFLG)
           ENDIF

           IF (.NOT. ANGINDOC) THEN
C             GET PROJECTION ANGLES FROM HEADER
              ANGBUFT(1) = INUMBRT(NWANT1)
              CALL LUNGETVALS(INPROJ,IAPLOC+1,8,ANGBUFT(2),IRTFLG)
              CALL BUILDM(INUMBRT,DM,1,ANGBUFT,.FALSE.,SSDUM,
     &                    .FALSE.,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999
           ENDIF

C          PAD: PROJ  TO SIZE: N
           CALL PADD2(PROJ,NX,BI,LSD,N)

C          FORWARD FOURIER TRANSFORM OF: BI
           INV = +1
           CALL FMRS_2(BI,N,N,INV)

c$omp      parallel do private(i,j)
           DO J=1,N
              DO I=0,NX
                 BI(I,J) = BI(I,J) * (-1)**(I+J+1)
              ENDDO
           ENDDO

           DO  ISYM=1,MAXSYM
              IF (MAXSYM. GT. 1)  THEN
C                SYMMETRIES, MULTIPLY MATRICES
                 DMS = MATMUL(SM(:,:,ISYM),DM(:,:,IMGNUM))
              ELSE
                 DMS = DM(:,:,IMGNUM)
              ENDIF
#ifdef SP_MP
              IF (RANDLIST(NWANT1)) THEN
                 DO JT=1,LN1
c$omp               parallel do private(j)
                    DO J=-NX+JT,NX,LN1
                       CALL ONELINE(J,N,NX,XE,WE,BI,DMS)
                    ENDDO
                ENDDO
             ELSE
                DO JT=1,LN1
c$omp              parallel do private(j)
                   DO J=-NX+JT,NX,LN1
                      CALL ONELINE(J,N,NX,XO,WO,BI,DMS)
                   ENDDO
                ENDDO
             ENDIF
#else
             IF (RANDLIST(NWANT1)) THEN
                 DO J=-NX+1,NX
                    CALL ONELINE(J,N,NX,XE,WE,BI,DMS)
                 ENDDO
             ELSE
                 DO J=-NX+1,NX
                    CALL ONELINE(J,N,NX,XO,WO,BI,DMS)
                 ENDDO
             ENDIF
#endif
          ENDDO                           ! END OF SYMMETRIES LOOP

          IF (NWANT1 .GE. NANG) EXIT      ! END OF INPUT LIST
              
C         OPEN NEXT SET OF I/O FILES 
          CALL NEXTFILES(NWANT1, NWANT2,        INUMBRT,INUMBROUT, 
     &                   .FALSE.,LUNXM1,LUNXM2,
     &                   NANG,NANGT,            MAXIM1,MAXIM2,   
     &                   INPROJ,INPROJ,LUNROT,  FILPAT,FILPATOUT,
     &                   IMGNUM,IMGNUMOUT,IRTFLG) 

          IF (IRTFLG .LT. 0) EXIT         ! END OF INPUT STACK
          IF (IRTFLG .NE. 0) GOTO 9999    ! ERROR
        ENDDO                             ! END OF PROJECTIONS LOOP 
#endif
C       --------------------  END OF NON-MPI VERSION ------------------    

        !write(6,*) '   ' 
        !write(6,*) ' xe:' , xe(0,1,1),xe(0,1,2),xe(75,75,75)
        !write(6,*) ' xo:' , xo(0,1,1),xo(0,1,2),xo(75,75,75)
        !write(6,*) ' we:' , we(0,1,1),we(0,1,2),we(75,75,75)
        !write(6,*) ' wo:' , wo(0,1,1),wo(0,1,2),wo(75,75,75)


C       SYMMETRIZE PLANE 0 OF BOTH VOLUMES
c$omp   parallel sections
c$omp   section
        CALL  SYMPLANE0(XE,WE,NX,N)
c$omp   section
        CALL  SYMPLANE0(XO,WO,NX,N)
c$omp   end parallel sections

9999    IF (ALLOCATED(PROJ))      DEALLOCATE(PROJ)
        IF (ALLOCATED(BI))        DEALLOCATE(BI)
        IF (ALLOCATED(RANDLIST))  DEALLOCATE(RANDLIST)
        IF (ALLOCATED(INUMBROUT)) DEALLOCATE(INUMBROUT)
        CLOSE(LUNROTT)
        CLOSE(LUNXM2)

        END

C       ------------------- WINDKB2A -------------------------------

C       DERIVED FROM:  WINDKB2 (now in var3d.f) 

        SUBROUTINE WINDKB2A(BI,R, L,LSD,N, ALPHA,AAAA,NNN)

        REAL      :: BI(LSD,N,N)      ! INPUT
        REAL      :: R(L,L,L)         ! OUTPUT

        PARAMETER (QUADPI = 3.14159265358979323846)
        PARAMETER (TWOPI  = 2*QUADPI)

        IP = (N-L ) / 2 + MOD(L,2)

ccc$omp   parallel do private(i,j,k)
        DO K=1,L
           DO J=1,L
              DO I=1,L
                 R(I,J,K) = BI(IP+I,IP+J,IP+K)
              ENDDO
           ENDDO
        ENDDO

        L2   = (L / 2)**2
        L2P  = (L / 2 - 1)**2
        IP   = L / 2 + 1
        XNU  = REAL(NNN) / 2.0
        
        RI   = RIBSL(ALPHA,XNU)

        WKB0 = ALPHA**XNU / RI
        QRT  = (TWOPI * AAAA)**2
        TNR  = 0.0
        M    = 0

        DO K=1,L
           KMIPSQ = (K-IP)**2

           DO J=1,L
              JMIPSQ = (J-IP)**2 + KMIPSQ

              DO I=1,L

C                SPHERE RADIUS
                 LR = JMIPSQ + (I-IP)**2

                 IF (LR <= L2) THEN
C                  INSIDE SPHERE RADIUS

                   SIGMA = QRT * LR - ALPHA * ALPHA

                   IF (ABS(SIGMA) < 1.0E-7) THEN
                     WKB = 1.0

                   ELSEIF(SIGMA > 0.0) THEN
C                    2PI A R > ALPHA
                     ART = SQRT(SIGMA)
                     RI  = RJBSL(ART, XNU)
                     WKB = WKB0 * RI / ART**XNU

                   ELSE
C                    2PI A R < ALPHA
                     ART = SQRT(ABS(SIGMA))
                     RI  = RIBSL(ART,XNU)
                     WKB = WKB0 * RI / ART**XNU
                   ENDIF

C                  DIVIDE BY WEIGHT?
                   R(I,J,K) = R(I,J,K) / ABS(WKB)

                   IF (LR >= L2P .AND. LR <= L2) THEN
C                     RUNNING SUM OF VALUES
                      TNR = TNR + R(I,J,K)
                      M   = M + 1
                   ENDIF
                 ENDIF
              ENDDO
           ENDDO
        ENDDO

        TNR = TNR / REAL(M)

c$omp   parallel do private(i,j,jmipsq,k,kmipsq,lr)
        DO K=1,L
           KMIPSQ = (K-IP)**2

           DO J=1,L
              JMIPSQ = (J-IP)**2 + KMIPSQ

              DO I=1,L

                 LR = JMIPSQ + (I-IP)**2

                 IF (LR <= L2) THEN
C                   INSIDE SPHERE
                    R(I,J,K) = R(I,J,K) - TNR
                 ELSE
C                   OUTSIDE SPHERE
                    R(I,J,K) = 0.0
                 ENDIF
              ENDDO
           ENDDO
        ENDDO

        END



C       ----------------SYMPLANE0 ---------------------------------------
 
        SUBROUTINE  SYMPLANE0(X,W,N2,N)

        REAL    :: W(0:N2,N,N)
        COMPLEX :: X(0:N2,N,N)

C       SYMMETRIZE PLANE 0
        DO IZA=2,N2
           DO IYA=2,N2
              X(0,IYA,IZA)         = X(0,IYA,IZA) + 
     &                               CONJG(X(0,N-IYA+2,N-IZA+2))
              W(0,IYA,IZA)         = W(0,IYA,IZA) +
     &                               W(0,N-IYA+2,N-IZA+2)
              X(0,N-IYA+2,N-IZA+2) = CONJG(X(0,IYA,IZA))
              W(0,N-IYA+2,N-IZA+2) = W(0,IYA,IZA)
              X(0,N-IYA+2,IZA)     = X(0,N-IYA+2,IZA) +
     &                               CONJG(X(0,IYA,N-IZA+2))
              W(0,N-IYA+2,IZA)     = W(0,N-IYA+2,IZA) +
     &                               W(0,IYA,N-IZA+2)
              X(0,IYA,N-IZA+2)     = CONJG(X(0,N-IYA+2,IZA))
              W(0,IYA,N-IZA+2)     = W(0,N-IYA+2,IZA)
           ENDDO
        ENDDO

        DO IYA=2,N2
           X(0,IYA,1)     = X(0,IYA,1)+CONJG(X(0,N-IYA+2,1))
           W(0,IYA,1)     = W(0,IYA,1)+W(0,N-IYA+2,1)
           X(0,N-IYA+2,1) = CONJG(X(0,IYA,1))
           W(0,N-IYA+2,1) = W(0,IYA,1)
        ENDDO

        DO IZA=2,N2
           X(0,1,IZA)     = X(0,1,IZA)+CONJG(X(0,1,N-IZA+2))
           W(0,1,IZA)     = W(0,1,IZA)+W(0,1,N-IZA+2)
           X(0,1,N-IZA+2) = CONJG(X(0,1,IZA))
           W(0,1,N-IZA+2) = W(0,1,IZA)
        ENDDO

        END

C       ---------ADDADA ---------------------------------------

        SUBROUTINE  ADDADA(X,Y,N)

        REAL ::  X(N),Y(N)

C       JUST A SIMPLE ELEMENT BY ELEMENT SUM OF THE TWO ARRAYS 

c$omp   parallel do private(k)
        DO K=1,N
           X(K) = X(K) + Y(K)
        ENDDO

        END

