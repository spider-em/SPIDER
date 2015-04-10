
C++*********************************************************************
C
C  BP3F.F                                          MAY 02 P.A. Penczek
C              REPLACED BESSEL FUNCTIONS           MAY 02 P.A. Penczek
C              OPFILEC                             FEB 03 ArDean Leith
C              BUILDM PARAMETERS                   JUL 03 ArDean Leith
C              BUILDM PARAMETERS                   SEP 03 ArDean Leith
C              ALLOCATION ERROR HANDLING           OCT 04 ArDean Leith
C              OMP PRIVATE XX                      NOV 04 ArDean Leith
C              MPI DECONSTRUCTED FROM BP32D        NOV 06 ArDean Leith 
C              REWRITE                             DEC 06 ArDean Leith
C              FILLBESSIL                          DEC 08 ArDean Leith
C              CLOSE LUNIN                         SEP 09 ArDean Leith
C              OPFILES PARAMETERS                  DEC 10 ArDean Leith
C              PREVIOUSLY NAMED WIW2D**            JAN 11 ArDean Leith
C              ROT2QS --> RTSQ RENAMED             DEC 11 ArDean Leith
C              RTSF SUPPORT                        JAN 12 ArDean Leith
C              NSAM --> NX, RTSQ PARAM.            JAN 12 ArDean Leith
C              ILIST ALLOC                         JAN 13 ArDean Leith
C              HALT ADDED                          FEB 13 ArDean Leith
C
C=**********************************************************************
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright (C)2002,2014 P. A. Penczek & ArDean Leith                *
C=* University of Texas - Houston Medical School                       *
C=* Email:  pawel.a.penczek@uth.tmc.edu                                *
C=*                                                                    *
C=* SPIDER is free software; you can redistribute it and/or            *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* SPIDER is distributed in the hope that it will be useful,          *
C=* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C=* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
C=* General Public License for more details.                           *
C=*                                                                    *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program. If not, see <http://www.gnu.org/licenses> *
C=*                                                                    *
C=**********************************************************************
C
C   BP3F(CALLRTSQ,HALT)
C
C   PURPOSE:  FOR: 'BP 3F' IMPROVED WITH LESS MEMORY REQUIRED AND
C             MORE STACK USAGE. WITH CALLRTSQ CAN OUTPUT THE
C             ROTATED/SHIFTED FILES USED FOR THE BACKPROJECTION
C             FOR: 'RB 3F' ALSO
C
C  OPERATION: 'BP 3F'
C  BP3F  ---->  BUILDM1  
C    |          BUILDS  
C    |          BP3FQ  -->  FMRS_PLAN
C    |                      PADD2
C    |                      FMRS_2
C    |                      ONELINE
C    |                      SYMPLANE0
C    |                      NRMW2 
C    |                      WINDKB2A 
C    v  MPI
C    | 
C  BP3F_MPI --> BUILDM1  
C                BUILDS  
C                BP3FQ_MPI  --> BUILDM1
C                    |          FMRS_2
C                    |          ONELINE
C                    |          SYMPLANE0
C                    |          NRMW2 
C                    |          WINDKB2 
C                   COP 
C
CC     
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE BP3F(CALLRTSQ,HALT)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'F90ALLOC.INC'

        LOGICAL               :: CALLRTSQ,HALT

C       DOC FILE POINTERS
        REAL, POINTER         :: ANGBUF(:,:), ANGSYM(:,:)

        REAL,    ALLOCATABLE  :: DM(:,:),SM(:,:)
        REAL,    ALLOCATABLE  :: W(:,:,:)
        COMPLEX, ALLOCATABLE  :: X(:,:,:)
        REAL,    ALLOCATABLE  :: TEMP(:)
        REAL,    ALLOCATABLE  :: ROTANG(:),SHX(:),SHY(:)
        INTEGER, ALLOCATABLE  :: ILIST(:)

        LOGICAL               :: ANGINDOC,WANTVOL
        CHARACTER(LEN=MAXNAM) :: FILPAT,VOLNAM,VOLWGH
        CHARACTER(LEN=MAXNAM) :: ANGDOC

        INTEGER, PARAMETER    :: LUNIN   = 16
        INTEGER, PARAMETER    :: LUNVOL  = 17
        INTEGER, PARAMETER    :: LUNVOLV = 18
        INTEGER, PARAMETER    :: LUNVOLW = 19
        INTEGER, PARAMETER    :: LUNROTT = 20
        INTEGER, PARAMETER    :: LUNDOC  = 80
        INTEGER, PARAMETER    :: LUNXM1  = 0  ! UNUSABLE NEED #s
        INTEGER, PARAMETER    :: LUNXM2  = 82

#ifdef USE_MPI
C       KLUDGE TO OVERCOME BUG REPORTED BY L. ALAMO, NOV 2006
        CALL BP3F_MPI()
        RETURN
#else
        MYPID = -1
#endif

        NILMAX = NIMAXPLUS
        ALLOCATE(ILIST(NILMAX),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'BP32F; ILIST....',NILMAX)
           RETURN
        ENDIF 

C       OPEN INPUT IMAGE FILES 
        CALL OPFILES(0,LUNIN,LUNDOC,LUNXM1, 
     &             .TRUE.,FILPAT,NLET, 'O',
     &             ITYPE,NX,NY,NZ,MAXIM1,
     &             'IMAGE FILE NAME OR TEMPLATE (E.G. STK@****)~~',
     &             .FALSE., ILIST,NILMAX, 
     &             NDUM,NANG,IMGNUM, IRTFLG) 
        IF (IRTFLG .NE. 0) GOTO 999

        IF (NX .NE. NY) THEN
           CALL ERRT(101,'ONLY WORKS ON SQUARE IMAGES',NDUM)
           GOTO 999
        ENDIF

        MAXNUM = MAXVAL(ILIST(1:NANG))
         
C       NANG - TOTAL NUMBER OF IMAGES
        IF (MYPID .LE. 0) WRITE(NOUT,2001) NANG
2001    FORMAT('  NUMBER OF IMAGES: ',I7,I7)

C       RETRIEVE ARRAY WITH ANGLES DATA IN IT
C       PSI, THE, PHI, REF#, EXP#, INPLANE, SX, SY  
        MAXXT = 8 + 1
        MAXYT = MAXNUM
        CALL GETDOCDAT('ANGLES DOC',.TRUE.,ANGDOC,LUNDOC,.FALSE.,MAXXT,
     &                       MAXYT,ANGBUF,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        !write(6,*)'gotdocdat:MAXXT,MAXYT,MAXNUM',MAXXT,MAXYT,MAXNUM

C
C       RETRIEVE ARRAY WITH SYMMETRIES DATA IN IT
        MAXXS  = 4
        MAXSYM = 0
        CALL GETDOCDAT('SYMMETRIES DOC',.TRUE.,ANGDOC,LUNDOC,
     &                 .TRUE.,MAXXS, MAXSYM,ANGSYM,IRTFLG)
        IF (IRTFLG .NE. 0) MAXSYM = 1

        N2      = 2 * NX
        LSD     = N2 + 2 - MOD(N2,2)
        NMAT    = LSD * N2 * N2
        NDIMSYM = MAX(1,MAXSYM)

        ALLOCATE(DM(9,MAXNUM), 
     &           SM(9,NDIMSYM), 
     &           W(0:NX,N2,N2), 
     &           X(0:NX,N2,N2), STAT=IRTFLG)

        IF (IRTFLG .NE. 0) THEN 
           MWANT = 9*MAXNUM + 9*NDIMSYM + 3*(NX+1)*N2*N2 
           CALL ERRT(46,'BP3F, DM..', MWANT)
           GOTO 999
        ENDIF

C       GET ANGLES FROM DOCUMENT FILE AND PLACE IN DM
        CALL BUILDM1(ILIST,NANG, DM,ANGBUF,9,MAXNUM, IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999
        IF (.NOT. CALLRTSQ) DEALLOCATE(ANGBUF)
 
        IF (MAXSYM > 1)  THEN
C          HAVE SYMMETRIES, CONSTRUCT SM ANGLES ARRAY
           CALL BUILDS(SM,MAXSYM,ANGSYM(1,1),IRTFLG)
           DEALLOCATE(ANGSYM)
        ENDIF

        IF (HALT) THEN
C          OPEN THE NEW FOURIER VOLUME AND REAL WEIGHTS FILES

           IFORM  = 3
           MAXIM  = 0
           NXP1   = NX + 1
           NXP1T2 = 2  * NXP1

           !write(6,*) 'nx,n2,nxp1,nxp1t2:',nx,n2,nxp1,nxp1t2

           CALL OPFILEC(0,.TRUE.,VOLNAM,LUNVOLV,'U',IFORM,NXP1T2,N2,N2,
     &                  MAXIM,'PARTIAL VOLUME',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 999

           LENAT   = INDEX(VOLNAM,'@')
           IF (LENAT > 1) THEN
              VOLWGH = VOLNAM(1:LENAT-1) // '_W' // VOLNAM(LENAT:NLET)
           ELSE              
              VOLWGH = VOLNAM(1:NLET) // '_W' 
           ENDIF

C          OPEN WEIGHTS VOLUMES 
           IFORM  = 3
           MAXIM  = 0
           !write(6,*) 'volwgh:',volwgh(1:22)
           CALL OPFILEC(0,.FALSE.,VOLWGH,LUNVOLW,'U',IFORM,NXP1,N2,N2,
     &                  MAXIM,'PARTIAL WEIGHT',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 999
            
        ELSE

C          OPEN REGULAR OUTPUT VOLUME
           IFORM = 3
           CALL OPFILEC(0,.TRUE.,VOLNAM,LUNVOL,'U',IFORM,NX,NX,NX,
     &               MAXIM,'RECONSTRUCTED VOLUME',.FALSE.,IRTFLG)
           WANTVOL = (IRTFLG == 0) 
        ENDIF

C       BACK PROJECTION
        CALL BP3FQ(NX,X,W,LSD,N2, CALLRTSQ,FILPAT,ANGBUF,
     &        ILIST,NANG, DM,IMGNUM,SM,MAXSYM, 
     &        LUNIN,LUNROTT,LUNDOC,LUNXM1,LUNXM2,
     &        MAXIM1,HALT, MAXNUM,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        IF (HALT) THEN
C          SAVE THE FOURIER VOLUME AND REAL WEIGHTS TO FILES

           !X(0:NX,N2,N2)
           CALL WRTVOL(LUNVOLV,NXP1T2,N2, 1,N2, X,IRTFLG)
           CLOSE(LUNVOLV)

           !W(0:NX,N2,N2) 
           CALL WRTVOL(LUNVOLW,NXP1,N2, 1,N2, W,IRTFLG)
           CLOSE(LUNVOLW)

           IF (.NOT. WANTVOL) GOTO 999   ! FINISHED
        ENDIF


C       ADDITIONAL SYMMETRIZATION OF THE VOLUME IN REAL SPACE 05/03/02
	IF (MAXSYM > 1)  THEN
           ALLOCATE(TEMP(NX*NX*NX), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN 
              CALL ERRT(46,'BP 3F, TEMP',NX*NX*NX)
              GOTO 999
           ENDIF

	   CALL COP(X,TEMP, NX*NX*NX)  ! JUST COPIES X INTO: TEMP

c$omp      parallel do private(i,j,k)
           DO K=1,N2
              DO J=1,N2
                 DO I=0,NX
                    X(I,J,K) = CMPLX(0.0,0.0)
                 ENDDO
              ENDDO
           ENDDO

           IF (MOD(NX,2) == 0)  THEN
              KNX = NX/2-1
           ELSE
              KNX = NX/2
           ENDIF
           KLX = -NX/2

           CALL SYMVOL(TEMP,X,KLX,KNX,KLX,KNX,KLX,KNX,SM,MAXSYM)
           DEALLOCATE(TEMP)
        ENDIF


C       NOTE: NX=NY=NZ 
        CALL WRTVOL(LUNVOL,NX,NX, 1,NX, X,IRTFLG)

999     CLOSE(LUNVOL)
        CLOSE(LUNIN)

        IF (ALLOCATED(DM)) DEALLOCATE(DM)
        IF (ALLOCATED(SM)) DEALLOCATE(SM)
        IF (ALLOCATED(X))  DEALLOCATE(X)
        IF (ALLOCATED(W))  DEALLOCATE(W)

        END

 
C       ------------------ BP3FQ ---------------------------------

        SUBROUTINE BP3FQ(NX,X,W,LSD,N2, CALLRTSQ,FILPAT,ANGBUF,
     &                   INUMBRT,NANG, DM,IMGNUM,SM,MAXSYM,
     &                   LUNIN,LUNROTT,LUNDOC,LUNXM1,LUNXM2,
     &                   MAXIM1,HALT,MAXNUM,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        INTEGER                  :: NX,LSD,N2
        COMPLEX                  :: X(0:NX, N2,N2)
        REAL                     :: W(0:NX, N2,N2)
        LOGICAL                  :: CALLRTSQ
        REAL                     :: ANGBUF(9,MAXNUM)
        INTEGER                  :: NANG,IMGNUM,MAXSYM
        REAL                     :: DM(3,3,MAXNUM)
        REAL                     :: SM(3,3,MAXSYM)
        INTEGER                  :: LUNIN,LUNROTT,LUNDOC,LUNXM1,LUNXM2
        INTEGER                  :: MAXIM1
        LOGICAL                  :: HALT
        INTEGER                  :: IRTFLG

        REAL                     :: DMS(3,3)

        CHARACTER(LEN=*)         :: FILPAT

        INTEGER                  :: INUMBRT(NANG)
        REAL,    ALLOCATABLE     :: PROJ(:,:),PROJTEMP(:,:)
        COMPLEX, ALLOCATABLE     :: BI(:,:)
        CHARACTER (LEN=MAXNAM)   :: FILPATOUT,DOCNAM
        INTEGER, ALLOCATABLE     :: INUMBROUT(:)

        LOGICAL                  :: ANGINDOC
        INTEGER                  :: NXLD

        DOUBLE PRECISION         :: PI

C       COMMON: /TABS/ IS USED IN ONELINE, EXTRACTLINE, PUTLINE3, ETC
        INTEGER, PARAMETER          :: LTAB=4999
        COMMON  /TABS/ LN2,FLTB,TABI(0:LTAB)

        MYPID = -1                ! NOT MPI

C       GENERALIZED KAISER-BESSEL WINDOW ACCORDING TO LEWITT
        CALL FILLBESSIL(N2,LN2,LTAB,FLTB,TABI,ALPHA,AAAA,NNN)

        NSIZE  = 1
        LUNROT = 0    ! IF ZERO --> NO RTSQ OUTPUT WANTED
        NXLD   = NX

        IF (CALLRTSQ) THEN
           IF ( USE_FBS_INTERP) NXLD = NX + 2 - MOD(NX,2)

C          READ IMAGE AND ROTATE, SCALE & SHIFT INTO: PROJ
           NSIZE = NX
           ALLOCATE(INUMBROUT(NANG), STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN 
              CALL ERRT(46,'INUMBROUT', NANG)
              GOTO 9999
           ENDIF

C          OPEN OUTPUT IMAGE(S) FOR RTSQ/RTSF 
           ITYPE = 1
           CALL OPFILES(LUNIN,LUNROTT,LUNDOC,LUNXM2, 
     &            .TRUE., FILPATOUT,NLETO,'U',
     &            ITYPE,NX,NX,1, MAXIM2,
     &            'TRANSFORMED OUTPUT IMAGES TEMPLATE (E.G. ROT@****)~',
     &            .FALSE., INUMBROUT,NANG, 
     &            NDUM,NANGT, IMGNUMOUT, IRTFLG) 

C          IF IRTFLG IS NEGATIVE THEN NO RTSQ OUTPUT WANTED
           IF (IRTFLG == 0) LUNROT = LUNROTT
           IF (IRTFLG > 0) GOTO 9999
        ENDIF

        ALLOCATE(PROJ(NX,NX), 
     &           PROJTEMP(NXLD,NSIZE),
     &           BI(0:NX,N2), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           NEEDED = NX*NX + NXLD*NSIZE + 2*(NX+1)*N2*N2 +
     &             (NX+1)*N2
           CALL ERRT(46,'BP3F, PROJ...',NEEDED)
           RETURN
        ENDIF

C       CREATE FFTW3 PLAN FOR 2D FFT ON BI USING ALL THREADS
        CALL FMRS_PLAN(.TRUE.,BI,N2,N2,1, 0,+1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
 
c$omp   parallel do private(i,j,k)
        DO K=1,N2
           DO J=1,N2
              DO I=0,NX
                 X(I,J,K) = CMPLX(0.0,0.0)
                 W(I,J,K) = 0.0
              ENDDO
           ENDDO
        ENDDO

#ifdef  SP_MP
C	LN1 = LN + 1           ! WHY?? LN ALWAYS EQUALS: 5
	LN1 = 6                ! WHY?? LN ALWAYS EQUALS: 5
#endif

        NWANT1 = 1
        NWANT2 = 1
        DO 
c          write(6,*)  'Projection #: ',IMGNUM

           IF (CALLRTSQ) THEN
C             READ IMAGE INTO: PROJTEMP (X PADDED IF FBS)
              CALL READV(LUNIN,PROJTEMP,NXLD,NX, NX,NX,1)

C             Reg. numbers for angle & shift =(6,7,8)
  
              IF (USE_FBS_INTERP) THEN
C                ROTATE & SHIFT FROM: PROJTEMP  INTO: PROJ
	         CALL RTSF(PROJTEMP,PROJ, 
     &                     NXLD,NX,NX,
     &                     ANGBUF(7,IMGNUM), 1.0,
     &                     ANGBUF(8,IMGNUM),ANGBUF(9,IMGNUM),IRTFLG)

              ELSE
C                ROTATE & SHIFT FROM:PROJTEMP  INTO: PROJ
	         CALL RTSQ(PROJTEMP,PROJ, 
     &                     NX,NX, NX,NX,
     &                     ANGBUF(7,IMGNUM), 1.0,
     &                     ANGBUF(8,IMGNUM),ANGBUF(9,IMGNUM),IRTFLG)
              ENDIF

              IF (LUNROT > 0) THEN
                 CALL WRTVOL(LUNROT, NX,NX,1,1, PROJ,IRTFLG)
                 CLOSE(LUNROT)
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
              CALL REDVOL(LUNIN,NX,NX, 1,1, PROJ,IRTFLG)
           ENDIF


C          PAD: PROJ TO SIZE: N2
           CALL PADD2(PROJ,NX,BI,LSD,N2)

C          FOURIER TRANSFORM OF: BI
           INV = +1
           CALL FMRS_2(BI,N2,N2,INV)  

c$omp      parallel do private(i,j)
           DO J=1,N2
              DO I=0,NX
                 BI(I,J) = BI(I,J)*(-1)**(I+J+1)
              ENDDO
           ENDDO

           DO ISYM=1,MAXSYM
             IF (MAXSYM > 1)  THEN
C               SYMMETRIES, MULTIPLY MATRICES
                DMS = MATMUL(SM(:,:,ISYM),DM(:,:,IMGNUM))
             ELSE
                DMS = DM(:,:,IMGNUM)
             ENDIF

#ifdef SP_MP
	     DO JT=1,LN1
c$omp           parallel do private(j)
                DO J=-NX+JT,NX,LN1
                   CALL ONELINE(J,N2,NX,X,W,BI,DMS)
                ENDDO
             ENDDO
#else
             DO J=-NX+1,NX
               CALL ONELINE(J,N2,NX,X,W,BI,DMS)
             ENDDO
#endif
           ENDDO                          ! END OF SYMMETRIES LOOP

           IF (NWANT1 >= NANG) EXIT      ! END OF LIST

C          OPEN NEXT SET OF I/O FILES 
           CALL NEXTFILES(NWANT1, NWANT2, INUMBRT,INUMBROUT, 
     &                     .FALSE.,LUNXM1,LUNXM2,
     &                     NANG,NANGT,  
     &                     MAXIM1,MAXIM2,   
     &                     LUNIN,LUNIN,LUNROT, FILPAT,FILPATOUT,
     &                     IMGNUM,IMGNUMOUT,IRTFLG) 
           IF (IRTFLG .LT. 0) EXIT         ! END OF INPUT STACK
           IF (IRTFLG .NE. 0) GOTO 9999    ! ERROR
        ENDDO                              ! END OF PROJECTIONS LOOP 

C       SYMMETRIZE PLANE 0
        CALL SYMPLANE0(X,W,NX,N2)

        IF ( .NOT. HALT) THEN
C          WEIGHT AND FOURIER TRANSFORM
           CALL NRMW2(X,W,NX,N2)

C          WINDOW
           CALL WINDKB2A(X,X,NX,LSD,N2,ALPHA,AAAA,NNN)
        ENDIF

        IRTFLG = 0

9999    IF (ALLOCATED(PROJ)) DEALLOCATE(PROJ)
        IF (ALLOCATED(BI))   DEALLOCATE(BI)

        CLOSE(LUNIN)
        CLOSE(LUNROTT)
        CLOSE(LUNXM2)

        END













C       ------------------ BP3F_MERGE ------------------------------

        SUBROUTINE BP3F_MERGE()

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'F90ALLOC.INC'
C       FROM CMBLOCK.INC --> INTEGER,DIMENSION(NIMAX) :: INUMBR

C       DOC FILE POINTERS
        REAL, POINTER         :: ANGBUF(:,:), ANGSYM(:,:)

        REAL,    ALLOCATABLE  :: SM(:,:)
        COMPLEX, ALLOCATABLE  :: X(:,:,:),XOUT(:,:,:)
        REAL,    ALLOCATABLE  :: W(:,:,:),WOUT(:,:,:)

        CHARACTER(LEN=MAXNAM) :: ANGDOC
        CHARACTER(LEN=MAXNAM) :: FILPATV,FILPATW,VOLNAM,FILNAMW

        INTEGER, PARAMETER    :: LUNINV  = 16
        INTEGER, PARAMETER    :: LUNINW  = 17
        INTEGER, PARAMETER    :: LUNVOL  = 18
        INTEGER, PARAMETER    :: LUNWAT  = 19
        INTEGER, PARAMETER    :: LUNDOC  = 80
        INTEGER, PARAMETER    :: LUNXMV  = 0  ! UNUSABLE NEED #s
        INTEGER, PARAMETER    :: LUNXMW  = 0  ! UNUSABLE NEED #s

        INTEGER, PARAMETER    :: LTAB=4999
        REAL                  :: TABI(0:LTAB)
        INTEGER               :: LN2,NNN
        REAL                  :: FLTB,ALPHA,AAAA

C       OPEN FIRST INPUT VOLUME FILE 
        NILMAX = NIMAX
        CALL OPFILES(0,LUNINV,LUNDOC,LUNXMV, 
     &             .TRUE.,FILPATV,NLET, 'O',
     &             ITYPEV,NXV,NYV,NZV,MAXIMV,
     &             'PARTIAL VOLUME FILES TEMPLATE (E.G. BPSTK@****)~~',
     &             .FALSE., INUMBR,NILMAX, 
     &             NDUM,NVOLV,IMGNUM, IRTFLG) 
        IF (IRTFLG .NE. 0) GOTO 999

        IF (MYPID <= 0) WRITE(NOUT,2001) NVOLV
2001    FORMAT('  NUMBER OF VOLUMES: ',I7,I7)


        LENAT   = INDEX(FILPATV,'@')
        LENSTAR = INDEX(FILPATV,'*')
        IF (LENAT > 1) THEN
           FILPATW = FILPATV(1:LENAT-1) // '_W' // FILPATV(LENAT:NLET)
        ELSEIF (LENSTAR > 1) THEN             
           FILPATW = FILPATV(1:LENAST-1) // '_W'//FILPATV(LENAST:NLET)
        ELSE
           CALL ERRT(101,'INVALID TEMPLATE',NDUM)
           GOTO 999
        ENDIF

C       OPEN FIRST INPUT WEIGHT FILE 
        FILNAMW = FILPATW
        CALL OPFILES(0,LUNINW,LUNDOC,LUNXMW, 
     &             .FALSE.,FILPATW,NLETW, 'O',
     &             ITYPEW,NXW,NYW,NZW,MAXIMW,
     &             FILNAMW,
     &             .TRUE., INUMBR,NILMAX, 
     &             IDUM,NVOLW,IMGNUMW, IRTFLG) 
        IF (IRTFLG .NE. 0) GOTO 999

        CALL SIZCHK(UNUSED,0,NYV,NZV,0,
     &                     0,NYW,NZW,0,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        NXP1T2 = NXV
        N2     = NYV  
        NY     = N2 / 2 
        NX     = NY
        NVOLW  = NVOLV
        MAXNUM = MAXVAL(INUMBR(1:NVOLV))
        NMAT   = NXV * NYV * NZV 
        LSD    = N2 + 2 - MOD(N2,2)

        !write(6,*)' maxsym:',maxsym,lsd,nmat,ndimsym

C       RETRIEVE ARRAY WITH SYMMETRIES DATA IN IT
        MAXXS  = 4
        MAXSYM = 0
        CALL GETDOCDAT('SYMMETRIES DOC',.TRUE.,ANGDOC,LUNDOC,
     &                 .TRUE.,MAXXS, MAXSYM,ANGSYM,IRTFLG)

        IF (IRTFLG .NE. 0) MAXSYM = 1
        NDIMSYM = MAX(1,MAXSYM)

        ALLOCATE(SM(9,NDIMSYM), 
     &           W(0:NX,N2,N2), WOUT(0:NX,N2,N2),
     &           X(0:NX,N2,N2), XOUT(0:NX,N2,N2),
     &          STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           MWANT = 9*MAXNUM + 6*(NX+1)*N2*N2
           CALL ERRT(46,'BP3F_MERGE, SM..', MWANT)
           GOTO 999
        ENDIF
 
        IF (MAXSYM > 1)  THEN
C          HAVE SYMMETRIES, CONSTRUCT SM ANGLES ARRAY FROM ANGDOC
           CALL BUILDS(SM,MAXSYM,ANGSYM(1,1),IRTFLG)
           DEALLOCATE(ANGSYM)
        ENDIF

C       OPEN OUTPUT VOLUME
        IFORM = 3
        CALL OPFILEC(LUNINV,.TRUE.,VOLNAM,LUNVOL,'U',IFORM,NX,NX,NX,
     &                  MAXIM,'RECONSTRUCTED VOLUME',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999
        

        NWANTV = 1
        NWANTW = 1
        DO 
           !write(6,*)  ' Volume #: ',IMGNUM,nwantv,nvol

C          LOAD VOLUME INTO: X
           CALL REDVOL(LUNINV,NXV,NYV, 1,NZV, X,IRTFLG)

C          LOAD WEIGHTS INTO: W
           CALL REDVOL(LUNINW,NXW,NYW, 1,NZW, W,IRTFLG)

C          ADD VOLUMES & WEIGHTS, PLACE IN XOUT AND WOUT
           CALL ADDADA(XOUT,X, NMAT)
           CALL ADDADA(WOUT,W, NMAT/2)

           !write(6,*)  ' Added: ',NXV,nyv,nxw,nyw,nmat

           IF (NWANTV >= NVOLV) EXIT      ! END OF LIST

C          OPEN NEXT INPUT VOLUME FILE
           CALL NEXTFILE(NWANTV,  INUMBR, 
     &                  .TRUE.,   LUNXMV,
     &                   NVOLV,   MAXIMV,   
     &                   LUNINV,  0, 
     &                   FILPATV, 'O',
     &                   IMGNUMV, IRTFLG) 
           IF (IRTFLG < 0)    EXIT        ! END OF INPUT STACK
           IF (IRTFLG .NE. 0) GOTO 999    ! ERROR

C          OPEN NEXT INPUT WEIGHT FILE
           CALL NEXTFILE(NWANTW,   INUMBR, 
     &                   .TRUE.,   LUNXMW,
     &                   NVOLW,    MAXIMW,   
     &                   LUNINW,   0, 
     &                   FILPATW,  'O',
     &                   IMGNUMW,  IRTFLG) 
           IF (IRTFLG < 0) THEN   ! END OF INPUT STACK
              CALL ERRT(101,'INSUFFICIENT WEIGHT VOLUMES', NE)
              GOTO 999    ! ERROR
           ENDIF
           IF (IRTFLG .NE. 0) GOTO 999    ! ERROR
        ENDDO                             ! END OF PROJECTIONS LOOP 


C       WEIGHT:XOUT  USING: WOUT & THEN INVERSE FOURIER TRANSFORM IT
        CALL NRMW2(XOUT,WOUT,NX,N2)

        !write(6,*)  'called nrmw2: NX,N2',NX,N2

C       GENERALIZED KAISER-BESSEL WINDOW ACCORDING TO LEWITT
        CALL FILLBESSIL(N2,LN2,LTAB,FLTB,TABI,ALPHA,AAAA,NNN)

C       WINDOWS FROM 2X PADDED NON-CENTERED SPACE TO REAL VOLUME
        CALL WINDKB2A(XOUT,XOUT,NX,LSD,N2, ALPHA,AAAA,NNN)

C       ADDITIONAL SYMMETRIZATION OF THE VOLUME IN REAL SPACE 
	IF (MAXSYM > 1)  THEN

           CALL COP(XOUT, X,NX*NX*NX)  ! COPIES: XOUT INTO: X

C          ZEROES XOUT
c$omp      parallel do private(i,j,k)
           DO K=1,N2
              DO J=1,N2
                 DO I=0,NX
                    X(I,J,K) = CMPLX(0.0,0.0)
                 ENDDO
              ENDDO
           ENDDO

           IF (MOD(NX,2) == 0)  THEN
              KNX = NX/2-1
           ELSE
              KNX = NX/2
           ENDIF
           KLX = -NX/2

C          INPUTS: X   OUTPUTS TO: XOUT
           CALL SYMVOL(X,XOUT,KLX,KNX,KLX,KNX,KLX,KNX,SM,MAXSYM)

        ENDIF

C       NOTE: NX=NY=NZ 
        CALL WRTVOL(LUNVOL,NX,NX, 1,NX, XOUT,IRTFLG)

999     CLOSE(LUNINV)
        CLOSE(LUNINW)
        CLOSE(LUNVOL)

        IF (ALLOCATED(SM))    DEALLOCATE(SM)
        IF (ALLOCATED(X))     DEALLOCATE(X)
        IF (ALLOCATED(W))     DEALLOCATE(W)
        IF (ALLOCATED(XOUT))  DEALLOCATE(XOUT)
        IF (ALLOCATED(WOUT))  DEALLOCATE(WOUT)

        END



	
C++*********************************************************************
C
C BP3F_MPI.F   DECONSTRUCTED FROM BP32D        NOV  06 ARDEAN LEITH
C                 ORIGINAL BP3F GAVE WRONG RESULTS
C                 WHEN USED UNDER MPI
C                 FFTW3 CORRECTIONS               SEP  08 ARDEAN LEITH
C **********************************************************************
C
C  BP3F_MPI
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE BP3F_MPI

C       NOTE: STUPID TRANSFORM OF N2-->N AND N2/2-->N2 !!!!al
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'F90ALLOC.INC'

C       DOC FILE POINTERS
        REAL, POINTER         :: ANGBUF(:,:), ANGSYM(:,:)

        REAL,    ALLOCATABLE  :: DM(:,:), SM(:,:) 
        COMPLEX, ALLOCATABLE  :: XE(:,:,:)
        REAL,    ALLOCATABLE  :: WE(:,:,:)
        INTEGER, ALLOCATABLE  :: ILIST(:)

        REAL,    ALLOCATABLE  :: TEMP(:)
        LOGICAL               :: ANGINDOC
        CHARACTER(LEN=1)      :: NULL = CHAR(0)
        CHARACTER(LEN=MAXNAM) :: ANGDOC
        CHARACTER(LEN=MAXNAM) :: FILPAT,FILNAM
        CHARACTER(LEN=MAXNAM) :: VOLNAM
 
        INTEGER, PARAMETER    :: LUNVOL   = 18
        INTEGER, PARAMETER    :: LUNIN   = 19
        INTEGER, PARAMETER    :: LUNDOC  = 77


#ifndef USE_MPI
        PRINT *, ' THIS ROUTINE for MPI ONLY'
        STOP
#else
        CALL SET_MPI(ICOMM,MYPID,MPIERR)  ! SETS ICOMM AND MYPID

        NILMAX = NIMAXPLUS
        ALLOCATE(ILIST(NILMAX),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'BP32F; ILIST....',NILMAX)
           RETURN
        ENDIF 

        CALL FILELIST(.TRUE.,LUNIN,FILPAT,NLET,ILIST,NILMAX,NANG,
     &                 'TEMPLATE FOR 2-D IMAGES',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CLOSE(LUNIN)     ! USED FOR SELECTION FILE??
        MAXNUM = MAXVAL(ILIST(1:NANG))

C       NANG - TOTAL NUMBER OF IMAGES
        IF (MYPID <= 0) WRITE(NOUT,2001) NANG
2001    FORMAT('  NUMBER OF IMAGES: ',I8)

C       RETRIEVE ARRAY WITH ANGLES DATA IN IT
        ANGINDOC = .TRUE.

C       PSI, THE, PHI, REF#, EXP#, INPLANE, SX, SY  
        MAXXT    = 8 + 1
        MAXYT    = MAXNUM
        CALL GETDOCDAT('ANGLES DOC',.TRUE.,ANGDOC,LUNDOC,.FALSE.,MAXXT,
     &                 MAXYT,ANGBUF,IRTFLG)
        IF (IRTFLG .NE. 0) ANGINDOC = .FALSE.

C       RETRIEVE ARRAY WITH SYMMETRIES DATA IN IT
        MAXXS  = 0
        MAXSYM = 0
        CALL GETDOCDAT('SYMMETRIES DOC',.TRUE.,ANGDOC,LUNDOC,
     &                 .TRUE.,MAXXS,MAXSYM,ANGSYM,IRTFLG)
        IF (IRTFLG .NE. 0)  MAXSYM = 1

C       OPEN FIRST IMAGE FILE TO DETERMINE NX, NY, NSL
        CALL FILGET(FILPAT,FILNAM,0,ILIST(1),IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FILNAM,LUNIN,'O',IFORM,NX,NY,NSL,
     &             MAXIM,'DUMMY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        IF (MYPID  .LE. 0) CLOSE(LUNIN)

        N2   = 2 * NX
        LSD  = N2 + 2 - MOD(N2,2)
        NMAT = LSD * N2 * N2

        IF (ANGINDOC) THEN
C          GET ANGLES FROM DOC. FILE 
           ALLOCATE(DM(9,MAXNUM), STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN
              MEMWANT = 9 * NANG
              CALL ERRT(46,'BP 3F, DM',MEMWANT)
              GOTO 9999
           ENDIF

           CALL BUILDM1(ILIST,NANG, DM,ANGBUF,9,MAXNUM,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
        
           IF (ASSOCIATED(ANGBUF)) DEALLOCATE(ANGBUF)
        ELSE
C          GET ANGLES FROM IMAGE FILE HEADER 
           ALLOCATE(DM(9,1), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'BP 3F, DM',IER)
              GOTO 9999
           ENDIF
        ENDIF

        IF (MAXSYM .GT. 1)  THEN
C          HAVE SYMMETRIES 
           ALLOCATE(SM(9,MAXSYM), STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN 
              MEMWANT = 9 * MAXSYM
              CALL ERRT(46,'BP 3F, SM',MEMWANT)
              GOTO 9999
           ENDIF
           CALL BUILDS(SM,MAXSYM,ANGSYM(1,1),IRTFLG)
           IF (ASSOCIATED(ANGSYM)) DEALLOCATE(ANGSYM)

       ELSE
           ALLOCATE(SM(1,1), STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN 
               CALL ERRT(46,'BP 3F, SM(1,1)',IER)
               GOTO 9999
           ENDIF
        ENDIF

        ALLOCATE(XE(0:NX,N2,N2),
     &           WE(0:NX,N2,N2), STAT=IRTFLG)

        IF (IRTFLG.NE.0) THEN 
C           X ARRAYS ARE COMPLEX SO 3 NOT 2
            MEMWANT = 3 * ((NX)+1) * N2 * N2
            CALL ERRT(46,'BP 3F; XE & WE ',MEMWANT)
            GOTO 9999
        ENDIF

        CALL BP3FQ_MPI(NX,XE,WE,LSD,N2, FILPAT,FILNAM,LUNIN,
     &                  ILIST,NANG, DM,NANG,SM,MAXSYM,ANGINDOC,IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9999

C       NOW XE IS READY, SYMMETRIZE IF NECESSARY
C       ADDITIONAL SYMMETRIZATION OF THE VOLUME XE IN REAL SPACE 05/03/02
	IF (MAXSYM .GT. 1)  THEN
           ALLOCATE (TEMP(NX*NX*NX), STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN 
              MEMWANT = NX * NX * NX 
              CALL ERRT(46,'BP 3F, TEMP',MEMWANT)
              GOTO 9999
           ENDIF

C          COPY XE to TEMP
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

C       OPEN OUTPUT VOLUME
        IFORM = 3
        CALL OPFILEC(0,.TRUE.,VOLNAM,LUNVOL,'U',IFORM,NX,NX,NX,
     &                  MAXIM,'RECONSTRUCTED 3-D',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       NOTE: NX=NY=NZ 
        CALL WRTVOL(LUNVOL,NX,NX, 1,NX, XE,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (MYPID .LE. 0) CLOSE(LUNIN)
        IF (MYPID .LE. 0) CLOSE(LUNVOL)

9999    IF (ALLOCATED(DM))      DEALLOCATE(DM)
        IF (ALLOCATED(SM))      DEALLOCATE(SM)
        IF (ALLOCATED(XE))      DEALLOCATE(XE)
        IF (ALLOCATED(WE))      DEALLOCATE(WE)
        IF (ALLOCATED(TEMP))    DEALLOCATE(TEMP)
        IF (ASSOCIATED(ANGBUF)) DEALLOCATE(ANGBUF)
        IF (ASSOCIATED(ANGSYM)) DEALLOCATE(ANGSYM)
#endif
 
        END

C       ---------------- BP3FQ_MPI -------------------------------------

        SUBROUTINE BP3FQ_MPI(NX,XE,WE,LSD,N2, 
     &                 FILPAT,FILNAM,LUNIN, 
     &                 INUMBRT,NGOT,DM, NANG,SM,MAXSYM,ANGINDOC,IRTFLG)

        INCLUDE 'CMLIMIT.INC'

        LOGICAL               :: ANGINDOC
        REAL,  ALLOCATABLE    :: PROJ(:,:)
        COMPLEX, ALLOCATABLE  :: BI(:,:)
        REAL                  :: ANGBUF(4)
        REAL                  :: WE(0:NX,N2,N2)
        COMPLEX               :: XE(0:NX,N2,N2)
        INTEGER               :: INUMBRT(NANG)
        REAL                  :: DM(3,3,NANG)
        REAL                  :: SM(3,3,MAXSYM)
        REAL                  :: DMS(3,3)

        CHARACTER(LEN=*)      :: FILPAT,FILNAM

        REAL   , ALLOCATABLE  :: PRJLOC(:,:,:)
        REAL   , ALLOCATABLE  :: PRJBUF(:,:,:)
        REAL   , ALLOCATABLE  :: WELOC(:,:,:)
        COMPLEX, ALLOCATABLE  :: XELOC(:,:,:)
        INTEGER, ALLOCATABLE  :: PSIZE(:), NBASE(:)

        LOGICAL, PARAMETER    :: SPIDER_SIGN  = .TRUE.
        LOGICAL, PARAMETER    :: SPIDER_SCALE = .TRUE.
 
C       COMMON: /TABS/ IS USED IN ONELINE, EXTRACTLINE, PUTLINE3, ETC
        INTEGER, PARAMETER    :: LTAB = 499
        COMMON  /TABS/ LN2,FLTB,TABI(0:LTAB)

C       THIS MPI VERSION IS MEMORY INTENSIVE.
C       2-D IMAGES ARE READ INTO MEMORY AND DISTRIBUTED.
C       EACH PROCESSOR WILL HOLD ROUGHLY NANG/NPROC  2-D IMAGES. 

#ifndef USE_MPI
        MYPID = -1
        PRINT *, ' THIS ROUTINE FOR MPI COMPILATION AND USE ONLY'
        STOP
#else
        INCLUDE 'mpif.h'
        INTEGER          :: ISTAT(MPI_STATUS_SIZE)
                          
        ICOMM = MPI_COMM_WORLD
        CALL MPI_COMM_RANK(ICOMM, MYPID, MPIERR)
        CALL MPI_COMM_SIZE(ICOMM, NPROCS, MPIERR)

C       GENERALIZED KAISER-BESSEL WINDOW ACCORDING TO LEWITT
        CALL FILLBESSIL(N2,LN2,LTAB,FLTB,TABI,ALPHA,AAAA,NNN)

        DO K=1,N2
           DO J=1,N2
              DO I=0,NX
                 XE(I,J,K) = CMPLX(0.0,0.0)
                 WE(I,J,K) = 0.0
              ENDDO
           ENDDO
        ENDDO

        ALLOCATE(PROJ(NX,NX),
     &           BI(0:NX,N2), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           MEMWANT = NX*NX + (NX + 1) * N2  
           CALL ERRT(46,'BP 3F, PROJ, ...',MEMWANT)
           RETURN
        ENDIF

C       DISTRIBUTE PARTICLES TO PROCESSORS.
C       NANGLOC IS THE NUMBER OF PARTICLES ASSIGNED TO EACH PROCESSOR.

        ALLOCATE(PSIZE(NPROCS),
     &           NBASE(NPROCS), 
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MEMWANT = 2 * NPROCS 
           CALL ERRT(46,'BP 3F, PSIZE, NBASE...',MEMWANT)
           RETURN
        ENDIF

C       SETPART RETURNS SIZE OF LOCAL PIECE AND THE GLOBAL OFFSET.

        CALL SETPART(NANG, PSIZE, NBASE)
        NANGLOC = PSIZE(MYPID+1)

C       2-D IMAGES ARE DISTRIBUTED AND HELD IN PRJLOC ON EACH PROCESSOR

        ALLOCATE(PRJBUF(NX,NX,PSIZE(1)),
     &           PRJLOC(NX,NX,NANGLOC),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MEMWANT = NX*NX*PSIZE(1) + NX*NX*NANGLOC
           CALL ERRT(46,'BP 3F, PRJBUF, PRJLOC',MEMWANT)
           RETURN
        ENDIF

C       PROCESSOR 0 READS IMAGE FILES AND DISTRIBUTE THEM.
C       (THIS VERSION ASSUMES THAT THERE IS SUFFICIENT
C        MEMORY TO HOLD NANG/NPROCS IMAGES)

        DO IPROC = 1, NPROCS
           NANGLOC = PSIZE(IPROC)

C          READ IMAGES INTO THE BUFFER FIRST, THEN DISTRIBUTE

           DO JLOC = 1, NANGLOC
              JGLB = NBASE(IPROC) + JLOC
              CALL FILGET(FILPAT,FILNAM,0,INUMBRT(JGLB),IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999
              MAXIM = 0
              CALL OPFILEC(0      ,.FALSE., FILNAM, LUNIN, 'O'  ,
     &                     IFORM  , NX  , NX  , NSL   , MAXIM,
     &                     'DUMMY',.FALSE., IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999
              CALL READV1P(LUNIN,PRJBUF(1,1,JLOC),
     &                     NX,NX,NX,NX,1)
            ENDDO

           IF (IPROC .GT. 1) THEN
               IF  (MYPID .EQ. 0) THEN
C                 SEND TO ANOTHER PROCESSOR
                  CALL MPI_SEND(PRJBUF , NX*NX*NANGLOC, MPI_REAL,
     &                          IPROC-1, IPROC-1      , ICOMM    ,
     &                          MPIERR)
                  IF (MPIERR .NE. 0) THEN
                     WRITE(6,*) ' BP32DQ: SEND ERROR!'
                     STOP
                  ENDIF
               ELSEIF (MYPID .EQ. IPROC-1) THEN
C                 RECEIVE PROJECTION IMAGES FROM PROCESSOR 0
                  CALL MPI_RECV(PRJLOC, NX*NX*NANGLOC, MPI_REAL,
     &                          0     , MPI_ANY_TAG  , ICOMM    ,
     &                          ISTAT , MPIERR)
                  IF (MPIERR .NE. 0) THEN
                     WRITE(6,*) ' BP32DQ: RECV FAILED'
                     STOP
                  ENDIF
              ENDIF  
           ELSEIF (MYPID .EQ. 0) THEN
C             KEEP FOR MYSELF
              DO JLOC = 1, NANGLOC
                 DO ISAM = 1, NX
                     DO JROW = 1, NX
                        PRJLOC(ISAM,JROW,JLOC) = PRJBUF(ISAM,JROW,JLOC)
                     ENDDO
                  ENDDO
               ENDDO
           ENDIF
        ENDDO
        IF (ALLOCATED(PRJBUF)) DEALLOCATE(PRJBUF)

        IF (.NOT. ANGINDOC) THEN
C           GET ANGLES FROM HEADER
            ANGBUF(1) = INUMBRT(K)
            CALL LUNGETVALS(LUNIN,IAPLOC+1,3,ANGBUF(2),IRTFLG)

            CALL BUILDM1(INUMBRT,1, DM,ANGBUF,4,1,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999
        ENDIF
        CLOSE(LUNIN)

C       PERFORM CALCULATIONS IN PARALLEL NOW

C       CREATE FFTW3 PLAN FOR THIS PROCESSOR
        INV       = +1            ! FORWARD FFT
        NUMTHWANT = 1             ! MAY RETURN THIS

C       CREATE FFTW3 PLAN FOR 2D FFT ON BI USING ONE THREAD
        CALL FMRS_PLAN(.TRUE.,BI, N2,N2,1, NUMTHWANT,INV,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        ALLOCATE (WELOC(0:NX,N2,N2), 
     &            XELOC(0:NX,N2,N2), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MEMWANT = 2 * (NX+1)*N2*N2 
           CALL ERRT(46,'BP 3F, WELOC, XELOC',MEMWANT)
           GOTO 9999
        ENDIF

        DO K=1,N2
           DO J=1,N2
              DO I=0,NX
                 XELOC(I,J,K) = CMPLX(0.0,0.0)
                 WELOC(I,J,K) = 0.0
              ENDDO
           ENDDO
        ENDDO

        NANGLOC = PSIZE(MYPID+1)
        DO JLOC = 1, NANGLOC
           JGLB = NBASE(MYPID+1) + JLOC

C          PAD: PRJLOC TO SIZE: N2
           CALL PADD2(PRJLOC(1,1,JLOC),NX,BI,LSD,N2)

C          FOURIER TRANSFOR OF: BI
           CALL FMRS_2(BI,N2,N2,INV)

           DO J=1,N2
              DO I=0,NX
                 BI(I,J) = BI(I,J)*(-1)**(I+J+1)
              ENDDO
           ENDDO

           DO ISYM=1,MAXSYM
              IF (MAXSYM .GT. 1)  THEN
C                SYMMETRIES, MULTIPLY MATRICES
                 DMS = MATMUL(SM(:,:,ISYM),DM(:,:,JGLB))
              ELSE
                 DMS = DM(:,:,JGLB)
              ENDIF

              DO J=-NX+1,NX
                 CALL ONELINE(J,N2,NX,XELOC,WELOC,BI,DMS)
              ENDDO
           ENDDO
        ENDDO

C       SUM UP X & W FROM LOCAL PIECES (XLOC, WLOC) RESIDING ON EACH PROCESSOR

        DO K3 = 1, N2
           CALL MPI_ALLREDUCE(XELOC(0,1,K3), XE(0,1,K3), (NX+1)*N2, 
     &                        MPI_COMPLEX  , MPI_SUM   , ICOMM    ,
     &                        MPIERR)
           IF (MPIERR. NE. 0) THEN
              WRITE(6,*) ' BP32DQ: FAILED TO ALLREDUCE XELOC'
              STOP
           ENDIF

           CALL MPI_ALLREDUCE(WELOC(0,1,K3), WE(0,1,K3), (NX+1)*N2, 
     &                        MPI_REAL     , MPI_SUM   , ICOMM    ,
     &                        MPIERR)
           IF (MPIERR. NE. 0) THEN
              WRITE(6,*) ' BP32DQ: FAILED TO ALLREDUCE WELOC'
              STOP
           ENDIF
        ENDDO 

C       SYMMETRIZE VOLUME
        CALL SYMPLANE0(XE,WE,NX,N2)

C       WEIGHT AND INVERSE FOURIER TRANSFORM
        CALL NRMW2(XE,WE,NX,N2)

C       WINDOW
        CALL WINDKB2A(XE,XE,NX,LSD,N2,ALPHA,AAAA,NNN)



9999    IF (ALLOCATED(PRJLOC)) DEALLOCATE(PRJLOC)
        IF (ALLOCATED(XELOC))  DEALLOCATE(XELOC)
        IF (ALLOCATED(WELOC))  DEALLOCATE(WELOC)
        IF (ALLOCATED(PROJ))   DEALLOCATE(PROJ)
        IF (ALLOCATED(BI))     DEALLOCATE(BI)
#endif

        END





C++*********************************************************************
C
C  BUILDM1.F        MERGED WITH REANG                JUL 03 ARDEAN LEITH
C                   BYLIST ADDED                     SEP 03 ARDEAN LEITH
C                   BYLIST FILLSS,SS REMOVED         FEB 13 ARDEAN LEITH
C **********************************************************************
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright (C) 1985-2013  Health Research Inc.                      *
C=*                                                                    *
C=* HEALTH RESEARCH INCORPORATED (HRI),                                *   
C=* ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455.                   *
C=*                                                                    *
C=* Email:  spider@wadsworth.org                                       *
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
C **********************************************************************
C
C  BUILDM1(ILIST,DM,NANG, ANGBUF,MAXX,MAXY, IRTFLG)
C
C  PURPOSE: BULID ROTATION MATRICES FROM THREE EULERIAN ANGLES
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE BUILDM1(ILIST,NANG, DM,ANGBUF,MAXX,MAXY, IRTFLG)

        INCLUDE 'CMBLOCK.INC' 

        INTEGER :: ILIST(NANG)
        INTEGER :: NANG
        REAL    :: DM(9,MAXY), ANGBUF(MAXX,MAXY)
        INTEGER :: MAXX,MAXY
        INTEGER :: IRTFLG

        INTEGER :: ICOMM,MYPID,MPIERR,K,ITMP,ICOUNT,J
        REAL    :: SSDUM

        LOGICAL, PARAMETER :: FILLSS = .FALSE.


        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C       READ ANGLES FROM THE DOCUMENT FILE.
C       ORDER IN THE DOCUMENT FILE IS PSI, THETA, PHI AND ANGLES 
C       ARE IN DEGREES! IN ANG ARRAY IT IS OTHER WAY AROUND
C       OUTPUT IS COMPACTED TO 1...NANG LINES (NOT BY SELECTOR)

        DO K=1,NANG

C          GET ANGLE SELECTOR INDEX FROM ILIST
           ITMP   = ILIST(K)

           !write(6,*) 'nang,k,itmp,maxy:',nang,k,itmp,maxy

           ICOUNT = ANGBUF(1,ITMP)
           IF (ICOUNT <= 0) THEN
C             MISSING KEY
              CALL ERRT(102,'MISSING ANGLE FOR IMAGE',ITMP)
              IRTFLG = 1
              RETURN
           ENDIF

           CALL CANG(ANGBUF(4,ITMP),ANGBUF(3,ITMP),ANGBUF(2,ITMP),
     &               FILLSS,SSDUM,DM(1,ITMP))

           IF (MYPID <= 0 .AND. VERBOSE) THEN
              WRITE(NOUT,333) K,(ANGBUF(J,ITMP),J=2,4)
333           FORMAT('  PROJECTION #',I7,
     &               '  PSI=',F6.1,' THETA=',F6.1,' PHI=',F6.1)
           ENDIF
         ENDDO

        IRTFLG = 0

        END


C       ------------------- PADD2 -------------------------------

        SUBROUTINE PADD2(PROJ,L,BI,LSD,N)

C       PADS: PROJ OF SIZE: L  INTO: BI  WITH SIZE: N

        REAL             :: PROJ(L,L),BI(LSD,N)

        DOUBLE PRECISION :: QS
        DOUBLE PRECISION :: DKLP

        KLP = 0
        R   = L/2
        QS  = 0.0D0

C       ESTIMATE AVERAGE OUTSIDE THE CIRCLE
        CALL ASTA_D(PROJ,L,R,QS,DKLP)

        QS = QS / REAL(DKLP)

C       ZEROS ALL OF: BI
c$omp   parallel do private(i,j)
        DO J=1,N
           DO I=1,N
              BI(I,J) = 0.0
           ENDDO
        ENDDO

C       FOR L ODD ADD ONE.  N IS ALWAYS EVEN
        IP = (N-L) / 2 + MOD(L,2)

c$omp   parallel do private(i,j)
        DO J=1,L
           DO I=1,L
              BI(IP+I,IP+J) = PROJ(I,J) - QS
           ENDDO
        ENDDO

        END

C       ------------------- NRMW2 -------------------------------

        SUBROUTINE NRMW2(R,W,N2,N)

        COMPLEX :: R(0:N2,N,N)
        REAL    :: W(0:N2,N,N)
        INTEGER :: N2,N

        INTEGER :: I,J,K,INV

c$omp   parallel do private(i,j,k)
        DO K=1,N
           DO J=1,N
              DO I=0,N2
                IF (W(I,J,K) > 0.1)  THEN
		   R(I,J,K) = R(I,J,K) * (-1)**(I+J+K)/W(I,J,K)
		ELSE
		   R(I,J,K) = (0.0,0.0)
		ENDIF
              ENDDO
           ENDDO
        ENDDO

        INV = -1
        CALL FMRS_3(R,N,N,N,INV) 

        END

	





	
