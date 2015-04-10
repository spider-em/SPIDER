
C++*********************************************************************
C                                                                      *
C  BPCG.F    CORRECTIONS APPLIED ON VOLUME SIDE  11/09/98              *
C            COMPRESSION OF ANGLES               08/14/96              *
C            USED PROJT FOR BCKPJ CALL           FEB 2000 ARDEAN LEITH *
C            OPFILEC                             FEB 2003 ARDEAN LEITH *
C            VERBOSE                             FEB 2006 ARDEAN LEITH *
C            MPI BUG FIXED                       OCT 2008 ARDEAN LEITH *
C            REFACTORED                          OCT 2008 ARDEAN LEITH *
C            REFACTORED                          DEC 2010 ARDEAN LEITH *
C            COMMON PAR REMOVED,                 DEC 2010 ARDEAN LEITH * 
C            THREE OUTPUTS                       JAN 2011 ARDEAN LEITH *
C            REPCG --> BPCG                      JAN 2011 ARDEAN LEITH *
C            MPI_RECV BUG                        APR 2011 ARDEAN LEITH *
C            FBS ADDED                           OCT 2011 G KISHCHENKO *
C            FBS2                                DEC 2011 ARDEAN LEITH *
C            RENAMED SUBROUTINES                 DEC 2011 ARDEAN LEITH *
C            FBS2 FOR ALL THREE                  JAN 2012 ARDEAN LEITH *
C            LIN RI TRAP BUG                     APR 2012 ARDEAN LEITH *
C            KLP_8                               APR 2012 ARDEAN LEITH *
C            D_KLP,D_KLS                           JAN 2013 ARDEAN LEITH *
C            gfort COMPILER BUG?                 MAY 2013 ARDEAN LEITH *                                                            *
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
C  BPCG()
C
C  PURPOSE:  CALCULATES 3D RECONSTRUCTION USING CONJUGATE GRADIENTS 
C            WITH REGULARIZATION.
C            RECONSTRUCTION IS CALCULATED INSIDE A SPHERE OF
C            SPECIFIED RADIUS
C            RECONSTRUCTION KEPT IN SQUARE TO INTRODUCE OTHER 
C            CONSTRAINTS. AVERAGE OUTSIDE WINDOW IS SUBTRACTED.
C            SYMMETRIES NOT IMPLEMENTED
C            CAN CALCULATE 2 SAMPLED VOLUMES ALSO
C
C  MEMORY DEMAND:  NON-MPI KEEPS 4 VOLUMES
C
C  CALLS:    BPCG  --->     PREPCUB-S
C                  --->     BPCG2      ---> ATASQ
C                  --->                -->  RPRQ
C                  --->                -->  BCKPJ_LIN or BCKPJ_FBS
C                  --->     HIANG
C                  --->     BPCG_3_LIN --> BCKPJ_LIN
C                  --->     BPCG_3_FBS --> BCKPJ_FBS
C          
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE BPCG

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'F90ALLOC.INC'

        REAL, POINTER         :: PANG(:,:)
        REAL, ALLOCATABLE     :: CB(:,:,:)
        REAL, ALLOCATABLE     :: ANG(:,:), DM(:,:)
        REAL, ALLOCATABLE     :: BCKN(:)

        INTEGER, ALLOCATABLE  :: ICUBE(:,:)
        INTEGER, ALLOCATABLE  :: LB(:) 

        INTEGER, ALLOCATABLE  :: ILIST(:),ILIST1(:),ILIST2(:)
 
        LOGICAL               :: WANT3,USELISTS,FBS_WANTED
        CHARACTER(LEN=MAXNAM) :: ANGDOC,FILPAT    ! MAXNAM FROM CMLIMIT
        CHARACTER(LEN=MAXNAM) :: FILVOL 
        CHARACTER(LEN=1)      :: ANSW
        CHARACTER(LEN=1)      :: NULL = CHAR(0)
        LOGICAL               :: ASKNAM = .TRUE.
        LOGICAL               :: FOUROK = .FALSE.

        INTEGER, PARAMETER    :: LUNPROJ = 20 
        INTEGER, PARAMETER    :: LUNDOC  = 80 
        INTEGER, PARAMETER    :: LUNXM   = 0    ! SELFILE NOT ALLOWED
        INTEGER, PARAMETER    :: LUNANG  = 81 
        INTEGER, PARAMETER    :: LUNVOL  = 21
        INTEGER, PARAMETER    :: LUNVOL1 = 22
        INTEGER, PARAMETER    :: LUNVOL2 = 23

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID
     
C       OPEN FIRST INPUT FILE
C       RETURNS: NANG = NUMBER OF ANGLES = NUMBER OF PROJECTIONS
        NILMAX  = NIMAXPLUS      ! FROM CMLIMIT
        ALLOCATE(ILIST1(NILMAX),
     &           ILIST2(NILMAX),
     &           ILIST(NILMAX),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'BPCG; ILIST....',3*NILMAX)
           RETURN
        ENDIF 

C       OPEN FIRST OVERALL INPUT FILE
        CALL OPFILES(0,LUNPROJ,LUNDOC,LUNXM,  
     &             ASKNAM,FILPAT,NLET, 'O',
     &             IFORM ,NSAM,NROW,NSLICE,NSTACK,
     &             'TEMPLATE FOR IMAGE FILES~',
     &             FOUROK, ILIST,NILMAX, 
     &             NDUM,NANG,IMG1, IRTFLG) 
        IF (IRTFLG .NE. 0) RETURN
        
        IF (FILPAT(NLET:NLET) .EQ. '@') THEN
           CALL ERRT(101,'OPERATION DOES NOT WORK ON WHOLE STACKS',NE)
           GOTO 9999
        ENDIF  

        MAXNUM = MAXVAL(ILIST(1:NANG))

C       NANG - TOTAL NUMBER OF IMAGES
        IF (MYPID <= 0) WRITE(NOUT,90) NANG
90      FORMAT('  NUMBER OF IMAGES:',I8)
        
        WANT3    = (FCHAR(4:7) .EQ. 'CG 3')   ! WANT THREE VOLUMES
        USELISTS = (FCHAR(4:8) .EQ. 'CG 3L')  ! WANT THREE LISTS

        !write(6,*)want3,uselists,fchar(4:8)
        IF (WANT3 .AND. USELISTS) THEN 
           CALL FILELIST(.FALSE.,LUNDOC,NULL,NLETP,
     &                   ILIST1,NILMAX,NANG1,
     &                   'IMAGES FOR FIRST SAMPLE VOLUME',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           CALL FILELIST(.FALSE.,LUNDOC,NULL,NLETP,
     &                   ILIST2,NILMAX,NANG2,
     &                   'IMAGES FOR SECOND SAMPLE VOLUME',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

        ELSEIF (WANT3) THEN
C          RANDOMLY DIVIDE ILIST INTO 2 LISTS: ILIST1,ILIST2 
           CALL MAKETWOLISTS(ILIST,NANG,ILIST1,NANG1,ILIST2,NANG2)

           WRITE(NOUT,*) ' Random lists:',NANG1, NANG2
        ENDIF

        RI = (NSAM / 2) - 2    ! DEFAULT VALUE
        CALL RDPRM1S(RI,NDUM,'RADIUS OF RECONSTRUCTED OBJECT',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       RETRIEVE ARRAY WITH ANGLES DATA IN IT
        MAXXT = 4
        MAXYT = MAXNUM
        CALL GETDOCDAT('ANGLES DOC',.TRUE.,ANGDOC,
     &                  LUNANG,.FALSE.,MAXXT,
     &                  MAXYT,PANG,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(101,'CAN NOT OPEN ANGLES FILE',NE)
           GOTO 9999
        ENDIF  

        CALL RDPRMC (ANSW, NLETI, .TRUE.,
     &     'LINEAR OR FBS INTERPOLATION (L,F)', NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        FBS_WANTED = (ANSW(1:1) == 'F')  

        N     = NSAM   ! LINEAR DIMENSION OF PROJECTIONS & RESTORATION

C       IDUM IS A DUMMY VARIABLE, VALUE OF NN IS DETERMINED
        LDP   = N / 2 + 1
        LDPNM = LDP

        CALL PREPCUB_S(N,NN,IDUM,RI,.FALSE.,LDP)     ! RETURNS: NN

C       USE NN TO ALLOCATE: ICUBE 
C       TOTAL MEMORY IS VOLUMES: CB, BCKN, BCKE...
C       PLUS TWO 2D PROJECTIONS.
C       CB   - BACK-PROJECTED ORIGINAL PROJECTIONS, READ FROM FILE
C       BCKE - WORKING VOLUME
C       BCKN - CURRENT RECONSTRUCTION
        ALLOCATE (ICUBE(5,NN), 
     &            CB(N,N,N), 
     &            ANG(3,NANG),
     &            DM(9,NANG),
     &            LB(NANG), 
     &            BCKN(N*N*N), 
     &            STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           MWANT = 5*NN + N*N*N + 3*NANG + 9*NANG + NANG + N*N*N 
           CALL ERRT(46,'BPCG, ICUBE...',MWANT)
           GOTO 9999
        ENDIF

C       MAKES LIST OF VOXEL LOCS ON EACH LINE VOLUME WITHIN RADIUS          *
        CALL PREPCUB_S(N,NN,ICUBE,RI,.TRUE.,LDP)   ! RETURNS: IPCUBE

C       OPEN OUTPUT VOLUME
        MAXIM = 0
        IFORM = 3
        CALL OPFILEC(0,.TRUE.,FILVOL,LUNVOL,'U',IFORM,N,N,N,
     &             MAXIM,'RECONSTRUCTED VOLUME',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (WANT3) THEN
           MAXIM = 0
           IFORM = 3
           CALL OPFILEC(0,.TRUE.,FILVOL,LUNVOL1,'U',IFORM,N,N,N,
     &             MAXIM,'FIRST SAMPLE VOLUME',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FILVOL,LUNVOL2,'U',IFORM,N,N,N,
     &             MAXIM,'SECOND SAMPLE VOLUME',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0)  GOTO 9999
        ENDIF

        ERRM = 0.00005
        CHIM = 0.0
        CALL RDPRM2S(ERRM,CHIM,NDUM,'ERROR LIMIT & CHI^2 LIMIT',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        MAXIT = 20
        MODE  = 1
        CALL RDPRIS(MAXIT,MODE,NOT_USED,'ITERATION LIMIT, MODE',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        ANOISE = 2000
        CALL RDPRM1S(ANOISE,NOT_USED,'LAMBDA',IRTFLG)  ! NOISE SUPRESSOR
        IF (IRTFLG .NE. 0) GOTO 9999

C       FIND NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)

        NXLD = N + 2 - MOD(N,2)

        WRITE(NOUT,92)'  CREATING VOLUME ------------------------------'
92      FORMAT(/,A)

C       BACK-PROJECTION AND BACKGROUND CORRECTION, RETURNS: CHI2
        CALL BPCG_2(N,NANG,CB,ANG,ILIST,ICUBE,NN,DM,RI,PANG,
     &               FILPAT,MAXIM,BNORM,CHI2, 
     &               LUNPROJ,LDP,LDPNM, FBS_WANTED,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       COMPRESS ANGLES, CONVERTS ANGLES TO 'DM' FORMAT
        CALL HIANG(ANG,NANG,DM,LB,LO)     ! ALTERS ANG,DM,LB,& LO !!
        NANG = LO
        IF (MYPID <= 0) WRITE(NOUT,91) NANG
91      FORMAT('  EFFECTIVE NUMBER OF ANGLES:    ',I6)

        DI      = N / 2
        IRADIUS = INT(RI)

        IF (FBS_WANTED) THEN
           IF (MYPID <= 0) WRITE(NOUT,*) ' USING FBS INTERPOLATION'

           IF (MODE == 1 .AND. (RI+1.0) > DI  .OR. 
     &         MODE == 2 .AND. (RI+2.0) > DI  .OR.
     &         MODE == 3 .AND. (RI+3.0) > DI)  THEN
              CALL ERRT(102,'RADIUS TOO LARGE',IRADIUS)
              GOTO 9999
           ENDIF   

           CALL BPCG_3_FBS(CB,BCKN,NXLD,N,ICUBE,NN,DM,LB,NANG,RI,
     &                 NUMTH,BNORM,CHI2,
     &                 ERRM,CHIM,MAXIT,MODE,ANOISE, LDP,LDPNM,IRTLFG)

        ELSE
           IF (MYPID <= 0) WRITE(NOUT,*) ' USING LINEAR INTERPOLATION'

           IF (MODE == 1 .AND. (RI+1.0+1.0) > DI  .OR. 
     &         MODE == 2 .AND. (RI+2.0+1.0) > DI .OR.
     &         MODE == 3 .AND. (RI+3.0+1.0) > DI)  THEN
              CALL ERRT(102,'RADIUS TOO LARGE',IRADIUS)
              GOTO 9999
           ENDIF   

           CALL BPCG_3_LIN(CB,BCKN,NXLD,N,ICUBE,NN,DM,LB,NANG,RI,
     &                 NUMTH,BNORM,CHI2,
     &                 ERRM,CHIM,MAXIT,MODE,ANOISE, LDP,LDPNM,IRTLFG)
        ENDIF
        IF (IRTFLG .NE. 0) GOTO 9999

C       SAVE OVERALL OUTPUT VOLUME
        CALL WRTVOL(LUNVOL,N,N,1,N,BCKN,IRTFLG)

        IF (WANT3) THEN

           WRITE(NOUT,*) 
     &          ' CREATING FIRST SAMPLED VOLUME ------------------'
           WRITE(NOUT,*) ' '
           !write(6,*) 'Creating  first output sampled volume -----'

C          OPEN NEXT INPUT PROJECTION
           K = 1
           CALL NEXTFILE(K,       ILIST1, 
     &                   FOUROK,  0,
     &                   NANG1,   MAXIM,   
     &                   LUNPROJ, 0, 
     &                   FILPAT,  'O',
     &                   IMGNUM,  IRTFLG)

C          BACK-PROJECTION AND BACKGROUND CORRECTION, RETURNS: CHI2
           CALL BPCG_2(N,NANG1,CB,ANG,ILIST1,ICUBE,NN,DM,RI,PANG,
     &               FILPAT,MAXIM,BNORM,CHI2, 
     &               LUNPROJ,LDP,LDPNM, FBS_WANTED,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          COMPRESS ANGLES, CONVERTS ANGLES TO 'DM' FORMAT
           CALL HIANG(ANG,NANG1,DM,LB,LO)     ! ALTERS NANG1!!
           NANG1 = LO
           IF (MYPID <= 0) WRITE(NOUT,91) NANG1

           IF (FBS_WANTED) THEN
              IF (MYPID <= 0) WRITE(NOUT,*) ' USING FBS INTERPOLATION'
              CALL BPCG_3_FBS(CB,BCKN,NXLD,N,ICUBE,NN,DM,LB,NANG1,RI,
     &                    NUMTH,BNORM,CHI2,
     &                    ERRM,CHIM,MAXIT,MODE,ANOISE, LDP,LDPNM,IRTLFG)

           ELSE
              IF (MYPID <= 0) WRITE(NOUT,*)' USING LINEAR INTERPOLATION'
              CALL BPCG_3_LIN(CB,BCKN,NXLD,N,ICUBE,NN,DM,LB,NANG1,RI,
     &                 NUMTH,BNORM,CHI2,
     &                 ERRM,CHIM,MAXIT,MODE,ANOISE, LDP,LDPNM,IRTLFG)
           ENDIF
           IF (IRTFLG .NE. 0) GOTO 9999

C          SAVE FIRST OUTPUT SAMPLED VOLUME
           CALL WRTVOL(LUNVOL1,N,N,1,N,BCKN,IRTFLG)



           WRITE(NOUT,*) 
     &         ' CREATING SECOND SAMPLED VOLUME ------------------'
           WRITE(NOUT,*) ' '
           !write(6,*) 'Creating second output sampled volume -----'

C          OPEN NEXT INPUT PROJECTION
           K = 1
           CALL NEXTFILE(K,       ILIST2, 
     &                   FOUROK,  0,
     &                   NANG2,   MAXIM,   
     &                   LUNPROJ, 0, 
     &                   FILPAT,  'O',
     &                   IMGNUM,  IRTFLG)

C          BACK-PROJECTION AND BACKGROUND CORRECTION. RETURNS: CHI2
           CALL BPCG_2(N,NANG2,CB,ANG,ILIST2,ICUBE,NN,DM,RI,PANG,
     &                  FILPAT,MAXIM,BNORM,CHI2, 
     &                  LUNPROJ,LDP,LDPNM, FBS_WANTED,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          COMPRESS ANGLES, CONVERTS ANGLES TO 'DM' FORMAT
           CALL HIANG(ANG,NANG2,DM,LB,LO)     
           NANG2 = LO
           IF (MYPID <= 0) WRITE(NOUT,91) NANG2 

           IF (FBS_WANTED) THEN
              IF (MYPID <= 0) WRITE(NOUT,*) ' USING FBS INTERPOLATION'
              CALL BPCG_3_FBS(CB,BCKN,NXLD,N,ICUBE,NN,DM,LB,NANG2,RI,
     &                    NUMTH,BNORM,CHI2,
     &                    ERRM,CHIM,MAXIT,MODE,ANOISE, LDP,LDPNM,IRTLFG)

           ELSE
              IF (MYPID <= 0)WRITE(NOUT,*)' USING LINEAR INTERPOLATION'
              CALL BPCG_3_LIN(CB,BCKN,NXLD,N,ICUBE,NN,DM,LB,NANG2,RI,
     &                 NUMTH,BNORM,CHI2,
     &                 ERRM,CHIM,MAXIT,MODE,ANOISE, LDP,LDPNM,IRTLFG)
           ENDIF
           IF (IRTFLG .NE. 0) GOTO 9999

C          SAVE SECOND OUTPUT SAMPLED VOLUME2
           CALL WRTVOL(LUNVOL2,N,N,1,N, BCKN,IRTFLG)

        ENDIF

9999    IF (ALLOCATED(BCKN))  DEALLOCATE(BCKN)
        IF (ALLOCATED(LB))    DEALLOCATE(LB)
        IF (ALLOCATED(DM))    DEALLOCATE(DM)
        IF (ALLOCATED(ANG))   DEALLOCATE(ANG)
        IF (ALLOCATED(CB))    DEALLOCATE(CB)
        IF (ALLOCATED(ICUBE)) DEALLOCATE(ICUBE)
        IF (ASSOCIATED(PANG)) DEALLOCATE(PANG)
        IF (ALLOCATED(ILIST)) DEALLOCATE(ILIST)
        IF (ALLOCATED(ILIST1))DEALLOCATE(ILIST1)
        IF (ALLOCATED(ILIST2))DEALLOCATE(ILIST2)
                
        CLOSE(LUNVOL) 
        CLOSE(LUNVOL1) 
        CLOSE(LUNVOL2) 
        CLOSE(LUNPROJ) 
  
        END


C ***************************** BPCG_2 *********************************
C
C  BPCG_2
C
C  PURPOSE:  LOADS PROJECTIONS AND ANGLES. FINDS BACKGROUND OUTSIDE
C            OF PROJECTED CIRCLE.  BACK-PROJECTS IMAGES INTO VOLUME.
C            DETERMINES AVERAGE OUTSIDE OF CIRCLE WHICH WILL LATER BE
C            SUBTRACTED FROM FINAL VOLUME.
C
C  *********************************************************************

        SUBROUTINE BPCG_2(N,NANG,CUBE,ANG,ILIST,IPCUBE,NN,DM,
     &                    RI,ANGBUF,FINPAT,MAXIM, BNORM,CHI2,
     &                    LUNPROJ,LDP,LDPNM,FBS_WANTED, IRTFLG)

        USE TYPE_KINDS
 
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        INTEGER               :: N,NANG
        REAL                  :: CUBE(N,N,N),ANG(3,NANG)
        INTEGER               :: ILIST(NANG),IPCUBE(5,NN)
        INTEGER               :: NN
        REAL                  :: DM(9,NANG)
        REAL                  :: RI
        REAL                  :: ANGBUF(4,NANG)
        CHARACTER(LEN=MAXNAM) :: FINPAT
        INTEGER               :: MAXIM
        REAL                  :: BNORM,CHI2
        INTEGER               :: LUNPROJ,LDP,LDPNM,IRTFLG
        LOGICAL               :: FBS_WANTED

        REAL                  :: PROJ(N*N)
        DOUBLE PRECISION      :: D_ABA,D_SUS,D_SSQ
        CHARACTER(LEN=MAXNAM) :: FILNAM

        LOGICAL               :: FOUROK = .FALSE.
        REAL                  :: ADUM
        INTEGER(KIND=I_8)     :: KLP_8
        DOUBLE PRECISION      :: D_KLP,D_KLP_LOC, D_KLS,D_KLS_LOC

#ifdef USE_MPI

        include 'mpif.h'
        REAL,    ALLOCATABLE  :: CUBE_LOC(:,:,:)
        REAL,    ALLOCATABLE  :: PRJBUF(:,:,:),PRJLOC(:,:,:)
        REAL,    ALLOCATABLE  :: ANG_LOC(:,:)
        INTEGER, ALLOCATABLE  :: PSIZE(:)
        INTEGER, ALLOCATABLE  :: NBASE(:)
        DOUBLE PRECISION      :: D_ABA_LOC,D_SUS_LOC,D_SSQ_LOC

        INTEGER               :: MPISTAT(MPI_STATUS_SIZE)

        ICOMM = MPI_COMM_WORLD
        CALL MPI_COMM_RANK(ICOMM, MYPID , MPIERR)
        CALL MPI_COMM_SIZE(ICOMM, NPROCS, MPIERR)

        ALLOCATE(PSIZE(NPROCS),
     &           NBASE(NPROCS), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,' BPCG_2; PSIZE & NBASE.',2*NPROCS)
           STOP
        ENDIF
#else
C       AUTOMATIC ARRAY
        REAL               :: PROJT(4,N*N)

        REAL, ALLOCATABLE  :: PROJPAD(:,:)  ! PADDED 2D PROJ ARRAY
        REAL, ALLOCATABLE  :: XDER   (:,:)  ! X  DERIVATIVE OF PROJ
        REAL, ALLOCATABLE  :: YDER   (:,:)  ! Y  DERIVATIVE OF PROJ
        REAL, ALLOCATABLE  :: XYDER  (:,:)  ! XY DERIVATIVE OF PROJ

        MYPID = -1
        IF (FBS_WANTED) THEN
           NXLD = N + 2 - MOD(N,2)
           ALLOCATE(PROJPAD(NXLD, N),
     &              XDER   (NXLD, N),
     &              YDER   (NXLD, N),
     &              XYDER  (NXLD, N),
     &              STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              MWANT = 4*NXLD*N
              CALL ERRT(46,'BPCG_2; PROJPAD,...',MWANT)
              GOTO 9999
           ENDIF
        ENDIF
#endif

        D_ABA  = 0.0D0
        D_SUS  = 0.0D0
        D_SSQ  = 0.0D0

        D_KLP  = 0.0D0
        D_KLS  = 0.0D0
        
c$omp   parallel do private(i,j,k)
        DO K=1,N
           DO J=1,N
              DO I=1,N
                 CUBE(I,J,K) = 0.0
              ENDDO
           ENDDO
        ENDDO

#ifdef USE_MPI
C ------------------------- MPI ONLY CODE -----------------------------
        D_KLP_LOC = 0.0D0
        D_KLS_LOC = 0.0D0
        D_ABA_LOC = 0.0D0 
        D_SUS_LOC = 0.0D0
        D_SSQ_LOC = 0.0D0

        CALL SETPART(NANG, PSIZE, NBASE)
        NANGLOC = PSIZE(MYPID+1)

        ALLOCATE(PRJBUF  (N,N,PSIZE(1)),
     &           PRJLOC  (N,N,NANGLOC),
     &           ANG_LOC (3,NANGLOC),
     &           CUBE_LOC(N,N,N), 
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MWANT = N*N*PSIZE(1) + N*N*NANGLOC + 3*NANGLOC + N*N*N
           CALL ERRT(46,'BPCG_2; PRJBUF, PRJLOC..',MWANT)
           STOP
        ENDIF

C       ZERO LOCAL ARRAYS
        CUBE_LOC = 0.0
        ANG_LOC  = 0.0

        DO IPROC = 1, NPROCS
           NLOC = PSIZE(IPROC)

C          READ A SUBSET OF IMAGES (ONLY ONE PROCESSOR READS)

           DO K=1,NLOC
              KGLB = K + NBASE(IPROC)
              NLET = 0
              CALL  FILGET(FINPAT,FILNAM,NLET,ILIST(KGLB),IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

              MAXIM = 0
              CALL OPFILEC(0,.FALSE.,FILNAM,LUNPROJ,'O',IFORM,
     &                     LSAM,LROW,NSL,
     &                     MAXIM,'DUMMY',.FALSE.,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN
C
              IF (MYPID <= 0) THEN
C                ONLY PROCESS WITH MYPID = 0 NEEDS TO READ
                 DO K2=1,N
                   CALL REDLIN1P(LUNPROJ,PRJBUF(1,K2,K),N,K2)
                 ENDDO
              ENDIF
              CLOSE(LUNPROJ)
           ENDDO

C          DISTRIBUTE IMAGES
           IF (IPROC > 1) THEN
              IF (MYPID == 0) THEN
                !write(6,*) ' send_mpi; iproc,mypid:',iproc,mypid,nloc

                 CALL SEND_MPI('BPCG_2','PRJBUF', PRJBUF, N*N*NLOC, 
     &                         'R',IPROC-1,IPROC-1, ICOMM)

              ELSE IF (MYPID .EQ. IPROC-1) THEN

                 !write(6,*) ' recv_mpi; iproc,mypid:',iproc,mypid,nloc
                 CALL RECV_MPI('BPCG_2','PRJLOC', PRJLOC, N*N*NLOC, 
     &                          'R', 0,MPI_ANY_TAG, ICOMM)              

c                call mpi_recv(prjloc,  n*n*nloc,    mpi_real,
c     &                         0,       mpi_any_tag, icomm,
c     &                         mpistat, mpierr)
                 !write(6,*) ' recv_mpi; ok:         ',iproc,mypid,nloc
              ENDIF

           ELSE IF (MYPID .EQ. 0) THEN
              CALL SCOPY(N*N*NLOC,PRJBUF,1,PRJLOC,1)
           ENDIF
        ENDDO

        DO K = 1, NANGLOC
           KGLB = K + NBASE(MYPID+1)

C          ORDER IN DOCUMENT FILE IS PSI, THETA, PHI AND ANGLES ARE IN
C          DEGREES!
C          IN ANG ARRAY IT IS OTHER WAY AROUND

           ITMP         = ILIST(KGLB)
           ANG_LOC(3,K) = ANGBUF(2,ITMP)
           ANG_LOC(2,K) = ANGBUF(3,ITMP)
           ANG_LOC(1,K) = ANGBUF(4,ITMP)

c          write(6,91) kglb,itmp,(ang_loc(j,k),j=3,1,-1)
 91        format('  kglb:',i5,'  proj:',i6,
     &               '; psi:',f6.1,' theta:',f6.1,' phi:',f6.1)

C          ESTIMATE AVERAGE OUTSIDE THE CIRCLE?
           CALL ASTASQ(PRJLOC(1,1,K), N, RI, D_ABA_LOC, D_KLP_LOC, 
     &                 D_SUS_LOC, D_SSQ_LOC, D_KLS_LOC)
 
C          CREATE ROTATION MATRIX: DM
           CALL CANG(ANG_LOC(1,K),ANG_LOC(2,K),ANG_LOC(3,K),
     &              .FALSE.,ADUM,DM(1,KGLB))

C          BACKPROJECT: PRJLOC  INTO: CUBE_LOC
           CALL RPRQ(N,PRJLOC(1,1,K),CUBE_LOC,IPCUBE,NN,
     &               DM(1,KGLB), LDP,LDPNM,IRTFLG)
        ENDDO

        IF (ALLOCATED(PRJBUF)) DEALLOCATE(PRJBUF) 
        IF (ALLOCATED(PRJLOC)) DEALLOCATE(PRJLOC) 

C       GATHER ANG_LOC INTO ANG

        DO J = 1, NPROCS
           NBASE(J) = 3 * NBASE(J)
           PSIZE(J) = 3 * PSIZE(J)
        ENDDO
        CALL MPI_ALLGATHERV(ANG_LOC,  3*NANGLOC, MPI_REAL,
     &                      ANG,      PSIZE,     NBASE,
     &                      MPI_REAL, ICOMM,     IERR)    

        DO K = 1, NANG
           IF (VERBOSE) THEN
              IF (MYPID == 0) WRITE(NOUT,3331) K,(ANG(J,K),J=3,1,-1)
3331          FORMAT('  IMAGE #:',I7,
     &               '  PSI:', F6.1,
     &               ' THETA:',F6.1,
     &               ' PHI:',  F6.1)
           ENDIF
        END DO

        NNN  = N*N*N
        IERR = 0
        CALL ALLREDUCE_MPI('REPRCG','CUBE', CUBE_LOC,CUBE,
     &                         NNN, 'R','S',ICOMM)
        CALL ALLREDUCE_MPI('REPRCG','D_ABA', D_ABA_LOC,D_ABA,
     &                           1, 'D','S',ICOMM)
        CALL ALLREDUCE_MPI('REPRCG','D_KLP', D_KLP_LOC,D_KLP,
     &                           1, 'D','S',ICOMM)
        CALL ALLREDUCE_MPI('REPRCG','D_SSQ', D_SSQ_LOC,D_SSQ,
     &                           1, 'D','S',ICOMM)
        CALL ALLREDUCE_MPI('REPRCG','D_SUS', D_SUS_LOC,D_SUS,
     &                           1, 'D','S',ICOMM)
        CALL ALLREDUCE_MPI('REPRCG','D_KLS', D_KLS_LOC,D_KLS,
     &                           1, 'D','S',ICOMM)

        IF (ALLOCATED(CUBE_LOC)) DEALLOCATE(CUBE_LOC)
        IF (ALLOCATED(ANG_LOC))  DEALLOCATE(ANG_LOC)
        IF (ALLOCATED(PSIZE))    DEALLOCATE(PSIZE)
        IF (ALLOCATED(NBASE))    DEALLOCATE(NBASE)


C       --------------- END OF MPI CODE -----------------------------
#else
C       --------------- START NON-MPI CODE --------------------------

        K   = 1
        NNN = N*N*N

        DO      ! LOOP OVER ALL INPUT IMAGES

           ! LOAD NEXT IMAGE (SQUARE)
           CALL REDVOL(LUNPROJ,N,N,1,1,PROJ,IRTFLG)

C          LOAD ANGLES FOR THIS IMAGE
C          ORDER IN DOCUMENT FILE IS PSI, THETA, PHI AND ANGLES ARE IN 
C          DEGREES! IN ANG ARRAY IT IS OTHER WAY AROUND

           ITMP     = ILIST(K)
           ICOUNT   = ANGBUF(1,ITMP)

           IF (ICOUNT <= 0) THEN
C             MISSING KEY
              CALL ERRT(102,'MISSING ANGLES FOR IMAGE',ITMP)
              IRTFLG = 1
              GOTO 9999
           ENDIF

           ANG(3,K) = ANGBUF(2,ITMP)  ! NOTE ORDER REVERSAL!
           ANG(2,K) = ANGBUF(3,ITMP)
           ANG(1,K) = ANGBUF(4,ITMP)
           IF (VERBOSE) THEN
              WRITE(NOUT,3331) K,(ANG(J,K),J=3,1,-1)
3331          FORMAT('  IMAGE #:',I7,
     &               '  PSI:',    F6.1,
     &               '  THETA:',  F6.1,
     &               '  PHI:',    F6.1)
           ENDIF

C          CREATE ROTATION MATRIX: DM
           CALL CANG(ANG(1,K),ANG(2,K),ANG(3,K),.FALSE.,ADUM,DM(1,K))

C          ESTIMATE AVERAGE OUTSIDE THE CIRCLE?
           CALL ASTASQ(PROJ,N,RI,D_ABA,D_KLP,D_SUS,D_SSQ,D_KLS)

           IF (FBS_WANTED) THEN

C             PAD PROJECTION TO: PROJPAD
c$omp         parallel do private(iy,iloc,ix)
              DO IY = 1,N
                 ILOC = (IY-1)*N + 1 
                 DO IX = 1,N 
                    PROJPAD(IX,IY) = PROJ(ILOC)
                    ILOC           = ILOC + 1
                 ENDDO
              ENDDO

C             CALCULATE PROJECTION DERIVATIVES USING FFT
C             IN FBS2_PREP PROJPAD DOES NOT RETURN USEFULL VALUE!
              CALL FBS2_PREP(PROJPAD, XDER,YDER,XYDER, 
     &                       NXLD,N,N, IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999
   
C             BACKPROJECT FROM PLANE: PROJ TO VOLUME: CUBE
              CALL BCKPJ_FBS(CUBE,NNN,DM(1,K),
     &                  PROJ, XDER,YDER,XYDER,
     &                  NXLD,N,IPCUBE,NN, LDP,LDPNM)
 
           ELSE
C             COMPUTE IMAGE DIFFERENCES FOR INTERP. PIXELS
c$omp         parallel do private(i,pt)
              DO I = 1,N*N - N - 1
                 PT         = PROJ(I)
                 PROJT(1,I) = PT
                 PROJT(2,I) = PROJ(I+N)   - PT
                 PROJT(3,I) = PROJ(I+1)   - PT
 
C                may 2013 unfathom-able gfort bug requires -fwrapv
                 PROJT(4,I) = PROJ(I+N+1) - PROJ(I+1) - PROJT(2,I)
 
              ENDDO

C             BACKPROJECTION FROM PROJT INTO: CUBE
              CALL BCKPJ_LIN(CUBE,NNN,DM(1,K),
     &                     PROJT,N,IPCUBE,NN,LDP,LDPNM)
           ENDIF
  
C          OPEN NEXT INPUT PROJECTION
           CALL NEXTFILE(K,       ILIST, 
     &                   FOUROK,  0,
     &                   NANG,    MAXIM,   
     &                   LUNPROJ, 0, 
     &                   FINPAT,  'O',
     &                   IMGNUM,  IRTFLG) 
           IF (IRTFLG ==  -1) EXIT      !  END OF INPUT STACK
           IF (IRTFLG .NE. 0) RETURN
 
        ENDDO
        IRTFLG = 0

9999    IF (ALLOCATED(PROJPAD)) DEALLOCATE (PROJPAD)
        IF (ALLOCATED(XDER))    DEALLOCATE (XDER)
        IF (ALLOCATED(YDER))    DEALLOCATE (YDER)
        IF (ALLOCATED(XYDER))   DEALLOCATE (XYDER)


C       --------------- END OF NON-MPI CODE -------------------------
#endif


C       BACKGROUND AVERAGE CORRECTION FOR VOLUME
        D_ABA = D_ABA / D_KLP
        BNORM = 0.0
        QT    = D_ABA * NANG

c$omp   parallel do private(kn,j,k,i),reduction(+:bnorm)
        DO KN=1,NN
           J = IPCUBE(4,KN)
           K = IPCUBE(5,KN)
           DO I=IPCUBE(3,KN), IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)
              CUBE(I,J,K) = CUBE(I,J,K) - QT
              BNORM     = BNORM + CUBE(I,J,K) * CUBE(I,J,K)
           ENDDO
        ENDDO

C       CALCULATE SUM OF SQUARES OF PROJECTIONS WITHIN CIRCLES AFTER
C       SUBTRACTION OF THE AVERAGE
        D_SSQ = D_SSQ - D_ABA * (2 * D_SUS - D_ABA * D_KLS)
        CHI2  = D_SSQ

C       PRINT STATISTICS

        KLP_8 = D_KLP

        IF (MYPID <= 0) THEN
           WRITE(NOUT,2044) KLP_8,D_ABA,BNORM,CHI2
2044       FORMAT(/,'  TOTAL POINTS IN PROJECTIONS:',I10,/,
     &              '  AVERAGE OUTSIDE THE WINDOW: ',1PE10.3,/,
     &              '  SQUARED BP:                 ',1PE10.3,/,
     &              '  CHI2:                       ',1PE10.3)
        ENDIF

        !write(6,*) ' k,sum-cube:',k,sum(cube)

        END

C**************************** BPCG_3_LIN ******************************

        SUBROUTINE BPCG_3_LIN(CB,BCKN,NXLD,N,IPCUBE,NN,DM,LB,NANG,RI,
     &                      NUMTH,BNORM,CHI2,
     &                      ERRM,CHIM,MAXIT,MODE,ANOISE, 
     &                      LDP,LDPNM,IRTFLG)

C       NUMTH = NUMTHREADS() FOR MP, OTHERWISE=1.

        INCLUDE 'CMBLOCK.INC'

        REAL                :: CB(N,N,N),BCKN(N,N,N)
        INTEGER             :: NXLD,N,NN
        INTEGER             :: IPCUBE(5,NN),LB(NANG)
        INTEGER             :: NANG
        REAL                :: RI
        INTEGER             :: NUMTH
        REAL                :: CHI2,ERRM,CHIM
        INTEGER             :: MAXIT,MODE
        REAL                :: ANOISE, DELSQ
        INTEGER             :: LDP,LDPNM,IRTFLG

        REAL                :: DM(9,NANG),PROJ(N*N,NUMTH)
        REAL                :: PROJT(4,N*N)
        REAL,ALLOCATABLE    :: BCKE(:,:,:), BCKP(:,:,:)

#ifdef USE_MPI
        include 'mpif.h'
        INTEGER             :: MPISTAT(MPI_STATUS_SIZE)
        INTEGER,ALLOCATABLE :: PSIZE(:), NBASE(:)
        INTEGER,ALLOCATABLE :: LBLOC(:)
        REAL,   ALLOCATABLE :: DMLOC(:,:)
        REAL,   ALLOCATABLE :: BCKE_SUM(:,:,:)
        REAL,   ALLOCATABLE :: RESID(:,:,:)

        ICOMM = MPI_COMM_WORLD
        CALL MPI_COMM_RANK(ICOMM, MYPID , MPIERR)
        CALL MPI_COMM_SIZE(ICOMM, NPROCS, MPIERR)

        ALLOCATE(PSIZE(NPROCS), 
     &           BCKE_SUM(N,N,N), 
     &           NBASE(NPROCS), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MWANT = 2 * NPROCS + N*N*N
           CALL ERRT(46,'BPCG_3_LIN; NBASE,...',MWANT)
           RETURN
        ENDIF
        CALL SETPART(NANG, PSIZE, NBASE)
        NANGLOC = PSIZE(MYPID+1)

#else
        MYPID = -1
#endif
        
        NNN = N*N*N
        ALLOCATE (BCKE(N,N,N), 
     &            BCKP(N,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           MWANT = 2 * NNN
           CALL ERRT(46,'BPCG_3_LIN; BCKE,BCKP',MWANT)
           GOTO 9999
        ENDIF

c$omp   parallel do private(k,j,i)
        DO K=1,N
           DO J=1,N
              DO I=1,N
                 BCKN(I,J,K)     = 0.0
                 BCKP(I,J,K)     = 0.0
#ifdef USE_MPI
                 BCKE_SUM(I,J,K) = 0.0
#endif
              ENDDO
           ENDDO
        ENDDO

        ERR  = 1.0
        ITER = 0
        Q    = 0.0

C       ------------------ START OF:  MPI CODE ------------------------
#ifdef USE_MPI

C       DISTRIBUTE LB AND DM

        ALLOCATE(LBLOC(NANGLOC),
     &           DMLOC(9, NANGLOC), 
     &           STAT = IRTFLAG)
        IF (IRTFLG .NE. 0) THEN
           MWANT = 10 * NANGLOC
           CALL ERRT(46,'BPCG_3_LIN; DMLOC',MWANT)
           STOP
        ENDIF

        IDMTAG  = 1
        LBTAG   = 2
        MASTER  = 0

        IF (MYPID .EQ. 0) THEN
C          PROCESS 0 DISTRIBUTES DM() AND LB() TO OTHER PROCESSORS

           DO IP = 2, NPROCS
              IBEGIN = NBASE(IP) + 1
              NLOC   = PSIZE(IP) 

              CALL SEND_MPI('BPCG_3_LIN','DM', DM(1,IBEGIN), 9*NLOC, 
     &                      'R',IP-1,IDMTAG,ICOMM)

              CALL SEND_MPI('BPCG_3_LIN','LB', LB(IBEGIN), NLOC, 
     &                      'I',IP-1,LBTAG,ICOMM)

           ENDDO

C          COPY DATA INTO LOCAL ARRAYS: DMLOC & LBLOC
           CALL SCOPY(9*NANGLOC, DM, 1, DMLOC, 1)
           CALL SCOPY(NANGLOC  , LB, 1, LBLOC, 1)
        ELSE

C          === SLAVES RECEIVE DATA FROM THE MASTER ===

 
           CALL RECV_MPI('BPCG_3_LIN','DMLOC', DMLOC, 9*NANGLOC, 'R',
     &                   MASTER,IDMTAG, ICOMM)
 
           CALL RECV_MPI('BPCG_3_LIN','LBLOC', LBLOC, NANGLOC, 'I',
     &                   MASTER,LBTAG, ICOMM)


        ENDIF

#ifdef MPI_DEBUG
        write(6,*)' BPCG_3_LIN: data distributed, mypid:',mypid
#endif

        DO ITER=1,MAXIT

           DELSQ = 0.0
c$omp      parallel do private(kn,j,k,i),reduction(+:delsq)
           DO KN=1,NN
              J = IPCUBE(4,KN)
              K = IPCUBE(5,KN)
              DO I=IPCUBE(3,KN),IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)
                 DELSQ = DELSQ + CB(I,J,K) * CB(I,J,K)
              ENDDO
           ENDDO
           Q = Q * DELSQ

c$omp      parallel do private(kn,j,k,i)
           DO KN=1,NN
              J = IPCUBE(4,KN)
              K = IPCUBE(5,KN)
              DO I=IPCUBE(3,KN),IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)
                 BCKP(I,J,K) = CB(I,J,K) - Q * BCKP(I,J,K)
              ENDDO
           ENDDO
           Q = -1.0 / DELSQ

c$omp      parallel do private(k,j,i)
           DO K=1,N
              DO J=1,N
                 DO I=1,N
                    BCKE(I,J,K)     = 0.0
                    BCKE_SUM(I,J,K) = 0.0
                 ENDDO
              ENDDO
           ENDDO

C          BCKP -> PROJ -> BCKE

C          LOOP OVER PROJECTIONS

           DO K=1,NANGLOC
              L_EN  = MIN(NANGLOC,K+NUMTH-1)
              L_NUM = MIN(NUMTH,NANGLOC-K+1)

c$omp         parallel do private(i,j)
              DO J=1,L_NUM
                 DO I=1,N*N
                    PROJ(I,J) = 0.0
                 ENDDO
              ENDDO

C             PROJECT
              CALL PRJCQ_N(BCKP,NNN,DMLOC(1,K),
     &                     PROJ,N,IPCUBE,NN, LDP,LDPNM)

C             HERE BCKPJ ITSELF IS MP  
              DO  L_TH=1,L_NUM
C                MULTIPLY PROJECTIONS BY THEIR WEIGHTS
                 IF (LBLOC(K+L_TH-1) .GT. 1)  THEN
c$omp               parallel do private(i)
                    DO  I=1,N*N
                       PROJ(I,L_TH) = PROJ(I,L_TH) * LBLOC(K+L_TH-1)
                    ENDDO
                 ENDIF

c$omp            parallel do private(i,pt)
                 DO I = 1,N*N - N - 1
                    PT         = PROJ(I, L_TH)
                    PROJT(1,I) = PT
                    PROJT(2,I) = PROJ(I+N,  L_TH) - PT
                    PROJT(3,I) = PROJ(I+1,  L_TH) - PT 

C                   may 2013 unfathom-able gfort bug requires -fwrapv
                    PROJT(4,I) = PROJ(I+N+1,L_TH) - PROJ(I+1,L_TH) - 
     &                           PROJT(2,I)

                 ENDDO

C                BACKPROJECT
                 CALL BCKPJ_LIN(BCKE,NNN,DMLOC(1,K+L_TH-1),
     &                        PROJT,N,IPCUBE,NN, LDP,LDPNM)
              ENDDO
           ENDDO
C
           CALL ALLREDUCE_MPI('BPCG_3_LIN','BCKE_SUM',BCKE,BCKE_SUM,
     &                           NNN, 'R','S',ICOMM)


C          === IGNORE MODE > 0 FOR THE TIME BEING ===

           IF (MODE .EQ. 1)  THEN
              CALL FIXEDGE1(BCKP,NNN,BCKP,N,IPCUBE,NN)      
c             CALL BFIRSTS(BCKE,BCKP,N,ANOISE,IPCUBE,NN,RI) sep 01 al
              CALL BFIRSTS(BCKE_SUM,BCKP,N,ANOISE,IPCUBE,NN)

           ELSEIF(MODE.EQ.2)  THEN
              CALL FIXEDGE2(BCKP,NNN,BCKP,N,IPCUBE,NN)      
              CALL BSECOND(BCKE_SUM,BCKP,N,ANOISE,IPCUBE,NN)

           ELSEIF(MODE.EQ.3)  THEN
              CALL FIXEDGE3(BCKP,NNN,BCKP,N,IPCUBE,NN)      
              CALL BTHIRD(BCKE_SUM,BCKP,N,ANOISE,IPCUBE,NN)
           ENDIF

           AKDEN = 0.0

c$omp      parallel do private(kn,j,k,i),reduction(+:akden)
           DO KN=1,NN
              J = IPCUBE(4,KN)
              K = IPCUBE(5,KN)
              DO I=IPCUBE(3,KN),
     &             IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)
                 AKDEN = AKDEN+BCKP(I,J,K) * BCKE_SUM(I,J,K)
              ENDDO
           ENDDO

           P=DELSQ / AKDEN

c$omp      parallel do private(kn,j,k,i)
           DO KN=1,NN
              J = IPCUBE(4,KN)
              K = IPCUBE(5,KN)
              DO I=IPCUBE(3,KN),
     &           IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)
                 BCKN(I,J,K) = BCKN(I,J,K)+P*BCKP(I,J,K)
                 CB(I,J,K)   = CB(I,J,K)-P*BCKE_sum(I,J,K)
              ENDDO
           ENDDO
 
           ERR  = DELSQ / BNORM
           CHI2 = CHI2 - P * DELSQ
           IF (MYPID .EQ. 0) WRITE(NOUT,2041) ITER,ERR,CHI2
2041       FORMAT('  ITERATION:', I3,
     &            '  DIFFERENCE:',1PE12.4,
     &            '  CHISQ:',1PE12.4)

C          CHECK STOPPING CRITERIA
           IF (ABS(ERR) <= ERRM .OR. ABS(CHI2) <= CHIM) GOTO 9999            
        ENDDO


C       END OF: -------------- USE_MPI ------------------------------

#else

        DO  ITER=1,MAXIT

           DELSQ = 0.0
c$omp      parallel do private(kn,j,k,i),reduction(+:delsq)
           DO KN=1,NN
              J = IPCUBE(4,KN)
              K = IPCUBE(5,KN)
              DO I=IPCUBE(3,KN),IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)
                 DELSQ = DELSQ+CB(I,J,K)*CB(I,J,K)
              ENDDO
           ENDDO
           Q = Q * DELSQ

c$omp      parallel do private(kn,j,k,i)
           DO KN=1,NN
              J = IPCUBE(4,KN)
              K = IPCUBE(5,KN)
              DO I=IPCUBE(3,KN),IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)
                 BCKP(I,J,K) = CB(I,J,K) - Q * BCKP(I,J,K)
              ENDDO
           ENDDO
           Q = -1.0 / DELSQ

c$omp      parallel do private(k,j,i)
           DO K=1,N
              DO J=1,N
                 DO I=1,N
                    BCKE(I,J,K) = 0.0
                 ENDDO
              ENDDO
           ENDDO

C          BCKP -> PROJ -> BCKE
C          LOOP OVER PROJECTIONS
           DO K=1,NANG,NUMTH
              L_EN  = MIN(NANG,K+NUMTH-1)
              L_NUM = MIN(NUMTH,NANG-K+1)

c$omp         parallel do private(i,j)
              DO J=1,L_NUM
                 DO I=1,N*N
                    PROJ(I,J) = 0.0
                 ENDDO
              ENDDO

c$omp      parallel do private(l_th),shared(nnn,n,nn),schedule(static)
              DO  L_TH=1,L_NUM
C                PROJECTS BCKP INTO PROJ  
                 CALL PRJCQ_N(BCKP,NNN,DM(1,K+L_TH-1),
     &                       PROJ(1,L_TH),N,IPCUBE,NN, LDP,LDPNM)
              ENDDO

C             HERE BCKPJ ITSELF IS MP  
              DO  L_TH=1,L_NUM
C                MULTIPLY PROJECTIONS BY THEIR WEIGHTS
                 IF (LB(K+L_TH-1) > 1)  THEN
c$omp               parallel do private(i)
                    DO  I=1,N*N
                       PROJ(I,L_TH) = PROJ(I,L_TH)*LB(K+L_TH-1)
                    ENDDO
                 ENDIF

c$omp            parallel do private(i,pt)
                 DO I = 1,N*N - N - 1
                    PT         = PROJ(I, L_TH)
                    PROJT(1,I) = PT
                    PROJT(2,I) = PROJ(I+N,  L_TH) - PT
                    PROJT(3,I) = PROJ(I+1,  L_TH) - PT 

C                   may 2013 unfathom-able gfort bug requires -fwrapv
                    PROJT(4,I) = PROJ(I+N+1,L_TH) - PROJ(I+1,L_TH) - 
     &                           PROJT(2,I)
                 ENDDO

C                BACKPROJECT FROM PLANE: PROJT TO VOLUME: BCKE
                 CALL BCKPJ_LIN(BCKE,NNN,DM(1,K+L_TH-1),PROJT,N,
     &                          IPCUBE,NN, LDP,LDPNM)
              ENDDO
           ENDDO

           IF (MODE == 1) THEN
              CALL FIXEDGE1(BCKP,NNN,BCKP,N,IPCUBE,NN)      
              CALL BFIRSTS(BCKE,BCKP,N,ANOISE,IPCUBE,NN)

           ELSEIF(MODE == 2) THEN
              CALL FIXEDGE2(BCKP,NNN,BCKP,N,IPCUBE,NN)      
              CALL BSECOND(BCKE,BCKP,N,ANOISE,IPCUBE,NN)

           ELSEIF (MODE == 3) THEN
              CALL FIXEDGE3(BCKP,NNN,BCKP,N,IPCUBE,NN)      
              CALL BTHIRD(BCKE,BCKP,N,ANOISE,IPCUBE,NN)
           ENDIF

           AKDEN = 0.0

c$omp      parallel do private(kn,j,k,i),reduction(+:akden)
           DO KN=1,NN
              J = IPCUBE(4,KN)
              K = IPCUBE(5,KN)
              DO I=IPCUBE(3,KN),
     &             IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)
                 AKDEN = AKDEN + BCKP(I,J,K) * BCKE(I,J,K)
              ENDDO
           ENDDO
           P = DELSQ / AKDEN

c$omp      parallel do private(kn,j,k,i)
           DO KN=1,NN
              J = IPCUBE(4,KN)
              K = IPCUBE(5,KN)
              DO I=IPCUBE(3,KN),
     &           IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)
                 BCKN(I,J,K) = BCKN(I,J,K)+P*BCKP(I,J,K)
                 CB(I,J,K)   = CB(I,J,K)-P*BCKE(I,J,K)
              ENDDO
           ENDDO

           ERR  = DELSQ / BNORM
           CHI2 = CHI2 - P * DELSQ

           WRITE(NOUT,2041) ITER,ERR,CHI2
2041       FORMAT('  ITERATION: 'I3,
     &            '  DIFFERENCE:',1PE12.4,
     &            '  CHISQ:',1PE12.4)

C          CHECK STOPPING CRITERIA
           IF (ABS(ERR) <= ERRM .OR. ABS(CHI2) <= CHIM) THEN  
              WRITE(NOUT,*) '  '  
              GOTO 9999            
           ENDIF
       ENDDO


C      END OF: -------------- NON-MPI ----------------------------
#endif

C      NUMBER OF ITERATIONS EXCEEDED
 
9999   IF (ALLOCATED(BCKE))     DEALLOCATE (BCKE)
       IF (ALLOCATED(BCKP))     DEALLOCATE (BCKP)
#ifdef USE_MPI
       IF (ALLOCATED(BCKE_SUM)) DEALLOCATE(BCKE_SUM)
#endif

       END



C*************************  ASTASQ  ********************************

        SUBROUTINE  ASTASQ(X,N,RI,D_ABA,D_KLP,D_SUS,D_SSQ,D_KLS)

        REAL              :: X(N,N)
        DOUBLE PRECISION  :: D_ABA,D_SUS,D_SSQ
        DOUBLE PRECISION  :: D_KLP,D_KLS

C       ESTIMATE AVERAGE OUTSIDE THE CIRCLE. 
C       RETURNS: D_ABA, D_KLP, D_SUS, D_SSQ, D_KLS

        R  = RI * RI
        NC = N / 2+1

        DO J=1,N
           T  = J - NC
           XX = T*T

           DO I=1,N
              T=I - NC
              IF ( (XX + T*T) > R)    THEN
C                OUTSIDE THE CIRCLE MASK. 
                 D_ABA  = D_ABA  + DBLE(X(I,J))
                 D_KLP = D_KLP + 1
              ELSE
C                INSIDE THE CIRCLE MASK. 
                 D_SSQ  = D_SSQ + X(I,J) * DBLE(X(I,J))
                 D_SUS  = D_SUS + X(I,J)
                 D_KLS = D_KLS + 1
              ENDIF
           ENDDO
        ENDDO
        END

C*************************  FIXEDGE1  ********************************

        SUBROUTINE  FIXEDGE1(BCKP,NNN,BCK3,N,IPCUBE,NN)

        REAL    :: BCKP(NNN),BCK3(N,N,N) 
        INTEGER :: IPCUBE(5,NN)

C       PUT ZEROS OUTSIDE
        NT = 1
        DO I=1,NNN
           IF (NT > NN)  THEN
              BCKP(I) = 0.0
           ELSEIF (I < IPCUBE(1,NT))  THEN
              BCKP(I)= 0.0
           ELSEIF(I == IPCUBE(2,NT))  THEN
              NT = NT+1
           ENDIF
        ENDDO

C       ADD PIXELS ON THE EDGE
C       FIX THE EDGES IN BCKP
        DO  K=1,N
           DO J=1,N
              DO  I=1,N-1
                 IF (BCK3(I+1,J,K) .NE. 0.0) THEN
                    BCK3(I,J,K) = BCK3(I+1,J,K)
                    EXIT
                 ENDIF
              ENDDO

              DO  I=N,2,-1
                 IF (BCK3(I-1,J,K) .NE. 0.0) THEN
                    BCK3(I,J,K) =BCK3(I-1,J,K)
                    EXIT
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        DO  K=1,N
           DO I=1,N
              DO  J=1,N-1
                 IF (BCK3(I,J+1,K) .NE. 0.0) THEN
                    BCK3(I,J,K) = BCK3(I,J+1,K)
                    EXIT
                 ENDIF
              ENDDO

              DO  J=N,2,-1
                 IF (BCK3(I,J-1,K) .NE. 0.0) THEN
                    BCK3(I,J,K) = BCK3(I,J-1,K)
                    EXIT
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        DO  J=1,N
           DO I=1,N
              DO  K=1,N-1
                 IF (BCK3(I,J,K+1) .NE. 0.0) THEN
                    BCK3(I,J,K) = BCK3(I,J,K+1)
                    EXIT
                 ENDIF
              ENDDO

              DO  K=N,2,-1
                 IF (BCK3(I,J,K-1) .NE. 0.0) THEN
                    BCK3(I,J,K) = BCK3(I,J,K-1)
                    EXIT
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        END

C************************  FIXEDGE2 *********************************

        SUBROUTINE  FIXEDGE2(BCKP,NNN,BCK3,N,IPCUBE,NN)

        REAL    :: BCKP(NNN),BCK3(N,N,N)
        INTEGER :: IPCUBE(5,NN)

C       PUT ZEROS OUTSIDE

        NT = 1
        DO    I=1,NNN
           IF (NT .GT. NN)  THEN
              BCKP(I) = 0.0
           ELSEIF (I .LT. IPCUBE(1,NT))  THEN
              BCKP(I) = 0.0
           ELSEIF (I .EQ. IPCUBE(2,NT))  THEN
              NT = NT+1
           ENDIF
        ENDDO

C       ADD PIXELS ON THE EDGE
C       FIX THE EDGES IN BCKP
        DO  K=1,N
           DO J=1,N
              DO  I=2,N-1
                 IF (BCK3(I+1,J,K) .NE. 0.0) THEN
                     BCK3(I,J,K)   = BCK3(I+1,J,K)
                     BCK3(I-1,J,K) = BCK3(I+1,J,K)
                     EXIT
                 ENDIF
              ENDDO

              DO  I=N-1,2,-1
                 IF (BCK3(I-1,J,K) .NE. 0.0) THEN
                     BCK3(I,J,K)   = BCK3(I-1,J,K)
                     BCK3(I+1,J,K) = BCK3(I-1,J,K)
                     EXIT
                 ENDIF
             ENDDO
           ENDDO
        ENDDO
        DO  K=1,N
           DO I=1,N
              DO  J=2,N-1
                 IF(BCK3(I,J+1,K).NE.0.0) THEN
                    BCK3(I,J,K)=BCK3(I,J+1,K)
                    BCK3(I,J-1,K)=BCK3(I,J+1,K)
                    EXIT
                 ENDIF
              ENDDO

              DO  J=N-1,2,-1
                 IF(BCK3(I,J-1,K).NE.0.0) THEN
                    BCK3(I,J,K)   = BCK3(I,J-1,K)
                    BCK3(I,J+1,K) = BCK3(I,J-1,K)
                    EXIT
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        DO  J=1,N
           DO I=1,N
              DO  K=2,N-1
                 IF(BCK3(I,J,K+1).NE.0.0) THEN
                    BCK3(I,J,K)   = BCK3(I,J,K+1)
                    BCK3(I,J,K-1) = BCK3(I,J,K+1)
                    EXIT
                 ENDIF
              ENDDO

              DO  K=N-1,2,-1
                 IF(BCK3(I,J,K-1).NE.0.0) THEN
                    BCK3(I,J,K)   = BCK3(I,J,K-1)
                    BCK3(I,J,K+1) = BCK3(I,J,K-1)
                   EXIT
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        END

C************************  FIXEDGE3 ********************************

        SUBROUTINE  FIXEDGE3(BCKP,NNN,BCK3,N,IPCUBE,NN)

        REAL    :: BCKP(NNN),BCK3(N,N,N)
        INTEGER :: IPCUBE(5,NN)

C       PUT ZEROS OUTSIDE

        NT = 1
        DO I=1,NNN
           IF (NT.GT.NN)  THEN
              BCKP(I) = 0.0
           ELSEIF (I .LT. IPCUBE(1,NT))  THEN
              BCKP(I) = 0.0
           ELSEIF (I .EQ. IPCUBE(2,NT))  THEN
              N T= NT+1
           ENDIF
        ENDDO

C       ADD PIXELS ON THE EDGE
C       FIX THE EDGES IN BCKP
        DO  K=1,N
           DO J=1,N
              DO  I=3,N-1
                 IF(BCK3(I+1,J,K) .NE. 0.0) THEN
                    BCK3(I,J,K)   = BCK3(I+1,J,K)
                    BCK3(I-1,J,K) = BCK3(I+1,J,K)
                    BCK3(I-2,J,K) = BCK3(I+1,J,K)
                    EXIT
                 ENDIF
              ENDDO

              DO  I=N-2,2,-1
                 IF(BCK3(I-1,J,K) .NE. 0.0) THEN
                    BCK3(I,J,K)   = BCK3(I-1,J,K)
                    BCK3(I+1,J,K) = BCK3(I-1,J,K)
                    BCK3(I+2,J,K) = BCK3(I-1,J,K)
                    EXIT
                 ENDIF
              ENDDO
           ENDDO
        ENDDO

        DO  K=1,N
           DO I=1,N
              DO  J=3,N-1
                 IF(BCK3(I,J+1,K) .NE. 0.0) THEN
                    BCK3(I,J,K)   = BCK3(I,J+1,K)
                    BCK3(I,J-1,K) = BCK3(I,J+1,K)
                    BCK3(I,J-2,K) = BCK3(I,J+1,K)
                    EXIT
                 ENDIF
              ENDDO

              DO  J=N-2,2,-1
                 IF(BCK3(I,J-1,K) .NE. 0.0) THEN
                    BCK3(I,J,K)   = BCK3(I,J-1,K)
                    BCK3(I,J+1,K) = BCK3(I,J-1,K)
                    BCK3(I,J+2,K) = BCK3(I,J-1,K)
                   EXIT
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        DO  J=1,N
           DO I=1,N
              DO  K=3,N-1
                 IF(BCK3(I,J,K+1) .NE. 0.0) THEN
                    BCK3(I,J,K)   = BCK3(I,J,K+1)
                    BCK3(I,J,K-1) = BCK3(I,J,K+1)
                    BCK3(I,J,K-2) = BCK3(I,J,K+1)
                    EXIT
                 ENDIF
              ENDDO

              DO  K=N-2,2,-1
                 IF(BCK3(I,J,K-1) .NE. 0.0) THEN
                    BCK3(I,J,K)   = BCK3(I,J,K-1)
                    BCK3(I,J,K+1) = BCK3(I,J,K-1)
                    BCK3(I,J,K+2) = BCK3(I,J,K-1)
                    EXIT
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        END

C       ********************  BFIRSTS   ******************************

        SUBROUTINE  BFIRSTS(BCKE,BCKP,N,ANOISE,IPCUBE,NN)

        REAL      ::  BCKE(N,N,N),BCKP(N,N,N)
        INTEGER   ::  IPCUBE(5,NN)

c$omp   parallel do private(kn,j,k,i)
        DO KN=1,NN
           J = IPCUBE(4,KN)
           K = IPCUBE(5,KN)

           DO I=IPCUBE(3,KN),IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)

              BCKE(I,J,K) = BCKE(I,J,K) + ANOISE*(6*BCKP(I,J,K)-(
     &           + BCKP(I+1,J,K) + BCKP(I,J+1,K) + BCKP(I,J,K+1)
     &           + BCKP(I-1,J,K) + BCKP(I,J-1,K) + BCKP(I,J,K-1)))
           ENDDO
        ENDDO
        END

C****************************  BSECOND ******************************

        SUBROUTINE  BSECOND(BCKE,BCKP,N,ANOISE,IPCUBE,NN)

        REAL      ::  BCKE(N,N,N),BCKP(N,N,N)
        INTEGER   ::  IPCUBE(5,NN)

c$omp   parallel do private(kn,j,k,i)
        DO KN=1,NN
           J = IPCUBE(4,KN)
           K = IPCUBE(5,KN)
           DO I=IPCUBE(3,KN),IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)

               BCKE(I,J,K) = BCKE(I,J,K) + ANOISE*(18*BCKP(I,J,K)
     &          - 4.0*BCKP(I+1,J,K)+BCKP(I+2,J,K)
     &          - 4.0*BCKP(I-1,J,K)+BCKP(I-2,J,K)
     &          - 4.0*BCKP(I,J+1,K)+BCKP(I,J+2,K)
     &          - 4.0*BCKP(I,J-1,K)+BCKP(I,J-2,K)
     &          - 4.0*BCKP(I,J,K+1)+BCKP(I,J,K+2)
     &          - 4.0*BCKP(I,J,K-1)+BCKP(I,J,K-2)
     &          )
           ENDDO
        ENDDO
        END

C       ********************** BTHIRD *******************************

        SUBROUTINE  BTHIRD(BCKE,BCKP,N,ANOISE,IPCUBE,NN)

        REAL      :: BCKE(N,N,N),BCKP(N,N,N)
        INTEGER   :: IPCUBE(5,NN)

c$omp   parallel do private(kn,j,k,i)
        DO KN=1,NN
           J = IPCUBE(4,KN)
           K = IPCUBE(5,KN)
           DO I=IPCUBE(3,KN),IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)

                 BCKE(I,J,K) = BCKE(I,J,K) + ANOISE*(60*BCKP(I,J,K)
     &          - 15.0*BCKP(I+1,J,K)+6.0*BCKP(I+2,J,K)-BCKP(I+3,J,K)
     &          - 15.0*BCKP(I-1,J,K)+6.0*BCKP(I-2,J,K)-BCKP(I-3,J,K)
     &          - 15.0*BCKP(I,J+1,K)+6.0*BCKP(I,J+2,K)-BCKP(I,J+3,K)
     &          - 15.0*BCKP(I,J-1,K)+6.0*BCKP(I,J-2,K)-BCKP(I,J-3,K)
     &          - 15.0*BCKP(I,J,K+1)+6.0*BCKP(I,J,K+2)-BCKP(I,J,K+3)
     &          - 15.0*BCKP(I,J,K-1)+6.0*BCKP(I,J,K-2)-BCKP(I,J,K-3)
     &          )
           ENDDO
        ENDDO
        END


C       *************************** MAKETWOLISTS xxxxxxxxxxxxxxxxxxxxx

        SUBROUTINE MAKETWOLISTS(ILIST,N,ILIST1,N1,ILIST2,N2)

        IMPLICIT NONE

        INTEGER :: N,N1,N2
        INTEGER :: ILIST(N),ILIST1(N),ILIST2(N)

        LOGICAL :: RANDLIST(N)
        INTEGER :: IORD,K
        LOGICAL :: LTMP
        REAL    :: X

C       PURPOSE: CREATE 2 LIST OF IMAGES FOR EACH RECONSTRUCTION
C                FROM ONE LIST BY RANDOM SELECTION OF HALF OF THE
C                NUMBERS IN THE ORIGINAL LIST

	RANDLIST(1:N/2)   = .TRUE.
	RANDLIST(N/2+1:N) = .FALSE.

	DO  K=1,N
           CALL RANDOM_NUMBER(HARVEST=X)
           IORD           = MIN(N,MAX(1,INT(X*N+0.5)))
	   LTMP           = RANDLIST(IORD)
	   RANDLIST(IORD) = RANDLIST(K)
	   RANDLIST(K)    = LTMP
	ENDDO

        N1 = 0
        N2 = 0

        DO K=1,N
           IF (RANDLIST(K)) THEN
              N1         = N1 + 1
              ILIST1(N1) = ILIST(K)
           ELSE
              N2         = N2 + 1
              ILIST2(N2) = ILIST(K)
           ENDIF
        ENDDO

        END

C       *********************** BPCG_3_FBS ***************************

        SUBROUTINE BPCG_3_FBS(CB,BCKN,NXLD,N,IPCUBE,NN,DM,LB,NANG,RI,
     &                        NUMTH,BNORM,CHI2,
     &                        ERRM,CHIM,MAXIT,MODE,ANOISE, 
     &                        LDP,LDPNM,IRTFLG)

C       NUMTH = NUMTHREADS() FOR MP, OTHERWISE=1.
C          ( MP - MULTIPLE PROCESSING )

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'

        REAL               :: CB(N,N,N),BCKN(N,N,N)
        INTEGER            :: NXLD,N,NN
        INTEGER            :: IPCUBE(5,NN)
        REAL               :: DM(9,NANG)           
        INTEGER            :: LB(NANG)
        INTEGER            :: NANG
        REAL               :: RI
        INTEGER            :: NUMTH
        REAL               :: BNORM,CHI2,ERRM,CHIM
        INTEGER            :: MAXIT,MODE
        REAL               :: ANOISE
        INTEGER            :: LDP,LDPNM,IRTFLG

        INTEGER            :: KN, ICYCLE
        REAL               :: DELSQ
        REAL               :: PROJ(N*N,NUMTH)
        INTEGER            :: K,J,I,L_TH,IY,ILOC,IX,MYPID,NE,MWANT
        REAL               :: AKDEN,ERR,Q,P
        INTEGER            :: NNN,ITER,L_EN,L_NUM

        REAL, ALLOCATABLE  :: BCKE_SUM(:,:,:)
        REAL, ALLOCATABLE  :: BCKE    (:,:,:)
        REAL, ALLOCATABLE  :: BCKP    (:,:,:)

        REAL, ALLOCATABLE  :: XYZ  (:,:,:)
        REAL, ALLOCATABLE  :: X1   (:,:,:)
        REAL, ALLOCATABLE  :: Y1   (:,:,:)
        REAL, ALLOCATABLE  :: Z1   (:,:,:)
        REAL, ALLOCATABLE  :: XY2  (:,:,:)
        REAL, ALLOCATABLE  :: XZ2  (:,:,:)
        REAL, ALLOCATABLE  :: YZ2  (:,:,:)

        REAL, ALLOCATABLE  :: PROJPAD(:,:)   ! PADDED 2D PROJ ARRAY
        REAL, ALLOCATABLE  :: XDER   (:,:)   ! X  DERIVATIVE OF PROJ
        REAL, ALLOCATABLE  :: YDER   (:,:)   ! Y  DERIVATIVE OF PROJ
        REAL, ALLOCATABLE  :: XYDER  (:,:)   ! XY DERIVATIVE OF PROJ

        ICYCLE = 1

        MYPID  = -1

        NXLD   = N + 2 - MOD(N,2)

        ALLOCATE( PROJPAD(NXLD, N),
     &            XDER   (NXLD, N),
     &            YDER   (NXLD, N),
     &            XYDER  (NXLD, N),
     &            BCKE   (N,    N, N), 
     &            BCKP   (N,    N, N),
     &            XYZ    (NXLD, N, N),
     &            X1     (NXLD, N, N),
     &            Y1     (NXLD, N, N),
     &            Z1     (NXLD, N, N),
     &            XY2    (NXLD, N, N),
     &            XZ2    (NXLD, N, N),
     &            YZ2    (NXLD, N, N),
     &            STAT=IRTFLG)

        IF (IRTFLG .NE. 0) THEN
           MWANT = 4 * NXLD*N + 2 * N*N*N + 7*NXLD*N*N
           CALL ERRT(46,'BPCG_3_FBS; PROJPAD,...',MWANT)
           GOTO 9999
        ENDIF


        NNN = N*N*N
c$omp   parallel do private(k,j,i)
        DO K=1,N
           DO J=1,N
              DO I=1,N
                 BCKN(I,J,K)     = 0.0
                 BCKP(I,J,K)     = 0.0
#ifdef USE_MPI
                 BCKE_SUM(I,J,K) = 0.0
#endif
              ENDDO
           ENDDO
        ENDDO

        ERR  = 1.0
        ITER = 0
        Q    = 0.0

        DO  ITER=1,MAXIT

           DELSQ = 0.0
c$omp      parallel do private(kn,j,k,i),reduction(+:delsq)
           DO KN=1,NN
              J = IPCUBE(4,KN)
              K = IPCUBE(5,KN)
              DO I=IPCUBE(3,KN),IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)
                 DELSQ = DELSQ + CB(I,J,K)*CB(I,J,K)
              ENDDO
           ENDDO
           Q = Q * DELSQ

c$omp      parallel do private(kn,j,k,i)
           DO KN=1,NN
              J = IPCUBE(4,KN)
              K = IPCUBE(5,KN)
              DO I=IPCUBE(3,KN),IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)
                 BCKP(I,J,K) = CB(I,J,K) - Q * BCKP(I,J,K)
              ENDDO
           ENDDO
           Q = -1.0 / DELSQ

C          BCKP -> PROJ -> BCKE
C          LOOP OVER PROJECTIONS

c$omp      parallel do private(k,j,i)
           DO K=1,N
              DO J=1,N
                 DO I=1,N
                    BCKE(I,J,K) = 0.0
                    XYZ(I,J,K)  = BCKP(I,J,K)
                 ENDDO
                 DO I = N+1,NXLD
                    XYZ(I,J,K) = 0
                 ENDDO
               ENDDO
           ENDDO

C          CALCULATION OF PROJECTIONS DERIVATIVES USING 3D FFT
           CALL FBS3_PREP(XYZ, NXLD,N,N,N,
     &                    X1,Y1,Z1,XY2,XZ2, YZ2)

           !WRITE(6,*), 'Iteration cycle:', ICYCLE
           ICYCLE = ICYCLE + 1

           DO K=1,NANG,NUMTH
              L_EN  = MIN(NANG,K+NUMTH-1)
              L_NUM = MIN(NUMTH,NANG-K+1)

c$omp         parallel do private(i,j)
              DO J=1,L_NUM
                 DO I=1,N*N
                    PROJ(I,J) = 0.0
                 ENDDO
              ENDDO

c$omp         parallel do private(l_th),schedule(static)
              DO  L_TH=1,L_NUM

C                PROJECTS BCKP INTO PROJ  (FROM VOLUME TO PLANE)
                 CALL PRJCQ_FBS3(BCKP,DM(1,K+L_TH-1),
     &                PROJ(1,L_TH),N,NXLD,  XYZ,
     &                X1, Y1, Z1,
     &                XY2,XZ2,YZ2)
              ENDDO

C             HERE BCKPJ ITSELF IS MP  (MULTIPLE PROCESSING)

              DO  L_TH=1,L_NUM
C                MULTIPLY PROJECTIONS BY THEIR WEIGHTS
                 IF (LB(K+L_TH-1) .GT. 1)  THEN
c$omp               parallel do private(i)
                    DO  I=1,N*N
                       PROJ(I,L_TH) = PROJ(I,L_TH) * LB(K+L_TH-1)
                    ENDDO
                 ENDIF

C                PAD PROJECTION TO: PROJPAD
c$omp            parallel do private(iy,iloc,ix)
                 DO IY = 1,N
                    ILOC = (IY-1)*N + 1 
                    DO IX = 1,N 
                       PROJPAD(IX,IY) = PROJ(ILOC, L_TH)
                       ILOC           = ILOC + 1
                    ENDDO
                 ENDDO

C                CALCULATE PROJECTION DERIVATIVES USING FFT
C                IN FBS2_PREP PROJPAD DOES NOT RETURN USEFULL VALUE!
                 CALL FBS2_PREP(PROJPAD, XDER,YDER,XYDER, 
     &                          NXLD,N,N, IRTFLG)
   
C                BACKPROJECT FROM PLANE: PROJ TO VOLUME: BCKE
                 CALL BCKPJ_FBS(BCKE,NNN,DM(1,K+L_TH-1),
     &                  PROJ(1, L_TH), XDER,YDER,XYDER,
     &                  NXLD,N,IPCUBE,NN, LDP,LDPNM)
              ENDDO
           ENDDO

           IF (MODE .EQ. 1) THEN
              CALL FIXEDGE1(BCKP,NNN,BCKP,N,IPCUBE,NN)      
              CALL BFIRSTS(BCKE,BCKP,N,ANOISE,IPCUBE,NN)

           ELSEIF(MODE .EQ. 2) THEN
              CALL FIXEDGE2(BCKP,NNN,BCKP,N,IPCUBE,NN)      
              CALL BSECOND(BCKE,BCKP,N,ANOISE,IPCUBE,NN)

           ELSEIF (MODE .EQ. 3) THEN
              CALL FIXEDGE3(BCKP,NNN,BCKP,N,IPCUBE,NN)      
              CALL BTHIRD(BCKE,BCKP,N,ANOISE,IPCUBE,NN)
           ENDIF

           AKDEN = 0.0

c$omp      parallel do private(kn,j,k,i),reduction(+:akden)
           DO KN=1,NN
              J = IPCUBE(4,KN)
              K = IPCUBE(5,KN)
              DO I=IPCUBE(3,KN),
     &             IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)
                 AKDEN = AKDEN + BCKP(I,J,K) * BCKE(I,J,K)
              ENDDO
           ENDDO
           P = DELSQ / AKDEN

c$omp      parallel do private(kn,j,k,i)
           DO KN=1,NN
              J = IPCUBE(4,KN)
              K = IPCUBE(5,KN)
              DO I=IPCUBE(3,KN),
     &           IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)

                 BCKN(I,J,K) = BCKN(I,J,K)+P*BCKP(I,J,K)
                 CB(I,J,K)   = CB(I,J,K)-P*BCKE(I,J,K)
              ENDDO
           ENDDO

           ERR  = DELSQ / BNORM
           CHI2 = CHI2 - P * DELSQ

C         SKIP THIS IF NOT VERBOSE
          WRITE(NOUT,2041) ITER,ERR,CHI2
2041      FORMAT('  ITERATION: 'I3,
     &           '   DIFFERENCE:',1PE12.4,
     &           '   CHISQ:',1PE12.4)

C          CHECK STOPPING CRITERIA
           IF (ABS(ERR)  .LE. ERRM .OR. 
     &         ABS(CHI2) .LE. CHIM) THEN  
C             WRITE FINAL ERROR 
              IF (.NOT. VERBOSE)  WRITE(NOUT,2041) ITER,ERR,CHI2
              WRITE(NOUT,*) '  '  
              GOTO 9999            
           ENDIF
       ENDDO

9999   IF (ALLOCATED(BCKE))    DEALLOCATE (BCKE)
       IF (ALLOCATED(BCKP))    DEALLOCATE (BCKP)
       IF (ALLOCATED(XYZ))     DEALLOCATE (XYZ)
       IF (ALLOCATED(X1))      DEALLOCATE (X1)
       IF (ALLOCATED(Y1))      DEALLOCATE (Y1)
       IF (ALLOCATED(Z1))      DEALLOCATE (Z1)
       IF (ALLOCATED(XY2))     DEALLOCATE (XY2)
       IF (ALLOCATED(XZ2))     DEALLOCATE (XZ2)
       IF (ALLOCATED(YZ2))     DEALLOCATE (YZ2)
       IF (ALLOCATED(PROJPAD)) DEALLOCATE (PROJPAD)
       IF (ALLOCATED(XDER))    DEALLOCATE (XDER)
       IF (ALLOCATED(YDER))    DEALLOCATE (YDER)
       IF (ALLOCATED(XYDER))   DEALLOCATE (XYDER)

       END
