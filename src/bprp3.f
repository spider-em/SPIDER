C **********************************************************************
c  BPRP3
C         CORRECTIONS APPLIED ON THE VOL. SIDE      01/10/94
C         COMPRESSION OF ANGLES                     08/14/96
C         SYMMETRIES CORRECTED                      01/2001
C         REAL SPACE SYM. AFTER ITERATIVE PROCESS   05/02
C         OPFILEC                                   FEB 03 ARDEAN LEITH
C         ADDED REG_SET FOR ITER                    AUG 00 ARDEAN LEITH
C         BUILDM PARAMETERS                         JUL 03 ARDEAN LEITH
C         MPI                                       FEB 04 Chao Yang
C         REFACTORED                                OCT 08 ARDEAN LEITH
C         MPI REFACTORED                            OCT 08 ARDEAN LEITH
C         OUTPUT SHORTENED                          NOV 09 ARDEAN LEITH
C         ORDER OF PSI,THETA,PHI LISTED WRONG       OCT 10 ARDEAN LEITH
C         REPS RENAMED BPRP                         MAY 11 ARDEAN LEITH
C         ILIST ALLOCATED                           MAY 11 ARDEAN LEITH
C         PREPCUB                                   NOV 11 ARDEAN LEITH
C         BCKPJ_LIN                                 DEC 12 ARDEAN LEITH
C         BP RP 3                                   APR 12 ARDEAN LEITH
C         BP RP 3 NO OVERALL VOL                    AUG 12 ARDEAN LEITH
C         DKLP                                      JAN 13 ARDEAN LEITH *
C         NEWFILE IMGNUM == -1 BUG on gfort         APR 13 ARDEAN LEITH *
C
C=**********************************************************************
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright (C)2002,2013 P. A. Penczek & ArDean Leith                *
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
C  BPRP3
C
C  PURPOSE: REPROJECTIONS 3D, RICHARDSON'S METHOD, 
C           RECONSTRUCTION KEPT IN SQUARE TO INTRODUCE OTHER CONSTRAINTS.
C           AVERAGE OUTSIDE THE WINDOW IS SUBTRACTED
C           MIN, MAX RELATE TO THE PROJECTIONS.  SYMMETRIES IMPOSED.
C
C  CALLS: BPRP_2  (internal)
C         RPRQ 
C         ASTA 
C         PREPCUB 
C         BCKPJ_LIN 
C         PRJCQ_N 
C         BPRP_3  (internal)
C         SMT3_Q 
C         DOMIN3_S 
C         DOMAX3_S 
C         DOCORS3_S 
C         BMAX_C 
C         BMIN_C 
C         FMAX_Q 
C         FMIN_Q 
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE BPRP3()

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'F90ALLOC.INC'

        REAL, POINTER         :: PANG(:,:)
        REAL, POINTER         :: ANGSYM(:,:)
        REAL, ALLOCATABLE     :: SM(:,:),ANG(:,:), DM(:,:)
        REAL, ALLOCATABLE     :: CB(:,:,:)
        REAL, ALLOCATABLE     :: BCKN(:)
        REAL, ALLOCATABLE     :: BCKE(:)

        INTEGER, ALLOCATABLE  :: IPCUBE(:,:)
        INTEGER, ALLOCATABLE  :: LB(:) 

        INTEGER, ALLOCATABLE  :: ILIST(:),ILIST1(:),ILIST2(:)
 
        LOGICAL               :: WANT3,USELISTS,FBS_WANTED
        LOGICAL               :: MD,SAYIT,DO_OVERALL
        CHARACTER(LEN=MAXNAM) :: ANGDOC,FILPAT    ! MAXNAM FROM CMLIMIT
        CHARACTER(LEN=MAXNAM) :: SYMDOC           ! MAXNAM FROM CMLIMIT
        CHARACTER(LEN=MAXNAM) :: FILVOL 
        CHARACTER(LEN=1)      :: ANSW
        CHARACTER(LEN=1)      :: NULL = CHAR(0)
        DOUBLE PRECISION      :: ABA
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
 
C       MEMORY FOR LISTS OF IMAGES    
        NILMAX  = NIMAXPLUS      ! FROM CMLIMIT
        ALLOCATE(ILIST1(NILMAX),
     &           ILIST2(NILMAX),
     &           ILIST(NILMAX),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'BPRP3; ILIST....',3*NILMAX)
           RETURN
        ENDIF 

C       OPEN FIRST INPUT FILE
C       RETURNS: NANG = NUMBER OF ANGLES = NUMBER OF PROJECTIONS
        CALL OPFILES(0,LUNPROJ,LUNDOC,LUNXM,  
     &             ASKNAM,FILPAT,NLET, 'O',
     &             IFORM ,NX,NY,NZ,NSTACK,
     &             'TEMPLATE FOR IMAGE FILES~',
     &             FOUROK, ILIST,NILMAX, 
     &             NDUM,NANG,IMG1, IRTFLG) 
        IF (IRTFLG .NE. 0) RETURN
        
        IF (FILPAT(NLET:NLET) == '@') THEN
           CALL ERRT(101,
     &        'OPERATION DOES NOT WORK ON WHOLE STACKS',NE)
           GOTO 9999
        ENDIF  

        MAXNUM = MAXVAL(ILIST(1:NANG))

C       NANG - TOTAL NUMBER OF IMAGES
        IF (MYPID <= 0) WRITE(NOUT,90) NANG
90      FORMAT('  NUMBER OF IMAGES:',I8)
        
        WANT3    = (FCHAR(4:7) .EQ. 'RP 3')   ! WANT THREE VOLUMES
        USELISTS = (FCHAR(4:8) .EQ. 'RP 3L')  ! WANT THREE LISTS

#ifdef USE_MPI
        IF (WANT3) THEN   ! WANT THREE VOLUMES
           CALL ERRT(101,
     &     'OPERATION <BP RP 3> NOT AVAILABLE FOR MPI, USE <BP RP>',NE)
           GOTO 9999
        ENDIF
#endif

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

        RI = (NX / 2) - 2    ! DEFAULT VALUE
        CALL RDPRM1S(RI,N_UNUSED,
     &              'RADIUS OF RECONSTRUCTED OBJECT',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        IRI = RI

        NDIM  = NX   ! LINEAR DIMENSION OF PROJECTIONS & RESTORATION
        LDP   = NDIM / 2 + 1
        LDPNM = LDP 

C       IDUM IS A DUMMY VARIABLE, VALUE OF ,NDIM2 IS DETERMINED HERE
        MD = .FALSE.
        CALL PREPCUB_S(NDIM,NDIM2,IDUM,RI,MD,LDP)
 
C       USE NDIM2 TO ALLOCATE: IPCUBE 
C       TOTAL MEMORY IS VOLUMES: CB, BCKN, BCKE...
C       PLUS TWO 2D PROJECTIONS.
C       CB   - BACK-PROJECTED ORIGINAL PROJECTIONS, READ FROM FILE
C       BCKE - WORKING VOLUME
C       BCKN - CURRENT RECONSTRUCTION
        ALLOCATE (IPCUBE(5,NDIM2), 
     &            BCKN(NDIM*NDIM*NDIM), 
     &            BCKE(NDIM*NDIM*NDIM), 
     &            CB  (NDIM,NDIM,NDIM), 
     &            ANG (3,NANG),
     &            DM  (9,NANG),
     &            LB  (NANG), 
     &            STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           MWANT = 5*NDIM2 + NDIM*NDIM*NDIM + 3*NANG + 9*NANG +
     &             NANG + NDIM*NDIM*NDIM 
           CALL ERRT(46,'BPRP, IPCUBE...',MWANT)
           GOTO 9999
        ENDIF

C       MAKES LIST OF VOXEL LOCS ON EACH LINE VOLUME WITHIN RADIUS          *
        MD = .TRUE.
        CALL PREPCUB_S(NDIM,NDIM2,IPCUBE,RI,MD,LDP)  ! RETURNS: IPCUBE

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

C       RETRIEVE ARRAY WITH SYMMETRIES DATA IN IT
        MAXXS = 0
        NSYM  = 0
        CALL GETDOCDAT('SYMMETRIES DOC',.TRUE.,SYMDOC,
     &                 LUNANG,.TRUE.,MAXXS,
     &                 NSYM,ANGSYM,IRTFLG)
        IF (IRTFLG .NE. 0) NSYM = 1

        ALLOCATE(SM(9,NSYM), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'BP RP; SM',IER)
           GOTO 9999
        ENDIF
 
        IF (NSYM > 1)  THEN
           CALL BUILDS(SM,NSYM,ANGSYM(1,1),IRTFLG)
           DEALLOCATE(ANGSYM)

           WRITE(NOUT,2021) NSYM
2021       FORMAT(/,'  NUMBER OF SYMMETRIES:',I7,/)
        ENDIF

C       OPEN OUTPUT VOLUME
        MAXIM = 0
        IFORM = 3
        CALL OPFILEC(0,.TRUE.,FILVOL,LUNVOL,'U',
     &               IFORM,NDIM,NDIM,NDIM,
     &               MAXIM,'RECONSTRUCTED VOLUME',
     &               .FALSE.,IRTFLG)
        DO_OVERALL = (FILVOL(1:1) .NE. '*')
        IF (IRTFLG .NE. 0 .AND. DO_OVERALL ) GOTO 9999

        IF (WANT3) THEN
           MAXIM = 0
           IFORM = 3
           CALL OPFILEC(0,.TRUE.,FILVOL,LUNVOL1,'U',
     &                  IFORM,NDIM,NDIM,NDIM,
     &                  MAXIM,'FIRST SAMPLE VOLUME',
     &                 .FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FILVOL,LUNVOL2,'U',
     &               IFORM,NDIM,NDIM,NDIM,
     &               MAXIM,'SECOND SAMPLE VOLUME',
     &                .FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0)  GOTO 9999
        ENDIF

        ALA = 1.0e-4
        AIM = 0
        CALL  RDPRM2S(ALA,AIM, NOT_USED,   
     &               'LAMBDA, CORRECTION LIMIT',IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9999

        MAXIT = 100
        MODE  = 8
        CALL  RDPRI2S(MAXIT,MODE,NOT_USED,
     &                'ITERATION LIMIT, MODE',IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9999

        TMIN = FMIN   ! DEFAULT FROM COMMON
        TMAX = FMAX
        CALL  RDPRM2S(TMIN,TMAX, NOT_USED, 
     &                'MINIMUM, MAXIMUM',IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9999

        FBSFLAG = 0
        CALL RDPRM2S(SMOOTH, FBSFLAG,NOT_USED,
     &              'SMOOTHING CONST (0.0-0.999)',IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9999

        FBS_WANTED = (FBSFLAG > 0)  ! NOT DOCUMENTED!!!!!  

C       FIND NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)

C       CB   - KEEPS BACK-PROJECTED ORIGINAL PROJECTIONS, 
C              READ FROM THE DISK
C       BCKE - WORKING VOLUME
C       BCKN - CURRENT RECONSTRUCTION

        IF (DO_OVERALL) THEN
           WRITE(NOUT,92)'  CREATING VOLUME ---------------------------'
92         FORMAT(/,A)

C          RETRIEVE ANGLES, ANG ARRAY, AND DM FROM PANG
           SAYIT = (MYPID <= 0) .AND. VERBOSE
           CALL BPRP_DM(NANG,PANG,ILIST,ANG,DM,SAYIT,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          READS IMAGES AND BACKPROJECTS INTO: CB
           CALL BPRP_2(NDIM,NANG,CB,ILIST,IPCUBE,NDIM2,DM,
     &              RI,FILPAT,MAXIM, 
     &              ABA,SM,NSYM,LUNPROJ,LUNVOL,
     &              BNORM, LDP,LDPNM,FBS_WANTED,  IRTFLG)


C          COMPRESS ANGLES INTO DM - ALTERS: ANG,NANG1,DM & LO!!
           CALL HIANG(ANG,NANG,DM,LB,LO)
           NANG = LO

           IF (MYPID <= 0) WRITE(NOUT,2027) NANG
2027       FORMAT('  EFFECTIVE NUMBER OF ANGLES:',I7)

#if defined ( USE_MPI) && defined (MPI_DEBUG) 
           T0 = MPI_WTIME()
#endif

           CALL BPRP_3(BCKN,BCKE,NDIM,IPCUBE,NDIM2,
     &              DM,LB,NANG,IRI,ABA,
     &              SM,NSYM,NUMTH,LUNVOL,ITERDONE,BNORM,
     &              LDP,LDPNM,ALA,AIM,MAXIT,MODE,TMIN,TMAX,SMOOTH)

#if defined ( USE_MPI) && defined (MPI_DEBUG) 
           T1 = MPI_WTIME()  - T0
           IF (MYPID == 0) THEN
              WRITE(6,222) T1
 222          FORMAT('  BPRP TIME:',ES11.3)
           ENDIF 
#endif

C          SAVE OVERALL OUTPUT VOLUME
           CALL WRTVOL(LUNVOL,NDIM,NDIM, 1,NDIM,BCKN,IRTFLG)
           CLOSE(LUNVOL)
        ENDIF

C       SET ITERDONE IN REG NSEL(1)
        CALL REG_SET_NSEL(1,1,REAL(ITERDONE),0.0, 0.0, 
     &                    0.0, 0.0,IRTFLG)

        IF (WANT3) THEN

           WRITE(NOUT,*) ' '
           WRITE(NOUT,*) 
     &          ' CREATING FIRST SAMPLED VOLUME ------------------'

C          OPEN NEXT INPUT PROJECTION
           K      = 1
           IMGNUM = 0
           CALL NEXTFILE(K,       ILIST1, 
     &                   FOUROK,  0,
     &                   NANG1,   MAXIM,   
     &                   LUNPROJ, 0, 
     &                   FILPAT,  'O',
     &                   IMGNUM,  IRTFLG)
 
C          RETRIEVE ANGLES, ANG ARRAY, AND DM FROM PANG

           CALL BPRP_DM(NANG1,PANG,ILIST1,ANG,DM,SAYIT,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          READS IMAGES AND BACKPROJECTS THEM INTO: CB
           CALL BPRP_2(NDIM,NANG1,CB,ILIST1,IPCUBE,NDIM2,DM,
     &                 RI,FILPAT,MAXIM, 
     &                 ABA,SM,NSYM,LUNPROJ,LUNVOL1,
     &                 BNORM, LDP,LDPNM,FBS_WANTED,  IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          COMPRESS ANGLES INTO DM - ALTERS: ANG,NANG1,DM & LO!!
           CALL HIANG(ANG,NANG1,DM,LB,LO)
           NANG1 = LO
           IF (MYPID <= 0) WRITE(NOUT,2027) NANG1

           CALL BPRP_3(BCKN,BCKE,NDIM,IPCUBE,NDIM2,
     &              DM,LB,NANG1,IRI,ABA,
     &              SM,NSYM,NUMTH,LUNVOL1,ITERDONE,BNORM,
     &              LDP,LDPNM,ALA,AIM,MAXIT,MODE,TMIN,TMAX,SMOOTH)


C          SAVE FIRST OUTPUT SAMPLED VOLUME
           CALL WRTVOL(LUNVOL1,NDIM,NDIM, 1,NDIM,BCKN,IRTFLG)

           WRITE(NOUT,*) ' '
           WRITE(NOUT,*) 
     &         ' CREATING SECOND SAMPLED VOLUME ------------------'

C          OPEN NEXT INPUT PROJECTION
           K      = 1
           IMGNUM = 0
           CALL NEXTFILE(K,       ILIST2, 
     &                   FOUROK,  0,
     &                   NANG2,   MAXIM,   
     &                   LUNPROJ, 0, 
     &                   FILPAT,  'O',
     &                   IMGNUM,  IRTFLG)

C          RETRIEVE ANGLES, ANG ARRAY, AND DM FROM PANG
           CALL BPRP_DM(NANG2,PANG,  ILIST2,ANG,DM,SAYIT,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          READS IMAGES AND BACKPROJECTS THEM INTO: CB
           CALL BPRP_2(NDIM,NANG2,CB,ILIST2,IPCUBE,NDIM2,DM,
     &                 RI,FILPAT,MAXIM, 
     &                 ABA,SM,NSYM,LUNPROJ,LUNVOL2,
     &                 BNORM, LDP,LDPNM,FBS_WANTED,  IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          COMPRESS ANGLES INTO DM - ALTERS:  ANG,NANG2,DM & LO!!
           CALL HIANG(ANG,NANG2,DM,LB,LO)
           NANG2 = LO
           IF (MYPID <= 0) WRITE(NOUT,2027) NANG2

           CALL BPRP_3(BCKN,BCKE,NDIM,IPCUBE,NDIM2,
     &              DM,LB,NANG2,IRI,ABA,
     &              SM,NSYM,NUMTH,LUNVOL2,ITERDONE,BNORM,
     &              LDP,LDPNM,ALA,AIM,MAXIT,MODE,TMIN,TMAX,SMOOTH)

C          SAVE 2ND OUTPUT SAMPLED VOLUME
           CALL WRTVOL(LUNVOL2,NDIM,NDIM, 1,NDIM,BCKN,IRTFLG)
        ENDIF
        WRITE(NOUT,*) ' '


9999    IF (ALLOCATED(BCKN))   DEALLOCATE(BCKN)
        IF (ALLOCATED(BCKE))   DEALLOCATE(BCKE)
        IF (ALLOCATED(CB))     DEALLOCATE(CB)
        IF (ALLOCATED(LB))     DEALLOCATE(LB)
        IF (ALLOCATED(DM))     DEALLOCATE(DM)
        IF (ALLOCATED(ANG))    DEALLOCATE(ANG)
        IF (ALLOCATED(IPCUBE)) DEALLOCATE(IPCUBE)
        IF (ALLOCATED(ILIST))  DEALLOCATE(ILIST)
        IF (ALLOCATED(ILIST1)) DEALLOCATE(ILIST1)
        IF (ALLOCATED(ILIST2)) DEALLOCATE(ILIST2)
        IF (ALLOCATED(SM))     DEALLOCATE(SM)
        IF (ASSOCIATED(PANG))  DEALLOCATE(PANG)
                
        CLOSE(LUNVOL) 
        CLOSE(LUNVOL1) 
        CLOSE(LUNVOL2) 
        CLOSE(LUNPROJ) 
  
        END

C++****************************** BPRP_DM *****************************

        SUBROUTINE BPRP_DM(NANG,ANGBUF,ILIST,ANG,DM,SAYIT,IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'

        INTEGER               :: NANG
        REAL                  :: ANGBUF(4,NANG)
        INTEGER               :: ILIST(NANG)
        REAL                  :: ANG(3,NANG)
        REAL                  :: DM(9,NANG)
        INTEGER               :: IRTFLG

        INTEGER               :: K,ITMP,ICOUNT,J
        LOGICAL               :: SAYIT
        REAL                  :: SSDUM

        IRTFLG = 1

        DO K=1,NANG

C          LOAD ANGLES FOR THIS IMAGE
C          ORDER IN DOCUMENT FILE IS PSI, THETA, PHI AND ANGLES ARE IN 
C          DEGREES! IN ANG ARRAY IT IS OTHER WAY AROUND

           ITMP   = ILIST(K)
           ICOUNT = ANGBUF(1,ITMP)
           IF (ICOUNT <= 0) THEN
C             MISSING KEY
              CALL ERRT(102,'MISSING ANGLES FOR IMAGE',ITMP)
              RETURN
           ENDIF

           ANG(3,K) = ANGBUF(2,ITMP)  ! NOTE ORDER REVERSAL!
           ANG(2,K) = ANGBUF(3,ITMP)
           ANG(1,K) = ANGBUF(4,ITMP)

C          CREATE ROTATION MATRIX: DM
           CALL CANG(ANG(1,K),ANG(2,K),ANG(3,K),
     &              .FALSE.,SSDUM,DM(1,K))
 
           IF (SAYIT) THEN
              WRITE(NOUT,333) K,(ANG(J,K),J=3,1,-1)
333           FORMAT('  IMAGE #:',I7,
     &               '  PSI:',F6.1,' THETA:',F6.1,' PHI:',F6.1)
           ENDIF
         ENDDO
         IRTFLG = 0
         END

C++*********************************************************************
C                                                                      
C  BPRP_2.F   ADDED REG_SET FOR ITER             AUG 2000 ARDEAN LEITH
C             REFACTORED REDPRS-->BPRP_2         APR 2012 ARDEAN LEITH
C                                                                      
C++*********************************************************************
C
C  PURPOSE:  LOADS PROJECTIONS. FINDS BACKGROUND OUTSIDE
C            OF PROJECTED CIRCLE.  BACK-PROJECTS IMAGES INTO VOLUME.
C            DETERMINES AVERAGE OUTSIDE OF CIRCLE WHICH WILL LATER BE
C            SUBTRACTED FROM FINAL VOLUME.  BILINEAR INTERPOLATION      
C                                                                      
C++*********************************************************************

        SUBROUTINE BPRP_2(N,NANG,CB,ILIST,IPCUBE,NN,DM,
     &                    RI,FILPAT,MAXIM, 
     &                    ABA,SM,NSYM,LUNPROJ,LUNVOL,
     &                    BNORM, LDP,LDPNM,FBS_WANTED,  IRTFLG)
 
        USE TYPE_KINDS
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        INTEGER               :: N,NANG
        REAL                  :: CB(N,N,N),ANG(3,NANG)
        INTEGER               :: ILIST(NANG),IPCUBE(5,NN)
        INTEGER               :: NN
        REAL                  :: DM(3,3,NANG)
        REAL                  :: RI
        REAL                  :: ANGBUF(4,NANG)
        CHARACTER(LEN=*)      :: FILPAT
        INTEGER               :: MAXIM
        DOUBLE PRECISION      :: ABA
        REAL                  :: SM(3,3,NSYM)
        INTEGER               :: NSYM
        INTEGER               :: LUNPROJ,LUNVOL 
        REAL                  :: BNORM
        INTEGER               :: LDP,LDPNM
        LOGICAL               :: FBS_WANTED
        INTEGER               :: IRTFLG

        REAL                  :: DMS(3,3)
        CHARACTER(LEN=MAXNAM) :: FILNAM
        DOUBLE PRECISION      :: SUS,SSQ
        LOGICAL               :: FOUROK = .FALSE.
        REAL                  :: ADUM

C       AUTOMATIC
	REAL                  :: PROJ(N,N)

        INTEGER(KIND=I_8)     :: KLP_8
        DOUBLE PRECISION      :: DKLP,DKLP_LOC

#ifndef USE_MPI
C       AUTOMATIC ARRAY
        REAL                  :: PROJT(4,N*N)

        REAL, ALLOCATABLE     :: PROJPAD(:,:) ! PADDED 2D PROJ ARRAY
        REAL, ALLOCATABLE     :: XDER   (:,:) ! X  DERIVATIVE OF PROJ
        REAL, ALLOCATABLE     :: YDER   (:,:) ! Y  DERIVATIVE OF PROJ
        REAL, ALLOCATABLE     :: XYDER  (:,:) ! XY DERIVATIVE OF PROJ
#endif

#ifdef USE_MPI
        INCLUDE 'mpif.h'
        INTEGER               :: ISTAT(MPI_STATUS_SIZE)
        REAL   , ALLOCATABLE  :: CB_LOC(:,:,:)
        REAL   , ALLOCATABLE  :: PRJLOC(:,:,:), PRJBUF(:,:,:)
        INTEGER, ALLOCATABLE  :: PSIZE(:)
        INTEGER, ALLOCATABLE  :: NBASE(:)
        DOUBLE PRECISION      :: ABA_LOC

        ICOMM = MPI_COMM_WORLD
        CALL MPI_COMM_RANK(ICOMM, MYPID, MPIERR)
        CALL MPI_COMM_SIZE(ICOMM, NPROCS, MPIERR)

        ALLOCATE(PSIZE(NPROCS),
     &           NBASE(NPROCS),
     &           CB_LOC(N,N,N), 
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MWANT = 2*NPROCS + N*N*N
           CALL ERRT(46,'BPRP_2; PSIZE...',MWANT)
           RETURN
        ENDIF

        CB_LOC = 0.0

C       DATA DISTRIBUTION
        CALL SETPART(NANG, PSIZE, NBASE)
        NANG_LOC = PSIZE(MYPID+1) 

#ifdef MPI_DEBUG
        WRITE(6,111) NBASE(MYPID+1), MYPID
 111    FORMAT('  BPRP_2: NBASE: ', I5, ' MYPID:', I5)
        CALL FLUSHFILE(6)
#endif

        ABA      = 0.0D0
        DKLP     = 0
	CB       = 0.0      ! ARRAY INITIALIZATION

        ABA_LOC  = 0.0D0
        DKLP_LOC = 0

        ALLOCATE(PRJLOC(N,N,NANG_LOC), PRJBUF(N,N,PSIZE(1)), 
     &            STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'BPRP_2; PRJLOC, PRJBUF',IER)
           RETURN
        ENDIF

        DO IPROC = 1, NPROCS
           NLOC = PSIZE(IPROC)

C          READ A SUBSET OF IMAGES (ONLY ONE PROCESSOR READS)
 
           DO K=1,NLOC
              KGLB = K + NBASE(IPROC)
              NLET = 0
              CALL FILGET(FILPAT,FILNAM,NLET,ILIST(KGLB),IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 999
                                     
              MAXIM = 0
              CALL OPFILEC(0,.FALSE.,FILNAM,LUNPROJ,'O',IFORM,
     &                     LSAM,LROW,NSL,
     &                     MAXIM,'DUMMY',.FALSE.,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 999
                                                                                
              DO K2=1,N
                 CALL  REDLIN1P(LUNPROJ,PRJBUF(1,K2,K),N,K2)
              ENDDO
              IF (MYPID == 0) CLOSE(LUNPROJ)
           ENDDO

C          DISTRIBUTE IMAGES

           IF (IPROC .GT. 1) THEN
              IF (MYPID == 0) THEN
                 CALL SEND_MPI('BPRP_2','PRJBUF', PRJBUF, N*N*NLOC, 
     &                        'R',IPROC-1,IPROC-1, ICOMM)

              ELSE IF (MYPID .EQ. IPROC-1) THEN
                 CALL MPI_RECV(PRJLOC, N*N*NLOC   , MPI_REAL,
     &                         0     , MPI_ANY_TAG, ICOMM    ,
     &                         ISTAT , MPIERR)
                 IF (MPIERR .NE. 0) THEN
                     WRITE(6,*) ' BPRP_2: RECV FAILED'
                     STOP
                 ENDIF
              ENDIF
           ELSE IF (MYPID .EQ. 0) THEN
               CALL  SCOPY(N*N*NLOC,PRJBUF,1,PRJLOC,1)
           ENDIF 
        ENDDO

        DO K = 1, NANG_LOC
           KGLB = K + NBASE(MYPID+1)

C          ESTIMATE AVERAGE OUTSIDE THE CIRCLE
           CALL ASTA_D(PRJLOC(1,1,K),N,RI,ABA_LOC,DKLP_LOC)

           DO  ISYM=1,NSYM
              IF (NSYM .GT. 1)  THEN
C                SYMMETRIES, MULTIPLY MATRICES
                 DMS = MATMUL(SM(:,:,ISYM), DM(:,:,KGLB))
              ELSE
                 DMS = DM(:,:,KGLB)
              ENDIF

C             BACKPROJECTS PRJLOC INTO CB_LOC. BILINEAR INTERPOLATION   
              CALL RPRQ(N,PRJLOC(1,1,K),CB_LOC,IPCUBE,NN,DMS, 
     &                  LDP,LDPNM,IRTFLG)
              IF (IRTFLG .NE. 0) THEN
                 WRITE(6,*) ' BPRP_2: RPRQ FAILED'
                 STOP
              ENDIF
           ENDDO
        ENDDO
        N3 = N*N*N
        IF (ALLOCATED(PRJLOC)) DEALLOCATE(PRJLOC)

#ifdef MPI_DEBUG
        write(6,*) '  redprs: mpi_allreduce on cb..., mypid: ', mypid
#endif
        CALL ALLREDUCE_MPI('BPRP_2','CB', CB_LOC,CB,
     &                          N3, 'R','S',ICOMM)
        CALL ALLREDUCE_MPI('BPRP_2','ABA', ABA_LOC,ABA,
     &                           1, 'D','S',ICOMM)
        CALL ALLREDUCE_MPI('BPRP_2','DKLP', DKLP_LOC,DKLP,
     &                           1, 'D','S',ICOMM)

#ifdef MPI_DEBUG
        WRITE(6,*) '  BPRP_2: DONE MPI_ALLREDUCE..., MYPID: ', MYPID
#endif

C       --------------- END OF: MPI CODE   --------------------------
#else
C       --------------- START NON-MPI CODE --------------------------

        MYPID = -1
        ABA   = 0.0D0
        DKLP  = 0
	CB    = 0.0      ! ARRAY INITIALIZATION

        K = 1
        DO                ! LOOP OVER ALL INPUT PROJECTIONS

           ! LOAD NEXT PROJECTION (SQUARE) ALREADY OPENED
           CALL REDVOL(LUNPROJ,N,N,1,1,PROJ,IRTFLG)
  
C          ESTIMATE AVERAGE OUTSIDE THE CIRCLE
           CALL ASTA_D(PROJ,N,RI,ABA,DKLP)

           DO  ISYM=1,NSYM
              IF (NSYM > 1 )  THEN
C                SYMMETRIES, MULTIPLY MATRICES
                 DMS = MATMUL(SM(:,:,ISYM),DM(:,:,K))
              ELSE
                 DMS = DM(:,:,K)
              ENDIF

C             BACKPROJECTS PROJ INTO CB USING BILINEAR INTERPOLATION   
              CALL RPRQ(N,PROJ,CB,IPCUBE,NN,DMS, LDP,LDPNM,IRTFLG)
	   ENDDO

C          OPEN NEXT INPUT PROJECTION
           IMGNUM = 0
           CALL NEXTFILE(K,       ILIST, 
     &                   FOUROK,  0,
     &                   NANG,    MAXIM,   
     &                   LUNPROJ, 0, 
     &                   FILPAT,  'O',
     &                   IMGNUM,  IRTFLG)
 
           IF (IRTFLG .EQ. -1) EXIT      !  END OF INPUT STACK
           IF (IRTFLG .NE. 0) RETURN
        ENDDO

#endif
C       --------------- END OF:   NON-MPI CODE --------------------


C       CLOSE DOCUMENT FILE (LUNANG??) ??????????
        IF (MYPID <= 0) CLOSE(77)

        ABA   = ABA / DKLP
        KLP_8 = DKLP

C       PRINT STATISTICS
        IF (MYPID <= 0) WRITE(NOUT,2044) KLP_8,ABA
2044    FORMAT ('  TOTAL POINTS IN PROJECTIONS:',I10,
     &          '  AVERAGE OUTSIDE THE WINDOW: ',  ES10.3)

C       SUBTRACT AVERAGE FROM BACKPROJECTED VOL.
        BNORM = 0.0
        QT    = ABA * NANG * NSYM

        DO KN=1,NN
           J = IPCUBE(4,KN)
           K = IPCUBE(5,KN)
           DO I=IPCUBE(3,KN),IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)
              CB(I,J,K) = CB(I,J,K) - QT
              BNORM     = BNORM + CB(I,J,K) * CB(I,J,K)
           ENDDO

           CALL WRTLIN(LUNVOL,CB(1,J,K),N,(K-1)*N+J)
        ENDDO

        IRTFLG = 0

999     CONTINUE
#ifdef USE_MPI
        IF (ALLOCATED(CB_LOC)) DEALLOCATE(CB_LOC)
        IF (ALLOCATED(PSIZE))  DEALLOCATE(PSIZE)
        IF (ALLOCATED(NBASE))  DEALLOCATE(NBASE)
#endif
        END


C++*********************************************************************
C
C  BPRP_3.F     SPEEDED UP                     FEB 2000 ARDEAN LEITH
C               LDP PARAMETERS                 FEB 2009 ARDEAN LEITH
C               ALA... PARAMETERS              APR 2012 ARDEAN LEITH
C
C--*********************************************************************
C
C   BPRP_3(BCKN,N,IPCUBE,NN,DM,LB,NANG,IRI,ABA,YM, 
C          NSYM,NUMTH,INPIC, LDP,LDPNM,
C          ALA,AIM,MAXIT,MODE,TMIN,TMAX,SMOOTHT)
C
C  PARAMETERS:
C       NUMTH - NUMTHREDS FOR MP  
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE BPRP_3(BCKN,BCKE,N,IPCUBE,NN,DM,LB,NANG,IRI,ABA,
     &                    SM,NSYM,NUMTH,INPIC,ITERDONE,BNORM, 
     &               LDP,LDPNM,ALA,AIM,MAXIT,MODE,TMIN,TMAX,SMOOTHT)


        INCLUDE 'CMBLOCK.INC'

        REAL                  :: BCKN(N,N,N)
        REAL                  :: BCKE(N,N,N)
        REAL                  :: SM(3,3,NSYM),DM(3,3,NANG)
        INTEGER               :: IPCUBE(5,NN),LB(NANG)

C       CB,PROJ & PROJT ARE AUTOMATIC ARRAYS
        REAL ::  CB(N),PROJT(4,N*N),PROJ(N*N,NUMTH),DMS(3,3,NUMTH)

C       MASK IS AN AUTOMATIC ARRAY
        LOGICAL*1             :: MASK(N,N)
        DOUBLE PRECISION      :: ABA

        LOGICAL*1             :: ACTIVE_MIN, ACTIVE_MAX

#ifdef USE_MPI
        INCLUDE 'mpif.h'
        INTEGER               :: MPISTAT(MPI_STATUS_SIZE)
        REAL                  :: ALA, SQ, SQOLD, QT
        INTEGER, ALLOCATABLE  :: PSIZE(:)
        INTEGER, ALLOCATABLE  :: NBASE(:)
        INTEGER, ALLOCATABLE  :: LB_LOC(:)
        REAL,    ALLOCATABLE  :: DM_LOC(:,:,:)
        REAL,    ALLOCATABLE  :: BCKE_SUM(:,:,:)
        DOUBLE PRECISION      :: TSUM, TSUM0, TSUM1

        ICOMM = MPI_COMM_WORLD
        CALL MPI_COMM_RANK(ICOMM, MYPID , MPIERR)
        CALL MPI_COMM_SIZE(ICOMM, NPROCS, MPIERR)
#else
        MYPID = -1
#endif

        TMIN = TMIN - ABA
        TMAX = TMAX - ABA

        IF (MYPID .LE. 0) WRITE(NOUT,2059) TMIN,TMAX
2059    FORMAT('  MIN & MAX AFTER AVG SUBTRACTION:',
     &            ES10.3,' ... ',ES10.3)

C       CHANGE SMOOTH TO (ZERO,INFINITY) RANGE
        SMOOTH = SMOOTHT / (1.0 - SMOOTHT)

C       PREPARE THE LOGICAL MASK FOR MIN-MAX

        R  = IRI * IRI
        NC = N / 2 + 1
c$omp   parallel do private(j,i,qt,xx)
        DO J=1,N
           QT = J-NC
           XX = QT * QT
           DO I=1,N
              QT = I - NC
              IF (XX+QT*QT .LT. R) THEN
                 MASK(I,J) = .TRUE.
              ELSE
                 MASK(I,J) = .FALSE.
              ENDIF
           ENDDO
        ENDDO

        NMAT = N*N*N
        LTB  = N*N

c$omp   parallel do private(k,j,i)
C       ZEROS BCKE AND BCKN
        DO K=1,N
           DO J=1,N
              DO I=1,N
                 BCKE(I,J,K) = 0.0
                 BCKN(I,J,K) = 0.0
              ENDDO
           ENDDO
        ENDDO
        
        ACTIVE_MIN = .FALSE.
        ACTIVE_MAX = .FALSE.
        SQOLD      = 1.E23

C       ------------------------ START OF: MPI CODE -----------------
#ifdef USE_MPI
        ALLOCATE (BCKE_SUM(N,N,N),
     &            PSIZE(NPROCS),
     &            NBASE(NPROCS),
     &            STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           MWANT = 2*NPROCS + N*N*N
           CALL ERRT(46,'BPRP_3, BCKE_SUM...',IER)
           RETURN
        ENDIF
        BCKE_SUM = 0.0

        CALL SETPART(NANG, PSIZE, NBASE)
        NANG_LOC = PSIZE(MYPID+1)

        ALLOCATE(LB_LOC(NANG_LOC),DM_LOC(3, 3, NANG_LOC), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           MWANT = NANG_LOC+ 3*3*NANG_LOC
           CALL ERRT(46,'BPRP_3, LB_LOC',MWANT)
           RETURN
        ENDIF

        IDMTAG  = 1
        LBTAG   = 2
        MASTER  = 0
        IF (MYPID .EQ. 0) THEN

c          === MASTER DISTRIBUTES DM() and LB() ===

           DO IP = 2, NPROCS
              IBEGIN = NBASE(IP) + 1
              CALL MPI_SEND(DM(1,1,IBEGIN), 9*PSIZE(IP), MPI_REAL,
     &                      IP-1          , IDMTAG      , ICOMM    ,
     &                      MPIERR)
              call MPI_SEND(LB(IBEGIN)    , PSIZE(IP)  , MPI_INTEGER,
     &                      IP-1          , LBTAG      , ICOMM       ,
     &                      MPIERR)
           ENDDO
           CALL SCOPY(9*NANG_LOC, DM, 1, DM_LOC, 1)
           CALL SCOPY(NANG_LOC  , LB, 1, LB_LOC, 1)
        ELSE
c
c          === SLAVES RECEIVE FROM THE MASTER ===
C
           CALL MPI_RECV(DM_LOC, 9*NANG_LOC, MPI_REAL, MASTER, 
     &                   IDMTAG , ICOMM      , MPISTAT , MPIERR)
           CALL MPI_RECV(LB_LOC, NANG_LOC  , MPI_INTEGER, MASTER,
     &                   LBTAG , ICOMM      , MPISTAT    , MPIERR)
        ENDIF

#ifdef MPI_DEBUG
        WRITE(6,444) MYPID
 444    FORMAT('  BPRP_3: DATA DISTRIBUTION COMPLETED, MYPID = ',I3) 
        TSUM = 0.0
#endif

        DO  ITER=1,MAXIT
           IF (ITER .GT. 1)  THEN
              DO K=1,N
                 DO J=1,N
                    DO I=1,N
                       BCKE(I,J,K)     = 0.0
                       BCKE_SUM(i,j,k) = 0.0
                    ENDDO
                 ENDDO
              ENDDO
               
C             BCKN -> PROJ -> BCKE

              SBQ  = 0.0
              GMIN = 1.0E23
              GMAX = -1.0E23

C             LOOP OVER SYMMETRIES
	      DO ISYM=1,NSYM
C                LOOP OVER PROJECTIONS (NANG TIMES!)
                 DO K=1,NANG_LOC
                    DO I=1,LTB
                       PROJ(I,1)=0.0
                    ENDDO
c
                    IF(NSYM.GT.1)  THEN
C                      SYMMETRIES, MULTIPLY MATRICES
                       DMS(:,:,1)
     &                 = MATMUL(SM(:,:,ISYM),DM_LOC(:,:,K))
                    ELSE
                       DMS(:,:,1)= DM_LOC(:,:,K)
                    ENDIF
C                   CREATES PROJ FROM BCKN
                    CALL PRJCQ_N(BCKN,NMAT,DMS,PROJ,N,IPCUBE,
     &                           NN,LDP,LDPNM)
C
c                   *** have considered MODE = 0 ONLY ***
c  
                    IF ((MODE.EQ.2.OR.MODE.EQ.3.OR.MODE.EQ.7.OR
     &                   .MODE.EQ.8).AND. .NOT.ACTIVE_MIN)  
     &              CALL FMIN_Q(PROJ,MASK,LTB,GMIN)
     
                    IF ((MODE.EQ.5.OR.MODE.EQ.6.OR.MODE.EQ.7.OR.
     &                   MODE.EQ.8) .AND. .NOT.ACTIVE_MAX)  
     &              CALL FMAX_Q(PROJ,MASK,LTB,GMAX)

C                   HERE BCKPJ_LIN ITSELF IS MP
C                   MULTIPLY PROJECTIONS BY THEIR WEIGHTS
                    IF (LB_LOC(K) .GT. 1)  THEN
                       DO I=1,N*N
                          PROJ(I,1) = PROJ(I,1)*LB_LOC(K)
                       ENDDO
                    ENDIF

                    DO I = 1,N*N - N - 1
                       PT         = PROJ(I,1)
                       PROJT(1,i) = PT
                       PROJT(2,i) = PROJ(I+N,  1) - PT
                       PROJT(3,i) = PROJ(I+1,  1) - PT
                       PROJT(4,i) = PROJ(I+N+1,1) - PROJ(I+1,1)
     &                            - PROJT(2,I)
                    ENDDO

C                   BACKPROJECT FROM PROJT INTO BCKE
                    IF (NSYM .GT. 1)  THEN
C                      SYMMETRIES, MULTIPLY MATRICES
                       DMS(:,:,1) =
     &                     MATMUL(SM(:,:,ISYM),DM_LOC(:,:,K))
                    ELSE
                       DMS(:,:,1) = DM_LOC(:,:,K)
                    ENDIF
                    CALL BCKPJ_LIN(BCKE, NMAT  , DMS,PROJT,
     &                             N   , IPCUBE, NN, LDP,LDPNM)
C                END LOOP OVER PROJECTIONS
                 ENDDO

c                === GLOBAL SUM ===

#ifdef MPI_DEBUG
                 TSUM0 = MPI_WTIME() 
#endif
                 CALL MPI_ALLREDUCE(BCKE   , BCKE_SUM,
     &                              NMAT   , MPI_REAL,
     &                              MPI_SUM, ICOMM    ,
     &                              MPIERR)
                 IF (MPIERR .NE. 0) THEN
                    WRITE(0,*)'REDPRS: FAILED TO ALLREDUCE BCKE_SUM'
                    STOP
                 ENDIF
#ifdef MPI_DEBUG
                 TSUM1 = MPI_WTIME() 
                 TSUM  = TSUM + (TSUM1-TSUM0) 
#endif
C             END LOOP OVER SYMMETRIES
              ENDDO

C             END OF SECTION DONE FOR ITERATIONS > 1 --------------
           ENDIF

C          BEGIN ITERATIONS HERE

C          ALTERS BCKN IN SMT3_Q
c          *** NOT working in MPI yet
           IF (MODE.EQ.1.OR.MODE.EQ.3.OR.MODE.EQ.6.OR.MODE.EQ.8) 
     &          CALL SMT3_Q(T,ALA,BCKN,BCKN,N,N,N,IPCUBE,NN)
           SQ = 0.0

C          ONLY PROCESSOR READS CB IN THE FOLLOWING           

           DO KN=1,NN
              J = IPCUBE(4,KN)
              K = IPCUBE(5,KN)
              CALL  REDLIN1P(INPIC,CB,N,(K-1)*N+J)

              DO I=IPCUBE(3,KN), IPCUBE(3,KN)+IPCUBE(2,KN)
     &            -IPCUBE(1,KN)
                 QT          = CB(I) - BCKE_SUM(I,J,K)
                 SQ          = SQ + QT* QT
                 BCKN(I,J,K) = BCKN(I,J,K) + ALA * QT
              ENDDO
           ENDDO

           CALL MPI_BCAST(BCKN, NMAT, MPI_REAL, 0, ICOMM, IERR)
           CALL MPI_BCAST(QT, 1, MPI_REAL, 0, ICOMM, IERR)
           CALL MPI_BCAST(SQ, 1, MPI_REAL, 0, ICOMM, IERR)

           IF (MYPID .LE. 0)  WRITE(NOUT,2041) SQ, SQ/BNORM
2041       FORMAT('  SQUARED STRUCTURE CORRECTION:',ES12.4,2X,ES12.4)

           IF (MODE .NE. 0)  THEN
C             MODE > 0
              IF (ITER .GT. 1)  THEN
                 IF ((MODE.EQ.2.OR.MODE.EQ.3.OR.MODE.EQ.7.OR.MODE.EQ.8)
     &              .AND. .NOT.ACTIVE_MIN)  THEN
                    WRITE(NOUT,2061) GMIN
2061                FORMAT('  MINIMUM IN PROJECTIONS:',ES10.3) 
                    IF (GMIN .LT. TMIN)  THEN
                       CALL BMIN_C(BCKN,NMAT,IPCUBE,NN,BMIN)
                       WRITE(NOUT,2051)  BMIN
2051                  FORMAT('  MIN CONSTRAINT ACTIVATED, VALUE IN 3D:',
     &                        ES10.3)
                       ACTIVE_MIN = .TRUE.
                    ENDIF
                 ENDIF
C
                 IF ((MODE.EQ.5.OR.MODE.EQ.6.OR.MODE.EQ.7.OR.MODE.EQ.8) 
     &               .AND. .NOT.ACTIVE_MAX)  THEN
                    WRITE(NOUT,2062) GMAX
2062                FORMAT('  MAXIMUM IN PROJECTIONS:',ES10.3) 
                    IF (GMAX .GT. TMAX)  THEN
                       CALL BMAX_C(BCKN,NMAT,IPCUBE,NN,BMAX)
                       WRITE(NOUT,2052)  BMAX
2052                  FORMAT('  MAX CONSTRAINT ACTIVATED, VALUE IN 3D:',
     &                          ES10.3)
                       ACTIVE_MAX = .TRUE.
                    ENDIF
                 ENDIF
              ENDIF

C             ENFORCE MIN CONSTRAINTS
              IF (ACTIVE_MIN) CALL DOMIN3_S(BCKN,NMAT,IPCUBE,NN,BMIN)

C             ENFORCE MAX CONSTRAINTS
              IF (ACTIVE_MAX) CALL DOMAX3_S(BCKN,NMAT,IPCUBE,NN,BMAX)
C             END OF MODE>0
           ENDIF

           ITERDONE = ITER

C          CHECK STOPPING CRITERIA
           IF (SQ .GT. AIM .AND. ITER .LT. MAXIT)  THEN
              IF (SQ .LT. SQOLD) THEN
                 SQOLD = SQ
              ELSE
C               PERFORM ADDITIONAL SYMMETRIZATION IN REAL SPACE
                IF (NSYM .GT. 1)  THEN
                   DO K=1,N
                      DO J=1,N
                         DO I=1,N
                            BCKE(I,J,K) = BCKN(I,J,K)
                            BCKN(I,J,K) = 0.0
                         ENDDO
                      ENDDO
                   ENDDO

                   IF (MOD(N,2) .EQ. 0)  THEN
                      KNX = N/2-1
                   ELSE
                      KNX = N/2
                   ENDIF
                   KLX = -N/2
                   CALL SYMVOL(BCKE,BCKN,KLX,KNX,KLX,KNX,KLX,KNX,SM,
     &                         NSYM)
                ENDIF
                GOTO 999
              ENDIF
           ELSE
              IF (NSYM .GT. 1)  THEN
                 DO K=1,N
                    DO J=1,N
                       DO I=1,N
                          BCKE(I,J,K) = BCKN(I,J,K)
                          BCKN(I,J,K) = 0.0
                       ENDDO
                    ENDDO
                 ENDDO
                 IF (MOD(N,2) .EQ. 0)  THEN
                    KNX = N/2-1
                 ELSE
                    KNX = N/2
                 ENDIF
                 KLX = -N/2
                 CALL SYMVOL(BCKE,BCKN,KLX,KNX,KLX,KNX,KLX,KNX,SM,NSYM)
              ENDIF
              GOTO 999
           ENDIF
C       END ITERATION
        ENDDO

C       ------------------------ END OF: MPI CODE --------------------
#else

        DO  ITER=1,MAXIT

           IF (ITER .GT. 1)  THEN
c$omp         parallel do private(k,j,i)
              DO K=1,N
                 DO J=1,N
                    DO I=1,N
                       BCKE(I,J,K) = 0.0
                    ENDDO
                 ENDDO
              ENDDO
               
C             BCKN -> PROJ -> BCKE
              SBQ  = 0.0
              GMIN = 1.0E23
              GMAX = -1.0E23

C         LOOP OVER SYMMETRIES
	  DO  ISYM=1,NSYM
C             LOOP OVER PROJECTIONS (NANG TIMES!)
              DO K=1,NANG,NUMTH
                 L_EN  = MIN0(NANG,K+NUMTH-1)
                 L_NUM = MIN0(NUMTH,NANG-K+1)

c$omp            parallel do private(i,j)
                 DO J=1,L_NUM
                     DO I=1,LTB
                        PROJ(I,J) = 0.0
                     ENDDO
                 ENDDO

             DO  L_TH=1,L_NUM
               IF (NSYM .GT. 1)  THEN
C                 SYMMETRIES, MULTIPLY MATRICES
                  DMS(:,:,L_TH) = MATMUL(SM(:,:,ISYM),DM(:,:,K+L_TH-1))
               ELSE
                  DMS(:,:,L_TH) = DM(:,:,K+L_TH-1)
               ENDIF
	     ENDDO
c$omp            parallel do private(l_th),shared(nmat,n,nn)
                 DO L_TH=1,L_NUM

C                  CREATES PROJ FROM BCKN
                   CALL PRJCQ_N(BCKN,NMAT,DMS(1,1,L_TH),
     &                       PROJ(1,L_TH),N,IPCUBE,NN,LDP,LDPNM)
                 ENDDO
C  
                 DO L_TH=1,L_NUM
C                   LOOP OVER ALL PROCESSORS

                    IF ((MODE.EQ.2.OR.MODE.EQ.3.OR.MODE.EQ.7.OR
     &                  .MODE.EQ.8).AND. .NOT.ACTIVE_MIN)  
     &              CALL FMIN_Q(PROJ(1,L_TH),MASK,LTB,GMIN)
     
                    IF ((MODE.EQ.5.OR.MODE.EQ.6.OR.MODE.EQ.7.OR.
     &                   MODE.EQ.8) .AND. .NOT.ACTIVE_MAX)  
     &              CALL FMAX_Q(PROJ(1,L_TH),MASK,LTB,GMAX)
                 ENDDO

C                HERE BCKPJ_LIN ITSELF IS MP  
                 DO L_TH=1,L_NUM
C                   LOOP OVER ALL PROCESSORS
C                   MULTIPLY PROJECTIONS BY THEIR WEIGHTS
                    IF (LB(K+L_TH-1) .GT. 1)  THEN
c$omp                  parallel do private(i)
                       DO I=1,N*N
                           PROJ(I,L_TH) = PROJ(I,L_TH) * LB(K+L_TH-1)
                       ENDDO
                    ENDIF

c$omp               parallel do private(i,pt)
                    DO I = 1,N*N - N - 1
                       PT         = PROJ(I, L_TH)
                       PROJT(1,I) = PT
                       PROJT(2,I) = PROJ(I+N,  L_TH) - PT
                       PROJT(3,I) = PROJ(I+1,  L_TH) - PT 
                       PROJT(4,I) = PROJ(I+N+1,L_TH) - PROJ(I+1,L_TH) - 
     &                              PROJT(2,I)
                    ENDDO

C                   BACKPROJECT FROM PROJT INTO BCKE
                    IF (NSYM.GT.1)  THEN
C                      SYMMETRIES, MULTIPLY MATRICES
                       DMS(:,:,1)=MATMUL(SM(:,:,ISYM),DM(:,:,K+L_TH-1))
                    ELSE
                       DMS(:,:,1)=DM(:,:,K+L_TH-1)
                    ENDIF
                    CALL BCKPJ_LIN(BCKE,NMAT,DMS(1,1,1),PROJT,
     &                             N,IPCUBE,NN, LDP,LDPNM)
                 ENDDO
C             END LOOP OVER PROJECTIONS
              ENDDO
C            END LOOP OVER SYMMETRIES
             ENDDO

C             END OF SECTION DONE FOR ITERATIONS > 1 --------------
           ENDIF


C          BEGIN ITERATIONS HERE

C          ALTERS BCKN IN SMT3_Q
           IF (MODE.EQ.1.OR.MODE.EQ.3.OR.MODE.EQ.6.OR.MODE.EQ.8) 
     &          CALL SMT3_Q(SMOOTH,ALA,BCKN,BCKN,N,N,N,IPCUBE,NN)
           SQ = 0.0
           DO KN=1,NN
              J = IPCUBE(4,KN)
              K = IPCUBE(5,KN)
              CALL  REDLIN(INPIC,CB,N,(K-1)*N+J)
              DO I=IPCUBE(3,KN), IPCUBE(3,KN)+IPCUBE(2,KN)
     &                -IPCUBE(1,KN)
                 QT          = CB(I) - BCKE(I,J,K)
                 SQ          = SQ + QT* QT
                 BCKN(I,J,K) = BCKN(I,J,K) + ALA * QT
              ENDDO
           ENDDO

           WRITE(NOUT,2041) ITER, SQ,SQ/BNORM
2041       FORMAT('  ITERATION:',I6,
     &        '   SQUARED STRUCTURE CORRECTION:',ES12.4,2X,ES12.4)

          IF (MODE .NE. 0)  THEN
C             MODE > 0
              IF (ITER .GT. 1)  THEN
                 IF((MODE.EQ.2.OR.MODE.EQ.3.OR.MODE.EQ.7.OR.MODE.EQ.8)
     &              .AND. .NOT.ACTIVE_MIN)  THEN
                    WRITE(NOUT,2061) GMIN
2061                FORMAT('  MINIMUM IN PROJECTIONS =',1PE10.3) 
                    IF (GMIN .LT. TMIN)  THEN
                       CALL BMIN_C(BCKN,NMAT,IPCUBE,NN,BMIN)
                       WRITE(NOUT,2051)  BMIN
2051                  FORMAT('  MIN CONSTRAINT ACTIVATED, VALUE IN 3D:',
     &                        ES10.3)
                       ACTIVE_MIN=.TRUE.
                    ENDIF
                 ENDIF
C
                 IF ((MODE.EQ.5.OR.MODE.EQ.6.OR.MODE.EQ.7.OR.MODE.EQ.8) 
     &               .AND. .NOT.ACTIVE_MAX)  THEN
                    WRITE(NOUT,2062) GMAX
2062                FORMAT('  MAXIMUM IN PROJECTIONS:',ES10.3) 
                    IF (GMAX .GT. TMAX)  THEN
                       CALL BMAX_C(BCKN,NMAT,IPCUBE,NN,BMAX)
                       WRITE(NOUT,2052)  BMAX
2052                  FORMAT('  MAX CONSTRAINT ACTIVATED, VALUE IN 3D:',
     &                          ES10.3)
                       ACTIVE_MAX = .TRUE.
                    ENDIF
                 ENDIF
              ENDIF

C             ENFORCE MIN CONSTRAINTS
              IF (ACTIVE_MIN) CALL DOMIN3_S(BCKN,NMAT,IPCUBE,NN,BMIN)

C             ENFORCE MAX CONSTRAINTS
              IF (ACTIVE_MAX) CALL DOMAX3_S(BCKN,NMAT,IPCUBE,NN,BMAX)
C             END OF MODE>0
           ENDIF

           ITERDONE = ITER

C          CHECK STOPPING CRITERIA
           IF (SQ .GT. AIM .AND. ITER .LT. MAXIT)  THEN
              IF (SQ .LT. SQOLD) THEN
                 SQOLD = SQ
              ELSE

C               Perform additional symmetrization in real space
		IF (NSYM .GT. 1)  THEN
c$omp             parallel do private(k,j,i)
                  DO K=1,N
                     DO J=1,N
                       DO I=1,N
                         BCKE(I,J,K) = BCKN(I,J,K)
                         BCKN(I,J,K) = 0.0
                       ENDDO
                     ENDDO
                  ENDDO
                 IF (MOD(N,2) .EQ. 0)  THEN
                    KNX = N/2-1
                 ELSE
                    KNX = N/2
                 ENDIF
                 KLX = -N/2
	        CALL SYMVOL(BCKE,BCKN,KLX,KNX,KLX,KNX,KLX,KNX,SM,NSYM)
	      ENDIF 
              GOTO 999
           ENDIF
        ELSE

	   IF (NSYM .GT. 1)  THEN
c$omp         parallel do private(k,j,i)
               DO K=1,N
                  DO J=1,N
                    DO I=1,N
                       BCKE(I,J,K) = BCKN(I,J,K)
                       BCKN(I,J,K) = 0.0
                    ENDDO
                  ENDDO
              ENDDO

              IF (MOD(N,2) == 0)  THEN
                 KNX = N/2-1
              ELSE
                 KNX = N/2
              ENDIF
              KLX = -N/2

	      CALL SYMVOL(BCKE,BCKN,KLX,KNX,KLX,KNX,KLX,KNX,SM,NSYM)

           ENDIF 
           GOTO 999
        ENDIF

C       END ITERATION
        ENDDO
#endif
C       ------------------------ END OF: NON-MPI CODE ----------------

        IF (NSYM .GT. 1)  THEN
c$omp      parallel do private(k,j,i)
           DO K=1,N
              DO J=1,N
                 DO I=1,N
                    BCKE(I,J,K) = BCKN(I,J,K)
                    BCKN(I,J,K) = 0.0
                 ENDDO
              ENDDO
           ENDDO

           IF (MOD(N,2) == 0)  THEN
              KNX = N/2-1
           ELSE
              KNX = N/2
           ENDIF

           KLX = -N/2
           CALL SYMVOL(BCKE,BCKN,KLX,KNX,KLX,KNX,KLX,KNX,SM,NSYM)
        ENDIF 

999   CONTINUE

#ifdef USE_MPI
      IF (ALLOCATED(BCKE_SUM)) DEALLOCATE(BCKE_SUM)
      IF (ALLOCATED(DM_LOC))   DEALLOCATE(DM_LOC)
      IF (ALLOCATED(LB_LOC))   DEALLOCATE(LB_LOC)
      IF (ALLOCATED(PSIZE))    DEALLOCATE(PSIZE)
      IF (ALLOCATED(NBASE))    DEALLOCATE(NBASE)
#ifdef MPI_DEBUG
      IF (MYPID == 0) WRITE(6,888) TSUM
 888  FORMAT('  REDUCTION TIME: ', ES11.3) 
#endif
#endif

      END
