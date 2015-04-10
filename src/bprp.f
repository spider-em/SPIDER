C **********************************************************************
c  BPRP.F
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
C         DKLP                                      JAN 13 ARDEAN LEITH *
C         VERBOSE                                   MAY 14 ARDEAN LEITH *
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
C  BPRP
C
C  PURPOSE: REPROJECTIONS 3D, RICHARDSON'S METHOD, 
C           RECONSTRUCTION KEPT IN SQUARE TO INTRODUCE OTHER CONSTRAINTS.
C           AVERAGE OUTSIDE THE WINDOW IS SUBTRACTED
C           MIN, MAX RELATE TO THE PROJECTIONS
C           SYMMETRIES IMPOSED ...
C
C  CALLS: REDPRS
C         RPRQ 
C         ASTA 
C         PREPCUB 
C         BCKPJ_LIN 
C         PRJCQ_N 
C         REPR3Q
C         SMT3_Q 
C-        DOMIN3_S 
C-        DOMAX3_S 
C         DOCORS3_S 
C-        BMAX_C 
C-        BMIN_C 
C         FMAX_Q 
C         FMIN_Q 
C
C  NOTE: NEEDS ALLOC ERROR LEAK TRAPS
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE BPRP

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'F90ALLOC.INC'

        REAL,    POINTER        :: ANGSYM(:,:)
        REAL,    ALLOCATABLE    :: SM(:,:),ANG(:,:),DM(:,:)
        REAL,    ALLOCATABLE    :: BCKN(:,:,:)
        INTEGER, ALLOCATABLE    :: LB(:)
        INTEGER, ALLOCATABLE    :: IPCUBE(:,:)
        INTEGER, ALLOCATABLE    :: ILIST(:)
        CHARACTER(LEN=MAXNAM)   :: ANGDOC
        CHARACTER(LEN=MAXNAM)   :: FILNAM,FILPAT
        DOUBLE PRECISION        :: ABA
        LOGICAL                 :: MD
        DOUBLE PRECISION        :: T0, T1

        INTEGER, PARAMETER      :: LUNANG  = 77 
        INTEGER, PARAMETER      :: LUNPROJ = 97 
        INTEGER, PARAMETER      :: LUNVOL  = 98 
        INTEGER, PARAMETER      :: INPIC   = 99

        CALL SET_MPI(ICOMM,MYPID,MPIERR)

        NILMAX = NIMAXPLUS
        ALLOCATE(ILIST(NILMAX),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'BPRP; ILIST....',NILMAX)
           RETURN
        ENDIF 

        CALL FILELIST(.TRUE.,INPIC,FILPAT,NLET,ILIST,NILMAX,NANG,
     &                 'TEMPLATE FOR IMAGE FILES~',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       NANG - TOTAL NUMBER OF IMAGES
        IF (MYPID <= 0) WRITE(NOUT,2001) NANG
2001    FORMAT('  NUMBER OF IMAGES: ',I7)

        CALL  RDPRM(RI,NOT_USED,'RADIUS OF RECONSTRUCTED OBJECT')
        IRI = RI

        ALLOCATE (ANG(3,NANG), DM(9,NANG),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'BPRP, ANG,DM',12*NANG)
           GOTO 9999
        ENDIF

C       RETRIEVE ANGLES FROM DOC FILE (IN REVERSED ORDER)
        CALL GETDOCLIST('ANGLES DOC',LUNANG,ILIST,NANG,
     &                  1,3,.TRUE.,ANG,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       FIND DM
        DO K=1,NANG
           CALL CANG(ANG(1,K),ANG(2,K),ANG(3,K),
     &              .FALSE.,SSDUM,DM(1,K))

           IF (VERBOSE .AND. MYPID <= 0) THEN
              WRITE(NOUT,333) K,(ANG(J,K),J=3,1,-1)
333           FORMAT('  PROJECTION #',I7,
     &               '; PSI=',F6.1,' THETA=',F6.1,' PHI=',F6.1)
           ENDIF
         ENDDO
C       OPEN FIRST IMAGE FILE TO DETERMINE NSAM, NROW, NSL
        NLET = 0
        CALL FILGET(FILPAT,FILNAM,NLET,ILIST,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FILNAM,LUNPROJ,'O',IFORM,NSAM,NROW,NSL,
     &             MAXIM,'DUMMY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        CLOSE(LUNPROJ)

        NDIM    = NSAM
        NM      = 0
        LDP     = NDIM / 2 + 1
        LDPNM   = LDP + NM

C       RETRIEVE ARRAY WITH SYMMETRIES DATA IN IT
        MAXXS = 0
        NSYM  = 0
        CALL GETDOCDAT('SYMMETRIES DOC',.TRUE.,FILNAM,INPIC,
     &                  .TRUE.,MAXXS,NSYM,ANGSYM,IRTFLG)
        IF (IRTFLG .NE. 0)  NSYM = 1

        IF (NSYM > 1)  THEN
           ALLOCATE(SM(9,NSYM), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN 
              CALL ERRT(46,'BP RP, SM',IER)
              GOTO 9999
           ENDIF
           CALL BUILDS(SM,NSYM,ANGSYM(1,1),IRTFLG)
           DEALLOCATE(ANGSYM)

           WRITE(NOUT,2021) NSYM
2021       FORMAT(/,'  NUMBER OF SYMMETRIES: ',I7,/)

        ELSE
           ALLOCATE(SM(1,1), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN 
              CALL ERRT(46,'BP RP, SM-2nd',IER)
              GOTO 9999
           ENDIF
        ENDIF

C       DUM IS A DUMMY VARIABLE
        MD = .FALSE.
        CALL PREPCUB_S(NDIM,NDIMSQ,IDUM,RI,MD,LDP)

C       USE NDIMSQ TO ALLOCATE (IPCUBE) 
        ALLOCATE(IPCUBE(5,NDIMSQ), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'BPRP, IPCUBE',5*NDIMSQ)
           GOTO 9999
        ENDIF

        MD = .TRUE.
        CALL PREPCUB_S(NDIM,NDIMSQ,IPCUBE,RI,MD,LDP)

C       IN THIS VERSION TOTAL MEMORY IS THREE VOLUMES BCKN, BCKE AND 
C       CUBE PLUS 2 2D PROJECTIONS.
C
C       CUBE - KEEPS BACK-PROJECTED ORIGINAL PROJECTIONS, READ FROM THE DISK
C       BCKE - WORKING VOLUME
C       BCKN - CURRENT RECONSTRUCTION

C       FIND NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)

        IF (MYPID <= 0)  WRITE(NOUT,1001)
1001    FORMAT(/,'  REPROJECTION PROGRAM FOR 3-D BACK-PROJECTION',/)
     
        IFORM = 3
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNVOL,'U',IFORM,NDIM,NDIM,NDIM,
     &             MAXIM,'RECONSTRUCTED 3-D',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        CALL REDPRS(NDIM,NANG,ANG,ILIST,IPCUBE,NDIMSQ,DM,
     &         RI,ABA,SM,NSYM,LUNPROJ,LUNVOL,BNORM,FILPAT,
     &         LDP,LDPNM, IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        ALLOCATE(LB(NANG), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'BPRP, LB',NANG)
           GOTO 9999
        ENDIF

C       COMPRESS ANGLES - IT CHANGES NANG !!
        CALL HIANG(ANG,NANG,DM,LB,LO)
        NANG = LO

        IF (MYPID <= 0)  WRITE(NOUT,2027) NANG
2027    FORMAT('  EFFECTIVE NUMBER OF ANGLES: ',I7)

        ALLOCATE(BCKN(NDIM,NDIM,NDIM), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           MWANT = NDIM*NDIM*NDIM
           CALL ERRT(46,'REPS, BCKN',MWANT)
           GOTO 9999
        ENDIF

#if defined ( USE_MPI) && defined (MPI_DEBUG) 
        T0 = MPI_WTIME()
#endif

        CALL REPR3Q(BCKN,NDIM,IPCUBE,NDIMSQ,
     &               DM,LB,NANG,IRI,ABA,
     &               SM,NSYM,NUMTH,LUNVOL,ITERDONE,BNORM, LDP,LDPNM)

#if defined ( USE_MPI) && defined (MPI_DEBUG) 
        T1 = MPI_WTIME()
        T1 = T1 - T0
        IF (MYPID == 0) THEN
            WRITE(6,222) T1
 222        FORMAT(' BPRP TIME: ', 1PE11.3)
        ENDIF 
#endif

        CALL WRITEV(LUNVOL,BCKN,NDIM,NDIM,NDIM,NDIM,NDIM)
        CLOSE(LUNVOL)

C       SET ITERDONE IN REG NSEL(1)
        CALL REG_SET_NSEL(1,1,REAL(ITERDONE),0.0, 0.0, 0.0, 0.0,IRTFLG)

9999    IF (ALLOCATED(SM))     DEALLOCATE(SM)
        IF (ALLOCATED(BCKN))   DEALLOCATE(BCKN)
        IF (ALLOCATED(LB))     DEALLOCATE(LB)
        IF (ALLOCATED(DM))     DEALLOCATE(DM)
        IF (ALLOCATED(IPCUBE)) DEALLOCATE(IPCUBE)
        IF (ALLOCATED(ANG))    DEALLOCATE(ANG)
        IF (ALLOCATED(ILIST))  DEALLOCATE(ILIST)

        CALL FLUSHRESULTS()

        END


C++*********************************************************************
C *  REDPRS.F                 ADDED REG_SET FOR ITER AUG 00 ARDEAN LEITH
C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2001, P. A. Penczek                                   *
C=*                                                                    *
C=* University of Texas - Houston Medical School                       *
C=*                                                                    *
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

        SUBROUTINE REDPRS(N,NANG,ANG,ILIST,IPCUBE,NN,DM,
     &              RI,ABA,SM,NSYM,LUNPROJ,LUNVOL,
     &              BNORM,FILPAT, LDP,LDPNM, IRTFLG)
 
        USE TYPE_KINDS

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

	REAL, ALLOCATABLE         :: CB(:,:,:)
        REAL                      :: ANG(3,NANG)
        REAL                      :: SM(3,3,NSYM)
        INTEGER                   :: ILIST(NANG),IPCUBE(5,NN)
        REAL                      :: DM(3,3,NANG),DMS(3,3)

        CHARACTER(LEN=*)          :: FILPAT
        CHARACTER(LEN=MAXNAM)     :: FILNAM
        DOUBLE PRECISION          :: ABA,SUS,SSQ
        INTEGER(KIND=I_8)         :: KLP_8
        DOUBLE PRECISION          :: DKLP,DKLP_LOC

C       AUTOMATIC
	REAL                      :: PROJ(N,N)


#ifdef USE_MPI
        INCLUDE 'mpif.h'
        INTEGER                   :: ISTAT(MPI_STATUS_SIZE)
        REAL   , ALLOCATABLE      :: CB_LOC(:,:,:)
        REAL   , ALLOCATABLE      :: PRJLOC(:,:,:), PRJBUF(:,:,:)
        INTEGER, ALLOCATABLE      :: PSIZE(:)
        INTEGER, ALLOCATABLE      :: NBASE(:)
        DOUBLE PRECISION          :: ABA_LOC

        ICOMM = MPI_COMM_WORLD
        CALL MPI_COMM_RANK(ICOMM, MYPID, MPIERR)
        CALL MPI_COMM_SIZE(ICOMM, NPROCS, MPIERR)

        ALLOCATE(PSIZE(NPROCS),
     &           NBASE(NPROCS),
     &           CB_LOC(N,N,N), 
     &           STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           MWANT = 2*NPROCS + N*N*N
           CALL ERRT(46,'REDPRS, PSIZE...',MWANT)
           RETURN
        ENDIF

        CB_LOC = 0.0

C       DATA DISTRIBUTION

        CALL SETPART(NANG, PSIZE, NBASE)
        NANG_LOC = PSIZE(MYPID+1) 

#ifdef MPI_DEBUG
        WRITE(6,111) NBASE(MYPID+1), MYPID
 111    FORMAT('  REDPRS: NBASE: ', I5, ' MYPID: ', I5)
        CALL FLUSHFILE(6)
#endif
#else
        MYPID = -1
#endif

        ALLOCATE(CB(N,N,N), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'REDPRS, CB',N*N*N)
           RETURN
        ENDIF

        ABA   = 0.0D0
        DKLP  = 0
	CB    = 0.0

C       ------------------------ START OF: MPI CODE --------------------
#ifdef USE_MPI
        ABA_LOC  = 0.0D0
        DKLP_LOC = 0
        ALLOCATE(PRJLOC(N,N,NANG_LOC), PRJBUF(N,N,PSIZE(1)), 
     &            STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'REDPRS, PRJLOC, PRJBUF',IER)
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

           IF (IPROC > 1) THEN
              IF (MYPID == 0) THEN
                 CALL SEND_MPI('REDPRS','PRJBUF', PRJBUF, N*N*NLOC, 
     &                        'R',IPROC-1,IPROC-1, ICOMM)

              ELSE IF (MYPID == IPROC-1) THEN
                 CALL MPI_RECV(PRJLOC, N*N*NLOC   , MPI_REAL,
     &                         0     , MPI_ANY_TAG, ICOMM    ,
     &                         ISTAT , MPIERR)
                 IF (MPIERR .NE. 0) THEN
                     WRITE(6,*) 'REDPRS: RECV FAILED'
                     STOP
                 ENDIF
              ENDIF
           ELSE IF (MYPID == 0) THEN
               CALL  SCOPY(N*N*NLOC,PRJBUF,1,PRJLOC,1)
           ENDIF 
        ENDDO

        DO K = 1, NANG_LOC
           KGLB = K + NBASE(MYPID+1)

           CALL ASTA_D(PRJLOC(1,1,K),N,RI,ABA_LOC,DKLP_LOC)

           DO  ISYM=1,NSYM
              IF (NSYM > 1)  THEN
C                SYMMETRIES, MULTIPLY MATRICES
                 DMS = MATMUL(SM(:,:,ISYM),DM(:,:,KGLB))
              ELSE
                 DMS = DM(:,:,KGLB)
              ENDIF
C              CALL RPRQD(N,PRJLOC(1,1,K),CB_LOC,IPCUBE,NN,DMS,RI, 
C                         LDP,LDPNM)
              CALL RPRQ(N,PRJLOC(1,1,K),CB_LOC,IPCUBE,NN,DMS, 
     &                  LDP,LDPNM,IRTFLG)
              IF (IRTFLG .NE. 0) THEN
                 WRITE(6,*) 'REDPRS: RPRQ FAILED'
                 STOP
              ENDIF
           ENDDO
        ENDDO
        N3 = N*N*N
        IF (ALLOCATED(PRJLOC)) DEALLOCATE(PRJLOC)

#ifdef MPI_DEBUG
        write(6,*) 'redprs: mpi_allreduce on cb..., mypid = ', mypid
#endif
        CALL ALLREDUCE_MPI('REDPRS','CB', CB_LOC,CB,
     &                           N3, 'R','S',ICOMM)
        CALL ALLREDUCE_MPI('REDPRS','ABA', ABA_LOC,ABA,
     &                           1, 'D','S',ICOMM)
        CALL ALLREDUCE_MPI('REDPRS','DKLP', DKLP_LOC,DKLP,
     &                           1, 'D','S',ICOMM)

#ifdef MPI_DEBUG
        WRITE(6,*) 'REDPRS: DONE MPI_ALLREDUCE..., MYPID = ', MYPID
#endif

C       ------------------------ END OF: MPI CODE --------------------
#else

        DO K=1,NANG
           NLET = 0
           CALL FILGET(FILPAT,FILNAM,NLET,ILIST(K),IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 999

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FILNAM,LUNPROJ,'O',IFORM,
     &                  LSAM,LROW,NSL,
     &                  MAXIM,'DUMMY',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 999

           DO K2=1,N
              CALL  REDLIN(LUNPROJ,PROJ(1,K2),N,K2)
           ENDDO
           CLOSE(LUNPROJ)
  
           CALL ASTA_D(PROJ,N,RI,ABA,DKLP)

           DO  ISYM=1,NSYM
              IF (NSYM > 1 )  THEN
C                SYMMETRIES, MULTIPLY MATRICES
                 DMS = MATMUL(SM(:,:,ISYM),DM(:,:,K))
              ELSE
                 DMS = DM(:,:,K)
              ENDIF
              CALL RPRQ (N,PROJ,CB,IPCUBE,NN,DMS, LDP,LDPNM,IRTFLG)

	   ENDDO
        ENDDO

#endif
C       --------------- END OF:   NON-MPI CODE --------------------
C       CLOSE DOCUMENT FILE (LUNANG??)
        IF (MYPID <= 0) CLOSE(77)

        ABA   = ABA / DKLP
        KLP_8 = DKLP

C       PRINT STATISTICS
        IF (MYPID <= 0) WRITE(NOUT,2044)  KLP_8,ABA
2044    FORMAT
     &      ('  TOTAL POINTS IN PROJECTIONS:',I10,
     &     /,'  AVERAGE OUTSIDE THE WINDOW :',1PE10.3,/)

C       SUBTRACT THE AVERAGE

        BNORM = 0.0
        QT    = ABA*NANG*NSYM
        DO KN=1,NN
           J = IPCUBE(4,KN)
           K = IPCUBE(5,KN)
           DO I=IPCUBE(3,KN),IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)
              CB(I,J,K) = CB(I,J,K)-QT
              BNORM     = BNORM+CB(I,J,K)*CB(I,J,K)
           ENDDO
           CALL WRTLIN(LUNVOL,CB(1,J,K),N,(K-1)*N+J)
        ENDDO
        IRTFLG = 0

999     IF (ALLOCATED(CB)) DEALLOCATE(CB)
#ifdef USE_MPI
        IF (ALLOCATED(CB_LOC)) DEALLOCATE(CB_LOC)
        IF (ALLOCATED(PSIZE))  DEALLOCATE(PSIZE)
        IF (ALLOCATED(NBASE))  DEALLOCATE(NBASE)
#endif
        END


C++*********************************************************************
C *  REPR3Q.F  
C                                    SPEEDED UP FEB. 2000 ARDEAN LEITH
C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2001, P. A. Penczek                                   *
C=*                                                                    *
C=* University of Texas - Houston Medical School                       *
C=*                                                                    *
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
C   REPR3Q(BCKN,N,IPCUBE,NN,DM,LB,NANG,IRI,ABA,YM, 
C          NSYM,NUMTH,INPIC, LDP,LDPNM)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE REPR3Q(BCKN,N,IPCUBE,NN,DM,LB,NANG,IRI,ABA,
     &                    SM,NSYM,NUMTH,INPIC,ITERDONE,BNORM, 
     &                    LDP,LDPNM)

C       NUMTH - NUMTHREDS FOR MP; =NUMTHREADS() FOR ONYX, OTHERWISE=1.

        INCLUDE 'CMBLOCK.INC'

        REAL                  :: BCKN(N,N,N)
        REAL                  :: SM(3,3,NSYM),DM(3,3,NANG)
        INTEGER               :: IPCUBE(5,NN),LB(NANG)
        REAL, ALLOCATABLE     :: BCKE(:,:,:)

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
#ifdef MPI_DEBUG
        DOUBLE PRECISION      :: TSUM, TSUM0, TSUM1
#endif

        ICOMM = MPI_COMM_WORLD
        CALL MPI_COMM_RANK(ICOMM, MYPID , MPIERR)
        CALL MPI_COMM_SIZE(ICOMM, NPROCS, MPIERR)
#else
        MYPID = -1
#endif

        CALL RDPRM2(ALA,AIM,NOT_USED,   'LAMBDA, CORRECTION LIMIT')

        CALL RDPRMI(MAXIT,MODE,NOT_USED,'ITERATION LIMIT, MODE')

        CALL RDPRM2(TMIN,TMAX,NOT_USED, 'MINIMUM, MAXIMUM')

        TMIN = TMIN - ABA
        TMAX = TMAX - ABA

        IF (MYPID <= 0) WRITE(NOUT,2059) TMIN,TMAX
2059    FORMAT('  MINIMUM AND MAXIMUM AFTER AVERAGE SUBTRACTION',/,
     &           2(5X,1PE10.3))

        CALL RDPRM(T,NOT_USED,'SMOOTHING CONST (0.0-0.999)')

C       CHANGE IT TO (ZERO,INFINITY) RANGE
        T = T / (1.0-T)

C       PREPARE THE LOGICAL MASK FOR MIN-MAX

        R  = IRI*IRI
        NC = N/2+1
c$omp   parallel do private(j,i,qt,xx)
        DO J=1,N
           QT = J-NC
           XX = QT*QT
           DO I=1,N
              QT = I - NC
              IF (XX+QT*QT < R) THEN
                 MASK(I,J) = .TRUE.
              ELSE
                 MASK(I,J) = .FALSE.
              ENDIF
           ENDDO
        ENDDO

        NMAT = N*N*N
        LTB  = N*N

        ALLOCATE(BCKE(N,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'REPR3Q, BCKE',N*N*N)
           RETURN
        ENDIF

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
        IF (IRTFLG .NE. 0) THEN
           MWANT = 2*NPROCS + N*N*N
           CALL ERRT(46,'REPR3Q, BCKE_SUM...',MWANT)
           RETURN
        ENDIF
        BCKE_SUM = 0.0

        CALL SETPART(NANG, PSIZE, NBASE)
        NANG_LOC = PSIZE(MYPID+1)

        ALLOCATE(LB_LOC(NANG_LOC),DM_LOC(3, 3, NANG_LOC), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           MWANT = NANG_LOC+ 3*3*NANG_LOC
           CALL ERRT(46,'REPR3Q; LB_LOC...',MWANT)
           RETURN
        ENDIF

        IDMTAG  = 1
        LBTAG   = 2
        MASTER  = 0
        IF (MYPID == 0) THEN

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
 444    FORMAT('  REPR3Q: DATA DISTRIBUTION COMPLETED, MYPID: ',I3) 
        TSUM = 0.0
#endif

        IF (MYPID == 0)  WRITE(NOUT,971)
971     FORMAT(/'  COMMENCING ITERATION; 1 ')

        DO  ITER=1,MAXIT
           IF (ITER > 1)  THEN
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
                    IF(NSYM > 1)  THEN
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
                    IF ((MODE == 2 .OR. MODE == 3 .OR. MODE == 7 .OR.
     &                   MODE == 8) .AND. .NOT. ACTIVE_MIN)  
     &              CALL FMIN_Q(PROJ,MASK,LTB,GMIN)
     
                    IF ((MODE == 5.OR.MODE == 6 .OR. MODE == 7 .OR.
     &                   MODE == 8) .AND. .NOT. ACTIVE_MAX)  
     &              CALL FMAX_Q(PROJ,MASK,LTB,GMAX)

C                   HERE BCKPJ_LIN ITSELF IS MP
C                   MULTIPLY PROJECTIONS BY THEIR WEIGHTS
                    IF (LB_LOC(K) > 1)  THEN
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
                    IF (NSYM > 1)  THEN
C                      SYMMETRIES, MULTIPLY MATRICES
                       DMS(:,:,1)
     &                 = MATMUL(SM(:,:,ISYM),DM_LOC(:,:,K))
                    ELSE
                       DMS(:,:,1)=DM_LOC(:,:,K)
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
                    WRITE(0,*) 'REDPRS: FAILED TO ALLREDUCE BCKE_SUM'
                    STOP
                 ENDIF
#ifdef MPI_DEBUG
                 TSUM1 = MPI_WTIME() 
                 TSUM = TSUM + (TSUM1-TSUM0) 
#endif
C             END LOOP OVER SYMMETRIES
              ENDDO

C             END OF SECTION DONE FOR ITERATIONS > 1 --------------
           ENDIF

C          BEGIN ITERATIONS HERE

C          ALTERS BCKN IN SMT3_Q
c          *** NOT working in MPI yet
           IF (MODE == 1 .OR. MODE == 3 .OR. MODE == 6 .OR. MODE == 8) 
     &        CALL SMT3_Q(T,ALA,BCKN,BCKN,N,N,N,IPCUBE,NN)
           SQ = 0.0

C          ONLY PROCESSOR READS CB IN THE FOLLOWING           

           DO KN=1,NN
              J = IPCUBE(4,KN)
              K = IPCUBE(5,KN)
              CALL  REDLIN1P(INPIC,CB,N,(K-1)*N+J)
              DO I=IPCUBE(3,KN),IPCUBE(3,KN)+IPCUBE(2,KN)
     &            -IPCUBE(1,KN)
                 QT          = CB(I) - BCKE_SUM(I,J,K)
                 SQ          = SQ + QT* QT
                 BCKN(I,J,K) = BCKN(I,J,K) + ALA * QT
              ENDDO
           ENDDO
           CALL MPI_BCAST(BCKN, NMAT, MPI_REAL, 0, ICOMM, IERR)
           CALL MPI_BCAST(QT, 1, MPI_REAL, 0, ICOMM, IERR)
           CALL MPI_BCAST(SQ, 1, MPI_REAL, 0, ICOMM, IERR)

           IF ((MYPID <= 0) .AND.
     &        (VERBOSE .AND. (ITER == 1 .OR. ITER > (MAXIT  -9))))  
     &           WRITE(NOUT,2041) SQ, SQ/BNORM
2041       FORMAT('  SQUARED CORRECTION OF THE STRUCTURE',2(2X,1PE12.4))

           IF (MODE .NE. 0)  THEN
C             MODE > 0
              IF (ITER > 1)  THEN
                 IF ((MODE == 2 .OR. MODE == 3 .OR.
     &                MODE == 7 .OR. MODE == 8)
     &               .AND. .NOT. ACTIVE_MIN)  THEN
                    WRITE(NOUT,2061) GMIN
2061                FORMAT('  MINIMUM IN PROJECTIONS =',1PE10.3) 
                    IF (GMIN < TMIN)  THEN
                       CALL BMIN_C(BCKN,NMAT,IPCUBE,NN,BMIN)
                       WRITE(NOUT,2051)  BMIN
2051                  FORMAT('  MIN CONSTRAINT ACTIVATED, VALUE IN 3D:',
     &                        1PE10.3)
                       ACTIVE_MIN = .TRUE.
                    ENDIF
                 ENDIF
C
                 IF ((MODE == 5 .OR. MODE == 6  .OR. 
     &                MODE == 7 .OR. MODE == 8) .AND. 
     &               .NOT. ACTIVE_MAX)THEN
                    WRITE(NOUT,2062) GMAX
2062                FORMAT('  MAXIMUM IN PROJECTIONS=',1PE10.3) 
                    IF (GMAX > TMAX)  THEN
                       CALL BMAX_C(BCKN,NMAT,IPCUBE,NN,BMAX)
                       WRITE(NOUT,2052)  BMAX
2052                  FORMAT('  MAX CONSTRAINT ACTIVATED, VALUE IN 3D=',
     &                          1PE10.3)
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
           IF (SQ > AIM .AND. ITER < MAXIT)  THEN
              IF (SQ < SQOLD) THEN
                 SQOLD = SQ
              ELSE
C               PERFORM ADDITIONAL SYMMETRIZATION IN REAL SPACE
                IF (NSYM > 1)  THEN
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
                   CALL SYMVOL(BCKE,BCKN,KLX,KNX,KLX,KNX,KLX,KNX,SM,
     &                         NSYM)
                ENDIF
                GOTO 999
              ENDIF
           ELSE
              IF (NSYM > 1)  THEN
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

C       ------------------------ END OF: MPI CODE --------------------
#else
        !WRITE(NOUT,975)
975     FORMAT('  COMMENCING ITERATION: 1 ')

        DO  ITER=1,MAXIT

           IF (ITER > 1)  THEN
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
               IF (NSYM > 1)  THEN
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
                    IF ((MODE == 2 .OR. MODE == 3 .OR. MODE == 7 .OR.
     &                   MODE == 8 ) .AND. .NOT. ACTIVE_MIN)  
     &              CALL FMIN_Q(PROJ(1,L_TH),MASK,LTB,GMIN)
     
                    IF ((MODE == 5 .OR. MODE == 6 .OR. MODE == 7 .OR.
     &                   MODE == 8) .AND. .NOT. ACTIVE_MAX)  
     &              CALL FMAX_Q(PROJ(1,L_TH),MASK,LTB,GMAX)
                 ENDDO

C                HERE BCKPJ_LIN ITSELF IS MP  
                 DO L_TH=1,L_NUM
C                   LOOP OVER ALL PROCESSORS
C                   MULTIPLY PROJECTIONS BY THEIR WEIGHTS
                    IF (LB(K+L_TH-1) > 1)  THEN
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
                    IF (NSYM > 1)  THEN
C                      SYMMETRIES, MULTIPLY MATRICES
                       DMS(:,:,1)=MATMUL(SM(:,:,ISYM),DM(:,:,K+L_TH-1))
                    ELSE
                       DMS(:,:,1) = DM(:,:,K+L_TH-1)
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
           IF (MODE == 1 .OR.  MODE == 3 .OR. MODE == 6 .OR. MODE == 8) 
     &          CALL SMT3_Q(T,ALA,BCKN,BCKN,N,N,N,IPCUBE,NN)
           SQ = 0.0
           DO KN=1,NN
              J = IPCUBE(4,KN)
              K = IPCUBE(5,KN)
              CALL  REDLIN(INPIC,CB,N,(K-1)*N+J)
              DO I=IPCUBE(3,KN),IPCUBE(3,KN)+IPCUBE(2,KN)
     &                -IPCUBE(1,KN)
                 QT          = CB(I) - BCKE(I,J,K)
                 SQ          = SQ + QT* QT
                 BCKN(I,J,K) = BCKN(I,J,K) + ALA * QT
              ENDDO
           ENDDO

           IF (VERBOSE .OR. (ITER == 1 .OR. ITER  > (MAXIT  -9)))  
     &        WRITE(NOUT,2041) ITER, SQ,SQ/BNORM
2041       FORMAT('  ITERATION: ',I5,
     &        '   SQUARED CORRECTION OF STRUCTURE:',2(2X,1PE12.4))

          IF (MODE .NE. 0)  THEN
C             MODE > 0
              IF (ITER > 1)  THEN
                 IF((MODE == 2 .OR. MODE == 3 .OR. 
     &               MODE == 7 .OR. MODE == 8) .AND. 
     &              .NOT. ACTIVE_MIN)THEN
                    WRITE(NOUT,2061) GMIN
2061                FORMAT('  MINIMUM IN PROJECTIONS =',1PE10.3) 
                    IF (GMIN < TMIN)  THEN
                       CALL BMIN_C(BCKN,NMAT,IPCUBE,NN,BMIN)
                       WRITE(NOUT,2051)  BMIN
2051                  FORMAT('  MIN CONSTRAINT ACTIVATED, VALUE IN 3D=',
     &                        1PE10.3)
                       ACTIVE_MIN = .TRUE.
                    ENDIF
                 ENDIF
C
                 IF ((MODE == 5 .OR. MODE == 6 .OR.
     &                MODE == 7 .OR. MODE == 8) .AND. 
     &               .NOT. ACTIVE_MAX)  THEN
                    WRITE(NOUT,2062) GMAX
2062                FORMAT('  MAXIMUM IN PROJECTIONS =',1PE10.3) 
                    IF (GMAX > TMAX)  THEN
                       CALL BMAX_C(BCKN,NMAT,IPCUBE,NN,BMAX)
                       WRITE(NOUT,2052)  BMAX
2052                  FORMAT('  MAX CONSTRAINT ACTIVATED, VALUE IN 3D=',
     &                          1PE10.3)
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
           IF (SQ > AIM .AND. ITER < MAXIT)  THEN
              IF (SQ < SQOLD) THEN
                 SQOLD = SQ
              ELSE

C               Perform additional symmetrization in real space
		IF (NSYM > 1)  THEN
c$omp             parallel do private(k,j,i)
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
        ELSE

	   IF (NSYM > 1)  THEN
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

        IF (NSYM > 1)  THEN
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

999   IF (ALLOCATED(BCKE))  DEALLOCATE (BCKE)
#ifdef USE_MPI
      if (ALLOCATED(BCKE_SUM)) DEALLOCATE(BCKE_SUM)
      if (ALLOCATED(DM_LOC))   DEALLOCATE(DM_LOC)
      if (ALLOCATED(LB_LOC))   DEALLOCATE(LB_LOC)
      if (ALLOCATED(PSIZE))    DEALLOCATE(PSIZE)
      if (ALLOCATED(NBASE))    DEALLOCATE(NBASE)
#ifdef MPI_DEBUG
      IF (MYPID == 0) WRITE(6,888) TSUM
 888  FORMAT('  REDUCTION TIME: ', 1PE11.3) 
#endif
#endif
      END




