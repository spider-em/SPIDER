C **********************************************************************
C *  NN4                  COSMETIC & ERROR TRAPS  JUNE 08 ARDEAN LEITH
C
C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2003,2008 P. A. Penczek                                   *
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
C PURPOSE??: BACK PROJECTION USING NEAREST NEIGHBOR INTERPOLATION IN
C            FOURIER SPACE???? AL
C
C***********************************************************************

        SUBROUTINE NN4

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'F90ALLOC.INC'

        REAL, DIMENSION(:,:), ALLOCATABLE      :: DM,SM
        COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: X
        REAL, DIMENSION(:), ALLOCATABLE        :: TEMP
C       DOC FILE POINTERS
        REAL, DIMENSION(:,:), POINTER          :: ANGBUF, ANGSYM

        CHARACTER(LEN=MAXNAM)                  :: FINPAT,FINPIC,FILNAM
        CHARACTER(LEN=MAXNAM)                  :: ANGDOC

        DATA  IOPIC/68/,INPIC/69/

        NILMAX = NIMAX

        CALL FILELIST(.TRUE.,INPIC,FINPAT,NLET,INUMBR,NILMAX,NANG,
     &                 'TEMPLATE FOR 2-D IMAGES',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CLOSE(INPIC)
        MAXNUM = MAXVAL(INUMBR(1:NANG))

C       N    - LINEAR DIMENSION OF PROJECTIONS AND RESTORED CUBE
C       NANG - TOTAL NUMBER OF IMAGES
        WRITE(NOUT,2001) NANG
2001    FORMAT(' NUMBER OF IMAGES =',I5)

C       RETRIEVE ARRAY WITH ANGLES DATA IN IT
        MAXXT = 4
        MAXYT = MAXNUM
        CALL GETDOCDAT('ANGLES DOC',.TRUE.,ANGDOC,77,.FALSE.,MAXXT,
     &                       MAXYT,ANGBUF,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

C       RETRIEVE ARRAY WITH SYMMETRIES DATA IN IT
        MAXXS  = 0
        MAXSYM = 0
        CALL GETDOCDAT('SYMMETRIES DOC',.TRUE.,ANGDOC,77,.TRUE.,MAXXS,
     &                   MAXSYM,ANGSYM,IRTFLG)
        IF (IRTFLG .NE. 0) MAXSYM = 1

C       OPEN FIRST IMAGE FILE TO DETERMINE NSAM, NROW, NSL
        CALL FILGET(FINPAT,FINPIC,NLET,INUMBR(1),INTFLG)

        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,NSAM,NROW,NSL,
     &             MAXIM,'DUMMY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 999
        CLOSE(INPIC)

        N2     = 4 * NSAM   ! removed jun 2008 al (seems wrong)
        N2     = 2 * NSAM
        LSD    = N2+2-MOD(N2,2)
        NMAT   = LSD*N2*N2

        ALLOCATE(DM(9,NANG), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'BP NF, DM',9*NANG)
           GOTO 999
        ENDIF

        CALL BUILDM(INUMBR,DM,NANG,ANGBUF(1,1),.FALSE.,SSDUM,
     &              .FALSE.,IRTFLG)
        DEALLOCATE(ANGBUF)
        IF (IRTFLG .NE. 0) GOTO 999

        IF (MAXSYM .GT. 1)  THEN
           ALLOCATE(SM(9,MAXSYM), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN 
              MWANT = 9 * MAXSYM
              CALL ERRT(46,'BP NF, SM',MWANT)
              GOTO 999
           ENDIF
           CALL  BUILDS(SM,MAXSYM,ANGSYM(1,1),IRTFLG)
           DEALLOCATE(ANGSYM)
        ELSE
           ALLOCATE(SM(1,1), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN 
              CALL ERRT(46,'BP NF, SM-2nd',1)
              GOTO 999
           ENDIF
        ENDIF

        ALLOCATE (X(0:N2/2,N2,N2), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           MWANT = (N2/2 + 1)*N2*N2
           write(nout,*) 'mwant: ',mwant

           CALL ERRT(46,'BP NF, X',MWANT)
           GOTO 999
        ENDIF

        CALL NN4Q(NSAM,X, LSD,N2,N2/2,INUMBR,DM,NANG,SM,
     &            MAXSYM,FINPAT,FINPIC,NLET)

C       ADDITIONAL SYMMETRIZATION OF THE VOLUME IN REAL SPACE 05/03/02
        IF (MAXSYM .GT. 1)  THEN
           ALLOCATE(TEMP(NSAM*NSAM*NSAM), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN 
              MWANT = NSAM*NSAM*NSAM
              CALL ERRT(46,'BP NF, TEMP',MWANT)
              GOTO 999
           ENDIF

           CALL COP(X,TEMP,NSAM*NSAM*NSAM)

c$omp      parallel do private(i,j,k)
           DO K=1,N2
              DO J=1,N2
                 DO I=0,N2/2
                    X(I,J,K) = CMPLX(0.0,0.0)
                 ENDDO
              ENDDO
           ENDDO

           IF (MOD(NSAM,2) .EQ. 0) THEN
              KNX = NSAM/2-1
           ELSE
              KNX = NSAM/2
           ENDIF
           KLX = -NSAM/2

           CALL SYMVOL(TEMP,X,KLX,KNX,KLX,KNX,KLX,KNX,SM,MAXSYM)
           DEALLOCATE(TEMP)
        ENDIF

        IFORM = 3
        CALL OPFILEC(0,.TRUE.,FILNAM,IOPIC,'U',IFORM,NSAM,NSAM,NSAM,
     &           MAXIM,'RECONSTRUCTED 3-D',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9997

C       NOTE: NSAM=NROW=NSLICE 
        CALL WRITEV(IOPIC,X,NSAM,NSAM,NSAM,NSAM,NSAM)

9997    CLOSE(IOPIC)

999     IF (ALLOCATED(DM)) DEALLOCATE(DM)
        IF (ALLOCATED(SM)) DEALLOCATE(SM)
        IF (ALLOCATED(X))  DEALLOCATE(X)

        END


C       ------------------ NN4Q ----------------------------------

        SUBROUTINE NN4Q(NS,X,LSD,N,N2,ILIST,DM,NANG,SM,MAXSYM,
     &                   FINPAT,FINPIC,NLET)
    
        INCLUDE 'CMLIMIT.INC'

        COMPLEX                                :: X(0:N2,N,N)
        DIMENSION                              :: ILIST(NANG)
        DIMENSION                              :: DM(3,3,NANG)
        DIMENSION                              :: SM(3,3,MAXSYM)
        DIMENSION                              :: DMS(3,3)
        CHARACTER(LEN=MAXNAM)                  :: FINPAT,FINPIC

        INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: NR
        REAL, ALLOCATABLE, DIMENSION(:,:)      :: PROJ
        COMPLEX, ALLOCATABLE, DIMENSION(:,:)   :: BI
        DOUBLE PRECISION                       :: PI

        DATA  INPROJ/69/

        ALLOCATE (NR(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           MWANT = (N2 + 1) * N * N
           CALL ERRT(46,'BP NF, NR',MWANT)
           RETURN
        ENDIF

c$omp   parallel do private(i,j,k)
        DO K=1,N
           DO J=1,N
              DO I=0,N2
                 X(I,J,K)  = CMPLX(0.0,0.0)
                 NR(I,J,K) = 0
              ENDDO
           ENDDO
        ENDDO

        ALLOCATE(PROJ(NS,NS),BI(0:N2,N), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           MWANT = NS*NS + (N2+1)*N
           CALL ERRT(46,'BP NF, PROJ, BI',MWANT)
           GOTO 999
        ENDIF

        DO    K=1,NANG
           !! write(6,*)  *,' Projection #',k,nlet

C          OPEN DESIRED FILE
           CALL FILGET(FINPAT,FINPIC,NLET,ILIST(K),IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPROJ,'O',IFORM,NSAM,NSAM,
     &                  NSL,MAXIM,'DUMMY',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 999

           CALL READV(INPROJ,PROJ,NS,NS,NS,NS,1)
           CLOSE(INPROJ)

C          PAD TO DOUBLE SIZE
           CALL PADD2(PROJ,NS,BI,LSD,N)

C          FOURIER TRANSFORM
           INV = +1
           CALL FMRS_2(BI,N,N,INV)

c$omp      parallel do private(i,j)
           DO J=1,N
              DO I=0,N2
                 BI(I,J) = BI(I,J) * (-1)**(I+J+1)
              ENDDO
           ENDDO
C
           DO ISYM=1,MAXSYM
              IF (MAXSYM .GT. 1)  THEN
C                SYMMETRIES, MULTIPLY MATRICES
                 DMS = MATMUL(SM(:,:,ISYM),DM(:,:,K))
              ELSE
                 DMS = DM(:,:,K)
              ENDIF
C,schedule(static)
c$omp         parallel do private(j),shared(N,N2,JT,X,NR,BI,DMS)
              DO J=-N2+1,N2
                 CALL ONELINENN(J,N,N2,X,NR,BI,DMS)
              ENDDO
           ENDDO   ! END OF SYMMETRIES LOOP
        ENDDO      ! END OF PROJECTIONS LOOP

C       SYMMETRIZE PLANE: 0
        CALL SYMPLANEI(X,NR,N2,N)

C       CALCULATE REAL SPACE VOLUME
        CALL NORMN4(X,NR,N2,N)

C       WINDOW?
        CALL WINDUM(X,X,NS,LSD,N)

999     IF (ALLOCATED(NR))   DEALLOCATE (NR)
        IF (ALLOCATED(PROJ)) DEALLOCATE (PROJ)
        IF (ALLOCATED(BI))   DEALLOCATE (BI)
        END


C       --------------------- ONELINENN ---------------------------------

        SUBROUTINE  ONELINENN(J,N,N2,X,NR,BI,DM)

        DIMENSION      :: NR(0:N2,N,N)
        COMPLEX        :: BI(0:N2,N),X(0:N2,N,N),BTQ
        DIMENSION      :: DM(6)

        IF (J .GE. 0) THEN
           JP = J+1
        ELSE
           JP = N+J+1
        ENDIF

        DO  I=0,N2
           IF ((I*I+J*J.LT.N*N/4) .AND. .NOT.(I.EQ.0.AND.J.LT.0)) THEN
              XNEW = I*DM(1)+J*DM(4)
              YNEW = I*DM(2)+J*DM(5)
              ZNEW = I*DM(3)+J*DM(6)
              IF (XNEW .LT. 0.0)  THEN
                 XNEW = -XNEW
                 YNEW = -YNEW
                 ZNEW = -ZNEW
                 BTQ  = CONJG(BI(I,JP))
              ELSE
                 BTQ = BI(I,JP)
              ENDIF

              IXN = IFIX(XNEW+0.5+N) - N
              IYN = IFIX(YNEW+0.5+N) - N
              IZN = IFIX(ZNEW+0.5+N) - N
              IF (IXN.LE.N2 .AND.
     &            IYN.GE.-N2.AND.IYN.LE.N2 .AND.
     &            IZN.GE.-N2.AND.IZN.LE.N2) THEN
                 IF (IXN .GE. 0) THEN
                    IF (IZN .GE. 0) THEN
                       IZA = IZN+1
                    ELSE
                       IZA = N+IZN+1
                    ENDIF

                    IF (IYN .GE. 0) THEN
                       IYA = IYN+1
                    ELSE
                       IYA = N+IYN+1
                    ENDIF

                    X(IXN,IYA,IZA)  = X(IXN,IYA,IZA)+BTQ
                    NR(IXN,IYA,IZA) = NR(IXN,IYA,IZA)+1
                 ELSE
                    IF (IZN .GT. 0)  THEN
                       IZT = N-IZN+1
                    ELSE
                       IZT = -IZN+1
                    ENDIF

                    IF (IYN .GT. 0) THEN
                       IYT = N- IYN + 1
                    ELSE
                       IYT = -IYN + 1
                    ENDIF

                    X(-IXN,IYT,IZT)  = X(-IXN,IYT,IZT)+CONJG(BTQ)
                    NR(-IXN,IYT,IZT) = NR(-IXN,IYT,IZT)+1
                 ENDIF
              ENDIF
           ENDIF
        ENDDO   !  END J-I LOOP

        END

C       ------------------- WINDUM -------------------------------

        SUBROUTINE WINDUM(BI,R,L,LSD,N)

        DIMENSION  R(L,L,L),BI(LSD,N,N)


        IP = (N-L) / 2 + MOD(L,2)
        DO K=1,L
           DO J=1,L
              DO I=1,L
                 R(I,J,K) = BI(IP+I,IP+J,IP+K)
              ENDDO
           ENDDO
        ENDDO

        L2  = (L/2)**2
        L2P = (L/2-1)**2
        IP  = L / 2 + 1
        TNR = 0.0
        M   = 0

        DO K=1,L
           DO J=1,L
              DO I=1,L
                 LR = (K-IP)**2+(J-IP)**2+(I-IP)**2
                 IF (LR .LE. L2) THEN
                    IF (LR.GE.L2P .AND. LR.LE.L2) THEN
                       TNR = TNR + R(I,J,K)
                       M   = M+1
                    ENDIF
                 ENDIF
              ENDDO
           ENDDO
        ENDDO

        TNR = TNR/REAL(M)
c$omp   parallel do private(i,j,k,lr)
        DO  K=1,L
           DO  J=1,L
              DO  I=1,L
                 LR = (K-IP)**2 + (J-IP)**2 + (I-IP)**2
                 IF (LR .LE. L2) THEN
                    R(I,J,K) = R(I,J,K) - TNR
                 ELSE
                    R(I,J,K) = 0.0
                 ENDIF
              ENDDO
           ENDDO
        ENDDO

        END


C       ------------------- NORMN4 -------------------------------

        SUBROUTINE  NORMN4(X,NR,N2,N)

        DIMENSION  :: NR(0:N2,N,N)
        COMPLEX    :: X(0:N2,N,N)

c$omp   parallel do private(i,j,k)
        DO K=1,N
           DO J=1,N
              DO I=0,N2
                 IF (NR(I,J,K).GT.0)  
     &               X(I,J,K) = X(I,J,K) * (-1)**(I+J+K) / NR(I,J,K)
              ENDDO
           ENDDO
        ENDDO

C       FOURIER BACK TRANSFORM
        INV = -1
        CALL FMRS_3(X,N,N,N,INV)

        END

C       ----------------SYMPLANEI ---------------------------------------
 
        SUBROUTINE  SYMPLANEI(X,W,N2,N)

C       POSSIBLE PURPOSE??: CaLCULATE WIENER SUMMATION FROM THE 
C       INSERTED 2D SLICE PUT THE SUMMATION INTO 3D GRIDS USING 
C       NEAREST NEIGHBOUR APPROXIMATION 

        INTEGER  :: W(0:N2,N,N)
        COMPLEX  :: X(0:N2,N,N)

C       SYMMETRIZE PLANE 0
        DO  IZA=2,N2
           DO  IYA=2,N2
              X(0,IYA,IZA) = X(0,IYA,IZA)+CONJG(X(0,N-IYA+2,N-IZA+2))
              W(0,IYA,IZA) = W(0,IYA,IZA)+W(0,N-IYA+2,N-IZA+2)
              X(0,N-IYA+2,N-IZA+2) = CONJG(X(0,IYA,IZA))
              W(0,N-IYA+2,N-IZA+2) = W(0,IYA,IZA)
              X(0,N-IYA+2,IZA)=X(0,N-IYA+2,IZA)+CONJG(X(0,IYA,N-IZA+2))
              W(0,N-IYA+2,IZA) = W(0,N-IYA+2,IZA)+W(0,IYA,N-IZA+2)
              X(0,IYA,N-IZA+2) = CONJG(X(0,N-IYA+2,IZA))
              W(0,IYA,N-IZA+2) = W(0,N-IYA+2,IZA)
           ENDDO
        ENDDO

        DO  IYA=2,N2
           X(0,IYA,1)     = X(0,IYA,1)+CONJG(X(0,N-IYA+2,1))
           W(0,IYA,1)     = W(0,IYA,1)+W(0,N-IYA+2,1)
           X(0,N-IYA+2,1) = CONJG(X(0,IYA,1))
           W(0,N-IYA+2,1) = W(0,IYA,1)
        ENDDO

        DO  IZA=2,N2
           X(0,1,IZA)     = X(0,1,IZA)+CONJG(X(0,1,N-IZA+2))
           W(0,1,IZA)     = W(0,1,IZA)+W(0,1,N-IZA+2)
           X(0,1,N-IZA+2) = CONJG(X(0,1,IZA))
           W(0,1,N-IZA+2) = W(0,1,IZA)
        ENDDO

        END
