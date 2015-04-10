C++*********************************************************************
C
C  VAR3D.F                                         05/20/02
C              FMRS_PLAN                           MAY  08 ARDEAN LEITH
C  
C **********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2002, P. A. Penczek                                   *
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
C  Selection without replacement (M out of N)
C   VAR3DQSWOR
C  Selection without replacement (M out of N), Version NN
C   VAR3DNNWOR
C  Selection with replacements
C   VAR3DQSWR
C  Selection with replacements, Version NN
C   VAR3DNNWR
C  Selection with replacements per angular direction
C   VAR3DQSWRA
C  Selection with replacements per angular direction, NN version
C   VAR3DNNWRA 
C  Jackknife method (leave one out)
C   VAR3DQJACK
C
C NOTE:    SLOPPILY WRITTEN.  CLONES ITSELF.
C          NOBODY HAS EVER MANAGED TO GET THIS OPERATION TO WORK. al
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************
C Replaced bessel functions

        SUBROUTINE VAR3D(MODE)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'F90ALLOC.INC'

        REAL, DIMENSION(:,:), ALLOCATABLE :: DM,SM
C       DOC FILE POINTERS
        REAL, DIMENSION(:,:), POINTER :: ANGBUF, ANGSYM
C        Auxiliary files for option A
        INTEGER, DIMENSION(:), ALLOCATABLE :: NUMANG,KEYNUMANG,LISTANG

        COMMON /F_SPEC/ FINPAT,NLET,FINPIC
        CHARACTER*80    FINPIC,FINPAT,FILNAM,ANGDOC
        CHARACTER*1     MODE
        DATA  INPIC/99/

        NILMAX = NIMAX

        CALL FILELIST(.TRUE.,INPIC,FINPAT,NLET,INUMBR,NILMAX,NANG,
     &                 'TEMPLATE FOR 2-D IMAGES',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        MAXNUM = MAXVAL(INUMBR(1:NANG))
        CLOSE(INPIC)

C       N    - LINEAR DIMENSION OF PROJECTIONS AND RESTORED CUBE
C       NANG - TOTAL NUMBER OF IMAGES
        WRITE(NOUT,2001) NANG
2001    FORMAT(' NUMBER OF IMAGES =',I5)

C       RETRIEVE ARRAY WITH ANGLES DATA IN IT
        MAXXT = 4
        MAXYT = MAXNUM
        CALL GETDOCDAT('ANGLES DOC',.TRUE.,ANGDOC,77,.FALSE.,MAXXT,
     &                       MAXYT,ANGBUF,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9998

C       RETRIEVE ARRAY WITH SYMMETRIES DATA IN IT
        MAXXS=0
        MAXSYM=0
        CALL GETDOCDAT('SYMMETRIES DOC',.TRUE.,ANGDOC,77,.TRUE.,MAXXS,
     &                   MAXSYM,ANGSYM,IRTFLG)
        IF(IRTFLG.NE.0)  MAXSYM=1

C       OPEN FIRST IMAGE FILE TO DETERMINE NSAM, NROW, NSL
        CALL FILGET(FINPAT,FINPIC,NLET,INUMBR(1),INTFLG)

        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,NSAM,NROW,NSL,
     &             MAXIM,'DUMMY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CLOSE(INPIC)


	IF(MODE.EQ.'N'.OR.MODE.EQ.'Q')  THEN
		N2 = 4*NSAM
	ELSE
		N2 = 2*NSAM
	ENDIF
        LSD    = N2+2-MOD(N2,2)
        NMAT   = LSD*N2*N2

        ALLOCATE(DM(9,NANG), STAT=IRTFLG)
        CALL BUILDM(INUMBR,DM,NANG,ANGBUF(1,1),.FALSE.,SSDUM,.FALSE.,
     &              IRTFLG)

        IF (MODE.EQ.'A' .OR. MODE.EQ.'N')  THEN
C          Build list of unique directions & how many projections per direction.
           ALLOCATE(NUMANG(NANG),KEYNUMANG(NANG), STAT=IRTFLG)

           CALL  TIEANG(INUMBR,NANG,NUMANG,ANGBUF(1,1),NUMDIR,IRTFLG)
		  IF (NUMDIR==NANG)  THEN
			PRINT  *,'NO MULTIPLE PROJECTIONS!'
C			CALL  ERRT(???)
			STOP
		  ENDIF
        	 ALLOCATE(LISTANG(NUMDIR), STAT=IRTFLG)
		 LISTANG=0
		 DO  I=1,NANG
		   KEYNUMANG(I)=I
		   LISTANG(NUMANG(I))=LISTANG(NUMANG(I))+1
		 ENDDO
C  Sort according to assigned numbers, the original numbers will be in KEYNUMANG
		 CALL ISORT2(NUMANG,KEYNUMANG,NANG)
		ENDIF
        DEALLOCATE(ANGBUF)
        IF (IRTFLG .NE. 0) GOTO 9998

        IF(MAXSYM.GT.1)  THEN
           ALLOCATE(SM(9,MAXSYM), STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN
              CALL ERRT(46,'VA 3, SM',IER)
              DEALLOCATE (DM)
           ENDIF
           CALL  BUILDS(SM,MAXSYM,ANGSYM(1,1),IRTFLG)
           DEALLOCATE(ANGSYM)
        ELSE
           ALLOCATE(SM(1,1), STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN
              CALL ERRT(46,'VA 3, SM-2nd',IER)
              DEALLOCATE (DM)
           ENDIF
        ENDIF

        IF(MODE.EQ.'O')  THEN
C  Selection without replacement (M out of N)
         CALL VAR3DQSWOR(NSAM,LSD,N2,N2/2,INUMBR,DM,NANG,SM,MAXSYM)
        ELSEIF(MODE.EQ.'M')  THEN
C  Selection without replacement (M out of N), Version NN
         CALL VAR3DNNWOR(NSAM,LSD,N2,N2/2,INUMBR,DM,NANG,SM,MAXSYM)
        ELSEIF(MODE.EQ.'W')  THEN
C  Selection with replacements
         CALL VAR3DQSWR(NSAM,LSD,N2,N2/2,INUMBR,DM,NANG,SM,MAXSYM)
        ELSEIF(MODE.EQ.'Q')  THEN
C  Selection with replacements, Version NN
         CALL VAR3DNNWR(NSAM,LSD,N2,N2/2,INUMBR,DM,NANG,SM,MAXSYM)
        ELSEIF(MODE.EQ.'B')  THEN
C  Selection with replacements, Version NN  Same as Q, but buffer on disk
         CALL VAR3DNNWB(NSAM,LSD,N2,N2/2,INUMBR,DM,NANG,SM,MAXSYM)
        ELSEIF(MODE.EQ.'A')  THEN
C  Selection with replacements per angular direction
         CALL VAR3DQSWRA(NSAM,LSD,N2,N2/2,INUMBR,DM,NANG,SM,MAXSYM,
     &		NUMANG,KEYNUMANG,LISTANG,NUMDIR)
	 DEALLOCATE(NUMANG,KEYNUMANG,LISTANG)
        ELSEIF(MODE.EQ.'N')  THEN
C  Selection with replacements per angular direction, NN version
         CALL VAR3DNNWRA(NSAM,LSD,N2,N2/2,INUMBR,DM,NANG,SM,MAXSYM,
     &		NUMANG,KEYNUMANG,LISTANG,NUMDIR)
	 DEALLOCATE(NUMANG,KEYNUMANG,LISTANG)
        ELSEIF(MODE.EQ.'J')  THEN
C  Jackknife method (leave one out)
         CALL VAR3DQJACK(NSAM,LSD,N2,N2/2,INUMBR,DM,NANG,SM,MAXSYM)
        ENDIF

9997    DEALLOCATE(DM)
9998    DEALLOCATE(SM)
        END

CPPPPPPPPP------------------ VAR3DQSWOR ----------------------------------
C                         without replacements (M out of N)
        SUBROUTINE  VAR3DQSWOR(NS,LSD,N,N2,ILIST,DM,NANG,SM,MAXSYM)

        INCLUDE 'CMBLOCK.INC'
        DIMENSION         ILIST(NANG)
        DIMENSION         DM(3,3,NANG),SM(3,3,MAXSYM),DMS(3,3)

C  2D Fourier transforms of the input data
        COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: X,BI,XX
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: W,WW,PROJ
C  Additional matrices
        INTEGER, DIMENSION(:), ALLOCATABLE :: NORD

        COMMON  /F_SPEC/  FINPAT,NLET,FINPIC
        CHARACTER*80      FINPIC,FINPAT,FILNAM
        CHARACTER*1  NULL
        DOUBLE PRECISION  PI
        PARAMETER         (LTAB=4999)
        COMMON  /TABS/    LN2,FLTB,TABI(0:LTAB)
C In this version the order of the Bessel function is mmm=1
        COMMON  /BESSEL_PARAM/  ALPHA,AAAA,NNN
C,mmm
        DATA  IOPIC/98/,INPROJ/99/
        PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
        PARAMETER (TWOPI = 2*QUADPI)

        MONO(K1,K2)=MIN0(K1,K2)+((MAX0(K1,K2)-1)*(MAX0(K1,K2)-2)/2)

        CALL RDPRMI(MREM,NITER,NOT_USED,
     &          'Select M, number of volumes to be created')
C
C K=6
        LN=5
        LN2=LN/2
C Generalized Kaiser-Bessel window according to Lewitt
C M=NS, N=N
        R=NS/2
        V=REAL(LN-1)/2.0/REAL(N)
        ALPHA=6.5
C       AAAA=0.0079
        AAAA=0.9*V
        NNN=3
C       mmm=1
C       GENERATE TABLE WITH INTERPOLANTS
C       B0=(SQRT(ALPHA)**mmm)*BESSI(mmm,ALPHA)
        B0=SQRT(ALPHA)*BESI1(ALPHA)
        FLTB=REAL(LTAB)/REAL(LN2+1)
C  Cannot be parallel as there are DATA satements in BESI1
cc$omp parallel do private(i,s,x),shared(mmm)
        DO  I=0,LTAB
         S=REAL(I)/FLTB/N
         IF(S.LE.AAAA)  THEN
          XXX=SQRT(1.0-(S/AAAA)**2)
          TABI(I)=
C(SQRT(ALPHA*XXX)**mmm)*BESSI(mmm,ALPHA*XXX)/B0
CTEMPO=
     &          SQRT(ALPHA*XXX)*BESI1(ALPHA*XXX)/B0
         ELSE
          TABI(I)=0.0
         ENDIF
        ENDDO

        ALLOCATE (X(0:N2,N,N),W(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, W',IER)
           RETURN
        ENDIF

C       DETERMINE NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)

        ALLOCATE (PROJ(NS,NS,NUMTH),NORD(NANG),BI(0:N2,N,NANG),
     &		STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, BI',IER)
        ENDIF

        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,+1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,-1,IRTFLG)  ! UNUSED??
        IF (IRTFLG .NE. 0) RETURN

C       Read all the projection data
        DO    K=1,NANG
C          PRINT  *,' PROJECTION #',K

C          OPEN DESIRED FILE
           CALL FILGET(FINPAT,FINPIC,NLET,ILIST(K),IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPROJ,'O',IFORM,NSAM,NSAM,NSL,
     &                   MAXIM,'DUMMY',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN


           CALL READV(INPROJ,BI(0,1,K),NS,NS,NS,NS,1)
           CLOSE(INPROJ)
        ENDDO
C       Prepare all the projection data
        NS2=NS*NS
        INV = +1
        DO  IMI=1,NANG,NUMTH
c$omp     parallel do private(K,J,I),SHARED(INV)
          DO  K=IMI,MIN(NANG,IMI+NUMTH-1)
           CALL COP(BI(0,1,K),PROJ(1,1,K-IMI+1),NS2)
           CALL PADD2(PROJ(1,1,K-IMI+1),NS,BI(0,1,K),LSD,N)
           CALL FMRS_2(BI(0,1,K),N,N,INV)
           DO  J=1,N
              DO  I=0,N2
                 BI(I,J,K)=BI(I,J,K)*(-1)**(I+J+1)
              ENDDO
           ENDDO
          ENDDO
        ENDDO
        DEALLOCATE (PROJ)

        NULL=CHAR(0)
        CALL  FILERD(FINPAT,NLET,NULL,
     &          'TEMPLATE FOR OUTPUT VOLUME',IRTFLG)
        CALL RDPRMI(IT,NOT_USED,NOT_USED,'First output volume number')

#ifdef SP_MP
        LN1=LN+1
#endif
C **********************************************************************
        IF(MREM.GT.NANG/2)  THEN
C  In case the user want to select more than half projections
C     run the deleting mode
        MREM=NANG-MREM

        ALLOCATE (XX(0:N2,N,N),WW(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, WW',IER)
           RETURN
        ENDIF

C  First, build the full structure
c$omp parallel sections
c$omp section
        X=CMPLX(0.0,0.0)
c$omp section
        W=0.0
c$omp end parallel sections
C
        DO    K=1,NANG
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,K))
            ELSE
             DMS=DM(:,:,K)
            ENDIF
#ifdef SP_MP
            DO  JT=1,LN1
c$omp parallel do private(j),shared(N,N2,JT,X,W,BI,DMS)
             DO J=-N2+JT,N2,LN1
              CALL ONELINE(J,N,N2,X,W,BI(0,1,K),DMS)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
              CALL ONELINE(J,N,N2,X,W,BI(0,1,K),DMS)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO
C
C          END OF PROJECTIONS LOOP
        ENDDO

C       DO NOT SYMMETRIZE PLANE 0, WILL BE NEEDED AS IS FOR PARTIAL VOLUMES
CCCCCCCC        CALL  SYMPLANE0(X,W,N2,N)

C  Have volume X and weighting W in Fourier space.
C  Sample without replacement
        DO NK=1,NITER
C  Create new copies of X and W.
c$omp parallel sections
c$omp section
        XX=X
c$omp section
        WW=W
c$omp section
C Sample MREM without replacements
         DO MM=1,NANG
          NORD(MM)=MM
         ENDDO
         DO MM=1,NANG
          CALL RANDOM_NUMBER(HARVEST=ERND)
          IMTEMP = MIN0(NANG,MAX0(1,INT(ERND*NANG+0.5)))
          INORD=NORD(MM)
          NORD(MM)=NORD(IMTEMP)
          NORD(IMTEMP)=INORD
         ENDDO
c$omp end parallel sections
C  Go through selected projections, remove from 3D M on the list
        DO  MM=1,MREM
        LK2=NORD(MM)
C Remove LK2
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,LK2))
            ELSE
             DMS=DM(:,:,LK2)
            ENDIF
#ifdef SP_MP
            DO  JT=1,LN1
c$omp parallel do private(j),shared(N,N2,JT,XX,WW,BI,DMS)
             DO J=-N2+JT,N2,LN1
              CALL ONELINEDEL(J,N,N2,XX,WW,BI(0,1,LK2),DMS)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
              CALL ONELINEDEL(J,N,N2,XX,WW,BI(0,1,LK2),DMS)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO
        ENDDO

C       SYMMETRIZE PLANE 0
        CALL  SYMPLANE0(XX,WW,N2,N)
C  Calculate real space volume
        CALL NRMW2(XX,WW,N2,N)
        CALL WINDKB2(XX,XX,NS,LSD,N)
C
        CALL OUTVOL(XX,NS,NS,NS,NK+IT-1)
        ENDDO
        DEALLOCATE (X, W, WW, XX, BI, NORD)
C **********************************************************************
        ELSE
C  user wants to select less than half of the projections
C    run the adding mode
C  Sample without replacement
        DO NK=1,NITER
c$omp parallel sections
c$omp section
        X=CMPLX(0.0,0.0)
c$omp section
        W=0.0
c$omp section
C Sample MREM without replacements
         DO MM=1,NANG
          NORD(MM)=MM
         ENDDO
         DO MM=1,NANG
          CALL RANDOM_NUMBER(HARVEST=ERND)
          IMTEMP = MIN0(NANG,MAX0(1,INT(ERND*NANG+0.5)))
          INORD=NORD(MM)
          NORD(MM)=NORD(IMTEMP)
          NORD(IMTEMP)=INORD
         ENDDO
c$omp end parallel sections
C  Go through selected projections, add to 3D M on the list
        DO  MM=1,MREM
        LK2=NORD(MM)
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,LK2))
            ELSE
             DMS=DM(:,:,LK2)
            ENDIF
#ifdef SP_MP
            DO  JT=1,LN1
c$omp parallel do private(j),shared(N,N2,JT,X,W,BI,DMS)
             DO J=-N2+JT,N2,LN1
              CALL ONELINE(J,N,N2,X,W,BI(0,1,LK2),DMS)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
              CALL ONELINE(J,N,N2,X,W,BI(0,1,LK2),DMS)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO
C
C          END OF PROJECTIONS LOOP
        ENDDO

C       SYMMETRIZE PLANE 0
        CALL  SYMPLANE0(X,W,N2,N)
C  Calculate real space volume
        CALL NRMW2(X,W,N2,N)
        CALL WINDKB2(X,X,NS,LSD,N)
C
        CALL OUTVOL(X,NS,NS,NS,NK+IT-1)
        ENDDO
        DEALLOCATE (X, W, BI, NORD)
        ENDIF

C
        END


CPPPPPPP------------------ VAR3DNNWOR ----------------------------------
C   NOT DONE YET
C                         without replacements, Version NN
        SUBROUTINE  VAR3DNNWOR(NS,LSD,N,N2,ILIST,DM,NANG,SM,MAXSYM)

        INCLUDE 'CMBLOCK.INC'
        DIMENSION         ILIST(NANG)
        DIMENSION         DM(3,3,NANG),SM(3,3,MAXSYM),DMS(3,3)

C  2D Fourier transforms of the input data
        COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: X,BI,XX
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: W,WW,PROJ
C  Additional matrices
        INTEGER, DIMENSION(:), ALLOCATABLE :: NORD

        COMMON  /F_SPEC/  FINPAT,NLET,FINPIC
        CHARACTER*80      FINPIC,FINPAT,FILNAM
        CHARACTER*1  NULL
        DOUBLE PRECISION  PI
        PARAMETER         (LTAB=4999)
        COMMON  /TABS/    LN2,FLTB,TABI(0:LTAB)
C In this version the order of the Bessel function is mmm=1
        COMMON  /BESSEL_PARAM/  ALPHA,AAAA,NNN
C,mmm
        DATA  IOPIC/98/,INPROJ/99/
        PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
        PARAMETER (TWOPI = 2*QUADPI)

        MONO(K1,K2)=MIN0(K1,K2)+((MAX0(K1,K2)-1)*(MAX0(K1,K2)-2)/2)

        CALL RDPRMI(MREM,NITER,NOT_USED,
     &          'Select M, number of volumes to be created')
C
C K=6
        LN=5
        LN2=LN/2
C Generalized Kaiser-Bessel window according to Lewitt
C M=NS, N=N
        R=NS/2
        V=REAL(LN-1)/2.0/REAL(N)
        ALPHA=6.5
C       AAAA=0.0079
        AAAA=0.9*V
        NNN=3
C       mmm=1
C       GENERATE TABLE WITH INTERPOLANTS
C       B0=(SQRT(ALPHA)**mmm)*BESSI(mmm,ALPHA)
        B0=SQRT(ALPHA)*BESI1(ALPHA)
        FLTB=REAL(LTAB)/REAL(LN2+1)
C  Cannot be parallel as there are DATA satements in BESI1
cc$omp parallel do private(i,s,x),shared(mmm)
        DO  I=0,LTAB
         S=REAL(I)/FLTB/N
         IF(S.LE.AAAA)  THEN
          XXX=SQRT(1.0-(S/AAAA)**2)
          TABI(I)=
C(SQRT(ALPHA*XXX)**mmm)*BESSI(mmm,ALPHA*XXX)/B0
CTEMPO=
     &          SQRT(ALPHA*XXX)*BESI1(ALPHA*XXX)/B0
         ELSE
          TABI(I)=0.0
         ENDIF
        ENDDO

        ALLOCATE (X(0:N2,N,N),W(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, W',IER)
           RETURN
        ENDIF

C       DETERMINE NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)

        ALLOCATE (PROJ(NS,NS,NUMTH),NORD(NANG),BI(0:N2,N,NANG),
     &		STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, BI',IER)
        ENDIF

        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,+1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,-1,IRTFLG)  ! UNUSED??
        IF (IRTFLG .NE. 0) RETURN

C
C  Read all the projection data
        DO    K=1,NANG
C          PRINT  *,' PROJECTION #',K

C          OPEN DESIRED FILE
           CALL FILGET(FINPAT,FINPIC,NLET,ILIST(K),IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPROJ,'O',IFORM,NSAM,NSAM,NSL,
     &                   MAXIM,'DUMMY',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN


           CALL READV(INPROJ,BI(0,1,K),NS,NS,NS,NS,1)
           CLOSE(INPROJ)
        ENDDO
C  Prepare all the projection data
        NS2=NS*NS
        INV = +1
        DO  IMI=1,NANG,NUMTH
c$omp parallel do private(K,J,I),SHARED(INV)
          DO  K=IMI,MIN(NANG,IMI+NUMTH-1)
           CALL COP(BI(0,1,K),PROJ(1,1,K-IMI+1),NS2)
           CALL PADD2(PROJ(1,1,K-IMI+1),NS,BI(0,1,K),LSD,N)
           CALL FMRS_2(BI(0,1,K),N,N,INV)
           DO  J=1,N
              DO  I=0,N2
                 BI(I,J,K)=BI(I,J,K)*(-1)**(I+J+1)
              ENDDO
           ENDDO
          ENDDO
        ENDDO
        DEALLOCATE (PROJ)

        NULL=CHAR(0)
        CALL  FILERD(FINPAT,NLET,NULL,
     &          'TEMPLATE FOR OUTPUT VOLUME',IRTFLG)
        CALL RDPRMI(IT,NOT_USED,NOT_USED,'First output volume number')

#ifdef SP_MP
        LN1=LN+1
#endif
C **********************************************************************
        IF(MREM.GT.NANG/2)  THEN
C  In case the user want to select more than half projections
C     run the deleting mode
        MREM=NANG-MREM

        ALLOCATE (XX(0:N2,N,N),WW(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, WW',IER)
           RETURN
        ENDIF

C  First, build the full structure
c$omp parallel sections
c$omp section
        X=CMPLX(0.0,0.0)
c$omp section
        W=0.0
c$omp end parallel sections
C
        DO    K=1,NANG
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,K))
            ELSE
             DMS=DM(:,:,K)
            ENDIF
#ifdef SP_MP
            DO  JT=1,LN1
c$omp parallel do private(j),shared(N,N2,JT,X,W,BI,DMS)
             DO J=-N2+JT,N2,LN1
              CALL ONELINE(J,N,N2,X,W,BI(0,1,K),DMS)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
              CALL ONELINE(J,N,N2,X,W,BI(0,1,K),DMS)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO
C
C          END OF PROJECTIONS LOOP
        ENDDO

C       DO NOT SYMMETRIZE PLANE 0, WILL BE NEEDED AS IS FOR PARTIAL VOLUMES
CCCCCCCC        CALL  SYMPLANE0(X,W,N2,N)

C  Have volume X and weighting W in Fourier space.
C  Sample without replacement
        DO NK=1,NITER
C  Create new copies of X and W.
c$omp parallel sections
c$omp section
        XX=X
c$omp section
        WW=W
c$omp section
C Sample MREM without replacements
         DO MM=1,NANG
          NORD(MM)=MM
         ENDDO
         DO MM=1,NANG
          CALL RANDOM_NUMBER(HARVEST=ERND)
          IMTEMP = MIN0(NANG,MAX0(1,INT(ERND*NANG+0.5)))
          INORD=NORD(MM)
          NORD(MM)=NORD(IMTEMP)
          NORD(IMTEMP)=INORD
         ENDDO
c$omp end parallel sections
C  Go through selected projections, remove from 3D M on the list
        DO  MM=1,MREM
        LK2=NORD(MM)
C Remove LK2
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,LK2))
            ELSE
             DMS=DM(:,:,LK2)
            ENDIF
#ifdef SP_MP
            DO  JT=1,LN1
c$omp parallel do private(j),shared(N,N2,JT,XX,WW,BI,DMS)
             DO J=-N2+JT,N2,LN1
              CALL ONELINEDEL(J,N,N2,XX,WW,BI(0,1,LK2),DMS)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
              CALL ONELINEDEL(J,N,N2,XX,WW,BI(0,1,LK2),DMS)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO
        ENDDO

C       SYMMETRIZE PLANE 0
        CALL  SYMPLANE0(XX,WW,N2,N)
C  Calculate real space volume
        CALL NRMW2(XX,WW,N2,N)
        CALL WINDKB2(XX,XX,NS,LSD,N)
C
        CALL OUTVOL(XX,NS,NS,NS,NK+IT-1)
        ENDDO
        DEALLOCATE (X, W, WW, XX, BI, NORD)
C **********************************************************************
        ELSE
C  user wants to select less than half of the projections
C    run the adding mode
C  Sample without replacement
        DO NK=1,NITER
c$omp parallel sections
c$omp section
        X=CMPLX(0.0,0.0)
c$omp section
        W=0.0
c$omp section
C Sample MREM without replacements
         DO MM=1,NANG
          NORD(MM)=MM
         ENDDO
         DO MM=1,NANG
          CALL RANDOM_NUMBER(HARVEST=ERND)
          IMTEMP = MIN0(NANG,MAX0(1,INT(ERND*NANG+0.5)))
          INORD=NORD(MM)
          NORD(MM)=NORD(IMTEMP)
          NORD(IMTEMP)=INORD
         ENDDO
c$omp end parallel sections
C  Go through selected projections, add to 3D M on the list
        DO  MM=1,MREM
        LK2=NORD(MM)
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,LK2))
            ELSE
             DMS=DM(:,:,LK2)
            ENDIF
#ifdef SP_MP
            DO  JT=1,LN1
c$omp parallel do private(j),shared(N,N2,JT,X,W,BI,DMS)
             DO J=-N2+JT,N2,LN1
              CALL ONELINE(J,N,N2,X,W,BI(0,1,LK2),DMS)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
              CALL ONELINE(J,N,N2,X,W,BI(0,1,LK2),DMS)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO
C
C          END OF PROJECTIONS LOOP
        ENDDO

C       SYMMETRIZE PLANE 0
        CALL  SYMPLANE0(X,W,N2,N)
C  Calculate real space volume
        CALL NRMW2(X,W,N2,N)
        CALL WINDKB2(X,X,NS,LSD,N)
C
        CALL OUTVOL(X,NS,NS,NS,NK+IT-1)
        ENDDO
        DEALLOCATE (X, W, BI, NORD)
        ENDIF

C
        END
CPPPPPPP------------------ VAR3DQSWR ----------------------------------
C  Selection with replacements

        SUBROUTINE  VAR3DQSWR(NS,LSD,N,N2,ILIST,DM,NANG,SM,MAXSYM)

        INCLUDE 'CMBLOCK.INC'
        DIMENSION         ILIST(NANG)
        DIMENSION         DM(3,3,NANG),SM(3,3,MAXSYM),DMS(3,3)

C  2D Fourier transforms of the input data
        COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: X,BI,XX
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: W,WW,PROJ
C  Additional matrices
        INTEGER, DIMENSION(:), ALLOCATABLE :: NORD
        REAL, DIMENSION(:), ALLOCATABLE :: ERND

        COMMON  /F_SPEC/  FINPAT,NLET,FINPIC
        CHARACTER*80      FINPIC,FINPAT,FILNAM
        CHARACTER*1  NULL
        DOUBLE PRECISION  PI
        PARAMETER         (LTAB=4999)
        COMMON  /TABS/    LN2,FLTB,TABI(0:LTAB)
C In this version the order of the Bessel function is mmm=1
        COMMON  /BESSEL_PARAM/  ALPHA,AAAA,NNN
C,mmm
        DATA  IOPIC/98/,INPROJ/99/
        PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
        PARAMETER (TWOPI = 2*QUADPI)

        MONO(K1,K2)=MIN0(K1,K2)+((MAX0(K1,K2)-1)*(MAX0(K1,K2)-2)/2)

        CALL RDPRMI(NITER,NOT_USED,NOT_USED,
     &          'Number of volumes to be created')
C
C K=6
        LN=5
        LN2=LN/2
C Generalized Kaiser-Bessel window according to Lewitt
C M=NS, N=N
        R=NS/2
        V=REAL(LN-1)/2.0/REAL(N)
        ALPHA=6.5
C       AAAA=0.0079
        AAAA=0.9*V
        NNN=3
C       mmm=1
C       GENERATE TABLE WITH INTERPOLANTS
C       B0=(SQRT(ALPHA)**mmm)*BESSI(mmm,ALPHA)
        B0=SQRT(ALPHA)*BESI1(ALPHA)
        FLTB=REAL(LTAB)/REAL(LN2+1)
C  Cannot be parallel as there are DATA satements in BESI1
cc$omp parallel do private(i,s,x),shared(mmm)
        DO  I=0,LTAB
         S=REAL(I)/FLTB/N
         IF(S.LE.AAAA)  THEN
          XXX=SQRT(1.0-(S/AAAA)**2)
          TABI(I)=
C(SQRT(ALPHA*XXX)**mmm)*BESSI(mmm,ALPHA*XXX)/B0
CTEMPO=
     &          SQRT(ALPHA*XXX)*BESI1(ALPHA*XXX)/B0
         ELSE
          TABI(I)=0.0
         ENDIF
        ENDDO

#ifdef SP_MP
        LN1=LN+1
#endif
        ALLOCATE (X(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, X',IER)
        ENDIF

       ALLOCATE (W(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, W',IER)
           RETURN
        ENDIF

C       DETERMINE NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)

        ALLOCATE (PROJ(NS,NS,NUMTH), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, PROJ',IER)
           DEALLOCATE (W)
           RETURN
        ENDIF

        ALLOCATE (BI(0:N2,N,NANG), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, BI',IER)
           DEALLOCATE (W,PROJ)
        ENDIF

        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,+1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,-1,IRTFLG)  ! UNUSED??
        IF (IRTFLG .NE. 0) RETURN

CC  Read all the projection data
        DO    K=1,NANG
C          PRINT  *,' PROJECTION #',K

C          OPEN DESIRED FILE
           CALL FILGET(FINPAT,FINPIC,NLET,ILIST(K),IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPROJ,'O',IFORM,NSAM,NSAM,NSL,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN


           CALL READV(INPROJ,BI(0,1,K),NS,NS,NS,NS,1)
           CLOSE(INPROJ)
        ENDDO
C  Prepare all the projection data
        NS2=NS*NS
        INV = +1
        DO  IMI=1,NANG,NUMTH
c$omp parallel do private(K,J,I),SHARED(INV)
          DO  K=IMI,MIN(NANG,IMI+NUMTH-1)
           CALL  COP(BI(0,1,K),PROJ(1,1,K-IMI+1),NS2)
           CALL PADD2(PROJ(1,1,K-IMI+1),NS,BI(0,1,K),LSD,N)
           CALL FMRS_2(BI(0,1,K),N,N,INV)
           DO  J=1,N
              DO  I=0,N2
                 BI(I,J,K)=BI(I,J,K)*(-1)**(I+J+1)
              ENDDO
           ENDDO
          ENDDO
        ENDDO
        DEALLOCATE (PROJ)

        NULL=CHAR(0)
        CALL  FILERD(FINPAT,NLET,NULL,
     &          'TEMPLATE FOR OUTPUT VOLUME',IRTFLG)
        CALL RDPRMI(IT,NOT_USED,NOT_USED,'First output volume number')

C  First, build the full structure
c$omp parallel sections
c$omp section
        X=CMPLX(0.0,0.0)
c$omp section
        W=0.0
c$omp end parallel sections
C
        DO    K=1,NANG
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,K))
            ELSE
             DMS=DM(:,:,K)
            ENDIF
#ifdef SP_MP
            DO  JT=1,LN1
c$omp parallel do private(j),shared(N,N2,JT,X,W,BI,DMS)
             DO J=-N2+JT,N2,LN1
              CALL ONELINE(J,N,N2,X,W,BI(0,1,K),DMS)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
              CALL ONELINE(J,N,N2,X,W,BI(0,1,K),DMS)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO
C
C          END OF PROJECTIONS LOOP
        ENDDO

C       DO NOT SYMMETRIZE PLANE 0, WILL BE NEEDED AS IS FOR PARTIAL VOLUMES
CCCCCCCC        CALL  SYMPLANE0(X,W,N2,N)

C  Have volume X and weighting W in Fourier space.
C  Sample with replacement
        ALLOCATE (ERND(NANG),NORD(NANG), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'BP 3F, NANG',IER)
           RETURN
        ENDIF

        ALLOCATE (XX(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, XX',IER)
           RETURN
        ENDIF

        ALLOCATE (WW(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, WW',IER)
           RETURN
        ENDIF

C Endless loop
        DO NK=1,NITER
c$omp parallel sections
c$omp section
         CALL RANDOM_NUMBER(HARVEST=ERND)
         NORD = MIN0(NANG,MAX0(1,INT(ERND*NANG+0.5)))
C  Sort the selection of images
        CALL ISORT(NORD,NANG)
c$omp section
C  Create new copies of X and W.
        XX=X
c$omp section
        WW=W
c$omp end parallel sections

C  Go through selected projections, remove from 3D missing,
C                                   add multiple...
        LK1=0
        LK2=0
100     LK1=LK1+1
101     LK2=LK2+1
        IF(NORD(LK1).EQ.LK2)  THEN
        MULT=1
102      IF(LK1.EQ.NANG)  THEN
          IF(MULT.GT.1) THEN
C  Add multiple LK2
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,LK2))
            ELSE
             DMS=DM(:,:,LK2)
            ENDIF
#ifdef SP_MP
            MULTI=MULT-1
            DO  JT=1,LN1
c$omp parallel do private(j),shared(N,N2,JT,XX,WW,BI,DMS)
             DO J=-N2+JT,N2,LN1
              CALL ONELINEM(J,N,N2,XX,WW,BI(0,1,LK2),DMS,MULTI)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
              CALL ONELINEM(J,N,N2,XX,WW,BI(0,1,LK2),DMS,MULT-1)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO
C

          ENDIF
C  End
C  Check how many from this series are missing at the end,
C  i.e.,  if the series is 1-12, what if 11 and 12 are missing?
	  IF(NORD(NANG).LT.NANG)  THEN
C  Remove the missing ones at the end
	   DO ID=NORD(NANG)+1,NANG
C Remove ID
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,ID))
            ELSE
             DMS=DM(:,:,ID)
            ENDIF
c$omp parallel do private(j),shared(N,N2,XX,BI,DMS)
            DO J=-N2+1,N2
             CALL ONELINEDEL(J,N,N2,XX,WW,BI(0,1,ID),DMS)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
	   
	   ENDDO
	  ENDIF
C  End and out
           GOTO  103
         ELSE
          IF(NORD(LK1).EQ.NORD(LK1+1))  THEN
          LK1=LK1+1
          MULT=MULT+1
          GOTO  102
          ELSE
           IF(MULT.GT.1) THEN
C  Add multiple LK2
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,LK2))
            ELSE
             DMS=DM(:,:,LK2)
            ENDIF
#ifdef SP_MP
            MULTI=MULT-1
            DO  JT=1,LN1
c$omp parallel do private(j),shared(N,N2,JT,XX,WW,BI,DMS)
             DO J=-N2+JT,N2,LN1
              CALL ONELINEM(J,N,N2,XX,WW,BI(0,1,LK2),DMS,MULTI)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
              CALL ONELINEM(J,N,N2,XX,WW,BI(0,1,LK2),DMS,MULT-1)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO
C

           ENDIF
C  Continue
           GOTO 100
          ENDIF
         ENDIF
        ELSE
C Remove LK2
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,LK2))
            ELSE
             DMS=DM(:,:,LK2)
            ENDIF
#ifdef SP_MP
            DO  JT=1,LN1
c$omp parallel do private(j),shared(N,N2,JT,XX,WW,BI,DMS)
             DO J=-N2+JT,N2,LN1
              CALL ONELINEDEL(J,N,N2,XX,WW,BI(0,1,LK2),DMS)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
              CALL ONELINEDEL(J,N,N2,XX,WW,BI(0,1,LK2),DMS)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO
C
        GOTO  101
        ENDIF
103     CONTINUE
C       SYMMETRIZE PLANE 0
        CALL  SYMPLANE0(XX,WW,N2,N)
C  Calculate real space volume
        CALL NRMW2(XX,WW,N2,N)
        CALL WINDKB2(XX,XX,NS,LSD,N)

        CALL OUTVOL(XX,NS,NS,NS,NK+IT-1)

        ENDDO
        DEALLOCATE (X, W, XX, WW, BI, NORD, ERND)

C
        END


CPPPPPPP------------------ VAR3DNNWR ----------------------------------
C  Selection with replacements, Version NN

        SUBROUTINE  VAR3DNNWR(NS,LSD,N,N2,ILIST,DM,NANG,SM,MAXSYM)

        INCLUDE 'CMBLOCK.INC'
        DIMENSION         ILIST(NANG)
        DIMENSION         DM(3,3,NANG),SM(3,3,MAXSYM),DMS(3,3)

C  2D Fourier transforms of the input data
        COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: BI,XX
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: PROJ
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: WW
C  Additional matrices
        INTEGER, DIMENSION(:), ALLOCATABLE :: NORD
        REAL, DIMENSION(:), ALLOCATABLE :: ERND

        COMMON  /F_SPEC/  FINPAT,NLET,FINPIC
        CHARACTER*80      FINPIC,FINPAT,FILNAM
        CHARACTER*1  NULL
        DATA  IOPIC/98/,INPROJ/99/
        CALL RDPRMI(NITER,NOT_USED,NOT_USED,
     &          'Number of volumes to be created')

C       DETERMINE NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)

        ALLOCATE (BI(0:N2,N,NANG), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, BI',IER)
        ENDIF

        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,+1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,-1,IRTFLG)  ! UNUSED??
        IF (IRTFLG .NE. 0) RETURN

C
        ALLOCATE (WW(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, WW',IER)
           DEALLOCATE (BI)
           RETURN
        ENDIF

        ALLOCATE (PROJ(NS,NS,NUMTH), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, PROJ',IER)
           DEALLOCATE (WW,BI)
           RETURN
        ENDIF

           DO  J=1,N
              DO  I=0,N2
                 WW(I,J,1)=(-1)**(I+J+1)
              ENDDO
           ENDDO
C  Read all the projection data
        DO    K=1,NANG
C          PRINT  *,' PROJECTION #',K

C          OPEN DESIRED FILE
           CALL FILGET(FINPAT,FINPIC,NLET,ILIST(K),IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPROJ,'O',IFORM,NSAM,NSAM,NSL,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN


           CALL READV(INPROJ,BI(0,1,K),NS,NS,NS,NS,1)
           CLOSE(INPROJ)
        ENDDO
C  Prepare all the projection data
        NS2=NS*NS
        INV = +1
        DO  IMI=1,NANG,NUMTH
c$omp parallel do private(K,J,I),SHARED(INV)
          DO  K=IMI,MIN(NANG,IMI+NUMTH-1)
           CALL  COP(BI(0,1,K),PROJ(1,1,K-IMI+1),NS2)
           CALL PADD2(PROJ(1,1,K-IMI+1),NS,BI(0,1,K),LSD,N)
           CALL FMRS_2(BI(0,1,K),N,N,INV)
           BI(:,:,K)=BI(:,:,K)*WW(:,:,1)
          ENDDO
        ENDDO
        DEALLOCATE (PROJ)

        NULL=CHAR(0)
        CALL  FILERD(FINPAT,NLET,NULL,
     &          'TEMPLATE FOR OUTPUT VOLUME',IRTFLG)
        CALL RDPRMI(IT,NOT_USED,NOT_USED,'First output volume number')

C  Sample with replacement
        ALLOCATE (ERND(NANG),NORD(NANG), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'BP 3F, NANG',IER)
           RETURN
        ENDIF

        ALLOCATE (XX(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, XX',IER)
           RETURN
        ENDIF

C Endless loop
        DO NIK=1,NITER
c$omp parallel sections
c$omp section
         CALL RANDOM_NUMBER(HARVEST=ERND)
         NORD = MIN0(NANG,MAX0(1,INT(ERND*NANG+0.5)))
C  Sort the selection of images
        CALL ISORT(NORD,NANG)
c$omp section
C  Reset XX and WW.
        XX=CMPLX(0.0,0.0)
c$omp section
        WW=0
c$omp end parallel sections
C  Go through selected projections, 
C          add multiple and single separately...
        LK1=0
	MULT=1
101     LK1=LK1+1
        IF(NORD(LK1).EQ.NORD(LK1+1))  THEN
C  The same entry repeated
         MULT=MULT+1
	 IF(LK1+1.EQ.NANG)  THEN
C  Add multiple NORD(LK1)
         NK=NORD(LK1)
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,NK))
            ELSE
             DMS=DM(:,:,NK)
            ENDIF
c$omp parallel do private(j),shared(N,N2,XX,WW,BI,DMS,MULT)
            DO J=-N2+1,N2
              CALL ONELINENWM(J,N,N2,XX,WW,BI(0,1,NK),DMS,MULT)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
C    Over and out
	 GOTO 102
	 ENDIF
	 GOTO 101
	ELSE
C  Different entry
	 NK=NORD(LK1)
	 IF(MULT.EQ.1)  THEN
C  Add single LK1
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,NK))
            ELSE
             DMS=DM(:,:,NK)
            ENDIF
c$omp parallel do private(j),shared(N,N2,XX,WW,BI,DMS)
            DO J=-N2+1,N2
              CALL ONELINENN(J,N,N2,XX,WW,BI(0,1,NK),DMS)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
	 ELSE
C  Add multiple LK1
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,NK))
            ELSE
             DMS=DM(:,:,NK)
            ENDIF
c$omp parallel do private(j),shared(N,N2,XX,WW,BI,DMS,MULT)
            DO J=-N2+1,N2
              CALL ONELINENWM(J,N,N2,XX,WW,BI(0,1,NK),DMS,MULT)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
	 ENDIF
	 IF(LK1+1.EQ.NANG)  THEN
C  Add single LK1+1
	 NK=NORD(LK1+1)
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,NK))
            ELSE
             DMS=DM(:,:,NK)
            ENDIF
c$omp parallel do private(j),shared(N,N2,XX,WW,BI,DMS)
            DO J=-N2+1,N2
              CALL ONELINENN(J,N,N2,XX,WW,BI(0,1,NK),DMS)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
          GOTO  102
         ENDIF
	 MULT=1
	 GOTO  101
        ENDIF
C  It should never get here....
	WRITE(NOUT,*)  'ERROR'
	STOP
C
102     CONTINUE
C  SYMMETRIZE PLANE 0
        CALL  SYMPLANEI(XX,WW,N2,N)
C  Calculate real space volume
        CALL NORMN4(XX,WW,N2,N)
        CALL WINDUM(XX,XX,NS,LSD,N)

        CALL OUTVOL(XX,NS,NS,NS,NIK+IT-1)

        ENDDO
        DEALLOCATE (XX, WW, BI, NORD, ERND)
C
        END

CPPPPPPP------------------ VAR3DNNWB ----------------------------------
C  Selection with replacements, Version NN, same as VAR3DNNWR, but buffer on disk
        SUBROUTINE  VAR3DNNWB(NS,LSD,N,N2,ILIST,DM,NANG,SM,MAXSYM)

        INCLUDE 'CMBLOCK.INC'
        DIMENSION         ILIST(NANG)
        DIMENSION         DM(3,3,NANG),SM(3,3,MAXSYM),DMS(3,3)

C  2D Fourier transforms of the input data
        COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: XX
        COMPLEX, DIMENSION(:,:),   ALLOCATABLE :: BI
        REAL,    DIMENSION(:,:),   ALLOCATABLE :: PROJ
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: WW
C  Additional matrices
        INTEGER, DIMENSION(:), ALLOCATABLE :: NORD
        REAL, DIMENSION(:), ALLOCATABLE :: ERND

        COMMON  /F_SPEC/  FINPAT,NLET,FINPIC
        CHARACTER*80      FINPIC,FINPAT,FILNAM
        CHARACTER*1  NULL
        DATA  IOPIC/98/,INPROJ/99/
        CALL RDPRMI(NITER,NOT_USED,NOT_USED,
     &          'Number of volumes to be created')

C       DETERMINE NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)

        ALLOCATE (BI(0:N2,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, BI',IER)
        ENDIF

        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,+1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,-1,IRTFLG)  ! UNUSED??
        IF (IRTFLG .NE. 0) RETURN

C        ALLOCATE (WW(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, WW',IER)
           DEALLOCATE (BI)
           RETURN
        ENDIF

        ALLOCATE (PROJ(NS,NS), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, PROJ',IER)
           DEALLOCATE (WW,BI)
           RETURN
        ENDIF
	
       DO  J=1,N
          DO  I=0,N2
             WW(I,J,1)=(-1)**(I+J+1)
          ENDDO
       ENDDO
C  Read all the projection data and store in the scratch file
       MAXIM = 0
       CALL OPFILEC(0,.TRUE.,FILNAM,IOPIC,'Z/E',IFORM,NSAMB,NROWB,NANGB,
     &                   MAXIM,'Scratch',.FALSE.,IRTFLG)
       IF (IRTFLG .EQ. 0) THEN
C  Scratch file exists, check the dimensions
         IF(NANGB.NE.NANG .OR. NSAMB.NE.2*(N2+1) .OR.NROWB.NE.N)  THEN
          CALL ERRT(1,'VA 3B, Scratch',IER)
	  CLOSE(IOPIC)
	  RETURN
	 ENDIF      
       ELSE
C Scratch file did not exist, open it
        MAXIM = 0
	IFORM=3
	NSAMB=2*(N2+1)
        CALL OPFILEC(0,.FALSE.,FILNAM,IOPIC,'N',IFORM,NSAMB,N,NANG,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
C Prepare all projections and write to the scratch file
         NS2=NS*NS
         DO    K=1,NANG
C          PRINT  *,' PROJECTION #',K

C          OPEN DESIRED FILE
           CALL FILGET(FINPAT,FINPIC,NLET,ILIST(K),IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPROJ,'O',IFORM,NSAM,NSAM,NSL,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN


           CALL READV(INPROJ,PROJ,NS,NS,NS,NS,1)
           CLOSE(INPROJ)
           CALL PADD2(PROJ,NS,BI,LSD,N)
           INV = +1
           CALL FMRS_2(BI,N,N,INV)
           BI=BI*WW(:,:,1)
	   CALL  WRTVOL(IOPIC,NSAMB,N,K,K,BI,IRTFLG)
         ENDDO
         DEALLOCATE (PROJ)
C end of creation of the scratch file, Close it and open again
	 CLOSE(IOPIC)
         MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FILNAM,IOPIC,'O',IFORM,NSAMB,NROWB,NANGB,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
     	 IF(IRTFLG.NE.0)  THEN
	  PRINT  *, 'SCRATCH FILE WRONG?'
	  RETURN
	 ENDIF
	ENDIF

        NULL=CHAR(0)
        CALL  FILERD(FINPAT,NLET,NULL,
     &          'TEMPLATE FOR OUTPUT VOLUME',IRTFLG)
        CALL RDPRMI(IT,NOT_USED,NOT_USED,'First output volume number')

C  Sample with replacement
        ALLOCATE (ERND(NANG),NORD(NANG), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'BP 3F, NANG',IER)
           RETURN
        ENDIF

        ALLOCATE (XX(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, XX',IER)
           RETURN
        ENDIF

C Endless loop
        DO NIK=1,NITER
c$omp parallel sections
c$omp section
         CALL RANDOM_NUMBER(HARVEST=ERND)
         NORD = MIN0(NANG,MAX0(1,INT(ERND*NANG+0.5)))
C  Sort the selection of images
        CALL ISORT(NORD,NANG)
c$omp section
C  Reset XX and WW.
        XX=CMPLX(0.0,0.0)
c$omp section
        WW=0
c$omp end parallel sections
C  Go through selected projections, 
C          add multiple and single separately...
        LK1=0
	MULT=1
101     LK1=LK1+1
        IF(NORD(LK1).EQ.NORD(LK1+1))  THEN
C  The same entry repeated
         MULT=MULT+1
	 IF(LK1+1.EQ.NANG)  THEN
C  Add multiple NORD(LK1)
         NK=NORD(LK1)
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,NK))
            ELSE
             DMS=DM(:,:,NK)
            ENDIF
	    CALL  REDVOL(IOPIC,NSAMB,N,NK,NK,BI,IRTFLG)
c$omp parallel do private(j),shared(N,N2,XX,WW,BI,DMS,MULT)
            DO J=-N2+1,N2
              CALL ONELINENWM(J,N,N2,XX,WW,BI,DMS,MULT)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
C    Over and out
	 GOTO 102
	 ENDIF
	 GOTO 101
	ELSE
C  Different entry
	 NK=NORD(LK1)
	 IF(MULT.EQ.1)  THEN
C  Add single LK1
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,NK))
            ELSE
             DMS=DM(:,:,NK)
            ENDIF
	    CALL  REDVOL(IOPIC,NSAMB,N,NK,NK,BI,IRTFLG)
c$omp parallel do private(j),shared(N,N2,XX,WW,BI,DMS)
            DO J=-N2+1,N2
              CALL ONELINENN(J,N,N2,XX,WW,BI,DMS)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
	 ELSE
C  Add multiple LK1
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,NK))
            ELSE
             DMS=DM(:,:,NK)
            ENDIF
	    CALL  REDVOL(IOPIC,NSAMB,N,NK,NK,BI,IRTFLG)
c$omp parallel do private(j),shared(N,N2,XX,WW,BI,DMS,MULT)
            DO J=-N2+1,N2
              CALL ONELINENWM(J,N,N2,XX,WW,BI,DMS,MULT)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
	 ENDIF
	 IF(LK1+1.EQ.NANG)  THEN
C  Add single LK1+1
	 NK=NORD(LK1+1)
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,NK))
            ELSE
             DMS=DM(:,:,NK)
            ENDIF
	    CALL  REDVOL(IOPIC,NSAMB,N,NK,NK,BI,IRTFLG)
c$omp parallel do private(j),shared(N,N2,XX,WW,BI,DMS)
            DO J=-N2+1,N2
              CALL ONELINENN(J,N,N2,XX,WW,BI,DMS)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
          GOTO  102
         ENDIF
	 MULT=1
	 GOTO  101
        ENDIF
C  It should never get here....
	WRITE(NOUT,*)  'ERROR'
	STOP
C
102     CONTINUE
C  SYMMETRIZE PLANE 0
        CALL  SYMPLANEI(XX,WW,N2,N)
C  Calculate real space volume
        CALL NORMN4(XX,WW,N2,N)
        CALL WINDUM(XX,XX,NS,LSD,N)

        CALL OUTVOL(XX,NS,NS,NS,NIK+IT-1)

        ENDDO
        DEALLOCATE (XX, WW, BI, NORD, ERND)
C  Close Scratch file
	CLOSE(IOPIC)
C
        END


#if 1
C  In this version we build averages and calculate the reconstruction
CPPPPPPP------------------ VAR3DQSWRA ----------------------------------
C  Selection with replacements per angular direction
        SUBROUTINE  VAR3DQSWRA(NS,LSD,N,N2,ILIST,DM,NANG,SM,MAXSYM,
     &		NUMANG,KEYNUMANG,LISTANG,NUMDIR)

        INCLUDE 'CMBLOCK.INC'
        DIMENSION         ILIST(NANG)
        DIMENSION         DM(3,3,NANG),SM(3,3,MAXSYM),DMS(3,3)
	DIMENSION         NUMANG(NANG),KEYNUMANG(NANG),LISTANG(NUMDIR)

C  2D Fourier transforms of the input data
        COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: BI,X
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: W
        COMPLEX, DIMENSION(:,:), ALLOCATABLE :: BT
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: PROJ
C  Additional matrices
        INTEGER, DIMENSION(:), ALLOCATABLE :: NORD
        REAL, DIMENSION(:), ALLOCATABLE :: ERND

        COMMON  /F_SPEC/  FINPAT,NLET,FINPIC
        CHARACTER*80      FINPIC,FINPAT,FILNAM
        CHARACTER*1  NULL
        DOUBLE PRECISION  PI
        PARAMETER         (LTAB=4999)
        COMMON  /TABS/    LN2,FLTB,TABI(0:LTAB)
C In this version the order of the Bessel function is mmm=1
        COMMON  /BESSEL_PARAM/  ALPHA,AAAA,NNN
C,mmm
        DATA  IOPIC/98/,INPROJ/99/
        PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
        PARAMETER (TWOPI = 2*QUADPI)

        CALL RDPRMI(NITER,NOT_USED,NOT_USED,
     &          'Number of volumes to be created')
C
C K=6
        LN=5
        LN2=LN/2
C Generalized Kaiser-Bessel window according to Lewitt
C M=NS, N=N
        R=NS/2
        V=REAL(LN-1)/2.0/REAL(N)
        ALPHA=6.5
C       AAAA=0.0079
        AAAA=0.9*V
        NNN=3
C       mmm=1
C       GENERATE TABLE WITH INTERPOLANTS
C       B0=(SQRT(ALPHA)**mmm)*BESSI(mmm,ALPHA)
        B0=SQRT(ALPHA)*BESI1(ALPHA)
        FLTB=REAL(LTAB)/REAL(LN2+1)
C  Cannot be parallel as there are DATA satements in BESI1
cc$omp parallel do private(i,s,x),shared(mmm)
        DO  I=0,LTAB
         S=REAL(I)/FLTB/N
         IF(S.LE.AAAA)  THEN
          XXX=SQRT(1.0-(S/AAAA)**2)
          TABI(I)=
C(SQRT(ALPHA*XXX)**mmm)*BESSI(mmm,ALPHA*XXX)/B0
CTEMPO=
     &          SQRT(ALPHA*XXX)*BESI1(ALPHA*XXX)/B0
         ELSE
          TABI(I)=0.0
         ENDIF
        ENDDO

        ALLOCATE (X(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, X',IER)
        ENDIF

       ALLOCATE (W(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, W',IER)
           RETURN
        ENDIF

C It will not have the parallel version
C       DETERMINE NUMBER OF OMP THREADS
C        CALL GETTHREADS(NUMTH)
	NUMTH=1

        ALLOCATE (PROJ(NS,NS,NUMTH), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, PROJ',IER)
           DEALLOCATE (W)
           RETURN
        ENDIF

        ALLOCATE (BI(0:N2,N,NANG), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, BI',IER)
           DEALLOCATE (W,PROJ)
        ENDIF

        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,+1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,-1,IRTFLG)  ! UNUSED??
        IF (IRTFLG .NE. 0) RETURN

CC  Read all the projection data
        DO    K=1,NANG
C          PRINT  *,' PROJECTION #',K

C          OPEN DESIRED FILE
           CALL FILGET(FINPAT,FINPIC,NLET,ILIST(K),IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPROJ,'O',IFORM,NSAM,NSAM,NSL,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN


           CALL READV(INPROJ,BI(0,1,K),NS,NS,NS,NS,1)
           CLOSE(INPROJ)
        ENDDO
C  Prepare all the projection data
        NS2=NS*NS
        INV = +1
        DO  IMI=1,NANG,NUMTH
c$omp parallel do private(K,J,I),SHARED(INV)
          DO  K=IMI,MIN(NANG,IMI+NUMTH-1)
           CALL  COP(BI(0,1,K),PROJ(1,1,K-IMI+1),NS2)
           CALL PADD2(PROJ(1,1,K-IMI+1),NS,BI(0,1,K),LSD,N)
           CALL FMRS_2(BI(0,1,K),N,N,INV)
           DO  J=1,N
              DO  I=0,N2
                 BI(I,J,K)=BI(I,J,K)*(-1)**(I+J+1)
              ENDDO
           ENDDO
          ENDDO
        ENDDO
        DEALLOCATE (PROJ)

        NULL=CHAR(0)
        CALL  FILERD(FINPAT,NLET,NULL,
     &          'TEMPLATE FOR OUTPUT VOLUME',IRTFLG)
        CALL RDPRMI(IT,NOT_USED,NOT_USED,'First output volume number')

C  First, build the volume of weights
        W=0.0
C
C  NPRJ is the beginning of the projections with the same direction
	NPRJ=0
	DO LNG=1,NUMDIR
C  Current number of projections to process
	LANG=LISTANG(LNG)
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,KEYNUMANG(1+NPRJ)))
            ELSE
             DMS=DM(:,:,KEYNUMANG(1+NPRJ))
            ENDIF
            DO J=-N2+1,N2
              CALL ONELINEM_W(J,N,N2,W,DMS,LANG)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
C
C
	NPRJ=NPRJ+LISTANG(LNG)
C  End of do loop over angular directions....
        ENDDO

C       SYMMETRIZE WEIGHTS W
        CALL  SYMPLANE0_WW(W,N2,N)
C
C  Have weights W in Fourier space.
C  Sample with replacement per angular direction
C  Array NORD used for resampling has length equal to the maximum number
C   of projections per direction

	INAGMX=MAXVAL(LISTANG)
        ALLOCATE (ERND(INAGMX),NORD(INAGMX), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'BP 3F, NANG',IER)
           RETURN
        ENDIF
C BT will contain current sum of projections at the same angular direction
        ALLOCATE (BT(0:N2,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           DEALLOCATE (W)
           CALL ERRT(46,'VA 3, BI',IER)
        ENDIF

C Endless loop
        DO NK=1,NITER
C Reset the new volume X.
C The weights should remain the same....            **************************
        X=CMPLX(0.0,0.0)
C Process each angular direction independently
C  NPRJ is the beginning of the projections with the same direction
	NPRJ=0
	DO LNG=1,NUMDIR
C  Current number of projections to process
	LANG=LISTANG(LNG)
	IF(LANG.GT.1)  THEN
         CALL RANDOM_NUMBER(HARVEST=ERND(1:LANG))
         NORD(1:LANG) = MIN0(LANG,MAX0(1,INT(ERND(1:LANG)*LANG+0.5)))
C  LK2 is the projection number within the current (i.e., LNG)
C   projection direction. 
C   The actual projection number is KEYNUMANG(LK2+NPRJ)
	BT=CMPLX(0.0,0.0)
	DO LK2=1,LANG
	 BT=BT+BI(:,:,KEYNUMANG(NORD(LK2)+NPRJ))
	ENDDO
	ELSE
C  There was just one projection (should not be the case...)
	 BT=BI(:,:,KEYNUMANG(1+NPRJ))
	ENDIF
C  Add BT corresponding to multiple projections at the direction LNG
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,KEYNUMANG(1+NPRJ)))
            ELSE
             DMS=DM(:,:,KEYNUMANG(1+NPRJ))
            ENDIF
            DO J=-N2+1,N2
              CALL ONELINE_X(J,N,N2,X,BT,DMS)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
C
	NPRJ=NPRJ+LISTANG(LNG)
C  End of do loop over angular directions....
	ENDDO
C       SYMMETRIZE PLANE 0
        CALL  SYMPLANE0_X(X,N2,N)
C  Calculate real space volume
        CALL NRMW2(X,W,N2,N)
        CALL WINDKB2(X,X,NS,LSD,N)

        CALL OUTVOL(X,NS,NS,NS,NK+IT-1)

        ENDDO
        DEALLOCATE (X, W, BI, BT, NORD, ERND)

C
        END

#else
C  In this version the projections are added individually
CPPPPPPP------------------ VAR3DQSWRA ----------------------------------
C  Selection with replacements per angular direction
        SUBROUTINE  VAR3DQSWRA(NS,LSD,N,N2,ILIST,DM,NANG,SM,MAXSYM,
     &		NUMANG,KEYNUMANG,LISTANG,NUMDIR)

        INCLUDE 'CMBLOCK.INC'
        DIMENSION         ILIST(NANG)
        DIMENSION         DM(3,3,NANG),SM(3,3,MAXSYM),DMS(3,3)
	DIMENSION         NUMANG(NANG),KEYNUMANG(NANG),LISTANG(NUMDIR)

C  2D Fourier transforms of the input data
        COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: X,XX
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: W
        COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: BI
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: PROJ
C  Additional matrices
        INTEGER, DIMENSION(:), ALLOCATABLE :: NORD
        REAL, DIMENSION(:), ALLOCATABLE :: ERND

        COMMON  /F_SPEC/  FINPAT,NLET,FINPIC
        CHARACTER*80      FINPIC,FINPAT,FILNAM
        CHARACTER*1  NULL
        DOUBLE PRECISION  PI
        PARAMETER         (LTAB=4999)
        COMMON  /TABS/    LN2,FLTB,TABI(0:LTAB)
C In this version the order of the Bessel function is mmm=1
        COMMON  /BESSEL_PARAM/  ALPHA,AAAA,NNN
C,mmm
        DATA  IOPIC/98/,INPROJ/99/
        PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
        PARAMETER (TWOPI = 2*QUADPI)

        MONO(K1,K2)=MIN0(K1,K2)+((MAX0(K1,K2)-1)*(MAX0(K1,K2)-2)/2)

        CALL RDPRMI(NITER,NOT_USED,NOT_USED,
     &          'Number of volumes to be created')
C
C K=6
        LN=5
        LN2=LN/2
C Generalized Kaiser-Bessel window according to Lewitt
C M=NS, N=N
        R=NS/2
        V=REAL(LN-1)/2.0/REAL(N)
        ALPHA=6.5
C       AAAA=0.0079
        AAAA=0.9*V
        NNN=3
C       mmm=1
C       GENERATE TABLE WITH INTERPOLANTS
C       B0=(SQRT(ALPHA)**mmm)*BESSI(mmm,ALPHA)
        B0=SQRT(ALPHA)*BESI1(ALPHA)
        FLTB=REAL(LTAB)/REAL(LN2+1)
C  Cannot be parallel as there are DATA satements in BESI1
cc$omp parallel do private(i,s,x),shared(mmm)
        DO  I=0,LTAB
         S=REAL(I)/FLTB/N
         IF(S.LE.AAAA)  THEN
          XXX=SQRT(1.0-(S/AAAA)**2)
          TABI(I)=
C(SQRT(ALPHA*XXX)**mmm)*BESSI(mmm,ALPHA*XXX)/B0
CTEMPO=
     &          SQRT(ALPHA*XXX)*BESI1(ALPHA*XXX)/B0
         ELSE
          TABI(I)=0.0
         ENDIF
        ENDDO

#ifdef SP_MP
        LN1=LN+1
#endif
        ALLOCATE (X(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, X',IER)
        ENDIF

       ALLOCATE (W(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, W',IER)
           RETURN
        ENDIF

C       DETERMINE NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)

        ALLOCATE (PROJ(NS,NS,NUMTH), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, PROJ',IER)
           DEALLOCATE (W)
           RETURN
        ENDIF

        ALLOCATE (BI(0:N2,N,NANG), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, BI',IER)
           DEALLOCATE (W,PROJ)
        ENDIF

        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,+1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,-1,IRTFLG)  ! UNUSED??
        IF (IRTFLG .NE. 0) RETURN

C       Read all the projection data
        DO    K=1,NANG
C          PRINT  *,' PROJECTION #',K

C          OPEN DESIRED FILE
           CALL FILGET(FINPAT,FINPIC,NLET,ILIST(K),IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPROJ,'O',IFORM,NSAM,NSAM,NSL,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN


           CALL READV(INPROJ,BI(0,1,K),NS,NS,NS,NS,1)
           CLOSE(INPROJ)
        ENDDO
C  Prepare all the projection data
        NS2=NS*NS
        INV = +1
        DO  IMI=1,NANG,NUMTH
c$omp parallel do private(K,J,I),SHARED(INV)
          DO  K=IMI,MIN(NANG,IMI+NUMTH-1)
           CALL  COP(BI(0,1,K),PROJ(1,1,K-IMI+1),NS2)
           CALL PADD2(PROJ(1,1,K-IMI+1),NS,BI(0,1,K),LSD,N)
           CALL FMRS_2(BI(0,1,K),N,N,INV)
           DO  J=1,N
              DO  I=0,N2
                 BI(I,J,K)=BI(I,J,K)*(-1)**(I+J+1)
              ENDDO
           ENDDO
          ENDDO
        ENDDO
        DEALLOCATE (PROJ)

        NULL=CHAR(0)
        CALL  FILERD(FINPAT,NLET,NULL,
     &          'TEMPLATE FOR OUTPUT VOLUME',IRTFLG)
        CALL RDPRMI(IT,NOT_USED,NOT_USED,'First output volume number')

C  First, build the full structure
c$omp parallel sections
c$omp section
        X=CMPLX(0.0,0.0)
c$omp section
        W=0.0
c$omp end parallel sections
C
        DO    K=1,NANG
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,K))
            ELSE
             DMS=DM(:,:,K)
            ENDIF
#ifdef SP_MP
            DO  JT=1,LN1
C  DOES NOT WORK IN PARALLEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Cc$omp parallel do private(j),shared(N,N2,JT,X,W,BI,DMS)
             DO J=-N2+JT,N2,LN1
              CALL ONELINE(J,N,N2,X,W,BI(0,1,K),DMS)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
              CALL ONELINE(J,N,N2,X,W,BI(0,1,K),DMS)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO
C
C          END OF PROJECTIONS LOOP
        ENDDO

C       DO NOT SYMMETRIZE PLANE 0 OF X, WILL BE NEEDED AS IS FOR PARTIAL VOLUMES
        CALL  SYMPLANE0_WW(W,N2,N)
C
C  Have volume X and weighting W in Fourier space.
C  Sample with replacement per angular direction
C  Array NORD used for resampling has length equal to the maximum number
C   of projections per direction

	INAGMX=MAXVAL(LISTANG)
        ALLOCATE (ERND(INAGMX),NORD(INAGMX), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'BP 3F, NANG',IER)
           RETURN
        ENDIF

        ALLOCATE (XX(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, XX',IER)
           RETURN
        ENDIF

C  The weights should remain the same....            **************************
C Endless loop
        DO NK=1,NITER
C  Create new copy of X.
        XX=X
C Process each angular direction independently
C  NPRJ is the beginning of the projections with the same direction
	NPRJ=0
	DO LNG=1,NUMDIR
C  Current number of projections to process
	LANG=LISTANG(LNG)
	IF(LANG.GT.1)  THEN
         CALL RANDOM_NUMBER(HARVEST=ERND(1:LANG))
         NORD(1:LANG) = MIN0(LANG,MAX0(1,INT(ERND(1:LANG)*LANG+0.5)))
C  Sort the selection of images
        CALL ISORT(NORD,LANG)
C  Go through selected projections, remove from 3D missing,
C                                   add multiple...
C  LK1 and LK2 number projections within the current (i.e., LNG)
C   projection direction. 
C   The actual projection number is KEYNUMANG(LK2+NPRJ)
        LK1=0
        LK2=0
100     LK1=LK1+1
101     LK2=LK2+1
        IF(NORD(LK1).EQ.LK2)  THEN
        MULT=1
102      IF(LK1.EQ.LANG)  THEN
          IF(MULT.GT.1) THEN
C  Add multiple LK2
            MULTI=MULT-1
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,KEYNUMANG(LK2+NPRJ)))
            ELSE
             DMS=DM(:,:,KEYNUMANG(LK2+NPRJ))
            ENDIF
#ifdef SP_MP
            DO  JT=1,LN1
Cc$omp parallel do private(j),shared(N,N2,JT,XX,BI,DMS)
             DO J=-N2+JT,N2,LN1
              CALL ONELINEM_DP(J,N,N2,XX,BI(0,1,KEYNUMANG(LK2+NPRJ)),
     &				DMS,MULTI)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
              CALL ONELINEM_DP(J,N,N2,XX,BI(0,1,KEYNUMANG(LK2+NPRJ)),
     &				DMS,MULTI)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO
C
          ENDIF
C  End
C  Check how many from this series are missing at the end,
C  i.e.,  if the series is 1-12, what if 11 and 12 are missing?
	  IF(NORD(LANG).LT.LANG)  THEN
C  Remove the missing ones at the end
	   DO ID=NORD(LANG)+1,LANG
C Remove ID
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,KEYNUMANG(ID+NPRJ)))
            ELSE
             DMS=DM(:,:,KEYNUMANG(ID+NPRJ))
            ENDIF
#ifdef SP_MP
            DO  JT=1,LN1
Cc$omp parallel do private(j),shared(N,N2,JT,XX,BI,DMS)
             DO J=-N2+JT,N2,LN1
           CALL ONELINEDEL_DP(J,N,N2,XX,BI(0,1,KEYNUMANG(ID+NPRJ)),DMS)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
           CALL ONELINEDEL_DP(J,N,N2,XX,BI(0,1,KEYNUMANG(ID+NPRJ)),DMS)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO
	   
	   ENDDO
	  ENDIF
C  End and out
           GOTO  103
         ELSE
          IF(NORD(LK1).EQ.NORD(LK1+1))  THEN
          LK1=LK1+1
          MULT=MULT+1
          GOTO  102
          ELSE
           IF(MULT.GT.1) THEN
C  Add multiple LK2
           MULTI=MULT-1
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,KEYNUMANG(LK2+NPRJ)))
            ELSE
             DMS=DM(:,:,KEYNUMANG(LK2+NPRJ))
            ENDIF
#ifdef SP_MP
            DO  JT=1,LN1
Cc$omp parallel do private(j),shared(N,N2,JT,XX,BI,DMS)
             DO J=-N2+JT,N2,LN1
              CALL ONELINEM_DP(J,N,N2,XX,BI(0,1,KEYNUMANG(LK2+NPRJ)),
     &				DMS,MULTI)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
              CALL ONELINEM_DP(J,N,N2,XX,BI(0,1,KEYNUMANG(LK2+NPRJ)),
     &				DMS,MULTI)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO
C

           ENDIF
C  Continue
           GOTO 100
          ENDIF
         ENDIF
        ELSE
C Remove LK2
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,KEYNUMANG(LK2+NPRJ)))
            ELSE
             DMS=DM(:,:,KEYNUMANG(LK2+NPRJ))
            ENDIF
#ifdef SP_MP
            DO  JT=1,LN1
Cc$omp parallel do private(j),shared(N,N2,JT,XX,BI,DMS)
             DO J=-N2+JT,N2,LN1
           CALL ONELINEDEL_DP(J,N,N2,XX,BI(0,1,KEYNUMANG(LK2+NPRJ)),DMS)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
           CALL ONELINEDEL_DP(J,N,N2,XX,BI(0,1,KEYNUMANG(LK2+NPRJ)),DMS)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO
C
        GOTO  101
        ENDIF
103     CONTINUE
	ENDIF
	NPRJ=NPRJ+LISTANG(LNG)
C  End of do loop over angular directions....
	ENDDO
C  Convert XX to single precision
C	CALL  CONV_TO_SP(XX,XX,N2,N)
C       SYMMETRIZE PLANE 0
        CALL  SYMPLANE0_X(XX,N2,N)
C  Calculate real space volume
        CALL NRMW2(XX,W,N2,N)
        CALL WINDKB2(XX,XX,NS,LSD,N)

        CALL OUTVOL(XX,NS,NS,NS,NK+IT-1)

        ENDDO
        DEALLOCATE (X, W, XX, BI, NORD, ERND)

C
        END

#endif

#if 1
CPPPPPPP------------------ VAR3DNNWRA ----------------------------------
C  Selection with replacements per angular direction, NN version, uses averages
        SUBROUTINE  VAR3DNNWRA(NS,LSD,N,N2,ILIST,DM,NANG,SM,MAXSYM,
     &		NUMANG,KEYNUMANG,LISTANG,NUMDIR)

        INCLUDE 'CMBLOCK.INC'
        DIMENSION         ILIST(NANG)
        DIMENSION         DM(3,3,NANG),SM(3,3,MAXSYM),DMS(3,3)
	DIMENSION         NUMANG(NANG),KEYNUMANG(NANG),LISTANG(NUMDIR)

        COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: X
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: W
        COMPLEX, DIMENSION(:,:), ALLOCATABLE :: BT
C  2D Fourier transforms of the input data
        COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: BI
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: PROJ
C  Additional matrices
        INTEGER, DIMENSION(:), ALLOCATABLE :: NORD
        REAL, DIMENSION(:), ALLOCATABLE :: ERND

        COMMON  /F_SPEC/  FINPAT,NLET,FINPIC
        CHARACTER*80      FINPIC,FINPAT,FILNAM
        CHARACTER*1  NULL
        DATA  IOPIC/98/,INPROJ/99/

        CALL RDPRMI(NITER,NOT_USED,NOT_USED,
     &          'Number of volumes to be created')
        ALLOCATE (X(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, X',IER)
        ENDIF

       ALLOCATE (W(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, W',IER)
           RETURN
        ENDIF

C It will not have the parallel version
C       DETERMINE NUMBER OF OMP THREADS
C        CALL GETTHREADS(NUMTH)
	NUMTH=1

        ALLOCATE (PROJ(NS,NS,NUMTH), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, PROJ',IER)
           DEALLOCATE (W)
           RETURN
        ENDIF

        ALLOCATE (BI(0:N2,N,NANG), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, BI',IER)
           DEALLOCATE (W,PROJ)
        ENDIF

        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,+1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,-1,IRTFLG)  ! UNUSED??
        IF (IRTFLG .NE. 0) RETURN

C
C  Read all the projection data
        DO    K=1,NANG
C          PRINT  *,' PROJECTION #',K

C          OPEN DESIRED FILE
           CALL FILGET(FINPAT,FINPIC,NLET,ILIST(K),IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPROJ,'O',IFORM,NSAM,NSAM,NSL,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN


           CALL READV(INPROJ,BI(0,1,K),NS,NS,NS,NS,1)
           CLOSE(INPROJ)
        ENDDO
C  Prepare all the projection data
        NS2=NS*NS
        INV = +1
        DO  IMI=1,NANG,NUMTH
c$omp parallel do private(K,J,I),SHARED(INV)
          DO  K=IMI,MIN(NANG,IMI+NUMTH-1)
           CALL  COP(BI(0,1,K),PROJ(1,1,K-IMI+1),NS2)
           CALL PADD2(PROJ(1,1,K-IMI+1),NS,BI(0,1,K),LSD,N)
           CALL FMRS_2(BI(0,1,K),N,N,INV)
           DO  J=1,N
              DO  I=0,N2
                 BI(I,J,K)=BI(I,J,K)*(-1)**(I+J+1)
              ENDDO
           ENDDO
          ENDDO
        ENDDO
        DEALLOCATE (PROJ)

        NULL=CHAR(0)
        CALL  FILERD(FINPAT,NLET,NULL,
     &          'TEMPLATE FOR OUTPUT VOLUME',IRTFLG)
        CALL RDPRMI(IT,NOT_USED,NOT_USED,'First output volume number')

C  First, build the volume of weights
        W=0.0
C
C  NPRJ is the beginning of the projections with the same direction
	NPRJ=0
	DO LNG=1,NUMDIR
C  Current number of projections to process
	LANG=LISTANG(LNG)
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,KEYNUMANG(1+NPRJ)))
            ELSE
             DMS=DM(:,:,KEYNUMANG(1+NPRJ))
            ENDIF
            DO J=-N2+1,N2
              CALL ONELINENW(J,N,N2,W,DMS)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
C
C
	NPRJ=NPRJ+LISTANG(LNG)
C  End of do loop over angular directions....
        ENDDO
C       SYMMETRIZE PLANE 0 FOR WEIGHTS
        CALL  SYMPLANE0_W(W,N2,N)

C  Have weights W in Fourier space.
C  Sample with replacement per angular direction
C  Array NORD used for resampling has length equal to the maximum number
C   of projections per direction

	INAGMX=MAXVAL(LISTANG)
        ALLOCATE (ERND(INAGMX),NORD(INAGMX), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'BP 3F, NANG',IER)
           RETURN
        ENDIF

C BT will contain current sum of projections at the same angular direction
        ALLOCATE (BT(0:N2,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           DEALLOCATE (W)
           CALL ERRT(46,'VA 3, BI',IER)
        ENDIF
C Endless loop
        DO NK=1,NITER
C Reset the new volume X.
C The weights should remain the same....            **************************
        X=CMPLX(0.0,0.0)
C Process each angular direction independently
C  NPRJ is the beginning of the projections with the same direction
	NPRJ=0
	DO LNG=1,NUMDIR
C  Current number of projections to process
	LANG=LISTANG(LNG)
	IF(LANG.GT.1)  THEN
         CALL RANDOM_NUMBER(HARVEST=ERND(1:LANG))
         NORD(1:LANG) = MIN0(LANG,MAX0(1,INT(ERND(1:LANG)*LANG+0.5)))
C  LK2 is the projection number within the current (i.e., LNG)
C   projection direction. 
C   The actual projection number is KEYNUMANG(LK2+NPRJ)
	BT=CMPLX(0.0,0.0)
	DO LK2=1,LANG
	 BT=BT+BI(:,:,KEYNUMANG(NORD(LK2)+NPRJ))
	ENDDO
	ELSE
C  There was just one projection (should not be the case...)
	 BT=BI(:,:,KEYNUMANG(1+NPRJ))
	ENDIF
C  Add BT corresponding to multiple projections at the direction LNG
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,KEYNUMANG(1+NPRJ)))
            ELSE
             DMS=DM(:,:,KEYNUMANG(1+NPRJ))
            ENDIF
            DO J=-N2+1,N2
              CALL ONELINEN(J,N,N2,X,BT,DMS)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
C
	NPRJ=NPRJ+LISTANG(LNG)
C  End of do loop over angular directions....
	ENDDO
C       SYMMETRIZE PLANE 0
        CALL  SYMPLANE0_X(X,N2,N)
C  Calculate real space volume
        CALL NORMN4(X,W,N2,N)
        CALL WINDUM(X,X,NS,LSD,N)

        CALL OUTVOL(X,NS,NS,NS,NK+IT-1)

        ENDDO
        DEALLOCATE (X, W, BT, BI, NORD, ERND)

C
        END

#else
CPPPPPPP------------------ VAR3DNNWRA ----------------------------------
C  Selection with replacements per angular direction, NN version
C  THIS VERSION DOES NOT WORK FOR SOME REASON
        SUBROUTINE  VAR3DNNWRA(NS,LSD,N,N2,ILIST,DM,NANG,SM,MAXSYM,
     &		NUMANG,KEYNUMANG,LISTANG,NUMDIR)

        INCLUDE 'CMBLOCK.INC'
        DIMENSION         ILIST(NANG)
        DIMENSION         DM(3,3,NANG),SM(3,3,MAXSYM),DMS(3,3)
	DIMENSION         NUMANG(NANG),KEYNUMANG(NANG),LISTANG(NUMDIR)

C  2D Fourier transforms of the input data
C        COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: X,XX
C        DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: W,WW
        COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: X,XX
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: W,WW
        COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: BI
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: PROJ
C  Additional matrices
        INTEGER, DIMENSION(:), ALLOCATABLE :: NORD
        REAL, DIMENSION(:), ALLOCATABLE :: ERND

        COMMON  /F_SPEC/  FINPAT,NLET,FINPIC
        CHARACTER*80      FINPIC,FINPAT,FILNAM
        CHARACTER*1  NULL
        DATA  IOPIC/98/,INPROJ/99/

        CALL RDPRMI(NITER,NOT_USED,NOT_USED,
     &          'Number of volumes to be created')
        ALLOCATE (X(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, X',IER)
        ENDIF

       ALLOCATE (W(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, W',IER)
           RETURN
        ENDIF

C       DETERMINE NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)

        ALLOCATE (PROJ(NS,NS,NUMTH), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, PROJ',IER)
           DEALLOCATE (W)
           RETURN
        ENDIF

        ALLOCATE (BI(0:N2,N,NANG), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, BI',IER)
           DEALLOCATE (W,PROJ)
        ENDIF

        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,+1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,-1,IRTFLG)  ! UNUSED??
        IF (IRTFLG .NE. 0) RETURN

C       Read all the projection data
        DO    K=1,NANG
C          PRINT  *,' PROJECTION #',K

C          OPEN DESIRED FILE
           CALL FILGET(FINPAT,FINPIC,NLET,ILIST(K),IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPROJ,'O',IFORM,NSAM,NSAM,NSL,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN


           CALL READV(INPROJ,BI(0,1,K),NS,NS,NS,NS,1)
           CLOSE(INPROJ)
        ENDDO
C  Prepare all the projection data
        NS2=NS*NS
        INV = +1
        DO  IMI=1,NANG,NUMTH
c$omp parallel do private(K,J,I),SHARED(INV)
          DO  K=IMI,MIN(NANG,IMI+NUMTH-1)
           CALL  COP(BI(0,1,K),PROJ(1,1,K-IMI+1),NS2)
           CALL PADD2(PROJ(1,1,K-IMI+1),NS,BI(0,1,K),LSD,N)
           CALL FMRS_2(BI(0,1,K),N,N,INV)
           DO  J=1,N
              DO  I=0,N2
                 BI(I,J,K)=BI(I,J,K)*(-1)**(I+J+1)
              ENDDO
           ENDDO
          ENDDO
        ENDDO
        DEALLOCATE (PROJ)

        NULL=CHAR(0)
        CALL  FILERD(FINPAT,NLET,NULL,
     &          'TEMPLATE FOR OUTPUT VOLUME',IRTFLG)
        CALL RDPRMI(IT,NOT_USED,NOT_USED,'First output volume number')

C  First, build the full structure
c$omp parallel sections
c$omp section
        X=CMPLX(0.0,0.0)
c$omp section
        W=0.0
c$omp end parallel sections
C
        DO    K=1,NANG
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,K))
            ELSE
             DMS=DM(:,:,K)
            ENDIF
c$omp parallel do private(j),shared(N,N2,X,W,BI,DMS)
            DO J=-N2+1,N2
              CALL ONELINENN(J,N,N2,X,W,BI(0,1,K),DMS)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
C
C          END OF PROJECTIONS LOOP
        ENDDO

C       SYMMETRIZE PLANE 0 FOR WEIGHTS
        CALL  SYMPLANE0_W(W,N2,N)

C  Have volume X and weighting W in Fourier space.
C  Sample with replacement per angular direction
C  Array NORD used for resampling has length equal to the maximum number
C   of projections per direction

	INAGMX=MAXVAL(LISTANG)
        ALLOCATE (ERND(INAGMX),NORD(INAGMX), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'BP 3F, NANG',IER)
           RETURN
        ENDIF

        ALLOCATE (XX(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, XX',IER)
           RETURN
        ENDIF

C Endless loop
        DO NK=1,NITER
C  Create new copy of X.
        XX=X
C Process each angular direction independently
C  NPRJ is the beginning of the projections with the same direction
	NPRJ=0
	DO LNG=1,NUMDIR
C  Current number of projections to process
	LANG=LISTANG(LNG)
	IF(LANG.GT.1)  THEN
         CALL RANDOM_NUMBER(HARVEST=ERND(1:LANG))
         NORD(1:LANG) = MIN0(LANG,MAX0(1,INT(ERND(1:LANG)*LANG+0.5)))
C  Sort the selection of images
        CALL ISORT(NORD,LANG)
C  Go through selected projections, remove from 3D missing,
C                                   add multiple...
C  LK1 and LK2 number projections within the current (i.e., LNG)
C   projection direction. 
C   The actual projection number is KEYNUMANG(LK2+NPRJ)
        LK1=0
        LK2=0
100     LK1=LK1+1
101     LK2=LK2+1
        IF(NORD(LK1).EQ.LK2)  THEN
        MULT=1
102      IF(LK1.EQ.LANG)  THEN
          IF(MULT.GT.1) THEN
C  Add multiple LK2
	  MULTI=MULT-1
            DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,KEYNUMANG(LK2+NPRJ)))
            ELSE
             DMS=DM(:,:,KEYNUMANG(LK2+NPRJ))
            ENDIF
c$omp parallel do private(j),shared(N,N2,XX,BI,DMS,MULTI)
            DO J=-N2+1,N2
              CALL ONELINENM(J,N,N2,XX,BI(0,1,KEYNUMANG(LK2+NPRJ)),
     &				DMS,MULTI)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
C
          ENDIF
C  End
C  Check how many from this series are missing at the end,
C  i.e.,  if the series is 1-12, what if 11 and 12 are missing?
	  IF(NORD(LANG).LT.LANG)  THEN
C  Remove the missing ones at the end
	   DO ID=NORD(LANG)+1,LANG
C Remove ID
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,KEYNUMANG(ID+NPRJ)))
            ELSE
             DMS=DM(:,:,KEYNUMANG(ID+NPRJ))
            ENDIF
c$omp parallel do private(j),shared(N,N2,XX,BI,DMS)
            DO J=-N2+1,N2
             CALL ONELINEND(J,N,N2,XX,BI(0,1,KEYNUMANG(ID+NPRJ)),DMS)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
	   
	   ENDDO
	  ENDIF
C  End and out
           GOTO  103
         ELSE
          IF(NORD(LK1).EQ.NORD(LK1+1))  THEN
          LK1=LK1+1
          MULT=MULT+1
          GOTO  102
          ELSE
           IF(MULT.GT.1) THEN
C  Add multiple LK2
	   MULTI=MULT-1
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,KEYNUMANG(LK2+NPRJ)))
            ELSE
             DMS=DM(:,:,KEYNUMANG(LK2+NPRJ))
            ENDIF
c$omp parallel do private(j),shared(N,N2,XX,BI,DMS,MULTI)
            DO J=-N2+1,N2
             CALL ONELINENM(J,N,N2,XX,BI(0,1,KEYNUMANG(LK2+NPRJ)),
     &				DMS,MULTI)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
C

           ENDIF
C  Continue
           GOTO 100
          ENDIF
         ENDIF
        ELSE
C Remove LK2
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,KEYNUMANG(LK2+NPRJ)))
            ELSE
             DMS=DM(:,:,KEYNUMANG(LK2+NPRJ))
            ENDIF
c$omp parallel do private(j),shared(N,N2,XX,BI,DMS)
            DO J=-N2+1,N2
             CALL ONELINEND(J,N,N2,XX,BI(0,1,KEYNUMANG(LK2+NPRJ)),DMS)
            ENDDO
C   END OF SYMMETRIES LOOP
           ENDDO
C
        GOTO  101
        ENDIF
103     CONTINUE
	ENDIF
	NPRJ=NPRJ+LISTANG(LNG)
C  End of do loop over angular directions....
	ENDDO
C  Convert XX and WW to single precision
C	CALL  CONV_TO_SP(XX,WW,XX,WW,N2,N)
C       SYMMETRIZE PLANE 0
        CALL  SYMPLANE0_X(XX,N2,N)
C  Calculate real space volume
        CALL NORMN4(XX,W,N2,N)
        CALL WINDUM(XX,XX,NS,LSD,N)

        CALL OUTVOL(XX,NS,NS,NS,NK+IT-1)

        ENDDO
        DEALLOCATE (X, W, XX, BI, NORD, ERND)

C
        END
#endif
        SUBROUTINE  CONV_TO_SP(X,XX,N2,N)
        COMPLEX  X(0:N2,N*N)
	COMPLEX*16  XX(0:N2,N*N)
	DO I=1,N*N
	DO J=0,N2
	X(J,I)=XX(J,I)
	ENDDO
	ENDDO
	END

        SUBROUTINE  CONV_TO_SP_W(W,WW,N2,N)
        DIMENSION  W(0:N2,N*N)
	DOUBLE PRECISION WW(0:N2,N*N)
	DO I=1,N*N
	DO J=0,N2
	W(J,I)=WW(J,I)
	ENDDO
	ENDDO
	END

        SUBROUTINE  CONV_TO_SP_(X,W,XX,WW,N2,N)
        DIMENSION  W(0:N2,N*N)
        COMPLEX  X(0:N2,N*N)
	DOUBLE PRECISION WW(0:N2,N*N)
	COMPLEX*16  XX(0:N2,N*N)
	DO I=1,N*N
	DO J=0,N2
	W(J,I)=WW(J,I)
	X(J,I)=XX(J,I)
	ENDDO
	ENDDO
	END

	SUBROUTINE  TIEANG(ILIST,NANG,NUMANG,ANGBUF,NUMDIR,IRTFLG)
        DIMENSION  NUMANG(NANG),ILIST(NANG)
        DIMENSION  ANGBUF(4,NANG)
	NUMDIR=1
	NUMANG(1)=NUMDIR
         DO K=2,NANG

C          READ ANGLES FROM THE DOCUMENT FILE.
C          ORDER IN THE DOCUMENT FILE IS PSI, THETA, PHI AND ANGLES 
C          ARE IN DEGREES! 
C     ** PSI IS ASSUMED TO BE ZERO **

           ITMP   = ILIST(K)
	   DO  L=1,K-1
	    IL = ILIST(L)
	    IF(ANGBUF(3,ITMP).EQ.ANGBUF(3,IL) .AND. 
     &         ANGBUF(4,ITMP).EQ.ANGBUF(4,IL)) THEN
	     NUMANG(K)=NUMANG(L)
	     GOTO  12
	    ENDIF
	   ENDDO
	 NUMDIR=NUMDIR+1
	 NUMANG(K)=NUMDIR
12	 CONTINUE
	 ENDDO
	END

CPPPPPPP------------------ VAR3DQJACK ----------------------------------
C  Jackknife method (leave one out)
        SUBROUTINE  VAR3DQJACK(NS,LSD,N,N2,ILIST,DM,NANG,SM,MAXSYM)

        INCLUDE 'CMBLOCK.INC'
        DIMENSION         ILIST(NANG)
        DIMENSION         DM(3,3,NANG),SM(3,3,MAXSYM),DMS(3,3)

C  2D Fourier transforms of the input data
        COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: X,BI,XX
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: W,WW,PROJ
C  Additional matrices
        INTEGER, DIMENSION(:), ALLOCATABLE :: NORD

        COMMON  /F_SPEC/  FINPAT,NLET,FINPIC
        CHARACTER*80      FINPIC,FINPAT,FILNAM
        CHARACTER*1  NULL
        DOUBLE PRECISION  PI
        PARAMETER         (LTAB=4999)
        COMMON  /TABS/    LN2,FLTB,TABI(0:LTAB)
C In this version the order of the Bessel function is mmm=1
        COMMON  /BESSEL_PARAM/  ALPHA,AAAA,NNN
C,mmm
        DATA  IOPIC/98/,INPROJ/99/
        PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
        PARAMETER (TWOPI = 2*QUADPI)

        MONO(K1,K2)=MIN0(K1,K2)+((MAX0(K1,K2)-1)*(MAX0(K1,K2)-2)/2)
C
C K=6
        LN=5
        LN2=LN/2
C Generalized Kaiser-Bessel window according to Lewitt
C M=NS, N=N
        R=NS/2
        V=REAL(LN-1)/2.0/REAL(N)
        ALPHA=6.5
C       AAAA=0.0079
        AAAA=0.9*V
        NNN=3
C       mmm=1
C       GENERATE TABLE WITH INTERPOLANTS
C       B0=(SQRT(ALPHA)**mmm)*BESSI(mmm,ALPHA)
        B0=SQRT(ALPHA)*BESI1(ALPHA)
        FLTB=REAL(LTAB)/REAL(LN2+1)
C  Cannot be parallel as there are DATA satements in BESI1
cc$omp parallel do private(i,s,x),shared(mmm)
        DO  I=0,LTAB
         S=REAL(I)/FLTB/N
         IF(S.LE.AAAA)  THEN
          XXX=SQRT(1.0-(S/AAAA)**2)
          TABI(I)=
C(SQRT(ALPHA*XXX)**mmm)*BESSI(mmm,ALPHA*XXX)/B0
CTEMPO=
     &          SQRT(ALPHA*XXX)*BESI1(ALPHA*XXX)/B0
         ELSE
          TABI(I)=0.0
         ENDIF
        ENDDO

        ALLOCATE (X(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, X',IER)
        ENDIF

       ALLOCATE (W(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, W',IER)
           RETURN
        ENDIF

C       DETERMINE NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)

        ALLOCATE (PROJ(NS,NS,NUMTH), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, PROJ',IER)
           DEALLOCATE (W)
           RETURN
        ENDIF

        ALLOCATE (NORD(NANG), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, NORD',IER)
           RETURN
        ENDIF

        ALLOCATE (BI(0:N2,N,NANG), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, BI',IER)
           DEALLOCATE (W)
           DEALLOCATE (PROJ)
        ENDIF

        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,+1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CALL FMRS_PLAN(.TRUE.,BI,N,N,1, 0,-1,IRTFLG)  ! UNUSED??
        IF (IRTFLG .NE. 0) RETURN

C       Read all the projection data
        DO    K=1,NANG
C          PRINT  *,' PROJECTION #',K

C          OPEN DESIRED FILE
           CALL FILGET(FINPAT,FINPIC,NLET,ILIST(K),IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPROJ,'O',IFORM,NSAM,NSAM,NSL,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN


           CALL READV(INPROJ,BI(0,1,K),NS,NS,NS,NS,1)
           CLOSE(INPROJ)
        ENDDO
C  Prepare all the projection data
        NS2=NS*NS
        INV = +1
        DO  IMI=1,NANG,NUMTH
c$omp parallel do private(K,J,I),SHARED(INV)
          DO  K=IMI,MIN(NANG,IMI+NUMTH-1)
           CALL  COP(BI(0,1,K),PROJ(1,1,K-IMI+1),NS2)
           CALL PADD2(PROJ(1,1,K-IMI+1),NS,BI(0,1,K),LSD,N)
           CALL FMRS_2(BI(0,1,K),N,N,INV)
           DO  J=1,N
              DO  I=0,N2
                 BI(I,J,K)=BI(I,J,K)*(-1)**(I+J+1)
              ENDDO
           ENDDO
          ENDDO
        ENDDO
        DEALLOCATE (PROJ)

        NULL=CHAR(0)
        CALL  FILERD(FINPAT,NLET,NULL,
     &          'TEMPLATE FOR OUTPUT VOLUME',IRTFLG)

#ifdef SP_MP
        LN1=LN+1
#endif
C **********************************************************************
        ALLOCATE (XX(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, XX',IER)
           RETURN
        ENDIF

        ALLOCATE (WW(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN
           CALL ERRT(46,'VA 3, WW',IER)
           RETURN
        ENDIF

C  First, build the full structure
c$omp parallel sections
c$omp section
        X=CMPLX(0.0,0.0)
c$omp section
        W=0.0
c$omp end parallel sections
C
        DO    K=1,NANG
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,K))
            ELSE
             DMS=DM(:,:,K)
            ENDIF
#ifdef SP_MP
            DO  JT=1,LN1
c$omp parallel do private(j),shared(N,N2,JT,X,W,BI,DMS)
             DO J=-N2+JT,N2,LN1
              CALL ONELINE(J,N,N2,X,W,BI(0,1,K),DMS)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
              CALL ONELINE(J,N,N2,X,W,BI(0,1,K),DMS)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO
C
C          END OF PROJECTIONS LOOP
        ENDDO

C       DO NOT SYMMETRIZE PLANE 0, WILL BE NEEDED AS IS FOR PARTIAL VOLUMES
CCCCCCCC        CALL  SYMPLANE0(X,W,N2,N)

C  Have volume X and weighting W in Fourier space.
C  Sample without replacement
        DO NK=1,NANG
C  Create new copies of X and W.
c$omp parallel sections
c$omp section
        XX=X
c$omp section
        WW=W
c$omp end parallel sections
C Sample MREM without replacements
C  REMOVE ONE PROJECTION
C Remove NK
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,NK))
            ELSE
             DMS=DM(:,:,NK)
            ENDIF
#ifdef SP_MP
            DO  JT=1,LN1
c$omp parallel do private(j),shared(N,N2,JT,XX,WW,BI,DMS)
             DO J=-N2+JT,N2,LN1
              CALL ONELINEDEL(J,N,N2,XX,WW,BI(0,1,NK),DMS)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
              CALL ONELINEDEL(J,N,N2,XX,WW,BI(0,1,NK),DMS)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO

C       SYMMETRIZE PLANE 0
        CALL  SYMPLANE0(XX,WW,N2,N)
C  Calculate real space volume
        CALL NRMW2(XX,WW,N2,N)
        CALL WINDKB2(XX,XX,NS,LSD,N)
C
        CALL OUTVOL(XX,NS,NS,NS,NK)
        ENDDO
        DEALLOCATE (X, W, WW, XX, BI, NORD)
C
        END

      SUBROUTINE ISORT(A,N)
      INTEGER A(N),T,TT
C
C     SINGLETON SORT PROGRAM TO ORDER B AND C USING A AS A KEY
C     AS OF THE PRESENT TIME (FEB. 1971) THIS IS THE FASTEST GENERAL
C     PURPOSE SORTING METHOD KNOWN.
C MODIFIED VERSION WITH REAL KEY ARRAY. J.FRANK, FEB. 1977
C
      INTEGER IL(16), IU(16)

C
      M = 1
      I = 1
      J = N
    5 IF (I .GE. J) GO TO 70
C
C     ORDER THE TWO ENDS AND THE MIDDLE
C
   10 K = I
      IJ = (I + J)/2
      T = A(IJ)
      IF (A(I) .LE. T) GO TO 20
      A(IJ) = A(I)
      A(I) = T
      T = A(IJ)
   20 L = J
      IF (A(J) .GE. T) GO TO 40
      IF (A(J) .LT. A(I)) GO TO 25
      A(IJ) = A(J)
      A(J) = T
      T = A(IJ)
      GO TO 40
   25 A(IJ) = A(I)
      A(I) = A(J)
      A(J) = T
      T = A(IJ)
      GO TO 40
C
C     SPLIT THE SEQUENCE BETWEEN I AND J INTO TWO SEQUENCES.  THAT
C     SEQUENCE BETWEEN I AND L WILL CONTAIN ONLY ELEMENTS LESS THAN OR
C     EQUAL TO T, WHILE THAT BETWEEN K AND J WILL CONTAIN ONLY ELEMENTS
C     GREATER THAN T.
C
   30 A(L) = A(K)
      A(K) = TT
   40 L = L - 1
      IF (A(L) .GT. T) GO TO 40
      TT = A(L)
   50 K = K + 1
      IF (A(K) .LT. T) GO TO 50
      IF (K .LE. L) GO TO 30
C
C     SAVE THE END POINTS OF THE LONGER SEQUENCE IN IL AND IU, AND SORT
C     THE SHORTER SEQUENCE.
C
      IF (L - I .LE. J - K) GO TO 60
      IL(M) = I
      IU(M) = L
      I = K
      M = M + 1
      GO TO 80
   60 IL(M) = K
      IU(M) = J
      J = L
      M = M + 1
      GO TO 80
C
C     RETRIEVE END POINTS PREVIOUSLY SAVED AND SORT BETWEEN THEM.
C
   70 M = M - 1
      IF (M .EQ. 0) RETURN
      I = IL(M)
      J = IU(M)
C
C     IF THE SEQUENCE IS LONGER THAN 11 OR IS THE FIRST SEQUENCE, SORT
C     BY SPLITTING RECURSIVELY.
C
   80 IF (J - I .GE. 11) GO TO 10
      IF (I .EQ. 1) GO TO 5
C
C     IF THE SEQUENCE IS 11 OR LESS LONG, SORT IT BY A SHELLSORT.
C
      I = I - 1
   90 I = I + 1
      IF (I .EQ. J) GO TO 70
      T = A(I + 1)
      IF (A(I) .LE. T) GO TO 90
      K = I
  100 A(K+1) = A(K)
      K = K - 1
      IF (T .LT. A(K)) GO TO 100
      A(K+1) = T
      GO TO 90
C
      END
C++*******************************************************************
C
C
C  PARAMETERS:
C       A    KEY ARRAY USED IN SORTING
C       B    ARRAY TO BE SORTED
C       N
C
C--*******************************************************************


      SUBROUTINE ISORT2 ( A, B, N)

      INTEGER  A(N),T,TT
      INTEGER  IL(16), IU(16)
      INTEGER  B(N),X

      M = 1
      I = 1
      J = N
    5 IF (I .GE. J) GO TO 70

C    ORDER THE TWO ENDS AND THE MIDDLE

   10 K = I
      IJ = (I + J)/2
      T = A(IJ)
      IF (A(I) .LE. T) GO TO 20
      A(IJ) = A(I)
      A(I) = T
      T = A(IJ)
      X = B(I)
      B(I) = B(IJ)
      B(IJ) = X
   20 L = J
      IF (A(J) .GE. T) GO TO 40
      IF (A(J) .LT. A(I)) GO TO 25
      A(IJ) = A(J)
      A(J) = T
      T = A(IJ)
      X = B(IJ)
      B(IJ) = B(J)
      B(J) = X
      GO TO 40
   25 A(IJ) = A(I)
      A(I) = A(J)
      A(J) = T
      T = A(IJ)
      X = B(J)
      B(J) = B(IJ)
      B(IJ) = B(I)
      B(I) = X
      GO TO 40

C     SPLIT THE SEQUENCE BETWEEN I AND J INTO TWO SEQUENCES.  THAT
C     SEQUENCE BETWEEN I AND L WILL CONTAIN ONLY ELEMENTS LESS THAN OR
C     EQUAL TO T, WHILE THAT BETWEEN K AND J WILL CONTAIN ONLY ELEMENTS
C     GREATER THAN T.

   30 A(L) = A(K)
      A(K) = TT
      X = B(L)
      B(L) = B(K)
      B(K) = X
   40 L = L - 1
      IF (A(L) .GT. T) GO TO 40
      TT = A(L)
   50 K = K + 1
      IF (A(K) .LT. T) GO TO 50
      IF (K .LE. L) GO TO 30

C     SAVE THE END POINTS OF THE LONGER SEQUENCE IN IL AND IU, AND SORT
C     THE SHORTER SEQUENCE.

      IF (L - I .LE. J - K) GO TO 60
      IL(M) = I
      IU(M) = L
      I = K
      M = M + 1
      GO TO 80
   60 IL(M) = K
      IU(M) = J
      J = L
      M = M + 1
      GO TO 80

C     RETRIEVE END POINTS PREVIOUSLY SAVED AND SORT BETWEEN THEM.

   70 M = M - 1
      IF (M .EQ. 0) RETURN
      I = IL(M)
      J = IU(M)

C     IF THE SEQUENCE IS LONGER THAN 11 OR IS THE FIRST SEQUENCE, SORT
C     BY SPLITTING RECURSIVELY.

   80 IF (J - I .GE. 11) GO TO 10
      IF (I .EQ. 1) GO TO 5

C    IF THE SEQUENCE IS 11 OR LESS LONG, SORT IT BY A SHELLSORT.

      I = I - 1
   90 I = I + 1
      IF (I .EQ. J) GO TO 70
      T = A(I + 1)
      IF (A(I) .LE. T) GO TO 90
      X = B(I+1)
      K = I
  100 A(K+1) = A(K)
      B(K+1) = B(K)
      K = K - 1
      IF (T .LT. A(K)) GO TO 100
      A(K+1) = T
      B(K+1) = X
      GO TO 90

      END
         SUBROUTINE  OUTVOL(X,NSAM,NROW,NSLICE,IT)

         INCLUDE 'CMBLOCK.INC'

        COMMON  /F_SPEC/  FINPAT,NLET,FINPIC
        CHARACTER*80      FINPIC,FINPAT,FILNAM
         DIMENSION     X(NSAM,NROW,NSLICE)
         DATA INPIC/69/

         CALL  FILGET(FINPAT,FILNAM,NLET,IT,IRTFLG)
         IF(IRTFLG .NE. 0) RETURN

         MAXIM  = 0
         IFORM  = 3
         CALL OPFILEC(0,.FALSE.,FILNAM,INPIC,'U',IFORM,NSAM,NROW,NSLICE,
     &               MAXIM,' ',.FALSE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         CALL WRITEV(INPIC,X,NSAM,NROW,NSAM,NROW,NSLICE)

         CLOSE(INPIC)
         END
C       --------------------- ONELINEDEL_DP ---------------------------------

        SUBROUTINE  ONELINEDEL_DP(J,N,N2,X,BI,DM)

        COMPLEX        BI(0:N2,N)
	COMPLEX     X(0:N2,N,N),BTQ
        DIMENSION      DM(6)
        PARAMETER      (LTAB=4999)
        COMMON  /TABS/ LN2,FLTB,TABI(0:LTAB)

        IF (J .GE. 0)  THEN
           JP=J+1
        ELSE
           JP=N+J+1
        ENDIF

        DO  I=0,N2
           IF ((I*I+J*J.LT.N*N/4).AND..NOT.(I.EQ.0.AND.J.LT.0))  THEN
              XNEW=I*DM(1)+J*DM(4)
              YNEW=I*DM(2)+J*DM(5)
              ZNEW=I*DM(3)+J*DM(6)
              IF (XNEW.LT.0.0)  THEN
                 XNEW=-XNEW
                 YNEW=-YNEW
                 ZNEW=-ZNEW
                 BTQ=CONJG(BI(I,JP))
              ELSE
                 BTQ=BI(I,JP)
              ENDIF
              IXN=IFIX(XNEW+0.5+N)-N
              IYN=IFIX(YNEW+0.5+N)-N
              IZN=IFIX(ZNEW+0.5+N)-N
              IF (IXN.LE.N2-LN2-1 .AND.
     &        IYN.GE.-N2+2+LN2.AND.IYN.LE.N2-LN2-1 .AND.
     &        IZN.GE.-N2+2+LN2.AND.IZN.LE.N2-LN2-1) THEN
                 IF (IXN.GE.0) THEN
C                   MAKE SURE THAT LOWER LIMIT FOR X DOES NOT GO BELOW 0
                    LB=-MIN0(IXN,LN2)
                    DO LZ=-LN2,LN2
                       IZP=IZN+LZ
                       IF(IZP.GE.0) THEN
                          IZA=IZP+1
                       ELSE
                          IZA=N+IZP+1
                       ENDIF
                       TZ = TABI(NINT(ABS(ZNEW-IZP)*FLTB))
           	       IF(TZ.NE.0.0)  THEN
                       DO  LY=-LN2,LN2
                          IYP=IYN+LY
                          IF (IYP.GE.0) THEN
                             IYA=IYP+1
                          ELSE
                             IYA=N+IYP+1
                          ENDIF
                          TY = TABI(NINT(ABS(YNEW-IYP)*FLTB))*TZ
           		  IF(TY.NE.0.0)  THEN
                          DO  IXP=LB+IXN,LN2+IXN
C                            GET THE WEIGHT
C                            WG=FINF(XNEW,YNEW,ZNEW,IXT,IYT,IZT)
                             WG=TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY
			     IF(WG.NE.0.0)  THEN
                              X(IXP,IYA,IZA)=X(IXP,IYA,IZA)-BTQ*WG
			     ENDIF
                          ENDDO
			  ENDIF
                       ENDDO
		       ENDIF
                   ENDDO
                ENDIF

C               ADD REFLECTED POINTS
                IF (IXN .LT. LN2) THEN
                   DO  LZ=-LN2,LN2
                      IZP=IZN+LZ
                      IF (IZP.GT.0)  THEN
                         IZT=N-IZP+1
                      ELSE
                         IZT=-IZP+1
                      ENDIF
                      TZ=TABI(NINT(ABS(ZNEW-IZP)*FLTB))
C
           	      IF(TZ.NE.0.0)  THEN
                      DO  LY=-LN2,LN2
                         IYP=IYN+LY
                         IF (IYP.GT.0) THEN
                            IYT=N-IYP+1
                         ELSE
                            IYT=-IYP+1
                         ENDIF
                         TY=TABI(NINT(ABS(YNEW-IYP)*FLTB))*TZ
			 IF(TY.NE.0.0)  THEN
                         DO  IXP=-LN2+IXN,-1
C                           GET THE WEIGHT
C                           WG=FINF(XNEW,YNEW,ZNEW,IXT,IYT,IZT)
                            WG=TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY
			    IF(WG.NE.0.0)  THEN
                          X(-IXP,IYT,IZT)=X(-IXP,IYT,IZT)-CONJG(BTQ)*WG
			    ENDIF
                         ENDDO
			 ENDIF
                      ENDDO
		      ENDIF
                   ENDDO
                ENDIF
              ENDIF
           ENDIF
C          END J-I LOOP
        ENDDO

        END
C       --------------------- ONELINEM_DP ---------------------------------

        SUBROUTINE  ONELINEM_DP(J,N,N2,X,BI,DM,MULT)

        COMPLEX        BI(0:N2,N)
	COMPLEX     X(0:N2,N,N),BTQ
        DIMENSION      DM(6)
        PARAMETER      (LTAB=4999)
        COMMON  /TABS/ LN2,FLTB,TABI(0:LTAB)

        IF (J .GE. 0)  THEN
           JP=J+1
        ELSE
           JP=N+J+1
        ENDIF

        DO  I=0,N2
           IF ((I*I+J*J.LT.N*N/4).AND..NOT.(I.EQ.0.AND.J.LT.0))  THEN
              XNEW=I*DM(1)+J*DM(4)
              YNEW=I*DM(2)+J*DM(5)
              ZNEW=I*DM(3)+J*DM(6)
              IF (XNEW.LT.0.0)  THEN
                 XNEW=-XNEW
                 YNEW=-YNEW
                 ZNEW=-ZNEW
                 BTQ=CONJG(BI(I,JP))
              ELSE
                 BTQ=BI(I,JP)
              ENDIF
              IXN=IFIX(XNEW+0.5+N)-N
              IYN=IFIX(YNEW+0.5+N)-N
              IZN=IFIX(ZNEW+0.5+N)-N
              IF (IXN.LE.N2-LN2-1 .AND.
     &        IYN.GE.-N2+2+LN2.AND.IYN.LE.N2-LN2-1 .AND.
     &        IZN.GE.-N2+2+LN2.AND.IZN.LE.N2-LN2-1) THEN
                 IF (IXN.GE.0) THEN
C                   MAKE SURE THAT LOWER LIMIT FOR X DOES NOT GO BELOW 0
                    LB=-MIN0(IXN,LN2)
                    DO LZ=-LN2,LN2
                       IZP=IZN+LZ
                       IF(IZP.GE.0) THEN
                          IZA=IZP+1
                       ELSE
                          IZA=N+IZP+1
                       ENDIF
                       TZ = TABI(NINT(ABS(ZNEW-IZP)*FLTB))
           	       IF(TZ.NE.0.0)  THEN
                       DO  LY=-LN2,LN2
                          IYP=IYN+LY
                          IF (IYP.GE.0) THEN
                             IYA=IYP+1
                          ELSE
                             IYA=N+IYP+1
                          ENDIF
                          TY = TABI(NINT(ABS(YNEW-IYP)*FLTB))*TZ
           		  IF(TY.NE.0.0)  THEN
                          DO  IXP=LB+IXN,LN2+IXN
C                            GET THE WEIGHT
C                            WG=FINF(XNEW,YNEW,ZNEW,IXT,IYT,IZT)
                             WG=TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY*MULT
			     IF(WG.NE.0.0)  THEN
                              X(IXP,IYA,IZA)=X(IXP,IYA,IZA)+BTQ*WG
			     ENDIF
                          ENDDO
			  ENDIF
                       ENDDO
		       ENDIF
                   ENDDO
                ENDIF

C               ADD REFLECTED POINTS
                IF (IXN .LT. LN2) THEN
                   DO  LZ=-LN2,LN2
                      IZP=IZN+LZ
                      IF (IZP.GT.0)  THEN
                         IZT=N-IZP+1
                      ELSE
                         IZT=-IZP+1
                      ENDIF
                      TZ=TABI(NINT(ABS(ZNEW-IZP)*FLTB))
C
           	      IF(TZ.NE.0.0)  THEN
                      DO  LY=-LN2,LN2
                         IYP=IYN+LY
                         IF (IYP.GT.0) THEN
                            IYT=N-IYP+1
                         ELSE
                            IYT=-IYP+1
                         ENDIF
                         TY=TABI(NINT(ABS(YNEW-IYP)*FLTB))*TZ
			 IF(TY.NE.0.0)  THEN
                         DO  IXP=IXN-LN2,-1
C                           GET THE WEIGHT
C                           WG=FINF(XNEW,YNEW,ZNEW,IXT,IYT,IZT)
                            WG=TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY*MULT
			    IF(WG.NE.0.0)  THEN
                          X(-IXP,IYT,IZT)=X(-IXP,IYT,IZT)+CONJG(BTQ)*WG
			    ENDIF
                         ENDDO
			 ENDIF
                      ENDDO
		      ENDIF
                   ENDDO
                ENDIF
              ENDIF
           ENDIF
C          END J-I LOOP
        ENDDO

        END
C***********************************************************************

        SUBROUTINE  ONELINE_DP(J,N,N2,X,W,BI,DM)

        DOUBLE PRECISION      W(0:N2,N,N)
        COMPLEX        BI(0:N2,N)
	COMPLEX*16     X(0:N2,N,N),BTQ
        DIMENSION      DM(6)

        PARAMETER      (LTAB=4999)
        COMMON  /TABS/ LN2,FLTB,TABI(0:LTAB)

        IF (J .GE. 0)  THEN
           JP = J+1
        ELSE
           JP = N+J+1
        ENDIF

        DO  I=0,N2
           IF (((I*I+J*J) .LT.  (N*N/4)) .AND..NOT. 
     &          (I.EQ. 0  .AND. J.LT.0)) THEN
              XNEW = I * DM(1) + J * DM(4)
              YNEW = I * DM(2) + J * DM(5)
              ZNEW = I * DM(3) + J * DM(6)

              IF (XNEW .LT. 0.0)  THEN
                 XNEW = -XNEW
                 YNEW = -YNEW
                 ZNEW = -ZNEW
                 BTQ  = CONJG(BI(I,JP))
              ELSE
                 BTQ  = BI(I,JP)
              ENDIF

              IXN = IFIX(XNEW+0.5+N) - N
              IYN = IFIX(YNEW+0.5+N) - N
              IZN = IFIX(ZNEW+0.5+N) - N

              IF (IXN .LE. (N2-LN2-1)  .AND.
     &            IYN .GE. (-N2+2+LN2) .AND. IYN .LE. (N2-LN2-1) .AND.
     &            IZN .GE. (-N2+2+LN2) .AND. IZN .LE. (N2-LN2-1)) THEN

                 IF (IXN .GE. 0) THEN
C                   MAKE SURE THAT LOWER LIMIT FOR X DOES NOT GO BELOW 0
                    LB = -MIN0(IXN,LN2)
                    DO LZ=-LN2,LN2
                       IZP = IZN + LZ
                       IF(IZP .GE. 0) THEN
                          IZA = IZP + 1
                       ELSE
                          IZA = N + IZP + 1
                       ENDIF
 
                       TZ  = TABI(NINT(ABS(ZNEW-IZP) * FLTB))

                       IF (TZ .NE. 0.0)  THEN
                          DO  LY=-LN2,LN2
                             IYP = IYN + LY
                             IF (IYP .GE .0) THEN
                                IYA = IYP + 1
                             ELSE
                                IYA = N + IYP + 1
                             ENDIF
 
                             TY  = TABI(NINT(ABS(YNEW-IYP) * FLTB)) * TZ
                             IF (TY .NE. 0.0)  THEN
                                DO  IXP=LB+IXN,LN2+IXN

C                                  GET THE WEIGHT
                                   WG=TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY
                                   IF (WG .NE. 0.0) THEN

                                      X(IXP,IYA,IZA) =
     &                                    X(IXP,IYA,IZA) + BTQ * WG
                                      W(IXP,IYA,IZA) = 
     &                                    W(IXP,IYA,IZA) + WG
                                   ENDIF
                                ENDDO
                             ENDIF
                          ENDDO
                       ENDIF
                   ENDDO
                ENDIF

C               ADD REFLECTED POINTS
                IF (IXN .LT. LN2) THEN
                   DO  LZ=-LN2,LN2
                      IZP = IZN + LZ
                      IZT =  - IZP + 1
                      IF (IZP .GT. 0)  IZT = N + IZT

                      TZ = TABI(NINT(ABS(ZNEW-IZP) * FLTB))

                      IF (TZ .NE. 0.0)  THEN
                         DO  LY=-LN2,LN2
                            IYP = IYN + LY
                            IYT = -IYP + 1
                            IF (IYP .GT. 0) IYT = IYT + N

                            TY = TABI(NINT(ABS(YNEW-IYP) * FLTB)) * TZ
                            IF (TY .NE. 0.0)  THEN
                               DO  IXP=IXN-LN2,-1

C                                 GET THE WEIGHT
                                  WG = TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY

                                  IF (WG .NE. 0.0)  THEN
                                     X(-IXP,IYT,IZT) = 
     &                                   X(-IXP,IYT,IZT) + CONJG(BTQ)*WG
                                     W(-IXP,IYT,IZT) =
     &                                   W(-IXP,IYT,IZT) + WG
                                  ENDIF
                               ENDDO
                            ENDIF
                         ENDDO
                      ENDIF
                   ENDDO
                ENDIF
              ENDIF
           ENDIF
C          END J-I LOOP
        ENDDO

        END

C       --------------------- ONELINENW ---------------------------------

        SUBROUTINE  ONELINENW(J,N,N2,NR,DM)

        DIMENSION      NR(0:N2,N,N)
        DIMENSION      DM(6)

        IF (J .GE. 0)  THEN
           JP=J+1
        ELSE
           JP=N+J+1
        ENDIF

        DO  I=0,N2
           IF ((I*I+J*J.LT.N*N/4).AND..NOT.(I.EQ.0.AND.J.LT.0))  THEN
              XNEW=I*DM(1)+J*DM(4)
              YNEW=I*DM(2)+J*DM(5)
              ZNEW=I*DM(3)+J*DM(6)
              IF (XNEW.LT.0.0)  THEN
                 XNEW=-XNEW
                 YNEW=-YNEW
                 ZNEW=-ZNEW
              ENDIF
              IXN=IFIX(XNEW+0.5+N)-N
              IYN=IFIX(YNEW+0.5+N)-N
              IZN=IFIX(ZNEW+0.5+N)-N
              IF (IXN.LE.N2 .AND.
     &        IYN.GE.-N2.AND.IYN.LE.N2 .AND.
     &        IZN.GE.-N2.AND.IZN.LE.N2) THEN
               IF(IXN.GE.0) THEN
                       IF(IZN.GE.0) THEN
                          IZA=IZN+1
                       ELSE
                          IZA=N+IZN+1
                       ENDIF
                       IF(IYN.GE.0) THEN
                          IYA=IYN+1
                       ELSE
                          IYA=N+IYN+1
                       ENDIF
	               NR(IXN,IYA,IZA)=NR(IXN,IYA,IZA)+1
C              
               ELSE
                      IF (IZN.GT.0)  THEN
                         IZT=N-IZN+1
                      ELSE
                         IZT=-IZN+1
                      ENDIF
                      IF(IYN.GT.0) THEN
                         IYT=N-IYN+1
                      ELSE
                         IYT=-IYN+1
                      ENDIF
                      NR(-IXN,IYT,IZT)=NR(-IXN,IYT,IZT)+1
                ENDIF
              ENDIF
           ENDIF
C          END J-I LOOP
        ENDDO

        END
C       --------------------- ONELINEND ---------------------------------

        SUBROUTINE  ONELINEND(J,N,N2,X,BI,DM)

        COMPLEX        BI(0:N2,N),X(0:N2,N,N),BTQ
        DIMENSION      DM(6)

        IF (J .GE. 0)  THEN
           JP=J+1
        ELSE
           JP=N+J+1
        ENDIF

        DO  I=0,N2
           IF ((I*I+J*J.LT.N*N/4).AND..NOT.(I.EQ.0.AND.J.LT.0))  THEN
              XNEW=I*DM(1)+J*DM(4)
              YNEW=I*DM(2)+J*DM(5)
              ZNEW=I*DM(3)+J*DM(6)
              IF (XNEW.LT.0.0)  THEN
                 XNEW=-XNEW
                 YNEW=-YNEW
                 ZNEW=-ZNEW
                 BTQ=CONJG(BI(I,JP))
              ELSE
                 BTQ=BI(I,JP)
              ENDIF
              IXN=IFIX(XNEW+0.5+N)-N
              IYN=IFIX(YNEW+0.5+N)-N
              IZN=IFIX(ZNEW+0.5+N)-N
              IF (IXN.LE.N2 .AND.
     &        IYN.GE.-N2.AND.IYN.LE.N2 .AND.
     &        IZN.GE.-N2.AND.IZN.LE.N2) THEN
               IF(IXN.GE.0) THEN
                       IF(IZN.GE.0) THEN
                          IZA=IZN+1
                       ELSE
                          IZA=N+IZN+1
                       ENDIF
                       IF(IYN.GE.0) THEN
                          IYA=IYN+1
                       ELSE
                          IYA=N+IYN+1
                       ENDIF
                 X(IXN,IYA,IZA)=X(IXN,IYA,IZA)-BTQ
C              
               ELSE
                      IF (IZN.GT.0)  THEN
                         IZT=N-IZN+1
                      ELSE
                         IZT=-IZN+1
                      ENDIF
                      IF(IYN.GT.0) THEN
                         IYT=N-IYN+1
                      ELSE
                         IYT=-IYN+1
                      ENDIF
                      X(-IXN,IYT,IZT)=X(-IXN,IYT,IZT)-CONJG(BTQ)
                ENDIF
              ENDIF
           ENDIF
C          END J-I LOOP
        ENDDO

        END
CC       --------------------- ONELINENWD ---------------------------------

        SUBROUTINE  ONELINENWD(J,N,N2,X,NR,BI,DM)

        COMPLEX        BI(0:N2,N),X(0:N2,N,N),BTQ
        DIMENSION      NR(0:N2,N,N)
        DIMENSION      DM(6)

        IF (J .GE. 0)  THEN
           JP=J+1
        ELSE
           JP=N+J+1
        ENDIF

        DO  I=0,N2
           IF ((I*I+J*J.LT.N*N/4).AND..NOT.(I.EQ.0.AND.J.LT.0))  THEN
              XNEW=I*DM(1)+J*DM(4)
              YNEW=I*DM(2)+J*DM(5)
              ZNEW=I*DM(3)+J*DM(6)
              IF (XNEW.LT.0.0)  THEN
                 XNEW=-XNEW
                 YNEW=-YNEW
                 ZNEW=-ZNEW
                 BTQ=CONJG(BI(I,JP))
              ELSE
                 BTQ=BI(I,JP)
              ENDIF
              IXN=IFIX(XNEW+0.5+N)-N
              IYN=IFIX(YNEW+0.5+N)-N
              IZN=IFIX(ZNEW+0.5+N)-N
              IF (IXN.LE.N2 .AND.
     &        IYN.GE.-N2.AND.IYN.LE.N2 .AND.
     &        IZN.GE.-N2.AND.IZN.LE.N2) THEN
               IF(IXN.GE.0) THEN
                       IF(IZN.GE.0) THEN
                          IZA=IZN+1
                       ELSE
                          IZA=N+IZN+1
                       ENDIF
                       IF(IYN.GE.0) THEN
                          IYA=IYN+1
                       ELSE
                          IYA=N+IYN+1
                       ENDIF
                       X(IXN,IYA,IZA)=X(IXN,IYA,IZA)-BTQ
	               NR(IXN,IYA,IZA)=NR(IXN,IYA,IZA)-1
C              
               ELSE
                      IF (IZN.GT.0)  THEN
                         IZT=N-IZN+1
                      ELSE
                         IZT=-IZN+1
                      ENDIF
                      IF(IYN.GT.0) THEN
                         IYT=N-IYN+1
                      ELSE
                         IYT=-IYN+1
                      ENDIF
                      X(-IXN,IYT,IZT)=X(-IXN,IYT,IZT)-CONJG(BTQ)
                      NR(-IXN,IYT,IZT)=NR(-IXN,IYT,IZT)-1
                ENDIF
              ENDIF
           ENDIF
C          END J-I LOOP
        ENDDO

        END
C       --------------------- ONELINENM ---------------------------------

        SUBROUTINE  ONELINENM(J,N,N2,X,BI,DM,MULT)

        COMPLEX        BI(0:N2,N),X(0:N2,N,N),BTQ
        DIMENSION      DM(6)

        IF (J .GE. 0)  THEN
           JP=J+1
        ELSE
           JP=N+J+1
        ENDIF

        DO  I=0,N2
           IF ((I*I+J*J.LT.N*N/4).AND..NOT.(I.EQ.0.AND.J.LT.0))  THEN
              XNEW=I*DM(1)+J*DM(4)
              YNEW=I*DM(2)+J*DM(5)
              ZNEW=I*DM(3)+J*DM(6)
              IF (XNEW.LT.0.0)  THEN
                 XNEW=-XNEW
                 YNEW=-YNEW
                 ZNEW=-ZNEW
                 BTQ=CONJG(BI(I,JP))
              ELSE
                 BTQ=BI(I,JP)
              ENDIF
              IXN=IFIX(XNEW+0.5+N)-N
              IYN=IFIX(YNEW+0.5+N)-N
              IZN=IFIX(ZNEW+0.5+N)-N
              IF (IXN.LE.N2 .AND.
     &        IYN.GE.-N2.AND.IYN.LE.N2 .AND.
     &        IZN.GE.-N2.AND.IZN.LE.N2) THEN
               IF(IXN.GE.0) THEN
                       IF(IZN.GE.0) THEN
                          IZA=IZN+1
                       ELSE
                          IZA=N+IZN+1
                       ENDIF
                       IF(IYN.GE.0) THEN
                          IYA=IYN+1
                       ELSE
                          IYA=N+IYN+1
                       ENDIF
                 X(IXN,IYA,IZA)=X(IXN,IYA,IZA)+BTQ*MULT
C              
               ELSE
                      IF (IZN.GT.0)  THEN
                         IZT=N-IZN+1
                      ELSE
                         IZT=-IZN+1
                      ENDIF
                      IF(IYN.GT.0) THEN
                         IYT=N-IYN+1
                      ELSE
                         IYT=-IYN+1
                      ENDIF
                      X(-IXN,IYT,IZT)=X(-IXN,IYT,IZT)+CONJG(BTQ)*MULT
                ENDIF
              ENDIF
           ENDIF
C          END J-I LOOP
        ENDDO

        END
C       --------------------- ONELINENWM ---------------------------------

        SUBROUTINE  ONELINENWM(J,N,N2,X,NR,BI,DM,MULT)

        COMPLEX        BI(0:N2,N),X(0:N2,N,N),BTQ
        DIMENSION      NR(0:N2,N,N)
        DIMENSION      DM(6)

        IF (J .GE. 0)  THEN
           JP=J+1
        ELSE
           JP=N+J+1
        ENDIF

        DO  I=0,N2
           IF ((I*I+J*J.LT.N*N/4).AND..NOT.(I.EQ.0.AND.J.LT.0))  THEN
              XNEW=I*DM(1)+J*DM(4)
              YNEW=I*DM(2)+J*DM(5)
              ZNEW=I*DM(3)+J*DM(6)
              IF (XNEW.LT.0.0)  THEN
                 XNEW=-XNEW
                 YNEW=-YNEW
                 ZNEW=-ZNEW
                 BTQ=CONJG(BI(I,JP))
              ELSE
                 BTQ=BI(I,JP)
              ENDIF
              IXN=IFIX(XNEW+0.5+N)-N
              IYN=IFIX(YNEW+0.5+N)-N
              IZN=IFIX(ZNEW+0.5+N)-N
              IF (IXN.LE.N2 .AND.
     &        IYN.GE.-N2.AND.IYN.LE.N2 .AND.
     &        IZN.GE.-N2.AND.IZN.LE.N2) THEN
               IF(IXN.GE.0) THEN
                       IF(IZN.GE.0) THEN
                          IZA=IZN+1
                       ELSE
                          IZA=N+IZN+1
                       ENDIF
                       IF(IYN.GE.0) THEN
                          IYA=IYN+1
                       ELSE
                          IYA=N+IYN+1
                       ENDIF
                       X(IXN,IYA,IZA)=X(IXN,IYA,IZA)+BTQ*MULT
	               NR(IXN,IYA,IZA)=NR(IXN,IYA,IZA)+MULT
C              
               ELSE
                      IF (IZN.GT.0)  THEN
                         IZT=N-IZN+1
                      ELSE
                         IZT=-IZN+1
                      ENDIF
                      IF(IYN.GT.0) THEN
                         IYT=N-IYN+1
                      ELSE
                         IYT=-IYN+1
                      ENDIF
                      X(-IXN,IYT,IZT)=X(-IXN,IYT,IZT)+CONJG(BTQ)*MULT
                      NR(-IXN,IYT,IZT)=NR(-IXN,IYT,IZT)+MULT
                ENDIF
              ENDIF
           ENDIF
C          END J-I LOOP
        ENDDO

        END
C       --------------------- ONELINEN ---------------------------------

        SUBROUTINE  ONELINEN(J,N,N2,X,BI,DM)

        COMPLEX        BI(0:N2,N),X(0:N2,N,N),BTQ
        DIMENSION      DM(6)

        IF (J .GE. 0)  THEN
           JP=J+1
        ELSE
           JP=N+J+1
        ENDIF

        DO  I=0,N2
           IF ((I*I+J*J.LT.N*N/4).AND..NOT.(I.EQ.0.AND.J.LT.0))  THEN
              XNEW=I*DM(1)+J*DM(4)
              YNEW=I*DM(2)+J*DM(5)
              ZNEW=I*DM(3)+J*DM(6)
              IF (XNEW.LT.0.0)  THEN
                 XNEW=-XNEW
                 YNEW=-YNEW
                 ZNEW=-ZNEW
                 BTQ=CONJG(BI(I,JP))
              ELSE
                 BTQ=BI(I,JP)
              ENDIF
              IXN=IFIX(XNEW+0.5+N)-N
              IYN=IFIX(YNEW+0.5+N)-N
              IZN=IFIX(ZNEW+0.5+N)-N
              IF (IXN.LE.N2 .AND.
     &        IYN.GE.-N2.AND.IYN.LE.N2 .AND.
     &        IZN.GE.-N2.AND.IZN.LE.N2) THEN
               IF(IXN.GE.0) THEN
                       IF(IZN.GE.0) THEN
                          IZA=IZN+1
                       ELSE
                          IZA=N+IZN+1
                       ENDIF
                       IF(IYN.GE.0) THEN
                          IYA=IYN+1
                       ELSE
                          IYA=N+IYN+1
                       ENDIF
                 X(IXN,IYA,IZA)=X(IXN,IYA,IZA)+BTQ
C              
               ELSE
                      IF (IZN.GT.0)  THEN
                         IZT=N-IZN+1
                      ELSE
                         IZT=-IZN+1
                      ENDIF
                      IF(IYN.GT.0) THEN
                         IYT=N-IYN+1
                      ELSE
                         IYT=-IYN+1
                      ENDIF
                      X(-IXN,IYT,IZT)=X(-IXN,IYT,IZT)+CONJG(BTQ)
                ENDIF
              ENDIF
           ENDIF
C          END J-I LOOP
        ENDDO

        END
C       ----------------SYMPLANE0_W ---------------------------------------
 
        SUBROUTINE  SYMPLANE0_W(W,N2,N)

        INTEGER  W(0:N2,N,N)

C       SYMMETRIZE PLANE 0
        DO  IZA=2,N2
           DO  IYA=2,N2
              W(0,IYA,IZA)=W(0,IYA,IZA)+W(0,N-IYA+2,N-IZA+2)
              W(0,N-IYA+2,N-IZA+2)=W(0,IYA,IZA)
              W(0,N-IYA+2,IZA)=W(0,N-IYA+2,IZA)+W(0,IYA,N-IZA+2)
              W(0,IYA,N-IZA+2)=W(0,N-IYA+2,IZA)
           ENDDO
        ENDDO

        DO  IYA=2,N2
           W(0,IYA,1)=W(0,IYA,1)+W(0,N-IYA+2,1)
           W(0,N-IYA+2,1)=W(0,IYA,1)
        ENDDO

        DO  IZA=2,N2
           W(0,1,IZA)=W(0,1,IZA)+W(0,1,N-IZA+2)
           W(0,1,N-IZA+2)=W(0,1,IZA)
        ENDDO

        END
C       ----------------SYMPLANE0_WW ---------------------------------------
 
        SUBROUTINE  SYMPLANE0_WW(W,N2,N)

        DIMENSION  W(0:N2,N,N)

C       SYMMETRIZE PLANE 0
        DO  IZA=2,N2
           DO  IYA=2,N2
              W(0,IYA,IZA)=W(0,IYA,IZA)+W(0,N-IYA+2,N-IZA+2)
              W(0,N-IYA+2,N-IZA+2)=W(0,IYA,IZA)
              W(0,N-IYA+2,IZA)=W(0,N-IYA+2,IZA)+W(0,IYA,N-IZA+2)
              W(0,IYA,N-IZA+2)=W(0,N-IYA+2,IZA)
           ENDDO
        ENDDO

        DO  IYA=2,N2
           W(0,IYA,1)=W(0,IYA,1)+W(0,N-IYA+2,1)
           W(0,N-IYA+2,1)=W(0,IYA,1)
        ENDDO

        DO  IZA=2,N2
           W(0,1,IZA)=W(0,1,IZA)+W(0,1,N-IZA+2)
           W(0,1,N-IZA+2)=W(0,1,IZA)
        ENDDO

        END
C       ----------------SYMPLANE0_X ---------------------------------------
 
        SUBROUTINE  SYMPLANE0_X(X,N2,N)

        COMPLEX  X(0:N2,N,N)

C       SYMMETRIZE PLANE 0
        DO  IZA=2,N2
           DO  IYA=2,N2
              X(0,IYA,IZA)=X(0,IYA,IZA)+CONJG(X(0,N-IYA+2,N-IZA+2))
              X(0,N-IYA+2,N-IZA+2)=CONJG(X(0,IYA,IZA))
              X(0,N-IYA+2,IZA)=X(0,N-IYA+2,IZA)+CONJG(X(0,IYA,N-IZA+2))
              X(0,IYA,N-IZA+2)=CONJG(X(0,N-IYA+2,IZA))
           ENDDO
        ENDDO

        DO  IYA=2,N2
           X(0,IYA,1)=X(0,IYA,1)+CONJG(X(0,N-IYA+2,1))
           X(0,N-IYA+2,1)=CONJG(X(0,IYA,1))
        ENDDO

        DO  IZA=2,N2
           X(0,1,IZA)=X(0,1,IZA)+CONJG(X(0,1,N-IZA+2))
           X(0,1,N-IZA+2)=CONJG(X(0,1,IZA))
        ENDDO

        END
C       --------------------- ONELINEM_W ---------------------------------

        SUBROUTINE  ONELINEM_W(J,N,N2,W,DM,MULT)

        DIMENSION      W(0:N2,N,N)
        DIMENSION      DM(6)
        PARAMETER      (LTAB=4999)
        COMMON  /TABS/ LN2,FLTB,TABI(0:LTAB)

        IF (J .GE. 0)  THEN
           JP=J+1
        ELSE
           JP=N+J+1
        ENDIF

        DO  I=0,N2
           IF ((I*I+J*J.LT.N*N/4).AND..NOT.(I.EQ.0.AND.J.LT.0))  THEN
              XNEW=I*DM(1)+J*DM(4)
              YNEW=I*DM(2)+J*DM(5)
              ZNEW=I*DM(3)+J*DM(6)
              IF (XNEW.LT.0.0)  THEN
                 XNEW=-XNEW
                 YNEW=-YNEW
                 ZNEW=-ZNEW
              ENDIF
              IXN=IFIX(XNEW+0.5+N)-N
              IYN=IFIX(YNEW+0.5+N)-N
              IZN=IFIX(ZNEW+0.5+N)-N
              IF (IXN.LE.N2-LN2-1 .AND.
     &        IYN.GE.-N2+2+LN2.AND.IYN.LE.N2-LN2-1 .AND.
     &        IZN.GE.-N2+2+LN2.AND.IZN.LE.N2-LN2-1) THEN
                 IF (IXN.GE.0) THEN
C                   MAKE SURE THAT LOWER LIMIT FOR X DOES NOT GO BELOW 0
                    LB=-MIN0(IXN,LN2)
                    DO LZ=-LN2,LN2
                       IZP=IZN+LZ
                       IF(IZP.GE.0) THEN
                          IZA=IZP+1
                       ELSE
                          IZA=N+IZP+1
                       ENDIF
                       TZ = TABI(NINT(ABS(ZNEW-IZP)*FLTB))
           	       IF(TZ.NE.0.0)  THEN
                       DO  LY=-LN2,LN2
                          IYP=IYN+LY
                          IF (IYP.GE.0) THEN
                             IYA=IYP+1
                          ELSE
                             IYA=N+IYP+1
                          ENDIF
                          TY = TABI(NINT(ABS(YNEW-IYP)*FLTB))*TZ
           		  IF(TY.NE.0.0)  THEN
                          DO  IXP=LB+IXN,LN2+IXN
C                            GET THE WEIGHT
C                            WG=FINF(XNEW,YNEW,ZNEW,IXT,IYT,IZT)
                             WG=TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY*MULT
			     IF(WG.NE.0.0)  THEN
                              W(IXP,IYA,IZA)=W(IXP,IYA,IZA)+WG
			     ENDIF
                          ENDDO
			  ENDIF
                       ENDDO
		       ENDIF
                   ENDDO
                ENDIF

C               ADD REFLECTED POINTS
                IF (IXN .LT. LN2) THEN
                   DO  LZ=-LN2,LN2
                      IZP=IZN+LZ
                      IF (IZP.GT.0)  THEN
                         IZT=N-IZP+1
                      ELSE
                         IZT=-IZP+1
                      ENDIF
                      TZ=TABI(NINT(ABS(ZNEW-IZP)*FLTB))
C
           	      IF(TZ.NE.0.0)  THEN
                      DO  LY=-LN2,LN2
                         IYP=IYN+LY
                         IF (IYP.GT.0) THEN
                            IYT=N-IYP+1
                         ELSE
                            IYT=-IYP+1
                         ENDIF
                         TY=TABI(NINT(ABS(YNEW-IYP)*FLTB))*TZ
			 IF(TY.NE.0.0)  THEN
                         DO  IXP=IXN-LN2,-1
C                           GET THE WEIGHT
C                           WG=FINF(XNEW,YNEW,ZNEW,IXT,IYT,IZT)
                            WG=TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY*MULT
			    IF(WG.NE.0.0)  THEN
                            W(-IXP,IYT,IZT)=W(-IXP,IYT,IZT)+WG
			    ENDIF
                         ENDDO
			 ENDIF
                      ENDDO
		      ENDIF
                   ENDDO
                ENDIF
              ENDIF
           ENDIF
C          END J-I LOOP
        ENDDO

        END
C       --------------------- ONELINE_X ---------------------------------

        SUBROUTINE  ONELINE_X(J,N,N2,X,BI,DM)

        COMPLEX        BI(0:N2,N)
	COMPLEX     X(0:N2,N,N),BTQ
        DIMENSION      DM(6)
        PARAMETER      (LTAB=4999)
        COMMON  /TABS/ LN2,FLTB,TABI(0:LTAB)

        IF (J .GE. 0)  THEN
           JP=J+1
        ELSE
           JP=N+J+1
        ENDIF

        DO  I=0,N2
           IF ((I*I+J*J.LT.N*N/4).AND..NOT.(I.EQ.0.AND.J.LT.0))  THEN
              XNEW=I*DM(1)+J*DM(4)
              YNEW=I*DM(2)+J*DM(5)
              ZNEW=I*DM(3)+J*DM(6)
              IF (XNEW.LT.0.0)  THEN
                 XNEW=-XNEW
                 YNEW=-YNEW
                 ZNEW=-ZNEW
                 BTQ=CONJG(BI(I,JP))
              ELSE
                 BTQ=BI(I,JP)
              ENDIF
              IXN=IFIX(XNEW+0.5+N)-N
              IYN=IFIX(YNEW+0.5+N)-N
              IZN=IFIX(ZNEW+0.5+N)-N
              IF (IXN.LE.N2-LN2-1 .AND.
     &        IYN.GE.-N2+2+LN2.AND.IYN.LE.N2-LN2-1 .AND.
     &        IZN.GE.-N2+2+LN2.AND.IZN.LE.N2-LN2-1) THEN
                 IF (IXN.GE.0) THEN
C                   MAKE SURE THAT LOWER LIMIT FOR X DOES NOT GO BELOW 0
                    LB=-MIN0(IXN,LN2)
                    DO LZ=-LN2,LN2
                       IZP=IZN+LZ
                       IF(IZP.GE.0) THEN
                          IZA=IZP+1
                       ELSE
                          IZA=N+IZP+1
                       ENDIF
                       TZ = TABI(NINT(ABS(ZNEW-IZP)*FLTB))
           	       IF(TZ.NE.0.0)  THEN
                       DO  LY=-LN2,LN2
                          IYP=IYN+LY
                          IF (IYP.GE.0) THEN
                             IYA=IYP+1
                          ELSE
                             IYA=N+IYP+1
                          ENDIF
                          TY = TABI(NINT(ABS(YNEW-IYP)*FLTB))*TZ
           		  IF(TY.NE.0.0)  THEN
                          DO  IXP=LB+IXN,LN2+IXN
C                            GET THE WEIGHT
C                            WG=FINF(XNEW,YNEW,ZNEW,IXT,IYT,IZT)
                             WG=TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY
			     IF(WG.NE.0.0)  THEN
                              X(IXP,IYA,IZA)=X(IXP,IYA,IZA)+BTQ*WG
			     ENDIF
                          ENDDO
			  ENDIF
                       ENDDO
		       ENDIF
                   ENDDO
                ENDIF

C               ADD REFLECTED POINTS
                IF (IXN .LT. LN2) THEN
                   DO  LZ=-LN2,LN2
                      IZP=IZN+LZ
                      IF (IZP.GT.0)  THEN
                         IZT=N-IZP+1
                      ELSE
                         IZT=-IZP+1
                      ENDIF
                      TZ=TABI(NINT(ABS(ZNEW-IZP)*FLTB))
C
           	      IF(TZ.NE.0.0)  THEN
                      DO  LY=-LN2,LN2
                         IYP=IYN+LY
                         IF (IYP.GT.0) THEN
                            IYT=N-IYP+1
                         ELSE
                            IYT=-IYP+1
                         ENDIF
                         TY=TABI(NINT(ABS(YNEW-IYP)*FLTB))*TZ
			 IF(TY.NE.0.0)  THEN
                         DO  IXP=IXN-LN2,-1
C                           GET THE WEIGHT
C                           WG=FINF(XNEW,YNEW,ZNEW,IXT,IYT,IZT)
                            WG=TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY
			    IF(WG.NE.0.0)  THEN
                          X(-IXP,IYT,IZT)=X(-IXP,IYT,IZT)+CONJG(BTQ)*WG
			    ENDIF
                         ENDDO
			 ENDIF
                      ENDDO
		      ENDIF
                   ENDDO
                ENDIF
              ENDIF
           ENDIF
C          END J-I LOOP
        ENDDO

        END




C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
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
C       --------------------- ONELINEM ---------------------------------

        SUBROUTINE  ONELINEM(J,N,N2,X,W,BI,DM,MULT)

        DIMENSION      W(0:N2,N,N)
        COMPLEX        BI(0:N2,N),X(0:N2,N,N),BTQ
        DIMENSION      DM(6)
        PARAMETER      (LTAB=4999)
        COMMON  /TABS/ LN2,FLTB,TABI(0:LTAB)

        IF (J .GE. 0)  THEN
           JP=J+1
        ELSE
           JP=N+J+1
        ENDIF

        DO  I=0,N2
           IF ((I*I+J*J.LT.N*N/4).AND..NOT.(I.EQ.0.AND.J.LT.0))  THEN
              XNEW=I*DM(1)+J*DM(4)
              YNEW=I*DM(2)+J*DM(5)
              ZNEW=I*DM(3)+J*DM(6)
              IF (XNEW.LT.0.0)  THEN
                 XNEW=-XNEW
                 YNEW=-YNEW
                 ZNEW=-ZNEW
                 BTQ=CONJG(BI(I,JP))
              ELSE
                 BTQ=BI(I,JP)
              ENDIF
              IXN=IFIX(XNEW+0.5+N)-N
              IYN=IFIX(YNEW+0.5+N)-N
              IZN=IFIX(ZNEW+0.5+N)-N
              IF (IXN.LE.N2-LN2-1 .AND.
     &        IYN.GE.-N2+2+LN2.AND.IYN.LE.N2-LN2-1 .AND.
     &        IZN.GE.-N2+2+LN2.AND.IZN.LE.N2-LN2-1) THEN
                 IF (IXN.GE.0) THEN
C                   MAKE SURE THAT LOWER LIMIT FOR X DOES NOT GO BELOW 0
                    LB=-MIN0(IXN,LN2)
                    DO LZ=-LN2,LN2
                       IZP=IZN+LZ
                       IF(IZP.GE.0) THEN
                          IZA=IZP+1
                       ELSE
                          IZA=N+IZP+1
                       ENDIF
                       TZ = TABI(NINT(ABS(ZNEW-IZP)*FLTB))
           	       IF(TZ.NE.0.0)  THEN
                       DO  LY=-LN2,LN2
                          IYP=IYN+LY
                          IF (IYP.GE.0) THEN
                             IYA=IYP+1
                          ELSE
                             IYA=N+IYP+1
                          ENDIF
                          TY = TABI(NINT(ABS(YNEW-IYP)*FLTB))*TZ
           		  IF(TY.NE.0.0)  THEN
                          DO  IXP=LB+IXN,LN2+IXN
C                            GET THE WEIGHT
C                            WG=FINF(XNEW,YNEW,ZNEW,IXT,IYT,IZT)
                             WG=TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY*MULT
			     IF(WG.NE.0.0)  THEN
                              X(IXP,IYA,IZA)=X(IXP,IYA,IZA)+BTQ*WG
                              W(IXP,IYA,IZA)=W(IXP,IYA,IZA)+WG
			     ENDIF
                          ENDDO
			  ENDIF
                       ENDDO
		       ENDIF
                   ENDDO
                ENDIF

C               ADD REFLECTED POINTS
                IF (IXN .LT. LN2) THEN
                   DO  LZ=-LN2,LN2
                      IZP=IZN+LZ
                      IF (IZP.GT.0)  THEN
                         IZT=N-IZP+1
                      ELSE
                         IZT=-IZP+1
                      ENDIF
                      TZ=TABI(NINT(ABS(ZNEW-IZP)*FLTB))
C
           	      IF(TZ.NE.0.0)  THEN
                      DO  LY=-LN2,LN2
                         IYP=IYN+LY
                         IF (IYP.GT.0) THEN
                            IYT=N-IYP+1
                         ELSE
                            IYT=-IYP+1
                         ENDIF
                         TY=TABI(NINT(ABS(YNEW-IYP)*FLTB))*TZ
			 IF(TY.NE.0.0)  THEN
                         DO  IXP=IXN-LN2,-1
C                           GET THE WEIGHT
C                           WG=FINF(XNEW,YNEW,ZNEW,IXT,IYT,IZT)
                            WG=TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY*MULT
			    IF(WG.NE.0.0)  THEN
                          X(-IXP,IYT,IZT)=X(-IXP,IYT,IZT)+CONJG(BTQ)*WG
                            W(-IXP,IYT,IZT)=W(-IXP,IYT,IZT)+WG
			    ENDIF
                         ENDDO
			 ENDIF
                      ENDDO
		      ENDIF
                   ENDDO
                ENDIF
              ENDIF
           ENDIF
C          END J-I LOOP
        ENDDO

        END



C **********************************************************************
C
C  ONELINEDEL
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
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
C ONELINEDEL(J,N,N2,X,W,BI,DM)
C
C NOTE: THIS WAS A FOOLISH MODIFICATION OF: ONELINE.F IN THAT IT ONLY
C       DIFFERS BY SIGNS IN THE LINES THAT APPLY THE WEIGHT (WG) TO
C       X.  ANOTHER SLOPPY ADDITION THAT MAKES MAINTENANCE OF SPIDER
C       DIFFICULT. al
C
C CALLED BY: VAR3D, VAR3DISK
C
C **********************************************************************

        SUBROUTINE  ONELINEDEL(J,N,N2,X,W,BI,DM)

        DIMENSION      W(0:N2,N,N)
        COMPLEX        BI(0:N2,N),X(0:N2,N,N),BTQ
        DIMENSION      DM(6)
        PARAMETER      (LTAB=4999)
        COMMON  /TABS/ LN2,FLTB,TABI(0:LTAB)

        IF (J .GE. 0)  THEN
           JP=J+1
        ELSE
           JP=N+J+1
        ENDIF

        DO  I=0,N2
           IF ((I*I+J*J.LT.N*N/4).AND..NOT.(I.EQ.0.AND.J.LT.0))  THEN
              XNEW=I*DM(1)+J*DM(4)
              YNEW=I*DM(2)+J*DM(5)
              ZNEW=I*DM(3)+J*DM(6)
              IF (XNEW.LT.0.0)  THEN
                 XNEW=-XNEW
                 YNEW=-YNEW
                 ZNEW=-ZNEW
                 BTQ=CONJG(BI(I,JP))
              ELSE
                 BTQ=BI(I,JP)
              ENDIF
              IXN=IFIX(XNEW+0.5+N)-N
              IYN=IFIX(YNEW+0.5+N)-N
              IZN=IFIX(ZNEW+0.5+N)-N
              IF (IXN.LE.N2-LN2-1 .AND.
     &        IYN.GE.-N2+2+LN2.AND.IYN.LE.N2-LN2-1 .AND.
     &        IZN.GE.-N2+2+LN2.AND.IZN.LE.N2-LN2-1) THEN
                 IF (IXN.GE.0) THEN
C                   MAKE SURE THAT LOWER LIMIT FOR X DOES NOT GO BELOW 0
                    LB=-MIN0(IXN,LN2)
                    DO LZ=-LN2,LN2
                       IZP=IZN+LZ
                       IF(IZP.GE.0) THEN
                          IZA=IZP+1
                       ELSE
                          IZA=N+IZP+1
                       ENDIF
                       TZ = TABI(NINT(ABS(ZNEW-IZP)*FLTB))
           	       IF(TZ.NE.0.0)  THEN
                       DO  LY=-LN2,LN2
                          IYP=IYN+LY
                          IF (IYP.GE.0) THEN
                             IYA=IYP+1
                          ELSE
                             IYA=N+IYP+1
                          ENDIF
                          TY = TABI(NINT(ABS(YNEW-IYP)*FLTB))*TZ
           		  IF(TY.NE.0.0)  THEN
                          DO  IXP=LB+IXN,LN2+IXN
C                            GET THE WEIGHT
C                            WG=FINF(XNEW,YNEW,ZNEW,IXT,IYT,IZT)
                             WG=TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY
			     IF(WG.NE.0.0)  THEN
                              X(IXP,IYA,IZA)=X(IXP,IYA,IZA)-BTQ*WG
                              W(IXP,IYA,IZA)=W(IXP,IYA,IZA)-WG
			     ENDIF
                          ENDDO
			  ENDIF
                       ENDDO
		       ENDIF
                   ENDDO
                ENDIF

C               ADD REFLECTED POINTS
                IF (IXN .LT. LN2) THEN
                   DO  LZ=-LN2,LN2
                      IZP=IZN+LZ
                      IF (IZP.GT.0)  THEN
                         IZT=N-IZP+1
                      ELSE
                         IZT=-IZP+1
                      ENDIF
                      TZ=TABI(NINT(ABS(ZNEW-IZP)*FLTB))
C
           	      IF(TZ.NE.0.0)  THEN
                      DO  LY=-LN2,LN2
                         IYP=IYN+LY
                         IF (IYP.GT.0) THEN
                            IYT=N-IYP+1
                         ELSE
                            IYT=-IYP+1
                         ENDIF
                         TY=TABI(NINT(ABS(YNEW-IYP)*FLTB))*TZ
			 IF(TY.NE.0.0)  THEN
                         DO  IXP=-LN2+IXN,-1
C                           GET THE WEIGHT
C                           WG=FINF(XNEW,YNEW,ZNEW,IXT,IYT,IZT)
                            WG=TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY
			    IF(WG.NE.0.0)  THEN
                          X(-IXP,IYT,IZT)=X(-IXP,IYT,IZT)-CONJG(BTQ)*WG
                            W(-IXP,IYT,IZT)=W(-IXP,IYT,IZT)-WG
			    ENDIF
                         ENDDO
			 ENDIF
                      ENDDO
		      ENDIF
                   ENDDO
                ENDIF
              ENDIF
           ENDIF
C          END J-I LOOP
        ENDDO

        END



C       ------------------- WINDKB2 -------------------------------

        SUBROUTINE WINDKB2(BI,R,L,LSD,N)

        DIMENSION  R(L,L,L),BI(LSD,N,N)
        COMMON  /BESSEL_PARAM/  ALPHA,AAAA,NNN

        PARAMETER (QUADPI = 3.14159265358979323846)
        PARAMETER (TWOPI = 2*QUADPI)

        IP = (N-L)/2+MOD(L,2)
        DO K=1,L
           DO J=1,L
              DO I=1,L
                 R(I,J,K) = BI(IP+I,IP+J,IP+K)
              ENDDO
           ENDDO
        ENDDO

        L2  = (L/2)**2
        L2P = (L/2-1)**2
        IP  = L/2+1
        XNU = REAL(NNN)/2.
        
        RI = RIBSL(ALPHA,XNU)
C       IF (ABS(RI-RIN).GT.1.E-5)  PRINT  *,'BESSIK'

        WKB0 = ALPHA**XNU/RI
        QRT  = (TWOPI*AAAA)**2
        TNR  = 0.0
        M    = 0
        DO K=1,L
           DO J=1,L
              DO I=1,L
                 LR = (K-IP)**2+(J-IP)**2+(I-IP)**2
                 IF (LR<=L2) THEN
                 SIGMA=QRT*LR-ALPHA*ALPHA
                  IF (ABS(SIGMA).LT.1.0E-7)  THEN
                     WKB=1.0
                  ELSEIF(SIGMA.GT.0.0)  THEN
C                    2PI A R > ALPHA
                     ART = SQRT(SIGMA)
                     RI = RJBSL(ART, XNU)
C       if(abs(ri-rin)/rin.gt.1.e-5)  print  *,'bessjy',i,j,k
                     WKB=WKB0*RI/ART**XNU
                  ELSE
C                    2PI A R < ALPHA
                     ART = SQRT(ABS(SIGMA))
                     RI = RIBSL(ART,XNU)
C       if(abs(ri-rin)/rin.gt.1.e-5)  print  *,'bessik',i,j,k,ri,rin
                     WKB=WKB0*RI/ART**XNU
                  ENDIF
                  R(I,J,K) = R(I,J,K)/ABS(WKB)
                  IF (LR .GE. L2P .AND. LR <= L2) THEN
                      TNR = TNR+R(I,J,K)
                      M   = M+1
                  ENDIF
                 ENDIF
              ENDDO
           ENDDO
        ENDDO

        TNR = TNR/REAL(M)
c$omp   parallel do private(i,j,k,lr)
        DO K=1,L
           DO J=1,L
              DO I=1,L
                 LR=(K-IP)**2+(J-IP)**2+(I-IP)**2
                 IF (LR<=L2) THEN
                    R(I,J,K)=R(I,J,K)-TNR
                 ELSE
                    R(I,J,K)=0.0
                 ENDIF
              ENDDO
           ENDDO
        ENDDO

        END


