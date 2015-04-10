
C **********************************************************************
C  SSNR3       BESSEL FUNCTION REPLACED
C              OPFILEC                             FEB  03 ARDEAN LEITH
C              BUILDM PARAMETERS                   SEP  03 ARDEAN LEITH
C                  MAXNAM                         JUL  14 ARDEAN LEITH *
C
C **********************************************************************
C *  SSNR3
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
C SSNR3                        
C                        
C **********************************************************************

        SUBROUTINE SSNR3

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'F90ALLOC.INC'

        !COMMON /F_SPEC/ FINPAT,NLET,FINPIC
        !CHARACTER*80    FINPIC,FINPAT,FILNAM,ANGDOC

        CHARACTER(LEN=MAXNAM) :: FINPIC,FINPAT,FILNAM,ANGDOC
        COMMON /F_SPEC/          FINPAT,NLET,FINPIC

        INTEGER, ALLOCATABLE, DIMENSION(:) :: ILIST
        REAL, DIMENSION(:,:), ALLOCATABLE :: DM,SM

C       DOC FILE POINTERS
        REAL, DIMENSION(:,:), POINTER :: ANGBUF, ANGSYM

        DATA  IOPIC/98/,INPIC/99/

        NILMAX = NIMAX

        ALLOCATE(ILIST(NILMAX), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'SSNR3; ILIST',IER)
           RETURN
        ENDIF

        CALL FILELIST(.TRUE.,INPIC,FINPAT,NLET,ILIST,NILMAX,NANG,
     &                 'TEMPLATE FOR 2-D IMAGES',IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9997

        MAXNUM = MAXVAL(ILIST(1:NANG))
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
        IF (IRTFLG .NE. 0) GOTO 9997

C       RETRIEVE ARRAY WITH SYMMETRIES DATA IN IT
        MAXXS  = 0
        MAXSYM = 0
        CALL GETDOCDAT('SYMMETRIES DOC',.TRUE.,ANGDOC,77,.TRUE.,MAXXS,
     &                  MAXSYM,ANGSYM,IRTFLG)
        IF (IRTFLG .NE. 0)  MAXSYM=1

C       OPEN FIRST IMAGE FILE TO DETERMINE NSAM, NROW, NSL
        CALL FILGET(FINPAT,FINPIC,NLET,ILIST(1),INTFLG)

        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,NSAM,NROW,NSL,
     &               MAXIM,'DUMMY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CLOSE(INPIC)

C       WINDOW IF ODD DIMENSION
        N2     = NSAM-MOD(NSAM,2)
        LSD    = N2+2-MOD(N2,2)
        NMAT   = LSD*N2*N2

        ALLOCATE(DM(9,NANG), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'BP 3F, DM',9*NANG)
           GOTO 9997
        ENDIF

        CALL BUILDM(ILIST,DM,NANG,ANGBUF(1,1),.FALSE.,SSDUM,
     &              .FALSE.,IRTFLG)
        DEALLOCATE(ANGBUF)
        IF (IRTFLG .NE. 0) GOTO 9997

        IF (MAXSYM .GT. 1)  THEN
           ALLOCATE(SM(9,MAXSYM), STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN 
              CALL ERRT(46,'BP 3F, SM..',N2*N2*N2+9*MAXSYM)
              GOTO 9997
           ENDIF

           CALL  BUILDS(SM,MAXSYM,ANGSYM(1,1),IRTFLG)
           DEALLOCATE(ANGSYM)

        ELSE
           ALLOCATE(SM(1,1), STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN 
              CALL ERRT(46,'BP 3F, SM..',N2*N2*N2)
              GOTO 9997
           ENDIF
        ENDIF

        CALL WIW3DP(NSAM,LSD,N2,N2/2,ILIST,DM,NANG,SM,MAXSYM)

9997    IF (ALLOCATED(ILIST)) DEALLOCATE (ILIST)        
        IF (ALLOCATED(DM))    DEALLOCATE(DM)
        IF (ALLOCATED(SM))    DEALLOCATE(SM)

        END



C       ------------------ WIW3DP ----------------------------------

        SUBROUTINE  WIW3DP(NS,LSD,N,N2,ILIST,DM,NANG,SM,MAXSYM)

        
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        DIMENSION         ILIST(NANG)
        DIMENSION         DM(3,3,NANG),SM(3,3,MAXSYM),DMS(3,3)

        COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: X,SSNR
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: X2
        INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: NR
        REAL, DIMENSION(:,:), ALLOCATABLE :: PROJ
        COMPLEX, DIMENSION(:,:), ALLOCATABLE :: BI
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: W,W2

        !CHARACTER*80      FINPIC,FINPAT
        !COMMON  /F_SPEC/  FINPAT,NLET,FINPIC

        CHARACTER(LEN=MAXNAM) :: FINPIC,FINPAT
        COMMON /F_SPEC/          FINPAT,NLET,FINPIC

        DOUBLE PRECISION  PI
        PARAMETER         (LTAB=5000)
        COMMON  /TABS/    LN2,FLTB,TABI(0:LTAB)
	COMMON  /BESSEL_PARAM/  ALPHA,AAAA,NNN
C,mmm
        DATA  INPROJ/99/
	PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
	PARAMETER (TWOPI = 2*QUADPI)
C K=6
        LN=5
        LN2=LN/2
C Generalized Kaiser-Bessel window according to Lewitt
C M=NS, N=N
	R=NS/2
	V=REAL(LN-1)/2.0/REAL(N)
	ALPHA=6.5
	AAAA=1.1*V
	NNN=3
C	mmm=1
C       GENERATE TABLE WITH INTERPOLANTS
	B0=SQRT(ALPHA)*BESI1(ALPHA)
        FLTB=REAL(LTAB)/REAL(LN2+1)
cc$omp parallel do private(i,s,x),shared(mmm)
        DO  I=0,LTAB
	 S=REAL(I)/FLTB/N
	 IF(S.LE.AAAA)  THEN
	  XX=SQRT(1.0-(S/AAAA)**2)
	  TABI(I)=SQRT(ALPHA*XX)*BESI1(ALPHA*XX)/B0
	 ELSE
	  TABI(I)=0.0
	 ENDIF
        ENDDO
c	do i=ltab,0,-1
c	 if(tabi(i).ne.0.0)  then
c	  print  *,i
c	  goto 111
c	 endif
c	enddo
c111	continue
        ALLOCATE (X(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'BP 3F, X2',IER)
           RETURN
        ENDIF

        ALLOCATE (X2(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'BP 3F, X2',IER)
           DEALLOCATE (X)
           RETURN
        ENDIF

        ALLOCATE (NR(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'BP 3F, NR',IER)
          DEALLOCATE (X, X2)
           RETURN
         ENDIF

        ALLOCATE (W(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'BP 3F, W',IER)
           DEALLOCATE (X, X2, NR)
           RETURN
        ENDIF

        ALLOCATE (W2(0:N2,N,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'BP 3F, W',IER)
           DEALLOCATE (X, X2, NR, W)
           RETURN
        ENDIF
	
        ALLOCATE (SSNR(-N2:N2-1,-N2:N2-1,-N2:N2-1), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'BP 3F, W',IER)
           DEALLOCATE (X2, NR, W, W2)
           RETURN
        ENDIF

c$omp parallel do private(i,j,k)
        DO    K=1,N
           DO    J=1,N
              DO    I=0,N2
                 X(I,J,K)=CMPLX(0.0,0.0)
                 X2(I,J,K)=0.0
                 W(I,J,K)=0.0
                 W2(I,J,K)=0.0
		 NR(I,J,K)=0
              ENDDO
           ENDDO
        ENDDO
c$omp parallel do private(i,j,k)
        DO    K=-N2,N2-1
           DO    J=-N2,N2-1
              DO    I=-N2,N2-1
                 SSNR(I,J,K)=0.0
              ENDDO
           ENDDO
        ENDDO



        ALLOCATE (PROJ(NS,NS), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'BP 3F, PROJ',IER)
           DEALLOCATE (X2, NR, W, W2, SSNR)
           RETURN
        ENDIF



        ALLOCATE (BI(0:N2,N), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'BP 3F, BI',IER)
          DEALLOCATE (X2, NR, W, W2, SSNR, PROJ)
        ENDIF

#ifdef SP_MP
	LN1=LN+1
#endif
        DO    K=1,NANG
C          PRINT  *,' PROJECTION #',K

C          OPEN DESIRED FILE
           CALL FILGET(FINPAT,FINPIC,NLET,ILIST(K),IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPROJ,'O',IFORM,NSAM,NSAM,NSL,
     &                  MAXIM,'DUMMY',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           DO J=1,NS
              CALL  REDLIN(INPROJ,PROJ(1,J),NS,J)
           ENDDO
           CLOSE(INPROJ)

           CALL WINDOW2(PROJ,NS,BI,LSD,N)
           INV = +1
           CALL FMRS_2(BI,N,N,INV)
c$omp      parallel do private(i,j)
           DO  J=1,N
              DO  I=0,N2
                 BI(I,J)=BI(I,J)*(-1)**(I+J+1)
              ENDDO
           ENDDO
C
           DO  ISYM=1,MAXSYM
            IF(MAXSYM.GT.1)  THEN
C  symmetries, multiply matrices
             DMS=MATMUL(SM(:,:,ISYM),DM(:,:,K))
            ELSE
             DMS=DM(:,:,K)
            ENDIF
#ifdef SP_MP
	    DO  JT=1,LN1
C,schedule(static)
c$omp       parallel do private(j)
c$omp& shared(N,N2,JT,X,X2,W,NR,BI,DM)
             DO J=-N2+JT,N2,LN1
              CALL ONELINE3(J,N,N2,X,X2,W,W2,NR,BI,DMS)
             ENDDO
            ENDDO
#else
            DO J=-N2+1,N2
              CALL ONELINE3(J,N,N2,X,X2,W,W2,NR,BI,DMS)
            ENDDO
#endif
C   END OF SYMMETRIES LOOP
           ENDDO

C          END OF PROJECTIONS LOOP
        ENDDO

C       SYMMETRIZE PLANE 0
        DO  IZA=2,N2
           DO  IYA=2,N2
              X(0,IYA,IZA)=X(0,IYA,IZA)+CONJG(X(0,N-IYA+2,N-IZA+2))
              X2(0,IYA,IZA)=X2(0,IYA,IZA)+X2(0,N-IYA+2,N-IZA+2)
              W(0,IYA,IZA)=W(0,IYA,IZA)+W(0,N-IYA+2,N-IZA+2)
              W2(0,IYA,IZA)=W2(0,IYA,IZA)+W2(0,N-IYA+2,N-IZA+2)
              NR(0,IYA,IZA)=NR(0,IYA,IZA)+NR(0,N-IYA+2,N-IZA+2)

              X(0,N-IYA+2,N-IZA+2)=CONJG(X(0,IYA,IZA))
              X2(0,N-IYA+2,N-IZA+2)=X2(0,IYA,IZA)
              W(0,N-IYA+2,N-IZA+2)=W(0,IYA,IZA)
              W2(0,N-IYA+2,N-IZA+2)=W2(0,IYA,IZA)
              NR(0,N-IYA+2,N-IZA+2)=NR(0,IYA,IZA)

              X(0,N-IYA+2,IZA)=X(0,N-IYA+2,IZA)+CONJG(X(0,IYA,N-IZA+2))
              X2(0,N-IYA+2,IZA)=X2(0,N-IYA+2,IZA)+X2(0,IYA,N-IZA+2)
              W(0,N-IYA+2,IZA)=W(0,N-IYA+2,IZA)+W(0,IYA,N-IZA+2)
              W2(0,N-IYA+2,IZA)=W2(0,N-IYA+2,IZA)+W2(0,IYA,N-IZA+2)
              NR(0,N-IYA+2,IZA)=NR(0,N-IYA+2,IZA)+NR(0,IYA,N-IZA+2)

              X(0,IYA,N-IZA+2)=CONJG(X(0,N-IYA+2,IZA))
              X2(0,IYA,N-IZA+2)=X2(0,N-IYA+2,IZA)
              W(0,IYA,N-IZA+2)=W(0,N-IYA+2,IZA)
              W2(0,IYA,N-IZA+2)=W2(0,N-IYA+2,IZA)
              NR(0,IYA,N-IZA+2)=NR(0,N-IYA+2,IZA)
           ENDDO
        ENDDO
        DO  IYA=2,N2
           X(0,IYA,1)=X(0,IYA,1)+CONJG(X(0,N-IYA+2,1))
           X2(0,IYA,1)=X2(0,IYA,1)+X2(0,N-IYA+2,1)
           W(0,IYA,1)=W(0,IYA,1)+W(0,N-IYA+2,1)
           W2(0,IYA,1)=W2(0,IYA,1)+W2(0,N-IYA+2,1)
           NR(0,IYA,1)=NR(0,IYA,1)+NR(0,N-IYA+2,1)

           X(0,N-IYA+2,1)=CONJG(X(0,IYA,1))
           X2(0,N-IYA+2,1)=X2(0,IYA,1)
           W(0,N-IYA+2,1)=W(0,IYA,1)
           W2(0,N-IYA+2,1)=W2(0,IYA,1)
           NR(0,N-IYA+2,1)=NR(0,IYA,1)
        ENDDO
        DO  IZA=2,N2
           X(0,1,IZA)=X(0,1,IZA)+CONJG(X(0,1,N-IZA+2))
           X2(0,1,IZA)=X2(0,1,IZA)+X2(0,1,N-IZA+2)
           W(0,1,IZA)=W(0,1,IZA)+W(0,1,N-IZA+2)
           W2(0,1,IZA)=W2(0,1,IZA)+W2(0,1,N-IZA+2)
           NR(0,1,IZA)=NR(0,1,IZA)+NR(0,1,N-IZA+2)

           X(0,1,N-IZA+2)=CONJG(X(0,1,IZA))
           X2(0,1,N-IZA+2)=X2(0,1,IZA)
           W(0,1,N-IZA+2)=W(0,1,IZA)
           W2(0,1,N-IZA+2)=W2(0,1,IZA)
           NR(0,1,N-IZA+2)=NR(0,1,IZA)
        ENDDO

        CALL SSNR3D(X,X2,W,W2,NR,SSNR,N2,N)
        IFORM = 3
	MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FINPIC,INPROJ,'U',IFORM,N,N,N,
     &           MAXIM,'3-D SSNR',.FALSE.,IRTFLG)

        IF (IRTFLG .EQ. 0) 
     &    CALL WRITEV(INPROJ,SSNR,N,N,N,N,N)

	CLOSE(INPROJ)
        DEALLOCATE (X,W,PROJ,BI,X2,W2,NR,SSNR)

        END


C       --------------------- ONELINE3 ---------------------------------

        SUBROUTINE  ONELINE3(J,N,N2,X,X2,W,W2,NR,BI,DM)

        DIMENSION      W(0:N2,N,N),W2(0:N2,N,N)
	DIMENSION      X2(0:N2,N,N),NR(0:N2,N,N)
        COMPLEX        BI(0:N2,N),X(0:N2,N,N),BTQ
        DIMENSION      DM(6)
        PARAMETER      (LTAB=5000)
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
		       IF(TZ.NE.0)  THEN
                       DO  LY=-LN2,LN2
                          IYP=IYN+LY
                          IF (IYP.GE.0) THEN
                             IYA=IYP+1
                          ELSE
                             IYA=N+IYP+1
                          ENDIF
                          TY = TABI(NINT(ABS(YNEW-IYP)*FLTB))*TZ
		          IF(TY.NE.0)  THEN
                          DO  IXP=LB+IXN,LN2+IXN
C                            GET THE WEIGHT
C                            WG=FINF(XNEW,YNEW,ZNEW,IXT,IYT,IZT)
                             WG=TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY
 			    IF(WG.NE.0.0)  THEN
                            X(IXP,IYA,IZA)=X(IXP,IYA,IZA)+BTQ*WG
       X2(IXP,IYA,IZA)=X2(IXP,IYA,IZA)+(REAL(BTQ)**2+AIMAG(BTQ)**2)*WG
                             W(IXP,IYA,IZA)=W(IXP,IYA,IZA)+WG
                             W2(IXP,IYA,IZA)=W2(IXP,IYA,IZA)+WG*WG
                             NR(IXP,IYA,IZA)=NR(IXP,IYA,IZA)+1
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
		      IF(TZ.NE.0)  THEN
                      DO  LY=-LN2,LN2
                         IYP=IYN+LY
                         IF (IYP.GT.0) THEN
                            IYT=N-IYP+1
                         ELSE
                            IYT=-IYP+1
                         ENDIF
                        TY=TABI(NINT(ABS(YNEW-IYP)*FLTB))*TZ
		        IF(TY.NE.0)  THEN
                        DO  IXP=IXN-LN2,-1
C                           GET THE WEIGHT
C                           WG=FINF(XNEW,YNEW,ZNEW,IXT,IYT,IZT)
                            WG=TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY
			  IF(WG.NE.0.0)  THEN
                          X(-IXP,IYT,IZT)=X(-IXP,IYT,IZT)+CONJG(BTQ)*WG
       X2(-IXP,IYT,IZT)=X2(-IXP,IYT,IZT)+(REAL(BTQ)**2+AIMAG(BTQ)**2)*WG
                            W(-IXP,IYT,IZT)=W(-IXP,IYT,IZT)+WG
                            W2(-IXP,IYT,IZT)=W2(-IXP,IYT,IZT)+WG*WG
                            NR(-IXP,IYT,IZT)=NR(-IXP,IYT,IZT)+1
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

C       ------------------- WINDOW2 -------------------------------

        SUBROUTINE WINDOW2(PROJ,L,BI,LSD,N)

        DIMENSION  PROJ(L,L),BI(LSD,N)
        DOUBLE     PRECISION QS

        KLP=0
        R=L/2
        QS=0.0D0

        CALL ASTA(PROJ,L,R,QS,KLP)
        QS = QS/REAL(KLP)
c$omp parallel do private(i,j)
        DO  J=1,N
           DO  I=1,N
              BI(I,J)=PROJ(I,J)-QS
           ENDDO
        ENDDO

        END

C       ------------------- SSNR3D -------------------------------

        SUBROUTINE  SSNR3D(X,X2,W,W2,NR,SSNR,N2,N)
        PARAMETER (NDLI=6)
        DIMENSION  W(0:N2,N,N),W2(0:N2,N,N),X2(0:N2,N,N),NR(0:N2,N,N)
        COMPLEX    X(0:N2,N,N)
	INCLUDE  'CMBLOCK.INC'
	DIMENSION  DLIST(NDLI),NRS(0:N2)
	DIMENSION  SSNR(-N2:N2-1,-N2:N2-1,-N2:N2-1)
	DOUBLE PRECISION  SIGNAL(0:N2),FR(0:N2),TMP,TMPS
C	DOUBLE PRECISION  FRS(0:N2)
	DOUBLE PRECISION  SW(0:N2)
        DATA  LUN9/99/
C  X(i,j,k) - weighted sum of Fs.  SUM(wF)/SUM(w)
C  X2(i,j,k) - weighted sum of F^2.  SUM(wF^2)
C  W(i,j,k) - sum of weights.  SUM(w)
C  W2(i,j,k) - sum of squares of weights.  SUM(w^2)
C  NR(i,j,k) - number of elements added at a given freq.
C  NRS: length - frequencies; NRS(i) total number of elements added within i'th shell. SUM(NR)
C  SIGNAL: length - frequencies; Sum over the shell of (SUM(wF)/SUM(w))^2
C  FR: length - frequencies; FR(i) -  
	LR=0
	NRS=0
	SIGNAL=0.0D0
	FR=0.0D0
C	FRS=0.0
	SW=0.0D0
C
       DO  K=1,N
        KK=K-1
        IF(KK.GE.N/2)  KK=KK-N
           DO  J=1,N
           JJ=J-1
           IF(JJ.GE.N/2)  JJ=JJ-N
              DO  I=0,N2
		IF(W(I,J,K).LT.1.0)  THEN
		 W(I,J,K)=0.0
		 X(I,J,K)=0.0
		 NR(I,J,K)=0
		 X2(I,J,K)=0.0
		ENDIF
              ENDDO
           ENDDO
        ENDDO
        write(nout,*)  sum(nr)
C
        DO  K=1,N
        KK=K-1
        IF(KK.GE.N/2)  KK=KK-N
           DO  J=1,N
           JJ=J-1
           IF(JJ.GE.N/2)  JJ=JJ-N
              DO  I=0,N2
               IF(NR(I,J,K).NE.0)  THEN
		X(I,J,K)=X(I,J,K)*(-1)**(I+J+K)/W(I,J,K)
	        IF(.NOT.(I.EQ.0.AND.JJ.LT.0)) THEN
                 PII=SQRT((REAL(KK)/REAL(N/2))**2+
     &           (REAL(JJ)/REAL(N/2))**2+(REAL(I)/REAL(N2))**2)
                 IF(PII.LE.1.0)  THEN
                  L=MIN0(MAX0(NINT(PII*N2),0),N2)
		  TMPS=W(I,J,K)/W2(I,J,K)*CABS(X(I,J,K))**2
                  SIGNAL(L)=SIGNAL(L)+TMPS
	        IF(NR(I,J,K).GT.1)  THEN
         	  TMP=
     &		( X2(I,J,K)-W(I,J,K)*
     &		( REAL(X(I,J,K))**2+AIMAG(X(I,J,K))**2 )  )
     &		/W(I,J,K)/(1.-W2(I,J,K)/W(I,J,K)**2)
                ELSE
		  TMP=TMPS
		ENDIF
         	  FR(L)=FR(L)+TMP/W(I,J,K)
		  SSNR(I,JJ,KK)=DMAX1(0.0D0,W(I,J,K)*TMPS/TMP-1.0D0)
	 	  SSNR(-I,-JJ,-KK)=SSNR(I,JJ,KK)
C		  FRS(L)=FRS(L)+TMP
		  SW(l)=SW(l)+W2(I,J,K)/W(I,J,K)
 		  NRS(L)=NRS(L)+1
 		 ENDIF
 		 ENDIF
	       ENDIF
              ENDDO
           ENDDO
        ENDDO
C SAVE RESULTS
	print  *,nrs
        DO   L=0,N2
           IF(NRS(L).NE.0)  THEN
              DLIST(1)=L+1
              DLIST(2)=(REAL(L)+0.5)/REAL(N2)*0.5
              DLIST(3)=SW(L)/NRS(L)*DMAX1(0.0D0,SIGNAL(L)/FR(L)-1.0D0)
 	      DLIST(4)=SIGNAL(L)
	      DLIST(5)=FR(L)
C 	      DLIST(4)=DLIST(3)/(1.0+DLIST(3))
C	      DLIST(5)=FRS(L)/NRS(L)
	      DLIST(6)=NRS(L)
              CALL  SAVD(LUN9,DLIST,NDLI,IRTFLG)
           ENDIF
        ENDDO
        CLOSE(LUN9)
        CALL  SAVDC
        END
