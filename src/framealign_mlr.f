C++*********************************************************************
C 
C  FRAMEALIGN_MLR                                 JAN 15 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2015  Health Research Inc.,                         *
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
C  FRAMEALIGN_MLR
C
C  PURPOSE: DETERMINE ALIGNMENT FOR FRAME SERIES USING 
C           MULTIPLE LINEAR REGRESSION
C
C  CALLS:   MINPACK ROUTINES LMDIF & LMDIF1
C                      
C--*******************************************************************

	SUBROUTINE FRAMEALIGN_MLR

        IMPLICIT  NONE

	INCLUDE 'CMBLOCK.INC'
      	INCLUDE 'F90ALLOC.INC'
	INCLUDE 'CMLIMIT.INC'

        INTEGER, PARAMETER            :: MAXFRAMES = 40
        INTEGER, PARAMETER            :: MAXFN     = 39 * 39 / 2
      	REAL                          :: ALIXY(MAXFN)
      	REAL                          :: COEF(MAXFRAMES,MAXFN)
        COMMON /COEFFS/COEF,ALIXY

        CHARACTER(LEN=MAXNAM)         :: DOCIN
	CHARACTER(LEN=MAXNAM)         :: DOCOUT
	CHARACTER(LEN=MAXNAM)         :: COMMENT

        INTEGER, PARAMETER            :: NLIST = 7
      	REAL	                      :: DLIST(NLIST)  
	REAL, POINTER                 :: DOCBUF(:,:)
      	DOUBLE PRECISION, ALLOCATABLE :: FVEC(:)
      	DOUBLE PRECISION, ALLOCATABLE :: X(:),     Y(:)
      	DOUBLE PRECISION, ALLOCATABLE :: X_OK(:),  Y_OK(:)
      	DOUBLE PRECISION, ALLOCATABLE :: X1(:),    Y1(:)
      	DOUBLE PRECISION, ALLOCATABLE :: SHX1(:),  SHY1(:)
      	DOUBLE PRECISION, ALLOCATABLE :: SHXC(:),  SHYC(:)
      	DOUBLE PRECISION, ALLOCATABLE :: WA(:)
      	INTEGER, ALLOCATABLE          :: IWA(:)
      	REAL                          :: ALIX(MAXFN),ALIY(MAXFN)
      	REAL, ALLOCATABLE             :: ERRX(:),ERRY(:)
      	REAL                          :: XT,YT
      	REAL                          :: AVG_ERR_X,AVG_ERR_Y

      	DOUBLE PRECISION              :: XVAL, YVAL
        LOGICAL                       :: ASKNAM, GETSIZE  
        INTEGER                       :: MAXXT,MAXYT,IRTFLG
        INTEGER                       :: LUNDOCOUTT,NLET 

        LOGICAL                       :: ADDEXT,ISOLD
        LOGICAL                       :: APPEND,MESSAGE,NEWFILE

	INTEGER                       :: M,N,N1,N2,I,MWANT,K,NF,NFC
        INTEGER                       :: MP,LWA,NFGO,NFEND,J,M_OK
        INTEGER                       :: NLETD,IKEY,INFOX,INFOY,ITER

	DOUBLE PRECISION              :: TOL

        INTEGER, PARAMETER            :: LUNDOC    = 80
        INTEGER, PARAMETER            :: LUNDOCOUT = 81

	DOUBLE PRECISION, PARAMETER   :: ONE_D = 1.0D0

	EXTERNAL                      :: fcn_mlr

        ! MAX ACCEPTATBLE RESIDUAL ERROR
        REAL, PARAMETER               :: FMIN_ERR = 1.2  

	MAXXT   = 0  ! MAXXT & MAXYT SET TO ZERO ON ENTRY, OTHER-WISE 
	MAXYT   = 0  ! THEY ARE MAX VALUES FOR DOCBUF ARRAY CREATION
	ASKNAM  = .TRUE.
	GETSIZE = .TRUE.

	CALL GETDOCDAT("FRAME ALIGNMENT DOC", ASKNAM,DOCIN,
     &                 LUNDOC,GETSIZE,MAXXT, MAXYT,DOCBUF,IRTFLG)
	IF (IRTFLG.NE.0) RETURN	

        IF (MAXYT < 4 ) THEN
           CALL ERRT(102,'TOO FEW FUNCTIONS',MAXYT)
           RETURN
         
        ELSEIF (MAXYT > MAXFN ) THEN
           CALL ERRT(102,'INCREASE MAXFRAMES + RECOMPILE',MAXYT)
           RETURN

        ELSEIF (MAXXT < 4 ) THEN
           CALL ERRT(102,'TOO FEW REGISTERS IN DOC FILE',MAXXT)
           RETURN
        ENDIF

        N1 = MAXVAL(DOCBUF(2,:))
        N2 = MAXVAL(DOCBUF(3,:))

C       GET NF = NUMBER OF SUBFRAMES == (NUMBER OF VARIABLES + 1)
        NF = MAX(N1,N2)

C       GET N = NUMBER OF VARIABLES
        N = NF - 1

C       GET M = NUMBER OF FUNCTIONS IN USE
        M = MAXYT

C       FIND MP = NUMBER OF POSSIBLE FUNCTIONS
        MP = (N-1) * N / 2

C       LWA IS A POSITIVE INTEGER VARIABLE >=  M*N+5*N+M.
C       WA  IS A WORK ARRAY OF LENGTH LWA.
        LWA = M * N + 5 * N + M

C       IWA IS AN INTEGER WORK ARRAY OF LENGTH N.
 
	ALLOCATE(X1(N),    Y1(N), 
     &           X(N),     Y(N),    
     &           X_OK(N),  Y_OK(N),  
     &           SHXC(NF), SHYC(NF),
     &           SHX1(NF), SHY1(NF),
     &           FVEC(M),  
     &           ERRX(M),  ERRY(M),
     &           WA(LWA),
     &           IWA(N),   STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           MWANT = 13*N + 4*M + 2*LWA + 4*NF + MAXXT*MAXYT
           CALL ERRT(46,'X1...',MWANT)
           RETURN
        ENDIF

C       GET COEFFICIENTS FOR THE FUNCTIONS
        COEF = 0.0   ! ARRAY OP
      
        DO I = 1,MAXYT
           NFGO  = DOCBUF(2,I) 
           NFEND = DOCBUF(3,I)
  
           !write(6,*) 'i,nfgo,nfend:',i,nfgo,nfend

C          COEF LISTS PRESENCE IN EACH FUNCTION
           DO J = NFGO,NFEND-1
              COEF(J,I) = 1.0
              COEF(J,I) = 1.0
           ENDDO

           ALIX(I) = DOCBUF(4,I)
           ALIY(I) = DOCBUF(5,I)

           IF (NFGO == 1 .AND. NFEND == NF) THEN
C             NEED INITIAL ESTIMATES OF UNKOWNS
              XVAL = ALIX(I) / (NFEND -1)
              YVAL = ALIY(I) / (NFEND -1)
              X1   = XVAL            ! ASSUME ALL X1 THE SAME
              Y1   = YVAL            ! ASSUME ALL Y1 THE SAME
              !write(6,*) 'i,nfgo,nfend:',i,nfgo,nfend
           ENDIF
        ENDDO

        TOL  = SQRT(EPSILON(ONE_D))  ! SQUARE ROOT OF MACHINE PRECISION
        !IF (VERBOSE) WRITE(NOUT,*) ' '

C       LOOP UNTIL RESIDUAL ERROR OF ESTIMATES IS SATISFACTORY

        ITER = 0
        DO WHILE (M > NF)           ! NEED AT LEAST NF FUNCTIONS

           ITER = ITER + 1

           X     = X1(:)            ! STARTING VALUES FOR X 
           ALIXY = ALIX             ! PASSED TO FCN IN COMMON
           CALL LMDIF1(FCN_MLR,M,N, X,FVEC, TOL, INFOX, IWA,WA,LWA)

           Y     = Y1(:)            ! STARTING VALUES FOR Y
           ALIXY = ALIY             ! PASSED TO FCN IN COMMON
           CALL LMDIF1(FCN_MLR,M,N, Y,FVEC, TOL, INFOY, IWA,WA,LWA)

           !DO I = 1,NF-1            ! LOOP OVER ALL FRAMES > FIRST
           !  write(6,95) i,x(i),y(i)
c95         !  format('Serial Offset',i2, f6.2,f6.2)
           !ENDDO

C          FIND RESIDUAL ERRORS FOR EACH FUNCTION
           AVG_ERR_X = 0
           AVG_ERR_Y = 0

           IF (VERBOSE) THEN   ! KLUDGE FOR OUTPUT FORMATTING
              DO I = 1,M
                 N1      = DOCBUF(2,I)
                 N2      = DOCBUF(3,I)

                 XT      = SUM(X(N1:N2-1))  
                 YT      = SUM(Y(N1:N2-1))  

                 ERRX(I) = ABS(XT - DOCBUF(4,I))
                 ERRY(I) = ABS(YT - DOCBUF(5,I))

                 AVG_ERR_X = AVG_ERR_X + ERRX(I)
                 AVG_ERR_Y = AVG_ERR_Y + ERRY(I)
              ENDDO

              AVG_ERR_X = AVG_ERR_X / M
              AVG_ERR_Y = AVG_ERR_Y / M
         
              WRITE(NOUT,95)ITER,M,INFOX,INFOY,AVG_ERR_X,AVG_ERR_Y

95            FORMAT('  ITERATION: ',I5,' FUNCTIONS: ',I5,
     &               '  LMDIF INFO X & Y: ',I0,', ',I0,
     &               '   AVG. RESIDUAL ERROR: ',F6.2,', ',F6.2)
           ENDIF

           DO I = 1,M
              N1      = DOCBUF(2,I)
              N2      = DOCBUF(3,I)

              XT      = SUM(X(N1:N2-1))  
              YT      = SUM(Y(N1:N2-1))  

              ERRX(I) = ABS(XT - DOCBUF(4,I))
              ERRY(I) = ABS(YT - DOCBUF(5,I))

              IF (VERBOSE .AND. ITER == 1) THEN
 
                 WRITE(NOUT,94) I,N1,N2, XT,YT,ERRX(I),ERRY(I)
94               FORMAT('  ',I4,':  ',I2,'..',I2,
     &                  '   Shift: ',         F6.2,' ,',F6.2,
     &                  '   Residual error: ',F6.2,' ,',F6.2)
              ENDIF
           ENDDO
           IF (VERBOSE) WRITE(NOUT,*) ' '

C          PRESERVE CURRENT VALUES IN CASE TOO MANY DISCARDED
           X_OK      = X   
           Y_OK      = Y
           M_OK      = M

C          DISCARD FUNCTIONS HAVING LARGE RESIDUAL ERRORS
           DO I = 1,M 

              IF (ERRX(I) > FMIN_ERR  .OR. 
     &            ERRY(I) > FMIN_ERR) THEN

C                HIGH ERROR, DISCARD THIS FUNCTION (EQUATION)
                 ALIX(1:I-1)     = ALIX(1:I-1)     ! USED IN FCN 
                 ALIX(I:M-1)     = ALIX(I+1:M) 

                 ALIY(1:I-1)     = ALIY(1:I-1)     ! USED IN FCN 
                 ALIY(I:M-1)     = ALIY(I+1:M) 

                 COEF(:,1:I-1)   = COEF(:,1:I-1)   ! USED IN FCN
                 COEF(:,I:M-1)   = COEF(:,I+1:M) 

                 DOCBUF(:,1:I-1) = DOCBUF(:,1:I-1) ! FOR LABELS
                 DOCBUF(:,I:M-1) = DOCBUF(:,I+1:M) 

                 M = M - 1
              ENDIF
           ENDDO
           IF (M == M_OK) EXIT    ! ALL FUNCTIONS HAVE OK ERRORS
        
        ENDDO

C       RETRIEVE LAST USABLE X, Y, ETC
        M = M_OK
        X = X_OK(1:M)
        Y = Y_OK(1:M)

        IF (VERBOSE) THEN
           DO I = 1,M
              N1      = DOCBUF(2,I)
              N2      = DOCBUF(3,I)

              XT      = SUM(X(N1:N2-1))  
              YT      = SUM(Y(N1:N2-1))  

              ERRX(I) = ABS(XT - DOCBUF(4,I))
              ERRY(I) = ABS(YT - DOCBUF(5,I))

              WRITE(NOUT,94) I,N1,N2, XT,YT,ERRX(I),ERRY(I)

           ENDDO
           WRITE(NOUT,*) ' '
        ENDIF


C       FIND REFINED SHIFTS VS INITAL AND CENTRAL FRAMES 
        
        NFC  = 1                 ! VS INITIAL FRAME

        SHX1 = 0.0               ! ARRAY OP
        SHY1 = 0.0               ! ARRAY OP

        DO I = 2,NF              ! LOOP OVER ALL FRAMES > FIRST
          SHX1(I) = SUM(X(1:I-1))
          SHY1(I) = SUM(Y(1:I-1))

          !write(6,92) i,shx1(i),shy1(i)
92        format('1..',i2, f6.2,f6.2)
        ENDDO

        NFC  = NF / 2            ! VS CENTRAL FRAME

        SHXC = 0.0               ! ARRAY OP
        SHYC = 0.0               ! ARRAY OP

        DO I = 1,NFC-1         ! LOOP OVER ALL FRAMES < CENTRAL
          SHXC(I) = -SUM(X(I:NFC-1))
          SHYC(I) = -SUM(Y(I:NFC-1))
        ENDDO

        DO I = NFC+1,NF         ! LOOP OVER ALL FRAMES > CENTRAL
          SHXC(I) = SUM(X(NFC:I-1))
          SHYC(I) = SUM(Y(NFC:I-1))

          !write(6,93) i,shxc(i),shyc(i)
93        format('1..',i2, f6.2,f6.2)
        ENDDO

C       SAVE REFINED SHIFTS IN OUTPUT DOC FILE
        ADDEXT  = .TRUE.                         
        ASKNAM  = .TRUE.                         
        ISOLD   = .FALSE.                        
        APPEND  = .TRUE.                         
        MESSAGE = .TRUE.                        
        IRTFLG  = -8         ! NO IC USE         

        CALL OPENDOC(DOCOUT,ADDEXT,NLETD,LUNDOCOUT,LUNDOCOUTT,ASKNAM,
     &        'SHIFT DOC',ISOLD,APPEND,MESSAGE,NEWFILE,IRTFLG)  
        IF (IRTFLG .NE. 0) GOTO 9999                              

C                  123456789 123456789 123456789 123456789 123456789 123456789
        COMMENT = '  SHIFTS vs FRAME: 1 AND FRAME: '
        CALL INTTOCHAR(NFC,COMMENT(33:35),NLET,1)
        CALL LUNDOCPUTCOM(LUNDOCOUTT,COMMENT(1:36),IRTFLG)                     

C                 123456789 123456789 123456789 123456789 123456789 123456789
        COMMENT ='        FRAME,' //
     &           '      SX vs 1,      SY vs 1,'  //
     &           '     SX vs   ,      SY vs   '                         
        CALL INTTOCHAR(NFC,COMMENT(54:55),NLET,2)
        CALL INTTOCHAR(NFC,COMMENT(69:70),NLET,2)
        CALL LUNDOCPUTCOM(LUNDOCOUTT,COMMENT(1:72),IRTFLG)                     

        DO I = 1, NF

           IKEY     = I 
           DLIST(1) = I
           DLIST(2) = SHX1(I)    ! VS INITIAL FRAME
           DLIST(3) = SHY1(I)    ! VS INITIAL FRAME
           DLIST(4) = SHXC(I)    ! VS CENTRAL FRAME
           DLIST(5) = SHYC(I)    ! VS CENTRAL FRAME

           CALL LUNDOCWRTDAT(LUNDOCOUTT,IKEY,DLIST,5,IRTFLG)

        ENDDO

9999    CLOSE(LUNDOCOUT)

	END

C       -------------------------- FCN_MLR -----------------------


C       TASK SPECIFIC SUBROUTINE FOR MULTIPLE LINEAR REGRESSION
C       FOR FINDING 'BEST' ALIGNMENT SHIFTS 

	SUBROUTINE FCN_MLR(M,N,X,FVEC,IFLAG)

        IMPLICIT NONE

        INTEGER             :: M,N,IFLAG
	DOUBLE PRECISION    :: FVEC(M)
	DOUBLE PRECISION    :: X(N)

        INTEGER, PARAMETER  :: MAXFRAMES = 40
        INTEGER, PARAMETER  :: MAXFN     = 39 * 39 / 2
      	REAL                :: COEF(MAXFRAMES,MAXFN)
      	REAL                :: ALIXY(MAXFN) 
        COMMON /COEFFS/COEF,ALIXY

        INTEGER             :: I,J

	DO I=1,M         ! LOOP OVER ALL FUNCTIONS

           FVEC(I) = ALIXY(I)

           DO J= 1,N     ! LOOP OVER ALL VARIABLES 
              FVEC(I) = FVEC(I) - COEF(J,I) * X(J)
           ENDDO
        ENDDO

	END 



#ifdef NEVER
c ---------------------- Unused below -------------------
           write(6,*) 'coef:'
           write(6,*)  coef(:,1)

           write(6,*) 'm,n,fvec:', m,n
           write(6,*)  fvec
 
           write(6,*) 'docbuf:'
           write(6,*) docbuf(:,1)
           write(6,*) 'coef:'
           write(6,*) coef(:,1)
           write(6,*) 'x1,y1:'
           write(6,*) x1(1),y1(1)
           write(6,*) 'alix:'
           write(6,*) alix(1),aliy(1)
           stop
#endif
