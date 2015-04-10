
C **********************************************************************
C
C   SPCA3     USED ALLOCATE                  JAN 2001 ARDEAN LEITH
C             ENLARGED MADAR                 MAR 2002 ARDEAN LEITH
C             TSIZE BUG                      AUG 2002 ARDEAN LEITH
C             REWRITTEN                      SEP 2003 ARDEAN LEITH
C             REWRITTEN                      MAR 2006 ARDEAN LEITH
C             LUNE HEADER SPACE              MAR 2009 ARDEAN LEITH
C             PIXEL FILE FDUM                JUN 2009 ARDEAN LEITH
C             PCASCOOR: WEIGHTP/NUMIM        FEB 2010 ARDEAN LEITH
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
c PURPOSE:  CARRIES OUT ITERATIVE PCA 
C
C   SPCA3(NUMIM,NFAC,NPIX, INUMBR,LUNS,LUNI,LUNP, LUNE, 
C         WEIGHTI, WEIGHTP, SUMW, VARP, LUNT)
C
C       NUMIM	  NUMBER OF IMAGES                               (INPUT)
C	NFAC 	  NUMBER OF EIGENVECTORS REQUESTED               (INPUT)
C	NPIX      NUMBER OF ACTIVE PIXELS PER IMAGE              (INPUT)
C       INUMBR	  IMAGE NUMBER LIST                              (INPUT)
C	LUNS	  SEQUENTIAL IMAGE FILE I/O UNIT (INPUT FILE)    (INPUT)
C	LUNI      IMAGE COORDINATE FILE I/O UNIT (OUTPUT FILE)   (INPUT)
C	LUNP      PIXEL COORDINATE FILE I/O UNIT (OUTPUT FILE)   (INPUT)
C	LUNE	  EIGENVALUE FILE I/O UNIT (OUTPUT FILE)         (INPUT)
C	WEIGHTP	  SUM OF PIXEL VALUES AT THIS PIXEL              (INPUT)
C	WEIGHTI	  SUM OF PIXEL VALUES IN THIS IMAGE              (INPUT)
C	SUMW      SUM OF ALL THE PIXEL VALUES IN ALL IMAGES      (INPUT)
C	VARP	  VARIANCE OF PIXEL VALUES AT THIS PIXEL         (INPUT)
C	LUNT	  TEMP FILE I/O UNIT (OUTPUT FILE)               (INPUT)
C
C       ALL FILES ARE FORMATTED EXCEPT FOR LUNS!
C
C  CALL TREE:
C     JPMSK1 ---> SPCA3 --->  PCA_STOCHA
C             iter            | PCA_GSMOD
C                             | PCA_ITMPOW
C                             | PCA_CPROJ
C
C                       --->  ORDERE
C                       --->  PCASCOOR

C **********************************************************************

        SUBROUTINE SPCA3(NUMIM,   NFAC,    NPIX, INUMBR, 
     &                   LUNS,    LUNI,    LUNP, LUNE, 
     &                   WEIGHTI, WEIGHTP, SUMW, VARP, LUNT)
 
        INCLUDE 'CMBLOCK.INC'

        INTEGER, DIMENSION(NUMIM)          :: INUMBR
        REAL, DIMENSION(NUMIM)             :: WEIGHTI
        REAL, DIMENSION(NPIX)              :: WEIGHTP,VARP

        REAL, ALLOCATABLE, DIMENSION(:,:)  :: EVECTS,BLBZ,BLAD
        REAL, ALLOCATABLE, DIMENSION(:)    :: BLU,EVALS,BLCO

        CHARACTER(LEN=1)                   :: NULL = CHAR(0)

C       AUTOMATIC ARRAYS
        REAL, DIMENSION(NFAC)              :: BLW

        WRITE(NOUT,*) ' ITERATIVE SOLUTION --- '

C 	COMPUTE WEIGHTP AND VARIANCE
        DO J = 1, NPIX
           WEIGHTP(J) = WEIGHTP(J) / NUMIM
           VARP(J)    = VARP(J) - (WEIGHTP(J) **2) * NUMIM
        ENDDO
        TRACE = SUM(VARP)

        JBASE = NFAC + 3

        ALLOCATE(
     &     EVECTS(NPIX,JBASE), 
     &     BLU(NPIX), 
     &     BLBZ(NPIX,JBASE),    
     &     EVALS(NPIX),               ! should be jbase ? al
     &     BLAD(JBASE,JBASE),   
     &     BLCO(NPIX),   
     &     STAT=IRTFLG)

	IF (IRTFLG .NE. 0)  THEN  
           MWANT =2*(NPIX*JBASE) + (JBASE*JBASE) + 3*(NPIX)
           WRITE(NOUT,*)
     &          '*** INSUFFICIENT MEMORY FOR ITERATIVE SOLUTION'
           CALL ERRT(46, 'SPCA3 NEEDED', MWANT)
           RETURN
        ENDIF

        CALL OPAUXFILE(.FALSE.,'TEMP_SCRATCH',NULL,-LUNT,0,
     &                  'S',' ',.TRUE.,IRTFLG)
     
        CALL PCA_STOCHA(JBASE,  NPIX, NUMIM, NFAC, 
     &	                EVECTS, BLU,  TRACE, BLBZ,
     &	                EVALS,  BLAD, LUNS,  LUNT,
     &                  WEIGHTP,VARP)

        CLOSE(LUNT,STATUS='DELETE')

C       PUT EIGENVECTORS IN DECENDING ORDER
        CALL ORDERE(NPIX,NFAC,EVECTS,EVALS)
  
        WRITE(NOUT,*) ' ORDER OF EIGENVALUES:'
        WRITE(NOUT,90) (EVALS(I),I=1,NFAC)  
90      FORMAT(8(1X,G10.3))

        CALL PCASCOOR(NFAC,   JBASE, NPIX,    NUMIM,   INUMBR,
     &                EVECTS, EVALS, WEIGHTI, WEIGHTP, SUMW,  
     &                BLCO,   BLU,   BLW,     LUNS,    LUNI,  LUNP)
        	
C       WRITE EIGEN FILE HEADER
        KIND_PCA = 1
        WRITE(LUNE,91) NFAC,SUMW,TRACE,KIND_PCA,NPIX
91      FORMAT(I10,' ',1PG12.5,' ',1PG12.5,' ',I10,' ',I10)

C       SAVE EIGENVALUES TO _EIG FILE (ONE FACTOR PER LINE)
	CUL = 0.0
	DO I = 1, NFAC
           PER = 100.0 * EVALS(I) / TRACE
           CUL = CUL + PER

           WRITE(LUNE,92) EVALS(I),PER,CUL
 92        FORMAT(1PE12.5,'  ',E12.5,'  ', E12.5)
        ENDDO

C       SAVE EIGENVECTORS ARRAY TO _EIG FILE 
	DO I = 1,NPIX
           WRITE(LUNE,93) (EVECTS(I,J),J=1,JBASE)
93         FORMAT(10(1PG12.5,' '))
        ENDDO

        IF (ALLOCATED(EVECTS)) DEALLOCATE(EVECTS) 
        IF (ALLOCATED(BLU))    DEALLOCATE(BLU) 
        IF (ALLOCATED(BLBZ))   DEALLOCATE(BLBZ)    
        IF (ALLOCATED(EVALS))  DEALLOCATE(EVALS) 
        IF (ALLOCATED(BLAD))   DEALLOCATE(BLAD)   
        IF (ALLOCATED(BLCO))   DEALLOCATE(BLCO)   
 
        RETURN
        END 

C *****************************  ORDERE ********************************

        SUBROUTINE ORDERE(NPIX,NFAC,S,D)

        DIMENSION  S(NPIX,NFAC),D(NFAC)

C       SORT THE EIGENSTUFF
C       OK. I KNOW IT'S NOT EFFICIENT, BUT SURE SAVES SPACE ...ml?

        LOGICAL :: L

1       L = .FALSE.

        DO I=2,NFAC
          IF (D(I-1) .LT. D(I))  THEN
            L      = .TRUE.
            T      = D(I)
            D(I)   = D(I-1)
            D(I-1) = T
            DO J=1,NPIX
              T        = S(J,I)
              S(J,I)   = S(J,I-1)
              S(J,I-1) = T
	    ENDDO
          ENDIF
	ENDDO

        IF (L) GOTO  1

	END

C++*********************************************************************
C
C PCASCOOR          CREATED, AUTHOR       11/01/1993 MAHIEDDINE LADJADJ
C                   REWRITTEN             03/06/2006 ArDean Leith
C
C--*********************************************************************
C
C PURPOSE: SAVES TO FILES
C
C PARAMETERS:
C	NFAC 	NUMBER OF EIGENVECTORS REQUESTED              (INPUT)
C       JBASE   NFAC  + 3                                     (INPUT)
C	NPIX    NUMBER OF PIXELS PER IMAGE                    (INPUT) 
C       NUMIM	NUMBER OF IMAGES                              (INPUT) 
C       EVECTS  EIGENVECTORS ARRAY                            (INPUT)                   (     )
C       EVALS   EIGENVALUES ARRAY                             (INPUT)
C	WEIGHTP	SUM OF PIXEL VALUES AT THIS PIXEL             (INPUT)
C       SUMW    SUM OF  WEIGHT PIXEL VALUES  ?????                 ()
C       CO                                                    (     )
C       U                                                     (     )
C       W                                                     (     )
C	LUNS    SEQUENTIAL IMAGE FILE I/O UNIT (INPUT FILE)   (INPUT)
C	LUNI    IMAGE COORDINATE FILE I/O UNIT (OUTPUT FILE)  (INPUT)
C	LUNP    PIXEL COORDINATE FILE I/O UNIT (OUTPUT FILE)  (INPUT)
C
C	S(,) HAS THE EIGENVECTORS (COLUMN) OF X'X AND D() HAS THE
C       EIGENVALUES. X(I,*)=U() J=1,JTOT W/ I=1,ITOT
C         
C--*********************************************************************

        SUBROUTINE PCASCOOR(NFAC,JBASE,NPIX,NUMIM,INUMBR,
     &             EVECTS,EVALS,WEIGHTI,WEIGHTP,SUMW,
     &             CO,BLU,BLW,LUNS,LUNI,LUNP)

        REAL    :: EVECTS(NPIX,JBASE), EVALS(NPIX)
        REAL    :: CO(NPIX), BLU(NPIX), BLW(NFAC)
        REAL    :: WEIGHTP(NPIX), WEIGHTI(NUMIM)
        INTEGER :: INUMBR(NUMIM)

C       POSITION _SEQ, _IMC, & _PIX FILES
        CALL REW(LUNS,1)
        REWIND(LUNI)
        READ(LUNI,*)  IDUM, IDUM, IDUM, IDUM, IDUM, IDUM

        REWIND(LUNP)
        READ(LUNP,*) IDUM, IDUM, IDUM, IDUM, IDUM
        
C       INITIALIZE CO ARRAY TO ZERO
        CO = 0.0

        ACT  = 0.0
        FDUM = 0.0

C       WRITE _IMC  DATA
        DO I  =  1, NUMIM
C         READ THE WHOLE IMAGE  IN BLU ARRAY FROM _SEQ FILE.
          READ(LUNS) (BLU(J),J=1,NPIX), FIM

C         COORDINATES  TO ORIGIN FOR THE ROWS
C         [BLW] = [X -<X>].[EIGENVECTORS]   MATIX NOTATION. 
C         <Xi> AVERAGE VALUE OF PIXELS = SUM (XIJ),J=1,NUMIM

          DO  K=1,NFAC                                                   
             BLW(K) =  0.0 
             DO J=1,NPIX
                BLW(K) = BLW(K) + (BLU(J) - WEIGHTP(J)) * EVECTS(J,K)
             ENDDO
          ENDDO

C         DISTANCES TO ORIGIN FOR THE ROWS                    
          DOR = 0.0   
          DO J=1,NPIX
             DOR = DOR + (BLU(J) - WEIGHTP(J)) **2
          ENDDO

          PIA = WEIGHTI(I) / SUMW   
          FIM = INUMBR(I)
          WRITE(LUNI,90) (BLW(K),K=1,NFAC), PIA, DOR, FIM, ACT
        ENDDO

C       WRITE _PIX  DATA,  COORDINATES  FOR THE COLS
        DO J=1, NPIX
          DO K=1, NFAC 
C            IF IMAGES HAVE SOME CROSS-CORRELATION (I.E; (PART OF) ONE 
C	     IMAGE IS THE SAME AS (PART OF) ANOTHER), ONE FACTOR
C	     WILL BE 0 (SO A NEGATIVE EIGENVALUE MAY EXIST).

             EVALS(K) = MAX(EVALS(K), 1.0E-9)
             BLW(K)   = EVECTS(J,K) * SQRT(EVALS(K))            
          ENDDO

C         WRITE  DATA TO PIX*** FILE 
          FPIX = J
C          PIAT = WEIGHTP(J) / NPIX feb 2010 al
          PIAT = WEIGHTP(J) * NUMIM  / NPIX
          WRITE(LUNP,90) (BLW(K),K=1,NFAC),PIAT, CO(J),FPIX,FDUM
90        FORMAT(10(1PG12.5,' '))
        ENDDO

        END


C++*********************************************************************
C
C PCA_STOCHA         CREATED, AUTHOR       11/01/1993 MAHIEDDINE LADJADJ      
C                    REWRITTEN             03/06/2006 ArDean Leith
C
C **********************************************************************
C
C PCA_STOCHA(JBASE, NPIX, NUMIM, NFAC, S, U, TRACE, BB, 
C            BLD, AD, LSAV, LUNT, WEIGHTP, VARP)
C
C PARAMETERS:
C       JBASE   NFAC  + 3                                      (INPUT)
C       NPIX    NUMBER OF PIXELS PER IMAGE                     (INPUT)
C       NUMIM   NUMBER OF IMAGES                               (INPUT)
C       NFAC    NUMBER OF EIGENVECTORS REQUESTED               (INPUT)
C       S                                                      (     )
C       U                                                      (     )
C       TRACE                                                  (INPUT)
C       BB                                                     (     )
C       BLD                                                    (OUTPUT)
C       AD                                                     (     )
C       LSAV    INPUT UNIT FOR TEMP FILE HOLDING PIXEL VALUES  (INPUT)
C       LUNT    INPUT AND OUTPUT UNIT FOR TEMP FILE USED HERE  (INPUT)
C       WEIGHTP SUM OF PIXEL VALUES (UNDER MASK) BY PIXEL      (INPUT)
C       VARP    VARIANCE OF PIXEL VALUES (UNDER MASK) BY PIXEL (INPUT)
C
C--*********************************************************************

        SUBROUTINE PCA_STOCHA(JBASE,   NPIX, NUMIM, NFAC, 
     &                        S,       U,    TRACE, BB, 
     &                        BLD,     AD,   LSAV,  LUNT,
     &                        WEIGHTP, VARP)

        INCLUDE 'CMBLOCK.INC'

        REAL :: S(NPIX, JBASE), U(NPIX),  
     &          BB(NPIX,JBASE), BLD(NPIX), AD(JBASE,JBASE),
     &          WEIGHTP(NPIX),  VARP(NPIX) 

C       AUTOMATIC ARRAYS
        LOGICAL :: INB(JBASE)

        LOGICAL :: S_ON_DISK, LDUM

C       NAR  =  CENTERING FREQUENCY
        DATA  EPS/1.0E-4/, NAR/10/, KITER/8/

        S_ON_DISK = .FALSE.

C	PUT RANDOM NUMBERS IN  'S' MATRIX.
        DO I=1,NPIX
           DO J=1,JBASE
              S(I,J) = SEN3A(BID)
           ENDDO
        ENDDO

C       ORTHONORMALIZATION OF THE JBASE (FIRST) COLUMNS OF S(NPIX,*)
        CALL PCA_GSMOD(NPIX,NPIX,JBASE,S,KR,VARP)
        
C       INITIALIZE INB ARRAY                    
        INB = .TRUE.

        NIT  = 0
	LDUM = .TRUE.
        DO WHILE (LDUM)
          NIT = NIT + 1
          WRITE(NOUT,*)  ' Iteration: ', NIT  
          CALL PCA_ITPOW(NUMIM,NPIX,JBASE,NAR,WEIGHTP,S,BB,U, 
     &                   VARP,TRACE,LSAV,INB)

          IF (NIT. GE. KITER - 1) THEN
            IF ((NIT .GT. KITER) .AND. (MOD(NIT-12,4) .NE. 0)) THEN
              CONTINUE
            ELSE
              JCASE = 0
              DO L=1,JBASE
                 IF (INB(L)) JCASE = JCASE + 1
              ENDDO

              CALL PCA_CPROJ(NUMIM,JBASE,NPIX,S,VARP,BB,U,AD, 
     &                      WEIGHTP,LSAV,JCASE,INB)
              JCASE = 0
              DO L=1,JBASE
                IF (INB(L)) THEN
                   JCASE  = JCASE + 1
                   BLD(L) = VARP(JCASE)  ! BLD DEFINED OVER 1...JBASE
                ENDIF
              ENDDO

              IF (S_ON_DISK)  THEN
                REWIND LUNT
                AVE = 1.0
                DO  L=1,NFAC
                  VARP(L)  =  0.0                          
                  READ(LUNT) (U(J),J=1,NPIX)

                  DO  J =1,NPIX
                     VARP(L) = VARP(L) + S(J,L) * U(J)
                  ENDDO

                  IF ((1.0 - ABS(VARP(L))) .LT. EPS)  THEN
                     INB(L) = .FALSE.
                  ELSE
                     AVE = MIN(AVE, ABS(VARP(L)))
                  ENDIF
                ENDDO

                WRITE(NOUT,501) AVE
501             FORMAT('  Cosines of eigenvectors  ----- ', /, 
     &	               '       Minimum cosine: ', F8.5)

                WRITE(NOUT,502) (VARP(L),L=1,NFAC)
502             FORMAT(8(2X, F8.5))

                WRITE(NOUT,505) (INB(L),L=1,JBASE)
505             FORMAT(8(9X, L1))
    
                IF ((1.0-AVE) .LT. EPS) EXIT
              ENDIF
              REWIND LUNT

              DO L =1,JBASE
                 WRITE(LUNT) (S(J,L),J=1,NPIX)
              ENDDO

              S_ON_DISK = .TRUE.
            ENDIF
          ENDIF
        ENDDO

        END


C **********************************************************************
C
C       AUTHOR: MAHIEDDINE LADJADJ          11/1/93
C
C--*********************************************************************

        SUBROUTINE PCA_ITPOW(NUMIM, NPIX, JBASE, NAR, WEIGHTP, 
     &                       S, BB, U, V, TRACE, LSAV, INB)
                                                                               
        REAL     :: WEIGHTP(NPIX), S(NPIX,JBASE), BB(NPIX,JBASE)
        REAL     :: U(NPIX),       V(NPIX)

        LOGICAL  :: INB(JBASE)
        LOGICAL  :: CENTER
                                                                               
        CALL REW(LSAV,1)
        DO L =1,JBASE
          IF (INB(L))  THEN
            DO  J = 1,NPIX
               BB(J, L) =  S(J, L)
               S(J, L)  =  0.0
            ENDDO
          ENDIF
        ENDDO
                                                                               
        DO  IA=1,NUMIM
            READ(LSAV) (U(K),K=1,NPIX),FIM
            CENTER = MOD(IA,NAR) .EQ. 0 .OR. IA .EQ. NUMIM
            DO L = 1,JBASE
              IF (INB(L)) THEN
                T1  = 0.0                                                
                DO K=1,NPIX
                   T1 = T1 + BB(K, L) * (U(K) - WEIGHTP(K))
                ENDDO
                                                                              
                DO K=1,NPIX
                   S(K,L) = S(K,L) + (U(K) - WEIGHTP(K)) * T1
                ENDDO

                IF (CENTER)  THEN
C                 PERIODIC CENTERING
                  T2  = 0.0                                                
                  DO JP=1,NPIX
                     T2 = T2 + S(JP,L)
                  ENDDO

                  T2 = T2 / NPIX
                  DO J=1,NPIX
                     S(J,L) = S(J,L) - T2                            
                  ENDDO
                ENDIF
              ENDIF
            ENDDO
        ENDDO

        DOP  =  TRACE/2/(NPIX - 1)
                                                                             
        DO  L=1,JBASE
          IF (INB(L)) THEN
            DO  J=1,NPIX
              S(J,L) =  S(J,L) - DOP * BB(J,L)
            ENDDO
          ENDIF
        ENDDO

        DO L=1,JBASE
          IF (.NOT.INB(L))  THEN
             CALL PCA_GSMODL(NPIX,JBASE,S,KRANG,V,INB)
             RETURN
          ENDIF
        ENDDO

        CALL PCA_GSMOD(NPIX,NPIX,JBASE,S,KRANG,V)

        END            


C **********************************************************************
C
C     AUTHOR: MAHIEDDINE LADJADJ          11/1/93
C
C     ORTHONORMALIZATION OF THE  JCARD (FIRST) COLUMNS  OF                      
C     S(ICARD,*) BY THE METHOD GRAM - SCHMIDT (MODIFIED) .
C     INPUT         1/ IDIM   RESERVED DIMENSION FOR S(IDIM,*)
C                   2/ ICARD  ACTUAL NUMBER OF ROWS FOR S(ICARD,*)
C                   3/ JCARD  NBR. OF COLUMNS TO PROCESS S(ICARD, JCARD)
C                   4/ P(*)   WEIGHT VECTOR , DIMENSION P(IDIM) = 1
C     INPUT - OUTPUT  5/ S(*,*) INPUT  THE MATRIX TO BE PROCESSED
C                             OUTPUT ORTHONORMALIZED MATRIX (P METRIC)          
C     OUTPUT        6/ KRANG  RANK OF MATRIX S(ICARD, JCARD)  
C                   7, 8/ T(*), V(*) WORKING ARRAYS T(IDIM), V(IDIM)
C     IF THERE IS COLINEARITY,  THE CORRESPONDING COLUMN IS SET TO ZERO 
C
C--*********************************************************************
        
        SUBROUTINE PCA_GSMOD(IDIM,ICARD,JCARD,S,KRANG,V)

        INCLUDE 'CMBLOCK.INC'

        DIMENSION  S(IDIM,JCARD), V(IDIM)

        DATA  EPS / 1.0E-10 /

        KRANG  =  JCARD 

C       INITIAL NORMS. ORTHONORMALIZATION OF S(* , 1)                           
        DO J=1,JCARD 
           V(J) = 0.0                          
           DO I=1,ICARD 
              V(J) = V(J) + S(I,J) * S(I,J)
           ENDDO
           V(J) = MAX(V(J),EPS)
        ENDDO

        C  = 1.0 / SQRT(V(1))                                                 
        DO I=1,ICARD                                                      
           S(I,1) =  C * S(I,1)                     
        ENDDO
        IF (JCARD .EQ. 1)  RETURN

C       ORTHOGONALISATION OF S(*, J1). MODIFICATION OF THE NEXT.               
        DO J=1,JCARD-1
          J1 = J + 1                      
          DO JJ=J1,JCARD                
            TJJ  = 0.0   
                     
            DO I=1,ICARD                
               TJJ = TJJ +  S(I,JJ) * S(I,J)
            ENDDO

            DO I=1,ICARD             
               S(I,JJ) = S(I,JJ) - TJJ * S(I,J)
            ENDDO
          ENDDO

C         COLINEARITY TEST. NORMALIZATION OF S(*,J1)
          C  = 0.0 
          DO I=1,ICARD             
             C = C + S(I,J1) * S(I,J1)
          ENDDO

          IF (C/V(J1) .LE. EPS) THEN
            KRANG = KRANG - 1

            DO I=1,ICARD
               S(I,J1) = 0.0
            ENDDO
          ELSE
            C = 1.0 / SQRT(C)

            DO I=1,ICARD
               S(I,J1) = C * S(I,J1)
            ENDDO
          ENDIF
        ENDDO

       END


C **********************************************************************
C
C     AUTHOR: MAHIEDDINE LADJADJ          11/1/93
C
C     ORTHONORMALIZATION OF THE CHOSEN COLUMNS  OF S(IDIM,*)
C     INPUT         1/ IDIM   RESERVED DIMENSION FOR S(IDIM,*)
C                   3/ JCARD  NBR. OF COLUMNS TO PROCESS S(IDIM, JCARD)         
C                   4/ P(*)   WEIGHT VECTOR , DIMENSION P(IDIM) = 1
C     INPUT - OUTPUT  5/ S(*,*) INPUT  THE MATRIX TO BE PROCESSED
C                             OUTPUT ORTHONORMALIZED MATRIX (P METRIC)
C     IF INB(J) = .TRUE.  THEN J'TH COLUMN HAS TO BE PROCESSED
C     IF INB(J) = .FALSE. THEN J'TH COLUMN IS ALREADY ORTHONORMALIZED
C
C--*********************************************************************

        SUBROUTINE PCA_GSMODL(IDIM,JCARD,S,KRANG,V,INB)

        INCLUDE 'CMBLOCK.INC'

        DIMENSION  S(IDIM, JCARD)
        LOGICAL  INB(JCARD), V(JCARD)

        DATA  EPS / 1.0E-10 /

        V(1) = .TRUE.
        DO J=1,JCARD
           IF (.NOT. INB(J)) V(1) = .FALSE.
        ENDDO

        IF (V(1)) THEN
           CALL PCA_GSMOD(IDIM, IDIM, JBASE, X, KRANG, V)
           RETURN
        ENDIF

        KRANG  = JCARD
        DO  J=1,JCARD
          V(J) = .NOT.INB(J)
        ENDDO
        DO  J=1,JCARD
          IF (INB(J))  THEN
            DO JJ=1,JCARD
              IF (V(JJ)) THEN
                TJJ =  0.0
                DO I=1,IDIM
                   TJJ = TJJ + S(I,JJ) * S(I, J)
                ENDDO

                DO I=1,IDIM
                   S(I,J)= S(I,J) - TJJ * S(I,JJ)
                ENDDO

                VV = 0.0
                DO I=1,IDIM
                  VV =  VV  +  S(I,J) * S(I,J)
                ENDDO
                VV  = MAX(VV, EPS)
                C  =  1.0 / SQRT(VV)

                DO I=1,IDIM
                   S(I,J) = C * S(I,J)
                ENDDO
              ENDIF
            ENDDO
            V(J) = .TRUE.
          ENDIF
        ENDDO

        END

C **********************************************************************
C
C	AUTHOR: MAHIEDDINE LADJADJ          11/1/93
C
C	OPERATION OF PROJECTION AND DIAGONALIZATION, 
C       FOR DIAGONALIZATION BY DIRECT READING.
C        
C--*********************************************************************

        SUBROUTINE PCA_CPROJ(NUMIM,  JBASE, NPIX,  S, D, BB, U, AD, 
     &                       WEIGHTP, LSAV, JCASE, INB)

        DIMENSION S(NPIX,JBASE), BB(NPIX,JBASE), AD(JCASE,JCASE), 
     &            U(NPIX),       WEIGHTP(NPIX),  D(NPIX)
        LOGICAL  INB(JBASE)

        CALL REW(LSAV,1)

C       ZERO AD ARRAY
        AD = 0.0

        DO  I=1,NUMIM                                     
            READ(LSAV) (U(KK),KK=1,NPIX),FIM
            MC = 0
            DO M=1,JBASE
              IF (INB(M))  THEN
                MC = MC + 1
                LC = 0
                DO L=1,M
                  IF (INB(L))  THEN
                    LC = LC + 1
                    CIL  = 0.0
                    CIM  = 0.0
                    DO K=1,NPIX        
                       UUU = U(K) -  WEIGHTP(K)
                       CIL = CIL  +  S(K,L) * UUU
                       CIM = CIM  +  S(K,M) * UUU
                    ENDDO
                    AD(LC,MC) = AD(LC,MC) + CIL * CIM
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
        ENDDO
                                                                               
        DO M=1,JCASE                                                       
          DO L=1, M                                                   
            AD(M,L) = AD(L,M)
          ENDDO
        ENDDO

        CALL VPROP(JCASE, JCASE, AD,  D,  U, KOD)

        MC = 0
        DO M=1,JBASE
          IF (INB(M))  THEN
            MC = MC + 1
            DO J=1, NPIX
              BB(J,M) = 0.0
              KC = 0
              DO  K=1,JBASE
                IF (INB(K)) THEN
                  KC      = KC + 1
                  BB(J,M) = BB(J,M) +  S(J,K) * AD(KC,MC)
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO

        DO L=1,JBASE
          IF (INB(L)) THEN
             DO J=1,NPIX
                S(J,L) = BB(J,L)
             ENDDO
          ENDIF
        ENDDO

        END
