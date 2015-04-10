C **********************************************************************
C
C  JPMSK3   NEW                                FEB 04 ARDEAN LEITH
C           PIXEL FILE FDUM                    JUN 09 ARDEAN LEITH
C           IMC FILE FDUM                      JUN 09 ARDEAN LEITH
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
C   JPMSK3()
C
C   PURPOSE:  PLOT BOTH ACTIVE AND INACTIVE IMAGE LOCATIONS 
C
C  OPERATIONS SUPPORTED : CA SMI  -  RUN CORAN
C
C  CALL TREE:
C     JPMSK3 ---> JPMSK4 --->  INCOR3
C                             | GETCOO 
C                             | GETCOOT (IF TRANSPOSED)
C
C  VARIABLES:           
C	NFAC 	  NUMBER OF EIGENVECTORS REQUESTED              (INPUT)
C	NPIX      NUMBER OF PIXELS PER IMAGE                    (INPUT)
C       NUMIM	  NUMBER OF IMAGES                              (INPUT)
C       INUMBR()  IMAGE NUMBER LIST                             (INPUT)
C	USE_PCA   CORAN VS PCA FLAG                             (INPUT)
C       EVECTS()  EIGENVECTORS (COLUMN)  OF X'X AND             (INPUT)
C                 X(I,*)= BLU() J=1,JTOT W/ I=1,ITOT
C       EVALS()   EIGENVALUE ARRAY                           (INPUT/OUTPUT)
C	WEIGHTP() SUM OF PIXEL VALUES AT THIS PIXEL             (INPUT)
C       CO()      WORKING ARRAY  
C	SUMW      SUM OF ALL THE PIXEL VALUES IN ALL IMAGES     (INPUT)
C       BLU()     WORKING ARRAY FOR INPUTS
C       BLW()     WORKING ARRAY FOR OUTPUTS
C	LUNS	  SEQUENTIAL IMAGE I/O UNIT ( FOR INPUT)        (INPUT)
C	LUNI      IMAGE COORDINATE I/O UNIT (FOR OUTPUT)        (INPUT)
C	LUNP      PIXEL COORDINATE I/O UNIT (FOR OUTPUT)        (INPUT)
C
C	S(,) HAS THE EIGENVECTORS (COLUMN)  OF X'X AND D() HAS THE
C       EIGENVALUES. X(I,*)= U() J=1,JTOT W/ I=1,ITOT
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C-----------------------------------------------------------------------------

        SUBROUTINE JPMSK3(LUNDOC,LUNM,LUNP,LUNE,LUNIN)
        
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        
        COMMON /IOBUF/ BUF(NBUFSIZ)

        CHARACTER(LEN=MAXNAM) :: FILPATIA,FILPATII,FILPATI,FILNMI,
     &                           FILPRE,FILNMP,FILNME,FILNMM
        COMMON /COMMUN/          FILPATIA,FILPATII,FILPATI,FILNMI,
     &                           FILPRE,FILNMP,FILNME,FILNMM
        
        REAL, ALLOCATABLE, DIMENSION(:)     :: BUFM,BUFI,BUFP,WEIGHTP
        REAL, ALLOCATABLE, DIMENSION(:,:)   :: EVECTS
        REAL, ALLOCATABLE, DIMENSION(:,:)   :: BLW

        CHARACTER(LEN=1)                    :: NULL, ANS
        CHARACTER(LEN=3)                    :: ANS1
        
        LOGICAL                             :: USE_PCA,ADANEG

        INTEGER, ALLOCATABLE                :: INUMBRI(:),INUMBRA(:)

        NULL = CHAR(0)
        LUNI = 10

        ALLOCATE(INUMBRI(NIMAX),
     &           INUMBRA(NIMAX),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'JPMSK3; INUMBRI....',2*NIMAX)
           RETURN
        ENDIF

        CALL FILERD(FILPRE,NLET,NULL,
     &              'CORAN/PCA FILE PREFIX (e.g. CORAN_)~',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       READ IN SELECTED ACTIVE IMAGE FILE NUMBERS; STORE IN INUMBRA().
        NUMIMA  =  NIMAX 
        CALL FILELIST(.TRUE.,LUNDOC,FILPATIA,NLET1,
     &         INUMBRA,NIMAX,NUMIMA,'ACTIVE IMAGE FILE TEMPLATE',
     &         IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        WRITE(NOUT, *) ' NUMBER OF ACTIVE IMAGES: ',NUMIMA,CHAR(12)

C       READ IN SELECTED IMAGE FILE NUMBERS; STORE IN INUMBRI().
        NUMIMI  =  NIMAX 
        CALL FILELIST(.TRUE.,LUNDOC,FILPATII,NLET1,
     &         INUMBRI,NIMAX,NUMIMI,
     &         'INACTIVE IMAGE FILE TEMPLATE',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        WRITE(NOUT, *) ' NUMBER OF INACTIVE IMAGES: ',NUMIMI

        FILNMP = FILPRE(1:NLET) // '_PIX'//NULL
        FILNME = FILPRE(1:NLET) // '_EIG'//NULL
        FILNMM = FILPRE(1:NLET) // '_MAS'//NULL

        NUMIM  = NUMIMA + NUMIMI

C       OPEN MASK FILE
        CALL OPFILEC(0,.FALSE.,FILNMM,LUNM,'O',IFORM,
     &              NSAMM,NROWM,NSLICEM, MAXIMT,' ',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
 
 
C       ALLOCATE MASK & IMAGE INPUT BUFFER
        NPIX = NSAMM*NROWM*NSLICEM 
        ALLOCATE (BUFM(NPIX),BUFI(NPIX),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'JPMSK2; BUFM & BUFI',2*NPIX)
           RETURN
        ENDIF

C       READ MASK IN BUFM
        CALL REDVOL(LUNM,NSAMM,NROWM,1,NSLICEM,BUFM,IRTFLG)
        NMASK = 0
        DO I=1,NSAMM*NROWM*NSLICEM
           IF (BUFM(I) .GT. 0.5) NMASK = NMASK + 1
        ENDDO

        WRITE(NOUT,*)' NUMBER OF PIXELS UNDER MASK: ',NMASK
        IF (NMASK .LE. 0) THEN
           CALL ERRT(101,'NO PIXELS UNDER MASK',NE)
           GOTO 9999
        ENDIF

C       OPEN EIGENVALUE FILE
        CALL OPAUXFILE(.FALSE.,FILNME,DATEXC,LUNE,0,
     &                       'O', ' ',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       READ _EIG FILE HEADER
	READ(LUNE,*) NFAC, SUMP, TRACE, KIND_PCA, NVAL
80      FORMAT(I10,G12.5,1X,G12.5,I10,I10)

        USE_PCA = (KIND_PCA .EQ. 1)

C       READ _EIG EIGEN VALUES
        DO I = 1,NFAC
	   READ(LUNE,*) EVALSDUM, PERDUM, CULDUM
        ENDDO

        WRITE(NOUT,*)' NUMBER OF FACTORS: ',NFAC

C       ALLOCATE EVECTS INPUT AND BLW WORKING BUFFER
        ALLOCATE (EVECTS(NVAL,NVAL), 
     &            BLW(NFAC,NUMIM), 
     &            STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MWANT = (NVAL * NVAL) + (NFAC * NUMIM)
           CALL ERRT(46,'JPMSK3; EVECTS & BLW',MWANT)
           GOTO 9999
        ENDIF

C       READ _EIG EIGENVECTORS 
	DO I = 1, NVAL
           READ(LUNE,*) (EVECTS(I,J),J=1,NVAL)
C          write(6,*)   (EVECTS(I,J),J=1,NVAL)
        ENDDO
      
C       ALLOCATE WEIGHTP BUFFER, END MAY BE UNUSED IF MASKED
        ALLOCATE (WEIGHTP(NPIX),
     &            STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'JPMSK3; WEIGHTP....',NPIX)
           GOTO 9999
        ENDIF

        IF (USE_PCA) THEN
C          OPEN PIXEL COORDINATE FILE TO GET WEIGHTP.
           CALL OPAUXFILE(.FALSE.,FILNMP,DATEXC,LUNP,0,
     &                       'O', ' ',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
 	   READ(LUNP,*) NMASKPT 

C          ALLOCATE PIXEL INPUT BUFFER
           ALLOCATE (BUFP(NMASKPT),
     &               STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'JPMSK3; BUFP....',NMASKPT)
              GOTO 9999
           ENDIF

C          READ THE _PIX  FILE TO GET WEIGHTP.
           DO I=1,NPIX
              READ(LUNP,*) BUFP,WEIGHTP(I),FDUMCO,FDUMPIX,FDUM 
           ENDDO
        ELSE
           WEIGHTP = 0.0
        ENDIF

        CALL JPMSK4(LUNIN,NUMIMA,NPIX,NFAC,NSAMM,NROWM,NSLICEM,
     &                    BUFM,BUFI,FILPATIA,INUMBRA,
     &                    WEIGHTP,USE_PCA,BLW,EVECTS,NVAL)

        CALL JPMSK4(LUNIN,NUMIMI,NPIX,NFAC,NSAMM,NROWM,NSLICEM,
     &                    BUFM,BUFI,FILPATII,INUMBRI,
     &                    WEIGHTP,USE_PCA,BLW(1,NUMIMA+1),EVECTS,NVAL)

C       OPEN FORMATTED IMAGE COORDINATE FILE (_IMC)
        CALL FILERD(FILNMI,NLET,NULL,'COORDINATE OUTPUT FILE PREFIX~',
     &                IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        FILPATI = FILNMI(1:NLET) // '_IMC'//NULL
        CALL OPAUXFILE(.FALSE.,FILPATI,DATEXC,LUNI,0,
     &                       'U', ' ',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       CREATE IMAGE COORDINATE FILE HEADER 
        WRITE(LUNI,95) NUMIM, NFAC, NSAMM, NROWM, NUMIM, KIND_PCA
95      FORMAT(10I10)

        FDUM = 0.0
        ACTA = 1.0
        ACTI = 0.0

        DO I = 1, NUMIMA
C         WRITE IMC DATA TO IMC FILE 
          FIM = INUMBRA(I)
          IF (USE_PCA) THEN
             WRITE(LUNI,90) (BLW(K,I),K=1,NFAC),FDUM,FDUM,FIM,ACTA
90           FORMAT(10(1PG12.5,' '))
          ELSE
             WRITE(LUNI,90) (BLW(K,I),K=1,NFAC),FDUM,FDUM,FIM,ACTA
          ENDIF
        ENDDO

        DO I = 1, NUMIMI
C         WRITE IMC DATA TO IMC FILE 
          FIM = INUMBRI(I)
          IF (USE_PCA) THEN
            WRITE(LUNI,90)(BLW(K,I+NUMIMA),K=1,NFAC),FDUM,FDUM,FIM,ACTI
          ELSE
            WRITE(LUNI,90)(BLW(K,I+NUMIMA),K=1,NFAC),FDUM,FDUM,FIM,ACTI
          ENDIF
        ENDDO

9999    IF (ALLOCATED(BUFM))    DEALLOCATE(BUFM)
        IF (ALLOCATED(BUFI))    DEALLOCATE(BUFI)
        IF (ALLOCATED(BUFP))    DEALLOCATE(BUFP)
        IF (ALLOCATED(WEIGHTP)) DEALLOCATE(WEIGHTP)
        IF (ALLOCATED(INUMBRI)) DEALLOCATE(INUMBRI)
        IF (ALLOCATED(INUMBRA)) DEALLOCATE(INUMBRA)


C       CLOSE ALL FILES THAT MIGHT BE OPEN
        CLOSE(LUNDOC)
        CLOSE(LUNP)
        CLOSE(LUNM)
        CLOSE(LUNE)
        CLOSE(LUNI)

        END

C       ------------------------- JPMSK4 ---------------------------


        SUBROUTINE JPMSK4(LUNIN,NUMIM,NPIX,NFAC,NSAMM,NROWM,NSLICEM,
     &                    BUFM,BUFI,FILPAT,ILIST,
     &                    WEIGHTP,USE_PCA,BLW,EVECTS,NVAL)
        
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        
        CHARACTER(LEN=*)                       :: FILPAT
        
        INTEGER,DIMENSION(NUMIM)               :: ILIST

        REAL, DIMENSION(NPIX)                  :: BUFM,BUFI
        REAL, DIMENSION(NPIX)                  :: WEIGHTP
        REAL, DIMENSION(NVAL,NVAL)             :: EVECTS
        REAL, DIMENSION(NFAC,NUMIM)            :: BLW
        REAL, DIMENSION(NSAMM)                 :: BUF

        LOGICAL                                :: USE_PCA,ADANEG
        LOGICAL                                :: TRANSPOSE

        REAL, DIMENSION(NUMIM)                 :: PDUM
        CHARACTER(LEN=MAXNAM)                  :: FILNAM

        ADANEG    = .FALSE.   ! NOT IMPLEMENTED YET
        TRANSPOSE = .FALSE.   ! NOT IMPLEMENTED YET

C       LOOP OVER ALL IMAGES
        DO IM = 1, NUMIM
              NLET  = 0
              CALL FILGET(FILPAT,FILNAM,NLET,ILIST(IM),IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999

              MAXIMT = 0
              CALL OPFILEC(0,.FALSE.,FILNAM,LUNIN,'O',IFORM,
     &                   NSAM,NROW,NSLICE,MAXIMT,' ',.FALSE.,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO  9999 

              IF (FMIN .LT. 0.0 .AND. .NOT. USE_PCA)  THEN
                 CALL ERRT(101,'CORAN DOES NOT ACCEPT NEGATIVE PIXELS',
     &                     NDUM)
                 GOTO 9999
               ENDIF

C             COMPARE MASK DIMENSIONS WITH IMAGE DIMENSIONS,ETC
              IF (IFORM .NE. 1 .AND. IFORM .NE. 3) THEN
                  CALL ERRT(2, 'JPMSK3', NDUM)
                  GOTO  9999

              ELSEIF ((NSAMM   .NE. NSAM)   .OR. 
     &                (NROWM   .NE. NROW)   .OR.
     &                (NSLICEM .NE. NSLICE)) THEN
                 WRITE(NOUT, 93) NSAM, NROW, NSAMM, NROWM, 
     &                            NSLICE,NSLICEM
93               FORMAT('*** IMAGE DIMENSION (',I4,',',I4,',',I4,
     &                  ')  NOT SAME AS MASK (',I4,',',I4,',',I4,')')
                 CALL ERRT(100,'JPMSK3',NE)
                 GOTO 9999
              ENDIF
    
C             READ IMAGE AREA WHERE MASK > 0.5 INTO CORE AND
C	      COMPUTE ITS MASK-RELATED,  PRECISE AVERAGE
              WEIGHTIT = 0.0
              ILOCM    = 0
              IPIX     = 0

              DO IROW = 1,NROW*NSLICE
                 CALL REDLIN(LUNIN,BUF,NSAM,IROW)
                 DO ISAM = 1,NSAM
                    ILOCM = ILOCM + 1
                    IF (BUFM(ILOCM) .GT. 0.5) THEN
C                      INSIDE MASK, USE THIS PIXEL
                       VAL = BUF(ISAM)
                       IF (ADANEG) VAL = VAL + ADA

C	               WEIGHTIT = SUM OF THE ELEMENTS (MASK > 0.5)    .
                       WEIGHTIT      = WEIGHTIT + VAL               

                       IPIX          = IPIX + 1
                       BUFI(IPIX)    = VAL 
                    ENDIF
                 ENDDO
              ENDDO
c           write(6,*) 'ipix:',ipix,npix,use_pca

C          WEIGHTIT = SUM OF ALL PIXEL VALUES UNDER MASK IN IMAGE.

C          FIND IMC COORDINATES
           IF (TRANSPOSE) THEN
C             NOT IMPLEMENTED YET!!! NEED EIG. CALCULATION
           ELSE
              DO K=1,NFAC                                                   
                 BLW(K,IM) = 0.0 
                 IF (USE_PCA) THEN
                    DO J=1,IPIX
                      BLW(K,IM) = BLW(K,IM) + 
     &                      (BUFI(J) - WEIGHTP(J)/NUMIM) * EVECTS(J, K)
                    ENDDO
                 ELSE
                    DO J=1,IPIX
                      BLW(K,IM) = BLW(K,IM) + 
     &                            (BUFI(J) * EVECTS(J, K)) / WEIGHTIT
                    ENDDO
                 ENDIF
              ENDDO
           ENDIF

C       END OF:  DO IM = 1, NUMIM
        ENDDO


9999    CONTINUE
        END



C       ------------------------- JPMSK5 ---------------------------


        SUBROUTINE JPMSK5(LUNIN,NUMIM,NFAC, FILPAT,
     &                    WEIGHTI,CO,ID,USE_PCA,IRTFLG)
        
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        
        CHARACTER(LEN=*)                       :: FILPAT
        
        REAL, DIMENSION(NUMIM)                 :: WEIGHTI,CO
        INTEGER, DIMENSION(NUMIM)              :: ID
        REAL, DIMENSION(NFAC)                  :: DUM

        LOGICAL                                :: USE_PCA

        CHARACTER(LEN=MAXNAM)                  :: FILNAM

C       OPEN _IMC COORDINATE FILE TO GET WEIGHTI, CO.
        CALL OPAUXFILE(.FALSE.,FILPAT,DATEXC,LUNIN,0,
     &                       'O', ' ',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       READ IMAGE COORDINATE FILE HEADER 
        READ(LUNIN,*) NUMIMT, NFACT, NSAMMT, NROWMT, NUMIMT,KIND_PCAT

        DO I=1,NUMIM

C         READ THE _IMC  FILE TO GET WEIGHTI....
          READ(LUNIN,*) (DUM(K),K=1,NFAC),WEIGHTI(I),CO(I),FPIX,FDUM
          ID(I) = FPIX
        ENDDO

        CLOSE(LUNIN)
        IRTFLG = 0

        END




