C **********************************************************************
C
C  JPMSK1          RANDOM SHUFFLING              JULY 87 
C                  NEW IN - CORE                 12/30/87 
C                  CO - PROJ. OF ARBITRARY DATA  5/03/89 JF
C                  LONG FILE NAMES               FEB 89  ArDean Leith
C                  CRAY COMPATIBLE FORTRAN       JAN 90  JF	
C                  MODIFIED                     9/1/93   M. LADJADJ   
C                  FORCE_INCORE ADDED           JUN 2000 ArDean Leith
C                  USED OPAUXFILE                 4/2/01 ArDean Leith
C                  INCREASED MASK SIZE           3/12/02 ArDean Leith
C                  ALLOWED VOLUME INPUT          8/23/02 ArDean Leith
C                  INCREASED MAXIM              10/22/02 ArDean Leith
C                  OPFILEC                        FEB 03 ARDEAN LEITH
C                  REWRITTEN                      SEP 03 ARDEAN LEITH
C                  _MAS FILE                      DEC 03 ARDEAN LEITH
C                  CLOSE LUN                      APR 04 ARDEAN LEITH
C                  MAX. NFAC LISTED               MAY 04 ARDEAN LEITH
C                  SKIPREST                       MAY 04 ARDEAN LEITH
C                  ITERATIVE PCA RECOVERED        MAR 06 ARDEAN LEITH
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
C   JPMSK1()
C
C   PURPOSE:  RUN CORAN & PCA PROGRAMS
C             USE MASK FILE TO DECIDE WHICH POINTS TO INCLUDE
C
C   OPERATIONS SUPPORTED : 'CA S'
C
C   CALL TREE:
C     JPMSK1 ---> SCORAN3 --->  INCOR3
C                             | GETCOO 
C                             | GETCOOT (IF TRANSPOSED)
C
C
C  JPMSK1
C    WEIGHTI = SUM OF PIXEL VALUES (UNDER MASK) BY IMAGE().
C    WEIGHTP = SUM OF PIXEL VALUES (UNDER MASK) BY PIXEL().
C    IN LUNS : LIST OF ALL PIXELS  (UNDER MASK) BY IMAGE
C    IN LUNT : TRANSPOSE OF LUNS
C
C  SCORAN3   
C    NATIVE
C       INCOR3
C          PIA = WEIGHTI(IMAGE)
C          if PCA:    MATS(J,JJ) = MATS(J,JJ) + (BLU(J) * BLU(JJ))
C          if CORAN:  MATS(J,JJ) = MATS(J,JJ) + (BLU(J) * BLU(JJ)) /PIA 
C
C          if PCA:    MATS(J,JJ) = (MATS(J, JJ) - 
C                             (WEIGHTP(J) * WEIGHTP(JJ)) / NUMIM)
C          if CORAN:  AAA        =  SQRT(WEIGHTP(J) * WEIGHTP(JJ))
C                     MATS(J,JJ) =  MATS(J, JJ) / AAA  -  AAA / SUMW
C
C          if CORAN:  MATS(K,L) = MATS(K,L) * SQRT(SUMW / WEIGHTP(K)) 
C       GETCOO
C          if PCA:    BLW(K) = BLW(K) + (BLU(J) - WEIGHTP(J) / NUMIM) * 
C                                    EVECTS(J, K) 
C          if CORAN:  BLW(K) = BLW(K) + (BLU(J) * EVECTS(J, K)) / PIA
C          IN LUNP :  PIAT   = WEIGHTP(J) / NPIX
C    TRANSPOSED
C       INCORT
C       GETCOOT
C
C    IN LUNE : 
C       LIST OF ALL EIGENVALUES BY FACTOR
C       if NATIVE:      LIST OF NPIX EIGENVECTORS 
C       if TRANSPOSED:  LIST OF NUMIM EIGENVECTORS
C
C
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C-----------------------------------------------------------------------------

        SUBROUTINE JPMSK1()
        
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        
        COMMON /IOBUF/ BUF(NBUFSIZ)

        CHARACTER(LEN=MAXNAM) :: FILPATI,FILNAM,FILPATP,FILPATE,FILPATT,
     &                           FILPATC,MSKNAM,FILPATS,FILPRE,FILPATM

        COMMON /COMMUN/          FILPATI,FILNAM,FILPATP,FILPATE,FILPATT,
     &                           FILPATC,MSKNAM,FILPATS,FILPRE,FILPATM
        
        REAL, ALLOCATABLE, DIMENSION(:,:)           :: TRANS
        REAL, ALLOCATABLE, DIMENSION(:)             :: BUFM,BUFI,BUFT
        REAL, ALLOCATABLE, DIMENSION(:)             :: WEIGHTI,WEIGHTP
        REAL, ALLOCATABLE, DIMENSION(:)             :: VARP

        CHARACTER(LEN=1) ::    NULL, ANS
        CHARACTER(LEN=2) ::    ANS1
        
        DOUBLE PRECISION :: DAV, VARIT
        LOGICAL          :: EX, FLAGR,USE_PCA,TRANSPOSE,ADANEG,SKIPREST
        LOGICAL          :: USE_ITERPCA

#ifndef SP_32
        INTEGER *8       SIZE1,SIZE2
#else
        INTEGER *4       SIZE1,SIZE2
#endif

        DATA LUN,LUNS,LUNI,LUNM,LUNP,LUNE,LUNDOC/70,71,72,73,74,75,76/
        DATA LUNT/77/

        NULL    = CHAR(0)

C         READ IN SELECTED IMAGE FILE NUMBERS; STORE IN INUMBR().
          NUMIM  =  NIMAX 
          CALL FILELIST(.TRUE.,LUNDOC,FILPATI,NLET1,
     &           INUMBR,NIMAX,NUMIM,'IMAGE FILE TEMPLATE',IRTFLG)
          IF (IRTFLG .NE. 0) RETURN
          WRITE(NOUT, *) ' NUMBER OF IMAGES: ',NUMIM

          MAXIMT = 0
          CALL FILERD(MSKNAM,NLET,NULL,'MASK',IRTFLG)
          IF (IRTFLG  .NE.  0) THEN
C             OPEN FIRST IMAGE TO GET SIZE
              NLET = 0
              CALL FILGET(FILPATI,FILNAM,NLET,INUMBR(1),IRTFLG)
              IF (IRTFLG .NE. 0) GOTO  9999
              CALL OPFILEC(0,.FALSE.,FILNAM,LUNM,'O',IFORM,
     &                   NSAMM,NROWM,NSLICEM,MAXIMT,' ',.FALSE.,IRTFLG)
             IF (IRTFLG .NE. 0) RETURN 

             NPIX = NSAMM * NROWM * NSLICEM          
             ALLOCATE (BUFM(NPIX),STAT=IRTFLG)
             IF (IRTFLG .NE. 0) THEN
                MWANT = NSAMM * NROWM * NSLICEM 
                CALL ERRT(46,'BUFM',MWANT)
                GOTO 9999
             ENDIF

C            FILL MASK
             BUFM = 1.0
         ELSE
C            OPEN MASK TO GET SIZE
             CALL OPFILEC(0,.FALSE.,MSKNAM,LUNM,'O',IFORM,
     &                 NSAMM,NROWM,NSLICEM, MAXIMT,'  ',.FALSE.,IRTFLG)
             IF (IRTFLG .NE. 0) GOTO  9999 

             NPIX = NSAMM * NROWM * NSLICEM          
             ALLOCATE (BUFM(NPIX),STAT=IRTFLG)
             IF (IRTFLG .NE. 0) THEN
                MWANT = NSAMM * NROWM * NSLICEM 
                CALL ERRT(46,'BUFM',MWANT)
                GOTO 9999
             ENDIF

C            READ IN MASK IN BUFM
             NPIX = 0
             ILOCM = 0
             DO I=1,NROWM*NSLICEM
                CALL REDLIN(LUNM,BUF,NSAMM,I)
                DO J = 1,NSAMM
                   ILOCM = ILOCM + 1
                   IF (BUF(J) .GT. 0.5) THEN
                      BUFM(ILOCM) = 1.0
                      NPIX       = NPIX + 1
                   ELSE
                      BUFM(ILOCM) = 0.0
                   ENDIF
                ENDDO
             ENDDO
             IF (IRTFLG .GT. 0) GOTO 9999

             IF (NPIX .LE. 0) THEN
                CALL ERRT(101,'NO PIXELS UNDER MASK',IDUM)
                GOTO 9999
             ENDIF
             WRITE(NOUT,*)' NUMBER OF PIXELS UNDER MASK: ',NPIX
          ENDIF

C         WEIGHTP = SUM OF PIXEL VALUES AT THIS PIXEL
          ALLOCATE(WEIGHTI(NUMIM),WEIGHTP(NPIX),STAT=IRTFLG)
          IF (IRTFLG .NE. 0) THEN
             MWANT = NPIX + NUMIM 
             CALL ERRT(46,'JPMSK1, WEIGHTI....',MWANT)
             RETURN
          ENDIF

C         ZERO WEIGHTP ARRAY
          WEIGHTP = 0.0

          FLAGR = .FALSE.
C         RANDOM PERMUTATIONS APPEARS NOT TO DO ANYTHING (LOST SOMETIME) 

          NFACMAXTP =  MIN(NUMIM,NPIX) 
          NFACMAXTC =  NFACMAXTP - 1
          WRITE(NOUT,97) NFACMAXTP,NFACMAXTC
97        FORMAT('  MAX NUMBER OF FACTORS USING PCA: ',I7,
     &           ',  USING CORAN: ',I7)

          NFAC = MIN(NFACMAXTP,NFACMAXTC,10)
          CALL RDPRI1S(NFAC,NOT_USED,'NUMBER OF FACTORS',IRTFLG)
          IF (IRTFLG .NE. 0 .OR. NFAC .LE. 0) THEN
             CALL ERRT(102,'INVALID NUMBER OF FACTORS',NFAC)
             RETURN
          ENDIF

          CALL RDPRMC(ANS1,NCHAR,.TRUE.,
     &       'CORAN, PCA, ITERATIVE PCA, OR SKIP ANALYSIS (C/P/I/S)',
     ^       NULL,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN
          USE_PCA     = (ANS1(1:1) .EQ. 'P' )

          USE_ITERPCA = (ANS1(1:1) .EQ. 'I')
          IF (USE_ITERPCA) THEN
             USE_PCA = .TRUE.

C            VARP = VARIANCE OF PIXEL VALUES AT THIS PIXEL
             ALLOCATE(VARP(NPIX),STAT=IRTFLG)
             IF (IRTFLG .NE. 0) THEN
                CALL ERRT(46,'JPMSK1, VARP',NPIX)
                RETURN
             ENDIF
C            ZERO VARP ARRAY
             VARP = 0.0
          ENDIF

          ADA = 0.0
          IF (USE_PCA) THEN
             NFACMAX  =  MIN(NFAC,NUMIM,NPIX) 
             IF (NFAC .NE. NFACMAX) WRITE(NOUT,96) NFACMAX
96           FORMAT(' WARNING: NUMBER OF FACTORS LIMITED TO: ',I7)
             NFAC     = NFACMAX
             KIND_PCA = 1

          ELSEIF (ANS1(1:1) .EQ. 'C') THEN
C            FOR CORAN
             CALL RDPRM1S(ADA, NOT_USED, 'ADDITIVE CONSTANT',IRTFLG)
             IF (IRTFLG .NE. 0) RETURN

             NFAC     =  MIN(NFAC + 1, NUMIM, NPIX) - 1
             KIND_PCA = 0
          ENDIF
          ADANEG = (ADA .NE. 0.0)

          CALL FILERD(FILPRE,NLET,NULL,'OUTPUT FILE PREFIX~',IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

          FILPATS = FILPRE(1:NLET) // '_SEQ'//NULL
          FILPATP = FILPRE(1:NLET) // '_PIX'//NULL
          FILPATE = FILPRE(1:NLET) // '_EIG'//NULL
          FILPATC = FILPRE(1:NLET) // '_IMC'//NULL
          FILPATM = FILPRE(1:NLET) // '_MAS'//NULL
          FILPATT = FILPRE(1:NLET) // '_SET'//NULL

C         OPEN ALL FILES NEEDED IN SUBSEQUENT ROUTINES ON IO UNITS:
C         LUNS  =  OUTPUT FILE FOR IMAGE DATA
C         LUNT  =  OUTPUT FILE FOR TRANSPOSED IMAGE DATA (DIRECT ACCESS)
C         LUNI  =  OUTPUT FILE FOR IMAGE COORDINATES
C         LUNP  =  OUTPUT FILE FOR PIXEL COORDINATES
C         LUNE  =  OUTPUT FILE FOR EIGENVALUES & EIGENVECTORS
        
          CLOSE(LUNM)
          MAXIMT = 0
          CALL OPFILEC(0,.FALSE.,FILPATM,LUNM,'U',IFORM,
     &                 NSAMM,NROWM,NSLICEM,MAXIMT,' ',.FALSE.,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO  9999 
          CALL WRTVOL(LUNM,NSAMM,NROWM,1,NSLICEM,BUFM,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO  9999 


C         OPEN SEQUENTIAL ACCESS UNFORMATED FILE (_SEQ)
          CALL OPAUXFILE(.FALSE.,FILPATS,DATEXC,-LUNS,0,
     &                       'U', ' ',.TRUE.,IRTFLG)
C         CREATE _SEQ FILE HEADER 
          WRITE(LUNS) NUMIM, NPIX


C         OPEN FORMATTED IMAGE COORDINATE FILE (_IMC)
          CALL OPAUXFILE(.FALSE.,FILPATC,DATEXC,LUNI,0,
     &                       'U', ' ',.TRUE.,IRTFLG)
C         CREATE IMAGE COORDINATE FILE HEADER 
     
          WRITE(LUNI,95) NUMIM, NFAC, NSAMM, NROWM, NUMIM, KIND_PCA
95        FORMAT(10I10)


C         OPEN FORMATTED PIXEL COORDINATE FILE (_PIX)
          CALL OPAUXFILE(.FALSE.,FILPATP,DATEXC,LUNP,0,
     &                       'U', ' ',.TRUE.,IRTFLG)
C         CREATE PIXEL COORDINATE FILE HEADER 
          WRITE(LUNP,95) NPIX, NFAC, NSAMM, NROWM, NUMIM, KIND_PCA
                           
C         OPEN FORMATTED EIGENVALUE FILE
          CALL OPAUXFILE(.FALSE.,FILPATE,DATEXC,LUNE,0,
     &                       'U', ' ',.TRUE.,IRTFLG)

C         ALLOCATE IMAGE INPUT BUFFER
          ALLOCATE (BUFI(NPIX), STAT=IRTFLG)
          IF (IRTFLG .NE. 0) THEN
             CALL ERRT(46,'JPMSK1, BUFI',NPIX)
             GOTO 9999
          ENDIF

C         CHOOSE QUICKEST METHOD ACCORDING TO PROBLEM SIZE
          SIZE1 = NPIX  
          SIZE1 = SIZE1 * NPIX + 4 * NPIX 

C         IBM "MAX" CAN NOT HANDLE INTEGER * 8
          IF (USE_PCA) THEN
             SIZE2 = NUMIM
             SIZE2 = SIZE2 * NUMIM + 4 * NUMIM 
          ELSE
             SIZE2 = NUMIM
             SIZE2 = SIZE2 * NUMIM + 4 * NUMIM
          ENDIF

          IF (.NOT. USE_ITERPCA) THEN
	     WRITE(NOUT,98)  SIZE1,SIZE2
98	     FORMAT(/,'  MEMORY USE. IN-CORE:',I12,'  TRANSPOSED:',I12)

             TRANSPOSE = (SIZE2 .LT. SIZE1)
             IF (ANS1(2:2) .EQ. 'N') TRANSPOSE = .FALSE.
          ELSE
             TRANSPOSE = .FALSE.
          ENDIF

          IF (TRANSPOSE) THEN
             ALLOCATE(BUFT(NUMIM), STAT=IRTFLG)
             IF (IRTFLG .NE. 0) THEN
                CALL ERRT(46,'JPMSK1, BUFT',NUMIM)
                GOTO 9999
             ENDIF

C            FIND A INCORE SIZE FOR TRANSPOSITION MATRIX
             IDIV = 1
             DO 
                IWIDE = NUMIM / IDIV + MOD(NUMIM,IDIV)
                ALLOCATE(TRANS(IWIDE,NPIX), STAT=IRTFLG)
                IF (IRTFLG .EQ. 0) EXIT
                IDIV = IDIV + 1
             ENDDO

C            OPEN DIRECT ACCESS UNFORMATED FILE (_SET) FOR TRANSPOSE
             CALL OPAUXFILE(.FALSE.,FILPATT,DATEXC,LUNT,NUMIM*4,
     &                       'U', ' ',.TRUE.,IRTFLG)

          ENDIF

C         LOOP OVER ALL IMAGES
          IGOT     = 0
          ICOLS    = 0
          SKIPREST = .FALSE.
          FMINALL  = HUGE(FMINALL)

          DO IM = 1, NUMIM
C           RANDOM SHUFFLE OPTION
            IF (IM .GT. 1 .AND. FLAGR) THEN
C              'R'  OPTION USED.
C              NSH = SH * FLOAT(NPIX) + 0.5
C              CALL PERMUT(BUFM, NPIX, NSH)

            ELSE
              NLET  = 0
              CALL FILGET(FILPATI,FILNAM,NLET,INUMBR(IM),IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999

              MAXIMT = 0
              CLOSE(LUN)
              CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'O',IFORM,
     &                   NSAM,NROW,NSLICE,MAXIMT,' ',.FALSE.,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO  9999 

              IF ((FMIN + ADA) .LT. 0.0 .AND. .NOT. USE_PCA)  THEN
                 FMINALL  = MIN(FMIN,FMINALL)
                 SKIPREST = .TRUE.
                 CYCLE
              ENDIF

C             COMPARE MASK DIMENSIONS WITH IMAGE DIMENSIONS,ETC
              IF (IFORM .NE. 1 .AND. IFORM .NE. 3) THEN
                  CALL ERRT(2, 'JPMSK1', NDUM)
                  GOTO  9999
              ELSEIF ((NSAMM   .NE. NSAM)   .OR. 
     &                (NROWM   .NE. NROW)   .OR.
     &                (NSLICEM .NE. NSLICE)) THEN
                 WRITE(NOUT, 93) NSAM, NROW, NSAMM, NROWM, 
     &                            NSLICE,NSLICEM
93               FORMAT('*** IMAGE DIMENSION (',I4,',',I4,',',I4,
     &                  ')  NOT SAME AS MASK (',I4,',',I4,',',I4,')')
                 CALL ERRT(100,'JPMSK1',NE)
                 GOTO 9999
              ENDIF
    
C             READ IMAGE AREA WHERE MASK > 0.5 INTO CORE AND
C	      COMPUTE ITS MASK-RELATED,  PRECISE AVERAGE
              WEIGHTIT = 0.0
              ILOCM    = 0
              IPIX     = 0

              DO IROW = 1,NROW*NSLICE
                 CALL REDLIN(LUN,BUF,NSAM,IROW)
                 DO ISAM = 1,NSAM
                    ILOCM = ILOCM + 1
                    IF (BUFM(ILOCM) .GT. 0.5) THEN
C                      INSIDE MASK, USE THIS PIXEL
                       VAL    = BUF(ISAM)
                       IF (ADANEG) VAL = VAL + ADA

C	               WEIGHTIT = SUM OF THE ELEMENTS(MASK > 0.5)    .
                       WEIGHTIT      = WEIGHTIT + VAL               

                       IPIX          = IPIX + 1

                       WEIGHTP(IPIX) = WEIGHTP(IPIX) + VAL
                       BUFI(IPIX)    = VAL 

                       IF (USE_ITERPCA) THEN
C                         CALCULATE VARIANCE ALSO
                          VARP(IPIX) = VARP(IPIX) + VAL ** 2
                       ENDIF
                    ENDIF
                 ENDDO
              ENDDO
           ENDIF

           IF (SKIPREST) THEN
              WRITE(NOUT,198) FMINALL
 198          FORMAT(' *** OVERALL MINIMUM: ',G13.7)  
              CALL ERRT(101,'CORAN CAN NOT ACCEPT NEGATIVE PIXELS',NDUM)
              GOTO 9999
           ENDIF

C          WEIGHTIT = SUM OF ALL PIXEL VALUES UNDER MASK IN IMAGE.
           WEIGHTI(IM) = WEIGHTIT

           FIM = INUMBR(IM)
           WRITE(LUNS,IOSTAT=IERR) BUFI,FIM

           ICOLS = ICOLS + 1
           IF (TRANSPOSE) THEN
C             INCORE, TRANSPOSED (DOES THE TRANSPOSING)
              TRANS(ICOLS, :) = BUFI

              IF (ICOLS .GE. IWIDE) THEN
C                MUST WRITE OUT TRANSPOSED MATRIX SECTION

                 DO I = 1,NPIX
                    IF (IGOT .GT. 0) THEN
C                      HAVE DIVISION ALREADY IN FILE
                       READ(LUNT,REC=I,IOSTAT=IERR) BUFT
                    ENDIF
	            BUFT(IGOT+1:) = TRANS(1:ICOLS,I)
                    
                    WRITE(LUNT,REC=I,IOSTAT=IERR) BUFT
                 ENDDO
                 IGOT  = IM
                 ICOLS = 0
              ENDIF
           ENDIF
        ENDDO
 
C       SUMW IS THE SUM OF ALL THE PIXEL VALUES IN ALL THE IMAGES.
        SUMW = SUM(WEIGHTP)

C       FREE UP ALLOCATIONS FOR USE IN NEXT STEP
        IF (ALLOCATED(BUFI))  DEALLOCATE(BUFI)
        IF (ALLOCATED(BUFM))  DEALLOCATE(BUFM)
        IF (ALLOCATED(BUFT))  DEALLOCATE(BUFT)
        IF (ALLOCATED(TRANS)) DEALLOCATE(TRANS)

        IF (USE_ITERPCA) THEN
C          ITERATIVE PCA
           CALL SPCA3(NUMIM,   NFAC,    NPIX, INUMBR, 
     &                LUNS,    LUNI,    LUNP, LUNE, 
     &                WEIGHTI, WEIGHTP, SUMW, VARP, LUNT)

         ELSEIF (ANS1(1:1) .NE. 'S') THEN
C          PCA OR CORAN
           LUNIN = LUNS
           IF (TRANSPOSE) LUNIN = LUNT
           CALL SCORAN3(NUMIM,    NFAC,    NPIX,      INUMBR, USE_PCA, 
     &                   LUNIN,   LUNI,    LUNP,      LUNE, 
     &                   WEIGHTI, WEIGHTP, TRANSPOSE, SUMW)
        ENDIF

9999    IF (ALLOCATED(BUFM))    DEALLOCATE(BUFM)
        IF (ALLOCATED(BUFI))    DEALLOCATE(BUFI)
        IF (ALLOCATED(TRANS))   DEALLOCATE(TRANS)
        IF (ALLOCATED(WEIGHTI)) DEALLOCATE(WEIGHTI)
        IF (ALLOCATED(WEIGHTP)) DEALLOCATE(WEIGHTP)
        IF (ALLOCATED(VARP))    DEALLOCATE(VARP)

C       CLOSE ALL FILES THAT MIGHT BE OPEN
        CLOSE(LUNS)
        CLOSE(LUNI)
        CLOSE(LUNP)
        CLOSE(LUNM)
        CLOSE(LUNE)
        CLOSE(LUNT)
        CLOSE(LUN)

	RETURN
        END

