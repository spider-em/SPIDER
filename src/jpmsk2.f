
C++************************************************************ 5/14/85   
C
C  JPMSK2   
C           RECONS RESTORED FROM LISTING       9/16/85 JF
C           LONG FILE NAMES                    FEB 89 ARDEAN LEITH
C           PCA REWRITTEN                      OCT 03 MAHIEDDINE LADJADJ
C           USED OPAUXFILE & OPFILE            NOV 00 ARDEAN LEITH
C           SRIPE PARAMETERS CHANGED           JAN 01 ARDEAN LEITH
C           FMIN --> FMINT & INSERTED RECDUM   APR 01 ARDEAN LEITH
C           REWRITTEN                          OCT 03 ARDEAN LEITH
C           MASK BUG                           NOV 03 ARDEAN LEITH
C           _EIG & _IMC FORMATS CHANGED        JAN 04 ARDEAN LEITH
C           _REMAKE J BUG                      APR 04 ARDEAN LEITH
C           SRE COORDINATE BUG                 JUL 04 ARDEAN LEITH
C           RDPRAF REMOVED                     DEC 05 ARDEAN LEITH 
C           CREATES VOLUMES OK                 MAR 06 ARDEAN LEITH 
C           PIXEL FILE FDUM                    JUN 09 ARDEAN LEITH
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
C JPMSK2(LUN,LUNI,LUNE,LUNP,LUNM,LUNDOC)
C
C PURPOSE:   RECONSTITUTE IMAGES USING SELECTED FACTORS
C            OF A CORRESPONDENCE ANALYSIS STORED IN
C            FILES: _IMC , _PIX , _MAS, AND _EIG 
C
C PARAMETERS:
C    LUN      LOGICAL UNIT NUMBER OF INPUT/OUTPUT FILE
C    LUNI     LOGICAL UNIT NUMBER OF IMAGE COORDINATES FILE
C    LUNE     LOGICAL UNIT NUMBER OF EIGEN FILE
C    LUNP     LOGICAL UNIT NUMBER OF PIXEL COORDINATES FILE
C    LUNM     LOGICAL UNIT NUMBER OF MASK FILE
C
C COMMANDS SUPPORTED :
C    CA SR  -- RECONSTITUTE FULL IMAGE INCLUDING 0-FACTOR
C    CA SRD -- RECONSTITUTE DIFFERENTIAL IMAGE, WITHOUT 0-FACTOR
C    CA SRI -- RECONSTITUTE IMPORTANCE IMAGE (AS CA_SRD BUT
C                 WITHOUT COORDINATE WEIGHTING)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

	SUBROUTINE JPMSK2(LUN,LUNI,LUNE,LUNP,LUNM)

	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'

        COMMON /IOBUF/ BUF(NBUFSIZ)

        COMMON  /COMMUN/         FILNAM,FILNMP,FILNMC,FILNME, FILNMM,
     &                           FILPRE
        CHARACTER(LEN=MAXNAM) :: FILNAM,FILNMP,FILNMC,FILNME, FILNMM,
     &                           FILPRE

        CHARACTER(LEN=1)                   :: NULL, ANS
        LOGICAL                            :: USE_PCA,ADD_AVRG

        REAL, ALLOCATABLE, DIMENSION(:)    :: BUFM
        INTEGER, ALLOCATABLE, DIMENSION(:) :: NUMFAC
        REAL, ALLOCATABLE, DIMENSION(:)    :: QBUF
        REAL, ALLOCATABLE, DIMENSION(:)    :: COO
        LOGICAL                            :: SRE, SRA, SR

#ifndef SP_32
        INTEGER * 8 :: JTOT,IASK8,IOK
#else
        INTEGER * 4 :: JTOT,IASK8,IOK
#endif

        NULL = CHAR(0)
       
        SRE  = (FCHAR(4:6) .EQ. 'SRE')
        SRA  = (FCHAR(4:6) .EQ. 'SRA')
        SR   = (.NOT. SRE .AND. .NOT. SRA) .OR. 
     &         (FCHAR(4:6) .EQ. 'SRI') .OR. 
     &         (FCHAR(4:6) .EQ. 'SRD')

        CALL FILERD(FILPRE,NLET,NULL,
     &              'CORAN/PCA FILE PREFIX (e.g. CORAN_)~',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (SR) THEN
C	   IMUSED IS THE NUMBER OF IMAGES WE WANT TO RECONSTRUCT.
           NILMAX = NIMAX	 
           CALL FILELIST(.FALSE.,LUNDOC,NULL,IDUM,
     &               INUMBR,NILMAX,IMUSED,NULL,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
        ELSE

        ENDIF

        FILNMC = FILPRE(1:NLET) // '_IMC'//NULL
        FILNMP = FILPRE(1:NLET) // '_PIX'//NULL
        FILNME = FILPRE(1:NLET) // '_EIG'//NULL
        FILNMM = FILPRE(1:NLET) // '_MAS'//NULL

C       OPEN IMAGE COORDINATE FILE
        CALL OPAUXFILE(.FALSE.,FILNMC,DATEXC,LUNI,0,
     &                       'O', ' ',.TRUE.,IRTFLG)
        READ(LUNI,*) NUMIM, NFAC, NSAM, NROW, NDUM, KIND_PCA

        IF (KIND_PCA .EQ. 1) THEN

          IF (FCHAR(4:6) .EQ. 'SRD' .OR. FCHAR(4:6) .EQ. 'SRI') THEN
             CALL ERRT(101,'OPERATION DOES NOT WORK FOR PCA FILES',NE)
             GOTO 9999
          ENDIF

          USE_PCA = .TRUE.

          CALL RDPRMC(ANS, NCHAR, .TRUE.,
     &      'SUBTRACT AVERAGE IN PCA RECONSTITUTION? (N/Y)',NULL,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN
          ADD_AVRG = (ANS .EQ. 'N') 

        ELSE
          USE_PCA = .FALSE.
     
C         IDIFF= SWITCH DETERMINING WHETHER THE AVERAGE IMAGE 
C         SHOULD BE SUBTRACTED.
C         IDIFF=1 -- SUBTRACT. IDIFF=-1 -- SUBTRACT 
C         AVERAGE IMAGE, AND THEN DIVIDE RESULT BY AVERAGE IMAGE 
C         TO OBTAIN  SO-CALLED "IMPORTANCE IMAGE" USED IN EARLIER
C         LITERATURE.

   	  IDIFF = 0
          IF (FCHAR(6:6) .EQ. 'D') THEN
C            'CA SRD'
             IDIFF = 1 
	  ELSEIF (FCHAR(6:6) .EQ. 'I') THEN
C            'CA SRI'
             IDIFF = -1
          ENDIF
        ENDIF
 
C       OPEN PIXEL COORDINATE FILE
        CALL OPAUXFILE(.FALSE.,FILNMP,DATEXC,LUNP,0,
     &                       'O', ' ',.TRUE.,IRTFLG)
	READ(LUNP,*) NMASKP,NDUM1,NDUM2,NDUM3,NDUM4 
        WRITE(NOUT,*)' NUMBER OF ACTIVE PIXELS: ',NMASKP
        IF (NMASKP .LE. 0) THEN
           CALL ERRT(101,'NO PIXELS',NE)
           GOTO 9999
        ENDIF

C       OPEN MASK FILE, RETURNS IFORM  AND NSLICE
        CALL OPFILEC(0,.FALSE.,FILNMM,LUNM,'O',IFORM,
     &              NSAMM,NROWM,NSLICE, MAXIMT,' ',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
 
        IF ((NSAMM .NE. NSAM) .OR. (NROWM .NE. NROW)) THEN
              WRITE(NOUT, 93) NSAM, NROW, NSAMM, NROWM 
93            FORMAT('*** IMAGE DIMENSION (',I4,',',I4,')', 
     &               '  NOT SAME AS MASK  (',I4,',',I4,')')
              CALL ERRT(100,'JPMSK2',NE)
              GOTO 9999
        ENDIF
 
        NPIX = NSAM*NROW*NSLICE 
        ALLOCATE (BUFM(NPIX),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'JPMSK2; BUFM',NPIX)
           RETURN
        ENDIF

C       READ MASK IN BUFM
        CALL REDVOL(LUNM,NSAM,NROW,1,NSLICE,BUFM,IRTFLG)
        NMASK = 0
        DO I=1,NSAM*NROW*NSLICE
           IF (BUFM(I) .GT. 0.5) NMASK = NMASK + 1
        ENDDO

        WRITE(NOUT,*)' NUMBER OF PIXELS UNDER MASK: ',NMASK
        IF (NMASK .LE. 0) THEN
           CALL ERRT(101,'NO PIXELS UNDER MASK',NE)
           GOTO 9999
        ENDIF

C       FACTOR SELECTION 
        ALLOCATE(NUMFAC(NFAC), COO(NFAC),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            MWANT = 2 * NFAC
            CALL ERRT(46,'NUMFAC & COO',MWANT)
            GOTO 9999
	ENDIF

        WRITE(NOUT,*)' NUMBER OF FACTORS AVAILABLE: ',NFAC
        NUSE   = NFAC
        CALL FILELISTB(LUNDOC,NUMFAC,NFAC,NUSE,
     &     'FACTOR NUMBERS OR DOC. FILE NAME FOR FACTOR LIST',
     &     IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
  
C       OPEN EIGENVALUE FILE
        CALL OPAUXFILE(.FALSE.,FILNME,DATEXC,LUNE,0,
     &                       'O', ' ',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       READ _EIG FILE HEADER
	READ(LUNE,*) NFACT, SUMP, TRACE, KIND_PCA

C       MEMORY ALLOCATION

C       JPIC  = OFFSET OF ARRAY CONTAINING RECONSTITUTED IMAGE

        IF (USE_PCA) THEN
           JPIC = 1
	   LCI	= JPIC  + NPIX
	   LIDI	= LCI 	+ NUMIM * NFAC
	   LCP	= LIDI	+ NUMIM
	   LD	= LCP 	+ NMASK  * NFAC
	   JIM  = LD	+ NFAC
	   LPJ  = JIM   + NMASK
	   JTOT = LPJ   + NMASK

        ELSE
           JPIC = 1
	   LCI	= JPIC  + NPIX
	   LWI	= LCI 	+ NUMIM * NFAC
	   LIDI	= LWI	+ NUMIM
	   LCP	= LIDI	+ NUMIM
	   LD	= LCP 	+ NMASK  * NFAC
	   JIM	= LD	+ NFAC
	   LWP  = JIM   + NMASK
	   JTOT = LWP   + NMASK
        ENDIF

C       COMPLAIN IF EXCESSIVE ALLOCATION
        IASK8 = JTOT * 4
        CALL BIGALLOC(IASK8,IOK,.FALSE.,.TRUE.,IRTFLG)

        ALLOCATE(QBUF(JTOT),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'QBUF',JTOT)
           GOTO 9999
        ENDIF

C       LOAD THE _IMC, _PIX, & _EIG DATA INTO LOCAL ARRAYS
        IF (USE_PCA) THEN
           CALL PCA_SRIP2(NUMIM, NMASK, NFAC, LUNI, LUNP, LUNE,
     &                    QBUF(LCI), QBUF(LIDI), QBUF(LCP), 
     &                    QBUF(LD),  QBUF(LPJ), IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
        ELSE
	   CALL SRIP2(NUMIM,NMASK,NFAC,LUNI,LUNP,LUNE,
     &		      QBUF(LCI),QBUF(LWI),QBUF(LIDI),
     &                QBUF(LCP),QBUF(LWP), QBUF(LD))

	ENDIF 

	CLOSE(LUNI)
	CLOSE(LUNP)

C       NOW THAT THE NUMBER OF FACTORS USED IN THE INITIAL ANALYSIS IS
C       KNOWN, SWAP FACTOR SELECTION ARRAY AND CHECK AGAINST THIS NUMBER

        DO  I=1,NUSE
           IF (NUMFAC(I) .GT. NFAC) THEN
              WRITE(NOUT,9801) NFAC
9801          FORMAT(' ** AXIS NUMBER OUT OF RANGE. MAXIMUM = ',I2)
              WRITE(NOUT,9802) NFAC, NUMFAC(I)
9802          FORMAT(' Number of factors (or eigenvectors) used to ',
     &               ' create data files is: ',I6,
     &               ' you have asked for factor: ',I4)
	     GOTO 9999
          ENDIF
        ENDDO

        IF (SRA) THEN
C          ----------------------------- ARBITRARY ----------- 'CA SRA'
C          RECONSTITUTION FOR DUMMY IMAGES USING ARBITRARY COORDINATES
           NVAL = NUSE
           CALL RDPRA('COORDINATES FOR EACH FACTOR SELECTED',
     &         NVAL,0,.FALSE.,QBUF,NUSE,IRTFLG)

C          ZERO COO ARRAY
	   COO = 0.0

           IMAX = 0
           DO K=1,NUSE
               DO  I=1,NFAC
                 IF (NUMFAC(K) .EQ. I) THEN
                    COO(I) = QBUF(K)
                    IF (I .GT. IMAX) IMAX = I
                    EXIT
                 ENDIF
	       ENDDO
          ENDDO

C         WRITE OUT COORDINATES TO BE USED
          WRITE(NOUT,*) ' COORDINATES USED:'
          WRITE(NOUT,281) (COO(I),I=1,IMAX)
281       FORMAT(1X,8F8.3)
          WRITE(NOUT,*) ' '

C         NOW CALL RECONSTITUTION ROUTINE FOR DUMMY IMAGE
          IF (USE_PCA) THEN
              CALL PCA_RECDUM(NUMIM, NMASK, NUSE, NUMFAC, 
     &              NFAC, AV,QBUF(JIM), QBUF(LCP),
     &              COO, QBUF(LD), QBUF(LPJ), ADD_AVRG)
	  ELSE
              CALL RECDUM(NUMIM,NMASK,NUSE,NUMFAC,NFAC,SUMP,
     &              FMAXT,FMINT,AV,SIG,IDIFF,QBUF(JIM),
     &              QBUF(LCP),COO,QBUF(LD),QBUF(LWI),QBUF(LWP))
          ENDIF

          MAXIM = 0
          CALL OPFILEC(0,.TRUE.,FILNAM,LUN,'U',IFORM,NSAM,NROW,NSLICE,
     &                   MAXIM,'RECONSTITUTED OUTPUT',.FALSE.,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9999

C         REPLACE PIXELS IN MASK AREA BY THEIR TRUE VALUES
C         COPY OUTPUT IMAGE QBUF(JIM:-) TO RECONSTITUTED IMAGE
C         QBUF(JPIC:-) SUBJECT TO MASK IN BUFM

          IF (USE_PCA) THEN
              FMINT = MINVAL(QBUF(JIM:JIM+NMASK-1)) 
          ENDIF

C         PUT FMINT IN THE RECONSTITUTED IMAGE: QBUF(JPIC:)
          QBUF(JPIC:JPIC+NPIX-1) = FMINT

          ILOC = JIM 
          DO  L = 1, NPIX
             IF (BUFM(L) .GT. 0.5) THEN
                QBUF(JPIC+L-1) = QBUF(ILOC)
                ILOC           = ILOC + 1
             ENDIF
          ENDDO

C         WRITE IMAGE ONTO DISK
          CALL WRTVOL(LUN,NSAM,NROW,1,NSLICE,QBUF(JPIC),IRTFLG)
          CLOSE(LUN)
	   
        ELSEIF (SRE) THEN
C         ------------- EIGENIMAGE RECONSTITUTION ----------- 'CA SRE'
 
C          ZERO COO ARRAY
	   COO = 0.0

C          SET COORDINATE TO 1.0 FOR EACH FACTOR IN USE
           DO K=1,NUSE
               DO  I=1,NFAC
                 IF (NUMFAC(K) .EQ. I) THEN
                    COO(I) = 1.0
                    EXIT
                 ENDIF
	       ENDDO
           ENDDO


C         CALL RECONSTITUTION ROUTINE FOR DUMMY IMAGE
          IF (USE_PCA) THEN
              CALL PCA_RECDUM(NUMIM, NMASK, NUSE, NUMFAC, 
     &              NFAC, AV,QBUF(JIM), QBUF(LCP),
     &              COO, QBUF(LD), QBUF(LPJ), ADD_AVRG)
	  ELSE
              CALL RECDUM(NUMIM,NMASK,NUSE,NUMFAC,NFAC,SUMP,
     &              FMAXT,FMINT,AV,SIG,IDIFF,QBUF(JIM),
     &              QBUF(LCP),COO,QBUF(LD),QBUF(LWI),QBUF(LWP))
          ENDIF

          MAXIM = 0
          CALL OPFILEC(0,.TRUE.,FILNAM,LUN,'U',IFORM,NSAM,NROW,NSLICE,
     &                 MAXIM,'RECONSTITUTED OUTPUT',.FALSE.,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9999

C         REPLACE PIXELS IN MASK AREA BY THEIR TRUE VALUES
C         COPY OUTPUT IMAGE QBUF(JIM:-) TO RECONSTITUTED IMAGE
C         QBUF(JPIC:-) SUBJECT TO MASK IN BUFM

          IF (USE_PCA) THEN
              FMINT = MINVAL(QBUF(JIM:JIM+NMASK-1)) 
          ENDIF

C         PUT FMINT IN THE RECONSTITUTED IMAGE: QBUF(JPIC:)
          QBUF(JPIC:JPIC+NPIX-1) = FMINT

          ILOC = JIM 
          DO  L = 1, NPIX
             IF (BUFM(L) .GT. 0.5) THEN
                QBUF(JPIC+L-1) = QBUF(ILOC)
                ILOC           = ILOC + 1
             ENDIF
          ENDDO

C         WRITE IMAGE ONTO DISK
          CALL WRTVOL(LUN,NSAM,NROW,1,NSLICE,QBUF(JPIC),IRTFLG)
          CLOSE(LUN)
	   
        ELSEIF (SR) THEN
C          ---------------- IMAGE RECONSTITUTION  ------------ 'CA SR'
           CALL FILERD(FILPRE,NLET,NULL,'OUTPUT FILE PREFIX~',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          LOOP OVER ALL REQUESTED IMAGES
           DO IM=1,IMUSED
             IMNO = 0
             DO I = 1,NUMIM 
	        IF (TRANSFER(QBUF(LIDI + I - 1),I) .EQ. INUMBR(IM))THEN
                   IMNO = I
                   EXIT
                ENDIF
             ENDDO
             IF (IMNO .EQ. 0) THEN
                WRITE(NOUT,*) ' IMAGE: ',INUMBR(IM),'  NOT FOUND'
                CYCLE
             ENDIF

C            RECONSTITUTION FOR CURRENT IMAGE
             IF (USE_PCA) THEN
               CALL PCA_REMAKE(IMNO, NUMIM, NMASK, NUSE, NUMFAC, 
     &              NFAC, AV,QBUF(JIM), QBUF(LCP), 
     &              QBUF(LCI), QBUF(LD), QBUF(LPJ), ADD_AVRG)
             ELSE
               CALL RECONS(IMNO,NUMIM,NMASK,NUSE,NUMFAC,NFAC,SUMP,
     &                  FMAXT,FMINT,AV,IDIFF,QBUF(JIM),QBUF(LCP),
     &                  QBUF(LCI),QBUF(LD),QBUF(LWI),QBUF(LWP))
             ENDIF

C            OPEN IMAGE OUTPUT FILE
             NLET1 = 0
             CALL FILGET(FILPRE,FILNAM,NLET1,IMNO,IRTFLG)
             IF (IRTFLG .NE. 0)  GOTO 9999

             MAXIM = 0
             CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'U',IFORM,
     &                   NSAM,NROW,NSLICE,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
             IF (IRTFLG .NE. 0) GOTO 9999

C            REPLACE PIXELS IN MASK AREA BY THEIR TRUE VALUES
C            COPY OUTPUT IMAGE QBUF(JIM:-) TO RECONSTITUTED IMAGE
C            QBUF(JPIC:-) SUBJECT TO MASK IN BUFM

             IF (USE_PCA) THEN
                 FMINT = MINVAL(QBUF(JIM:JIM+NMASK-1)) 
             ENDIF

C	     PUT FMINT IN THE RECONSTITUTED IMAGE: QBUF(JPIC:)
             QBUF(JPIC:JPIC+NPIX-1) = FMINT

             ILOC = JIM 
             DO  L = 1, NPIX
                IF (BUFM(L) .GT. 0.5) THEN
                   QBUF(JPIC+L-1) = QBUF(ILOC)
                   ILOC           = ILOC + 1
                ENDIF
	     ENDDO

C            WRITE IMAGE ONTO DISK
             CALL WRTVOL(LUN,NSAM,NROW,1,NSLICE,QBUF(JPIC),IRTFLG)
             CLOSE(LUN)
	   ENDDO 
        ENDIF

9999	CLOSE(LUN)
    	CLOSE(LUNI)
    	CLOSE(LUNP)
        IF (ALLOCATED(COO))    DEALLOCATE(COO)
	IF (ALLOCATED(NUMFAC)) DEALLOCATE(NUMFAC)
	IF (ALLOCATED(QBUF))   DEALLOCATE(QBUF)
	IF (ALLOCATED(BUFM))   DEALLOCATE(BUFM)

	END


C ++********************************************************************
C                                                                      *
C  SRIP2(NIMA,NPIX,NFAC,LUNI,LUNP,LUNE,CI,WI,IDI,CP,WP,D)
C
C  REMOVED FROM SRIPE.F    DEC 90 AL
C
C  PURPOSE:  LOADS CORAN SPECIFIC FILES                                                          *
C                                                                      *
C  PARAMETERS:  CI                                            RETURNED                                                         *
C               WI                                            RETURNED                                                         *
C               IDI                                           RETURNED                                                         *
C               CP                                            RETURNED                                                         *
C               WP                                            RETURNED                                                         *
C               D                                             RETURNED                                                         *
C                                                                      *
C***********************************************************************

	SUBROUTINE SRIP2(NIMA,NPIX,NFAC,LUNI,LUNP,LUNE,
     &                   CI,WI,IDI,CP,WP,D) 

	DIMENSION CI(NIMA,NFAC), WI(NIMA), IDI(NIMA)
	DIMENSION CP(NPIX,NFAC), WP(NPIX), D(NFAC)

	DO  I=1,NIMA
	   READ(LUNI,*) (CI(I,N),N=1,NFAC), WI(I),CODUM,FIM,FACTDUM
           IDI(I) = FIM
	ENDDO

	DO  I=1,NPIX
	   READ(LUNP,*) (CP(I,N),N=1,NFAC), WP(I),CODUM,FPIXDUM,FDUM
	ENDDO

	DO  I=1,NFAC
	   READ(LUNE,*) D(I)
        ENDDO

	END

C **********************************************************************
C
C *  AUTHOR :  MAHIEDDINE LADJADJ                                          *
C
C  PURPOSE: LOADS PCA INPUT FILES                                                          *
C                                                                      *
C  PARAMETERS:  CI                                            RETURNED                                                         *
C               IDI                                           RETURNED                                                         *
C               CP                                            RETURNED                                                         *
C               D                                             RETURNED                                                         *
C               PJ                                            RETURNED                                                         *
C                                                                      *
C **********************************************************************

        SUBROUTINE PCA_SRIP2(NIMA,NPIX,NFAC,LUNI,LUNP,LUNE,
     &             CI,IDI,CP,D,PJ,IRTFLG)

        INCLUDE 'CMBLOCK.INC'

        DIMENSION CI(NIMA,NFAC), IDI(NIMA), CP(NPIX,NFAC)
        DIMENSION PJ(NPIX), D(NFAC)

        IRTFLG = 0

c       READ THE _IMC  FILE.
        DO I=1,NIMA
           READ(LUNI,*) (CI(I,N),N=1,NFAC),WDUM,DORDUM,FIM,ACTDUM
           IDI(I) = FIM
        ENDDO

C       READ THE _PIX  FILE.
        DO I=1,NPIX
           READ(LUNP,*) (CP(I,N),N=1,NFAC),PJ(I),CODUM,FPIXDUM,FDUM
        ENDDO

C       READ THE PCA _EIG  FILE TO GET THE EIGENVALUES.
        DO I=1,NFAC
           READ(LUNE,*) D(I)
        ENDDO

        WRITE(NOUT,*) ' DATA IS FROM PCA USING COVARIANCE MATRIX'
        RETURN

        END

C ++********************************************************************
C  RECONS                                                                     *
C                                                                      *
C  RESTORED FROM LISTING 9/16/85 JF
C                                                                      *
C **********************************************************************
C **********************************************************************
C
C PURPOSE: RECONSTITUTES IMAGES USING COORDINATES FROM _IMC & _PIX FILES
C                                                                      *
C  PARAMETERS:  IMNO                                          SENT
C               NIMA                                          SENT
C               NPIX                                          SENT
C               NF                                            SENT
C               NV                                            SENT                                                         *
C               NFAC                                          SENT
C               SUMP                                          SENT
C               FMAX1,FMIN1                                   RETURNED
C               PIA                                           RETURNED                                                         *
C               IDIFF                                         SENT
C               RIM                                           RETURNED                                                         *
C               CP                                            SENT                                                         *
C               CI                                            SENT                                                         *
C               D                                             SENT                                                         *
C               WI                                            SENT                                                        *
C               WP                                            SENT                                                                                                   *
C                                                                      *
C***********************************************************************

        SUBROUTINE RECONS(IMNO,NIMA,NPIX,NF,NV,NFAC,SUMP,
     &              FMAXT,FMINT,PIA,IDIFF,RIM,CP,CI,D,WI,WP)

        DIMENSION RIM(NPIX), CP(NPIX,NFAC), CI(NIMA,NFAC), D(NFAC)
        DIMENSION WP(NPIX), WI(NIMA), NV(NF)

        COMMON BUF(1)

        INCLUDE 'CMBLOCK.INC'

C       ZERO BUF
        BUF(1:NFAC) = 0.0

        DO  L = 1,NF
           K      = NV(L)
           BUF(K) = CI(IMNO,K)
	ENDDO

C       WRITE OUT COORDINATES EFFECTIVE IN THIS RECONSTITUTION
        WRITE(NOUT,*) ' COORDINATES EFFECTIVE IN THIS RECONSTITUTION:'

        WRITE(NOUT,281) (BUF(L),L=1,NFAC)
281     FORMAT(2X,8F9.5)

        WRITE(NOUT,*) ' '

        DELTA = 1.0 - FLOAT(IDIFF)
        SUM   = 0

C       CHANGED FROM 10000 APR 01 al
        FMINT = HUGE(FMINT)
        FMAXT = -FMIN

C       LOOP OVER ALL PIXELS

        DO  I=1,NPIX
           RIM(I) = 0.0
           DO  L=1,NF
              K      = NV(L)
              RIM(I) = RIM(I) + CP(I,K) * CI(IMNO,K) / SQRT(D(K))
           ENDDO

           IF (IDIFF .GE. 0) THEN
              RIM(I) = WP(I) * WI(IMNO) * (DELTA + RIM(I)) * SUMP
           ENDIF

           FMAXT = MAX (RIM(I),FMAXT)
           FMINT = MIN (RIM(I),FMINT)

           SUM = SUM + RIM(I)
        ENDDO

        PIA = SUM / FLOAT(NPIX)

        END


C ++********************************************************************
C                                                                      *
C   RECDUM                                                             *
C                                                                      *
C **********************************************************************
C **********************************************************************
C                                                                      *
C  RECDUM(NIMA,NPIX,NF,NV,NFAC,SUMP,
C         FMAX,FMIN,PIA,SIG,IDIFF, RIM,CP,CI,D,WI,WP) 
C                                                                      *
C  PURPOSE:  RECONSTITUTES DUMMY IMAGES USING USER-SUPPLIED COORDINATES
C
C  PARAMETERS:
C   NIMA  		NUMBER OF ACTIVE IMAGES
C   NPIX   		NUMBER OF ACTIVE PIXELS
C   NFAC    		NUMBER OF FACTORS USED IN THE ANALYSIS
C   SUMP     		TOTAL WEIGHT
C   FMAX, FMIN 		MAXIMUM, MINIMUM OF RECONSTITUTED IMAGE
C   PIA			AVERAGE OF RECONSTITUTED IMAGE (SUPPLIED)
C   SIG			STANDARD DEVIATION OF RECONSTITUTED IMAGE
C   IDIFF    		FLAG INDICATING THE TYPE OF IMAGE COMPUTED
C 	=1    		   DIFFERENTIAL IMAGE (0-FACTOR NOT ADDED)
C       =0		   TOTAL IMAGE        (0-FACTOR IMAGE ADDED)
C       =-1    		   IMPORTANCE IMAGE (SAME AS DIFFERENTIAL BUT
C			   NO WEIGHTING USED)
C   RIM(NPIX)		OUTPUT IMAGE
C   CP(NPIX,NFAC)	PIXEL COOS AS NEEDED BY SRIPE
C   CI(NFAC)		IMAGE COOS
C   D(NFAC)		EIGENVALUES AS NEEDED BY SRIPE
C   WP(NPIX)		PIXEL WEIGHTS
C   WI(NIMA)		IMAGE WEIGHTS 
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

         SUBROUTINE RECDUM(NIMA,NPIX,NF,NV,NFAC,SUMP,
     &             FMAX,FMIN,PIA,SIG,IDIFF, RIM,CP,CI,D,WI,WP)

         DIMENSION RIM(NPIX), CP(NPIX,NFAC), CI(NFAC), D(NFAC)
         DIMENSION WP(NPIX), WI(NIMA), NV(NF)

         IF (IDIFF .GE. 0) THEN
            WIAVG = SUM(WI) / FLOAT(NIMA) 
         ENDIF

         FMIN  = HUGE(FMIN)
         FMAX  = -FMIN
         PIA   = 0
         DELTA = 1. - FLOAT(IDIFF)

         DO I = 1, NPIX
           RIM(I) = 0.0
           DO L = 1,NF
              K      = NV(L)
              RIM(I) = RIM(I) + CP(I,K) * CI(K) / SQRT(D(K))
           ENDDO

           IF (IDIFF .GE. 0) THEN
              RIM(I) = WP(I) * WIAVG * (DELTA + RIM(I)) * SUMP
           ENDIF

           IF (RIM(I) .GT .FMAX) FMAX = RIM(I)
           IF (RIM(I) .LT .FMIN) FMIN = RIM(I)
           PIA = PIA + RIM(I)
         ENDDO

         PIA = PIA / FLOAT(NPIX)

         RETURN
         END


C++*********************************************************************
C
C PCA_RECDUM(NIMA, NPIX, NF, NV, NFAC,
C            PIA, RIM, CP, CI, D, PJ, ADD_AVRG)
C
C PURPOSE: RECONSTITUTES DUMMY IMAGES USING USER-SUPPLIED COORDINATES
C
C PARAMETERS:
C    NIMA	        NUMBER OF ACTIVE IMAGES
C    NPIX	        NUMBER OF ACTIVE PIXELS
C    NFAC	        NUMBER OF FACTORS USED IN THE ANALYSIS
C    SUMP	        TOTAL WEIGHT
C    PIA	        AVERAGE OF RECONSTITUTED IMAGE (SUPPLIED)
C    RIM(NPIX)          OUTPUT IMAGE
C    CP(NPIX,NFAC)	PIXEL COOS AS NEEDED BY SRIPE
C    CI(NFAC)           IMAGE COOS
C    D(NFAC)            EIGENVALUES AS NEEDED BY SRIPE
C    PJ(I)		AVERAGE OF INDIVIDUAL IMAGE.
C    SDV(I)             STANDARD DEVIATION
C
C++*********************************************************************

       SUBROUTINE PCA_RECDUM(NIMA, NPIX, NF, NV, NFAC,
     &                       PIA, RIM, CP, CI, D, PJ, ADD_AVRG)

        DIMENSION :: RIM(NPIX), CP(NPIX,NFAC), CI(NFAC), D(NFAC)
        DIMENSION :: PJ(NPIX), NV(NF)

	LOGICAL   :: ADD_AVRG

        PIA = 0.0

        DO I = 1, NPIX
	  IF (ADD_AVRG) THEN
             RIM(I) = PJ(I)
          ELSE
             RIM(I) = 0
          ENDIF

          DO L = 1, NF
             K      = NV(L)
             RIM(I) = RIM(I) + CP(I,K) * CI(K) / SQRT(D(K))
          ENDDO
                                                
          PIA = PIA + RIM(I)

        ENDDO

        END

 
C++*********************************************************************
C
C  PCA_REMAKE(IMNO, NIMA, NPIX, NF, NV, NFAC,
C                 PIA, RIM, CP, CI, D, PJ, ADD_AVRG)
C
C PURPOSE: RECONSTITUTES IMAGES USING COORDINATES FROM _IMC & _PIX FILES
C
C PARAMETERS:
C    IMNO          IMAGE NUMBER
C    NIMA          NUMBER OF ACTIVE IMAGES
C    NPIX          NUMBER OF ACTIVE PIXELS
C    NFAC          NUMBER OF FACTORS USED IN THE ANALYSIS
C    PIA           AVERAGE OF RECONSTITUTED IMAGE (SUPPLIED)
C    RIM(NPIX)     OUTPUT IMAGE
C    CP(NPIX,NFAC) PIXEL COOS AS NEEDED BY SRIPE
C    CI(NFAC)      IMAGE COOS
C    D(NFAC)       EIGENVALUES AS NEEDED BY SRIPE
C    PJ(I)         AVERAGE OF INDIVIDUAL IMAGE.
C    ADD_AVRG      ADD THE AVERAGE TO IMAGE
C
C++*********************************************************************
  
        SUBROUTINE PCA_REMAKE(IMNO, NIMA, NPIX, NF, NV, NFAC,
     &                     PIA, RIM, CP, CI, D, PJ, ADD_AVRG)

        DIMENSION :: RIM(NPIX), CP(NPIX,NFAC), CI(NIMA,NFAC), D(NFAC)
        DIMENSION :: PJ(NPIX),  NV(NF)
        REAL      :: SUM
        LOGICAL   :: ADD_AVRG

        SUM = 0.0

C       LOOP OVER ALL PIXELS
        DO I = 1, NPIX

C	  REPLACE THE AVERAGES TO THE ORIGINAL DATA IF(PCA_AVRG)=TRUE
          IF (ADD_AVRG) THEN
             RIM(I) = PJ(I)
          ELSE
             RIM(I) = 0
          ENDIF

          DO L = 1, NF
             K      = NV(L)
#ifdef DEBUG
             if (i .lt. 1 .or. i .gt. npix) write(6,*) 'bad i:',i,l,k
             if (k .lt. 1 .or. k .gt. nfac) write(6,*) 'bad k:',k,l,i
#endif
             RIM(I) = RIM(I) + CP(I, K) * CI(IMNO, K) / SQRT(D(K))
          ENDDO

          SUM = SUM + RIM(I)
        ENDDO

        PIA = SUM / FLOAT(NPIX)

        END
 
 

