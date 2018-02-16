C++*********************************************************************
C
C  FOUR_FQ.F
C             OPFILEC                              FEB 03 ARDEAN LEITH    
C             COS                                  JUL 12 G. KISHCHENKO    
C             FQ Q PARAMS                          OCT 12 ARDEAN LEITH    
C             RAISED SINC                          OCT 13 ARDEAN LEITH    
C             CREATED FROM FOUR1A                  NOV 14 ARDEAN LEITH    
C             VOLUME IOPT 12 & 13 TRAP for ifort   JUN 15 ARDEAN LEITH    
C             VOLUME FILTER BUG                    MAY 17 ARDEAN LEITH    
C             SIMPLIFIED                           JAN 18 ARDEAN LEITH    
C             SQRT2M1 in BUTTERWORTH               JAN 18 ARDEAN LEITH
C       
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2018  Health Research Inc.,                         *
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
C FOUR_FQ 
c
C PURPOSE: APPLIES FOURIER FILTERS TO 2-D OR 3-D REAL PICTURES
C            'FQ'    : QUICK FILTERING (IN CORE, 2-D OR 3-D)
C            'FQ NP' : QUICK FILTERING (IN CORE, 2-D OR 3-D) NO PAD
C
C NOTE:  BUTTERWORTH FILTERS APPEAR TO BE OFFSET FROM MY EXPECTED
C        CENTRAL LINE BETWEEN PASS BAND AND STOP BAND BUT THIS GOES
C        WAY BACK BEFOR SPIDER 17.0 AND MAYBE FOREVER.  I SUGGEST
C        NOT USING IT??
C
C NOTE:  APPEARS TO HAVE UNDOCUMENTED AND UNTESTED ELLIPTICAL 
C        FILTRATION FOR OPTIONS: 1,2,3,4.   THIS WILL GIVE BAD 
C        ERRORS IF PERSON ENTERS MORE THAN ONE VALUE ON THE INPUT
C        PARAMETER LINE AS IT IS USED FOR ELLIPSES. al nov 2014
C
C        SUBSTITUTION OF PRECALCULATED INTERMEDIATE VARIABLES NOT
C        USEFULL FOR SPEED SINCE TIME IS DOMINATED BY FFT. SEE RCS
C        VERSIONS FOR TESTING THIS
C
C        'FF' (ffilts.f)     SETS XD2SQ TO: (NX  / 2)**2 BUT
C        'FQ' (four_fq.f)    SETS XD2SQ TO: (NXF / 2)**2 WHERE
C             NX  IS X DIMENSION OF POSSIBLY PADDED IMAGE
C             NXF IS SLIGHTLY LARGER DUE TO MIXED RADIX FOURIER PAD
C        SO THEY GIVE SLIGHTLY DIFFERENT RESULTS.  I SUSPECT THAT
C        'FF' IS ACTUALLY CORRECT?
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE FOUR_FQ

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC'

        CHARACTER (LEN=MAXNAM) :: FILNAM

        REAL, ALLOCATABLE      :: BB(:,:,:)
       
        INTEGER                :: N2X,N2Y,N2Z,LSD
        INTEGER                :: NX,NY,NZ
        INTEGER                :: MAXIM,IOPT,IRTFLG,NOT_USED
        INTEGER                :: IDUM,NE,IPADTYPE
        INTEGER                :: NMAX,NGOT,INV,I,J,K,IRECT
        REAL                   :: ORD,FP,FS,PARM1T,PARM2T,AVET
        REAL                   :: PARM1,PARM2,PARM3,PADVAL
        REAL                   :: BFPS(4)
        INTEGER                :: IWANTFOU,NLETF,ITYPE
        LOGICAL                :: WANTFOU
        CHARACTER (LEN=1)      :: NULL = CHAR(0)

        REAL, PARAMETER        :: EPS = 0.882
        REAL, PARAMETER        :: AA  = 10.624

        INTEGER, PARAMETER     :: LUN1 = 21
        INTEGER, PARAMETER     :: LUN2 = 22

C       OPEN REAL OR FOURIER SPACE INPUT FILE
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NX,NY,NZ,
     &             MAXIM,'INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        AVET = AV      ! FROM CMBLOCK

        IF (IFORM .NE. 1 .AND. IFORM .NE. 3) THEN
           WRITE(NOUT,*) '  OPERATION PROCESSING FOURIER INPUT'
           CALL ERRT(101,'OPERATION DOES NOT ACCEPT FOURIER INPUT',IDUM)
           GOTO 999

        ELSEIF (IFORM == 1 .AND. NY == 1) THEN
           CALL ERRT(101,'OPERATION DOES NOT WORK ON 1D INPUT',IDUM)
           GOTO 999
        ENDIF
        
        CALL FILERD(FILNAM,NLETF,NULL,'OUTPUT',IRTFLG)
        IF (IRTFLG. NE. 0) GOTO 999
 
 1000   WRITE(NOUT,90)
 90     FORMAT
     &  ('  1 - LOW-PASS,                2 - HIGH-PASS'         ,/,
     &   '  3 - GAUSS.  LOW-PASS,        4 - GAUSS.  HIGH-PASS' ,/,
     &   '  5 - FERMI   LOW-PASS         6 - FERMI   HIGH-PASS' ,/,
     &   '  7 - BUTTER. LOW-PASS,        8 - BUTTER. HIGH-PASS' ,/,
     &   '  9 - RAISED COS. LOW-PASS,   10 - RAISED COS. HIGH-PASS' ,/,
     &   ' 13 - RAISED SINC WINDOW,     14 - B FACTOR')

        IWANTFOU = 0    ! HIDDEN FOURIER OUTPUT OPTION
        CALL RDPRI2S(IOPT,IWANTFOU,NOT_USED,
     &               'FILTER TYPE (1-10,13,14)',IRTFLG)
        IF (IRTFLG. NE. 0) GOTO 999
        
        IF (IOPT < 1 .OR. IOPT > 14) THEN
           CALL ERRT(102,'ILLEGAL VALUE FOR FILTER TYPE',IOPT)
           GOTO 1000
        ENDIF

       IF (NZ > 1 .AND. IOPT == 13 ) THEN
           CALL ERRT(101,'OPTION NOT IMPLEMENTED FOR VOLUMES',IOPT)
           GOTO 1000
        ENDIF
         
        WANTFOU = (IWANTFOU == 1)

        PARM2 = 0.0
        PARM3 = 0.0
        BFPS  = 0.0    ! ARRAY OP

        IF (IOPT == 5 .OR. IOPT == 6)  THEN
C          FERMI DISTRIBUTION FILTER ****************************

           PARM1 = 0.25
           CALL RDPRM1S(PARM1,NOT_USED,
     &        'FILTER RADIUS (IN FREQUENCY OR PIXEL UNITS)',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (PARM1 < 0.0 .OR. PARM1 > 1.0) 
     &         PARM1 = PARM1 / NX     ! RADIUS INPUT

           CALL RDPRM1S(PARM3,NOT_USED,
     &                     'TEMPERATURE (0=CUTOFF)',IRTFLG)

C          EXPONENTIAL FOR HIGH-PASS OPTION
           IF (IOPT == 6) PARM3 = -PARM3
 
        ELSEIF (IOPT == 7  .OR. IOPT == 8) THEN
C          BUTTERWORTH  FILTER ********************************

           NMAX = 4
           CALL RDPRA(
     &        'LOWER & UPPER LIMITING FREQ. (IN FREQ OR PIXEL UNITS)',
     &        NMAX,0,.FALSE.,BFPS,NGOT,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (BFPS(1) > 1.0) THEN
C             RADIUS INPUT
              PARM1 = BFPS(1) / NX
              PARM2 = BFPS(2) / NX
           ELSE
C             FREQUENCY INPUT
              PARM1 = BFPS(1)
              PARM2 = BFPS(2)
           ENDIF

           IF (BFPS(3) == 0.0 .AND. BFPS(4) == 0.0) THEN
C             BUTTERWORTH CIRCULAR FILTER

              ORD   = 2.0 * ALOG10(EPS / SQRT(AA**2 - 1.0) )

              ORD   = ORD   / ALOG10(PARM1 / PARM2)
              PARM1 = PARM1 / EPS**(2. / ORD)

              PARM3 = ORD

              WRITE(NOUT,*)
     &         ' BUTTERWORTH CIRCULAR FILTER, ORD: ',ORD

           ELSE
C             BUTTERWORTH ELLIPTIC FILTER
C             LOW-PASS  IOPT=11,  HIGH-PASS IOPT=12

              IOPT = IOPT + 4
              IF (NZ > 1 ) THEN
                 CALL ERRT(101,
     &              'OPTION NOT IMPLEMENTED FOR VOLUMES',IOPT)
                 GOTO 1000
              ENDIF

              ORD   = 2.0 * ALOG10(EPS / SQRT(AA**2 - 1.0) )
              WRITE(NOUT,*)
     &         ' BUTTERWORTH ELIPTICAL FILTER, ORD: ',ORD
           ENDIF

        ELSEIF ( IOPT == 9  .OR. IOPT == 10)  THEN
C          RAISED COSINE FILTER ********************************

           PARM1 = 0.25
           PARM2 = -9999999
           CALL RDPRM2S(PARM1,PARM2,NOT_USED,
     &       'LOWER & UPPER LIMITING FREQ. (IN FREQ OR PIXEL UNITS)',
     &        IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (PARM1 < 0.0 .OR. PARM1 > 1.0) 
     &         PARM1 = PARM1 / NX     ! RADIUS INPUT

           IF (PARM2 == -9999999) THEN
              ! SAME AS PARM1
              PARM2 = PARM1

           ELSEIF  (PARM2 < 0.0 .OR. PARM2 > 1.0) THEN
              ! RADIUS INPUT
              PARM2 = PARM2 / NY   
           ENDIF

        ELSEIF (IOPT == 14)  THEN
C          B FACTOR FILTER  ************************************

           WRITE(NOUT,*)
     &         ' NORMALIZES AMPLITUDES BY A "B" TEMPERATURE FACTOR'
           WRITE(NOUT,*)' AMP = AMP*D*(EXP(B*F**2))'

           PARM1 = -200
           CALL RDPRM1S(PARM1, NOT_USED,'B FACTOR (PIXEL**2)',   IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           PARM2 = 1.0
           CALL RDPRM1S(PARM2, NOT_USED,'D MULTIPLIER', IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           CALL RDPRM1S(PARM3, NOT_USED,'FREQUENCY CUTOFF',  IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
           IF (PARM3 < 0.0 .OR. PARM3 > 0.5) THEN
              WRITE(NOUT,'(A)') '  Allowed Cutoff Range: 0 .. 0.5'
           ENDIF

           PARM3 = PARM3**2
         

        ELSE
C          ALL OTHER FILTERS *************************************
C          LOW-PASS,  HIGH-PASS,  GAUSS.LOW-PASS,  GAUSS.HIGH-PASS
C          RAISED SINC WINDOW 

           PARM1T = 0.25
           PARM2T = -9999999
           CALL RDPRM2S(PARM1T,PARM2T,NOT_USED,
     &        'FILTER RADIUS (IN FREQUENCY OR PIXEL UNITS)',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           PARM1 = PARM1T
           IF (PARM1T < 0.0 .OR. PARM1T > 1.0) THEN
              PARM1 = PARM1 / NX     ! RADIUS INPUT
              WRITE(NOUT,'(A,F8.4)') '  Frequency: ',PARM1
           ENDIF

           IF (PARM2T == -9999999) THEN
              ! DEFAULTS PARM2 SAME AS PARM1
              PARM2 = PARM1

           ELSEIF  (PARM2T < 0.0 .OR. PARM2T > 1.0) THEN
C             RADIUS INPUT FOR PARM2
              PARM2 = PARM2T / NY   
           ENDIF

        ENDIF


C       IPADTYPE: 0 == NONE,  1== ZEROPAD, 2==BORDERPAD, 3==AVG PAD
        IPADTYPE = 2    ! PAD WITH ZEROES

        IF ( INDEX(FCHAR,'NP') > 0)  THEN
C          NO PADDING
           IPADTYPE = 0    ! NO PADDING
           N2X      = NX
           N2Y      = NY
           N2Z      = NZ

        ELSE
C          2x PADDING
           N2X = 2 * NX
           IF (NY > 1)  THEN
              N2Y = 2 * NY
           ELSE
              N2Y = 1
           ENDIF
           IF (NZ > 1)  THEN
              N2Z = 2 * NZ
           ELSE
              N2Z = 1
           ENDIF
        ENDIF

        LSD = N2X + 2 - MOD(N2X,2)
         
        ALLOCATE (BB(LSD,N2Y,N2Z), STAT=IRTFLG)           
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'FOUR_FQ; BB',LSD*N2Y*N2Z)
           GOTO 999
        ENDIF


        IF (IPADTYPE == 0) THEN

C          NO PADDING AT ALL, EXCEPT FOR FFT PAD IN BUFFER
           PADVAL = 0

C          LOAD IMAGE/VOLUME INTO FFT PADDED BUFFER
           CALL REDNPADVOL(LUN1,PADVAL, 
     &                  NX, NY, NZ, LSD,N2Y,N2Z,
     &                  BB,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 999 

        ELSEIF (IPADTYPE == 1) THEN
C          USE ZERO AS PADDING VALUE
           PADVAL = 0.0

C          LOAD IMAGE/VOLUME INTO FFT PADDED BUFFER
           CALL REDNPADVOL(LUN1,PADVAL, 
     &                  NX, NY, NZ, LSD,N2Y,N2Z,
     &                  BB,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 999 

        ELSEIF (IPADTYPE == 3) THEN
C          USE AVE FROM IMAGE AS PADDING VALUE
           PADVAL = AVET

C          LOAD IMAGE/VOLUME INTO FFT PADDED BUFFER
           CALL REDNPADVOL(LUN1,PADVAL, 
     &                  NX, NY, NZ, LSD,N2Y,N2Z,
     &                  BB,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 999 
 
        ELSEIF (IPADTYPE == 2) THEN
C          USE AVERAGE OF BORDER VALUES AS PADDING VALUE

C          LOAD IMAGE/VOLUME INTO FFT PADDED BUFFER
           CALL READV(LUN1,BB, LSD,N2Y, NX,NY,NZ)

C          DETERMINE PADDING VALUE & PAD IMAGE/VOLUME
           IF (IFORM == 1) THEN
              PADVAL = 
     &          (SUM(BB(1:NX,1,1))   + SUM(BB(1:NX,NY,1)) +
     &           SUM(BB(1,2:NY-1,1)) + SUM(BB(NX,2:NY-1,1)) ) / 
     &           REAL(2*(NX+NY)-4)

C             PAD IMAGE
c$omp         parallel do private(j)
	      DO J=1,N2Y
	         BB(NX+1:N2X,J,1) = PADVAL
	      ENDDO

c$omp         parallel do private(j)
	      DO J=NY+1,N2Y
	         BB(1:NX,J,1) = PADVAL
	      ENDDO

           ELSE
              PADVAL =
     &          (SUM(BB(1:NX,1:NY,1))     + SUM(BB(1:NX,1:NY,NZ))   +
     &           SUM(BB(1:NX,1,2:NZ-1))   + SUM(BB(1:NX,NY,2:NZ-1)) +
     &           SUM(BB(1,2:NY-1,2:NZ-1)) + SUM(BB(NX,2:NY-1,2:NZ-1)))
     &         / REAL( 4 * (NX+NY+NZ) - 16)

C             PAD VOLUME
c$omp         parallel do private(k)
	      DO K = 1,NZ
 	         BB(NX+1:N2X, 1:NY,     K) = PADVAL
 	         BB(1:N2X,    NY+1:N2Y, K) = PADVAL
              ENDDO

c$omp         parallel do private(k)
	      DO K=NZ+1,N2Z
	         BB(1:N2X, 1:N2Y, K) = PADVAL
	      ENDDO
	   ENDIF

        ENDIF

        !write(6,*) ' pad:',padval

C       FORWARD FFT, FFT PADDED
        INV = 1
        CALL FMRS(BB,N2X,N2Y,N2Z, 0.0D0, .TRUE.,.TRUE.,INV,IRTFLG)
        IF (INV == 0) THEN
           IRTFLG = 1
           GOTO 999
        ENDIF

        IF (IFORM == 1)  THEN
C          IMAGE
           WRITE(NOUT,"(A,I0,A,I0)")
     &                '  Dimensions used: ',N2X,' x ',N2Y

           CALL FQ2(IOPT, PARM1,PARM2,PARM3,BFPS, BB, 
     &              LSD,N2X,N2Y, NX,NY, IRTFLG)
        ELSE
C          VOLUME
           WRITE(NOUT,"(A,I0,A,I0,A,I0)") 
     &           '  Dimensions used: ',N2X,' x ',N2Y,' x ',N2Z

           CALL FQ3(IOPT, PARM1,PARM2,PARM3,BFPS, BB, 
     &              LSD,N2X,N2Y,N2Z, NX,NY,NZ, IRTFLG)
        ENDIF

C       OPEN OUTPUT FILE
        MAXIM = 0
        ITYPE = IFORM

        IF (WANTFOU) THEN
           IF (IFORM == 1)  THEN
C             NON-FOURIER IMAGE
              ITYPE = -11
              IF (MOD(NX,2) == 0)   ITYPE = -12
           ELSEIF (IFORM == 3)  THEN
C             NON-FOURIER VOLUME
              ITYPE = -21
              IF (MOD(NX,2) == 0)  ITYPE = -22
           ENDIF

           CALL OPFILEC(LUN1,.FALSE.,FILNAM,LUN2,'U',ITYPE,N2X,N2Y,N2Z,
     &                  MAXIM,' ',.TRUE.,IRTFLG)
           IF (IRTFLG. NE. 0) GOTO 999

C          STORE FILTERED FFT
           CALL WRTVOL(LUN2,N2X,N2Y,1,N2Z,BB,IRTFLG) 
           GOTO 999

        ENDIF

        CALL OPFILEC(LUN1,.FALSE.,FILNAM,LUN2,'U',ITYPE,NX,NY,NZ,
     &               MAXIM,' ',.FALSE.,IRTFLG)
        IF (IRTFLG. NE. 0) GOTO 999


C       INVERSE FFT,  POSSIBLY PADDED
        INV = -1         ! INVERSE FFT
        CALL FMRS(BB, N2X,N2Y,N2Z, 0.0D0, .TRUE.,.TRUE.,INV,IRTFLG)
        IF (INV == 0) THEN
           CALL ERRT(101,'USING FFTW',NE)
           GOTO 999 
        ENDIF

        !write(6,*) ' sum3:',sum(bb(1:nx, 1:ny, 1))

C       STORE FILTERED IMAGE
        IRECT = 1
        DO K=1,NZ
           DO J=1,NY
 	      CALL WRTLIN(LUN2,BB(1,J,K),NX,IRECT)
              IRECT = IRECT + 1
           ENDDO
	ENDDO

999     IF (ALLOCATED(BB)) DEALLOCATE(BB)

        CLOSE(LUN1)
        CLOSE(LUN2)
       
        END



C++*********************************************************************
C
C FQ2.F                                         12/22/94  
C                RDPRAF REMOVED                 DEC 2005 ARDEAN LEITH 
C                COSINE FILTER                  JUL 2012 G. KISHCHENKO 
C                INCORE WITHOUT LUNS            OCT 2012 ARDEAN LEITH 
C                FREQ + PIXELS                  NOV 2012 G. KISHCHENKO 
C                PARM2 BUG                      DEC 2012 ARDEAN LEITH 
C                FREQ UNIT CUTOFF = 1           AUG 2013 ARDEAN LEITH 
C                RAISED SINC                    FEB 2014 ARDEAN LEITH 
C                FROM FQ_Q                      NOV 2014 ARDEAN LEITH 
C        
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2018, Health Research Inc.,                         *
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
C FQ2(IOPT,PARM1,PARM2,PARM3,BFPS,B, LSD,N2X,N2Y, NX,NY, IRTFLG)
C
C PURPOSE: QUICK FILTERING OF REAL-SPACE IMAGE FILE BY FFT
C
C PARAMETERS:
C        IOPT          TYPE OF FILTER
C        PARM1..       FILTER PARAMETERS
C        BFPS          FILTER PARAMETERS
C        B             BUFFER 
C        LSD           X DIMENSION OF FOURIER FILE
C        NX,NY         DIMENSIONS OF REAL-SPACE FILE
C        N2X,N2Y       DIMENSIONS OF PADDED FILE
C
C23456789012345678901234567890123456789012345678901234567890123456789012  
C--*******************************************************************

        SUBROUTINE FQ2(IOPT, PARM1,PARM2,PARM3,BFPS, 
     &                  B, LSD,N2X,N2Y, NX,NY, IRTFLG)
        
        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'

        INTEGER            :: IOPT      
        REAL               :: PARM1,PARM2,PARM3
        REAL               :: BFPS(4)
        REAL               :: B(LSD,N2Y)
        INTEGER            :: LSD,N2X,N2Y
        INTEGER            :: NX,NY,IRTFLG

        REAL               :: ZEROTERM
        REAL               :: ORD,PARM1SQ,PARM2SQ,SQRT2M1

        REAL               :: XD2SQ,YD2SQ,FSQDP,FSQ,FSQD4
        REAL               :: F,F2,FPE,FSE,ORDT,PARMT,EXPCORR
        INTEGER            :: J,I,IX,IY, NR2
        
        REAL, PARAMETER    :: PI  = 3.14159265358979323846
        REAL, PARAMETER    :: EPS = 0.882

        real       :: bsavi(lsd),bsavv(lsd),bsavf(lsd)
        integer    :: it
        real       :: xd2sqdp1sqinv,xd2sqdp1sq

        XD2SQ   = FLOAT(N2X/2)**2

        NR2     = N2Y / 2
        YD2SQ   = FLOAT(NR2)**2

        ORD     = PARM3       ! FOR BUTTERWORTH FILTER
        PARM1SQ = PARM1**2
        PARM2SQ = PARM2**2

        SQRT2M1 = (SQRT(2.0)-1)

C       KEEP ZERO TERM FOR HIGH PASS OPTIONS
        ZEROTERM = B(1,1)
            
        !xd2sqdp1sq    =   xd2sq/parm1sq
        !xd2sqdp1sqinv = 1/xd2sq/parm1sq

        !write(6,*) ' p1,p2,p3:        ',        parm1,parm2,parm3
        !write(6,*) ' p1sq,p2sq:       ',        parm1sq,parm2sq
        !write(6,*) ' nx,xd2sq,n2x:    ',        nx,xd2sq,n2x
        !write(6,*) ' xd2sqdp1sqinv,zeroterm: ', xd2sqdp1sqinv,zeroterm

        !write(6,*) ' sum1:',sum(b(1:n2x, 1:n2y)),zeroterm

c$omp  parallel do private(j,i,iy,ix,it,f,f2,fpe,fse,ordt,parmt,
c$omp&                     fsq,fsqdp,fsqd4,expcorr)
        DO J=1,N2Y

           IY = (J-1)
           IF (IY > NR2) IY = IY - N2Y


           DO I=1,LSD,2

              IX    = (I-1) / 2
              FSQ   = (FLOAT(IX*IX)/XD2SQ +
     &                 FLOAT(IY*IY)/YD2SQ)

              FSQDP = (FLOAT(IX*IX)/XD2SQ/PARM1SQ +
     &                 FLOAT(IY*IY)/YD2SQ/PARM2SQ)

              !it        = (i -1)/2 + 1    ! 1,1,2
              !f         = float(ix)/n2x + float(iy)/n2y
              !bsavf(it) = f
              !write(6,*) ' it,ix,iy,f:',it,ix,iy,f

              IF (IOPT == 1) THEN
C                LOWPASS *************************************

                 !bsavv(it) =  1.0
                 IF (0.25 * FSQDP > 1.0) THEN

                    B(I,  J)  = 0.0
                    B(I+1,J)  = 0.0

                   !bsavv(it) = 0.0
                 ENDIF

              ELSEIF (IOPT == 2) THEN   
C                HIGH PASS ***********************************

                 !bsavv(it) =  1.0
 
                 IF ((IX.NE.0 .OR. IY.NE.0) .AND. 
     &               (0.25*FSQDP) <= 1.0) THEN

                    B(I,  J)  = 0.0
                    B(I+1,J)  = 0.0
 
                   !bsavv(it) = 0.0
                 ENDIF

              ELSEIF(IOPT == 3)  THEN
C                GAUSSIAN LOW PASS ***************************

                 F = 0.125 * FSQDP

                 !bsavi(it) = F

                 IF (F < 16.0)  THEN
                    F         = EXP(-F)
                    B(I,  J)  = B(I,  J) * F
                    B(I+1,J)  = B(I+1,J) * F
                   !bsavv(it) = F
                 ELSE
                    B(I,  J)  = 0.0
                    B(I+1,J)  = 0.0
                   !bsavv(it) = 0.0
                 ENDIF

              ELSEIF (IOPT==4)  THEN    
C                GAUSSIAN HIGH PASS **************************

                 IF (IX .NE. 0 .OR. IY .NE. 0)  THEN
                    F = 0.125 * FSQDP

                   !bsavv(it) = 1.0
                   !bsavi(it) = F

                    IF (F < 16.0)  THEN
                       F         = 1.0 - EXP(-F)
                       B(I,  J)  = B(I,  J) * F
                       B(I+1,J)  = B(I+1,J) * F

                      !bsavv(it) = F
                    ENDIF
                 ENDIF

              ELSEIF (IOPT == 5 .OR. IOPT == 6)  THEN
C                FERMI DISTRIBUTION FILTER *******************
              
                 F         = (0.5 * SQRT(FSQ) - PARM1) / PARM3

                 F         = MIN(MAX(F,-10.0), 10.0)
                 F         = (1.0/(1.0 + EXP(F)))

                 B(I,  J)  = B(I,  J) * F
                 B(I+1,J)  = B(I+1,J) * F

                 !bsavi(it) = (0.5 * sqrt(fsq) - parm1) / parm3
                 !bsavv(it) = f

              ELSEIF (IOPT == 7) THEN
C                BUTTERWORTH LOWPASS FILTER ******************

                 F = 0.5 * SQRT(FSQ)

                !F = SQRT(1.0 / (1.0 +           (F/PARM1)**ORD)) ! al jan 2018
                 F = SQRT(1.0 / (1.0 + SQRT2M1 * (F/PARM1)**ORD)) 

                 B(I,J)    = B(I,  J) * F
                 B(I+1,J)  = B(I+1,J) * F

                 !bsavi(it) = 0.5 * sqrt(fsq)
                 !bsavv(it) = f

              ELSEIF (IOPT == 8) THEN
C                BUTTERWORTH HIGHPASS FILTER *****************
        
                 !bsavi(it) = 0

                 IF (IX.NE.0 .OR. IY.NE. 0) THEN
                    F = 0.5 * SQRT(FSQ)

                   !F  = (1.0 - SQRT(1.0 / (1.0 +(F/PARM1)**ORD)))
                    F  = (1.0 - SQRT(1.0 / 
     &                    (1.0 + SQRT2M1 *(F/PARM1)**ORD))) ! al jan 2018

                    B(I,  J)  = B(I,  J) * F
                    B(I+1,J)  = B(I+1,J) * F

                   !bsavi(it) = 0.5 * sqrt(fsq)
                   !bsavv(it) = f
                 ENDIF


              ELSEIF (IOPT == 9) THEN
C                RAISED COSINE LOWPASS FILTER ******************

                 F = 0.5 * SQRT(FSQ)
 
                 F = (F - PARM1) / (PARM2 - PARM1)

                 IF (F < 0) THEN
                    F2 = 1
                 ELSEIF (F > 1) THEN
                    F2 = 0
                 ELSE
                    F2 = 0.5 * (1 + COS(PI*F))
                 ENDIF

                 B(I,  J)  = B(I,  J) * F2
                 B(I+1,J)  = B(I+1,J) * F2
                 !bsavv(it) = F2

              ELSEIF (IOPT == 10) THEN
C                RAISED COSINE HIGHPASS FILTER ******************

                 F = 0.5 * SQRT(FSQ)

                 F = (F - PARM1) / (PARM2 - PARM1)

                 IF (F < 0) THEN
                    F2 = 0
                 ELSEIF (F > 1) THEN
                    F2 = 1
                 ELSE
                    F2 = 0.5 * (1 - COS(PI*F))
                 ENDIF

                 B(I,  J)  = B(I,  J) * F2
                 B(I+1,J)  = B(I+1,J) * F2
                 !bsavv(it) = F2

              ELSEIF (IOPT == 11) THEN 
C                BUTTERWORTH ELLIPTIC LOWPASS FILTER *********  
C                CALCULATE EFFECTIVE FP AND FS IN A GIVEN 
C                DIRECTION ON THE PLANE

                 IF (IX.NE.0 .OR. IY.NE.0) THEN

                    FPE = ATAN2(BFPS(1)*SQRT(FLOAT(IY*IY)/YD2SQ),
     &                          BFPS(3)*SQRT(FLOAT(IX*IX)/XD2SQ))
                    FPE = SQRT((BFPS(1)*COS(FPE))**2 + 
     &                         (BFPS(3)*SIN(FPE))**2)

                    FSE = ATAN2(BFPS(2)*SQRT(FLOAT(IY*IY)/YD2SQ),
     &                          BFPS(4)*SQRT(FLOAT(IX*IX)/XD2SQ))
                    FSE = SQRT((BFPS(2)*COS(FSE))**2 + 
     &                         (BFPS(4)*SIN(FSE))**2)

                    ORDT      = ORD  /ALOG10(FPE / FSE)
                    PARMT     = FPE/ EPS**(2. / ORDT)

                    F         = 0.5*SQRT(FSQ)
                    F         = SQRT(1.0 / (1.0 + (F/PARMT)**ORDT))

                    B(I,  J)  = B(I,  J) * F
                    B(I+1,J)  = B(I+1,J) * F
                    !bsavv(it) = F
                ENDIF

              ELSEIF (IOPT == 12) THEN
C                BUTTERWORTH ELLIPTIC HIGHPASS FILTER ********* 

                 IF (IX .NE. 0 .OR. IY.NE. 0) THEN

                    FPE = ATAN2(BFPS(1)*SQRT(FLOAT(IY*IY)/YD2SQ),
     &                          BFPS(3)*SQRT(FLOAT(IX*IX)/XD2SQ))
                    FPE = SQRT((BFPS(1)*COS(FPE))**2 +
     &                         (BFPS(3)*SIN(FPE))**2)

                    FSE = ATAN2(BFPS(2)*SQRT(FLOAT(IY*IY)/YD2SQ),
     &                          BFPS(4)*SQRT(FLOAT(IX*IX)/XD2SQ))
                    FSE = SQRT((BFPS(2)*COS(FSE))**2 + 
     &                         (BFPS(4)*SIN(FSE))**2)

                    ORDT     = ORD / ALOG10(FPE / FSE)
                    PARMT    = FPE / (EPS)**(2. / ORDT)

                    F        = 0.5*SQRT(FSQ)
                    F        = (1.0-SQRT(1.0/(1.0+(F/PARMT)**ORDT)))

                    B(I,  J)  = B(I,  J) * F
                    B(I+1,J)  = B(I+1,J) * F

                    !bsavv(it) = F
                 ENDIF

              ELSEIF (IOPT == 13) THEN
C                RAISED SINC WINDOW **************************
                 F = 0.5 * SQRT(FSQDP)

                 IF (F <= 0.0001) THEN
                    F2 = 1
                 ELSEIF (F >= 1.0) THEN
                    F2 = 0
                 ELSE
                    F2 = SIN(PI*F) / (PI*F)
                 ENDIF

                 B(I,  J)  = B(I,  J) * (1 + 9 * F2)
                 B(I+1,J)  = B(I+1,J) * (1 + 9 * F2)

                !bsavv(it) = (1 + 9 * F2)

             ELSEIF (IOPT == 14) THEN
C                B FACTOR FILTER **************************

                 FSQD4     = FSQ * 0.25

                 !bsavi(it) = fsqd4
                 !bsavv(it) = 1.0

                 IF (FSQD4 < PARM3)  THEN
                    EXPCORR  = EXP(FSQD4 * PARM1)

                    B(I,  J)  = B(I,  J) * PARM2 * EXPCORR
                    B(I+1,J)  = B(I+1,J) * PARM2 * EXPCORR
                   !bsavv(it) = PARM2 * EXPCORR

                 ENDIF
              ENDIF
           ENDDO

#ifdef NEVER
           if (j ==1) then
             do it = 1,33
               write(6,888) bsavi(it), bsavf(it), bsavv(it)
888            format (f12.4, f12.4, f12.4)
             enddo
             stop
           endif
#endif

        ENDDO

C       RESTORE ZERO TERM FOR HIGH PASS OPTIONS
        IF (IOPT == 2 .OR. IOPT == 4 .OR. 
     &      IOPT == 6 .OR. IOPT == 8)     B(1,1) = ZEROTERM

        !write(6,*) ' sum2:',sum(b(1:n2x, 1:n2y))

        IRTFLG = 0

        END








C++*********************************************************************
C
C FQ_3.F                                        12/22/94  
C                RDPRAF REMOVED                 DEC 2005 ARDEAN LEITH 
C                COSINE FILTER                  JUL 2012 G. KISHCHENKO 
C                INCORE WITHOUT LUNS            OCT 2012 ARDEAN LEITH 
C                FREQ + PIXELS                  NOV 2012 G. KISHCHENKO 
C                PARM2 BUG                      DEC 2012 ARDEAN LEITH 
C                FREQ UNIT CUTOFF = 1           AUG 2013 ARDEAN LEITH 
C                RAISED SINC                    FEB 2014 ARDEAN LEITH 
C                FROM FQ_Q                      NOV 2014 ARDEAN LEITH 
C                VOLUME ZD2SQ FILTER BUG           MAY 2017 ARDEAN LEITH    
C        
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2018, Health Research Inc.,                         *
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
C FQ_3(IOPT,PARM1,PARM2,PARM3,BFPS, B, LSD,N2X,N2Y, NX,NY, IRTFLG)
C
C PURPOSE: QUICK FILTERING OF REAL-SPACE IMAGE FILE BY FFT
C
C PARAMETERS:
C        IOPT          TYPE OF FILTER
C        PARM1..       FILTER PARAMETERS
C        BFPS          FILTER PARAMETERS
C        B             BUFFER 
C        LSD           X DIMENSION OF FOURIER FILE
C        NX,NY,NZ      DIMENSIONS OF REAL-SPACE FILE
C        N2X,N2Y,N2Z   DIMENSIONS OF PADDED FILE
C
C23456789012345678901234567890123456789012345678901234567890123456789012  
C--*******************************************************************

        SUBROUTINE FQ3(IOPT, PARM1,PARM2,PARM3,BFPS, 
     &                  B, LSD,N2X,N2Y,N2Z, NX,NY,NZ, IRTFLG)
        
        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'

        INTEGER            :: IOPT       
        REAL               :: PARM1,PARM2,PARM3
        REAL               :: BFPS(4)
        REAL               :: B(LSD,N2Y,N2Z)
        INTEGER            :: LSD,N2X,N2Y,N2Z
        INTEGER            :: NX,NY,NZ,IRTFLG

        REAL               :: ORD,PARM1SQ,PARM2SQ,SQRT2M1

        REAL               :: XD2SQ,YD2SQ,ZD2SQ,FSQDP,FSQD4
        REAL               :: F,F2,FPE,FSE,ORDT,PARMT,FSQ,EXPCORR
        REAL               :: ZEROTERM
        INTEGER            :: K,J,I, IX,IY,IZ, NR2,NL2
        
        REAL, PARAMETER    :: PI  = 3.14159265358979323846
        REAL, PARAMETER    :: EPS = 0.882

        NR2     = N2Y / 2
        NL2     = N2Z / 2

        XD2SQ      = FLOAT(N2X/2)**2
        YD2SQ      = FLOAT(NR2)  **2
        ZD2SQ      = FLOAT(NL2)  **2

        ORD     = PARM3
        PARM1SQ = PARM1**2
        PARM2SQ = PARM2**2
      
        SQRT2M1 = (SQRT(2.0)-1)


        IF (IOPT == 11 .OR. IOPT == 12) THEN 
C          BUTTERWORTH ELLIPTIC LOWPASS FILTER *********       
C          BUTTERWORTH ELLIPTIC HIGHPASS FILTER *********      
           CALL ERRT(101,
     &        'OPTION NOT IMPLEMENTED FOR VOLUMES',IOPT)
           RETURN
        ENDIF

C       KEEP ZERO TERM FOR HIGH PASS OPTIONS
        ZEROTERM = B(1,1,1)

c$omp   parallel do private(i,j,k,ix,iy,iz,f,f2,fpe,fse,ordt,parmt,
c$omp&                      fsq,fsqdp,fsqd4,expcorr)
      
        DO K=1,N2Z

           IZ = (K-1)
           IF (IZ > NR2) IZ = IZ - N2Z

           DO J=1,N2Y
              IY = (J-1)
              IF (IY > NR2) IY = IY - N2Y

              DO I=1,LSD,2
                 IX = (I-1) / 2
                 FSQDP = (FLOAT(IX*IX)/XD2SQ/PARM1SQ +
     &                    FLOAT(IY*IY)/YD2SQ/PARM2SQ +
     &                    FLOAT(IZ*IZ)/ZD2SQ/PARM2SQ)

                 FSQ   = (FLOAT(IX*IX)/XD2SQ +
     &                    FLOAT(IY*IY)/YD2SQ +
     &                    FLOAT(IZ*IZ)/ZD2SQ)

                 IF (IOPT == 1) THEN
C                   LOWPASS *************************************
                    IF (0.25 * FSQDP > 1.0) THEN

                       B(I,  J,K) = 0.0
                       B(I+1,J,K) = 0.0
                    ENDIF

                 ELSEIF (IOPT == 2) THEN        
C                   HIGH PASS ***********************************
                    IF ( (IX.NE.0 .OR. IY .NE.0 .OR. IZ .NE. 0) .AND.
     &                 0.25 * FSQDP <= 1.0) THEN

                       B(I,  J,K) = 0.0
                       B(I+1,J,K) = 0.0
                    ENDIF

                 ELSEIF(IOPT == 3)  THEN
C                   GAUSSIAN LOW PASS ***************************
                    F = 0.125 * FSQDP

                    IF (F < 16.0)  THEN
                       F          = EXP(-F)
                       B(I,  J,K) = B(I,  J,K) * F
                       B(I+1,J,K) = B(I+1,J,K) * F
                    ELSE
                       B(I,  J,K) = 0.0
                       B(I+1,J,K) = 0.0
                    ENDIF

                 ELSEIF (IOPT==4)  THEN 
C                   GAUSSIAN HIGH PASS **************************

                    IF (IX .NE. 0 .OR. IY .NE. 0 .OR. IZ .NE. 0)THEN
                       F = 0.125 * FSQDP

                       IF (F < 16.0)  THEN
                          F          = 1.0 - EXP(-F)
                          B(I,  J,K) = B(I,  J,K) * F
                          B(I+1,J,K) = B(I+1,J,K) * F
                       ENDIF
                    ENDIF

                 ELSEIF (IOPT == 5 .OR. IOPT == 6)  THEN
C                   FERMI DISTRIBUTION FILTER *******************
              
                    F = (0.5 * SQRT(FSQ)-PARM1) / PARM3

                    F          = MIN(MAX(F,-10.0), 10.0)
                    F          = 1.0 / (1.0 + EXP(F))

                    B(I,  J,K) = B(I,  J,K) * F
                    B(I+1,J,K) = B(I+1,J,K) * F

                 ELSEIF (IOPT == 7) THEN
C                   BUTTERWORTH LOWPASS FILTER ******************

                    F          = 0.5 * SQRT(FSQ)

                    F          = SQRT(1.0 / 
     &                           (1.0 + SQRT2M1 * (F / PARM1)**ORD) )

                    B(I,  J,K) = B(I,  J,K) * F
                    B(I+1,J,K) = B(I+1,J,K) * F

                 ELSEIF (IOPT == 8) THEN
C                   BUTTERWORTH HIGHPASS FILTER *****************
        
                    IF (IX.NE.0 .OR. IY.NE. 0 .OR. IZ.NE. 0) THEN
                       F = 0.5 * SQRT(FSQ)
   
                       F = 1.0 - SQRT(1.0 / 
     &                          (1.0 + SQRT2M1 * (F / PARM1)**ORD))

                       B(I,  J,K) = B(I,  J,K) * F
                       B(I+1,J,K) = B(I+1,J,K) * F
                    ENDIF


                 ELSEIF (IOPT == 9) THEN
C                   RAISED COSINE LOWPASS FILTER ******************

                    F = 0.5 * SQRT(FSQ)

C                   F = (F-FP)    / (FS-FP)
                    F = (F-PARM1) / (PARM2 - PARM1)

                    IF (F < 0) THEN
                       F2 = 1
                    ELSEIF (F > 1) THEN
                       F2 = 0
                    ELSE
                       F2 = 0.5 * (1 + COS(PI*F))
                    ENDIF

                    B(I,  J,K) = B(I,  J,K) * F2
                    B(I+1,J,K) = B(I+1,J,K) * F2

                 ELSEIF (IOPT == 10) THEN
C                   RAISED COSINE HIGHPASS FILTER ******************

                    F = 0.5 * SQRT(FSQ)

                    F = (F-PARM1) / (PARM2 - PARM1)

                    IF (F < 0) THEN
                       F2 = 0
                    ELSEIF (F > 1) THEN
                       F2 = 1
                    ELSE
                       F2 = 0.5 * ( 1 - COS(PI * F))
                    ENDIF

                    B(I,  J,K) = B(I,  J,K) * F2
                    B(I+1,J,K) = B(I+1,J,K) * F2

                 ELSEIF (IOPT == 13) THEN
C                   RAISED SINC WINDOW **************************

                    F = 0.5 * SQRT(FSQDP)

                    IF (F <= 0.0001) THEN
                       F2 = 1
                    ELSEIF (F >= 1.0) THEN
                       F2 = 0
                    ELSE
                       F2 = SIN(PI * F) / (PI * F)
                    ENDIF
   
                    B(I,  J,K) = B(I,  J,K) * (1 + 9 * F2)
                    B(I+1,J,K) = B(I+1,J,K) * (1 + 9 * F2)
   
                 ELSEIF (IOPT == 14) THEN
C                   B FACTOR FILTER *****************************

                    FSQD4 = FSQ * 0.25
      
                    IF (FSQD4 < PARM3)  THEN
                       EXPCORR    = EXP(F * PARM1)

                       B(I,J,K)   = B(I,J,K)  * PARM2 * EXPCORR
                       B(I+1,J,K) = B(I+1,J,K)* PARM2 * EXPCORR
                    ENDIF

                 ENDIF
              ENDDO
           ENDDO
        ENDDO

C       RESTORE ZERO TERM FOR HIGH PASS OPTIONS
        IF (IOPT == 2 .OR. IOPT == 4 .OR. 
     &      IOPT == 6 .OR. IOPT == 8) 
     &     B(1,1,1) = ZEROTERM

        IRTFLG = 0

        END




