C **********************************************************************
C  TFED.F  
C            IMPROPER PS BUG FIXED               SEP   10 ArDean Leith
C            PROMPTS IMPROVED                    MAY   12 ArDean Leith
C            CTF ED                              JUN   12 ArDean Leith
C            DEF. DOC FILE COMMENTS CHANGED      APR   13 ArDean Leith
C                                                            
C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2002, Z. Huang & P. A. Penczek                        *
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
C PURPOSE: DETERMINE DEFOCUS & ASTIGMATISM OF  
C          MICROGRAPH FROM 2D POWER SPECTRUM 
C            
C NOTE:    FULL OF BUGS. DEFOCUS IS OK FOR NON-ASTIGMATIC IMAGES
C          ASTIGMATISM IS TOTALLY BUGGY! THIS IS GARBAGE al mar 2013
C
C OPERATIONS:  TFED(), CTFED()
C
C **********************************************************************

        SUBROUTINE  CTFED()

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'F90ALLOC.INC'

        INTEGER, PARAMETER                :: NT=4
        REAL                              :: DLIST(6)
        INTEGER                           :: NUMBER
        REAL                              :: GETDEFOCUS,LAMBDA,KEV
        CHARACTER(LEN=1)                  :: NULL = CHAR(0)
        REAL, DIMENSION(:,:), ALLOCATABLE :: POW2,FNOI_OUT,BUFIN,BUFWIN
        REAL, DIMENSION(:,:), ALLOCATABLE :: BUFFFT,BUFPWS,BUFPOW
        CHARACTER(LEN=MAXNAM)             :: FILNAM,DOCNAM,FILPAT
        CHARACTER(LEN=90)                 :: COMMENT
        CHARACTER(LEN=60)                 :: FORMOUT

        LOGICAL                           :: ADDEXT,GETNAME,ISOLD
        LOGICAL                           :: APPEND,MESSAGE,NEWFILE

        LOGICAL                           :: ASKNAM = .TRUE.
        LOGICAL                           :: FOUROK = .FALSE.
        LOGICAL                           :: MASKIT = .TRUE.

        character(len=maxnam)             :: filjnk = 'jnkwin'
        character(len=maxnam)             :: filpow = 'jnkpow'

        COMMON /SEARCH_RANGE/ STRT_SRCH,STOP_SRCH,NMB

        REAL, PARAMETER    :: QUADPI = 3.141592653589793238462

        INTEGER, PARAMETER   :: LUNM    = 20
        INTEGER, PARAMETER   :: LUNPOW  = 21
        INTEGER, PARAMETER   :: LUNDOCM = 81
        INTEGER, PARAMETER   :: LUNDOCN = 82
        INTEGER, PARAMETER   :: LUNDOCD = 83
        INTEGER, PARAMETER   :: LUNXM   = 0
        CHARACTER, PARAMETER :: MODE    = ' '
        INTEGER                           :: IMG1

C       DEFOCUSGUESSING3 FOR WEAK DEFOCUS ASTIGMATISM ESTIMATION  
C       DEFOCUSGUESSING2 FOR STRONG ASTIGMATISM ESTIMATION

        STRT_SRCH = 3000.0
        STOP_SRCH = 200000.0

C       OPEN FIRST MICROGRAPH FILE, NOT FOURIER
        MAXIMM  = 0                       
        NILMAX  = NIMAX         ! INUMBR FROM CMLIMIT
        IMG1 = 0    ! Passing uninitialized variable may result in 'INVALID IMAGE NUMBER' error
        CALL OPFILES(0,LUNM,LUNDOCM,LUNXM,
     &             ASKNAM,FILPAT,NLETM, 'O',         
     &             IFORMM,NXM,NYM,NZM,MAXIMM,         
     &             'MICROGRAPH IMAGE',                              
     &             FOUROK, INUMBR,NILMAX,             
     &             NDUM,NGOT1,IMG1, IRTFLG)           
        IF (IRTFLG .NE. 0) RETURN                    
        IF (NZM > 1) THEN
           CALL ERRT(101,'OPERATION DOES NOT WORK ON VOLUMES',NE)
           RETURN
        ENDIF

        NXT    = 500
        IXOVER =  50
        IYOVER = -999
        CALL RDPRI3S(NXT,IXOVER,IYOVER,NOT_USED,
     &             'TILE SIZE, X & Y TILE PERCENT OVERLAP',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        IF (IYOVER < 0) IYOVER = IXOVER

        IF (IXOVER < 0 .OR. IXOVER > 99) THEN
           CALL ERRT(102,'X OVERLAP OUT OF RANGE (0..99)',IXOVER)
           RETURN
        ELSEIF (IYOVER < 0 .OR. IYOVER > 99) THEN
           CALL ERRT(102,'Y OVERLAP OUT OF RANGE (0..99)',IYOVER)
           RETURN
        ELSEIF (NXT < 5 .OR. NXT > NXM .OR. NXT > NYM) THEN
           CALL ERRT(102,'TILE SIZE OUT OF RANGE',NXT)
           RETURN
        ENDIF

        IXB =  500
        IYB = -999
        CALL RDPRI2S(IXB,IYB,NOT_USED,'X & Y TILING BORDER',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        IF (IYB < 0) IYB = IXB

        IF (IXB < 0 .OR. IXB > NXM) THEN
           CALL ERRT(102,'X BORDER OUT OF RANGE' ,IXB)
           RETURN
        ELSEIF (IYB < 0 .OR. IYB > NYM) THEN
           CALL ERRT(102,'Y BORDER OUT OF RANGE',IYB)
           RETURN
        ENDIF

        ! NORMALIZATION OF % OVERLAP IN X & Y
        ROX  = 100.0 / (100.0 - IXOVER)                  
        ROY  = 100.0 / (100.0 - IYOVER)                  

        ! NO. OF TILES HORIZONTAL DIM.(X)
        XOFF  = FLOAT(IXOVER) * FLOAT(NXT) / 100.0
        IXOFF = NXT - INT(XOFF)
        NPX   = 0
        IXGO  = IXB + 1
        !WRITE(6,*) '  ixover:',ixover,'xoff:',xoff,'ixoff:',ixoff
        !write(6,*)  ' npx,ixgo,ixb:',npx,ixgo,ixb
        DO
           IF ((IXGO + NXT) > (NXM - IXB)) EXIT
           IXGO = IXGO + IXOFF
           NPX  = NPX + 1
        ENDDO

        YOFF  = FLOAT(IYOVER) * FLOAT(NXT) / 100.0
        IYOFF = NXT - INT(YOFF)
        NPY   = 0
        IYGO  = IYB + 1
        DO
           IF ((IYGO + NXT) > (NYM - IYB)) EXIT
           IYGO = IYGO + IYOFF
           NPY  = NPY + 1
        ENDDO
        WRITE(NOUT,'(A,I4, A,I4, A,I4)') 
     &          '  Tiles:',NPX,' x',NPY,'  Size:',NXT

        PS  = 1
        CSM = 2.0
        CALL RDPRM2S(PS,CSM,NOT_USED,
     &      'PIXEL SIZE [A] & SPHERICAL ABERRATION CS [MM]',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
           IF (CS < 0.0001)    CS = 0.0001

C       CONVERT CS TO [A]
        CS = CSM * 1.0E07

        KEV = 200.0
        !LAMBDA = 2.508E-02
        CALL RDPRM1S(KEV,NOT_USED,'ELECTON VOLTAGE [KEV]',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        LAMBDA = 12.398 / SQRT(KEV*(1022+KEV))
          
        CONTRAST = 0.100 
        CALL RDPRM1S(CONTRAST,NOT_USED,
     &           'AMPLITUDE CONTRAST RATIO [0-1]',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999



        NXLD = NXT + 2 - MOD(NXT,2)
        NMB  = NXT/2  + MOD(NXT,2)

        ALLOCATE(FNOI_OUT(NMB, NT),
     &           BUFIN   (NXM, NXT),
     &           BUFWIN  (NXT, NXT),
     &           BUFFFT  (NXLD,NXT),
     &           BUFPOW  (NXT, NXT),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0)  THEN
           MWANT = NXM*NXT + 2*NXT*NXT + NXLD*NXT + NMB*NT
           CALL ERRT(46,'TFED; BUFIN..',MWANT)
           RETURN
        ENDIF

        IYWINGO = IYB + 1         ! FIRST Y
        IGOT    = 0               ! TILE COUNTER

        DO KY = 1, NPY            ! LOOP OVER ALL TILE ROWS

C          READ IN NXT ROWS
           IYW = 1
           DO IYF = IYWINGO, IYWINGO+NXT-1
              CALL REDLIN(LUNM,BUFIN(1,IYW),NXM, IYF)
              IYW = IYW + 1
           ENDDO

           !call chkfile(filjnk(1:8),66,1, nxm,nxt,1, bufin,irtflg)

           IXW = IXB + 1          ! FIRST X
           DO KX = 1, NPX         ! LOOP OVER ALL HORIZ. TILE ROWS

              BUFWIN(1:NXT, 1:NXT) = BUFIN(IXW:IXW+NXT-1, 1:NXT)
              IXW = IXW + IXOFF

              !call inttochar(igot,filjnk(7:9),nlet,3)
              !call chkfile(filjnk(1:9),66,1,nxt,nxt,1,bufwin,irtflg)

C             CORRECT RAMP INTENSITIES
              CALL RAMP_PB(BUFWIN,NXT,NXT,.FALSE.,NOUT)

C             PAD INTO BUFFFT
              BUFFFT(1:NXT, 1:NXT) = BUFWIN(1:NXT,1:NXT)

              INV = +1   ! FORWARD FFT 
              CALL FMRS_2(BUFFFT,NXT,NXT,INV)
              CALL PW2SR (BUFFFT,NXT,NXT,MODE)

              IGOT = IGOT + 1

              IF (IGOT > 1) THEN
                 BUFPOW(1:NXT, 1:NXT) = BUFPOW(1:NXT, 1:NXT) + 
     &                                  BUFFFT(1:NXT, 1:NXT)
              ELSE
                 BUFPOW(1:NXT, 1:NXT) = BUFFFT(1:NXT, 1:NXT)
              ENDIF
              !call inttochar(igot,filjnk(7:9),nlet,3)
              !call chkfile(filjnk(1:9),66,1,nxt,nxt,1,bufwin,irtflg)

           ENDDO

           IYWINGO = IYWINGO + IYOFF
        ENDDO

        SQRTNT = SQRT(FLOAT(IGOT))     ! Why?? Legacy continued
        BUFPOW = BUFPOW / SQRTNT
        !call chkfile('jnkpow',66,1,nxt,nxt,1,bufpow,irtflg)

        IF (MASKIT) THEN  ! NOT ACTIVATED
           ICX      = NXT/2 + 1          ! X,Y CENTER
           BACKPIX  = BUFPOW(ICX+10,5)   ! PERIPHERAL PIXEL VALUE

           RADMASKI = .95 * NXT          ! POWER SPECTRUM MASK RADIUS (A)  (0 = DON'T MASK)
           RADMASKI = .5 * NXT          ! POWER SPECTRUM MASK RADIUS (A)  (0 = DON'T MASK)

           CIRCF    = 2 * PS / RADMASKI  ! MASK RADIUS (NOW IN 1/PX)
           CIRCI    = CIRCF * NXT        ! COMPUTE MASK RADIUS

C          CENTRAL MASKING
           DO  IY = 1,NXT
              FI1 =  FLOAT(IY-ICX)**2 
              DO  IX = 1,NXT
                 CRAD2 = FI1 + FLOAT(IX-ICX)**2
                 IF (CRAD2 < CIRCI) BUFPOW(IX,IY) = BACKPIX
              ENDDO
           ENDDO
        ENDIF
        write(NOUT,'(A,F5.0)') '  Central PS mask radius:',CIRCI

C       ESTIMATING DEFOCUS BASED ON OVERALL POWER SPECTRUM
C       SET THE START ANGLE AND STEP TO BE 0.0 AND 180.0.
        BETA0 = 0.0
        BETA1 = 180.0

        !write(6,*) 'CALLING DEFOCUSGUESSING1 -----------:'
        CALL DEFOCUSGUESSING1(BUFPOW,NXT,NXT,NMB,BETA0,BETA1,
     &           PS,CS,LAMBDA,CONTRAST,AV_DEFO,
     &           ICUT_LOW_FRQ,ICN_SND,NDGREE,FNOI_OUT,NT,XX_CC,NUMBER)
        IF (NUMBER <= 0) GOTO 9999

        FNOI_OUT(1:NUMBER,2:4)= EXP(FNOI_OUT(1:NUMBER,2:4))
        FNOI_OUT(1:NUMBER,1)  = FNOI_OUT(1:NUMBER,1)/2.0/PS/NMB

        FNOI_OUT(1:NUMBER,3) = FNOI_OUT(1:NUMBER,3)-FNOI_OUT(1:NUMBER,2)
        FNOI_OUT(1:NUMBER,4) = FNOI_OUT(1:NUMBER,4)-FNOI_OUT(1:NUMBER,2)

        ADDEXT  = .TRUE.
        GETNAME = .TRUE.
        ISOLD   = .FALSE.
        APPEND  = .FALSE.
        MESSAGE = .TRUE.
        IRTFLG  = -8         ! NO IC USE

        CALL OPENDOC(DOCNAM,ADDEXT,NLET,LUNDOCN,LUNDOCNO,GETNAME,
     &           'DEFOCUS NOISE DOC',ISOLD,APPEND,MESSAGE,
     &            NEWFILE,IRTFLG)

C           123456789 123456789 123456789 123456789 123456789 123456789 
        COMMENT = 
     &     '  Spa. freq.(1/A), Back. noise, Back. subtr. PS,  Env(f)**2'
        CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:69),IRTFLG)

        DO KEY=1,NUMBER
            DLIST(1) = FNOI_OUT(KEY,1)
            DLIST(2) = FNOI_OUT(KEY,2)
            DLIST(3) = FNOI_OUT(KEY,3)
            DLIST(4) = FNOI_OUT(KEY,4)
            CALL LUNDOCWRTDAT(LUNDOCNO,KEY,DLIST,4,IRTFLG)
        ENDDO
        CLOSE(LUNDOCN)

        !write(6,*) 'DEFOCUS FOUND:',AV_DEFO 

        IF ( XX_CC == 9999 .OR. AV_DEFO <= STRT_SRCH) THEN
C           ESTIMATE DEFOCUS OF DIFFICULT CASES
            XX_CC   = 9999.
            AST_AGL = 0.0
            TMP_AMP = 0.0
            TMP_SUM = 0.0
        ELSE

           WRITE(NOUT,*) ' ASTIGMATISM CALCULATION NO LONGER DONE'
c          calculation is buggy often crashes. removed mar 2013 al           

        ENDIF

C       ROUGH ESTIMATION OF DEFOCUS IN DIFFICULT CASE AND NO ASTIGMATISM ESTIMATION IN SUCH CASE!
        IF (XX_CC > 1) THEN
            CALL DEFOCUSGUESSING3(BUFPOW,NXT,NXT,NMB,
     &           PS,CS,LAMBDA,CONTRAST,AV_DEFO,XX_CC)
        ENDIF
        
        WRITE(NOUT,'(A,F8.1)')'  OVERALL DEFOCUS: ',AV_DEFO

C       OPEN OUTPUT DOC FILE (FOR APPENDING)
        ADDEXT  = .TRUE.
        GETNAME = .TRUE.
        ISOLD   = .FALSE.
        APPEND  = .TRUE.
        MESSAGE = .TRUE.
        IRTFLG  = -8         ! NO IC USE

        LUNRET = LUNDOCD
        CALL OPENDOC(DOCNAM,addext,NLET,LUNDOCD,LUNRET,GETNAME,
     &             'OUTPUT DEFOCUS DOCUMENT',ISOLD,APPEND,MESSAGE,
     &             NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (NEWFILE) THEN

C                     123456789 123456789 123456789 123456789 123456789 123456789 
           COMMENT = 'CONTENTS:   MICROGRAPH DEFOCUS VALUES' 

           CALL LUNDOCPUTCOM(LUNRET,COMMENT(1:37),IRTFLG)
C                     123456789 123456789 123456789 123456789 123456789 123456789012
           COMMENT = '          MICR.  DEFOCUS  CUTOFF'
           CALL LUNDOCPUTCOM(LUNRET,COMMENT(1:35),IRTFLG)
        ENDIF

        KEY = 1
        CALL RDPRI1S(KEY, NOT_USED,
     &     'KEY/IMAGE NUMBER FOR DOC FILE',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
      
        FORMOUT  = '(I7,1X,I2,1X,F7.0,1X,F8.1,1X,F6.3)'
        DLIST(1) = KEY 
        DLIST(2) = AV_DEFO 
        DLIST(3) = XX_CC 
        CALL LUNDOCWRTDATF(LUNRET,KEY,DLIST,3,FORMOUT,IRTFLG)

        CLOSE(LUNDOCD)

C       SAVE AVERAGE POWER SPECTRUM
        IFORM = 1
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNPOW,'U',IFORM,NXT,NXT,1,
     &               MAXIM,'AVERAGE POWER SPECTRUM',FOUROK,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        CALL WRTVOL(LUNPOW,NXT,NXT,1,1,BUFPOW,IRTFLG)

        CALL REG_SET_NSEL(1,3, AV_DEFO, XX_CC,FLOAT(IGOT),0,0,IRTFLG)

9999    IF (ALLOCATED(POW2))     DEALLOCATE(POW2)
        IF (ALLOCATED(FNOI_OUT)) DEALLOCATE(FNOI_OUT)
        IF (ALLOCATED(BUFIN))    DEALLOCATE(BUFIN)
        IF (ALLOCATED(BUFWIN))   DEALLOCATE(BUFWIN)
        IF (ALLOCATED(BUFFFT))   DEALLOCATE(BUFFFT)
        IF (ALLOCATED(BUFPWS))   DEALLOCATE(BUFPWS)
        IF (ALLOCATED(BUFPOW))   DEALLOCATE(BUFPOW)

        END







C ****************************** TFED ********************************

        SUBROUTINE  TFED

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'F90ALLOC.INC'

        PARAMETER  (NLIST=5,NT=4)
        REAL                              :: DLIST(NLIST)
        INTEGER                           :: NUMBER
        LOGICAL                           :: NEWFILE
        REAL                              :: ASTIGFLAG
        LOGICAL                           :: DOASTIG  ! for debugging
        REAL                              :: GETDEFOCUS,LAMBDA
        CHARACTER(LEN=1)                  :: NULL = CHAR(0)
        REAL, DIMENSION(:,:), ALLOCATABLE :: POW2,TMP_OUT
        REAL, DIMENSION(:),ALLOCATABLE    :: XDEFO
        CHARACTER(LEN=MAXNAM)             :: FILNAM,DOCNAM
        CHARACTER(LEN=80)                 :: COMMENT

        COMMON /SEARCH_RANGE/ STRT_SRCH,STOP_SRCH,NMB

        PARAMETER (QUADPI = 3.141592653589793238462)

        DATA      LUN1/88/
        DATA      NDOC/88/


C       DEFOCUSGUESSING3 FOR WEAK DEFOCUS ASTIGMATISM ESTIMATION  
C       DEFOCUSGUESSING2 FOR STRONG ASTIGMATISM ESTIMATION

        STRT_SRCH = 3000.0
        STOP_SRCH = 200000.0
        MAXIM1    = 0

        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYP,NSAM1,NROW1,NSLIC1,
     &             MAXIM1,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(101,',OPENING INPUT FILE',NE)
            RETURN
        ENDIF

        IF (NSAM1 .NE. NROW1) THEN
            CLOSE(LUN1)
            CALL ERRT(101,'IMAGE MUST BE SQUARE',NE)
            RETURN
        ENDIF

        CALL RDPRM2(PS,CS,NOT_USED,
     &      'PIXEL SIZE[A] & SPHERICAL ABERRATION CS [MM]')
           IF (CS < 0.0001)    CS = 0.0001

C       CONVERT CS TO [A]
        CS = CS*1.0E07

        CALL RDPRM(LAMBDA,NOT_USED,'WAVELENGTH LAMBDA [A]')

        ASTIGFLAG = 0
        CALL RDPRM2S(CONTRAST,ASTIGFLAG,NOT_USED,
     &           'AMPLITUDE CONTRAST RATIO [0-1]',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        DOASTIG = (ASTIGFLAG > 0)

        NMB  = NSAM1/2  + MOD(NSAM1,2)

        ALLOCATE(POW2(NSAM1,NROW1),TMP_OUT(NMB,NT),STAT=IRTFLG)
        IF (IRTFLG .NE. 0)  THEN
           MWANT = NSAM1*NROW1 + NMB*NT
           CALL ERRT(46,'TFED; POW2..',MWANT)
           RETURN
        ENDIF

        CALL READV(LUN1,POW2,NSAM1,NROW1,NSAM1,NROW1,NSLIC1)
        CLOSE(LUN1)

C       ESTIMATING DEFOCUS BASED ON OVERALL POWER SPECTRUM
C       SET THE START ANGLE AND STEP TO BE 0.0 AND 180.0.
        BETA0 = 0.0
        BETA1 = 180.

        CALL DEFOCUSGUESSING1(POW2,NSAM1,NROW1,NMB,BETA0,BETA1,
     &           PS,CS,LAMBDA,CONTRAST,AV_DEFO,
     &           ICUT_LOW_FRQ,ICN_SND,NDGREE,TMP_OUT,NT,XX_CC,NUMBER)

        IF (NUMBER .LE. 0) THEN   ! ERROR RETURN al
           GOTO 9999
        ENDIF

        TMP_OUT(1:NUMBER,2:4)= EXP(TMP_OUT(1:NUMBER,2:4))
        TMP_OUT(1:NUMBER,1)  = TMP_OUT(1:NUMBER,1)/2.0/PS/NMB
        TMP_OUT(1:NUMBER,3)  = TMP_OUT(1:NUMBER,3)-TMP_OUT(1:NUMBER,2)
        TMP_OUT(1:NUMBER,4)  = TMP_OUT(1:NUMBER,4)-TMP_OUT(1:NUMBER,2)

        CALL OPENDOC(DOCNAM,.TRUE.,NLET,NDOC,LUNDOC,.TRUE.,
     &           'DEFOCUS DOCUMENT',.FALSE.,.FALSE.,.TRUE.,
     &            NEWFILE,IRTFLG)

C           123456789 123456789 123456789 123456789 123456789 123456789 
        COMMENT = 
     &     '  Spa. freq.(1/A), Back. noise, Back. subtr. PS,  Env(f)**2'
        CALL LUNDOCPUTCOM(LUNDOC,COMMENT(1:69),IRTFLG)

        DO II=1,NUMBER
            DLIST(1) = II
            DLIST(2) = TMP_OUT(II,1)
            DLIST(3) = TMP_OUT(II,2)
            DLIST(4) = TMP_OUT(II,3)
            DLIST(5) = TMP_OUT(II,4)
            CALL LUNDOCWRTDAT(LUNDOC,II,DLIST(2),NLIST-1,IRTFLG)
        ENDDO
        CLOSE(NDOC)

        AST_AGL = 0.0
        AST_AMP = 0.0
        AST_DEF = AV_DEFO

        IF (DOASTIG .AND.
     &      XX_CC .NE. 9999 .AND.
     &      AV_DEFO > STRT_SRCH ) THEN
           
C          ESTIMATION OF ASTIGMATISM; STEP IS 18 DEGREES  BUGGY!!!!!!
           X_AV_DEFO = AV_DEFO
           IDGREE    = NDGREE
           CALL AST_CALC(POW2,NSAM1,NROW1,NMB,NT,PS,CS,LAMBDA,
     &           CONTRAST,X_AV_DEFO,IDGREE,AST_AGL,AST_AMP,AST_DEF)

        ENDIF

C       ROUGH ESTIMATION OF DEFOCUS IN DIFFICULT CASE AND NO ASTIGMATISM ESTIMATION IN SUCH CASE!
        IF (XX_CC > 1) THEN
            CALL DEFOCUSGUESSING3(POW2,NSAM1,NROW1,NMB,
     &           PS,CS,LAMBDA,CONTRAST,AV_DEFO,XX_CC)
        ENDIF

        CALL REG_SET_NSEL(1,5,AST_AGL,AST_AMP,AST_DEF,AV_DEFO,
     &                       XX_CC,IRTFLG)

9999    DEALLOCATE(POW2,TMP_OUT)

        END



CC======================================================================

        SUBROUTINE CTF_SIGNAL1(VALUE, X,LENGTH,
     &           PS, CS, LAMBDA, DEFOCUS, CONTRAST)

C       X IS X COORDINATES!
        IMPLICIT NONE

        REAL        QUADPI
        PARAMETER   (QUADPI = 3.1415926535897932384626)
        REAL        IPS
        COMMON /SEARCH_RANGE/ DZLOW1,DZHIGH1,NMB
        INTEGER     LENGTH, I,NMB
        REAL        VALUE(LENGTH),X(LENGTH)
        REAL        PS, CS, LAMBDA, DEFOCUS, CONTRAST, FF
        REAL        AA, BB, CC,DZLOW1,DZHIGH1

        AA = 0.5*QUADPI*CS*LAMBDA**3
        BB = QUADPI*LAMBDA
        CC = ATAN(CONTRAST/(1.0-CONTRAST))
        IPS = 1.0/(2.0*PS*NMB)

        DO I=1, LENGTH
            FF = X(I)*IPS
            VALUE(I) = SIN((AA*FF**2-BB*DEFOCUS)*FF**2-CC)
        ENDDO
        END

C=========================================================================

         SUBROUTINE GET_CTF_ZERO(NUMBER,PS, CS,LAMBDA,
     &           DEFOCUS, CONTRAST,IFRQ1,IFRQ2,NZERO,IFIRST,ISECND)

        COMMON /SEARCH_RANGE/ DZLOW1,DZHIGH1,NMB
        PARAMETER (QUADPI = 3.1415926535897932384626)
         REAL LAMBDA,DZLOW1,DZHIGH1

         AA = -0.5*QUADPI*CS*LAMBDA**3
         BB = QUADPI*LAMBDA
         BB=BB*DEFOCUS
         CC = -ATAN(CONTRAST/(1.0-CONTRAST))
        FRQ1=REAL(IFRQ1)/2.0/PS/REAL(NMB)
        FRQ2=REAL(IFRQ2)/2.0/PS/REAL(NMB)
         TMPN1=AA*FRQ1**4+BB*FRQ1**2+CC
         TMPN1=TMPN1/QUADPI
         N1=INT(TMPN1)
        TMPN2=AA*FRQ2**4+BB*FRQ2**2+CC
         TMPN2=TMPN2/QUADPI
         N2=INT(TMPN2)
        N=N2-N1+1
        TMP1=4.*AA

        DO II=1,N2
            TMP2=CC-REAL(II-1)*QUADPI
            TMP2=SQRT(BB**2-TMP1*TMP2)
            TMP2=TMP2-BB
            TMP2=SQRT(TMP2/2.0/AA)
            NFRQ=INT(TMP2*2.0*PS*NMB)
        ENDDO

        IF(N1.GT.1) THEN
            DO II=1,N1
          TMP2=CC-REAL(II-1)*QUADPI
          TMP2=SQRT(BB**2-TMP1*TMP2)
          TMP2=TMP2-BB
          TMP2=SQRT(TMP2/2.0/AA)
          NFRQ=INT(TMP2*2.0*PS*NMB)

          IF(II.EQ.1) THEN
              IFIRST=NFRQ
          ELSEIF(II.EQ.2) THEN
              ISECND=NFRQ
          ENDIF
            ENDDO
        ELSE
CC==FURTHER SEARCH CTF ZERO, FOR THOSE N.LE.0
            IFIRST=IFRQ1
            ISECND=0
        ENDIF

        NZERO=N
3000    RETURN
        END

CC===================================================================

        SUBROUTINE ZH_CRCSE2(BUF,SEC,NSAM,NROW,IR,BETA0,BETA1)

CC==    ROTATIONAL AVERAGE FROM BETA0 TO BETA1!

        DIMENSION BUF(NSAM,NROW), SEC(IR), SNO(IR)

        PARAMETER (QUADPI= 3.1415926535897932384626)
        PARAMETER (DGR_TO_RAD = (QUADPI/180.))

CC==    BEGIN
        XBETA0 = BETA0*DGR_TO_RAD
        XBETA1 = BETA1*DGR_TO_RAD+XBETA0
CC==       
         SEC = 0.0
         SNO = 0.0        
         DO J=1,NROW
            KJ = J-NROW/2-1

            IF (IABS(KJ) .LE. IR-1)  THEN

               DO I=1,NSAM
                  KI = I-NSAM/2-1
                  R  = SQRT(FLOAT(KJ*KJ)+FLOAT(KI*KI))+1.0
                  L  = R

                  IF (L.LE.IR-1) THEN
                     ANGLE = ATAN2(REAL(KI),REAL(KJ))
                     IF (ANGLE.LT.0.0) ANGLE=ANGLE+2.*QUADPI
                     IF(ANGLE.GE.XBETA0.AND.ANGLE.LT.XBETA1) THEN
                           XD       = R-L
                           SEC(L)   = SEC(L)  + BUF(I,J)*(1.0-XD)
                           SEC(L+1) = SEC(L+1)+ BUF(I,J)*XD
                           SNO(L)   = SNO(L)  + 1.0-XD
                           SNO(L+1) = SNO(L+1)+ XD
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

         DO I=1,IR
             SEC(I)= SEC(I) / AMAX1(1.0,SNO(I))
         ENDDO

         END

CC=============================================================

        SUBROUTINE PARTI8 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT1,N,ICUT1,ICUT2,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ
CC==LOW FREQUENCIES REGION RE-FITTING!!!
CC==INTERMDEDIATE REGION FITTING
        ISTART=ICUT1
        ISTOP=ICUT2
        N_3=N
CC==
        KS=ISTOP-ISTART+1
        L=2
        ISWI=4
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
        TT1_Q2(KS+1,1:N-1)=Q2(ISTART,1:N-1)
        TT1_Q2(KS+2,1:N-1)=Q2(ISTOP,1:N-1)
        TT1_Q2(KS+2,N)=1.0
        TT1_Q2(KS+1,N)=1.0
        TT1_Q2(KS+2,N+1)=PLOT1(ISTOP,3)
        TT1_Q2(KS+1,N+1)=PLOT1(ISTART,3)

        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-3)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
CC==
        KS=ISTART
        L=1
        ISWI=5
        TT1_Q2(1:KS,1:N+1)=Q2(1:ISTART,1:N+1)
        TT1_Q2(KS+1,1:N-1)=Q2(ISTART,1:N-1)
        TT1_Q2(KS+1,N)=1.0
        TT1_Q2(KS+1,N+1)=PLOT1(ISTART,3)

        CALL LSFIT(KS,L,N2,PLOT1(1,2+(ISWI-4)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
        END

CC========================================================

        SUBROUTINE XFIT_ZH(K,N2,N,PLOT,ICUT1,ICUT2)

CCC==INEQUALITY CONSTRAINED LINEAR SQUARE MINIMIZATION
CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
CC==ENVELOPE FITTING 

        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA

        L=2
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF
CC==
        DO I=1,N-1
            Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
        Q2(K+1:2*K,1:N+1)=Q2(1:K,1:N+1)

        CALL PARTI8 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICUT1,ICUT2,KP)
        DEALLOCATE(Q2)

        END

CCC=================================

        SUBROUTINE PARTI9(KLMD,NKLMD,KLM2D,N2D,Q2,
     &           K,N2,N,PLOT1,IFIT1,IFIT2,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ
CC==LOW FREQUENCIES REGION RE-FITTING!!!
CC==THIS SUBROUTINE IS FOR THE PURPOSE OF SEARCHING POINT TO CUT LOW FREQUENCY REGION OFF
        ISTART=IFIT1
        ISTOP=IFIT2
        N_3=N
        KS=ISTOP-ISTART+1
        L=0
CC==ENV
        ISWI=1
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)

        CALL LSFIT(KS,L,N2,PLOT1(ISTART,(ISWI-1)*2+2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
CC==BACKGROUND
        ISWI=2
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
        CALL LSFIT(KS,L,N2,PLOT1(ISTART,(ISWI-1)*2+2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
        END

CC===============================================================

        SUBROUTINE XFITLINE1(K,N2,PLOT,IFIT1,IFIT2,ICUT)

CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)

CC==IN THIS CASE, K=NMB
        X_THR=0.000001
        N=2
        L=0
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF
CC==
        Q2(1:K,1)=PLOT(1:K,1)
        Q2(1:K,2)=1.0
        Q2(1:K,3)=PLOT(1:K,3)
CC==
        CALL PARTI9(KLMD,NKLMD,KLM2D,N2D,Q2,
     &      K,N2,N,PLOT,IFIT1,IFIT2,KP)
CC==SEARCHING THE CUTTING POINT IN THE LOW FREQUENCIES
        PLOT(IFIT1:IFIT2,4)=PLOT(IFIT1:IFIT2,3)-PLOT(IFIT1:IFIT2,2)
        PLOT(IFIT1:IFIT2,4)=PLOT(IFIT1:IFIT2,4)/PLOT(IFIT1:IFIT2,2)

        DO I=IFIT1,IFIT2
            XX=ABS(PLOT(I,4))
            IF(XX.LE.X_THR) GOTO 1100
        ENDDO
1100    ICUT=I

        DEALLOCATE(Q2)
        END

CC==============================================================================
        SUBROUTINE XFITLINE2(K,N2,PLOT,IFIT1,IFIT2,ICUT)

CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
        REAL, DIMENSION(:), ALLOCATABLE :: PLOT1
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
 
CC==IN THIS CASE, K=NMB
        X_THR=0.000001
        N=2
        L=0
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(PLOT1(K), Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED, Q2..',KLM2D*N2D)
            RETURN
        ENDIF

CC==
        Q2(1:K,1)=PLOT(1:K,1)
        Q2(1:K,2)=1.0
        Q2(1:K,3)=PLOT(1:K,3)
CC==
        CALL PARTI9(KLMD,NKLMD,KLM2D,N2D,Q2,
     &               K,N2,N,PLOT,IFIT1,IFIT2,KP)

CC==SEARCHING THE CUTTING POINT IN THE LOW FREQUENCIES
        PLOT1(IFIT1:IFIT2)=PLOT(IFIT1:IFIT2,4)-PLOT(IFIT1:IFIT2,3)

        DO II=IFIT1,IFIT2
            IF(PLOT(II,3).NE.0.0) THEN
          PLOT1(II)=PLOT1(II)/PLOT(II,4)
          XX=ABS(PLOT1(II))

          IF(XX.LE.X_THR) GOTO 1100
            ENDIF
        ENDDO
1100    ICUT=II
        DEALLOCATE(Q2,PLOT1)
        END

CC================================================================

        SUBROUTINE PARTI11 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT1,N,ICUT1,ICUT2,ICNSTRNT,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ

CC==LOW FREQUENCIES REGION RE-FITTING!!!
        ISTART=ICUT1
        ISTOP=ICUT2
        N_3=N
        L=2
        KS=ISTOP-ISTART+1
CC==EVN with 2 EQ CONSTRAINTS
        ISWI=3
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
CC== EQUALITY CONSTRAINT CONDITION IN TT1_Q2
        TT1_Q2(KS+1,1:N-1)=Q2(ISTART,1:N-1)
        TT1_Q2(KS+2,1:N-1)=Q2(ICNSTRNT,1:N-1)
        TT1_Q2(KS+2,N)=1.0
        TT1_Q2(KS+1,N)=1.0
        TT1_Q2(KS+1,N+1)=PLOT1(ISTART,3)
        TT1_Q2(KS+2,N+1)=PLOT1(ICNSTRNT,3)

        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-3)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)

        END

CC=======================================================================

        SUBROUTINE PARTI10 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT1,N,ICUT1,ICUT2,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ

CC==LOW FREQUENCIES REGION RE-FITTING!!!
        ISTART=ICUT1
        ISTOP=ICUT2
        N_3=N
        L=1
        KS=ISTOP-ISTART+1
CC==ENV WITH TWO CONSTRAINTS
        ISWI=6
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
CC== EQUALITY CONSTRAINT CONDITION IN TT1_Q2
        TT1_Q2(KS+1,1:N-1)=Q2(ISTART,1:N-1)
        TT1_Q2(KS+1,N)=1.0
        TT1_Q2(KS+1,N+1)=PLOT1(ISTART,3)

        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-6)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)

        END

CC======================================================================

        SUBROUTINE XFIT_BKGND1(K,N2,N,PLOT,ICUT1,ICUT2,
     &                   ICNSTRNT,ICN_SND,
     &                   PS,CS,LAMBDA,CONTRAST,DEFOCUS,IILOOP,NZERO)

CCC==INEQUALITY CONSTRAINED LINEAR SQUARE MINIMIZATION
CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA

        L=1
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF
CC==
        DO I=1,N-1
           Q2(ICUT1:ICUT2,I)=PLOT(ICUT1:ICUT2,1)**I
        ENDDO
        Q2(ICUT1:ICUT2,N)=1.0
        Q2(ICUT1:ICUT2,N+1)=PLOT(ICUT1:ICUT2,3)

        CALL PARTI15(KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICUT1,ICUT2,ICNSTRNT,KP)

        PLOT(ICUT1:ICUT2,3)=PLOT(ICUT1:ICUT2,3)-PLOT(ICUT1:ICUT2,2)
        PLOT(ICUT1:ICUT2,4)=PLOT(ICUT1:ICUT2,4)-PLOT(ICUT1:ICUT2,2)

        ISWI=1
        DZMAX=GETDEFOCUS(PLOT(1,3),PLOT(1,1),K,PS,CS,LAMBDA,
     &           CONTRAST,XSCORE,PLOT(1,4),ICUT1,ICUT2,ISWI)

         DEFOCUS=DZMAX

        CALL GET_CTF_ZERO(K,PS, CS,LAMBDA,
     &           DEFOCUS, CONTRAST,ICUT1,ICUT2,NZERO,IFIRST,ISECND)

        CALL CTF_SIGNAL1(PLOT(1,2), PLOT(1,1),K, PS, CS,
     &                LAMBDA, DEFOCUS, CONTRAST)

        PLOT(ICUT1:ICUT2,2)=PLOT(ICUT1:ICUT2,2)**2*PLOT(ICUT1:ICUT2,4)

        DEALLOCATE(Q2)
        END

CC======================================================================

        SUBROUTINE DEFOCUSGUESSING2(POW2,
     &           NSAM1,NROW1,NMB,XSTART,XSTEP,PS,CS,LAMBDA,CONTRAST,
     &           DEFOCUS,IILOOP,ICUT_LOW_FRQ,ICN_SND,NDGREE)

        DIMENSION POW2(NSAM1,NROW1)
        REAL, DIMENSION(:,:), ALLOCATABLE :: PLOT,TMP
        REAL LAMBDA

CC==ASTIGMATISM ESTIMATION!!!
        ALLOCATE(TMP(NMB,2),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
        ENDIF
        CALL ZH_CRCSE2(POW2,TMP(1,2),NSAM1,NROW1,NMB,XSTART,XSTEP)
        IPW_ZERO=0

        DO I=1,NMB
           TMP(I,1)=I
        ENDDO
CC==
        IF(NMB.LT.5) THEN
            CALL ERRT(36,'TFED',NE)
            RETURN
        ENDIF

        DO I=1,NMB
            IF(TMP(I,2).EQ.0.0) THEN
            IPW_ZERO=I
            ENDIF
        ENDDO

CC==WE FIT A STRAIGHT LINE BELOW PW TO LOCATE POINT CUT OFF LOW FREQUENCY REGION
        N2=4
        NUMBER=NMB-IPW_ZERO

        ALLOCATE(PLOT(NUMBER,N2),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
         ENDIF
CC==USING GIVEN POLYNOMIAL TO FIT BACKGROUND AND OVERALL ENVELOPE! 
1100    III=NDGREE
        PLOT(:,1)=TMP(IPW_ZERO+1:NMB,1)-1.0
        PLOT(:,3)=ALOG(TMP(IPW_ZERO+1:NMB,2))

CCC==LOCATE THE FIRST ZERO OF CTF AS CUT-OFF POINT: ICUTB
        IFIT11=ICUT_LOW_FRQ-IPW_ZERO
        IF(IFIT11.LE.0) IFIT11=1

        DO II=1,999
            IFIT11=IFIT11+II-1
            IFIT2=NUMBER
            CALL XFITLINE1(NUMBER,N2,PLOT,IFIT11,IFIT2,ICUTB_1)

            IF(ICUTB_1.GT.IFIT11) GOTO 1200

        ENDDO

1200    IFIT1=IFIT11
        IFIT2=NUMBER

CC==LOCATE THE HIGHEST PEAK OF CTF AS FIRST CONTROL POINT(EQUALITY CONSTRAINT POINT)

        CALL XFITLINE2(NUMBER,N2,PLOT,IFIT1,IFIT2,ICUTA)
CC=LOCATE THE SECOND CTF ZERO
        IFIT1=ICUTA
        IFIT2=NUMBER
        CALL XFITLINE1(NUMBER,N2,PLOT,IFIT1,IFIT2,ICUTC)

        IFIT1=ICUTC
        IFIT2=NUMBER
        CALL XFITLINE2(NUMBER,N2,PLOT,IFIT1,IFIT2,ICUTD)

        IFIT1=ICUTD
        IFIT2=NUMBER
        CALL XFITLINE1(NUMBER,N2,PLOT,IFIT1,IFIT2,ICUTE)

        N=III
        ICUT1=ICUTA
        ICUT2=NUMBER
        CALL XFIT_ZH_AST(NUMBER,N2,PLOT,ICUT1,ICUT2)

        ICUT1=IFIT11
        N=III
        ICNSTRNT=ICUTC
        IBACK=1
CC==FIT OVERALL EVELOPE FROM ICUTA TO ICUT2
        CALL XFIT_LOW4(NUMBER,N2,PLOT,IFIT11,ICUTA)

        CALL XFIT_LOW5(NUMBER,N2,N,PLOT,ICUTA,ICUT2)

CC==CC==FIT BACKGROUND FROM ICUTB TO ICUT2
        CALL XFIT_LOW2(NUMBER,N2,PLOT,IFIT11,ICUT2,N)

CC==CC==FIT BACKGROUND FROM IBACK TO ICUTB USING SPECIFIC POLYNOMIAL DEGREE 4
        CALL XFIT_LOW10(NUMBER,N2,PLOT,IBACK,IFIT11)

        ICUT1=IFIT11
        CALL XFIT_BKGND1(NUMBER,N2,N,PLOT,ICUT1,ICUT2,ICNSTRNT,ICN_SND,
     &           PS,CS,LAMBDA,CONTRAST,DEFOCUS,IILOOP,NZERO)
        DEALLOCATE(PLOT,TMP)
        END
CC==================================================================================

        SUBROUTINE XFIT_BKGND2(K,N2,N,PLOT,ICUT1,ICUT2,
     &           PS,CS,LAMBDA,CONTRAST,DEFOCUS,IILOOP)

CCC==INEQUALITY CONSTRAINED LINEAR SQUARE MINIMIZATION
CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA

        L=0
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF

CC==
        DO I=1,N-1
            Q2(ICUT1:ICUT2,I)=PLOT(ICUT1:ICUT2,1)**I
        ENDDO

        Q2(ICUT1:ICUT2,N)=1.0
        Q2(ICUT1:ICUT2,N+1)=PLOT(ICUT1:ICUT2,3)

        CALL PARTI10(KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICUT1,ICUT2,KP)

CC==OUTPUT
        PLOT(ICUT1:ICUT2,2:4)=EXP(PLOT(ICUT1:ICUT2,2:4))
        PLOT(ICUT1:ICUT2,3)=PLOT(ICUT1:ICUT2,3)-PLOT(ICUT1:ICUT2,2)
        PLOT(ICUT1:ICUT2,4)=PLOT(ICUT1:ICUT2,4)-PLOT(ICUT1:ICUT2,2)

        ISWI=1

        DZMAX=GETDEFOCUS(PLOT(1,3),PLOT(1,1),K,PS,CS,LAMBDA,
     &           CONTRAST,XSCORE,PLOT(1,4),ICUT1,ICUT2,ISWI)

         DEFOCUS=DZMAX

        CALL CTF_SIGNAL1(PLOT(1,2), PLOT(1,1),K, PS, CS,
     F           LAMBDA, DEFOCUS, CONTRAST)

        PLOT(1:K,2)=PLOT(1:K,2)**2*PLOT(1:K,4)
        DEALLOCATE(Q2)
        END

CC=============================================================

        SUBROUTINE XFIT_LOW1(K,N2,PLOT,ICUT1,ICUT2,ICUTD,N)

CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA
CC==BEGIN: FIRST ALLOCATE ARRAYS FOR CALCULATION
        L=3
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)

        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF

CC==
        DO I=1,N-1
            Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
CC==
        Q2(K+1:2*K,1:N+1)=Q2(1:K,1:N+1)

      !write(6,*) 'in xfit_low1, calling parti13, icutd:',icutd
        CALL PARTI13 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICUT1,ICUT2,ICUTD,KP)
      !write(6,*) 'in xfit_low1, after parti13, icutd:',icutd

        DEALLOCATE(Q2)
        END

CC===========================================================
C       alters plot1

        SUBROUTINE PARTI13 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT1,N,ICUT1,ICUT2,ICUTD,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ

CC==LOW FREQUENCIES REGION RE-FITTING!!!
CC==INTERMDEDIATE REGION FITTING
        ISTART=ICUT1
        ISTOP=ICUT2
CC====ZHONG HUANG, APR,16,03
        N_3=N
        KS=ISTOP-ISTART+1
        L=3
        ISWI=7
CC=BACKGROUND FITTING WITH THREE EQ CONSTRAINTS
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
CC== EQUALITY CONSTRAINT CONDITION IN TT1_Q2

        TT1_Q2(KS+1,1:N-1)=Q2(ISTART,1:N-1)
        TT1_Q2(KS+2,1:N-1)=Q2(ISTOP,1:N-1)
        TT1_Q2(KS+3,1:N-1)=Q2(ICUTD,1:N-1)
        TT1_Q2(KS+1,N)=1.0
        TT1_Q2(KS+2,N)=1.0
        TT1_Q2(KS+3,N)=1.0
        TT1_Q2(KS+1,N+1)=PLOT1(ISTART,3)
        TT1_Q2(KS+2,N+1)=PLOT1(ISTOP,3)
        TT1_Q2(KS+3,N+1)=PLOT1(ICUTD,3)

        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-6)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)

        END

CC===============================================================

        SUBROUTINE XFIT_LOW2(K,N2,PLOT,ICUT1,ICUT2,N)

CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA

CC==BEGIN: FIRST ALLOCATE ARRAYS FOR CALCULATION
        L=2
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF

CC==
        DO I=1,N-1
            Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)

        CALL PARTI12 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &                K,N2,PLOT,N,ICUT1,ICUT2,KP)
        DEALLOCATE(Q2)
        END

CC===========================================================

        SUBROUTINE PARTI12 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT1,N,ICUT1,ICUT2,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ

CC==LOW FREQUENCIES REGION RE-FITTING!!!
CC==INTERMDEDIATE REGION FITTING
        ISTART=ICUT1
        ISTOP=ICUT2
CC==
        N_3=N
        KS=ISTOP-ISTART+1
        L=1
        ISWI=6
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
CC== EQUALITY CONSTRAINT CONDITION IN TT1_Q2
        TT1_Q2(KS+1,1:N-1)=Q2(ISTART,1:N-1)
        TT1_Q2(KS+1,N)=1.0
        TT1_Q2(KS+1,N+1)=PLOT1(ISTART,3)

        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-6)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
        END

CC=================================================================

        SUBROUTINE XFIT_LOW3(K,N2,PLOT,ICUT1,ICUT2)

CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA

CC==BEGIN: FIRST ALLOCATE ARRAYS FOR CALCULATION
        N=6
        L=2
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
        RETURN
        ENDIF
CC==
        DO I=1,N-1
            Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
CC==
        CALL PARTI14 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,ICUT1,ICUT2,KP)
        DEALLOCATE(Q2)
        END

CC===========================================================

        SUBROUTINE PARTI14 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT1,ICUT1,ICUT2,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ

CC==LOW FREQUENCIES REGION RE-FITTING!!!
CC==INTERMDEDIATE REGION FITTING
        N=4
        ISTART=ICUT1
        ISTOP=ICUT2
CC==
        N_3=N
        KS=ISTOP-ISTART+1
        L=1
        ISWI=6
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
CC== ONE EQUALITY CONSTRAINT CONDITION IN TT1_Q2
        TT1_Q2(KS+1,1:N-1)=Q2(ISTART,1:N-1)
        TT1_Q2(KS+1,N)=1.0
        TT1_Q2(KS+1,N+1)=PLOT1(ISTART,3)

        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-6)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
        END

CC===============================================================================

        SUBROUTINE DEFOCUSGUESSING1(POW2,
     &           NSAM1,NROW1,NMB,XSTART,XSTEP,PS,CS,
     &           LAMBDA,CONTRAST,AV_DEFO,ICUT_LOW_FRQ,
     &           ICN_SND,NDGREE,TMP_OUT,NT,XX_CC,NUMBER)

        INCLUDE 'CMBLOCK.INC'

        PARAMETER (NMRANK=6)
        DIMENSION POW2(NSAM1,NROW1)
        DIMENSION TMP_OUT(NMB,NT)
        DIMENSION X_DEFO(NMRANK),N_SELECT(NMRANK)

        dimension prev_x_defo(nmrank),prev_n_select(nmrank) !feb 2013 al

        INTEGER   TMP_FIRST(NMRANK),TMP_SECND(NMRANK)
        REAL, DIMENSION(:,:), ALLOCATABLE :: PLOT,TMP
        REAL                              :: LAMBDA,QUADPI

        PARAMETER (QUADPI = 3.1415926535897932384626)

        integer, save :: idone = 0

        NUMBER = 0  ! INITIALIZE FOR DEFAULT RETURN AL

        ALLOCATE(TMP(NMB,2),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED, TMP',NMB*2)
            RETURN
         ENDIF

CC==    ICUT_LOW_FRQ
CC==    N_ZERO IS MAXIMUM NZEROS OF CTF CONTROL POINTS FOR ASTIGMATISM ESTIMATION

C       rotational averging, returns tmp(1,2)  1...nmb
        CALL ZH_CRCSE2(POW2,TMP(1,2),NSAM1,NROW1,NMB,XSTART,XSTEP)

        !call chkreal('tmp(1,2)', tmp(1,2),nmb, nmb,1, 1)
        !call chkmaxloc('rotp', tmp(1,2),nmb)

        IPW_ZERO = 0
        DO I=1,NMB
           TMP(I,1) = I
        ENDDO
CC==
        IF (NMB .LT. 5) THEN
            CALL ERRT(101,'TFED, NMB MUST BE > 4',NE)
            RETURN
        ENDIF

        DO I=1,NMB
            IF (TMP(I,2) .EQ. 0.0) IPW_ZERO = I ! Finds last zero
        ENDDO

        !call chkreal('tmp', tmp(1,2),nmb, nmb,1, 1)
        !call chkminloc('tmp', tmp(1,2),nmb)
        !call chkmaxloc('tmp', tmp(1,2),nmb)

CC==    ESTIMATE MINIMUM DEFOCUS BASED ON GIVEN VOLTAGE, CS, AND PIXEL SIZE
        N2        = 4
        TMP1_DEFO = (2.*PS)**2/LAMBDA
        TMP2_DEFO = CONTRAST/(1.0-CONTRAST)
        TMP2_DEFO = ATAN(TMP2_DEFO)*TMP1_DEFO/QUADPI
        TMP3_DEFO = .5*CS*LAMBDA**2/(2.*PS)**2
        DEFO_MIN  = TMP2_DEFO+TMP3_DEFO

        TMP1_DEFO = .5*REAL(NMB-1)/PS/REAL(NMB)
        TMP1_DEFO = TMP1_DEFO*TMP1_DEFO
        TMP2_DEFO = .5/PS
        TMP2_DEFO = TMP2_DEFO*TMP2_DEFO
        DEFO_MAX  = 1.0/(TMP2_DEFO-TMP1_DEFO)/LAMBDA
        DEFO_MAX  = DEFO_MAX+.5*CS*LAMBDA**2*(TMP2_DEFO+TMP1_DEFO)

        !write(6,*) ' guessing1; def range:',DEFO_MIN,defo_max

        NUMBER    = NMB - IPW_ZERO
        IF (NUMBER .LE. 0) THEN
           CALL ERRT(101,
     &       'ABNORMAL POWER SPECTRUM, UNABLE TO DETERMINE DEFOCUS',NE)
           RETURN
        ENDIF

        TMP_XX    = NUMBER/3.
        IPT_LOW   = INT(TMP_XX)
        TMP_XX    = NUMBER*2./3.
        IPT_HIGH  = INT(TMP_XX)
        TMP_XX    = NUMBER/2.
        IPT_HALF  = INT(TMP_XX)

        NDGREE    = 0

        ALLOCATE(PLOT(NUMBER,N2),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'TFED, PLOT',NUMBER*N2)
           NUMBER = 0  ! INITIALZE FOR ERROR RETURN al
           RETURN
        ENDIF

        DO III=2,5
            PLOT(:,1) = TMP(IPW_ZERO+1:NMB,1)-1.0
            PLOT(:,3) = ALOG(TMP(IPW_ZERO+1:NMB,2))
            IFIT1     = 1
            IFIT2     = NUMBER

C==         LOCATE THE HIGHEST PEAK IN THE POWER SPECTRUM
            CALL XFITLINE2(NUMBER,N2,PLOT,IFIT1,IFIT2,ICUTL)

            IFIT11 = ICUTL
            IFIT2  = NUMBER

CC==        LOCATE THE FIRST ZERO OF CTF AS CUT-OFF POINT: ICUTB
            CALL XFITLINE1(NUMBER,N2,PLOT,IFIT11,IFIT2,ICUTB)

           ! write(6,'(a,i3,2i6,f10.1)')'gu1; cutlb:',iii,icutl,icutb
           !write(6,*)'GUESSING1 high peak:',iii,icutl,IFIT11,IFIT2,ICUTB

CC==        TO POWER SPECTRUM HAVING WEAK CTF EFFECT IN MEDIUM AND HIGH FREQUENCY
CC==        IPT_HALF IS TO AVOID INCLUDING TOO MUCH USELESS INFORMATION WHEN PIXEL SIZE IS SAMLL OR
CC==        THE MICRPGRAPHS HAVE TOO POOR RESOLUTION

            IF (ICUTB .GT. IPT_HIGH) THEN
               IFIT2 = IPT_HIGH
               CALL XFITLINE1(NUMBER,N2,PLOT,IFIT11,IFIT2,ICUTB)
            ENDIF

            IF (ICUTB.GT.IPT_LOW.AND.ICUTB.LT.IPT_HIGH) THEN
               IFIT2 = IPT_HALF
               CALL XFITLINE1(NUMBER,N2,PLOT,IFIT11,IFIT2,ICUTB)
            ENDIF

CC==        LOCATE HIGHEST PEAK OF CTF AS 1'ST CONTROL POINT(EQUALITY CONSTRAINT POINT)
            IFIT1 = ICUTB
            IFIT2 = NUMBER
            CALL XFITLINE2(NUMBER,N2,PLOT,IFIT1,IFIT2,ICUTA)

CC=         LOCATE THE SECOND CTF ZERO.
            IFIT1 = ICUTA
            IFIT2 = NUMBER
            CALL XFITLINE1(NUMBER,N2,PLOT,IFIT1,IFIT2,ICUTC)

            IFIT1 = ICUTC
            IFIT2 = NUMBER
            CALL XFITLINE2(NUMBER,N2,PLOT,IFIT1,IFIT2,ICUTD)

            IFIT1 = ICUTD
            IFIT2 = NUMBER
            CALL XFITLINE1(NUMBER,N2,PLOT,IFIT1,IFIT2,ICUTE)

            ICUT1 = ICUTA
            ICUT2 = NUMBER
            N     = III
            CALL XFIT_ZH(NUMBER,N2,N,PLOT,ICUT1,ICUT2)

            !write(6,*) 'in GUESSING1 2nd peak:',iii,ICUTa,icutc,icutd,icute

CC==        FIT BACKGROUND FROM ICUTB TO ICUT2
            ICUT1    = ICUTB
            ICNSTRNT = ICUTB
            IBACK    = 1
            CALL XFIT_LOW2(NUMBER,N2,PLOT,ICUTB,ICUT2,N)

CC==CC==    FIT BACKGROUND FROM IBACK TO ICUTB USING SPECIFIC POLYNOMIAL DEGREE 4
            IPT_START = 1
            IPT_CUT   = ICUTB
            CALL XFIT_LOW10(NUMBER,N2,PLOT,IPT_START,IPT_CUT)

            CALL XFIT_BKGND3(NUMBER,N2,N,PLOT,ICUT1,ICUT2,ICNSTRNT)

            CALL DEFOCUS_EST(NUMBER,N2,N,PLOT,ICUT1,ICUT2,
     &           PS,CS,LAMBDA,CONTRAST,DEFOX,NZERO,
     &           IFIRST,ISECND)

           !write(6,'(a,i3,2i6,f10.1)')'gu1; cuts:',iii,icut1,icut2,defox

            N_SELECT (III-1) = NZERO
            TMP_FIRST(III-1) = IFIRST
            TMP_SECND(III-1) = ISECND
            X_DEFO   (III-1) = DEFOX
        ENDDO

       !write(6,*) 'in GUESSING1 X_DEFO:',X_DEFO
       !write(6,*) 'in GUESSING1 DEFO_MIN,DEFO_MAX:',DEFO_MIN,DEFO_MAX

CC=====GUESSING DEFOCUS FROM THE FOUR TIMES CALCULATION!
C       if (x_defo(ii) .gt. defo_max  processing may
c       result in random segmentation
c       faults later in the run.  I have altered
c       it by guessing at what the author may have meant
c       It no longer crashes but am uncertain it is 
c       correct.  Feb 2013 al

        prev_x_defo   = x_defo  ! added feb 2013 al
        prev_n_select = x_defo  ! added feb 2013 al

        !write(6,'(a,5f10.1)') 'x_defo: ', x_defo(1:4)

C       can not go below min
        DO II=1,4
           IF(X_DEFO  (II) .LT. DEFO_MIN) THEN
              X_DEFO  (II) = DEFO_MIN
              N_SELECT(II) = 1
           ENDIF

C         if (x_defo(ii) .gt. defo_max  processing often 
c         result in random segmentation faults later in the run.  
C         I have altered
c         it by guessing at what the author may have meant
c         It no longer crashes but am uncertain it is 
c         correct.  Feb 2013 al
c         IF (X_DEFO(II) .GT. DEFO_MAX) THEN
c            IF (II .EQ. 1) THEN
c               X_DEFO  (II) = DEFO_MIN 
c               N_SELECT(II) = 1
c            ELSE
c               X_DEFO  (II) = X_DEFO  (II-1)  ! removed feb 2013 al
c               N_SELECT(II) = N_SELECT(II-1)  ! removed feb 2013 al
c            ENDIF
c         ENDIF
 
          IF (X_DEFO(II) .GT. DEFO_MAX) THEN
             !IF (II .EQ. 1) THEN
                X_DEFO  (II) = DEFO_MIN 
                N_SELECT(II) = 1
             !ELSE
             !   X_DEFO  (II) = prev_X_DEFO  (II-1)  
             !   N_SELECT(II) = prev_N_SELECT(II-1)  
             !ENDIF
          ENDIF

           !write(6,*) 'in GUESSING1 X_DEFO(II):',ii,X_DEFO(II)
        ENDDO

        !write(6,'(a,4f10.1)') 'x_defo1:', x_defo(1:4)

CCC=    SELECT POLYNOMIAL DEGREE ACCORDING TO CTF NZEROS!
        IF (X_DEFO(1) .GT. X_DEFO(2)) THEN
            TMP_FIRST(2) = TMP_FIRST(1)
            TMP_SECND(2) = TMP_SECND(1)
            N_SELECT(2)  = N_SELECT(1)
            X_DEFO(2)    = X_DEFO(1)
        ENDIF

        IF (X_DEFO(2) .GT. X_DEFO(3)) THEN
            TMP_FIRST(3)= TMP_FIRST(2)
            TMP_SECND(3)= TMP_SECND(2)
            N_SELECT(3) = N_SELECT(2)
            X_DEFO(3)   = X_DEFO(2)
        ENDIF

        IF (X_DEFO(3) .GT. X_DEFO(4)) THEN
            AV_DEFO      = X_DEFO(3)
            ICUT_LOW_FRQ = TMP_FIRST(3)+IPW_ZERO
            ICN_SND      = TMP_SECND(3)+IPW_ZERO
            N_TMP        = N_SELECT(3)
            NDGREE       = N_SELECT(3)-1
        ELSE
            AV_DEFO      = X_DEFO(4)
            ICUT_LOW_FRQ = TMP_FIRST(4)+IPW_ZERO
            ICN_SND      = TMP_SECND(4)+IPW_ZERO
            N_TMP        = N_SELECT(4)
            NDGREE       = N_SELECT(4)-1
        ENDIF

CC==    JUDGING WHETHER IT IS TRUE NEAR DEFOCUS:
        X_FST_ITVL = REAL(ICUTA-ICUTB)/REAL(NUMBER)

        IF (N_TMP.LE.2 .AND. X_FST_ITVL.GT..2) THEN

            !write(6,*) 'in GUESSING1 calling ctf_2_zero,icutd:',icutd
            !write(6,*) CS,        LAMBDA, CONTRAST,  PS, NDGREE, AV_DEFO 
c                 2.26E+7   1.967E-2   0.1    1.584   0     454.1

            CALL CTF_2_ZERO(TMP,NMB,NT,NUMBER,N2,TMP_OUT,IPW_ZERO,
     &           CS,LAMBDA,CONTRAST,PS,NDGREE,AV_DEFO
     &           ,ICUTA,ICUTB,ICUTD,DEFO_MIN,DEFO_MAX,XX_CC )
            RETURN

        ELSEIF (N_TMP.LE.2 .AND. X_FST_ITVL.LT..2) THEN

            !write(6,*) 'in GUESSING1 calling ctf_CORRECT:',N_TMP,X_FST_ITVL
            CALL CTF_CORRECT(TMP,NMB,NT,NUMBER,N2,TMP_OUT,IPW_ZERO,
     &          CS,LAMBDA,CONTRAST,PS,AV_DEFO
     &          ,ICUTA,ICUTB,ICUTD,DEFO_MIN,DEFO_MAX,XX_CC )

            IF (AV_DEFO.LT.DEFO_MIN) THEN
               CALL DEFOCUS_GUESSING4(TMP,NMB,NT,NUMBER,N2,TMP_OUT,
     &                          IPW_ZERO,CS,LAMBDA,CONTRAST,
     &                          PS,AV_DEFO,DEFO_MIN,DEFO_MAX,XX_CC)
            ENDIF
            RETURN
        ENDIF

CC==    WE DO NOT WISH TOO HIGH POLYNOMIAL DEGREE !
        IF (NDGREE .GE. 6) THEN
           NDGREE=6
        ENDIF

      tmp_av_defo = av_defo  ! al feb 2013


CC==    RESET PLOT AND DO DEFOCUS REFINEMENT
        PLOT(:,1) = TMP(IPW_ZERO+1:NMB,1)-1.0
        PLOT(:,3) = ALOG(TMP(IPW_ZERO+1:NMB,2))

CC==    DEFOCUS REFINEMENT
       !write(6,*) 'in guessing1 calling defo_refine:',av_defo

        CALL DEFO_REFINE(PLOT,TMP,TMP_OUT,NMB,
     &    NUMBER,N2,NT,NDGREE,ICUTA,ICUTB,PS,CS,LAMBDA,CONTRAST,
     &    IPW_ZERO,DEFO_MIN,DEFO_MAX,AV_DEFO,XX_CC)
        DEALLOCATE(PLOT,TMP)

        !write(6,'(a,f10.1," --",f10.1)') 'Defocus:',tmp_av_defo, av_defo

        END

CC=====================================================================

        SUBROUTINE XFIT_BKGND3(K,N2,N,PLOT,ICUT1,ICUT2,ICNSTRNT)

CCC==INEQUALITY CONSTRAINED LINEAR SQUARE MINIMIZATION
CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4

        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)

        L=2
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF
CC==
        DO I=1,N-1
            Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO
        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
CC==GENERATE Q4 FOR UDR
        CALL PARTI11(KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICUT1,ICUT2,ICNSTRNT,KP)
        DEALLOCATE(Q2)
        END

CC=============

        SUBROUTINE DEFOCUS_EST(K,N2,N,PLOT,ICUT1,ICUT2,
     &      PS,CS,LAMBDA,CONTRAST,DEFOCUS,NZERO,
     &      IFIRST,ISECND)

C       TMP2 is automatic array
        REAL TMP2(K)
        REAL PLOT(K,N2)
        REAL LAMBDA
        REAL, DIMENSION(:), ALLOCATABLE :: TMP1

        ALLOCATE(TMP1(K),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
         ENDIF

        DO KKK=1,ICUT2
            XXX   = PLOT(KKK,1)/2./PS/REAL(K)
            XXXXX = EXP(PLOT(KKK,4))-EXP(PLOT(KKK,2))

            IF (PLOT(KKK,2) .NE. 0.0) THEN
          PLOT(KKK,2) = EXP(PLOT(KKK,2))
          PLOT(KKK,3) = EXP(PLOT(KKK,3))
          PLOT(KKK,4) = EXP(PLOT(KKK,4))
            ELSE
          PLOT(KKK,2) = 0.0
          PLOT(KKK,3) = 0.0
          PLOT(KKK,4) = 0.0
            ENDIF
        ENDDO

        PLOT(1:ICUT2,3) = PLOT(1:ICUT2,3)-PLOT(1:ICUT2,2)
        PLOT(1:ICUT2,4) = PLOT(1:ICUT2,4)-PLOT(1:ICUT2,2)
        TMP2(1:ICUT2)   = PLOT(1:ICUT2,4)

        ISWI = 1
        DZMAX = GETDEFOCUS(PLOT(1,3),PLOT(1,1),K,PS,CS,LAMBDA,
     &         CONTRAST,XSCORE,PLOT(1,4),ICUT1,ICUT2,ISWI)

CCCC==  OUTPUT RESULT 
         DEFOCUS = DZMAX

        CALL GET_CTF_ZERO(K,PS, CS,LAMBDA,
     &           DEFOCUS, CONTRAST,ICUT1,ICUT2,NZERO,IFIRST,ISECND)

        CALL CTF_SIGNAL1(TMP1, PLOT(1,1),K, PS, CS,
     &                  LAMBDA, DEFOCUS, CONTRAST)

        DEALLOCATE(TMP1)
        END

CCC=====================================================================

        SUBROUTINE PARTI15(KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT1,N,ICUT1,ICUT2,ICNSTRNT,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ

CC==LOW FREQUENCIES REGION RE-FITTING!!!
CC==INTERMDEDIATE REGION FITTING
CC====ZHONG HUANG, APR,16,03
        ISTART=ICUT1
        ISTOP=ICUT2
        N_3=N
        L=1
        KS=ISTOP-ISTART+1
        ISWI=1
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
CC== EQUALITY CONSTRAINT CONDITION IN TT1_Q2
        TT1_Q2(KS+1,1:N-1)=Q2(ICNSTRNT,1:N-1)
        TT1_Q2(KS+1,N)=1.0
        TT1_Q2(KS+1,N+1)=PLOT1(ICNSTRNT,3)

        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-1)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     & Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)

        END

CC=====================================================================

        SUBROUTINE XFIT_LOW4(K,N2,PLOT,ICUT1,ICUT2)

CCC==INEQUALITY CONSTRAINED LINEAR SQUARE MINIMIZATION
CC==FOR LOW FREQUENCIES FITTING,IT GIVES HIGH ORDER FITTING FROM 1 TO ICUT1
CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA

CC==BEGIN: FIRST ALLOCATE ARRAYS FOR CALCULATION
        N=5
        L=1
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF


CC==
        DO I=1,N-1
          Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
CC==
        Q2(K+1:2*K,1:N+1)=Q2(1:K,1:N+1)
CC==
        CALL PARTI17 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICUT1,ICUT2,KP)
        DEALLOCATE(Q2)
        END

CCC=======================================================

        SUBROUTINE XFIT_ZH_AST(K,N2,PLOT,ICUT1,ICUT2)

CCC==INEQUALITY CONSTRAINED LINEAR SQUARE MINIMIZATION
CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4

        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA

        N=7
        L=2
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF
CC==
        DO I=1,N-1
            Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
CC==
        Q2(K+1:2*K,1:N+1)=Q2(1:K,1:N+1)
CC==GENERATE Q4 FOR UDR
        CALL PARTI8 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICUT1,ICUT2,KP)
        DEALLOCATE(Q2)

        END

CC===================================================================

        SUBROUTINE PARTI17 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT1,N,ICUT1,ICUT2,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ

CC==LOW FREQUENCIES REGION RE-FITTING!!!
CC==INTERMDEDIATE REGION FITTING
        ISTART=ICUT1
        ISTOP=ICUT2
        N_3=N
        KS=ISTOP-ISTART+1
        L=2
        ISWI=4
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
CC== EQUALITY CONSTRAINT CONDITION IN TT1_Q2
        TT1_Q2(KS+1,1:N-1)=Q2(ISTART,1:N-1)
        TT1_Q2(KS+2,1:N-1)=Q2(ISTOP,1:N-1)
        TT1_Q2(KS+1,N)=1.0
        TT1_Q2(KS+2,N)=1.0
        TT1_Q2(KS+1,N+1)=PLOT1(ISTART,3)
        TT1_Q2(KS+2,N+1)=PLOT1(ISTOP,4)

        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-3)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
        END

CCC=================================================================

        SUBROUTINE XFIT_LOW5(K,N2,N,PLOT,ICUT1,ICUT2)

CCC==INEQUALITY CONSTRAINED LINEAR SQUARE MINIMIZATION
CC==FOR LOW FREQUENCIES FITTING,IT GIVES HIGH ORDER FITTING FROM 1 TO ICUT1
CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA

CC==BEGIN: FIRST ALLOCATE ARRAYS FOR CALCULATION
        L=1
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF
CC==
        DO I=1,N-1
            Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
CC==
        Q2(K+1:2*K,1:N+1)=Q2(1:K,1:N+1)
CC==GENERATE Q4 FOR UDR
        CALL PARTI18 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICUT1,ICUT2,KP)
        DEALLOCATE(Q2)
        END

CC==========================================================

        SUBROUTINE PARTI18 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT1,N,ICUT1,ICUT2,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ

CC==LOW FREQUENCIES REGION RE-FITTING!!!
CC==INTERMDEDIATE REGION FITTING
        ISTART=ICUT1
        ISTOP=ICUT2
CC====ZHONG HUANG, APR,16,03
        N_3=N
        KS=ISTOP-ISTART+1
        L=1
        ISWI=5
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
CC== EQUALITY CONSTRAINT CONDITION IN TT1_Q2
        TT1_Q2(KS+1,1:N-1)=Q2(ISTART,1:N-1)
        TT1_Q2(KS+1,N)=1.0
        TT1_Q2(KS+1,N+1)=PLOT1(ISTART,3)
        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-4)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
        END

CC============================================================

        SUBROUTINE XFIT_EN(NUMBER,N2,PLOT,
     &           PS,CS,LAMBDA,CONTRAST,DEFOCUS,IPT1,IPT2,IPT3,IPT4)

CC==RE-FIT CTF ZEROS AND PEAKS AFTER OBTAINING DEFOCUS, AND USE AS CONSTRAINT
CC==POINTS TO GET MORE SMOOTH EVELOPE FUCNTION
CC==ESTIMATE CTF ZERO
        INTEGER, DIMENSION(:), ALLOCATABLE :: FOUR_ZERO,NCON
        REAL PLOT(NUMBER,N2)
        REAL LAMBDA

        N  = 5
        N4 = 4
        ALLOCATE(FOUR_ZERO(N4),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',N4)
            RETURN
        ENDIF
        N0=0

       !! defocus bad!!!!!!!!!!

       !write(6,*) '88',CS,    LAMBDA, CONTRAST,  PS,  DEFOCUS,  NZERO 
c                     2.26E+7 1.96E-2  0.1     1.58   454.1191  32767        

        CALL GET_CTF_ZERO4(NUMBER,N0,PS, CS,
     &           LAMBDA, DEFOCUS, CONTRAST,NZERO,N4,FOUR_ZERO)

       !write(6,*) '89',CS,  LAMBDA, CONTRAST,  PS,  DEFOCUS,  NZERO 
c89               2.26E+7   1.96E-2  0.1     1.58  22639.0    44
c89               2.26      1.96     0.1     1.58  454.1    1

        IF (NZERO.GE.5) THEN
            N3 = N4-1
        ELSE
            N3=NZERO-1
        ENDIF

c       BUG zero allocation here for ncon  feb 2013 al
c       n3 = 0  ,  n4 = 4,  nzero=1
        ntoalloc = n4-1

      !write(6,*) 'in xfit_en, nzero,n3:',nzero,n3,n4,four_zero
      !write(6,*) 'in xfit_en, allocating ncon :',ntoalloc,'==',n3

C	ALLOCATE(NCON(N3),STAT=IRTFLG)
        ALLOCATE(NCON(NTOALLOC),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED, NCON',NTOALLOC)
            RETURN
        ENDIF
        DO II=1,N3
            III=II+1
            IFIT1=FOUR_ZERO(II)
            IFIT2=FOUR_ZERO(III)

            IF(II.EQ.1) THEN
          X_FIT=REAL(IFIT2-IFIT1)
          X_FIT=X_FIT/4.0
          IFIT1=INT(X_FIT)+IFIT1
            ENDIF
          CALL XFITLINE2(NUMBER,N2,PLOT,IFIT1,IFIT2,IFIT3)
          NCON(II)=IFIT3
        ENDDO
CC==
        IPT1=FOUR_ZERO(2)
        IPT2=FOUR_ZERO(3)
        IPT3=NCON(1)
        IPT4=NCON(2)

        CALL XFIT_EN4(NUMBER,N2,N4,PLOT,NCON,N,NZERO)
        DEALLOCATE(FOUR_ZERO,NCON)
        END

CC=====================

        SUBROUTINE GET_CTF_ZERO4(NUMBER,N0,PS, CS,
     &           LAMBDA, DEFOCUS, CONTRAST,NZERO,N4,FOUR_ZERO)

        COMMON /SEARCH_RANGE/ DZLOW1,DZHIGH1,NMB
        PARAMETER (QUADPI = 3.1415926535897932384626)
         REAL LAMBDA,DZLOW1,DZHIGH1
        INTEGER FOUR_ZERO(N4)

       !write(6,*) '04',CS,  LAMBDA, CONTRAST,  PS,  DEFOCUS,  NZERO 
       !write(6,*) '04',DZLOW1,DZHIGH1,NMB 
C                  04  3000.0  200000.0  250   defocus bad!!!!

       !write(6,*) '04',CS,  LAMBDA, CONTRAST,  PS,  DEFOCUS,  NZERO 

CC==GIVE FIRST FOUR CTF ZEROS ACCORDING THE GIVEN DEFOCUS
         AA    = -0.5*QUADPI*CS*LAMBDA**3
         BB    = QUADPI*LAMBDA
         BB    = BB*DEFOCUS             !bad?????
         CC    = -ATAN(CONTRAST/(1.0-CONTRAST))

        FRQ1  = 0.0
        FRQ2  = REAL(NUMBER+N0) / 2.0 / PS / REAL(NMB)

        !write(6,*) '04, aa,bb,cc,defocus',aa,bb,cc,defocus
        !write(6,*) '04, frq1,2',frq1,frq2

         TMPN1 = AA*FRQ1**4+BB*FRQ1**2+CC
         TMPN1 = TMPN1 / QUADPI
         N1    = INT(TMPN1)

        TMPN2 = AA*FRQ2**4+BB*FRQ2**2+CC   ! bad!!!!!!!!
         TMPN2 = TMPN2 / QUADPI
         N2    = INT(TMPN2)

        !write(6,*) '04, TMPN2,1',TMPN2,TMPN1

        N     = N2 - N1 + 1
        TMP1  = 4.*AA

        !write(6,*) '04, n2,n1,n',n2,n1,n 


        DO II=1,N
            TMP2 = CC - REAL(II-1)*QUADPI
            TMP2 = SQRT(BB**2-TMP1*TMP2)
            TMP2 = TMP2-BB
            TMP2 = SQRT(TMP2/2.0/AA)
            NFRQ = INT(TMP2*2.0*PS*NMB)

            IF(II .LE. 4) THEN
          FOUR_ZERO(II) = NFRQ
            ENDIF

        ENDDO

        IF(N .LT. 4) THEN
            NN = 4-N
            DO II=1,NN
          NNN            = 4-II+1
          FOUR_ZERO(NNN) = 0
            ENDDO
        ENDIF
        NZERO = N

3000    RETURN
        END
CC============================================================================

        SUBROUTINE XFIT_EN4(K,N2,N3,PLOT,NCON,N,NZERO)

CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        INTEGER NCON(N3)
CC==BEGIN: FIRST ALLOCATE ARRAYS FOR CALCULATION
        IF(NZERO.GE.5) THEN
            CALL EN_FIT7(PLOT,N,NCON,K,N2,N3)
        ELSEIF(NZERO.EQ.4) THEN
            L=3
            NT_2=3
            CALL EN_FIT6(PLOT,NT_2,NCON,K,N2,N3,L)
            NT_1=3
            CALL EN_FIT5(PLOT,NT_1,NCON,K,N2,N3,L)
        ELSEIF(NZERO.EQ.3) THEN
            L=2
            NT_3=3
            CALL EN_FIT5(PLOT,NT_3,NCON,K,N2,N3,L)
            IFIT1=1
            IFIT2=NCON(1)
            CALL XFITLINE3(K,N2,PLOT,IFIT1,IFIT2)
        ELSEIF(NZERO.LE.2) THEN
            IF(NCON(1).EQ.0) RETURN
            CALL EN_FIT8(PLOT,NCON,K,N2,N4)
        ENDIF
CC--EXTEND ENVELOPE TO LOW FREQUENCIES
          CALL EN_LOW(PLOT,N2,K,NCON,N4)
        END

CC===========================================================

        SUBROUTINE PARTI19 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT1,N,ICUT1,ICUT2,ICUTD,ICUTC,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ

CC==LOW FREQUENCIES REGION RE-FITTING!!!
CC==INTERMDEDIATE REGION FITTING
        ISTART=ICUT1
        ISTOP=ICUT2
CC====ZHONG HUANG, APR,16,03
        N_3=N
        KS=ISTOP-ISTART+1
        L=4
        ISWI=8
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
CC== EQUALITY CONSTRAINT CONDITION IN TT1_Q2
        TT1_Q2(KS+1,1:N-1)=Q2(ISTART,1:N-1)
        TT1_Q2(KS+2,1:N-1)=Q2(ISTOP,1:N-1)
        TT1_Q2(KS+3,1:N-1)=Q2(ICUTD,1:N-1)
        TT1_Q2(KS+4,1:N-1)=Q2(ICUTC,1:N-1)
        TT1_Q2(KS+1,N)=1.0
        TT1_Q2(KS+2,N)=1.0
        TT1_Q2(KS+3,N)=1.0
        TT1_Q2(KS+4,N)=1.0
        TT1_Q2(KS+1,N+1)=PLOT1(ISTART,3)
        TT1_Q2(KS+2,N+1)=PLOT1(ISTOP,3)
        TT1_Q2(KS+3,N+1)=PLOT1(ICUTD,3)
        TT1_Q2(KS+4,N+1)=PLOT1(ICUTC,3)
CC==
        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-7)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
        END
CC===================================================================

        SUBROUTINE EN_FIT5(PLOT,N,NCON,K,N2,N4,L)

        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
        REAL PLOT(K,N2)
        INTEGER NCON(N4)

        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF

CC==
        DO I=1,N-1
          Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
CC==
        Q2(K+1:2*K,1:N+1)=Q2(1:K,1:N+1)
CC==GENERATE Q4 FOR UDR
        ICUT1=NCON(1)
        ICUT2=K-1
        ICON3=NCON(2)
        ICON4=NCON(3)
        IF(L.EQ.4) THEN
            CALL PARTI19 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICUT1,ICUT2,ICON3,ICON4,KP)
        ELSEIF(L.EQ.3) THEN
            CALL PARTI17 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICUT1,ICON3,KP)
        ELSEIF(L.EQ.2) THEN
            CALL PARTI20 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICUT1,ICUT2,KP)
        ENDIF
        DEALLOCATE(Q2)
        END
CC============================================================

        SUBROUTINE EN_FIT6(PLOT,N,NCON,K,N2,N4,L)

        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
        REAL PLOT(K,N2)
        INTEGER NCON(N4)

        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF

CC==
        DO I=1,N-1
            Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
CC==
        Q2(K+1:2*K,1:N+1)=Q2(1:K,1:N+1)
CC==GENERATE Q4 FOR UDR
        ICUT1=NCON(1)
        ICUT2=K-1
        ICON3=NCON(2)
        ICON4=NCON(3)
        IF(L.EQ.4) THEN
            CALL PARTI19 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICUT1,ICUT2,ICON3,ICON4,KP)
        ELSEIF(L.EQ.3) THEN
            CALL PARTI20 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICON3,ICUT2,KP)
        ELSEIF(L.EQ.2) THEN
            CALL PARTI20 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICUT1,ICUT2,KP)
        ENDIF
        DEALLOCATE(Q2)
        END

CCC====================================================================

        SUBROUTINE PARTI20 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT1,N,ICUT1,ICUT2,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ

CC==LOW FREQUENCIES REGION RE-FITTING!!!
CC==INTERMDEDIATE REGION FITTING
        ISTART=ICUT1
        ISTOP=ICUT2
CC====ZHONG HUANG, APR,16,03
        N_3=N
        KS=ISTOP-ISTART+1
        L=0
        ISWI=2
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-1)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
        END

CC=======================================================================

        SUBROUTINE EN_FIT7(PLOT,N,NCON,K,N2,N4)

        REAL PLOT(K,N2)
        INTEGER NCON(N4)

        ICUT1=NCON(1)
        ICUT2=K-1
        ICON3=NCON(2)
        ICON4=NCON(3)
        N_2=11
        CALL EN_FIT_THREE(PLOT,N_2,K,N2,ICUT1,ICUT2,ICON4)

        END

CC=======================================================

        SUBROUTINE EN_FIT_TWO(PLOT,N,K,N2,IPT1,IPT2)

        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
        REAL PLOT(K,N2)
CC==
        L=1
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF
CC==
        DO I=1,N-1
            Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO
        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
CC==
        Q2(K+1:2*K,1:N+1)=Q2(1:K,1:N+1)
CC==GENERATE Q4 FOR UDR
        CALL PARTI21 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,IPT1,IPT2,KP)
        DEALLOCATE(Q2)
        END

CC=========================================================

        SUBROUTINE EN_FIT_THREE(PLOT,N,K,N2,IPT1,IPT2,IPT3)

        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
        REAL PLOT(K,N2)

        L=0
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF
CC==
        DO I=1,N-1
            Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO
        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
CC==
        Q2(K+1:2*K,1:N+1)=Q2(1:K,1:N+1)
CC==GENERATE Q4 FOR UDR
        CALL PARTI20 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,IPT1,IPT2,KP)
        DEALLOCATE(Q2)
        END

CC========================================================

        SUBROUTINE PARTI21 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &                      K,N2,PLOT1,N,ICUT1,ICUT2,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ

CC==LOW FREQUENCIES REGION RE-FITTING!!!
CC==INTERMDEDIATE REGION FITTING
        ISTART=ICUT1
        ISTOP=ICUT2
CC====ZHONG HUANG, APR,16,03
        N_3=N
        KS=ISTOP-ISTART+1
        L=1
        ISWI=5
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
CC== EQUALITY CONSTRAINT CONDITION IN TT1_Q2
        TT1_Q2(KS+1,1:N-1)=Q2(ISTOP,1:N-1)
        TT1_Q2(KS+1,N)=1.0
        TT1_Q2(KS+1,N+1)=PLOT1(ISTOP,4)
        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-4)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
        END

CC===============================================================

        SUBROUTINE EN_LOW(PLOT,N2,K,NCON,N4)

        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
        REAL    PLOT(K,N2)
        INTEGER NCON(N4)
        N=3
        L=1
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',KLM2D*N2D)
         RETURN
         ENDIF
CC==
        DO I=1,N-1
            Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO
        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
CC==
        Q2(K+1:2*K,1:N+1)=Q2(1:K,1:N+1)
CC==GENERATE Q4 FOR UDR
        IPT1=1
        IPT2=NCON(1)   ! is undefined here!!!!
        CALL PARTI22 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &                K,N2,PLOT,N,IPT1,IPT2,KP)
        DEALLOCATE(Q2)
        END

CCC====================================================================

        SUBROUTINE PARTI22(KLMD,NKLMD,KLM2D,N2D,Q2,
     &                     K,N2,PLOT1,N,IPT1,IPT2,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC-TT1_Q2, Q,X,RES,CU,S,and IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ

CC==LOW FREQUENCIES REGION RE-FITTING!!!
CC==INTERMDEDIATE REGION FITTING
        ISTART = IPT1
        ISTOP  = IPT2
CC====ZHONG HUANG, APR,16,03
        N_3=N
        KS=ISTOP-ISTART+1
        L=1
        ISWI=5

      !write(6,*)'in parti22, istart,istop,bnd:',istart,istop,k

      !write(6,*)'in parti22, ks,n,klm2d:',ks,n,klm2d

        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)

CC== EQUALITY CONSTRAINT CONDITION IN TT1_Q2
        TT1_Q2(KS+1,1:N-1)= Q2(ISTOP,1:N-1)
        TT1_Q2(KS+1,N)    = 1.0
        TT1_Q2(KS+1,N+1)  = PLOT1(ISTOP,4)

        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-4)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
        END

CC===============================================================

        SUBROUTINE PARTI23(KLMD,NKLMD,KLM2D,N2D,Q2,
     &                     K,N2,PLOT1,N,IPT1,IPT2,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ

CC==LOW FREQUENCIES REGION RE-FITTING!!!
CC==INTERMDEDIATE REGION FITTING
        ISTART=IPT1
        ISTOP=IPT2
        N_3=N
        KS=ISTOP-ISTART+1
        L=1
        ISWI=5
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
CC== EQUALITY CONSTRAINT CONDITION IN TT1_Q2
        TT1_Q2(KS+1,1:N-1)=Q2(ISTART,1:N-1)
        TT1_Q2(KS+1,N)=1.0
        TT1_Q2(KS+1,N+1)=PLOT1(ISTART,3)
        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-4)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
        END

CC===============================================================

        SUBROUTINE XFIT_LOW6(K,N2,PLOT,IPT,IPT1,IPT2,IPT3)

CCC==INEQUALITY CONSTRAINED LINEAR SQUARE MINIMIZATION
CC==FOR LOW FREQUENCIES FITTING,IT GIVES HIGH ORDER FITTING FROM 1 TO ICUT1
CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA

CC==BEGIN: FIRST ALLOCATE ARRAYS FOR CALCULATION
        ICUT1=IPT1
        ICUT2=IPT2
        N=4
        L=0
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF
CC==
        DO I=1,N-1
          Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
CC==
C        Q2(K+1:2*K,1:N+1)=Q2(1:K,1:N+1)
CC==GENERATE Q4 FOR UDR
        CALL PARTI24 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,ICUT1,ICUT2,KP,N)
        DEALLOCATE(Q2)
        END

CC======================================================

        SUBROUTINE PARTI24 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT1,ICUT1,ICUT2,KP,N)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ

CC==LOW FREQUENCIES REGION RE-FITTING!!!
CC==INTERMDEDIATE REGION FITTING
        ISTART=ICUT1
        ISTOP=ICUT2
CC==
        N_3=N
        KS=ICUT2-ICUT1+1
        L=0
        ISWI=1
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
CC== EQUALITY CONSTRAINT CONDITION IN TT1_Q2
        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-1)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
        END
CC==============================================================

        SUBROUTINE XFIT_LOW7(K,N2,PLOT,IPT1,IPT2)

CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA

CC==BEGIN: FIRST ALLOCATE ARRAYS FOR CALCULATION
        ICUT1=IPT1
        ICUT2=IPT2
        N=3
        L=0
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF
CC==
        DO I=1,N-1
          Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO
        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
CC==GENERATE Q4 FOR UDR

        CALL PARTI24 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,ICUT1,ICUT2,KP,N)
        DEALLOCATE(Q2)
        END

CC======================================================

        SUBROUTINE XFITLINE3(K,N2,PLOT,IFIT1,IFIT2)
CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)

CC==IN THIS CASE, K=NMB
        N=2
        L=0
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG

        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF

        Q2(1:K,1)=PLOT(1:K,1)
        Q2(1:K,2)=1.0
        Q2(1:K,3)=PLOT(1:K,3)
CC==
        CALL PARTI9(KLMD,NKLMD,KLM2D,N2D,Q2,
     &              K,N2,N,PLOT,IFIT1,IFIT2,KP)
        DEALLOCATE(Q2)
        END

CC=================================================================

        SUBROUTINE XFIT_LOW8(K,N2,PLOT,IPT1,IPT2)
CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA

CC==BEGIN: FIRST ALLOCATE ARRAYS FOR CALCULATION
        ICUT1=IPT1
        ICUT2=IPT2
        N=9
        L=0
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG

        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF

        DO I=1,N-1
            Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
CC==GENERATE Q4 FOR UDR
        CALL PARTI25 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,ICUT1,ICUT2,KP)
        DEALLOCATE(Q2)
        END

CCC========================================================

        SUBROUTINE PARTI25 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT1,ICUT1,ICUT2,KP)
        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S, and IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ

CC==LOW FREQUENCIES REGION RE-FITTING!!!
CC==INTERMDEDIATE REGION FITTING
        N=9
        ISTART=ICUT1
        ISTOP=ICUT2
CC==
        N_3=N
        KS=ICUT2-ICUT1+1
        L=0
        ISWI=1
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
CC== EQUALITY CONSTRAINT CONDITION IN TT1_Q2
        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-1)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
        END

CC==================================================

        SUBROUTINE RESOLUTION_ANALYZE(NUMBER,N2,N0,PS, CS,
     &           LAMBDA, DEFOCUS, CONTRAST,PLOT1,XX_CC)
        COMMON /SEARCH_RANGE/ DZLOW1,DZHIGH1,NMB
        PARAMETER (QUADPI = 3.1415926535897932384626)
         REAL LAMBDA,DZLOW1,DZHIGH1
        REAL PLOT1(NUMBER,N2)
CC==GIVE FIRST FOUR CTF ZEROS ACCORDING THE GIVEN DEFOCUS
        REAL,DIMENSION (:,:),ALLOCATABLE :: TMP_1,X_ITGL
        REAL,DIMENSION (:),ALLOCATABLE :: TMP_2,TMP_3,TMP_NZERO_M,
     &                                    X_CC,X_PEAK
        INTEGER, DIMENSION (:), ALLOCATABLE :: NZERO
        X_CCTHR=0.8
        XX_CC=9999.
        AA = -0.5*QUADPI*CS*LAMBDA**3
         BB = QUADPI*LAMBDA
         BB=BB*DEFOCUS
         CC = -ATAN(CONTRAST/(1.0-CONTRAST))
        FRQ1=0.0
        FRQ2=REAL(NUMBER+N0)/2.0/PS/REAL(NMB)
         TMPN1=AA*FRQ1**4+BB*FRQ1**2+CC
         TMPN1=TMPN1/QUADPI
         N_1=INT(TMPN1)
        TMPN2=AA*FRQ2**4+BB*FRQ2**2+CC
         TMPN2=TMPN2/QUADPI
         N_2=INT(TMPN2)
        N_X=N_2-N_1+1

        IF(N_X.LE.1) THEN
            XX_CC=9999
            RETURN
        ENDIF

        N_XXX=N_X-1
        TMP1=4.*AA

        ALLOCATE(NZERO(N_X),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
        ENDIF
        ALLOCATE(TMP_NZERO_M(N_X),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
        ENDIF
        ALLOCATE(X_CC(N_X),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
        ENDIF
        ALLOCATE(X_PEAK(N_X),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
        ENDIF

        DO II=1,N_X
            TMP2=CC-REAL(II-1)*QUADPI
            TMP2=SQRT(BB**2-TMP1*TMP2)
            TMP2=TMP2-BB
            TMP2=SQRT(TMP2/2.0/AA)
            NFRQ=INT(TMP2*2.0*PS*NMB)
            NZERO(II)=NFRQ
            TMP3=CC-REAL(II-1)*QUADPI-.5*QUADPI
            TMP_XX=abs(BB**2-TMP1*TMP3)
            TMP3=SQRT(TMP_XX)
            TMP3=TMP3-BB
            TMP3=SQRT(TMP3/2.0/AA)
            TMP_NZERO_M(II)=TMP3
        ENDDO
CC
        ALLOCATE(TMP_1(NUMBER,4),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
        ENDIF
        ALLOCATE(TMP_2(NUMBER),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
        ENDIF
        ALLOCATE(X_ITGL(N_X,3),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
        ENDIF
        ALLOCATE(TMP_3(NUMBER),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
        ENDIF
        TMP_1(:,:)=PLOT1(:,:)
        DO KKK=1,NUMBER
            XXXXX=EXP(TMP_1(KKK,4))-EXP(TMP_1(KKK,2))
            IF(TMP_1(KKK,2).NE.0.0) THEN
          TMP_1(KKK,2)=EXP(TMP_1(KKK,2))
          TMP_1(KKK,3)=EXP(TMP_1(KKK,3))
          TMP_1(KKK,4)=EXP(TMP_1(KKK,4))
            ELSE
          TMP_1(KKK,2)=0.0
          TMP_1(KKK,3)=0.0
          TMP_1(KKK,4)=0.0
            ENDIF
        ENDDO

        TMP_1(:,3)=TMP_1(:,3)-TMP_1(:,2)
        TMP_1(:,4)=TMP_1(:,4)-TMP_1(:,2)
CC==SEARCHING NZEROS!
        CALL CTF_SIGNAL1(TMP_2,TMP_1(1,1),NUMBER, PS, CS,
     &                LAMBDA, DEFOCUS, CONTRAST)
        TMP_2(1:NUMBER)=TMP_2(1:NUMBER)**2*TMP_1(1:NUMBER,4)

        DO II=1,NUMBER
            TMP_3(II)=SQRT(ABS(TMP_1(II,3)))*SQRT(ABS(TMP_2(II)))
        ENDDO

        IF(N_XXX.LE.1) THEN
CC==WE CAN EXTRACT INFORMATION TO RESOLUTION LIMIT!
           XX_CC=1.0/(2.*PS)
           RETURN
        ENDIF

        DO III=1,N_XXX
            XX=TMP_NZERO_M(III)*PS*2.0*NUMBER
            I_XX=INT(XX)
            X_PEAK(III)=TMP_1(I_XX,3)
            IIII=III+1
            III_1=NZERO(III)
            III_2=NZERO(IIII)
            N_XX=III_2-III_1+1
            CALL POWER_SUM (TMP_2(III_1),N_XX,X_ITGL(III,1))
            CALL POWER_SUM (TMP_1(III_1,3),N_XX,X_ITGL(III,2))
            CALL POWER_SUM (TMP_3(III_1),N_XX,X_ITGL(III,3))
            CALL CC_TWO_CURVES(TMP_1(III_1,3),
     &           TMP_2(III_1),N_XX,X_CC(III))

            IF(X_CC(III).GE.X_CCTHR) THEN
CC==WE SET THRESHOLD AS 0.8
          XX_CC=(TMP_NZERO_M(III)+TMP_NZERO_M(IIII))/2.
            ENDIF

        ENDDO
        DEALLOCATE(TMP_1,NZERO,TMP_2,TMP_3,X_ITGL,TMP_NZERO_M,
     &              X_CC,X_PEAK)
        END

CC===================================================================

        SUBROUTINE CC_TWO_CURVES(TMP_X,TMP_Y,N_XX,X_CC)

        REAL TMP_X(N_XX),TMP_Y(N_XX)

        TMP1=0.0
        TMP2=0.0
        TMP3=0.0
        TMP4=0.0
        TMP5=0.0
CC==
        TMPX=0.0
        TMPY=0.0
        DO II=1,N_XX
            TMPX=TMPX+TMP_X(II)
            TMPY=TMPY+TMP_Y(II)
        ENDDO

        IF(TMPX.EQ.0.0.OR.TMPY.EQ.0.0) THEN
            X_CC=0
            RETURN
        ENDIF

        DO II=1,N_XX
            TMP1=TMP1+TMP_X(II)
            TMP2=TMP2+TMP_Y(II)
            TMP3=TMP3+TMP_X(II)*TMP_X(II)
            TMP4=TMP4+TMP_Y(II)*TMP_Y(II)
            TMP5=TMP5+SQRT(ABS(TMP_X(II)*TMP_Y(II)))
        ENDDO

        TMP6=REAL(N_XX)
        TMP1=TMP1/TMP6
        TMP2=TMP2/TMP6
        X_CC=TMP5/TMP6
        X_CC=X_CC/SQRT(ABS(TMP1))/SQRT(ABS(TMP2))
        END

CC=========================================================================

        SUBROUTINE CTF_2_ZERO(TMP,NMB,NT,NUMBER,N2,TMP_OUT,
     &      IPW_ZERO,CS,LAMBDA,CONTRAST,PS,NDGREE,AV_DEFO,
     &      ICUTA,ICUTB,ICUTD,DEFO_MIN,DEFO_MAX,XX_CC)

        REAL LAMBDA
        DIMENSION TMP(NMB,2),TMP_OUT(NMB,NT)
        REAL, DIMENSION(:,:), ALLOCATABLE :: PLOT,TMP_PLOT

        ALLOCATE(PLOT(NUMBER,N2),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
        ENDIF

        ALLOCATE(TMP_PLOT(NUMBER,N2),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
        ENDIF

        PLOT(:,1)=TMP(IPW_ZERO+1:NMB,1)-1.0
        PLOT(:,3)=ALOG(TMP(IPW_ZERO+1:NMB,2))
        N=NDGREE+1
CC==CC==FIT BACKGROUND FROM ICUTB TO ICUT2
        ICUT1=ICUTA
        ICUT2=NUMBER
        CALL XFIT_ZH(NUMBER,N2,N,PLOT,ICUT1,ICUT2)

        N_TEST=NDGREE+1

       !write(6,*) 'in CTF_2_ZERO calling XFIT_LOW1,icutd:',icutd

        CALL XFIT_LOW1(NUMBER,N2,PLOT,ICUTA,ICUT2,ICUTD,N_TEST)

        ICUT1=ICUTB
        ICNSTRNT=ICUTB
        IBACK=1
        CALL XFIT_LOW2(NUMBER,N2,PLOT,ICUTB,ICUT2,N)
        CALL XFIT_LOW10(NUMBER,N2,PLOT,IBACK,ICUTB)

        ICUT_LOW=1
        ICUT_HIGH=ICUTB
        CALL XFIT_LOW9(NUMBER,N2,PLOT,ICUT_LOW,ICUT_HIGH)

       !write(6,*) '2zero',CS,       LAMBDA, CONTRAST,  PS,  AV_DEFO 
c                 2.26E+7  1.96E-2   0.1    1.584   454.1


        CALL XFIT_EN(NUMBER,N2,PLOT,
     &           PS,CS,LAMBDA,CONTRAST,AV_DEFO,IPT1,IPT2,IPT3,IPT4)

        IF (IPT1 .EQ. 0) RETURN
CC====
        CALL XFIT_BKGND3(NUMBER,N2,N,PLOT,ICUT_LOW,ICUT2,ICNSTRNT)
        IFIT_B_STOP=ICUT2
        TMP_PLOT(:,:)=PLOT(:,:)
        TMP_OUT(1:NUMBER,:)=PLOT(:,:)

        CALL DEFOCUS_EST1(NUMBER,N2,N,PLOT,ICUT1,ICUT2,
     &           PS,CS,LAMBDA,CONTRAST,DEFOCUS)

        IF(DEFOCUS.GT.DEFO_MAX) THEN
            AV_DEFO=DEFO_MAX
            XX_CC=9999
            GO TO 1122
        ELSE
            AV_DEFO=DEFOCUS
        ENDIF
        N0=0
        CALL RESOLUTION_ANALYZE(NUMBER,N2,N0,PS, CS,
     &           LAMBDA, AV_DEFO, CONTRAST,TMP_PLOT,XX_CC)
1122    DEALLOCATE(PLOT,TMP_PLOT)
        END

CC============================================================================

        SUBROUTINE EN_FIT8(PLOT,NCON,K,N2,N4)

        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
        REAL PLOT(K,N2)
        INTEGER NCON(N4)

        L=0
        N=5
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG

        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF

        DO I=1,N-1
          Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
        Q2(K+1:2*K,1:N+1)=Q2(1:K,1:N+1)
CC==GENERATE Q4 FOR UDR
        IF(ICUT1.EQ.0) RETURN

        ICUT1=NCON(1)
        ICUT2=K
        ICON3=NCON(2)
        ICON4=NCON(3)

        CALL PARTI20 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICUT1,ICUT2,KP)
        DEALLOCATE(Q2)
        END

CCC=============================================================

        SUBROUTINE XFIT_LOW9(K,N2,PLOT,IPT1,IPT2)

CCC==INEQUALITY CONSTRAINED LINEAR SQUARE MINIMIZATION
CC==FOR LOW FREQUENCIES FITTING,IT GIVES HIGH ORDER FITTING FROM 1 TO ICUT1
CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA

CC==BEGIN: FIRST ALLOCATE ARRAYS FOR CALCULATION
        ICUT1=IPT1
        ICUT2=IPT2
        N=6
        L=0
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF
CC==
        DO I=1,N-1
          Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)

        CALL PARTI24 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,ICUT1,ICUT2,KP,N)
        DEALLOCATE(Q2)
        END

CC=========================================================================

        SUBROUTINE XFIT_HIGH_EN(K,N2,PLOT,IPT1,IPT2)

CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA

CC==BEGIN: FIRST ALLOCATE ARRAYS FOR CALCULATION
        ICUT1=IPT1
        ICUT2=IPT2
        N=7
        L=0
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG

        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF
        DO I=1,N-1
            Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
CC==GENERATE Q4 FOR UDR
        CALL PARTI20 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICUT1,ICUT2,KP)
        DEALLOCATE(Q2)
        END

CC===============================================================================

        SUBROUTINE CTF_CORRECT(TMP,NMB,NT,NUMBER,N2,TMP_OUT,
     &           IPW_ZERO,CS,LAMBDA,CONTRAST,PS,AV_DEFO,
     &           ICUTA,ICUTB,ICUTD,DEFO_MIN,DEFO_MAX ,XX_CC)

CC==THIS SUBROUTINE CORRECT MIS_ESTIMATED DUE TO LONG AND FLAT HIGH FREQUENCY TAIL
        REAL LAMBDA
        DIMENSION TMP(NMB,2),TMP_OUT(NMB,NT)
        REAL, DIMENSION(:,:), ALLOCATABLE :: PLOT,TMP_PLOT

        ALLOCATE(PLOT(NUMBER,N2),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
         ENDIF

        ALLOCATE(TMP_PLOT(NUMBER,N2),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
         ENDIF
        PLOT(:,1)=TMP(IPW_ZERO+1:NMB,1)-1.0
        PLOT(:,3)=ALOG(TMP(IPW_ZERO+1:NMB,2))
        N=NDGREE+1
CC==CC==FIT BACKGROUND FROM ICUTB TO ICUT2
        ICUT1=ICUTB
        X_NUM=REAL(NUMBER)/2.
        ICUT2=INT(X_NUM)
CC==ICUTB is too far!
        IF(ICUT2.LE.ICUT1) THEN
            XX_CC=9999
            GO TO 1122
        ENDIF

        CALL FIT_ALL(NUMBER,N2,PLOT)

        TMP_OUT(1:NUMBER,:)=PLOT(:,:)
        TMP_PLOT(:,:)=PLOT(:,:)

        CALL DEFOCUS_EST(NUMBER,N2,N,PLOT,ICUT1,ICUT2,
     &           PS,CS,LAMBDA,CONTRAST,DEFOCUS,NZERO,
     &           IFIRST,ISECND)

        IF(DEFOCUS.GT.DEFO_MAX) THEN
            XX_CC=9999
            GO TO 1122
        ELSE
            AV_DEFO=DEFOCUS
        ENDIF
        N0=0
        CALL RESOLUTION_ANALYZE(NUMBER,N2,N0,PS,CS,
     &           LAMBDA,AV_DEFO,CONTRAST,TMP_PLOT,XX_CC)

1122    DEALLOCATE(PLOT,TMP_PLOT)
        END

CC===================================================================

        SUBROUTINE FIT_ALL(K,N2,PLOT)
CCC==INEQUALITY CONSTRAINED LINEAR SQUARE MINIMIZATION
CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4

        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA

        IFIT1=1
        IFIT2=K
        N=6
        L=0
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF
CC==
        DO I=1,N-1
          Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
        Q2(K+1:2*K,1:N+1)=Q2(1:K,1:N+1)

        CALL PARTI9(KLMD,NKLMD,KLM2D,N2D,Q2,
     &              K,N2,N,PLOT,IFIT1,IFIT2,KP)
        DEALLOCATE(Q2)
        END

CC====================================================================

        SUBROUTINE DEFO_REFINE(PLOT,TMP,TMP_OUT,NMB,NUMBER,
     &       N2,NT,NDGREE,ICUTA,ICUTB,PS,CS,LAMBDA,CONTRAST,IPW_ZERO,
     &       DEFO_MIN,DEFO_MAX,AV_DEFO,XX_CC)

        REAL LAMBDA
        DIMENSION TMP(NMB,2),PLOT(NUMBER,N2),TMP_OUT(NMB,NT)
        REAL, DIMENSION(:,:), ALLOCATABLE :: TMP_PLOT

        ALLOCATE(TMP_PLOT(NUMBER,N2),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
        ENDIF

        N=NDGREE+1
        ICUT1=ICUTB
        ICUT2=NUMBER

c       ICUTD was undefined until bb set it to 0 in 2009,
c       but that is also wrong since it is 'out of bounds'.
c       So in feb 2013  I set it arbitrarily to 1 since I
c       can not figure out what it should be.  al

C	ICUTD=0             !!!!!!!!!!!!!!!!!!!!!!!!!
        ICUTD=(icut2 - icuta)/3             !al feb 2013

       !write(6,*) 'in DEFO_REFINE calling XFIT_ZH,icut:', icut1,icut2
        CALL XFIT_ZH(NUMBER,N2,N,PLOT,ICUT1,ICUT2)

c       write(6,*) 'in DEFO_REFINE calling XFIT_LOW1,icut:',
c     &             icut1,icut2,icuta,icutd

        CALL XFIT_LOW1(NUMBER,N2,PLOT,ICUTA,ICUT2,ICUTD,N)

        ICUT1=ICUTB
        ICNSTRNT=ICUTB
        IBACK=1

        CALL XFIT_LOW2(NUMBER,N2,PLOT,ICUTB,ICUT2,N)
        CALL XFIT_LOW10(NUMBER,N2,PLOT,IBACK,ICUTB)
        CALL XFIT_EN(NUMBER,N2,PLOT,
     &           PS,CS,LAMBDA,CONTRAST,AV_DEFO,IPT1,IPT2,IPT3,IPT4)

        IF(ICUTB.GT.IPT2) THEN
            CALL XFIT_LOW6(NUMBER,N2,PLOT,ICUTB,IPT1,IPT2,IPT3)
        ENDIF

        IF(IPT2.NE.0) THEN
            CALL XFIT_LOW7(NUMBER,N2,PLOT,IBACK,IPT2)
            CALL XFIT_LOW7(NUMBER,N2,PLOT,IPT2,ICUT2)
        ENDIF

        CALL XFIT_BKGND3(NUMBER,N2,N,PLOT,ICUT1,ICUT2,ICNSTRNT)

        IF(N_TEST.GT.4) THEN
            CALL XFIT_LOW8(NUMBER,N2,PLOT,IPT2,ICUT2)
        ENDIF

        IF(IPT4.NE.0) THEN
            IFIT_B_STRT=IPT4
            IPT_END=NUMBER
             CALL XFIT_HIGH_EN(NUMBER,N2,PLOT,IPT4,IPT_END)
        ELSE
            IFIT_B_STRT=IPT3
        ENDIF

            IFIT_B_STOP=ICUT2
            TMP_OUT(1:NUMBER,:)=PLOT(:,:)
            TMP_PLOT(:,:)=PLOT(:,:)
            ICUT1=ICUTB
        IF(ICUT1.GT.IPT1) THEN
            ICUT1=IPT1
            ICUT2=NUMBER
        ENDIF

        IF(ICUTA.LT.ICUT2/2.0) THEN
            ICUT1=ICUTA
        ENDIF
CC==TO SMALL PIXEL SIZE MICORGRAPH
        IF(PS.LT.2.) THEN
            ICUT1=ICUTB
            TMP_XX=NUMBER/2.
            ICUT2=INT(TMP_XX)
        ENDIF

        CALL DEFOCUS_EST1(NUMBER,N2,N,PLOT,ICUT1,ICUT2,
     &           PS,CS,LAMBDA,CONTRAST,DEFOCUS)

        IF(DEFOCUS.GT.DEFO_MAX.OR.DEFOCUS.LT.AV_DEFO*2./3.) THEN
            PLOT(:,:)=TMP_PLOT(:,:)
            ICUT1=ICUTB
            CALL DEFOCUS_EST1(NUMBER,N2,N,PLOT,ICUT1,ICUT2,
     &           PS,CS,LAMBDA,CONTRAST,DEFOCUS)

          IF(DEFOCUS.GT.DEFO_MAX.OR.DEFOCUS.LT.AV_DEFO*2./3.)THEN
              XX_CC=9999
              GO TO 1122
          ENDIF
        ENDIF

        IF(DEFOCUS.GT.DEFO_MIN) THEN
            AV_DEFO=DEFOCUS
        ELSE
            CALL DEFOCUS_GUESSING4(TMP,NMB,NT,NUMBER,N2,TMP_OUT,
     &                 IPW_ZERO,CS,LAMBDA,CONTRAST,PS,AV_DEFO,
     &                 DEFO_MIN,DEFO_MAX ,XX_CC)

        ENDIF

        N0=0
        CALL RESOLUTION_ANALYZE(NUMBER,N2,N0,PS, CS,
     &           LAMBDA, AV_DEFO, CONTRAST,TMP_PLOT,XX_CC)
1122    DEALLOCATE(TMP_PLOT)
        END

CC======================================================================================

        SUBROUTINE XFIT_LOW10(K,N2,PLOT,ICUT1,ICUT2)
CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA

CC==BEGIN: FIRST ALLOCATE ARRAYS FOR CALCULATION
        N=6
        L=0
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG
CC==
        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
         RETURN
         ENDIF
CC==
        DO I=1,N-1
            Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
CC==
        CALL PARTI26 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,ICUT1,ICUT2,KP,N)
        DEALLOCATE(Q2)
        END
CC==================================================================
        SUBROUTINE PARTI26 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT1,ICUT1,ICUT2,KP,N)
        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)

        INTEGER  TMPI,TMPJ

        ISTART=ICUT1
        ISTOP=ICUT2

        N_3=N
        KS=ISTOP-ISTART+1
        L=0
        ISWI=1
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)

        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-1)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
        END

CC===========================================================================

        SUBROUTINE DEFOCUSGUESSING3(POW2,
     &           NSAM1,NROW1,NMB,PS,CS,LAMBDA,CONTRAST,AV_DEFO,XX_CC)
        DIMENSION POW2(NSAM1,NROW1)
        REAL, DIMENSION(:,:), ALLOCATABLE :: PLOT,TMP,TMP_X
        REAL, DIMENSION(:), ALLOCATABLE :: X_DEFO
        REAL LAMBDA
        REAL QUADPI
        PARAMETER (QUADPI = 3.1415926535897932384626)

        NMRANK=6

        ALLOCATE(TMP(NMB,2),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
        ENDIF

        ALLOCATE(X_DEFO(NMRANK),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
        ENDIF

CC=======ICUT_LOW_FRQ
CC==N_ZERO IS MAXIMUM NZEROS OF CTF CONTROL POINTS FOR ASTIGMATISM ESTIMATION

        XSTART = 0.0
        XSTEP  = 180.0
        CALL ZH_CRCSE2(POW2,TMP(1,2),NSAM1,NROW1,NMB,XSTART,XSTEP)
        IPW_ZERO=0

        DO I=1,NMB
           TMP(I,1)=I
        ENDDO
CC==
        IF(NMB.LT.5) THEN
            CALL ERRT(36,'TFED',NE)
            RETURN
        ENDIF

        DO I=1,NMB
            IF(TMP(I,2).EQ.0.0) THEN
            IPW_ZERO=I
            ENDIF
        ENDDO
CC==ESTIMATE THE MINIMUM DEFOCUS BASED ON GIVEN VOLTAGE, CS, AND PIXEL SIZE
        TMP1_DEFO=(2.*PS)**2/LAMBDA
        TMP2_DEFO=CONTRAST/(1.0-CONTRAST)
        TMP2_DEFO=ATAN(TMP2_DEFO)*TMP1_DEFO/QUADPI
        TMP3_DEFO=.5*CS*LAMBDA**2/(2.*PS)**2
        DEFO_MIN=TMP2_DEFO+TMP3_DEFO

        TMP1_DEFO=.5*REAL(NMB-1)/PS/REAL(NMB)
        TMP1_DEFO=TMP1_DEFO*TMP1_DEFO
        TMP2_DEFO=.5/PS
        TMP2_DEFO=TMP2_DEFO*TMP2_DEFO
        DEFO_MAX=1.0/(TMP2_DEFO-TMP1_DEFO)/LAMBDA
        DEFO_MAX=DEFO_MAX+.5*CS*LAMBDA**2*(TMP2_DEFO+TMP1_DEFO)
CC==WE FIT A STRAIGHT LINE BELOW PW TO LOCATE POINT CUT OFF LOW FREQUENCY REGION
        N2=4
        NUMBER=NMB-IPW_ZERO
        NDGREE=0
        IBACK=1
CCC=========THRESHOLD OF POWER
        X_TRH1=0.45
        X_TRH2=0.90

        ALLOCATE(PLOT(NUMBER,N2),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
        ENDIF
         ALLOCATE(TMP_X(NUMBER,N2),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
         ENDIF
CC==Using power of PSP to cut off low frequencies region (threshold 45%)
        PLOT(:,3)=TMP(IPW_ZERO+1:NMB,2)
        XX_SUM=SUM(PLOT(:,3))
        XXX_SUM=0.0

        DO II=1,NUMBER
          XXX_SUM=XXX_SUM+PLOT(II,3)/XX_SUM
            IF(XXX_SUM.GT.X_TRH1) THEN
          ICUT1=II
          GO TO 1101
            ENDIF
        ENDDO

1101    XXX_SUM=0.0
CC==Using power of PSP to cut off high frequencies region (threshold 90%)
        DO II=1,NUMBER
           XXX_SUM=XXX_SUM+PLOT(II,3)/XX_SUM

           IF(XXX_SUM.GT.X_TRH2) THEN
              ICUT2=II
              GO TO 1100
           ENDIF
        ENDDO

1100    AV_DEFO=DEFO_MIN
        DO III=2,5
            PLOT(:,1)=TMP(IPW_ZERO+1:NMB,1)-1.0
            PLOT(:,3)=ALOG(TMP(IPW_ZERO+1:NMB,2))
            IFIT1=ICUTB
            IFIT2=NUMBER
            N=III
            CALL XFIT_ZH1(NUMBER,N2,N,PLOT,ICUT1,ICUT2)

            TMP_X(:,:)=PLOT(:,:)

            CALL DEFOCUS_EST1(NUMBER,N2,N,PLOT,ICUT1,ICUT2,
     &           PS,CS,LAMBDA,CONTRAST,DEFOX)
            X_DEFO(III-1)=DEFOX

            N0=0
            CALL RESOLUTION_ANALYZE(NUMBER,N2,N0,PS, CS,
     &           LAMBDA, AV_DEFO, CONTRAST,TMP_X,X_CC)
            IF(AV_DEFO.LT.DEFOX.AND.DEFOX.LT.DEFO_MAX) THEN
               AV_DEFO=DEFOX
               XX_CC=X_CC
            ENDIF

        ENDDO
        DEALLOCATE(PLOT,TMP,X_DEFO,TMP_X)
        END

cc================================================================

        SUBROUTINE XFIT_ZH1(K,N2,N,PLOT,ICUT1,ICUT2)

CCC==INEQUALITY CONSTRAINED LINEAR SQUARE MINIMIZATION
CCC==Q1, Q2=KLM2D*N2D
CC==PLOT=K*4
CC==ENVELOPE FITTING  NO EQUALITY CONTRAINTS 

        REAL, DIMENSION(:,:), ALLOCATABLE :: Q2
CC==INCREASE ARRAY Q4 FOR DIFFERENT RANK FITTING TO ABV OR UDR 
        REAL PLOT(K,N2)
        REAL LAMBDA

        L=0
        M=K
        N1=N+1
CC==DEFINE SIZE OF WORKING ARRAYS!
        KLMD=K+L+M
        KLM2D=K+L+M+2
        NKLMD=K+L+M+N
        N2D=N+2
        KP=K
CC==FOR N.GT.2 CACULATION, ARRAY SIZE NEED ENLARGED!
        N_LARG=KLMD*2
        KLM2D=KLM2D+N_LARG
        KLMD=KLMD+N_LARG
        N2D=N2D+N_LARG
        NKLMD=NKLMD+N_LARG

        ALLOCATE(Q2(KLM2D,N2D),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
         ENDIF

        DO I=1,N-1
            Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
        Q2(K+1:2*K,1:N+1)=Q2(1:K,1:N+1)

        CALL PARTI28 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &           K,N2,PLOT,N,ICUT1,ICUT2,KP)

        DO I=1,N-1
            Q2(1:K,I)=PLOT(1:K,1)**I
        ENDDO

        Q2(1:K,N)=1.0
        Q2(1:K,N+1)=PLOT(1:K,3)
        Q2(K+1:2*K,1:N+1)=Q2(1:K,1:N+1)

        CALL PARTI29 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT,N,ICUT1,ICUT2,KP)
        DEALLOCATE(Q2)
        END

CC=============================================================

        SUBROUTINE PARTI28 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &           K,N2,PLOT1,N,ICUT1,ICUT2,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ
CC==LOW FREQUENCIES REGION RE-FITTING!!!
CC==INTERMDEDIATE REGION FITTING
        ISTART=ICUT1
        ISTOP=ICUT2
        N_3=N
        KS=ISTOP-ISTART+1
        L=0
        ISWI=2
CC==TT1_02 is destroyed after calling lsfit

        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-1)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)

        KS=ISTART
        TT1_Q2(1:KS,1:N+1)=Q2(1:ISTART,1:N+1)

        CALL LSFIT(KS,L,N2,PLOT1(1,2+(ISWI-1)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)

        KS=K-ISTOP+1
        TT1_Q2(1:KS,1:N+1)=Q2(ISTOP:K,1:N+1)

        CALL LSFIT(KS,L,N2,PLOT1(ICUT2,2+(ISWI-1)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
        END

CC=============================================================

        SUBROUTINE PARTI29 (KLMD,NKLMD,KLM2D,N2D,Q2,
     &  K,N2,PLOT1,N,ICUT1,ICUT2,KP)

        REAL Q2(KLM2D,N2D),PLOT1(K,N2)
CC==TT1_Q2,Q,X,RES,CU,S,IU are automatic arrays
        REAL TT1_Q2(KLM2D,N2D)
        DOUBLE PRECISION Q(KLM2D,N2D),X(N2D),RES(KLMD),CU(2,NKLMD),
     &                   S(KLMD)
        INTEGER IU(2,NKLMD)
        INTEGER  TMPI,TMPJ

        ISTART=ICUT1
        ISTOP=ICUT2
        N_3=N
        ISWI=1
        L=0

        KS=ISTOP-ISTART+1
        TT1_Q2(1:KS,1:N+1)=Q2(ISTART:ISTOP,1:N+1)
        CALL LSFIT(KS,L,N2,PLOT1(ISTART,2+(ISWI-1)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)


        KS=ISTART
        TT1_Q2(1:KS,1:N+1)=Q2(1:ISTART,1:N+1)
        CALL LSFIT(KS,L,N2,PLOT1(1,2+(ISWI-1)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)

        KS=K-ISTOP+1
        TT1_Q2(1:KS,1:N+1)=Q2(ISTOP:K,1:N+1)
        CALL LSFIT(KS,L,N2,PLOT1(ISTOP,2+(ISWI-1)*2),
     &         KLM2D,N2D,TT1_Q2,N_3,ISWI,
     &         Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)
        END

CC=============================================================

        SUBROUTINE DEFOCUS_EST1(K,N2,N,PLOT,ICUT1,ICUT2,
     &           PS,CS,LAMBDA,CONTRAST,DEFOCUS)
        REAL PLOT(K,N2)
        REAL LAMBDA

        DO KKK=1,ICUT2
            XXX=PLOT(KKK,1)/2./PS/REAL(K)
            XXXXX=EXP(PLOT(KKK,4))-EXP(PLOT(KKK,2))

            IF(PLOT(KKK,2).NE.0.0) THEN
          PLOT(KKK,2)=EXP(PLOT(KKK,2))
          PLOT(KKK,3)=EXP(PLOT(KKK,3))
          PLOT(KKK,4)=EXP(PLOT(KKK,4))
            ELSE
          PLOT(KKK,2)=0.0
          PLOT(KKK,3)=0.0
          PLOT(KKK,4)=0.0
            ENDIF

        ENDDO

        PLOT(1:ICUT2,3)=PLOT(1:ICUT2,3)-PLOT(1:ICUT2,2)
        PLOT(1:ICUT2,4)=PLOT(1:ICUT2,4)-PLOT(1:ICUT2,2)

        ISWI=1
        DZMAX=GETDEFOCUS(PLOT(1,3),PLOT(1,1),K,PS,CS,LAMBDA,
     &           CONTRAST,XSCORE,PLOT(1,4),ICUT1,ICUT2,ISWI)
CCCC=======OUTPUT RESULT======================================
         DEFOCUS=DZMAX
 
        END

CC================================================================

        SUBROUTINE POWER_SUM (TMP_X,N_XX,X_ITGL)

        REAL TMP_X(N_XX)
        X_ITGL=SUM(TMP_X(1:N_XX))
        END

CC====================================================================

        SUBROUTINE DEFOCUS_GUESSING4(TMP,NMB,NT,NUMBER,N2,TMP_OUT,
     &              IPW_ZERO,CS,LAMBDA,CONTRAST,PS,AV_DEFO,
     &              DEFO_MIN,DEFO_MAX ,XX_CC)

CC==THIS SUBROUTINE CORRECT DEFOCUS MIS_ESTIMATED DUE TO LONG AND FLAT HIGH FREQUENCY TAIL
        REAL LAMBDA
        DIMENSION TMP(NMB,2),TMP_OUT(NMB,NT)
        REAL, DIMENSION(:,:), ALLOCATABLE :: PLOT,TMP_PLOT

        X_TRH1=0.45
        X_TRH2=0.90

        ALLOCATE(PLOT(NUMBER,N2),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
         ENDIF
        ALLOCATE(TMP_PLOT(NUMBER,N2),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED',NE)
            RETURN
         ENDIF
        PLOT(:,3)=TMP(IPW_ZERO+1:NMB,2)
        XX_SUM=SUM(PLOT(:,3))
        XXX_SUM=0.0

        DO II=1,NUMBER
            XXX_SUM=XXX_SUM+PLOT(II,3)/XX_SUM
            IF(XXX_SUM.GT.X_TRH1) THEN
          ICUT1=II
          GO TO 1101
            ENDIF

        ENDDO

1101    XXX_SUM=0.0
CC==Using power of PSP to cut off high frequencies region (threshold 90%)
        DO II=1,NUMBER
            XXX_SUM=XXX_SUM+PLOT(II,3)/XX_SUM
            IF(XXX_SUM.GT.X_TRH2) THEN
              ICUT2=II
              GO TO 1100
        ENDIF

        ENDDO
 
1100    PLOT(:,1)=TMP(IPW_ZERO+1:NMB,1)-1.0
        PLOT(:,3)=ALOG(TMP(IPW_ZERO+1:NMB,2))
        N=NDGREE+1
CC==CC==FIT BACKGROUND FROM ICUTB TO ICUT2

        CALL FIT_ALL(NUMBER,N2,PLOT)

        TMP_OUT(1:NUMBER,:)=PLOT(:,:)
        TMP_PLOT(:,:)=PLOT(:,:)

        CALL DEFOCUS_EST(NUMBER,N2,N,PLOT,ICUT1,ICUT2,
     &           PS,CS,LAMBDA,CONTRAST,DEFOCUS,NZERO,
     &           IFIRST,ISECND)

        IF(DEFOCUS.GT.DEFO_MAX) THEN
            XX_CC=9999
            GO TO 1122
        ELSE
            AV_DEFO=DEFOCUS
        ENDIF
        N0=0
        CALL RESOLUTION_ANALYZE(NUMBER,N2,N0,PS, CS,
     &           LAMBDA, AV_DEFO, CONTRAST,TMP_PLOT,XX_CC)
1122    DEALLOCATE(PLOT,TMP_PLOT)
        END

CC==================================================================================

        SUBROUTINE AST_CALC1(NNLOOP,BETA,XDEFO,AV_DEFO,
     &           AST_AGL,TMP_AMP,TMP_SUM,RATIO)
        REAL XDEFO(NNLOOP)
        PARAMETER (QUADPI = 3.1415926535897932384626)
        TMP_SIN=0.0
        TMP_COS=0.0
        TMP1=BETA/180.0*QUADPI
        TMP_SUM=SUM(XDEFO)/REAL(NNLOOP)

        DO I=1,NNLOOP
            FREQ=2.0*(I-1+.5)*TMP1
            TMP_SIN=TMP_SIN+SIN(FREQ)*(XDEFO(I)-TMP_SUM)
            TMP_COS=TMP_COS+COS(FREQ)*(XDEFO(I)-TMP_SUM)
        ENDDO

        TMP_SIN=2.0*TMP_SIN/REAL(NNLOOP)
        TMP_COS=2.0*TMP_COS/REAL(NNLOOP)
        TMP_AMP=SQRT(TMP_SIN**2+TMP_COS**2)*2
        AST_AGL=ATAN2(TMP_COS,TMP_SIN)

        IF(AST_AGL.LT.0.0) AST_AGL=2.0*QUADPI+AST_AGL

        AST_AGL=AST_AGL*180.0/QUADPI/2.0-90.0
        TMP_SUM=0.0

        DO I=1,NNLOOP
            TMP_SUM=TMP_SUM+XDEFO(I)-
     &      TMP_AMP/2.0*SIN(2*(REAL(I+.5)*BETA-AST_AGL))
        ENDDO

        TMP_SUM=TMP_SUM/REAL(NNLOOP)
        RATIO=TMP_AMP/AV_DEFO
        END

CC===========================================================================================

        SUBROUTINE AST_CALC(POW2,NSAM1,NROW1,NMB,NT,PS,CS,
     &           LAMBDA,CONTRAST,AV_DEFO,NDGREE,AST_AGL,TMP_AMP,TMP_SUM)

        REAL POW2(NSAM1,NROW1)
        REAL LAMBDA
        REAL, DIMENSION(:,:), ALLOCATABLE :: TMP_OUT
        REAL, DIMENSION(:),   ALLOCATABLE :: XDEFO

        BETA   = 18.
        X_LOOP = 180. / BETA

1122    NNLOOP = INT(X_LOOP)     ! number of angular determinations

        DEFO_CUT1 = .5  * AV_DEFO
        DEFO_CUT2 = 2.0 * AV_DEFO

        ALLOCATE(XDEFO(NNLOOP),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED; XDEFO',NNLOOP)
            RETURN
        ENDIF

        ALLOCATE(TMP_OUT(NMB,NT),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'TFED; TMP_OUT',NMB*NT)
            RETURN
        ENDIF

        DO IILOOP=1,NNLOOP
            TMP1       = REAL(IILOOP-1)
            TMP2       = REAL(IILOOP)
            BETA0      = TMP1 * BETA
            BETA1      = TMP2 * BETA
            NDGREE_AST = NDGREE

            !write(6,*) 'ast_calc; call GUESSING1 number:',beta1

            CALL DEFOCUSGUESSING1(POW2,NSAM1,NROW1,NMB,BETA0,BETA1,
     &           PS,CS,LAMBDA,CONTRAST,DEFOCUS,
     &           ICUT_LOW_FRQ,ICN_SND,NDGREE,TMP_OUT,NT,AST_X_CC,NUMBER)

            !write(6,*) 'after ast_calc; :',number,iiloop

CC==        IN CASE OF VERY STRONG ASTIGMATISM; ADJUST POLYNOMIAL DEGREE
            IF (BETA .LT. 18.) THEN
C               this can never be used!!!!! al feb 2013

          TMP1       = DEFOCUS / FST_DEFO * REAL(NDGREE)
          NDGREE_AST = INT(TMP1)

              IF (NDGREE_AST .LT. 2) NDGREE_AST = 2

          CALL DEFOCUSGUESSING2(POW2,NSAM1,NROW1,NMB,BETA0,BETA1
     &                 ,PS,CS,LAMBDA,CONTRAST,DEFOCUS,IILOOP,
     &                 ICUT_LOW_FRQ,ICN_SND,NDGREE_AST)

            ENDIF

            write(6,'(A,f6.2,f10.1,/)') 'ast_calc :',beta0,defocus

            XDEFO(IILOOP) = DEFOCUS
        ENDDO
CC==    CHECK

        DO I=1,NNLOOP
           IF(XDEFO(I) .LT. DEFO_CUT1) THEN
              XDEFO(I) = AV_DEFO
           ELSEIF(XDEFO(I).GT.DEFO_CUT2) THEN
              XDEFO(I) = AV_DEFO
           ENDIF
        ENDDO

        CALL AST_CALC1(NNLOOP,BETA,XDEFO,
     &                 AV_DEFO,AST_AGL,TMP_AMP,TMP_SUM,RATIO)

CC==    FOR STRONG ATIGMATISM ONLY

        IF (RATIO.GT.2 .AND. BETA.GE.18) THEN
           BETA=3.
           GO TO 1122
        ENDIF

        DEALLOCATE(XDEFO,TMP_OUT)
        END

CC===========================================================================

c**************************************************************************
C *  LSFIT.F
C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2002, Z. Huang & P. A. Penczek                        *
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
C **********************************************************************

       SUBROUTINE LSFIT(KS,L,N2,P,KLM2D,N2D,Q1,N,ISWI,
     &           Q,X,RES,CU,S,IU,KLMD,NKLMD,KP)

C  DRIVER PROGRAM FOR SUBROUTINE CL1.
C                                                                       
C THIS PROGRAM SOLVES A K BY N OVERDETERMINED SYSTEM                    
C                                                                       
C    AX=B                                                               
C                                                                       
C IN THE L1 SENSE SUBJECT TO L EQUALITY CONSTRAINTS                     
C                                                                       
C    CX=D                                                               
C                                                                       
C AND M INEQUALITY CONSTRAINTS                                          
C                                                                       
C    EX.LE.F.                                                           
C                                                                       
C COMPLETE DETAILS OF THE PARAMETERS MAY BE                             
C FOUND IN THE DOCUMENTATION OF THE SUBROUTINE.                         
C                                                                       
C THE ARRAYS ARE CURRENTLY DIMENSIONED TO ALLOW PROBLEMS                
C FOR WHICH K+L+M .LE. 100, N .LE. 10.                                  
C                                                                       
C THE PROGRAM MAY BE TESTED ON THE FOLLOWING DATA.  
                                                                                 
        DOUBLE PRECISION Q(KLM2D,N2D),TOLER,X(N2D),RES(KLMD),
     &           CU(2,NKLMD),S(KLMD),TMP
        DIMENSION IU(2,NKLMD)
        REAL Q1(KLM2D,N2D),P(KP)
C modified by Zhong Huang, July,12,02                                   
cc==TOLER, ITER, KODE are set in the program;TOLER and ITER are set according to tests
        L=0
        KODE=0
        TOLER=.000001
        ITER=500
CC==ZHONG HUANG,JULY,12,02;L=0,1,2,3,4,5,6 correspond to different equality constraints

        M=KS
        N1=N+1
        IF(ISWI.EQ.1) THEN
           Q(1:KS,1:N1)=DBLE(Q1(1:KS,1:N1))
           Q(KS+1:2*KS,1:N1)=DBLE(Q1(1:KS,1:N1))
        ELSEIF(ISWI.EQ.2) THEN
           Q(1:KS,1:N1)=DBLE(Q1(1:KS,1:N1))
           Q(KS+1:2*KS,1:N1)=-DBLE(Q1(1:KS,1:N1))
        ELSEIF(ISWI.EQ.3) THEN
           L=2
           Q(1:KS+2,1:N1)=DBLE(Q1(1:KS+2,1:N1))
           Q(KS+3:2+2*KS,1:N1)=DBLE(Q1(1:KS,1:N1))
        ELSEIF(ISWI.EQ.4) THEN
           L=2
           Q(1:KS+2,1:N1)=DBLE(Q1(1:KS+2,1:N1))
           Q(KS+3:2+2*KS,1:N1)=-DBLE(Q1(1:KS,1:N1))
        ELSEIF(ISWI.EQ.5) THEN
           L=1
           Q(1:KS+1,1:N1)=DBLE(Q1(1:KS+1,1:N1))
           Q(KS+2:2*KS+1,1:N1)=-DBLE(Q1(1:KS,1:N1))
        ELSEIF(ISWI.EQ.6) THEN
           L=1
           Q(1:KS+1,1:N1)=DBLE(Q1(1:KS+1,1:N1))
           Q(KS+2:2*KS+1,1:N1)=DBLE(Q1(1:KS,1:N1))
        ELSEIF(ISWI.EQ.7) THEN
           L=3
           Q(1:KS+3,1:N1)=DBLE(Q1(1:KS+3,1:N1))
           Q(KS+4:2*KS+3,1:N1)=-DBLE(Q1(1:KS,1:N1))
        ELSEIF(ISWI.EQ.8) THEN
           L=4
           Q(1:KS+4,1:N1)=DBLE(Q1(1:KS+4,1:N1))
           Q(KS+5:2*KS+4,1:N1)=-DBLE(Q1(1:KS,1:N1))
        ENDIF
        CALL CL1(KS, L, M, N, KLMD, KLM2D, NKLMD, N2D, Q,
     * KODE, TOLER, ITER, X, RES, ERROR, CU, IU, S) 
                        
CC==CALCULATING THE RESTRAINED RESULTS
        IF (KODE .GT. 0) THEN
            !write(6,*) 'Bad kode:',kode

           RETURN
        ENDIF
        DO I=1,KS
           TMP=0.0

           DO J=1,N-1
              TMP=TMP+Q1(I,1)**J*X(J)
           ENDDO

           TMP=TMP+X(N)
           P(I)=SNGL(TMP)
        ENDDO

        END
                                                         
cc=========================================================================

C **********************************************************************
C *  CL1.F
C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2002, Z. Huang & P. A. Penczek                        *
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
C **********************************************************************

CC==Inequality&Equality Constrained Least Square Fitting Program

      SUBROUTINE CL1(K, L, M, N, KLMD, KLM2D, NKLMD, N2D,      
     * Q, KODE, TOLER, ITER, X, RES, ERROR, CU, IU, S)

C THIS SUBROUTINE USES A MODIFICATION OF THE SIMPLEX
C METHOD OF LINEAR PROGRAMMING TO CALCULATE AN L1 SOLUTION
C TO A K BY N SYSTEM OF LINEAR EQUATIONS
C             AX=B
C SUBJECT TO L LINEAR EQUALITY CONSTRAINTS
C             CX=D
C AND M LINEAR INEQUALITY CONSTRAINTS
C             EX.LE.F.
C DESCRIPTION OF PARAMETERS
C K      NUMBER OF ROWS OF THE MATRIX A (K.GE.1).
C L      NUMBER OF ROWS OF THE MATRIX C (L.GE.0).
C M      NUMBER OF ROWS OF THE MATRIX E (M.GE.0).
C N      NUMBER OF COLUMNS OF THE MATRICES A,C,E (N.GE.1).
C KLMD   SET TO AT LEAST K+L+M FOR ADJUSTABLE DIMENSIONS.
C KLM2D  SET TO AT LEAST K+L+M+2 FOR ADJUSTABLE DIMENSIONS.
C NKLMD  SET TO AT LEAST N+K+L+M FOR ADJUSTABLE DIMENSIONS.
C N2D    SET TO AT LEAST N+2 FOR ADJUSTABLE DIMENSIONS
C Q      TWO DIMENSIONAL REAL ARRAY WITH KLM2D ROWS AND
C        AT LEAST N2D COLUMNS.
C        ON ENTRY THE MATRICES A,C AND E, AND THE VECTORS
C        B,D AND F MUST BE STORED IN THE FIRST K+L+M ROWS
C        AND N+1 COLUMNS OF Q AS FOLLOWS
C             A B
C         Q = C D
C             E F
C        THESE VALUES ARE DESTROYED BY THE SUBROUTINE.
C KODE   A CODE USED ON ENTRY TO, AND EXIT
C        FROM, THE SUBROUTINE.
C        ON ENTRY, THIS SHOULD NORMALLY BE SET TO 0.
C        HOWEVER, IF CERTAIN NONNEGATIVITY CONSTRAINTS
C        ARE TO BE INCLUDED IMPLICITLY, RATHER THAN
C        EXPLICITLY IN THE CONSTRAINTS EX.LE.F, THEN KODE
C        SHOULD BE SET TO 1, AND THE NONNEGATIVITY
C        CONSTRAINTS INCLUDED IN THE ARRAYS X AND
C        RES (SEE BELOW).
C        ON EXIT, KODE HAS ONE OF THE
C        FOLLOWING VALUES
C             0- OPTIMAL SOLUTION FOUND,
C             1- NO FEASIBLE SOLUTION TO THE
C                CONSTRAINTS,
C             2- CALCULATIONS TERMINATED
C                PREMATURELY DUE TO ROUNDING ERRORS,
C             3- MAXIMUM NUMBER OF ITERATIONS REACHED.
C TOLER  A SMALL POSITIVE TOLERANCE. EMPIRICAL
C        EVIDENCE SUGGESTS TOLER = 10**(-D*2/3),
C        WHERE D REPRESENTS THE NUMBER OF DECIMAL
C        DIGITS OF ACCURACY AVAILABLE. ESSENTIALLY,
C        THE SUBROUTINE CANNOT DISTINGUISH BETWEEN ZERO
C        AND ANY QUANTITY WHOSE MAGNITUDE DOES NOT EXCEED
C        TOLER. IN PARTICULAR, IT WILL NOT PIVOT ON ANY
C        NUMBER WHOSE MAGNITUDE DOES NOT EXCEED TOLER.
C ITER   ON ENTRY ITER MUST CONTAIN AN UPPER BOUND ON
C        THE MAXIMUM NUMBER OF ITERATIONS ALLOWED.
C        A SUGGESTED VALUE IS 10*(K+L+M). ON EXIT ITER
C        GIVES THE NUMBER OF SIMPLEX ITERATIONS.
C X      ONE DIMENSIONAL REAL ARRAY OF SIZE AT LEAST N2D.
C        ON EXIT THIS ARRAY CONTAINS A
C        SOLUTION TO THE L1 PROBLEM. IF KODE=1
C        ON ENTRY, THIS ARRAY IS ALSO USED TO INCLUDE
C        SIMPLE NONNEGATIVITY CONSTRAINTS ON THE
C        VARIABLES. THE VALUES -1, 0, OR 1
C        FOR X(J) INDICATE THAT THE J-TH VARIABLE
C        IS RESTRICTED TO BE .LE.0, UNRESTRICTED,
C        OR .GE.0 RESPECTIVELY.
C RES    ONE DIMENSIONAL REAL ARRAY OF SIZE AT LEAST KLMD.
C        ON EXIT THIS CONTAINS THE RESIDUALS B-AX
C        IN THE FIRST K COMPONENTS, D-CX IN THE
C        NEXT L COMPONENTS (THESE WILL BE =0),AND
C        F-EX IN THE NEXT M COMPONENTS. IF KODE=1 ON
C        ENTRY, THIS ARRAY IS ALSO USED TO INCLUDE SIMPLE
C        NONNEGATIVITY CONSTRAINTS ON THE RESIDUALS
C        B-AX. THE VALUES -1, 0, OR 1 FOR RES(I)
C        INDICATE THAT THE I-TH RESIDUAL (1.LE.I.LE.K) IS
C        RESTRICTED TO BE .LE.0, UNRESTRICTED, OR .GE.0
C        RESPECTIVELY.
C ERROR  ON EXIT, THIS GIVES THE MINIMUM SUM OF
C        ABSOLUTE VALUES OF THE RESIDUALS.
C CU     A TWO DIMENSIONAL REAL ARRAY WITH TWO ROWS AND
C        AT LEAST NKLMD COLUMNS USED FOR WORKSPACE.
C IU     A TWO DIMENSIONAL INTEGER ARRAY WITH TWO ROWS AND
C        AT LEAST NKLMD COLUMNS USED FOR WORKSPACE.
C S      INTEGER ARRAY OF SIZE AT LEAST KLMD, USED FOR
C        WORKSPACE.
C IF YOUR FORTRAN COMPILER PERMITS A SINGLE COLUMN OF A TWO
C DIMENSIONAL ARRAY TO BE PASSED TO A ONE DIMENSIONAL ARRAY
C THROUGH A SUBROUTINE CALL, CONSIDERABLE SAVINGS IN
C EXECUTION TIME MAY BE ACHIEVED THROUGH THE USE OF THE
C FOLLOWING SUBROUTINE, WHICH OPERATES ON COLUMN VECTORS.
C     SUBROUTINE COL(V1, V2, XMLT, NOTROW, K)
C THIS SUBROUTINE ADDS TO THE VECTOR V1 A MULTIPLE OF THE
C VECTOR V2 (ELEMENTS 1 THROUGH K EXCLUDING NOTROW).
C     DIMENSION V1(K), V2(K)
C     KEND = NOTROW - 1
C     KSTART = NOTROW + 1
C     IF (KEND .LT. 1) GO TO 20
C     DO 10 I=1,KEND
C        V1(I) = V1(I) + XMLT*V2(I)
C  10 CONTINUE
C     IF(KSTART .GT. K) GO TO 40
C  20 DO 30 I=KSTART,K
C       V1(I) = V1(I) + XMLT*V2(I)
C  30 CONTINUE
C  40 RETURN
C     END
C SEE COMMENTS FOLLOWING STATEMENT LABELLED 440 FOR
C INSTRUCTIONS ON THE IMPLEMENTATION OF THIS MODIFICATION.
      DOUBLE PRECISION XSUM
C      DOUBLE PRECISION DBLE
C      REAL Q, X, Z, CU, SN, ZU, ZV, CUV, RES, XMAX, XMIN,
CC     * ERROR, PIVOT, TOLER, TPIVOT
      DOUBLE PRECISION Q,X,RES,CU,S,Z, SN, ZU, ZV, CUV, XMAX, XMIN,
     * PIVOT, TPIVOT,TOLER
C      REAL ABS
C      INTEGER I, J, K, L, M, N, S, IA, II, IN, IU, JS, KK,
C     * NK, N1, N2, JMN, JPN, KLM, NKL, NK1, N2D, IIMN,
C     * IOUT, ITER, KLMD, KLM1, KLM2, KODE, NKLM, NKL1,
C     * KLM2D, MAXIT, NKLMD, IPHASE, KFORCE, IINEG
C      INTEGER IABS
      DIMENSION Q(KLM2D,N2D), X(N2D), RES(KLMD),
     * CU(2,NKLMD), IU(2,NKLMD), S(KLMD)
C INITIALIZATION.
      MAXIT = 500
      N1 = N + 1
      N2 = N + 2
      NK = N + K
      NK1 = NK + 1
      NKL = NK + L
      NKL1 = NKL + 1
      KLM = K + L + M
      KLM1 = KLM + 1
      KLM2 = KLM + 2
      NKLM = N + KLM
      KFORCE = 1
      ITER = 0
      JS = 1
      IA = 0
C SET UP LABELS IN Q.
      DO 10 J=1,N
         Q(KLM2,J) = J
   10 CONTINUE
      DO 30 I=1,KLM
         Q(I,N2) = N + I
         IF (Q(I,N1).GE.0.) GO TO 30
         DO 20 J=1,N2
            Q(I,J) = -Q(I,J)
   20    CONTINUE
   30 CONTINUE
C SET UP PHASE 1 COSTS.
      IPHASE = 2
      DO 40 J=1,NKLM
         CU(1,J) = 0.
         CU(2,J) = 0.
         IU(1,J) = 0
         IU(2,J) = 0
   40 CONTINUE
      IF (L.EQ.0) GO TO 60
      DO 50 J=NK1,NKL
         CU(1,J) = 1.
         CU(2,J) = 1.
         IU(1,J) = 1
         IU(2,J) = 1
   50 CONTINUE
      IPHASE = 1
   60 IF (M.EQ.0) GO TO 80
      DO 70 J=NKL1,NKLM
         CU(2,J) = 1.
         IU(2,J) = 1
         JMN = J - N
         IF (Q(JMN,N2).LT.0.) IPHASE = 1
   70 CONTINUE
   80 IF (KODE.EQ.0) GO TO 150
      DO 110 J=1,N
         IF (X(J)) 90, 110, 100
   90    CU(1,J) = 1.
         IU(1,J) = 1
         GO TO 110
  100    CU(2,J) = 1.
         IU(2,J) = 1
  110 CONTINUE
      DO 140 J=1,K
         JPN = J + N
         IF (RES(J)) 120, 140, 130
  120    CU(1,JPN) = 1.
         IU(1,JPN) = 1
         IF (Q(J,N2).GT.0.0) IPHASE = 1
         GO TO 140
  130    CU(2,JPN) = 1.
         IU(2,JPN) = 1
         IF (Q(J,N2).LT.0.0) IPHASE = 1
  140 CONTINUE
  150 IF (IPHASE.EQ.2) GO TO 500
C COMPUTE THE MARGINAL COSTS.
  160 DO 200 J=JS,N1
         XSUM = 0.D0
         DO 190 I=1,KLM
            II = Q(I,N2)
            IF (II.LT.0) GO TO 170
            Z = CU(1,II)
            GO TO 180
  170       IINEG = -II
            Z = CU(2,IINEG)
  180       XSUM = XSUM + DBLE(Q(I,J))*DBLE(Z)
C  180       XSUM = XSUM + Q(I,J)*Z
  190    CONTINUE
         Q(KLM1,J) = XSUM
  200 CONTINUE
      DO 230 J=JS,N
         II = Q(KLM2,J)
         IF (II.LT.0) GO TO 210
         Z = CU(1,II)
         GO TO 220
  210    IINEG = -II
         Z = CU(2,IINEG)
  220    Q(KLM1,J) = Q(KLM1,J) - Z
  230 CONTINUE
C DETERMINE THE VECTOR TO ENTER THE BASIS.
  240 XMAX = 0.
      IF (JS.GT.N) GO TO 490
      DO 280 J=JS,N
         ZU = Q(KLM1,J)
         II = Q(KLM2,J)
         IF (II.GT.0) GO TO 250
         II = -II
         ZV = ZU
         ZU = -ZU - CU(1,II) - CU(2,II)
         GO TO 260
  250    ZV = -ZU - CU(1,II) - CU(2,II)
  260    IF (KFORCE.EQ.1 .AND. II.GT.N) GO TO 280
         IF (IU(1,II).EQ.1) GO TO 270
         IF (ZU.LE.XMAX) GO TO 270
         XMAX = ZU
         IN = J
  270    IF (IU(2,II).EQ.1) GO TO 280
         IF (ZV.LE.XMAX) GO TO 280
         XMAX = ZV
         IN = J
  280 CONTINUE
      IF (XMAX.LE.TOLER) GO TO 490
      IF (Q(KLM1,IN).EQ.XMAX) GO TO 300
      DO 290 I=1,KLM2
         Q(I,IN) = -Q(I,IN)
  290 CONTINUE
      Q(KLM1,IN) = XMAX
C DETERMINE THE VECTOR TO LEAVE THE BASIS.
  300 IF (IPHASE.EQ.1 .OR. IA.EQ.0) GO TO 330
      XMAX = 0.
      DO 310 I=1,IA
         Z = ABS(Q(I,IN))
         IF (Z.LE.XMAX) GO TO 310
         XMAX = Z
         IOUT = I
  310 CONTINUE
      IF (XMAX.LE.TOLER) GO TO 330
      DO 320 J=1,N2
         Z = Q(IA,J)
         Q(IA,J) = Q(IOUT,J)
         Q(IOUT,J) = Z
  320 CONTINUE
      IOUT = IA
      IA = IA - 1
      PIVOT = Q(IOUT,IN)
      GO TO 420
  330 KK = 0
      DO 340 I=1,KLM
         Z = Q(I,IN)
         IF (Z.LE.TOLER) GO TO 340
         KK = KK + 1
         RES(KK) = Q(I,N1)/Z
         S(KK) = I
  340 CONTINUE
  350 IF (KK.GT.0) GO TO 360
      KODE = 2
      GO TO 590
  360 XMIN = RES(1)
      IOUT = S(1)
      J = 1
      IF (KK.EQ.1) GO TO 380
      DO 370 I=2,KK
         IF (RES(I).GE.XMIN) GO TO 370
         J = I
         XMIN = RES(I)
         IOUT = S(I)
  370 CONTINUE
      RES(J) = RES(KK)
      S(J) = S(KK)
  380 KK = KK - 1
      PIVOT = Q(IOUT,IN)
      II = Q(IOUT,N2)
      IF (IPHASE.EQ.1) GO TO 400
      IF (II.LT.0) GO TO 390
      IF (IU(2,II).EQ.1) GO TO 420
      GO TO 400
  390 IINEG = -II
      IF (IU(1,IINEG).EQ.1) GO TO 420
c 400 II = IABS(II)
  400 II =  ABS(II) 
      CUV = CU(1,II) + CU(2,II)
      IF (Q(KLM1,IN)-PIVOT*CUV.LE.TOLER) GO TO 420
C BYPASS INTERMEDIATE VERTICES.
      DO 410 J=JS,N1
         Z = Q(IOUT,J)
         Q(KLM1,J) = Q(KLM1,J) - Z*CUV
         Q(IOUT,J) = -Z
  410 CONTINUE
      Q(IOUT,N2) = -Q(IOUT,N2)
      GO TO 350
C GAUSS-JORDAN ELIMINATION.
  420 IF (ITER.LT.MAXIT) GO TO 430
      KODE = 3
      GO TO 590
  430 ITER = ITER + 1
      DO 440 J=JS,N1
         IF (J.NE.IN) Q(IOUT,J) = Q(IOUT,J)/PIVOT
  440 CONTINUE
C IF PERMITTED, USE SUBROUTINE COL OF THE DESCRIPTION
C SECTION AND REPLACE THE FOLLOWING SEVEN STATEMENTS DOWN
C TO AND INCLUDING STATEMENT NUMBER 460 BY..
C     DO 460 J=JS,N1
C        IF(J .EQ. IN) GO TO 460
C        Z = -Q(IOUT,J)
C        CALL COL(Q(1,J), Q(1,IN), Z, IOUT, KLM1)
C 460 CONTINUE
      DO 460 J=JS,N1
         IF (J.EQ.IN) GO TO 460
         Z = -Q(IOUT,J)
         DO 450 I=1,KLM1
            IF (I.NE.IOUT) Q(I,J) = Q(I,J) + Z*Q(I,IN)
  450    CONTINUE
  460 CONTINUE
      TPIVOT = -PIVOT
      DO 470 I=1,KLM1
         IF (I.NE.IOUT) Q(I,IN) = Q(I,IN)/TPIVOT
  470 CONTINUE
      Q(IOUT,IN) = 1./PIVOT
      Z = Q(IOUT,N2)
      Q(IOUT,N2) = Q(KLM2,IN)
      Q(KLM2,IN) = Z
      II = ABS(Z)
      IF (IU(1,II).EQ.0 .OR. IU(2,II).EQ.0) GO TO 240
      DO 480 I=1,KLM2
         Z = Q(I,IN)
         Q(I,IN) = Q(I,JS)
         Q(I,JS) = Z
  480 CONTINUE
      JS = JS + 1
      GO TO 240
C TEST FOR OPTIMALITY.
  490 IF (KFORCE.EQ.0) GO TO 580
      IF (IPHASE.EQ.1 .AND. Q(KLM1,N1).LE.TOLER) GO TO 500
      KFORCE = 0
      GO TO 240
C SET UP PHASE 2 COSTS.
  500 IPHASE = 2
      DO 510 J=1,NKLM
         CU(1,J) = 0.
         CU(2,J) = 0.
  510 CONTINUE
      DO 520 J=N1,NK
         CU(1,J) = 1.
         CU(2,J) = 1.
  520 CONTINUE
      DO 560 I=1,KLM
         II = Q(I,N2)
         IF (II.GT.0) GO TO 530
         II = -II
         IF (IU(2,II).EQ.0) GO TO 560
         CU(2,II) = 0.
         GO TO 540
  530    IF (IU(1,II).EQ.0) GO TO 560
         CU(1,II) = 0.
  540    IA = IA + 1
         DO 550 J=1,N2
            Z = Q(IA,J)
            Q(IA,J) = Q(I,J)
            Q(I,J) = Z
  550    CONTINUE
  560 CONTINUE
      GO TO 160
  570 IF (Q(KLM1,N1).LE.TOLER) GO TO 500
      KODE = 1
      GO TO 590
  580 IF (IPHASE.EQ.1) GO TO 570
C PREPARE OUTPUT.
      KODE = 0
  590 XSUM = 0.D0
      DO 600 J=1,N
         X(J) = 0.
  600 CONTINUE
      DO 610 I=1,KLM
         RES(I) = 0.
  610 CONTINUE
      DO 640 I=1,KLM
         II = Q(I,N2)
         SN = 1.
         IF (II.GT.0) GO TO 620
         II = -II
         SN = -1.
  620    IF (II.GT.N) GO TO 630
         X(II) = SN*Q(I,N1)
         GO TO 640
  630    IIMN = II - N
         RES(IIMN) = SN*Q(I,N1)
         IF (II.GE.N1 .AND. II.LE.NK) XSUM = XSUM +
     *    DBLE(Q(I,N1))
C     *    Q(I,N1)
  640 CONTINUE
      ERROR = XSUM
      RETURN
      END
