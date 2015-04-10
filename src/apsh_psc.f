C++*********************************************************************
C
C   APSH_PSC.F      RTSQ, RENAMED                 JAN 11 ARDEAN LEITH
C                   MAKE_CLOSE_LIST, GETANGAS     FEB 11 ARDEAN LEITH
C                   QU,QM ALLOCATED               FEB 11 ARDEAN LEITH
C                   IRAY1,IRAY2                   NOV 11 ARDEAN LEITH
C                   FBS_WANTED                    NOV 11 ARDEAN LEITH
C                   IRAY CALCULATION              FEB 12 ARDEAN LEITH
C                   AP_ENDS                       MAR 12 ARDEAN LEITH
C                   DENOISE                       SEP 12 ARDEAN LEITH
C                   FOURIER LOWPASS DENOISE       OCT 12 ARDEAN LEITH
C
C **********************************************************************
C=* This file is part of:                                              * 
C=* SPIDER - Modular Image Processing System.                          *
C=* Authors: J. Frank & A. Leith                                       *
C=* Copyright 1985-2012  Health Research Inc.                          *
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
C=* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
C=* General Public License for more details.                           *
C=*                                                                    *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program. If not, see <http://www.gnu.org/licenses> *       
C=*                                                                    *
C **********************************************************************
C
C  APSH_PSC 
C  
C  PURPOSE: FIND ROTATIONAL AND SHIFT PARAMETERS TO ALIGN A SERIES OF
C           REFERENCE IMAGES WITH SAMPLE IMAGES.
C           USED IF:  CIRCREF_IN_CORE .AND. NUMEXP >= NUMTH. 
C           BY DEFAULT USES COEFFICIENTS, NON-TRANSPOSED RINGS AND
C           COMPLEX VARIABLES TO SPEED UP CALCULATIONS.
C
C SOME PARAMETERS:
C       IREFLIST            LIST OF REF. IMAGE FILE NUMBERS      (SENT)
C       NUMREF              NO. OF IMAGES                        (SENT)
C       IEXPLIST            LIST OF EXP. IMAGE FILE NUMBERS      (SENT)
C       NUMEXP              NO. OF IMAGES                        (SENT)
C       NX,NY               ACTUAL (NOT-WINDOWED) IMAGE SIZE     (SENT)
C       ISHRANGEX,ISHRANGEY ALLOWED SHIFT RANGE                  (SENT)
C       ISTEP               SHIFT STEP WITHIN RANGE              (SENT)
C       NRING               NUMBER OF CIRCULAR RINGS             (SENT)
C       LCIRC               LENGTH OF CIRCULAR RING ARRAY        (SENT)
C       NUMR                CIRCULAR RING ARRAY INDICES          (SENT)
C       CIRCREF             CIRCULAR RINGS ARRAY                 (SENT)
C       CIRCREF_IN_CORE     INCORE MEMORY USED                   (SENT)
C       MODE                                                     (SENT)
C       REFANGDOC           REF. ANGLES FILE NAME                (SENT)
C       EXPANGDOC           EXP. ANGLES FILE NAME                (SENT)
C       SCRFILE                                                  (SENT)
C       FFTW_PLANS          FFTW PLAN POINTERS                   (SENT)
C       REFPAT              REF. IMAGE SERIES FILE TEMPLATE      (SENT)
C       EXPPAT              EXP. IMAGE SERIES FILE TEMPLATE      (SENT)
C       RANGE                                                    (SENT)
C       ROTFIRST            ROTATE/SHIFT EXP INPUT IMAGES        (SENT)
C       DENOISE             DENOISE EXP INPUT IMAGES             (SENT)
C       CKMIRROR            LOGICAL FLAG TO CHECK MIRRORING      (SENT)
C       CTYPE               SUB-OPERATION                        (SENT)
C       LUNDOC              DOC FILE I/O UNIT                    (SENT)
C       RAY1,RAY2           RAY RANGE                            (SENT)
C       FBS_WANTED          WANT TO USE RTSF IF ROTFIRST         (SENT)
C
C NOTES:  NON-MPI CODE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

cpgi$g opt=3

        SUBROUTINE APSH_PSC(IREFLIST,NUMREF,IEXPLIST,NUMEXP, 
     &               NX,NY,ISHRANGEX,ISHRANGEY,ISTEP,
     &               NRING,LCIRC,NUMR, CIRCREF,CIRCREF_IN_CORE,
     &               MODE, REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS,
     &               REFPAT,EXPPAT,RANGE,ROTFIRST,DENOISE,
     &               CKMIRROR,CTYPE,LUNDOC,RAY1,RAY2,FBS_WANTED)

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

	INTEGER              :: IREFLIST(NUMREF) 
	INTEGER              :: NUMREF 
	INTEGER              :: IEXPLIST(NUMEXP)
	INTEGER              :: NUMEXP,NX,NY,ISHRANGEX,ISHRANGEY
        INTEGER              :: ISTEP,NRING,LCIRC 
        INTEGER              :: NUMR(3,NRING)
	REAL                 :: CIRCREF(LCIRC,NUMREF)
	LOGICAL              :: CIRCREF_IN_CORE
        CHARACTER (LEN=1)    :: MODE
        CHARACTER (LEN=*)    :: REFANGDOC,EXPANGDOC,SCRFILE
        INTEGER *8           :: FFTW_PLANS(*)   ! STRUCTURE POINTERS
        CHARACTER (LEN=*)    :: REFPAT,EXPPAT 
        REAL                 :: RANGE
        LOGICAL              :: ROTFIRST,DENOISE,CKMIRROR
        CHARACTER (LEN=*)    :: CTYPE
        INTEGER              :: LUNDOC
        REAL                 :: RAY1,RAY2 
        LOGICAL              :: FBS_WANTED,CHNG_ORDER,SAY_RAW
  

	LOGICAL              :: MIRRORNEW,GOTREFANG,LIMITRANGE
        LOGICAL              :: WEIGHT
        LOGICAL              :: TRANS        ! FLAG FOR REFORMED RINGS
        LOGICAL              :: CPLX         ! FLAG FOR COMPLEX CROSRNG
        LOGICAL              :: USECOEF      ! FOR TESTING
        LOGICAL              :: ANGINHEADER,GOTEXPANG   
        INTEGER              :: NXT

C       AUTOMATIC ARRAYS
	REAL                 :: ANGOUT(3)

C       ALLOCATED ARRAYS
	REAL,    ALLOCATABLE :: EXPBUF(:,:,:)
	REAL,    ALLOCATABLE :: TMPBUF(:,:)
	INTEGER, ALLOCATABLE :: NPROJA(:)
	REAL,    ALLOCATABLE :: DLIST(:,:) 
	REAL,    ALLOCATABLE :: REFDIR(:,:), EXPDIR(:,:) 
	REAL,    ALLOCATABLE :: ANGREF(:,:), ANGEXP(:,:)
        REAL,    ALLOCATABLE :: COEFFS(:,:)
        INTEGER, ALLOCATABLE :: IXY(:,:)
        INTEGER, ALLOCATABLE :: NLOCS(:,:)

        INTEGER, PARAMETER   :: NLISTMAX = 15
        REAL                 :: PARLIST(NLISTMAX)

        REAL,    PARAMETER   :: QUADPI = 3.1415926535897932384626
        REAL,    PARAMETER   :: DGR_TO_RAD = (QUADPI/180)

        INTEGER, PARAMETER   :: LUNT    = 77
        INTEGER, PARAMETER   :: INANG   = 78
        INTEGER, PARAMETER   :: LUNRING = 50

        INTEGER              :: NBORDER = 0       ! # BORDER PIXELS
        INTEGER              :: NSUBPIX = 0       ! # SUBPIX PIXELS
        INTEGER              :: MYPID   = -1      ! NOT FOR MPI

        LOGICAL, PARAMETER   :: MPIBCAST = .TRUE.
 
        LOGICAL              :: erri2
 
C       SET TYPE OF OUTPUT DOC FILES WANTED
        NWANTOUT = 15

C       INITIALIZE CCROT STATISTICS COUNTERS
        ANGDIFTHR   = 0.0
        CALL AP_STAT_ADD(-1,CCROT,ANGDIF,ANGDIFTHR,CCROTLAS,
     &                   CCROTAVG,IBIGANGDIF,ANGDIFAVG,IMPROVCCROT,
     &                   CCROTIMPROV,IWORSECCROT,CCROTWORSE)

        LIMITRANGE  = (RANGE > 0.0)
        RANGECOS    = COS(RANGE * DGR_TO_RAD)  
        !write(6,*) ' range,rangecos:',range,rangecos

C       FIND NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)
        
        NRAYSC = NUMR(3,NRING) / 2     ! # OF RAYS   (FOURIER)
        NRAYS  = NUMR(3,NRING) -2      ! # OF RAYS   (REAL)

        ! CONVERT 0...360 ANGLES TO 1...NRAYS  RAY POSITIONS
        IF (RAY1< 0 .OR. RAY1 >360) THEN
            WRITE(NOUT,*) ' Ray Range:',RAY1,RAY2
            CALL ERRT(101,'INTEGER RAY LIMITS = 0...360',NDUM)
            GOTO 9999
        ELSEIF (RAY2 < 0 .OR. RAY2 >360) THEN
            WRITE(NOUT,*) ' Ray Range:',RAY1,RAY2
            CALL ERRT(101,'INTEGER RAY RANGE = 0...360',NDUM)
            GOTO 9999
        ENDIF

        !IRAY1 = ((RAY1+1) / FLOAT(359)) * NRAYS 
        !IRAY2 = ((RAY2+1) / FLOAT(359)) * NRAYS 
        !IF (RAY1 == 0)    IRAY1 = 1  
        !IF (RAY2 <= 0)    IRAY2 = NRAYS  
        !IF (IRAY2  > 359) IRAY2 = NRAYS  

        IRAY1 = 1 + RAY1 * FLOAT(NRAYS - 1) / 360
        IRAY2 = 1 + RAY2 * FLOAT(NRAYS - 1) / 360
        !IF (RAY2 <= 0)    IRAY2 = NRAYS  
        !IF (IRAY2  > 359) IRAY2 = NRAYS  

        !rangnew1 = ang_n(float(iray1),mode,nrays)
        !rangnew2 = ang_n(float(iray2),mode,nrays)
        !write(6,*) ' Nrays:',nrays,nring,nraysc,ray1,ray2
        !write(6,*) ' Ray range:',ray1,ray2,'-->',iray1,iray2
        !write(6,*) ' raynum:  ',iray1,iray2,'-->',rangnew1,rangnew2

	ALLOCATE(EXPBUF(NX,NY,NUMTH), 
     &           DLIST(5,NUMTH), 
     &           NPROJA(NUMTH),    STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           MWANT = NX*NY*NUMTH + 6*NUMTH    
           CALL ERRT(46,'APSH_PSC; EXPBUF & DLIST...',MWANT)
           GOTO 9999
        ENDIF
 
        IF (DENOISE) THEN

           N2X   = NX * 2
           N2Y   = NY * 2
           N2XLD = N2X + 2 - MOD(N2X,2)
           NXLD  = NX + 2 - MOD(NX,2) ! FFT PAD
           NXT   = NX
           IF (FBS_WANTED) NXT = NXLD ! FOR BUFFER IN CALLEE

	   ALLOCATE(TMPBUF(N2XLD,N2Y), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'APSH_PSC; TMPBUF',N2XLD*N2Y)
              GOTO 9999
           ENDIF

        ELSEIF (ROTFIRST) THEN
	   ALLOCATE(TMPBUF(NX,NY), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'APSH_PSC; TMPBUF',NX*NY)
              GOTO 9999
           ENDIF

           anginheader = (expangdoc .eq. '*') ! unfinished !!!!!!!!!!!
        ENDIF

        DLIST = 0.0        ! ZEROS WHOLE DLIST

C       CHOICES: SHC, SHCT,SHCU, SHCUC, SHCTC

        ILENC      = LNBLNKN(CTYPE) 
        CHNG_ORDER = .TRUE.    ! ORDER OF SHIFT AND ROTATE
        SAY_RAW    = .FALSE.   ! REPORT SHIFT AFTER ROTATION COMPENSATION
        IF (ILENC > 3) SAY_RAW = (INDEX(CTYPE(:ILENC),'R') > 0)
        IF (SAY_RAW) THEN
            WRITE(NOUT,*)' REPORTING RAW SHIFT VALUES IN REG: 13 & 14'
        ENDIF

        USECOEF = .TRUE. 
        WEIGHT  = .TRUE.    ! REF. IMAGES HAVE WEIGHTED FFT'S
        TRANS   = .FALSE.   ! DO NOT USE REFORMED RINGS/RAYS
        CPLX    = .TRUE.    ! USE COMPLEX CROSRNG VARIABLES
        IF (ILENC > 3) TRANS   = (CTYPE(4:4) .EQ. 'T')
        IF (ILENC > 4) CPLX    = (CTYPE(5:5) .EQ. 'C')
        IF (ILENC > 5) USECOEF = .NOT. (CTYPE(5:5) .EQ. 'N')
 
        !write(6,*) 'trans:',trans,'  cplx:',cplx,'  usecoef:',usecoef

        IF (TRANS) THEN
C          SET # OF POINTS ON EACH RAY AND RAY STARTING INDEX IN CIRC
C          THIS IS ONLY USED FOR TRANSFORMED RING ORDERS
	   ALLOCATE( NLOCS(2,NRAYSC+1), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
              MWANT = 2*NRAYSC+2  
              CALL ERRT(46,'APSH_PSC; NLOCS',MWANT)
              GOTO 9999
           ENDIF 
  
           CALL APRINGS_TRANS_LOCS(NUMR,NRING, NLOCS,NRAYSC)
        ELSE
	   ALLOCATE(NLOCS(1,1),  STAT=IRTFLG)
        ENDIF

        IF (USECOEF) THEN
	   ALLOCATE(COEFFS(6,LCIRC), IXY(2,LCIRC), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
              MWANT = 6*LCIRCR +2*LCIRCR 
              CALL ERRT(46,'APSH_PSC; COEFFS...',MWANT)
              GOTO 9999
           ENDIF 
           IXY     = -100       ! FOR HANDLING CIRC PADS
           COEFFS  = 0.0

c                     write(6,*) ' aprings  coef Threads:',numth
c          if (trans) write(6,*) ' aprings  transformed rings' 
c          if (cplx)  write(6,*) ' aprings  complx crosrngs'

C          READ REFERENCE IMAGES INTO REFERENCE RINGS ARRAY (CIRCREF) 
C          OR CREATE REFERENCE RINGS FILE FOR LATER READING
C          SAVES COEFFS FOR LATER USE 
           CALL APRINGS_NEW_COEF(IREFLIST,NUMREF,  NX,NY,
     &                       NRING,LCIRC,NUMR, NLOCS,NRAYSC,
     &                       COEFFS,IXY,
     &                       MODE,FFTW_PLANS,
     &                       REFPAT,LUNT, CIRCREF,CIRCREF_IN_CORE,
     &                       LUNRING,SCRFILE, WEIGHT, TRANS,IRTFLG)
        ELSE
C          ONLY USED DURING TESTING!
C          READ REFERENCE IMAGES INTO REFERENCE RINGS (CIRCREF) ARRAY 
           CALL APRINGS_NEW(IREFLIST,NUMREF, NX,NY,
     &               NRING,LCIRC,NUMR, MODE,FFTW_PLANS,
     &               REFPAT,LUNT,CIRCREF,CIRCREF_IN_CORE,
     &               LUNRING,SCRFILE,IRTFLG)
        ENDIF

       !call chkring(2,.false., circref,lcirc, numr,nring, ndum,ndum)
       !call chkray (2,.false., circref,lcirc, numr,nring, ndum,ndum)

	ALLOCATE(REFDIR(3,NUMREF),
     &           ANGREF(3,NUMREF), STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'APSH_PSC; REFDIR, ANGREF',6*NUMREF)
           GOTO 9999
        ENDIF 

C       READ REF. ANGLES INTO ANGREF FROM HEADER OR DOC FILE
C       CONVERT REF. ANGLES TO UNITARY DIRECTIONAL VECTORS (REFDIR).
        CALL AP_GETANGAS(IREFLIST,NUMREF,0,REFANGDOC,REFPAT,
     &                   LUNT,INANG,3,ANGREF,GOTREFANG,NGOTREF,
     &                   .TRUE.,REFDIR,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        NGOTPAR = 0
        IF (EXPANGDOC .NE. CHAR(0)) THEN
	   ALLOCATE(ANGEXP(8,NUMEXP), 
     &              EXPDIR(3,NUMEXP), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'APSH_PSC; ANGEXP....',11*NUMEXP)
              GOTO 9999
           ENDIF 

C          READ EXP. ANGLES INTO ANGEXP
C          CONVERT EXP. ANGLES TO UNITARY DIRECTIONAL VECTORS (EXPDIR).
	   CALL AP_GETANGAS(IEXPLIST,NUMEXP,0,EXPANGDOC,EXPPAT,
     &                      LUNT,INANG,8,ANGEXP,GOTEXPANG,NGOTPAR,
     &                      .TRUE.,EXPDIR,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

        ELSE
C          DUMMY ALLOCATE TO AVOID BUS ERROR ON SOME SYSTEMS
	   ALLOCATE(ANGEXP(8,1), EXPDIR(3,1), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'APSH_PSC; ANGEXP....',11)
              GOTO 9999
           ENDIF 
        ENDIF

        !print *,'gotexpang:',gotexpang,ngotpar
        !print *,'angexp:',angexp(:,1)
        !print *,'expdir:',expdir(:,1)
        !write(6,*) ' Rotfirst,angin:',ROTFIRST,denoise,ANGINHEADER

C       LOOP OVER ALL SETS OF EXPERIMENTAL (SAMPLE) IMAGES -------------


 	DO IEXPT=1,NUMEXP,NUMTH

C          LOAD EXP. IMAGE DATA FOR THIS MP SET OF IMAGES INTO EXPBUF
           IEND = MIN(NUMEXP,IEXPT+NUMTH-1)

           IF (ROTFIRST .AND. DENOISE) THEN
C             ROTATE/SHIFT EXP IMAGES WHEN READING THEM, NEEDS 2XFFT PAD
              !write(6,*) ' calling getdata_den'

	      CALL AP_GETDATA_DEN(IEXPLIST,NUMEXP, 
     &                    NXT,NX,NY, N2XLD,N2X,N2Y, 0.0,
     &                    NUMTH,EXPPAT,LUNT, IEXPT,IEND,
     &                    ANGINHEADER, ANGEXP, 
     &                    MPIBCAST,TMPBUF,EXPBUF,
     &                    FBS_WANTED,IRTFLG)

           ELSEIF (ROTFIRST .AND. .NOT. DENOISE) THEN
C             ROTATE/SHIFT EXP IMAGES WHEN READING THEM, NO PADDING,

	      CALL AP_GETDATA_RTSQ(IEXPLIST,NUMEXP, 
     &                    NX,NY, NX,NY,0.0,
     &                    NUMTH,EXPPAT,LUNT, IEXPT,IEND,
     &                    ANGINHEADER, ANGEXP, 
     &                    MPIBCAST,TMPBUF,EXPBUF,
     &                    .FALSE., ADUM, ADUM, FBS_WANTED,IRTFLG)

          ELSE
C             NO ROTATE, NO PADDING,
	      CALL AP_GETDATA(IEXPLIST,NUMEXP,
     &                    NX,NY, NX,NY,0.0,
     &                    NUMTH,EXPPAT,LUNT, IEXPT,IEND,
     &                    MPIBCAST,  EXPBUF,  
     &                    .FALSE., ADUM, ADUM, IRTFLG)
           ENDIF
           IF (IRTFLG .NE. 0) GOTO 9999

C          NUMTH EXP. IMAGES READY TO BE ALIGNED


C          NUMTH EXP. IMAGES READY TO BE ALIGNED

c$omp      parallel do private(iexp,it)
	   DO IEXP=IEXPT,MIN(NUMEXP,IEXPT+NUMTH-1)
              IT = IEXP - IEXPT + 1

	      CALL APRQ2DC(EXPBUF(1,1,IT),CIRCREF,NUMR,
     &	            NX,NY,ISHRANGEX,ISHRANGEY,ISTEP,
     &	            IRAY1,IRAY2,LCIRC,NRING,NUMREF,MODE,
     &              REFDIR,EXPDIR(1,IEXP),RANGECOS,
     &              DLIST(1,IT),DLIST(2,IT),
     &              DLIST(3,IT),DLIST(4,IT),
     &              DLIST(5,IT),NPROJA(IT),
     &              CKMIRROR,LIMITRANGE,FFTW_PLANS,
     &              COEFFS,IXY,NLOCS,NRAYSC,TRANS,CPLX,USECOEF,
     &              NBORDER,NSUBPIX)
	   ENDDO
c$omp      end parallel do

C          OUTPUT FROM DLIST POSITION 
C          1 - NUMBER OF MOST SIMILAR REFERENCE PROJECTION.
C          2 - NOT-NORMALIZED CORRELATION COEFFICIENT.
C          3 - PSI ANGLE. (IN=PLANE ROTATION) RANGNEW
C          4 - SX SHIFT
C          5 - SY SHIFT
         
           DO IEXP=IEXPT,MIN(NUMEXP,IEXPT+NUMTH-1)
              IT     = IEXP - IEXPT + 1
              IMGEXP = IEXPLIST(IEXP)

C             DLIST(1,IT) IS LIST NUMBER OF MOST SIMILAR REF. IMAGE 
C                 (<0 IF MIRRORED, 0 IF NO SIMILAR IMAGE )

              IREF = INT(DLIST(1,IT))
              IF (IREF < 0) THEN
C                MIRRORED REFERENCE IMAGE
                 IMGREF = IREFLIST(-IREF)

C                IREFT IS FOR REFDIR INDEX
                 IREFT     = -IREF
                 MIRRORNEW = .TRUE.

              ELSEIF (IREF == 0) THEN
C                NO NEARBY REFERENCE IMAGE
                 IMGREF = 0

C                IREFT IS FOR REFDIR INDEX
                 IREFT     = 1
                 MIRRORNEW = .FALSE.

              ELSE
                 IMGREF = IREFLIST(IREF)
C                IREFT IS FOR REFDIR INDEX
                 IREFT     = IREF
                 MIRRORNEW = .FALSE.
              ENDIF
 
              CCROT    = DLIST(2,IT)
              RANGNEW  = DLIST(3,IT)
              XSHRAW   = DLIST(4,IT)
              YSHRAW   = DLIST(5,IT)
              NPROJ    = NPROJA(IT)
              PEAKV    = 1.0

C             WRITE DATA TO DOC FILE
              CALL AP_ENDS(IEXP,IMGEXP,IMGREF,
     &                ANGREF(1,IREFT),REFDIR(1,IREFT),
     &                ANGEXP(1,IEXP), EXPDIR(1,IEXP),ISHRANGEX,
     &                GOTREFANG, NGOTPAR, CCROT,PEAKV,
     &                RANGNEW,XSHRAW,YSHRAW, MIRRORNEW,REFPAT,
     &                NPROJ, CTYPE, LUNDOC,
     &                CHNG_ORDER,SAY_RAW,PARLIST)

C             WRITE DATA TO IMAGE HEADER
              CALL AP_END_HEAD(IMGEXP,EXPPAT,LUNT,PARLIST,8,IRTFLG)

C             UPDATE CCROT & ANGULAR DISPLACEMENT STATISTICS COUNTERS
              CALL AP_STAT_ADD(NGOTPAR,CCROT,PARLIST(10),
     &                       ANGDIFTHR,ANGEXP(8,IEXP),
     &                       CCROTAVG,IBIGANGDIF,ANGDIFAVG,IMPROVCCROT,
     &                       CCROTIMPROV,IWORSECCROT,CCROTWORSE)

	   ENDDO
	ENDDO

        IF (LUNDOC > 0) THEN
C         SAVE CCROT & ANGULAR DISPLACEMENT STATISTICS
          CALL AP_STAT(NUMEXP,ANGDIFTHR,IBIGANGDIF,
     &                 ANGDIFAVG, CCROTAVG,
     &                 IMPROVCCROT,CCROTIMPROV,
     &                 IWORSECCROT,CCROTWORSE,
     &                 NBORDER,NSUBPIX,LUNDOC)
        ENDIF

9999    CONTINUE

C       DEALLOCATE  ARRAYS
        IF (ALLOCATED(DLIST))      DEALLOCATE(DLIST)
        IF (ALLOCATED(NPROJA))     DEALLOCATE(NPROJA)
	IF (ALLOCATED(EXPBUF))     DEALLOCATE(EXPBUF)
	IF (ALLOCATED(REFDIR))     DEALLOCATE(REFDIR)
	IF (ALLOCATED(EXPDIR))     DEALLOCATE(EXPDIR)
	IF (ALLOCATED(ANGEXP))     DEALLOCATE(ANGEXP)
	IF (ALLOCATED(ANGREF))     DEALLOCATE(ANGREF)

        IF (ALLOCATED(NLOCS))      DEALLOCATE(NLOCS)
        IF (ALLOCATED(COEFFS))     DEALLOCATE(COEFFS)
        IF (ALLOCATED(IXY))        DEALLOCATE(IXY)
        IF (ALLOCATED(TMPBUF))     DEALLOCATE(TMPBUF)

       
       END



C+**********************************************************************
C
C APRQ2DC.F
C 
C  PARAMETERS:
C                DIREF    NUMBER OF  MOST SIMILAR REF. PROJ.  (OUTPUT)
C                            (NEGATIVE IF MIRRORED)
C                CCROT    CORR COEFF.                         (OUTPUT)
C                RANGNEW  INPLANE ANGLE                       (OUTPUT)
C                XSHRAW   SHIFT                               (OUTPUT)
C                YSHRAW   SHIFT                               (OUTPUT)
C                NPROJ                                        (OUTPUT)
C
C-**********************************************************************

	SUBROUTINE APRQ2DC(EXPBUF,CIRCREF,NUMR,
     &	             NX,NY,ISHRANGEX,ISHRANGEY,ISTEP,
     &	             IRAY1,IRAY2,LCIRC,NRING,NUMREF,MODE,
     &               REFDIR,EXPDIR,RANGECOS,
     &               DIREF,CCROT,RANGNEW,XSHRAW,YSHRAW,NPROJ,
     &               CKMIRROR,LIMITRANGE,FFTW_PLANS,
     &               COEFFS,IXY,NLOCS,NRAYSC, TRANS,CPLX,USECOEF,
     &               NBORDER,NSUBPIX)

C       NOTE: RUNS WITHIN OMP PARALLEL SECTION OF CODE!
        INCLUDE 'MAKE_CLOSE_LIST.INC'  

	REAL                   :: EXPBUF(NX,NY)
        COMPLEX                :: CIRCREF(LCIRC/2,NUMREF)
        INTEGER                :: NUMR(3,NRING) 
        INTEGER                :: NX,NY,ISHRANGEX,ISHRANGEY,ISTEP 
        INTEGER                :: IRAY1,IRAY2,LCIRC,NRING,NUMREF 
        CHARACTER (LEN=1)      :: MODE
	REAL                   :: REFDIR(3,NUMREF), EXPDIR(3)
	REAL                   :: RANGECOS,DIREF,CCROT,RANGNEW
        REAL                   :: XSHRAW,YSHRAW
        INTEGER                :: NPROJ 
        LOGICAL                :: CKMIRROR,LIMITRANGE
        INTEGER *8             :: FFTW_PLANS(*)  ! STRUCTURE POINTERS
        REAL                   :: COEFFS(6,LCIRC)
        INTEGER                :: IXY(2,LCIRC)
        INTEGER                :: NLOCS(2,*)
        INTEGER                :: NRAYSC
        LOGICAL                :: TRANS          ! TRANSFORMED RAYS
        LOGICAL                :: CPLX           ! COMPLEX CROSRNG
        LOGICAL                :: USECOEF        ! FOR TESTING
        INTEGER                :: NBORDER,NSUBPIX

C       AUTOMATIC ARRAYS
	REAL                   :: CC(-ISTEP:ISTEP,-ISTEP:ISTEP)
	REAL                   :: CCP(-1:1,-1:1)

C       ALLOCATABLE ARRAYS
        REAL, ALLOCATABLE      :: CIRCT(:)
        COMPLEX, ALLOCATABLE   :: CIRCEXP(:)
        COMPLEX, ALLOCATABLE   :: QU(:),QM(:)
        INTEGER, POINTER       :: LCG(:)

	REAL                   :: PEAK
	REAL                   :: CCROT_INTERP,CCNOW
        LOGICAL                :: MIRRORED,ALLRAYS
        LOGICAL                :: ISMIRRORED,USE_UN,USE_MIR
        LOGICAL                :: ISMIRROREDT

        REAL                   :: WR(1)
 	REAL                   :: POS_MAX

        LOGICAL, PARAMETER     :: USE_OMP = .FALSE.

	REAL                   :: ang_n

        PEAK    = 0.0    ! UNUSED?
        WR(1)   = 0.0    ! DUMMY VALUE FLAG FOR APRINGS CALL

        MAXRIN  = NUMR(3,NRING) - 2 ! ACTUAL LENGTH OF LONGEST RING
        LCIRCD2 = LCIRC / 2

C       IF LIMITRANGE, MAKE LIST OF CLOSE REF. IMAGES, RETURNS: NPROJ
        CALL MAKE_CLOSE_LIST(NUMREF,LIMITRANGE,
     &                       REFDIR,EXPDIR,
     &                       RANGECOS, CKMIRROR, 
     &                       LCG, NPROJ, IRTFLG)

        ALLRAYS = (IRAY1 == 1 .AND. IRAY2 == MAXRIN) 

 	IF (NPROJ <= 0) THEN
C          THERE IS NO REFERENCE WITHIN SEARCH RANGE
           XSHRAW  = 0
	   YSHRAW  = 0
           DIREF   = 0     ! FEEDS INTO DLIST
           RANGNEW = 0
           CCROT   = -1.0 
           RETURN	
        ENDIF
	  
	ALLOCATE(CIRCEXP(LCIRCD2),
     &           CIRCT(LCIRC), 
     &           QU(LCIRCD2+1), 
     &           QM(LCIRCD2+1), 
     &           STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'APSH_PSC; CIRCEXP...',4*LCIRC+4)
           RETURN 
        ENDIF
 
	CCROT    = -HUGE(CCROT)    ! HIGHEST CC VALUE

C       SEARCH BOTH MIRRORED & NON-MIRRORED IF CHKMIR
        USE_UN  = .TRUE.
        USE_MIR = CKMIRROR

c       GO THROUGH CENTERS FOR SHIFT ALIGNMENT ------------------------
	DO JT= -ISHRANGEY,ISHRANGEY,ISTEP
	   CNR2 = NY / 2 + 1 + JT

	   DO IT= -ISHRANGEX,ISHRANGEX,ISTEP
	      CNS2 = NX / 2 + 1 + IT

C             NORMALIZE EXP VALUES UNDER MASK OVER VARIANCE RANGE.
C             INTERPOLATE TO POLAR COORD. AROUND: CNS2,CNR2.
C             CREATE FOURIER OF RINGS, NO WEIGHTING OF RINGS        

              CALL APRINGS_SATU(EXPBUF,NX,NY,CNS2,CNR2, 
     &                          MODE,NUMR,NRING,LCIRC,USE_OMP,
     &                          WR,FFTW_PLANS, NLOCS,NRAYSC,
     &                          COEFFS,IXY, USECOEF,TRANS,CPLX,
     &                          CIRCT,  CIRCEXP, IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

  
C             COMPARE EXP. RING WITH ALL REFERENCE RINGS ------------

	      DO IRR=1,NPROJ  ! LOOP OVER ALL REF. IMAGES
                 IR = IRR
                 IF (LIMITRANGE) THEN
C                    ONLY CHECKING NEARBY PROJECTIONS
                     IR = ABS(LCG(IRR))
                     IF (CKMIRROR ) THEN
C                       ONLY SEARCH EITHER MIRRORED OR NON-MIRRORED
                        USE_UN  = (LCG(IRR) >= 0)
                        USE_MIR = (LCG(IRR) <  0)
                     ENDIF
                 ENDIF

C                CROSS CORRELATION PEAK DETERMINATION
#ifdef DEBUGT  
                 CALL CROSRNG_SATU(CIRCEXP,CIRCREF(1,IR),
     &                        NUMR,NRING,LCIRC,
     &                        FFTW_PLANS, NLOCS,NRAYSC,
     &                        TRANS,CPLX, USE_UN,USE_MIR,
     &                        ISMIRRORED,CCNOW,POS_MAX, IRTFLG)
#else     
                 IF (ALLRAYS) THEN
                    CALL CROSRNG_2CQR(CIRCEXP,CIRCREF(1,IR),LCIRC/2,
     &                        NRING,MAXRIN,NUMR,
     &                        .FALSE., FFTW_PLANS, 
     &                        USE_UN,USE_MIR,
     &                        QU,QM,
     &                        ISMIRRORED,CCNOW,POS_MAX, MAXL)
                 ELSE
                    CALL CROSRNG_2CQRR(CIRCEXP,CIRCREF(1,IR),LCIRC/2,
     &                        NRING,MAXRIN,NUMR,IRAY1,IRAY2,
     &                        .FALSE., FFTW_PLANS, 
     &                        USE_UN,USE_MIR,
     &                        QU,QM,
     &                        ISMIRRORED,CCNOW,POS_MAX, MAXL)
                 ENDIF
#endif

                 IF (CCNOW .GE. CCROT)  THEN
C                   THIS POSITION IS BEST SO FAR OVER ALL REFS & SHIFTS
	            RANGNEW  = ANG_N(POS_MAX,MODE,MAXRIN)! BEST ROTATION
C                   write(6,'(a,i5,3f9.2,3i4)')' rangnew0:',
C     &                      maxl,pos_max,rangnew,ccnow,it,jt
	            CCROT    = CCNOW  ! BEST CC VALUE
	            IBE      = IR     ! BEST REFERENCE IMAGE
	            ISX      = IT     ! BEST X SHIFT
	            ISY      = JT     ! BEST Y SHIFT
	            IDIS     = IR     ! BEST MIRROR FLAG & REF. #
                    IF (ISMIRRORED) IDIS = -IR
	         ENDIF

	      ENDDO ! END OF:  DO IRR=1,NPROJ --------------------- REFS.
           ENDDO    ! END OF:  DO IT=-ISHRANGEX,ISHRANGEX,ISTEP
	ENDDO       ! END OF:  DO JT=-ISHRANGEY,ISHRANGEY,ISTEP

        !write(6,'(a,2i4,3x,f8.2)') ' 0 X,Y:',isx,isy,ccrot 
        !write(6,*) ' '
        !write(6,'(a,2f9.2,i5)') ' rangnew1:',rangnew,CCROT

        ISTEPX = ISTEP
        IF (ISHRANGEX <= 0) ISTEPX = 0
        ISTEPY = ISTEP
        IF (ISHRANGEY <= 0) ISTEPY = 0

        IF (ISTEP > 1) THEN
C          CHECK NEIGHBORING LOCATIONS OF BEST REF's BEST CC'S

           IR   = IBE     ! BEST REF IMAGE
	   ISXT = ISX     ! BEST X SHIFT FOR THIS REF IMAGE
           ISYT = ISY     ! BEST Y SHIFT FOR THIS REF IMAGE

           !write(6,'a,3i4,3x,f8.2') ' Ref:',ir,isxt,isyt, ccrot 

C          WHAT IF X OR Y RANGE IS ZERO?
           JT1 = -ISTEP+1
           JT2 =  ISTEP-1
           IF (ISHRANGEY <= 0) THEN
               JT1 = 0
               JT2 = 0
           ENDIF
           IT1 = -ISTEP+1
           IT2 =  ISTEP-1
           IF (ISHRANGEX <= 0) THEN
               IT1 = 0
               IT2 = 0
           ENDIF

           DO JT= JT1,JT2                          ! OVER Y
	      CNR2 = NY / 2 + 1 + JT + ISYT         ! CENTER

	      DO IT= IT1,IT2                        ! OVER X
	         CNS2 = NX / 2 + 1 + IT + ISXT      ! CENTER

   	         IF (IT == 0 .AND. JT == 0) CYCLE     ! ALREADY HAVE 
                 IF (((ISXT + IT) < -ISHRANGEX) .OR.
     &               ((ISXT + IT) >  ISHRANGEX) .OR.
     &               ((ISYT + JT) < -ISHRANGEY) .OR.
     &               ((ISYT + JT) >  ISHRANGEY)) CYCLE ! OUT OF SHIFT 

                 !write(6,'a,2i5') ' Checking: ',isxt+it,isyt+jt

C                NORMALIZE EXP VALUES UNDER MASK OVER VARIANCE RANGE.
C                INTERPOLATE TO POLAR COORD. AROUND: CNS2,CNR2.
C                CREATE FOURIER OF RINGS, NO WEIGHTING OF RINGS        
                 CALL APRINGS_SATU(EXPBUF,NX,NY,CNS2,CNR2, 
     &                     MODE,NUMR,NRING,LCIRC,USE_OMP,
     &                     WR,FFTW_PLANS, NLOCS,NRAYSC,
     &                     COEFFS,IXY, USECOEF,TRANS,CPLX,
     &                     CIRCT,  CIRCEXP, IRTFLG)
                 IF (IRTFLG .NE. 0) RETURN
 
#ifdef DEBUGT  
C                CROSS CORRELATION PEAK DETERMINATION
                 CALL CROSRNG_SATU(CIRCEXP,CIRCREF(1,IR),
     &                             NUMR,NRING,LCIRC,
     &                             FFTW_PLANS, NLOCS,NRAYSC,
     &                             TRANS,CPLX, USE_UN,USE_MIR,
     &                             ISMIRROREDT,CC(IT,JT),POS_MAX,IRTFLG)
#else
                 IF (ALLRAYS) THEN
                    CALL CROSRNG_2CQR(CIRCEXP,CIRCREF(1,IR),LCIRC/2,
     &                             NRING,MAXRIN,NUMR,
     &                             .FALSE., FFTW_PLANS, 
     &                             USE_UN,USE_MIR,
     &                             QU,QM,
     &                             ISMIRROREDT,CC(IT,JT),
     &                             POS_MAX, MAXL)
                 ELSE
                    CALL CROSRNG_2CQRR(CIRCEXP,CIRCREF(1,IR),LCIRC/2,
     &                             NRING,MAXRIN,NUMR,IRAY1,IRAY2,
     &                             .FALSE., FFTW_PLANS, 
     &                             USE_UN,USE_MIR,
     &                             QU,QM,
     &                             ISMIRROREDT,CC(IT,JT),
     &                             POS_MAX, MAXL)
                 ENDIF
#endif
                 IF (IRTFLG .NE. 0) RETURN

	         IF (CC(IT,JT) > CCROT)  THEN
C                   CC VALUE HIGHER THAN CURRENT BEST
	            RANGNEW = ANG_N(POS_MAX ,MODE,MAXRIN)
	            CCROT   = CC(IT,JT) ! BEST CC VALUE
                    IBE     = IR        ! BEST REF. # 
	            ISX     = IT + ISXT ! BEST X SHIFT
	            ISY     = JT + ISYT ! BEST Y SHIFT
                    IBE     = IR        ! BEST REF. # 
 	            IDIS    = IR        ! BEST MIRROR FLAG & REF. #
                    IF (ISMIRROREDT) IDIS = -IR
	         ENDIF

	      ENDDO     ! END OF:  DO IT=-ISTEP,ISTEP
           ENDDO        ! END OF:  DO JT=-ISTEP,ISTEP
      
        ELSE
C          NO NEED TO CHECK STEPPED SHIFT LOCATIONS  
        ENDIF
        !write(6,'(a,2f8.2)') ' rangnew2:',rangnew,pos_max

        SX       = ISX              ! BEST X SHIFT
        SY       = ISY              ! BEST Y SHIFT
        DIREF    = IDIS             ! REAL AS FEEDS INTO DLIST
        MIRRORED = (IDIS .LT. 0 )   ! CHECK MIRRORED ONLY IF(IDIS.LT.0) 

        !write(6,'a,2i4,3x,f8.2') ' 1 X,Y:',isx,isy,ccrot 

        IF (ABS(ISX) == ISHRANGEX .OR. ABS(ISY) == ISHRANGEY) THEN

           NBORDER = NBORDER + 1    ! ON BORDER, SKIP SUB-PIXEL
        ELSE
C          NOTE: RECOVER ARRAY FOR SUBPIXEL CCROT AROUND BEST SHIFT
C          NOTE: NO EASY WAY TO RECOVER SUBPIXEL ARRAY EVEN IF ISTEP=1!!

C          LOAD CC VALUES FOR SHIFTS ADJACENT TO MAX CC SHIFT 
	   CCP(0,0) = CCROT

           DO JT= -1,1
	      CNR2 = NY / 2 + 1 + JT + ISY     ! CENTER
	      DO IT= -1,1
	         CNS2 = NX / 2 + 1 + IT + ISX  ! CENTER

   	         IF (IT == 0 .AND. JT == 0) CYCLE ! ALREADY HAVE THIS

C                NORMALIZE EXP VALUES UNDER MASK OVER VARIANCE RANGE.
C                INTERPOLATE TO POLAR COORD. AROUND: CNS2,CNR2.
C                CREATE FOURIER OF RINGS, NO WEIGHTING OF RINGS        
                 CALL APRINGS_SATU(EXPBUF,NX,NY,CNS2,CNR2, 
     &                        MODE,NUMR,NRING,LCIRC,USE_OMP,
     &                        WR,FFTW_PLANS, NLOCS,NRAYSC,
     &                        COEFFS,IXY, USECOEF,TRANS,CPLX,
     &                        CIRCT, CIRCEXP,   IRTFLG)
 
C                CROSS CORRELATION PEAK DETERMINATION
#ifdef DEBUG  
                 CALL CROSRNG_SATU(CIRCEXP,CIRCREF(1,IBE),
     &                        NUMR,NRING,LCIRC,
     &                        FFTW_PLANS, NLOCS,NRAYSC,
     &                        TRANS,CPLX, .NOT. MIRRORED,MIRRORED,
     &                        ISMIRRORED,CCP(IT,JT), POS_MAX,IRTFLG)
#else
                IF (ALLRAYS ) THEN
                    CALL CROSRNG_2CQR(CIRCEXP,CIRCREF(1,IBE),LCIRC/2,
     &                        NRING,MAXRIN,NUMR,
     &                        .FALSE., FFTW_PLANS, 
     &                        .NOT.MIRRORED,MIRRORED,
     &                        QU,QM,
     &                        ISMIRROREDT,CCP(IT,JT),
     &                        POS_MAX, MAXL)
                ELSE
                   CALL CROSRNG_2CQRR(CIRCEXP,CIRCREF(1,IBE),LCIRC/2,
     &                        NRING,MAXRIN,NUMR,IRAY1,IRAY2,
     &                        .FALSE., FFTW_PLANS, 
     &                        .NOT.MIRRORED,MIRRORED,
     &                        QU,QM,
     &                        ISMIRROREDT,CCP(IT,JT),
     &                        POS_MAX, MAXL)
                 ENDIF
#endif

	      ENDDO     ! END OF:  DO IT=-ISTEP,ISTEP
           ENDDO        ! END OF:  DO JT=-ISTEP,ISTEP

777        continue


C          SUB-PIXEL INTERPOLATION -----------------------------------

           !write(6,'(a,2f8.2)') ' rangnew3:',rangnew,pos_max 
           !write(6,*) 'cc' write(6,910) cc
910        format(3f14.4,/,3f14.4,/,3f14.4)

C          SUB-PIXEL INTERPOLATION 
           CALL PARABL(CCP,SSX,SSY,PEAK)  ! RETURNS: SSX,SSY,PEAK 

	   CNS2 = NX/2+1 + SX + SSX
	   CNR2 = NY/2+1 + SY + SSY

           !write(6,'a,2f8.2,3x,2f8.2')' 3 X,Y:',sx+ssx,sy+ssy,CNR2,CNS2 
        !write(6,'(a,f8.2)') ' rangnew4:',rangnew 

C          NORMALIZE EXP VALUES UNDER MASK OVER VARIANCE RANGE.
C          INTERPOLATE TO POLAR COORD. AROUND: CNS2,CNR2.
C          CREATE FOURIER OF RINGS, NO WEIGHTING OF RINGS        
C          CAN NOT USE: APRINGS_ONE_COEF AS NOT INTEGRAL SHIFT.

           CALL APRINGS_SATU(EXPBUF,NX,NY,CNS2,CNR2, 
     &                     MODE,NUMR,NRING,LCIRC,USE_OMP,
     &                     WR,FFTW_PLANS, NLOCS,NRAYSC,
     &                     COEFFS,IXY, .FALSE.,TRANS,CPLX,
     &                     CIRCT,CIRCEXP,   IRTFLG)

C          CROSS CORRELATION PEAK DETERMINATION
#ifdef DEBUGT  
           CALL CROSRNG_SATU(CIRCEXP,CIRCREF(1,IBE),
     &                       NUMR,NRING,LCIRC,
     &                       FFTW_PLANS, NLOCS,NRAYSC,
     &                       TRANS,CPLX, .NOT. MIRRORED,MIRRORED,
     &                       ISMIRRORED,CCROT_INTERP,POS_MAX_INTERP, 
     &                       IRTFLG)
#else
           IF ( ALLRAYS) THEN
               CALL CROSRNG_2CQR(CIRCEXP,CIRCREF(1,IBE),LCIRC/2,
     &                       NRING,MAXRIN,NUMR,
     &                       .FALSE., FFTW_PLANS, 
     &                       .NOT. MIRRORED,MIRRORED,
     &                       QU,QM,
     &                       ISMIRRORED,CCROT_INTERP,
     &                       POS_MAX_INTERP, MAXL)
           ELSE
              CALL CROSRNG_2CQRR(CIRCEXP,CIRCREF(1,IBE),LCIRC/2,
     &                       NRING,MAXRIN,NUMR,iray1,iray2,
     &                       .FALSE., FFTW_PLANS, 
     &                       .NOT. MIRRORED,MIRRORED,
     &                       QU,QM,
     &                       ISMIRRORED,CCROT_INTERP,
     &                       POS_MAX_INTERP, MAXL)
           ENDIF
#endif
 
           IF (CCROT_INTERP > CCROT) THEN
C                USE SUB-PIXEL LOCATION INSTEAD OF INTEGRAL LOCATION

              !write(6,961) ssx,ssy, ccrot_interp,ccrot 
961           format(' Subpixel:',2f8.2,'    ',f9.2,' >> ',f9.2)

	      RANGNEW = ANG_N(POS_MAX_INTERP,MODE,MAXRIN) ! BEST ROT.
              CCROT   = CCROT_INTERP   ! HIGHEST CC 
	      SX      = SX + SSX       ! BEST X SHIFT
              SY      = SY + SSY       ! BEST Y SHIFT
              NSUBPIX = NSUBPIX + 1    ! SUB-PIXEL COUNTER

           ENDIF  ! END OF: (CCROT_INTERP > CCROT)
        ENDIF     ! END OF: IF (ABS(ISX).NE. ISHRANGEX.......

        XSHRAW = SX
        YSHRAW = SY

C       CHANGE ORDER OF SHIFT & ROTATION IN: AP_ENDS
        !write(6,'(a,2i4,f8.2)') 'X,Y:',isx,isy,ccrot 

9999    IF (ASSOCIATED(LCG))    DEALLOCATE(LCG)
        IF (ALLOCATED(CIRCEXP)) DEALLOCATE(CIRCEXP)
        IF (ALLOCATED(QM))      DEALLOCATE(QM)
        IF (ALLOCATED(QU))      DEALLOCATE(QU)
        IF (ALLOCATED(CIRCT))   DEALLOCATE(CIRCT)
        NULLIFY(LCG)

	END


#ifdef NEVER
           IF (.FALSE.) THEN !!!!!!!!!!!!!!!!!!!!!!!!!!
C             DENOISE THE EXP IMAGES

c$omp         parallel do private(iexp,it)
	      DO IEXP=IEXPT,MIN(NUMEXP,IEXPT+NUMTH-1)
                 IT = IEXP - IEXPT + 1

C                MEDIAN FILTER OVER 3x3 SQUARE BOX
                 CALL MD2_NOLUN(EXPBUF(1,1,IT),TMPBUF(1,1,IT), IT!!!!
     &                          NX,NY, 3,9,'B')
                 EXPBUF(1:NX,1:NY,IT) = TMPBUF(1:NX,1:NY,IT)

              ENDDO
           ENDIF
#endif
