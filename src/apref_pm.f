 
C++********************************************************************
C
C    APREF_PM.F
C                   ADDED APSHIFT                 NOV 03 ARDEAN LEITH
C                   AP_END                        FEB 04 ARDEAN LEITH
C                   CROSRNG_E SPEEDS UP           AUG 04 ARDEAN LEITH
C                   AP_END CALL HAS PARLIST       OCT 04 ARDEAN LEITH
C                   CCROT STATISTICS              FEB 05 ARDEAN LEITH
C                   APRINGS_ONE                   MAR 05 ARDEAN LEITH
C                   CROSRNG_E REWRITE             MAR 08 ARDEAN LEITH
C                   IREFLIST VAR. NAME            MAR 08 ARDEAN LEITH
C                   FMRS_PLAN                     MAR 08 ARDEAN LEITH
C                   APRINGS_ONE_NEW               MAR 08 ARDEAN LEITH
C                   AP_END REDO                   NOV 08 ARDEAN LEITH
C                   AP_STAT_ADD                   NOV 08 ARDEAN LEITH
C                   ANGDIF EXPDIR BUG             MAR 09 ARDEAN LEITH
C                   CROSRNG_2, TT REMOVED         JUN 10 ARDEAN LEITH
C                   AP_STAT NBORDER               OCT 10 ARDEAN LEITH
C                   RENAMED FROM DSGR_PM          JAN 11 ARDEAN LEITH
C                   ROTFIRST                      FEB 11 ARDEAN LEITH
C                   MAKE_CLOSE_LIST, GETANGAS     FEB 11 ARDEAN LEITH
C                   AP_GETDATA USED               NOV 11 ARDEAN LEITH
C                   FBS_WANTED                    JAN 12 ARDEAN LEITH
C                   ANGINHEADER BUG               APR 15 ARDEAN LEITH
C                   AP_STAT_R USED                APR 15 ARDEAN LEITH
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
C  APREF_PM
C
C  PURPOSE: FIND ROTATIONAL AND SHIFT PARAMETERS TO ALIGN A SERIES OF
C           REFERENCE IMAGES WITH SAMPLE IMAGES  'AP REF'
C           USED WHEN: CTYPE == 'T' OR
C           (CIRCREF_IN_CORE AND NUMBER OF THREADS > 1 AND
C           NUMBER OF EXP IMAGES > NUMBER OF THREADS).  THIS OCCURS
C           ALMOST ALWAYS ON MULTIPORCESSOR MACHINES EXCEPT WHEN 
C           USING A SINGLE EXP. IMAGE. 
C
C  ORDER OF PROCESSING:
C 1 - CONVERT EACH REFERENCE IMAGE TO RINGS, DO THE FFTS
C     FOR ALL THE RINGS, APPLY WEIGHTS TO THE RINGS, STORE IT IN CIRCREF.
C     IN ADDITION, HIGHEST FREQUENCY FOR ALL THE RINGS EXCEPT
C     MAXRING ARE DIVIDED BY 2.
C 2 - CONVERT EACH INPUT IMAGE TO RINGS, DO THE FFTS,
C     COMPARE EACH INPUT IMAGE WITH ALL THE REFERENCE IMAGES
C     USING CROSRNG_E. 
C     SINCE Y WERE PRE-WEIGHTED THE RESULT IS ALREADY CORRECT.
C
C PARAMETERS:
C       IREFLIST            LIST OF REF. IMAGE FILE NUMBERS   (INPUT)
C       NUMREF              NO. OF IMAGES                     (INPUT)
C       IEXPLIST            LIST OF EXP. IMAGE FILE NUMBERS   (INPUT)
C       NUMEXP              NO. OF IMAGES                     (INPUT)
C       NX,NY               ACTUAL (NOT-WINDOWED) IMAGE SIZE  (INPUT)
C       REFANGDOC           REF. ANGLES FILE NAME             (INPUT)
C       EXPANGDOC           EXP. ANGLES FILE NAME             (INPUT)
C       REFPAT              REF. IMAGE SERIES FILE TEMPLATE   (INPUT)
C       EXPPAT              EXP. IMAGE SERIES FILE TEMPLATE   (INPUT)
C
C  OPERATIONS:  'AP REF'
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE APREF_PM(IREFLIST,NUMREF,IEXPLIST,NUMEXP,
     &             NX,NY,ANGDIFTHR,RANGE,
     &             NRING,LCIRC,NUMR,CIRCREF,
     &             MODE, REFANGDOC,EXPANGDOC,SCRFILE,FFTW_PLANS,
     &             REFPAT,EXPPAT,CKMIRROR,CTYPE,
     &             ROTFIRST,ISHRANGE,LUNDOC,FBS_WANTED)

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        INTEGER                     :: NUMREF,NUMEXP
	INTEGER                     :: IREFLIST(NUMREF) 
	INTEGER                     :: IEXPLIST(NUMEXP) 
        INTEGER                     :: NX,NY
	REAL                        :: ANGDIFTHR,RANGE
        INTEGER                     :: NRING,LCIRC
        INTEGER                     :: NUMR(3,NRING)
	REAL                        :: CIRCREF(LCIRC,NUMREF)
        CHARACTER (LEN=1)           :: MODE
        CHARACTER (LEN=*)           :: REFANGDOC,EXPANGDOC,SCRFILE
        INTEGER*8, INTENT(IN)       :: FFTW_PLANS(*)
        CHARACTER (LEN=*)           :: REFPAT,EXPPAT
        LOGICAL                     :: CKMIRROR
        CHARACTER (LEN=*)           :: CTYPE
        LOGICAL                     :: ROTFIRST
        INTEGER                     :: ISHRANGE,LUNDOC
        LOGICAL                     :: FBS_WANTED,FIRST,FIRSTCC

C       ALLOCATABLE ARRAYS
	REAL,    ALLOCATABLE        :: EXPBUF(:,:,:)
	REAL,    ALLOCATABLE        :: AVI(:),SIGI(:)
	REAL,    ALLOCATABLE        :: CCLIST(:) 
	REAL,    ALLOCATABLE        :: RANGNLIST(:) 
	INTEGER, ALLOCATABLE        :: NPROJ(:)
	INTEGER, ALLOCATABLE        :: ILIST(:) 
	LOGICAL, ALLOCATABLE        :: MIRLIST(:) 
	REAL,    ALLOCATABLE        :: TMPBUF(:,:)

        INTEGER *8                  :: IPLAN ! STRUCTURE POINTER

        CHARACTER (LEN=1)           :: NULL = CHAR(0)
        LOGICAL                     :: CIRCREF_IN_CORE
        LOGICAL                     :: MIRRORNEW
        LOGICAL                     :: GOTREFANG,GOTEXPANG
        LOGICAL,PARAMETER           :: ANGINHEADER = .FALSE.

C       AUTOMATIC ARRAYS
	REAL                        :: DLIST(6)
	REAL                        :: ANGEXP(8,NUMEXP)
	REAL                        :: EXPDIR(3,NUMEXP)
	REAL                        :: REFDIR(3,NUMREF) 
	REAL                        :: ANGREF(3,NUMREF)

        INTEGER, PARAMETER          :: NLISTMAX = 15
        REAL                        :: PARLIST(NLISTMAX) 

	REAL, PARAMETER             :: QUADPI = 3.14159265358979323846
	REAL, PARAMETER             :: DGR_TO_RAD = (QUADPI/180)

        INTEGER                     :: NBORDER = 0       ! UNUSED
        INTEGER                     :: NSUBPIX = 0       ! UNUSED

        INTEGER, PARAMETER          :: INPIC   = 77
        INTEGER, PARAMETER          :: INANG   = 78
        INTEGER, PARAMETER          :: LUNRING = 50

C       FIND NUMBER OF OMP THREADS 
        CALL GETTHREADS(NUMTH) 

	RANGECOS = COS(RANGE*DGR_TO_RAD)

C       MAKE EXPBUF BIG ENOUGH FOR ALL THREADS
        MWANTX = NX * NY * NUMTH

	ALLOCATE(EXPBUF(NX,NY,NUMTH), 
     &           AVI(NUMTH),
     &           SIGI(NUMTH),
     &           ILIST(NUMTH),
     &           CCLIST(NUMTH),
     &           RANGNLIST(NUMTH),
     &           NPROJ(NUMTH),
     &           MIRLIST(NUMTH),
     &           STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           MWANT = MWANTX + 7*NUMTH 
           CALL ERRT(46,'EXPBUF...',MWANT)
           GOTO 9999
        ENDIF
        IF (ROTFIRST) THEN
	   ALLOCATE(TMPBUF(NX,NY), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'APREF_PM; TMPBUF',NX*NY)
              GOTO 9999
           ENDIF
           IF (FBS_WANTED) THEN
            WRITE(NOUT,*)' ALIGNING INPUT IMAGES WITH FBS INTERPOLATION'
           ELSE
           WRITE(NOUT,*)' ALIGNING INPUT IMAGES WITH QUAD INTERPOLATION'
           ENDIF
        ENDIF 

C       DUMMY CALL TO INITIALIZE REV. FFTW3 PLAN FOR USE WITHIN OMP ||
 	CALL FMRS_PLAN(.FALSE.,DUM,NX,NY,1, 1, -1, IPLAN, IRTFLG)
 	CALL FMRS_PLAN(.FALSE.,DUM,NX,NY,1, 1, +1, IPLAN, IRTFLG)

C       ZERO DLIST ARRAY
	DLIST = 0.0

C       GET REF PROJ. ANGLES FROM DOC FILE  OR IMAGE HEADER (IF WANTED)
C       CONVERT REF ANGLES TO UNITARY DIRECTIONAL VECTORS (REFDIR).
        CALL AP_GETANGAS(IREFLIST,NUMREF,0,REFANGDOC,REFPAT,
     &                   INPIC,INANG,3,ANGREF,GOTREFANG,NGOTREF,
     &                   .TRUE.,REFDIR,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       GET EXP PROJ. ANGLES FROM DOC FILE  OR IMAGE HEADER (IF WANTED)
C       CONVERT EXP. ANGLES TO UNITARY DIRECTIONAL VECTORS (EXPDIR).
        IRTFLG = -8999   ! WANT CC MIR NOT CCROT
        CALL AP_GETANGAS(IEXPLIST,NUMEXP,0,EXPANGDOC,EXPPAT,
     &                   INPIC,INANG,8,ANGEXP,GOTEXPANG,NGOTPAR,
     &                  .TRUE.,EXPDIR,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        ! NEED EXISTING ALIGN FILE FOR CC DIFFERENCES
        FIRST   = (NGOTPAR == 0)        ! FIRST FILE, NO DIFFERENCES

        ! NEED EXISTING CC CONTAINING ALIGN FILE FOR CC DIFFERENCES
        FIRSTCC = ANGEXP(8,1) ==  1.0 .OR. 
     &            ANGEXP(8,1) == -1.0 

        !write(6,*) 'angexp:',angexp(:,1), firstcc


C       READ REF IMAGES INTO REF RINGS (CIRCREF) ARRAY 
        CIRCREF_IN_CORE = .TRUE.

        CALL APRINGS_NEW(IREFLIST,NUMREF, NX,NY,
     &               NRING,LCIRC,NUMR,MODE,FFTW_PLANS, 
     &               REFPAT,INPIC,CIRCREF,CIRCREF_IN_CORE,
     &               LUNRING,SCRFILE,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       INITIALIZE CC STATISTICS
        CALL AP_STAT_ADD_R(-1,ANGDIF,ANGDIFTHR, IBIGANGDIF,ANGDIFAVG,
     &               CC,CCLAS,CCAVG,IMPROVCC,CCIMPROV,IWORSECC,CCWORSE)


        DO IEXPT=1,NUMEXP,NUMTH
C         LOOP OVER ALL EXPERIMENTAL (SAMPLE) IMAGES

C         LOAD NUMTH EXPERIMENTAL IMAGES INTO ARRAY EXPBUF
          IEND = MIN(NUMEXP,IEXPT+NUMTH-1)

          IF (ROTFIRST) THEN
C            WANT TO ROTATE/SHIFT EXP IMAGES WHEN READING THEM
	     CALL AP_GETDATA_RTSQ(IEXPLIST,NUMEXP,
     &                        NX,NY, NX,NY, 0.0,
     &                        NUMTH,EXPPAT,INPIC, IEXPT,IEND,
     &                        ANGINHEADER, ANGEXP, 
     &                        .TRUE.,TMPBUF,EXPBUF,
     &                        .TRUE.,AVI,SIGI,FBS_WANTED,IRTFLG)
          ELSE
	     CALL AP_GETDATA(IEXPLIST,NUMEXP,
     &                       NX,NY, NX,NY, 0.0,
     &                       NUMTH,EXPPAT,INPIC, IEXPT,IEND,
     &                       .TRUE., EXPBUF,
     &                       .TRUE.,AVI,SIGI, IRTFLG)
          ENDIF
          IF (IRTFLG .NE. 0) GOTO 9999

C         FIND CLOSEST MATCHING REFERENCE IMAGE FOR EACH EXP. IMAGE
c$omp     parallel do private(IEXP,IT)
	  DO IEXP=IEXPT,IEND
             IT  = IEXP-IEXPT+1
	     CALL APREF_2(EXPBUF(1,1,IT),NX,NY,NUMR,NRING,
     &              MODE,CIRCREF,LCIRC,NUMREF,
     &              REFDIR,EXPDIR(1,IEXP),RANGECOS,
     &	            ILIST(IT),CCLIST(IT),RANGNLIST(IT),MIRLIST(IT),
     &              NPROJ(IT),CKMIRROR,FFTW_PLANS,IRTFLG)
	  ENDDO

C         LOOP OVER CURRENT SET OF EXPERIMENTAL (SAMPLE) IMAGES
          DO IEXP=IEXPT,MIN(NUMEXP,IEXPT+NUMTH-1)

C            IT POINTS TO CURRENT EXP. IMAGE IN THE SUBLIST
             IT = IEXP - IEXPT + 1

C            ILIST(IT) IS LIST NUMBER OF MOST SIMILAR REF. IMAGE 
             IREF      = ILIST(IT)

             IF (IREF <= 0) THEN
C               NO NEARBY REFERENCE IMAGE
                IMGREF = 0
C               IREFT IS FOR REFDIR INDEX
                IREFT  = 1
             ELSE
                IMGREF = IREFLIST(IREF)
C               IREFT IS FOR REFDIR INDEX
                IREFT  = IREF
             ENDIF
  
             IMGEXP    = IEXPLIST(IEXP)
             CCROT     = CCLIST(IT)
             RANGNEW   = RANGNLIST(IT)
             MIRRORNEW = MIRLIST(IT) 
             PEAKV     = 0.0
             XSHNEW    = 0.0
             YSHNEW    = 0.0

             IF (IMGREF > 0 .AND. ISHRANGE > 0) THEN
C               DETERMINE SHIFT PARAMETERS FOR EXP. IMAGE IN EXPBUF 
                NXP = 2*NX+2
                NYP = 2*NY

                CALL APSHIFT(INPIC,  REFPAT,IMGREF,
     &                 NX,NY, NXP,NYP,
     &                 EXPBUF(1,1,IT),AVI(IT),SIGI(IT), ISHRANGE,
     &                 RANGNEW,XSHNEW,YSHNEW,MIRRORNEW,PEAKV,IRTFLG)
                IF (IRTFLG .NE. 0) RETURN
             ENDIF

C            AP_END WRITES ALIGNMENT PARAMETERS TO DOC FILE 
             CALL AP_END(IEXP,IMGEXP,IMGREF, 
     &            ANGREF(1,IREFT),REFDIR(1,IREFT),
     &            ANGEXP(1,IEXP),EXPDIR(1,IEXP),ISHRANGE,
     &            GOTREFANG, NGOTPAR, CCROT,PEAKV,
     &            RANGNEW,XSHNEW,YSHNEW,MIRRORNEW,REFPAT,
     &            NPROJ(IT), CTYPE, LUNDOC,PARLIST)

C           WRITE ALIGNMENT DATA TO IMAGE HEADER
            CALL AP_END_HEAD(IMGEXP,EXPPAT,INPIC,PARLIST,8,IRTFLG)

            CALL AP_STAT_ADD_R(NGOTPAR,
     &            PARLIST(10),ANGDIFTHR,IBIGANGDIF,ANGDIFAVG,
     &            PEAKV,ANGEXP(8,IEXP),CCAVG,
     &            IMPROVCC,CCIMPROV,IWORSECC,CCWORSE ,FIRSTCC)
          ENDDO
       ENDDO

      IF (LUNDOC > 0) THEN
C         SAVE CC & ANGULAR DISPLACEMENT STATISTICS IN DOC FILE

          CALL AP_STAT_R(NUMEXP,ANGDIFTHR,IBIGANGDIF,ANGDIFAVG, 
     &                   CCAVG,IMPROVCC,CCIMPROV, IWORSECC,CCWORSE,
     &                   FIRST,FIRSTCC, LUNDOC)
 
       ENDIF

9999   IF (ALLOCATED(EXPBUF))    DEALLOCATE(EXPBUF)
       IF (ALLOCATED(AVI))       DEALLOCATE(AVI)
       IF (ALLOCATED(SIGI))      DEALLOCATE(SIGI)
       IF (ALLOCATED(ILIST))     DEALLOCATE(ILIST)
       IF (ALLOCATED(CCLIST))    DEALLOCATE(CCLIST)
       IF (ALLOCATED(RANGNLIST)) DEALLOCATE(RANGNLIST)
       IF (ALLOCATED(MIRLIST))   DEALLOCATE(MIRLIST)
       IF (ALLOCATED(NPROJ))     DEALLOCATE(NPROJ)
       IF (ALLOCATED(TMPBUF))    DEALLOCATE(TMPBUF)

       END



C++************************** APREF_2 **********************************


	SUBROUTINE APREF_2(XIM,NX,NY, NUMR,NRING,MODE,
     &		         CIRCREF,LCIRC,NUMREF,
     &                   REFDIR,EXPDIR,RANGECOS,
     &                   IMGREFL,CCROT,RANGNEW,MIRNEW,NPROJ,
     &                   CKMIRROR,FFTW_PLANS,IRTFLG)

C       NOTE: RUNS WITHIN OMP PARALLEL SECTION OF CODE!
        INCLUDE 'MAKE_CLOSE_LIST.INC'  

        REAL                    :: XIM(NX,NY)
        INTEGER                 :: NUMR(3,NRING)
	CHARACTER (LEN=1)       :: MODE
	REAL                    :: CIRCREF(LCIRC,NUMREF)
	REAL                    :: REFDIR(3,NUMREF)
	REAL                    :: EXPDIR(3)
	REAL                    :: RANGECOS 
	LOGICAL                 :: CKMIRROR,ONLYMIRROR
	LOGICAL                 :: MIRNEW,LIMITRANGE
	LOGICAL                 :: ISMIRRORED
        LOGICAL                 :: USE_UN,USE_MIR
        INTEGER *8              :: FFTW_PLANS(*)

	DOUBLE PRECISION        :: CCROTD,TOTMIN
        INTEGER, POINTER        :: LCG(:)

C       ALLOCATABLE ARRAYS
        REAL, ALLOCATABLE       :: CIRCEXP(:)


        MAXRIN     = NUMR(3,NRING) - 2

        ALLOCATE(CIRCEXP(LCIRC),  STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'APREF_2; CIRCEXP',LCIRC)
           RETURN
        ENDIF

        WR         = 0.0
        NPROJ      = NUMREF
        LIMITRANGE = (RANGECOS < 1.0)
        NULLIFY(LCG)  ! FOR INTEL COMPILER

C       IF LIMITRANGE, LIST CLOSE REF. IMAGES, RETURNS: NPROJ
        CALL MAKE_CLOSE_LIST(NUMREF,LIMITRANGE,
     &                       REFDIR,EXPDIR,
     &                       RANGECOS, .TRUE., 
     &                       LCG, NPROJ, IRTFLG)

        IF (NPROJ <= 0) THEN
C          NO REF. IMAGE WITHIN COMPARISON ANGLE
           IMGREFL = 0 
           CCROT   = -1.0
           RANGNEW = 0.0
           MIRNEW  = .FALSE.
           RETURN
        ENDIF
    
C       CALCULATE DIMENSIONS FOR NORMALIZING IN APRINGS_ONE
	CNS2 = NX/2+1
	CNR2 = NY/2+1

C       EXTRACT EXP. IMAGE RINGS, NORMALIZE & FFT THEM
        CALL APRINGS_ONE_NEW(NX,NY, CNS2,CNR2, XIM, .FALSE.,
     &                       MODE,NUMR,NRING,LCIRC, WR,FFTW_PLANS,
     &                       CIRCEXP,IRTFLG)
    
        CCROTD  = -1.0D20
        IMGREFL = 0
      
        DO IMIL=1,NPROJ
           IMI = IMIL
           IF (LIMITRANGE) IMI = ABS(LCG(IMIL))

           IF (CKMIRROR .AND. LIMITRANGE) THEN
C              ONLY SEARCH EITHER MIRRORED OR NON-MIRRORED
               USE_UN  = (LCG(IMIL) >= 0)
               USE_MIR = (LCG(IMIL) < 0)
           ELSE
C              SEARCH BOTH MIRRORED & NON-MIRRORED IF CHKMIR
               USE_UN  = .TRUE.
               USE_MIR = CKMIRROR
           ENDIF

C          CHECK MIRRORED/ NON-MIRRORED POSITION
           CALL CROSRNG_2(CIRCREF(1,IMI),CIRCEXP,
     &                    LCIRC,NRING, MAXRIN,NUMR,
     &                    .FALSE.,FFTW_PLANS(1),
     &                    USE_UN,USE_MIR,
     &                    ISMIRRORED,  TOTMIN,TOT)

           IF (TOTMIN >= CCROTD)  THEN
C             GOOD MATCH WITH TOTA (MIRRORED OR NOT)  POSITION 
              CCROTD  = TOTMIN
              RANGNEW = TOT
              MIRNEW  =  ISMIRRORED
              IMGREFL = IMI
           ENDIF

       ENDDO !END OF: DO IMIL=1,NPROJ
          
       IF (MODE == 'F')  THEN
           RANGNEW = (RANGNEW-1.0) / MAXRIN*360.0
       ELSE
           RANGNEW = (RANGNEW-1.0) / MAXRIN*180.0
       ENDIF

       CCROT = CCROTD
       IF (ASSOCIATED(LCG))    DEALLOCATE(LCG)
       NULLIFY(LCG)
       IF (ALLOCATED(CIRCEXP)) DEALLOCATE(CIRCEXP)

       END

#ifdef DEBUGD
        NR     = NUMR(1,NRING)     ! OUTER RING NUMBER
        MAXRIN = NUMR(3,NRING)

        write(88,*) ' DEBUG OUTPUT FROM APREF_PM'
        write(88,902) NR,MAXRIN
902     format(' OUTER RING NUMBER: ',I10,' RING LENGTH:',I10)
        do i = 1, nring
           write(88,901) NUMR(1,I),NUMR(2,I),NUMR(3,I)
901        format(3i10)
        enddo
#endif

