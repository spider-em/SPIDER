
C++*********************************************************************
C
C    AP_GETDATA.F     EXTRACTED                 OCT. 2003 ARDEAN LEITH
C                     LS1... INTERNALIZED       DEC. 2003 ARDEAN LEITH
C                     CALLS NORM3               JUN  2004 ARDEAN LEITH
C                     REMOVED WINDOWING         MAR  2005 ARDEAN LEITH
C                     AVI & SIGI PARAMETERS     NOV  2008 ARDEAN LEITH
C                     ADDED GETDATS_RTSQ        FEB  2011 ARDEAN LEITH
C                     MERGED WITH AP_GETDAT     JUN  2011 ARDEAN LEITH
C                     PAD                       AUG  2011 ARDEAN LEITH
C                     AP_GETDATA_MASK           NOV  2011 ARDEAN LEITH
C                     AP_GETDATA_MASK_RTSQ      DEC  2011 ARDEAN LEITH
C                     RTSQ CALL                 DEC  2011 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C  AP_GETDATA           (ILIST,NUMIMG, NX,NY,NXP,NYP,
C             NUMTH,EXPPAT,LUNIN,IGO,IEND,
C             MPIBCAST,BUFOUT,
C             WANTSTATS, AVI,SIGI, IRTFLG)
C
C  AP_GETDATA_RTSQ      (ILIST,NUMIMG, NX,NY, NXP,NYP, PADVAL,
C             NUMTH,EXPPAT,LUNIN, IGO,IEND,
C             ANGINHEADER,ANGEXP, 
C             MPIBCAST,BUFEXP,BUFOUT, 
C             WANTSTATS,AVI,SIGI,IRTFLG)
C
C  AP_GETDATA_MASK     (ILIST,NUMIMG, NX,NY, NXP,NYP, RADI,
C             NUMTH,EXPPAT,LUNIN, IGO,IEND,
C             MPIBCAST, BUFOUT, 
C             WANTSTATS,AVI,SIGI,  IRTFLG)
C
C  AP_GETDATA_MASK_RTSQ(ILIST,NUMIMG, NX,NY, NXP,NYP, RADI,
C             NUMTH,EXPPAT,LUNIN, IGO,IEND,
C             ANGINHEADER,ANGEXP, 
C             MPIBCAST,BUFEXP,BUFOUT,
C             WANTSTATS,AVI,SIGI,IRTFLG)
C
C  AP_GETDATA_DEN( ILIST,NUMIMG, NXT,NX,NY, N2XLD,N2X,N2Y, RADI,
C             NUMTH,EXPPAT,LUNIN, IGO,IEND,
C             ANGINHEADER,ANGEXP, MPIBCAST,
C             BUFEXP,BUFOUT,
C             FBS_WANTED,IRTFLG)
C
C PURPOSE:  RETURN SERIES OF IMAGE DATA IN ARRAY BUFOUT FOR 'AP' OPS.
C           CALLS NORM3 IF NEEDED. CAN LOAD MUTIPLE IMAGES.
C           RETURNS ARRAY OF IMAGE STATISTICS
C
C PARAMETERS:
C       ILIST               LIST OF IMAGE FILE NUMBERS        (INPUT)
C       NUMIMG              NO. OF IMAGES                     (INPUT)
C       NX,NY               IMAGE DIMENSIONS                  (INPUT)
C       NXP,NYP             PADDED IMAGE DIMENSIONS           (INPUT)
C       PADVAL              PAD VALUE                         (INPUT)
C       RADI                MASK RADIUS( OUTSIDE USE PADVAL)  (INPUT)
C       NUMTH               # THREADS                         (INPUT)
C       EXPPAT              IMAGE SERIES FILE TEMPLATE        (INPUT)
C       LUNIN               IMAGE FILE IO UNIT                (INPUT)
C       IGO,IEND            IMAGE INDEX RANGE                 (INPUT)
C       ANGINHEADER         ANGLES IN IMAGE HEADER            (INPUT)
C       MPIBCAST            USE MPI BCAST READ                (INPUT)
C       ANGEXP              ANGLES, ETC                       (IN/OUT)
C       BUFEXP              UNROTATED IMAGE ARRAY             (OUTPUT)
C       BUFOUT              FINAL IMAGE ARRAY                 (OUTPUT)
C       WANTSTATS           WANT STATISTICS                   (INPUT)
C       AVI,SIGI            STATISTICS ARRAYS                 (OUTPUT)
C       IRTFLG              ERROR FLAG                        (OUTPUT)
C
C--*********************************************************************

C       --------------------- AP_GETDATA ------------------------------

	SUBROUTINE AP_GETDATA(ILIST,NUMIMG, 
     &                       NX,NY, NXP,NYP, PADVAL,
     &                       NUMTH,EXPPAT,LUNIN, IGO,IEND,
     &                       MPIBCAST, BUFOUT, 
     &                       WANTSTATS,AVI,SIGI,  IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        INTEGER, INTENT(IN)          :: ILIST(NUMIMG)
        INTEGER, INTENT(IN)          :: NUMIMG
        INTEGER, INTENT(IN)          :: NX,NY,NUMTH, NXP,NYP
	REAL                         :: PADVAL
        CHARACTER (LEN=*),INTENT(IN) :: EXPPAT
        INTEGER, INTENT(IN)          :: LUNIN,IGO,IEND
	REAL                         :: BUFOUT(NXP,NYP,NUMTH)
        LOGICAL                      :: MPIBCAST
        LOGICAL                      :: WANTSTATS
	REAL,    INTENT(OUT)         :: AVI(NUMTH),SIGI(NUMTH)
        INTEGER, INTENT(OUT)         :: IRTFLG

        CHARACTER(LEN=MAXNAM)        :: FILNAM
        LOGICAL                      :: ONEIMAGE,PADIT
        INTEGER                      :: ITI,NLET,MAXIM,LX,LY
        INTEGER                      :: IFORMT,LZ,IT
        DOUBLE PRECISION             :: DAV,DSIG
        REAL                         :: UNUSED

        INTEGER                      :: INV
        INTEGER                      :: lnblnk
 
        !write(6,*) 'Reading images: ',igo,'...',iend,exppat(1:20)
        NLET   = lnblnk(EXPPAT)

        ONEIMAGE = (IGO <= 0 .OR. ILIST(IGO) <= 0)
        IF (ONEIMAGE) THEN
           FILNAM = EXPPAT
        ENDIF

        PADIT  = (NXP > NX .OR. NYP > NY)

        DO ITI=IGO,IEND

C          BUFOUT STARTING LOCATION
           IT = ITI - IGO + 1

           IF (.NOT. ONEIMAGE) THEN
              NLET = 0
              CALL FILGET(EXPPAT,FILNAM,NLET,ILIST(ITI),IRTFLG)
              IF (IRTFLG .NE. 0)  RETURN
           ENDIF

           !write(6,*) 'Reading image: ',filnam(1:30),' ',oneimage,igo

C          OPEN EXISTING IMAGE FILE
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FILNAM,LUNIN,'O',IFORMT,
     &               LX,LY,LZ,MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           CALL SIZCHK(UNUSED,LX,LY,  0,0,
     &                        NX, NY, 0,0,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
          
C          LOAD THE WHOLE IMAGE
           IF (PADIT) THEN
              CALL REDNPADVOL_SEL(LUNIN,PADVAL, 
     &                        NX, NY, 1,  NXP,NYP,1,MPIBCAST,
     &                        BUFOUT(1,1,IT), IRTFLG)
           ELSE
              CALL REDVOL_SEL(LUNIN,NX,NY,1,1,MPIBCAST, 
     &                        BUFOUT(1,1,IT), IRTFLG)
           ENDIF

           CLOSE(LUNIN)
           IF (IRTFLG .NE. 0) RETURN

           IF (WANTSTATS) THEN
              IF (IMAMI .NE. 1) THEN
C                FIND IMAGE STATISTICS, USE OPEN MP
                 CALL NORMVALSP(BUFOUT(1,1,IT), NX,NY,1, 
     &                          NXP,NYP,1,
     &                          DAV,DSIG,.TRUE.)
                 AV  = DAV     ! INTO COMMON
                 SIG = DSIG
              ENDIF
C             RECORD THE AVERAGE AND SD (FROM COMMON) FOR RETURN
              AVI(IT)  = AV
              SIGI(IT) = SIG
           ENDIF
        ENDDO

        IRTFLG = 0
        END

C       ************************** AP_GETDATA_RTSQ *******************
          
	SUBROUTINE AP_GETDATA_RTSQ(ILIST,NUMIMG, 
     &                          NX,NY, NXP,NYP, PADVAL,
     &                          NUMTH,EXPPAT,LUNIN, IGO,IEND,
     &                          ANGINHEADER,ANGEXP, 
     &                          MPIBCAST,BUFEXP,BUFOUT,
     &                          WANTSTATS,AVI,SIGI,FBS_WANTED,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        INTEGER, INTENT(IN)   :: ILIST(NUMIMG)
        INTEGER, INTENT(IN)   :: NUMIMG
        INTEGER, INTENT(IN)   :: NX,NY,NUMTH,NXP,NYP
        CHARACTER (LEN=*)     :: EXPPAT
        INTEGER, INTENT(IN)   :: LUNIN, IGO,IEND
        LOGICAL, INTENT(IN)   :: ANGINHEADER
        REAL,    INTENT(IN)   :: ANGEXP(8,NUMIMG) 
        LOGICAL, INTENT(IN)   :: MPIBCAST
        REAL                  :: BUFEXP(NX,NY), PADVAL
	REAL,    INTENT(OUT)  :: BUFOUT(NXP,NYP,NUMTH)
        LOGICAL, INTENT(IN)   :: WANTSTATS
	REAL,    INTENT(OUT)  :: AVI(NUMTH),SIGI(NUMTH)
        LOGICAL, INTENT(IN)   :: FBS_WANTED
        INTEGER, INTENT(OUT)  :: IRTFLG

        INTEGER               :: ITI,NLET,MAXIM,IFORMT,LX,LY,LZ,IDUM,IT
        REAL                  :: SCALE,UNUSED,SHXI,SHYI,THETA
        CHARACTER(LEN=MAXNAM) :: FILNAM
        LOGICAL               :: ONEIMAGE,PADIT
        DOUBLE PRECISION      :: DAV,DSIG

        SCALE    = 1.0
        ONEIMAGE = (ILIST(IGO) .LE. 0)
        PADIT     = (NXP > NX .OR. NYP > NY)
        IF (PADIT) THEN
            CALL ERRT(101,'PGM ERROR, PADDING NOT IMPLEMENTED',IDUM)
            STOP
        ENDIF
C       write(6,*) 'Reading images: ',igo,'...',iend

        DO ITI=IGO,IEND
C          BUFOUT STARTING LOCATION
           IT = ITI-IGO+1

           IF (ONEIMAGE) THEN
              FILNAM = EXPPAT
           ELSE
              NLET = 0
              CALL FILGET(EXPPAT,FILNAM,NLET,ILIST(ITI),IRTFLG)
              IF (IRTFLG .NE. 0)  RETURN
           ENDIF

C          OPEN EXISTING IMAGE FILE
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FILNAM,LUNIN,'O',IFORMT,
     &                  LX,LY,LZ, MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0)  RETURN

           CALL SIZCHK(UNUSED,LX,LY,0,0,
     &                        NX, NY, 0,0,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (ANGINHEADER) THEN
C             GET ANGLES FROM IMAGE HEADER
              CALL LUNGETVALS(LUNIN,IAPLOC + 1,8,ANGEXP(1,ITI),IRTFLG)
              IF (IRTFLG .NE. 0) RETURN
           ENDIF

C          READ IN THE WHOLE IMAGE INTO BUFEXP
           CALL REDVOL_SEL(LUNIN,NX,NY,1,1, MPIBCAST,BUFEXP,IRTFLG)

           CLOSE(LUNIN)
           IF (IRTFLG .NE. 0) RETURN

C          SHIFT & ROTATE THE INPUT IMAGE --> BUFOUT  
           THETA = ANGEXP(4,ITI)
           SHXI  = ANGEXP(5,ITI)
           SHYI  = ANGEXP(6,ITI)
           !print *,'rtsq; thta,xy',THETA,SHXI,SHYI

C          ROTATE AND SHIFT
           IF (FBS_WANTED) THEN
              CALL RTSF_PAD(BUFEXP,BUFOUT(1,1,IT),
     &                      NX,NY, NXP,NYP,
     &                      THETA,SCALE,SHXI,SHYI, IRTFLG)
           ELSE
              CALL RTSQ(BUFEXP,BUFOUT(1,1,IT),
     &                  NX,NY, NXP,NYP,
     &                  THETA,SCALE,SHXI,SHYI, IRTFLG)
           ENDIF

           IF (WANTSTATS) THEN
C             FIND IMAGE STATISTICS, (HAS MOVED SO NOT SAME)
              CALL NORMVALSP(BUFOUT(1,1,IT),NX,NY,1, 
     &                       NXP,NYP,1,
     &                       DAV,DSIG,.TRUE.)

C             RECORD THE AVERAGE AND SD 
              AVI(IT)  = DAV
              SIGI(IT) = DSIG
           ENDIF
        ENDDO

        IRTFLG = 0
        END

C       ************************** AP_GETDATA_MASK ***************

	SUBROUTINE AP_GETDATA_MASK(ILIST,NUMIMG, 
     &                       NX,NY, NXP,NYP, RADI,
     &                       NUMTH,EXPPAT,LUNIN, IGO,IEND,
     &                       MPIBCAST, BUFOUT, 
     &                       WANTSTATS,AVI,SIGI,  IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        INTEGER, INTENT(IN)          :: ILIST(NUMIMG)
        INTEGER, INTENT(IN)          :: NUMIMG
        INTEGER, INTENT(IN)          :: NX,NY,NUMTH, NXP,NYP
	REAL,    INTENT(IN)          :: RADI
        CHARACTER (LEN=*),INTENT(IN) :: EXPPAT
        INTEGER, INTENT(IN)          :: LUNIN,IGO,IEND
	REAL                         :: BUFOUT(NXP,NYP,NUMTH)
        LOGICAL                      :: MPIBCAST
        LOGICAL                      :: WANTSTATS
	REAL,    INTENT(OUT)         :: AVI(NUMTH),SIGI(NUMTH)
        INTEGER, INTENT(OUT)         :: IRTFLG

        CHARACTER(LEN=MAXNAM)        :: FILNAM
        LOGICAL                      :: ONEIMAGE,PADIT
        INTEGER                      :: ITI,NLET,MAXIM,LX,LY
        INTEGER                      :: IFORMT,LZ,IT
        DOUBLE PRECISION             :: DAV,DSIG
        REAL                         :: UNUSED

        INTEGER                      :: INV
        REAL                         :: XCEN,YCEN,RADISQ,RSQ,YSQTMP
        INTEGER                      :: IX,IY,IZ
 
        integer :: itype
        real :: avtmp

c       write(6,*) 'Reading images: ',igo,'...',iend

        ONEIMAGE = (IGO .LE. 0 .OR. ILIST(IGO) .LE. 0)
        IF (ONEIMAGE) THEN
           FILNAM = EXPPAT
        ENDIF

        PADIT  = (NXP > NX .OR. NYP > NY)

        IF (RADI > 0) THEN
           RADISQ = RADI **2

           XCEN   = NX/2 + 1
           YCEN   = NY/2 + 1
        ENDIF
 
        DO ITI=IGO,IEND

C          BUFOUT STARTING LOCATION
           IT = ITI - IGO + 1

           IF (.NOT. ONEIMAGE) THEN
              NLET = 0
              CALL FILGET(EXPPAT,FILNAM,NLET,ILIST(ITI),IRTFLG)
              IF (IRTFLG .NE. 0)  RETURN
           ENDIF

C          OPEN EXISTING IMAGE FILE
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FILNAM,LUNIN,'O',IFORMT,
     &                  LX,LY,LZ,MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           CALL SIZCHK(UNUSED,LX,LY,  0,0,
     &                        NX,NY,  0,0,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
          
           IF (WANTSTATS) THEN
              IF (IMAMI .NE. 1) THEN
C                FIND IMAGE STATISTICS
                 CALL NORM3(LUNIN,LX,LY,LZ, FMAX,FMIN,AV)
              ENDIF

C             RECORD THE AVERAGE AND SD (FROM COMMON)
              AVI(IT)  = AV
              SIGI(IT) = SIG ! IN COMMON
           ENDIF

C          LOAD THE WHOLE IMAGE
           IF (PADIT) THEN
              !write(6,*) ' padit; radi,av:',radi,av
              CALL REDMASKNPADVOL(LUNIN,AV,RADI, 
     &                            NX, NY, 1, NXP,NYP,1,MPIBCAST,
     &                            BUFOUT(1,1,IT),IRTFLG)
           ELSE
C             NOT PADDED
              CALL REDVOL_SEL(LUNIN,NX,NY,1,1,MPIBCAST, 
     &                        BUFOUT(1,1,IT),IRTFLG)
           ENDIF

#ifdef NEVER
           !write(6,*) ' radi,av:',radi,av,bufout(1,1,it)
           maxim = 0
           itype = 1
           call opfilec(0,.false.,'jnk0',98,'U',itype,
     &               nxp,nyp,1,maxim,' ',.true.,irtflg)
           call wrtvol(98,nxp,nyp, 1,1, bufout,irtflg)
           close(98)
#endif

           CLOSE(LUNIN)
           IF (IRTFLG .NE. 0) RETURN

           IF (.NOT. PADIT .AND. RADI > 0 .AND. WANTSTATS) THEN
C              MASK OUTSIDE CIRCLE USING AV
               AVTMP = AVI(IT)

               DO IY = 1,NY
                  YSQTMP = (FLOAT(IY) - YCEN) **2 
 
                  DO IX = 1,NX
                     RSQ = (FLOAT(IX) - XCEN) **2 + YSQTMP  
                     IF (RSQ > RADISQ) BUFOUT(IX,IY,IT) = AVTMP 
                  ENDDO
               ENDDO
            ENDIF

#ifdef NEVER
           maxim = 0
           itype = 1
           call opfilec(0,.false.,'jnk1',98,'U',itype,
     &               nxp,nyp,1,maxim,' ',.true.,irtflg)
           call wrtvol(98,nxp,nyp, 1,1, bufout,irtflg)
           close(98)
#endif

        ENDDO

        IRTFLG = 0
        END


C       ************************** AP_GETDATA_MASK_RTSQ ***************
    
	SUBROUTINE AP_GETDATA_MASK_RTSQ(ILIST,NUMIMG, 
     &                          NX,NY, NXP,NYP, RADI,
     &                          NUMTH,EXPPAT,LUNIN, IGO,IEND,
     &                          ANGINHEADER,ANGEXP, 
     &                          MPIBCAST,BUFEXP,BUFOUT,
     &                          WANTSTATS,AVI,SIGI,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        INTEGER, INTENT(IN)   :: ILIST(NUMIMG)
        INTEGER, INTENT(IN)   :: NUMIMG
        INTEGER, INTENT(IN)   :: NX,NY,NUMTH,NXP,NYP
        REAL,    INTENT(IN)   :: RADI
        CHARACTER (LEN=*)     :: EXPPAT
        INTEGER, INTENT(IN)   :: LUNIN, IGO,IEND
        LOGICAL, INTENT(IN)   :: ANGINHEADER
        REAL,    INTENT(IN)   :: ANGEXP(8,NUMIMG) 
        LOGICAL, INTENT(IN)   :: MPIBCAST
        REAL                  :: BUFEXP(NX,NY)
	REAL,    INTENT(OUT)  :: BUFOUT(NXP,NYP,NUMTH)
        LOGICAL, INTENT(IN)   :: WANTSTATS
	REAL,    INTENT(OUT)  :: AVI(NUMTH),SIGI(NUMTH)
        INTEGER, INTENT(OUT)  :: IRTFLG

        INTEGER               :: ITI,NLET,MAXIM,IFORMT,LX,LY
        INTEGER               :: LZ,IT,IDUM
        REAL                  :: UNUSED,SHXI,SHYI,THETA,YSQTMP,AVTMP
        CHARACTER(LEN=MAXNAM) :: FILNAM
        LOGICAL               :: ONEIMAGE,PADIT
        DOUBLE PRECISION      :: DAV,DSIG
        REAL                  :: XCEN,YCEN,RADISQ,RSQ
        INTEGER               :: IX,IY,IZ

        REAL, PARAMETER       :: SCALE = 1.0

        ONEIMAGE = (ILIST(IGO) <= 0)
        IF (ONEIMAGE) THEN
           FILNAM = EXPPAT
        ENDIF
        PADIT    = (NXP > NX .OR. NYP > NY)

        IF (RADI > 0) THEN
           RADISQ = RADI **2

           XCEN   = NX/2 + 1
           YCEN   = NY/2 + 1
        ENDIF
  
C       write(6,*) 'Reading images: ',igo,'...',iend

        DO ITI=IGO,IEND

C          BUFOUT STARTING LOCATION
           IT = ITI - IGO + 1

           IF (.NOT. ONEIMAGE) THEN
              NLET = 0
              CALL FILGET(EXPPAT,FILNAM,NLET,ILIST(ITI),IRTFLG)
              IF (IRTFLG .NE. 0)  RETURN
           ENDIF

C          OPEN EXISTING IMAGE FILE
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FILNAM,LUNIN,'O',IFORMT,
     &                  LX,LY,LZ, MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0)  RETURN

           CALL SIZCHK(UNUSED, LX,LY,0,0,  NX,NY,0,0, IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (ANGINHEADER) THEN
C             GET ANGLES FROM IMAGE HEADER
              CALL LUNGETVALS(LUNIN,IAPLOC + 1,8,ANGEXP(1,ITI),IRTFLG)
              IF (IRTFLG .NE. 0) RETURN
           ENDIF

C          READ THE WHOLE IMAGE INTO: BUFEXP
           CALL REDVOL_SEL(LUNIN,NX,NY,1,1, MPIBCAST, BUFEXP,IRTFLG)

           CLOSE(LUNIN)
           IF (IRTFLG .NE. 0) RETURN

C          SHIFT & ROTATE THE INPUT IMAGE --> BUFOUT  
           THETA = ANGEXP(4,ITI)
           SHXI  = ANGEXP(5,ITI)
           SHYI  = ANGEXP(6,ITI)
           !print *,'rqs; thta,xy',THETA,SHXI,SHYI

C          ROTATE AND SHIFT
           CALL RTSQ(BUFEXP,BUFOUT(1,1,IT),
     &               NX,NY,NXP,NYP,
     &               THETA,SCALE,SHXI,SHYI,IDUM,0)

           IF (WANTSTATS) THEN
C             FIND IMAGE STATISTICS, (HAS MOVED SO NOT SAME)
              CALL NORMVALSP(BUFOUT(1,1,IT),NX,NY,1, 
     &                       NXP,NYP,1,
     &                       DAV,DSIG,.TRUE.)

C             RECORD THE AVERAGE AND SD 
              AVI(IT)  = DAV
              SIGI(IT) = DSIG

              IF (RADI > 0 ) THEN
C                MASK OUTSIDE CIRCLE USING: AV

                 AVTMP = AVI(IT)
                 DO IY = 1,NY
                    YSQTMP = (FLOAT(IY) - YCEN) **2 
 
                    DO IX = 1,NX
                       RSQ = (FLOAT(IX) - XCEN) **2 + YSQTMP
                       IF (RSQ > RADISQ) BUFEXP(IX,IY) = AVTMP 
                    ENDDO 
                 ENDDO     
              ENDIF
           ENDIF                
        ENDDO

        IRTFLG = 0
        END

          
C       ************************** AP_GETDATA_DEN *******************
          
C       INPUT: 2XFFT PAD, OUTPUT: UNPADDED    
	SUBROUTINE AP_GETDATA_DEN(ILIST,NUMIMG, 
     &                          NXT,NX,NY, N2XLD,N2X,N2Y, RADI,
     &                          NUMTH,EXPPAT,LUNIN, IGO,IEND,
     &                          ANGINHEADER,ANGEXP, MPIBCAST,
     &                          BUFEXP,BUFOUT,
     &                          FBS_WANTED,IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        INTEGER, INTENT(IN)          :: ILIST(NUMIMG)
        INTEGER, INTENT(IN)          :: NUMIMG
        INTEGER, INTENT(IN)          :: NXT,NX,NY, N2XLD,N2X,N2Y
	REAL,    INTENT(IN)          :: RADI
        INTEGER, INTENT(IN)          :: NUMTH
        CHARACTER (LEN=*),INTENT(IN) :: EXPPAT
        INTEGER, INTENT(IN)          :: LUNIN,IGO,IEND
        LOGICAL, INTENT(IN)          :: ANGINHEADER
        REAL,    INTENT(IN)          :: ANGEXP(8,NUMIMG) 
        LOGICAL                      :: MPIBCAST
	REAL                         :: BUFEXP(N2XLD,N2Y) ! 2XFFT PAD
	REAL                         :: BUFOUT(NX,NY,NUMTH)
        LOGICAL, INTENT(IN)          :: FBS_WANTED
        INTEGER, INTENT(OUT)         :: IRTFLG

        CHARACTER(LEN=MAXNAM)        :: FILNAM
        LOGICAL                      :: ONEIMAGE,PADIT,MASKIT
        INTEGER                      :: ITI,NLET,MAXIM,LX,LY
        INTEGER                      :: IFORMT,LZ,IT,IDUM
        REAL                         :: UNUSED
        REAL                         :: BFPS(4),FDUM
        REAL                         :: THETA,SHXI,SHYI,SCALE
        INTEGER                      :: INV,IOPT
        INTEGER                      :: IX,IY,IZ

	REAL                         :: BUFTMP(NXT,NY)  ! UNPADDED
 
        !write(6,*) 'Reading images: ',igo,'...',iend,fbs_wanted

        ONEIMAGE = (IGO .LE. 0 .OR. ILIST(IGO) .LE. 0)
        IF (ONEIMAGE) FILNAM = EXPPAT

        MASKIT  = (RADI > 0)
        IOPT    = 7          ! BUTTERWORTH LOW PASS FILTER
        BFPS(1) = 0.25       ! PASS-BAND FREQUENCY
        BFPS(2) = 0.25 +0.15 ! STOP BAND FREQUENCY
        BFPS(3) = 0.0 
        BFPS(4) = 0.0 
        SCALE   = 1.0

        DO ITI=IGO,IEND

C          BUFOUT STARTING LOCATION
           IT = ITI - IGO + 1

           IF (.NOT. ONEIMAGE) THEN
              NLET = 0
              CALL FILGET(EXPPAT,FILNAM,NLET,ILIST(ITI),IRTFLG)
              IF (IRTFLG .NE. 0)  RETURN
           ENDIF

C          OPEN EXISTING IMAGE FILE
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FILNAM,LUNIN,'O',IFORMT,
     &                  LX,LY,LZ,MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           CALL SIZCHK(UNUSED,LX,LY,  0,0,
     &                        NX,NY,  0,0,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
          
           IF (ANGINHEADER) THEN
C             GET ANGLES FROM IMAGE HEADER
              CALL LUNGETVALS(LUNIN,IAPLOC + 1,8,ANGEXP(1,ITI),IRTFLG)
              IF (IRTFLG .NE. 0) RETURN
           ENDIF

           IF (MASKIT) THEN
C             PUTPUT: CIRCULAR MASKED, 2XFFT PADDED IMAGE
              !write(6,*) ' padit; radi,av:',radi,av
              CALL REDMASKNPADVOL(LUNIN,AV,RADI, 
     &                            NX, NY, 1, N2XLD,N2Y,1,MPIBCAST,
     &                            BUFEXP,IRTFLG)
           ELSE
C             OUTPUT: 2XFFT PADDED IMAGE
              CALL REDNPADVOL(LUNIN,AV, 
     &                        NX, NY, 1, N2XLD,N2Y,1,
     &                        BUFEXP,IRTFLG)
           ENDIF
           CLOSE(LUNIN)
           IF (IRTFLG .NE. 0) RETURN

           !write(6,*) ' radi,av:',radi,av,bufout(1,1,it)
           !call chkpadfile('jnkld',66,1, N2XLD,n2y,1, nx,ny,1,bufexp,irtflg)

C          FOURIER LOWPASS BUTTERWORTH, 2XFFT PAD IN/OUT
           CALL FQ_BUF(IOPT,BFPS,FDUM,FDUM,FDUM,
     &                 BUFEXP, N2XLD,N2X,N2Y, NX,NY, IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           !call chkpadfile('jnkfq',66,1, N2XLD,n2y,1, nx,ny,1, bufexp,irtflg)
           !call chkfile('jnktmp',66,1,nx,ny,1, buftmp,irtflg)

C          SHIFT & ROTATE THE INPUT IMAGE --> BUFOUT  
           THETA = ANGEXP(4,ITI)
           SHXI  = ANGEXP(5,ITI)
           SHYI  = ANGEXP(6,ITI)
           !print *,'theta,shifts',theta,shxi,shyi,fbs_wanted

C          ROTATE AND SHIFT, FROM: 2xFFT PADDED  TO: UNPADDED
           IF (FBS_WANTED) THEN
              CALL RTSF_PADIN(BUFEXP,BUFOUT(1,1,IT),
     &                        NXT,NX,NY, N2XLD, 
     &                        THETA,SCALE, SHXI,SHYI, IRTFLG)
           ELSE
              CALL RTSQ_PADIN(BUFEXP,BUFOUT(1,1,IT), 
     &                        N2XLD,N2Y, NX,NY,
     &                        THETA,SCALE, SHXI,SHYI, IRTFLG)
           ENDIF

           !call chkfile('jnkout',66,1,nx,ny,1, BUFOUT(1,1,it),irtflg)

        ENDDO

        IRTFLG = 0

        END

