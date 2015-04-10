
C++*********************************************************************
C
C APSHIFT.F         NEW                           NOV 03 ARDEAN LEITH
C                   AVI ERROR IF NOT NORM3        JUN 04 ARDEAN LEITH
C                   GETDATS HAD UNDEFINED NUMTH   FEB 05 ARDEAN LEITH
C                   USED APCC                     FEB 08 ARDEAN LEITH
C                   REWRITE                       NOV 08 ARDEAN LEITH
C                   MIRRORED SHIFT BUG            MAR 11 ARDEAN LEITH
C                   APCC DOES FFTW SCALIING NOW   AUG 11 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2011  Health Research Inc.,                         *
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
C  APSHIFT(LUNIN, REFPAT,IMGREF,  NSAM,NROW, NSAMP,NROWP, 
C          EXPBUF,AVI,SIGI, ISHRANGE,
C          RANGNEW,XSHNEW,YSHNEW,MIRRORNEW,PEAKV,IRTFLG)
C                                  
C  PURPOSE:  MASTER SHIFT DETERMINATION ROUTINE FOR MANY 'AP ..' 
C            OPERATIONS
C
C            ROTATE EXPBUF IMAGE BY RANGNEW, MIRROR IF NECESSARY, PAD 
C            TO 2X SIZE, CROSS-CORRELATE WITH REFERENCE IMAGE, FIND
C            CC PEAK.
C
C  PARAMETERS: LUNIN       IO UNIT                                SENT
C              EXPPAT      EXP. IMAGE TEMPLATE                    SENT
C              IMGEXP      EXP. IMAGE NUMBER                      SENT
C              REFPAT      REF. IMAGE TEMPLATE                    SENT
C              IMGREF      REF. IMAGE NUMBER                      SENT
C              NSAM,NROW   ACTUAL INPUT IMAGE DIMENSIONS          SENT
C              NSAMP,NROWP NSAMP=2*NSAM+2,NROWP=2*NSAM            SENT
C              EXPBUF      EXP. IMAGE BUFFER                      SENT  
C              AVI,SIGI    EXP. IMAGE STATISTICS                  SENT
C              ISHRANGE    POSSIBLE IMAGE SHIFT                   SENT
C              RANGNEW     INPLANE ROTATION ANGLE                 SENT
C              XSHNEW..    SHIFT                                  RET.
C              MIRRORNEW   LOGICAL FLAG THAT REF. NEEDS MIRROR    SENT
C              PEAKV       PEAK HEIGHT                            RET. 
C              IRTFLG      ERROR FLAG                             RET.  
C
C  NOTE:   NSAMP INCLUDES THE EXTRA SPACE FOR FOURIER XFORM
C
C--*********************************************************************

        SUBROUTINE APSHIFT(LUNIN, REFPAT,IMGREF,
     &                     NSAM,NROW, NSAMP,NROWP, 
     &                     EXPBUF,AVI,SIGI, ISHRANGE,
     &                     RANGNEW,XSHNEW,YSHNEW,
     &                     MIRRORNEW,PEAKV,IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

         
        INTEGER,INTENT(IN)    :: LUNIN
        CHARACTER(LEN=*)      :: REFPAT
        INTEGER,INTENT(IN)    :: IMGREF
        INTEGER,INTENT(IN)    :: NSAM,NROW,NSAMP,NROWP
        REAL,INTENT(IN)       :: EXPBUF(NSAM,NROW)
        REAL,INTENT(INOUT)    :: AVI,SIGI
        INTEGER,INTENT(IN)    :: ISHRANGE
        REAL,INTENT(OUT)      :: RANGNEW,XSHNEW,YSHNEW
        LOGICAL,INTENT(IN)    :: MIRRORNEW
        REAL,INTENT(OUT)      :: PEAKV
        INTEGER,INTENT(OUT)   :: IRTFLG

        LOGICAL               :: SKIP_PEAK,NORMIT,SPIDER_SIGN
        CHARACTER(LEN=MAXNAM) :: FILNAM

C       AUTOMATIC ARRAYS
        REAL                  :: BUFPADI(NSAMP,NROWP)   ! NSAMP=2*NSAM+2
        REAL                  :: BUFPADR(NSAMP,NROWP)   ! NROWP=2*NROW

        INTEGER               :: ISAM,ICENT,KCENT,NSAMP1,IROW
        INTEGER               :: IGO,NLET,ITYPE,NSLICE,MAXIM
        INTEGER               :: NE,IODD

        REAL                  :: XI,XOLD,YOLD,RN2,SN2,RW2,RS2
        REAL                  :: COD,SID,SIGR,PADVAL,YI,YCOD
        REAL                  :: YSID,ZDUM

        REAL                  :: quadri

        REAL, PARAMETER       :: PI = 3.1415926535897932384626
        REAL, PARAMETER       :: DGR_TO_RAD = PI/180

#ifdef DEBUG
         write(6,*) '----- in apshift ---'
         write(6,*) 'IMGREF,RANGNEW, MIRRORNEW'
         write(6,*)  IMGREF,RANGNEW, MIRRORNEW 
#endif

C        ROTATE IMAGE BY RANGNEW & PAD TO 2X SIZE

         ICENT = NROW/2+1
         KCENT = NSAM/2+1
         RN2   = -NROW/2
         SN2   = -NSAM/2
         RW2   = -RN2
         RS2   = -SN2

         IF (MOD(NSAM,2) .EQ. 0) RW2 = RW2 - 1.0
         IF (MOD(NROW,2) .EQ. 0) RS2 = RS2 - 1.0

         NSAMP1 = NSAM + 1

C        ROTATE THE IMAGE ------------------------------------- ROTATE

         COD = COS(RANGNEW * DGR_TO_RAD)
         SID = SIN(RANGNEW * DGR_TO_RAD)

         DO IROW=1,NROW
           YI = IROW - ICENT
           IF (YI .LT. RN2)  YI = MIN(RW2+YI-RN2+1.0, RW2)
           IF (YI .GT. RW2)  YI = MAX(RN2+YI-RW2-1.0, RN2)
           YCOD   =  YI * COD + ICENT
           YSID   = -YI * SID + KCENT
           IGO    = (IROW - 1) * NSAM

c$omp      parallel do private(isam,xi,xold,yold)
           DO ISAM=1,NSAM
              XI = ISAM - KCENT 
              IF (XI .LT. SN2) XI = MIN(RS2+XI-SN2+1.0,RS2)
              IF (XI .GT. RS2) XI = MAX(SN2+XI-RS2-1.0,SN2)
              YOLD = XI * SID  + YCOD
              XOLD = XI * COD  + YSID

              IF (.NOT. MIRRORNEW) THEN
C                NO MIRROR 
                 BUFPADI(ISAM,IROW) = 
     &                   QUADRI(XOLD,YOLD,NSAM,NROW,EXPBUF)
              ELSE
C                MIRROR THE IMAGE IN PLACE  
                 BUFPADI(NSAMP1-ISAM,IROW) = 
     &                   QUADRI(XOLD,YOLD,NSAM,NROW,EXPBUF)
              ENDIF
	   ENDDO

C          FILL REMAINING EMPTY COLS IN PADDED IMAGE
           BUFPADI(NSAMP1:NSAMP, IROW) = AVI
	 ENDDO

C        FILL REMAINING EMPTY ROWS IN PADDED IMAGE
         BUFPADI(1:NSAMP, NROW+1:NROWP) = AVI

C        OPEN REFERENCE IMAGE INPUT FILE
         NLET = 0
         CALL FILGET(REFPAT,FILNAM,NLET,IMGREF,IRTFLG)
         IF (IRTFLG .NE. 0)  RETURN

         CALL OPFILEC(0,.FALSE.,FILNAM,LUNIN,'O',ITYPE,
     &                NSAM,NROW,NSLICE,MAXIM,'INPUT',.FALSE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
         IF (IMAMI .NE. 1) 
     &       CALL NORM3(LUNIN,NSAM,NROW,NSLICE,FMAX,FMIN,AV)
         SIGR = SIG

C        LOAD & PAD REFERENCE IMAGE TO DOUBLE SIZE (WITH ZEROS!)
         PADVAL = 0.0
         CALL REDNPADVOL(LUNIN,PADVAL, 
     &                   NSAM,NROW,1,  2*NSAM+2,2*NROW,1,
     &                   BUFPADR, IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         CLOSE(LUNIN)

#ifdef DEBUGNEVER
c-----------------debug
         write(6,*) ' nsamp,nrowp: ', nsamp,nrowp

         call opfilec(0,.false.,'jnkpadi',lunin,'U',itype,
     &                nsamp,nrowp,1,maxim,' ',.false.,irtflg)
         call wrtvol(lunin,nsamp,nrowp,1,1,bufpadi,irtflg)
         close(lunin)
         call opfilec(0,.false.,'jnkpadr',lunin,'U',itype,
     &                nsamp,nrowp,1,maxim,' ',.false.,irtflg)
         call wrtvol(lunin,nsamp,nrowp,1,1,bufpadr,irtflg)
         close(lunin)
c---------------------------
#endif


C        CROSS CORRELATION BUFPADI vs BUFPADR ------------------- CC
         SKIP_PEAK   = .FALSE.
         NORMIT      = .TRUE.
         SPIDER_SIGN = .FALSE.

         CALL APCC(NSAMP, 2*NSAM,2*NROW,1, BUFPADI,BUFPADR,
     &             SKIP_PEAK,NORMIT,SPIDER_SIGN, 
     &             ISHRANGE,ISHRANGE,0,
     &             XSHNEW,YSHNEW,ZDUM, PEAKV,IRTFLG)

        IF (IRTFLG .NE. 0)  THEN
           CALL ERRT(101,'APSHIFT CC ERROR',NE)
           RETURN
        ENDIF

C       NORMALIZATION (PROPERLY WOULD NEED SIG AFTER PADDING!) al

        PEAKV = PEAKV / FLOAT(NSAMP*NROWP-1) / SIGI / SIGR

C       ADJUST FOR MIRROR OF IMAGE NOT MIRROR OF REF. IMAGE
        IF (MIRRORNEW) THEN              ! BUG FIX MAR 2011
           IODD = MOD(NSAM,2)
           IF (IODD == 0) THEN
              XSHNEW = 1 - XSHNEW
           ELSE
              XSHNEW = -XSHNEW
           ENDIF
        ENDIF            
        
        END

#ifdef NEVER
        write(6,*) 'peakv0:',peakv,NSAMP*NROWP,sigi,sigr
        write(6,*) 'peakv1:',peakv
        write(6,*) 'peakv1/npixm1:',peakv/FLOAT(NSAMP*NROWP-1)
        write(6,*) 'peakv2:',peakv
        write(6,*) '1/peakv2:',1.0/peakv
#endif

