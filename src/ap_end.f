
C++*********************************************************************
C
C    AP_END
C           PEAKV REMOVED                           MAR 04 ArDean Leith
C           IKEY                                    MAR 04 ArDean Leith
C           ANGDIF IF CHANGE IN MIRROR              AUG 04 ArDean Leith
C           ANGLES IN REG OUTPUT WRONG              AUG 04 ArDean Leith
C           APEND2                                  SEP 04 Chao Yang
C           REWRITE FOR END OF 'AP MD'..            SEP 04 ArDean Leith
C           REWRITE FOR UN-MIRROR                   OCT 04 ArDean Leith
C           RANGOUT <= 360                          DEC 04 ArDean Leith
C           MPI WRITES TO REGISTERS NOW             OCT 08 ArDean Leith
C           REMOVED APSHIFT CALL                    NOV 08 ArDean Leith
c           IF (NPROJ.EQ.0)NPROJ=1                  AUG 09 ArDean Leith
C           IF(ABS(ANGDIF).LT.ALMOSTZERO)           AUG 09 ArDean Leith
C           REMOVED OBSLT                           OCT 10 ArDean Leith
C           AP_ENDS FOR REPORTING RAW SHIFTS        MAR 12 ArDean Leith
C           ABS(PEAKV)                              MAY 12 ArDean Leith
C           KEEP OLD PARAM ALWAYS                   JUN 12 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C   AP_END(IMGEXP,IMGREF, 
C          ANGREF,REFDIR, ANGEXP,EXPDIR,ISHRANGE,
C          GOTREFANG, NGOTPAR,CCROTT,PEAKVT,
C          RANGNEW,MIRRORNEW,REFPAT,
C          NPROJT,LUNDOC,DLIST)
C   
C PURPOSE : WRITE ALIGNMENT PARAMETERS TO OUTPUT FILE
C
C PARAMETERS:
C       KEY               IMAGE KEY NUMBER                 (INPUT)
C       IMGEXP            EXP IMAGE NUMBER                 (INPUT)
C       IMGREF            REF IMAGE NUMBER                 (INPUT)
C       ANGREF            REF IMAGE DIRECTION              (INPUT)
C       REFDIR            REF IMAGE DIRECTION              (INPUT)
C       ANGEXP            EXP IMAGE DIRECTION              (INPUT)
C       EXPDIR            EXP IMAGE DIRECTION              (INPUT)
C       ISHRANGE          SHIFT RANGE                      (INPUT)
C       GOTREFANG         FLAG FOR REF. ANGLES AVAIL.      (INPUT)
C       NGOTPAR           NUMBER OF PARAMETERS IN DOC      (INPUT)
C       CCROTT            ROTATIONAL CC VALUE              (INPUT)
C       PEAKVT            FINAL CC VALUE                   (INPUT)            
C       RANGNEW           ROTATION                         (INPUT)
C       XSHNEW,YSHNEW     SHIFT VALUES                     (INPUT)
C       MIRRORNEW         MIRRORED FLAG                    (INPUT)
C       REFPAT            REF. IMAGE SERIES FILE TEMPLATE  (INPUT)
C       NPROJT            # OF PROJECTIONS                 (INPUT)
C       CTYPE             'SH', 'FOU', or 'REF' FLAG       (INPUT)
C       LUNDOC            DOC. FILE OUTPUT UNIT            (INPUT)  
C       CHNG_ORDER        CHANGE ORDER OF SHIFT/ROTATE     (INPUT)
C       SAY_RAW           SAY RAW SHIFTS                   (INPUT)
C       DLIST             PARAMETERS                       (OUTPUT)
C
C  OPERATIONS:  'AP REF', 'AP SH'
C
C--*********************************************************************

         SUBROUTINE AP_ENDS(KEY,IMGEXP,IMGREF, 
     &                  ANGREF,REFDIR, 
     &                  ANGEXP,EXPDIR,ISHRANGE,
     &                  GOTREFANG, NGOTPAR, CCROTT,PEAKVT,
     &                  RANGNEW,XSHRAW,YSHRAW, MIRRORNEW,REFPAT,
     &                  NPROJT, CTYPE, LUNDOC,
     &                  CHNG_ORDER,SAY_RAW,DLIST)

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        INTEGER             :: KEY,IMGEXP,IMGREF
	REAL                :: ANGREF(3),REFDIR(3)
	REAL                :: ANGEXP(7),EXPDIR(3)
        INTEGER             :: ISHRANGE
	LOGICAL             :: GOTREFANG
        INTEGER             :: NGOTPAR
	REAL                :: CCROTT,PEAKVT,RANGNEW,XSHRAW,YSHRAW
        LOGICAL             :: MIRRORNEW
        CHARACTER (LEN=*)   :: REFPAT 
        INTEGER             :: NPROJT
        CHARACTER (LEN=*)   :: CTYPE 
        INTEGER             :: LUNDOC
        REAL                :: DLIST(*)
        LOGICAL             :: CHNG_ORDER,SAY_RAW

        LOGICAL             :: MIRROROLD
	REAL                :: ANGOUT(3)

        REAL, PARAMETER     :: QUADPI = 3.14159265358979323846
        REAL, PARAMETER     :: DGR_TO_RAD = (QUADPI/180)
        REAL, PARAMETER     :: ALMOST_ZERO = 0.05

         NPROJ   = NPROJT
c         IF (NPROJ .EQ. 0) NPROJ = 1 removed aug 09 al

         PEAKV = PEAKVT

         IF (IMGREF > 0) THEN
C           IMGREF IS NUMBER OF MOST SIMILAR REF. IMAGE 
            CCROT   = CCROTT
         ELSE
C           NO REFERENCE IMAGE SELECTED
            CCROT = -1.0
            PEAKV = 0.0
         ENDIF

C        SET NEW PROJECTION ANGLES
         ANGOUT = 0.0            ! DEFAULT VALUE
         IF (IMGREF > 0 .AND. GOTREFANG) THEN
C           USE REF. ANGLES AS NEW PROJECTION ANGLES
            ANGOUT = ANGREF(1:3)

            IF (MIRRORNEW) THEN
C              REF. PROJECTION MUST BE MIRRORED
               ANGOUT(1) = -ANGOUT(1)
               ANGOUT(2) = 180+ANGOUT(2)
            ENDIF
         ELSEIF (NGOTPAR >= 3) THEN
C           KEEP OLD EXP. PROJ. ANGLES 
            ANGOUT = ANGEXP(1:3)
         ENDIF

         RANGOLD   = 0.0
         XSHOLD    = 0.0
         YSHOLD    = 0.0

         IF (NGOTPAR >= 7 ) THEN
C           USE OLD INPLANE ROT. & SHIFT  
            RANGOLD   = ANGEXP(4)
            XSHOLD    = ANGEXP(5)
            YSHOLD    = ANGEXP(6)

            MIRROROLD = (ANGEXP(7) > 0)
            IF (MIRROROLD) THEN
               WRITE(NOUT,90)
90             FORMAT(
     &           ' *** MIRRORED PROJECTION INPUT NO LONGER ALLOWED.',/,
     &           ' *** CONVERT YOUR ALIGNMENT PARAMETER FILES TO ',
     &           ' NON-MIRRORED FORM OR USE OLDER SPIDER RELEASE.')
               CALL ERRT(101,
     &          'MIRRORED PROJECTION INPUT NO LONGER ALLOWED',NE)
            ENDIF

         ENDIF

         IF (CHNG_ORDER) THEN
C           HAVE TO CHANGE ORDER OF SHIFT & ROTATION.
C           IN 'AP **' IMAGE IS SHIFTED FIRST, ROTATED SECOND.
C           IN 'RT SQ' IT IS ROTATED FIRST, SHIFTED SECOND.
C           THIS CODE CORRESPONDS TO OLD OPERATION: 'SA P'.
            SX     = -XSHRAW              ! BEST X SHIFT
            SY     = -YSHRAW              ! BEST Y SHIFT
	    CO     =  COS(RANGNEW * DGR_TO_RAD)
	    SO     = -SIN(RANGNEW * DGR_TO_RAD)

	    XSHNEW = SX*CO - SY*SO
	    YSHNEW = SX*SO + SY*CO
         ELSE
	    XSHNEW = XSHRAW
	    YSHNEW = YSHRAW
         ENDIF
   
C        COMBINE ROT. & SHIFT WITH PREVIOUS TRANSFORMATIONS
         C       =  COS(RANGNEW * DGR_TO_RAD)
         S       = -SIN(RANGNEW * DGR_TO_RAD)

         XSHOUT  = XSHNEW  + XSHOLD*C - YSHOLD*S
         YSHOUT  = YSHNEW  + XSHOLD*S + YSHOLD*C
         RANGOUT = RANGOLD + RANGNEW

C        LIST ANGLES IN RANGE 0...360
         DO WHILE(RANGOUT < 0.0)
            RANGOUT = RANGOUT + 360.0
         ENDDO
         DO WHILE(RANGOUT >= 360.0)
            RANGOUT = RANGOUT - 360.0
         ENDDO

C        CONVERT ~360 TO: 0
         IF ( ABS(RANGOUT - 360.0) <= ALMOST_ZERO) RANGOUT = 0.0


C        SET FLAG FOR NO ANGDIF DETERMINED
         ANGDIF = -1.0 

         IF (IMGREF <= 0)  THEN
C            NO RELEVANT REF. IMAGE FOUND
             ANGDIF = 0.0

         ELSEIF (GOTREFANG .AND. NGOTPAR >= 3) THEN
C            CAN FIND ANGDIF
             ANGDIF = ABS(EXPDIR(1) * REFDIR(1) + 
     &                    EXPDIR(2) * REFDIR(2) + 
     &                    EXPDIR(3) * REFDIR(3))
             ANGDIF = MIN(1.0,ANGDIF)
             ANGDIF = ACOS(ANGDIF) / DGR_TO_RAD
         ENDIF

         IKEY     = KEY
         IF (CTYPE(1:2) == 'SH'  .OR. 
     &       CTYPE(1:3) == 'FOU' .OR.
     &       CTYPE(1:3) == 'REF') THEN
             IKEY     = IMGEXP
         ENDIF

         IF (SAY_RAW) THEN
C           WANT TO REPORT SHIFT BEFORE ROTATE
            XSHT = XSHRAW
            YSHT = YSHRAW
         ELSE
            XSHT = XSHNEW
            YSHT = YSHNEW
         ENDIF


C        SAVE 15 ALIGNMENT PARAMETERS
         CALL AP_OUT(IKEY,IMGEXP,
     &               ANGOUT, IMGREF,CCROT,
     &               RANGNEW,XSHT,YSHT,MIRRORNEW,
     &               RANGOUT,XSHOUT,YSHOUT,
     &               NPROJ,ANGDIF,PEAKV,
     &               LUNDOC,DLIST,15,IRTFLG)

          END




C++******************************  AP_END *****************************



        SUBROUTINE AP_END(KEY,IMGEXP,IMGREF, 
     &                  ANGREF,REFDIR, 
     &                  ANGEXP,EXPDIR,ISHRANGE,
     &                  GOTREFANG, NGOTPAR, CCROTT,PEAKVT,
     &                  RANGNEW,XSHNEW,YSHNEW, MIRRORNEW,REFPAT,
     &                  NPROJT, CTYPE, LUNDOC,DLIST)

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        CHARACTER (LEN=*)                           :: REFPAT 

        LOGICAL                                     :: MIRROROLD
        LOGICAL                                     :: MIRRORNEW
	LOGICAL                                     :: GOTREFANG

	REAL, DIMENSION(7)                          :: ANGEXP
	REAL, DIMENSION(3)                          :: ANGREF
	REAL, DIMENSION(3)                          :: ANGOUT
	REAL, DIMENSION(3)                          :: EXPDIR,REFDIR
        CHARACTER (LEN=*)                           :: CTYPE 

        REAL, DIMENSION(*)                          :: DLIST

        PARAMETER (QUADPI = 3.14159265358979323846)
        PARAMETER (DGR_TO_RAD = (QUADPI/180))

         NPROJ   = NPROJT
c         IF (NPROJ .EQ. 0) NPROJ = 1 removed aug 09 al

         PEAKV = PEAKVT  

         IF (IMGREF > 0) THEN
C           IMGREF IS NUMBER OF MOST SIMILAR REF. IMAGE 
            CCROT   = CCROTT
         ELSE
C           NO REFERENCE IMAGE SELECTED
            CCROT = -1.0
            PEAKV = 0.0
         ENDIF

C        SET NEW PROJECTION ANGLES
         ANGOUT = 0.0            ! DEFAULT VALUE
         IF (IMGREF > 0 .AND. GOTREFANG) THEN
C           USE REF. ANGLES AS NEW PROJECTION ANGLES
            ANGOUT = ANGREF(1:3)

            IF (MIRRORNEW) THEN
C              REF. PROJECTION MUST BE MIRRORED
               ANGOUT(1) = -ANGOUT(1)
               ANGOUT(2) = 180+ANGOUT(2)
            ENDIF
         ELSEIF (NGOTPAR >= 3) THEN
C           KEEP OLD EXP. PROJ. ANGLES 
            ANGOUT = ANGEXP(1:3)
         ENDIF

         RANGOLD   = 0.0
         XSHOLD    = 0.0
         YSHOLD    = 0.0

         IF (NGOTPAR >= 7 .AND. ISHRANGE > 0) THEN
C           USE OLD INPLANE ROT. & SHIFT  
            RANGOLD   = ANGEXP(4)
            XSHOLD    = ANGEXP(5)
            YSHOLD    = ANGEXP(6)

            MIRROROLD = (ANGEXP(7) > 0)
            IF (MIRROROLD) THEN
               WRITE(NOUT,90)
90             FORMAT(
     &           ' *** MIRRORED PROJECTION INPUT NO LONGER ALLOWED.',/,
     &           ' *** CONVERT YOUR ALIGNMENT PARAMETER FILES TO ',
     &           ' NON-MIRRORED FORM OR USE OLDER SPIDER RELEASE.')
               CALL ERRT(101,
     &          'MIRRORED PROJECTION INPUT NO LONGER ALLOWED',NE)
            ENDIF

         ENDIF
    
C        COMBINE ROT. & SHIFT WITH PREVIOUS TRANSFORMATION
         C       =  COS(RANGNEW * DGR_TO_RAD)
         S       = -SIN(RANGNEW * DGR_TO_RAD)

         XSHOUT  = XSHNEW  + XSHOLD*C - YSHOLD*S
         YSHOUT  = YSHNEW  + XSHOLD*S + YSHOLD*C
         RANGOUT = RANGOLD + RANGNEW

C        LIST ANGLES IN RANGE 0...360
         DO WHILE(RANGOUT .LT. 0.0)
            RANGOUT = RANGOUT + 360.0
         ENDDO
         DO WHILE(RANGOUT .GE. 360.0)
            RANGOUT = RANGOUT - 360.0
         ENDDO

C        CONVERT ~360 TO: 0
         IF ( ABS(RANGOUT - 360.0) <= 0.05) RANGOUT = 0.0


C        SET FLAG FOR NO ANGDIF DETERMINED
         ANGDIF = -1.0 

         IF (IMGREF <= 0)  THEN
C            NO RELEVANT REF. IMAGE FOUND
             ANGDIF = 0.0

         ELSEIF (GOTREFANG .AND. NGOTPAR >= 3) THEN
C            CAN FIND ANGDIF
             ANGDIF = ABS(EXPDIR(1) * REFDIR(1) + 
     &                    EXPDIR(2) * REFDIR(2) + 
     &                    EXPDIR(3) * REFDIR(3))
             ANGDIF = MIN(1.0,ANGDIF)
             ANGDIF = ACOS(ANGDIF) / DGR_TO_RAD
         ENDIF

         IKEY     = KEY
         IF (CTYPE(1:2) == 'SH'  .OR. 
     &       CTYPE(1:3) == 'FOU' .OR.
     &       CTYPE(1:3) == 'REF') THEN
             IKEY     = IMGEXP
         ENDIF

C        SAVE 15 ALIGNMENT PARAMETERS
         CALL AP_OUT(IKEY,IMGEXP,
     &               ANGOUT, IMGREF,CCROT,
     &               RANGNEW,XSHNEW,YSHNEW,MIRRORNEW,
     &               RANGOUT,XSHOUT,YSHOUT,
     &               NPROJ,ANGDIF,PEAKV,
     &               LUNDOC,DLIST,15,IRTFLG)

          END


C++****************************** AP_OUT *****************************

        SUBROUTINE AP_OUT(KEY,IMGEXP,
     &                    ANGOUT, IMGREF,CCROT,
     &                    RANGNEW,XSHNEW,YSHNEW,MIRRORNEW,
     &                    RANGOUT,XSHOUT,YSHOUT,
     &                    NPROJ,ANGDIF,PEAKV,
     &                    LUNDOC,DLIST,NLIST,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'

        INTEGER           :: KEY,IMGEXP
        REAL              :: ANGOUT(3)
        INTEGER           :: IMGREF
        REAL              :: CCROT,RANGNEW,XSHNEW,YSHNEW
        LOGICAL           :: MIRRORNEW
        REAL              :: RANGOUT,XSHOUT,YSHOUT 
        INTEGER           :: NPROJ
        REAL              :: ANGDIF,PEAKV
        INTEGER           :: LUNDOC
        REAL              :: DLIST(NLIST)
        INTEGER           :: NLIST,IRTFLG

        REAL, PARAMETER   :: ALMOSTZERO = 0.05
        INTEGER           :: NSEL_USED

        CHARACTER(LEN=150) :: FORMOUT

        !FORMOUT  = '(I10,1X,I4,10000(1X,1PG13.6))'  old format 3/30/12

C       <,<,<, REF#,IMG#,
C       ROT<, SX,SY,NPROJ, <DIF,CCROT,
C       ROT<,SX,SY,PEAK
        FORMOUT  = 
     &      '(I7,1X,I2,1X,F7.2,1X,F7.2,1X,F7.2,2X,F6.0,2X,F7.0,2X,'// 
     &      'F7.2,1X,F7.2,1X,F7.2,2X,F6.0,1X, F7.2,1X, F11.2,1X,'  //
     &      'F7.2,1X,F7.2,1X,F7.2,2X,F7.4)'              
C            123456789 123456789 123456789 123456789 123456789 1234567890123456789012

C       ZERO DLIST ARRAY
        DLIST     = 0.0

        DLIST(1)  = ANGOUT(1)
        DLIST(2)  = ANGOUT(2)
        DLIST(3)  = ANGOUT(3)
        DLIST(4)  = IMGREF
        DLIST(5)  = IMGEXP

        DLIST(6)  = RANGOUT
        IF (ABS(RANGOUT) < ALMOSTZERO) DLIST(6) = 0.0
        DLIST(7)  = XSHOUT
        IF (ABS(XSHOUT)  < ALMOSTZERO) DLIST(7) = 0.0
        DLIST(8)  = YSHOUT
        IF (ABS(YSHOUT)  < ALMOSTZERO) DLIST(8) = 0.0

        DLIST(9)  = NPROJ

        DLIST(10) = ANGDIF
        IF (ABS(ANGDIF) < ALMOSTZERO) DLIST(10) = 0.0

        DLIST(11) = CCROT

        DLIST(12) = RANGNEW
        IF (ABS(RANGNEW) < ALMOSTZERO) DLIST(12) = 0.0
        DLIST(13) = XSHNEW
        IF (ABS(XSHNEW)  < ALMOSTZERO) DLIST(13) = 0.0
        DLIST(14) = YSHNEW
        IF (ABS(YSHNEW)  < ALMOSTZERO) DLIST(14) = 0.0

        DLIST(15) = ABS(PEAKV)   ! IN CASE OF NEGATIVE
        IF (MIRRORNEW) DLIST(15) = -(ABS(PEAKV)) 

#ifndef USE_MPI
         IF (LUNDOC .GT. 0) THEN
C           SAVE IN ALIGNMENT DOC FILE IF NOT USING MPI
C       <,<,<, REF#,IMG#,ROT<, SX,SY,NPROJ, <DIF,CCROT,ROT<,SX,SY,PEAK

            IF (USE_LONGCOL) THEN
              CALL LUNDOCWRTDAT(LUNDOC,KEY,DLIST,NLIST,IRTFLG)

            ELSE
              CALL LUNDOCWRTDATF(LUNDOC,KEY,DLIST,NLIST,FORMOUT,IRTFLG)
            ENDIF
         ENDIF
#endif
 
         CALL REG_GET_USED(NSEL_USED)

         IF (NSEL_USED .GT. 0) THEN
C            OUTPUT TO SPIDER'S REGISTERS
             CALL REG_SET_NSEL(1,5, DLIST(1),DLIST(2),DLIST(3),
     &                              DLIST(4),DLIST(5),IRTFLG)
             CALL REG_SET_NSEL(6,5, DLIST(6),DLIST(7),DLIST(8),
     &                              DLIST(9),DLIST(10),IRTFLG)
             CALL REG_SET_NSEL(11,5,DLIST(11),DLIST(12),DLIST(13),
     &                              DLIST(14),DLIST(15),IRTFLG)
          ENDIF

          END

C+**********************************************************************

       SUBROUTINE AP_END_HEAD(IMNUM,FILPAT,LUNIN,PARLIST,NVALS,IRTFLG)

C      PURPOSE: PUT ANGLES, ETC AS HEADER VALUES IN IMAGE

       IMPLICIT NONE
       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC'

       INTEGER               :: IMNUM
       CHARACTER(LEN=*)      :: FILPAT
       INTEGER               :: LUNIN
       REAL                  :: PARLIST(NVALS)
       INTEGER               :: NVALS
       INTEGER               :: IRTFLG

       CHARACTER(LEN=MAXNAM) :: FILNAM
       INTEGER               :: MAXIM,ITYPE,NX,NY,NZ,NLET
       LOGICAL               :: FOUROK = .FALSE.

C      OPEN EXISTING IMAGE FILE FOR OUTPUT
       IF (IMNUM >= 0) THEN
C         CREATE IMAGE FILE NAME FIRST
          NLET = 0
          CALL FILGET(FILPAT,FILNAM,NLET,IMNUM,IRTFLG)
          IF (IRTFLG .NE. 0)  RETURN
       ELSE
          FILNAM = FILPAT
       ENDIF
 
       MAXIM = 0
       CALL OPFILEC(0,.FALSE.,FILNAM,LUNIN,'O',ITYPE,
     &              NX,NY,NZ,MAXIM,' ',FOUROK,IRTFLG)
       IF (IRTFLG .NE. 0)  RETURN

C      PUT ANGLES, ETC IN HEADER  (MRC OK)
       CALL LUNSETVALS(LUNIN,IAPLOC+1,NVALS,PARLIST,IRTFLG)

       CLOSE(LUNIN)

       END
