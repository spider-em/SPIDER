
C++********************************************************************
C
C SETPRMB.F      MODIFIED                       4/22/96  ArDean Leith
C                LUNRED                         FEB 03   ArDean Leith
C                REMOVED UNUSED, ADDED SIGD     MAY 09   ArDean Leith
C                MRC FRIENDLY                   AUG 19   ArDean Leith
C               
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-20190  Health Research Inc.,                        *
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
C SETPRMB(LUN, FMAXD,FMIND, AVD,SIGD)
C
C PURPOSE:
C    SET HEADER PARAMETERS FOR  NORMALIZATION STATUS OF FILE, AND 
C    WRITE HEADER LABEL INTO FILE.
C
C    SETPRMB(LUN, FMAXD,FMIND, AVD,SIGD)
C      LUN            LOGICAL UNIT NUMBER OF FILE TO BE LABELED  (SENT)
C      FMIND,FMAXD    IF AVAILABLE, THE MINIMUM & MAXIMUM OF THE
C                     IMAGE IS STORED IN THE FILE                (SENT)
C      AVD            AVERAGE VALUE                              (SENT)
C      SIGD           S.D. VALUE                                 (SENT)
C
C    NOTES:  THE HEADER RECORD OF SPIDER FILE CONTAINS FOLLOWING INFO:
C            LOCATION # 5  FLAG INDICATING DATA TYPE (=IFORM)
C                       6  IMAMI = FLAG INDICATING IF THE IMAGE HAS 
C				BEEN SEARCHED FOR MAX AND MIN. 
C				IF FMAXD=FMIND IN SETPRM ARGUMENTS THEN
C				IMAMI IS SET TO 0, OTHERWISE, TO +1.
C                       7  FMAXD = IMAGE MAXIMUM
C                       8  FMIND = IMAGE MINIMUM
C                       9  AVD   = IMAGE AVERAGE
C                      10  SIG   = STANDARD DEVIATION (SQ. ROOT OF VARIANCE)
C
C       NOTE: ROUTINE REPLACES SETPRM WHICH SHOULD NO LONGER BE USED!
C
C--*******************************************************************

      SUBROUTINE SETPRMB(LUN, FMAXD,FMIND, AVD,SIGD)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
 
      INTEGER          :: LUN
      REAL             :: FMAXD,FMIND, AVD,SIGD

      CHARACTER(LEN=2) :: TYPE
      INTEGER          :: IRTFLG

C     UPDATE THE HEADER VALUES (?? THIS ALSO SETS COMMON: IMAMI !!)
      IMAMI = 1
      IF (FMAXD <= FMIND) THEN
C        SET IMAMI TO UNDETERMINED (ZERO)
         IMAMI = 0
      ENDIF

      !write(6,*) ' In setprmb - imami,fmin:',imami,fmind

C     UPDATE THE INCORE HEADER STATISTICS    (MRC OK)
      CALL LUNSETSTAT(LUN,IMAMI, FMIND,FMAXD, AVD,SIGD,IRTFLG)

C     WRITE UPDATED HEADER BACK IN THE FILE  (MRC OK)
      CALL LUNWRTCURHED(LUN,IRTFLG)  
      
      END
                           
C----------------------------- SETPRMS --------------------

      SUBROUTINE SETPRMS(LUN, SCALE,IRTFLG)

C     PURPOSE: RESET PIXSIZ HEADER PARAMETERS FOR SCALING
      IMPLICIT NONE

      INTEGER          :: LUN,IRTFLG
      REAL             :: SCALE

      REAL             :: PIXSIZ,PIXSIZOLD

C     GET CURRENT FILE HEADER VALUE FOR PIXSIZ  (MRC OK)
      CALL LUNGETPIXSIZ(LUN,PIXSIZOLD,IRTFLG)

C     RESCALE PIXSIZ
      PIXSIZ = PIXSIZOLD * SCALE

C     SET CURRENT FILE HEADER VALUE FOR PIXSIZ  (MRC OK)
      CALL LUNSETPIXSIZ(LUN,PIXSIZ,IRTFLG)

C     WRITE UPDATED HEADER BACK IN THE FILE     (MRC OK)
      CALL LUNWRTCURHED(LUN,IRTFLG)  
      
      END
