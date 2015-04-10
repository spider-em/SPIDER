
C++********************************************************************
C
C SETPRMB.F      MODIFIED                       4/22/96 ARDEAN LEITH
C                LUNRED                         FEB 03  ARDEAN LEITH
C                REMOVED UNUSED, ADDED SIGD     MAY 09  ARDEAN LEITH
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
C SETPRMB(LUN, FMAXD,FMIND, AVD,SIGD)
C
C PURPOSE:
C    WILL SET HEADER PARAMETERS FOR  NORMALIZATION STATUS OF FILE, AND 
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
C       NOTE: THIS ROUTINE REPLACES SETPRM.  SETPRM SHOULD NO LONGER
C             BE USED!
C
C--*******************************************************************

      SUBROUTINE SETPRMB(LUN, FMAXD,FMIND, AVD,SIGD)

      INCLUDE 'CMBLOCK.INC'
 
      CHARACTER(LEN=2) :: TYPE

C     UPDATE THE HEADER VALUES
      IMAMI = 1
      IF (FMAXD .EQ. FMIND) THEN
C        SET IMAMI TO UNDETERMINED (ZERO)
         IMAMI  = 0
      ENDIF

C     UPDATE THE INCORE HEADER STATISTICS 
      CALL LUNSETSTAT(LUN,IMAMI, FMIND,FMAXD, AVD,SIGD,IRTFLG)

C     WRITE UPDATED HEADER BACK IN THE FILE
      CALL LUNWRTCURHED(LUN,IRTFLG) 
      
      RETURN
      END
                           

      SUBROUTINE SETPRMS(LUN, SCALE,IRTFLG)

C     UPDATE THE INCORE HEADER VALUE AND FILE HEADER FOR PIXSIZ
      CALL LUNGETVALS(LUN,38,1,PIXSIZOLD,IRTFLG)
      PIXSIZ = PIXSIZOLD * SCALE
      CALL LUNSETVALS(LUN,38,1,PIXSIZ,IRTFLG)
      
      RETURN
      END
