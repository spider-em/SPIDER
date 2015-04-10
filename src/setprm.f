
C++********************************************************************
C
C SETPRM.FOR     MODIFIED:                      4/22/96 ARDEAN LEITH
C                REMOVED SETPRMB PARAMETERS     MAY 09  ARDEAN LEITH
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
C    SETPRM(LUN,IDUM,IDUM, FMAXD,FMIND,AVD, CDUM)
C
C    PURPOSE: WILL SET LABEL PARAMETERS IDENTIFYING TYPE, NORMALIZATION 
C             STATUS OF FILE, AND WRITE HEADER LABEL INTO FILE.
C
C    PARAMETERS:
C      LUN            LOGICAL UNIT NUMBER OF FILE TO BE LABELED
C      NSAM           NUMBER OF SAMPLES
C      IDUM           UNUSED
C      FMIND,FMAXD    IF AVAILABLE, THE MINIMUM AND MAXIMUM OF THE
C                     IMAGE STORED IN THE FILE (OR ZERO IF UNKNOWN)
C      AVD            AVERAGE VALUE
C      CDUM           UNUSED
C
C    NOTE:            DO NOT USE THIS ROUTINE ANY MORE.  USE
C                     SETPRMB INSTEAD!!!!!!!!!!!!!!!!!!!!!
C
C--*******************************************************************

      SUBROUTINE SETPRM(LUN, IDUM1,IDUM2, FMAXD,FMIND,AVD, CDUM)

C     'CMBLOCK HOLDS SIG!!
      INCLUDE 'CMBLOCK.INC'
 
      CHARACTER(LEN=2) :: CDUM

      CALL SETPRMB(LUN, FMAXD,FMIND, AVD,SIG)
      
      END
                           
