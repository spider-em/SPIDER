
C++*********************************************************************
C
C CSVMS.FOR  -- CREATED                             JAN 85 R. BANERJEE
C               SAYIT                               APR 13 ARDEAN LEITH                            
C **********************************************************************
C *  AUTHOR: R. BANERJEE                                                   *
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2013  Health Research Inc.,                         *
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
C    CSVMS(COMLIN,SAYIT,IRET)
C
C    PURPOSE:  EXECUTE SYSTEM COMMANDS LIKE COPY, RENAME, 
C              DELETE FILES, ETC.
C
C    PARAMETERS:
C        SAYIT    FOR ECHO TO SCREEN & RESULTS
C        COMLIN   COMMAND LINE FOR THE SYSTEM COMMAND 
C                 CALLER PUTS THE COMMAND IN THE COMLIN
C
C--********************************************************************

	SUBROUTINE CSVMS(COMLIN,SAYIT,IRET)

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'

	CHARACTER (LEN=*)  :: COMLIN
        LOGICAL            :: SAYIT
        INTEGER            :: IRET

        INTEGER            :: NCHAR
        INTEGER            :: system,lnblnk

        NCHAR = lnblnk(COMLIN)
        IF (NCHAR <= 0) NCHAR = LEN(COMLIN)

        IF (SAYIT) THEN
           WRITE(NOUT,102)COMLIN(1:NCHAR)
           IF (NDAT .NE. NOUT) WRITE(NDAT,102)COMLIN(1:NCHAR)

  102      FORMAT(/,1X,A)
        ENDIF

 	IRET = system(COMLIN(1:NCHAR))

	END

