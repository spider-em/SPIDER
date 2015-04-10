
C **********************************************************************
C
C   SHOSLI.FOR  -- CREATED OCT 90
C **********************************************************************
C *  AUTHOR: ArDean Leith 
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
C      SHOSLI(LUNOUT,BUF,NSAM,NREC1,NREC2,SLICE)  
C
C      PURPOSE:     PRINTS SLICE ARRAY 
C
C      PARAMETERS:  
C
C      CALLS:       WRTLIN 
C
C      CALLED BY:   CONINT
C
C--********************************************************************

         SUBROUTINE SHOSLI(LUNOUT,BUF,NSAM,NREC1,NREC2,SLICE)
 
         INTEGER   * 2  SLICE(*)
         DIMENSION      BUF(*)

         DO  I = NREC1,NREC2
            IPTR = (I-1) * NSAM + 1

            WRITE(LUNOUT,90) (SLICE(J),J=IPTR,IPTR+NSAM-1)
90          FORMAT(4I3)

         END DO

         WRITE(LUNOUT,*) ' '

         END
    
