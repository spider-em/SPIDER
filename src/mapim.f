
C **********************************************************************
C
C  MAPIM.FOR  -- CREATED OCT 90
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
C      MAPIM(LUNIN,LUNOUT,NSAM,NREC1,NREC2,TABLE,NMAX,BUF,IRTFLG) 
C
C      PURPOSE:  APPLIES A LOOKUP TABLE TO EACH COLUMN OF AN IMAGE 
C                FOR ALL RECORDS FROM NREC1 TO NREC2.  THE 
C                LOOK-UP-TABLE IS A FLOATING POINT TABLE.  IF OLD-
C                VALUE IS 1 THEN NEW-VALUE IS SET TO TABLE(1), ETC.
C                NO ERROR CHECKING IS DONE (TO SAVE TIME) !!!!!!! 
C
C      PARAMETERS:  
C
C      CALLED BY:  CONINT 
C
C      CALLS:      REDLIN      WRTLIN   
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

       SUBROUTINE MAPIM(LUNIN,LUNOUT,NSAM,NREC1,NREC2,TABLE,NMAX,
     &                  BUF,IRTFLG)

 

       DIMENSION TABLE(*),BUF(*)

C      GO THRU THE STACK SLICE BY SLICE AND MAKE REPLACEMENTS FROM TABLE
         
       DO IREC = NREC1,NREC2
          CALL REDLIN(LUNIN,BUF,NSAM,IREC)

          DO ICOL = 1,NSAM
C            CHECK THIS VOXEL FOR REPLACEMENT
             IT = BUF(ICOL)
             IF (IT .GT. 0)  BUF(ICOL) = TABLE(IT) 
          END DO

          CALL WRTLIN(LUNOUT,BUF,NSAM,IREC)
       END DO

       RETURN
       END
