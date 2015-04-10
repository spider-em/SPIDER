
C **********************************************************************
C
C   FILSLI.FOR  -- CREATED OCT 90
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
C      FILSLI(LUNIM,BUF,NSAM,NREC1,NREC2,THLEV,SLICE)  
C
C      PURPOSE:     READS SPIDER PICTURE FILES SLICES INTO SLICE ARRAY 
C
C      PARAMETERS:  
C
C      CALLED BY:   CONINT
C
C--********************************************************************

         SUBROUTINE FILSLI(LUNIM,BUF,NSAM,NREC1,NREC2,THRESH,
     &                     THLEV,SLICE)
 
 

         INTEGER   * 2  SLICE(*)
         DIMENSION      BUF(*)
         LOGICAL        THRESH

C        READ THE SPIDER FILE INTO SLICE ARRAY

         IPTR = 0 

         IF (THRESH) THEN
           DO  I = NREC1,NREC2
             CALL REDLIN(LUNIM,BUF,NSAM,I)

             DO  J = 1,NSAM
               IPTR = IPTR + 1

               IF (BUF(J) .GT. THLEV) THEN
C                 VOXEL IS ABOVE THRESHOLD
                  SLICE(IPTR) = -1
               ELSE
                  SLICE(IPTR) = 0
               ENDIF
             END DO
           END DO

         ELSE

           DO  I = NREC1,NREC2
             CALL REDLIN(LUNIM,BUF,NSAM,I)

             DO  J = 1,NSAM
               IPTR = IPTR + 1
               SLICE(IPTR) = BUF(J)
             END DO
           END DO

         ENDIF

         RETURN
         END
    
