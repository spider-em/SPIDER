
C++*********************************************************************
C
C IMSTAT.FOR -- CREATED NOV 90 
C **********************************************************************
C *  AUTHOR:  ArDean Leith
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
C    IMSTAT(IRTFLG)
C
C    PURPOSE:       CALLS MAPSTAT TO GET IMAGE STATISTICS.
C
C    PARAMETERS     IRTFLG       ERROR RETURN
C
C    CALLS:               
C
C    CALLED BY:        
C
C        0         2         3         4         5         6         7     
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

       SUBROUTINE IMSTAT(IRTFLG)

       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC' 
 
       CHARACTER(LEN=MAXNAM)   ::  IMFILE

       DATA LUNIM,LUND/20,21/

       IRTFLG = 1

       MAXIM  = 0
20     CALL OPFILEC(0,.TRUE.,IMFILE,LUNIM,'OLD',IFORM,NSAM,NROW,NSLICE,
     &             MAXIM,'CLUSTER INPUT',.FALSE.,IRTFLG)
       IF (IRTFLG .EQ. -1) RETURN

       IF (FMIN .LT. 0.0) THEN
C         WILL NOT WORK RIGHT
          WRITE(NOUT,*) ' *** CAN NOT USE FILE WITH FMIN < 0.0'
          CLOSE(LUNIM)
          CALL ERRT(100,'IMSTAT',NE)
          RETURN
       ENDIF

       NREC2 = NROW * NSLICE
       CALL MAPSTAT(LUNIM,LUND,NSAM,NROW,NSLICE,1,NREC2,NVOX,IRTFLG)            
       CLOSE(LUNIM)
         
       RETURN
       END




