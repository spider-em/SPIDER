
C++*********************************************************************
C
C GETLAB.FOR -- CREATED NOV 87  BY ArDean Leith
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
C  GETLAB(LUN,NSAM,BUF,IBUF1,NBUF,VALUE,IRTFLG)
C
C  PURPOSE:    THIS SUBROUTINE RETRIEVES VARIABLES FROM THE HEADER 
C              BY BUFFER POSITION NUMBER
C
C  PARAMETERS:
C             LUN          LOGICAL UNIT NUMBER OF FILE 
C             NSAM         NUMBER OF SAMPLES IN FILE
C             UNUSED       UNUSED 
C             IBUF1        FIRST BUFFER POSITION TO BE RETRIEVED
C             NBUF         NUMBER OF BUFFER POSITIONS TO BE RETRIEVED
C             VALUES       ARRAY FOR BUFFER VALUES RETRIEVED
C             IRTFLG       ERROR FLAG
C 
C       NOTE:   THE HEADER RECORD(S) OF THE FILE CONTAINS THE FOLLOWING 
C               BUFFER POSITIONS WHICH CAN BE RETRIEVED (AMONG OTHERS)
C                   7  FMAXD = IMAGE MAXIMUM
C                   8  FMIND = IMAGE MINIMUM
C                   9  AVD   = IMAGE AVERAGE
C                  10  SIG   = STANDARD DEVIATION (SQ. ROOT OF VARIANCE)
C                  14  IANGLE= FLAG INDIACATING THAT THERE ARE TILT
C                              ANGLES ARE STORED
C                  15  PHI   = TILT ANGLE
C                  16  THETA = TILT ANGLE
C                  17  GAMMA = TILT ANGLE
C                  18  XOFF  = X OFFSET
C                  19  YOFF  = Y OFFSET
C                  20  ZOFF  = Z OFFSET
C                  21  SCALE = SCALE
C
C--*******************************************************************
  
        SUBROUTINE GETLAB(LUN,NSAM,UNUSED,IBUF1,NBUF,VALUES,IRTFLG)

        DIMENSION VALUES(*), UNUSED(1)   

        ISTOP = IBUF1 + NBUF - 1
        IF (ISTOP .GT. 256) THEN
           CALL ERRT(102,'MAXIMUM HEADER LOCATION MUST BE < 257',ISTOP)
           IRTFLG = 1
           RETURN
        ENDIF

        NVAL = ISTOP - IBUF1 + 1

C       COPY THE HEADER INTO VALUES
        CALL LUNGETVALS(LUN,IBUF1,NVAL,VALUES,IRTFLG)
                 
        RETURN
        END

