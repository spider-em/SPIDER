
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
C  GETLAB(LUN,IGO,NVAL,VALUE,IRTFLG)
C
C  PURPOSE:   RETRIEVES VARIABLES FROM THE HEADER 
C             BY BUFFER POSITION NUMBER
C
C  PARAMETERS:
C             LUN        LOGICAL UNIT NUMBER OF FILE 
C             IGO        FIRST BUFFER POSITION TO BE RETRIEVED
C             NVAL       NUMBER OF BUFFER POSITIONS TO BE RETRIEVED
C             VALUES     ARRAY FOR BUFFER VALUES RETRIEVED
C             IRTFLG     ERROR FLAG
C 
C       NOTE:   THE HEADER RECORD(S) OF SPIDER FILES CONTAINS FOLLOWING 
C               POSITIONS WHICH CAN BE RETRIEVED (AMONG OTHERS)
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
  
C       ----------- GETLAB -------------------------------------

        SUBROUTINE GETLAB(LUN,IGO,NVAL,VALUES,IRTFLG)

        IMPLICIT NONE

        INTEGER   :: LUN,IGO,NVAL,IRTFLG
        REAL      :: VALUES(*)   

        INTEGER   :: ISTOP

        ISTOP = IGO + NVAL - 1
        IF (ISTOP > 256) THEN
           CALL ERRT(102,'MAX HEADER LOCATION MUST BE < 257',ISTOP)
           IRTFLG = 1
           RETURN
        ENDIF

C       COPY THE HEADER LOCATIONS INTO: VALUES
        CALL LUNGETVALS(LUN,IGO,NVAL,VALUES,IRTFLG)
                 
        END


C       ----------- GETLAB_R_MRC -------------------------------------

        SUBROUTINE GETLAB_R_MRC(LUN,IGO,NVAL,VALUES,IRTFLG)

        IMPLICIT NONE

        REAL      :: VALUES(*)   
        INTEGER   :: LUN,IGO,NVAL,IRTFLG

        INTEGER   :: ISTOP

        ISTOP = IGO + NVAL - 1
        IF (ISTOP > 56) THEN
           CALL ERRT(102,'MAXIMUM HEADER LOCATION MUST BE < 57',ISTOP)
           IRTFLG = 1
           RETURN
        ENDIF

C       COPY THE HEADER LOCATIONS INTO: VALUES (AS FLOATS)
        CALL LUNGETVALS_R_MRC(LUN,IGO,NVAL,VALUES,IRTFLG)
                         
        END


