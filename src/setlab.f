
C++*********************************************************************
C
C SETLAB.F      CREATED                           NOV 87  ArDean Leith
C               LUNRED                            FEB 03  ARDEAN LEITH
C               1PG FORMAT                        NOV 10  ARDEAN LEITH
C               LINE FORMAT, MPISET               NOV 13  ARDEAN LEITH
C
C **********************************************************************
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
C    SETLAB(LUN,BUF,IGO,NBUF,VALUE,TYPE,IRTFLG)
C
C    PURPOSE:    THIS SUBROUTINE SETS HEADER PARAMETERS BY BUFFER NUMBER
C                WRITES THE HEADER INTO THE FILE.
C
C    PARAMETERS:
C         
C             LUN          LOGICAL UNIT NUMBER OF FILE 
C             BUF          WORK SPACE FOR READ/WRITE BUFFER
C             IGO          FIRST BUFFER POSITION TO BE SET
C             NBUF         NUMBER OF BUFFER POSITIONS TO BE SET
C             VALUES       ARRAY FOR BUFFER VALUES TO BE SET
C             TYPE         CHARACTER VARIABLE CONTAINING FLAG FOR IFORM
C			   TYPE SYMBOL	DATA TYPE	        IFORM
C                           R     2-D REAL                       +1
C                           R3    3-D REAL                       +3
C                           P     2-D POLAR                      +2
C                           D     NON-IMAGE DATA                  0
C                           F     2-D FOURIER                    -1
C                           F3    3-D FOURIER                    -3
C                           U     UNCHANGED	 UNCHANGED	
C                           O2    2-D FOURIER, MIXED RADIX ODD   -11  
C                           E2    2-D FOURIER, MIXED RADIX EVEN  -12
C                           O3    3-D FOURIER, MIXED RADIX ODD   -21
C                           E3    3-D FOURIER, MIXED RADIX EVEN  -22
C
C             IRTFLG       ERROR FLAG  (-1 ON ENTRY SUPRESSES PRINT-OUT)
C
C 
C       NOTE:   THE HEADER RECORD(S) OF THE FILE CONTAIN THE FOLLOWING 
C               BUFFER POSITIONS WHICH CAN BE ALTERED (AMONG OTHERS):
C               POSITION   1  NSLICE = NUMBER OF SLICES (PLANES) IN VOLUME
C                                (=1 FOR AN IMAGE)  ON VAX  LONG HEADER
C                                FORMAT FILES THE VALUE OF NSLICE STORED IN
C                                THE FILE IS NEGATIVE.
C                          2  NROW   = NUMBER OF ROWS PER SLICE
C                          3  IREC   = (UNUSED)
C                          4  NHISTREC = (UNUSED)
C                          5  FLAG INDICATING DATA TYPE (=IFORM)
C                          6  IMAMI = FLAG INDICATING IF THE IMAGE HAS 
C				      BEEN SEARCHED FOR MAX AND MIN. 
C				      IMAMI IS SET TO +1 IF SEARCHED.
C                          7  FMAXD = IMAGE MAXIMUM
C                          8  FMIND = IMAGE MINIMUM
C			   9  AVD   = IMAGE AVERAGE
C                         10  SIG   = STANDARD DEVIATION (SQ. ROOT OF VARIANCE)
C                         11  IHIST = UNUSED
C                         13  LABLN = NUMBER OF FLOATING POINT HEADER VARIABLES
C                         14  IANGLE= ANGLE FILL FLAG
C                         15  PHI   = TILT ANGLE
C                         16  THETA = TILT ANGLE
C                         17  PSI   = TILT ANGLE
C                         18  XOFF  = X TRANSLATION
C                         19  YOFF  = Y TRANSLATION
C                         20  ZOFF  = Z TRANSLATION
C
C--*********************************************************************
  
        SUBROUTINE SETLAB(LUN,NSAM,UNUSED,IGO,NBUF,VALUES,TYPE,IRTFLG)
  
        INCLUDE 'CMBLOCK.INC' 

        REAL             ::  UNUSED(*),VALUES(*)
        CHARACTER   *(*) ::  TYPE
        LOGICAL          ::  PRNT

C       AUTOMATIC ARRAYS    
        REAL             :: OLDVALUES(NBUF)

        INTEGER          :: ICOMM,MYPID,MPIERR
        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

        PRNT = .TRUE.
        IF (IRTFLG == -1) PRNT = .FALSE.
        IRTFLG = 1

                 
C       UPDATE THE HEADER BUFFER
        IF (TYPE(1:1) .NE. 'U') THEN
           LENC = LEN(TYPE)
           IF     (LENC == 2 .AND. TYPE(1:2) == 'R3') THEN 
                ITYPE = 3
           ELSEIF (LENC == 2 .AND. TYPE(1:2) == 'O2') THEN
                ITYPE = -11
           ELSEIF (LENC == 2 .AND. TYPE(1:2) == 'E2') THEN
                ITYPE = -12
           ELSEIF (LENC == 2 .AND. TYPE(1:2) == 'O3') THEN
                ITYPE = -21
           ELSEIF (LENC == 2 .AND. TYPE(1:2) == 'E3') THEN
                ITYPE = -22
           ELSEIF (TYPE(1:1) == 'P') THEN
                ITYPE = 2
           ELSEIF (TYPE(1:1) == 'D') THEN
                ITYPE = 0
           ELSE
                ITYPE = 1
           ENDIF
           CALL LUNSETTYPE(LUN,ITYPE,IRTFLG)

           IF (PRNT .AND. MYPID <= 0) THEN 
              WRITE(NOUT,9993) ITYPE,TYPE
9993          FORMAT('  NEW IFORM:',I5,' TYPE:',A2)
           ENDIF
        ENDIF

        ISTOP = MIN(256,IGO+NBUF-1)
        NVAL  = ISTOP - IGO + 1

        CALL LUNGETVALS(LUN,IGO,NVAL,OLDVALUES,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CALL LUNSETVALS(LUN,IGO,NVAL,VALUES,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (PRNT .AND. VERBOSE) THEN
           J = 0
           DO I = IGO,ISTOP
              J = J + 1
              IF (MYPID <= 0) THEN
                 WRITE(NOUT,9999) I,OLDVALUES(J),VALUES(J)
9999             FORMAT('  HEADER LOCATION: ',I3,' CHANGED FROM: ',
     &                    1PG10.3,' TO: ',1PG10.3)
              ENDIF
           ENDDO
        ENDIF

C       WRITE ALTERED HEADER BACK IN THE FILE
        CALL LUNWRTCURHED(LUN,IRTFLG) 

        END

