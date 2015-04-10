C++*********************************************************************
C
C FOUR1A_FP.F                                  
C                 ADAPTED TO NEW FOURIER FORMAT    7/24/00 BIMAL                  
C                 OPFILE                            NOV 00 ARDEAN LEITH
C                 OPFILEC                           FEB 03 ARDEAN LEITH    
C                 SETPRMS FOR PIXSIZ                NOV 10 ARDEAN LEITH    
C                 DEFAULT RATIO                     NOV 10 ARDEAN LEITH    
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
C  FOUR1A_FP
C
C  PURPOSE: 2D OR 3D IMAGES OF ANY (EVEN/ODD) DIMENSION  IS 
C           INTERPOLATED TO ANY INTGER DIMENSION. 
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE FOUR1A_FP

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC'

        CHARACTER (LEN=MAXNAM) :: FILNAM
        CHARACTER (LEN=1)      :: NULL = CHAR(0)

        REAL, ALLOCATABLE      :: X(:,:)  
        REAL, ALLOCATABLE      :: Y(:,:)     
        REAL, ALLOCATABLE      :: X3(:,:,:) 
        REAL, ALLOCATABLE      :: Y3(:,:,:) 
 
	INTEGER, PARAMETER     :: LUN1 = 50
	INTEGER, PARAMETER     :: LUN2 = 51

C       OPEN OLD FILE
        MAXIM   = 0       
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NSAM1,NROW1,NSLICE1,
     &             MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C	GET OUTPUT FILE  NAME
        CALL FILERD(FILNAM,NLETO,NULL,'OUTPUT',IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 9000

        NSAM2   = 0
        NROW2   = 0
        NSLICE2 = 0

        IF (NSLICE1 > 1) THEN
           CALL RDPRI3S(NSAM2,NROW2,NSLICE2,NOT_USED,
     &                'DIMENSIONS, NSAM, NROW, & NSLICE',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

C	   THE USER IS ALLOWED TO ENTER ONLY
C	   ONE DIMENSION, THE OTHER DIMS ARE COMPUTED TO KEEP THE SAME 
C	   SIZE RELATION AS INPUT.
           IF (NROW2 .LE. 0) THEN
              NROW2 = (FLOAT(NSAM2)   / FLOAT(NSAM1)) * FLOAT(NROW1)
           ENDIF
           IF (NSLICE2 .LE. 0) THEN
              NSLICE2 = (FLOAT(NSAM2) / FLOAT(NSAM1)) * FLOAT(NSLICE1)
           ENDIF
        ELSE
           NSLICE2 = 1
           CALL RDPRIS(NSAM2,NROW2,NOT_USED,
     &                'DIMENSIONS, NSAM & NROW',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

C	   THE USER IS ALLOWED TO ENTER ONLY
C	   ONE DIMENSION, THE OTHER DIMS ARE COMPUTED TO KEEP THE SAME 
C	   SIZE RELATION AS INPUT.
           IF (NROW2 .EQ. 0) THEN
              NROW2 = (FLOAT(NSAM2) / FLOAT(NSAM1)) * FLOAT(NROW1)
           ENDIF

        ENDIF

C	OPEN THE OUTPUT FILE
        MAXIM2 = 0
	CALL OPFILEC(LUN1,.FALSE.,FILNAM,LUN2,'U',IFORM,
     &             NSAM2,NROW2,NSLICE2, MAXIM2,'OUTPUT',.TRUE.,IRTFLG)
	IF (IRTFLG .NE. 0) GOTO 9000

        IF (NSLICE1 .GT. 1) THEN
C          3D CASE

           LSD  = NSAM1+2-MOD(NSAM1,2)
           LSDN = NSAM2+2-MOD(NSAM2,2)

           ALLOCATE (X3(LSD, NROW1,NSLICE1),
     &               Y3(LSDN,NROW2,NSLICE2),STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN 
              MWANT = LSD*NROW1*NSLICE1 + LSDN*NROW2*NSLICE2
              CALL ERRT(46,'FOUR1A_FP, X3 & Y3',MWANT)
              GOTO 9001
           ENDIF 

           CALL READV(LUN1,X3,LSD,NROW1,NSAM1,NROW1,NSLICE1)

           CALL FINT3(X3,Y3,NSAM1,NROW1,NSLICE1,NSAM2,NROW2,
     &                NSLICE2,LSD,LSDN)

           CALL WRITEV(LUN2,Y3,LSDN,NROW2,NSAM2,NROW2,NSLICE2)

C          VOLUME FILE HEADER FOR PIXSIZ HAS CHANGED
           SCALEX = FLOAT(NSAM1)   / FLOAT(NSAM2)
           SCALEY = FLOAT(NROW1)   / FLOAT(NROW2)
           SCALEZ = FLOAT(NSLICE1) / FLOAT(NSLICE2) 

           SCALET  = SCALEX
           IF (SCALEY .NE. SCALEX) SCALET = 0.0   ! X NOT SAME AS Y
           IF (SCALEZ .NE. SCALEX) SCALET = 0.0   ! Z NOT SAME AS Y
  
        ELSE
C          2D CASE
           NSLICE2 = 1
           LSD     = NSAM1+2-MOD(NSAM1,2)
           LSDN    = NSAM2+2-MOD(NSAM2,2)

           ALLOCATE (X(LSD,NROW1), Y(LSDN,NROW2), STAT=IRTFLG)
           IF (IRTFLG.NE.0) THEN 
              MWANT = LSD*NROW1 + LSDN*NROW2
              CALL ERRT(46,'FOUR1A_FP, X & Y',MWANT)
              GOTO 9000
          ENDIF 

           CALL READV(LUN1,X,LSD,NROW1,NSAM1,NROW1,NSLICE2)

           CALL FINT(X,Y,NSAM1,NROW1,NSAM2,NROW2,LSD,LSDN)

           CALL WRITEV(LUN2,Y,LSDN,NROW2,NSAM2,NROW2,NSLICE2)

C          IMAGE FILE HEADER FOR PIXSIZ HAS CHANGED
           SCALEX = FLOAT(NSAM1)   / FLOAT(NSAM2)
           SCALEY = FLOAT(NROW1)   / FLOAT(NROW2)
           SCALET = SCALEX
           IF (SCALEY .NE. SCALEX) SCALET = 0.0   ! X NOT SAME AS Y

        ENDIF

C       UPDATE THE INCORE HEADER VALUE & FILE HEADER FOR PIXSIZ
        CALL SETPRMS(LUN2, SCALET,IRTFLG)


9000   IF (NSLICE .EQ. 1) THEN
           IF (ALLOCATED(X)) DEALLOCATE (X)
           IF (ALLOCATED(Y)) DEALLOCATE (Y)
        ELSE
           IF (ALLOCATED(X3)) DEALLOCATE (X3)
           IF (ALLOCATED(Y3)) DEALLOCATE (Y3)
        ENDIF

9001    CLOSE (LUN2)
        CLOSE (LUN1)

        END
