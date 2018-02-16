C++*********************************************************************
C
C    FOUR1B           OPFILEC                    FEB 2003 ARDEAN LEITH    
C                     REFACTORED                 FEB 2018 ARDEAN LEITH    
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2018  Health Research Inc.,                         *
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
C  PURPOSE: CALCULATES POWER SPECTRUM ------------------------ 'PW'
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

	SUBROUTINE FOUR1B

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: QA
        REAL, ALLOCATABLE, DIMENSION(:,:)   :: qt

        CHARACTER (LEN=MAXNAM) ::  FILNAM

        INTEGER                :: NX,NY,NZ,MAXIM,IRTFLG,IFORMIN,LSD,IRL
        INTEGER                :: NE,INV

        INTEGER, PARAMETER     :: LUN1 = 21
        INTEGER, PARAMETER     :: LUN2 = 22

C       FOURIER MODULI FROM COMPLEX FOURIER TRANSFORM 
C       ON DISK FOR 2-D OR 3-D PICTURES.

C       OPEN INPUT IMAGE/VOLUME, MAY BE FOURIER FORMAT
        MAXIM = 0
	CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NX,NY,NZ,
     &		   MAXIM,'INPUT',.TRUE.,IRTFLG)
	IF (IRTFLG .NE. 0) RETURN

C       FOURIER INPUT 2D & 3D
        NZ      = MAX(1,NZ)
	IFORMIN = IFORM

        IF (IFORM == 1)  THEN
C          REAL IMAGE
           IF(MOD(NX,2) == 0)  THEN
              LSD = NX + 2
           ELSE
              LSD = NX + 1
           ENDIF
           IRL = NX
          
        ELSEIF (IFORM == 3)  THEN
C          REAL VOLUME
           IF (MOD(NX,2) == 0)  THEN
              LSD = NX + 2
           ELSE
              LSD = NX + 1
           ENDIF
           IRL = NX
          
        ELSEIF (IFORM == -11)  THEN
C          FOURIER IMAGE, ODD X
           IFORM = 1
           LSD   = NX
           IRL   = NX
           NX    = NX - 1
          
        ELSEIF (IFORM == -12)  THEN
C          FOURIER IMAGE, EVEN X
           IFORM = 1
           LSD   = NX
           IRL   = NX
           NX    = NX - 2
           
        ELSEIF (IFORM == -21)  THEN
C         FOURIER VOLUME, ODD X
           IFORM = 3
           LSD   = NX
           IRL   = NX
           NX    = NX - 1
           
        ELSEIF(IFORM == -22)  THEN
C          FOURIER VOLUME, EVEN X
           IFORM = 3
           LSD   = NX
           IRL   = NX
           NX    = NX - 2
           
        ELSE
           CALL ERRT(102,'UNKNOWN FILE FORMAT',IFORM)
           CLOSE(LUN1)
           RETURN
        ENDIF

C       OPEN REAL OUTPUT FILE FOR POWER SPECTRUM
        MAXIM = 0
	CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'N',IFORM,NX,NY,NZ,
     &		       MAXIM,'OUTPUT',.FALSE.,IRTFLG)
	IF (IRTFLG .NE. 0)  RETURN

       !ALLOCATE (qt(nx,ny), QA(LSD,NY,NZ), STAT=IRTFLG)
 	ALLOCATE (QA(LSD,NY,NZ), STAT=IRTFLG)

        IF (IRTFLG  ==  0) THEN 
C          ADEQUATE SPACE IN MEMORY, LOAD INPUT IMAGE/VOLUME          
	   CALL READV(LUN1,QA,LSD,NY,IRL,NY,NZ)

	   IF (IFORMIN  >  0) THEN
C             REAL INPUT 2D & 3D, CONVERT TO FOURIER
	      INV = +1
              IF (NZ == 1)  THEN
	         CALL FMRS_2(QA,NX,NY,INV)
	      ELSE
	         CALL FMRS_3(QA,NX,NY,NZ,INV)
	      ENDIF
	      IF (INV == 0)  THEN
                 CALL  ERRT(101,'IN FFT',NE)
                 GOTO 999
	      ENDIF
	   ENDIF
           
           IF (NZ == 1)  THEN
C             IMAGE

              IF (INDEX(FCHAR,'T') > 1 ) THEN
                 CALL PW2SR  (QA,NX,NY,FCHAR(4:4))
c                OUTPUT REAL IMAGE	   
                 CALL WRITEV(LUN2,QA,LSD,NY,NX,NY,1)
                !qt(1:nx,1:ny) = qa(1:nx,1:ny,1)
                !call chkmaxloc2d(' max loc qt:',qt,nx,ny)

              ELSE
                CALL PW2SR_A(QA,NX,LSD,NY,FCHAR(4:4))

C                OUTPUT REAL IMAGE	   
                 CALL WRITEV(LUN2,QA,LSD,NY,NX,NY,1)
              ENDIF 

	   ELSE
C             VOLUME
              CALL PW3SR(QA,NX,NY,NZ,FCHAR(4:4))

C             OUTPUT REAL VOLUME	   
              CALL WRITEV(LUN2,QA,LSD,NY,NX,NY,NZ)
           ENDIF


        ELSEIF (IFORMIN > 0)  THEN
C          REAL INPUT FOR ON DISK VERSION NOT SUPPORTED
           CALL ERRT(101,'DISK VERSION NEEDS FOURIER INPUT',NE)
           GOTO 999

	ELSE
C          ON DISK VERSION, ANCIENT LEGACY VERSION
           WRITE(NOUT,*)' ** WARNING: SLOW ON-DISK VERSION USED.'

	   IF (IFORM == 1)  THEN
C             IMAGE
	      CALL PW2SDR(LUN1,LUN2,NX,NY,FCHAR(4:4))
	   ELSE
C             VOLUME
	      CALL PW3SDR(LUN1,LUN2,NX,NY,NZ,FCHAR(4:4))
	   ENDIF

	ENDIF

999     CONTINUE
	CLOSE(LUN1)
	CLOSE(LUN2)
        DEALLOCATE(QA)  
      
        END
