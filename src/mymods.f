
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
C++************************************************************************
C
C   MYMODS.F
C
C   A SHELL FOR USER SUPLIED SUBROUTINES, INCLUDES A COMMENTED SAMPLE
C 
C--************************************************************************

	SUBROUTINE MYMODS(MAXMEM)

C       INCLUDE COMMONS WITH SPIDER SYSTEM INFORMATION      
	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        CHARACTER(LEN=MAXNAM)   ::  FILNAM

        REAL, ALLOCATABLE, DIMENSION(:) ::  VOL1,VOL2

C       DEFINE IO UNITS, FOR IMAGE FILES USE IOUNITS 8...99 ONLY
	DATA            LUN1,LUN2/8,9/

C       FILE OPENING PART --------------------------------------------

C       OPEN THE INPUT FILE
C       NSAM - X-DIMENSION, NROW - Y-DIMENSION, NSLICE - Z-DIMENSION
C       'O' - INDICATES EXISTING FILE, .FALSE. IS FOURIER INPUT NOT OK
C       MAXIM > 0 - INDICATES THAT OPERATION WAORKS ON WHOLE STACKS
         
        MAXIM = 0
 	CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,
     &              NSAM,NROW,NSLICE,MAXIM,'INPUT',.FALSE.,IRTFLG)

C       CHECK WHETHER THE FILE OPENED OK, IRTFLG<>0 INDICATES IT WAS NOT
	IF (IRTFLG .NE. 0)  GOTO 999

C       OPEN THE OUTPUT FILE, SAME DIMENSIONS AS INPUT FILE
C       'UNKNOWN' TELLS THE PROGRAM TO OVERWRITE ANY EXISTING FILE

        MAXIM = 0
 	CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',IFORM,
     &              NSAM,NROW,NSLICE,MAXIM,'OUTPUT',.FALSE.,IRTFLG)

C       CHECK WHETHER FILE OPENED OK, IRTFLG<>0 INDICATES IT WAS NOT
	IF (IRTFLG .NE. 0)  GOTO 999

C       MEMORY ALLOCATION PART ---------------------------------------

C       THIS OPERATION NEEDS TWO VOLUMES, EACH NSAM*NROW*NSLICE LARGE.
	MEMWANT = NSAM * NROW * NSLICE

        ALLOCATE (VOL1(MEMWANT),VOL2(MEMWANT), STAT=IRTFLG)

C       CHECK WHETHER MEMORY OK, IRTFLG<>0 INDICATES IT WAS NOT
        IF (IRTFLG .NE. 0) THEN
C          ERRT WILL PRINT ERROR MESSAGE AND STOP IF IN BATCH MODE
           CALL ERRT(46,'MYMODS, VOL1 & VOL2',NDUM)
           GO TO 999
        ENDIF

C       SPECIFIC PROCESSING PART -------------------------------------

C       SOLICIT FLOATING POINT NUMERICAL PARAMETERS USED BY USONE 
        CALL  RDPRM2(CW,SW,NOT_USED,'CENTRAL WEIGHT, SIDE WEIGHT')

C       CALL THE SUBROUTINE THAT WILL DO THE ACTUAL PROCESSING
        CALL USONE(VOL1,LUN1,VOL2,LUN2,NSAM,NROW,NSLICE,CW,SW)

C       CLEANUP ------------------------------------------------------

C       DEALLOCATE RUN-TIME MEMORY
        IF(ALLOCATED(VOL1)) DEALLOCATE(VOL1)
        IF(ALLOCATED(VOL2)) DEALLOCATE(VOL2)

C       CLOSE LOGICAL IO UNITS.
999     CLOSE(LUN1)
        CLOSE(LUN2)
        RETURN
	END




C       ----------------------------------------------------------


	SUBROUTINE  USONE(BIN,LUN1,BOU,LUN2,NSAM,NROW,NSLICE,CW,SW)

	DIMENSION  BIN(NSAM,NROW,NSLICE),BOU(NSAM,NROW,NSLICE)
 
C       THIS PROGRAM APPLIES 1D MOVING AVERAGE IN X-DIRECTION (NSAM)
C       THE RESULT IS "DIRECTIONAL" FILTRATION IN ONE DIMENTION

C       BOU(I,*,*) = (CW*BIN(I,*,*) + SW * (BIN(I-1,*,*) +
C                     BIN(I+1,*,*))) / (CW+SW+SW)

C       THIS PROGRAM WORKS FOR BOTH 2D AND 3D FILES.  THIS IS DUE TO 
C       THE WAY 3D FILES ARE ORGANIZED IN SPIDER - THEY HAVE NSLICE 
C       2D SLICES. A 2D FILE IS A 3D FILE WITH NSLICE=1.

C       READ INPUT VOLUME
        DO K=1,NSLICE
           DO J=1,NROW
              IREC=J+(K-1)*NROW
              CALL  REDLIN(LUN1,BIN(1,J,K),NSAM,IREC)
           ENDDO
        ENDDO


C       DEFINE CENTER OF THE VOLUME ACCORDING TO SPIDER CONVENTION, NOT
C       USED IN THIS PROGRAM, JUST AN EXAMPLE.
C	NXC = NSAM   / 2 + 1
C	NYC = NROW   / 2 + 1
C	NZC = NSLICE / 2 + 1

	Q = 1.0 / (CW + SW + SW)
C       DO THE PROCESSING

C       FOLLOWING LINE IS A COMPILER PARALLEL DIRECTIVE FOR OPEN-MP.
c$omp parallel do private(I,J,K)
	DO K=1,NSLICE
	   DO J=1,NROW
C            NEXT TWO LINES TAKE CARE OF BORDER EFFECT.  IF THE IMAGE 
C            IS NOT TO BE TREATED AS CIRCULARLY CLOSED, THEY SHOULD 
C            BE COMMENTED OUT.
	     BOU(1,J,K) = (CW * BIN(1,J,K) + SW * (BIN(NSAM,J,K)+
     &                     BIN(2,J,K))) * Q
	     BOU(NSAM,J,K) = (CW*BIN(NSAM,J,K) + SW * (BIN(NSAM-1,J,K)+
     &                       BIN(1,J,K)))*Q

C            HERE IT IS REQUIRED THAT NSAM IS AT LEAST 3.  IT IS 
C            NOWHERE CHECKED, BUT IT CAN BE EASILY ADDED.
	     DO I=2,NSAM-1
                BOU(I,J,K) = (CW * BIN(I,J,K) + SW * (BIN(I-1,J,K) +
     &                     BIN(I+1,J,K))) * Q
             ENDDO
          ENDDO
        ENDDO

C       WRITE OUTPUT VOLUME
        DO K=1,NSLICE
           DO J=1,NROW
              IREC = J + (K - 1) * NROW
              CALL  WRTLIN(LUN2,BOU(1,J,K),NSAM,IREC)
           ENDDO
        ENDDO

C       IO FILES WILL BE CLOSED IN THE CALLING PROGRAM 
        RETURN
	END
