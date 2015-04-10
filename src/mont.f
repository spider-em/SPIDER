
C++************************************************************** 3/4/81
C
C  MONT.F             LONG FILNAMES               JAN  89 ARDEAN LEITH
C                     REWRITTEN                   MAY  01 ARDEAN LEITH
C                     ALLOC                       FEB  12 ARDEAN LEITH
C                     DEALLOC                     MAR  15 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2015  Health Research Inc.,                         *
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
C PURPOSE: MAKE MONTAGE FROM FILE SERIES
C
C***********************************************************************

	SUBROUTINE MONT(MAXDIM)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

        REAL, ALLOCATABLE       :: BUF(:)

        LOGICAL                 :: NORMLZ, PUTOUT
        CHARACTER(LEN=MAXNAM)   :: FILNAM,FILOUT,FILPAT
        CHARACTER(LEN=1)        :: NULL = CHAR(0)

        INTEGER,PARAMETER       :: LUN  = 20
        INTEGER,PARAMETER       :: LUNO = 21

C       GET FILENAME PATTERN, NUMBERS
        NILMAX = NIMAX
        CALL FILSEQP(FILPAT,NLET,INUMBR,NILMAX,NUMTOT, 
     &     'FILE PREFIX OR TEMPLATE (EG. PIC****)' ,IRTFLG)
     
        IF (IRTFLG .NE. 0) RETURN

C       CREATE FIRST FILENAME
        CALL FILGET(FILPAT,FILNAM,NLET,INUMBR(1),IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       OPEN FIRST FILE
        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'O',IFORM,NX,NY,NZ,
     &               MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        IF (IFORM == 3) THEN
           WRITE(NOUT,*) 'WARNING: CAN ONLY USE FIRST SLICE OF VOLUME'
        ENDIF

        IF (IMAMI.NE.1) CALL NORM3(LUN,NX,NY,NZ,FMAX,FMIN,AV)
        WRITE(NOUT,97) FMIN,FMAX
97      FORMAT('  FIRST IMAGE RANGE: ',G11.3,' ... ',G11.3,/)

	IF (FCHAR(4:4) == 'S') THEN
           NORMLZ = .TRUE.
           IF (IMAMI.NE.1) CALL NORM3(LUN,NX,NY,NZ,FMAX,FMIN,AV)
           FMINT = FMIN
           SCT   = 0.0
           IF (FMAX .NE. FMIN) SCT = 2.0 / (FMAX-FMIN)
        ELSE
           NORMLZ = .FALSE.
        ENDIF

	CALL RDPRI2S(NPR,MAR,NOT_USED,
     &      'NUMBER OF IMAGES PER ROW, MARGIN WIDTH',IRTFLG)

        BACK = 1.0
10      CALL RDPRM1S(BACK,NOT_USED,'MARGIN VALUE',IRTFLG)
	IF (FCHAR(4:4) == 'S' .AND. 
     &     (BACK < 0.0 .OR. BACK > 2.0)) THEN
           CALL ERRT(101,'MARGIN VALUE MUST BE IN RANGE 0 ... 2',NE)
           GOTO 10
        ENDIF

C       FIND SIZE OF OUTPUT MONTAGE (MARGINS ALL AROUND)
	NXO    = (NX * NPR) + MAR *(NPR + 1)
	NUMROW = NUMTOT / NPR + 1
	IF ((NUMTOT / NPR) * NPR == NUMTOT) NUMROW = NUMTOT / NPR
	NYO     = NY * NUMROW + MAR * ( NUMROW + 1)
        NUMMARS = NUMROW + 1

        MWANT = NXO * (NYO + MAR)
        ALLOCATE(BUF(MWANT), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'MONT; BUF',MWANT)
           GOTO 999
        ENDIF
                  
C       OPEN OUTPUT FILE
        IFORM = 1
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILOUT,LUNO,'N',IFORM,NXO,NYO,1,
     &             MAXIM,'OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999
          
C       FILL HORIZONTAL MARGIN ROWS IN FILE
        IROW = 1
        IF (MAR > 0) THEN
C          FILL BUFFER WITH BACKGROUND FIRST
           BUF = BACK

C          PLACE MARGIN ROWS IN FILE
           DO M = 1, NUMMARS
C             PLACE MAR ROWS OF BUFFER IN FILE
              DO I = 1,MAR
                 CALL WRTLIN(LUNO,BUF(1),NXO,IROW)
                 IROW = IROW + 1
              ENDDO
             IROW = IROW + NY
           ENDDO
        ENDIF

C       READ IN IMAGE ROWS               
        IPICONY = NPR
        IROWO     = MAR 

        DO IFIL = 1,NUMTOT
C          FIRST FILE ALREADY OPEN!
	   IF (IFIL > 1) THEN
C             OPEN THE NEW INPUT FILE
              CLOSE(LUN)
              CALL FILGET(FILPAT,FILNAM,NLET,INUMBR(IFIL),IRTFLG)
              MAXIM = 0
              CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'O',IFORM,
     &             NXN,NYN,NZ,
     &             MAXIM,'UNUSED',.FALSE.,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 998

              CALL SIZCHK(UNUSED,NX, NY, 1,1,
     &                           NXN,NYN,1,1,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 998

              IF (NORMLZ) THEN
                 IF (IMAMI .NE. 1)CALL NORM3(LUN,NX,NY,NZ,FMAX,FMIN,AV)
                 FMINT = FMIN
                 SCT   = 0.0
                 IF (FMAX .NE. FMIN) SCT = 2./(FMAX-FMIN)
              ENDIF
	   ENDIF

	   IPICONY = IPICONY + 1
C          START NEW ROW IF NECESSARY
	   IF (IPICONY > NPR)  IPICONY = 1

           PUTOUT = .TRUE.

C          PUT IN VERTICAL LEFT MARGIN
	   NSTART = 1
	   DO IROW = 1, NY
              DO ISAM = NSTART, NSTART + MAR - 1
                 BUF(ISAM) = BACK
              ENDDO
              NSTART = NSTART + NXO
	   ENDDO
      
C          READ ROWS FROM FILE
	   DO  IROW = 1,NY
C             COMPUTE CURRENT BUFFER POSITION
              NSTART = (IPICONY - 1) * (NX + MAR) + MAR + 1 +
     &                 (IROW -1) * NXO
              CALL REDLIN(LUN,BUF(NSTART),NX,IROW)
              IF (NORMLZ) THEN
                 DO  K = NSTART ,NSTART+NX-1
                    BUF(K) = (BUF(K) - FMINT) * SCT
                 ENDDO
              ENDIF

C             ADD VERTICAL MARGIN TO RIGHT
              DO J = NSTART+NX, NSTART+NX+MAR
                 BUF(J) = BACK
              ENDDO
	   ENDDO

           IF (IPICONY == NPR) THEN
C             OUTPUT THIS ROW OF IMAGES
              ISTART =  1
              DO IROWT = 1, NY
	         CALL WRTLIN(LUNO,BUF(ISTART),NXO,IROWO + IROWT)
                 ISTART = ISTART + NXO
              ENDDO
              IROWO  = IROWO + NY + MAR 
              PUTOUT = .FALSE.
           ENDIF
        ENDDO

        IF (PUTOUT) THEN
C          OUTPUT REMAINING PARTIAL ROW OF IMAGES

C          MUST FILL REMAINING PART OF BUFFER WITH BACKGROUND
           ISTART =  IPICONY  * (NX + MAR) + MAR + 1
           NVALS  = NXO - ISTART
           DO IROWT = 1, NY 
              DO ISAMO = ISTART, ISTART + NVALS
	         BUF(ISAMO) = BACK
              ENDDO
              ISTART = ISTART + NXO
           ENDDO

C          PUT OUT THE ROW OF IMAGES 
           ISTART = 1
           DO IROWT = 1, NY
	      CALL WRTLIN(LUNO,BUF(ISTART),NXO,IROWO + IROWT)
              ISTART = ISTART + NXO
           ENDDO
        ENDIF

998     CLOSE(LUNO)
999     CLOSE(LUN)

        IF (ALLOCATED(BUF))    DEALLOCATE(BUF)

	END

