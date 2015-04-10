

C++*********************************************************************
C
C STACK.F      LONG FILENAMES                      JAN 89 ARDEAN LEITH
C              REWRITE                             MAY 13 ArDean Leith
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
C    STACK()
C
C    PURPOSE:  STACK 2-D SLICES INTO 3-D IMAGE
C              CAN OPERATE ON IMAGE SERIES
C
C    NOTE:     OBSOLETE   USE 'CP TO VOL' INSTEAD!!
C
C--*******************************************************************

 	SUBROUTINE STACK()

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 

        CHARACTER(LEN=MAXNAM) :: FILNM1,FILNM2,FILENM,FILPAT,FILDUM
        CHARACTER(LEN=1)      :: NULL = CHAR(0)

        INTEGER,PARAMETER     :: LUNI = 21 
        INTEGER,PARAMETER     :: LUNO = 23

        COMMON /IOBUF/ BUF(NBUFSIZ)

C       USE FILE OPEN TO FIND VALUES TO USE FOR NX, NY
        MAXIM  = 0
	IFOUND = -4
        CALL OPFILEC(0,.TRUE.,FILNM1,LUNI,'O',IFORM,NX,NY,NZ,
     &               MAXIM,'FIRST',.FALSE.,IRTFLG)
	IFOUND = 0
        IF (IRTFLG .NE. 0) RETURN
	CLOSE(LUNI)
        NLET = lnblnk(FILNM1)
        IF (NLET .LE. 0) NLET = LEN(FILNM1)

C       FIND NUMBER OF FIRST FILE
	CALL FILCAD(FILNM1,FILPAT,N1,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       FIND NUMBER OF LAST FILE
	CALL FILERD(FILNM2,NLET,NULL,'LAST',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
	CALL FILCAD(FILNM2,FILDUM,N2,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        NFILES = N2 - N1 + 1

        IF (NZ .LE. 1) THEN
           NZS = NFILES
        ELSE
C          FIND TOTAL NUMBER OF SLICES
           NZS = 0
	   DO  IFILE=1,NFILES

C            CREATE CURRENT INPUT FILE NAME
             CALL FILGET(FILPAT,FILNM2,NLET,N1+IFILE-1,IRTFLG)
             IF (IRTFLG .NE. 0) THEN
                CALL ERRT(2,'STACK',NE)
                GOTO 9900
             ENDIF

C            OPEN INPUT FILE TO GET NZT
             MAXIM  = 0
	     IFOUND = -4
             CALL OPFILEC(0,.FALSE.,FILNM2,LUNI,'O',IFORM,
     &                NXT,NYT,NZT,MAXIM,'DUMMY',.FALSE.,IRTFLG)
	     IFOUND = 0
             IF (IRTFLG .NE. 0) GOTO 9900
	     IF (NXT.NE.NX .OR. NYT.NE.NY) THEN
                CALL ERRT(1,'STACK ',NE)
                GOTO 9900
             ENDIF
             CLOSE(LUNI)

             NZS = NZS + NZT
           ENDDO
        ENDIF

C       OPEN OUTPUT VOLUME
	IFORM  = 3
        MAXIM  = 0
        CALL OPFILEC(0,.TRUE.,FILENM,LUNO,'U',IFORM,NX,NY,NZS,
     &                   MAXIM,'OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        WRITE(NOUT,*) ' '

        IRECOUT = 0
	DO  IFILE=1,NFILES

C         CREATE CURRENT INPUT FILE NAME
          CALL FILGET(FILPAT,FILNM2,NLET,N1+IFILE-1,IRTFLG)
          IF (IRTFLG .NE. 0) THEN
             CALL ERRT(2,'STACK',NE)
             GOTO 9900
          ENDIF

C         OPEN INPUT FILE
          CALL OPFILEC(0,.FALSE.,FILNM2,LUNI,'O',IFORM,
     &                NXT,NYT,NZT,MAXIM,'DUMMY',.FALSE.,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9900
	  IF (NXT.NE.NX .OR. NYT.NE.NY) THEN
             CALL ERRT(1,'STACK ',NE)
             GOTO 9900
          ENDIF

          DO  IRECIN=1,NY*NZT
             CALL REDLIN(LUNI,BUF,NX,IRECIN)
             IRECOUT = IRECOUT + 1
	     IF (IRECOUT .GT. NY * NZS) THEN
                CALL ERRT(102,'RECORD LIMIT (NY*NZ)',NY*NZ)
                GOTO 9900
             ENDIF
             CALL WRTLIN(LUNO,BUF,NX,IRECOUT)
	  ENDDO
          CLOSE(LUNI)
	ENDDO


9900    CLOSE(LUNI)
        CLOSE(LUNO)

	END



C     --------------------- STACK_REPLACE ----------------------------

      SUBROUTINE STACK_REPLACE()

      IMPLICIT NONE
      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      REAL,    ALLOCATABLE   :: BUFIMG(:)

      CHARACTER (LEN=MAXNAM) :: FILOUT,FILIN

      INTEGER                :: ITYPE,NX,NY,NZ,NXI,NYI,NZI
      INTEGER                :: IZ,NPIX,NDUM,NOT_USED
      INTEGER                :: ISLICE,MAXIM,IRTFLG

      INTEGER,PARAMETER      :: LUNIN     = 21 
      INTEGER,PARAMETER      :: LUNOUT    = 23

C     OPEN EXISTING OUTPUT (3D) FILE
      MAXIM = 0
      CALL OPFILEC(0,.TRUE.,FILOUT,LUNOUT,'O',
     &               ITYPE,NX,NY,NZ,
     &               MAXIM,'OUTPUT VOLUME',.FALSE.,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      NPIX = NX * NY
      ALLOCATE(BUFIMG(NPIX),STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN 
         CALL ERRT(46,'STACK, BUFIMG',NPIX)
         GOTO 9999
      ENDIF

      IF (VERBOSE) WRITE(NOUT,*) ' '

      ISLICE = 0

      DO                   ! ENDLESS LOOP

C        OPEN FIRST OR NEXT INPUT FILE
         MAXIM = 0
         CALL OPFILEC(0,.TRUE.,FILIN,LUNIN,'O',ITYPE,NXI,NYI,NZI,
     &                MAXIM,'INPUT',.FALSE.,IRTFLG)

C        CONTINUE UNTIL ASTERICK RECEIVED FOR INPUT FILENAME
         IF (IRTFLG .NE. 0) GOTO 9999

         IF (NXI.NE.NX .OR. NYI.NE.NY) THEN
            CALL ERRT(101,'INCOMPATIBLE SIZES',NDUM)
            GOTO 9999
         ENDIF

        ISLICE = ISLICE + 1
        CALL RDPRI1S(ISLICE,NOT_USED,'SLICE NUMBER',IRTFLG)

        IF (ISLICE <= 0 .OR. ISLICE > NZ) THEN
           CALL ERRT(102,'SLICE OUT OF RANGE',NZ)
           GOTO 9999
        ENDIF

        CALL REDVOL(LUNIN,NX,NY,1,1,BUFIMG,IRTFLG)

        CALL WRTVOL(LUNOUT,NX,NY,ISLICE,ISLICE,BUFIMG,IRTFLG)

        CLOSE(LUNIN)

C       CONTINUE UNTIL ASTERICK RECEIVED FOR INPUT FILENAME
      ENDDO

9999  IF (ALLOCATED(BUFIMG)) DEALLOCATE(BUFIMG)

      IF (VERBOSE) WRITE(NOUT,*) ' '

      CLOSE(LUNIN)
      CLOSE(LUNOUT)
      
      END






















