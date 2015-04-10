
C **********************************************************************
C
C ADS.F                            FOR SPEED     APR 03 ARDEAN LEITH
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
C  PURPOSE:   ADD A SERIES OF IMAGES (USES LESS MEMORY THAN 'AS R'
C             FASTER THAN 'AD'
C
C--*********************************************************************

      SUBROUTINE ADS(LUNIN,LUNSUM,LUNDOC)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=MAXNAM) ::    FILNAM,FILA,FILPAT
      COMMON /COMMUN/             FILNAM,FILA,FILPAT

      COMMON /IOBUF/ BUF(NBUFSIZ)

      REAL, ALLOCATABLE, DIMENSION(:) :: SUMARAY
#ifndef SP_32
      INTEGER *8       NVOX,IOK8
#else
      INTEGER *4       NVOX,IOK8
#endif
      
      CALL FILELIST(.TRUE.,LUNDOC,FILPAT,NLETP,INUMBR,NIMAX,NUMT,
     &                 'INPUT FILE TEMPLATE (E.G. PIC****)',IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     OPEN FIRST PICTURE TO DETERMINE DIMS
      CALL FILGET(FILPAT,FILNAM,NLETP,INUMBR(1),IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      MAXIM = 0
      CALL OPFILEC(0,.FALSE.,FILNAM,LUNIN,'O',IFORM,NSAM,NROW,NSLICE,
     &               MAXIM,' ',.TRUE.,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      NVOX = NSAM * NROW * NSLICE
C     COMPLAIN IF EXCESSIVE ALLOCATION
      CALL BIGALLOC(NVOX,IOK8,.FALSE.,.TRUE.,IRTFLG)

      ALLOCATE(SUMARAY(NVOX),STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
          CALL ERRT(102,'FAILED TO ALLOCATE',NVOX3)
          GOTO 999
      ENDIF

C     FILL THE SUM ARRAY WITH FIRST FILE
      CALL REDVOL(LUNIN,NSAM,NROW,1,NSLICE,SUMARAY,IRTFLG)
      CLOSE(LUNIN)

      NREC = NROW * NSLICE

      DO IFIL=2,NUMT

         CALL FILGET(FILPAT,FILNAM,NLETP,INUMBR(IFIL),IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(3,'ADS',NE)
            GOTO 997
         ENDIF

         MAXIM = 0
         CALL OPFILEC(0,.FALSE.,FILNAM,LUNIN,'Z',IFORM,
     &                NSAM1,NROW1,NSLICE1,MAXIM,' ',.TRUE.,IRTFLG)
         IF (IRTFLG .NE. 0) THEN
C           ALLOW GAPS IN FILE SERIES
            WRITE(NOUT,100) FILNAM
100         FORMAT(' FILE: ',A,' NOT FOUND, FILE SKIPPED')
         
         ELSE
C           INPUT FILE OPENED OK
            IF (NSAM1.NE.NSAM.OR.NROW1.NE.NROW.OR.NSLICE1.NE.NSLICE)THEN
               CALL ERRT(1,'ADS',NE)
               GOTO 997
            ENDIF

            ILOC = 0            
            DO I=1,NREC
               CALL  REDLIN(LUNIN,BUF,NSAM,I)
               DO J = 1,NSAM
                  SUMARAY(ILOC+J) = SUMARAY(ILOC+J) + BUF(J)
               ENDDO
               ILOC = ILOC + NSAM
            ENDDO
            CLOSE(LUNIN)
         ENDIF
      ENDDO

C     FINISHED, OPEN OUTPUT FILE
      MAXIM = 0
      CALL OPFILEC(LUNIN,.TRUE.,FILA,LUNSUM,'U',IFORM,NSAM,NROW,NSLICE,
     &                      MAXIM,'OUTPUT',.TRUE.,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 997

      CALL WRTVOL(LUNSUM,NSAM,NROW,1,NSLICE,SUMARAY,IRTFLG)

997   CLOSE(LUNSUM)
      CLOSE(LUNIN)

999   IF (ALLOCATED(SUMARAY)) DEALLOCATE(SUMARAY)

      RETURN
      END

