C++*********************************************************************
C
C  VAR3R.F                                        06/03/02
C                                         
C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2002, P. A. Penczek                                   *
C=* University of Texas - Houston Medical School                       *
C=* Email:  pawel.a.penczek@uth.tmc.edu                                *
C=*                                                                    *
C=* This program is free software; you can redistribute it and/or      *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* This program is distributed in the hope that it will be useful,    *
C=* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C=* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
C=* General Public License for more details.                           *
C=*                                                                    *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program; if not, write to the                      *
C=* Free Software Foundation, Inc.,                                    *
C=* 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.      *
C=*                                                                    *
C=**********************************************************************
C
C VAR3R  
C
C PURPOSE:  CALCULATE REL VARIANCE FOR A SET OF 3D VOLUMES
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE VAR3R

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER (LEN=MAXNAM) :: FILNAM,FILA,FILV,FILPAT
#ifndef SP_32
      INTEGER *8             :: NVOX,IOK8
#else
      INTEGER *4             :: NVOX,IOK8
#endif
      REAL, ALLOCATABLE      :: AVGARAY(:)
      REAL, ALLOCATABLE      :: VARARAY(:)
      REAL, ALLOCATABLE      :: RELVAR(:)
      REAL, ALLOCATABLE      :: INPARAY(:)

      INTEGER, PARAMETER     :: LUN1 = 21
         
      CALL FILELIST(.TRUE.,NDOC,FILPAT,NLETP,INUMBR,NIMAX,NUMT,
     &                 'INPUT FILE TEMPLATE',IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET FIRST PICTURE TO DETERMINE DIMS
      CALL FILGET(FILPAT,FILNAM,NLETP,INUMBR(1),IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      MAXIM = 0
      CALL OPFILEC(0,.FALSE.,FILNAM,LUN1,'O',IFORM,NSAM,NROW,NSLICE,
     &               MAXIM,' ',.TRUE.,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      CLOSE(LUN1)

      NVOX = NSAM * NROW * NSLICE
      CALL BIGALLOC(NVOX,IOK8,.FALSE.,.TRUE.,IRTFLG)

      ALLOCATE(AVGARAY(NVOX),VARARAY(NVOX),
     &	       RELVAR(NVOX), INPARAY(NSAM),STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
          MWANT = 3*NVOX + NSAM
          CALL ERRT(46,'VAR3R, AVGARAY...',MWANT)
          RETURN
      ENDIF

      AVGARAY = 0.0
      NUMA    = 0

      DO IFIL=1,NUMT

         CALL FILGET(FILPAT,FILNAM,NLETP,INUMBR(IFIL),IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 997

         CALL OPFILEC(0,.FALSE.,FILNAM,LUN1,'Z',IFORM,
     &                   NSAM,NROW,NSLICE,MAXIM,' ',.TRUE.,IRTFLG)
         IF (IRTFLG .EQ. 0) THEN
            NUMA = NUMA + 1
	    DO  I=1,NROW*NSLICE
		CALL  REDLIN(LUN1,INPARAY,NSAM,I)
	        AVGARAY(1+(I-1)*NSAM:I*NSAM) =
     &			AVGARAY(1+(I-1)*NSAM:I*NSAM) + INPARAY
	    ENDDO 
	    CLOSE(LUN1)
	 ENDIF
       ENDDO

C     COMPUTE AVERAGE

      AVGARAY = AVGARAY / NUMA
C     OPEN NEW AVERAGE OUTPUT FILE
      CALL OPFILEC(0,.TRUE.,FILA,LUN1,'U',IFORM,NSAM,NROW,NSLICE,
     &                      MAXIM,'AVERAGE',.TRUE.,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 997

      CALL WRITEV(LUN1,AVGARAY,NSAM,NROW,NSAM,NROW,NSLICE)
      CLOSE(LUN1)

      VARARAY = 0.0
      RELVAR  = 0.0

      DO IFIL=1,NUMT

          CALL FILGET(FILPAT,FILNAM,NLETP,INUMBR(IFIL),IRTFLG)
          IF (IRTFLG .NE. 0) THEN
             CALL ERRT(3,'VA3R',NE)
             GOTO 997
          ENDIF

          CALL OPFILEC(0,.FALSE.,FILNAM,LUN1,'Z',IFORM,
     &                   NSAM,NROW,NSLICE,MAXIM,' ',.TRUE.,IRTFLG)
          IF (IRTFLG .EQ. 0) THEN
	     DO  I=1,NROW*NSLICE
	        CALL  REDLIN(LUN1,INPARAY,NSAM,I)

	        VARARAY(1+(I-1)*NSAM:I*NSAM) =
     &                VARARAY(1+(I-1)*NSAM:I*NSAM)+
     &		      (AVGARAY(1+(I-1)*NSAM:I*NSAM)-INPARAY)**2

	        RELVAR(1+(I-1)*NSAM:I*NSAM) =
     &                RELVAR(1+(I-1)*NSAM:I*NSAM)+
     &		      (AVGARAY(1+(I-1)*NSAM:I*NSAM)-INPARAY)**4
	     ENDDO 
	     CLOSE(LUN1)
	  ENDIF
      ENDDO

      VARARAY = VARARAY / NUMA

C     OPEN NEW VARIANCE OUTPUT FILE
      CALL OPFILEC(0,.TRUE.,FILV,LUN1,'U',IFORM,NSAM,NROW,NSLICE,
     &                   MAXIM,'VARIANCE',.TRUE.,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 997
      CALL WRITEV(LUN1,VARARAY,NSAM,NROW,NSAM,NROW,NSLICE)
      CLOSE(LUN1)

      RELVAR = RELVAR/NUMA

      DO  I=1,NSAM*NROW*NSLICE
	 IF (VARARAY(I) .GT. 0.0)  THEN
	     RELVAR(I) = 
     &            (RELVAR(I)/VARARAY(I)**2-(NUMA-3)/REAL(NUMA-1))/NUMA
	 ELSE
	     RELVAR(I) = 0.0
	 ENDIF
      ENDDO

C     OPEN NEW REL VARIANCE OUTPUT FILE
      CALL OPFILEC(0,.TRUE.,FILV,LUN1,'U',IFORM,NSAM,NROW,NSLICE,
     &               MAXIM,'RELVARIANCE',.TRUE.,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 997

      CALL WRITEV(LUN1,RELVAR,NSAM,NROW,NSAM,NROW,NSLICE)
      CLOSE(LUN1)

997   IF (ALLOCATED(INPARAY)) DEALLOCATE(INPARAY)
      IF (ALLOCATED(RELVAR))  DEALLOCATE(RELVAR)
      IF (ALLOCATED(VARARAY)) DEALLOCATE(VARARAY)
      IF (ALLOCATED(AVGARAY)) DEALLOCATE(AVGARAY)

      END
