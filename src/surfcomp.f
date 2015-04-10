C **********************************************************************
C
C   SURFCOMP.F   -- CREATED                        MAR 00
C                   MAXNAM                         JUL 14 ARDEAN LEITH
C **********************************************************************
C=*  AUTHOR: ArDean Leith                                              *
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C  SURFCOMP()
C
C  PURPOSE:     READS TWO SPIDER BINARY VOLUMES CONTAINING
C               SURFACE VOXELS FROM A SPIDER VOLUME, (THE SURFACE
C               VOXELS HAVE VALUES GREATEER THAN ZERO).
C               FINDS DIFFERENCE STATISTICS FOR THE TWO SURFACES.     
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE SURFCOMP()

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL,    ALLOCATABLE, DIMENSION(:)   :: BUFV1,BUFV2,BUFS1,BUFS2
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISURF1,ISURF2
        REAL,    ALLOCATABLE, DIMENSION(:)   :: RSURF1,RSURF2

        DOUBLE PRECISION :: DISTMINAV,DISTMINAV2
        DOUBLE PRECISION :: DSIG
        REAL, PARAMETER  :: BACKLEVEL = 0.0
        REAL, PARAMETER  :: SURFLEVEL = 1.0
        REAL, PARAMETER  :: FLTZER    = 10E-10

        CHARACTER(LEN=MAXNAM) :: FILNAM

        LUNIMV1   = 11
        LUNIMS1   = 12
        LUNIMV2   = 13
        LUNIMS2   = 14

C       OPEN FIRST SPIDER VOLUME AS INPUT
        MAXIM    = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNIMV1,'O',IFORM,NSAM,NROW,NSLICE,
     &             MAXIM,'FIRST INPUT VOLUME',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       OPEN FIRST SPIDER SURFACE VOLUME AS INPUT
        MAXIM    = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNIMS1,'O',IFORM,
     &             NSAM1,NROW1,NSLICE1,
     &             MAXIM,'FIRST INPUT SURFACE',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        IF (NSAM .NE. NSAM1 .OR. NROW .NE. NROW1 .OR. 
     &      NSLICE .NE. NSLICE1) THEN
            CALL ERRT(101,'VOLUMES MUST HAVE SAME DIMENSIONS',NDUM)
            GOTO 9999
        ENDIF

C       OPEN SECOND SPIDER VOLUME AS INPUT
        MAXIM    = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNIMV2,'O',IFORM,
     &              NSAM1,NROW1,NSLICE1,
     &             MAXIM,'SECOND INPUT VOLUME',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        IF (NSAM .NE. NSAM1 .OR. NROW .NE. NROW1 .OR. 
     &      NSLICE .NE. NSLICE1) THEN
            CALL ERRT(101,'VOLUMES MUST HAVE SAME DIMENSIONS',NDUM)
            GOTO 9999
        ENDIF

C       OPEN SECOND SPIDER SURFACE VOLUME AS INPUT
        MAXIM    = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNIMS2,'O',IFORM,
     &             NSAM1,NROW1,NSLICE1,
     &             MAXIM,'SECOND INPUT SURFACE',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        IF (NSAM .NE. NSAM1 .OR. NROW .NE. NROW1 .OR. 
     &      NSLICE .NE. NSLICE1) THEN
            CALL ERRT(101,'VOLUMES MUST HAVE SAME DIMENSIONS',NDUM)
            GOTO 9999
        ENDIF

        ALLOCATE (BUFV1(NSAM),BUFV2(NSAM),
     &            BUFS1(NSAM),BUFS2(NSAM),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(101,'ALLOCATION FAILED',NDUM)
            GOTO 9999
        ENDIF

        ALLOCATE (ISURF1(3,NSAM * NROW * NSLICE), 
     &            ISURF2(3,NSAM * NROW * NSLICE),
     &            RSURF1(NSAM * NROW * NSLICE), 
     &            RSURF2(NSAM * NROW * NSLICE), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(101,'ALLOCATION FAILED',NDUM)
            GOTO 9999
        ENDIF

C       INCASE SURFLEVEL IS NOT EXACT DIGIT
        THLEV  = SURFLEVEL - FLTZER

        NSURF1 = 0
        NSURF2 = 0
        DO ISLICE = 1,NSLICE
          DO IROW = 1, NROW
             IREC = (ISLICE -1) * NROW + IROW

             CALL REDLIN(LUNIMV1,BUFV1,NSAM,IREC)
             CALL REDLIN(LUNIMS1,BUFS1,NSAM,IREC)
             CALL REDLIN(LUNIMV2,BUFV2,NSAM,IREC)
             CALL REDLIN(LUNIMS2,BUFS2,NSAM,IREC)

             DO ISAM = 1,NSAM
                IF (BUFS1(ISAM) .GE. THLEV) THEN
C                  SURFACE VOXEL
                   NSURF1           = NSURF1 + 1
                   ISURF1(1,NSURF1) = ISAM
                   ISURF1(2,NSURF1) = IROW
                   ISURF1(3,NSURF1) = ISLICE
                   RSURF1(NSURF1)   = BUFV2(ISAM)
                ENDIF

                IF (BUFS2(ISAM) .GE. THLEV) THEN
C                  SURFACE VOXEL
                   NSURF2           = NSURF2 + 1
                   ISURF2(1,NSURF2) = ISAM
                   ISURF2(2,NSURF2) = IROW
                   ISURF2(3,NSURF2) = ISLICE
                   RSURF2(NSURF2)   = BUFV1(ISAM)
                ENDIF

             ENDDO
           ENDDO
        ENDDO
        CLOSE(LUNIMV1)
        CLOSE(LUNIMV2)
        CLOSE(LUNIMS1)
        CLOSE(LUNIMS2)

C       FIND INTER-SURFACE DISTANCES

        NUMDIST    = 0
        DISTMINAV  = 0.0
        DISTMINAV2 = 0.0
        IDIF       = 0.0

        DO N1 = 1,NSURF1
           IX1 = ISURF1(1,N1)
           IY1 = ISURF1(2,N1)
           IZ1 = ISURF1(3,N1)

           IDISTMIN = MAX(NSAM,NROW,NSLICE) ** 2

           DO N2 = 1, NSURF2
              IX2 = ISURF2(1,N2)
              IY2 = ISURF2(2,N2)
              IZ2 = ISURF2(3,N2)

              IDIST = (IX2 -IX1)**2 + (IY2 - IY1)**2 + (IZ2 - IZ1)**2
              IF (IDIST .LT. IDISTMIN) THEN
                 IDISTMIN = IDIST
                 NODE2    = N2
              ENDIF
           ENDDO
           
C          FOUND MINIMAL DISTANCE FROM SURFACE 1 TO VOXEL ON SURFACE 2
           DIF        = RSURF2(NODE2) - RSURF1(N1)
           IF (DIF .NE. 0) IDIF = IDIF + SIGN(1.0,DIF)

           NUMDIST    = NUMDIST + 1
           DISTMINAV  = DISTMINAV  + SQRT(REAL(IDISTMIN))
           DISTMINAV2 = DISTMINAV2 + IABS(IDISTMIN)
                    
        ENDDO

        WRITE(NOUT,*) ' '

        WRITE(NOUT,90) NSURF1
90      FORMAT(' Total Voxels on First Surface:     ',I10)
        WRITE(NOUT,91) NSURF2
91      FORMAT(' Total Voxels on 2nd Surface:       ',I10)

        WRITE(NOUT,92) NUMDIST
92      FORMAT(' Total Number of Distances:         ',I10)

        AVDIST2  = DISTMINAV / NUMDIST
        WRITE(NOUT,93) AVDIST2
93      FORMAT(' Mean Inter-surface Distance:       ',G14.3)

        DSIG   = SQRT((DISTMINAV2 - DISTMINAV**2/NUMDIST) / 
     &                 DBLE(NUMDIST-1.0))

        WRITE(NOUT,94) DSIG
94      FORMAT(' S.D. of Distances:                 ',G14.3)

        WRITE(NOUT,95) IDIF
95      FORMAT(' Summed Sign of Density Difference: ',I10)

9999    CONTINUE
        IF (ALLOCATED(BUFV1))  DEALLOCATE(BUFV1)
        IF (ALLOCATED(BUFV2))  DEALLOCATE(BUFV2)
        IF (ALLOCATED(BUFS1))  DEALLOCATE(BUFS1)
        IF (ALLOCATED(BUFS2))  DEALLOCATE(BUFS2)
        IF (ALLOCATED(ISURF1)) DEALLOCATE(ISURF1)
        IF (ALLOCATED(ISURF2)) DEALLOCATE(ISURF2)
        IF (ALLOCATED(RSURF1)) DEALLOCATE(RSURF1)
        IF (ALLOCATED(RSURF2)) DEALLOCATE(RSURF2)

        RETURN
        END
    

