
C **********************************************************************
C
C   SURFFIT.F   -- CREATED                         MAR 00 ARDEAN LEITH
C                  MAXNAM                          JUL 14 ARDEAN LEITH
C
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
C SURFFIT()
C
C PURPOSE:     READS SPIDER 3-D PICTURE FILE, CREATES 3-D IMAGE 
C              FILE CONTAINING NUMBERS FOR CONNECTED CLUSTERS
C
C PARAMETERS:  
C
C CALLS:       
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE SURFFIT()

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL, ALLOCATABLE, DIMENSION(:) :: SLICES
        REAL, ALLOCATABLE, DIMENSION(:) :: SURFACE
        REAL, ALLOCATABLE, DIMENSION(:) :: BUF

        REAL, PARAMETER                 :: BACKLEVEL = 0.0
        REAL, PARAMETER                 :: SURFLEVEL = 1.0

        CHARACTER(LEN=MAXNAM)           :: FILNAM
        LOGICAL                         :: COMPARE

        LUNIM    = 11
        LUNOUT   = 12

C       OPEN SPIDER VOLUME AS INPUT
20      MAXIM    = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNIM,'O',IFORM,NSAM,NROW,NSLICE,
     &             MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (IMAMI.NE.1) CALL NORM3(LUNIM,NSAM,NROW,NSLICE,FMAX,FMIN,AV)
 
C       DISPLAY MAX AND MIN VALUE OF PICTURE, ASK FOR THE LEVEL
        WRITE(NOUT,91) FMIN,FMAX
91      FORMAT(' IMAGE RANGE: ',1PG11.3,'....',1PG11.3)

C       FIND THRESHOLD LEVEL FOR SURFACE            
22      CALL RDPRM1S(THLEV,NOT_USED,'THRESHOLD LEVEL',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 20

C       OPEN NEW SPIDER VOLUME AS OUTPUT
23      MAXIM = 0
        CALL OPFILEC(LUNIM,.TRUE.,FILNAM,LUNOUT,'N',IFORM,
     &             NSAM,NROW,NSLICE,
     &             MAXIM,'SURFACE OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 22
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (ALLOCATED(SLICES)) DEALLOCATE(SLICES)
        ALLOCATE (SLICES(NSAM * NROW * 2), SURFACE(NSAM * NROW * 2), 
     &            STAT=IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (ALLOCATED(BUF)) DEALLOCATE(BUF)
        ALLOCATE (BUF(NSAM), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        NUMSURFS = 0
 
        DO ISLICE = 1,NSLICE

           IF (MOD(ISLICE,2) .NE. 0) THEN
C            CURRENT SLICE IS IN SLICE1
             INOW  = 1
             INEXT = NSAM*NROW+1
           ELSE
C            NEXT SLICE GOES INTO SLICE1 
             INOW  = NSAM*NROW+1
             INEXT = 1
           ENDIF 

           NREC1  = (ISLICE - 1) * NROW + 1
           NREC2  = NREC1 + NROW - 1

           IF (ISLICE .EQ. 1) THEN
C             MUST LOAD FIRST SLICE 
              CALL SURFFIT_RED(LUNIM,BUF,NSAM,NROW,NSLICE,
     &           1,1,THLEV,SURFLEVEL,BACKLEVEL,SLICES(1),IRTFLG)
           ENDIF

           IF (ISLICE .LT. NSLICE) THEN
C             LOAD NEW NEXT SLICE
              ISLICEN = ISLICE + 1
              CALL SURFFIT_RED(LUNIM,BUF,NSAM,NROW,NSLICE,ISLICEN,
     &          ISLICEN, THLEV,SURFLEVEL,BACKLEVEL,SLICES(INEXT),IRTFLG)
           ENDIF


C          PROCESS CURRENT SLICE FOR SURFACE 
           CALL SURFFIT_SET(NSAM,NROW,NSLICE,
     &                  SLICES(INOW), SLICES(INEXT),
     &                  ISLICE,SURFACE(INOW),SURFACE(INEXT),
     &                  NUMSURFST,SURFLEVEL,BACKLEVEL,IRTFLG)


           IF (IRTFLG .NE. 0) RETURN
           NUMSURFS = NUMSURFS + NUMSURFST

C          STORE CURRENT SLICE IN OUTPUT FILE
           IPTR = 1
           DO  IREC = (ISLICE-1) * NROW + 1,(ISLICE-1) * NROW + NROW 
              CALL WRTLIN(LUNOUT,SURFACE(INOW + IPTR),NSAM,IREC)
              IPTR = IPTR + NSAM
           END DO

cc           WRITE(NOUT,96) ISLICE,NUMSURFS
96         FORMAT(' After slice:',I4,',  Surface voxels=',I10)

        ENDDO

        WRITE(NOUT,97) NUMSURFS
97      FORMAT(' Total Surface voxels=',I10)

9999    CONTINUE
C       CLOSE THE FILES
        CLOSE(LUNOUT)
        CLOSE(LUNIM)
        IF (ALLOCATED(BUF))     DEALLOCATE(BUF)
        IF (ALLOCATED(SURFACE)) DEALLOCATE(SURFACE)
        IF (ALLOCATED(SLICES))  DEALLOCATE(SLICES)

        RETURN
        END
    

C      ----------------------- SURFFIT_SET ------------------------

       SUBROUTINE SURFFIT_SET(NSAM,NROW,NSLICE,SLICE1,SLICE2,
     &                        ISLICE,SURFACE1,SURFACE2,NSURF,
     &                        SURFLEVEL,BACKLEVEL,IRTFLG)

 
       COMMON /UNITS/LUN,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT

       REAL, DIMENSION(NSAM*NROW) :: SLICE1,SLICE2
       REAL, DIMENSION(NSAM*NROW) :: SURFACE1,SURFACE2
       REAL, DIMENSION(NSAM) ::      BUF
       LOGICAL ::                    ATSURF

       DATA FLTZER/10E-10/
       
       THLEV  = BACKLEVEL

       IRTFLG = 0

       NMT    = 0
       NOCC   = 0
       NSURF  = 0

C      ZERO SURFACE FOR NEXT SLICE DOWN
       DO I = 1,NSAM * NROW 
          SURFACE2(I) = 0.0
       ENDDO

       DO  IROW = 1, NROW
          IPTR0 = (IROW-1) * NSAM

       DO  ICOL = 1, NSAM
          IPTR1 = IPTR0 + ICOL

          IF (SLICE1(IPTR1) .LE. THLEV) THEN
C            EMPTY VOXEL ON THIS SLICE
             ATSURF = .FALSE.
             NMT    = NMT + 1
       
             IF (ISLICE .LT. NSLICE .AND. SLICE2(IPTR1) .GT. THLEV)THEN
C                WILL HAVE A VOXEL ON SLICE2 AT SURFACE
                 SURFACE2(IPTR1) = SURFLEVEL
                 NSURF           = NSURF + 1
             ENDIF
          ELSE
C            OCCUPIED VOXEL
             NOCC   = NOCC + 1

             ATSURF = .FALSE.

             IF (ICOL .GT. 1 .AND. SLICE1(IPTR1-1) .LE. THLEV) THEN
                ATSURF = .TRUE.

             ELSEIF (ICOL .LT. NSAM .AND. SLICE1(IPTR1+1) .LE. THLEV)
     &          THEN
                ATSURF = .TRUE.

             ELSEIF (IROW .GT. 1 .AND. SLICE1(IPTR1-NSAM) .LE. THLEV)
     &          THEN
                ATSURF = .TRUE.

             ELSEIF(IROW .LT. NSAM .AND. SLICE1(IPTR1+NSAM).LE. THLEV)
     &          THEN
                ATSURF = .TRUE.

             ELSEIF (ISLICE .LT. NSLICE .AND. SLICE2(IPTR1) .LE. THLEV)
     &          THEN
                ATSURF = .TRUE.
             ENDIF

             IF (ATSURF) THEN
C               VOXEL ON THIS SLICE1 IS AT SURFACE
             
                IF (SURFACE1(IPTR1) .NE. SLICE1(IPTR1)) THEN
                   SURFACE1(IPTR1) = SURFLEVEL
                   NSURF           = NSURF + 1
                ENDIF
             ENDIF
          ENDIF

         ENDDO
         ENDDo

         END


C      ----------------------- SURFFIT_RED ------------------------

         SUBROUTINE SURFFIT_RED(LUNIM,BUF,NSAM,NROW,NSLICE,
     &                    ISLICE1,ISLICE2,THLEV,SURFLEVEL,BACKLEVEL,
     &                     VOLUME,IRTFLG)
 
         REAL, DIMENSION(NSAM*NROW) :: VOLUME
         REAL, DIMENSION(NSAM) ::      BUF

C        READ THE SPIDER FILE INTO SLICE ARRAY
         NREC1 = (ISLICE1 - 1) * NROW + 1
         NREC2 = (ISLICE2 - 1) * NROW + NROW 

         IPTR = 0 

           DO  I = NREC1,NREC2
             CALL REDLIN(LUNIM,BUF,NSAM,I)

             DO  J = 1,NSAM
               IPTR = IPTR + 1

               IF (BUF(J) .GT. THLEV) THEN
C                 VOXEL IS ABOVE THRESHOLD
                  VOLUME(IPTR) = SURFLEVEL
               ELSE
                  VOLUME(IPTR) = BACKLEVEL
               ENDIF
             END DO
           END DO


         RETURN
         END
    
