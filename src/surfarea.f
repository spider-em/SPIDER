
C **********************************************************************
C
C   SURFAREA.F   -- CREATED SEPT. 98 al
C
C **********************************************************************
C *  AUTHOR: ArDean Leith 
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
C      SURFAREA(MAXDIM)
C
C      PURPOSE:     READS SPIDER  PICTURE FILE, FINDS INTERFACIAL AREA
C
C      PARAMETERS:  MAXDIM      UNLABELED COMMON BLOCK LENGTH
C
C      CALLED BY:   
C
C      CALLS:  
C
C      NOTES:       NOT OPTIMIZED FOR SPEED I JUST WANTED TO GET IT 
C                   RUNNING QUICK, DO NOT EXPECT MUCH USAGE     
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE SURFAREA(MAXDIM)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        COMMON /IOBUF/ BUF(NBUFSIZ)

        INTEGER, ALLOCATABLE, DIMENSION(:) :: IMGS
        CHARACTER(LEN=MAXNAM)              :: FILNAM
 
        DATA    FLTZER/10E-30/

        LUNIM    = 11
        
C       OPEN SPIDER FILE AS INPUT
20      MAXIM    = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNIM,'O',IFORM,NSAM,NROW,NSLICE,
     &             MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        NSLICE1 = NSLICE
        NPIXP1  = NSAM * NROW + 1

        IF (NSAM .GT. NBUFSIZ) THEN
            CALL ERRT(6,'SURFAREA',NE)
            GOTO 999
        ENDIF

        IF (IMAMI .NE. 1) THEN
            CALL NORM3(LUNIM,NSAM,NROW,NSLICE,FMAX,FMIN,AV)
        ENDIF

        IF ((FMAX - FMIN) .LT. FLTZER) THEN
            CALL ERRT(101,'BLANK FILE SKIPPED',NE)
            GOTO 999
        ENDIF

C       DISPLAY MAX AND MIN VALUE OF PICTURE , ASK FOR THE LEVEL
        WRITE(NOUT,91) FMIN,FMAX
91      FORMAT('  IMAGE RANGE: ',1PG11.3,'....',1PG11.3)

C       FIND THRESHOLD LEVEL FOR CLUSTERS            
22      CALL RDPRM2S(THLEV,FDUM,NOT_USED,'THRESHOLD LEVEL',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 20

        WRITE(NOUT,*) ' '

        MEMWANT = NSAM * NROW * 2
        ALLOCATE (IMGS(MEMWANT), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'SURFAREA, IMGS',MEMWANT)
           GO TO 999
        ENDIF

        INOW  = 1
        INEXT = NPIXP1
        IVOL  = 0
        IVOX  = 0
        ISURF = 0

        DO ISLICE = 1,NSLICE

           NREC1  = (ISLICE - 1) * NROW + 1
           NREC2  = NREC1 + NROW - 1

C          LOAD FIRST SLICE 
           CALL SURFAREA2(LUNIM,BUF,NSAM,NREC1,NREC2,
     &                THLEV,IMGS(INOW))

           IF (ISLICE .LT. NSLICE) THEN
C             MUST LOAD FOLLOWING SLICE ALSO
              NREC1N  = ISLICE * NROW + 1
              NREC2N  = NREC1N + NROW - 1

              CALL SURFAREA2(LUNIM,BUF,NSAM,NREC1N,NREC2N,
     &           THLEV,IMGS(INEXT))
           ENDIF

C          PROCESS CURRENT SLICE FOR AREA 
           CALL SURFAREA1(NSAM,NROW,NSLICE,ISLICE,
     &                IMGS(INOW),IMGS(INEXT),IVOLT,IVOXT,ISURFT,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (NSLICE .GT. 1) WRITE(NOUT,96) ISLICE,ISURFT
96         FORMAT(' SLICE:',I4,',  INTERFACE AREA=',I8)

           IVOL  = IVOL  + IVOLT
           IVOX  = IVOX  + IVOXT
           ISURF = ISURF + ISURFT
           
        ENDDO

        WRITE(NOUT,92) IVOX
92      FORMAT(/' TOTAL VOXELS:           ',I10)

        WRITE(NOUT,93) IVOL
93      FORMAT(' VOXELS ABOVE THRESHOLD: ',I10)

        WRITE(NOUT,95) ISURF
95      FORMAT(' INTERFACE AREA:         ',I10/)

C       IF (NSEL(1).NE. 0.) PARAM(NSEL(1)) = IVOX
C       IF (NSEL(2).NE. 0.) PARAM(NSEL(2)) = IVOL
C       IF (NSEL(3).NE. 0.) PARAM(NSEL(3)) = ISURF

        CALL REG_SET_NSEL(1,3,FLOAT(IVOX),FLOAT(IVOL),
     &                      FLOAT(ISURF), 0.0,0.0, IRTFLG)
       
999     CONTINUE
C       CLOSE THE FILE
        CLOSE(LUNIM)

C       DEALLOCATE RUN-TIME MEMORY
        IF (ALLOCATED(IMGS)) DEALLOCATE(IMGS)

        RETURN
        END
    
C **********************************************************************
C
C      SURFAREA1SURFAREA1(NSAM,NROW,NSLICE,ISLICE,
C                         MSLICE1,MSLICE2,ISURF,IRTFLG)
C
C      PURPOSE:  DETERMINES SURFACE AREA USING 2 SLICES AT A TIME    
C
C      PARAMETERS:  
C
C      CALLED BY:  SURFAREA1  
C
C      CALLS:       
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--********************************************************************

       SUBROUTINE SURFAREA1(NSAM,NROW,NSLICE,ISLICE,
     &                      MSLICE1,MSLICE2,IVOL,IVOX,ISURF,IRTFLG)

       INTEGER      MSLICE1(*),MSLICE2(*)

       ISURF  = 0
       IVOL   = 0
       IVOX   = 0

       DO  IROW = 1, NROW 
          IPTR0 = (IROW-1) * NSAM

          DO  ICOL = 1, NSAM 
             IPTR1 = IPTR0 + ICOL
             NOW   = MSLICE1(IPTR1)
            
             IVOX = IVOX + 1
             IF (NOW .LT. 0) THEN
C               VOXEL OCCUPIED

                IVOL = IVOL + 1

C               COUNT SURFACE  UP & DOWN (IF BORDER) 
                IF (IROW .EQ. 1 .OR. IROW .EQ. NROW) ISURF = ISURF + 1

C               COUNT SURFACE TO LEFT & RIGHT (IF BORDER) 
                IF (ICOL .EQ. 1 .OR. ICOL .EQ. NSAM) ISURF = ISURF + 1

C               CHECK VOXEL TO LEFT FROM THIS VOXEL FOR SURFACE
                IF (ICOL .GT. 1 .AND.
     &              MSLICE1(IPTR1 - 1) .EQ. 0)  ISURF = ISURF + 1

C               CHECK VOXEL TO RIGHT FROM THIS VOXEL FOR SURFACE
                IF (ICOL .LT. NSAM .AND.
     &              MSLICE1(IPTR1 + 1) .EQ. 0)  ISURF = ISURF + 1

C               CHECK VOXEL UP FROM THIS VOXEL FOR SURFACE
                IF (IROW .GT. 1 .AND.
     &              MSLICE1(IPTR1 - NSAM) .EQ. 0) ISURF = ISURF + 1

C               CHECK VOXEL DOWN FROM THIS VOXEL FOR SURFACE
                IF (IROW .LT. NROW .AND.
     &              MSLICE1(IPTR1 + NSAM) .EQ. 0) ISURF = ISURF + 1

C               COUNT SURFACE  ABOVE (IF BORDER) 
                IF (NSLICE .GT. 1 .AND. 
     &             (ISLICE .EQ. 1 .OR. ISLICE .EQ. NSLICE)) 
     &              ISURF = ISURF + 1

C               CHECK VOXEL BELOW (SLICE) FROM THIS VOXEL FOR SURFACE
                IF (ISLICE .LT. NSLICE .AND.
     &              MSLICE2(IPTR1) .EQ. 0) ISURF = ISURF + 1

             ELSE
C               CHECK NEXT VOXEL DOWN FROM THIS VOXEL FOR SURFACE
                IF (ISLICE .LT. NSLICE .AND.
     &              MSLICE2(IPTR1) .LT. 0)  ISURF = ISURF + 1
                
             ENDIF
          END DO
       END DO

       RETURN
       END

C **********************************************************************
C
C   SURFAREA2.FOR  -- CREATED OCT 98
C **********************************************************************
C *  AUTHOR: ArDean Leith 
C **********************************************************************
C
C      SURFAREA2(LUNIM,BUF,NSAM,NREC1,NREC2,THLEV,SLICE)  
C
C      PURPOSE:     READS SPIDER PICTURE FILES SLICES INTO SLICE ARRAY 
C
C      PARAMETERS:  
C
C      CALLED BY:   
C
C      CALLS:       REDLIN 
C
C--********************************************************************

         SUBROUTINE SURFAREA2(LUNIM,BUF,NSAM,NREC1,NREC2,
     &                     THLEV,SLICE)
 
         INTEGER        SLICE(*)
         DIMENSION      BUF(*)

C        READ THE SPIDER FILE INTO SLICE ARRAY

         IPTR = 0 

           DO  I = NREC1,NREC2
             CALL REDLIN(LUNIM,BUF,NSAM,I)

             DO  J = 1,NSAM
               IPTR = IPTR + 1

               IF (BUF(J) .GE. THLEV) THEN
C                 VOXEL IS ABOVE THRESHOLD
                  SLICE(IPTR) = -1
               ELSE
                  SLICE(IPTR) = 0
               ENDIF
             END DO
           END DO

         RETURN
         END



