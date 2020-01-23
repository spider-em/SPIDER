
C **********************************************************************
C
C   CONINT.FOR  -- CREATED                       OCT 1990 ArDean Leith
C
C **********************************************************************
C=* AUTHOR: ArDean Leith                                               *
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2019  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C  CONINT(IRTFLG)
C
C  PURPOSE:     READS SPIDER 3-D PICTURE FILE, CREATES 3-D IMAGE 
C               FILE CONTAINING NUMBERS FOR CONNECTED CLUSTERS
C
C  PARAMETERS:  
C
C  CALLS:       CCONECT     FILSLI
C               MAPIM
C               EMPSLI      EMPSLI
C               MAKTAB      SHOSLI
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE CONINT(IRTFLG)

        INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC' 

C       NEQMAX IS MAXIMUM NUMBER OF PLACES IN BRANCHING EQUIV. TABLE
        PARAMETER  (NEQMAX = 16000)
		                                                                                        
        INTEGER                              :: IEQUIV(2,NEQMAX)

        INTEGER*2, ALLOCATABLE               :: SLICES(:)
        REAL, ALLOCATABLE                    :: BUF(:)
        INTEGER, ALLOCATABLE                 :: MSLICES(:)
        REAL, ALLOCATABLE                    :: TABLE(:)

C       REAL                                 :: VALUES(3)
        CHARACTER (LEN=1)                    :: NULL = CHAR(0)
        LOGICAL                              :: LASTSLI,DEBUGING
        CHARACTER (LEN=MAXNAM)               :: FILNAM
        INTEGER, PARAMETER                   :: LUNIM  = 11
        INTEGER, PARAMETER                   :: LUNOUT = 12
        REAL, PARAMETER                      :: FLTZER = 10E-30

        DEBUGING = .FALSE.

        NEQUIV   = 0
        LASTCLUS = 0

C       OPEN SPIDER FILE AS INPUT
20      MAXIM    = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNIM,'O',IFORM,NX,NY,NZ,
     &             MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        NZ1    = NZ
        NPIXP1 = NX * NY + 1

        NLEN    = NX * NY * 2        ! INTEGER*2
        ALLOCATE(SLICES(NLEN),MSLICES(NZ),BUF(NX), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           MWANT = NLEN/2 + NZ + NX
           CALL ERRT(46,'CONINT; SLICES...',MWANT)
           GOTO 9999
        ENDIF

        IF (IMAMI .NE. 1) THEN
             IF (NZ > 1) 
     &          WRITE(NOUT,*) ' NORMALIZING 3D FILE, PLEASE WAIT.'
             CALL NORM3(LUNIM,NX,NY,NZ1,FMAX,FMIN,AV)
        ENDIF

        IF ((FMAX - FMIN) < FLTZER) THEN
            WRITE(NOUT,*) ' *** ERROR: BLANK FILE SKIPPED '
            CALL ERRT(100,'CONINT',NE)
            GOTO 9999
        ENDIF

21      NUMSLI = NZ
        CALL RDPRAI(MSLICES,NZ,NUMSLI,1,NZ1,
     &       'SLICE NUMBERS',NULL,IRTFLG)
        IF (IRTFLG .EQ. -1) THEN
           CLOSE(LUNIM)
           GOTO 20
        ENDIF
        NZ2 = NUMSLI
 
C       DISPLAY MAX AND MIN VALUE OF PICTURE , ASK FOR THE LEVEL
        WRITE(NOUT,91) FMIN,FMAX
91      FORMAT(' IMAGE RANGE: ',1PG11.3,'....',1PG11.3)

C       FIND THRESHOLD LEVEL FOR CLUSTERS            
22      CALL RDPRM1S(THLEV,NOT_USED,'THRESHOLD LEVEL',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 21

23      MAXIM = 0
        IFORM = 3          
        CALL OPFILEC(LUNIN,.TRUE.,FILNAM,LUNOUT,'N',IFORM,
     &         NX,NY,NZ,MAXIM,'CLUSTER OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 22
        IF (IRTFLG .NE. 0) GOTO 9999

        LASTSLI = .FALSE.
 
        DO IPTR = 1,NUMSLI

           IF (MOD(IPTR,2) .NE. 0) THEN
C            CURRENT SLICE IS IN SLICE1
             INOW  = 1
             INEXT = NPIXP1
           ELSE
C            NEXT SLICE GOES INTO SLICE1 
             INOW  = NPIXP1
             INEXT = 1
           ENDIF 

           ISLICE = MSLICES(IPTR)
           NREC1  = (ISLICE - 1) * NY + 1
           NREC2  = NREC1 + NY - 1

           IF (IPTR .EQ. 1) THEN
C             MUST LOAD CURRENT SLICE ALSO
              CALL FILSLI(LUNIM,BUF,NX,NREC1,NREC2,.TRUE.,
     &                    THLEV,SLICES(1))
           ENDIF

           IF (IPTR < NUMSLI) THEN
C             LOAD NEW NEXT SLICE
              ISLICEN = MSLICES(IPTR+1)
              IF (ISLICEN .GT. NZ1) THEN
C                NO SUCH INPUT SLICE
                 WRITE(NOUT,*) ' *** SLICE NOT AVAILABLE:',ISLICEN
                 GOTO 40
              ENDIF
              NREC1N  = (ISLICEN - 1) * NY + 1
              NREC2N  = NREC1N + NY - 1
              CALL FILSLI(LUNIM,BUF,NX,NREC1N,NREC2N,.TRUE.,
     &          THLEV,SLICES(INEXT))
           ELSE
              LASTSLI = .TRUE.
           ENDIF

C          PROCESS CURRENT SLICE FOR CONNECTIVITY 
           CALL CCONECT(NX,NY,LUNOUT,SLICES(INOW),
     &          SLICES(INEXT),LASTSLI,IEQUIV,NEQUIV,
     &          NEQMAX,LASTCLUS,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          STORE CURRENT SLICE IN OUTPUT FILE
           NREC1N = (IPTR - 1) * NY + 1
           NREC2N = NREC1N + NY - 1
           CALL EMPSLI(LUNOUT,BUF,NX,NREC1N,NREC2N,SLICES(INOW))
           IF (DEBUGING) CALL SHOSLI(NOUT,BUF,NX,1,NX,SLICES(INOW))

           WRITE(NOUT,96) ISLICE,NEQUIV,LASTCLUS
96         FORMAT(' After slice:',I4,',  Branches=',I6,'  Clusters=',I6)

40         CONTINUE           
        ENDDO

C       ALL SLICES PROCESSED, START SECOND PASS THRU DATA FOR BRANCHES

C****************** DEBUGING
        IF (DEBUGING) THEN
        WRITE(10,*) ' IEQUIV '
        WRITE(10,793) ((IEQUIV(I,J),I=1,2),J=1,NEQUIV)
793     FORMAT(6(2I5,2X)) 
        ISLICE = 10          
        WRITE(10,*) ' AFTER FIRST PASS STACK SLICE:',ISLICE

        NREC1  = (ISLICE - 1) * NY + 1
        NREC2  = NREC1 + NY - 1
        CALL FILSLI(LUNOUT,BUF,NX,NREC1,NREC2,.FALSE.,0.0,SLICES)

        WRITE(10,9099)
        WRITE(10,*) ' window one  (200,90)...(240,200)'
        NREC1 = 90
        NREC2 = 200

        NX1 = 200
        NX2 = 340

        NX1 = 210
        NX2 = 239
        CALL SHOSLI2(10,BUF,NX,NX1,NX2,NREC1,NREC2,SLICES)
        NX1 = 240
        NX2 = 269

        WRITE(10,9099)
        CALL SHOSLI2(10,BUF,NX,NX1,NX2,NREC1,NREC2,SLICES)

        WRITE(10,9099)
9099    FORMAT('1')

C        WRITE(10,*) ' window two (415,190)...(475,325)'
C        NREC1 = 190
C        NREC2 = 325
C        NX1 = 415
C        NX2 = 444
C        CALL SHOSLI2(10,BUF,NX,NX1,NX2,NREC1,NREC2,SLICES)
C        WRITE(10,9099)
C        NX1 = 445
C        NX2 = 475
C        CALL SHOSLI2(10,BUF,NX,NX1,NX2,NREC1,NREC2,SLICES)
        ENDIF
C*******************************************
        ALLOCATE(TABLE(LASTCLUS),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           MWANT = LASTCLUS
           CALL ERRT(46,'EC CL: , TABLE,',MWANT)
           GOTO 9999
        ENDIF

        WRITE(NOUT,*) ' '
        WRITE(NOUT,*) ' Constructing mapping table, please wait.'
        CALL MAKTAB(IEQUIV,NEQUIV,TABLE,LASTCLUS,NTAB,IRTFLG)

c************DEBUGING
        IF (DEBUGING) THEN
           WRITE(10,*) ' TABLE '
           WRITE(10,6993) (IT,TABLE(IT),IT=1,LASTCLUS)
6993       FORMAT(7(I5,1X,F5.0))
        ENDIF
C***************************************

        NREC1 = 1
        NREC2 = NY * NZ2
        WRITE(NOUT,*) ' Remapping cluster numbers, please wait.'
        CALL MAPIM(LUNOUT,LUNOUT,NX,NREC1,NREC2,TABLE,LASTCLUS,
     &       BUF,IRTFLG)

ccc        values(1) = 1.0       
ccc        values(2) = nlab
ccc        values(3) = 0.0
ccc        negative irtflg supresses label change output
ccc        irtflg = -1         
ccc        call setlab(lunout,nx,buf,6,3,values,'u',irtflg)
C****************** DEBUGING

        IF (DEBUGING) THEN
          ISLICE = 10          
          WRITE(10,*) ' FINAL STACK SLICE:',ISLICE
          NREC1  = (ISLICE - 1) * NY + 1
          NREC2  = NREC1 + NY - 1
          CALL FILSLI(LUNOUT,BUF,NX,NREC1,NREC2,.FALSE.,0.0,SLICES)

          WRITE(10,9099)
          WRITE(10,*) ' window one  (200,90)...(240,200)'
          NREC1 = 90
          NREC2 = 200

          NX1 = 200
          NX2 = 340

          NX1 = 210
          NX2 = 239
          CALL SHOSLI2(10,BUF,NX,NX1,NX2,NREC1,NREC2,SLICES)
          NX1 = 240
          NX2 = 269

          WRITE(10,9099)
          CALL SHOSLI2(10,BUF,NX,NX1,NX2,NREC1,NREC2,SLICES)
C        WRITE(10,9099)
C        WRITE(10,*) ' window two (415,190)...(475,325)'
C        NREC1 = 190
C        NREC2 = 325
C        NX1 = 415
C        NX2 = 444
C        CALL SHOSLI2(10,BUF,NX,NX1,NX2,NREC1,NREC2,SLICES)
C        WRITE(10,9099)
C        NX1 = 445
C        NX2 = 475
C       CALL SHOSLI2(10,BUF,NX,NX1,NX2,NREC1,NREC2,SLICES)
         ENDIF
C***********************************************
       
9999    CONTINUE
C       CLOSE THE FILES
        CLOSE(LUNOUT)
        CLOSE(LUNIM)

        END
    


