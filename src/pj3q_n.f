C++*********************************************************************
C                                                                      *
C PJ3Q_N.F      SPEEDED UP                       FEB 2000 ARDEAN LEITH *
C               LUNDOCREDANG PARAMETERS CHANGED  DEC 2000 ARDEAN LEITH *
C               OPENDOC PARAMETERS               DEC 2000 ARDEAN LEITH *
C               REWRITTEN                        SEP 2003 PAWEL        *
C               REFRINGS OPTION                  FEB 2005 ARDEAN LEITH *
C               FBS_WANTED                       JAN 2012 ARDEAN LEITH *                                                      *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C PJ3Q_N(FBS_WANTED,REFRINGS)
C
C PURPOSE:  COMPUTES PROJECTIONS OF A 3D VOLUME ACCORDING TO 
C           THREE EULERIAN ANGLES. DOES A WHOLE IMAGE SERIES. CAN
C           CREATE 'REFERENCE RINGS' FILE(S) ALSO 
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE PJ3Q_N(FBS_WANTED,REFRINGS) 
 
         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC'

         LOGICAL                  :: FBS_WANTED,REFRINGS
 
         CHARACTER(LEN=MAXNAM)    :: FINPAT,FINPIC
         CHARACTER(LEN=1)         :: MODE
         LOGICAL                  :: MD

         REAL, ALLOCATABLE        :: BCKE(:)
         INTEGER, ALLOCATABLE     :: IPCUBE(:,:)
         REAL, ALLOCATABLE        :: PROJ(:,:,:)
         REAL, ALLOCATABLE        :: ANGBUF(:,:)

         INTEGER                  :: ICOMM,MYPID,MPIERR,MAXIM,IRTFLG
         REAL                     :: R1

         INTEGER, PARAMETER       :: INVOL     = 98
         INTEGER, PARAMETER       :: INPRJ     = 97
         INTEGER, PARAMETER       :: INDOCAT   = 96
         INTEGER, PARAMETER       :: INDOCS    = 95 
         INTEGER, PARAMETER       :: LUNRINGST = 94 

         CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C        OPEN INPUT VOLUME
         MAXIM  = 0
         IRTFLG = 0
         CALL OPFILEC(0,.TRUE.,FINPAT,INVOL,'O',
     &               IFORM,NX,NY,NZ,
     &               MAXIM,'INPUT VOLUME~',.FALSE.,IRTFLG)
         IF (IRTFLG .NE. 0)  RETURN

         IF (NX .NE. NY .OR. NZ .NE. NX) THEN
            WRITE(NOUT,91) 
91          FORMAT(' WARNING NON-CUBIC VOLUME MAY GIVE INVALID RESULTS')
         ENDIF

         RI   = (MIN(NX,NY,NZ)/2) - 2
         CALL RDPRM1S(RI,NOT_USED,'OBJECT RADIUS',IRTFLG)
         IF (IRTFLG .NE. 0)  RETURN

C        FIND NUMBER OF OMP THREADS
         CALL GETTHREADS(NUMTH)

         NVOX = NX*NY*NZ

C        INITIALIZE NN  AND FIND NN 

         LDPX = NX  /2+1
         LDPY = NY  /2+1
         LDPZ = NZ/2+1

         NN   = 1
         MD   = .FALSE.
         CALL PREPCUB(NX,NY,NZ, NN,IDUM,
     &                RI,MD,LDPX,LDPY,LDPZ) 

C        USE NVOX TO ALLOCATE (BCKE) 
         MAXKEY = NIMAX
         ALLOCATE(BCKE(NVOX),
     &            IPCUBE(5,NN),
     &            ANGBUF(3,MAXKEY),STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN 
            MWANT = NVOX + 5 * NN + 3 * MAXKEY
            CALL ERRT(46,'PJ 3Q, BCKE...',MWANT)
            GOTO 9999   
         ENDIF

C        READ BCKE
	 CALL READV(INVOL,BCKE,NX,NY,NX,NY,NZ)
         CLOSE(INVOL)

C        PREPARE IPCUBE
         MD = .TRUE.
         CALL PREPCUB(NX,NY,NZ, NN,IPCUBE,
     &                RI,MD,LDPX,LDPY,LDPZ) 

         ALLOCATE(PROJ(NX,NY,NUMTH),STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN 
            MWANT = NX * NY  * NUMTH
            CALL ERRT(46,'PJ 3Q, PROJ',MWANT)
            GOTO 9999   
         ENDIF

C        READ SELECTION DOC FILE 
C        NANG - NUMBER OF ANGLES (PROJECTIONS)
         CALL FILELIST(.FALSE.,INDOCS,FINPAT,NLETA,INUMBR,MAXKEY,NANG,
     &        'ANGLE NUMBERS OR SELECTION DOC. FILE NAME',IRTFLG)
         IF (IRTFLG .NE. 0)  GOTO 9999         

C        OPEN ANGLES DOC FILE
         CALL OPENDOC(FINPIC,.TRUE.,NLETI,INDOCAT,INDOCA,.TRUE.,
     &     'ANGLES DOC',.TRUE.,.FALSE.,.FALSE.,NEWFILE,IRTFLG)
         IF (IRTFLG .NE. 0)  GOTO 9999

C        READ ANGLES FROM ANGLES DOC FILE.
C        ORDER IN THE DOCUMENT FILE IS PSI, THETA, PHI AND ANGLES ARE 
C        IN DEGREES! 

         CALL LUNDOCREDANG(INDOCA,ANGBUF,MAXKEY,NGOTY,MAXGOTY,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999
         IF (NGOTY .LT. NANG) THEN 
            CALL ERRT(102,'ONLY GOT ANGLES FOR IMAGES',NGOTY)
            GOTO 9999   
         ENDIF

         REFRINGS =  (FCHAR(7:7) .EQ. 'R') 
         LUNRINGS = 0
         IF (REFRINGS) THEN
C           WANT TO CREATE REFERENCE RINGS FILE
                     
            MR    = 5
            NR    = MIN(NX,NY)
            ISKIP = 0
            CALL RDPRI3S(MR,NR,ISKIP,NOT_USED,
     &                'FIRST, LAST RING, & RING SKIP',IRTFLG)
            IF (IRTFLG .NE. 0) GO TO 9999
            LUNRINGS = LUNRINGST

            MODE     = 'F'
         ENDIF

C        PROJECT NOW
         CALL WRITPRO_N(PROJ,INPRJ,NX,NY,NZ,NUMTH,BCKE,NVOX,
     &                 IPCUBE,NN,RI,INUMBR,NANG,MAXKEY,ANGBUF,
     &                 LUNRINGS,MODE,MR,NR,ISKIP,LDPX,LDPY,LDPZ,
     &                 FBS_WANTED,IRTFLG)

         IF (MYPID .LE. 0) THEN
            IF (FBS_WANTED) THEN
               WRITE(NOUT,92) NANG
92          FORMAT('  PJ 3F FINISHED FOR: ',I7,'   PROJECTIONS -----',/)
            ELSE
               WRITE(NOUT,93) NANG
93          FORMAT('  PJ 3Q FINISHED FOR: ',I7,'   PROJECTIONS -----',/)
            ENDIF

            CALL FLUSHRESULTS()
         ENDIF

9999     IF(ALLOCATED(ANGBUF)) DEALLOCATE(ANGBUF)
         IF(ALLOCATED(PROJ))   DEALLOCATE(PROJ)
         IF(ALLOCATED(IPCUBE)) DEALLOCATE(IPCUBE)
         IF(ALLOCATED(BCKE))   DEALLOCATE(BCKE)

         CLOSE(INDOCAT)
         CLOSE(INDOCS)
         IF (REFRINGS) CLOSE(LUNRINGS)

         RETURN
         END

