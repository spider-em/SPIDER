C **********************************************************************
C
C ADDS.F                      USED FILELIST    JULY 1999  ArDean Leith
C                             USED F90 J       JULY 1999  ArDean Leith
C                             SQRT(NEG) TRAP   APR  2000  ArDean Leith
C                             ALLOC TRAP       APR  2003  ArDean Leith
C                             REDVOL           MAY  2017  ArDean Leith
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2017  Health Research Inc.,                         *
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
C  PURPOSE:   CREATE AVERAGE AND VARIANCE FILES FROM A SET OF INPUT
C             IMAGES/VOLUMES OR EVEN/ODD SUBSETS FROM THE FILES 
C
C             OBSOLETE OPERATION 'AS' GIVES DIFFERING RESULTS ON
C             MODERN XEONS DUE TO SUMMATION OF VALUES NEAR ZERO
C             IN DIFFERING PRECISIONS.  al may 2017
C--*********************************************************************

      SUBROUTINE ADDS(LUNA,LUNIN,LUNV,NDOC,IDUM)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      INTEGER                :: LUNA,LUNIN,LUNV,NDOC,IDUM
      CHARACTER (LEN=MAXNAM) :: FILNAM,FILA,FILV,FILPAT

      DOUBLE PRECISION       :: VARAV, AV2, VV
      CHARACTER (LEN=1)      :: SER
      CHARACTER (LEN=1)      :: NULL   = CHAR(0)

      INTEGER                :: NILMAX,MAXIM,NLETP,NX,NY,NZ,NDUM,NUMT
      INTEGER                :: IMGNUM,IRTFLG,NVOX3,NUMC,LUP,NLETA
      INTEGER                :: ILUP,NLETV,NSEL_USED,NE,NUMPR,NOT_USED
      REAL                   :: OFFOLD,FDUM,AV1,OFF,FNUMT,VAV,SAV,SAVT
      INTEGER                :: NUMA,NMT,IFIL,I

      INTEGER,PARAMETER      :: LUNXM  = 0  ! SELFILE NOT ALLOWED
      LOGICAL,PARAMETER      :: FOUROK = .TRUE.

      REAL, ALLOCATABLE      :: AVGARAY(:),VARARAY(:),BUFIN(:)

#ifndef SP_32
      INTEGER *8             :: NVOX,IOK8
#else
      INTEGER *4             :: NVOX,IOK8
#endif
         
      NILMAX = NIMAX

C     OPEN INPUT IMAGE(S)
      MAXIM = 0
      CALL OPFILES(0,LUNIN,NDOC,LUNXM, 
     &             .TRUE.,FILPAT,NLETP, 'O',
     &             IFORM,NX,NY,NZ,MAXIM,
     &             'INPUT FILE TEMPLATE (E.G. PIC****)~',
     &             FOUROK,INUMBR,NILMAX, 
     &             NDUM,NUMT,IMGNUM, IRTFLG) 
      IF (IRTFLG .NE. 0) RETURN

      CLOSE(LUNIN)

      NVOX = NX * NY * NZ
C     COMPLAIN IF EXCESSIVE ALLOCATION
      NVOX3 = NVOX * 3
      CALL BIGALLOC(NVOX,IOK8,.FALSE.,.TRUE.,IRTFLG)

      ALLOCATE(AVGARAY(NVOX),
     &         VARARAY(NVOX),
     &         BUFIN(NVOX),STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
          CALL ERRT(46,'ADDS; AVGARAY..',NVOX3)
          GOTO 999
      ENDIF

      CALL  RDPRMC(SER,NUMC,.TRUE., 
     &             'ALL, or ODD-EVEN FILES (A/O)',NULL,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 999

      IF (SER.NE.'O' .AND. SER.NE.'E') THEN
         SER = 'A'
         LUP = 1
      ELSE
         LUP = 2
      ENDIF

      IF (VERBOSE) WRITE(NOUT,*) ' '
      DO ILUP=1,LUP
         IF (LUP == 2)  THEN
            IF (ILUP == 1)  THEN
               WRITE(NOUT,*) ' FOR ODD-NUMBERED IMAGES'
               SER = 'O'
               CALL FILERD(FILA,NLETA,NULL,
     &             'AVERAGE  FILE FOR ODD-NUMBERED  IMAGES~',IRTFLG)
               IF (IRTFLG .NE. 0)  GOTO 997
               CALL FILERD(FILV,NLETV,NULL,
     &             'VARIANCE FILE FOR ODD-NUMBERED  IMAGES~',IRTFLG)
               IF (IRTFLG .NE. 0)  GOTO 997
            ELSE
               WRITE(NOUT,*) ' FOR EVEN-NUMBERED IMAGES'
               SER = 'E'
               CALL FILERD(FILA,NLETA,NULL,
     &             'AVERAGE  FILE FOR EVEN-NUMBERED IMAGES~',IRTFLG)
               IF (IRTFLG .NE. 0)  GOTO 997
               CALL FILERD(FILV,NLETA,NULL,
     &             'VARIANCE FILE FOR EVEN-NUMBERED IMAGES~',IRTFLG)
               IF (IRTFLG .NE. 0)  GOTO 997
            ENDIF

         ELSEIF (FCHAR(4:5) .NE. 'AD') THEN
            CALL FILERD(FILA,NLETA,NULL,'AVERAGE ',IRTFLG)
            IF (IRTFLG .NE. 0)  GOTO 997
            CALL FILERD(FILV,NLETA,NULL,'VARIANCE',IRTFLG)
            IF (IRTFLG .NE. 0)  GOTO 997
         ENDIF

C        COMMAND  'AS AD'  MEANS THAT AVERAGE FILE ALREADY EXISTS
         IF (FCHAR(4:5) == 'AD')  THEN

C           FILES ALREADY EXIST,FIND OUT HOW MANY PICTURES ALREADY AVERAGED

5           CALL REG_GET_USED(NSEL_USED)
            IF (NSEL_USED <= 0)  THEN
               CALL ERRT(101,
     &                  'MISSING OFFSET VALUE IN INPUT REGISTER',NE)
               GOTO 997
            ENDIF

C           OPEN EXISTING AVERAGE FILE
            CALL OPFILEC(0,.TRUE.,FILA,LUNA,'Z',
     &                   IFORM,NX,NY,NZ,
     &                   MAXIM,'EXISTING AVERAGE',.TRUE.,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 997

            NUMPR = 1
            CALL  RDPRI1S(NUMPR,NOT_USED,
     &             'NO. OF IMAGES ALREADY AVERAGED',IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 997

C           OPEN EXISTING VARIANCE FILE
            CALL OPFILEC(0,.TRUE.,FILV,LUNV,'Z',IFORM,NX,NY,NZ,
     &                   MAXIM,'EXISTING VARIANCE',.TRUE.,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 997

C           INPUT EXISTING AVERAGE
            CALL REDVOL(LUNA,NX,NY, 1,NZ, AVGARAY,IRTFLG)

C           INPUT EXISTING VARIANCE
            CALL REDVOL(LUNV,NX,NY, 1,NZ, VARARAY,IRTFLG)

C           INITIALIZE BOTH AVERAGE AND VARIANCE ARRAYS
            VARARAY = VARARAY * (NUMPR-1) + AVGARAY * AVGARAY * NUMPR

            AVGARAY = AVGARAY * NUMPR

C           GET OLD OFFSET FROM REGISTER LINE
            CALL REG_GET_NSEL(1,OFFOLD,FDUM,FDUM,FDUM,FDUM,IRTFLG)
         ELSE

C           FILE DOES NOT EXIST - OPEN WITH DIMS OF FIRST FILE & INITIALIZE

C           OPEN NEW AVERAGE OUTPUT FILE
            MAXIM = 0
            CALL OPFILEC(0,.FALSE.,FILA,LUNA,'U',IFORM,NX,NY,NZ,
     &                      MAXIM,' ',.TRUE.,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 997

C           OPEN NEW VARIANCE OUTPUT FILE
            MAXIM = 0
            CALL OPFILEC(0,.FALSE.,FILV,LUNV,'U',IFORM,NX,NY,NZ,
     &                   MAXIM,' ',.TRUE.,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 997

C           CLEAR BOTH AVERAGE AND VARIANCE ARRAYS
            AVGARAY = 0.0
            VARARAY = 0.0

            NUMPR   = 0   ! NUMBER OF PREVIOUS PICTURES IS 0
         ENDIF

C        CHANGED TO ALLOW GAPS IN FILE SERIES 7/21/89 MR.
         AV1  = 0.0
         AV2  = 0.0
         NUMA = 0
         NMT  = 0

         DO IFIL=1,NUMT

            CALL FILGET(FILPAT,FILNAM,NLETP,INUMBR(IFIL),IRTFLG)
            IF (IRTFLG .NE. 0) THEN
               CALL ERRT(3,'ADDS',NE)
               GOTO 997
            ENDIF
           
            CALL OPFILEC(0,.FALSE.,FILNAM,LUNIN,'Z',IFORM,
     &                   NX,NY,NZ,
     &                   MAXIM,' ',.TRUE.,IRTFLG)
            IF (IRTFLG .NE. 0) THEN
C              TO ALLOW GAPS IN FILE SERIES 7/21/89 MR.:
               WRITE(NOUT,100) FILNAM(1:NLETP)
100            FORMAT(' FILE: ',A,' NOT FOUND, FILE SKIPPED')
               CYCLE
            ENDIF
            MAXIM = 0

            NMT = NMT + 1
            IF ( SER == 'O' .AND. MOD(NMT,2) == 0 )  GOTO  707
            IF ( SER == 'E' .AND. MOD(NMT,2) == 1 )  GOTO  707

            NUMA = NUMA + 1
 
            CALL  REDVOL(LUNIN,NX,NY,1,NZ,BUFIN,IRTFLG)

	    AV2 = SUM(BUFIN) / NVOX

            IF (FCHAR(4:4) .NE. 'R')  THEN
	       BUFIN = BUFIN - AV2
	    ENDIF

            AVGARAY = AVGARAY + BUFIN
            VARARAY = VARARAY + BUFIN * BUFIN
            AV1     = AV1 + AV2   ! OVERALL SUM (FOR AVERAGE)

707         CLOSE(LUNIN)
         ENDDO

C        NOW COMPUTE NORMALIZED AVERAGE AND VARIANCE IMAGE

         FNUMT   = NUMA + NUMPR
         AVGARAY = AVGARAY / FNUMT  ! DIVIDE BY # OF IMAGES

         VV      = DOT_PRODUCT(AVGARAY,AVGARAY)
         VV      = VV / REAL(NVOX)

         VARARAY = (VARARAY - AVGARAY*AVGARAY*FNUMT) /(FNUMT-1.0)
         VARAV   = SUM(VARARAY)
         VAV     = VARAV /FLOAT(NVOX)

         OFF     = AV1 / FLOAT(NUMA)

         CALL WRTVOL(LUNA,NX,NY, 1,NZ, AVGARAY,IRTFLG)
         
         CALL WRTVOL(LUNV,NX,NY, 1,NZ, VARARAY,IRTFLG)

         IF (VARAV < 0.0) THEN
            WRITE(NOUT,*)   ' WARNING: NEGATIVE TOTAL VARIANCE: ',VARAV,
     &                      '  SD SET TO ZERO!!!'
            IF (NOUT .NE. NDAT)
     &        WRITE(NDAT,*) ' WARNING: NEGATIVE TOTAL VARIANCE: ',VARAV,
     &                      '  SD SET TO ZERO!!!'
            SAVT = 0.0
         ELSE
            SAVT = SQRT(VARAV)
         ENDIF

         IF (VAV < 0.0) THEN
            WRITE(NOUT,*)   ' WARNING: NEGATIVE  VARIANCE P.P. : ',VAV,
     &                      '  SD SET TO ZERO!!!'
            IF (NOUT .NE. NDAT)
     &        WRITE(NDAT,*) ' WARNING: NEGATIVE TOTAL VARIANCE: ',VAV,
     &                      '  SD SET TO ZERO!!!'
            SAV = 0.0
         ELSE
            SAV  = SQRT(VAV)
         ENDIF

         CALL REG_GET_USED(NSEL_USED)
         IF (NUMPR .NE. 0 .AND. NSEL_USED .NE. 0) THEN
C           TRAP FOR NaN
            IF (OFFOLD .NE. 0.0)
     &         OFF = (OFFOLD*FLOAT(NUMPR)+AV1)/FLOAT(NUMA+NUMPR)
         ENDIF

         CALL REG_SET_NSEL(1,1, OFF,0.0, 0.0, 0.0, 0.0,IRTFLG)

         IF (NOUT .NE. NDAT)
     &      WRITE(NDAT,7001) NUMA+NUMPR,NVOX,VARAV,SAVT,
     &                       VAV,SAV,OFF,VV

         WRITE(NOUT,7001) NUMA+NUMPR, NVOX, OFF,VARAV,SAVT,
     &                    VAV,SAV,OFF,VV
7001     FORMAT(/,'  VARIANCE COMPUTATION FOR:',I7,' IMAGES ',
     &          '    CONTAINING: ',I0,' ELEMENTS',/,
     &      '   AVERAGE:        ',1PG12.5,/,
     &      '   TOTAL VARIANCE: ',1PG12.5,'  TOTAL S.D.: ',1PG12.5,/,
     &      '   P.P.  VARIANCE: ',1PG12.5,'  P.P.  S.D.: ',1PG12.5,/,
     &      '   AVERAGE OFFSET: ',1PG12.5,/,
     &      '   VARIANCE OF AVERAGE IMAGE: ',1PG12.5,/)

	ENDDO

997     CLOSE(LUNV)
        CLOSE(LUNA)
        CLOSE(LUNIN)

999     IF (ALLOCATED(BUFIN))   DEALLOCATE(BUFIN)
        IF (ALLOCATED(VARARAY)) DEALLOCATE(VARARAY)
        IF (ALLOCATED(AVGARAY)) DEALLOCATE(AVGARAY)

        END
