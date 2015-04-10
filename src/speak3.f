C++*********************************************************************
C
C SPEAK3.F
C               RESTRICTED DISTANCE                 NOV 04 ARDEAN LEITH
C               DOC FILE OPEN FOR 'PK 3' BUG        DEC 04 ARDEAN LEITH
C               NMAX LIMIT FOR DOC FILE OUTPUT      AUG 05 ARDEAN LEITH
C               LUNDOC CLOSURE BUG                  OCT 10 ARDEAN LEITH
C               PROMPTS, * FOR DOC                  JUN 13 ARDEAN LEITH
C               LUNDOCRET BUG                       NOV 13 ARDEAN LEITH
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
C  SPEAK3(LUN,NX,NY,NZ,MAXDIM,OPT,LUNDOC)
C
C  PURPOSE:    SEARCHES FOR THE ML HIGHEST PEAKS IN A 
C              SPIDER VOLUME AND PRINTS OUT POSITIONS AND VALUES OF 
C              THESE PEAKS.
C
C  PARAMETERS:
C       LUN         LOGICAL UNIT NUMBER OF VOLUME
C       NX,NY,NZ    DIMENSIONS OF VOLUME
C       OPT         OUTPUT OPTION
C          OPT=' '  DEFAULT: NO DOCUMENT OUTPUT
C          OPT='D'  DOCUMENT OUTPUT: NUMBER,POSITION, AND VALUE
C                      OF PEAKS ARE WRITTEN INTO A DOCUMENT FILE
C          OPT='R'  DOCUMENT OUTPUT: NUMBER,POSITION, AND VALUE
C                      OF PEAKS WHICH ARE GREATED THAN RESTRICTED 
C                      DISTANCE FROM PREVIOUS PEAKS ARE WRITTEN INTO 
C                      DOCUMENT FILE
C       LUNDOC      LOGICAL UNIT NUMBER FOR DOCUMENT FILE
C
C      REGISTER POSITIONS 1= INTEGER X-SHIFT
C                         2= INTEGER Y-SHIFT
C                         3= INTEGER Z-SHIFT
C                         4= FLOATING X-SHIFT (CENTER OF GRAVITY)
C                         5= FLOATING Y-SHIFT (CENTER OF GRAVITY)
C                         6= FLOATING Z-SHIFT (CENTER OF GRAVITY)
C                         7= ABSOLUTE PEAK HEIGHT
C                         8= RATIO
C
C--*********************************************************************

         SUBROUTINE SPEAK3(LUN,NX,NY,NZ,OPTT,LUNDOC)

         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC'

         CHARACTER(LEN=MAXNAM)      :: DOCNAM
         CHARACTER                  :: OPTT,OPT,OPTB,OPTC
         LOGICAL                    :: CGR,RESTRICT
         INTEGER                    :: NCTR(3)

C        RUN TIME ARRAYS
         LOGICAL, ALLOCATABLE       :: KEEP(:)
         REAL, ALLOCATABLE          :: QBUF(:,:,:)
         REAL, ALLOCATABLE          :: RPC(:,:)
         REAL, ALLOCATABLE          :: PEAK(:)
         INTEGER, ALLOCATABLE       :: NPC(:,:)
         CHARACTER                  :: NULL = CHAR(0)

         OPT      = OPTT
         RESTRICT =  (OPT == 'R')

         CALL RDPRMC(OPTB,NC,.TRUE.,'MAXIMA(+) OR MINIMA(-)',NULL,IRT)
         IF (IRT .NE. 0) RETURN

         IF (OPTB == '-') THEN
            SIGN = -1.0
         ELSE
            SIGN = +1.0
         ENDIF

         ML  = 1
         NOR = 0
         CALL RDPRIS(ML,NOR,NOT_USED,
     &      'NUMBER OF PEAKS, CENTER ORIGIN OVERRIDE (0/1)',IRT)
         IF (IRT .NE. 0) RETURN

         ML = MAX(ML,1)

         ELIPX = 1.0
         ELIPY = 1.0
         ELIPZ = 0.0

         IF (RESTRICT) THEN
           CALL RDPRM3S(ELIPX,ELIPY,ELIPZ,NOT_USED,
     &       'X, Y, & Z RADII OF OF EXCLUDED NEIGHBORHOOD',IRT)
           IF (IRT .NE. 0) GOTO 9000
           CGR = .FALSE. 
         ELSE
            CALL RDPRMC(OPTC,NC,.TRUE.,
     &              'CENTER OF GRAVITY CALCULATION? (Y/N)',NULL,IRT)
            IF (IRT .NE. 0) GOTO 9000

            CGR = (OPTC == 'Y') 
            IF (CGR) THEN
               CALL RDPRM3S(ELIPX,ELIPY,ELIPZ,NOT_USED,
     &                'X, Y, & Z RADII OF ELLIPSES',IRTFLG)
               IF (IRTFLG .NE. 0) GOTO 9000

               IF (ELIPZ <= 0) THEN
                  CALL RDPRM1S(ELIPZ,NOT_USED,
     &                     'Z-RADIUS OF ELLIPSES',IRTFLG)
                  IF (IRTFLG .NE. 0) GOTO 9000
               ENDIF
            ENDIF
         ENDIF

         IF (NOR .NE. 0)  THEN
            NCTR(1) =  0.0
            NCTR(2) =  0.0
            NCTR(3) = -1.0
            CALL RDPRI3S(NCTR(1),NCTR(2),NCTR(3),NOT_USED,
     &                'X, Y, & Z ORIGIN COORDINATES',IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9000

            IF (NCTR(3) < 0) THEN
              CALL RDPRI1S(NCTR(3),NOT_USED,
     &                     'Z ORIGIN COORDINATE',IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9000
            ENDIF

            CALL RDPRI1S(NTAB,NOT_USED,
     &                 'PEAK NUMBER FOR RATIO',IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9000

            NTAB = MAX(NTAB,1)
            IF (NTAB > ML)   THEN
               CALL ERRT(102,'PEAK NOT CONTAINED IN TABLE',NTAB)
               GOTO 9000
            ENDIF
         ELSE
            NTAB  = 1
            NCTR(1) = NX   / 2 + 1
            NCTR(2) = NY   / 2 + 1
            NCTR(3) = NZ / 2 + 1
         ENDIF

         CALL RDPRMC(OPTB,NC,.TRUE.,'BOX SELECTION?(Y/N)',NULL,IRT)
         IF (IRT .NE. 0) GOTO 9000

         IF (OPTB == 'Y') THEN
            CALL RDPRMI(NSA1,NSA2,NOT_USED,'LOWER, UPPER SAMPLE')
            CALL RDPRMI(NRO1,NRO2,NOT_USED,'LOWER, UPPER ROW')
            CALL RDPRMI(NSL1,NSL2,NOT_USED,'LOWER, UPPER SLICE')
         ELSE
            NSL1 = 1
            NSL2 = NZ
            NRO1 = 1
            NRO2 = NY
            NSA1 = 1
            NSA2 = NX
         ENDIF

         ALLOCATE (QBUF(NX,NY,3),PEAK(ML), 
     &             NPC(3,ML),RPC(3,ML),KEEP(ML), STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN 
            MWANT = NX*NY*3 + 8*ML
            CALL ERRT(46,'KEEP,...',MWANT)
            GOTO 9000
         ENDIF

C        FIND PEAKS
         CALL PKSR3(LUN,QBUF,NX,NY,NZ,NSA1,NSA2,NRO1,NRO2,
     &              NSL1,NSL2,SIGN, PEAK, NPC, RPC, ML,NMAX)

         IF (NMAX .LE. 0)  THEN
            WRITE(NOUT,*)  ' *** No peaks found'
            GOTO 9000
         ENDIF


C        RESTRICT CALCULATION:
         IF (RESTRICT) THEN

C            INITIALIZE KEEP ARRAY
             KEEP = .TRUE.

             DO I = 1,NMAX-1
                XI   = RPC(1,I) 
                YI   = RPC(2,I) 
                ZI   = RPC(3,I) 

                DO J = I+1,NMAX
                   IF (KEEP(J)) THEN
C                     DETERMINE IF PEAK: J IS WITHIN ELIPSOID AROUND: I
                      XJ   = RPC(1,J) 
                      YJ   = RPC(2,J) 
                      ZJ   = RPC(3,J)

                      RZI = ((ZI - ZJ) / ELIPZ)**2
                      RYI = ((YI - YJ) / ELIPY)**2 + RZI
                      REL = ((XI - XJ) / ELIPX)**2 + RYI
 
C                     write(6,*) '(',i,',',j,')',xi,yi,zi, xj,yj,zj, rzi,ryi,rel
                      IF (REL <= 1.0) THEN
C                        INSIDE ELIPSOID, DISCARD J PEAK
                         KEEP(J) = .FALSE.
                         CYCLE
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO

         ENDIF

C        CENTER OF GRAVITY CALCULATION:
         IF (CGR) THEN
            CALL CGR_3(LUN,QBUF,NX,NY,NZ,
     &                 ELIPX,ELIPY,ELIPZ,NPC,
     &                 RXCEN,RYCEN,RZCEN,RSUM)

            IF (RSUM .NE. 0.0)  THEN
               WRITE(NOUT,298)
298            FORMAT(/,'  The x,y,z coordinates of first peak ',
     &                  'replaced by center of gravity ',
     &                  'approximation within area specified.')
               RPC(1,1) = RXCEN
               RPC(2,1) = RYCEN
               RPC(3,1) = RZCEN
            ELSE
               WRITE(NOUT,297)
297            FORMAT(/,'  Center of gravity approximation within ',
     &                  'area specified cannot be calculated - ',
     &                  'negative values encountered.')
            ENDIF
         ENDIF

         LUNDOCRET = 0
         IF (OPT == 'R' .OR. OPT == 'D') THEN
C           OPEN DOC FILE HERE TO AVOID MSG IN RESULTS OUTPUT
            CALL OPENDOC(DOCNAM,.TRUE.,NLET,LUNDOC,LUNDOCRET,.TRUE.,
     &              'OUTPUT DOC.',.FALSE.,.TRUE.,.TRUE.,NEWFILE,IRTFLG)

            IF (IRTFLG .NE. 0) OPT  = ' '  ! NO DOC FILE WANTED
          ENDIF

         MLOUT = MIN(ML,NMAX)
         IF (MLOUT < ML) THEN
            WRITE(NOUT,90)MLOUT
90          FORMAT(/,'  ONLY: ',I5,' PEAKS FOUND')
         ENDIF

C        WRITE TO RESULTS AND DOC FILE

         CALL NEWROUT(LUNDOCRET,OPT,PEAK,NPC,RPC,NCTR,MLOUT,NTAB,KEEP)
         CLOSE(LUNDOC)

C        SET REGISTER VALUES (IF DESIRED) FOR FIRST PEAK
         CALL REG_SET_NSEL(1,5,FLOAT(NPC(1,1)),FLOAT(NPC(2,1)),
     &                         FLOAT(NPC(3,1)),RPC(1,1),RPC(2,1),
     &                         IRTFLG)

         CALL REG_SET_NSEL(6,2,RPC(3,1),PEAK(1),0.0, 0.0, 0.0, IRTFLG)

9000     IF (ALLOCATED(KEEP)) DEALLOCATE(KEEP)
         IF (ALLOCATED(PEAK)) DEALLOCATE(PEAK)
         IF (ALLOCATED(NPC))  DEALLOCATE(NPC)
         IF (ALLOCATED(RPC))  DEALLOCATE(RPC)
         IF (ALLOCATED(QBUF)) DEALLOCATE(QBUF)

         RETURN

         END

C++**************************** NEWROUT *******************************
C
C  PURPOSE:  WRITE PEAK PARAMETERS TO RESULTS AND DOC FILE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE NEWROUT(LUNDOC,OPT,PEAK,NPC,RPC,NCTR,ML,NTAB,KEEP)

         INCLUDE 'CMBLOCK.INC'

         INTEGER                 :: NPC(3,ML)  
         INTEGER                 :: NCTR (3)
         REAL                    :: RPC(3,ML) 
         REAL                    :: PEAK(ML) 
         LOGICAL                 :: KEEP(ML)
         CHARACTER               :: OPT

         REAL                    :: DLIST(9)
         LOGICAL                 :: RESTRICT
         CHARACTER               :: CDUM
         CHARACTER(LEN=120)      :: COMMENT
         INTEGER                 :: IRTFLG

         RESTRICT = (OPT == 'R')

C                   123456789 123456789 123456789 123456789 12345678
C                   9 123456789 123456789 123456789 123456789 1
C                   123456789 123456789 

         COMMENT = 'PEAK      PK-X          PK-Y          PK-Z     '//
     &             'SUB_PIX-X     SUB_PIX-Y     SUB_PIX-Z     ' //
     &             'PK-HEIGHT    HEIGHT-RATIO'

         IF (LUNDOC > 0) THEN
            CALL LUNDOCPUTCOM(LUNDOC,COMMENT,IRTFLG)
         ENDIF

         IF (VERBOSE) THEN

            WRITE(NOUT,*) ' '

C                      123456789 123456789 123456789 12345
C                      6789 123456789 123456789 123456789 123456789
C                      123456789 1
            COMMENT = ' NO  NX-O NY-O  NZ-O NX   NY   NZ '       //
     &                '     X        Y        Z    PEAK     RATIO'

            WRITE(NOUT,*)COMMENT(1:77)
         ENDIF

C        PEAK HEIGHT RATIO
         RTA = PEAK(NTAB)
         IF (RTA == 0.0)  RTA = 1.0

         NEWN = 0
         DO N=1,ML
            IF (RESTRICT) THEN
C             DO NOT MERGE IF, ALLOC ERROR ON SOME SYSTEMS
              IF (.NOT. KEEP(N)) CYCLE
            ENDIF
            NEWN = NEWN + 1

            DO M=1,3
               NPC(M,N) = NPC(M,N) - NCTR(M)
               RPC(M,N) = RPC(M,N) - NCTR(M)
            ENDDO

            IF (VERBOSE) THEN
               WRITE(NOUT,701) NEWN,
     &            (NPC(M,N)+NCTR(M),M=1,3),
     &            (NPC(M,N),M=1,3),
     &            (RPC(M,N),M=1,3), 
     &             PEAK(N), 
     &             PEAK(N) / RTA
701            FORMAT(1X,I3,6I5,3(1X,F8.2),2(1X,G9.2))
            ENDIF

            IF (OPT == 'D' .OR. RESTRICT) THEN
               DLIST(1) = NEWN
               DLIST(2) = NPC(1,N)
               DLIST(3) = NPC(2,N)
               DLIST(4) = NPC(3,N)
               DLIST(5) = RPC(1,N)
               DLIST(6) = RPC(2,N)
               DLIST(7) = RPC(3,N)
               DLIST(8) = PEAK(N)
               DLIST(9) = PEAK(N)/RTA

               CALL LUNDOCWRTDAT(LUNDOC,NEWN,DLIST(2),8,IRTFLG)
            ENDIF
         ENDDO

         IF (RESTRICT .AND. VERBOSE) THEN
            WRITE(NOUT,90) NEWN,ML
90          FORMAT(/,'  Retained: ',I8,' peaks out of: ',I8)
         ENDIF

         IF (VERBOSE) WRITE(NOUT,*) ' '

         END

