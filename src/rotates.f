C++*********************************************************************
C                                                                      *
C  ROTATES.F   NEW                               SEP 2011 ARDEAN LEITH *  
C              FBS VOLUME INTERP                 OCT 2011 ARDEAN LEITH *
C              ROT2QS --> RTSQ RENAMED           DEC 2011 ARDEAN LEITH *
C              NSAM --> NX                       JAN 2012 ARDEAN LEITH *
C              ROT3L CALL ADDED                  JAN 2012 ARDEAN LEITH *
C              ROT 2D BUG                        JUN 2012 ARDEAN LEITH *
C              ROTATESZ FBS BUF BUG              SEP 2012 ARDEAN LEITH *
C              ROT32 IOFF BUG                    DEC 2013 ARDEAN LEITH *
C ********************************************************************** 
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
C  ROTATES ()
C
C  CALLS:
C
C  For images (Linear Interp)
C    ROT32(LUN1,LUN2,NX,IREC1,IREC2,1,TH,BACK,SHX,SHY) 
C  For images (Quadratic Interp)
C    RTSQ(BUF,BUFLIN, NX,NY, THETA,SCLI,SHXI,SHYI, IREC1,LUN2)
C  For images (FBS Interp)
C    RTSF(BUF,BUF2, NXLD,NX,NY, THETA,SCLI,SHXI,SHYI, IRTFLG)
C
C  For volumes (All)
C    ROTATES3(LUN1,LUN2, NX,NY,NZ, ANGLE, THETA,PHI,PSI, 
C           MODE, NCX,NCY,NCZ,  P1,P2,INTERPC, BACKC,BACK,IRTFLG)
C
C  For volumes around point (Linear Interp)
C    ROTATES3L(LUN2,BUF,KLX,KNX,KLY,KNY,KLZ,KNZ, DRM, BACKC,BACK)
C  For volumes around point (Quadratic Interp)
C    ROTATES3Q(LUN2,Q1,KLX,KNX,KLY,KNY,KLZ,KNZ, DRM, BACKC,BACK) 
C  For volumes around point (FBS Interp)
C    ROTATES3FBS(LUN2,BUF,KLX,KNX,KLY,KNY,KLZ,KNZ,
C             DRM,USEBACK,BACK,IRTFLG)
C
C  For volumes around line (Quadratic Interp)
C    ROTATESL3Q(LUN2,BUF,NX,NY,NZ,P1,P2,ANGLE, BACKC,BACK)
C  For volumes around line (Linear Interp)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12 
C--********************************************************************* 

      SUBROUTINE ROTATES() 

      IMPLICIT NONE
      INCLUDE 'CMBLOCK.INC' 
      INCLUDE 'CMLIMIT.INC' 

      CHARACTER (LEN=MAXNAM) :: FILNAM,BACKC,BACKTIL 
      CHARACTER (LEN=1)      :: CDUM,NULL,MODE,INTERPC 
      CHARACTER (LEN=40 )    :: MSG 

      REAL                   :: AV1,FMAX1,FMIN1
      REAL                   :: BACK,SHX,SHY
      REAL                   :: PHI,THETA,PSI,ANGLE
      REAL                   :: P1(3),P2(3)

      INTEGER                :: ICOMM,MYPID,MPIERR
      INTEGER                :: MAXIM1,MAXIM2,IRTFLG
      INTEGER                :: NX,NY,NZ,ITYPE,NLET,NOT_USED
      INTEGER                :: IREC1,IREC2,NXLD,IRECOFF
      INTEGER                :: IMAMI1,ISLICE,NCX,NCY,NCZ,IDUM,MWANT
      REAL                   :: TH
      DOUBLE PRECISION       :: DRM(3,3)
      LOGICAL                :: USEBACK

      REAL,    ALLOCATABLE   :: BUFLIN(:),BUF(:,:),BUF2(:,:)

      INTEGER,PARAMETER      :: LUN1 = 50 
      INTEGER,PARAMETER      :: LUN2 = 51 

C     DATA FUNC/ '14'/ ROT

      CALL SET_MPI(icomm,mypid,mpierr)  ! SET MYPID

      NULL   = CHAR(0)
      MODE   = FCHAR(4:4)

C     OPEN INPUT FILE 
      MAXIM1 = 0 
      CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &             NX,NY,NZ, 
     & 		   MAXIM1,'INPUT',.FALSE.,IRTFLG) 
      IF (IRTFLG .NE. 0) GOTO 9999 
 
      IF (IMAMI == 0) CALL NORM3(LUN1,NX,NY,NZ,FMAX,FMIN,AV) 
C     RECORD INPUT IMAGE  AVERAGE IN CASE NEEDED FOR BACKGROUND 
      AV1    = AV 
      FMAX1  = FMAX 
      FMIN1  = FMIN 
      IMAMI1 = IMAMI 

C     OPEN OUTPUT FILE WITH SAME DIMENSIONS AS INPUT FILE 
      MAXIM2 = 0 
      CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &             NX,NY,NZ, 
     &             MAXIM2,'OUTPUT',.FALSE.,IRTFLG) 
      IF (IRTFLG .NE. 0) GOTO 9999 

 
      IF (ITYPE == 1) THEN 
C        INPUT IMAGE ------------------------------------------ IMAGE

C        SET ROTATION ANGLE 
         CALL RDPRM1S(ANGLE,NOT_USED,'ROTATION ANGLE',IRTFLG) 
         TH = ANGLE * DATAN(1.0D0) / 45.0D0

C        SET ROTATION CENTER (ARBITRARY)
         SHX = 0.0 
         SHY = 0.0 

         IF (FCHAR(4:4) == 'A') THEN 
C           ROTATE AROUND AN ARBITRARY CENTER  
            CALL RDPRM2S(SHX,SHY,NOT_USED,
     &                   'X & Y CENTER OFFSETS',IRTFLG) 
            IF (IRTFLG .NE. 0) GOTO 9999 
         ENDIF 

C        SET INTERPOLATION TYPE 
         CALL RDPRMC(INTERPC,NLET,.TRUE.,
     &       'LINEAR, QUADRATIC, OR FBS INTERPOLATION (L,Q,F)',
     &               NULL,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999 
 


         IF ((INTERPC .NE. 'L') .AND. 
     &       (INTERPC .NE. 'Q') .AND. 
     &       (INTERPC .NE. 'F')) THEN
C                      123456789 1123456789 23456789 
                MSG = 'UNKNOWN INTERPOLATION TYPE: ' // INTERPC
                CALL ERRT(101,MSG,IDUM)
            GOTO 9999
         ENDIF

C        SET BACKGROUND VALUE 
         CALL RDPRMC(BACKC,NLET,.TRUE.,
     &        'UNROTATED, AVG, MIN, OR SPECIFIED CORNERS (U,A,M,2.5)',
     &         NULL,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999 

         BACK = HUGE(BACK)
         IF (INDEX(BACKC(1:NLET),'A') > 0) BACK = AV1
         IF (INDEX(BACKC(1:NLET),'M') > 0) BACK = FMIN1
         IF (INDEX(BACKC(1:NLET),'U') > 0) BACK = 0.0 ! UNUSED
        
         IF (BACK == HUGE(BACK)) THEN
            BACKTIL = '~' // BACKC(1:NLET)
            CALL RDPRM1S(BACK,NOT_USED,BACKTIL,IRTFLG) 
            IF (IRTFLG .NE. 0) GOTO 9999 
         ENDIF
         USEBACK = (BACKC(1:1) .NE. 'U')

         IF (INTERPC == 'Q' ) THEN
            ALLOCATE(BUF(NX,NY),BUFLIN(NX), STAT=IRTFLG)
            IF (IRTFLG .NE. 0) THEN 
               MWANT = NX + NX*NY*NZ
               CALL ERRT(46,'ROTATES; BUF',MWANT)
               GOTO 9999
            ENDIF  
         ELSEIF (INTERPC == 'F' ) THEN
            NXLD   = NX + 2 - MOD(NX,2)
            ALLOCATE (BUF(NXLD, NY),BUF2(NX,NY), STAT=IRTFLG)
            IF (IRTFLG .NE. 0) THEN 
               MWANT = NXLD*NY + NX*NY 
               CALL ERRT(46,'ROTATES; BUF..',MWANT)
               GOTO 9999
            ENDIF  
         ENDIF

C       ROTATE IMAGE IN-CORE 

        IF (INTERPC == 'L')     THEN
C          LINEAR INTERP ('RT') 

C          READS OWN IMAGE IN ROT32 WRITES OUTPUT IMAGE
           CALL ROT32(LUN1,LUN2,NX, 1,NY,1,
     &                TH,BACK,SHX,SHY,0)

        ELSEIF (INTERPC == 'Q') THEN
C          QUADRATIC INTERP ('RT SQ') 

C          LOAD IMAGE INTO BUF
           CALL REDVOL(LUN1,NX,NY, 1,1, BUF,IRTFLG)

C          ROTATE WITH RTSQ AND WRITE OUTPUT IMAGE
     	   CALL RTSQ_BACK(BUF,BUFLIN, NX,NY,
     &                    ANGLE, 1.0,0.0,0.0, 
     &                    USEBACK,BACK, 0,LUN2)

        ELSEIF (INTERPC == 'F') THEN
C          FBS INTERP ('RT SF') 

C          LOAD INTO FFT X PADDED BUF
           CALL READV(LUN1,BUF,NXLD,NY, NX,NY,1)

C          ROTATE WITH RTSF 
           CALL RTSF_BACK(BUF,BUF2, NXLD, NX, NY, 
     &                   ANGLE, 1.0,0.0,0.0, 
     &                   USEBACK,BACK, IRTFLG)

C          WRITE UNPADDED OUTPUT IMAGE
           CALL WRITEV(LUN2,BUF2,NX,NY,NX,NY,1)
         ENDIF 

      ELSE

C        INPUT VOLUME ------------------------------------------ VOLUME
 
         IF (MODE == 'L' ) THEN
C           ROTATE VOLUME AROUND LINE, GET ANGLE

            ANGLE = 0.0
            CALL RDPRM1S(ANGLE,NOT_USED,
     &                  'ROTATION ANGLE',IRTFLG)
            IF (IRTFLG .NE. 0) RETURN

C           GET LINE 
            P1(1) = NX   / 2 + 1
            P1(2) = NY   / 2 + 1
            P1(3) = 1
            CALL  RDPRM3S(P1(1),P1(2),P1(3),NOT_USED,
     &           'X, Y, & Z FOR FIRST  POINT ON ROTATION AXIS',
     &           IRTFLG)
            IF (IRTFLG .NE. 0)  GOTO 9999

            P2(1) = NX   / 2 + 1
            P2(2) = NY   / 2 + 1
            P2(3) = NZ  
            CALL  RDPRM3S(P2(1),P2(2),P2(3),NOT_USED,
     &           'X, Y, & Z FOR SECOND POINT ON ROTATION AXIS',
     &           IRTFLG)
            IF (IRTFLG .NE. 0)  GOTO 9999

            CALL RDPRMC(INTERPC,NLET,.TRUE.,
     &           'LINEAR OR QUADRATIC INTERPOLATION (L,Q)',
     &            NULL,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999
 
            IF (INTERPC .NE. 'L' .AND. 
     &          INTERPC .NE. 'Q') THEN
C                      123456789 1123456789 23456789 
                MSG = 'UNKNOWN INTERPOLATION TYPE: ' // INTERPC
                CALL ERRT(101,MSG,IDUM)
               GOTO 9999
            ENDIF

         ELSE
C           ROTATE VOLUME BY EULER ANGLES  

C           GET EULER ANGLES
            PHI   = 0.0
            THETA = 0.0
            PSI   = 0.0
            ANGLE = 0.0

            CALL RDPRM3S(PHI,THETA,PSI,NOT_USED,
     &                  'ROTATION ANGLES (PHI, THETA, & PSI)',IRTFLG)
            IF (IRTFLG .NE. 0) RETURN

            IF (MODE == 'A' ) THEN
               CALL RDPRI3S(NCX,NCY,NCZ,NOT_USED,
     &                 'X, Y & Z FOR CENTER OF ROTATION',IRTFLG)
               IF (IRTFLG .NE. 0) GOTO 9999
            ENDIF

            CALL RDPRMC(INTERPC,NLET,.TRUE.,
     &           'TRI- LINEAR, QUADRATIC  OR FBS INTERPOLATION (L,Q,F)',
     &            NULL,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999 
            IF (INTERPC .NE. 'L' .AND. 
     &          INTERPC .NE. 'Q' .AND.
     &          INTERPC .NE. 'F') THEN
C                      123456789 1123456789 23456789 
                MSG = 'UNKNOWN INTERPOLATION TYPE: ' // INTERPC
                CALL ERRT(101,MSG,IDUM)
               GOTO 9999
            ENDIF
         ENDIF 

         !write(6,*) ' mode & fcahr:',MODE,INTERPC

         IF (MODE .NE. 'L' .OR. INTERPC .NE. 'L') THEN
C           SET BACKGROUND VALUE 
            CALL RDPRMC(BACKC,NLET,.TRUE.,
     &         'UNROTATED, AVG, MIN, OR SPECIFIED CORNERS (U,A,M,2.5)',
     &            NULL,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999 

            BACK = HUGE(BACK)
            IF (BACKC(1:1) == 'U')            BACK = 0.0 !(unused)
            IF (INDEX(BACKC(1:NLET),'A') > 0) BACK = AV1
            IF (INDEX(BACKC(1:NLET),'M') > 0) BACK = FMIN1

            IF (BACK == HUGE(BACK)) THEN
               BACKTIL = '~' // BACKC(1:NLET)
               CALL RDPRM1S(BACK,NOT_USED,BACKTIL,IRTFLG) 
               IF (IRTFLG .NE. 0) GOTO 9999 
            ENDIF
         ENDIF

         IF ((MODE .NE. 'L') .AND.
     &        PHI .NE. 0 .AND. THETA == 0.0 .AND. PSI == 0) THEN
            !write(6,*) '  Using rotatesz:',interpc,useback,back

            SHX     = 0.0 
            SHY     = 0.0 
            ANGLE   = PHI
            USEBACK = (BACKC(1:1) .NE. 'U')

            CALL ROTATESZ(LUN1,LUN2, NX,NY,NZ,
     &                    ANGLE, SHX,SHY,
     &                    INTERPC,USEBACK, BACK,   
     &                    IRTFLG) 
         ELSE

            CALL ROTATES3(LUN1,LUN2, NX,NY,NZ,
     &                    ANGLE, THETA,PHI,PSI, 
     &                    MODE,  NCX,NCY,NCZ,  P1,P2,
     &                    INTERPC, BACKC,BACK,IRTFLG)
         ENDIF
      ENDIF

9999  CONTINUE

C     SET HEADER FOR ALTERATIONS IN IMAGE DUE TO OPERATIONS 
      CALL SETPRMB(LUN2, 0.0,0.0, 0.0,0.0) 

      IF (ALLOCATED(BUF))    DEALLOCATE(BUF)
      IF (ALLOCATED(BUF2))   DEALLOCATE(BUF2)
      IF (ALLOCATED(BUFLIN)) DEALLOCATE(BUFLIN)

9998  CLOSE(LUN1)
      CLOSE(LUN2)

      END


C     -------------  ROTATES3 ----------------------- ROTATES3

        SUBROUTINE ROTATES3(LUN1,LUN2,NX,NY,NZ,
     &                       ANGLE, THETA,PHI,PSI, 
     &                       MODE,   NCX,NCY,NCZ,  P1,P2,
     &                       INTERPC, BACKC,BACK,IRTFLG) 

        IMPLICIT NONE    
        INCLUDE 'CMBLOCK.INC'
  
        INTEGER           :: LUN1,LUN2,NX,NY,NZ
        REAL              :: ANGLE, THETA,PHI,PSI 
        CHARACTER(LEN=1)  :: MODE,INTERPC,BACKC
        INTEGER           :: NCX,NCY,NCZ 
        REAL              :: P1(3),P2(3)
        REAL              :: BACK 
        INTEGER           :: IRTFLG 

        REAL, ALLOCATABLE :: BUF(:,:,:)
        INTEGER           :: KNX,KNY,KNZ,KLX,KLY,KLZ
        INTEGER           :: I
        DOUBLE PRECISION  :: DRM(3,3)
        LOGICAL           :: USEBACK


C       3-D  ROTATION
        USEBACK = (BACKC(1:1) .NE. 'U')

        !write(6,*) ' In rotates3a; mode & interpc:',mode,interpc

        ALLOCATE (BUF(NX,NY,NZ), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'ROTATES3A, BUF',NX*NY*NZ)
           RETURN
        ENDIF

        CALL REDVOL(LUN1,NX,NY,1,NZ,BUF,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (MODE == 'A')  THEN
C          USER SPECIFIES CENTER OF ROTATION
           KNX = NX - NCX
           KNY = NY - NCY
           KNZ = NZ - NCZ

           KLX = -NCX + 1
           KLY = -NCY + 1
           KLZ = -NCZ + 1

C          CALCULATE ROTATION MATRIX
           CALL BLDR(DRM,PSI,THETA,PHI)

           !write(6,*) ' specified rotation center:',PSI,THETA,PHI

        ELSEIF (MODE .NE. 'L') THEN
C          DEFAULT IS SPIDER CENTER OF ROTATION
           KNX = NX   / 2 - 1 + MOD(NX,2)
           KNY = NY   / 2 - 1 + MOD(NY,2)
           KNZ = NZ / 2 - 1 + MOD(NZ,2)

           KLX = -NX   / 2
           KLY = -NY   / 2
           KLZ = -NZ / 2

C          CALCULATE ROTATION MATRIX
           CALL BLDR(DRM,PSI,THETA,PHI)
           !write(6,*) ' spider rotation center:',PSI,THETA,PHI
        ENDIF

        IF ((MODE == 'L' .AND. ANGLE == 0.0) .OR.
     &     (THETA == 0.0 .AND. PHI   == 0.0 .AND. PSI == 0.0)) THEN
C          NO ROTATION NEEDED
           !write(6,*) ' Calling wrtvol:',psi,theta,phi
           CALL WRTVOL(LUN2,NX,NY,1,NZ,BUF,IRTFLG)
           IRTFLG = 0

        ELSEIF (MODE == 'L')  THEN
C          USER SPECIFIES LINE OF ROTATION, USE LIN/QUAD INTERP

           IF (INTERPC == 'L') THEN
C             ROTATE ABOUT LINE, TRI-LINEAR INTERPOLATION 
              !write(6,*) ' Linear about line:',angle
              CALL ROTL3(LUN2,BUF, NX,NY,NZ, P1,P2,ANGLE)

           ELSEIF (INTERPC == 'Q') THEN
C             ROTATE ABOUT LINE, TRI-QUADRATIC INTERPOLATION 
              !write(6,*) ' Quad about line:',angle
              CALL ROTATESL3Q(LUN2,BUF, NX,NY,NZ,
     &                     P1,P2,ANGLE, BACKC,BACK)
           ENDIF
           IRTFLG = 0


        ELSEIF (INTERPC == 'L') THEN
C          ROTATE ABOUT POINT, TRI-LINEAR INTERPOLATION 

           !write(6,*) ' Linear about point:',PSI,THETA,PHI
           CALL ROTATES3L(LUN2,BUF,KLX,KNX,KLY,KNY,KLZ,KNZ,
     &                    DRM, BACKC,BACK)
           IRTFLG = 0

        ELSEIF (INTERPC == 'Q') THEN
C          ROTATE ABOUT POINT, TRI-QUADRATIC INTERPOLATION
           !write(6,*) ' Quad about point:',PSI,THETA,PHI

           CALL ROTATES3Q(LUN2,BUF,KLX,KNX,KLY,KNY,KLZ,KNZ,
     &                    DRM, BACKC,BACK)
           IRTFLG = 0
        ELSE
C          ROTATE ABOUT POINT, FBS INTERPOLATION
           !write(6,*) ' FBS about point:',PSI,THETA,PHI
           CALL ROTATES3FBS(LUN2,BUF,KLX,KNX,KLY,KNY,KLZ,KNZ,
     &                      DRM,USEBACK,BACK,IRTFLG)
        ENDIF

9999    IF (ALLOCATED(BUF)) DEALLOCATE(BUF)
  
        END



C       -------------  ROTATESZ ----------------------- ROTATESZ

        SUBROUTINE ROTATESZ(LUN1,LUN2, NX,NY,NZ,
     &                      ANGLE, SHX,SHY,
     &                      INTERPC,USEBACK, BACK,   
     &                      IRTFLG) 

        IMPLICIT NONE    
        INCLUDE 'CMBLOCK.INC'
  
        INTEGER           :: LUN1,LUN2, NX,NY,NZ
        REAL              :: ANGLE, SHX,SHY
        CHARACTER(LEN=1)  :: INTERPC
        LOGICAL           :: USEBACK
        REAL              :: BACK
        INTEGER           :: IRTFLG 

        REAL, ALLOCATABLE :: BUF(:,:),BUFLIN(:),BUF2(:,:)
        REAL              :: TH 

        INTEGER           :: ISLICE,IREC1,IREC2,NXLD,IRECOFF,MWANT,IROW
        INTEGER           :: IOFF

C       ROTATION OF VOLUME AROUND Z AXIS
        TH = ANGLE * DATAN(1.0D0) / 45.0D0

        IF (INTERPC == 'Q') THEN
C       ROTATION OF VOLUME AROUND Z AXIS
 
           ALLOCATE (BUF(NX,NY), BUFLIN(NX),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN 
              MWANT = NX*NY + NX
              CALL ERRT(46,'ROTATESZ, BUF..',MWANT)
              RETURN
           ENDIF

        ELSEIF (INTERPC == 'F') THEN
           NXLD = NX + 2 - MOD(NX,2)
           ALLOCATE (BUF(NXLD,NY), BUF2(NX,NY),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN 
              MWANT = NXLD*NY + NX*NY
              CALL ERRT(46,'ROTATESZ, BUF..',MWANT)
              RETURN
           ENDIF
         ENDIF

         !write(6,*) 'linear interp :',interpc,shx,shy
         !write(6,*) 'angle,back,nz:',angle,back,nz

         DO ISLICE = 1,NZ 
C           ROTATE SLICE IN-CORE 

            IF (INTERPC == 'L')     THEN
C              LINEAR INTERP ('RT')

C              READS IMAGE IN ROT32, WRITES OUTPUT SLICE
               IREC1 = (ISLICE-1)*NY+1 
               IREC2 = IREC1 + NY - 1 
               IOFF  = IREC1 - 1 
               !write(6,*) 'irec1,irec2,ioff:',irec1,irec2,ioff

               CALL ROT32(LUN1,LUN2,NX,IREC1,IREC2,1,
     &                    TH,BACK,SHX,SHY,IOFF)

            ELSEIF (INTERPC == 'Q') THEN
C              QUADRATIC INTERP ('RT SQ')

C              LOAD IMAGE INTO BUF
               CALL REDVOL(LUN1,NX,NY, ISLICE,ISLICE, BUF,IRTFLG)

C              ROTATE WITH RTSQ AND WRITE OUTPUT SLICE
               IRECOFF = (ISLICE-1) * NY 
   	       CALL RTSQ_BACK(BUF,BUFLIN, NX,NY,
     &                     ANGLE, 1.0,0.0,0.0, 
     &                     USEBACK,BACK,  IRECOFF,LUN2)

            ELSEIF (INTERPC == 'F') THEN
C              FBS INTERP ('RT SF')

C              LOAD CURRENT SLICE INTO FFT X PADDED BUF
               DO IROW=1,NY
                  CALL REDLIN(LUN1,BUF(1,IROW),
     &                        NX,IROW+(ISLICE-1)*NX)
               ENDDO

C              ROTATE WITH RTSF 
               CALL RTSF_BACK(BUF,BUF2, NXLD, NX, NY, 
     &                        ANGLE, 1.0,0.0,0.0, 
     &                        USEBACK,BACK, IRTFLG)

C              WRITE OUTPUT SLICE
               CALL WRTVOL(LUN2,NX,NY, ISLICE,ISLICE,BUF2,IRTFLG)
            ENDIF
         ENDDO 

9999     IF (ALLOCATED(BUF))    DEALLOCATE(BUF)
         IF (ALLOCATED(BUFLIN)) DEALLOCATE(BUFLIN)
         IF (ALLOCATED(BUF2))   DEALLOCATE(BUF2)

         END
