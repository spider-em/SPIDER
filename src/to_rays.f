C **********************************************************************
C
C  TO_RAYS     NEW                                 AUG 13 ARDEAN LEITH
C
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
C TO_RAYS
C
C PURPOSE:  CREATE POLAR REPRESENTATION OF IMAGE WITH RAYS ALONG THE
C           X DIMENSION.
C
C **********************************************************************

        SUBROUTINE TO_RAYS

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        CHARACTER(LEN=MAXNAM) :: FILNAM
        CHARACTER(LEN=72)     :: FORMOUT
        CHARACTER(LEN=96)     :: COMMENT

        REAL, ALLOCATABLE     :: BUF(:,:)
        REAL, ALLOCATABLE     :: OUT(:,:)

        CHARACTER(LEN=3)      :: MODE
        LOGICAL               :: HALFCIRC
        INTEGER               :: MAXIM,ITYPE,NX,NY,NZ,IRTFLG,IER,NLET
        INTEGER               :: NA,NRAYS,NR,IXC,IYC,I,J
        INTEGER               :: NRAD,NXP,NYP,IRAY,IRAD
        REAL                  :: PI,DFI,FI,XS,YS,VTMP,VMIN
        REAL                  :: UNUSED 
 
        INTEGER               :: NMIN,NMAX

        REAL                  :: quadri

        INTEGER, PARAMETER    :: LUNIN  = 20
        INTEGER, PARAMETER    :: LUNOUT = 21

        CHARACTER(LEN=1)      :: NULL   = CHAR(0)


C       GET NAME FOR INPUT FILE
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNIN,'O',ITYPE,NX,NY,
     &              NZ,MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (NZ > 1)  THEN
           CALL ERRT(101,'DOES NOT WORK ON VOLUMES',IER)
           RETURN
        ENDIF

C       GET NAME FOR POLAR OUTPUT FILE
        CALL FILERD(FILNAM,NLET,NULL,'OUTPUT',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        CALL RDPRMC(MODE,NA,.TRUE.,
     &        'FULL OR HALF CIRCLE (F/H)',NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        HALFCIRC  = (INDEX(MODE,'H') > 0)
        PI        = 4*DATAN(1.0D0)
        
C       FIND DEFAULT LENGTH OF RAYS = RADIUS = X
        NRAD = MIN(((NX-1)/2), ((NY-1)/2))  ! LENGTH OF RAY = X SIZE

C       FIND DEFAULT NUMBER OF RAYS = CIRCUMFERANCE = Y
        NRAYS = INT(2*PI*NRAD)                ! LEN. OF CIR = # OF RAYS
        IF (HALFCIRC) NRAYS = INT(PI*NRAD)    ! USED FOR POWER SPECTRA

        WRITE(NOUT,99)'  DEFAULT, RADIUS: ',NRAD,'   RAYS:',NRAYS
 99     FORMAT(A,I0,A,I0)

        CALL RDPRI2S(NRAD,NRAYS,UNUSED,
     &              'RADIUS & NUMBER OF RAYS',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        NXP  = NRAD + 1                       ! X SIZE
        NYP  = NRAYS                          ! Y SIZE

        DFI  = 2 * PI / NYP                   ! RADIANS 
        IF (HALFCIRC) DFI = PI / NYP

C       OPEN OUTPUT FILE
        MAXIM = 0
        ITYPE = 1
        CALL OPFILEC(0,.FALSE.,FILNAM,LUNOUT,'N',ITYPE,NXP,NYP,
     &               1,MAXIM,'RAY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        ALLOCATE (BUF(NX,NY),
     &            OUT(NXP,NYP),         ! NRAD + 1 x NRAYS
     &            STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'TO_RAYS; BUF, OUT',NX*NY + NXP*NYP + 2*NRAYS)
           GOTO 9999
        ENDIF

C       LOAD INPUT IMAGE
        CALL REDVOL(LUNIN,NX,NY, 1,1, BUF,IRTFLG)
           
        IXC     = NX/2+1     ! X CENTER
        IYC     = NY/2+1

        !write(6,*) ' Center:(',ixc,iyc,')  rad,rays:',nrad,nrays

C       CREATE POLAR REPRESENTATION
C       NOTE: ALGORITHM IS DEPENDENT ON CIRCULAR CLOSURE IN QUADRI!
        DO  J=0,NRAD       ! OVER ALL RADII
          DO I=1,NRAYS     ! OVER ALL POINTS ON CIRCLE 

             FI         = (I-1)   * DFI
             XS         = COS(FI) * J
             YS         = SIN(FI) * J

             !VTMP      = QUADRI(XS+IXC,YS+IYC,NX,NY,BUF)!*SQRT(REAL(J))
             VTMP       = QUADRI(XS+IXC,YS+IYC,NX,NY,BUF)

             OUT(J+1,I) = VTMP

             !if(i == 1) write(6,*) j,' : ',out(j+1,1)     
          ENDDO
        ENDDO


C       OUTPUT IMAGE -------------------------------------
        CALL WRTVOL(LUNOUT,NXP,NYP, 1,1, OUT,IRTFLG)
         
9999    IF (ALLOCATED(BUF))     DEALLOCATE (BUF)
        IF (ALLOCATED(OUT))     DEALLOCATE (OUT)

        CLOSE(LUNIN)
        CLOSE(LUNOUT) 
 
        END

C       --------------------- FINDRIDGES ------------------------

        SUBROUTINE FINDRIDGES

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        CHARACTER(LEN=MAXNAM) :: FILNAM,DOCNAM
        CHARACTER(LEN=72)     :: FORMOUT
        CHARACTER(LEN=96)     :: COMMENT

        REAL, ALLOCATABLE     :: BUF(:,:)
        REAL, ALLOCATABLE     :: CURVE(:)

        REAL                  :: DLIST(8)
        INTEGER               :: IRAY_SHORT,IRAY_LONG
        INTEGER               :: MAXIM,ITYPE,NX,NY,NZ,IRTFLG,IER,NLET
        INTEGER               :: NA,NR,IROW
        REAL                  :: ANGLE_L,ANGLE_S
        INTEGER               :: IRADRIDGESHORT,IRADRIDGELONG
        INTEGER               :: IASTIG,NMIN,NMAX
        INTEGER               :: IRADRIDGEMAX,IRADRIDGEMIN,IRADRIDGE
        REAL                  :: VALRIDG 

        LOGICAL               :: ADDEXT,GETNAME,ISOLD
        LOGICAL               :: APPEND,MESSAGE,NEWFILE
        INTEGER               :: LUNDOCN,LUNDOCNO,NLETD,KEY,NDIGITS
 
        INTEGER, PARAMETER    :: NUMMINS = 3 
        INTEGER, PARAMETER    :: NUMMAXS = 3

        INTEGER               :: LOCVALLEY(NUMMINS),LOCRIDGE(NUMMAXS)
        REAL                  :: VALVALLEY(NUMMINS),VALRIDGE(NUMMAXS)

        INTEGER, PARAMETER    :: LUNIN  = 20

        CHARACTER(LEN=1)      :: NULL   = CHAR(0)


C       OPEN INPUT FILE
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNIN,'O',ITYPE,NX,NY,
     &              NZ,MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (NZ > 1)  THEN
           CALL ERRT(101,'DOES NOT WORK ON VOLUMES',IER)
           RETURN
        ENDIF

        ALLOCATE (BUF(NX,NY),
     &            CURVE(NX), 
     &            STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'FINDRIDGES; BUF, CURVE',NX*NY + NX)
           GOTO 9999
        ENDIF

C       LOAD INPUT IMAGE
        CALL REDVOL(LUNIN,NX,NY, 1,1, BUF,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
           
        ADDEXT  = .TRUE.
        GETNAME = .TRUE.
        ISOLD   = .FALSE.
        APPEND  = .TRUE.
        MESSAGE = .FALSE. 
        IRTFLG  = -8         ! NO IC USE

        CALL OPENDOC(DOCNAM,ADDEXT,NLETD,LUNDOCN,LUNDOCNO,GETNAME,
     &         'OUTPUT DOC',ISOLD,APPEND,MESSAGE,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C                  123456789 123456789 123456789 123456789 123456789 123456789
        COMMENT = '  NUM,     S,          ANG_S,        V_S,         '//
     &               'L,           ANG_L,        V_L'
        !CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:80),IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       FIND RIDGES & VALLEY LOCATIONS ALONG ALL RAYS ------------- 

        IRAY_LONG    = 0
        IRAY_SHORT   = 0
        IRADRIDGEMAX = 0
        IRADRIDGEMIN = HUGE(IRADRIDGEMIN)

        !DO IROW=1,NY     ! LOOP OVER ALL ROWS 
        DO IROW=1+2,NY-2     ! LOOP OVER ALL ROWS (SKIP TOP & BOTTEM)

           CURVE = BUF(1:NX,IROW)

           CALL FINDMINMAX(CURVE,NX,
     &             LOCVALLEY,LOCRIDGE, VALVALLEY,VALRIDGE,
     &             3,3, NMIN,NMAX,IRTFLG)

C          FIND LONGEST AXIS OF CIRCLE RAY ---------------------------
           IRADRIDGE = LOCRIDGE(2)   ! RADIUS AT FIRST RIDGE
           VALRIDG   = VALRIDGE(2)   ! VALUE  AT FIRST RIDGE 

C          FIND LONGEST DISTANCE TO RIDGE ---------------------------
           IF (IRADRIDGE > IRADRIDGEMAX) THEN
              IRAY_LONG    = IROW
              IRADRIDGEMAX = IRADRIDGE
           ENDIF
   
C          FIND SHORTEST DISTANCE TO RIDGE ---------------------------
           IF (IRADRIDGE < IRADRIDGEMIN) THEN
              IRAY_SHORT   = IROW
              IRADRIDGEMIN = IRADRIDGE
           ENDIF
        ENDDO

        ANGLE_S = 180.0 * FLOAT(IRAY_SHORT) / FLOAT(NY)
        ANGLE_L = 180.0 * FLOAT(IRAY_LONG)  / FLOAT(NY)

        WRITE(6,90) ' SHORT RIDGE Y:',IRAY_SHORT,
     &              '   ANG:',        ANGLE_S,
     &              '  DISTANCE: ',   IRADRIDGEMIN
        WRITE(6,90) ' LONG  RIDGE Y:',IRAY_LONG,
     &              '   ANG:',        ANGLE_L,
     &              '  DISTANCE: ',   IRADRIDGEMAX
90      FORMAT(1X, A,I6, A,F7.2, A,I5) 

        IASTIG = IRADRIDGEMAX - IRADRIDGEMIN 

        WRITE(6,91) ' ASTIGMATISM:', IASTIG
91      FORMAT(1X, A,I5, A,ES10.3, A,I5) 

        !CALL GETFILENUM(FILNAM,KEY,NDIGITS,.TRUE.,IRTFLG)
        KEY = 9999

C       SAVE PARAMETERS
        DLIST(1) = IRAY_SHORT 
        DLIST(2) = ANGLE_S
        DLIST(3) = IRADRIDGEMIN

        DLIST(4) = IRAY_LONG
        DLIST(5) = ANGLE_L
        DLIST(6) = IRADRIDGEMAX

        DLIST(7) = IASTIG
        CALL REG_SET_NSELA(7,DLIST,.TRUE.,IRTFLG)

C        123456789 123456789 123456789 123456789 123456789 123456789 
        FORMOUT  = 
     &  '(I7,1X,I2,1X,F5.0,3X,F6.2,1X,F8.3,1X,F5.0,3X,F6.2,1X,' //
     &  'F8.3,F5.0,F5.0)'
C           123456789 123456789 123456789 123456789 123456789 123456789 1234
        CALL LUNDOCWRTDATF(LUNDOCNO,KEY,DLIST,8,FORMOUT,IRTFLG)

9999    IF (ALLOCATED(BUF))     DEALLOCATE (BUF)
        IF (ALLOCATED(CURVE))   DEALLOCATE (CURVE)

        CLOSE(LUNIN)
        CLOSE(LUNDOCN) 
 
        END



C       --------------------- FINDMINMAX ------------------------
C       FIND VALLEYS & RIDGES LOCATIONS ALONG SPECIFIED CURVE


        SUBROUTINE FINDMINMAX(CURVE,N,
     &             LOCVALLEY,LOCRIDGE, VALVALLEY,VALRIDGE,
     &             NUMMINS,NUMMAXS, NMIN,NMAX,IRTFLG)

        INTEGER             :: N
        REAL                :: CURVE(N)
        INTEGER             :: NUMMINS,NUMMAXS

        INTEGER             :: LOCVALLEY(NUMMINS),LOCRIDGE(NUMMAXS)
        REAL                :: VALVALLEY(NUMMINS),VALRIDGE(NUMMAXS)
        INTEGER             :: NMIN,NMAX
        INTEGER             :: IRTFLG

        REAL                :: VAL,VTMP
        LOGICAL             :: WANTMIN,WANTMAX

        IRTFLG  = 0
        VAL     = CURVE(1)         ! VALUE AT ORIGIN

        NMIN    = 0
        NMAX    = 0

        WANTMIN = .TRUE.
        WANTMAX = .TRUE.

         DO ILOC = 2,N          ! LOOP ALONG THE CURVE
           VTMP = CURVE(ILOC)

           IF (VTMP < VAL) THEN
C             GOING DOWN
              IF (WANTMAX) THEN
C                FOUND A MAX (MAY BE AT INITIAL POINT)

                 IF (NMAX >= NUMMAXS) CYCLE  
                 NMAX         = NMAX + 1
   
                 LOCRIDGE(NMAX) = ILOC - 1  ! FOUND A RIDGE
                 VALRIDGE(NMAX) = VAL       ! VALUE AT RIDGE
                 WANTMAX      = .FALSE.
                 WANTMIN      = .TRUE.
              ENDIF
           ELSE
C             GOING UP
              IF (WANTMIN) THEN
C                FOUND A MIN (MAY BE AT INITIAL POINT)
                 IF (NMIN >= NUMMINS) CYCLE
                 NMIN         = NMIN + 1
                 LOCVALLEY(NMIN) = ILOC - 1  ! FOUND A VALLEY
                 VALVALLEY(NMIN) = VAL    ! VALUE AT VALLEY
                 WANTMIN      = .FALSE.
                 WANTMAX      = .TRUE.
              ENDIF
           ENDIF

           IF (NMIN >= NUMMINS .AND. NMAX >= NUMMAXS) RETURN
                 
           VAL = VTMP  

           !write(6,'(a,i5,a,f7.2,a,L)')' R:',irad,' ',raymin,' ',gotmin
        ENDDO

        END









