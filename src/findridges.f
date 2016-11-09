C **********************************************************************
C
C  FINDRIDGES     NEW                              AUG 13 ArDean Leith
C                 ADDED RIDGE SEPERATION           OCT 16 ArDean Leith
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2016 Health Research Inc.,                          *
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
C FINDRIDGES(RIDGESONLY)
C
C PURPOSE:   FIND RIDGES AND VALLEYS RUNNING UP/DOWN AN IMAGE.
C            USUALLY USED WITH POWER SPECTRA CONVERTED TO POLAR VIEW
C            WITH 'PO R'
C
C PARAMETERS:  RIDGESONLY    FINDS FIRST 2 RIDGES ACROSS IMAGE FOR
C                            RAY HAVING GREATEST FIRST RIDGE DISTANCE
C
C OPERATIONS: 'RI R'  & 'RI RV'
C
C **********************************************************************

        SUBROUTINE FINDRIDGES(RIDGESONLY)

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        LOGICAL               :: RIDGESONLY

        CHARACTER(LEN=MAXNAM) :: FILNAM,DOCNAM
        CHARACTER(LEN=84)     :: FORMOUT
        CHARACTER(LEN=96)     :: COMMENT

        REAL, ALLOCATABLE     :: BUF(:,:)
        REAL, ALLOCATABLE     :: CURVE(:)

        REAL                  :: DLIST(9)
        INTEGER               :: IRAY_SHORT,IRAY_LONG
        INTEGER               :: MAXIM,ITYPE,NX,NY,NZ,IRTFLG,IER,NLET
        INTEGER               :: NA,NR,IROW,ISEL,NSUM,N1,N2,J,I,MIND
        REAL                  :: ANGLE_L,ANGLE_S
        INTEGER               :: IRADRIDGESHORT,IRADRIDGELONG
        INTEGER               :: IASTIG,NMIN,NMAX,NSEPER,ISEP
        INTEGER               :: IRADRIDGEMAX,IRADRIDGEMIN,IRADRIDGE
        INTEGER               :: IRADRIDGEMAX2,IRADRIDGE2,NOT_USED
        INTEGER               :: IRADRIDGEMAX3,IRADRIDGE3
        REAL                  :: VALRIDG 

        LOGICAL               :: ADDEXT,GETNAME,ISOLD
        LOGICAL               :: APPEND,MESSAGE,NEWFILE
        INTEGER               :: LUNDOCNO,NLETD,KEY,NDIGITS
 
        INTEGER, PARAMETER    :: NUMMINS = 5 ! MAX NUMBER OF VALLEYS
        INTEGER, PARAMETER    :: NUMMAXS = 5 ! MAX NUMBER OF RIDGES

        INTEGER               :: LOCVALLEY(NUMMINS),LOCRIDGE(NUMMAXS)
        REAL                  :: VALVALLEY(NUMMINS),VALRIDGE(NUMMAXS)

        INTEGER, PARAMETER    :: LUNIN   = 20
        INTEGER, PARAMETER    :: LUNDOCN = 80

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


        NSEPER = 1
        NSUM   = 1
        CALL RDPRI2S(NSEPER,NSUM,NOT_USED,
     &             'SEPARATION, NUMBER OF SUMMED LINES',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        MIND   = 1
        KEY    = 1
        CALL RDPRI2S(MIND,KEY,NOT_USED,
     &             'MINIMUM DISTANCE, DOC FILE KEY',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999


        IF (RIDGESONLY) THEN
C                    123456789 123456789 123456789 123456789 123456789 123456789
         COMMENT= 
     &   '  NUM,     Y,     ANG,   R1,    R2,    R3,    ASTIG   SEP'

        ELSE
C                   123456789 123456789 123456789 123456789 123456789 123456789
          COMMENT= '  NUM,    YS,       ANG_S,        R_S,        '//
     &             'YL,       ANG_L,        R_L'
        ENDIF


        IF (KEY < 2) THEN
           !write(6,*) ' lundocno, comment:',lundocno,comment(1:40)
           CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:86),IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
        ENDIF

C       FIND RIDGES & VALLEY LOCATIONS ALONG ALL RAYS ------------- 

        IRAY_LONG    = 0
        IRAY_SHORT   = 0
        IRADRIDGEMAX = 0
        IRADRIDGEMIN = HUGE(IRADRIDGEMIN)

        DO IROW=1+2,NY-2  ! LOOP OVER ALL ROWS (SKIP 2 AT TOP & BOTTEM)

           IF (NSUM < 2) THEN
              CURVE = BUF(1:NX,IROW)
           ELSE
              N1 = IROW - NSUM/2
              IF (N1 < 3) THEN
                 N1 = 3
              ENDIF

              N2 = N1 + NSUM - 1
              IF (N2 > (NY - 2)) THEN
                 N2 = NY - 2
                 N1 = N2 - NSUM
              ENDIF
              CURVE = BUF(:, N1) 
              DO J = N1+1,N2
                 CURVE = CURVE + BUF(:, J)
              ENDDO
           ENDIF

           CALL FINDMINMAX(CURVE,NX,
     &             LOCVALLEY,LOCRIDGE, VALVALLEY,VALRIDGE,
     &             NUMMINS,NUMMAXS, NSEPER, NMIN,NMAX,IRTFLG)

           DO I = 1,NUMMAXS -1
C             FIND LONGEST AXIS OF CIRCLE RAY ---------------------------
              IRADRIDGE  = LOCRIDGE(I)   ! RADIUS AT FIRST  RIDGE
              IRADRIDGE2 = LOCRIDGE(I+1) ! RADIUS AT SECOND RIDGE
              IRADRIDGE3 = LOCRIDGE(I+2) ! RADIUS AT THIRD  RIDGE
              VALRIDG    = VALRIDGE(I)   ! VALUE  AT FIRST  RIDGE 

              IF ( IRADRIDGE >= MIND ) EXIT
           ENDDO

C          FIND LONGEST DISTANCE TO RIDGE ---------------------------
           IF (IRADRIDGE > IRADRIDGEMAX) THEN
              IRAY_LONG     = IROW
              IRADRIDGEMAX  = IRADRIDGE
              IRADRIDGEMAX2 = IRADRIDGE2
              IRADRIDGEMAX3 = IRADRIDGE3
           ENDIF
   
C          FIND SHORTEST DISTANCE TO RIDGE ---------------------------
           IF (IRADRIDGE < IRADRIDGEMIN) THEN
              IRAY_SHORT   = IROW
              IRADRIDGEMIN = IRADRIDGE
           ENDIF
        ENDDO

        ANGLE_S = 180.0 * FLOAT(IRAY_SHORT) / FLOAT(NY)
        ANGLE_L = 180.0 * FLOAT(IRAY_LONG)  / FLOAT(NY)
        IASTIG  = IRADRIDGEMAX  - IRADRIDGEMIN 
        ISEP    = IRADRIDGEMAX2 - IRADRIDGEMAX

        IF (RIDGESONLY) THEN

           WRITE(6,92) ' LONG  RIDGE Y:', IRAY_LONG,
     &                 '  ANG:',          ANGLE_L,
     &                 '  R1: ',          IRADRIDGEMAX,
     &                 '  R2: ',          IRADRIDGEMAX2,
     &                 '  R3: ',          IRADRIDGEMAX3,
     &                 '  ASTIG: ',       IASTIG,
     &                 '  SEP: ',         ISEP
92         FORMAT(1X, A,I6, A,F7.2, A,I5, A,I5, A,I5, A,I5, A,I5) 

        ELSE
           WRITE(6,90) ' SHORT RIDGE Y:',IRAY_SHORT,
     &                 '   ANG:',        ANGLE_S,
     &                 '  DISTANCE: ',   IRADRIDGEMIN

           WRITE(6,90) ' LONG  RIDGE Y:',IRAY_LONG,
     &                 '   ANG:',        ANGLE_L,
     &                 '  DISTANCE: ',   IRADRIDGEMAX
90         FORMAT(1X, A,I6, A,F7.2, A,I5) 

           WRITE(6,91) ' ASTIGMATISM:', IASTIG
91         FORMAT(1X, A,I5) 
        ENDIF


        !CALL GETFILENUM(FILNAM,KEY,NDIGITS,.TRUE.,IRTFLG)

C       SAVE PARAMETERS
        DLIST(1) = IRAY_SHORT 
        DLIST(2) = ANGLE_S
        DLIST(3) = IRADRIDGEMIN
        IF (RIDGESONLY) THEN
           DLIST(1) = IRAY_LONG 
           DLIST(2) = ANGLE_L
           DLIST(3) = IRADRIDGEMAX
           DLIST(4) = IRADRIDGEMAX2
           DLIST(5) = IRADRIDGEMAX3

           DLIST(6) = IASTIG
           DLIST(7) = ISEP
           ISEL     = 7

C            123456789 123456789 123456789 123456789 123456789 123456789 
           FORMOUT  = 
     &     '(I7,1X,I2,1X,F5.0,2X,F6.2,3(1X,F6.0),3X,F6.0,1X,F6.0,)'
         ELSE
           DLIST(1) = IRAY_SHORT 
           DLIST(2) = ANGLE_S
           DLIST(3) = IRADRIDGEMIN
           DLIST(4) = IRAY_LONG
           DLIST(5) = ANGLE_L
           DLIST(6) = IRADRIDGEMAX

           DLIST(7) = IASTIG
           ISEL     = 7

C            123456789 123456789 123456789 123456789 123456789 123456789 
           FORMOUT  = 
     &     '(I7,1X,I2,1X,F5.0,3X,F6.2,1X,F8.3,1X,F5.0,3X,F6.2,1X,' //
     &     'F8.3,F5.0,F5.0)'
        ENDIF

        CALL REG_SET_NSELA(ISEL,DLIST,.TRUE.,IRTFLG)

        CALL LUNDOCWRTDATF(LUNDOCNO,KEY,DLIST,ISEL,FORMOUT,IRTFLG)
      !write(6,*) ' key,dlist,isel:',key,dlist,isel

9999    IF (ALLOCATED(BUF))     DEALLOCATE (BUF)
        IF (ALLOCATED(CURVE))   DEALLOCATE (CURVE)

        CLOSE(LUNIN)
        CLOSE(LUNDOCN) 
 
        END



C       --------------------- FINDMINMAX ------------------------
C       FIND VALLEYS & RIDGES LOCATIONS ALONG SPECIFIED CURVE


        SUBROUTINE FINDMINMAX(CURVE,N,
     &             LOCVALLEY,LOCRIDGE, VALVALLEY,VALRIDGE,
     &             NUMMINS,NUMMAXS,NSEPER, NMIN,NMAX,IRTFLG)

        IMPLICIT NONE
        INTEGER             :: N
        REAL                :: CURVE(N)
        INTEGER             :: NUMMINS,NUMMAXS,NSEPER

        INTEGER             :: LOCVALLEY(NUMMINS),LOCRIDGE(NUMMAXS)
        REAL                :: VALVALLEY(NUMMINS),VALRIDGE(NUMMAXS)
        INTEGER             :: NMIN,NMAX
        INTEGER             :: IRTFLG

        INTEGER             :: ILOC,ISEPER,N1,N2
        REAL                :: VAL,VTMP,CMAX,CMIN
        LOGICAL             :: WANTMIN,WANTMAX
        INTEGER             :: MAXL, MAXL_ARRAY(1)

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
 
                 IF (NSEPER > 1 .AND. NMAX > 0 .AND. ILOC < N) THEN
                    N1   = ILOC + 1
                    N2   = MIN(N, (ILOC + NSEPER))
                    CMAX = MAXVAL(CURVE(N1:N2)) 
                    IF (CMAX > VTMP) THEN
                       !VAL = VTMP  
                       CYCLE
                    ENDIF
                 ENDIF
 
                 NMAX           = NMAX + 1
                 LOCRIDGE(NMAX) = ILOC - 1  ! FOUND A RIDGE
                 VALRIDGE(NMAX) = VAL       ! VALUE AT RIDGE
                 WANTMAX        = .FALSE.
                 WANTMIN        = .TRUE.
              ENDIF
           ELSE
C             GOING UP
              IF (WANTMIN) THEN
C                FOUND A MIN (MAY BE AT INITIAL POINT)
                 IF (NMIN >= NUMMINS) CYCLE

                 IF (NSEPER > 1 .AND. NMIN > 0 .AND. ILOC < N) THEN
                    N1   = ILOC + 1
                    N2   = MIN(N, (ILOC + NSEPER))
                    CMIN = MINVAL(CURVE(N1:N2)) 
                    IF (CMIN < VTMP) THEN
                       !VAL = VTMP  
                       CYCLE
                    ENDIF
                 ENDIF

                 NMIN            = NMIN + 1
                 LOCVALLEY(NMIN) = ILOC - 1  ! FOUND A VALLEY
                 VALVALLEY(NMIN) = VAL    ! VALUE AT VALLEY
                 WANTMIN         = .FALSE.
                 WANTMAX         = .TRUE.
              ENDIF
           ENDIF

           IF (NMIN >= NUMMINS .AND. NMAX >= NUMMAXS) RETURN
                 
           VAL = VTMP  

           !write(6,'(a,i5,a,f7.2,a,L)')' R:',irad,' ',raymin,' ',gotmin
        ENDDO

        END









