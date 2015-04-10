C **********************************************************************
C
C  FINDMINSORMAXS     NEW                         NOV 14 ARDEAN LEITH
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
C FINDMINSORMAXS
C
C PURPOSE:  FIND RIDGES / VALLEYS ACROSS AN IMAGE. PUT RESULTS IN DOC
C           FILE.
C
C **********************************************************************

        SUBROUTINE FINDMINSORMAXS(LUNIN,LUNDOC)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        INTEGER               :: LUNIN,LUNDOC 

C       MAXNAM, INUMBR & NIMAX IS FROM: CMLIMIT.INC
C       MAXNAM, INUMBR & NIMAX IS FROM: CMLIMIT.INC
 
        CHARACTER(LEN=MAXNAM) :: FILNAM,DOCNAM
        CHARACTER(LEN=92)     :: COMMENT
        CHARACTER(LEN=3)      :: OPTC

        REAL, ALLOCATABLE     :: BUF(:,:)
        REAL, ALLOCATABLE     :: CURVE(:)
        REAL, ALLOCATABLE     :: DLIST(:)
        REAL, ALLOCATABLE     :: SIZVAL(:)
        INTEGER, ALLOCATABLE  :: LOCVAL(:)

        INTEGER               :: MAXIM,ITYPE,NX,NY,NZ,IRTFLG,NLET
        INTEGER               :: IROW
        INTEGER               :: NMIN,NMAX,NEIB,NILMAX

        LOGICAL               :: WANTMINS
        LOGICAL               :: ADDEXT,GETNAME,ISOLD
        LOGICAL               :: APPEND,MESSAGE,NEWFILE
        INTEGER               :: LUNDOCO,NLETD,KEY,J
        INTEGER               :: NCHAR,NVALS,NOT_USED,NVAL,NGOT,I
 
        CHARACTER(LEN=1)      :: NULL = CHAR(0)

C       OPEN INPUT FILE
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNIN,'O',ITYPE,NX,NY,
     &              NZ,MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (NZ > 1)  THEN
           CALL ERRT(101,'DOES NOT WORK ON VOLUMES',IRTFLG)
           RETURN
        ENDIF

        ALLOCATE (BUF(NX,NY),
     &            CURVE(NX), 
     &            STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'FINDMINSORMAXS; BUF, CURVE',NX*NY + NX)
           GOTO 9999
        ENDIF

C       LOAD INPUT IMAGE
        CALL REDVOL(LUNIN,NX,NY, 1,1, BUF,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
           
        CALL RDPRMC(OPTC,NCHAR,.TRUE.,'VALLEYS OR RIDGES? (VAL/RID)',
     &              NULL,IRTFLG)
        WANTMINS = (INDEX(OPTC(:NCHAR),'R') == 0)

        NVALS  = 1        
        IF (WANTMINS) THEN
           CALL RDPRI1S(NVALS,NOT_USED,'VALLEYS WANTED',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           IF (NVALS <= 0 ) THEN
              CALL ERRT(102,'INVALID NUMBER OF VALLEYS',NVALS)
              GOTO 9999
           ENDIF
        ELSE
           CALL RDPRI1S(NVALS,NOT_USED,'RIDGES WANTED',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           IF (NVALS <= 0 ) THEN
              CALL ERRT(102,'INVALID NUMBER OF RIDGES',NVALS)
              GOTO 9999
           ENDIF
        ENDIF

        ALLOCATE (SIZVAL(NVALS),
     &            LOCVAL(NVALS), 
     &            DLIST(NVALS*2),   STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'FINDMINSORMAXS; SIZVAL...',NVALS*4)
           GOTO 9999
        ENDIF

        CALL RDPRI1S(NEIB,NOT_USED,
     &              'SEARCH NEIGHBORHOOD DISTANCE',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        ADDEXT  = .TRUE.
        GETNAME = .TRUE.
        ISOLD   = .FALSE.
        APPEND  = .FALSE.
        MESSAGE = .TRUE. 
        IRTFLG  = -8         ! NO IC USE

        CALL OPENDOC(DOCNAM,ADDEXT,NLETD,LUNDOC,LUNDOCO,GETNAME,
     &         'OUTPUT DOC',ISOLD,APPEND,MESSAGE,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        COMMENT = 
     &     'KEY=ROW, LOC1,        VALUE1,        LOC2,        VALUE2,'//
     &      '        LOC3,        VALUE3 ...'
C           123456789 123456789 123456789 123456789 123456789 123456789

        CALL LUNDOCPUTCOM(LUNDOCO,COMMENT(1:92),IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       FIND RIDGES & VALLEY LOCATIONS ALONG ALL RAYS ------------- 

        KEY = 0
       
        DO IROW=1,NY     ! LOOP OVER ALL ROWS 

           CURVE = BUF(1:NX,IROW)

           IF (WANTMINS) THEN
              CALL FINDMINS(CURVE,NX, LOCVAL,SIZVAL,
     &                      NVALS, NGOT, NEIB,IRTFLG)
           ELSE
              CALL FINDMAXS(CURVE,NX, LOCVAL,SIZVAL,
     &                      NVALS, NGOT, NEIB,IRTFLG)
           ENDIF

           KEY = IROW 
           J   = 0
          
           DO I = 1,NVALS     ! LOOP OVER ALL VALLEYS/RIDGES WANTED 
 
              DLIST(J+1) = LOCVAL(I)
              DLIST(J+2) = SIZVAL(I)
              J          = J + 2
           ENDDO

           CALL LUNDOCWRTDAT(LUNDOCO,KEY,DLIST,J,IRTFLG)
         ENDDO

9999    IF (ALLOCATED(BUF))     DEALLOCATE (BUF)
        IF (ALLOCATED(CURVE))   DEALLOCATE (CURVE)
        IF (ALLOCATED(SIZVAL))  DEALLOCATE (SIZVAL)
        IF (ALLOCATED(LOCVAL))  DEALLOCATE (LOCVAL)


        CLOSE(LUNIN)
        CLOSE(LUNDOC) 
 
        END

C       --------------------- FINDMINS ------------------------

C       FIND VALLEYS & RIDGES LOCATIONS ALONG SPECIFIED CURVE

        SUBROUTINE FINDMINS(CURVE,N, LOCVALLEY, VALVALLEY,
     &                      NUMMIN, MINGOT,  NEIB, IRTFLG)

        IMPLICIT NONE
        INTEGER             :: N
        REAL                :: CURVE(N)

        INTEGER             :: NUMMIN,MINGOT
        INTEGER             :: LOCVALLEY(NUMMIN)
        REAL                :: VALVALLEY(NUMMIN)
        INTEGER             :: NEIB
        INTEGER             :: IRTFLG

        REAL                :: VAL,VTMP
        LOGICAL             :: LOCMIN,WANTMIN,WANTMAX
        INTEGER             :: I,J

        IRTFLG  = 0
        VAL     = CURVE(1)    ! VALUE AT ORIGIN


        LOCVALLEY = 0
        VALVALLEY = 0.0
        MINGOT    = 0

        WANTMIN   = .TRUE.
        WANTMAX   = .TRUE.

        DO I = 2,N          ! LOOP ALONG THE CURVE
           VTMP = CURVE(I)

           !if (i<10) write(6,*) 'i,vt,v:', i,vtmp, val,wantmin,wantmax

           IF (VTMP < VAL) THEN
C             GOING DOWN NOW
              IF (WANTMAX) THEN
C                FOUND A MAX (MAY BE AT INITIAL POINT)

                 WANTMAX        = .FALSE.
                 WANTMIN        = .TRUE.
              ENDIF
           !if (i<10) write(6,*) 'down,vt,v:',vtmp,val,wantmin,wantmax
           ELSE
C             GOING UP NOW
              IF (WANTMIN) THEN
C                FOUND A MIN (MAY BE AT INITIAL POINT)

                 LOCMIN = .TRUE.
                 IF (NEIB > 0) THEN
C                   
C                   CHECK IF A LOCAL MIN OVER CURRENT NEIGHBORHOOD
                    DO J=I,I-NEIB,-1
                       IF (J > 0  .AND. CURVE(J) < VAL) LOCMIN = .FALSE.
                    ENDDO

                    DO J=I,I+NEIB
                       IF (J <= N .AND. CURVE(J) < VAL) LOCMIN = .FALSE.
                    ENDDO
                 ENDIF
 
                 IF (LOCMIN) THEN
                    MINGOT            = MINGOT + 1
                    LOCVALLEY(MINGOT) = I - 1  ! FOUND A VALLEY
                    VALVALLEY(MINGOT) = VAL    ! VALUE AT VALLEY
                 ENDIF

                 WANTMIN      = .FALSE.
                 WANTMAX      = .TRUE.
              ENDIF
         !if(i<10)write(6,*)'up  ,vt,v:',vtmp,val,wantmin,wantmax,locmin
           ENDIF

           VAL = VTMP
  
           IF (MINGOT >= NUMMIN) EXIT
          
        ENDDO

        END

C       --------------------- FINDMAXS ------------------------

C       FIND RIDGES LOCATIONS ALONG SPECIFIED CURVE

        SUBROUTINE FINDMAXS(CURVE,N, LOCRIDGE, VALRIDGE,
     &                      NUMMAX,MAXGOT, NEIB, IRTFLG)

        IMPLICIT NONE

        REAL                :: CURVE(N)
        INTEGER             :: N
        INTEGER             :: LOCRIDGE(NUMMAX)
        REAL                :: VALRIDGE(NUMMAX)
        INTEGER             :: NUMMAX,MAXGOT
        INTEGER             :: NEIB
        INTEGER             :: IRTFLG

        REAL                :: VAL,VTMP
        LOGICAL             :: WANTMIN,WANTMAX,LOCMAX
        INTEGER             :: I,J

        IRTFLG  = 0
        VAL     = CURVE(1)    ! VALUE AT ORIGIN

        LOCRIDGE = 0
        VALRIDGE = 0.0
        MAXGOT   = 0

        WANTMIN  = .TRUE.
        WANTMAX  = .TRUE.

         DO I = 2,N          ! LOOP ALONG THE CURVE
           VTMP = CURVE(I)

           IF (VTMP < VAL) THEN
C             GOING DOWN
              IF (WANTMAX) THEN
C                FOUND A MAX (MAY BE AT INITIAL POINT)

                 LOCMAX = .TRUE.
                 IF (NEIB > 0) THEN
C                   
C                   CHECK IF A LOCAL MAX OVER CURRENT NEIGHBORHOOD
                    DO J=I,I-NEIB,-1
                       IF (J > 0  .AND. CURVE(J) > VAL) LOCMAX = .FALSE.
                    ENDDO

                    DO J=I,I+NEIB
                       IF (J <= N .AND. CURVE(J) > VAL) LOCMAX = .FALSE.
                    ENDDO
                 ENDIF
 
                 IF (LOCMAX) THEN
                    MAXGOT           = MAXGOT + 1
                    LOCRIDGE(MAXGOT) = I - 1  ! FOUND A RIDGE
                    VALRIDGE(MAXGOT) = VAL    ! VALUE AT RIDGE
                 ENDIF

                 WANTMAX        = .FALSE.
                 WANTMIN        = .TRUE.
              ENDIF
           ELSE
C             GOING UP
              IF (WANTMIN) THEN
C                FOUND A MIN (MAY BE AT INITIAL POINT)

                 WANTMIN      = .FALSE.
                 WANTMAX      = .TRUE.
              ENDIF
           ENDIF

           VAL = VTMP  

           IF (MAXGOT >= NUMMAX) EXIT
                 
        ENDDO

        END









