C **********************************************************************
C
C  TO_RAYS     NEW                                 AUG 13 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2016  Health Research Inc.,                         *
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









