
C++*********************************************************************
C
C  UTIL_1010.F  FROM UTIL_1110                   SEP 2012 ARDEAN LEITH 

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
C   UTIL_1010()
C
C   PURPOSE:  NORMALIZE IMAGES/VOLUMES BY FINDING AVG AND STD.
C             DEVIATION OUTSIDE MASK (OR OVERALL)
C             SUBTRACTS AVERAGE AND DIVIDES BY STD. DEV.
C             CAN OPERATE ON IMAGE SERIES
C             CAN USE A SELECTION DOC FILE FOR IMAGE NUMBERS
C
C--*********************************************************************

      SUBROUTINE UTIL_1010()

      IMPLICIT NONE
      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

C     DOC FILE POINTER
      INCLUDE 'F90ALLOC.INC'
      REAL, POINTER          :: PBUF(:,:)

      REAL,    ALLOCATABLE   :: BUFIMG(:),BUFIN(:),PLIST(:,:)
      LOGICAL, ALLOCATABLE   :: BUFMSK(:,:,:)
      INTEGER, ALLOCATABLE   :: ILIST1(:),ILIST2(:)

      CHARACTER (LEN=MAXNAM) :: FILPATIN,FILPATOUT,FILNAM,FILMASK
      CHARACTER (LEN=MAXNAM) :: DOCNAM,DOCNAME,MSG

      CHARACTER (LEN=1)      :: NULL = CHAR(0)
      CHARACTER (LEN=1)      :: ANS,CDUM
      CHARACTER (LEN=83)     :: PROMPT

      INTEGER                :: NILMAX,ITYPE,NX,NY,NZ
      INTEGER                :: NUMB,MAXX,MAXY,NLET,NDUM,NIMGOUT
      INTEGER                :: NIMG,IMGNUM,IMGNUMOUT,MWANT
      INTEGER                :: NINDX1,NINDX2,NC
      INTEGER                :: MAXIMIN,MAXIMOUT
      INTEGER                :: LOCAT,LOCAST,NIMGT,NLETF
      INTEGER                :: IX,IY,IZ,IRECIN,IMAMIT
      INTEGER                :: MAXIM,IRTFLG
      INTEGER                :: LX,LY,LZ,NI,NO,IRADI
      INTEGER                :: PERCENTIN,PERCENTOUT,NPIX
      DOUBLE PRECISION       :: DAVI,DSIGI,DAVO,DSIGO
      LOGICAL                :: USE_OMP,OVERALL,RAMPIT
      REAL                   :: UNUSED, dgm_1, dgm_2, dgm_3, dgm_4

      INTEGER,PARAMETER      :: LUNIN     = 21 
      INTEGER,PARAMETER      :: LUNMSK    = 22
      INTEGER,PARAMETER      :: LUNOUT    = 23
      INTEGER,PARAMETER      :: LUNDOCSEL = 81
      INTEGER,PARAMETER      :: LUNXM1    = 82
      INTEGER,PARAMETER      :: LUNXM2    = 83
      INTEGER,PARAMETER      :: LUNDOC    = 84

      INTEGER                :: lnblnkn
 
      NILMAX  = NIMAXPLUS      ! FROM CMLIMIT
      ALLOCATE(ILIST1(NILMAX),
     &         ILIST2(NILMAX),
     &         STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
          CALL ERRT(46,'UTIL_1010; ILIST1....',2*NILMAX)
          RETURN
      ENDIF

C     OPEN INPUT IMAGE(S) (NOT FOURIER)
      CALL OPFILES(0,LUNIN,LUNDOCSEL,LUNXM1, 
     &             .TRUE.,FILPATIN,NLET, 'O',
     &             ITYPE,NX,NY,NZ,MAXIMIN,
     &             NULL,
     &             .FALSE.,ILIST1,NILMAX, 
     &             NDUM,NIMG,IMGNUM, IRTFLG) 
      IF (IRTFLG .NE. 0) RETURN

      LOCAT  = INDEX(FILPATIN,'@')
      LOCAST = INDEX(FILPATIN,'*')
      IF (NIMG > 0) ILIST2 = ILIST1
  
      IMAMIT  = IMAMI
      IF (IMAMI > 0) THEN
         DSIGO   = SIG
         DAVO    = AV
      ENDIF
      OVERALL = .FALSE.

C     SINGLE IMAGE OPERATION 
      IF (NIMG <= 1) IMGNUM = 0
 
C     OPEN MASK INPUT IMAGE (USED FOR ALL INPUT IMAGES) 
      MAXIM = 0
      CALL OPFILEC(0,.TRUE.,FILMASK,LUNMSK,'O',ITYPE,
     &             LX,LY,LZ,
     &             MAXIM,'MASK (* IF NO MASK)',.FALSE.,IRTFLG)

      IF (IRTFLG .EQ. -1 .AND. FILMASK(1:1) == '*') THEN
         OVERALL = .TRUE.
      ELSE
         IF (IRTFLG .NE. 0) GOTO 9999

         CALL SIZCHK(UNUSED, LX,LY,LZ,0,  NX,NY,NZ,0, IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999

         MWANT = NX*NY*NZ + NX 
         ALLOCATE(BUFMSK(NX,NY,NZ),
     &            BUFIN(NX), STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN 
            CALL ERRT(46,'UTIL_1010, BUFMSK...',MWANT)
            GOTO 9999
         ENDIF

C        READ MASK IMAGE.  EXTRACT LOGICAL MASK IMAGE
         DO IZ=1,NZ
            DO IY=1,NY
               IRECIN = (IZ -1) * NY + IY
               CALL REDLIN(LUNMSK,BUFIN, NX, IRECIN)
               DO IX = 1,NX
                  BUFMSK(IY,IX,IZ) = (BUFIN(IX) > 0.5)
                  !if (iy == 32)write(6,*) 'mask:',ix,iy,iz,bufmsk(IY,IX,IZ)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

C     OPEN OUTPUT IMAGE(S)
      PROMPT    = 'OUTPUT FILE NAME OR TEMPLATE (E.G. IMG@****)~'
      IF (LOCAST == 0)  PROMPT = 'OUTPUT'
      IF (LOCAST == 0 .AND. LOCAT > 0)  PROMPT = 'OUTPUT STACK'
 
      MAXIMOUT  = -1           ! ALLOW BARE STACK
      IMGNUMOUT = IMGNUM       ! IMAGE # WANTED
      NIMGT     = -NIMG        ! USE FIRST LIST FOR OUTPUT

      CALL OPFILES(LUNIN,LUNOUT,LUNDOCSEL,LUNXM2,
     &            .TRUE.,FILPATOUT,NLET, 'U',
     &            ITYPE,NX,NY,NZ,MAXIMOUT,
     &            PROMPT,
     &            .FALSE., ILIST2,NIMGT, 
     &            NDUM,NIMGOUT, IMGNUMOUT, IRTFLG) 

      IF (IRTFLG .NE. 0) GOTO 9999
      !write(6,*)'output Image #:',IMGNUMOUT,NIMG,maximout,filpatout(:11)
      !write(6,*) 'fchar:',fchar(1:10)

      RAMPIT = .FALSE.
      IF (NZ == 1) THEN
         CALL RDPRMC(ANS,NC,.TRUE.,
     &            'APPLY FLAT-FIELD RAMP (Y,N)',CDUM,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999
         RAMPIT  = (ANS(1:1) == 'Y')

         IF (RAMPIT) THEN
            MWANT = NX*NY 
            ALLOCATE(PLIST(3,MWANT),STAT=IRTFLG)
            IF (IRTFLG .NE. 0) THEN 
               CALL ERRT(46,'UTIL_1010, PLIST',MWANT*3)
               GOTO 9999
            ENDIF
         ENDIF
      ENDIF

      MWANT = NX*NY*NZ 
      ALLOCATE(BUFIMG(MWANT),STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN 
         CALL ERRT(46,'UTIL_1010, BUFIMG',MWANT)
         GOTO 9999
      ENDIF

      NINDX1  = 1
      NINDX2  = 1
      USE_OMP = .TRUE.
      NPIX    = NX*NY*NZ
      NO      = NPIX

      DO 
         !write(6,*)  'Image #:',IMGNUM

C        READ INPUT IMAGE
         CALL REDVOL(LUNIN,NX,NY,1,NZ,BUFIMG,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999 

         IF (RAMPIT .AND. OVERALL) THEN
C           RAMP WITH ALL IMAGE VALUES
            CALL RAMP_PB(BUFIMG,NX,NY,.TRUE.,NOUT)
            
         ELSEIF (RAMPIT) THEN
C          RAMP WITH BACKGROUND IMAGE VALUES ONLY
           CALL RAMP_PB_MASK(BUFIMG,BUFMSK,PLIST,
     &                       NX,NY,.TRUE.,NOUT)
         ENDIF

         IF (OVERALL) THEN
            IF (IMAMIT == 0) THEN
C              FIND OVERALL STATS
               CALL NORMVALS(BUFIMG, NX, NY, NZ,
     &                       DAVO,DSIGO, USE_OMP)
            ENDIF

C           USE OVERALL AVG AND STD WITHOUT ANY MASK
            IF (VERBOSE .AND. IMGNUM > 0) THEN
               WRITE(NOUT,94) IMGNUM,NO,DAVO,DSIGO 
94             FORMAT('  IMAGE:',I6,'  PIXELS: ',I7,
     &                '  AVG:',1PG11.3,'  SIG:',1PG11.3)

            ELSEIF (VERBOSE) THEN
               WRITE(NOUT,95) NO,DAVO,DSIGO 
95             FORMAT('  PIXELS: ',I7,
     &                '  AVG:',1PG11.3,'  SIG:',1PG11.3)
            ENDIF
           
         ELSE
            !call chkfile('jnk111',66,1,nx,ny,1,bufimg,irtflg)

C           FIND AVG AND STD USING MASK
            CALL NORMVALS_LMASKED(BUFIMG, BUFMSK,
     &                            NX,NY,NZ, USE_OMP,
     &                            NI,DAVI,DSIGI, NO,DAVO,DSIGO)

            PERCENTIN  = INT(100.0 * FLOAT(NI) / FLOAT(NPIX))
            PERCENTOUT = INT(100.0 * FLOAT(NO) / FLOAT(NPIX))

            IF (VERBOSE .AND. IMGNUM > 0) THEN
               WRITE(NOUT,90) IMGNUM,NI,PERCENTIN,DAVI,DSIGI 
90             FORMAT('  IMAGE:',I6,
     &                '  Pixels inside: ',I7,'  % Inside:  ',I3,'%',
     &                '  Avg:',1PG11.3,'  Sig:',1PG11.3)

               WRITE(NOUT,91) IMGNUM,NO,PERCENTOUT,DAVO,DSIGO 
91             FORMAT('  IMAGE:',I6,
     &                '  Pixels outside:',I7,'  % Outside: ',I3,'%',
     &                '  AVG:',1PG11.3,'  SIG:',1PG11.3)

            ELSEIF (VERBOSE) THEN
               WRITE(NOUT,92) NI,PERCENTIN,DAVI,DSIGI 
92             FORMAT(
     &                '  Pixels inside: ',I7,'  % Inside:  ',I3,'%',
     &                '  Avg:',1PG11.3,'  Sig:',1PG11.3)
               WRITE(NOUT,93) NO,PERCENTOUT,DAVO,DSIGO 
93             FORMAT(
     &                '  Pixels outside:',I7,'  % Outside: ',I3,'%',
     &                '  Avg:',1PG11.3,'  Sig:',1PG11.3)
            ENDIF
         ENDIF

C        NORMALIZE IMAGE/VOL USING: AVG AND STD 
         CALL ARITH_NORMIT(BUFIMG,NX,NY,NZ,
     &                     DAVO,DSIGO,USE_OMP,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999 

C        SAVE OUTPUT IMAGE/VOL
         CALL WRTVOL(LUNOUT,NX,NY, 1,NZ,BUFIMG,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999 

C dgm mod - dgm_N are REAL while the others are DOUBLE PRECISION
C           and the Intel ifort compiler from 13.0.1 chokes on the
C           FLOAT(DOUBLE PRECISION) call...

         dgm_1 = DAVO
         dgm_2 = DSIGO
         dgm_3 = DAVI
         dgm_4 = DSIGI

C end dgm mod
	 
C        SET OPERTION LINE REGISTERS
         CALL REG_SET_NSEL(1,5, 
     &        FLOAT(NO),dgm_1,dgm_2,
     &                  dgm_3,dgm_4, IRTFLG)
C dgm     &        FLOAT(NO),FLOAT(DAVO),FLOAT(DSIGO),
C dgm     &                  FLOAT(DAVI),FLOAT(DSIGI), IRTFLG)

         IF (NINDX1 >= NIMG) EXIT      ! END OF INPUT LIST
             
C        OPEN NEXT SET OF I/O FILES 
         !write(6,*) NINDX1, NINDX2,NIMG,NIMGOUT,MAXIMIN,MAXIMOUT
         !write(6,*) LUNIN,LUNIN,LUNOUT, FILPATIN,FILPATOUT,IMGNUM,IMGNUMOUT
         CALL NEXTFILES(NINDX1, NINDX2, ILIST1,ILIST2, 
     &                  .FALSE.,LUNDOC,LUNXM2,
     &                  NIMG,NIMGOUT,   
     &                  MAXIMIN,MAXIMOUT,   
     &                  LUNIN,LUNIN,LUNOUT, FILPATIN,FILPATOUT,
     &                  IMGNUM,IMGNUMOUT, IRTFLG) 

          !write(6,'(A,5i6)')
      !&       ' Nextfiles imgnum,imgnumout,nindx1,nindx2,irtflg:',
      !&                   imgnum,imgnumout,nindx1,nindx2,irtflg

         IF (IRTFLG == -99) THEN
             CALL ERRT(102,'INSUFFICIENT OUTPUT FILE NAMES',NINDX2)
             EXIT         
         ELSEIF (IRTFLG < 0) THEN
             EXIT         ! END OF INPUT FILES
         ENDIF
         IF (IRTFLG .NE. 0) GOTO 9999    ! ERROR

         IMAMIT = IMAMI
         IF (IMAMI > 0) THEN
            DSIGO   = SIG
            DAVO    = AV
         ENDIF
      ENDDO

9999  IF (ALLOCATED(BUFIMG))  DEALLOCATE(BUFIMG)
      IF (ALLOCATED(BUFMSK))  DEALLOCATE(BUFMSK)
      IF (ALLOCATED(BUFIN))   DEALLOCATE(BUFIN)
      IF (ALLOCATED(ILIST1))  DEALLOCATE(ILIST1)
      IF (ALLOCATED(ILIST2))  DEALLOCATE(ILIST2)
      IF (ALLOCATED(PLIST))   DEALLOCATE(PLIST)

      IF (VERBOSE) WRITE(NOUT,*) ' '

      CLOSE(LUNIN)
      CLOSE(LUNMSK)
      CLOSE(LUNOUT)
      CLOSE(LUNXM1)
      CLOSE(LUNXM2)
      
      END


C    ************************** ARITH_NORMIT **************************


      SUBROUTINE ARITH_NORMIT(BUF,NX,NY,NZ, DAVG,DSIG,USE_OMP, IRTFLG)

      IMPLICIT NONE

      REAL             :: BUF(NX,NY,NZ)
      INTEGER          :: NX,NY,NZ 
      DOUBLE PRECISION :: DAVG,DSIG 
      LOGICAL          :: USE_OMP
      INTEGER          :: IRTFLG

      INTEGER          :: IX,IY,IZ,NE
      DOUBLE PRECISION :: CON

      IF (DSIG == 0) THEN
         CALL ERRT(101,'CAN NOT NORMALIZE, STD. DEV. = 0',NE)
         IRTFLG = 1
         RETURN
      ENDIF
  
      CON = 1.0 / DSIG          ! FASTER

      !write(6,*)' davg,dsig,nx,con:',davg,dsig,nx,con

      IF (USE_OMP) THEN
c$omp    parallel do private(iz,iy,ix)
         DO IZ=1,NZ
           DO IY = 1,NY
              DO IX = 1,NX
                 BUF(IX,IY,IZ) = (BUF(IX,IY,IZ) - DAVG) * CON
	      ENDDO
           ENDDO
         ENDDO

      ELSE
         DO IZ=1,NZ
           DO IY = 1,NY
              DO IX = 1,NX
                 BUF(IX,IY,IZ) = (BUF(IX,IY,IZ) - DAVG) * CON
	      ENDDO
           ENDDO
         ENDDO
      ENDIF

      IRTFLG = 0

      END










