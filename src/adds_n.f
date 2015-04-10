C **********************************************************************
C
C ADDS_N.F                           USED FILELIST JULY 99 ARDEAN LEITH
C                                    USED F90 JULY 99 ARDEAN LEITH
C                                    SQRT(NEG) TRAP APR 00 ARDEAN LEITH
C                                    ALLOC TRAP     APR 03 ARDEAN LEITH
C
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
C  PURPOSE:  CREATE AVERAGE AND VARIANCE FILES FROM A SET OF INPUT
C            IMAGES/VOLUMES AND/OR EVEN/ODD SUBSETS FROM THE FILES
C
C--*********************************************************************

      SUBROUTINE ADDS_N()

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER (LEN=MAXNAM)        :: FILNAM,FILPAT

      DOUBLE PRECISION              :: VARAV,VV,FNUMT,FNUMT1,FNUMT2
      DOUBLE PRECISION              :: AVALL,DTEMP
      INTEGER                       :: NILMAX,IRTFLG,K,NLET,ITYPE,NCHAR
      INTEGER                       :: MAXIM
      INTEGER                       :: NUMA,NUMA1,NUMA2
      INTEGER                       :: NX,NY,NZ,NDUM,NIMG,IMGNUM
      LOGICAL                       :: WANTA,WANTV,WANT3
      REAL                          :: SAV,SAVT,VAV

      INTEGER, ALLOCATABLE          :: ILIST(:)

      REAL, ALLOCATABLE             :: BUFIN(:)
      DOUBLE PRECISION, ALLOCATABLE :: AVGARAY (:),VARARAY(:) 
      DOUBLE PRECISION, ALLOCATABLE :: AVGARAY1(:),VARARAY1(:) 
      DOUBLE PRECISION, ALLOCATABLE :: AVGARAY2(:),VARARAY2(:) 
      CHARACTER (LEN=3)             :: ANSW

#ifndef SP_32
      INTEGER *8                    :: NVOX,IOK8,MWANT
#else
      INTEGER *4                    :: NVOX,IOK8,MWANT
#endif

      INTEGER,PARAMETER             :: LUNIN     = 21 
      INTEGER,PARAMETER             :: LUNA      = 22
      INTEGER,PARAMETER             :: LUNA1     = 23
      INTEGER,PARAMETER             :: LUNA2     = 24
      INTEGER,PARAMETER             :: LUNV      = 25
      INTEGER,PARAMETER             :: LUNV1     = 26
      INTEGER,PARAMETER             :: LUNV2     = 27

      INTEGER,PARAMETER             :: LUNDOCSEL = 81
      INTEGER,PARAMETER             :: LUNXM     = 0  ! SELFILE NOT ALLOWED
         
      LOGICAL,PARAMETER             :: ASKNAM    = .TRUE.
      LOGICAL,PARAMETER             :: FOUROK    = .TRUE.

      CHARACTER (LEN=1)             :: NULL = CHAR(0)

      NILMAX  = NIMAXPLUS      ! FROM CMLIMIT
      ALLOCATE(ILIST(NILMAX),  STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
          CALL ERRT(46,'ADDS; ILIST',NILMAX)
          RETURN
      ENDIF

C     OPEN INPUT IMAGE(S)
      MAXIM = 0
      CALL OPFILES(0,LUNIN,LUNDOCSEL,LUNXM, 
     &             ASKNAM,FILPAT,NLET, 'O',
     &             ITYPE,NX,NY,NZ,MAXIM,
     &             'TEMPLATE FOR INPUT FILES~',
     &             FOUROK,ILIST,NILMAX, 
     &             NDUM,NIMG,IMGNUM, IRTFLG) 
      IF (IRTFLG .NE. 0) GOTO 9999

      IF (COPT == 'I') WRITE(NOUT,90)
  90  FORMAT(
     &    ' .MENU: A   -- AVERAGE'/
     &    '        AS  -- AVERAGE  WITH SUBSET AVERAGES'/
     &    '        V   -- VARIANCE'/
     &    '        VS  -- VARIANCE WITH SUBSET VARIANCES'/
     &    '        AV  -- AVERAGE  AND VARIANCE'/
     &    '        AVS -- AVERAGE  AND VARIANCE WITH SUBSETS')

      CALL RDPRMC(ANSW,NCHAR,.TRUE.,
     &      'MENU OPTION (A/AS/V/VS/AV/AVS)',NULL,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

      WANTA = ( INDEX(ANSW(:NCHAR),'A') >0)  
      WANTV = ( INDEX(ANSW(:NCHAR),'V') >0) 
      WANT3 = ( INDEX(ANSW(:NCHAR),'S') >0)  

      NVOX = NX * NY * NZ
C     COMPLAIN IF EXCESSIVE ALLOCATION
      CALL BIGALLOC(NVOX,IOK8,.FALSE.,.TRUE.,IRTFLG)

      ALLOCATE(BUFIN  (NVOX), 
     &         AVGARAY(NVOX), STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
          CALL ERRT(46,'ADDS_N; BUFIN..',2*NVOX)
          GOTO 9999
      ENDIF

      AVGARAY = 0.0  ! ZERO THE ARRAY

      IF (WANTA) THEN
C        OPEN NEW AVERAGE OUTPUT FILE
         MAXIM = 0
         CALL OPFILEC(LUNIN,ASKNAM,FILNAM,LUNA,'U',
     &                ITYPE,NX,NY,NZ, MAXIM,
     &                'AVERAGE',.TRUE.,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999

         IF (WANT3) THEN

C           OPEN NEW AVERAGE OUTPUT FILE SUBSET 1
            MAXIM = 0
            CALL OPFILEC(LUNIN,ASKNAM,FILNAM,LUNA1,'U',
     &                   ITYPE,NX,NY,NZ, MAXIM,
     &                   'AVERAGE  FILE FOR FIRST  SUBSET~',
     &                   .TRUE.,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999

            MAXIM = 0
            CALL OPFILEC(LUNIN,ASKNAM,FILNAM,LUNA2,'U',
     &                   ITYPE,NX,NY,NZ, MAXIM,
     &                   'AVERAGE  FILE FOR SECOND SUBSET~',
     &                   .TRUE.,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999

            ALLOCATE(AVGARAY1(NVOX),
     &               AVGARAY2(NVOX), STAT=IRTFLG)
            IF (IRTFLG .NE. 0) THEN
               MWANT = 2 * NVOX
               CALL ERRT(46,'ADDS; AVGARAY1..',MWANT)
               GOTO 9999
            ENDIF

            AVGARAY1 = 0.0  ! ZERO THE AVERAGE ARRAY 1
            AVGARAY2 = 0.0  ! ZERO THE AVERAGE ARRAY 2
         ENDIF
      ENDIF

      IF (WANTV) THEN
C        OPEN NEW VARIANCE OUTPUT FILE
         MAXIM = 0
         CALL OPFILEC(LUNIN,ASKNAM,FILNAM,LUNV,'U',
     &                ITYPE,NX,NY,NZ, MAXIM,
     &                'VARIANCE',.TRUE.,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999

         ALLOCATE(VARARAY(NVOX), STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'ADDS; VARARAY',NVOX)
            GOTO 9999
         ENDIF

         VARARAY = 0.0  ! ZERO THE VARIANCE ARRAY

         IF (WANT3) THEN

C           OPEN NEW VARIANCE OUTPUT FILE SUBSET 1
            MAXIM = 0
            CALL OPFILEC(LUNIN,ASKNAM,FILNAM,LUNV1,'U',
     &                   ITYPE,NX,NY,NZ, MAXIM,
     &                   'VARIANCE FILE FOR FIRST  SUBSET~',
     &                   .TRUE.,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999

            MAXIM = 0
            CALL OPFILEC(LUNIN,ASKNAM,FILNAM,LUNV2,'U',
     &                   ITYPE,NX,NY,NZ, MAXIM,
     &                   'VARIANCE FILE FOR SECOND SUBSET~',
     &                   .TRUE.,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999

            ALLOCATE(VARARAY1(NVOX),
     &               VARARAY2(NVOX), STAT=IRTFLG)
            IF (IRTFLG .NE. 0) THEN
               MWANT = 2 * NVOX
               CALL ERRT(46,'ADDS; VARARAY1..',MWANT)
               GOTO 9999
            ENDIF

            VARARAY1 = 0.0  ! ZERO THE VARIANCE  ARRAY 1
            VARARAY2 = 0.0  ! ZERO THE VARIANCE  ARRAY 2
         ENDIF
      ENDIF



      K      = 1
      NUMA   = 0
      NUMA1  = 0
      NUMA2  = 0

      DO 
         NUMA = NUMA + 1    ! TOTAL NUMBER OF IMAGES
 
C        READ INPUT IMAGE/VOLUME
         CALL REDVOL(LUNIN,NX,NY, 1,NZ, BUFIN,IRTFLG)

         AVGARAY = AVGARAY + BUFIN
         IF (WANTV) THEN
            VARARAY = VARARAY + BUFIN * BUFIN
         ENDIF

         IF (WANT3) THEN 
             IF (MOD(NUMA,2) .NE. 0) THEN  
                AVGARAY1 = AVGARAY1 + BUFIN 
                NUMA1    = NUMA1 + 1
                IF (WANTV) THEN
                   VARARAY1 = VARARAY1 + BUFIN * BUFIN
                ENDIF

             ELSE
                AVGARAY2 = AVGARAY2 + BUFIN 
                NUMA2    = NUMA2 + 1
                IF (WANTV) THEN
                   VARARAY2 = VARARAY2 + BUFIN * BUFIN
                ENDIF
             ENDIF
         ENDIF 

C        OPEN NEXT INPUT IMAGE
         CALL NEXTFILE(K,       ILIST, 
     &                 FOUROK,  0,
     &                 NIMG,    MAXIM,   
     &                 LUNIN, 0, 
     &                 FILPAT,  'O',
     &                 IMGNUM,  IRTFLG)
 
         IF (IRTFLG .EQ. -1) EXIT      !  END OF INPUT STACK
         IF (IRTFLG .NE. 0) RETURN

      ENDDO


      FNUMT   = NUMA 
      FNUMT1  = NUMA1 
      FNUMT2  = NUMA2 
      VV      = 0.0     ! FOR NO VAR. CASE

C     COMPUTE AVERAGE IMAGE 
      DTEMP   = 1.0 / FNUMT
      AVGARAY = AVGARAY * DTEMP   ! ARRAY DIVISION
      AVALL   = SUM(AVGARAY) / NVOX

      !WRITE(NOUT,*) 'numa: ',nvox,numa,avall

      IF (WANTA) THEN
C        SAVE OVERALL AVG. OUTPUT IMAGE
         BUFIN = AVGARAY    ! CONVERT TO REAL
         CALL WRTVOL(LUNA,NX,NY, 1,NZ, BUFIN,IRTFLG)

         IF (.NOT. WANTV .AND. .NOT. WANT3) THEN

            IF (NOUT .NE. NDAT) WRITE(NDAT,91) NUMA,NVOX, AVALL 
            WRITE(NOUT,91) NUMA,NVOX, AVALL 

91          FORMAT(/,
     &           '  FOR:',I7,' IMAGES   CONTAINING: ',I0,' ELEMENTS',/,
     &           '  AVERAGE  OF AVERAGE IMAGE:',1PG12.5,/)
         ENDIF
 
      ENDIF

      VV    = 0.0
      VARAV = 0.0
      SAVT  = 0.0
      VAV   = 0.0
      SAV   = 0.0

      IF (WANTV) THEN
C        COMPUTE OVERALL VARIANCE IMAGE
         VV      = DOT_PRODUCT(AVGARAY,AVGARAY)
         VV      = VV / REAL(NVOX)

         IF (NUMA > 1 ) THEN

           VARARAY = (VARARAY - AVGARAY*AVGARAY*FNUMT) / (FNUMT-1.0)
           VARAV   = SUM(VARARAY)
           VAV     = VARAV / (FLOAT(NVOX))

           IF (VARAV < 0.0) THEN
               WRITE(NOUT,*)   
     &            ' WARNING: NEGATIVE TOTAL VARIANCE: ',VARAV,
     &            '  SD SET TO ZERO!!!'
               IF (NOUT .NE. NDAT)
     &            WRITE(NDAT,*)
     &            ' WARNING: NEGATIVE TOTAL VARIANCE: ',VARAV,
     &                      '  SD SET TO ZERO!!!'
            ELSE
               SAVT = SQRT(VARAV)
            ENDIF

            IF (VAV < 0.0) THEN
               WRITE(NOUT,*)' WARNING: NEGATIVE  VARIANCE P.P.: ',VAV,
     &                      '   SD SET TO ZERO!!!'
               IF (NOUT .NE. NDAT)
     &            WRITE(NDAT,*) ' WARNING: NEGATIVE VARIANCE P.P:',VAV,
     &                          '   SD SET TO ZERO!!!'
            ELSE
               SAV = SQRT(VAV)
            ENDIF
         ENDIF

C        SAVE OVERALL VAR. OUTPUT IMAGE
         BUFIN = VARARAY    ! CONVERT TO REAL
         CALL WRTVOL(LUNV,NX,NY, 1,NZ, BUFIN,IRTFLG)

         IF (NOUT .NE. NDAT)
     &      WRITE(NDAT,7001) NUMA,NVOX, AVALL,VV,VARAV,SAVT,VAV,SAV
         WRITE(NOUT,7001)    NUMA,NVOX, AVALL,VV,VARAV,SAVT,VAV,SAV

7001     FORMAT(/,'  FOR:',I7,' IMAGES   CONTAINING: ',I0,' ELEMENTS',/,
     &        '  AVERAGE  OF AVERAGE IMAGE:',1PG12.5,/,
     &        '  VARIANCE OF AVERAGE IMAGE:',1PG12.5,/,
     &        '  TOTAL VARIANCE: ',1PG12.5,'  TOTAL S.D.: ',1PG12.5,/,
     &        '  P.P.  VARIANCE: ',1PG12.5,'  P.P.  S.D.: ',1PG12.5,/)

      ENDIF

      CALL REG_SET_NSEL(1,2, AVALL,VV,  0.0, 0.0, 0.0,IRTFLG)

      VV    = 0.0
      VARAV = 0.0
      SAVT  = 0.0
      VAV   = 0.0
      SAV   = 0.0

      IF (WANT3) THEN
C        NOW COMPUTE SUBSET AVERAGE IMAGES

         DTEMP    = 1.0 / FNUMT1
         AVGARAY1 = AVGARAY1 * DTEMP
         DTEMP    = 1.0 / FNUMT2
         AVGARAY2 = AVGARAY2  * DTEMP

         IF (WANTA) THEN
C           SAVE SUBSET AVG. OUTPUT IMAGES

C           COMPUTE AVERAGE IMAGE 1 (CONVERT TO REAL)
            BUFIN  = AVGARAY1  
            CALL WRTVOL(LUNA1,NX,NY, 1,NZ, BUFIN,IRTFLG)

C           COMPUTE AVERAGE IMAGE 2 (CONVERT TO REAL)
            BUFIN  = AVGARAY2      
            CALL WRTVOL(LUNA2,NX,NY, 1,NZ, BUFIN,IRTFLG)
         ENDIF
      ENDIF

      IF (WANT3 .AND. WANTV) THEN

         WRITE(NOUT,*)   ' FOR SUBSET 1 ----------------------'

         VV    = 0.0
         VARAV = 0.0
         SAVT  = 0.0
         VAV   = 0.0
         SAV   = 0.0

         VV    = DOT_PRODUCT(AVGARAY1,AVGARAY1)
         VV    = VV / REAL(NVOX)
         AVALL = SUM(AVGARAY1) / NVOX

         IF (NUMA1 > 1) THEN
C           COMPUTE VARIANCE IMAGE 1

            VARARAY1 = (VARARAY1 -AVGARAY1*AVGARAY1*FNUMT1)/(FNUMT1-1.0)
            VARAV    = SUM(VARARAY1)
            VAV      = VARAV / FLOAT(NVOX)

            !WRITE(nout,*) 'vararay: ',minval(vararay1), maxval(vararay1)

            IF (VARAV < 0.0) THEN
            WRITE(NOUT,*)   ' WARNING: NEGATIVE TOTAL VARIANCE: ',VARAV,
     &                      '  SD SET TO ZERO!!!'
            IF (NOUT .NE. NDAT)
     &         WRITE(NDAT,*)' WARNING: NEGATIVE TOTAL VARIANCE: ',VARAV,
     &                      '  SD SET TO ZERO!!!'
            ELSE
               SAVT = SQRT(VARAV)
            ENDIF

            IF (VAV < 0.0) THEN
            WRITE(NOUT,*)   ' WARNING: NEGATIVE  P.P.VARIANCE:',VAV,
     &                      '  SD SET TO ZERO!!!'
            IF (NOUT .NE. NDAT)
     &        WRITE(NDAT,*) ' WARNING: NEGATIVE P.P. VARIANCE:',VAV,
     &                      '  SD SET TO ZERO!!!'
            ELSE
               SAV  = SQRT(VAV)
            ENDIF
         ENDIF

C        SAVE SUBSET VAR. OUTPUT IMAGE
         BUFIN = VARARAY1    ! CONVERT TO REAL
         CALL WRTVOL(LUNV1,NX,NY, 1,NZ, BUFIN,IRTFLG)

         IF (NOUT .NE. NDAT)
     &      WRITE(NDAT,7001) NUMA1,NVOX, AVALL,VV,VARAV,SAVT,VAV,SAV
         WRITE(NOUT,7001)    NUMA1,NVOX, AVALL,VV,VARAV,SAVT,VAV,SAV


         WRITE(NOUT,*)   ' FOR SUBSET 2 ----------------------'

         VV    = 0.0
         VARAV = 0.0
         SAVT  = 0.0
         VAV   = 0.0
         SAV   = 0.0

C        COMPUTE VARIANCE IMAGE 2
         VV    = DOT_PRODUCT(AVGARAY1,AVGARAY1)
         VV    = VV / REAL(NVOX)
         AVALL = SUM(AVGARAY2) / NVOX

         IF (NUMA2 > 1) THEN

            VARARAY2 = (VARARAY2 -AVGARAY2*AVGARAY2*FNUMT2)/(FNUMT2-1.0)
            VARAV    = SUM(VARARAY1)
            VAV      = VARAV / FLOAT(NVOX)

            !WRITE(nout,*) 'vararay: ',minval(vararay2), maxval(vararay2)

            IF (VARAV < 0.0) THEN
            WRITE(NOUT,*)    'WARNING: NEGATIVE TOTAL VARIANCE:',VARAV,
     &                       '  SD SET TO ZERO!!!'
            IF (NOUT .NE. NDAT)
     &         WRITE(NDAT,*) 'WARNING: NEGATIVE TOTAL VARIANCE: ',VARAV,
     &                       '  SD SET TO ZERO!!!'
            ELSE
               SAVT = SQRT(VARAV)
            ENDIF

            IF (VAV < 0.0) THEN
            WRITE(NOUT,*)   ' WARNING: NEGATIVE  P.P. VARIANCE:',VAV,
     &                      '  SD SET TO ZERO!!!'
            IF (NOUT .NE. NDAT)
     &        WRITE(NDAT,*) ' WARNING: NEGATIVE P.P. VARIANCE:',VAV,
     &                      '  SD SET TO ZERO!!!'
            ELSE
               SAV  = SQRT(VAV)
            ENDIF
         ENDIF

C        SAVE SUBSET 2 VAR. OUTPUT IMAGE
         BUFIN = VARARAY2    ! CONVERT TO REAL
         CALL WRTVOL(LUNV2,NX,NY, 1,NZ, BUFIN,IRTFLG)

         IF (NOUT .NE. NDAT)
     &      WRITE(NDAT,7001) NUMA2,NVOX, AVALL,VV,VARAV,SAVT,VAV,SAV
         WRITE(NOUT,7001)    NUMA2,NVOX, AVALL,VV,VARAV,SAVT,VAV,SAV

         WRITE(NOUT,*)   ' '
      ENDIF

9999  CLOSE(LUNA)
      CLOSE(LUNA1)
      CLOSE(LUNA2)
      CLOSE(LUNV)
      CLOSE(LUNV1)
      CLOSE(LUNV2)

      IF (ALLOCATED(BUFIN))    DEALLOCATE(BUFIN)
      IF (ALLOCATED(VARARAY))  DEALLOCATE(VARARAY)
      IF (ALLOCATED(VARARAY1)) DEALLOCATE(VARARAY1)
      IF (ALLOCATED(VARARAY2)) DEALLOCATE(VARARAY2)
      IF (ALLOCATED(AVGARAY))  DEALLOCATE(AVGARAY)
      IF (ALLOCATED(AVGARAY1)) DEALLOCATE(AVGARAY1)
      IF (ALLOCATED(AVGARAY2)) DEALLOCATE(AVGARAY2)
      IF (ALLOCATED(ILIST))    DEALLOCATE(ILIST)

      END
