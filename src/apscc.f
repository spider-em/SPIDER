C++*********************************************************************
C
C APSCC            NEW                             FEB 08 ArDean Leith
C                  FLIPPED EXP & REF BUG           AUG 09 ArDean Leith
C                  REMOVED FROM APMASTER           AUG 09 ArDean Leith
C                  ISHX = 0 ALLOWED                FEB 10 ArDean Leith
C                  NORMIT BUG                      SEP 11 ArDean Leith
C                  NZP SIZE BUG                    SEP 14 ArDean Leith
C                  NEXTFILES USED, PADVAL BUG      SEP 14 ArDean Leith
C                  SIGR BUG                        DEC 14 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C APSCC
C
C PURPOSE:  2D & 3D PADDED, CROSS CORRELATION MULTI-REFERENCE 
C           SHIFT ALIGNMENT
C  
C VARIABLES:
C       IREFLST               LIST OF REF. IMAGE FILE NUMBERS   (INPUT)
C       NUMREF                NO. OF IMAGES                     (INPUT)
C       IEXPLST               LIST OF EXP. IMAGE FILE NUMBERS   (INPUT)
C       NUMEXP                NO. OF IMAGES                     (INPUT)
C       NX,NY,NZ              ACTUAL (NOT-WINDOWED) IMAGE SIZE  (INPUT)
C       ISHX,ISHY,ISHZ        SHIFT SEARCH RANGE                (INPUT)
C       NORMIT                NORMALIZATION WANTED FLAG         (INPUT)
C       REFPAT                REF. IMAGE SERIES FILE TEMPLATE   (INPUT)
C       EXPPAT                EXP. IMAGE SERIES FILE TEMPLATE   (INPUT)
C       LUNREF,LUNEXP,LUNDOC  I/O UNITS                         (INPUT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE APSCC()

        IMPLICIT NONE

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        CHARACTER (LEN=MAXNAM) :: REFPAT,EXPPAT,OUTANG
        INTEGER                :: IEXPLST(NIMAX)

	CHARACTER(LEN=1)       :: NULL = CHAR(0)
	CHARACTER(LEN=80)      :: YN
	CHARACTER(LEN=160)     :: COMMEN
        LOGICAL                :: NEWFILE,NORMIT,PADIT
        INTEGER                :: NILMAX,NX,NY,NZ,NDUM,NPADVAL
        INTEGER                :: IDUM,ITYPER,ITYPEE,NUMREF,IREF
        INTEGER                :: NUMEXP,IEXP,IRTFLG,NLET,NOT_USED
        INTEGER                :: ITYPE,MAXIME,MAXIMR,NSEL_USED,ipad
        REAL                   :: UNUSED
        REAL                   :: AVR,SIGR,AVE,SIGE

        INTEGER                :: ICOMM,MYPID,MPIERR
        INTEGER                :: ISHX,ISHY,ISHZ,NX2,NY2,NZ2,NOUTANG

        LOGICAL                :: FOUROK    = .TRUE.

        INTEGER,PARAMETER      :: LUNREF    = 21 
        INTEGER,PARAMETER      :: LUNEXP    = 22
        INTEGER,PARAMETER      :: LUNDOCSEL = 81
        INTEGER,PARAMETER      :: LUNXMR    = 82
        INTEGER,PARAMETER      :: LUNXME    = 83
        INTEGER,PARAMETER      :: LUNDOC    = 84

	!DATA  INPIC,INANG,NDOC,NSCF/77,78,55,50/ !USED IN CALLED ROUTINE

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID 

        NILMAX = NIMAX  ! FROM CMLIMIT.INC

C       ASK FOR TEMPLATE AND NUMBERS FOR REFERENCE IMAGES
C       OPEN FIRST REFERENCE IMAGE TO BE ALIGNED
        CALL OPFILES(0,LUNREF,LUNDOCSEL,LUNXMR, 
     &        .TRUE.,REFPAT,NLET, 'O',
     &        ITYPER,NX,NY,NZ,MAXIMR,
     &        'REFERENCE IMAGE FILE NAME OR TEMPLATE (E.G. STK@****)~',
     &        FOUROK,INUMBR,NILMAX, 
     &        NDUM,NUMREF,IREF, IRTFLG) 
        IF (IRTFLG .NE. 0) RETURN

C       NUMREF - TOTAL NUMBER OF REF. IMAGES
        IF (NUMREF == 0 .AND. IREF > 0) NUMREF = 1
        IF (NUMREF <= 0)  THEN
           CALL ERRT(101,'No reference images',IDUM)
           GOTO 9999
        ELSEIF (MYPID <= 0) THEN
           WRITE(NOUT,'(A,I0)')'  Number of reference images: ',NUMREF
        ENDIF


C       FIND SEARCH RANGE
        ISHX = 0
        ISHY = 0
        ISHZ = 0 
        CALL RDPRI3S(ISHX,ISHY,ISHZ,NOT_USED,
     &     'SEARCH RANGE IN X, Y, & Z (ZERO FOR ALL)' ,IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9999

C       POSSIBLE SHIFT AMOUNTS FOR SEARCH
        IF (ISHX == 0) ISHX = NX / 2
        IF (ISHY == 0) ISHY = NY / 2
        IF (ISHZ == 0) ISHZ = NZ / 2
        IF (ISHX <  0) ISHX = 0
        IF (ISHY <  0) ISHY = 0
        IF (ISHZ <  0) ISHZ = 0

        
        IRTFLG = 0
        IF (ITYPER > 0)  THEN
C          REAL IMAGE INPUT
           CALL RDPRMC(YN,NLET,.TRUE.,
     &       'NORMALIZE PEAK HEIGHT, 2x PAD IMAGES (Y/N)', NULL,IRTFLG)
        ELSE
C          FOURIER INPUT (MUST ALREADY BE PADDED IF DESIRED)
           CALL RDPRMC(YN,NLET,.TRUE.,
     &        'NORMALIZE PEAK HEIGHT, 2x PADDED IMAGES (Y/N)', 
     &         NULL,IRTFLG)
        ENDIF
        IF (IRTFLG .NE. 0)  GOTO 9999
        NORMIT = ( YN(1:1) .NE. 'N' )

        PADIT  = (INDEX(YN(2:NLET),'N') == 0 )

        AVR  = AV
        SIGR = SIG
        IF (IMAMI .NE. 1 .AND. ITYPER > 0 ) THEN
           CALL NORM3(LUNREF,NX,NY,NZ,FMAX,FMIN,AV)
           AVR  = AV
           SIGR = SIG
        ENDIF


C       OPEN FIRST EXPERIMENTAL IMAGE TO BE ALIGNED
        CALL OPFILES(0,LUNEXP,LUNDOCSEL,LUNXME, 
     &      .TRUE.,EXPPAT,NLET, 'O',
     &      ITYPEE,NX2,NY2,NZ2,MAXIME,
     &      'EXPERIMENTAL IMAGE FILE NAME OR TEMPLATE (E.G. STK@****)~',
     &      FOUROK,IEXPLST,NILMAX, 
     &      NDUM,NUMEXP,IEXP, IRTFLG) 
        IF (IRTFLG .NE. 0)  GOTO 9999

C       IMAGES MUST BE SAME SIZE AND IFORM 
        CALL SIZCHK(UNUSED, NX,NY,NZ,ITYPER, NX2,NY2,NZ2,ITYPEE, IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (NUMEXP == 0 .AND. IEXP > 0) NUMEXP = 1
        IF (MYPID <= 0) THEN
          WRITE(NOUT,'(A,I0)')'  Number of experimental images: ',NUMEXP
        ENDIF

        AVE  = AV
        SIGE = SIG
        IF (IMAMI .NE. 1 .AND. ITYPER > 0) THEN 
           CALL NORM3(LUNEXP,NX,NY,NZ,FMAX,FMIN,AV)
           AVE  = AV
           SIGE = SIG
        ENDIF

C       GET NAME FOR OUTPUT DOC FILE
        CALL REG_GET_USED(NSEL_USED)

C       OPEN OUTPUT DOC FILE (FOR APPENDING)
        NOUTANG = LUNDOC
        CALL OPENDOC(OUTANG,.TRUE.,NLET,LUNDOC,NOUTANG,.TRUE.,
     &           'OUTPUT ALIGNMENT DOCUMENT',.FALSE.,.TRUE.,.TRUE.,
     &            NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        IF (IRTFLG == -1) THEN
C          DO NOT WANT OUTPUT DOC FILE
           NOUTANG = 0
        ELSE
C          WANT OUTPUT DOC FILE
           COMMEN = '       ' //
     &         'EXP#,         REF#,         SX,           SY,       '//
     &         '    SZ,           PEAK'
           CALL LUNDOCPUTCOM(NOUTANG,COMMEN,IRTFLG)
        ENDIF


        CALL APSCC_DO(IREF,INUMBR,NUMREF ,IEXP,IEXPLST,NUMEXP,
     &              NX,NY,NZ, ISHX,ISHY,ISHZ, NORMIT, PADIT,
     &              AVE,SIGE,AVR,SIGR, MAXIMR,MAXIME,ITYPER,
     &              REFPAT,EXPPAT,LUNREF,LUNEXP,NOUTANG,LUNXME,LUNXMR)

9999    CLOSE(LUNDOC)
        CLOSE(LUNEXP)
        CLOSE(LUNREF)

        END

C       ------------------- APSCC_DO ----------------------------------

        SUBROUTINE APSCC_DO(IREF,IREFLST,NUMREF, IEXP,IEXPLST,NUMEXP,
     &             NXT,NYT,NZT,  ISHX,ISHY,ISHZ, NORMIT, PADIT,
     &             AVE,SIGE,AVR,SIGR, MAXIMR,MAXIME,ITYPE, 
     &             REFPAT,EXPPAT, LUNREF,LUNEXP,LUNDOC,LUNXME,LUNXMR)

        IMPLICIT NONE

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        INTEGER                :: IEXP,IREF
	INTEGER                :: IREFLST(NUMREF)
	INTEGER                :: NUMREF
	INTEGER                :: IEXPLST(NUMEXP) 
	INTEGER                :: NUMEXP 
	INTEGER                :: NXT,NYT,NZT, ISHX,ISHY,ISHZ
        LOGICAL                :: NORMIT,PADIT
        REAL                   :: AVR,SIGR,AVE,SIGE
        INTEGER                :: MAXIMR,MAXIME,ITYPE
        CHARACTER (LEN=*)      :: REFPAT,EXPPAT
	INTEGER                :: LUNREF,LUNEXP,LUNDOC,LUNXME,LUNXMR  

C       ALLOCATABLE ARRAYS
	REAL, ALLOCATABLE      :: BUFI(:),BUFR(:),BUFE(:)

C       AUTOMATIC ARRAYS
        REAL                   :: DLIST(6)

	INTEGER                :: NSEL_USED,NXP,NYP,NZP,LSE,IRTFLG
	INTEGER                :: NINDXR,NINDXE,IKEY
	INTEGER                :: MWANT,NLET,MAXIM,LX,LY,LZ  
	INTEGER                :: INV,NPADVAL  
	INTEGER                :: NX,NY,NZ

        REAL                   :: UNUSED,PADVAL
        REAL                   :: XSHNEW,YSHNEW,ZSHNEW,PEAKV

        LOGICAL                :: SPIDER_SIGN  = .FALSE.
        LOGICAL                :: SPIDER_SCALE = .FALSE.

        LOGICAL                :: APCC_NORM    = .TRUE.
        LOGICAL                :: DO_FFT_I     = .TRUE.
        LOGICAL                :: DO_FFT_R     = .TRUE.
        LOGICAL                :: FOUROK       = .TRUE.

        DOUBLE PRECISION       :: DAV,DSIG 

        CALL REG_GET_USED(NSEL_USED)

        NPADVAL = 1
        IF (PADIT) NPADVAL = 2 

        IF ( ITYPE >= 0 ) THEN
C          REAL SPACE INPUT IMAGES, MAY WANT PADDED x 2
           NX  = NXT
           NXP = NX * NPADVAL
           LSE = NXP + 2 - MOD(NXP,2)

           NY  = NYT
           NYP = NYT * NPADVAL
           NZ  = NZT
           NZP = NZT * NPADVAL

        ELSEIF ( ITYPE == -11 .OR. ITYPE == -21 ) THEN
C          FOURIER INPUT IMAGES, SINCE ODD WAS NOT PADDED x 2
           NX  = NXT - 1
           NXP = NX 
           LSE = NXT

           NY  = NYT 
           NYP = NYT  
           NZ  = NZT 
           NZP = NZT  

        ELSEIF ( ITYPE == -12 .OR. ITYPE == -22 ) THEN
C          FOURIER INPUT IMAGES, MAY HAVE BEEN PADDED x 2
           NXP = NXT - 2
           NX  = NXP / NPADVAL
           LSE = NXT

           NY  = NYT / NPADVAL
           NYP = NYT  
           NZ  = NZT / NPADVAL
           NZP = NZT

        ELSE
           CALL ERRT(101,'UNKNOWN IMAGE FOURIER TYPE',ITYPE)
           GOTO 9999
        ENDIF

        IF (NY == 1) NYP = 1
        IF (NZ == 1) NZP = 1

        ! write(6,*) 'nx,ny,nz:',nx,ny,nz, nxp,nyp,nzp
        ! write(6,*) 'lse,itype:',lse,itype

C       MAKE BUFI & BUFR FOR PADDED IMAGES WITH FOURIER ROW LENGTH
	ALLOCATE(BUFR(LSE*NYP*NZP), 
     &           BUFI(LSE*NYP*NZP),
     &           STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           MWANT = 2 * LSE * NYP * NZP
           CALL ERRT(46,'BUFI & BUFR...',MWANT)
           GOTO 9999
        ENDIF

        IF (NUMREF > 1) THEN
C          MUST SAVE BUFI SINCE IT IS OVERWRITTEN BY APCC
	   ALLOCATE(BUFE(LSE*NYP*NZP),   STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
              MWANT = LSE * NYP * NZP
              CALL ERRT(46,'BUFE',MWANT)
              GOTO 9999
           ENDIF
        ENDIF

        NINDXE = 1
        NINDXR = 1
        IKEY   = 0

        DO 

C          LOAD & PAD EXP. IMAGE INTO FOURIER SIZE BUFFER
           IF (ITYPE >= 0) THEN

              PADVAL = AVE        ! CAN NOT USE PADVAL==ZERO
              CALL REDNPADVOL(LUNEXP,PADVAL,
     &                        NX,NY,NZ,
     &                        LSE, NYP,NZP,
     &                        BUFI, IRTFLG)
           ELSE
              CALL REDVOL(LUNEXP,LSE,NY,1,NZ, BUFI,IRTFLG)
           ENDIF

           !call chkfile('jnkpaded',66,1,lse,nyp,nzp, bufi,irtflg)

           IF (NUMREF > 1) THEN
C             MUST SAVE BUFI SINCE IT IS OVERWRITTEN BY APCC
              BUFE = BUFI
           ENDIF

           NINDXR = 1
           DO     ! LOOP OVER REFERENCES

C             LOAD & PAD REF. IMAGE INTO FOURIER SIZE BUFFER
              IF (ITYPE >= 0) THEN

                 PADVAL = AVR        ! CAN NOT USE PADVAL==ZERO
                 CALL REDNPADVOL(LUNREF,PADVAL,
     &                        NX,NY,NZ,
     &                        LSE, NYP,NZP,
     &                        BUFR, IRTFLG)
              ELSE
                 CALL REDVOL(LUNREF,LSE,NY,1,NZ, BUFR,IRTFLG)
              ENDIF

   
C             CROSS CORRELATION --------------------------------------- CC

C             APCC RETURNS PEAK IMAGE IN: BUFR, BUFI(0,0) SET TO 0,0
              CALL APCC_NEW(LSE, NXP,NYP,NZP, BUFI,BUFR,
     &                      DO_FFT_I,DO_FFT_R, 
     &                      .FALSE.,APCC_NORM,SPIDER_SIGN, 
     &                      ISHX,ISHY,ISHZ,
     &                      XSHNEW,YSHNEW,ZSHNEW, PEAKV,IRTFLG)

              IF (IRTFLG .NE. 0)  GOTO 9999
 
              IF (NORMIT) THEN
C                NORMALIZATION 
C                USING PADDED SIGS DOES NOT GIVE BETTER CC VALUE!! al
                 PEAKV = PEAKV /FLOAT(NX*NY*NZ-1)/SIGE/SIGR
              ENDIF

              IKEY     = IKEY + 1
	      DLIST(1) = IEXPLST(IEXP)
              DLIST(2) = IREFLST(IREF)
              DLIST(3) = XSHNEW
              DLIST(4) = YSHNEW
              DLIST(5) = ZSHNEW
              DLIST(6) = PEAKV

              IF (LUNDOC > 0) THEN
C                 SAVE IN ALIGNMENT DOC FILE
C                 IMG#, REF#, SX,SY,SZ, PEAK-HEIGHT
                  CALL LUNDOCWRTDAT(LUNDOC,IKEY,DLIST,6,IRTFLG)
                  IF (IRTFLG .NE. 0) GOTO 9999
              ENDIF
 
              IF (NSEL_USED > 0) THEN
C                 OUTPUT TO SPIDER'S REGISTERS
                  CALL REG_SET_NSEL(1,5, DLIST(1),DLIST(2),DLIST(3),
     &                                   DLIST(4),DLIST(5),IRTFLG)
                  CALL REG_SET_NSEL(6,1, DLIST(6),0.0,0.0,
     &                                   0.0,0.0,IRTFLG)
              ENDIF

              IF (NUMREF > 1) THEN
C                MUST RETRIEVE BUFI SINCE IT IS OVERWRITTEN BY APCC
                 BUFI = BUFE
              ENDIF

C             OPEN NEXT REFERENCE INPUT FILE, UPDATE NINDXR 
              CALL NEXTFILE(NINDXR, IREFLST, 
     &                      FOUROK,  LUNXMR,
     &                      NUMREF,  MAXIMR,  
     &                      LUNREF,  0,  
     &                      REFPAT,  'O',
     &                      IREF,    IRTFLG)

              !write(6,*) 'iref, irtflg:',iref, irtflg,nindxe


             IF (IRTFLG .NE. 0) EXIT      ! ERROR / END OF INPUT STACK

             IF (IMAMI .NE. 1 .AND. ITYPE > 0 .AND. ITYPE > 0) THEN
                CALL NORM3(LUNREF,NX,NY,NZ,FMAX,FMIN,AV)
                AVR  = AV
                SIGR = SIG
             ENDIF

           ENDDO

           IF (IRTFLG < 0) THEN 
C             OPEN NEXT REFERENCE INPUT FILE, RESET NINDXR 
              NINDXR = 0
              CALL NEXTFILE(NINDXR, IREFLST, 
     &                      FOUROK,  LUNXMR,
     &                      NUMREF,  MAXIMR,  
     &                      LUNREF,  0,  
     &                      REFPAT,  'O',
     &                      IREF,    IRTFLG)
             IF (IRTFLG .NE. 0) EXIT      ! ERROR / END OF INPUT STACK
          ENDIF

          !write(6,*) 'iref, irtflg:',iref, irtflg,nindxe
          !write(6,*) 'calling nextfile exp, iexp, irtflg:',iexp, irtflg

C          OPEN NEXT EXP. INPUT FILE, UPDATE NINDXR 

           CALL NEXTFILE(NINDXE, IEXPLST, 
     &                   FOUROK,  LUNXME,
     &                   NUMEXP,  MAXIME,  
     &                   LUNEXP,  0,  
     &                   EXPPAT,  'O',
     &                   IEXP,    IRTFLG)

           IF (IRTFLG .NE. 0) EXIT      ! ERROR / END OF INPUT STACK

           AVE  = AV
           SIGE = SIG
           IF (IMAMI .NE. 1 .AND. ITYPE > 0 .AND. ITYPE > 0) THEN
              CALL NORM3(LUNEXP,NX,NY,NZ,FMAX,FMIN,AV)
              AVE  = AV
              SIGE = SIG
           ENDIF

        ENDDO

9999    IF (ALLOCATED(BUFI))   DEALLOCATE(BUFI)
        IF (ALLOCATED(BUFR))   DEALLOCATE(BUFR)
        IF (ALLOCATED(BUFE))   DEALLOCATE(BUFE)

        END








C        ----------------- UNUSED BELOW HERE --------------------------

#ifdef NEVER
#ifdef DEBUG
           write(6,*)'  '
           write(6,*)'peakv: ',peakv

           !peakv = peakv / float(nxp*nyp*nzp)
           !write(6,*)'  '
           !write(6,*)'peakv /(nxp*nyp*nzp): ',peakv

           write(6,*)'peakv/sigef/sigrf: ',peakv/sigef/sigrf
           write(6,*)'peakv/sige/sigr:   ',peakv/sige/sigr
           nnn  = nx*ny*nz
           nnnp = nxp*nyp*nzp

           write(6,*)'peakv/float(nnn)/sigef/sigrf: ',
     &                peakv/float(nnn)/sigef/sigrf

           write(6,*)'peakv/float(nnnp)/sigef/sigrf: ',
     &                peakv/float(nnnp)/sigef/sigrf

           write(6,*)'peakv/float(nnn)/sige/sigr: ',
     &                peakv/float(nnn)/sige/sigr

           write(6,*)'peakv/float(nnnp)/sige/sigr: ',
     &                peakv/float(nnnp)/sige/sigr
           write(6,*)'  '
           write(6,*)'final peakv: ',peakv
#endif

c----------------------
        maxim = 0
        itype = 3
        call opfilec(0,.false.,'jnkexppad',98,'u',itype,
     &                lse,nyp,nzp,maxim,' ',.false.,irtflg)
        call wrtvol(98,lse,nyp, 1,nzp, bufi,irtflg)
        close(98)

        maxim = 0
        itype = 3
        call opfilec(0,.false.,'jnkrefpad',98,'u',itype,
     &                lse,nyp,nzp,maxim,' ',.false.,irtflg)
        call wrtvol(98,lse,nyp,1,nzp, bufr,irtflg)
        close(98)
c-----------------------
#endif


