
C++*********************************************************************
C
C PR3D_NEW.F        CREATED FROM PR3D(FSCOP)        APR 16 ARDEAN LEITH
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
C PR3D_NEW()
C
C PURPOSE:  CALCULATE 3-D PHASE RESIDUE OUTSIDE MISSING CONE, 
C           PR OF FOURIER RINGS(RADIUS, DIRECTION RELATIVE TO Z) AND 
C           OF FOURIER SHELLS (RADIUS). CAN USE A MASK
C
C OPERATIONS:  'FSC NEW'
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE PR3D_NEW()

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        CHARACTER(LEN=MAXNAM) :: FILNAM1,FILNAM2 ,FILNAMM1,FILNAMM2

        REAL, ALLOCATABLE     :: AIMG(:,:,:),BIMG(:,:,:),CIMG(:,:,:)
        REAL, ALLOCATABLE     :: DIMG(:,:,:)
        REAL, ALLOCATABLE     :: CSUM1(:), CSUM(:)
        REAL, ALLOCATABLE     :: PR(:,:), AMP(:,:), AVSUM(:,:)
        REAL, ALLOCATABLE     :: CSUM2(:)
        INTEGER, ALLOCATABLE  :: LR(:)

        REAL                  :: SSANG, WI, WIP, VOXSIZ, FMAXSPFREQ
        REAL                  :: RADMASK
        LOGICAL               :: MASKOUT, WANTMASK  
        LOGICAL, PARAMETER    :: FSCOP     = .TRUE.
        LOGICAL, PARAMETER    :: WANTSQRTS = .TRUE.

        CHARACTER(LEN=1)      :: SER
        CHARACTER(LEN=1)      :: NULL = CHAR(0)

        INTEGER               :: ICOMM,MYPID,MPIERR

        INTEGER               :: ITYPE1,ITYPE2,ITYPE3,ITYPE4 
        INTEGER               :: NX1,NY1,NZ1,IRTFLG,NE,INC 
        INTEGER               :: NX,NX2,NY2,NZ2 
        REAL                  :: UNUSED,SCALE1,SCALE2,FACT,DSCALE,Y1  
        INTEGER               :: LSD, NY, NZ, NOT_USED, NUMC 
        INTEGER               :: MWANT,K,J,I,INV,NLEN
        REAL                  :: XCEN, YCEN, ZCEN, TNM, RADMASKSQI,EEE 
  
        INTEGER               :: MAXIM
        INTEGER, PARAMETER    :: NSCALE     = 20
        INTEGER, PARAMETER    :: LUN1       = 21
        INTEGER, PARAMETER    :: LUN2       = 22
        INTEGER, PARAMETER    :: LUNMASKIN  = 23
        INTEGER, PARAMETER    :: LUNMASKOUT = 24
        INTEGER, PARAMETER    :: LUNGP      = 88
        INTEGER, PARAMETER    :: LUNDOC     = 89

        INTEGER               :: lnblnkn

        CALL SET_MPI(ICOMM,MYPID,MPIERR)

C       INPUT FIRST INPUT VOLUME
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM1,LUN1,'O',ITYPE1,
     &               NX1,NY1,NZ1,MAXIM,
     &               'FIRST INPUT VOLUME~',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (NZ1 < 2) THEN
C          ONLY FOR VOLUMES
           CALL ERRT(101,'OPERATION ONLY FOR VOLUMES',NE)
           GOTO 9998
        ENDIF

C       INPUT SECOND INPUT VOLUME
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM2,LUN2,'O',ITYPE2,
     &               NX2,NY2,NZ2,MAXIM,
     &               'SECOND INPUT VOLUME~',.FALSE.,IRTFLG)

        IF (IRTFLG .NE. 0) GOTO 9998

C       FILES MUST HAVE SAME DIMENSIONS
        CALL SIZCHK(UNUSED,NX1,NY1,NZ1,0,
     &                     NX2,NY2,NZ2,0,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        NX         = NX1
        NY         = NY1  
        NZ         = NZ1
        LSD        = NX1 + 2 - MOD(NX1,2) ! FOURIER X SIZE

        SCALE1     = 0.2     ! SCALE FACTOR 1
        SCALE2     = 2.0     ! SCALE FACTOR 2
        DSCALE     = (SCALE2-SCALE1) / FLOAT(NSCALE-1)
        SER        = 'C'     ! CONE
        SSANG      = 90.0    ! MAX. TILT ANGLE
        FACT       = 3       ! FACTOR FOR NOISE COMPARISON

        WIP        = 2.0     ! SHELL WIDTH, DEFAULT VALUE
        RADMASK    = 0       ! USE FILE FOR MASK, DEFAULT VALUE

!     'SHELL WIDTH (VOXELS), MASKING RADIUS (0 == USE FILE)', hidden
        CALL RDPRM2S(WIP, RADMASK, NOT_USED,
     &               'SHELL WIDTH (VOXELS)', IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        WI      = 1.0 / WIP
        MASKOUT = (RADMASK > 0)

        VOXSIZ = 1      ! VOXEL SIZE, DEFAULT VALUE
        CALL RDPRM1S(VOXSIZ, NOT_USED,'VOXEL SIZE (A)',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        FMAXSPFREQ = 0.5 / VOXSIZ

        Y1         = FLOAT(MAX(NX,NY,NZ))
        INC        = INT(Y1/WI) / 2+1

        ALLOCATE (AIMG(LSD,NY,NZ),
     &            BIMG(LSD,NY,NZ), 
     &            CIMG(LSD,NY,NZ), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MWANT = 2*LSD*NY*NZ + NX*NY*NZ  
           CALL ERRT(46,'PR3D_MA; AIMG...',MWANT)
           GOTO 9999
        ENDIF

C       LOAD THE TWO INPUT VOLUMES
        CALL READV(LUN1,AIMG, LSD,NY,NX,NY,NZ)
        CALL READV(LUN2,BIMG, LSD,NY,NX,NY,NZ)

        MASKOUT = .FALSE. 
        IF (RADMASK <= 0) THEN
C          WANT MASK INPUT
           CALL OPFILEC(0,.TRUE.,FILNAMM1,LUNMASKIN,'O',ITYPE3,
     &                  NX2,NY2,NZ2,MAXIM,
     &                  'INPUT MASK VOLUME (* IF NONE)~',.TRUE.,IRTFLG)
           IF (IRTFLG > 0) GOTO 9999
           NLEN     = lnblnkn(FILNAMM1)
           WANTMASK = (NLEN > 1 .AND. FILNAMM1(1:1) .NE. '*')

           IF (WANTMASK) THEN


C             FILES MUST HAVE SAME DIMENSIONS
              CALL SIZCHK(UNUSED,NX, NY, NZ, 0,
     &                           NX2,NY2,NZ2,0,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999

C             LOAD THE INPUT MASK VOLUME
              CALL READV(LUNMASKIN,CIMG, LSD,NY,NX,NY,NZ)

C             MASK THE INPUT VOLUMES, ARRAY OPERATIONS
              AIMG = AIMG * CIMG
              BIMG = BIMG * CIMG
           ENDIF

       ELSEIF (RADMASK > 0) THEN
C          SUPERGAUSSIAN MASKING WANTED

           MAXIM  = 0
           ITYPE4 = 3
           CALL OPFILEC(0,.TRUE.,FILNAMM2,LUNMASKOUT,'U',ITYPE4,
     &                  NX,NY,NZ,MAXIM,
     &                  'OUTPUT MASK VOLUME (* IF NONE)~',
     &                  .FALSE.,IRTFLG)

           IF (IRTFLG > 0) GOTO 9999

           NLEN    = lnblnkn(FILNAMM2)
           MASKOUT = (NLEN > 1 .AND. FILNAMM2(1:1) .NE. '*')

           XCEN       = (NX/2) + 1
           YCEN       = (NY/2) + 1
           ZCEN       = (NZ/2) + 1
	   TNM        = ALOG(1.0 / TINY(TNM))
           RADMASKSQI = 1.0 / (RADMASK**2)

           IF ( MASKOUT ) THEN
C            WANT MASK VOLUME OUTPUT
 
             ALLOCATE (DIMG(NX,NY,NZ), STAT=IRTFLG)
             IF (IRTFLG .NE. 0) THEN
                 MWANT = NX*NY*NZ  
                 CALL ERRT(46,'PR3D_MA; DIMG',MWANT)
                 GOTO 9999
             ENDIF
           ENDIF

           DO K = 1,NZ
              DO J = 1,NY
                 DO I = 1,NX
	           EEE = 0.5 * ((I - XCEN) **2 * RADMASKSQI +
     &                          (J - YCEN) **2 * RADMASKSQI +
     &                          (K - ZCEN) **2 * RADMASKSQI)

	           IF (EEE  >= TNM) THEN
C                     OUTSIDE OF POSITIVE MASK AREA
	              AIMG(I,J,K) = 0.0
	              BIMG(I,J,K) = 0.0
                      IF ( MASKOUT ) DIMG(I,J,K) = 0.0
	           ELSE
C                     INSIDE OF POSITIVE MASK AREA 
	              EEE         = -0.5 * (2*EEE)**2
                      AIMG(I,J,K) = EXP(EEE) * AIMG(I,J,K)
                      BIMG(I,J,K) = EXP(EEE) * BIMG(I,J,K)

C                     FILL OUTPUT MASK
                      IF ( MASKOUT ) DIMG(I,J,K) = EXP(EEE)
	           ENDIF
                ENDDO
             ENDDO
           ENDDO

           IF ( MASKOUT ) THEN
C             SAVE MASK VOLUME 
              CALL WRTVOL(LUNMASKOUT,NX,NY,1,NZ,DIMG,IRTFLG)
	      CLOSE(LUNMASKOUT)
           ENDIF
        ENDIF

 
C       INPUT FILES ARE REAL SPACE, CONVERT TO  FOURIER

        INV = 1
        CALL FMRS_3(AIMG,NX1,NY,NZ,INV)
        IF (INV == 0) THEN
           CALL ERRT(101,'FFTW ERROR',NE)
           GOTO 9999
        ENDIF

        INV = 1
        CALL FMRS_3(BIMG,NX2,NY,NZ,INV)
        IF (INV == 0)THEN
           CALL ERRT(101,'FFTW ERROR',NE)
           GOTO 9999
        ENDIF 

        ALLOCATE(PR(NSCALE,INC), AMP(NSCALE,INC), CSUM1(INC),
     &           LR(INC),CSUM(INC),CSUM2(INC), AVSUM(NSCALE,INC),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           MWANT = 3*NSCALE*INC + 4*INC
           CALL ERRT(46,'PR3D; PR...',MWANT)
           GOTO 9999
        ENDIF

C       CALCULATIONS
        CALL PR3DB(AIMG,BIMG,PR,AMP,CSUM1,LR,CSUM,CSUM2,
     &      AVSUM,LSD,NX,NY,NZ,DSCALE,NSCALE,SCALE1,
     &      SSANG,INC,Y1,WI,SER)

C       WRITE RESULT INTO DOC FILE AND RESULT FILE       
        CALL RFACTSD2_NEW(PR,AMP,CSUM1,LR,CSUM,CSUM2,AVSUM,
     &                   NSCALE,INC,FACT,
     &                   LUNDOC,FMAXSPFREQ,LUNGP)


9999	IF (ALLOCATED(AIMG))  DEALLOCATE (AIMG)
        IF (ALLOCATED(BIMG))  DEALLOCATE (BIMG)
        IF (ALLOCATED(CIMG))  DEALLOCATE (CIMG)
        IF (ALLOCATED(DIMG))  DEALLOCATE (DIMG)
        IF (ALLOCATED(PR))    DEALLOCATE (PR)
        IF (ALLOCATED(AMP))   DEALLOCATE (AMP)
        IF (ALLOCATED(CSUM1)) DEALLOCATE (CSUM1)
        IF (ALLOCATED(LR))    DEALLOCATE (LR)
        IF (ALLOCATED(CSUM))  DEALLOCATE (CSUM)
        IF (ALLOCATED(CSUM2)) DEALLOCATE (CSUM2)
        IF (ALLOCATED(AVSUM)) DEALLOCATE (AVSUM)

9998    CLOSE(LUN1)
        CLOSE(LUN2)
        CLOSE(LUNDOC)
        CLOSE(LUNGP)
        CLOSE(LUNMASKIN)

        END



C        ------------  RFACTSD2_NEW----------------------------------

         SUBROUTINE RFACTSD2_NEW(PR,AMP,CSUM1,LR,CSUM,CSUM2,AVSUM,
     &                          NSCALE,INC,FACT,
     &                          LUNDOC,FMAXSPFREQ,LUNGP )

         IMPLICIT NONE

         INCLUDE 'CMLIMIT.INC' 
         INCLUDE 'CMBLOCK.INC'
 
         REAL     :: PR(NSCALE,INC),AMP(NSCALE,INC),CSUM1(INC)
         INTEGER  :: LR(INC)
         REAL     :: CSUM(INC),CSUM2(INC),AVSUM(NSCALE,INC)
         REAL     :: RVAL(INC),FSCVAL(INC)
         INTEGER  :: NSCALE,INC
         REAL     :: WI,FACT
         INTEGER  :: LUNDOC
         REAL     :: FMAXSPFREQ
         INTEGER  :: LUNGP,ILOCMIN
 
         CHARACTER (LEN=MAXNAM)   :: DOCNAM,GPLOTFILE,PROMPT
         CHARACTER (LEN=2*MAXNAM) :: MSG

         INTEGER, PARAMETER :: NLIST = 9
         REAL     :: DLIST(NLIST)

         INTEGER  :: ICOMM,MYPID,MPIERR
         INTEGER  :: NSEL_USED, NSCM, NSC, NPIXLEND2, NLEN, K 
         INTEGER  :: NLET, NICDOC, IRTFLG, LNBLNKN, NMSG, L, NVOX 

         REAL     :: FSCMIN,SPFMIN
         REAL     :: XPREVL, XPREVG, FSCNOWG, FSCNOWL,FSCPREVG,FSCPREVL  
         REAL     :: SPFINTERPG, SPFINTERPL, RESOLG, RESOLL 
         REAL     :: SPFPREVG, SPFPREVL, SPFNOWG, SPFNOWL
         REAL     :: XINTERPG,XINTERPL 
         REAL     :: FSCLAST, WIP  
         REAL     :: RFM, BK1, BK2, BK3, BK4, FSCZ
         REAL     :: SPFLAST, RFMIN, FINTERP 

         LOGICAL  :: NEWFILE
         LOGICAL  :: WANTGPLOT 

         REAL, PARAMETER :: FSCCUTG = 0.143
         REAL, PARAMETER :: FSCCUTL = 0.50


         CALL SET_MPI(ICOMM,MYPID,MPIERR)

C        GET NUMBER OF OPERATION LINE REGISTERS SPECIFIED
         CALL REG_GET_USED(NSEL_USED)

         XPREVL = 0
         XPREVG = 0
         DLIST  = HUGE(FSCLAST)

C        OPEN OUTPUT DOC FILE
         CALL OPENDOC(DOCNAM,.TRUE.,NLET,LUNDOC,NICDOC,.TRUE.,
     &                'FSC OUTPUT DOCUMENT',.FALSE.,.FALSE.,.TRUE.,
     &                NEWFILE,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

C        OPEN FORMATTED, SEQUENTIAL FILE FOR GNUPLOT COMMANDS
         CALL OPAUXFILE(.TRUE.,GPLOTFILE,DATEXC,LUNGP,0,'N',
     &                       'GNUPLOT',.TRUE.,IRTFLG)
         NLET      = lnblnkn(GPLOTFILE)
         WANTGPLOT = (NLET > 0) .AND. (IRTFLG == 0)  
         IF (IRTFLG > 0) RETURN

C               1234567890123456789012345678901234567890 
         MSG = '       NORM-FREQ      ANGST         FSC'//
     &         '      |SQRT(FSC)| |SQRT(2FSC/(FSC+1))|'    
C                          10        20        30        40
         NMSG = 79
         CALL LUNDOCPUTCOM(LUNDOC,MSG(:NMSG),IRTFLG)

         FSCNOWG = 0.0
         FSCNOWL = 0.0
         FSCMIN  = HUGE(FSCMIN)      ! MIN OF FSC CURVE

         DO  L=1,INC
            NVOX = LR(L)

            IF (NVOX .NE. 0) THEN
               DLIST(1) = L         ! NUMBER

               SPFLAST  = DLIST(2)          ! WAS HUGE AT START
               DLIST(2) = FLOAT(L-1) / FLOAT(INC-1)*0.5  ! NORM FREQ

               IF (DLIST(2) == 0.0) THEN
C                 TRAP FOR DIV BY ZERO
                  DLIST(3) = 999999  ! ANGSTROMS
               ELSE
                  DLIST(3) = 0.5 / (FMAXSPFREQ * DLIST(2))  ! ANGSTROMS
               ENDIF

               RVAL(L)  = DLIST(3) 

               DLIST(5) = MIN(1.0, FACT/SQRT(FLOAT(NVOX))) ! FSCCRIT
               DLIST(6) = NVOX                             ! # VOXELS

               RFMIN    = -HUGE(RFMIN)
               NSCM     = 1

               DO NSC=1,NSCALE
		 IF (AMP(NSC,L) > TINY(RFMIN))  THEN
                    RFM = AVSUM(NSC,L) / MAX(1.0,AMP(NSC,L))
                    IF (RFM  <  RFMIN) THEN
                       NSCM     = NSC
                       RFMIN    = RFM
                    ENDIF
		 ENDIF
               ENDDO

C              NSCM IS THE NUMBER OF THE ELEMENT IN EACH ARRAY WITH THE
C              CORRECT SCALING. SCALE IS THE CORRECT SCALING.

               BK1     = AMP(NSCM,L)
               BK2     = PR(NSCM,L)
               BK3     = CSUM2(L)
               BK4     = CSUM1(L)

               FSCLAST = DLIST(4)  ! PREVIOUS FSC (STARTS AS HUGE)

               IF (BK3 > TINY(BK3) .AND. 
     &             BK4 > TINY(BK3)) THEN
		  DLIST(4) = CSUM(L) / SQRT(BK4 * BK3)   ! FSC
               ELSE
		  DLIST(4) = 0.0
               ENDIF
               FSCVAL(L)  = DLIST(4) 

C              WRITE FSC FILTER REGISTER COLUMNS
               FSCZ   = MAX(DLIST(4), 0.0)
               IF (FSCZ == 0.0) THEN
                  DLIST(5) = 0.0
		  DLIST(6) = 0.0
               ELSE
                  DLIST(5) = SQRT(FSCZ)
                  DLIST(6) = SQRT( 2*FSCZ/ (FSCZ+1) )
               ENDIF

C              WRITE TO DOC FILE
               NLEN = 5
               CALL LUNDOCWRTDAT(LUNDOC,L,DLIST(2),NLEN,IRTFLG)

               !IF (VERBOSE .AND. MYPID <= 0) THEN
               !   WRITE(NOUT,98) L,(DLIST(K),K=2,4)
               !  98  FORMAT (1X,I4,4(2X,F12.5),4X,I6)
               !ENDIF

               IF (DLIST(4) < FSCMIN) THEN
C                 MIN VALUE ON FSC CURVE
                  SPFMIN  = DLIST(2)
                  FSCMIN  = DLIST(4)
                  ILOCMIN = L
               ENDIF
               IF (L  >=  3 .AND. 
     &             FSCLAST  >= FSCCUTG .AND.
     &             DLIST(4) <  FSCCUTG .AND.
     &             FSCNOWG  <= 0  ) THEN

C                  CROSSED FSCCUTG FIRST TIME GOING DOWN
                   XPREVG   = L - 1     ! LAST INDEX ABOVE CUTOFF

                   FSCPREVG = FSCLAST
                   FSCNOWG  = DLIST(4) ! RESOLUTION FOR GOLD FSC 

                   SPFPREVG = SPFLAST
                   SPFNOWG  = DLIST(2)
               ENDIF

               IF (L  >=  3 .AND. 
     &             FSCLAST  >= FSCCUTL .AND.
     &             DLIST(4) <  FSCCUTL .AND.
     &             FSCNOWL  <= 0  ) THEN

C                  CROSSED FSCCUTL FIRST TIME GOING DOWN
                   XPREVL  = L - 1     ! LAST INDEX ABOVE CUTOFF

                   FSCPREVL = FSCLAST
                   FSCNOWL  = DLIST(4) ! RESOLUTION FOR LEGACY FSC 

                   SPFPREVL = SPFLAST
                   SPFNOWL  = DLIST(2)
               ENDIF

            ENDIF
         ENDDO

C        RESOLUTION NEEDED
         IF (XPREVG > 0) THEN
            FINTERP    = (FSCCUTG - FSCPREVG) / (FSCNOWG - FSCPREVG)  

            XINTERPG   = XPREVG   + FINTERP * (1)
            SPFINTERPG = SPFPREVG + FINTERP * (SPFNOWG - SPFPREVG)          
         ELSE
            XINTERPG   = ILOCMIN
            SPFINTERPG = SPFMIN    ! NOT INTERPOLATED?
         ENDIF

         RESOLG = 0.5 / (FMAXSPFREQ * SPFINTERPG )
         !write(6,*) 'resol:',resolg,fmaxspfreq,spfinterp,xinterpg 

         IF (XPREVL > 0) THEN
            FINTERP    = (FSCCUTL - FSCPREVL) / (FSCNOWL - FSCPREVL)  

            XINTERPL   = XPREVL   + FINTERP * (1)
            SPFINTERPL = SPFPREVL + FINTERP * (SPFNOWL - SPFPREVL)          
         ELSE
            XINTERPL   = ILOCMIN
            SPFINTERPL = SPFMIN    ! NOT INTERPOLATED?
         ENDIF

         RESOLL = 0.5 / (FMAXSPFREQ * SPFINTERPL )
         !write(6,*) 'resol:',resoll,fmaxspfreq,spfinterpl,xinterpl 

         IF (NSEL_USED > 0) THEN
C           OUTPUT TO SPIDER'S REGISTERS NEEDED
            CALL REG_SET_NSEL(1,5,XINTERPG,SPFINTERPG,RESOLG,
     &                            XINTERPL,SPFINTERPL,IRTFLG)
            CALL REG_SET_NSEL(6,1,RESOLL,0.0,0.0,0.0,0.0,IRTFLG)
         ENDIF

         IF (WANTGPLOT) THEN
C           WRITE GNUPLOT FILE OUTPUT

            WRITE(LUNGP,'(A)') 'set xlabel "Angstroms"' 
            WRITE(LUNGP,191)   'set title " At FSC: 0.143  Resolution:', 
     &                          RESOLG,   '   At FSC: 0.5  Resolution:',
     &                          RESOLL ,' (A)"'
191         FORMAT(A, F7.2, A,F7.2,A)

            WRITE(LUNGP,'(A)') 'set yrange [0:1.0]' 

            NPIXLEND2 = NINT(INC * 0.5)

            WRITE(LUNGP,192) 'set xrange [0:',NPIXLEND2,'] reverse' 
192         FORMAT(A, I6, A)

            WRITE(LUNGP,'(A,F6.2,A,F6.2,A)') 'plot ', FSCCUTG,
     &                  ', ', FSCCUTL,  
     &                  ', "-" using 2:3 with line'

            DO  L=1,INC
               IF (LR(L) .NE. 0) THEN   

                  WRITE(LUNGP,195) ,L,RVAL(L),FSCVAL(L)
195               FORMAT(' ',I5,'  ',F7.2,' ',F10.3)
               ENDIF
           ENDDO

         ENDIF

         END




