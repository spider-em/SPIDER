
C++*********************************************************************
C
C TRAFL.F
C           ADDED 'TF L FLIP'                      OCT 15 ArDean Leith                           
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2015  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@wadsworth.org                                        *
C=*                                                                    *
C=* SPIDER is free software! you can redistribute it and/or            *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation! either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* SPIDER is distributed in the hope that it will be useful,          *
C=* but WITHOUT ANY WARRANTY! without even the implied warranty of     *
C=* merchantability or fitness for a particular purpose.  See the GNU  *
C=* General Public License for more details.                           *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program. If not, see <http://www.gnu.org/licenses> *
C=*                                                                    *
C***********************************************************************
C
C  TRAFL
C
C  PURPOSE: GENERATE THE PHASE CONTRAST TRANSFER FUNCTION  FOR
C           BRIGHT-FIELD ELECTRON MICROSCOPY. THIS OPERATION WRITES 
C           THE 1-DIMENSIONAL TRANSFER FUNCTION (OR ITS SQUARE,
C           THE ENVELOPE FUNCTION) IN REAL, DISPLAYABLE FORM TO 
C           A DOCUMENT FILE.
C
C           ADDED 'TF L FLIP' WHICH CREATES PHASE FLIPPING DOC
C           FILE FOR USE BY 'FD' 
C
C                          
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

       SUBROUTINE TRAFL

       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC'

       CHARACTER             :: ANS
       REAL                  :: LAMBDA,KM

       INTEGER               :: NOT_USED,NX
       REAL                  :: CS

       CHARACTER(LEN=MAXNAM) :: DOCNAM
       REAL                  :: DLIST(3)
       CHARACTER(LEN=80)     :: COMMENT
          
       LOGICAL               :: ADDEXT,GETNAME,ISOLD
       LOGICAL               :: APPEND,MESSAGE,NEWFILE

       CHARACTER             :: NULL = CHAR(0)

       INTEGER,PARAMETER     :: LUNDOC = 81


       IF (FCHAR(6:6) == 'F') THEN    
          CALL TRAFL_FLIP
          RETURN
       ENDIF


       CALL RDPRM1S(CS,NOT_USED,
     &                  'SPHERICAL ABERRATION CS [MM]',IRTFLG)

       IF (CS < 0.0001) CS = 0.0001

       CALL RDPRM2S(DZ,LAMBDA,NOT_USED,
     &              'DEFOCUS [A], WAVELENGTH LAMBDA [A]',IRTFLG)

       CALL RDPRI1S(NX,NOT_USED,
     &              'NUMBER OF SPATIAL FREQUENCY POINTS',IRTFLG)

       CALL RDPRM1S(KM,NOT_USED,
     &              'MAX SPATIAL FREQUENCY [1/A]',IRTFLG)

       CALL RDPRM2S(Q,DS,NOT_USED,
     &              'SOURCE SIZE [1/A], DEFOCUS SPREAD [A]',IRTFLG)

       CALL RDPRM2S(WGH,ENV,NOT_USED,
     &'AMPL. CONTRAST RATIO [0-1], GAUSSIAN ENV. HALFW. [1/A]',
     &   IRTFLG)

       CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &       'DIFFRACTOGRAM, ENVELOPE, OR STRAIGHT (D/E/S)',
     &        NULL,IRTFLG)

       ENV = 1.0 / ENV**2
       SC  = KM / FLOAT(NX / 2)

C      OPEN OUTPUT DOC FILE
       ADDEXT  = .TRUE.
       GETNAME = .TRUE.
       ISOLD   = .FALSE.
       APPEND  = .FALSE.
       MESSAGE = .TRUE.
       IRTFLG  = -8         ! NO IC USE

       CALL OPENDOC(DOCNAM,ADDEXT,NLET,LUNDOC,LUNDOCNO,GETNAME,
     &           'OUTPUT DOC FILE',ISOLD,APPEND,MESSAGE,
     &            NEWFILE,IRTFLG)

C                123456789 123456789 123456789 123456789 123456789 123456789 
       COMMENT ='KEY=RAD. TRANSFER   RAD(PIX^-1}'    
       CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:34),IRTFLG)

       IE = 0
       IF (ANS == 'E') IE = 1

       WGH = ATAN(WGH / (1.0 - WGH))
       CS  = CS * 1.E7
       NS1 = NX / 2 + 1

       DO K=NS1,NX+1

          AK = (NS1 - K) * SC
          CALL TFD(B,CS,DZ,LAMBDA,Q,DS,IE,AK,WGH,ENV)

          IF (ANS .NE. 'S') B = B * B

          IKEY     = K - NS1 + 1
          DLIST(1) = B
          DLIST(2) = REAL(K - NS1) / NX

C         WRITE TO CTF DOC
          CALL LUNDOCWRTDAT(LUNDOCNO,IKEY,DLIST,2,IRTFLG)
       ENDDO

       CLOSE(LUNDOC)

       END



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



       SUBROUTINE TRAFL_FLIP

       IMPLICIT NONE
       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC'

       CHARACTER(LEN=1)      :: ANS
       CHARACTER(LEN=1)      :: NULL = CHAR(0)
       CHARACTER(LEN=100)    :: COMMENT
       REAL                  :: LAMBDA,KM

       REAL, ALLOCATABLE     :: DLIST(:,:)
       REAL                  :: DLIST1(3)
          
       CHARACTER(LEN=MAXNAM) :: DOCNAM
       INTEGER,PARAMETER     :: LUNDOC = 81

       LOGICAL               :: ADDEXT,GETNAME,ISOLD
       LOGICAL               :: APPEND,MESSAGE,NEWFILE
 
       INTEGER               :: NOT_USED,IRTFLG,NCHAR,NX,IE,NS1,IKEY
       INTEGER               :: NUM_BINS,K,IRAD,NLET,LUNDOCNO
       REAL                  :: CS,DZ,Q,DS,WGH
       REAL                  :: ENV,SP_PIXSIZ,SC,AK,B,RAD_PX,RAD_ANGST
       REAL                  :: PREV_CTF,FIRST_MIN
       REAL                  :: FIRSTMIN_RAD,FIRSTZERO_RAD,CURR_CTF
       REAL                  :: CTF_VALUE,STRAIGHT_CTF,TRAPPED_CTF
       REAL                  :: FLIPPED_CTF


       CALL RDPRM1S(CS,NOT_USED,
     &                  'SPHERICAL ABERRATION CS[MM]',IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

       IF (CS < 0.0001)  CS = 0.0001

       CALL RDPRM2S(DZ,LAMBDA,NOT_USED,
     &              'DEFOCUS[A], LAMBDA[A]',IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

       CALL RDPRI1S(NX,NOT_USED,
     &              'NUMBER OF SPATIAL FREQUENCY PTS',IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

       CALL RDPRM1S(KM,NOT_USED,
     &              'MAXIMUM SPATIAL FREQUENCY[A-1]',IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

       CALL RDPRM2S(Q,DS,NOT_USED,
     &              'SOURCE SIZE[A-1], DEFOCUS SPREAD[A]',IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

       CALL RDPRM2S(WGH,ENV,NOT_USED,
     &   'AMPL CONTRAST RATIO[0-1], GAUSSIAN ENV. HALFW.[FOU. UNITS]',
     &   IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

       CALL RDPRM1S(SP_PIXSIZ,NOT_USED, 'PIXEL SIZE[A]', IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

       CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &       'DIFFRACTOGRAM, ENVELOPE, OR STRAIGHT (D/E/S)',
     &        NULL,IRTFLG)

       ENV = 1.0 / ENV**2
       SC  = KM / FLOAT(NX / 2)
       IE  = 0
       IF (ANS == 'E') IE = 1

       WGH = ATAN(WGH / (1.0 - WGH))
       CS  = CS * 1.E7
       NS1 = NX / 2 + 1

C      GET NUMBER OF FOURIER BINS
       NUM_BINS = NX - NS1 + 2

       ALLOCATE (DLIST(6,NUM_BINS), STAT=IRTFLG)
       IF (IRTFLG .NE. 0) THEN
          CALL ERRT(46,'DLIST',6*NUM_BINS)
          RETURN
       ENDIF

       DO K=NS1,NX+1

          AK = (NS1-K) * SC
          CALL TFD(B,CS,DZ,LAMBDA,Q,DS,IE,AK,WGH,ENV)

          IF (ANS .NE. 'S') B = B * B     ! ENVELOPE

          IKEY          = K - NS1 + 1

          RAD_PX        = REAL(K-NS1) / NX   ! RAD_PX 
          DLIST(4,IKEY) = RAD_PX             ! RAD_PX

          RAD_ANGST     = RAD_PX / SP_PIXSIZ 
          DLIST(5,IKEY) = RAD_ANGST

          DLIST(6,IKEY) = B

       ENDDO

C      INITIALIZE FIRST MIN, ABSOLUTE MIN
       PREV_CTF      = DLIST(1,1)

       FIRST_MIN     = 1        ! RADIUS FOR FIRST MIN.
       FIRSTMIN_RAD  = -1       ! INITIALIZE FIRST MIN. RADIUS
       FIRSTZERO_RAD = -1       ! FIRST ZERO RADIUS

C      LOOP THROUGH FOURIER RADII TO FIND FIRST MINIMUM, FIRST ZERO
       DO IRAD=2,NUM_BINS       

          CURR_CTF  = DLIST(1,IRAD)
          RAD_PX    = DLIST(4,IRAD)         

C         CHECK FOR FIRST LOCAL MIN
          IF ( FIRSTMIN_RAD < 0 ) THEN
             IF ( CURR_CTF  > PREV_CTF ) THEN
                FIRST_MIN    = IRAD - 1     ! RADIUS TO END TRAP 
                FIRSTMIN_RAD = SP_PIXSIZ / RAD_PX
             ENDIF
          ENDIF

C         FIND FIRST ZERO
          IF ( FIRSTZERO_RAD < 0 ) THEN
C            LOOK FOR WHEN CTF CROSSES ORIGIN
             IF ( (CURR_CTF * PREV_CTF) <= 0 ) 
     &          FIRSTZERO_RAD = SP_PIXSIZ / RAD_PX
          ENDIF

          PREV_CTF = CURR_CTF    ! NEW, PREVIOUS CTF VALUE==CURRENT CTF 
       ENDDO                     ! END RADIUS-LOOP

C      OPEN OUTPUT DOC FILE
       ADDEXT  = .TRUE.
       GETNAME = .TRUE.
       ISOLD   = .FALSE.
       APPEND  = .FALSE.
       MESSAGE = .TRUE.
       IRTFLG  = -8         ! NO IC USE

       CALL OPENDOC(DOCNAM,ADDEXT,NLET,LUNDOC,LUNDOCNO,GETNAME,
     &           'OUTPUT DOC FILE',ISOLD,APPEND,MESSAGE,
     &            NEWFILE,IRTFLG)

C      WRITE RADII (IN ANGSTROMS) TO DOC FILE COMMENT KEY
C      NOTE: IT WOULD BE MORE ACCURATE TO INTERPOLATE, BI-LINEARLY 
C      PERHAPS, SO THESE VALUES WILL BE ON AVERAGE 1/2 FOURIER PIXEL OFF

C            123456789 123456789 123456789 123456789 123456789 123456789 
       COMMENT=
     &       '            DEFOCUS  RAD_FIRST_MIN(A)  RAD_FIRST_ZERO(A)'
       CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:60),IRTFLG)
       DLIST1(1) = DZ 
       DLIST1(2) = FIRSTMIN_RAD 
       DLIST1(3) = FIRSTZERO_RAD 
       CALL LUNDOCWRTDAT(LUNDOCNO,-999,DLIST1,3,IRTFLG)
       COMMENT = ' ------------------- '
       CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:20),IRTFLG)

C                123456789 123456789 123456789 123456789 123456789 123456789 
       COMMENT = '  TRANSFER '
       CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:16),IRTFLG)
       COMMENT = 'RAD    FLIPPED,     STRAIGHT,      TRAPPED, ' //
     &           '    RAD(PIX^-1)    RAD(A**-1)      RAW'    
       CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:94),IRTFLG)

C      LOOP THROUGH FOURIER RADII
       DO IRAD=1,NUM_BINS       

        ! GET ORIGINAL VALUES
        RAD_PX    = DLIST(4,IRAD)
        CTF_VALUE = DLIST(6,IRAD)    ! RAW VALUE

        ! STRAIGHT SIGN
        STRAIGHT_CTF = -CTF_VALUE    ! FOR UNTRAPPED CTF

        ! FLIP SIGN
        TRAPPED_CTF  = -CTF_VALUE    ! FOR TRAPPED CTF

        ! TRAP FOR LOW RESOLUTION FLIP SIGN
        IF ( IRAD < FIRST_MIN ) TRAPPED_CTF = 1

        IF ( STRAIGHT_CTF == 0  ) FLIPPED_CTF  = 0
        IF ( STRAIGHT_CTF .NE. 0 ) 
     &       FLIPPED_CTF = ABS( STRAIGHT_CTF ) / STRAIGHT_CTF 

        DLIST(1,IRAD) = FLIPPED_CTF 
        DLIST(2,IRAD) = STRAIGHT_CTF 
        DLIST(3,IRAD) = TRAPPED_CTF 
  
C       WRITE TO CTF DOC FILE
        CALL LUNDOCWRTDAT(LUNDOCNO,IRAD,DLIST(1,IRAD),6,IRTFLG)
      ENDDO

      CALL REG_SET_NSEL(1,2,FIRSTMIN_RAD,FIRSTZERO_RAD,
     &                  0.0,0.0,0.0,IRTFLG)

C     CLOSE DOCUMENT FILE
      CLOSE(LUNDOC)

      IF (ALLOCATED(DLIST))   DEALLOCATE(DLIST)

      END

      
