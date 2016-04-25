
C++*********************************************************************
C
C TRAFL.F
C           ADDED 'TF L FLIP'                      OCT 15 ArDean Leith                           
C           REWORKED 'TF L FLIP'                   NOV 15 ArDean Leith                           
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
C  TRAFL_LIS
C  PURPOSE:  'TF LIS' WHICH CREATES A DOC FILE CONTAINING: 
C             RADII, STRAIGHT CTF, NEGATIVE STRAIGHT CTF,
C             PHASE FLIPPING CTF, TRAPPED-CTF, CTF ENVELOPE FUNCTION, 
C             AND DIFFRACTOGRAM  COLUMNS.
C             THE DOC FILE CAN BE USED BY 'FD C'  FOR CTF CORRECTION
C             DO NOT USE FLIPPED CTF WITH 'FD C'
C                 
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

       SUBROUTINE TRAFL

       IMPLICIT NONE
       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC'

       CHARACTER             :: ANS
       REAL                  :: LAMBDA,FMAXSPFREQ

       INTEGER               :: NOT_USED,NX,NY
       REAL                  :: CS,FVAL,PIXSIZ
       REAL                  :: DZ,Q,DS,FDUM,ACR,GEH,SIGN,SC,AK,B
       INTEGER               :: NDIM,IRTFLG,NCHAR,NLET,LUNDOCNO,IE,NS1
       INTEGER               :: K,IKEY

       CHARACTER(LEN=MAXNAM) :: DOCNAM
       REAL                  :: DLIST(3)
       CHARACTER(LEN=80)     :: COMMENT
          
       LOGICAL               :: ADDEXT,GETNAME,ISOLD
       LOGICAL               :: APPEND,MESSAGE,NEWFILE
       LOGICAL               :: WANT_AST,WANT_GEH,WANT_SIGN
       LOGICAL               :: WANT_SPFREQ,WANT_PIXSIZ

       CHARACTER             :: NULL = CHAR(0)

       INTEGER, PARAMETER    :: LUNDOC = 81


       IF (FCHAR(5:5) == 'I') THEN    
          CALL TRAFL_LIS
          RETURN
       ENDIF

C      GET COMMON TF INPUTS
       NDIM        =  1         ! SQUARE ONLY
       WANT_AST    = .FALSE.    ! DO NOT ASK FOR ASTIG
       WANT_GEH    = .TRUE.     ! ASK FOR GEH
       WANT_SIGN   = .FALSE.    ! DO NOT ASK FOR SIGN
       WANT_SPFREQ = .TRUE.     ! ASK FOR SPFREQ
       WANT_PIXSIZ = .FALSE.    ! DO NOT ASK FOR PIXEL SIZE

       CALL GET_TF_INPUT(CS,DZ,LAMBDA,
     &                   NDIM, NX,NY,
     &                   WANT_SPFREQ,FMAXSPFREQ,
     &                   WANT_PIXSIZ,PIXSIZ,
     &                   Q,DS,
     &                   WANT_AST,  FDUM,FDUM,
     &                   WANT_GEH,  ACR,GEH,   
     &                   WANT_SIGN, FDUM,
     &                   IRTFLG)   

       CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &       'DIFFRACTOGRAM, ENVELOPE, OR STRAIGHT (D/E/S)',
     &        NULL,IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

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

C                   123456789 123456789 123456789 123456789 123456789 123456789 
       COMMENT =   '      TRANSFER'    

       CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:38),IRTFLG)

       IF (ANS == 'D') THEN
          COMMENT ='    DIFFRACTOGRAM    RAD(1/PIX}'
       ELSEIF  (ANS == 'C') THEN    
          COMMENT ='        ENVELOPE     RAD(1/PIX}'    
       ELSEIF  (ANS == 'S') THEN    
          COMMENT ='      STRAIGHT       RAD(1/PIX}'
       ENDIF    
       CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:34),IRTFLG)

       IE  = 0
       IF (ANS == 'E') IE = 1

       GEH = 1.0 / GEH**2
       SC  = FMAXSPFREQ / FLOAT(NX / 2)

       ACR = ATAN(ACR / (1.0 - ACR))
       CS  = CS * 1.E7
       NS1 = NX / 2 + 1

       DO K=NS1,NX+1

          AK = (NS1 - K) * SC
          CALL TFD(B,CS,DZ,LAMBDA,Q,DS,IE,AK,ACR,GEH)

          IF (ANS .NE. 'S') B = B * B   ! ENVELOPE OR DIFFRACTOGRAM

          IKEY     = K - NS1 + 1
          DLIST(1) = B
          DLIST(2) = REAL(K - NS1) / NX

C         WRITE TO CTF DOC
          CALL LUNDOCWRTDAT(LUNDOCNO,IKEY,DLIST,2,IRTFLG)
       ENDDO

       CLOSE(LUNDOC)

       END


C      ------------------- TRAFL_LIS -----------------------------


       SUBROUTINE TRAFL_LIS

       IMPLICIT NONE
       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC'

       CHARACTER(LEN=1)      :: ANS
       CHARACTER(LEN=1)      :: NULL = CHAR(0)
       CHARACTER(LEN=130)    :: COMMENT
       REAL                  :: LAMBDA,FMAXSPFREQ

       REAL, ALLOCATABLE     :: DLIST(:,:)
       REAL                  :: DLIST1(3)
          
       CHARACTER(LEN=MAXNAM) :: DOCNAM
       INTEGER,PARAMETER     :: LUNDOC = 81

       LOGICAL               :: ADDEXT,GETNAME,ISOLD
       LOGICAL               :: APPEND,MESSAGE,NEWFILE
       LOGICAL               :: WANT_AST,WANT_GEH,WANT_SIGN
       LOGICAL               :: WANT_SPFREQ,WANT_PIXSIZ
 
       INTEGER               :: NOT_USED,IRTFLG,NCHAR,NX,NS1,IKEY
       INTEGER               :: NUM_BINS,K,IRAD,NLET,LUNDOCNO,NDIM,NY
       REAL                  :: CS,DZ,Q,DS,ACR, B_ENV,RAD_ANGST
       REAL                  :: GEH,SP_PIXSIZ,SC,AK,B_STR,RAD_PX
       REAL                  :: PREV_CTF,FIRST_MIN,FIRST_ZERO
       REAL                  :: FIRSTMIN_RAD,FIRSTZERO_RAD,CURR_CTF
       REAL                  :: STRAIGHT_CTF,STRAIGHT_CTF_NEG
       REAL                  :: TRAPPED_CTF,FLIPPED_CTF,FLIPPED_CTF_TRI 
       REAL                  :: FDUM,SIGN,FVAL,FKEV,PIXSIZ

C      OPEN OUTPUT DOC FILE  (first for campatibility with 'TF C' ops
       ADDEXT  = .TRUE.
       GETNAME = .TRUE.
       ISOLD   = .FALSE.
       APPEND  = .FALSE.
       MESSAGE = .TRUE.
       IRTFLG  = -8         ! NO IC USE

       CALL OPENDOC(DOCNAM,ADDEXT,NLET,LUNDOC,LUNDOCNO,GETNAME,
     &           'OUTPUT DOC FILE',ISOLD,APPEND,MESSAGE,
     &            NEWFILE,IRTFLG)
       IF (IRTFLG .NE. 0) RETURN
C      GET COMMON TF INPUTS
       NDIM        =  1         ! SQUARE ONLY
       WANT_AST    = .FALSE.    ! DO NOT ASK FOR ASTIG
       WANT_GEH    = .TRUE.     ! ASK FOR GEH
       WANT_SIGN   = .FALSE.    ! DO NOT ASK FOR SIGN
       WANT_SPFREQ = .FALSE.    ! DO NOT ASK FOR SPFREQ
       WANT_PIXSIZ = .TRUE.     ! ASK FOR PIXEL SIZE

       CALL GET_TF_INPUT(CS,DZ,LAMBDA,
     &                NDIM, NX, NY,
     &                WANT_SPFREQ,FMAXSPFREQ,
     &                WANT_PIXSIZ,PIXSIZ,
     &                Q, DS,
     &                WANT_AST, FDUM, FDUM,
     &                WANT_GEH, ACR, GEH,
     &                WANT_SIGN, FDUM,
     &                IRTFLG) 


       IF (IRTFLG .NE. 0) RETURN

C      GET NUMBER OF FOURIER BINS
       NS1 = NX / 2 + 1
       NUM_BINS = NX - NS1 + 2

       ALLOCATE (DLIST(9,NUM_BINS), STAT=IRTFLG)
       IF (IRTFLG .NE. 0) THEN
          CALL ERRT(46,'DLIST',9*NUM_BINS)
          RETURN
       ENDIF

       IF (GEH .NE. 0) GEH = 1.0 / GEH**2
       SC  = FMAXSPFREQ / FLOAT(NX / 2)

       ACR = ATAN(ACR / (1.0 - ACR))
       CS  = CS * 1.E7

       DO K=NS1,NX+1

          AK = (NS1 - K) * SC
          CALL TFD_PLUS(B_STR,CS,DZ,LAMBDA,Q,DS,AK,ACR,GEH, B_ENV)

          IKEY          = K - NS1 + 1

          RAD_PX        = REAL(K-NS1) / NX   ! RAD_PX 
          DLIST(1,IKEY) = RAD_PX             ! RAD_PX

          RAD_ANGST     = RAD_PX / PIXSIZ 
          DLIST(2,IKEY) = RAD_ANGST

          DLIST(3,IKEY) = B_STR           ! STRAIGHT = RAW

          DLIST(7,IKEY) = B_ENV * B_ENV   ! ENVELOPE 

          DLIST(8,IKEY) = B_STR * B_STR   ! DIFFRACTOGRAM

       ENDDO

       !write(6,*) ' filled bins: ',ns1,'...',nx+1,'  numbins:',NUM_BINS

C      INITIALIZE FIRST MIN, ABSOLUTE MIN  FROM STRAIGHT CTF
       PREV_CTF   = DLIST(3,1)

       FIRST_MIN  = -1       ! RADIUS FOR FIRST MIN (PIXELS)
       FIRST_ZERO = -1       ! FIRST ZERO RADIUS (A)

C      LOOP THROUGH FOURIER RADII TO FIND FIRST MINIMUM, FIRST ZERO
       DO IRAD=2,NUM_BINS       

          CURR_CTF  = DLIST(3,IRAD)
          RAD_PX    = DLIST(1,IRAD)         

C         CHECK FOR FIRST LOCAL MIN
          IF ( FIRST_MIN < 0 ) THEN
             IF ( CURR_CTF  > PREV_CTF ) THEN
                FIRST_MIN    = IRAD - 1     ! RADIUS TO END TRAP 
                FIRSTMIN_RAD = PIXSIZ / RAD_PX
             ENDIF
          ENDIF

C         FIND FIRST ZERO
          IF ( FIRST_ZERO  < 0 ) THEN
C            LOOK FOR WHEN CTF CROSSES ORIGIN
             IF ( (CURR_CTF * PREV_CTF) <= 0 ) THEN
                FIRST_ZERO    = IRAD
                FIRSTZERO_RAD = PIXSIZ / RAD_PX
             ENDIF
          ENDIF

          PREV_CTF = CURR_CTF    ! NEW, PREVIOUS CTF VALUE==CURRENT CTF 
       ENDDO                     ! END RADIUS-LOOP


C      WRITE RADII (IN ANGSTROMS) TO DOC FILE COMMENT KEY
C      NOTE: IT WOULD BE MORE ACCURATE TO INTERPOLATE, BI-LINEARLY 
C      PERHAPS, SO THESE VALUES WILL BE ON AVERAGE 1/2 FOURIER PIXEL OFF

C            123456789 123456789 123456789 123456789 123456789 123456789 
       COMMENT=
     &   '            DEFOCUS  RAD_FIRST_MIN(PIX)  RAD_FIRST_ZERO(PIX)'
       CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:60),IRTFLG)
       DLIST1(1) = DZ 
       DLIST1(2) = FIRST_MIN 
       DLIST1(3) = FIRST_ZERO 
       CALL LUNDOCWRTDAT(LUNDOCNO,-998,DLIST1,3,IRTFLG)

       COMMENT = ' ------------------- '
       CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:20),IRTFLG)
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
       COMMENT = '  TRANSFER: '
       CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:16),IRTFLG)

       COMMENT = '       RAD(1/PIX),    RAD(1/A),  STRAIGHT,' // 
     &           '     -STRAIGHT,      FLIPPED,      TRAPPED,'   //   
     &           '     ENVELOPE,    DIFFRACTOGRAM'

       CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:126),IRTFLG)

C      LOOP THROUGH FOURIER RADII
       DO IRAD=1,NUM_BINS       

        ! GET ORIGINAL VALUES
        RAD_PX       = DLIST(1,IRAD)
        STRAIGHT_CTF = DLIST(3,IRAD)        ! STRAIGHT = RAW CTF VALUE

        ! STRAIGHT_CTF_NEG HAS NEGATIVE SIGN FOR UNDERFOCUS CONTRAST REVERSAL?
        STRAIGHT_CTF_NEG = -STRAIGHT_CTF     

        IF ( IRAD < FIRST_MIN ) THEN
           TRAPPED_CTF = 1
        ELSE
           TRAPPED_CTF = STRAIGHT_CTF_NEG  ! FOR TRAPPED CTF  
        ENDIF

C       FOR PHASE FLIPPING, GIVES BINARY CTF
        IF (STRAIGHT_CTF_NEG >= 0.0) THEN
            FLIPPED_CTF =  1.0
        ELSE
            FLIPPED_CTF = -1.0
        ENDIF

!C      FOR PHASE FLIPPING, GIVES TRINARY CTF
!       IF ( STRAIGHT_CTF_NEG == 0 ) THEN
!          FLIPPED_CTF_TRI = 0
!       ELSE 
!          FLIPPED_CTF_TRI = ABS( STRAIGHT_CTF_NEG ) / STRAIGHT_CTF_NEG 
!       ENDIF

       !DLIST(1,IRAD) = RAD(1/PIX)          ! ALREADY SET
       !DLIST(2,IRAD) = RAD(1/A)            ! ALREADY SET
       !DLIST(8,IRAD) = STRAIGHT_CTF        ! ALREADY SET
        DLIST(4,IRAD) = STRAIGHT_CTF_NEG 

        DLIST(5,IRAD) = FLIPPED_CTF 
        DLIST(6,IRAD) = TRAPPED_CTF 

       !DLIST(7,IRAD) = ENVELOPE            ! ALREADY SET
       !DLIST(8,IRAD) = DIFFRACTOGRAM       ! ALREADY SET 
       !DLIST(9,IRAD) = FLIPPED_CTF_TRI     ! ABANDONED 
  
C       WRITE TO CTF DOC FILE
        CALL LUNDOCWRTDAT(LUNDOCNO,IRAD,DLIST(1,IRAD),8,IRTFLG)

      ENDDO

      CALL REG_SET_NSEL(1,2,FIRSTMIN_RAD,FIRSTZERO_RAD,
     &                  0.0,0.0,0.0,IRTFLG)

C     CLOSE DOCUMENT FILE
      CLOSE(LUNDOC)

      IF (ALLOCATED(DLIST))   DEALLOCATE(DLIST)

      END

