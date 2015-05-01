
C++*********************************************************************
C
C    AP_GETANGAS.F
C                          ADDED CCROT RETRIEVAL  FEB 05 ARDEAN LEITH
C                          '-'                    FEB 11 ARDEAN LEITH
C                          ADDED DIR. VECTORS     FEB 11 ARDEAN LEITH
C                          ADDED CC RETURN        APR 15 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2015  Health Research Inc.,                         *
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
C AP_GETANGAS(ILIST,NIMA,IMIT,ANGDOCNAM,IMGPAT, LUNIN,SA,ANGS,IRTFLG)
C
C PURPOSE:       READ PROJ. ANGLES AND CONVERT TO UNITARY 
C                DIRECTIONAL VECTORS.
C
C PARAMETERS:
C       ILIST               LIST OF IMAGE FILE NUMBERS        (INPUT)
C       NIMA                NUMBER OF IMAGES                  (INPUT)
C       IMIT                IMAGE NUMBER (0 IF ALL)           (INPUT)
C       ANGDOCNAM           ANGLES FILE NAME                  (INPUT)
C       IMGPAT              IMAGE SERIES FILE TEMPLATE        (INPUT)
C       LUNIMG              IMAGE FILE IO UNIT                (INPUT)
C                               < 0 IMAGE ALREADY OPEN
C       LUNANG              PROJ. ANGLE FILE IO UNIT          (INPUT)
C       NWANT               NUMBER OF PARAMS TO RETRIEVE      (INPUT)
C       ANGS                PROJ. ANGLES                      (OUTPUT)
C       GOTANGS             FLAG FOR SUCCESS READING ANGLES   (OUTPUT)
C       NGOTX               # OF ALIGNMENT PAR. RETURNED      (OUTPUT)                      (OUTPUT)
C       WANTDIRS            WANT DIRECTIONAL VECTORS          (INPUT)                      (OUTPUT)
C       DIRS                DIRECTIONAL VECTORS               (OUTPUT)                      (OUTPUT)
C       IRTFLG              ERROR FLAG                        (OUTPUT)
C                              IF == -8999 WANT MIRCC NOT ROTCC 
C
C--*********************************************************************

	SUBROUTINE AP_GETANGAS(ILIST,NIMA,IMIT,ANGDOCNAM,IMGPAT,
     &                         LUNIMG,LUNANG,NWANT,ANGS,GOTANGS,NGOTX,
     &                         WANTDIRS,DIRS,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'F90ALLOC.INC'

        INTEGER                 :: ILIST(NIMA)
        INTEGER                 :: NIMA,IMIT
        CHARACTER (LEN=*)       :: ANGDOCNAM
        CHARACTER (LEN=*)       :: IMGPAT
        INTEGER                 :: LUNIMG,LUNANG,NWANT 
        REAL                    :: ANGS(NWANT,NIMA)
        LOGICAL                 :: GOTANGS
        INTEGER                 :: NGOTX
        LOGICAL                 :: WANTDIRS
        REAL                    :: DIRS(3,NIMA)
        INTEGER                 :: IRTFLG

        CHARACTER (LEN=MAXNAM)  :: ANGDOCNAMPE
        CHARACTER (LEN=MAXNAM)  :: FILNAM
        REAL                    :: BUFIN(16)

C       DOC FILE POINTER
        REAL, POINTER           :: ANGBUF(:,:)

        LOGICAL                 :: ANGINHEADER,WANTCCMIR
        INTEGER                 :: IGO,IEND,NLET,MAXX,MAXY,IMI,IV
        INTEGER                 :: MAXIM,LSAM,LROW,NSLICE,KEY,NLIST
        INTEGER                 :: INTFLAG
        REAL                    :: ZT,CCMIR

        GOTANGS     = .FALSE.
        WANTCCMIR   = ( IRTFLG == -8999)    !WANT MIRCC NOT ROTCC 
        ANGINHEADER = (ANGDOCNAM(1:1) == '-')
        IRTFLG      = 0

        IF (ANGDOCNAM(1:1) == '*' .OR.
     &      ANGDOCNAM(1:1) == CHAR(0) ) THEN
C          NO ANGLES AVAILABLE IN DOC FILE OR HEADER
           IF (NWANT > 0) ANGS = 0.0    ! ARRAY ZERO
           IF (WANTDIRS)  DIRS = 0.0    ! ARRAY ZERO
           NGOTX  = 0
           IRTFLG = 0
           RETURN
        ENDIF

        NGOTX = 3

C       READ  ANGLES 
        IGO  = 1
        IEND = NIMA
        IF (IMIT > 0) THEN
C          ONLY WANT ONE IMAGE'S ANGLES
           IGO  = IMIT
           IEND = IMIT
        ENDIF

        IF (.NOT. ANGINHEADER) THEN
C          ANGLES ARE IN DOC FILE, GET THE FILE NAME
           CALL FILNAMANDEXT(ANGDOCNAM,DATEXC,ANGDOCNAMPE,NLET,
     &                        .TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          RETRIEVE ARRAY WITH ANGLES DATA IN IT, CAN NOT USE
C          LUNDOCGETANG BECAUSE OF ILIST!
           MAXX = 8 + 1
           IF (NWANT > 7) MAXX = 11 + 1      ! GET CCROT ALSO
           IF (WANTCCMIR) MAXX = 15 + 1      ! GET CCMIR 
           !write(6,*) ' nwant,maxx:',nwant,maxx


           MAXY = MAXVAL(ILIST(1:NIMA))
           CALL GETDOCDAT(' ',.FALSE.,ANGDOCNAMPE,LUNANG,.FALSE.,
     &                 MAXX,MAXY,ANGBUF,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

        ENDIF

	DO  IMI=IGO,IEND
           IV = IMI - IGO + 1
           IF (IMIT > 0) IV = 1

           IF (ANGINHEADER) THEN
C             ANGLES IN IMAGE HEADER

              MAXX = NWANT + 1
              IF (NWANT >= 8) MAXX = 11 +1

              IF (LUNIMG > 0) THEN
C                MUST OPEN IMAGE ON LUNIMG FIRST
                 NLET = 0
                 CALL FILGET(IMGPAT,FILNAM,NLET,ILIST(IMI),INTFLAG)
	         MAXIM = 0
	         CALL OPFILEC(0,.FALSE.,FILNAM,LUNIMG,'O',IFORM,
     &                  LSAM,LROW,NSLICE,MAXIM,' ',.FALSE.,IRTFLG)
                 IF (IRTFLG .NE. 0)  GOTO 9999

                 CALL LUNGETVALS(LUNIMG,IAPLOC + 1,MAXX,BUFIN,IRTFLG)
                 CLOSE(LUNIMG)

              ELSEIF (LUNIMG < 0) THEN
C                IMAGE ALREADY OPEN ON LUNIMG, DO NOT CLOSE
                 CALL LUNGETVALS(-LUNIMG,IAPLOC + 1,MAXX,BUFIN,IRTFLG)
              ENDIF
              IF (IRTFLG .NE. 0)  GOTO 9999

              IF (BUFIN(4) <= 0) THEN
	         CALL ERRT(102,
     &              'NO ANGLES IN HEADER OF IMAGE NUMBER',IMI)
                 ANGS(:,IV) = 0.0
                 IRTFLG     = 1
	         GOTO 9999
              ENDIF

C             GET PROJECTION ANGLES FROM HEADER BUFFER
              ANGS(1:3,IV) = BUFIN(1:3)

              IF (NWANT >= 7) THEN
C                WANT OTHER ALIGNMENT PARAMETERS ALSO
                 ANGS(4:7,IV) = BUFIN(6:9)
                 ZT           = MAXVAL(BUFIN(6:9))
                 IF (ZT > 0.0) NGOTX = 7
              ENDIF
              IF (NWANT == 8) THEN
C                WANT CCROT PARAMETER ALSO
                 ANGS(8,IV) = BUFIN(11)
                 NGOTX      = 8
              ENDIF

           ELSE
C             READ ANGLES FROM DOC FILE ON LUNANG

C             EXTRACT ALIGN. PARAM. FROM ANGBUF
              KEY   = ILIST(IMI)
              NLIST = 3
              IF (NWANT >= 7) NLIST = 8
              IF (NWANT == 8) NLIST = 11
              IF (WANTCCMIR)  NLIST = 17
              CALL LUNDOCGETKEY(LUNANG,ANGBUF(1,1),MAXX,MAXY,KEY,
     &                          BUFIN,NLIST,.TRUE.,IRTFLG)
	      IF (IRTFLG .NE. 0) THEN
	         CALL ERRT(102,'MISSING ANGLE FOR IMAGE',KEY)
                 GOTO 9999
              ENDIF

              ANGS(1:3,IV) = BUFIN(1:3)
              IF (NWANT >= 7) THEN
C                WANT OTHER ALIGNMENT PARAMETERS ALSO
                 ANGS(4:6,IV) = BUFIN(6:8)
                 ZT           = MAXVAL(BUFIN(6:8))
                 IF (ZT > 0.0) NGOTX = 7
                 ANGS(7,IV)   = 0.0
                 IF (BUFIN(4) < 0) ANGS(7,IV) = 1.0
              ENDIF
              IF (NWANT == 8) THEN
C                WANT CCROT PARAMETER ALSO
                 ANGS(8,IV) = BUFIN(11)
                 NGOTX      = 8
                 IF (WANTCCMIR) THEN
C                   WANT MIR-CC PARAMETER INSTEAD OF CCROT
                    ANGS(8,IV) = BUFIN(15)
                 ENDIF
              ENDIF
           ENDIF
	ENDDO

        IRTFLG = 0

        IF (WANTDIRS) THEN
C          CONVERT ANGLES TO UNITARY DIRECTIONAL VECTORS
	   CALL AP_GETSATA(ANGS,DIRS,NWANT,NIMA,IRTFLG)
        ENDIF

        GOTANGS = .TRUE.

9999    IF (IMIT <= 0) CLOSE(LUNANG)

C       DEALLOCATE DOC. FILE MEMORY
        IF (ASSOCIATED(ANGBUF)) DEALLOCATE(ANGBUF)

        END 

C       **************************** AP_GETSATA  ***********************

	SUBROUTINE AP_GETSATA(ANGIN,SATA,NCOL,NDIM,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        INTEGER         :: NCOL,NDIM,IRTFLG

        REAL            :: ANGIN(NCOL,NDIM)
        REAL            :: SATA(3,NDIM)

	REAL, PARAMETER :: QUADPI = 3.1415926535897932384626
	REAL, PARAMETER :: DGR_TO_RAD = (QUADPI/180)

        INTEGER         :: I

        DO I=1,NDIM
           SATA(1,I) = COS(ANGIN(3,I)*DGR_TO_RAD) * 
     &                 SIN(ANGIN(2,I)*DGR_TO_RAD)

           SATA(2,I) = SIN(ANGIN(3,I)*DGR_TO_RAD) * 
     &                 SIN(ANGIN(2,I)*DGR_TO_RAD)

           SATA(3,I) = COS(ANGIN(2,I)*DGR_TO_RAD)
        ENDDO
        IRTFLG = 0

        END


#ifdef NEVER
C       **************************** AP_TEST  ***********************

	SUBROUTINE AP_TEST(ARRAY,nsam,nrow,IRTFLG)

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        REAL, ALLOCATABLE     :: ARRAY(:,:)
        REAL, POINTER     :: ARRAYT(:,:)

        write(6,*) ' nsam:',nsam,nrow
        ALLOCATE(ARRAY(nsam,nrow), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            MWANT = nsam* nrow 
            CALL ERRT(46,'; ARRAY...',MWANT)
            GOTO 999
        ENDIF 
        write(6,*) ' setting:'
        array = 2.45

        write(6,*) ' array in tst:',array(1,5)

        !arrayt =>  array
        !write(6,*) ' arrayt in tst:',arrayt 

        IRTFLG = 0

 999    return

        END
#endif
