
C++*********************************************************************
C
C FILTDOC.F     ALLOCATION                        MAY 2000 ARDEAN LEITH
C               OPENDOC PARAMETERS                DEC 2000 ARDEAN LEITH
C               LUNDOCREDSEQ RETURNS MAXY         APR 2003 ARDEAN LEITH
C               INCORE OPENDOC                    JUL 2003 ARDEAN LEITH
C               FOURIER INPUT AND OUTPUT FILE     OCT 2003 BIMAL RATH
C               NGOTY BUG                         FEB 2004 ARDEAN LEITH
C               REGISTER SELECTION                MAY 2024 ARDEAN LEITH
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
C  FILTDOC(LUNINT,LUNOUT,NX,NY,NZ,IFORMIN)
C
C  PURPOSE:      FILTERS AN IMAGE USING COEEFICIENTS FROM A
C                DOC. FILE
C
C  PARAMETERS:   LUNINT         INPUT UNIT FOR IMAGE & DOC     (SENT)
C                LUNOUT         OUTPUT UNIT FOR IMAGE          (SENT)
C                NX,NY & NZ                                    (SENT)
C
C--********************************************************************

	SUBROUTINE FILTDOC(LUNINT,LUNOUT,NX,NY,NZ,IFORMIN)

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 

        INTEGER               :: LUNINT,LUNOUT,NX,NY,NZ,IFORMIN

        COMPLEX, ALLOCATABLE  :: G(:,:,:)
        REAL, ALLOCATABLE     :: COEFFS(:,:)
        CHARACTER(LEN=MAXNAM) :: DOCFIL
        LOGICAL               :: LDUM
  
C       SET VALUE OF ISPACE AS PER REAL OR FT INPUT IMAGE
        IF (IFORMIN > 0) THEN
	   ISPACE = (NX + 2 - MOD(NX,2))/2
        ELSE
           ISPACE = NX
        ENDIF

	NZ = MAX(1,NZ)
        ALLOCATE (G(ISPACE,NY,NZ), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'FILTDOC, G',NDUM)
           RETURN
        ENDIF

C       LOAD IMAGE/VOLUME DATA (CAN NOT USE REDVOL HERE)
        DO K=1,NZ
           DO J=1,NY
              NR = J+(K-1)*NY
              CALL REDLIN(LUNINT,G(1,J,K),NX,NR)
	   ENDDO
        ENDDO
        CLOSE(LUNINT)

C       OPEN DOC FILE
        CALL OPENDOC(DOCFIL,.TRUE.,NLET,LUNINT,LUNIN,.TRUE.,'DOCUMENT',
     &            .TRUE.,.FALSE.,.TRUE.,LDUM,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        MAXX = 1
        IF (FCHAR(4:4) == 'C') THEN
C          READ DOC FILE REGISTER TO BE USED

           CALL RDPRI1S(MAXX,NOT_USED, 'DOC REGISTER NUMBER',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
        ENDIF

C       FOR 2D   (SQRT(2)/2) =~ .75
	LBS  = 0.75 * MAX(ISPACE,NY,NZ)

C       FOR 3D   (SQRT(3)/2) =~0.9
        IF (NZ > 1) LBS  = 0.9 * MAX(ISPACE,NY,NZ)

        MAXY = LBS

        ALLOCATE (COEFFS(MAXX,MAXY), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'FILTDOC, COEFFS',MAXY*MAXY)
           GOTO 9999
        ENDIF

        IF (IRTFLG .NE. 0) GOTO 9999

C       LOAD COEFFS DATA USING SEQUENTIAL READ FROM FIRST REGISTER
        CALL LUNDOCREDSEQ(LUNIN,COEFFS,MAXX,MAXY,MAXY,NGOTY,IRTFLG)
        CLOSE(LUNINT)
        IF (IRTFLG == 1) GOTO 9999

C       ZERO REMAINING COEFFS IN THE DESIRED COLUMN
	IF (NGOTY < LBS) COEFFS(MAXX, NGOTY+1:LBS) = 0.0

        !write(6,*) 'coeffs:',coeffs(maxx,:)

C       DO FT FOR REAL INPUT IMAGES
        IF (IFORMIN > 0) THEN
 	   INV = +1
 	   IF (NZ <= 1)  THEN
              CALL FMRS_2(G,NX,NY,INV)
 	   ELSE
              CALL FMRS_3(G,NX,NY,NZ,INV)
	   ENDIF
        ENDIF

c$omp   parallel do private(k,j,i,iz,iy,ir,rr,dd,fr)
	DO K=1,NZ
           IZ = K-1
           IF (IZ > NZ/2)  IZ = IZ - NZ

	   DO J=1,NY
              IY = J-1
              IF (IY > NY/2) IY = IY - NY

	      DO I=1,ISPACE
	         RR       = SQRT(FLOAT((I-1)**2 + IY*IY + IZ*IZ))
	         IR       = IFIX(RR)
	         DD       = RR - FLOAT(IR)
	         FR       = (1.0-DD) * COEFFS(MAXX,IR+1) +
     &                            DD * COEFFS(MAXX,IR+2)
	         G(I,J,K) = G(I,J,K) * FR	
	      ENDDO
	   ENDDO
	ENDDO

C       DO INVERSE FT FOR REAL INPUT IMAGES
        IF (IFORMIN > 0) THEN
	   INV = -1
 	   IF (NZ <= 1)  THEN
              CALL FMRS_2(G,NX,NY,INV)
 	   ELSE
              CALL FMRS_3(G,NX,NY,NZ,INV)
	   ENDIF
        ENDIF

        DO K=1,NZ
          DO J=1,NY
             NR = J+(K-1)*NY
             CALL WRTLIN(LUNOUT,G(1,J,K),NX,NR)
	  ENDDO
	ENDDO

C       CLOSE OUTPUT FILE
9999    CLOSE(LUNOUT)

C       DEALLOCATE RUN-TIME MEMORY
        IF (ALLOCATED(COEFFS)) DEALLOCATE(COEFFS)
        IF (ALLOCATED(G))      DEALLOCATE(G)

	END	
