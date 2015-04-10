
C++*********************************************************************
C
C POLAR_CC.F
C             OPFILEC                              FEB 03 ARDEAN LEITH
C             DE-OMP ||   LOOP FOR FFTW3 USE       MAR 08 ARDEAN LEITH
C             COSMETIC & ERROR TRAP REWRITE        DEC 12 ARDEAN LEITH
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

        SUBROUTINE POLAR_CC

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        CHARACTER(LEN=MAXNAM) :: FILNAM

        REAL, ALLOCATABLE     :: X(:,:), REF(:,:)
        CHARACTER(LEN=1)      :: MODE,ASK
        DOUBLE PRECISION      :: T7(-3:3)

        CHARACTER(LEN=1)      :: NULL = CHAR(0)

        INTEGER, PARAMETER    :: INPIC  = 20
        INTEGER, PARAMETER    :: INREF  = 21
        INTEGER, PARAMETER    :: LUNO   = 22

C       ASK FOR DATA FILES
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,INPIC,'O',ITYPE,NX,NY,
     &             NZ,MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,INREF,'O',ITYPE,NXP,NYP,
     &              NZ,MAXIM,'REFERENCE',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (NX .NE. NXP .OR. NY .NE .NYP)  THEN
           CALL ERRT(101,'IMAGES MUST HAVE SAME DIMENSIONS',IER)
           RETURN

        ELSEIF (NZ > 1)  THEN
           CALL ERRT(101,'DOES NOT WORK ON VOLUMES',IER)
           RETURN
        ENDIF

C       ASK FOR OUTPUT FILE
        MAXIM  = 0
	ITYPE  = 1
	NYP    = 1
	NZP    = 1
        CALL OPFILEC(INPIC,.TRUE.,FILNAM,LUNO,'N',ITYPE,NX,NYP,
     &             NZP,MAXIM,'OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CALL  RDPRMC(ASK,NA,.TRUE.,
     &              'FULL OR HALF CIRCLE (F/H)',NULL,IRTFLG)
        IF (ASK == 'H')  THEN
           MODE = 'H'
        ELSE
           MODE = 'F'
        ENDIF

        NXLD = NX+2-MOD(NX,2)

        ALLOCATE (X  (NXLD,NY),
     &            REF(NXLD,NY), 
     &            STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'POLAR_CC, X,REF',2*NXLD*NY)
           GOTO 9999
        ENDIF

        DO J=1,NY
           CALL REDLIN(INPIC,X(1,J),NX,J)
        ENDDO

        DO J=1,NY
           CALL REDLIN(INREF,REF(1,J),NX,J)
        ENDDO

C       parallel do private(j),shared(invt) REMOVED FOR FFTW3 al

C       FORWARD FFT, LINE BY LINE
        DO J=1,NY
           INV = +1
           CALL FMRS_1(X(1,J),  NX,INV)

           INV = +1
           CALL FMRS_1(REF(1,J),NX,INV)
        ENDDO

C       CROSS CORRELATION IN FOURIER SPACE
        CALL MF1D(X,REF,NXLD/2,NY)

C       REVERSE FFT ON FIRST LINE
        INV = -1
        CALL FMRS_1(X(1,1),NX,INV)

C       FIND MAXIMUM & LOCATION IN FIRST LINE
        QM = -1.0E20
        DO I=1,NX
           IF (X(I,1) >= QM)  THEN
              QM   = X(I,1)
              JTOT = I
           ENDIF
        ENDDO

        !write(6,*) 'jtot:',jtot, x(jtot,1)

C       SUB-PIXEL MAX.
        DO K=-3,3
           J    = MOD(JTOT+K+NX-1,NX)+1
           T7(K)= X(J,1)
        ENDDO

C       FIT PARABOLA
        CALL PRB1D(T7,7,POS)

        RANG = REAL(JTOT) + POS
        IF (MODE == 'F')  THEN
           RANG = (RANG-1.0) / NX * 360.0
        ELSE
           RANG = (RANG-1.0) / NX * 180.0
        ENDIF 
  
        WRITE(NOUT,90)  RANG,QM
90      FORMAT('    Angle: ',F10.4,'    Peak: ',G12.5)

        CALL REG_SET_NSEL(1,2,RANG,QM,0.0,0.0,0.0,IRTFLG)

C       SHIFT 1D CCF TO HAVE THE ORIGIN AT INT(NX/2)+1
	DO I=1,NX
	  REF(MOD(I+INT(NX/2)-1,NX)+1, 1) = X(I,1)
	ENDDO

	CALL WRTLIN(LUNO,REF,NX,1)

9999    IF (ALLOCATED(X))   DEALLOCATE (X)
        IF (ALLOCATED(REF)) DEALLOCATE (REF)

        CLOSE(INPIC)
        CLOSE(INREF)
	CLOSE(LUNO)

        END


C       -------------------- MF1D -------------------------------------

C       CROSS CORRELATION IN FOURIER SPACE

        SUBROUTINE MF1D(X,REF,L,NY)

        COMPLEX  :: X(L,NY), REF(L,NY)

        X(:,1) = X(:,1) * CONJG(REF(:,1))

        IF (NY == 1)  RETURN

        DO J=2,NY
           X(:,1) = X(:,1) + X(:,J) * CONJG(REF(:,J))
        ENDDO

        END
