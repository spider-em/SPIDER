C++*********************************************************************
C
C    MULTISHIFT
C                  OPFILEC                         FEB 03 ARDEAN LEITH
C                  PARALLEL FFTW BUG               JAN 08 ARDEAN LEITH
C                  REWRITE                         JAN 08 ARDEAN LEITH
C                  FMRS, AND OLD SUB PIXEL BUG     FEB 08 ARDEAN LEITH
C                  DOCALC ON FMRS CALL             FEB 08 ARDEAN LEITH
C
C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2001, P. A. Penczek                                   *
C=*                                                                    *
C=* University of Texas - Houston Medical School                       *
C=*                                                                    *
C=* Email:  pawel.a.penczek@uth.tmc.edu                                *
C=*                                                                    *
C=* This program is free software; you can redistribute it and/or      *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* This program is distributed in the hope that it will be useful,    *
C=* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C=* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
C=* General Public License for more details.                           *
C=*                                                                    *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program; if not, write to the                      *
C=* Free Software Foundation, Inc.,                                    *
C=* 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.      *
C=*                                                                    *
C=**********************************************************************
C                  
C  MULTISHIFT 
C
C  PURPOSE:  MULTI-REFERENCE SHIFT (WITH 180 DEGEES CHECK) ALIGNMENT
C            HAS BUG IN SUB_PIXEL RESULTS.
C  CALLS:
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

	SUBROUTINE MULTISHIFT

	INCLUDE 'CMBLOCK.INC' 
	INCLUDE 'CMLIMIT.INC' 

        CHARACTER (LEN=MAXNAM)               :: FINPAT,FINPIC,FILTOA

	INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NUMR
	INTEGER, ALLOCATABLE, DIMENSION(:)   :: ILIST,IRIST
	REAL, ALLOCATABLE, DIMENSION(:,:,:)  :: REFER_A

	DATA  INPIC/77/

C       ALLOCATE SPACE FOR REFERENCE IMAGE FILE LIST
	NILMAX = NIMAX
	ALLOCATE(ILIST(NILMAX),STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'AP MS, NUMR',NILMAX)
           RETURN
        ENDIF

C       ASK FOR REFERENCE IMAGE FILE LIST
	CALL FILELIST(.TRUE.,INPIC,FINPAT,NLETR,ILIST,NILMAX,NIMA,
     &      'TEMPLATE FOR REFERENCE IMAGE SERIES',IRTFLG)
	IF (IRTFLG .NE. 0) GOTO 9999

        IF (NIMA .GT. 0)  THEN
           WRITE(NOUT,2001) NIMA
2001       FORMAT('  Number of reference images: ',I6)
        ELSE
           CALL ERRT(101,'No reference images!',IDUM)
           GOTO 9999
        ENDIF

C       GET FIRST REFERENCE IMAGE TO DETERMINE DIMS
        CALL  FILGET(FINPAT,FINPIC,NLETR,ILIST(1),INTFLG)
        MAXIM = 0
	CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,NSAM,NROW,NSLICE,
     &               MAXIM,' ',.FALSE.,IRTFLG)
        IF (IRTFLG.NE.0)  GOTO 9999
        CLOSE(INPIC)

	CALL RDPRMI(NSIX,NSIY,NOT_USED,
     &             'TRANSLATION SEARCH RANGE X and Y')
	NSIX   = MAX(NSIX,1)
	NSIY   = MAX(NSIY,1)

C       DETERMINE NUMBER OF OMP THREADS
        CALL GETTHREADS(NUMTH)
        WRITE(NOUT,*) ' NUMBER OF OMP THREADS: ',NUMTH
 
C       ALLOCATE SPACE FOR IMAGES TO BE ALIGNED
	ALLOCATE(IRIST(NILMAX),STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'AP MS, NUMR',IER)
           GOTO 9999
        ENDIF

C       GET LIST OF SAMPLE IMAGES TO BE ALIGNED
C       NTOTAL IS NUMBER OF THE SAMPLE IMAGES
	CALL FILELIST(.TRUE.,INPIC,FILTOA,NLETI,IRIST,NILMAX,NTOTAL,
     &     'TEMPLATE FOR IMAGE SERIES TO BE ALIGNED',IRTFLG)
	IF (IRTFLG.NE.0) GOTO 9999

        IF (NTOTAL .GT. 0)  THEN
           WRITE(NOUT,2002) NTOTAL
2002       FORMAT('  Number of experimental images: ',I6/)
        ELSE
           CALL ERRT(101,'No experimental images!',IDUM)
           GOTO 9999
        ENDIF

C       NIMA IS NUMBER OF THE REFERENCE IMAGES
C       ALLOCATE SPACE FOR REFERENCE IMAGES
	ALLOCATE(REFER_A(2*(NSAM+1),2*NROW,NIMA),STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           MWANT = 2*(NSAM+1) * 2*NROW * NIMA
           CALL  ERRT(46,'AP MS, REFER_A',MWANT)
           RETURN
        ENDIF 

C       READ REFERENCE IMAGES
        CALL READIMS(FINPAT,NLETR,ILIST,NIMA,NSAM,NROW,REFER_A,IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           CALL ERRT(47,'AP MS',IER)
           RETURN
        ENDIF 
	
        IF (NTOTAL .GE. NUMTH)  THEN
C          WORKFLOW STRATEGY FOR LARGE NUMBER OF IMAGES TO BE ALIGNED
           write(nout,*) ' Using: mrshift_ps '

           CALL MRSHIFT_PS(ILIST,NIMA,IRIST,NTOTAL,NSAM,NROW,NSIX,NSIY,
     &                     NUMTH,REFER_A,FILTOA,NLETI)
	ELSE
C          FOR FEW IMAGES TO BE ALIGNED USE DIFFERENT STRATEGY.
           write(nout,*) ' Using: mrshift_ss '
           CALL MRSHIFT_SS(ILIST,NIMA,IRIST,NTOTAL,NSAM,NROW,NSIX,NSIY,
     &                     NUMTH,REFER_A,FILTOA,NLETI)
 	ENDIF

        WRITE (NOUT,2600)
2600    FORMAT (/,' ',7('-'),
     &      ' Multi-reference shift alignment, end of computation ',
     &      7('-'),/)

9999    IF (ALLOCATED(IRIST))      DEALLOCATE(IRIST)
        IF (ALLOCATED(ILIST))      DEALLOCATE(ILIST)
        IF (ALLOCATED(NUMR))       DEALLOCATE(NUMR)
        IF (ALLOCATED(REFER_A))    DEALLOCATE(REFER_A)

         END

C *************************** MRSHIFT_PS ****************************

        SUBROUTINE MRSHIFT_PS(ILIST,NIMA,IRIST,NTOTAL,
     &                NSAM,NROW,NSIX,NSIY,NUMTH,REFER_A,FILTOA,NLETI)

        USE TYPE_KINDS      

C       FFTW_PLANR IS A POINTER TO A STRUCTURE 
        INTEGER(KIND=I_8) :: FFTW_PLANF=0, FFTW_PLANR=0

	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'

	DIMENSION  REFER_A(2*(NSAM+1),2*NROW,NIMA)

	INTEGER, DIMENSION(NIMA)   :: ILIST
	INTEGER, DIMENSION(NTOTAL) :: IRIST
        CHARACTER (LEN=*)          :: FILTOA

C       AUTOMATIC ARRAYS
        PARAMETER (NDLI=7)
        DIMENSION  DLIST(NDLI,NUMTH)
	INTEGER    :: NASSIG(NUMTH)

	REAL, ALLOCATABLE, DIMENSION(:,:,:)    :: A
	COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: OUT

        CHARACTER(LEN=MAXNAM)   :: FINPIC 
        CHARACTER(LEN=1)        :: MODE

        DATA  NDOC/56/,INPIC/58/

	ALLOCATE(A(2*(NSAM+1),2*NROW,NUMTH),
     &		 OUT((NSAM+1),2*NROW,NUMTH),STAT=IRTFLG)
	IF (IRTFLG.NE.0) THEN
            MWANT = 3*(NSAM+1) * 2*NROW * NUMTH
            CALL ERRT(46,'AM MS, A and OUT ',MWANT)
            RETURN
        ENDIF 

C       NTOTAL INPUT IMAGES READY TO BE ALIGNED
        LSD       = 2 * (NSAM+1)
        IAX       = 2 * NSAM
        LSR       = 2 * NROW
        NUMTHFFTW = 1

        CALL FFTW3_MAKEPLAN(IAX,LSR,1,NUMTHFFTW,FFTW_PLANF, 1,IRTFLG)
        CALL FFTW3_MAKEPLAN(IAX,LSR,1,NUMTHFFTW,FFTW_PLANR,-1,IRTFLG)

C       LOOP OVER IMAGES TO BE ALIGNED
 	DO IMIT=1,NTOTAL,NUMTH
	   DO IMI=IMIT,MIN(NTOTAL,IMIT+NUMTH-1)

	      CALL FILGET(FILTOA,FINPIC,NLETI,IRIST(IMI),IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999

	      MAXIM = 0
	      CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,NSAMT,NROWT,
     &                     NSLICE, MAXIM,' ',.FALSE.,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999

              IF (NSAMT.NE.NSAM .OR. NROWT.NE.NROW)  THEN
                 CALL ERRT(1,'AP MS',NE)
                 CLOSE(INPIC)
	         GOTO 9999
	      ENDIF
	      DO J=1,NROW
	         CALL REDLIN(INPIC,A(1,J,IMI-IMIT+1),NSAM,J)
	      ENDDO
	      CLOSE(INPIC)
	   ENDDO

C          NUMTH INPUT IMAGES READY TO BE ALIGNED

c$omp      parallel do private(imi,irtflg) !!!!!!!!!!!!!!!
	   DO IMI=IMIT,MIN(NTOTAL,IMIT+NUMTH-1)
C             LOOP OVER ALL SAMPLE IMAGES

              CALL RAMP_C(A(1,1,IMI-IMIT+1),NSAM,NROW)      ! FLAT FIELDING
              CALL PD2(A(1,1,IMI-IMIT+1),LSD,LSR,NSAM,NROW) ! PAD 2X

              INV = +1                                  ! FORWARD FOURIER
              CALL FMRS(A(1,1,IMI-IMIT+1),IAX,LSR,1, FFTW_PLANF,
     &                  .FALSE.,.FALSE.,INV,IRTFLG)

              CALL SHFRA(A(1,1,IMI-IMIT+1),
     &           OUT(1,1,IMI-IMIT+1),
     &	         NSAM,NROW,
     &	         REFER_A,NIMA,NSIX,NSIY,
     &	         DLIST(3,IMI-IMIT+1),
     &	         DLIST(5,IMI-IMIT+1),
     &           DLIST(6,IMI-IMIT+1),
     &           NASSIG(IMI-IMIT+1),FFTW_PLANR)
	   ENDDO

C          OUTPUT (IN DLIST POSITION IS INCREASED BY 1, NO.1 IS THE KEY).
C          1 - NUMBER OF THE MOST SIMILAR REFERENCE PROJECTION
C          2 - CORRELATION COEFFICIENT
C          3 - SX
C          4 - SY
C          5 - IMAGE NUMBER

           DO  IMI=IMIT,MIN(NTOTAL,IMIT+NUMTH-1)
	      DLIST(2,IMI-IMIT+1) = ILIST(IABS(NASSIG(IMI-IMIT+1)))
	      IF (NASSIG(IMI-IMIT+1).GT.0)  THEN
	         DLIST(4,IMI-IMIT+1) = 0.0
	      ELSE
	         DLIST(4,IMI-IMIT+1) = 180.0
	      ENDIF

              DLIST(7,IMI-IMIT+1) = IRIST(IMI)
              DLIST(1,IMI-IMIT+1) = IMI
              CALL SAVD(NDOC,DLIST(1,IMI-IMIT+1),NDLI,IRTFLG)
	   ENDDO
	ENDDO

	CLOSE(NDOC)
        CALL SAVDC

C       DEALLOCATE LOCAL ARRAYS
9999    DEALLOCATE(A,OUT)
        CALL FFTW3_KILLPLAN(FFTW_PLANF,IRTFLG)
        CALL FFTW3_KILLPLAN(FFTW_PLANR,IRTFLG)

	END


C **********************************************************************
C
C  MRSHIFT_SS   : FOR SMALL NUMBER OF IMAGES TO BE ALIGNED
C
C--*********************************************************************

        SUBROUTINE MRSHIFT_SS(ILIST,NIMA,IRIST,NTOTAL,
     &                NSAM,NROW,NSIX,NSIY,NUMTH,REFER_A,FILTOA,NLETI)

        USE TYPE_KINDS      

C       FFTW_PLANR IS A POINTER TO A STRUCTURE 
        INTEGER(KIND=I_8) :: FFTW_PLANR=0

	INCLUDE 'CMBLOCK.INC'
	DIMENSION  REFER_A(2*(NSAM+1),2*NROW,NIMA)

	INTEGER, DIMENSION(NIMA)   :: ILIST
	INTEGER, DIMENSION(NTOTAL) :: IRIST
        CHARACTER (LEN=*)          :: FILTOA

C       AUTOMATIC ARRAYS
        PARAMETER (NDLI=7)
        DIMENSION  DLIST(NDLI)

	INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: NASSIG
	REAL,    ALLOCATABLE, DIMENSION(:,:,:) :: A,RESI
	COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: OUT

        CHARACTER(LEN=1)  :: MODE

        DATA  NDOC/56/,INPIC/58/

C       ALIGNMENT

	ALLOCATE(A(2*(NSAM+1),2*NROW,NTOTAL),
     &		 OUT((NSAM+1),2*NROW,NUMTH),
     &		 RESI(3,NIMA,NTOTAL),
     &		 NASSIG(NIMA,NTOTAL),
     &		 STAT=IRTFLG)

	IF (IRTFLG .NE. 0) THEN
            MWANT = 2*(NSAM+1) * 2*NROW * NTOTAL +
     &              1*(NSAM+1) * 2*NROW * NUMTH +
     &              4*NIMA*NTOTAL 

            CALL ERRT(46,'AM MS, A, RESI, NASSIG, and OUT',MWANT)
            RETURN
        ENDIF

C       READ ALL SAMPLE IMAGES TO BE ALIGNED
        CALL READIMS(FILTOA,NLETI,IRIST,NTOTAL,NSAM,NROW,A,IRTFLG)
	IF (IRTFLG .NE. 0) THEN
            CALL  ERRT(47,'AP MS',IER)
            GOTO 9999
        ENDIF

C       LOOP OVER SAMPLE IMAGES TO BE ALIGNED
C       NIMA   - NUMBER OF REFERENCE IMAGES
C       NTOTAL - NUMBER OF IMAGES TO BE ALIGNED

        NUMTHFFTW = 1
 	CALL FFTW3_MAKEPLAN(NSAM*2,NROW*2,1,NUMTHFFTW,
     &                      FFTW_PLANR,-1,IRTFLG)

 	DO  IMITT=1,NTOTAL*NIMA,NUMTH

c$omp      parallel do private(IMIT,NREF,IMI) !!!!!!!!!!!!!!!!
	   DO  IMIT=IMITT,MIN(NTOTAL*NIMA,IMITT+NUMTH-1)
		NREF = MOD(IMIT-1,NIMA)+1
		IMI  = (IMIT-1)/NIMA+1

		CALL SHFRS(A(1,1,IMI),OUT(1,1,IMIT-IMITT+1),
     &	           NSAM,NROW,
     &	           REFER_A(1,1,NREF),NSIX,NSIY,
     &	           RESI(1,NREF,IMI),
     &	           RESI(2,NREF,IMI),
     &             RESI(3,NREF,IMI),
     &             NASSIG(NREF,IMI),FFTW_PLANR)
	   ENDDO
	ENDDO

C       OUTPUT (IN DLIST POSITION IS INCREASED BY 1, NO.1 IS THE KEY).
C          1 - NUMBER OF THE MOST SIMILAR REFERENCE PROJECTION.
C          2 - CORRELATION COEFFICIENT.
C          3 - ANGLE (0 or 180)
C          4 - SX
C          5 - SY
C          6 - IMAGE NUMBER.

	DO IMI=1,NTOTAL
	  CM = -HUGE(CM)
	  DO NREF=1,NIMA
             IF (RESI(1,NREF,IMI) .GT. CM)  THEN
                CM  = RESI(1,NREF,IMI)
                NRI = NREF
             ENDIF
	  ENDDO

	  DLIST(1)=IMI
	  DLIST(2)   = ILIST(NRI)
	  DLIST(3)   = RESI(1,NRI,IMI)
	  DLIST(4)   = 180*NASSIG(NRI,IMI)
	  DLIST(5:6) = RESI(2:3,NRI,IMI)
	  DLIST(7)   = IRIST(IMI)
	  CALL SAVD(NDOC,DLIST,NDLI,IRTFLG)
	ENDDO

	CLOSE(NDOC)
        CALL  SAVDC

C       DEALLOCATE LOCAL ARRAYS
9999    DEALLOCATE(A,OUT,NASSIG,RESI)

        CALL FFTW3_KILLPLAN(FFTW_PLANR,IRTFLG)

	END


C       ------------------ SHFRS --------------------------------------

	SUBROUTINE SHFRS(A,OUT,NSAM,NROW,REFER_A,
     &		         NSIX,NSIY,CM,SX,SY,IDIM,FFTW_PLANR)

	COMPLEX :: A((NSAM+1),2*NROW), REFER_A((NSAM+1),2*NROW)
	COMPLEX :: OUT((NSAM+1),2*NROW)

	INV = -1
	LSD = 2 * (NSAM+1)
	IAX = 2 * NSAM
	LSR = 2 * NROW

        OUT = REFER_A * CONJG(A)         ! FOURIER MULTIPLICATION

        CALL FMRS(OUT,IAX,LSR,1, FFTW_PLANR, .FALSE.,.FALSE.,
     &            INV,IRTFLG)                 ! FOURIER BACK TRANSFORM

	CALL FINDRMX(OUT,LSD,LSR,NSIX,NSIY,CM,SX,SY)  ! MAX. LOCATION

	OUT = REFER_A * A                    ! FOURIER MULTIPLICATION

        CALL FMRS(OUT,IAX,LSR,1, FFTW_PLANR, .FALSE.,.FALSE.,
     &            INV,IRTFLG)                  ! FOURIER BACK TRANSFORM

	CALL FINDRMX(OUT,LSD,LSR,NSIX,NSIY,CMX,SXS,SYS) ! MAX. LOCATION

	IF (CMX .GT. CM)  THEN
	    CM   = CMX / FLOAT(IAX * LSR)  ! NO SCALING OF FFT
	    SX   = SXS
	    SY   = SYS
	    IDIM = 1
	ELSE
	    IDIM = 0
	ENDIF

	END


C       ------------------ SHFRA --------------------------------------

       SUBROUTINE SHFRA(A,OUT,NSAM,NROW,REFER_A,NIMA,
     &		        NSIX,NSIY,CM,SX,SY,IDIM,FFTW_PLANR)

       COMPLEX :: A((NSAM+1),2*NROW)
       COMPLEX :: REFER_A((NSAM+1),2*NROW,NIMA)
       COMPLEX :: OUT((NSAM+1),2*NROW)
       

       LSD = 2 * (NSAM+1)
       IAX = 2 * NSAM
       LSR = 2 * NROW
       CM  = -HUGE(SX)

       DO  IR=1,NIMA

         OUT = REFER_A(:,:,IR) * CONJG(A)     ! FOURIER MULTIPLICATION

         INV = -1                              ! BACK TRANSFORM
 	 CALL FMRS(OUT,IAX,LSR,1, FFTW_PLANR, .FALSE.,.FALSE.,
     &             INV,IRTFLG)                 ! FOURIER BACK TRANSFORM

         CALL FINDRMX(OUT,LSD,LSR,NSIX,NSIY,CMX,SXS,SYS) ! PEAK LOCATION

         IF (CMX .GT. CM)  THEN
	    CM   = CMX
	    SX   = SXS
	    SY   = SYS
	    IDIM = IR
         ENDIF

         OUT = REFER_A(:,:,IR) * A             ! FOURIER MULTIPLICATION

         INV = -1                              ! BACK TRANSFORM
 	 CALL FMRS(OUT,IAX,LSR,1, FFTW_PLANR, .FALSE.,.FALSE.,
     &             INV,IRTFLG)                 ! FOURIER BACK TRANSFORM

         CALL FINDRMX(OUT,LSD,LSR,NSIX,NSIY,CMX,SXS,SYS) ! PEAK LOCATION

         IF (CMX .GT. CM)  THEN
	    CM   = CMX
	    SX   = SXS
	    SY   = SYS
	    IDIM = -IR
         ENDIF
       ENDDO

       CM = CM / FLOAT(IAX * LSR) ! NO SCALING OF FFT

       END


C       ------------------ FINDRMX --------------------------------------

C        FIND  MAX. PEAK LOCATION

         SUBROUTINE  FINDRMX(D,LSD,NROW,NSIX,NSIY,CMX,SX,SY)

         DIMENSION  :: D(LSD,NROW),Z(-1:1,-1:1)
         LOGICAL    :: FOUND

         FOUND = .FALSE.
	 NSAM  = LSD - 2
	 NS21  = NSAM / 2+ 1
	 NR21  = NROW / 2+ 1
	
         CMX = -HUGE(SX)
         SX  = 0.0
         SY  = 0.0
         DO J=1,NSIY+1
            DO I=1,NSIX+1
               TMX = D(I,J) / FLOAT(NS21-I) / 
     &                        FLOAT(NR21-J)
               IF (CMX .LT. TMX)  THEN
		  CMX   = TMX
                  IX    = I-1
                  IY    = J-1
                  FOUND = .TRUE.
               ENDIF
            ENDDO         
         ENDDO
         DO J=1,NSIY+1
            DO I=NSAM-NSIX+1,NSAM
               TMX = D(I,J) / FLOAT(NS21-(NSAM-I+2)) /
     &                        FLOAT(NR21-J)
               IF ( CMX .LT. TMX)  THEN
		  CMX   = TMX
                  IX    = -(NSAM-I+1)
                  IY    = J-1
                  FOUND = .TRUE.
               ENDIF
            ENDDO         
         ENDDO

         DO J=NROW-NSIY+1,NROW
            DO I=1,NSIX+1
               TMX = D(I,J) / FLOAT(NS21-I) /
     &                        FLOAT(NR21-(NROW-J+2))
               IF (CMX .LT. TMX)  THEN
		  CMX   = TMX
                  IX    = I-1
                  IY    = -(NROW-J+1)
                  FOUND = .TRUE.
               ENDIF
            ENDDO         
         ENDDO
         DO J=NROW-NSIY+1,NROW
            DO I=NSAM-NSIX+1,NSAM
               TMX = D(I,J) / FLOAT(NS21-(NSAM-I+2)) /
     &                        FLOAT(NR21-(NROW-J+2))
               IF (CMX .LT. TMX)  THEN
		  CMX   = TMX
                  IX    = -(NSAM-I+1)
                  IY    = -(NROW-J+1)
                  FOUND = .TRUE.
               ENDIF
            ENDDO         
         ENDDO
	 
         IF (.NOT.FOUND)  RETURN
         SX = IX
         SY = IY

C       DO NOT INTERPOLATE IF THE MAX IS ON THE EDGE
	IF (NSIX.EQ.IABS(IX) .OR. NSIY.EQ.IABS(IY))  RETURN

	IF (IX .LE. 0)  THEN
	   IXS = NSAM + IX + 1
	ELSE
	   IXS = IX + 1
	ENDIF

	IF (IY .LE. 0)  THEN
	   IYS = NROW + IY + 1
	ELSE
	   IYS = IY + 1
	ENDIF

         DO J=-1,1
	   IYT = MOD(NROW+IYS+J-1,NROW) + 1
            DO I=-1,1
	      IXT    = MOD(NSAM+IXS+I-1, NSAM) + 1
              Z(I,J) = D(IXT,IYT) / FLOAT(NSAM/2-IABS(IX+I)) / 
     &                              FLOAT(NROW/2-IABS(IY+J))

            ENDDO
         ENDDO

         CALL PARABL(Z, XSH,YSH,CMX)

         SX = SX + XSH
         SY = SY + YSH

         END

C        -------------------- READIMS --------------------------------

        SUBROUTINE READIMS(FINPAT,NLET,ILIST,NIMA,NSAM,NROW,A,IRTFLG)

        USE TYPE_KINDS      

	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'

C       EXTERNAL ARRAYS
	DIMENSION  A(2*(NSAM+1),2*NROW,NIMA),ILIST(NIMA)
        CHARACTER(LEN=*)      :: FINPAT

        CHARACTER(LEN=MAXNAM) :: FINPIC

C       PLAN IS A POINTER TO A STRUCTURE 
        INTEGER(KIND=I_8) :: FFTW_PLANF=0

	DATA  INPIC/34/

C       READ ALL REFERENCE IMAGES
        DO  IMI=1,NIMA
           CALL FILGET(FINPAT,FINPIC,NLET,ILIST(IMI),INTFLAG)

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,NSAMT,NROWT,
     &                  NSLICE, MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0)  RETURN

           DO J=1,NROW
              CALL REDLIN(INPIC,A(1,J,IMI),NSAM,J)
           ENDDO
           CLOSE(INPIC)
	ENDDO

	LSD       = 2 * (NSAM+1)
        IAX       = 2 * NSAM
	LSR       = 2 * NROW
	INV       = +1
        NUMTHFFTW = 1
 
	CALL FFTW3_MAKEPLAN(IAX,LSR,1,NUMTHFFTW,FFTW_PLANF,INV,IRTFLG)

c$omp   parallel do private(IMI,IRTFLG) !!!!!!!!!!!!!!!!!!
	DO  IMI=1,NIMA
           CALL RAMP_C(A(1,1,IMI),NSAM,NROW)        ! FLAT FIELDING
           CALL PD2(A(1,1,IMI),LSD,LSR,NSAM,NROW)   ! PAD 2X

	   CALL FMRS(A(1,1,IMI),IAX,LSR,1, FFTW_PLANF,
     &              .FALSE.,.FALSE., INV,IRTFLG)
	ENDDO

	CALL FFTW3_KILLPLAN(FFTW_PLANF,IRTFLG)

        IRTFLG = 0
	IF (INV .GT. 1)  THEN
C          SIGNAL INCORRECT FFT RESULT
           IRTFLG = 1
	ENDIF

	END


C        -------------------- RAMP_C --------------------------------

	 SUBROUTINE  RAMP_C(X,NSAM,NROW)

C        PURPOSE: SUBTRACT THE RAMP AND NORMALIZE THE IMAGE

         DIMENSION        X(2*(NSAM+1),2*NROW)
         EXTERNAL         BETAI
         DOUBLE PRECISION BETAI
         DOUBLE PRECISION C,D,EPS,B1,B2,A,F,R2,DN1,DN2
         DOUBLE PRECISION QYX1,QYX2,QX1X2
     &                   ,QX1,QX2,QY,SYX1,SYX2,SX1X2,SX1
     &                   ,SX2,SY,SX1Q,SX2Q,SYQ

         PARAMETER (EPS=1.0D-5)
	 	 
         SYX1 = 0.0
         SYX2 = 0.0
         SY   = 0.0
         SX1Q = 0.0
         SX2Q = 0.0
         SYQ  = 0.0

         N1   = NSAM / 2
         N2   = NROW / 2
         SX1  = FLOAT(N1) * FLOAT(NSAM + 1)
         IF (MOD(NSAM,2) .EQ. 1)   SX1 = SX1 + 1 + N1
         SX2  = FLOAT(N2) * FLOAT(NROW + 1)
         IF (MOD(NROW,2) .EQ. 1)   SX2 = SX2 + 1 + N2
         SX1   = SX1 * NROW
         SX2   = SX2 * NSAM
         SX1X2 = 0.0D0

         DO  J = 1, NROW
           DO I = 1, NSAM
             SYX1 = SYX1 + X(I,J) * I
             SYX2 = SYX2 + X(I,J) * J
             SY   = SY   + X(I,J)
             SX1Q = SX1Q + I * I
             SX2Q = SX2Q + J * J
             SYQ  = SYQ  + X(I,J) * DBLE(X(I,J))
           END DO
         END DO

         DN    = FLOAT(NSAM) * FLOAT(NROW)
         QYX1  = SYX1 - SX1 * SY / DN
         QYX2  = SYX2 - SX2 * SY / DN
         QX1X2 = 0.0
         QX1   = SX1Q - SX1 * SX1 / DN
         QX2   = SX2Q - SX2 * SX2 / DN
         QY    = SYQ  - SY  * SY  / DN
         C     = QX1  * QX2 - QX1X2 * QX1X2

         IF (C .GT. EPS) THEN
           B1  = (QYX1 * QX2 - QYX2 * QX1X2) / C
           B2  = (QYX2 * QX1 - QYX1 * QX1X2) / C
           A   = (SY - B1 * SX1 - B2 * SX2)  / DN

           D   = A + B1 + B2
           SY  = 0.0
           SYQ = 0.0
           DO I = 1, NROW
             QY = D
             DO K = 1, NSAM
                X(K,I) = X(K,I) - QY
                SY     = SY   + X(K,I)
                SYQ    = SYQ  + X(K,I) * DBLE(X(K,I))
                QY     = QY + B1
             ENDDO
             D = D + B2
           ENDDO

C         ELSE
C           WRITE(NOUT,3030)
C3030       FORMAT(/,' No solution - image is not modified !')
         ENDIF

	 SY   = SY / NSAM / NROW
	 SYQ  = 1.0 / DSQRT((SYQ-NSAM*NROW*SY*SY )/ (NSAM*NROW-1))
	 X(1:NSAM,1:NROW) = SYQ * (X(1:NSAM,1:NROW) - SY)
         END

C       --------------------------- PD2 ------------------------------

	SUBROUTINE  PD2(A,LSD,LSR,NSAM,NROW)
	
	DIMENSION  A(LSD,LSR)

	IPA = NSAM/2 + MOD(NSAM,2)
	IPB = NROW/2 + MOD(NROW,2)

        DO J=NROW,1,-1
           DO I=NSAM,1,-1
              A(IPA+I,IPB+J) = A(I,J)
           ENDDO
        ENDDO

	DO J=1,IPB
	  DO I=1,LSD
              A(I,J) = 0.0
	  ENDDO
	ENDDO

	DO J=IPB+NROW+1,LSR
	  DO I=1,LSD
              A(I,J) = 0.0
	  ENDDO
	ENDDO

	DO J=IPB+1,IPB+NROW
	  DO I=1,IPA
              A(I,J) = 0.0
	  ENDDO

	  DO  I=IPA+NSAM+1,LSD
              A(I,J) = 0.0
	  ENDDO
	ENDDO

	END
