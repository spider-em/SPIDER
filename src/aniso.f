C ++********************************************************************
C                                                                      *
C ANISO                   CREATED             APRIL 2002 ARDEAN LEITH  * 
C                         STACKS SUPPORT        OCT 2002 ARDEAN LEITH  *
C                                                                      *
C  ANISO(LUN1,LUN2,NSAM,NROW,NSLICE,IRTFLG)
C
C  PARAMETERS: LUN1,LUN2   IO UNITS                             (INPUT)
C              NSAM        X DIMENSIONS                         (INPUT)
C              NROW        Y DIMENSIONS                         (INPUT)
C              NSLICE      Z DIMENSIONS                         (INPUT)
C
C  PURPOSE: ALTER CONTRAST IN AN IMAGE OR VOLUME USING ANISOTROPIC
C           DIFFUSION
C                                                                      *
C **********************************************************************

	SUBROUTINE ANISO(LUN1,LUN2,NSAM,NROW,NSLICE,MAXIM,IRTFLG)

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL, ALLOCATABLE, DIMENSION(:) :: BUF
        CHARACTER (LEN=1) ::               ANS,NULL
        LOGICAL ::                         LOADIT,NORMIT

        CALL RDPRMC(ANS,NC,.TRUE.,'CPF, MCD, OR HEG?',NULL,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        ITER = 10
        IF (ANS .EQ. 'H') ITER = 60
        CALL RDPRI1S(ITER,NOT_USED,'ITERATIONS',IRTFLG)
	IF (IRTFLG .NE. 0) RETURN

        IF (ANS .EQ. 'H') THEN
C          USE HEGERL & FRANGAKIS FORMULATION

C          HT IS TIME STEP
           HT = 0.01
           CALL RDPRM1S(HT,NOT_USED,'TIME STEP(0...0.25)',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          FLAMBDA IS A CONTRAST PARAMETER,  SIGMA IS A NOISE SCALE
           SIGMA   = 1.0
           FLAMBDA = 0.01
           CALL RDPRM2S(SIGMA,FLAMBDA,NOT_USED,'SIGMA & LAMBDA',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
           LOADIT = .FALSE.
        ELSE
C          MEAN CURVATURE OR EED FORMULATION
           IF (ANS .EQ. 'M') THEN
C             MEAN CURVATURE FORMULATION
              W = 1.0
              CALL RDPRM1S(W,NOT_USED,'WEIGHTING FACTOR',IRTFLG)
	      IF (IRTFLG .NE. 0) RETURN
           ENDIF

           ALLOCATE(BUF(NSAM*NROW*NSLICE), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN 
              CALL ERRT(46,'ANISO, BUF...',IER)
              RETURN
           ENDIF
           LOADIT = .TRUE.
        ENDIF
 
        IMGNUM = -3
        DO WHILE (IMGNUM .LT. MAXIM) 
C          INPUT VOLUME/IMAGE

           CALL GETSTACK(LUN1,LUN2,IMGNUM,MAXIM,VERBOSE,LOADIT,BUF,
     &                   LOADIT,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           IF (ANS .EQ. 'H') THEN
C             USE HEGERL & FRANGAKIS FORMULATION

              IF (NSLICE .GT. 3) THEN
C                VOLUME OR VOLUME STACK 
	         CALL ANISOF3(LUN1,LUN2,NSAM,NROW,NSLICE,ITER,HT,
     &                     FLAMBDA,SIGMA,IRTFLG)
              ELSE
C                2D IMAGE OR IMAGE STACK
	         CALL ANISOF(LUN1,LUN2,NSAM,NROW,NSLICE,ITER,HT,
     &                       FLAMBDA,SIGMA,IRTFLG)
               ENDIF
              IF (IRTFLG .NE. 0) GOTO 9999

           ELSE
C             NORMALIZE VOLUME OVER 0....1
              FCON = 1.0 / (FMAX - FMIN) 
              BUF  = (BUF - FMIN) * FCON

              ILOC = 1
  
C             GO THRU VOLUME SLICE-BY-SLICE 
              DO I = 1,NSLICE
                 IF (ANS .EQ. 'M') THEN
C                   MEAN CURVATURE FORMULATION
                    CALL ANISOE_M(W,BUF(ILOC),NSAM,NROW,ITER,IRTFLG)
                 ELSE
C                   CORNER PRESERVING FORMULATION
                    CALL ANISOE(BUF(ILOC),NSAM,NROW,ITER,IRTFLG)
                 ENDIF 
                 IF (IRTFLG .NE. 0) GOTO 9999
                 ILOC = ILOC + NSAM * NROW
              ENDDO

C             UNNORMALIZE VOLUME
              RANGE = (FMAX - FMIN) 
              BUF   = (RANGE * BUF) + FMIN

C             OUTPUT VOLUME
              CALL WRTVOL(LUN2,NSAM,NROW,1,NSLICE,BUF,IRTFLG)
           ENDIF
        ENDDO

9999    IF (ALLOCATED(BUF)) DEALLOCATE(BUF)

        RETURN
        END

