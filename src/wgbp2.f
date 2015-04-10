
C ++********************************************************************
C                                                                      *
C  WGBP2        SIMPLIFIED VERSION               01/04/94 PAWEL PENCZEK
c               OPFILEC                          FEB  03 ARDEAN LEITH
C               SPEEDED UP (36%) & USED ALLOC    APR  03 ARDEAN LEITH
C                                                                      *
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
C  WGBP2(MAXMEM) 
C                                                               
C  PURPOSE:   WEIGHTED BACK-PROJECTION WITH PARZEN FILTER.
C             RECONSTRUCTION KEPT IN THE SQUARE
C             RECONSTRUCTION FROM NROWL TO NROWH.
C             AVERAGE OUTSIDE THE WINDOW IS SUBTRACTED
C             GEOMETRY  CYLINDRICAL
C
C  CALL TREE:  WGBP2  ---->   PREPSL_2L
C                             RDPA
C                             FFTR_Q
C                             BCKC0_L

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE WGBP2(MAXMEM)

C       THIS ROUTINE OPENS FILES ON LUNS: IOFF....IOFF+NILMAX-1
#ifdef SP_MP
	PARAMETER  (NILMAX=25 )
#else
	PARAMETER  (NILMAX=91 )
#endif

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 
        INCLUDE 'F90ALLOC.INC'

	CHARACTER(LEN=MAXNAM) :: FINPAT,FINPIC,FINFO,ANGDOC
	COMMON  /COMMUN/         FINPAT,FINPIC,FINFO,ANGDOC

        REAL, DIMENSION(:,:), POINTER     :: PANG
        REAL, ALLOCATABLE, DIMENSION(:)   :: Q_ANG,Q_PROJ,Q_WGH,Q_BCK
        REAL, ALLOCATABLE, DIMENSION(:,:) :: DM,Q_PRJE
	DOUBLE PRECISION                  :: ABA
	LOGICAL                           :: OPENED

	COMMON        IPCUBE(5,1)

        DATA  LUNDOC/70/,LUNOUT/99/,IOFF/7/

        NUMMAX = NIMAX

        WRITE(NOUT,*)  ' SINGLE-TILT 3D-WEIGHTED BACKPROJECTION'

C       GET THE INPUT IMAGE FILE LISTING
        CALL FILELIST(.TRUE.,LUNOUT,FINPAT,NLET,
     &     INUMBR,NUMMAX,NANG,'TEMPLATE FOR 2-D IMAGES',IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(101,'WITH FILE TEMPLATE',NE)
           GOTO 9999
        ENDIF  

C       NANG - NUMBER OF ANGLES (PROJECTIONS)
	WRITE(NOUT,2001) NANG
2001	FORMAT('  NUMBER OF PROJECTION IMAGES: ',I0)

C       RETRIEVE ARRAY WITH ANGLES DATA IN IT
        MAXNUM = MAXVAL(INUMBR(1:NANG))
        MAXXT  = 4
        MAXYT  = MAXNUM
        CALL GETDOCDAT('ANGLES DOCUMENT',.TRUE.,ANGDOC,LUNDOC,.FALSE.,
     &                 MAXXT,MAXYT,PANG,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(4,'WGBP2 ',NE)
           GOTO 9999
        ENDIF       

        ALLOCATE (Q_ANG(NANG),DM(9,NANG), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'Q_ANG...',IER)
           GOTO 9999
        ENDIF

        DO I = 1, NANG
            Q_ANG(I) = PANG(3,INUMBR(I))
        ENDDO
        IF (ASSOCIATED(PANG))  DEALLOCATE(PANG)

C       OPEN ALL THE PROJECTION FILES,IF THE NUMBER IS LESS THAN NILMAX
	IF (NANG .LE. NILMAX)  THEN
	   OPENED = .TRUE.
           KE     = NANG
        ELSE
           OPENED = .FALSE.
           KE     = 1
        ENDIF

        MAXIM = 0
	DO K=1,KE
 	   CALL  FILGET(FINPAT,FINPIC,NLET,INUMBR(K),INTFLG)
           LUN = IOFF + K
           CALL OPFILEC(0,.FALSE.,FINPIC,LUN,'O',ITYPE,NSAM,NROW,NSL,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
 	   IF (IRTFLG .NE. 0) GOTO 9999
	ENDDO
        IF (.NOT. OPENED) CLOSE(IOFF+1)

        IRI    = NSAM / 2
        NSLICE = NSAM
        CALL RDPRIS(IRI,NSLICE,NOT_USED,
     &	   'RADIUS & HEIGHT OF RECONSTRUCTED OBJECT',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        RI = IRI

        NROWL = 1
        NROWH = NROW
        CALL  RDPRIS(NROWL,NROWH,NOT_USED,
     &	            'RECONSTRUCTION RANGE IN Y',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

	IF (NROWL.LT.1 .OR. NROWL.GT.NROW .OR.
     &      NROWH.LT.1 .OR. NROWH.GT.NROW .OR.
     &	    NROWL.GT.NROWH)  THEN
	   NROWL = 1
	   NROWH = NROW
        ENDIF
	LCYL = NROWH - NROWL + 1

        FM = 0.3
        CALL RDPRM1S(FM,NOT_USED,
     &		   'FREQUENCY CUT-OFF FOR PARZEN FILTER',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

	MK = 1
	M2 = 0
5	MK = MK*2
	M2 = M2+1
	IF (MK .LT. NSAM-1)  GOTO  5

 	MK     = MK*2
	M2     = M2+1

	LDPX   = NSAM/2+1
	LDPY   = NROW/2+1
	LDPZ   = NSLICE/2+1
	LDPNMX = NSAM/2+1
	LDPNMY = NROW/2+1

C       CREATES DM MATRIX HERE
	CALL RDPA(NANG,Q_ANG,DM)

C       RETURNS VALUE OF NN AND FILLS THE IPCUBE ARRAY
	CALL PREPSL_2L(NSAM,NSLICE,NN,NMAT,IPCUBE,RI,LDPX,LDPZ)
	NMAT = NSAM * NSLICE

	MEMTOT  = NANG + 9*NANG + 2*MK + NSAM*NANG + NSAM*NSLICE
	WRITE(NOUT,90)  MEMTOT *4 
90	FORMAT(/,'  MEMORY ALLOCATION (BYTES): ',I8,/)

        ALLOCATE(Q_PROJ(MK),Q_PRJE(NSAM,NANG),Q_WGH(MK),
     &           Q_BCK(NSAM*NSLICE),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'Q_PROJ...',IER)
           GOTO 9999
        ENDIF

C       OPEN OUTPUT VOLUME
	IFORM = 3
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FINFO,LUNOUT,'U',IFORM,NSAM,LCYL,NSLICE,
     &                   MAXIM,'3-D OUTPUT',.FALSE.,IRTFLG)
	IF (IRTFLG .NE. 0)  GOTO 9999

C   On Linux when the file is created and written in a nmon-consecutive order,
C    the resulting file is fragmented and access is very slow.  Thus, first
C    we have to write the whole file in a consecutive order, close it,
C    and open again.   Hopefuly, this will fix the problem.  PAP  08/03/05.

	Q_BCK(1:NSAM)=0.0
	DO  K=1,NSLICE
	   DO  J=1,LCYL
	      CALL  WRTLIN(LUNOUT,Q_BCK,NSAM,(K-1)*LCYL+J)
	   ENDDO
	ENDDO
	CLOSE(LUNOUT)

C       Open again!
        CALL OPFILEC(0,.FALSE.,FINFO,LUNOUT,'O',IFORM,NSAM,LCYL,NSLICE,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
	IF (IRTFLG .NE. 0)  GOTO 9999
C	
	DO I=3,MK,2
	   Q_WGH(I)   = I/2
	   Q_WGH(I+1) = Q_WGH(I)
	ENDDO

	Q_WGH(1) = 1
	Q_WGH(2) = MK/2

	IF (FM .NE. 0.0)  THEN
	   DO I=3,MK,2
	      J = I/2
	      FQ   = FLOAT(J)/FLOAT(MK)
	      PARZ = 0.0
	      IF (FQ .LE. FM) THEN
                 IF (FQ .LE. (FM/2.0)) THEN
  	            PARZ = 1.0-6.0*(FQ/FM)**2*(1.0-FQ/FM)
	         ELSE
  	            PARZ = 2.0*(1.0-FQ/FM)**3
                 ENDIF
              ENDIF
  	      Q_WGH(I)   = Q_WGH(I)   * PARZ
	      Q_WGH(I+1) = Q_WGH(I+1) * PARZ
	   ENDDO

	   FQ   = 0.5
	   PARZ = 0.0
	   IF (FQ .LE. FM) THEN
              IF (FQ .LE. (FM/2.0)) THEN
                 PARZ = 1.0-6.0*(FQ/FM)**2*(1.0-FQ/FM)
              ELSE
                 PARZ = 2.0*(1.0-FQ/FM)**3
              ENDIF
           ENDIF
 	   Q_WGH(2) = Q_WGH(2)*PARZ
       	ENDIF

	NMAT = NSAM * NSLICE

C       LOOP OVER ALL ROWS IN THE OUTPUT VOLUME
	DO LS=NROWL,NROWH

C          GET PROJECTION DATA FROM INPUT FILE
	   IF (OPENED) THEN 
	      DO J=1,NANG
	         CALL REDLIN(IOFF+J,Q_PRJE(1,J),NSAM,LS)
	      ENDDO
	   ELSE
	      DO J=1,NANG
                 MAXIM = 0
	         CALL FILGET(FINPAT,FINPIC,NLET,INUMBR(J),INTFLG)
                 CALL OPFILEC(0,.FALSE.,FINPIC,IOFF+1,'O',ITYPE,
     &                   NSAM,NROW,NSL, MAXIM,' ',.FALSE.,IRTFLG)
 	         IF (IRTFLG .NE. 0)  RETURN

	         CALL REDLIN(IOFF+1,Q_PRJE(1,J),NSAM,LS)
	         CLOSE(IOFF+1)
	      ENDDO
	   ENDIF

C          CALCULATE AVERAGE
	   IF (LS .EQ. NROWL) THEN
	      ABA = 0.0D0
c$omp      parallel do private(i,j),reduction(+:aba)
	      DO J=1,NANG
	        DO I=1,NSAM
	           ABA = Q_PRJE(I,J) + ABA
	        ENDDO
	      ENDDO
	      ABA = ABA / NANG / NSAM
	   ENDIF

C          REMOVE AVERAGE
	   DO J=1,NANG
c$omp parallel sections
c$omp section
	      DO I=1,NSAM
	        Q_PROJ(I) = Q_PRJE(I,J) - ABA
	      ENDDO
c$omp section
	      DO I=NSAM+1,MK
	         Q_PROJ(I) = 0.0
	      ENDDO
c$omp end parallel sections
	      CALL  FFTR_Q(Q_PROJ,M2)
c$omp      parallel do private(i)
	      DO I=1,MK
	         Q_PROJ(I) = Q_PROJ(I)*Q_WGH(I)
	      ENDDO
	      CALL FFTR_Q(Q_PROJ,-M2)
c$omp      parallel do private(i)
	      DO I=1,NSAM
	         Q_PRJE(I,J) = Q_PROJ(I)
	      ENDDO
	   ENDDO

C          ZERO THE OUTPUT SLAB
c$omp      parallel do private(i)
           DO I=1,NMAT
              Q_BCK(I) = 0.0
           ENDDO

C          BACK-PROJECT OVER ALL THE ANGLES
           DO I=1,NANG
              CALL BCKC0_L(Q_BCK,NMAT,DM(1,I),Q_PRJE(1,I),NSAM,
     &                    IPCUBE,NN,LDPNMX)
           ENDDO

C          DIVIDE BY NUMBER OF PIXELS PROJECTED
           FCON = 1.0  / NANG / NSAM
c$omp      parallel do private(i)
           DO I=1,NMAT
              Q_BCK(I) = Q_BCK(I) * FCON
	   ENDDO

C          SAVE THIS SLAB IN OUTPUT
	   DO IRC=1,NSLICE
	     IRB = 1+(IRC-1 )* NSAM
	     CALL WRTLIN(LUNOUT,Q_BCK(IRB),NSAM,(IRC-1)*LCYL+LS-NROWL+1)
	   ENDDO

	ENDDO

C       CLOSE OUTPUT FILE
	CLOSE(LUNOUT)
 
C       CLOSE ALL PROJECTION FILES ,IF THE NUMBER IS LESS THAN NILMAX
	IF (NANG .LE. NILMAX)  THEN
	   DO K=1,NANG
 	     CLOSE(IOFF+K)
	   ENDDO
	ENDIF

9999	IF (ALLOCATED(Q_PROJ))  DEALLOCATE (Q_PROJ)
   	IF (ALLOCATED(Q_PRJE))  DEALLOCATE (Q_PRJE)
   	IF (ALLOCATED(Q_WGH))   DEALLOCATE (Q_WGH)
   	IF (ALLOCATED(Q_BCK))   DEALLOCATE (Q_BCK)
   	IF (ALLOCATED(Q_ANG))   DEALLOCATE (Q_ANG)
   	IF (ALLOCATED(DM))      DEALLOCATE (DM)

        RETURN
	END



C ++********************************************************************
C                                                                      *
C                                                                      *
C **********************************************************************
C                                                                      *
C  PREPSL_2(NSAM,NSLICE,NN,NMAT,IPCUBE,RI)                                                                    *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C  IPCUBE: 1 - BEGINNING
C          2 - END
C          3 - IX
C          4 - IY
C          5 - IZ
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE PREPSL_2L(NSAM,NSLICE,NN,NMAT,IPCUBE,RI,LDPX,LDPZ)

	INTEGER       IPCUBE(5,*)
	LOGICAL       FIRST

	R    = RI * RI
	NN   = 0
	NMAT = 0

	DO I1=1,NSLICE
	   T     = I1 - LDPZ
	   XX    = T * T
	   FIRST = .TRUE.
	   DO I3=1,NSAM
	      NMAT = NMAT + 1
	      T    = I3 - LDPX
	      RC   = T * T + XX
	      IF (FIRST) THEN
	         IF (RC .LE. R) THEN
                    FIRST        = .FALSE.
	            NN           = NN + 1
	            IPCUBE(1,NN) = NMAT
	            IPCUBE(2,NN) = NMAT
	            IPCUBE(3,NN) = I3 - LDPX
	            IPCUBE(4,NN) = 1
	            IPCUBE(5,NN) = I1 - LDPZ
                 ENDIF
	      ELSE
	         IF (RC .LE. R) IPCUBE(2,NN) = NMAT
	      ENDIF
	   ENDDO
	ENDDO

	END


C ++********************************************************************
C                                                                      *
C                                                                      *
C **********************************************************************
C                                                                      *
C  BCKC0_L(CUBE,LTC,DM,B,NSAM,IPCUBE,NN,LDPX,LDPZ,LDPNMX)                                                                    *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

       SUBROUTINE BCKC0_L(CUBE,LTC,DM,B,NSAM,IPCUBE,NN,LDPNMX)

       DIMENSION  DM(9),CUBE(LTC),B(NSAM)
       INTEGER    IPCUBE(5,NN)

       DM1 = DM(1)
       DM3 = DM(3)

c$omp  parallel do  private(i,xb,xbb,j,iqx,dipx)
       DO I=1,NN
C          XB = (IPCUBE(3,I)-LDPX)*DM(1)+(IPCUBE(5,I)-LDPZ)*DM(3)
C	   XB = IPCUBE(3,I) *DM1 + IPCUBE(5,I) * DM3
	   XB = IPCUBE(3,I) *DM1 + IPCUBE(5,I) * DM3 - 
     &          IPCUBE(1,I) *DM1 + LDPNMX

	   DO J=IPCUBE(1,I),IPCUBE(2,I)
C	      XBB     = (J - IPCUBE(1,I)) * DM1 + XB
C	      XBB     =  J * DM1 - IPCUBE(1,I) * DM1 + XB
	      XBB     =  J * DM1 + XB

C	      IQX     = IFIX(XBB + FLOAT(LDPNMX))
	      IQX     = IFIX(XBB)

C	      DIPX    = XBB + LDPNMX - IQX
C	      DIPX    = XBB - IQX

C	      CUBE(J) = CUBE(J) + B(IQX) + DIPX * (B(IQX+1) - B(IQX))
	      CUBE(J) = CUBE(J) + B(IQX) +(XBB - IQX)*(B(IQX+1)-B(IQX))
	   ENDDO
       ENDDO

       END

C     1                 +(1.0-DIPX)*B(IQX)
C     2                 +     DIPX *B(IQX+1)

