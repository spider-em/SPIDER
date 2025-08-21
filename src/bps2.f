
C ++********************************************************************
C                                                                      *
C  UNIX-SPIDER VERSION
C  SIMPLIFIED VERSION                              01/05/94
C                  OPFILEC                         FEB 03 ARDEAN LEITH
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
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
C                                                                      *
C  BPS2(MAXMEM)                                                                    *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C  REPROJECTIONS 3D - SLICES, RICHARDSONS METHOD, 
C  RECONSTRUCTION KEPT IN THE SQUARE TO INTRODUCE OTHER CONSTRAINTS.
C  RECONSTRUCTION FROM NROWL TO NROWH.
C  AVERAGE OUTSIDE THE WINDOW IS SUBTRACTED
C  MIN, MAX RELATE TO THE PROJECTIONS
C  GEOMETRY  CYLINDRICAL
C  
CC23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE BPS2(MAXMEM)

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 
        INCLUDE 'PAR.INC'
C	    PAR includes INTEGER LDPX,LDPY,LDPZ,LDPNMX,LDPNMY,NZ1,LDP,NM,LDPNM

	CHARACTER*1  NULL
 
	CHARACTER(LEN=MAXNAM) :: FINPAT,FINPIC

	PARAMETER  (NILMAX=93,NILMXX=560)
	COMMON     DUMMY(80),BUF(1024),ILIST(NILMXX),
     &    NSAM,NROW,INANG,NN,NMAT,
     &	  LTB,LTBN,K_ANG,K_DM,K_LB,K_MAP,K_IPCUBE,
     &	  K_BCKE,K_PROJ,K_BCKN,K_PRJE,K_SIGMA,
     &    KDM(7),
     &	  IUNIT,Q(1)

	DOUBLE PRECISION  ABA

        DATA              INPIC/100/,IOFF/6/

	NULL = CHAR(0)

C       N - LINEAR DIMENSION OF PROJECTIONS AND RESTORED CUBE
C       NANG - NUMBER OF ANGLES (PROJECTIONS)

        WRITE(NOUT,*)' SINGLE-TILT ITERATIVE 3D RECONSTRUCTION'
	IUNIT = NOUT

 	CALL  FILERD(FINPAT,NLET,NULL,
     &	     'TEMPLATE FOR 2-D PROJECTIONS',IRTFLG)

	CALL  FILERD(FINPIC,NLETI,NULL,'SELECTION DOC',IRTFLG)

	K    = 0
	K2   = 1
	NANG = 0
778	LERR = -1
	IF (NANG .EQ. NILMAX)  THEN
            WRITE(NOUT,*) ' TOO MANY IMAGES, LIST TRUNCATED'
            GOTO  779
	ENDIF

	KP1 = K+1
	CALL  UNSAV(FINPIC,K,INPIC,KP1,Q,1,LERR,K2)
	IF (LERR .EQ. 0)  THEN
	   NANG        = NANG+1
	   ILIST(NANG) = Q(1)
	   K           = K+1
	   GOTO  778
	ENDIF
779	CLOSE(INPIC)

C       NANG - TOTAL NUMBER OF IMAGES

	WRITE(NOUT,2001) NANG
2001	FORMAT('  NUMBER OF IMAGES: ',I0)

C       GET THE ANGLES
 	K_ANG = 1
 	CALL  FILERD(FINPIC,NLETI,NULL,'ANGLES DOC',IRTFLG)

	K2 = 1
        DO K=0,NANG-1
	  LERR=-1
	  CALL  UNSAV(FINPIC,K,INPIC,ILIST(K+1),BUF,2,LERR,K2)
	  IF (LERR .EQ. 0)  THEN
	     Q(K_ANG+K) = BUF(2)
          ELSE
             CALL  ERRT(101,'READING THE ANGLES DOC FILE',NE)
             CLOSE(INPIC)
             RETURN
          ENDIF
	ENDDO
	CLOSE(INPIC)

        CALL RDPRMI(IRI,NSLICE,NOT_USED,
     &	   'RADIUS OF RECONSTRUCTED OBJECT, SLICE HEIGHT')
        RI = IRI

        CALL  RDPRMI(NROWL,NROWH,NOT_USED,'RECONSTRUCTION RANGE IN Y')

C       OPEN ALL THE PROJECTION FILES ....
	DO  K=1,NANG
 	   CALL  FILGET(FINPAT,FINPIC,NLET,ILIST(K),INTFLG)
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,IOFF+K,'O',IFORM,NSAM,NROW,NSL,
     &               MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
	ENDDO

	IF(NROWL.LT.1.OR.NROWL.GT.NROW.OR.NROWH.LT.1.OR.NROWH.GT.NROW
     &		     .OR.NROWL.GT.NROWH)  THEN
	   NROWL = 1
	   NROWH = NROW
	ENDIF
	LCYL = NROWH-NROWL+1

	INANG = NANG

	LDPX   = NSAM/2+1
	LDPY   = NROW/2+1
	LDPZ   = NSLICE/2+1
	LDPNMX = NSAM/2+1
	LDPNMY = NROW/2+1

 	K_DM     = IPALIGN64(K_ANG+NANG)
	K_LB     = IPALIGN64(K_DM+9*NANG)
	K_MAP    = K_LB
	K_IPCUBE = K_MAP

	CALL PREPSL_2(NSAM,NSLICE,NN,NMAT,Q(K_IPCUBE),RI)

	NMAT    = NSAM*NSLICE

	K_PROJ  = IPALIGN64(K_IPCUBE+5*NN)

        K_X     = IPALIGN64(K_PROJ+NSAM*NANG)
        K_IBN   = K_X

	MEMTOT  = K_IBN
C	MEMTOT  = IPALIGN64(K_IBN+NANG*N*N)

	IF (MEMTOT .GT. MAXMEM)  THEN
	   WRITE(NOUT,1001)  MEMTOT
	   WRITE(NOUT,1002)  MAXMEM
	   GOTO 9999
	ENDIF
C       LTB  will be found in READPRO = NANG*Nsam*Nrow

	CALL REDPRO2(NSAM,NROWL,NROWH,NANG,
     &	  Q(K_PROJ),Q(K_ANG),LTB,LTBN,ILIST,Q(K_IPCUBE),NN,Q(K_DM),
     &	  RI,ABA,NOUT)

	K_PRJE = IPALIGN64(K_PROJ+LTB)

C       IN THIS VERSION SIGMA IS ASSUMED TO BE PROPROTIONAL TO PROJ
C       AND THE CORRESPONDING ARRAY IS NOT USED ANYWHERE.

	K_SIGMA = K_PRJE
C	K_SIGMA = K_PRJE+LTB 
C       LTBN    = NSAM*NANG
	K_BCKE  = IPALIGN64(K_PRJE+LTBN)
	K_BCKN  = IPALIGN64(K_BCKE+NMAT*3)
        K_CB    = K_BCKE
	MEMTOT  = IPALIGN64(K_BCKN+NMAT)

	WRITE(NOUT,1001)  MEMTOT
1001	FORMAT(/,'  REPROJECTION MEMORY NEEDED: ',I0,/)
	IF (MEMTOT > MAXMEM)  THEN
	   WRITE(NOUT,1002) MAXMEM
1002	   FORMAT('  YOUR BUFFER LENGTH IS ONLY: ',I0,/,
     &            '  PROGRAM CANNOT RUN')
	   GOTO 9999
	ENDIF

	IFORM = 3
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FINPAT,INPIC,'U',IFORM,NSAM,LCYL,NSLICE,
     &               MAXIM,'OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

	CALL  REPR2_S
     &	  (Q(K_BCKE),Q(K_BCKN),NSAM,LCYL,NSLICE,NROWL,NROWH,NANG,
     &     Q(K_IPCUBE),NN,Q(K_PROJ),Q(K_PRJE),
     &		IRI,LTB,LTBN,ABA,INPIC)
	CLOSE(INPIC)
 
C       CLOSE ALL THE PROJECTION FILES ...
	DO K=1,NANG
           CLOSE(IOFF+K)
	ENDDO
C
9999	CONTINUE
	END
