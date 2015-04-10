C++*********************************************************************
C
C  BCQ.F                                                   02/06/97
C                             USED ALLOCATE NOT CHKMEM  DEC 2000 al
C                             USED OPFILEC              FEB 2003 al
C                             INPUT X,Y,Z TOGETHER      MAY 2003 al
C                             ALLOCATE & PARTITION      MAY 2003 al
C                             REANG --> BUILDM          JUL 2003 al
C                             BUILDM BUG                SEP 2003 al
C                             NSLICE2 BUG               APR 2004 al
C                             CW ALLOCATION BUG         JAN 2005 al
C                             BETTER ERROR MSG          AUG 2006 al
C
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
C
C   BCQ(UNUSED)
C
C   PURPOSE: CALCULATES BACK-PROJECTION STEP OF 3D RECONSTRUCTION 
C            USING THREE EULERIAN ANGLES.  OPTIONALLY ONE OF TWO 
C            POSSIBLE WEIGHTING FUNCTIONS IS APPLIED. 
C
C   PARAMETERS:    UNUSED                                    (UNUSED)
C
C   CALL TREE:      BCQ --------> BUILDM -----> CANG
C                    |            BPCQP -----> WTF --> FMRS_2
C                    |                         WTM
C                    |                         BPCQ
C                    |
C                    -----------> BUILDM -----> CANG
C                                 BPCMP -----> WTF --> FMRS_2
C                                              WTM
C                                              BPCM
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

	SUBROUTINE BCQ(UNUSED)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

	COMMON /PAR/    LDPX,LDPY,LDPZ,LDPNMX,LDPNMY,NZ1


        INCLUDE 'F90ALLOC.INC'
        REAL, DIMENSION(:,:), POINTER       :: ANGBUF
        REAL, ALLOCATABLE, DIMENSION(:,:)   :: DM,SS,CW
        REAL, ALLOCATABLE, DIMENSION(: )    :: CB,ILISTP
        REAL, ALLOCATABLE, DIMENSION(:,:)   :: PROJ
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: PROJS
        LOGICAL                             :: PARTITION

        CHARACTER(LEN=MAXNAM)               :: ANGDOC,FINPIC,FINPAT

	DATA  LUNDOC/97/,IOPIC/98/,INPIC/99/
	
C       READ INPUT TEMPLATE AND SELECTION DOC FILE CONTAINING IMAGE NO. 
        NILMAX = NIMAX
        CALL FILELIST(.TRUE.,LUNDOC,FINPAT,NLET,INUMBR,NILMAX,NANG,
     &                'TEMPLATE FOR INPUT IMAGES~',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       NANG - NUMBER OF ANGLES (PROJECTIONS)
	WRITE(NOUT,2001) NANG
2001	FORMAT(' TOTAL NUMBER OF IMAGES: ',I6)

        MAXXT = 4
        MAXYT = 0
        DO I = 1,NANG
           MAXYT = MAX(INUMBR(I),MAXYT)
        ENDDO

        CALL GETDOCDAT('ANGLES DOC',.TRUE.,ANGDOC,LUNDOC,.FALSE.,MAXXT,
     &                 MAXYT,ANGBUF,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        NZ3D = -1
        CALL RDPRI3S(NX3D,NY3D,NZ3D,NOT_USED,
     &                'OUTPUT VOLUME: X, Y & Z  DIMENSIONS',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (NZ3D .LE. -1) THEN
           CALL RDPRI1S(NZ3D,NOT_USED,
     &                  'OUTPUT VOLUME:  Z  DIMENSION',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ENDIF

	NZ1 = 1
	NZ2 = NZ3D
	CALL  RDPRIS(NZ1,NZ2,NOT_USED,
     &		'FIRST, LAST SLICE TO BE RECONSTRUCTED',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (NZ1 .LT. 1   .OR. NZ1 .GT.(NZ3D-1) .OR. 
     &      NZ2 .LE. NZ1 .OR. NZ2 .GT. NZ3D) THEN
            CALL ERRT(14,'A SLICE IS OUTSIDE VOLUME',NE)
            RETURN
        ENDIF 
	NZC = NZ2 - NZ1 + 1

	CALL  RDPRM(SNR,NOT_USED,'SNR/DIAMETER')
	IF (SNR .GT. 0.0)  SNR = 1.0 / SNR

        IFORM = 3
        CALL OPFILEC(0,.TRUE.,FINPIC,IOPIC,'U',IFORM,NX3D,NY3D,NZC,
     &               MAXIM,'RECONSTRUCTED 3-D OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       OPEN FIRST IMAGE FILE TO DETERMINE NSAM, NROW, & NSL
 	CALL FILGET(FINPAT,FINPIC,NLET,INUMBR(1),INTFLG)
        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,NSAM,NROW,NSL,
     &             MAXIM,'DUMMY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
	CLOSE(INPIC)

	LDPX   = NX3D/2+1
	LDPY   = NY3D/2+1
	LDPZ   = NZ3D/2+1
	LDPNMX = NSAM/2+1
	LDPNMY = NROW/2+1
	NMAT   = NX3D*NY3D*NZC
	NNNN   = NSAM+2-MOD(NSAM,2)

        ALLOCATE(DM(9,MAXYT),SS(6,MAXYT),CW(NNNN/2,NROW),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'DM,SS,&CW',15*MAXYT+NNNN/2*NROW) 
           GOTO 9999
        ENDIF

        PARTITION = (FCHAR(6:6) .EQ. 'P')
C       IF PARTITION BUILDM RETURNS DM & SS FOR ANGLES KEYED BY INUMBR
        CALL BUILDM(INUMBR,DM,NANG,ANGBUF,.TRUE.,SS,PARTITION,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        WRITE(NOUT,91) NANG
91      FORMAT(' PROJECTION ANGLES CREATED:',I8)

        NANGP = NANG
        IF (PARTITION) THEN
           ALLOCATE(ILISTP(NANG),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'ILISTP',NANG) 
              GOTO 9999
           ENDIF

C          READ INPUT TEMPLATE AND SELECTION DOC FILE CONTAINING IMAGE NO. 
           CALL FILELIST(.FALSE.,LUNDOC,FINPAT,NLET,ILISTP,NANG,NANGP,
     &       'FILE NUMBERS OR SELECTION DOC. FILE FOR THIS PARTITION~',
     &       IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          NANGP - NUMBER OF ANGLES (PROJECTIONS) IN THIS PARTITION
	   WRITE(NOUT,92) NANGP
92	   FORMAT(' NUMBER OF IMAGES IN THIS PARTITION: ',I6)
        ENDIF

C       VOLUME CB MAY BE TOO LARGE TO ALLOCATE??
        ALLOCATE(PROJ(NNNN,NROW),CB(NMAT),STAT=IRTFLG)

	IF (IRTFLG .EQ. 0)  THEN
C          VOLUME ALLOCATION SUCCESSFUL
C          3-D BACK-PROJECTION WITH VOLUME & ALL PROJECTIONS IN MEMORY

           IF (PARTITION) THEN
   	      WRITE(NOUT,93) NANGP
93            FORMAT(/,' 3-D BACK-PROJECTION WITH VOLUME AND:, ',I5,
     &                 ' PROJECTIONS IN MEMORY',/)

              CALL BPCQP(PROJ,CW,NNNN,NSAM,NROW,CB,NX3D,NY3D,NZC,
     &	         INUMBR,ILISTP,DM,SS,NANG,NANGP,SNR,FINPAT(1:NLET),
     &           FINPIC,INPIC)
           ELSE
   	      WRITE(NOUT,94) NANG
94            FORMAT(/,' 3-D BACK-PROJECTION WITH VOLUME & ALL: ',I5,
     &                 ' PROJECTIONS IN MEMORY',/)

              CALL BPCQP(PROJ,CW,NNNN,NSAM,NROW,CB,NX3D,NY3D,NZC,
     &	         INUMBR,INUMBR,DM,SS,NANG,NANG,SNR,FINPAT(1:NLET),
     &           FINPIC,INPIC)
           ENDIF
           CALL WRTVOL(IOPIC,NX3D,NY3D,1,NZC,CB,IRTFLG)

	ELSE
C          WHOLE VOLUME ALLOCATION NOT SUCCESSFUL
C          3-D BACK-PROJECTION WITH VOLUME & SOME PROJECTIONS IN MEMORY

           IF (PARTITION) THEN
              MWANT = NNNN*NROW + NMAT 
              CALL ERRT(46,'PROJ & CB (PROJECTIONS & OUTPUT VOLUME)',
     &                  MWANT) 
              GOTO 9999
           ENDIF

           ALLOCATE(CB(NMAT),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
   	     WRITE(NOUT,*) ' *** TRY PARTITIONING YOUR VOLUME'

             CALL ERRT(46,'OUTPUT VOLUME TOO LARGE',NMAT) 
             GOTO 9999
           ENDIF

C          FIND HOW MANY PROJECTIONS CAN FIT IN MEMORY
C          al 2006 THIS IS LIKELY TO CRASH LATER DUE TO NEED FOR
C          ADDITIONAL STACK LOCATED MEMORY???

	   DO LPRJ=NANG,1,-1
              ALLOCATE(PROJS(LPRJ,NNNN,NROW),STAT=IRTFLG)
              IF (IRTFLG .EQ. 0) EXIT
	   ENDDO

           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'PROJS (OUTPUT VOLUME & PROJECTIONS)',
     &                  NNNN*NROW) 
              GOTO 9999
           ENDIF

   	   WRITE(NOUT,93) LPRJ

           CALL BPCMP(PROJS,CW,NNNN,NSAM,NROW,LPRJ,CB,
     &	         NX3D,NY3D,NZC,INUMBR,DM,SS,NANG,SNR,
     &           IOPIC,FINPAT(1:NLET),FINPIC,INPIC)
        ENDIF

C       DEALLOCATE  ARRAYS
9999    IF (ALLOCATED(PROJ))    DEALLOCATE(PROJ)
        IF (ALLOCATED(PROJS))   DEALLOCATE(PROJS)
        IF (ALLOCATED(DM))      DEALLOCATE(DM)
        IF (ALLOCATED(SS))      DEALLOCATE(SS)
        IF (ALLOCATED(CW))      DEALLOCATE(CW)
        IF (ALLOCATED(CB))      DEALLOCATE(CB)
        IF (ALLOCATED(ILISTP))  DEALLOCATE(ILISTP)
        IF (ASSOCIATED(ANGBUF)) DEALLOCATE(ANGBUF)

	CLOSE(IOPIC)

        RETURN
	END
