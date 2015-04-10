
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

        SUBROUTINE ATMC
	
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
        INCLUDE 'F90ALLOC.INC'
	
	INTEGER  X42

        CHARACTER (LEN=MAXNAM)              :: FINPIC,FINPAT
        REAL, ALLOCATABLE, DIMENSION(:,:)   :: XYZ
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: CCF
 
        DIMENSION     DLIST(5)

        DATA NDOC/55/,INPIC/56/


        CALL RDPRMI(NGH,NPEAK, NOT_USED,
     &            'NEIGHBOURHOOD PIXELS FOR SEARCH')

C       NUMBER OF PEAKS')
        IF (NGH .LT. 3) THEN
           CALL ERRT(31,'ATMC',NE)
           RETURN
        ENDIF
        NQ = NGH/2
        DD=NGH*NGH

        NILMAX = NIMAX

        CALL FILELIST(.TRUE.,INPIC,FINPAT,NLET,INUMBR,NILMAX,NCCF,
     &                 'TEMPLATE FOR 2-D CCFs',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        MAXNUM = MAXVAL(INUMBR(1:NCCF))
        CLOSE(INPIC)
C       OPEN FIRST IMAGE FILE TO DETERMINE NSAM, NROW, NSL
        CALL FILGET(FINPAT,FINPIC,NLET,INUMBR(1),INTFLG)

        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,NSAM,NROW,NSL,
     &               MAXIM,'DUMMY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CLOSE(INPIC)

        ALLOCATE (CCF(NSAM,NROW,NCCF), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'AT MC, XYZ',IER)
           RETURN
        ENDIF

C       MAXIMUM NUMBER OF POSSIBLE PEAKS 
        ITMP = (NSAM/NGH +1)*(NROW/NGH +1)

        ALLOCATE (XYZ(4,ITMP), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'AT MC, XYZ',IER)
           RETURN
        ENDIF


	DO  I=1,NCCF
        CALL FILGET(FINPAT,FINPIC,NLET,INUMBR(I),INTFLG)

        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,NSAM,NROW,NSL,
     &               MAXIM,'DUMMY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
	CALL  READV(INPIC,CCF(1,1,I),NSAM,NROW,NSAM,NROW,1)
        CLOSE(INPIC)
	ENDDO

        CALL PKDMC(NSAM,NROW,NCCF,CCF,NQ,DD,XYZ,ITMP,L)

        NLIST    = 2
        DLIST(1) = -1
        DLIST(2) = L
        CALL SAVD(NDOC,DLIST,NLIST,IRTFLG)

        NLIST = 5
	
        DO I1=1,L
           DLIST(1) = I1
           DLIST(2) = XYZ(1,I1)
           DLIST(3) = XYZ(2,I1)
           DLIST(4) = XYZ(3,I1)
           DLIST(5) = XYZ(4,I1)
           CALL SAVD(NDOC,DLIST,NLIST,IRTFLG)
        ENDDO

        CALL SAVDC

        CLOSE(NDOC)

	ZM=MAXVAL(XYZ(3,1:L))

        NMAX = 0
        CALL  FILSEQP(FINPAT,NLET,INUMBR,NMAX,NIMA,
     &                'TEMPLATE FILENAME (ex: PROFILE***)',IRTFLG)
C
        NLS=2
        DO X42 = 1, L

C	  CREATE A FILENAME FOR CLASS X42
          CALL FILGET(FINPAT,FINPIC,NLET,X42,IRTFLG)
	  I=XYZ(1,X42)
	  J=XYZ(2,X42)

	  DO IC = 1, NCCF
              DLIST(1) = IC
              DLIST(2) = CCF(I,J,IC)/ZM
              IAP = 0
              CALL SAVDN1(NDOC,FINPIC,DLIST,NLS,IC-1,IAP)
	  END DO

          CLOSE(NDOC)
	END DO

#ifdef NEVER
        DO I1=1,L
	WRITE(88,100)  L
100	FORMAT(I6)
	WRITE(88,101)  (CCF(I,J,IC)/ZM,IC=1,NCCF)
101	FORMAT(5X,F12.5)
	ENDDO
	CLOSE(88)
#endif
        DEALLOCATE(XYZ,CCF)

        END




        SUBROUTINE PKDMC(NSAM,NROW,NCCF,CCF,NQ,DD,XYZ,ITMP,L)

        DIMENSION CCF(NSAM,NROW,NCCF)
        DIMENSION XYZ(4,ITMP) 
        LOGICAL   T
C
        L=0
        DO    J=NQ+1,NROW-NQ
           DO    I=NQ+1,NSAM-NQ
	      Z=-HUGE(Z)
	      DO  IC=1,NCCF
	        IF(CCF(I,J,IC).GT.Z)  THEN
		 NIC=IC
                 Z=CCF(I,J,IC)
	        ENDIF
	      ENDDO
	      TMP=CCF(I,J,NIC)
	      IF (.NOT. (
     &		TMP>CCF(I-1,J-1,NIC) .AND.
     &		TMP>CCF(I-1,J,NIC) .AND.
     &		TMP>CCF(I-1,J+1,NIC) .AND.
     &		TMP>CCF(I,J-1,NIC) .AND.
     &		TMP>CCF(I,J+1,NIC) .AND.
     &		TMP>CCF(I+1,J-1,NIC) .AND.
     &		TMP>CCF(I+1,J,NIC) .AND.
     &		TMP>CCF(I+1,J+1,NIC) ) ) GOTO  7

                 XX=FLOAT(I)
                 YY=FLOAT(J)
 
		 LB=1
	         NEWL=L
 6               L=NEWL
                 IF(L>0)  THEN
C                   Check whether there were any other peaks in the vicinity
                    DO  LT=1,L
                       IF( (XYZ(1,LT)-XX)**2+
     &                    (XYZ(2,LT)-YY)**2<=DD) THEN
                          IF(XYZ(3,LT)<Z)  THEN
C                               Found a smaller peak, remove the smaller one
				IF(L>1)  THEN
				XYZ(:,LT)=XYZ(:,L)
				NEWL=L-1
				LB=LT
				ELSE
				NEWL=0
				ENDIF
				GOTO  6
                          ENDIF
C                         Found a larger peak, skip the current one
                          GOTO 7
                       ENDIF
                    ENDDO
                 ENDIF
                 L=L+1
                 XYZ(1,L)=XX
                 XYZ(2,L)=YY
                 XYZ(3,L)=Z
                 XYZ(4,L)=NIC

7             CONTINUE
	
           ENDDO

        ENDDO

        END
