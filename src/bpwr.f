
C++*********************************************************************
C
C BPWR.F 
C              R**2 weighting                              03/04/92
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
C BPWR 
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE BPWR(MAXMEM)

         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC'
         INCLUDE 'F90ALLOC.INC'

         COMMON /IOBUF/ BUF(NBUFSIZ)

         CHARACTER (LEN=MAXNAM) :: FINPAT,FINPIC,FOUT
         CHARACTER (LEN=1)      :: NULL = CHAR(0)

         REAL, ALLOCATABLE      :: Q1(:), Q2(:)
         REAL                   :: PLIST(20)

         INTEGER, PARAMETER     :: NILMAX = 2000 
         INTEGER                :: ILIST(NILMAX)
 
         DATA                    INPIC/99/

         CALL FILERD(FINPAT,NLET,NULL,
     &      'TEMPLATE FOR 2-D IMAGES~',IRTFLG)

         CALL  FILERD(FOUT,NLOT,NULL,
     &      'TEMPLATE FOR 2-D OUTPUT IMAGES~',IRTFLG)

         CALL FILERD(FINPIC,NLETI,NULL,'SELECTION DOC',IRTFLG)

         CALL RDPRM(FM,NOT_USED,
     &      'CUT-OFF FREQUENCY FOR PARZEN FILTER')
         K    = 0
         K2   = 1
         NANG = 0
778      LERR = -1
         IF (NANG == NILMAX)  THEN
            WRITE(NOUT,*) '*** TOO MANY IMAGES, LIST TRUNCATED'
            GOTO  779
         ENDIF

         KP1 = K+1
         CALL  UNSAV(FINPIC,K,INPIC,KP1,PLIST,1,LERR,K2)
         IF (LERR == 0)  THEN
            NANG        = NANG+1
            ILIST(NANG) = PLIST(1)
            K           = K+1
            GOTO 778
         ENDIF

779      CLOSE(INPIC)
         WRITE(NOUT,2001) NANG
2001     FORMAT('  NUMBER OF IMAGES:',I0)

         CALL  FILGET(FINPAT,FINPIC,NLET,ILIST(1),INTFLG)
         MAXIM = 0
         CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,NSAM,NROW,NSL,
     &               MAXIM,' ',.FALSE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
         CLOSE(INPIC)

         N2S = 2*NSAM
         N2R = 2*NROW

         ALLOCATE(Q1(N2R*(N2S+2)),Q2(N2R*(N2S+2)), STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            MWANT = 2*N2R*(N2S+2)  
            CALL ERRT(46,'BPRW; Q1...',MWANT)
         ENDIF

	 IF (FM .GE. 0.0)  THEN
            WRITE(NOUT,1001)  N2S,N2R
1001        FORMAT(//'  R**2 Weighting of 2D images',/,
     &               '  Dimensions used:',2I8,/)
	 ELSE
            WRITE(NOUT,1003)  N2S,N2R
1003        FORMAT(//'  R*   Weighting of 2D images',/,
     &               '  Dimensions used:',2I8,/)
	 ENDIF

         IF (MEMTOT .GT. MAXMEM)  THEN
            CALL ERRT(102,'BUFFER LENGTH IS ONLY',MAXMEM)
            RETURN
         ENDIF

         CALL  BPWR_Q(BUF,Q1,Q2,FM,ILIST,NANG,N2S,N2R,
     &                NSAM,NROW,FINPAT,NLET,FOUT,NLOT)

         IF (ALLOCATED(Q1))  DEALLOCATE(Q1)
         IF (ALLOCATED(Q2))  DEALLOCATE(Q2)

         END
