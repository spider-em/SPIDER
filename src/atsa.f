C ++********************************************************************
C                                                                      *
C ATSA                                                                 *
C                  OPFILEC                         FEB 03 ARDEAN LEITH
C                                                                      *
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
C***********************************************************************

	SUBROUTINE ATSA(MAXMEM)

	PARAMETER  (NILMAX=4000)
	PARAMETER (MAXREG=3)
	PARAMETER (MAXKEY=2000)
	PARAMETER (MVAR=8)

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 

	COMMON  BUF(1024),Q(1)
	COMMON /DIMSS1/ K_Q,K_2,K_3,K_4,K_5,K_6,K_7,K_8
	DIMENSION IR(MAXREG,MAXKEY),VART(MVAR)
        CHARACTER(LEN=MAXNAM)   ::   FINPAT,FINPIC,FIL6,FIL5,DOCM
	CHARACTER*1 NULL

	DATA  INPIC/77/,LUN50/51/,LUN51/28/

	NULL=CHAR(0)
	NMAX=NILMAX

        CALL FILERD(FINPAT,NLET,NULL,
     &     'TEMPLATE FOR WINDOW SERIES',IRTFLG)
                                                       
        CALL RDPRMI(NNSAM,NNROW,NOT_USED,'SIZE OF MINI WINDOW')
        OPEN(UNIT=4,FILE='SCR',STATUS='UNKNOWN')

        CALL FILERD(DOCM,NLETDM,DATEXC,
     &             'DOCUMENT WITH CATEGORIES',IRTFLG)
	IKEY=1
        ISW=0
        CALL UNSDAL(DOCM,ISW,44,IKEY,RLIST,1,Q,MAXKEY,
     &              MAXREG,NMG,LERR)

        DO  II=1,NMG
           JI=(II-1)*MAXREG+2
           IR(1,II)=Q(JI)
           IR(2,II)=Q(JI+1)
	ENDDO
 
	DO  I=1,NMG
 	   CALL  FILGET(FINPAT,FINPIC,NLET,IR(1,I),INTFLG)
	   IF (INTFLG.NE.0) RETURN

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,NSAM,NROW,NSL,
     &               MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

	   GR=IR(2,I)
	   N2=NSAM/4
           K_Q=1
           K_2=IPALIGN64(K_Q+NSAM*NROW)
           K_3=IPALIGN64(K_2+NSAM*NROW)
           K_4=IPALIGN64(K_3+NSAM)
           K_5=IPALIGN64(K_4+NROW)
           K_6=IPALIGN64(K_5+NSAM)
           K_7=IPALIGN64(K_6+N2)
           K_8=IPALIGN64(K_7+N2)
           K_9=IPALIGN64(K_8+NSAM)

C          MEMORY REQUIREMENT
           MEMREQ=K_9
           IF (MEMREQ.GT.MAXMEM) THEN
	       CALL ERRT(6,'AT SA',NE)
               RETURN
           ENDIF

           DO IJK = 1,NROW
              CALL REDLIN(INPIC,Q(1+((IJK-1))*NSAM),NSAM,IJK)
           ENDDO

           LENGTH=7
           LENGTH=LENGTH/2+1
           K=LENGTH*LENGTH
C          MAKE SURE THAT LENGTH*LENGTH<NSAM ! OTHERWISE CHANGE XX.
           CALL MEED(Q(K_Q),Q(K_5),NSAM,LENGTH,Q(K_2),K)
C
           CALL TOMA(NSAM,NROW,NNSAM,NNROW,Q(K_2),VART,MVAR)
           CALL POJ(NSAM,Q(K_2),Q(K_3),Q(K_4),Q(K_5),Q(K_6),
     &           Q(K_7),Q(K_8),N2,VV)

           WRITE(4,299)GR,(VART(LJK),LJK=1,MVAR),VV
299        FORMAT(2X,F3.1,4(1X,1PE14.7),/,4(1X,1PE14.7),/,1(1X,1PE14.7))

           CLOSE(INPIC)
        ENDDO
c100	CONTINUE
	CLOSE(4)

        CALL FILERD(FIL6,NLET6,DATEXC,'ANALYSIS RESULTS',IRTFLG)
        IF (FIL6(1:1).EQ.'*') RETURN
	OPEN(LUN51,FILE=FIL6,STATUS='UNKNOWN')

        CALL FILERD(FIL5,NLET5,DATEXC,'DISCRIM FUNCTION',IRTFLG)
        IF (FIL5(1:1).EQ.'*') RETURN
  	OPEN(LUN50,FILE=FIL5,STATUS='UNKNOWN',FORM='UNFORMATTED')

        OPEN(4,FILE='SCR',STATUS='OLD')

	CALL DISC(LUN50,LUN51)

	CLOSE(4,STATUS='DELETE')
	END
