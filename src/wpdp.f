
C ++********************************************************************
C                                                             
C  WPDP                                                              
C                  LEN=MAXNAM                      JUL 14 ARDEAN LEITH  
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
C                                                                      *
C  WPDP(ILIST,QQ,WWN,VART,NSAM1,NROW1,NW,N2,MVAR,NIMA)                 *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE WPDP(ILIST,NQQ,WWN,VART,NSAM1,NROW1,NW,N2,MVAR,NIMA)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

	DIMENSION  ILIST(*)
        COMMON     Q(1)

        !CHARACTER*80  FINPAT,FINPIC,FILD
        !COMMON /FPIC/ FINPAT,FINPIC,NLETW

        CHARACTER(LEN=MAXNAM) ::  FINPAT,FINPIC,FILD
        COMMON /FPIC/             FINPAT,FINPIC,NLETW

        COMMON /DIMSS1/ K_Q,K_P,K_D,K_3,K_4,K_5,K_6,K_7,K_8,K_9

        INTEGER*2        M1,MD,NMAX,JV(9)
        DOUBLE PRECISION TMEAN(9),D(9,2),AR(3,2),E(3)
	PARAMETER        (NLIST=2)
	DIMENSION        DLIST(NLIST)
        DIMENSION        WWN(NSAM1,NROW1),QQ(NSAM1,NW),VART(MVAR)
	LOGICAL          FOUND

        DATA INPIC/55/,NDOC/38/,LUN50/50/

        CALL RDPRMI(NNSAM,NNROW,NOT_USED,'SIZE OF MINI WINDOW')

        CALL FILERD(FILD,NLET5,DATEXC,'DISCRIMINANT FUNCTION',IRTFLG)
        IF (FILD(1:1) .EQ. '*') RETURN

        OPEN(LUN50,FILE=FILD,STATUS='OLD',FORM='UNFORMATTED')  
 
C       READ THE DISCRIMINANT FUNCTION
        READ(LUN50)M1,MD,NMAX
        READ(LUN50)(TMEAN(K),K=1,M1)
        DO  JJ=1,M1
           READ(LUN50) (D(JJ,K),K=1,MD)
        ENDDO
        DO  IJ=1,NMAX
           READ(LUN50) (AR(IJ,K),K=1,MD)
        ENDDO
        READ(LUN50) (E(IJ),IJ=1,NMAX)
        READ(LUN50) (JV(IJ),IJ=1,M1)
	CLOSE(LUN50)

	FOUND = .FALSE.
	NI    = 0
	DO  K1=1,NIMA
           CALL  FILGET(FINPAT,FINPIC,NLETW,ILIST(K1),INTFLAG)
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,
     &               LSAM,LROW,NSLICE,
     &               MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           DO K2=1,LROW
              CALL REDLIN(INPIC,WWN(1,K2),LSAM,K2)
           ENDDO
           CLOSE(INPIC)

C          DO A MEDIAN FILTRATION AND THEN CALCULATE STATISTICS

           LENGTH = 7
           LENGTH = LENGTH/2+1
           K      = LENGTH*LENGTH
           CALL MEED(WWN,Q(K_3),NW,LENGTH,Q(K_D),K)

           CALL TOMA(NW,NW,NNSAM,NNROW,Q(K_D),VART,MVAR)      
           GR = 0.0        

           CALL POJ(NW,Q(K_D),Q(K_7),Q(K_8),Q(K_3),Q(K_5),Q(K_6),
     &           Q(K_4),N2,VV)

           CALL CLSS(M1,MD,NMAX,TMEAN,D,AR,E,JV,GR,VART,MVAR,VV,LL)

           IF (LL .EQ. 1)THEN
              NI       = NI+1 
              DLIST(1) = NI
              DLIST(2) = ILIST(K1)
              CALL SAVD(NDOC,DLIST,NLIST,IRTFLG)
           ENDIF
           FOUND = .TRUE.
	ENDDO

	IF (FOUND) THEN
	   CALL SAVDC
	   CLOSE(NDOC)
	ELSE
	   WRITE(NOUT,*)' No good particles found.'
	ENDIF

        END
