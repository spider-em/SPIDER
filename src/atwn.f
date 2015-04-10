C++*********************************************************************
C
C  ATWN.F                               
C                  OPFILEC                         FEB 03 ARDEAN LEITH
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
C   ATWN(MAXMEM)
C
C--*********************************************************************

        SUBROUTINE ATWN(MAXMEM)

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'
C       IN CMLIMIT::  COMMON /IMGMAX/ INUMBR(NIMAX)

        PARAMETER (MVAR=8)
        COMMON     Q(1)

        CHARACTER(LEN=MAXNAM) ::  FINPAT,FINPIC
        COMMON /FPIC/             FINPAT,FINPIC,NLETW

        COMMON /DIMSS1/ K_Q,K_P,K_D,K_3,K_4,K_5,K_6,K_7,K_8,K_9

        DATA INPIC/55/

C       OPEN SAMPLE INPUT FILE TO GET SIZING INFO
        NMAX = NIMAX
	CALL  FILSEQP(FINPAT,NLETW,INUMBR,NMAX,NIMA,
     &    'TEMPLATE FOR 2-D IMAGE NAME',IRTFLG)

C       NIMA - TOTAL NUMBER OF IMAGES
        CALL  FILGET(FINPAT,FINPIC,NLETW,INUMBR(1),INTFLG)

        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,NSAM1,NROW1,NSL,
     &               MAXIM,' ',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CLOSE(INPIC)   
	NW = NROW1

        N2=NW/4
        K_Q=1
        K_P=IPALIGN64(K_Q+NSAM1*NW)
        K_D=IPALIGN64(K_P+NW*NW)
        K_3=IPALIGN64(K_D+NW*NW)
        K_4=IPALIGN64(K_3+NW)
        K_5=IPALIGN64(K_4+NW)
        K_6=IPALIGN64(K_5+N2)
        K_7=IPALIGN64(K_6+N2) 
        K_8=IPALIGN64(K_7+NW)
        K_9=IPALIGN64(K_8+NW)

        MEMREQ = IPALIGN64(K_9+MVAR)
        IF (MEMREQ .GT. MAXMEM) THEN
           CALL ERRT(6,'AT WN',NE)
           RETURN
        ENDIF

	CALL WPDP(INUMBR,Q(K_Q),Q(K_P),Q(K_9),
     &            NSAM1,NROW1,NW,N2,MVAR,NIMA)

        END
