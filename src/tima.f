C ++********************************************************************
C                                                                      *
C TIMA                                               *
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

        SUBROUTINE TIMA(NSAM,Y,AVG,NNSAM,NNROW,AVAV,AVVR,SDAV,SDVR)

C	NSAM & NROW ARE THE DIMENSIONS OF THE ORIGINAL WINDOW
C	AND NNSAM AND NNROW ARE THE DIMENSIONS OF THE MINI WINDOW) 

        DIMENSION Y(NSAM,NSAM)


	ND   = NNSAM/2
	NT   = 0
	SUMT = 0
	VART = 0
	UUT  = 0
	UVT  = 0
	DO  I=1,(NSAM-NNSAM+1),ND
	   DO  J=1,(NSAM-NNROW+1),ND
              SUME = 0
	      SSUM = 0
	      AAD2 = 0
	      NT   = NT+1
	      DO  IW=I,I+NNSAM-1
	         DO  JW=J,J+NNROW-1
	            T    = Y(IW,JW)
	            R    = T-AVG
	            SUME = SUME+T
	            SSUM = SSUM+R
	            U    = T*T
	            AAD2 = AAD2+U
	         ENDDO
	      ENDDO
	      XX   = FLOAT(NNROW)*FLOAT(NNSAM)
	      AVV  = SUME/XX
	      SQQ  = AVV*AVV                  
	      VAR  = AAD2-(XX*SQQ)
	      AAVG = ABS(SSUM)/XX
	      SUMT = SUMT+AAVG
	      VART = VART+VAR
	      UU   = AAVG*AAVG
	      UV   = VAR*VAR
	      UUT  = UUT+UU
	      UVT  = UVT+UV	   
	   ENDDO
        ENDDO

C	VAR GIVES THE VARIANCE OF THE MINI WINDOW
C	SUMT AND VART GIVE THE SUMS OF THE LOCAL DISTANCE AVERAGES(LDA) 
C	AND VARIANCES OF ALL THE WINDOWS RESPECTIVELY
C	UU AND UV GIVE THE SQUARES OF LDA AND VARIANCE OF EACH WINDOW
C	UUT AND UVT GIVE THE SUMS OF SQUARES OF UU AND UV

        AVAV = SUMT/FLOAT(NT)
	AVVR = VART/FLOAT(NT)
        PS1  = AVAV*AVAV
	PS2  = AVVR*AVVR
	VRAV = UUT-(FLOAT(NT)*PS1)
	VRVR = UVT-(FLOAT(NT)*PS2)
	SDAV = SQRT(VRAV)
	SDVR = SQRT(VRVR)

        RETURN
        END
