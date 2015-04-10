
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
        SUBROUTINE prjs2
C       --------------------
C
C THIS ROUTINE PERFORMS A PROJECTION FROM IMAGE TO DATA SPACE...
	PARAMETER  (NILMAX=560)
	COMMON     DUMMY(80),BUF(1024),ILIST(NILMAX),
     A    NSAM,NROW,NANG,NN,NMAT,
     1	  LTB,LTBN,K_ANG,K_DM,K_LB,K_MAP,K_IPCUBE,
     2	  K_BCKE,K_PROJ,K_BCN,K_PRN,k_sigma,
     3    KDM(7),
     4	  IUNIT,Q(1)

	CALL  PRJC2
     &  (Q(K_BCKE),NMAT,Q(K_DM),NANG,Q(K_PRN),q(k_ipcube),NN,NSAM)
        END
