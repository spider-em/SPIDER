C++*******************************************************************
C
C TF.F
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
C  TF(B,CS,DZ,LAMBDA,KM,NSAM,Q,DS,IE)
C
C  CHANGED 9/5/94 TO INCLUDE INTELLIGIBLE GAUSSIAN PARAMETER JF
C
C--*******************************************************************

	SUBROUTINE TF(B,CS,DZ,LAMBDA,KM,NSAM,Q,DS,IE,WGH,ENV)

C       I DO NOT KNOW IF SAVE IS NEEDED FEB 99 al
        SAVE

	REAL B(1),LAMBDA,KM,KAPPA,KM1

	DATA PI/3.14159/

           IF (CS < 0.0001)    CS = 0.0001
	CS1   = CS*1.E7
	F1    = 1./SQRT(CS1*LAMBDA)
	F2    = SQRT(SQRT(CS1*LAMBDA**3))
	KM1   = KM*F2
	DK    = KM1/FLOAT(NSAM-1)
	Q1    = (Q*F2)**2
        ENV1  = ENV/F2**2
	F     = -PI**2
	DS1   = DS*F1
	KAPPA = (DS1**2)*F/(4.0)
  	DZ1   = DZ*F1

	DO K=1,NSAM
	   AK   = FLOAT(K-1)*DK
	   P    = AK**3-DZ1*AK
	   CH   = EXP(AK**4*KAPPA)
           B(K) = (EXP(F*Q1*P**2)*CH)*2*EXP(-ENV1*AK**2)
           IF (IE .NE. 1) THEN
              QQT  = 2.*PI*(.25*AK**4-.5*DZ1*AK**2)
              B(K) = B(K)*((1.0-WGH)*SIN(QQT)-WGH*COS(QQT))
           ENDIF
        END DO

        RETURN
	END

