C++*********************************************************************
C
C $$ TFD.FOR
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
C
C CHANGED 9/5/94 TO INCLUDE INTELLIGIBLE
C GAUSSIAN PARAMETER JF
C
C CHANGED LAST LINE TO CORRECT NORMALIZATION 4/2/03
C    WGH passed in by calling routine as ATAN(WGH/(1.0-WGH)),
C    CS1 passed in by calling routine as CS*1.E7,
C    (16.*ALOG(2.)) changed to parameter KPP,
C    F = -QUADPI**2 changed to a parameter
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE TFD(B,CS1,DZ,LAMBDA,Q,DS,IE,AK,WGH,ENV)

         REAL LAMBDA,KAPPA
         PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
         PARAMETER (KPP = 11.09035491943359380000)
         PARAMETER (F = -9.8696050643920898400000)

C         CS1 = CS*1.E7
         F1=1./SQRT(CS1*LAMBDA)
         F2=SQRT(SQRT(CS1*LAMBDA**3))
         Q1 = (Q*F2)**2
         ENV1=ENV/F2**2
         AKK=AK*F2
C         F=-QUADPI**2
         DS1=DS*DS*F1
         KAPPA=DS1*F/KPP
         DZ1=DZ*F1
         P=AKK**3-DZ1*AKK
         CH=EXP(AKK**4*KAPPA)
C         B=(EXP(F*Q1*P**2)*CH)*2*EXP(-ENV1*AKK**2)
C Factor "2" was removed to have the CTF=W at zero.
         B=(EXP(F*Q1*P**2)*CH)*EXP(-ENV1*AKK**2)
         IF(IE.EQ.1) RETURN
         QQT=2.*QUADPI*(.25*AKK**4-.5*DZ1*AKK**2)
         B=B*SIN(QQT-WGH)
C         B=B*SIN(QQT-ATAN(WGH/(1.0-WGH)))       ; old
C         B=B*((1.0-WGH)*SIN(QQT)-WGH*COS(QQT))  ; older
         END

