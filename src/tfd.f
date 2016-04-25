C++*********************************************************************
C
C TFD.F 
C        CHANGED TO INCLUDE INTELLIGIBLE GAUSSIAN PARAMETER JF 9/5/94
C        CHANGED LAST LINE TO CORRECT NORMALIZATION            4/2/03
C        ACR PASSED IN BY CALLING ROUTINE AS ATAN(ACR/(1.0-ACR)),
C        CS1 PASSED IN BY CALLING ROUTINE AS CS*1.E7,
C           (16.*ALOG(2.)) CHANGED TO PARAMETER KPP,
C        F = -QUADPI**2 CHANGED TO A PARAMETER
C        ADDED TFD_PLUS                                     AL 11/9/15
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2015  Health Research Inc.,                         *
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
C TFD     (B,CS1,DZ,LAMBDA,Q,DS,IE,AK,WGH,ENV)
C TFD_PLUS(B_RAW,CS1,DZ,LAMBDA,Q,DS,AK,WGH,ENV, B_ENV)
C
c PARAMETERS:
C          B      = TRANSFER                                  (OUTPUT)
C          CS1    = ADJUSTED SPHERICAL ABBERATION              (INPUT) 
C          DZ     = ADJUSTED DEFOCUS                           (INPUT)
C          LAMBDA = ADJUSTED ELECTRON WAVELENGHT               (INPUT) 
C          Q      = ADJUSTED SOURCE SIZE                       (INPUT)
C          DS     = ADJUSTED DEFOCUS SPREAD                    (INPUT)
C          AK     = ADJUSTED LOCATION                          (INPUT)
C          WGH    = ADJUSTED AMPLITUDE CONTRAST RATIO   == ACR (INPUT)
C          ENV    = ADJUSTED GAUSSIAN ENVELOP HALFWIDTH == GEH (INPUT)
C
C          B_RAW  = TRANSFER                                  (OUTPUT)
C          B_ENV  = ENVELOPE TRANSFER                         (OUTPUT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE TFD(B,CS1,DZ,LAMBDA,Q,DS,IE,AK,WGH,ENV)

         IMPLICIT  NONE
         REAL      :: B,CS1,DZ,LAMBDA,Q,DS,AK,WGH,ENV
         INTEGER   :: IE

         REAL      :: F1,F2,Q1,ENV1,AKK,DS1,KAPPA,DZ1,P,CH,QQT

         DOUBLE PRECISION, PARAMETER :: 
     &            QUADPI = 3.141592653589793238462643383279502884197
         DOUBLE PRECISION, PARAMETER :: KPP = 11.09035491943359380000
         DOUBLE PRECISION, PARAMETER :: F = -9.8696050643920898400000
C                                       F = -QUADPI**2

C        cs1   = cs * 1.e7   ! DONE IN CALLER
C        env   = 1 / geh**2  ! DONE IN CALLER

         F1    = 1. / SQRT(CS1 * LAMBDA)
         F2    = SQRT(SQRT(CS1 * LAMBDA**3))
         Q1    = (Q * F2)**2
         ENV1  = ENV / F2**2
         AKK   = AK * F2

         DS1   = DS  * DS * F1
         KAPPA = DS1 * F / KPP
         DZ1   = DZ  * F1
         P     = AKK**3 - DZ1 * AKK
         CH    = EXP(AKK**4 * KAPPA)

C        b     = (exp(f*q1*p**2)*ch)*2*exp(-ENV1*akk**2)
c        factor "2" was removed to have the ctf = w at zero.

         B     = (EXP(F * Q1 * P**2) * CH) * EXP(-ENV1 * AKK**2)
         IF (IE == 1) RETURN   ! FOR ENVELOPE ONLY

         QQT   = 2. * QUADPI*(.25 * AKK**4 - .5 * DZ1 * AKK**2)
         B     = B * SIN(QQT - WGH)

c        b     = b * sin(qqt - atan(WGH / (1.0 - WGH)))       ; old
c        b     = b * ((1.0 - WGH) * sin(qqt) - WGH*cos(qqt))  ; older

         END

C        ------------------ TFD_PLUS --------------------------------


         SUBROUTINE TFD_PLUS(B_STRAIGHT,CS1,DZ,LAMBDA,Q,DS,AK,WGH,ENV,
     &                       B_ENV)

         REAL      :: LAMBDA,KAPPA
         REAL      :: B,CS1,DZ,Q,DS,AK,WGH,ENV,B_ENV,B_STRAIGHT

         REAL      :: F1,F2,Q1,ENV1,AKK,DS1,DZ1,P,CH

         PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
         PARAMETER (KPP = 11.09035491943359380000)
         PARAMETER (F = -9.8696050643920898400000)

         F1    = 1. / SQRT(CS1 * LAMBDA)
         F2    = SQRT(SQRT(CS1 * LAMBDA**3))
         Q1    = (Q * F2)**2
         ENV1  = ENV / F2**2
         AKK   = AK  * F2

         DS1   = DS  * DS * F1
         KAPPA = DS1 * F / KPP
         DZ1   = DZ  * F1
         P     = AKK**3 - DZ1 * AKK
         CH    = EXP(AKK**4 * KAPPA)

         ! Use source size (q and q1) of zero for b_env? 

         B_ENV = (EXP(F * Q1 * P**2) * CH) * EXP(-ENV1 * AKK**2)

         QQT   = 2. * QUADPI * (.25 * AKK**4 - .5 * DZ1 * AKK**2)

         B_STRAIGHT = B_ENV  * SIN(QQT - WGH)

         END

         !write(6,*) 'env,f2,-env1,expenv: ',env,f2,-env1,exp(-env1*akk**2)



