
C************************************************  6/23/80 *** VAX 9/15/81
C
C   VTIL3.F                          
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
C    VTIL3
C    
C    DRIVER FOR SOME ANCIENT STUFF AND FOR SPHERICAL DECONVOLUTION
C               MULTIPLE LINEAR REGRESSION FOR FRAME ALIGNMENT
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*******************************************************************

        SUBROUTINE VTIL3(IDUM)

        INCLUDE 'CMBLOCK.INC'

C	DATA FUNC/'TA SP 19 ML'/   ! 19 is SPH

        INTEGER, PARAMETER  :: LUN1 = 20

        SELECT CASE(FCHAR(1:2))
       
          CASE ('TA')
	     CALL TILT(LUN1)
 
          CASE ('19')                   ! 'SPH = 19'
             CALL SPHDECON()    

          CASE ('ML')                  ! MULTIPLE LINEAR REGRESSION
   	     CALL FRAMEALIGN_MLR()

          CASE ('SP')                  ! SPOT DETECTION
   	     CALL DIFF1O(IRTFLG)

          CASE DEFAULT
             CALL  ERRT(101,'UNNOWN OPERATION',NE)

        END SELECT

        CLOSE(LUN1)

        END
