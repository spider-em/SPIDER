
C ++********************************************************************
C                                                                      
C   MRQUATER                                                            
C          PROMPTS & DOC FILE HEADERS IMPROVED   FEB 2014 ARDEAN LEITH                                                                *
C                                                                    
C                                                                      
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C   MRQUATER                                                            
C                                                                
C   PURPOSE: COMPUTES ROTATION MATRIX WHICH TRANSFORMS
C            THE POINTS OF IPT ONTO THE POINTS IN RPT.
C
C   USES UNIT QUATERIONS. THEORY FOUND IN APRIL 1987 J.OPT.SOC.AM.A
C     PP 629-642 BY BERNOLD K.P. HORN
C
C   ALL THE POINTS TO BE USED MUST BE FOUND IN BOTH VIEWS.
C
C   PARAMETERS:
C           RPT(3,LS) = COORDS OF POINTS IN REFERENCE VIEW
C           VPT(3,LS) = COORDS OF POINTS IN VIEW TO BE ALIGNED
C                       COORDS SHOULD BE CENTERED ABOUT THE CENTER 
C                       OF MASS OF THE POINTS USED.
C           NTPT      = NUMBER OF POINTS
C           LS        = RPT & VPT DIMENSION
C
C   OUTPUT: ROT(3,3) = ROTATIONAL MATRIX TO MATCH VPT TO RPT
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE MRQUATER(RPT, VPT, ROT, NTPT, LS)

      INCLUDE     'CMBLOCK.INC'

      DIMENSION ROT(3,3), QUAT(4,1), EVAL(1),
     &         S(3,3), TRIX(4,4), RPT(3,LS), VPT(3,LS)

C     TRIX IS MATRIX TO FIND QUATERION (QUAT)

C     SUM UP PRODUCTS
      DO   J = 1,3
        DO   M = 1,3
          S(M,J) = 0.0
          DO   I = 1,NTPT
  	      S(M,J) = S(M,J)+VPT(M,I)*RPT(J,I)
	  ENDDO
	ENDDO
      ENDDO
C
      TRIX(1,1) = S(1,1)+S(2,2)+S(3,3)
      TRIX(1,2) = S(2,3)-S(3,2)
      TRIX(1,3) = S(3,1)-S(1,3)
      TRIX(1,4) = S(1,2)-S(2,1)
      TRIX(2,1) = TRIX(1,2)
      TRIX(2,2) = S(1,1)-S(2,2)-S(3,3)
      TRIX(2,3) = S(1,2)+S(2,1)
      TRIX(2,4) = S(3,1)+S(1,3)
      TRIX(3,1) = TRIX(1,3)
      TRIX(3,2) = TRIX(2,3)
      TRIX(3,3) = -S(1,1)+S(2,2)-S(3,3)
      TRIX(3,4) = S(2,3)+S(3,2)
      TRIX(4,1) = TRIX(1,4)
      TRIX(4,2) = TRIX(2,4)
      TRIX(4,3) = TRIX(3,4)
      TRIX(4,4) = -S(1,1)-S(2,2)+S(3,3)

      CALL MREIGEN (TRIX,EVAL,QUAT,KODE)

      IF (KODE .NE. 0) THEN
         WRITE(NOUT,250)KODE
 250     FORMAT('  POOR FIT: IVIEW# ',I3,',    KODE= ',I3)

         WRITE(NDAT,260)(QUAT(J,1),J=1,4)
 260     FORMAT('  QUAT:  ( ',F10.7,2X,F10.7,2X,F10.7,2X,F10.7,' )' )
      ENDIF

C     INSURE IT IS A UNIT QUATERION:
      QLNG1=(QUAT(1,1)**2 + QUAT(2,1)**2 + QUAT(3,1)**2 + QUAT(4,1)**2)
      QLNG2     = SQRT(QLNG1)
      QUAT(1,1) = QUAT(1,1)/QLNG2
      QUAT(2,1) = QUAT(2,1)/QLNG2
      QUAT(3,1) = QUAT(3,1)/QLNG2
      QUAT(4,1) = QUAT(4,1)/QLNG2
C
      ROT(1,1) = QUAT(1,1)**2+QUAT(2,1)**2-QUAT(3,1)**2
     &           -QUAT(4,1)**2
      ROT(1,2) = 2*(QUAT(2,1)*QUAT(3,1)-QUAT(1,1)*QUAT(4,1))
      ROT(1,3) = 2*(QUAT(2,1)*QUAT(4,1)+QUAT(1,1)*QUAT(3,1))
      ROT(2,1) = 2*(QUAT(3,1)*QUAT(2,1)+QUAT(1,1)*QUAT(4,1))
      ROT(2,2) = QUAT(1,1)**2-QUAT(2,1)**2+QUAT(3,1)**2
     &           -QUAT(4,1)**2
      ROT(2,3) = 2*(QUAT(3,1)*QUAT(4,1)-QUAT(1,1)*QUAT(2,1))
      ROT(3,1) = 2*(QUAT(4,1)*QUAT(2,1)-QUAT(1,1)*QUAT(3,1))
      ROT(3,2) = 2*(QUAT(4,1)*QUAT(3,1)+QUAT(1,1)*QUAT(2,1))
      ROT(3,3) = QUAT(1,1)**2-QUAT(2,1)**2-QUAT(3,1)**2
     &            +QUAT(4,1)**2
      END
