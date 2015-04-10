C++*********************************************************************
C
C MRRT.FOR                                      06/12/96
C                  MAXNAM                        JUL 14   ARDEAN LEITH
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
C MRRT
C
C PURPOSE:
C   ROTATE AND SCALE 3D COORDINATES STORED IN A DOCUMENT FILE  06/12/96
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE MRRT(MUNUSED)

         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC'

         COMMON           DUMMY(1024),FI1(3),B(4)
         DOUBLE PRECISION R(3,3)

         CHARACTER(LEN=MAXNAM):: DOCFIL

         CHARACTER*1  NULL

         DATA  NDOC,NDOUT/55,56/

         NULL = CHAR(0)

         CALL FILERD(DOCFIL,NLETI,NULL,'3D COORDINATES DOCUMENT',IRTFLG)
         IF (IRTFLG .EQ. -1)  RETURN

         CALL  RDPRM2
     &      (FI1(1),FI1(2),NOT_USED,'ROTATION ANGLES PSI, THETA')

         CALL RDPRM(FI1(3),NOT_USED,'ROTATION ANGLE PHI')

         CALL RDPRM(SCALE,NOT_USED,'SCALE')

         CALL BLDR(R,FI1(1),FI1(2),FI1(3))

         NLIST = 4
         K     = 0
         K2    = 1

778      LERR = -1
         KP1  = K+1

         CALL UNSAV(DOCFIL,K,NDOC,KP1,FI1,3,LERR,K2)

         IF (LERR .EQ. 0)  THEN
            K = K+1
            DO KKK=1,3
               B(KKK+1)=0.0
               DO L=1,3
                  B(KKK+1) = B(KKK+1)+FI1(L)*R(KKK,L)
	       ENDDO
	    ENDDO
            DO KKK=1,3
               B(KKK+1) = B(KKK+1)*SCALE
	    ENDDO
            B(1) = K
            CALL  SAVD(NDOUT,B,NLIST,IRTFLG)
            IF (IRTFLG .EQ. -1) GOTO 5
            GOTO  778

         ENDIF

5        CALL  SAVDC
         CLOSE(NDOUT)
         CLOSE(NDOC)

         END
