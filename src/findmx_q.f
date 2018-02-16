C++*********************************************************************
C
C FINDMX_Q.F    COSMETIC                        NOV  2017 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2017  Health Research Inc.,                         *
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
C FINDMX_Q(D,LSD,NX,NY,NSI,PEAKV,SX,SY)
C
C PURPOSE: FIND MAXIMUM LOCATION AND VALUE IN IMAGE: 
C              D (NXxNY) ACTUAL SIZE: (LSDxNY).  ONLY SEARCHES NSI
C              PIXELS FROM CENTER OF IMAGE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE FINDMX_Q(D,LSD,NX,NY,NSI,PEAKV,SX,SY)

         REAL   :: D(LSD,NY), Z(-1:1,-1:1)

         JC    = NY / 2+1
         IC    = NX / 2+1
         PEAKV = D(IC,JC)

         DO JT= -NSI,NSI
            J = JT + JC
            DO IT = -NSI,NSI
               I = IT + IC
               IF (PEAKV <= D(I,J))  THEN
                  PEAKV = D(I,J)
                  IX    = I
                  IY    = J
               ENDIF
            ENDDO         
         ENDDO

         SX = IX - IC
         SY = IY - JC

         IF(IY < 2 .OR. IY > NY-1 .OR. IX < 2 .OR. IX > NX-1) RETURN

         DO J=-1,1
            DO I = -1,1
               Z(I,J) = D(IX+I,IY+J)
            ENDDO
         ENDDO

         CALL PARABL(Z,XSH,YSH,PEAKV)

         SX = SX + XSH
         SY = SY + YSH

         END
