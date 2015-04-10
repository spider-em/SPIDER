
      SUBROUTINE RMREPL(RADVOL1,RADVOL2,NSAM,NROW,NSLICE,FUZZ)
      
C     TRANSFER EACH LINE OF A 3D (OR 2D) RADON, IF COUNTER IS AT LEAST 1
C     TRANSFORM (FORMAT2000) FROM RADVOL1 TO RADVOL2.

C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2000  M. Radermacher                                  *
C=*                                                                    *
C=* Email:  spider@wadsworth.org                                       *
C=*                                                                    *
C=* This program is free software; you can redistribute it and/or      *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* This program is distributed in the hope that it will be useful,    *
C=* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C=* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
C=* General Public License for more details.                           *
C=*                                                                    *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program; if not, write to the                      *
C=* Free Software Foundation, Inc.,                                    *
C=* 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.      *
C=*                                                                    *
C=**********************************************************************


      DIMENSION RADVOL1(NSAM,NROW,NSLICE),RADVOL2(NSAM,NROW,NSLICE)
      
      F1=1.-FUZZ
      F2=FUZZ
      DO I=1,NSLICE
         DO K=1,NROW
           IF(RADVOL1(NSAM,K,I).GE.1) THEN
           DO L=1,NSAM
             RADVOL2(L,K,I)=RADVOL2(L,K,I)*F1+RADVOL1(L,K,I)*F2
           ENDDO
           ENDIF
         ENDDO
       ENDDO
       RETURN
       END
 
