C++*********************************************************************
C
C    VOMD.F  NEW                                               07/12/93
C            NEW FORMAT OF 'AP MD' OUTPUT DOC FILE             08/05/99
C            PSI INVERTED
C            MIRRORED IMAGES HAVE NEGATIVE SIGN                10/21/99
C            UNASSIGNED IMAGES MARKED BY 0
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C  PURPOSE: CREATE ANGULAR DOCUMENT FILE FROM 'AP MD' OUTPUT 
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE VOMD

         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC'

         COMMON  AAA(6,1)

         REAL                   :: X(6),Y(0:3)
         CHARACTER (LEN=MAXNAM) :: DOCFIL
         CHARACTER (LEN=1)      :: NULL = CHAR(0)

         INTEGER,PARAMETER      :: NDOC  = 81
         INTEGER,PARAMETER      :: NDOUT = 82

         CALL SET_MPI(ICOMM,MYPID,MPIERR)  

         CALL FILERD(DOCFIL,NLETI,NULL,'ANGULAR DOCUMENT',IRTFLG)
         IF (IRTFLG.EQ.-1)  RETURN

C        FILLS AAA WITH 3 ANGLES PER LINE
         K    = 0
         K2   = 1
778      LERR = -1
         KP1  = K+1
         CALL  UNSAV(DOCFIL,K,NDOC,KP1,AAA(1,KP1),3,LERR,K2)
         IF (LERR == 0)  THEN
            K = K+1
            GOTO 778
         ENDIF
         CLOSE(NDOC)

         CALL FILERD
     &      (DOCFIL,NLETI,NULL,'OUTPUT FROM AP ** DOCUMENT',IRTFLG)
         IF (IRTFLG == -1)  RETURN

         NLIST = 4
         K     = 0
         K2    = 1

779      LERR  = -1
         KP1   = K+1

         CALL  UNSAV(DOCFIL,K,NDOC,KP1,X,6,LERR,K2)

         IF (LERR == 0)  THEN
	    K  = K+1
            MM = X(1)
            IF (MM .NE. 0)  THEN
               Y(1) = -X(3) + 360.0

               IF (MM < 0) THEN
                  MM = -MM
                  Y(1) = Y(1)+180.0
                  Y(2) = 180.0 - AAA(2,MM)
                  Y(3) = AAA(3,MM)+180.0
               ELSE
                  Y(2) = AAA(2,MM)
                  Y(3) = AAA(3,MM)
               ENDIF

               IF (Y(1) >= 360.0)  Y(1) = Y(1)-360.0
               IF (Y(3) >= 360.0)  Y(3) = Y(3)-360.0
               Y(0) = X(6)

               CALL SAVD(NDOUT,Y,NLIST,IRTFLG)
               IF (IRTFLG == -1) GOTO  5

            ENDIF
            GOTO  779
         ENDIF

5        CALL  SAVDC
         IF (MYPID <= 0) CLOSE(NDOUT)
         IF (MYPID <= 0) CLOSE(NDOC)
         END
