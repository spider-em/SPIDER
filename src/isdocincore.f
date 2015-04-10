

C++*********************************************************************
C
C ISDOCINCORE.F      NEW                            JUNE 03 ARDEAN LEITH
C                    STRIP EXTENSION IF PRESENT     JUNE 04 ARDEAN LEITH
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
C   ISDOCINCORE(DOCNAM,NIC,MT,IRTFLG)
C
C   PURPOSE:  FINDS OUT IF DOCNAM IS CURRENTLY IN-CORE
C   
C   PARAMETERS:  DOCNAM    FILE NAME WITH EXTENSION              (SENT)
C                NIC       INCORE INDEX FOR DOCNAM               (RET)
C                MT        NEXT EMPTY IC INDEX                   (RET)
C                IRTFLG    ERROR FLAG                            (RET)
C
C--*********************************************************************

        SUBROUTINE ISDOCINCORE(DOCNAM,NIC,MT,IRTFLG)

        USE DOCIC_INFO

        INCLUDE 'CMBLOCK.INC' 

        CHARACTER(LEN=*) ::  DOCNAM

        DATA NICLAST/1/

C       STRIP EXTENSION IF PRESENT
        NDOT = INDEX(DOCNAM,'.',.TRUE.)

C       COMPARE WITH OLD NAMES
	NIC    = 0
        MT     = 0
        IRTFLG = 0

C       NAME IS MOST-LIKELY STILL THE SAME SO CHECK OLD ONE FIRST
        NLEND = LNBLNKN(DOCNAM)

C       STRIP EXTENSION IF PRESENT
        NDOT = INDEX(DOCNAM,'.',.TRUE.)
        IF (NDOT .GT. 1 .AND. NDOT .LT. NLEND) NLEND = NDOT - 1

        NLENO = NLETOLDNAM(NICLAST)
        IF (DOCNAM(1:NLEND) .EQ. OLDNAM(NICLAST)(1:NLENO) .AND.
     &      NLENO  .EQ. NLEND) THEN
C          DOCNAM IS ALREADY IN-CORE
           NLEND = LNBLNKN(DOCNAM)
           NIC  = NICLAST
           RETURN   
        ENDIF

        NLEND = LNBLNKN(DOCNAM)
C       WANT TO USE LOWEST NIC FIRST NOW al
	DO ICORE = MAXICDOCS,1,-1
           NLENO = NLETOLDNAM(ICORE)
           IF (DOCNAM(1:NLEND) .EQ. OLDNAM(ICORE)(1:NLENO) .AND.
     &         NLENO  .EQ. NLEND) THEN
C             DOCNAM IS ALREADY IN-CORE
              NIC     = ICORE
              NICLAST = NIC
              RETURN     
           ENDIF

C          REMEMBER WHICH LOCATIONS ARE EMPTY IN CASE DOCNAM NOT IN-CORE
           IF (OLDNAM(ICORE)(1:1) .EQ. CHAR(0)) MT = ICORE
        ENDDO

	RETURN
	END
