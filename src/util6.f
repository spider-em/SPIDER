
C **********************************************************************
C
C   UTIL6.F                                  AUTHOR: ArDean Leith
C                           ADDED SURFFIT      MAR  00 ARDEAN LEITH
C                           'LA' ADDED         OCT  02 ARDEAN LEITH
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
C   UTIL6(MAXDIM)
C
C   PURPOSE: HANDLES OPERATIONS FOR EXTRACT CLUSTERS (EC) AND IA
C	     THEY DEAL WITH CLUSTER CONNECTIVITY VOLUMES. AND SO
C            WHICH ANALYZES SURFACE FITTING BETWEEN TWO VOLUMES.
C
C--*********************************************************************
             
        SUBROUTINE UTIL6(MAXDIM)
		                                                                                        
        INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=MAXNAM)::    FILNAM
        CHARACTER(LEN=1)     ::    NULL
        LOGICAL              ::    BOTTEM

C	MAXIMUM NUMBER OF REGISTERS PER KEY IN 'UD IC'
	PARAMETER (MAXREG=7)
    
C	MAXIMUM NUMBER OF KEYS IN 'UD IC'
	PARAMETER (MAXKEY=4000)
 

        MAXIM1 = 0
        MAXIM2 = 0
        MAXIM3 = 0
	LUN1   = 8
	LUN2   = 12
	LUN3   = 7

        IF (FCHAR(1:2) .EQ. 'IA')  THEN
C          SURFACE AREA FROM VOLUME ------------------------------- IA

           CALL SURFAREA(MAXDIM)

	ELSEIF (FCHAR(1:2) .EQ. 'EC') THEN

	 IF (FCHAR(4:8) .EQ. 'CLOSE') THEN
C          EXTRACT CLUSTERS FROM VOLUME ---------------------- EC CLOSE

            CALL MAPDIST(.TRUE.,.TRUE.,IRTFLG)

         ELSEIF (FCHAR(4:5) .EQ. 'SE') THEN
C          EXTRACT CLUSTERS FROM VOLUME ------------------------- EC SE

           CALL MAPDOC(IRTFLG)

         ELSEIF (FCHAR(4:5) .EQ. 'ST') THEN
C          EXTRACT CLUSTERS FROM VOLUME ------------------------- EC ST

           CALL IMSTAT(IRTFLG)

         ELSEIF (FCHAR(4:5) .EQ. 'FI') THEN
C          EXTRACT CLUSTERS FROM VOLUME ------------------------- EC FI

           CALL MAPFILT(IRTFLG)

         ELSEIF (FCHAR(4:5) .EQ. 'CL')  THEN
C          EXTRACT CLUSTERS FROM VOLUME ------------------------- EC CL
           CALL CONINT(IRTFLG)
	 ENDIF

	ELSEIF (FCHAR(1:2) .EQ. 'FV') THEN
C          FIND THRESHOLD RESULTING IN A SPECIFIED VOLUME --------- FV
           CALL  FV

	ELSEIF (FCHAR(1:2) .EQ. 'SO' .AND. FCHAR(4:4) .EQ. 'C') THEN
C          FIND SURFACE DIFFERENCE ------------------------------ SO C
           CALL  SURFCOMP()

	ELSEIF (FCHAR(1:2) .EQ. 'SO') THEN
C          FIND SURFACE DIFFERENCE -------------------------------- SO
           CALL  SURFFIT()

	ELSEIF (FCHAR(1:2) .EQ. 'LA') THEN
C          LABEL SPIDER IMAGE WITH A LABEL ------------------------- LA

C          OPEN INPUT FILE AND FIND MIN & MAX
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,
     &             NSAM1,NROW1,NSLICE1, MAXIM1,'INPUT',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

           IF (IMAMI .NE. 1)
     &         CALL NORM3(LUN1,NSAM1,NROW1,NSLICE1,FMAX,FMIN,AV)
           FMIN1 = FMIN
           FMAX1 = FMAX

C          OPEN NEW OUTPUT FILE
C          NROWF IS DEPENDENT ON YOUR FONT IMAGE (KLUDGE)
           NROWF = 36

           BOTTEM = (FCHAR(4:4) .NE. 'L')
           NROW2  = NROW1
           IF (BOTTEM) NROW2 = NROW1 + NROWF

           CALL OPFILEC(0,.TRUE.,FILNAM,LUN3,'U',IFORM,
     &              NSAM1,NROW2,NSLICE1,MAXIM3,'OUTPUT',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

           CALL LABELSPI(LUN1,LUN2,LUN3,NSAM1,NROW1,NSLICE1,
     &                   FMIN1,FMAX1,BOTTEM)

9000       CLOSE(LUN1)
           CLOSE(LUN2)
           CLOSE(LUN3)

        ENDIF

	END

