cp
C ++********************************************************************
C                                                                      *
C VODA.F                                                               *
C                   OPENDOC PARAMETERS CHANGED   DEC 2000 ARDEAN LEITH *
C                   OPENDOC LUN &LUNDOCWRTDAT    JUL 2003 ARDEAN LEITH *                                                  *
C                   REMOVED SORTI                NOV 2011 ARDEAN LEITH *                                                  *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2011  Health Research Inc.,                         *
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
C                                                                      *
C  VODA   UNDOCUMENTED OPERATION!! al oct 2012                                                             *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE VODA

	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'

	PARAMETER (QUADPI = 3.14159265358979323846264)
	PARAMETER (DGR_TO_RAD = (QUADPI/180))
	PARAMETER (NDLI = 3)

	DIMENSION     DLIST(NDLI)

	REAL, ALLOCATABLE, DIMENSION(:,:) ::     REFA,VECTOR,X
	INTEGER, ALLOCATABLE, DIMENSION(:,:) ::  LISTA
	INTEGER, ALLOCATABLE, DIMENSION(:) ::    KALC,INGROUP
	INTEGER    ::                            X89
	LOGICAL    ::                            NEWFILE
	CHARACTER(LEN=MAXNAM) ::                 FINP,FINPAT,DOCFIL

	DATA  NDOCT/55/

C       REFERENCE ANGLES FILE
        CALL OPENDOC(DOCFIL,.TRUE.,NLETI,NDOCT,NDOC,.TRUE.,
     &              'REFERENCE ANGLES',
     &              .TRUE.,.FALSE.,.FALSE.,NEWFILE,IRTFLG)
	IF (IRTFLG .NE. 0) RETURN

	CALL LUNDOCINFO(NDOC,NUMREF,MAXCOLS,KEYUSED,.FALSE.,IRTFLG)
	IF (IRTFLG .NE. 0) RETURN

C       1 - THETA
C       2 - PHI
	ALLOCATE(REFA(2,NUMREF))

C       READS COL 1 & 2 INTO REFA
        CALL LUNDOCREDSLC(NDOC,.FALSE.,IDUM,REFA, 2,NUMREF,
     &                  .TRUE.,.FALSE., 1,2, 1,NUMREF, 
     &                  NGOTY, MAXGOTY,IRTFLG)
	CLOSE(NDOCT)
        IF (IRTFLG .NE. 0) GOTO 9999

C       CHANGE ALL ANGLES TO RADIANS
	REFA = REFA*DGR_TO_RAD

C       NEW ANGLES FILE
	CALL OPENDOC(DOCFIL,.TRUE.,NLETI,NDOCT,NDOC,.FALSE.,
     &       'ANGLES TO BE ASSIGNED',
     &       .TRUE.,.FALSE.,.FALSE.,NEWFILE,IRTFLG)
	IF (IRTFLG .NE. 0) GOTO 9999

	CALL LUNDOCINFO(NDOC,NUMA,MAXCOLS,KEYUSED,.FALSE.,IRTFLG)
	IF (IRTFLG .NE. 0) GOTO 9999

	ALLOCATE(X(3,NUMA))
        CALL LUNDOCREDSLC(NDOC,.FALSE.,IDUM,X,2,NUMA,
     &                  .TRUE., .FALSE., 1,2, 1,NUMA, 
     &                   NGOTY,MAXGOTY, IRTFLG)
	CLOSE(NDOCT)
        IF (IRTFLG .NE. 0) GOTO 9999

C       CHANGE ALL ANGLES TO RADIANS
	X = X*DGR_TO_RAD

	ALLOCATE(VECTOR(3,NUMA))
	ALLOCATE(KALC(NUMA))
	ALLOCATE(INGROUP(NUMA))

C       LENGTH OF AUXILIARY LIST IS SET TO NUMA
	LENI = NUMA
	ALLOCATE(LISTA(2,LENI))

C       CHANGE TO NORMAL VECTOR
        VECTOR(1,1:NUMA) = COS(X(2,1:NUMA)) * SIN(X(1,1:NUMA)) 
        VECTOR(2,1:NUMA) = SIN(X(2,1:NUMA)) * SIN(X(1,1:NUMA)) 
        VECTOR(3,1:NUMA) = COS(X(1,1:NUMA)) 

	CALL  RDPRM(OVERLAP,NOT_USED,
     &	     'OVERLAP BETWEEN ANGULAR ZONES (IN DEGREES)')

C       CHANGE TO RADIANS, 
C       DIVIDE BY TWO, AS THE BORDER IS +/- OVERLAP/2 IN EACH DIRECTION
	OVERLAP=COS(OVERLAP*DGR_TO_RAD/2.0)

	DO I=1,NUMA
	  X51 = 0.0
	  DO J=1,NUMREF 
	     X84  = COS(REFA(2,J))*SIN(REFA(1,J)) 
	     X85  = SIN(REFA(2,J))*SIN(REFA(1,J))
	     X86  = COS(REFA(1,J)) 
	     DISD = VECTOR(1,I)*X84+VECTOR(2,I)*X85+VECTOR(3,I)*X86
	     IF (DISD .GT. X51) THEN
	        X51 = DISD
	        X89 = J 
	     ENDIF
	  ENDDO
C         BELONGS TO THIS ONE: SMALLEST ANGLE
	  KALC(I) = X89
        ENDDO

	JODI = 0
	IF (OVERLAP .LT. 1.0)  THEN
C          FIND BORDER CASES, THEY HAVE TO BE WITHIN OVERLAP OF AN 
C          OBJECT IN THE CLASS
	   DO I=1,NUMA-1
	      DO J=I+1,NUMA
	         IF (KALC(I) .NE. KALC(J)) THEN
	            IF (DOT_PRODUCT(VECTOR(:,I),VECTOR(:,J)) .GE. 
     &                 OVERLAP) THEN
C                      PUT OBJECT I FROM GROUP KALC(I) ON A BORDER LIST OF GROUP KALC(J)
C                      AND VICE VERSA
C                      CHECK WHETHER ASSIGNEMENT HAS ALREADY BEEN MADE.
                       IF (JODI .GT. 0)  THEN
                          DO K=1,JODI
                            IF (LISTA(1,K).EQ.I .AND. 
     &                        LISTA(2,K).EQ.KALC(J)) GOTO 701
                          ENDDO
                       ENDIF

                      JODI = JODI+1
C                     CHECK THE LENGTH OF THE LIST FOR THE NEXT TWO
	              IF (JODI+1 .GT. LENI) THEN
	                 CALL ERRT(101,'AUXILIARY LIST TOO SHORT!',NE)
                         GOTO 9999
                      ENDIF
                      LISTA(1,JODI) = I
                      LISTA(2,JODI) = KALC(J)

701                   CONTINUE
	              DO K=1,JODI
	                 IF (LISTA(1,K).EQ.J .AND. 
     &                       LISTA(2,K).EQ.KALC(I)) GOTO 702
	              ENDDO
	              JODI=JODI+1
                      LISTA(1,JODI) = J
                      LISTA(2,JODI) = KALC(I)

702	              CONTINUE
                   ENDIF
                ENDIF
	     ENDDO
	   ENDDO  
	ENDIF

C       GET OUTPUT FILE NAME TEMPLATE 
	NMAX = 0
	CALL  FILSEQP(FINPAT,NLET,IDUML,NMAX,IDUM,
     &               'TEMPLATE FOR SELECTION DOC',IRTFLG)

C       SELECT OBJECTS FOR EACH GROUP
	DO J=1,NUMREF
	   LOBJ = 0
	   DO  I=1,NUMA
C             PUT OBJECT FROM GROUP J ON THE LIST
	      IF (KALC(I) .EQ. J)  THEN
	         LOBJ          = LOBJ+1
	         INGROUP(LOBJ) = I
	     ENDIF
	   ENDDO

C          ADD BORDER OBJECTS TO THE GROUP.
	   IF (JODI .GT. 0)  THEN
	     DO  K=1,JODI
	       IF (LISTA(2,K) .EQ. J) THEN
	          LOBJ          = LOBJ+1
	          INGROUP(LOBJ) = LISTA(1,K)
	       ENDIF
	     ENDDO
	   ENDIF

	   IF (LOBJ .GT. 0)  THEN
C             SAVE IN DOC FILE
	      IF (LOBJ.GT.1) CALL SORTI(INGROUP,LOBJ)
	      CALL FILGET(FINPAT,FINP,NLET,J,INTFLAG)
 
              CALL OPENDOC(FINP,.TRUE.,NLET,NDOCT,NDOC,.FALSE.,' ',
     &              .FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
	      IF (IRTFLG .NE. 0) GOTO 9999

	      DO I=1,LOBJ
	         IAP      = MIN0(I-1,1)
	         DLIST(1) = I
	         DLIST(2) = INGROUP(I)
                 CALL LUNDOCWRTDAT(NDOC,IKEY,DLIST,2,IRTFLG)
	     ENDDO
	     CLOSE(NDOCT)
	  ENDIF
	ENDDO 

9999    IF (ALLOCATED(REFA))   DEALLOCATE(REFA)
        IF (ALLOCATED(VECTOR)) DEALLOCATE(VECTOR)
        IF (ALLOCATED(X))      DEALLOCATE(X)
        IF (ALLOCATED(LISTA))  DEALLOCATE(LISTA)
        IF (ALLOCATED(KALC))   DEALLOCATE(KALC)
        IF (ALLOCATED(INGROUP))DEALLOCATE(INGROUP)
        CLOSE(NDOCT)

	END
