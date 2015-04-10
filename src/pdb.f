C **********************************************************************
C *  PDB
C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2004, 2014 P. A. Penczek                                   *
C=*                                                                    *
C=* University of Texas - Houston Medical School                       *
C=* Email:  pawel.a.penczek@uth.tmc.edu                                *
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

	SUBROUTINE PDB

        INCLUDE 'CMBLOCK.INC'

        IF (FCHAR(4:6) == 'CG3') THEN
C          ------------------------- 'PDB GRAVITY'
           CALL PDBCG3

        ELSEIF (FCHAR(4:5) == 'IF') THEN
C          ------------------------- 'GET PDB FILE PARAMETERS'
	   CALL PDBIF

	ELSEIF (FCHAR(4:5) == 'SH') THEN
C          ------------------------- 'PDB SHIFT
	   CALL PDBSH

	ELSEIF (FCHAR(4:7) == 'RT3A') THEN
C          ------------------------- 'PDB ROTATE AROUND AN ARBITRARY POINT
	   CALL PDBRT3A

	ELSEIF (FCHAR(4:7) == 'RT3L') THEN
C          ------------------------- 'PDB ROTATE AROUND A LINE DEFINED by Two points
	   CALL PDBRT3L

	ELSEIF (FCHAR(4:6) == 'RT3') THEN
C          ------------------------- 'PDB ROTATE AROUND AN ARBITRARY POINT
	   CALL PDBRT3

	ELSEIF (FCHAR(4:6) == 'STP') THEN
C          ------------------------- 'COPY SPIDER COORDINATES TO PDB FORMAT'
	   CALL SPTOPDB

	ELSEIF (FCHAR(4:6) == 'PTS') THEN
C          ------------------------- 'COPY PDB FILE TO SPIDER FORMAT'
	   CALL PDBTOSP
	ENDIF

	END

C **********************************************************************

	SUBROUTINE PDBSH

        INCLUDE 'CMLIMIT.INC'
	INCLUDE 'CMBLOCK.INC'
        CHARACTER (LEN=MAXNAM) :: PDBFILE,TRFILE
	CHARACTER *80  RECLIN,HEAD*10, NULL*1,ATOM*4,RESIDUE*3
	LOGICAL        EX,FLAGCELLDIM,FLAGATOM,ZNUM,COMBIND

	DATA  LUN2,LUN3/25,26/

	NATOM  = 1
	LENREC = 0
        CALL OPAUXFILE(.TRUE.,PDBFILE,NULL,LUN2,LENREC,'O',
     &                 'PDB INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CALL OPAUXFILE(.TRUE.,PDBFILE,NULL,LUN3,LENREC,'N',
     &                 'PDB OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
	
	CALL RDPRM3S(SXT,SYT,SZT,NOT_USED,
     &               'X, Y, & Z SHIFT IN ANGSTROMS:',IRTFLG)

C       CHANGE COORDINATE SYSTEM 
 	SX =  SYT
 	SY =  SXT
 	SZ = -SZT

70	READ(LUN2,30) RECLIN
30	FORMAT(A80)

	IF (RECLIN(1:4) .NE. 'ATOM' .AND.
     &      RECLIN(1:4) .NE. 'END'  .AND.
     & 	    RECLIN(1:4) .NE. 'TER' ) THEN
	    WRITE(LUN3,30) RECLIN
	    GOTO 70 

	ELSEIF (RECLIN(1:4) == 'ATOM') THEN
	   READ(RECLIN,75) HEAD,N,ATOM,NULL,RESIDUE,NULL,NR2,
     &                     NULL,XO,YO,ZO,
     &                     OCCUPANCY,TEMPERATURE,N
75	   FORMAT(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,1X,I3)	

	   X = XO + SX
	   Y = YO + SY
	   Z = ZO + SZ

C          REMOVED ALTERATION DEC 07 al
C          ALTERED NOV 07 al AFTER RECEIVING BUG REPORT
C           TEMP = X
C           X    = Y
C           Y    = TEMP
C           Z    = -Z 

	   WRITE(RECLIN(7:11),80) NATOM
80	   FORMAT(I5)

	   NATOM = NATOM + 1
	   WRITE(RECLIN(31:54),85) X,Y,Z
85	   FORMAT(3F8.3)

	   WRITE(LUN3,30) RECLIN
	   GOTO 70

	ELSEIF (RECLIN(1:3) == 'TER') THEN
	   WRITE(RECLIN(7:11),80)NATOM
	   WRITE(LUN3,30)RECLIN
	   NATOM = NATOM+1
	   GOTO 70

        ELSEIF (RECLIN(1:3) == 'END') THEN
	   WRITE(LUN3,30)RECLIN
	   GOTO 999
	ELSE
	   GOTO 70
	ENDIF

999	CLOSE(LUN2)
	CLOSE(LUN3)
	END


C **********************************************************************
C*
C*Rotate PDB file around an arbitrary center
C*
C********************************************

	SUBROUTINE PDBRT3A

	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'
        CHARACTER (LEN=MAXNAM) :: PDBFILE,TRFILE  
      
        CHARACTER *80 RECLIN,HEAD*10, NULL*1,ATOM*4,RESIDUE*3
	LOGICAL       EX,FLAGCELLDIM,FLAGATOM,ZNUM,COMBIND
	DOUBLE PRECISION RM(3,3)

	DATA  LUN2,LUN3/25,26/

	NATOM=1
	LENREC=0
	CALL OPAUXFILE(.TRUE.,PDBFILE,NULL,LUN2,LENREC,'O',
     &                 'PDB INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
 
	CALL OPAUXFILE(.TRUE.,PDBFILE,NULL,LUN3,LENREC,'N',
     &                 'PDB OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

	CALL RDPRM3S(PHI,THETA,PSI,NOT_USED,
     &              'PHI, THETA, PSI:',IRTFLG)

	CALL  BLDR(RM,PSI,THETA,PHI)

	CALL RDPRM3S(XC,YC,ZC,NOT_USED,
     &     'CENTER OF ROTATION IN ANGSTROMS X,Y,Z: ',IRTFLG)

70	READ(LUN2,30) RECLIN
30	FORMAT(A80)

	IF (RECLIN(1:4) .NE. 'ATOM' .AND.
     &      RECLIN(1:4).NE.'END'.AND.
     & 	    RECLIN(1:4).NE.'TER' ) THEN
	   WRITE(LUN3,30) RECLIN
	   GOTO 70 

	ELSEIF (RECLIN(1:4) == 'ATOM') THEN
	   READ(RECLIN,75) HEAD,N,ATOM,NULL,RESIDUE,NULL,NR2,
     &                     NULL,XO,YO,ZO,
     &                     OCCUPANCY,TEMPERATURE,N
75	   FORMAT(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,1X,I3)	
C	   WRITE(*,*) XO,YO,ZO

C          CHANGE COORDINATE SYSTEM
	   W  = XO
	   XO = YO
	   YO = W
	   ZO = -ZO
	   XO = XO-XC
	   YO = YO-YC
	   ZO = ZO-ZC

C          AFTER ROTATION CHANGE IT BACK, THAT'S WHY THE ORDER OF 
C          X,Y IS CHANGED AND Z HAS INVERTED SIGN

	   Y = RM(1,1)*XO+RM(1,2)*YO+RM(1,3)*ZO+XC
	   X = RM(2,1)*XO+RM(2,2)*YO+RM(2,3)*ZO+YC
	   Z = -(RM(3,1)*XO+RM(3,2)*YO+RM(3,3)*ZO+ZC)
	   WRITE(RECLIN(7:11),80) NATOM
80	   FORMAT(I5)
	   NATOM = NATOM+1
	   WRITE(RECLIN(31:54),85)  X,Y,Z
85	   FORMAT(3F8.3)
	   WRITE(LUN3,30)RECLIN
	   GOTO 70

	ELSEIF (RECLIN(1:3) == 'TER') THEN
	   WRITE(RECLIN(7:11),80)NATOM
	   WRITE(LUN3,30)RECLIN
	   NATOM = NATOM + 1
	   GOTO 70

        ELSEIF (RECLIN(1:3) == 'END') THEN
	   WRITE(LUN3,30)RECLIN
	   GOTO 999

	ELSE
	   GOTO 70
	ENDIF

999	CLOSE(LUN2)
	CLOSE(LUN3)
	END


C **********************************************************************
C*
C*Rotate PDB file around a line defined by two points in 3D space
C*
C********************************************

	SUBROUTINE PDBRT3L

	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'

        CHARACTER (LEN=MAXNAM) :: PDBFILE,TRFILE
        CHARACTER *80 RECLIN,HEAD*10, NULL*1,ATOM*4,RESIDUE*3
	LOGICAL          :: EX,FLAGCELLDIM,FLAGATOM,ZNUM,COMBIND
	DOUBLE PRECISION :: RM(3,3),R1(3,3),R2(3,3),R3(3,3)

        REAL, PARAMETER  :: QUADPI = 3.14159265358979323846 
        REAL, PARAMETER  :: DGR_TO_RAD = (QUADPI/180)

	DATA  LUN2,LUN3/25,26/

	NATOM  = 1
	LENREC = 0
	CALL OPAUXFILE(.TRUE.,PDBFILE,NULL,LUN2,LENREC,'O',
     &                 'PDB INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

30	FORMAT(A80)
 
	CALL OPAUXFILE(.TRUE.,PDBFILE,NULL,LUN3,LENREC,'N',
     &                 'PDB OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

	CALL RDPRM1S(ALPHA,NOT_USED,'ROTATION ANGLE ALPHA',IRTFLG)

	CALL RDPRM3S(X1,Y1,Z1,NOT_USED,
     &        'FIRST POINT DEFINING ROTATION AXIS IN ANGSTROMS X,Y,Z',
     &        IRTFLG)

	CALL RDPRM3S(X2,Y2,Z2,NOT_USED,
     &       'SECOND POINT DEFINING ROTATION AXIS IN ANGSTROMS X,Y,Z',
     &       IRTFLG)

	XX    = X2-X1
	YY    = Y2-Y1
	ZZ    = Z2-Z1

#if defined (SP_GFORTRAN)
	PSI   = -ATAN2(YY, XX)                / DGR_TO_RAD
	THETA =  ATAN2(ZZ, SQRT(XX*XX+YY*YY)) / DGR_TO_RAD
#else
	PSI   = -ATAN2D(YY, XX)
	THETA =  ATAN2D(ZZ, SQRT(XX*XX+YY*YY))
#endif

	CALL  BLDR(R1,PSI,THETA,90.0)
	CALL  BLDR(R2,0.0,ALPHA,0.0)

C       -1
C       R   U(ALPHA) R

         DO  11  I=1,3
            DO  11  J=1,3
               R3(J,I)=0.0
               DO  11  K=1,3
11       R3(J,I)=R3(J,I)+R2(K,I)*R1(J,K)

         DO  12  I=1,3
            DO  12  J=1,3
               RM(J,I)=0.0
               DO  12  K=1,3
12       RM(J,I)=RM(J,I)+R1(I,K)*R3(J,K)

70	READ(LUN2,30) RECLIN
	IF (RECLIN(1:4) .NE. 'ATOM' .AND.
     &      RECLIN(1:4) .NE. 'END'  .AND.
     &      RECLIN(1:4) .NE. 'TER' ) THEN
	    WRITE(LUN3,30) RECLIN
	    GOTO 70 

	ELSEIF (RECLIN(1:4) == 'ATOM') THEN
	   READ(RECLIN,75) HEAD,N,ATOM,NULL,RESIDUE,NULL,NR2,
     &                     NULL,XO,YO,ZO,
     &                     OCCUPANCY,TEMPERATURE,N
75	   FORMAT(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,1X,I3)	
C	   WRITE(*,*) XO,YO,ZO

C          CHANGE COORDINATE SYSTEM
	   W  = XO
	   XO = YO
	   YO = W
	   ZO = -ZO
	   XO = XO-X1
	   YO = YO-Y1
	   ZO = ZO-Z1

C          AFTER ROTATION CHANGE IT BACK, THAT'S WHY THE ORDER OF X,Y IS CHANGED
C          AND Z HAS INVERTED SIGN

	   Y=RM(1,1)*XO+RM(1,2)*YO+RM(1,3)*ZO+X1
	   X=RM(2,1)*XO+RM(2,2)*YO+RM(2,3)*ZO+Y1
	   Z=-(RM(3,1)*XO+RM(3,2)*YO+RM(3,3)*ZO+Z1)
	   WRITE(RECLIN(7:11),80) NATOM
80	   FORMAT(I5)
	   NATOM = NATOM+1
	   WRITE(RECLIN(31:54),85)  X,Y,Z
85	   FORMAT(3F8.3)
	   WRITE(LUN3,30)RECLIN
	   GOTO 70

	ELSEIF (RECLIN(1:3) == 'TER') THEN
	   WRITE(RECLIN(7:11),80)NATOM
	   WRITE(LUN3,30)RECLIN
	   NATOM = NATOM + 1
	   GOTO 70

        ELSEIF (RECLIN(1:3) == 'END') THEN
	   WRITE(LUN3,30)RECLIN
	   GOTO 999

	ELSE
	   GOTO 70
	ENDIF

999	CLOSE(LUN2)
	CLOSE(LUN3)
	END

C **********************************************************************
C*
C* Copy SPIDER document file to a PDB file
C*
C********************************************

	SUBROUTINE SPTOPDB

	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'

        REAL     :: DLIST(7)

        CHARACTER (LEN=MAXNAM) :: PDBFILE,TRFILE
 	CHARACTER *95 RECLIN,HEAD*10, NULL*1,ATOM*4
	CHARACTER *95 ANAM*4,RESIDUE*3,PHEAD*5,SENQ*1
	LOGICAL   EX,FLAGCELLDIM,FLAGATOM,ZNUM,COMBIND

	DATA  LUN2,LUN3/25,26/

	OPEN(99,FILE='ttt')

	PHEAD  = 'ATOM '
	LENREC = 0

30	FORMAT(A80)
        CALL OPAUXFILE(.TRUE.,PDBFILE,NULL,LUN2,LENREC,'O',
     &                 'SPIDER DOCUMENT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CALL OPAUXFILE(.TRUE.,PDBFILE,NULL,LUN3,LENREC,'N',
     &                 'PDB OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

	NCOUNT = 0
	DO

  	!READ(LUN2,30,END=701) TRFILE
	!IF (TRFILE(1:2) ==' ;') CYCLE
c33	!FORMAT(I5,1X,I1,7G12.6)
	!READ(TRFILE,33)NATOM,MMM,XO,YO,ZO, CATOM,TEMPERATURE,XSEN,XNUL

        CALL LUNDOCREDNXT(LUN2,IKEY,DLIST,7, UNUSED,ICOUNT,IRTFLG)
        write(6,*) ikey,icount, dlist(1)
	IF (IRTFLG .NE. 0) EXIT

	IF (ICOUNT < 1) CYCLE

        NATOM       = IKEY
        MMM         = ICOUNT
        XO          = DLIST(1)
        YO          = DLIST(2)
        ZO          = DLIST(3)
        CATOM       = DLIST(4)
        TEMPERATURE = DLIST(5)
        XSEN        = DLIST(6)
        XNUL        = DLIST(7)

	INUL        = INT(XNUL)
	IRES        = INT(CATOM/10000.)
	NR2         = INT(XSEN/1000.)
	ISEN        = INT(REAL(MOD(INT(XSEN),1000))/10.)
	OCCUPANCY   = REAL(MOD(INT(XSEN),2))
	
	    IF (IRES == 1) THEN
		RESIDUE ='LEU'
	    ElSEIF (IRES == 2) THEN
		RESIDUE ='THR'
	    ELSEIF (IRES == 3) THEN
		RESIDUE ='GLY'
	    ELSEIF (IRES == 4) THEN
		RESIDUE ='SER'
	    ELSEIF (IRES == 5) THEN
		RESIDUE ='VAL'
	    ELSEIF (IRES == 6) THEN
		RESIDUE ='PHE'
	    ELSEIF (IRES == 7) THEN
		RESIDUE ='ALA'
	    ELSEIF (IRES == 8) THEN
		RESIDUE ='LYS'
       	    ELSEIF (IRES == 9) THEN
		RESIDUE ='ARG'
	    ELSEIF (IRES == 10) THEN
		RESIDUE ='PRO'
            ELSEIF (IRES == 11) THEN
		RESIDUE ='TAU'
	    ELSEIF (IRES == 12) THEN
		RESIDUE ='HIS'
    	    ELSEIF (IRES == 13) THEN
		RESIDUE ='GLU'
	    ELSEIF (IRES == 14) THEN
		RESIDUE ='GLY'
	    ELSEIF (IRES == 15) THEN
		RESIDUE ='TYR'
	    ELSEIF (IRES == 16) THEN
		RESIDUE ='CYS'
  	    ELSEIF (IRES == 17) THEN
		RESIDUE ='TRY'
	    ELSEIF (IRES == 18) THEN
		RESIDUE ='ISO'
	    ELSEIF (IRES == 19) THEN
		RESIDUE ='MET'
	    ELSEIF (IRES == 20) THEN
		RESIDUE ='ILE'
	    ENDIF
cc
 	    IF (ISEN == 1) THEN
		SENQ = 'A'
	    ELSEIF (ISEN == 2) THEN
		SENQ ='B'
	    ELSEIF (ISEN == 3) THEN
		SENQ = 'C'
	    ELSEIF (ISEN == 4) THEN
		SENQ = 'D'
	    ELSEIF (ISEN == 5) THEN
		SENQ = 'E'
	    ELSEIF (ISEN == 6) THEN
		SENQ = 'F'
	    ELSEIF (ISEN == 7) THEN
		SENQ = 'G'
	    ELSEIF (ISEN == 8) THEN
		SENQ = 'H'
     	    ELSEIF (ISEN == 9) THEN
		SENQ = 'I'
	    ELSEIF (ISEN == 10) THEN
		SENQ = 'J'
	    ELSEIF (ISEN == 11) THEN
		SENQ = 'K'
	    ELSEIF (ISEN == 12) THEN
		SENQ = 'L'
	    ELSEIF (ISEN == 13) THEN
		SENQ = 'M'
	    ELSEIF (ISEN == 14) THEN
		SENQ = 'N'
	    ELSEIF (ISEN == 15) THEN
		SENQ = 'O'
	    ELSEIF (ISEN == 16) THEN
		SENQ = 'P'
	    ELSEIF (ISEN == 17) THEN
		SENQ = 'Q'
	    ENDIF
C
	    IF (INUL == 0) THEN
		NULL =' '
	    ELSEIF (INUL == 1) THEN
		NULL ='A'
	    ELSEIF (INUL == 2) THEN
		NULL ='B'
	    ELSEIF (INUL == 3) THEN
		NULL ='C'
	    ELSEIF (INUL == 4) THEN
		NULL ='D'
	    ELSEIF (INUL == 5) THEN
		NULL ='E'
	    ELSEIF (INUL == 6) THEN
		NULL ='F'
	    ELSEIF (INUL == 7) THEN
		NULL ='G'
	    ELSEIF (INUL == 8) THEN
		NULL ='H'
     	    ELSEIF (INUL == 9) THEN
		NULL ='I'
	    ELSEIF (INUL == 10) THEN
		NULL ='J'
	    ELSEIF (INUL == 11) THEN
		NULL ='K'
	    ELSEIF (INUL == 12) THEN
		NULL ='L'
	    ELSEIF (INUL == 13) THEN
		NULL ='M'
	    ELSEIF (INUL == 14) THEN
		NULL ='N'
	    ELSEIF (INUL == 15) THEN
		NULL ='O'
	    ELSEIF (INUL == 16) THEN
		NULL ='P'
	    ELSEIF (INUL == 17) THEN
		NULL ='Q'
	    ELSEIF (INUL == 18) THEN
		NULL ='R'
	    ELSEIF (INUL == 19) THEN
		NULL ='S'
	    ELSEIF (INUL == 20) THEN
		NULL ='T'
	    ELSEIF (INUL == 21) THEN
		NULL ='U'
	    ELSEIF (INUL == 22) THEN
		NULL ='V'
	   ELSE
	        NULL =' '
	   ENDIF
	   NCATOM = INT(REAL(MOD(INT(CATOM),10000))/100.)	
	   NPOS   = MOD(INT(CATOM),100)
	
	   IF (NCATOM == 1) THEN
              ANAM =' H  '
	   ELSEIF (NCATOM == 6) THEN
		
		IF (NPOS == 0) THEN
		    ANAM =' C  '
		ELSEIF (NPOS == 1) THEN
		    ANAM =' CA '
		ELSEIF (NPOS == 2) THEN
		    ANAM =' CB '
		ELSEIF (NPOS == 3) THEN
		    ANAM =' CE1'
		ELSEIF (NPOS == 4) THEN
		    ANAM =' CE2'
		ELSEIF (NPOS == 5) THEN
		    ANAM =' CD '
		ELSEIF (NPOS == 6) THEN
		    ANAM =' CD1'
		ELSEIF (NPOS == 7) THEN
		    ANAM =' CD2'
		ELSEIF (NPOS == 8) THEN
		    ANAM =' CG '
		ELSEIF (NPOS == 9) THEN
		    ANAM =' CG1'
		ELSEIF (NPOS == 10) THEN
		    ANAM =' CG2'
		ENDIF
	   ELSEIF (NCATOM == 7) THEN
		IF (NPOS == 0) THEN
		    ANAM =' N  '
		ELSEIF (NPOS == 1) THEN
		    ANAM =' N1'
		ELSEIF (NPOS == 2) THEN
		    ANAM =' NH2'
		ELSEIF (NPOS == 3) THEN
		    ANAM =' NE2'
		ELSEIF (NPOS == 4) THEN
		    ANAM =' NZ '
		ENDIF
	   ELSEIF (NCATOM == 8) THEN
		IF (NPOS == 0) THEN
		    ANAM =' O  '
		ELSEIF (NPOS == 1) THEN
		    ANAM =' OH'
		ELSEIF (NPOS == 2) THEN
		    ANAM =' OE1'
		ELSEIF (NPOS == 3) THEN
		    ANAM =' OE2'
		ELSEIF (NPOS == 4) THEN
		    ANAM =' OG1'
		ELSEIF (NPOS == 5) THEN
		    ANAM =' OD1'
		ELSEIF (NPOS == 6) THEN
		    ANAM =' OXT'
		ENDIF
	   ELSEIF (NCATOM == 16) THEN
		IF (NPOS == 0) THEN
		    ANAM =' S  '
		ELSEIF (NPOS == 1) THEN
		    ANAM =' SG '
		ELSEIF (NPOS == 2) THEN
		    ANAM =' SD '
		ENDIF
	   ENDIF
	   NCOUNT = NCOUNT + 1

	   WRITE(98,*) XNUL,INUL
	   WRITE(LUN3,76)PHEAD,NCOUNT,ANAM,RESIDUE,SENQ,NR2,
     &                   NULL,XO,YO,ZO,
     &                   OCCUPANCY, TEMPERATURE,ANAM(2:2)
76	   FORMAT(A5,1X,I5,1X,A4,1X,A3,1X,A1,I4,1A,3X,3F8.3,2F6.2,
     &        '           ',1A)
	   IF (NCATOM == 8 .AND. NPOS == 6) THEN
	      NCOUNT = NCOUNT+1
	      WRITE(LUN3,77)NCOUNT,RESIDUE,SENQ,NR2
	   ENDIF
77	   FORMAT('TER   ',I5,6X,A3,1X,A1,I4)	
	ENDDO	

701     WRITE(LUN3,62)
62	FORMAT('END')

	CLOSE(LUN2)
	CLOSE(LUN3)
	CLOSE(99)

	END

C **********************************************************************
C*
C* Copy PDB file to a SPIDER document file
C*
C********************************************

	SUBROUTINE PDBTOSP

	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'

	PARAMETER  (NLIST=8)
      	REAL	   DLIST(NLIST)
        CHARACTER  (LEN=MAXNAM) :: PDBFILE,TRFILE   
        CHARACTER  *95 RECLIN,HEAD*10, NULL*1,ATOM*4,RESIDUE*3
	CHARACTER  *95 NULL1*1,NULL2*1
	LOGICAL    EX,FLAGCELLDIM,FLAGATOM,ZNUM,COMBIND

	DATA  LUN2,LUN3/25,26/

30	FORMAT(A95)
	LENREC = 0
	CALL OPAUXFILE(.TRUE.,PDBFILE,NULL,LUN2,LENREC,'O',
     &                 'PDB INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

31	FORMAT(' ;  ',A60)
	NATOM = 0

70	READ(LUN2,30) RECLIN
	IF (RECLIN(1:4) == 'ATOM') THEN
	   READ(RECLIN,75)HEAD,N,ATOM,NULL,RESIDUE,NULL1,NR2,
     &          NULL2,XO,YO,ZO,
     &         OCCUPANCY,TEMPERATURE,N
75	   FORMAT(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,1X,I3)	

C I AM NOT SURE ABOUT THE NEXT LINE....
C	     IF (OCCUPANCY.GE.1) THEN

		IF (ATOM(2:2) == 'H') THEN
		   VATOM =100.
CCarbon
		ELSEIF (ATOM(2:2) == 'C') THEN
		   VATOM =600.
			IF (ATOM(3:3) == 'A') THEN
				 VATOM =601
			ELSEIF (ATOM(3:3) == 'B') THEN
				 VATOM =602
			ELSEIF (ATOM(3:4) == 'E1')THEN
			  	VATOM =603
			ELSEIF (ATOM(3:4) == 'E2') THEN
				VATOM =604
			ELSEIF (ATOM(3:3) == 'D') THEN	
				VATOM =605
			ELSEIF (ATOM(3:4) == 'D1') THEN
				VATOM =606
			ELSEIF (ATOM(3:4) == 'D2') THEN
				VATOM =607				
			ELSEIF (ATOM(3:3) == 'G') THEN
				VATOM =608
			ELSEIF (ATOM(3:4) == 'G1') THEN
				VATOM =609	
			ELSEIF (ATOM(3:4) == 'G2') THEN
				VATOM =610
			ELSEIF (ATOM(3:3) == 'Z') THEN
				VATOM =611
			ENDIF
CNitrogen
		ELSEIF (ATOM(2:2) == 'N') THEN
		   VATOM = 700.
			IF (ATOM(3:4) == 'H1') THEN
			  	VATOM =701
			    ELSEIF (ATOM(3:4) == 'H2') THEN
				VATOM =702.
			    ELSEIF (ATOM(3:3) == 'E')THEN
				VATOM =703.
			    ELSEIF (ATOM(3:4) == 'E2') THEN
				VATOM =704.
			    ELSEIF (ATOM(3:3) == 'Z') THEN
				VATOM =705.
			ENDIF
CS	
	        ELSEIF (ATOM(2:2) == 'S') THEN
		   VATOM = 1600.
			IF (ATOM(3:3) == 'G') THEN
				VATOM =1601.
			ELSEIF (ATOM(3:3) == 'D') THEN
				VATOM =1602.
			ENDIF
CP	
		ELSEIF (ATOM(2:2) == 'P') THEN
		      VATOM = 1500.
COxygen
		ELSEIF (ATOM(2:2) == 'O') THEN
		   VATOM = 800.
		   IF (ATOM(3:3) == 'H') THEN
			 VATOM = 801.
		   ELSEIF (ATOM(3:4) == 'E1') THEN
			 VATOM = 802.
		   ELSEIF (ATOM(3:4) == 'E2') THEN
			VATOM = 803.
		   ELSEIF (ATOM(3:4) == 'G1') THEN
			VATOM = 804.
		   ELSEIF (ATOM(3:4) == 'D1') THEN
			VATOM = 805.
		   ELSEIF (ATOM(3:4) == 'XT') THEN
			VATOM = 806.
		   ENDIF
		ENDIF
CResidue
	    IF (RESIDUE(1:3) == 'LEU') THEN
		IRES=1
	    ElSEIF (RESIDUE(1:3) == 'THR') THEN
		IRES=2
	    ELSEIF (RESIDUE(1:3) == 'GLY') THEN
		IRES=3
	    ELSEIF (RESIDUE(1:3) == 'SER') THEN
		IRES=4
	    ELSEIF (RESIDUE(1:3) == 'VAL') THEN
		IRES=5
	    ELSEIF (RESIDUE(1:3) == 'PHE') THEN
		IRES=6
	    ELSEIF (RESIDUE(1:3) == 'ALA') THEN
		IRES=7
	    ELSEIF (RESIDUE(1:3) == 'LYS') THEN
		IRES=8
       	    ELSEIF (RESIDUE(1:3) == 'ARG') THEN
		IRES=9
	    ELSEIF (RESIDUE(1:3) == 'PRO') THEN
		IRES=10
            ELSEIF (RESIDUE(1:3) == 'TAU') THEN
		IRES=11
	    ELSEIF (RESIDUE(1:3) == 'HIS') THEN
		IRES=12
    	    ELSEIF (RESIDUE(1:3) == 'GLU') THEN
		IRES=13
	    ELSEIF (RESIDUE(1:3) == 'GLY') THEN
		IRES=14
	    ELSEIF (RESIDUE(1:3) == 'TYR') THEN
		IRES=15
	    ELSEIF (RESIDUE(1:3) == 'CYS') THEN
		IRES=16
  	    ELSEIF (RESIDUE(1:3) == 'TRY') THEN
		IRES=17
	    ELSEIF (RESIDUE(1:3) == 'ISO') THEN
		IRES=18
	    ELSEIF (RESIDUE(1:3) == 'MET') THEN
		IRES=19
	    ELSEIF (RESIDUE(1:3) == 'ILE') THEN
		IRES=20
	    ENDIF
CCSequence
	    IF (NULL1(1:1) == 'A') THEN
		ISEN=1
	    ELSEIF (NULL1(1:1) == 'B') THEN
		ISEN=2
	    ELSEIF (NULL1(1:1) == 'C') THEN
		ISEN=3
	    ELSEIF (NULL1(1:1) == 'D') THEN
		ISEN=4
	    ELSEIF (NULL1(1:1) == 'E') THEN
		ISEN=5
	    ELSEIF (NULL1(1:1) == 'F') THEN
		ISEN=6
	    ELSEIF (NULL1(1:1) == 'G') THEN
		ISEN=7
	    ELSEIF (NULL1(1:1) == 'H') THEN
		ISEN=8
     	    ELSEIF (NULL1(1:1) == 'I') THEN
		ISEN=9
	    ELSEIF (NULL1(1:1) == 'J') THEN
		ISEN=10
	    ELSEIF (NULL1(1:1) == 'K') THEN
		ISEN=11
	    ELSEIF (NULL1(1:1) == 'L') THEN
		ISEN=12
	    ELSEIF (NULL1(1:1) == 'M') THEN
		ISEN=13
	    ELSEIF (NULL1(1:1) == 'N') THEN
		ISEN=14
	    ELSEIF (NULL1(1:1) == 'O') THEN
		ISEN=15
	    ELSEIF (NULL1(1:1) == 'P') THEN
		ISEN=16
	    ELSEIF (NULL1(1:1) == 'Q') THEN
		ISEN=17
	   ENDIF
	   IF (NULL2 == ' ') THEN
		INUL=0
	   ELSEIF (NULL2 == 'A') THEN
	        INUL=1	 	
	   ELSEIF (NULL2 == 'B') THEN
		INUL=2
	   ELSEIF (NULL2 == 'C') THEN
	        INUL=3	 	
	   ELSEIF (NULL2 == 'D') THEN
	        INUL=4	
	   ELSEIF (NULL2 == 'E') THEN
	        INUL=5	 	
	   ELSEIF (NULL2 == 'F') THEN
		INUL=6
	   ELSEIF (NULL2 == 'G') THEN
	        INUL=7	 	
	   ELSEIF (NULL2 == 'H') THEN
	        INUL=8	
	   ELSEIF (NULL2 == 'I') THEN
	        INUL=9	
	   ELSEIF (NULL2 == 'J') THEN
	        INUL=10	 	
	   ELSEIF (NULL2 == 'K') THEN
		INUL=11
	   ELSEIF (NULL2 == 'L') THEN
	        INUL=12	 	
	   ELSEIF (NULL2 == 'M') THEN
	        INUL=13
	   ELSEIF (NULL2 == 'N') THEN
	        INUL=14	 	
	   ELSEIF (NULL2 == 'O') THEN
		INUL=15
	   ELSEIF (NULL2 == 'P') THEN
	        INUL=16	 	
	   ELSEIF (NULL2 == 'Q') THEN
	        INUL=17	
	   ELSEIF (NULL2 == 'R') THEN
	        INUL=18	
	   ELSEIF (NULL2 == 'S') THEN
	        INUL=19	 	
	   ELSEIF (NULL2 == 'T') THEN
		INUL=20
	   ELSEIF (NULL2 == 'U') THEN
	        INUL=21	 	
	   ELSEIF (NULL2 == 'V') THEN
		INUL=22
	   ENDIF  
	    

           NATOM    = NATOM + 1

	   DLIST(1) = NATOM
	   DLIST(2) = XO
	   DLIST(3) = YO
           DLIST(4) = ZO
           DLIST(5) = IRES * 10000. + VATOM
	   DLIST(6) = TEMPERATURE
           DLIST(7) = NR2*1000. + ISEN*10. + OCCUPANCY
	   DLIST(8) = INUL
	   CALL SAVD(LUN3,DLIST,NLIST,IRTFLG)
	   GOTO 70

        ELSEIF (RECLIN(1:3) == 'END') THEN
	   CLOSE(LUN2)
	   CALL SAVDC
	   CLOSE(LUN3)
	   RETURN

	ELSE
	   GOTO 70
	ENDIF

	END

C **********************************************************************
C*
C* Calculate center of gravity of a PDB file
C*
C********************************************

	SUBROUTINE PDBCG3

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
	REAL, ALLOCATABLE, DIMENSION (:) :: R_TMP
        CHARACTER (LEN=MAXNAM) :: PDBFILE,TRFILE        
        CHARACTER *80 RECLIN,HEAD*10, NULL*1,ATOM*4,RESIDUE*3
	LOGICAL       EX,FLAGCELLDIM,FLAGATOM,ZNUM,COMBIND

	DATA  LUN2,LUN3/25,26/

	LENREC = 0
	CALL OPAUXFILE(.TRUE.,PDBFILE,NULL,LUN2,LENREC,'O',
     &                 'PDB INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
 
	NATOM=0
	UX=0.0
	UY=0.0
	UZ=0.0
	TMASS=0.0

70	READ(LUN2,30) RECLIN
30	FORMAT(A80)

	IF (RECLIN(1:6) == 'HEADER') THEN
           GOTO 70

	ELSEIF (RECLIN(1:4) == 'ATOM') THEN
	   READ(RECLIN,75) HEAD,N,ATOM,NULL,RESIDUE,NULL,NR2,
     &                     NULL,XO,YO,ZO,
     &                     OCCUPANCY,TEMPERATURE,N
75	   FORMAT(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,1X,I3)	
C	   WRITE(*,*) XO,YO,ZO

C            I AM NOT SURE ABOUT THE NEXT LINE....
	     IF (OCCUPANCY .GE. 1) THEN
		IF (ATOM(2:2) == 'H') THEN
		   VATOM = 1.0
		ELSEIF (ATOM(2:2) == 'C') THEN
		   VATOM = 6.0
		ELSEIF (ATOM(2:2) == 'N') THEN
		   VATOM = 7.0
		ELSEIF (ATOM(2:2) == 'O') THEN
		   VATOM = 8.0
	        ELSEIF (ATOM(2:2) == 'S') THEN
		   VATOM = 16.0
		ELSEIF (ATOM(2:2) == 'P') THEN
		   VATOM = 15.0
		ELSEIF (ATOM(2:2) == 'X' .AND. ATOM(1:3) == 'OXT')
     &             THEN
C                  TRNA THREE PRIME TERMINAL
		   VATOM = 8.0
		ELSE
		   WRITE(*, *) 'SPECIAL ATOM ENCOUNTERED IN:'
		   WRITE(*,*) RECLIN 
	           STOP
		ENDIF
C	        VATOM = VATOM*OCCUPANCY
	     ELSE 
	        VATOM = OCCUPANCY
	     ENDIF
	   TMASS = TMASS+VATOM
           UX    = UX+XO*VATOM
           UY    = UY+YO*VATOM
           UZ    = UZ+ZO*VATOM
           NATOM = NATOM+1
	   GOTO 70

        ELSEIF (RECLIN(1:3) == 'END') THEN
	   UX = UX/TMASS
	   UY = UY/TMASS
	   UZ = UZ/TMASS

	   WRITE(*,780)  NATOM
780	   FORMAT(' Number of atoms encountered:',I0)

	   WRITE(*,781)  UX,UY,UZ
781	   FORMAT('  Center of gravity of PDB file',/,3(2x,g12.4))

C          CHANGE COORDINATE SYSTEM
	   SX = UY
	   SY = UX
	   SZ = -UZ

	   WRITE(*,782)  SX,SY,SZ
782	   FORMAT('  Center of gravity of PDB file in SPIDER coords'
     &		  ,/,3(2x,g12.4))

	   ALLOCATE (R_TMP(6))
	   NREG     = 7
	   R_TMP(1) = UX
	   R_TMP(2) = UY
	   R_TMP(3) = UZ
	   R_TMP(4) = SX
	   R_TMP(5) = SY
	   R_TMP(6) = SZ
	   R_TMP(7) = NATOM
	   CALL REG_SET_NSELA(NREG,R_TMP,.FALSE.,IRTFLG)
	   DEALLOCATE(R_TMP)
	   CLOSE(LUN2) 
	   RETURN	
	ELSE
	   GOTO 70
	ENDIF

	CLOSE(LUN2)
	END

C **********************************************************************
C*
C* Copy PDB file to a SPIDER document file
C*
C********************************************

	SUBROUTINE PDBIF

	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'
 
	PARAMETER  (NLIST=4,NRG=6)
      	REAL	               :: DLIST(NLIST)
	REAL	               :: DLIST1(NLIST)
	REAL                   :: R_TMP(NRG)   
        CHARACTER (LEN=MAXNAM) :: PDBFILE,DOCNAM,COMMENT
        CHARACTER (LEN=80)     :: RECLIN
        CHARACTER (LEN=10)     :: HEAD
	CHARACTER (LEN=6 )     :: NULL
	CHARACTER (LEN=4)      :: ATOM
	CHARACTER (LEN=3)      :: RESIDUE
	LOGICAL                :: EX,FLAGCELLDIM,FLAGATOM,ZNUM,COMBIND

        LOGICAL                :: ADDEXT,GETNAME,ISOLD
        LOGICAL                :: APPEND,MESSAGE,NEWFILE
        INTEGER                :: NLET,LUNDOCNO,IRTFLG

        INTEGER, PARAMETER     :: LUN2    = 80
        INTEGER, PARAMETER     :: LUNDOCN = 81

	LENREC = 0
	CALL OPAUXFILE(.TRUE.,PDBFILE,NULL,LUN2,LENREC,'O',
     &                 'PDB INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        ADDEXT  = .TRUE.
        GETNAME = .TRUE.
        ISOLD   = .FALSE.
        APPEND  = .FALSE.
        MESSAGE = .TRUE. 
        IRTFLG  = -8         ! NO IC USE

        CALL OPENDOC(DOCNAM,ADDEXT,NLET,LUNDOCN,LUNDOCNO,GETNAME,
     &      'PARAMETER DOC',ISOLD,APPEND,MESSAGE,NEWFILE,IRTFLG)

C         123456789 123456789 123456789 123456789 123456789 123456789
        COMMENT =
     &   '  Cell constants, Cell angles, Origins, Scales, Atom numbers'
        CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:60),IRTFLG)

80	FORMAT(A80)
90	FORMAT(' ',A,A)

35	FORMAT(A6,3(2X,F7.4),3F7.4)
36	FORMAT(A6,4X,3F10.4,5X,F10.4)

	NCOUNT = 0
	NATOM  = 0
	NTERM  = 0
	NHET   = 0
	NPR    = 0

        WRITE(6,*) ' '
70	READ(LUN2,80) RECLIN

	IF ( RECLIN(1:4) == 'JRNL') THEN
	   WRITE(6,90) RECLIN         
	   GOTO 70

	ELSEIF (RECLIN(13:23) == 'FIT TO DATA') THEN
            WRITE(*,90) RECLIN(12:60)
            DO NPR=1,6
               READ(LUN2,80) RECLIN
               WRITE(6,90) RECLIN(12:60)
            ENDDO
            GOTO 70

	ELSEIF (RECLIN(13:21) == 'DATA USED') THEN
            WRITE(6,90) RECLIN(12:60)
            DO NPR=1,6
               READ(LUN2,80) RECLIN
               WRITE(6,90)RECLIN(12:60)
            ENDDO
            GOTO 70

	ELSEIF (RECLIN(14:20) == 'PROGRAM') THEN
            WRITE(6,90) 'USED PROGRAM',RECLIN(21:60)         
            GOTO 70

	ELSEIF (RECLIN(1:6) == 'CRYST1') THEN
            NLEN = lnblnk(RECLIN)
	    WRITE(6,90) RECLIN(1:NLEN)	

	    READ(RECLIN(1:54),35)NULL,DLIST(2:4),DLIST1(2:4)
	    NCOUNT     = NCOUNT+1
	    DLIST(1)   = NCOUNT
	    R_TMP(1:3) = DLIST(2:4)	
            CALL LUNDOCWRTDAT(LUNDOCNO,NCOUNT,DLIST(2),3,IRTFLG)
	    !CALL SAVD(LUNDOCNO,DLIST,NLIST,IRTFLG)

	    NCOUNT    = NCOUNT+1
	    DLIST1(1) = NCOUNT
	    DLIST(:)  = DLIST1(:)		
            CALL LUNDOCWRTDAT(LUNDOCNO,NCOUNT,DLIST(2),3,IRTFLG)
	    !CALL SAVD(LUNDOCNO,DLIST,NLIST,IRTFLG)
	    GOTO 70

	ELSEIF (RECLIN(1:3) == 'ORI') THEN	
	    NCOUNT   = NCOUNT+1
	    DLIST(1) = NCOUNT	
	    READ(RECLIN,36) NULL,DLIST(2:4),XX
            CALL LUNDOCWRTDAT(LUNDOCNO,NCOUNT,DLIST(2),3,IRTFLG)
	    !CALL SAVD(LUNDOCNO,DLIST,NLIST,IRTFLG)
	    GOTO 70

	ELSEIF (RECLIN(1:3) == 'SCA') THEN
	    NCOUNT   = NCOUNT+1
	    DLIST(1) = NCOUNT	
	    READ(RECLIN,36)NULL,DLIST(2:4),XX
            CALL LUNDOCWRTDAT(LUNDOCNO,NCOUNT,DLIST(2),3,IRTFLG)
	    !CALL SAVD(LUNDOCNO,DLIST,NLIST,IRTFLG)
	    GOTO 70

	ELSEIF (RECLIN(1:4) == 'ATOM') THEN
	    NATOM = NATOM+1
	    GOTO 70

	ELSEIF (RECLIN(1:3) == 'TER') THEN
	    NTERM = NTERM+1
	    GOTO 70

	ELSEIF (RECLIN(1:3) == 'HET') THEN
	    NHET = NHET+1
	    GOTO 70

        ELSEIF (RECLIN(1:3) == 'END') THEN
            CLOSE(LUN2)
            NCOUNT   = NCOUNT+1
            DLIST(1) = NCOUNT
            DLIST(2) = NATOM
            DLIST(3) = NTERM
            DLIST(4) = NHET

            CALL LUNDOCWRTDAT(LUNDOCNO,NCOUNT,DLIST(2),3,IRTFLG)
            !CALL SAVD(LUNDOCNO,DLIST,NLIST,IRTFLG)
            CALL SAVDC
            CLOSE(LUNDOCNO)

            XNATOM   = REAL(NATOM)
            XNTERM   = REAL(NTERM)
            XNHET    = REAL(NHET)
            R_TMP(4) = XNATOM
            R_TMP(5) = XNTERM
            R_TMP(6) = XNHET
            CALL REG_SET_NSELA(NRG,R_TMP,.FALSE.,IRTFLG)

            WRITE(6,*) ' '
            RETURN

	ELSE
            GOTO 70 
	ENDIF

	END


C **********************************************************************

	SUBROUTINE PDBRT3

	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'
        CHARACTER (LEN=MAXNAM) :: PDBFILE,TRFILE
        CHARACTER *80  RECLIN,HEAD*10, NULL*1,ATOM*4,RESIDUE*3
	LOGICAL        EX,FLAGCELLDIM,FLAGATOM,ZNUM,COMBIND
	DOUBLE PRECISION RM(3,3)

	DATA  LUN2,LUN3/25,26/

	NATOM  = 1
	LENREC = 0
	CALL OPAUXFILE(.TRUE.,PDBFILE,NULL,LUN2,LENREC,'O',
     &                 'PDB INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

	CALL OPAUXFILE(.TRUE.,PDBFILE,NULL,LUN3,LENREC,'N',
     &                 'PDB OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

   	CALL RDPRM3S(PHI, THETA, PSI,NOT_USED,
     &              'PHI, THETA, PSI: ',IRTFLG)

	CALL BLDR(RM,PSI,THETA,PHI)

70	READ(LUN2,30) RECLIN
30	FORMAT(A80)

	IF (RECLIN(1:4) .NE. 'ATOM' .AND. 
     &          RECLIN(1:4) .NE. 'END' .AND.
     & 	        RECLIN(1:4) .NE. 'TER' ) THEN
	   WRITE(LUN3,30) RECLIN
	   GOTO 70 

	ELSEIF (RECLIN(1:4) == 'ATOM') THEN
	   READ(RECLIN,75) HEAD,N,ATOM,NULL,RESIDUE,NULL,NR2,
     &                     NULL,XO,YO,ZO,
     &                     OCCUPANCY,TEMPERATURE,N
75	   FORMAT(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,1X,I3)	
C	   WRITE(*,*) XO,YO,ZO

c          CHANGE COORDINATE SYSTEM
	   W  = XO
	   XO = YO
	   YO = W
	   ZO = -ZO

C          AFTER ROTATION CHANGE IT BACK, THAT'S WHY THE ORDER OF X,Y 
C          IS CHANGED AND Z HAS INVERTED SIGN
	   Y=RM(1,1)*XO+RM(1,2)*YO+RM(1,3)*ZO
	   X=RM(2,1)*XO+RM(2,2)*YO+RM(2,3)*ZO
	   Z=-(RM(3,1)*XO+RM(3,2)*YO+RM(3,3)*ZO)
	   WRITE(RECLIN(7:11),80) NATOM
80	   FORMAT(I5)
	   NATOM=NATOM+1
	   WRITE(RECLIN(31:54),85)  X,Y,Z
85	   FORMAT(3F8.3)
	   WRITE(LUN3,30)RECLIN
	   GOTO 70

	ELSEIF (RECLIN(1:3) == 'TER') THEN
	   WRITE(RECLIN(7:11),80)NATOM
	   WRITE(LUN3,30)RECLIN
	   NATOM=NATOM+1
	   GOTO 70

        ELSEIF (RECLIN(1:3) == 'END') THEN
	   WRITE(LUN3,30)RECLIN
	   GOTO 999

	ELSE
	   GOTO 70
	ENDIF

999	CLOSE(LUN2)
	CLOSE(LUN3)
	END
