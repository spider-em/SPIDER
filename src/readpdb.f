
C **********************************************************************
C  READPDB.F       ADDED TEMPERATURE               MAR 02 ARDEAN LEITH
C                  HG11 ATOM                       APR 12 ARDEAN LEITH
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
C  READPDB
C
C  PURPOSE:  READ PDB FILE, GET COORDINATES,
C            GENERATE A SPIDER FORMAT VOLUME
C 
C **********************************************************************

	SUBROUTINE READPDB
	
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

	REAL, ALLOCATABLE     :: BUF(:,:,:)
        CHARACTER(LEN=MAXNAM) :: PDBFILE,SPIDERFILE,RECLIN
        CHARACTER(LEN=10)     :: HEAD 
        CHARACTER(LEN=4)      :: ATOM 
        CHARACTER(LEN=1)      :: NULL,CENTER,CTEMP
	LOGICAL               :: ZNUM,TEMP

	ZNUM = .FALSE.

	NATOM = 0
	TMASS = 0.0
   	NULL  = CHAR(0)
        LUN1  = 8
	LUN2  = 10
	XMAX  = -1.0E23
	XMIN  = 1.0E23
	YMAX  = -1.0E23
	YMIN  = 1.0E23
	ZMAX  = -1.0E23
	ZMIN  = 1.0E23
	UX    = 0.0
	UY    = 0.0
	UZ    = 0.0
	
C       OPEN PDB FILE (WHICH IS READABLE ASCII) AS FORMATTED SEQ.       
        LENREC = 0
        CALL OPAUXFILE(.TRUE.,PDBFILE,NULL,LUN1,LENREC,'O',
     &                 'PDB INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

15	READ(LUN1,20) RECLIN
20	FORMAT(A80)

C       GET MIN AND MAX OF X,Y,Z  AND CALCULATE OFFSET

35	IF (RECLIN(1:4) .EQ. 'ATOM') THEN

	   READ (RECLIN,50) HEAD,ATOM,X,Y,Z,OCCUPANCY,TEMPERATURE
	   IF (OCCUPANCY > 1.0 ) ZNUM = .TRUE.

C          COORDINATE SYSTEM IS DIFFERENT IN O AND SPIDER 
           W = X
           X = Y
           Y = W
           Z = -Z

	   XMAX  = MAX(XMAX,X)
	   XMIN  = MIN(XMIN,X)
	   YMAX  = MAX(YMAX,Y)
	   YMIN  = MIN(YMIN,Y)
	   ZMAX  = MAX(ZMAX,Z)
	   ZMIN  = MIN(ZMIN,Z)

           UX    = UX + X
           UY    = UY + Y
           UZ    = UZ + Z
           NATOM = NATOM+1

	   GOTO 15

	ELSEIF (RECLIN(1:3) .NE. 'END') THEN
	   GOTO 15
	ENDIF

        IF (NATOM <= 0)  THEN
          CALL ERRT(101,'NO ATOMS IN THE PDB FILE',NE)
          GOTO 9999
        ENDIF

C       SET CELL DIMENSION

	SX = (XMAX-XMIN) + 3.0
	SY = (YMAX-YMIN) + 3.0
	SZ = (ZMAX-ZMIN) + 3.0

	WRITE(NOUT,91) SX,SY,SZ
91      FORMAT('  Cell size: ',F7.1,',',F7.1,',',F7.1)

	REWIND (LUN1)

C       INPUT SAMPLE SIZE
        PIXEL = 1.0
	CALL RDPRM1S(PIXEL,NOT_USED,'VOXEL SIZE [A]',IRTFLG)

        CALL RDPRMC(CENTER,NCHAR,.TRUE.,'CENTER? (Y/N)',
     &               NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

	IF (CENTER .NE. 'Y')  THEN
          UX = 0.0
          UY = 0.0
          UZ = 0.0
        ELSE
          UX = UX / NATOM
          UY = UY / NATOM
          UZ = UZ / NATOM
        ENDIF

        CALL RDPRMC(CTEMP,NCHAR,.TRUE.,
     &      'ATOMS OR TEMPERATURE? (A/T)',NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

	TEMP = (CTEMP == 'T')
 
        XNEED = (UX-XMIN-PIXEL) * 2
        YNEED = (UY-YMIN-PIXEL) * 2
        ZNEED = (UZ-ZMIN-PIXEL) * 2

	NX    = NINT(XNEED) + 1
	NY    = NINT(YNEED) + 1
	NZ    = NINT(ZNEED) + 1

	WRITE(NOUT,506) NX,NY,NZ
506	FORMAT('  Minimum size needed for this volume: ',3I6)

        CALL RDPRI3S(NX,NY,NZ,NOT_USED,
     &              'SPIDER VOLUME SIZE; NX, NY & NZ',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

	ALLOCATE(BUF(NX,NY,NZ),STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           CALL  ERRT(46,'READPDB; BUF',NX*NY*NZ)
           GOTO 9999
        ENDIF
	
	I       = NX/2+1
	XOFFSET = PIXEL*I
	I       = NY/2+1
	YOFFSET = PIXEL*I
	I       = NZ/2+1
	ZOFFSET = PIXEL*I

#ifdef NEVER
        YMIN-UY < (-1.  * YOFFSET)
        YMIN-UY < (-1.  * PIXEL*I)
        YMIN-UY < (-1.  * PIXEL*(NY/2+1))
        YMIN-UY <  -1.  * (PIXEL*NY/2 + PIXEL) 
       -YMIN+UY <          PIXEL*NY/2 + PIXEL
       (-YMIN+UY-PIXEL) * 2 / PIXEL <  NY 
       NY >= (-YMIN+UY-PIXEL) * 2
       NY >= (UY-YMIN-PIXEL) * 2
        TMPY = (UY-YMIN-PIXEL) * 2
        WRITE(6,*)' NY, TMPY:',XMIN-UX,' >= ',TMPY
#endif
C       CHECK IF THE VOLUME IS LARGE ENOUGH TO CONTAIN THE OBJECT
        IF (XMIN-UX < (-1. * XOFFSET) .OR. 
     &      YMIN-UY < (-1. * YOFFSET) .OR.
     &      ZMIN-UZ < (-1. * ZOFFSET)) THEN

            WRITE(6,*)' XMIN-UX < -XOFFSET:',XMIN-UX,' < ',-XOFFSET
            WRITE(6,*)' YMIN-UY < -YOFFSET:',YMIN-UY,' < ',-YOFFSET
            WRITE(6,*)' ZMIN-UZ < -ZOFFSET:',ZMIN-UZ,' < ',-ZOFFSET

            CALL ERRT(101,'VOLUME TOO SMALL TO CONTAIN OBJECT',NE)
            GOTO 9999
        ENDIF
 
        TMPX = FLOAT(NX-1) * PIXEL - XOFFSET
        TMPY = FLOAT(NY-1) * PIXEL - YOFFSET
        TMPZ = FLOAT(NZ-1) * PIXEL - ZOFFSET

        IF (XMAX-UX > TMPX .OR. 
     &      YMAX-UY > TMPY .OR.
     &      ZMAX-UZ > TMPZ) THEN

            WRITE(6,*)' XMAX-UX > TMPX:',XMAX-UX,' > ',TMPX
            WRITE(6,*)' YMAX-UY > TMPY:',YMAX-UY,' > ',TMPY
            WRITE(6,*)' ZMAX-UZ > TMPZ:',ZMAX-UZ,' > ',TMPZ

            CALL ERRT(101,'VOLUME TOO SMALL TO CONTAIN OBJECT.',NE)
            GOTO 9999
        ENDIF

C       ZERO WHOLE BUF
        BUF = 0.0

        IFORM = 3
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,SPIDERFILE,LUN2,'U',IFORM,
     &              NX,NY,NZ, MAXIM,'SPIDER OUTPUT', .FALSE.,IRTFLG)
	IF (IRTFLG == 1) GOTO 9999
	
C       READ IN ATOMS
        NATOM = 0
        TMAX  = -100

40	CONTINUE
	READ(LUN1,20) RECLIN

	IF (RECLIN(1:4) .EQ. 'ATOM' ) THEN

!HEAD--NNNNxATOMxAAAx?????????? X      Y       Z        OCC   TEMP

 	   READ(RECLIN,50)HEAD,ATOM,X,Y,Z,OCCUPANCY,TEMPERATURE
50	   FORMAT(A6,5X,1X,A3,2X,3X,1X,1X,4X,1X,3X,3F8.3,2F6.2)
	
C            WRITE(NOUT,*)'ATOM #',NATOM,'  ',ATOM,' AT ',X,Y,Z,
C    &                    OCCUPANCY,TEMPERATURE
C            WRITE(6, *) ' N ATOM:',ATOM

	     IF (TEMP) THEN
                VATOM = TEMPERATURE
                TMAX  = MAX(TMAX,TEMPERATURE)

	     ELSEIF (.NOT. ZNUM) THEN
		IF (ATOM(2:2) .EQ. 'H' .OR.
     &              ATOM(1:1) .EQ. 'H'  ) THEN
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
     &             THEN   !TRNA THREE PRIME TERMINAL
		   VATOM = 8.0

		ELSE
                   WRITE(NOUT, *) ' SPECIAL ATOM:',ATOM,
     &                ':   ENCOUNTERED IN:'
                   WRITE(NOUT,*) RECLIN 
                   CALL ERRT(101,'UNSUPPORTED SPECIAL ATOM',NDUM)
	           GOTO 9999
		ENDIF

C	        VATOM = VATOM*OCCUPANCY
	     ELSE 
	        VATOM = OCCUPANCY
	     ENDIF

C            COORDINATE SYSTEM IS DIFFERENT IN O AND SPIDER 
	     W = X
	     X = Y
	     Y = W
	     Z = -Z

             X = X-UX
             Y = Y-UY
             Z = Z-UZ

C            TRILINEAR INTERPOLATION
	     IX = IFIX((X+XOFFSET)/PIXEL)
	     AX = (X+XOFFSET)/PIXEL-FLOAT(IX)
	     IY = IFIX((Y+YOFFSET)/PIXEL)
	     AY = (Y+YOFFSET)/PIXEL-FLOAT(IY)
	     IZ = IFIX((Z+ZOFFSET)/PIXEL)
	     AZ = (Z+ZOFFSET)/PIXEL-FLOAT(IZ)
             IF (TEMP) THEN
                BUF(IX,IY,IZ) = VATOM
             ELSE

                DO I3=1,0,-1
                   DO I2=1,0,-1
                     DO I1=1,0,-1
                        BUF(IX-I1+1,IY-I2+1,IZ-I3+1) =
     &                    BUF(IX-I1+1,IY-I2+1,IZ-I3+1) +
     &                    VATOM*(REAL(I1)-(2*I1-1)*AX)*
     &                    (REAL(I2)-(2*I2-1)*AY)*
     &                    (REAL(I3)-(2*I3-1)*AZ)    
                     ENDDO 
                   ENDDO 
                 ENDDO 
             ENDIF

             NATOM = NATOM + 1
	     TMASS = TMASS + VATOM
	     GOTO 40

	 ELSEIF (RECLIN(1:3) .NE. 'END') THEN
	     GOTO 40
	 ENDIF

         CALL WRTVOL(LUN2,NX,NY, 1,NZ, BUF,IRTFLG)

         IF (TEMP) WRITE(NOUT,*) ' Max. temperature:',TMAX

         WRITE(NOUT,90) NATOM,TMASS
90       FORMAT('  Number of atoms placed:',I8,
     &           '   Total mass:',ES10.3)

9999     CLOSE(LUN2)
         CLOSE(LUN1)
    	 IF (ALLOCATED(BUF)) DEALLOCATE(BUF) 

	 END
