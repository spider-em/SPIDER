
C++*********************************************************************
C
C    INTEGER FUNCTION NAREA
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
C    NAREA(IX,NP,IXWIN,IXPWIN)
C
C RETURNS THE NUMBER (1..NP) OF THE PATCH FOR A GIVEN COORDINATE
C PAIR SPECIFIED IN IX(1), IX(2). THE CORNER COOS (=TOP LEFT
C COOS IN SPIDER CONVENTION, BOTTOM LEFT IN MRC CONVENTION) OF ALL NP 
C PATCHES ARE DEFINED IN THE COMMON /AREA/. THE DIMENSIONS
C OF THE AVERAGING WINDOW ARE CONTAINED IN IXWIN, AND THOSE 
C OF A PATCH IN IXPWIN.
C 
C CALLED BY: WINAVE2
C
C--*******************************************************************

	INTEGER FUNCTION NAREA(IX,NP,IXWIN,IXPWIN)
	
 

	COMMON/AREA/ IXP(2,1)
	INTEGER   NP
	INTEGER   NXE(2),IX(2),IXWIN(2),IXPWIN(2)

C	NXE(1)=(IXPWIN(1)-IXWIN(1))/2 ! CHANGED 8/13/87
C	NXE(2)=(IXPWIN(2)-IXWIN(2))/2 ! CHANGED 8/13/87
	NXE(1)=(IXPWIN(1)-IXWIN(1))
	NXE(2)=(IXPWIN(2)-IXWIN(2))

	DO 20 IP=1,NP
	DO  I=1,2
C	IF(IX(I).LT.IXP(I,IP)-NXE(I)) GOTO 20 ! CHANGED 8/13/87
C	IF(IX(I).GT.IXP(I,IP)+NXE(I)) GOTO 20 ! CHANGED 8/13/87

	IF(IX(I).LT.IXP(I,IP))GOTO 20
	IF(IX(I).GT.IXP(I,IP)+NXE(I)) GOTO 20
	ENDDO

	NAREA=IP
	RETURN

20	CONTINUE
	NAREA=0

	RETURN
	END
