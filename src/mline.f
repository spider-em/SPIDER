C ++********************************************************************
C
C MLINE
C          LINE INTENSITY                        FEB 2014 ARDEAN LEITH
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
C MLINE(LUN,NX,NY,RP) 
C
C  PURPOSE: PUTS LINES IN SPIDER IMAGE. THIS SHOULD USE   
C           BRESENHAMS ALGORITHM FOR SPEED BUT I AM IN A HURRY SO I      
C           JUST BORROWED SOME EXISTING CODE!!!!!                       
C           CREATE LINE DEFINED BY TWO END POINTS
C
C PARAMETERS:
C
C	LUN      : LOGICAL UNIT NUMBER
C	NX,NY    : FILE DIMENSIONS
C	RP	 : INTENSITY VALUE FOR LINE
C
C CALLED BY:  pttern.f          
C
C **********************************************************************

C       -------------  MLINE2  ----------------------------------

	SUBROUTINE MLINE2(LUN,NX,NY,RP,IRTFLG)

        IMPLICIT NONE
	INCLUDE 'CMLIMIT.INC' 
	INCLUDE 'CMBLOCK.INC' 

        INTEGER           :: LUN,NX,NY,IRTFLG
        REAL              :: RP

        REAL              :: BUF(NX,NY)
        REAL              :: X1,Y1,X2,Y2
        INTEGER           :: IX1,IY1,IX2,IY2,NOT_USED
        REAL              :: FOREGR
 
        
C       LOAD IMAGE
        CALL REDVOL(LUN,NX,NY,1,1,BUF,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        FOREGR = RP

	CALL RDPRM3S(X1,Y1,FOREGR,NOT_USED,
     &               'COORDINATES & INTENSITY OF FIRST  POINT',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IX1 = IFIX(X1)
        IY1 = IFIX(Y1)
        IF (IX1 <= 0 .OR. IY1 <= 0) RETURN
           
        CALL RDPRM2S(X2,Y2,NOT_USED,
     &                  'COORDINATES OF SECOND POINT',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IX2 = IFIX(X2)
        IY2 = IFIX(Y2)
        IF (IX2 <= 0 .OR. IY2 <= 0) RETURN
             
        CALL MLINEB(BUF,NX,NY,IX1,IY1,IX2,IY2,FOREGR,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
	   
C       SAVE IMAGE
        CALL WRTVOL(LUN,NX,NY,1,1,BUF,IRTFLG)

        IRTFLG = 0
	END



C       -------------------  MLINEB  ---------------------------


        SUBROUTINE MLINEB(BUF,NX,NY,IX1,IY1,IX2,IY2,FOREGR,IRTFLG)

        IMPLICIT NONE
	INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC'

        REAL              :: BUF(NX,NY)
        INTEGER           :: NX,NY,IX1,IY1,IX2,IY2,IRTFLG
        REAL              :: FOREGR

        INTEGER           :: IGO,IEND,IX,IXL,IY,IGOX,IENDX,IHALF
        INTEGER           :: IYT,IYL,IYEND,IXT
        REAL              :: FACT,FCON

        IRTFLG = 1
        IF     ((IX1 > NX .OR. IX1 <= 0) .OR.
     &          (IY1 > NY .OR. IY1 <= 0)) THEN

           WRITE(NOUT,91) IX1,IY1
91         FORMAT('*** POINT: (',I0,',',I0,') OUTSIDE IMAGE')
           RETURN

        ELSEIF ((IX2 > NX .OR. IX2 <= 0) .OR.
     &          (IY2 > NY .OR. IY2 <= 0)) THEN

           WRITE(NOUT,91) IX2,IY2
           RETURN
        ENDIF

C       POINT IS WITHIN IMAGE
        IF (IY1 == IY2) THEN
C          HORIZONTAL LINE WOULD CAUSE DIVISION BY ZERO
           IGO  = MIN(IX1,IX2)
           IEND = MAX(IX1,IX2)

           BUF(IGO:IEND,IY1) = FOREGR
           !WRITE(6,*) 'horizontal',iy1,iy2

        ELSE
           FACT =  FLOAT(IX2-IX1) /  FLOAT(IY2-IY1)
           FCON =  - FACT * IY1 + 0.5
           IF (IX1 > IX2) FCON = - FACT * IY1 - 0.5
           IGO  = MIN(IY1,IY2)
           IEND = MAX(IY1,IY2)

           !WRITE(6,*) 'go,end:',igo,iend

           IXL = IX1
           IF (IY2 < IY1) IXL = IX2

           DO IY = IGO,IEND
C             FIND X VALUE FOR THIS Y COORDINATE
              IX = IX1 + IFIX(FACT * FLOAT(IY) + FCON)

C             SET BUFFER AT THIS LOCATION TO FOREGROUND
              BUF(IX,IY) = FOREGR

              !write(6,'(A,I0," , ",I0)'),' 1: ',ix,iy

             IF (IABS(IX - IXL) > 1) THEN
C               MUST ADD IN INTERPOLATED X VALUES 
                IGOX  = MIN(IX,IXL)
                IENDX = MAX(IX,IXL)
                IHALF = IGOX + (IENDX - IGOX) / 2

                IYT   = IYL
                IYEND = IY
                IF (IX < IXL) THEN
                    IYT   = IY
                    IYEND = IYL
                ENDIF

                DO IXT = IGOX,IENDX
                   BUF(IXT,IYT) = FOREGR
                   !write(6,'(A,I0," , ",I0)')' 2: ',ixt,iyt
                   IF (IXT == IHALF) IYT = IYEND                       
                ENDDO
             ENDIF
             IXL = IX
             IYL = IY         

           ENDDO 
        ENDIF

C********************************
        !WRITE(6,977) IX1,IY1,ix2,iy2
977      FORMAT(' (',I4,',',I4,')-->(',I4,',',I4,')')
C*********************************

        IRTFLG = 0

	END

C       ---------  MLINE  --------- (ABANDONED FEB 2014) ---------------


	SUBROUTINE MLINE(LUN,NX,NY,RP)

	INCLUDE 'CMLIMIT.INC' 
	INCLUDE 'CMBLOCK.INC' 

	COMMON /IOBUF/ BUF(NBUFSIZ)

	CALL RDPRI2S(IX1,IY1,NOT_USED,'COORDINATES OF FIRST  POINT',IRT)
        IF (IRT .NE. 0) RETURN

	CALL RDPRI2S(IX2,IY2,NOT_USED,'COORDINATES OF SECOND POINT',IRT)
        IF (IRT .NE. 0) RETURN

	IF (IY1  >  IY2) THEN
	  IG  = IX1
	  IX1 = IX2
	  IX2 = IG
	  IG  = IY1
	  IY1 = IY2
	  IY2 = IG
        ENDIF

	ISTART = MAX(1,IY1)
	IEND   = MIN(NY,IY2)

	IF (ISTART > NY .OR. IEND <= 0) THEN
	   CALL ERRT(101,'LINE GOES OUTSIDE IMAGE',NF)
           RETURN
        ENDIF
             
        !write(6,*) 'start..end:',istart,iend

	IF (IY1 .NE. IY2) THEN
	   DO IROW=ISTART,IEND
	      IF (IX1 > IX2) THEN
                IX0 = IX1 + IFIX((FLOAT(IX2-IX1) / 
     &                            FLOAT(IY2-IY1)) * FLOAT(IROW-IY1)-0.5)
              ELSE
	        IX0 = IX1 + IFIX((FLOAT(IX2-IX1) /
     &                            FLOAT(IY2-IY1)) * FLOAT(IROW-IY1)+0.5)
              ENDIF
            !write(6,*) 'i,ixo:',IROW,IX0

              IF (IX0 > 0 .AND. IX0 <= NX) THEN
                 CALL REDLIN(LUN,BUF,NX,IROW)

                 BUF(IX0) = RP

                 CALL WRTLIN(LUN,BUF,NX,IROW)
              ENDIF
           ENDDO

	   RETURN
        ENDIF


	IF (IX1 > IX2) THEN
	   IG  = IX1
	   IX1 = IX2
	   IX2 = IG
        ENDIF

	ISTART = MAX(1,IX1)
	IEND   = MIN(NX,IX2)

	IF (ISTART > NX .OR. IEND <= 0) THEN
	   CALL ERRT(101,'LINE GOES OUTSIDE IMAGE',NF)
           RETURN
        ENDIF

	CALL REDLIN(LUN,BUF,NX,IY1)

	DO  I=ISTART,IEND
	   BUF(I) = RP
	ENDDO

	CALL WRTLIN(LUN,BUF,NX,IY1)

	END



