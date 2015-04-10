C++*******************************************************************
C
C ROT32.F      USED ALLOCATE                       AUG 00 ARDEAN LEITH
C              VOLUME RECORD BUG                   DEC 13 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2013  Health Research Inc.,                         *
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
C    PURPOSE: ROTATE AN IMAGE BY AN ARBITRARY ANGLE OF
C     DEGREE BETWEEN 0 AND 360.  CALLED WHEN THE IMAGE IS SMALL
C     ENOUGH TO FIT IN THE BUFFER.
C
C    ROT32(LUNI,LUNO,NX,NROWS,NROWE,NROWSK,BUF,THETA,BACK,SHX,SHY)
C
C    PARAMETERS:
C
C       LUNI         LOGICAL UNIT NUMBER OF INPUT IMAGE
C       LUNO         LOGICAL UNIT NUMBER OF OUTPUT IMAGE
C	NX	     ROW LENGTH
C	NROWS,NROWE  STARTING AND ENDING ROW
C       NROWSK       SKIPPING FACTOR FOR ROWS
C                    (VALUE <0 MEANS NON-SEQUENTIAL INPUT AND OUTPUT) 
C       THETA        ROTATION ANGLE IN RADIANS
C       BACK         AVERAGE OF INPUT IMAGE
C       SHX,SHY      ORIGIN SHIFT
C       IOFF         OUTPUT OFFSET
C
C       POSITIVE THETA: COUNTER-CLOCKWISE ROTATION
C                       (MATHEMATICALLY POSITIVE DIRECTION)
C
C       AN IMAGE CAN BE A SLICE OF A THREE-DIMENSIONAL DENSITY
C       DISTRIBUTION. FOR THIS REASON, A STARTING ROW UNEQUAL
C       TO 1 AND AN ENDING ROW UNEQUAL TO NROW IS PERMITTED.
C
C--*******************************************************************

	SUBROUTINE ROT32(LUNI,LUNO,NX,NROWS,NROWE,NROWSK,
     &                   THETA,BACK,SHX,SHY,IOFF)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

	INTEGER  :: LUNI,LUNO,NX,NROWS,NROWE,NROWSK,IOFF
        REAL     :: THETA,BACK,SHX,SHY

        COMMON /IOBUF/ RBUF(NBUFSIZ)
 
C       DIMENSIONS OF BUF ARE NX + NX * NROW,  WHERE
C       NROW = ((NROWE-NROWS) / NROWSK) + 1

        REAL, ALLOCATABLE  :: BUF(:)
  
	REAL, PARAMETER    ::  PI = 3.14159

        IFLAG1 = 0
        IF (NROWSK < 0) IFLAG1 = 1
        IF (NROWSK < 0) NROWSK = -NROWSK

        IF ( THETA == 0.0 ) THEN
C          ZERO DEGREE ROTATION
           II = 0

	   DO  I=NROWS,NROWE,NROWSK
              IF (IFLAG1 == 0) II = II + 1
              IF (IFLAG1 == 1) II = I

              CALL REDLIN(LUNI,RBUF,NX,I)
	      CALL WRTLIN(LUNO,RBUF,NX,II)
	   ENDDO
	   RETURN
        ENDIF

 	NXH   = NX/2
	NROW  = ((NROWE-NROWS) / NROWSK) + 1
	NROWH = NROW/2
	KCENT = NXH+1
	ICENT = NROWH+1

        ALLOCATE (BUF(NX*NROW), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
            CALL ERRT(46,'ROT32, BUF',IER)
            RETURN
        ENDIF     

  	IF (THETA >  PI) THETA = -2.0 * PI + THETA
	IF (THETA < -PI) THETA =  2.0 * PI + THETA
	COD = COS(THETA)
	SID = SIN(THETA)

C       READ IN WHOLE INPUT IMAGE
        J = 0
	DO  I = NROWS,NROWE,NROWSK
           J = J + 1
           L = (J-1)*(NX)+1
           CALL REDLIN(LUNI,BUF(L),NX,I)        
	ENDDO

        !write(6,*) 'nrows,nrowe,nrowsk:',nrows,nrowe,nrowsk,irecoff

C       NOW GO THROUGH OUTPUT COO SYSTEM; COMPUTE, FOR EACH POINT IN
C       ROW, THE POSITION IN THE OLD COO SYSTEM.
C       THEN CALCULATE POINT FROM FOUR SURROUNDING POINTS USING BILINEAR
C       INTERPOLATION. WRITE OUT EACH LINE AS YOU GO ALONG.

C       JUST TO ALLOW A CHANGE IN THE ROTATIONAL CENTER
C       (RELATIVE TO THE CENTRAL PIXEL)

        RICENT = ICENT + SHY
        RKCENT = KCENT + SHX  

        JJ     = 0
	DO  I = 1,NROW
           JJ = JJ+1
           IF (IFLAG1 == 1) II = NROWS + (I-1) * NROWSK
           IF (IFLAG1 == 0) II = I

           Y    = I - RICENT
           YCOD =  Y * COD + RICENT
           YSID = -Y * SID + RKCENT

           DO K = 1,NX
              RBUF(K) = BACK
              X       = K - RKCENT
              XOLD    = X * COD + YSID
              YOLD    = X * SID + YCOD
              IYOLD   = YOLD
              YDIF    = YOLD - IYOLD
              YREM    = 1.   - YDIF
              IXOLD   = XOLD

              IF ((IYOLD >= 1 .AND. IYOLD <= NROW-1) .AND.
     &            (IXOLD >= 1 .AND. IXOLD <= NX-1)) THEN

c                INSIDE BOUNDARIES OF OUTPUT IMAGE
                 XDIF    = XOLD - IXOLD
                 XREM    = 1.   - XDIF
                 NADDR   = (IYOLD-1) * NX + IXOLD 
                 RBUF(K) = YDIF * (BUF(NADDR+NX) * XREM +
     &                     BUF(NADDR+NX+1)*XDIF) +
     &                     YREM * (BUF(NADDR)*XREM + BUF(NADDR+1)*XDIF)

              ENDIF
           ENDDO

           IRECT = II + IOFF
           CALL WRTLIN(LUNO,RBUF,NX,IRECT)
	ENDDO

        IF (ALLOCATED(BUF)) DEALLOCATE(BUF)

	END
