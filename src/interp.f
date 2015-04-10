
C++********************************************************************
C
C INTERP.FOR
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
C  INTERP(LUNO,LUNN,BUF,NSAMO,NROWO,NSAMN,NROWN)
C
C  PURPOSE:
C       THIS SUBROUTINE INTERPOLATES A GIVEN IMAGE 1 IN A RASTER OF
C       ARBITRARY DIMENSIONS TO GIVE IMAGE 2.
C
C  PARAMETERS:
C         LUNO          LOGICAL UNIT NUMBER OF INPUT IMAGE
C         LUNN          LOGICAL UNIT NUMBER OF OUTPUT IMAGE
C         BUF           BUFFER ARRAY OF SIZE 2*NSAMO+NSAMN
C         NSAMO,NROWO   DIMENSIONS OF INPUT PICTURE
C         NSAMN,NROWN   DIMENSIONS OF OUTPUT PICTURE
C
C       CODING:  BILINEAR INTERPOLATION AS DESCRIBED BY AEBI ET AL.,
C                ULTRASTR. RES.  IS APPLIED TO COMPUTE IMAGE ON
C                NEW RASTER.  PARAMETER LABEL IS COPIED OVER.  OUTPUT
C                PICTURE CAN BE LARGER OR SMALLER THAN INPUT PICTURE.
C
C*RPTX,RPTY SCALES I AND J
C*SUBSCRIPTS OF INTERPOLATED
C*PICTURE DOWN TO SCALE OF
C*ORIGINAL PICTURE,
C*X IS DISTANCE IN X DIR. FROM
C*POINT TO BE INTERPOLATED AT
C*TO OLDIMAGE(NPTY1,NPTX1)
C*Y IS Y DISTANCE
C*     OLDIMAGE(NPTY1,NPTX1)----------OLDIMAGE(NPTY1,NPTX2)
C           -              NEWIMAGE(J,I)             -
C*     OLDIMAGE(NPTY2,NPTX1)----------OLDIMAGE(NPTY2,NPTX2)
C*USED FOR ITS INTERPOLATION.
C*SEE J. OF SUPRAMOLECULAR STRUCTURE
C*PAGE 498 1973
C***************
C--*******************************************************************

      SUBROUTINE INTERP(LUNO,LUNN,NSAMO,NROWO,NSAMN,NROWN,IDUM)

      COMMON BUF(1)

	IF (NSAMO.EQ.2*NSAMN .AND. NROWO.EQ.2*NROWN)  THEN
C          TWO-FOLD DECIMATION BY SUMMATION OF NEIGHBOURING PIXELS.
           K1 = 1
           K2 = K1+NSAMO
           DO J=1,NROWO,2
              DO I=1,NSAMN
                 BUF(K2+I-1)=0.0
              ENDDO
              DO K=1,2
                 CALL  REDLIN(LUNO,BUF,NSAMO,J+K-1)
                 DO I=1,NSAMO,2
                    BUF(K2-1+(I+1)/2)=BUF(K2-1+(I+1)/2)+BUF(I)+BUF(I+1)
                 ENDDO
              ENDDO
              DO I=1,NSAMN
                 BUF(K2+I-1)=BUF(K2+I-1)/4
              ENDDO
              CALL  WRTLIN(LUNN,BUF(K2),NSAMN,(J+1)/2)
           ENDDO
           RETURN

	ELSEIF (NSAMO.EQ.4*NSAMN .AND. NROWO.EQ.4*NROWN)  THEN
C          Four-fold decimation by summation of neighbouring pixels.
           K1=1
           K2=K1+NSAMO
           DO    J=1,NROWO,4
              DO    I=1,NSAMN
                 BUF(K2+I-1)=0.0
              ENDDO
              DO    K=1,4
                 CALL  REDLIN(LUNO,BUF,NSAMO,J+K-1)
                 DO    I=1,NSAMO,4
                    BUF(K2-1+(I+3)/4)=BUF(K2-1+(I+3)/4)+
     &               BUF(I)+BUF(I+1)+BUF(I+2)+BUF(I+3)
                 ENDDO
              ENDDO
              DO    I=1,NSAMN
                 BUF(K2+I-1)=BUF(K2+I-1)/16
              ENDDO
              CALL  WRTLIN(LUNN,BUF(K2),NSAMN,(J+3)/4)
           ENDDO
           RETURN

	ELSEIF(NSAMO.EQ.6*NSAMN .AND. NROWO.EQ.6*NROWN)  THEN
C          SIX-FOLD DECIMATION BY SUMMATION OF NEIGHBOURING PIXELS.
           K1=1
           K2=K1+NSAMO
           DO J=1,NROWO,6
              DO I=1,NSAMN
                 BUF(K2+I-1)=0.0
              ENDDO
              DO K=1,6
                 CALL  REDLIN(LUNO,BUF,NSAMO,J+K-1)
                 DO I=1,NSAMO,6
                    BUF(K2-1+(I+5)/6)=BUF(K2-1+(I+5)/6)+
     &	             BUF(I)+BUF(I+1)+BUF(I+2)+BUF(I+3)+BUF(I+4)+BUF(I+5)
                 ENDDO
              ENDDO
              DO I=1,NSAMN
                 BUF(K2+I-1)=BUF(K2+I-1)/36
              ENDDO
              CALL  WRTLIN(LUNN,BUF(K2),NSAMN,(J+5)/6)
           ENDDO
           RETURN
	ENDIF
C---------------------------------------------------------------
      SIZX=FLOAT(NSAMN)/FLOAT(NSAMO)
      SIZY = FLOAT(NROWN)/FLOAT(NROWO)
      CORX = (NSAMO-1.)/(NSAMN-1.)
      CORY = (NROWO-1.)/(NROWN-1.)

C     WRITE(4,7777)SIZX,CORX,SIZY,CORY
C7777 FORMAT(1X,4F8.3)
      NSAMO2 = NSAMO*2
      NSAMO3 = NSAMO2+NSAMN

      NE  = NROWN - 1
      NES = NSAMN-1

C     DO FIRST LINE
         CALL REDLIN(LUNO,BUF,NSAMO,1)
	 RPTY = 1.
         DO  J = 2,NES
            RPTY = RPTY+CORX
            NPTY1 = INT(RPTY)
            Y = RPTY-NPTY1
            YREM = 1.-Y
            BUF(NSAMO2+J) = YREM*BUF(NPTY1)+Y*BUF(NPTY1+1)
         ENDDO
         BUF(NSAMO2+1) = BUF(1)
         BUF(NSAMO3)   = BUF(NSAMO)
         CALL WRTLIN(LUNN,BUF(NSAMO2+1),NSAMN,1)
         NS1 = 1
         NS2 = NSAMO+1
         RPTX = 1.
         DO  I = 2,NE
            RPTX=RPTX+CORY
            NPTX1=INT(RPTX)
            IF (NPTX1.EQ.NPREV) GOTO 800
            NPTX2=NPTX1+1
	    IF (SIZY .LT. 1) CALL REDLIN(LUNO,BUF(NS1),NSAMO,NPTX1)
            CALL REDLIN(LUNO,BUF(NS2),NSAMO,NPTX2)
            NPREV = NPTX1
C           ALTERNATE BUFFER ADDRESS
            NS1P = NS1
            NS1 = NS2
            NS2 = NS1P

800         X= RPTX-NPTX1
            XREM=1.-X
C           INTERPOLATE  LINE (I) 
	    RPTY = 1.
            DO  J = 2,NES
               RPTY=RPTY+CORX
               NPTY1=INT(RPTY)
               NPTY2=NPTY1+1
               Y=(RPTY-NPTY1)
               YREM=1.-Y
               BUF(NSAMO2+J)=X*(YREM*BUF(NPTY1+NS1-1)+Y*BUF(NPTY1+NS1))
     &             + XREM*(YREM*BUF(NPTY1+NS2-1) + Y*BUF(NPTY1+NS2))
            ENDDO
            BUF(NSAMO2+1) = X*BUF(NS1) + XREM*BUF(NS2)
            BUF(NSAMO3)=X*BUF(NSAMO+NS1-1) +XREM*BUF(NSAMO+NS2-1)
            CALL WRTLIN(LUNN,BUF(NSAMO2+1),NSAMN,I)
         ENDDO

C        DO LAST LINE
	 RPTY = 1.
	 IF(SIZY .LT. 1) CALL REDLIN(LUNO,BUF(NS1),NSAMO,NROWO)
         DO  J=2,NES
            RPTY=RPTY+CORX
            NPTY1=INT(RPTY)
            Y=RPTY-NPTY1
            YREM=1.-Y
            BUF(NSAMO2+J)=YREM*BUF(NPTY1+NS1-1)+Y*BUF(NPTY1+NS1)
         ENDDO
         BUF(NSAMO2+1) = BUF(NS1)
         BUF(NSAMO3)   = BUF(NS1+NSAMO-1)
         CALL WRTLIN(LUNN,BUF(NSAMO2+1),NSAMN,NROWN)
         END
