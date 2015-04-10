
C ++********************************************************************
C                                                                      
C  MRPUTINFO                                                           
C          PROMPTS IMPROVED                       FEB 2014 ARDEAN LEITH
C                                                                      
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C  MRPUTINFO(XYPTS,PRJ,ANGLE,P3D,NUMBER,SHIFT,SCALE,ERVW,ERPT,FULL)    *
C                                                                      *
C  PURPOSE:  RECORDS CALCULATED VALUES TO A FILE                       *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE MRPUTINFO(XYPTS,PRJ,ANGLE,P3D,NUMBER,SHIFT,SCALE,ERVW,
     &            ERPT,FULL,PTACTIVE,NTVW,NTPT,CIR)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      PARAMETER (LV=300)
      PARAMETER (LS=256)
      PARAMETER (NLIST=6)

      LOGICAL*1         FULL
      LOGICAL           PTACTIVE(LS,LV) 
      REAL              CIR(2)

      DIMENSION         XYPTS(2,LS,LV),PRJ(2,LS,LV),SHIFT(2,LV)
      DIMENSION         ANGLE(3,LV),P3D(3,LS)
      DIMENSION         SCALE(LV),ERVW(LV),ERPT(LS)
      INTEGER           NUMBER(LV)
      DIMENSION         DLIST(NLIST)

      CHARACTER(LEN=MAXNAM) :: XYIN1,XYIN2

      DATA  NDOC/55/
      DATA PI/3.141592654/
      DATA DTOR/0.017453292/

      IF (FULL)  THEN
c        WILL PUT ALL OBSERVATIONS INTO .DAT FILE.
         WRITE(NOUT,*)' VIEW POINTS:'

         DO  I=1,NTVW
           DTHETA = ANGLE(2,I)*180.0/PI
           DPSI   = ANGLE(1,I)*180.0/PI

           WRITE(NOUT,720)NUMBER(I),DTHETA,DPSI,SHIFT(1,I),SHIFT(2,I),
     &                    SCALE(I),ERVW(I)
 720       FORMAT('  VIEW:',I3,'  TILT=',F6.1,'  PSI=',F5.1,
     &            '  SHFT(',F5.1,1X,F5.1,')  SCL:',F5.3,'  ER:',E10.4)

           WRITE(NOUT,721)
721        FORMAT
     & ('  PT#        XYPTS                     PROJ            ERROR')

           WRITE(NOUT,*) ' ',('-',J=1,60)

           DO  J=1,NTPT
             WRITE(NOUT,740)J,XYPTS(1,J,I)+CIR(1),XYPTS(2,J,I)+CIR(2),
     &                      PRJ(1,J,I)+CIR(1),PRJ(2,J,I)+CIR(2)
     &  ,SQRT((XYPTS(1,J,I)-PRJ(1,J,I))**2+(XYPTS(2,J,I)-PRJ(2,J,I))**2)
 740         FORMAT(2X,I2,5(2X,1PE10.3))
	   ENDDO
           WRITE(NOUT,*)
         ENDDO
      ENDIF

      WRITE(NOUT,*)'  ANALYSIS PER POINT:'
      WRITE(NOUT,*)'  PT#      X      Y      Z     ERR OF POINT'

      DO  I=1,NTPT
         WRITE(NOUT,710)I,P3D(1,I),P3D(2,I),P3D(3,I),ERPT(I)
 710     FORMAT(2X,I2,5X,F6.1,1X,F6.1,1X,F6.1,1X,F10.6)
      ENDDO
      WRITE(NOUT,*)

      IF (FULL) THEN
C        FILES WITH PREDICTED MARKER POSITIONS ...
         NA   = 80
         NMAX = 0
         CALL  FILSEQP(XYIN1,NA,NUMBER,NMAX,NVW,
     &	    'PREFIX OF X,Y PREDICTED COORD OUTPUT FILES',IRTFLG)

         DO  IVIEWT=1,NTVW
C          GO THROUGH AND PUT ALL X,Y COORDS INTO DOC FILES
           NLETP = NA
           CALL FILGET(XYIN1,XYIN2,NLETP,NUMBER(IVIEWT),IRT2)
           IF (IRT2 .NE. 0) CYCLE

	   IAP  = 0
	   NRUN = 0
	   DO J=1,NTPT
	      IF (PTACTIVE(J,IVIEWT))  THEN
	         DLIST(1) = J
	         DLIST(2) = XYPTS(1,J,IVIEWT)+ CIR(1)
	         DLIST(3) = XYPTS(2,J,IVIEWT)+ CIR(2)
	         DLIST(4) = PRJ(1,J,IVIEWT)  + CIR(1)
	         DLIST(5) = PRJ(2,J,IVIEWT)  + CIR(2)
	         DLIST(6) =
     &             SQRT((XYPTS(1,J,IVIEWT) - PRJ(1,J,IVIEWT))**2 +
     &             (XYPTS(2,J,IVIEWT)-PRJ(2,J,IVIEWT))**2)

	         CALL  SAVDN1(NDOC,XYIN2,DLIST,NLIST,NRUN,IAP)
	         NRUN = 1
	      ENDIF
	   ENDDO

C          STORE ANGLES
           DO J=1,3
              DLIST(J+1) = ANGLE(J,IVIEWT) / DTOR
           ENDDO

           DLIST(1) = -1.0
           NLST     = 4
           CALL SAVDN1(NDOC,XYIN2,DLIST,NLST,NRUN,IAP)
           CLOSE(NDOC)
        ENDDO

      ENDIF

      END
