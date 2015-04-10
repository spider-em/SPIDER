
C++*********************************************************************
C
C SPEAKM.F    NEW                               MAY 2001 ARDEAN LEITH
C             REWRITE                           NOV 2012 ARDEAN LEITH
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
C  PURPOSE:  SEARCHES FOR THE  LOCATION OF MAXIMUM IN THE (REAL) IMAGE
C            AND PRINTS OUT POSITION AND VALUE OF THIS PEAK.
C
C  SPEAKM(LUN,NX,NY,NZ,FMAXT)
C
C  PARAMETERS: LUN      I/O UNIT NUMBER
C              NX,NY    DIMENSIONS OF IMAGE
C              NZ       DIMENSIONS OF IMAGE
C              FMAXT    MAX. VALUE OF IMAGE
C
C        REGISTER POSITIONS 
C          1 = X-PEAK LOCATION
C          2 = Y-PEAK LOCATION
C          3 = X-PEAK LOCATION RELATIVE TO ORIGIN
C          4 = X-PEAK LOCATION RELATIVE TO ORIGIN
C          5 = VALUE OF MAXIMUM
C          6 = NUMBER OF MAXIMUM LOCATIONS
C          7 = RADIUS OF BOUNDING BOX OF MAXIMA
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE SPEAKM(LUN,NX,NY,NZ,FMAXT)

         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC'

         COMMON /IOBUF/ BUF(NBUFSIZ)

         INTEGER :: LUN,NX,NY,NZ
         REAL    :: FMAXT
         REAL    :: RSQ(9)

         CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

         IF (NZ > 1) THEN
            CALL ERRT(101,'NOT IMPLEMENTED FOR VOLUMES',NE)
            RETURN           
         ENDIF

C        ORIGIN COORDINATES
         NXCTR = NX/2+1
         NYCTR = NY/2+1
         NZCTR = NZ/2+1
   
         CALL RDPRI3S(NXCTR,NYCTR,NZCTR,NOT_USED,
     &        'ENTER ORIGIN COORDINATES OR <CR> FOR CENTER',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         FX    = 0.0
         FY    = 0.0
         FZ    = 0.0
         NMAX  = 0

C        INITIALIZE BOUNDING BOX CORNERS
         IXMIN = NX
         IXMAX = 1
         IYMIN = NY
         IYMAX = 1
         IZMIN = NZ
         IZMAX = 1

         DO IZ = 1,NZ
            DO IY = 1,NY

               CALL REDLIN(LUN,BUF,NX,(NZ-1)*NY +IY)

               DO IX = 1,NX
                  IF (BUF(IX) >= FMAXT) THEN
                     NMAX  = NMAX + 1
                     FX    = FX + IX
                     FY    = FY + IY
                     FZ    = FZ + IZ

C                    FIND BOUNDING BOX FOR DISTANCES
                     IXMIN = MIN(IXMIN,IX)
                     IYMIN = MIN(IYMIN,IY)
                     IZMIN = MIN(IZMIN,IZ)
                     IXMAX = MIN(IXMAX,IX)
                     IYMAX = MIN(IYMAX,IY)
                     IZMAX = MIN(IZMAX,IZ)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         IF (VERBOSE .AND. MYPID <= 0) WRITE(NOUT,*) ' '
         IF (NMAX == 0)  THEN
            WRITE(NOUT,*) ' NO PEAK FOUND'
            IF (NDAT .NE. NOUT) WRITE(NDAT,*) ' NO PEAK FOUND'
            CALL REG_SET_NSEL(1, 5,0.0, 0.0, 0.0, 0.0, 0.0,IRTFLG)
            CALL REG_SET_NSEL(6, 2,0.0, 0.0, 0.0, 0.0, 0.0,IRTFLG)
            RETURN
         ENDIF

         FX  = FX / NMAX
         FY  = FY / NMAX
         FZ  = FZ / NMAX

         IF (NZ == 1) THEN
            IF (NMAX == 1) THEN
C              SINGLE MAXIMUM PEAK

C              9/25/81 PARABOLIC FIT TO 3X3 NEIGHBORHOOD OF PEAK POINT
C              PROGRAM SENT BY M.VAN HEEL, MODIFIED FOR SPIDER. JF
            
               KL = FX
               DO I=1,3
                  IY = FY + I - 2
                  IF (IY < 1)  IY = IY + NY
                  IF (IY > NY) IY = IY - NY
                  CALL REDLIN(LUN,BUF,NX,IY)

                  I1 = (I-1) * 3

                  DO K=1,3
                     IX = KL+K-2
                     IF (IX < 1)  IX = IX + NX
                     IF (IX > NX) IX = IX - NX

                     RSQ(I1+K) = BUF(IX)
                  ENDDO
               ENDDO

C              PARABL RETURNS PEAKV

               CALL PARABL(RSQ,XSH,YSH,PEAKV)
               FX  = XSH + FX
               FY  = YSH + FY
               RAD = 0.0

C              HACK FOR ALMOST BINARY PEAK ERROR OR NON-ADJACENT MAXS
               IF (PEAKV < MAXVAL(RSQ))  PEAKV = MAXVAL(RSQ)

            ELSE
C              MULTIPLE MAXIMUM PEAKS

C              FIND MAX DISTANCE TO CORNER OF BOUNDING BOX
               XMAX  = MAX((IXMAX-FX),(FX-IXMIN))
               YMAX  = MAX((IYMAX-FY),(FY-IYMIN))
               RAD   = SQRT((XMAX**2 + YMAX**2))

               PEAKV = FMAXT
            ENDIF

C           ADJUST FOR ORIGIN
            FXC = FX - NXCTR
            FYC = FY - NYCTR

            WRITE(NOUT,90)
            IF (NDAT .NE. NOUT) WRITE(NDAT,90)

90          FORMAT('    NMAX      X        Y         XO       YO',
     &             '       HEIGHT      RADIUS')

            WRITE(NOUT,91) NMAX, FX,FY, FXC,FYC, PEAKV,RAD
91          FORMAT(I7,2(F8.2,1X),2X,2(F8.2,1X),G16.5,3X,F5.2)

            IF (NDAT .NE. NOUT) WRITE(NDAT,91)
     &                    NMAX, FX,FY, FXC,FYC, PEAKV,RAD

            IF (VERBOSE .AND. MYPID <= 0) WRITE(NOUT,*) ' '


            CALL REG_SET_NSEL(1,5,FX,FY,FXC,FYC,PEAKV,IRTFLG)
            CALL REG_SET_NSEL(6,2,FLOAT(NMAX),RAD,0.0,0.0,0.0,IRTFLG)
 
         ELSE
C           VOLUME, NOT IMPLEMENTED YET
         ENDIF

         RETURN
         END

