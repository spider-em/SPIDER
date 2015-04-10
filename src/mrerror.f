
C ++********************************************************************
C                                                                     
C  MRERROR
C          PROMPTS & DOC FILE HEADERS IMPROVED   FEB 2014 ARDEAN LEITH                                                                *
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
C                                                                      
C MRERROR(P3DREF,P3D1,ERRORP,NTPT,PHI,THETA,PSI,SHIFT,SCAL,ERG,NOUT)   
C                                                                      
C PURPOSE: MARKER BASED ALIGNMENT - DOUBLE TILTED IMAGES   
C 
C INPUT:
C     P3DREF(3,NTPT)= POINTS IN REFERENCE IMAGE
C     P3D1(3,NTPT)  = POINTS IN VIEWS TO BE ALIGNED, on output
C                     aligned points and the errors
C     ERRORP(NTPT)  = error per point (output)
C     PHI           = ANGLE TO ROTATE VIEW ABOUT Z (IN PLANE)
C     THETA         = TILT ANGLE OF VIEW (ABOUT Y AXIS)
C     PSI           = ANGLE TO TURN ABOUT NEW Z AXIS
C     SCAL          = MULTIPLICATIVE SCALING FACTOR
C  PARAMETERS:                                                         
C                                                                      
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE MRERROR (P3DREF,P3D1,ERRORP,NTPT,
     &                    PHI,THETA,PSI,SHIFT,SCAL,ERG,NOUT)

      REAL       :: P3DREF(3,NTPT),P3D1(3,NTPT),ERRORP(NTPT)
      REAL       :: SHIFT(3),CNEW(3)

      CZ  = COS(PHI)
      SZ  = SIN(PHI)
      CY  = COS(THETA)
      SY  = SIN(THETA)
      CNZ = COS(PSI)
      SNZ = SIN(PSI)

      WRITE(NOUT,*)  ' '
      WRITE(NOUT,*)  ' REFERENCE POINT COORDINATES:'
      WRITE(NOUT,401)  (M,(P3DREF(I,M),I=1,3),M=1,NTPT)
401   FORMAT(1X,I5,2X,3F12.3)

      WRITE(NOUT,*)  ' '
      WRITE(NOUT,*)  ' ALIGNED POINT COORDINATES AND THE ERRORS:'

      ERG = 0.0
      DO M=1,NTPT
         CNEW(1) = SCAL*((CNZ*CZ*CY - SZ*SNZ)*P3D1(1,M)
     &    + (SZ*CY*CNZ + SNZ*CZ)*P3D1(2,M) - CNZ*SY*P3D1(3,M))

         CNEW(2) = SCAL*(-(SZ*CNZ + CZ*SNZ*CY)*P3D1(1,M)
     &    + (CZ*CNZ - SNZ*CY*SZ)*P3D1(2,M) +SNZ*SY*P3D1(3,M))

         CNEW(3) =
     &     SCAL*(SY*CZ*P3D1(1,M) + SY*SZ*P3D1(2,M) + CY*P3D1(3,M))

         ERT = 0.0

         DO J=1,3
            ERT = ERT+(CNEW(J)-P3DREF(J,M))**2

C           OVERWRITE INPUT POINTS WITH ALIGNED POINTS
            P3D1(J,M) = CNEW(J)
         ENDDO

         ERRORP(M) = SQRT(ERT)

         WRITE(NOUT,402) M, CNEW, ERRORP(M)
402      FORMAT(1X,I5,2X,4F12.3)

         ERG = ERG + ERT
      ENDDO

      ERG = SQRT(ERG/NTPT)

      WRITE(NOUT,*)  ' '
      WRITE(NOUT,*)  ' TOTAL ERROR:',ERG
      WRITE(NOUT,*)  ' '

      END
