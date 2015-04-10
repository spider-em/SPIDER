C ++********************************************************************
C                                                                      *
C   VONEA                                                              *
C                                                                      *
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
C   SUBROUTINE VONEA
C
C   PURPOSE: 
C
C   VARIABLES: 
C               
C **********************************************************************

      SUBROUTINE VONEA()

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      LOGICAL                         :: HALF,NEWFILE,DEBUGGING
      REAL, ALLOCATABLE               :: THETA(:),PHI(:)
      REAL, DIMENSION(6)              :: DLIST
      CHARACTER (LEN=MAXNAM)          :: DOCNAM
      CHARACTER (LEN=1)               :: ANSW

      REAL, PARAMETER        :: QUADPI = 3.14159265358979323846 
      REAL, PARAMETER        :: RAD_TO_DGR = (180.0/QUADPI)

      DEBUGGING = .FALSE.

      CALL  RDPRI1S(NPTS,NOT_USED,'NUMBER OF PROJECTIONS',IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      CALL RDPRMC(ANSW,NLET,.TRUE.,'HALF SPHERE (Y/N)',NULL,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      HALF = (ANSW .EQ. 'Y')
      IF (HALF) NPTS = 2 * NPTS
 
      ALLOCATE(THETA(NPTS),PHI(NPTS),STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
         CALL  ERRT(46,'VONEA',NPTS*2)
         GOTO 9999
      ENDIF

      NDOCT = 99
      CALL OPENDOC(DOCNAM,.TRUE.,NLET,NDOCT,NDOC,.TRUE.,
     &           'ANGULAR DOCUMENT',.FALSE.,.TRUE.,.TRUE.,
     &            NEWFILE,IRTFLG)
      IF (IRTFLG .EQ. -1) GOTO 9999

      CALL LUNDOCPUTCOM(NDOC,' PSI, THETA, AND PHI',IRTFLG)

      PSI  = 0.0

      SVAL =  3.6 / SQRT(FLOAT(NPTS))
      PHIT =  0.0

      PHI(1)   =  0.0
      THETA(1) =  0.0

      DO K = 2 , NPTS-1
         ZVAL = -1 + 2 * FLOAT(K) / FLOAT(NPTS)
         RVAL = SQRT(1 - ZVAL * ZVAL)
         PHIT = PHIT + SVAL / RVAL

         XD   = COS(PHIT) * RVAL
         YD   = SIN(PHIT) * RVAL
         ZD   = ZVAL

#if defined (SP_GFORTRAN)
         PHI  (NPTS - K + 1) = ATAN2(XD,YD) * RAD_TO_DGR
         THETA(NPTS - K + 1) = ACOS(ZVAL)   * RAD_TO_DGR
#else
         PHI  (NPTS - K + 1) = ATAN2D(XD,YD)
         THETA(NPTS - K + 1) = ACOSD(ZVAL)
#endif
         IF (DEBUGGING)WRITE(6,90) ZVAL,RVAL,PHIT,XD,YD,PHI(K),THETA(K)
90       FORMAT(3(F8.2,' '),' - ',2(F8.2,' ')' - ',3(F8.2,' '))
      ENDDO

      PHI(NPTS)   =  0.0
      THETA(NPTS) =  180

       NGOT = NPTS

       IF (HALF) THEN
C         REMOVE BOTTEM HALF OF SPHERE
          NGOT = 0
          DO I = 1,NPTS
             IF ( THETA(I) .GT. 90.0) CYCLE
             NGOT        = NGOT + 1
             PHI(NGOT)   = PHI(I)
             THETA(NGOT) = THETA(I)
          ENDDO
       ENDIF

      DO KEY = 1,NGOT
C         WRITE(6,91)PHI(KEY), THETA(KEY)
91        FORMAT(2(F8.2,' '),' - ',3(F8.2,' ')' ',3(F8.2,' '))

          DLIST(1) = PSI
          DLIST(2) = THETA(KEY)
          DLIST(3) = PHI(KEY)

         CALL LUNDOCWRTDAT(NDOC,KEY,DLIST,3,IRTFLG)

          IF (DEBUGGING) THEN
c             WRITE(6,92) KEY, PSI, THETA(KEY), PHI(KEY)
92           FORMAT(' Psi, Theta & Phi:',i3,'): ',
     &                F8.2,', ',F8.2,', ',F8.2)
          ENDIF
       ENDDO

       IF (HALF) THEN
          WRITE(NOUT,93) NGOT
93        FORMAT(' ACTUAL PROJECTIONS IN HEMISPHERE: ',I7)
       ENDIF

       CALL REG_SET_NSEL(1,1,FLOAT(NGOT),0.0,0.0,0.0,0.0,IRTFLG)

       IRTFLG = 0
9999   IF (ALLOCATED(THETA)) DEALLOCATE(THETA)
       IF (ALLOCATED(PHI))   DEALLOCATE(PHI)
       CLOSE(NDOCT)

       END



