C++*********************************************************************
C
C PJ3.F               REWRITTEN                 PAWEL PENCZYK SEPT 2003
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
C IMAGE_PROCESSING_ROUTINE
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE PJ3_N
 
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
        
        CHARACTER(LEN=MAXNAM)               :: FINPIC
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: CUBE
        REAL, ALLOCATABLE, DIMENSION(:,:)   :: B

        DATA  INPIC/60/,IOPIC/61/

        MAXIM=0
        CALL OPFILEC(0,.TRUE.,FINPIC,INPIC,'O',ITYPE,NSAM,NROW,NSLICE,
     &             MAXIM,'INPUT',.FALSE.,IRTFLG)

        IF (IRTFLG .NE. 0)    RETURN

        IF (ITYPE.NE.3)  THEN
           CALL  ERRT(2,'THIS OPERATION ONLY WORKS ON VOLUMES',NE)
           GOTO 9999
        ENDIF

        ALLOCATE(CUBE(NSAM,NROW,NSLICE), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           MWANT = NSAM*NROW*NSLICE
           CALL ERRT(46,'PJ 3, CUBE',MWANT)
           GOTO 9999
        ENDIF
        CALL READV(INPIC,CUBE,NSAM,NROW,NSAM,NROW,NSLICE)

        CALL RDPRMI(NSAMP,NROWP,NOT_USED,'PROJECTION DIMENSIONS X,Y')
        IF (NROWP.EQ.0) NROWP = NINT(NSAMP*(FLOAT(NROW)/FLOAT(NSAM)))
        IF (NROWP.EQ.0) THEN
           CALL  ERRT(31,'PJ 3',NE)
           GOTO 9999
        ENDIF

        MAXIM   = 0
        NSLICEP = 0
        ITYPE   = 1
        CALL OPFILEC(0,.TRUE.,FINPIC,IOPIC,'U',ITYPE,
     &             NSAMP,NROWP,NSLICEP,
     &             MAXIM,'OUTPUT',.FALSE.,IRTFLG)

        PSI = HUGE(PSI)
        CALL RDPRM3S(PHI,THETA,PSI,NOT_USED,
     &                 'PHI, THETA & PSI',IRTFLG)
        IF (IRTFLG.NE.0) GOTO 9999 

        TVAL = HUGE(PSI)
        IF (PSI .EQ. TVAL) THEN 
           CALL RDPRM1S(PSI,NOT_USED,'PSI',IRTFLG)
           IF (IRTFLG.NE.0) GOTO 9999 
        ENDIF

        ALLOCATE (B(NSAMP,NROWP), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           MWANT = NSAMP * NROWP
           CALL ERRT(46,'PJ 3, B',MWANT)
           GOTO 9999
        ENDIF
	B=0.
        CALL PJ33_N(CUBE,B,NSAM,NROW,NSLICE,NSAMP,NROWP,PSI,THETA,PHI)
 
        NSLICEP = 1
        CALL WRITEV(IOPIC,B,NSAMP,NROWP,NSAMP,NROWP,NSLICEP)

9999    IF (ALLOCATED(B))    DEALLOCATE (B)
        IF (ALLOCATED(CUBE)) DEALLOCATE (CUBE)

        CLOSE(IOPIC)
        CLOSE(INPIC)

        END





        SUBROUTINE  PJ33_N(CUBE,B,NSAM,NROW,NSLICE,NSAMP,NROWP,
     &                   PSI,THETA,PHI)

         DIMENSION         CUBE(NSAM,NROW,NSLICE),B(NSAMP,NROWP)
         DIMENSION         DM(9)
         DOUBLE PRECISION  CPHI,SPHI,CTHE,STHE,CPSI,SPSI
         DOUBLE PRECISION  QUADPI,DGR_TO_RAD,RAD_TO_DGR

         PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
         PARAMETER (DGR_TO_RAD = (QUADPI/180))


         CPHI = DCOS(DBLE(PHI)*DGR_TO_RAD)
         SPHI = DSIN(DBLE(PHI)*DGR_TO_RAD)
         CTHE = DCOS(DBLE(THETA)*DGR_TO_RAD)
         STHE = DSIN(DBLE(THETA)*DGR_TO_RAD)
         CPSI = DCOS(DBLE(PSI)*DGR_TO_RAD)
         SPSI = DSIN(DBLE(PSI)*DGR_TO_RAD)

         DM(1) = CPHI*CTHE*CPSI-SPHI*SPSI
         DM(2) = SPHI*CTHE*CPSI+CPHI*SPSI
         DM(3) = -STHE*CPSI
         DM(4) = -CPHI*CTHE*SPSI-SPHI*CPSI
         DM(5) = -SPHI*CTHE*SPSI+CPHI*CPSI
         DM(6) = STHE*SPSI

         DM(7) = STHE*CPHI
         DM(8) = STHE*SPHI
         DM(9) = CTHE

         DM1 = DM(1)
         DM2 = DM(2)
         DM3 = DM(3)

C        ZERO THE WHOLE B ARRAY
c        B = 0.0

        LDPX = NSAM/2+1
        LDPY = NROW/2+1
        LDPZ = NSLICE/2+1
        MPPX = NSAMP/2+1
        MPPY = NROWP/2+1

        DO K=1,NSLICE
           DO  J=1,NROW
               XB = (1-LDPX)*DM(1) + (J-LDPY)*DM(4) +
     &              (K-LDPZ)*DM(7) + LDPX

               YB = (1-LDPX)*DM(2) + (J-LDPY)*DM(5) +
     &              (K-LDPZ)*DM(8) + LDPY

               ZB = (1-LDPX)*DM(3) + (J-LDPY)*DM(6) +
     &              (K-LDPZ)*DM(9) + LDPZ
               J2 = J-LDPY+MPPY
               IF (.NOT.(J2.LT.1 .OR. J2.GT.NROWP))  THEN
                   DO  I=1,NSAM
                     IOX    = IFIX(XB)
                     IF(.NOT.(IOX.LT.1 .OR. IOX.GE.NSAM))  THEN
                        IOY    = IFIX(YB)
                        IF (.NOT.(IOY.LT.1 .OR. IOY.GE.NROW))  THEN
                            IOZ    = IFIX(ZB)
                            IF (.NOT.(IOZ.LT.1 .OR. IOZ.GE.NSLICE)) THEN
                               I2=I-LDPX+MPPX
                               IF (.NOT.(I2.LT.1 .OR. I2.GT.NSAMP))THEN

                                  DX = XB-IOX
                                  DY = YB-IOY
                                  DZ = ZB-IOZ

                                  A1 = CUBE(IOX,IOY,IOZ)
                                  A2 = CUBE(IOX+1,IOY,IOZ) - A1
                                  A3 = CUBE(IOX,IOY+1,IOZ) - A1
                                  A4 = CUBE(IOX,IOY,IOZ+1) - A1
                                  A5 = -A2 - CUBE(IOX,IOY+1,IOZ) + 
     &                                       CUBE(IOX+1,IOY+1,IOZ)
                                  A61= - CUBE(IOX,IOY,IOZ+1) + 
     &                                   CUBE(IOX+1,IOY,IOZ+1)
                                  A6 = -A2 + A61
                                  A7 = -A3 - CUBE(IOX,IOY,IOZ+1) + 
     &                                       CUBE(IOX,IOY+1,IOZ+1)
                                  A8 = -A5 - A61 -
     &                                  CUBE(IOX,IOY+1,IOZ+1) + 
     &                                  CUBE(IOX+1,IOY+1,IOZ+1)
                                  B(I2,J2) = B(I2,J2)+
     &                               A1 + DZ*(A4 + A6*DX + 
     &                              (A7 + A8*DX)*DY) + A3*DY +
     &                              DX*(A2 + A5*DY)

                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                        XB = XB + DM1
                        YB = YB + DM2
                        ZB = ZB + DM3
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
	END

