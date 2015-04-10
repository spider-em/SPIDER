
C ++********************************************************************
C                                                                      *
C SEEDS.F                                                              *
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
C                                                                      *
C  PURPOSE: CHOOSE NK EQUIDISTANT OBJECTS FOR SEEDS PROTOTYPES                                                           *
C                                                                      *
C  PARAMETERS:                                                         *
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE SEEDS(CIRSEED,CIRC,DIST,NK,LCIRC,IP,NIMA,NOUT)

        DIMENSION  CIRC(LCIRC,NIMA),CIRSEED(LCIRC,NK),DIST(NIMA)
        INTEGER*2  IP(NK)


C       CHOOSE NK EQUIDISTANT OBJECTS FOR SEEDS PROTOTYPES

        DMAX = DIST(1)
        DMIN = DMAX
        MAX  = 1
        MIN  = 1

        DO I=2,NIMA
           IF (DIST(I) .GT. DMAX)  THEN
              DMAX = DIST(I)
              MAX  = I
           ENDIF
           IF (DIST(I) .LT. DMIN)  THEN
              DMIN = DIST(I)
              MIN  = I
           ENDIF
        ENDDO

        LK     = 1
        IP(LK) = MIN
        WRITE(NOUT,4994) LK,IP(LK)

        DO    I=1,LCIRC
           CIRSEED(I,1 )= CIRC(I,MIN)
        ENDDO

C2      CIRSEED(I,2) = CIRC(I,MAX)

        LS = NIMA/NK
        DM = DMIN
        DO LK=2,NK
           DB = DMAX
           DO    J=1,LS
              DO    L=1,NIMA
                 IF (DIST(L).GT. DM .AND. DIST(L) .LT. DB)  THEN
                    DB    = DIST(L)
                    IP(LK) = L
                 ENDIF
              ENDDO
              DM=DB
           ENDDO

           WRITE(NOUT,4994)  LK,IP(LK)
4994       FORMAT('  Seed #:',I4,' = object #:',I5)

           DO I=1,LCIRC
              CIRSEED(I,LK) = CIRC(I,IP(LK))
           ENDDO
        ENDDO
        END
