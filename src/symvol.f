C **********************************************************************
c  REPS.F
C                 OMP FAILS ON SARK (ONLY)        AUG 04 ArDean Leith
C++*********************************************************************
C *
C *  SYMVOL
C *
C *  Symmetrize volume in real space
C *  Q1- input volume
C *  Q2 - output volume
C *  RM - matrices with symmetry operations (as calculated using BLDR or BUILDS)
C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2002, P. A. Penczek                                   *
C=*                                                                    *
C=* University of Texas - Houston Medical School                       *
C=*                                                                    *
C=* Email:  pawel.a.penczek@uth.tmc.edu                                *
C=*                                                                    *
C=* This program is free software; you can redistribute it and/or      *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* This program is distributed in the hope that it will be useful,    *
C=* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C=* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
C=* General Public License for more details.                           *
C=*                                                                    *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program; if not, write to the                      *
C=* Free Software Foundation, Inc.,                                    *
C=* 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.      *
C=*                                                                    *
C=**********************************************************************
C *
C **********************************************************************

         SUBROUTINE SYMVOL(Q1,Q2, KLX,KNX,KLY,KNY,KLZ,KNZ, RM,NSYM)

         DIMENSION  Q1(KLX:KNX,KLY:KNY,KLZ:KNZ)
         DIMENSION  Q2(KLX:KNX,KLY:KNY,KLZ:KNZ)
	 DIMENSION  RM(3,3,NSYM)

         DOUBLE     PRECISION  QR(3),QRT(3), DX,DY,DZ

         INTEGER, PARAMETER :: NSIZE = 27

         INTEGER, DIMENSION(27) :: X,Y,Z

#ifndef SP_MP
         REAL, DIMENSION(27) :: F
#endif


C        SET THE KNOWN COORDINATE GRID
C  Replaced by loops below, data does not agree with openmp.
c         DATA X/  
c     &          -1, 0, 1, -1, 0, 1, -1, 0, 1, 
c     &          -1, 0, 1, -1, 0, 1, -1, 0, 1, 
c     &          -1, 0, 1, -1, 0, 1, -1, 0, 1/ 

c         DATA Y/ 
c     &          -1,-1,-1,  0, 0, 0,  1, 1, 1, 
c     &          -1,-1,-1,  0, 0, 0,  1, 1, 1, 
c     &          -1,-1,-1,  0, 0, 0,  1, 1, 1/ 
 
c         DATA Z/  
c     &          -1,-1,-1, -1,-1,-1, -1,-1,-1, 
c     &           0, 0, 0,  0, 0, 0,  0, 0, 0,
c     &           1, 1, 1,  1, 1, 1,  1, 1, 1/

C      write(6,*) ' symvol -1'

C        SET THE KNOWN COORDINATE GRID
	 DO  L=1,NSIZE,3
	    X(L)   = -1
	    X(L+1) = 0
	    X(L+2) = 1
	    Y(L)   = MOD(L/3,3)-1
	    Y(L+1) = MOD(L/3,3)-1
	    Y(L+2) = MOD(L/3,3)-1
	 ENDDO

	 DO  L=1,NSIZE
	    Z(L) = (L-1) / 9-1
	 ENDDO

C        Calculate radius within which the rotation is valid
	 IRADI = (MIN(IABS(KLX),IABS(KNX),IABS(KLY),IABS(KNY),
     &		      IABS(KLZ),IABS(KNZ))-1)**2

c        write(6,*) ' symvol 0',klz,knz,kly,kny,klx,knx

#if defined (__ia64)
C        CAN NOT RUN IN PARALLEL ON ALTIX (CRASHES) 2004
#else
c$omp    parallel do private(IZ,IY,IX)
#endif
         DO IZ=KLZ,KNZ
           DO IY=KLY,KNY
              DO IX=KLX,KNX
		Q2(IX,IY,IZ) = 0.0
              ENDDO
           ENDDO
         ENDDO

C       write(6,*) ' symvol 111',nsym

	DO  ISYM=1,NSYM

           IF (RM(1,1,ISYM) .EQ. 1.0 .AND.
     &	       RM(2,2,ISYM) .EQ. 1.0 .AND. RM(3,3,ISYM) .EQ. 1.0)  THEN


C             NO ROTATION NEEDED
#if defined (__ia64)
C             CAN NOT RUN IN PARALLEL ON ALTIX (CRASHES)
#else
c$omp         parallel do private(IZ,IY,IX)
#endif
              DO IZ=KLZ,KNZ
                 DO IY=KLY,KNY
                    DO IX=KLX,KNX
                      Q2(IX,IY,IZ) = Q2(IX,IY,IZ) + Q1(IX,IY,IZ)
                    ENDDO
                 ENDDO
              ENDDO

	   ELSE   ! END OF:  IF (RM(1,1,ISYM) .EQ. 1.0 .AND. true
C             ROTATION NEEDED

#ifdef SP_MP
#if defined (__ia64)
C             CAN NOT RUN IN PARALLEL ON ALTIX (CRASHES)
#else
c$omp        parallel do private(IZ,IY,IX, QR,QRT, ICRD)
#endif

             DO IZ=KLZ,KNZ
                DO IY=KLY,KNY
                   QRT(1) = RM(2,1,ISYM)*IY + RM(3,1,ISYM)*IZ
                   QRT(2) = RM(2,2,ISYM)*IY + RM(3,2,ISYM)*IZ
                   QRT(3) = RM(2,3,ISYM)*IY + RM(3,3,ISYM)*IZ

                   DO IX=KLX,KNX
		      ICRD = IX*IX + IY*IY + IZ*IZ
		      IF (ICRD .LE. IRADI)  THEN
                         QR(1) = QRT(1) + RM(1,1,ISYM)*IX
                         QR(2) = QRT(2) + RM(1,2,ISYM)*IX
                         QR(3) = QRT(3) + RM(1,3,ISYM)*IX

C                        EVALUATE INTENSITY AT PX,PY,PZ
		         Q2(IX,IY,IZ) = Q2(IX,IY,IZ)+ 
     &	                   TRIQDF(QR,Q1,KLX,KNX,KLY,KNY,KLZ,KNZ,X,Y,Z)
		      ELSE
C                        ROTATED POSITION IS OUTSIDE VOLUME
                         Q2(IX,IY,IZ) = Q2(IX,IY,IZ)+ Q1(IX,IY,IZ)
		      ENDIF
                   ENDDO  ! END OF: DO IX=KLX,KNX
                ENDDO     ! END OF: DO IY=KLY,KNY
             ENDDO        ! END OF: DO IZ=KLZ,KNZ
           ENDIF

C          write(6,*) ' symvol 3'
#else
             DO IZ=KLZ,KNZ
                DO IY=KLY,KNY

                   QRT(1) = RM(2,1,ISYM)*IY+RM(3,1,ISYM)*IZ
                   QRT(2) = RM(2,2,ISYM)*IY+RM(3,2,ISYM)*IZ
                   QRT(3) = RM(2,3,ISYM)*IY+RM(3,3,ISYM)*IZ

                   DO IX=KLX,KNX
		      ICRD = IX*IX+IY*IY+IZ*IZ
		      IF (ICRD.LE.IRADI)  THEN
                         QR(1) = QRT(1) + RM(1,1,ISYM)*IX
                         QR(2) = QRT(2) + RM(1,2,ISYM)*IX
                         QR(3) = QRT(3) + RM(1,3,ISYM)*IX

C                        IOX..  INTEGER LOCATION IN -NSAM/2...NSAM/2 ARRAY
                         IOX = FLOOR(QR(1))   
                         IOY = FLOOR(QR(2))   
                         IOZ = FLOOR(QR(3))   

C                        DX.. OFFSET FROM INTEGER ARRAY
                         DX  = QR(1) - IOX
                         DY  = QR(2) - IOY
                         DZ  = QR(3) - IOZ
C                        ROTATED POSITION IS INSIDE OF VOLUME

C                        FIND INTENSITIES ON 3x3x3 COORDINATE GRID
                         DO L = 1,NSIZE
                            F(L) = Q1(IOX + X(L),IOY + Y(L),IOZ + Z(L))
                         ENDDO

C                        EVALUATE INTENSITY AT PX,PY,PZ
                       Q2(IX,IY,IZ) = Q2(IX,IY,IZ) + TRIQUAD(DX,DY,DZ,F)

                      ELSE
C                        ROTATED POSITION IS OUTSIDE VOLUME
                         Q2(IX,IY,IZ) = Q2(IX,IY,IZ) + Q1(IX,IY,IZ)
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
#endif
       ENDDO


#if defined (__ia64)
C             CAN NOT RUN IN PARALLEL ON ALTIX (CRASHES)
#else
c$omp  parallel do private(IZ,IY,IX)
#endif

         DO IZ=KLZ,KNZ
           DO IY=KLY,KNY
              DO IX=KLX,KNX
		Q2(IX,IY,IZ)=Q2(IX,IY,IZ)/NSYM
              ENDDO
           ENDDO
         ENDDO
        END


#ifdef SP_MP

	 FUNCTION  TRIQDF(QR,Q1,KLX,KNX,KLY,KNY,KLZ,KNZ,X,Y,Z)

         INTEGER, PARAMETER :: NSIZE = 27
         DOUBLE     PRECISION  QR(3),DX,DY,DZ
         DIMENSION  Q1(KLX:KNX,KLY:KNY,KLZ:KNZ)
         INTEGER, DIMENSION(27) :: X,Y,Z
	 DIMENSION  F(NSIZE)

C        IOX..  INTEGER LOCATION IN -NSAM/2...NSAM/2 ARRAY
         IOX = FLOOR(QR(1))
         IOY = FLOOR(QR(2))
         IOZ = FLOOR(QR(3))   

C        DX.. OFFSET FROM INTEGER ARRAY
         DX  = QR(1) - IOX
         DY  = QR(2) - IOY
         DZ  = QR(3) - IOZ

C        FIND INTENSITIES ON 3x3x3 COORDINATE GRID
         DO L = 1,NSIZE
            F(L) = Q1(IOX + X(L),IOY + Y(L),IOZ + Z(L))
         ENDDO

         TRIQDF = TRIQUAD(DX,DY,DZ,F)

	 END
#endif
