
C ++********************************************************************
C                                                    
C  SHIFT_PF                                     
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
C  SHIFT_PF(X,NNNN,NX,NY,NZ,SX,SY,SZ)                                    
C          
C  PURPOSE:    SHIFT (2)3-D IN FOURIER SPACE
C              IF SX AND SY AND SZ EQUAL ZERO THEN NO SHIFT
C
C  PARAMETERS:                       
C                                                                      
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE  SHIFT_PF(X,NNNN,NX,NY,NZ,SX,SY,SZ)

        IMPLICIT NONE

	COMPLEX           :: X(NNNN,NY,NZ)
	INTEGER           :: NNNN,NX,NY,NZ
	REAL              :: SX,SY,SZ
         
	DOUBLE PRECISION  :: PI2
        INTEGER           :: K,IZ,J,IY,I,IX
        REAL              :: ARGZ, ARGY, ARG

	IF (SX == 0.0 .AND. SY == 0.0 .AND. SZ == 0.0)  RETURN

	PI2 = 8.0 * DATAN(1.0D0)

c$omp   parallel do private(i,j,k,arg,argy,argz,iz,iy,ix)
	DO K=1,NZ

	   IZ = K - 1
	   IF (IZ > NZ/2)    IZ = IZ - NZ
	   ARGZ = PI2 * SZ * IZ / FLOAT(NZ)

	   DO J=1,NY
	      IY = J - 1
	      IF (IY > NY/2)  IY = IY - NY
	      ARGY = PI2 * SY * IY / FLOAT(NY) + ARGZ

	      DO I=1,NNNN
	         IX       = I - 1
	         ARG      = PI2 * SX * IX / FLOAT(NX) + ARGY
	         X(I,J,K) = X(I,J,K) * CMPLX(COS(ARG),SIN(ARG))
	      ENDDO
	   ENDDO
	ENDDO

	END
