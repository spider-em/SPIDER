C++*********************************************************************
C
C  CCRS_2.F                      
C      PGI BUG ON: O(I,J) = CTEMP*CONJG(Y(I,J))*CMPLX(COS(ARG),SIN(ARG))
C      PGI BUG                                  FEB 10 2006 ArDean Leith
C      PGI BUG                                  FEB 10 2006 ArDean Leith
C      MOD PGI COMPILER BUG                     FEB 19 2008 ArDean Leith
C      CCRS_2I                                  APR 27 2009 ArDean Leith
C      CCRS_2_PH                                JUN 30 2011 ArDean Leith
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
C  CCRS_2(X,Y,O,LS,NX,NY)
C
C  PURPOSE: CALCULATES CIRCULAR CROSCORRELATION, NON-POWER-OF-TWO 
C           DIMENSIONS
C
C  PARAMETERS:  X,Y  FOURIER TRANSFORMS                    (SENT)
C               O    F(X*CONJG(Y))                         (RET.)
C               LS   NX+2-MOD(NX,2)                    (SENT)
C
C NOTE:  USE CCRS_2I WHEN: O IS SAME ADDRESS AS X!!!!!
C        June 2013 CCRS_2I was buggy when used in oracfmskm.f so
C        reverted that use to: CCRS_2 
C        CCRS_2I is still in use in: corr1.f -> mccf.f -> mccf_p
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************
 
        SUBROUTINE CCRS_2(X,Y,O, LS,NX,NY)

C       REAL   X((NX+2-MOD(NX,2)),NY)
C       ABOVE ON PGI COMPILER 7.1 FAILS TO COMPILE PROPERLY SOMETIMES

        REAL  :: X(LS,NY)
        REAL  :: Y(LS,NY)
        REAL  :: O(LS,NY)
 
	PARAMETER (QUADPI = 3.1415926535897932384)
	PARAMETER (PI2    = 2*QUADPI)

        ITMP = NX / 2
        SX   = PI2 * FLOAT(ITMP)/ FLOAT(NX)
        ITMP = NY / 2
        SY   = PI2 * FLOAT(ITMP)/ FLOAT(NY)

C$omp   parallel do private(i,j,iy,argy,arg,tmpr,tmpi)
        DO J=1,NY
           IY = J - 1
           IF (IY > (NY/2)) IY = IY - NY
           ARGY = SY * IY

           DO I=1,LS,2
              ARG      = SX * (I-1)/2 + ARGY

    	      TMPR     = X(I,J)   * Y(I,J)  + X(I+1,J) * Y(I+1,J)
	      TMPI     = X(I+1,J) * Y(I,J)  - X(I,J)   * Y(I+1,J)

	      O(I,J)   = TMPR * COS(ARG) - TMPI * SIN(ARG)
	      O(I+1,J) = TMPI * COS(ARG) + TMPR * SIN(ARG)
           ENDDO
        ENDDO

        INS = -1
        CALL FMRS_2(O,NX,NY,INS)

        END

C       --------------- CCRS_PH_2 -------------------------------------

        SUBROUTINE CCRS_PH_2(X,Y, LS,NX,NY,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'

        REAL             :: X(LS,NY)
        REAL             :: Y(LS,NY)
        INTEGER          :: LS,NX,NY,IRTFLG
 
        REAL             :: SX,SY,ARGY,ARG,TMPR,TMPI,FNRM,BTM
        INTEGER          :: ITMP,K,J,IY,I,INS,NE

        LOGICAL          :: BADBTM 

	REAL, PARAMETER  :: QUADPI = 3.1415926535897932384
	REAL, PARAMETER  :: PI2    = 2*QUADPI

        ITMP   = NX / 2
        SX     = PI2 * FLOAT(ITMP) / FLOAT(NX)
        ITMP   = NY / 2
        SY     = PI2 * FLOAT(ITMP ) / FLOAT(NY)
        BADBTM = .FALSE.

C$omp   parallel do private(i,j,iy,argy,arg,tmpr,tmpi,fnrm,btm,badbtm)
        DO J=1,NY
           IY = J - 1
           IF (IY .GT. (NY/2)) IY = IY - NY
           ARGY = SY * IY

           DO I=1,LS,2
              ARG      = SX * (I-1) / 2 + ARGY

    	      TMPR     = X(I,J)   * Y(I,J)  + X(I+1,J) * Y(I+1,J)
	      TMPI     = X(I+1,J) * Y(I,J)  - X(I,J)   * Y(I+1,J)

	      X(I,  J) = TMPR * COS(ARG) - TMPI * SIN(ARG)
	      X(I+1,J) = TMPI * COS(ARG) + TMPR * SIN(ARG)

              BTM = SQRT(X(I,J)**2 + X(I+1,J)**2)
              IF (BTM .LE. 0.0) THEN
                 BADBTM = .TRUE.
                 EXIT
              ENDIF
                  
	      FNRM = 1.0 / BTM

	      X(I,  J) = X(I,  J) * FNRM
	      X(I+1,J) = X(I+1,J) * FNRM
           ENDDO
        ENDDO

        IF (BADBTM) THEN
           WRITE(NOUT,90) X(I,J), X(I+1,J)
90         FORMAT(' COMPLEX VALUES GIVE DIV. BY ZERO:',1PG14.3)
           CALL ERRT(101,'AVOIDED DIV BY ZERO',NE)
           IRTFLG = 1
           RETURN
        ENDIF
  
        INS = -1
        CALL FMRS_2(X,NX,NY,INS)

        IRTFLG = 0
        END

C      --------------- CCRS_2I -------------------------------------

       SUBROUTINE CCRS_2I(X,Y, LS,NX,NY)

C       REAL   X((NX+2-MOD(NX,2)),NY)

        REAL  :: X(LS,NY)
        REAL  :: Y(LS,NY)
 
	PARAMETER (QUADPI = 3.1415926535897932384)
	PARAMETER (PI2    = 2*QUADPI)

        ITMP = NX / 2
        SX   = PI2 * FLOAT(ITMP) / FLOAT(NX)

        ITMP = NY / 2
        SY   = PI2 * FLOAT(ITMP) / FLOAT(NY)

C$omp   parallel do private(i,j,iy,argy,arg,tmpr,tmpi)
        DO J=1,NY
           IY = J - 1
           IF (IY > (NY/2)) IY = IY - NY
           ARGY = SY * IY

           DO I=1,LS,2
              ARG      = SX * (I-1)/2 + ARGY

    	      TMPR     = X(I,J)   * Y(I,J)  + X(I+1,J) * Y(I+1,J)
	      TMPI     = X(I+1,J) * Y(I,J)  - X(I,J)   * Y(I+1,J)

	      X(I,J)   = TMPR * COS(ARG) - TMPI * SIN(ARG)
	      X(I+1,J) = TMPI * COS(ARG) + TMPR * SIN(ARG)
           ENDDO
        ENDDO

        INV = -1
        CALL FMRS_2(X,NX,NY,INV)

        END


#ifdef NEVER
	      OREAL    = TMPR * COS(ARG) - TMPI * SIN(ARG)
	      OIMG     = TMPI * COS(ARG) + TMPR * SIN(ARG)
	      O(I,J)   = OREAL
              O(I+1,J) = OIMG		 
#endif




#ifdef NEVER
        SUBROUTINE CCRS_PH_2_CMPLX(X,Y, LSC,NX,NY)

        IMPLICIT NONE

        COMPLEX         :: X(LSC, NY)
        COMPLEX         :: Y(LSC, NY)
        COMPLEX         :: CTEMPA,CTEMPB

	REAL, PARAMETER :: QUADPI = 3.141592653589793238462
	REAL, PARAMETER :: PI2=2*QUADPI

        INTEGER         :: LSC,NX,NY,ITMP,J,IY,IX,INS,I
        REAL            :: SX,SY,ARGY,ARG,FVAL
      
        ITMP = NX/2
        SX   = PI2*FLOAT(ITMP)/FLOAT(NX)
        ITMP = NY/2
        SY   = PI2*FLOAT(ITMP)/FLOAT(NY)

        !write(6,*) ' NX:',lsc,NX,NY

C       PHASE CORRELATION
        DO J=1,NY
           IY = J-1
           IF (IY .GT. (NY/2)) IY = IY - NY
           ARGY = SY*IY

           DO I=1,LSC
              IX     = I - 1
              ARG    = SX * IX + ARGY

              CTEMPA = CMPLX(COS(ARG),SIN(ARG))
              CTEMPB = X(I,J) * CONJG(Y(I,J)) * CTEMPA

              X(I,J) = CTEMPB /
     &                 (SQRT(REAL(CTEMPB)**2 + AIMAG(CTEMPB)**2))
           ENDDO
           !write(6,*) ' jfjf:',x(nnnn,j)
        ENDDO

        !call CHKMAX('x',  x,lsc*2*NY)
        !call CHKMin('x',  x,lsc*2*NY)
        !call chkcmplx('x',x,lsc*NY,10,1,lsc)

        INS = -1
        CALL FMRS_2(X,NX,NY,INS)

        !call chkreal('x',x,NX*NY,10,1,NX*NY)

        END
#endif




