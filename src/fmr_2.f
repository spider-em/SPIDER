C++*********************************************************************
C
C FMR_2.F  should be rewritten!!
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
C  For order of elements see fmr_1.
C
C CALLERS: crsm_2.f
C             SHOULD BE REPLACED!! BUT NOT WORTH THE EFFORT??
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

#ifdef SP_MP
        SUBROUTINE  FMR_2(X,NX,NY,DUMMY,INV)

        PARAMETER (LENBUF=8200)
        DIMENSION  X(NX,NY),WORK(LENBUF),DUMMY(*)
        LOGICAL*1  IFND

        IF (MAX(NX,2*NY) > LENBUF)  THEN
           CALL  ERRT(6,' FMR_2 ',NE)
           INV = 0
           RETURN
        ENDIF
#else

        SUBROUTINE  FMR_2(X,NX,NY,WORK,INV)
        DIMENSION   X(NX,NY),WORK(*)
        LOGICAL*1   IFND
#endif

        IFND = MOD(NX,2).EQ.0

        IF (IFND)  THEN
           LBD=2
        ELSE
           LBD=1
        ENDIF

        INS=INV*NX

C       work(max0(NX,2*NY))

        IF (INV.GE.0)  THEN
c$omp      parallel do private(i,work),shared(invt) 
           DO I=1,NY
              INVT=INV
              CALL  FMR_1(X(1,I),NX,WORK,INVT)
           ENDDO

           IF (INVT.LE.0)  THEN
              INV=INVT
              RETURN
           ENDIF
        ENDIF

c$omp   parallel do  private(i,j,work),shared(invt) 
        DO I=1,LBD
           DO J=1,NY
              WORK(NY+J)=X(I,J)
           ENDDO

           INVT = INV
           CALL  FMR_1(WORK(NY+1),NY,WORK,INVT)

           DO J=1,NY
              X(I,J)=WORK(NY+J)
           ENDDO
        ENDDO

        IF (INVT .EQ. 0)  THEN
           INV = 0
           RETURN
        ENDIF
c$omp   parallel do  private(i),shared(invt) 
        DO I=3,NX-1,2
           INVT = INS
           CALL FFTMCF(X(I,1),X(I+1,1),NY,NY,NY,INVT)
        ENDDO

        IF (INVT .EQ. 0)  THEN
           INV = 0
           RETURN
        ENDIF

        IF (.NOT.IFND)  CALL FFTMCF (X(NX,1),X(2,1),NY,NY,NY,INS)
        IF (INV .EQ. 1) RETURN

C       NORMALIZE FOR INVERSE
        Q = 1.0 / FLOAT(NY)

c$omp   parallel do  private(j,i) 
        DO J=1,NY
           DO I=LBD+1,NX
              X(I,J)=X(I,J)*Q
           ENDDO
        ENDDO

c$omp   parallel do  private(i,work),shared(invt) 
        DO I=1,NY
           INVT = INV
           CALL FMR_1(X(1,I),NX,WORK,INVT)
        ENDDO

        IF (INVT .LE. 0)  INV = INVT

        END
