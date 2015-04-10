
C ++********************************************************************
C                                                                      *
C                                                                      *
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
C                                                                      *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE  DISCRIM(LUN51,MXM,NMAX,NSUM,MGR,IP,N,JG,XNG,XNV,
     1  M1,B,C,XMEAN,TMEAN,JV,MD,A,AR,BSIG,D,
     2  X,T,IHISTI,E,IQ,XQ,VV,MDT,LEST)

        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*2 (I-N)
	INTEGER*4 LUN51
       DIMENSION A(M1,M1),AR(NMAX,MD),B(M1,M1),C(M1,M1),D(M1,M1)
       DIMENSION N(NMAX),JG(NMAX),XMEAN(M1,MGR),VV(MXM)
       DIMENSION TMEAN(M1),T(M1),X(M1)
       DIMENSION JV(M1),BSIG(M1,M1)
       DIMENSION E(MGR),IHISTI(MGR,MGR),XQ(M1),IQ(M1)
        CHARACTER*4  XNG(MGR),XNV(M1)
	CHARACTER*1 NULL

	NULL=CHAR(0)
       DO  I=1,MGR
       DO  J=1,MGR
       IHISTI(J,I)=0
       ENDDO
       ENDDO
       DF=NSUM-MGR
       DO  I=1,NMAX
       IF(E(I).EQ.0.0) E(I)=REAL(N(I))/REAL(NSUM)
       E(I)=2.*DLOG(E(I))
       ENDDO
C       WRITE(LUN51,98)
 98    FORMAT(//,1X,29('-'),'DISCRIMINANT ANALYSIS',29('-'))
 99    FORMAT(//,1X,40('- '))
       DO  I=1,M1
       DO  J=I,M1
       A(I,J)=C(I,J)-B(I,J)
       BSIG(I,J)=B(I,J)
       BSIG(J,I)=B(J,I)
       A(J,I)=A(I,J)
       ENDDO
       ENDDO
       CALL DIRNM(A,M1,B,D,T,XQ,IQ)
       TRACE=0.0
       DO  J=1,MD
       TRACE=TRACE+T(J)
       ENDDO
C       WRITE(LUN51,101)
C 101   FORMAT(/,
C     & 48H CHI-SQUARE TESTS WITH  SUCCESSIVE ROOTS REMOVED /)
C       WRITE(LUN51,102)
C 102   FORMAT(//,1X,13HROOTS REMOVED,2X,10HEIGENVALUE,2X,
C     & 11HCANONICAL R,2X,6HLAMBDA,2X,13HPERCENT TRACE,2X,
C     & 5HCHI^2,2X,2HDF,2X,4HSIGN /)
       J=MD
       VV(J+1)=1.
 22    VV(J)=VV(J+1)/(1.+T(J))
       J=J-1
       IF(J.GT.0) GOTO 22
       DO  J=1,MD
       JT=J-1
       RC=SQRT(T(J)/(1.0+T(J)))
       PR=100.*T(J)/TRACE
       CHI=NSUM-1-(M1+MGR)/2
       CHI=-CHI*DLOG(VV(J))
       NDF=(M1-JT)*(MGR-JT-1)
       ALPH=CHISQ(CHI,NDF)
C 14    WRITE(LUN51,103) JT,T(J),RC,VV(J),PR,CHI,NDF,ALPH
C 103   FORMAT(I8,7X,G10.4,3X,G10.4,2X,F6.2,2X,
C     & F10.2,4X,F7.2,1X,I2,1X,F5.3)
C       WRITE(LUN51,99)
       ENDDO
       EN=NSUM-1
       DO  J=1,MD
       DO  K=1,M1
       A(J,K)=0.
       DO  L=1,M1
       A(J,K)=A(J,K)+(D(L,J)*C(L,K))/EN
       ENDDO
       ENDDO
       ENDDO
       DO  J=1,MD
       T(J)=0.
       DO  L=1,M1
       T(J)=T(J)+A(J,L)*D(L,J)
       ENDDO
       ENDDO
       WRITE(LUN51,305)
 305   FORMAT(//' UNSTANDARDIZED DISCRIMINANT',
     & ' FUNCTION COEFFICIENTS'/)
       CALL PRT(LUN51,M1,MD,JV,XNV,D)
       DO  I=1,MD
       SUM=0.
       DO  J=1,M1
       SUM=SUM+D(J,I)*TMEAN(J)
       ENDDO
       X(I)=-SUM
       ENDDO
C       WRITE(LUN51,252)(X(K),K=1,MD)
C 252   FORMAT(10H CONSTANT  ,1X,6(1X,G10.4))
       WRITE(LUN51,306)
 306   FORMAT(//' STANDARDIZED DISCRIMINANT ',
     & 'FUNCTION COEFFICIENTS'/)
       DO    I=1,MD
       XX=SQRT(T(I))
       DO    J=1,M1
       D(J,I)=D(J,I)/XX
       ENDDO
       ENDDO
       DO  J=1,M1
       XX=SQRT(C(J,J)/EN)
       DO  I=1,MD
       A(J,I)=D(J,I)*XX
       ENDDO
       ENDDO
       CALL PRT(LUN51,M1,MD,JV,XNV,A)
       WRITE(LUN51,99)
C       WRITE(LUN51,280)
C 280  FORMAT(//' CENTROIDS OF GROUPS IN REDUCED SPACE'/)
       DO 270 J=1,NMAX
       IF(N(J).EQ.0) GOTO 271
       DO  L=1,MD
       SUM=0.
       DO  I=1,M1
       SUM=SUM+D(I,L)*(XMEAN(I,J)-TMEAN(I))
       ENDDO
       AR(J,L)=SUM
       ENDDO
       GOTO 270
 271   DO  I=1,MD
       AR(J,I)=0.
       ENDDO
 270   CONTINUE
C PRT1 had a printout removed, so it is commented out
C       CALL PRT1(MD)
C       WRITE(LUN51,99)
       IF(MD.EQ.MDT)  GOTO  203
       WRITE(LUN51,204)  MDT
 204   FORMAT(/,' ATTENTION ! ONLY ',I2,' DISCRIMINANT FUNCTION(S)'/
     * '     WILL BE USED DURING CLASSIFICATION PROCESS')
       MD=MDT
 203   CONTINUE
       WRITE(LUN51,200)
 200   FORMAT(//' TERITORIAL MAP OF STANDARDIZED ' ,
     & 'DISCRIMINANT SPACE  '   /)
       IF (MD.LE.1) GOTO 300
       WRITE(LUN51,201)
 201   FORMAT(//' HORIZONTAL 1-FUNCTION'    ,
     & ' VERTICAL 2-FUNCTION  '   /)
       CALL FGR(LUN51,IP,M1,D,XNG,NMAX,TMEAN,JG,N,
     & MD,D,AR,JV,VV,MXM,E,IHISTI,XQ,MDT,LEST)
       GOTO 401
  300  CONTINUE
       WRITE(LUN51,301)
 301   FORMAT(/' HISTOGRAM  OF  SUBJECTS  COUNTS   '/)
       CALL HGR(LUN51,IP,M1,D,XNG,NMAX,TMEAN,JG,
     & AR,JV,VV,MXM,E,IHISTI,XQ,LEST)
 401   CONTINUE
       WRITE(LUN51,99)
       EN=NSUM
       WRITE(LUN51,610)
 610   FORMAT(//' PREDICTON RESULTS'/)
       WRITE(LUN51,611)
 611   FORMAT(/' ACTUAL GRP.  N OF CASES',5X,'PREDICTION',
     & ' GROUP MEMBERSHIP '/)
       WRITE(LUN51,631)     (JG(K),XNG(K),K=1,NMAX)
       DO   I=1,NMAX
       T(I)=0.0
       DO   J=1,NMAX
       T(I)=T(I)+DBLE(IHISTI(I,J))
       ENDDO
       ENDDO
       WRITE(LUN51,634)     (T(I),I=1,NMAX)
       DO    I=1,NMAX
       T(I)=T(I)*100.0/DBLE(NSUM)
       ENDDO
       WRITE(LUN51,632)     (T(I),I=1,NMAX)
       DO  J=1,NMAX
       WRITE(LUN51,630)   JG(J),XNG(J),N(J),(IHISTI(I,J),I=1,NMAX)
       DO  I=1,NMAX
       T(I)=REAL(IHISTI(I,J))*100.0/EN
       ENDDO 
       WRITE(LUN51,632)   (T(K),K=1,NMAX)
       ENDDO
       CHI=0.0
       DO  652  I=1,NMAX
       DO  652  J=1,NMAX
       IF(I-J)  653,654,653
 654   CHI=CHI+DBLE((IHISTI(I,J)-N(I))**2)/DBLE(N(I))/2.0
       GOTO  652
 653   CHI=CHI+(DBLE(IHISTI(I,J)-N(I))/2.0/DBLE(NMAX-1))**2/DBLE(N(I))
     * /2.0/DBLE(NMAX-1)
 652   CONTINUE
       IF(NMAX.EQ.2)  NDF=2
       IF(NMAX.GT.2)  NDF=(NMAX-2)*NMAX
       ALPH=CHISQ(CHI,NDF)
       IHISTI(1,NMAX)=0
       DO   I=1,NMAX
       IHISTI(1,NMAX)=IHISTI(1,NMAX)+IHISTI(I,I)
       ENDDO
       WRITE(LUN51,635)     CHI,NDF,ALPH
       WRITE(LUN51,636)     IHISTI(1,NMAX)
       SUM=DBLE(IHISTI(1,NMAX))*100.0/EN
       WRITE(LUN51,637)     SUM,NSUM
 630   FORMAT(//I4,1X,A4,I12,4X,10(4X,I5))
 631   FORMAT(//25X,10(I4,1X,A4))
 632   FORMAT(25X,10(F5.1,' %  '))
 634   FORMAT(25X,10(1X,F5.0,4X))
 635   FORMAT(//,' CHI-SQUARE TEST OF THE CLASSIFICATION PROBABILITY'/
     *  ,' CHI-SQR.=',F10.3,'  WITH ',I5,'  NDF,  SGFN=',F6.4)
 636   FORMAT(' NUMBER OF PROPERLY CLASSIFIED CASES IS ',I5)
 637   FORMAT(' THIS IS ',F5.1,'% OF THE TOTAL NUMBER (',I5,')  OF',
     *                  ' CASES')
       WRITE(LUN51,98)
        ENDFILE  2
         RETURN
       END
