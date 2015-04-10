
C **********************************************************************
C
C                        NOUT USED JUNE 2000 ARDEAN LEITH
C
C **********************************************************************
C *  AUTHOR :                                                              *
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
C **********************************************************************

        SUBROUTINE STOCHA(JBASE,  NPIX,  NUMIM,  KFAC,  NITER, 
     &                    LUV, PJ,  S,  U,  TRACE,  SOMP,  V,  BB, 
     &                    D,  AD,  INB, NB, NTMP)

        INCLUDE 'CMBLOCK.INC'

        DIMENSION  S(NPIX, JBASE), PJ(NPIX), U(NPIX), LUV(NUMIM), 
     &           BB(NPIX, JBASE),  V(NPIX),  D(NPIX),  AD(JBASE, JBASE)
        LOGICAL  S_ON_DISK, INB(JBASE),LDUM

C       NAR  =  CENTERING FREQUENCY

        DATA  EPS/1.0E - 4/, NAR/10/, KITER/8/

        S_ON_DISK = .FALSE.
        TRACE     =  0 
 
        DO J  =  1, NPIX                                                        
          V(J)   =  0.0
          PJ(J)  =  0.0
        END DO
        CALL REW(NB,  1)
        SOMP  =  0.0 

C       CALCULATION OF WEIGHTS         
        DO I  =  1, NUMIM                                        
          READ(NB) (U(J), J = 1, NPIX),  PIA,  VAR,  INO

          IF(LUV(I) .NE. 0) THEN
            SOMP  =   SOMP  +  PIA                                              
            DO  J  =  1, NPIX
              PJ(J)  =  PJ(J)  +  U(J)                                    
            END DO

            DO J  =  1, NPIX

C             SET PROTECTION AGAINST UNDERFLOW
c             PJ(J) = AMAX1(PJ(J), 1.0E - 6)

              V(J)  =  V(J)  +  U(J) * U(J) / PIA
            END DO
          END IF
        END DO

        DO J  =  1, NPIX
          TRACE  =  TRACE  +  V(J) / PJ(J)
        END DO
        TRACE  =  TRACE  -  1.0 
        DO J  =  1, NPIX
           PJ(J)  =  PJ(J)/SOMP
        END DO


C	RANDOMIZE S MATRIX.
        DO J  =  1, NPIX
           DO L  =  1, JBASE
              S(J, L)  =  SEN3A(BID)
           END DO
        END DO
                                              
        CALL GSMOD (NPIX,  NPIX,  JBASE,  PJ,  S,  KR,  V)
        
        DO  L = 1, JBASE                    
           INB(L) = .TRUE.
        END DO
        NIT = 0
	LDUM=.TRUE.
        DO WHILE (LDUM)
88        NIT = NIT + 1
          WRITE(NOUT,*) ' ITERATION #:', NIT  
          CALL ITPOW (NUMIM, NPIX, JBASE, NAR, PJ, S, BB, U, V, 
     &                SOMP, TRACE, NB, LUV, INB)

          IF(NIT. GE. KITER - 1) THEN
            IF ((NIT.GT.KITER) .AND. (MOD(NIT - 12, 4) .NE. 0)) THEN
              CONTINUE
            ELSE
              JCASE = 0
              DO   L = 1, JBASE
                IF(INB(L))  JCASE = JCASE + 1
              END DO
              CALL CPROJ (NUMIM,  JBASE,  NPIX,  S,  V,  BB,  U,  AD, 
     &                      SOMP,  PJ,  NB,  LUV, JCASE, INB)
              JCASE = 0
              DO L = 1, JBASE
                IF(INB(L))  THEN
                   JCASE = JCASE + 1
                   D(L) = V(JCASE)
                ENDIF
              END DO
              IF (S_ON_DISK)  THEN
                REWIND NTMP
                AVE = 1.0
                DO  L  =  1, KFAC
                  V(L)  =  0.0                             

C	          U HAS TO HAVE THE AVERAGE SUBSTRATED AND PJ()=1

                  READ(NTMP) (U(J), J = 1, NPIX)
                  DO  J  =  1, NPIX
                     V(L)  =  V(L)  +  PJ(J) * S(J, L) * U(J)
                  END DO
                  IF((1.0 - ABS(V(L))).LT.EPS)  THEN
                     INB(L) = .FALSE.
                  ELSE
                     AVE = AMIN1(AVE, ABS(V(L)))
                  ENDIF
                END DO

                WRITE (NOUT,501) AVE
501             FORMAT('  ***   COSINES OF EIGENVECTORS   *** ', /, 
     &	               '      MINIMUM COSINE            =  ', F8.5)

                WRITE (NOUT,  502) (V(L), L = 1, KFAC)
502             FORMAT(8(2X, F8.5))

                WRITE (NOUT,  505) (INB(L), L = 1, JBASE)
505             FORMAT(8(9X, L1))
    
                IF ((1.0 - AVE) .LT. EPS)  GOTO  80
              ENDIF

              REWIND NTMP
              DO L  =  1, JBASE
                 WRITE(NTMP) (S(J, L), J = 1, NPIX)
              END DO
              S_ON_DISK = .TRUE.
            END IF
          END IF
        END DO

80      CONTINUE
        DO J = 1, NPIX                                      
           PJ(J)  =  PJ(J) * SOMP
        END DO

        END
