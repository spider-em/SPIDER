
C **********************************************************************
C
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
C
C **********************************************************************

        SUBROUTINE CPROJ(NUMIM, JBASE, NPIX, S, D, BB, U, AD, 
     $                    SOMP, PJ, NB, LUV,JCASE, INB)

C	OPERATION OF PROJECTION AND DIAGONALIZATION, FOR DIAGO-
C       NALIZATION BY DIRECT READING.

        DIMENSION S(NPIX, JBASE),  BB(NPIX, JBASE),  AD(JCASE, JCASE), 
     &              U(NPIX),  PJ(NPIX),  D(NPIX),  LUV(NUMIM)
        LOGICAL  INB(JBASE)

        CALL REW(NB,  1)
        DO  M  =  1, JCASE
          DO L  =  1, JCASE
            AD(L, M)  =  0.0
          END DO
        END DO

        DO  I  =  1, NUMIM                                     
          IF(LUV(I) .EQ. 0) THEN
            READ(NB)
          ELSE
C
C	U90 HAS TO HAVE AVERAGE SUBTRACTED
C
            READ (NB) (U(KK),  KK = 1, NPIX) , PIA,  VAR,  INO
C
C	PQ HAS TO BE = 1
C
            PQ  =  1.0 / (SOMP * PIA)
            MC = 0
            DO  M  =  1, JBASE
              IF(INB(M))  THEN
                MC = MC + 1
                LC = 0
                DO L  =  1, M
                  IF(INB(L))  THEN
                    LC = LC + 1
                    CIL  =  0.0
                    CIM  =  0.0
                    DO K  =  1, NPIX
C
C	IF PJ() CONTAINS AVERAGES , PIA HAS TOP BE =1 
C
C
                      UUU  =  U(K)  -  PJ(K) * PIA
                      CIL  =  CIL  +  S(K, L) * UUU
                      CIM  =  CIM  +  S(K, M) * UUU
                    END DO
                    AD(LC, MC)  =  AD(LC, MC)  +  PQ * CIL * CIM
                  ENDIF
                END DO
              ENDIF
            END DO
          ENDIF
        END DO
                                                                               
        DO M  =  1, JCASE                                                       
          DO L  =  1, M                                                   
            AD(M, L)  =  AD(L, M)
          END DO
        END DO
        
        CALL VPROP (JCASE ,  JCASE ,  AD ,  D ,  U ,  KOD)

        MC = 0
        DO  M  =  1, JBASE
          IF(INB(M))  THEN
            MC = MC + 1
            DO J  =  1, NPIX
              BB(J, M)  =  0.0
              KC = 0
              DO  K  =  1, JBASE
                IF(INB(K))  THEN
                  KC = KC + 1
                  BB(J, M)  =  BB(J, M)  +  S(J, K) * AD(KC, MC)
                ENDIF
              END DO
            END DO
          ENDIF
        END DO

        DO  L  =  1, JBASE
          IF(INB(L))  THEN
            DO J  =  1, NPIX
              S(J, L)  =  BB(J, L)
            END DO
          ENDIF
        END DO
        RETURN
        END
