
C
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
C
        
        SUBROUTINE ITPOW (NUMIM, NPIX, JBASE, NAR, PJ, S, BB, U, V,
     &                  SOMP, TRACE, NB, LUV, inb)
                                                                               
        DIMENSION PJ(NPIX),  S(NPIX, JBASE),  BB(NPIX, JBASE)
        DIMENSION U(NPIX),  V(NPIX),  LUV(NUMIM)
        logical  inb(jbase)
        LOGICAL CENTER
                                                                               
        CALL REW(NB,  1)
        DO  L  =  1, JBASE
          if(inb(l))  then
            DO  J  =  1, NPIX
              BB(J, L)  =  S(J, L)
              S(J, L)  =  0.0
            END DO
          endif
        END DO
                                                                               
        DO  IA  =  1, NUMIM
          IF(LUV(IA) .EQ. 0)  THEN
            READ(NB)
          ELSE
            READ(NB) (U(K),  K = 1, NPIX),  PIA,  VAR, INO
            CENTER = MOD(IA, NAR) .EQ. 0 .OR. IA .EQ. NUMIM
            DO L  =  1, JBASE
              if(inb(l))  then
                T1  =  0.0                                                
                DO K  =  1, NPIX
                  T1  =  T1  +  BB(K, L) * U(K) / PIA
                END DO
                                                                              
                DO K  =  1, NPIX
                  S(K, L)  =  S(K, L)  +  U(K) * T1/PJ(K) /SOMP
                END DO
                IF(CENTER)  THEN

C........ PERIODIC CENTERING
                  T2  =  0.0                                                
                  DO JP  =  1, NPIX
                    T2  =  T2  +  PJ(JP) * S(JP, L)
                  END DO
                  DO J  =  1, NPIX
                    S(J, L)  =  S(J, L)  -  T2                            
                  END DO
                ENDIF
              endif
            END DO
          ENDIF
        END DO
        DOP  =  TRACE/2/(NPIX  -  1)
                                                                             
        DO  L  =  1, JBASE
          if(inb(l))  then
            DO  J  =  1, NPIX
              S(J, L)  =  S(J, L)  -  DOP * BB(J, L)
            END DO
          endif
        END DO
        do l = 1, jbase
          if(.not.inb(l))  then
            CALL GSMODl(NPIX,  JBASE,  PJ,  S,  KRANG, v, inb)
            return
          endif
        enddo
        CALL GSMOD(NPIX, NPIX,  JBASE,  PJ,  S,  KRANG, v)
        RETURN
        END            
