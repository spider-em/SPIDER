C+*********************************************************************
C
C BOXX.F 
C         NORM3 CALL NOT INVOKED FOR MODE 2 BUG  MAR 01 ARDEAN LEITH
C         SETPRM REMOVED                         MAY 09 ARDEAN LEITH
C         MENU                                   FEB 12 ARDEAN LEITH
C         ANCIENT 3D SHIFT BUG & LOWPASS F BUG   AUG 12 ARDEAN LEITH
C           
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C  BOXX(LUN1,LUN2,NX,NY,NZ,NDUM)
C
C  PARAMETERS:
C        LUN1     LOGICAL UNIT NUMBER OF FILE
C        LUN2     LOGICAL UNIT NUMBER OF FILE
C        NX       NUMBER OF SAMPLES
C        NY       NUMBER OF ROWS
C        NZ       NUMBER OF SLICES
C        NDUM     UNUSED
C
C  PURPOSE: DOES OPERATIONS INVOLVING  LOCAL BOX
C         MODE 1 > SUBTRACTIVE BOX FILTERING        (HIGH PASS)
C         MODE 2 > RETURNS FMIN IF PT <  XAVG
C                          FMAX IF PT >= XAVG
C         MODE 3 > RETURNS LOCAL AVERAGE FOR POINT  (LOW PASS)
C         MODE 4 > DIVISIVE CONTRAST CORRECTION
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE BOXX(LUN1,LUN2, NX,NY,NZ, NDUM)

        INCLUDE 'CMBLOCK.INC'

        INTEGER            :: LUN1,LUN2,NX,NY,NZ,NDUM

        INTEGER            :: DS,DR,DL
        CHARACTER (LEN=1)  :: ANS
        CHARACTER (LEN=1)  :: NULL = CHAR(0)
        REAL, ALLOCATABLE  :: X(:),XOUT(:)


        WRITE(NOUT,92)
 92     FORMAT(
     &      ' .MENU: H   -- HIGH PASS FILTER'/
     &      '        L   -- LOW  PASS FILTER'/
     &      '        T   -- THRESHOLD FILTER'/
     &      '        D   -- DIVISIVE LOCAL FILTER'/)

        CALL RDPRMC(ANS,NC,.TRUE.,'FILTER OPTION (H,L,T,D)',
     &              NULL,IRTFLG)

        MODE = INDEX('HTLD',ANS)
        IF (MODE <= 0)  THEN
           CALL ERRT(101,'UNKNOWN FILTER OPTION',NERROR)
           RETURN
        ENDIF

        IF (MODE == 4) THEN 
C          DIVISIVE LOCAL FILTER  
           IF (IMAMI.NE.1) CALL NORM3(LUN1,NX,NY,NZ, FMAX,FMIN,AV)
           AV1 = AV * 0.1   ! why??
        ENDIF

        FMINN = FMIN
        FMAXX = FMAX
        IF (FMIN >= FMAX) THEN
           FMINN = 0.0
           FMAXX = 2.0
        ENDIF

        IF (NZ == 1 ) THEN
           DR = 0
           CALL RDPRI2S(DS,DR, NOT_USED,
     &               'BOX COLUMNS & ROWS',IRTFLG)
           IF (DS <= 0) THEN
             CALL ERRT(101,'BOX COLUMNS MUST BE > 0 AND ODD',NERROR)
             RETURN
           ENDIF
           IF (DR == 0) DR = DS

           IF     (MOD(DS,2) == 0) THEN
             CALL ERRT(101,'BOX COLUMNS MUST BE ODD',DS)
             RETURN
           ELSEIF (MOD(DR,2) == 0) THEN
             CALL ERRT(101,'BOX ROWS MUST BE ODD',   DR)
             RETURN
           ENDIF

           WRITE(NOUT,91) DS,DR
 91        FORMAT('  BOX SIZE: ',I3,' x',I3)

           ALLOCATE(X(NX*NY),XOUT(NX),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              MWANT = NX*NY + NX
              CALL ERRT(46,'BOXX; X & XOUT',MWANT)
              GOTO 9999
           ENDIF
  
        ELSE

           DR = 0
           DL = -999
           CALL RDPRI3S(DS,DR,DL, NOT_USED,
     &                  'BOX COLUMNS, ROWS, & SLICES',IRTFLG)
           IF (DS <= 0) THEN
              CALL ERRT(101,'BOX COLUMNS MUST BE > 0 AND ODD',NERROR)
              RETURN
           ENDIF
           IF (DR == 0) DR = DS
           IF (DL == 0) DL = DS

           IF (DL == -999) THEN
C             MUST FIND SLICE DIMENSION OF BOX  
              CALL RDPRI1S(DL, NOT_USED,'BOX SLICES',IRTLFG)
              IF (DL == -999) DL = DS
           ENDIF

           IF     (MOD(DS,2) == 0) THEN
             CALL ERRT(101,'BOX COLUMNS MUST BE ODD',DS)
             RETURN
           ELSEIF (MOD(DR,2) == 0) THEN
             CALL ERRT(101,'BOX ROWS MUST BE ODD',   DR)
             RETURN
           ELSEIF (MOD(DL,2) == 0) THEN
             CALL ERRT(101,'BOX SLICES MUST BE ODD', DL)
             RETURN
           ENDIF

           WRITE(NOUT,90) DS,DR,DL
 90        FORMAT('  BOX SIZE: ',I3,' x',I3,' x',I3)

           ALLOCATE(X(NX*NY*DL),XOUT(NX),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              MWANT = NX*NY*DL + NX
              CALL ERRT(46,'BOXX; X & XOUT',MWANT)
              GOTO 9999
           ENDIF
        ENDIF

        F = 1.0
        CALL RDPRM1S(F,NOT_USED,'FILTER WEIGHT (0.0 ... 1.0)',IRTLFG)
        FC = 1.0 - F

        IF (NZ .LE. 1)  THEN
C------------------------------------------------------ 2-D -----------

           KT     = 1
           KR     = 1
           MOVWAY = 2
           I      = -1

           DO J=1,NY
              JOFF = MOD(J-1,DR) * NX
              KT   = 1
              IF (MOVWAY == 3) KT = NX
              NO   = KT

C             *** START OFF OR GO DOWN ***
              IF (MODE == 1) THEN
                 CALL AVERG(LUN1,XAVG, J,KT, NX,NY, DS,DR,KR,X)

                 POINT    = X(JOFF+KT)
                 XOUT(NO) = F * (POINT-XAVG) + POINT * FC
 

              ELSEIF (MODE == 2) THEN
                 CALL AVERG(LUN1,XAVG, J,KT, NX,NY, DS,DR,KR,X)

                 POINT    = X(JOFF+KT)
                 XOUT(NO) = F * FMINN + FC * POINT
                 IF (POINT >= XAVG) XOUT(NO) = F * FMAXX + FC * POINT
 
              ELSEIF (MODE == 3) THEN
                 CALL AVERG(LUN1,XAVG, J,KT, NX,NY,DS,DR,KR,X)

                 POINT    = X(JOFF+KT)
                 XOUT(NO) = XAVG * F + FC * POINT
             
              ELSEIF (MODE == 4) THEN
                 CALL AVERG(LUN1,XAVG, J,KT, NX,NY, DS,DR,KR,X)

                 POINT    = X(JOFF+KT)
                 XOUT(NO) = F * POINT / (XAVG+AV1) + POINT * FC
              ENDIF

              DO  K=2,NX
                 KT = K
                 IF (MOVWAY == 3) KT = NX + 1 - K
                 NO = KT

C                *** MOVE RIGHT OR LEFT ***
                 IF (MODE == 1) THEN
                    CALL AVERG(LUN1,XAVG, J,KT, NX,NY,DS,DR,MOVWAY,X)

                    POINT    = X(JOFF+KT)
                    XOUT(NO) = F * (POINT - XAVG) + POINT * FC

                 ELSEIF (MODE == 2) THEN
                    CALL AVERG(LUN1,XAVG,J, KT, NX,NY,DS,DR,MOVWAY,X)

                    POINT    = X(JOFF+KT)
                    XOUT(NO) = F * FMINN + FC * POINT
                    IF (POINT >= XAVG) XOUT(NO) = F*FMAXX + FC*POINT

                 ELSEIF (MODE == 3) THEN
                    CALL AVERG(LUN1,XAVG, J,KT, NX,NY,DS,DR,MOVWAY,X)

                    POINT    = X(JOFF+KT)
                    XOUT(NO) = XAVG * F + FC * POINT

                 ELSEIF (MODE == 4) THEN
                    CALL AVERG(LUN1,XAVG, J,KT, NX,NY,DS,DR,MOVWAY,X)

                    POINT    = X(JOFF+KT)
                    XOUT(NO) = F*POINT / (XAVG+AV1) + POINT*FC
                 ENDIF
              ENDDO

              I      = -I
              MOVWAY = MOVWAY+I
              KT     = 1
              KR     = 4
              CALL WRTLIN(LUN2,XOUT,NX,J)

           ENDDO

C------------------------------------------------------ 3D ------------
        ELSE

           I      = -1
           MOVWAY = 2
           IS     = -1
           MOVSID = 4
           KR     = 1
           LREAD  = 0

           DO L=1,NZ    ! LOOP OVER SLICES
              LOFF = MOD(L-1,DL) * NX * NY    ! 0, 16k, 32k 48k

              DO J=1,NY     ! LOOP OVER ROWS
                 IF (MOVWAY == 3)  THEN
                    KT = NX
                 ELSE
                    KT = 1
                 ENDIF

                 IF (MOVSID == 5) THEN
                    KJ = NY + 1 - J
                 ELSE
                    KJ = J
                 ENDIF

                 NI    = LOFF + (KJ-1)*NX + KT
                 NO    = KT

                 IF (J == 1)  THEN
C                   FIRST LINE OF THIS SLICE

                    IF (MODE == 1)  THEN                   
                       CALL AVERG3(LUN1,XAVG, KJ,KT,L, NX,NY,NZ,
     &                             DS,DR,DL,KR,LREAD,X)
                       POINT    = X(NI)
                       XOUT(NO) = F*(POINT-XAVG) + POINT*FC

                    ELSEIF (MODE == 2)  THEN 
                       CALL AVERG3(LUN1,XAVG, KJ,KT,L, NX,NY,NZ,
     &                             DS,DR,DL,KR,LREAD,X)

                       POINT = X(NI)
                       IF (POINT >= XAVG) XOUT(NO) = F*FMAXX + FC*POINT

                    ELSEIF (MODE == 3)  THEN 
                       CALL AVERG3(LUN1,XAVG, KJ,KT,L, NX,NY,NZ,
     &                             DS,DR,DL,KR,LREAD,X)

                       POINT    = X(NI)
                       XOUT(NO) = XAVG*F + FC*POINT
                       
                    ELSEIF (MODE == 4)  THEN 
                       CALL AVERG3(LUN1,XAVG, KJ,KT,L, NX,NY,NZ,
     &                             DS,DR,DL,KR,LREAD,X)

                       POINT    = X(NI)
                       XOUT(NO) = F*POINT/(XAVG+AV1) + POINT*FC
                    ENDIF 
                                                                               
                 ELSE  
                 
C                   SECOND OR GREATER LINE OF THIS SLICE

                     IF (MODE == 1) THEN   
                       CALL AVERG3(LUN1,XAVG, KJ,KT,L, NX,NY,NZ,
     &                             DS,DR,DL,MOVSID,LREAD,X)

                       POINT    = X(NI)
                       XOUT(NO) = F*(POINT-XAVG) + POINT*FC

                    ELSEIF (MODE == 2) THEN                       
                       CALL AVERG3(LUN1,XAVG, KJ,KT,L, NX,NY,NZ,
     &                             DS,DR,DL,MOVSID,LREAD,X)

                       POINT    = X(NI)
                       XOUT(NO) = F*FMINN + FC*POINT
                       IF (POINT >= XAVG) XOUT(NO) = F*FMAXX + FC*POINT

                    ELSEIF (MODE == 3) THEN 
                       CALL AVERG3(LUN1,XAVG, KJ,KT,L, NX,NY,NZ,
     &                             DS,DR,DL,MOVSID,LREAD,X)

                       POINT    = X(NI)
                       XOUT(NO) = XAVG*F + FC*POINT

                    ELSEIF (MODE == 4) THEN                     
                       CALL AVERG3(LUN1,XAVG, KJ,KT,L, NX,NY,NZ,
     &                            DS,DR,DL,MOVSID,LREAD,X)

                       POINT    = X(NI)
                       XOUT(NO) = F*POINT/(XAVG+AV1) + POINT*FC               
                    ENDIF                    
                 ENDIF
             
                 DO  K=2,NX
                    IF (MOVWAY == 3) THEN
                       KT = NX + 1 - K    
                    ELSE
                       KT = K
                    ENDIF

                    NI    = LOFF + (KJ-1)*NX+KT
                    NO    = KT
                    
                    IF (MODE == 1) THEN  
                       CALL AVERG3(LUN1,XAVG, KJ,KT,L, NX,NY,NZ,
     &                             DS,DR,DL,MOVWAY,LREAD,X)

                       POINT    = X(NI)
                       XOUT(NO) = F*(POINT-XAVG) + POINT*FC
                       
                    ELSEIF (MODE == 2) THEN                                                               
                       CALL AVERG3(LUN1,XAVG, KJ,KT,L, NX,NY,NZ,
     &                             DS,DR,DL,MOVWAY,LREAD,X)

                       POINT    = X(NI)
                       XOUT(NO)  = F*FMINN + FC*POINT
                       IF (POINT >= XAVG) 
     &                     XOUT(NO) = F*FMAXX + FC*POINT

                    ELSEIF (MODE == 3) THEN
                       CALL AVERG3(LUN1,XAVG, KJ,KT,L, NX,NY,NZ,
     &                             DS,DR,DL,MOVWAY,LREAD,X)

                       POINT    = X(NI)
                       XOUT(NO) = XAVG*F + FC*POINT
                    
                    ELSEIF (MODE == 4) THEN                    
                       CALL AVERG3(LUN1,XAVG, KJ,KT,L, NX,NY,NZ,
     &                             DS,DR,DL,MOVWAY,LREAD,X)

                       POINT    = X(NI)
                       XOUT(NO) = F*POINT/(XAVG+AV1) + POINT*FC                    
                    ENDIF
                 ENDDO

                 I      = -I
                 MOVWAY = MOVWAY + I
                 NLIN   = (L-1) * NY + KJ

                 CALL WRTLIN(LUN2,XOUT,NX,NLIN)

              ENDDO

              IS     = -IS
              MOVSID = MOVSID + IS
              KR     = 6
           ENDDO
  
        ENDIF	

9999    IF (ALLOCATED(X))    DEALLOCATE(X)
        IF (ALLOCATED(XOUT)) DEALLOCATE(XOUT)
        END








































