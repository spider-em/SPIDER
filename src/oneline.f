C ++********************************************************************
C
C  ONELINE.F                                        05/03/02
C                                         
C **********************************************************************
C *  WIW3D.f
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
C **********************************************************************
C                                                                      *
C  ONELINE(J,N,N2,X,W,BI,DM)                                           *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE  ONELINE(J,N,N2,X,W,BI,DM)

        DIMENSION      W(0:N2,N,N)
        COMPLEX        BI(0:N2,N),X(0:N2,N,N),BTQ
        DIMENSION      DM(6)

        PARAMETER      (LTAB=4999)
        COMMON  /TABS/ LN2,FLTB,TABI(0:LTAB)

        IF (J .GE. 0)  THEN
           JP = J+1
        ELSE
           JP = N+J+1
        ENDIF

        DO  I=0,N2
           IF (((I*I+J*J) .LT.  (N*N/4)) .AND..NOT. 
     &          (I.EQ. 0  .AND. J.LT.0)) THEN
              XNEW = I * DM(1) + J * DM(4)
              YNEW = I * DM(2) + J * DM(5)
              ZNEW = I * DM(3) + J * DM(6)

              IF (XNEW .LT. 0.0)  THEN
                 XNEW = -XNEW
                 YNEW = -YNEW
                 ZNEW = -ZNEW
                 BTQ  = CONJG(BI(I,JP))
              ELSE
                 BTQ  = BI(I,JP)
              ENDIF

              IXN = IFIX(XNEW+0.5+N) - N
              IYN = IFIX(YNEW+0.5+N) - N
              IZN = IFIX(ZNEW+0.5+N) - N

              IF (IXN .LE. (N2-LN2-1)  .AND.
     &            IYN .GE. (-N2+2+LN2) .AND. IYN .LE. (N2-LN2-1) .AND.
     &            IZN .GE. (-N2+2+LN2) .AND. IZN .LE. (N2-LN2-1)) THEN

                 IF (IXN .GE. 0) THEN
C                   MAKE SURE THAT LOWER LIMIT FOR X DOES NOT GO BELOW 0
                    LB = -MIN0(IXN,LN2)
                    DO LZ=-LN2,LN2
                       IZP = IZN + LZ
                       IF(IZP .GE. 0) THEN
                          IZA = IZP + 1
                       ELSE
                          IZA = N + IZP + 1
                       ENDIF
 
                       TZ  = TABI(NINT(ABS(ZNEW-IZP) * FLTB))

                       IF (TZ .NE. 0.0)  THEN
                          DO  LY=-LN2,LN2
                             IYP = IYN + LY
                             IF (IYP .GE .0) THEN
                                IYA = IYP + 1
                             ELSE
                                IYA = N + IYP + 1
                             ENDIF
 
                             TY  = TABI(NINT(ABS(YNEW-IYP) * FLTB)) * TZ
                             IF (TY .NE. 0.0)  THEN
                                DO  IXP=LB+IXN,LN2+IXN

C                                  GET THE WEIGHT
                                   WG=TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY
                                   IF (WG .NE. 0.0) THEN

                                      X(IXP,IYA,IZA) =
     &                                    X(IXP,IYA,IZA) + BTQ * WG
                                      W(IXP,IYA,IZA) = 
     &                                    W(IXP,IYA,IZA) + WG
                                   ENDIF
                                ENDDO
                             ENDIF
                          ENDDO
                       ENDIF
                   ENDDO
                ENDIF

C               ADD REFLECTED POINTS
                IF (IXN .LT. LN2) THEN
                   DO  LZ=-LN2,LN2
                      IZP = IZN + LZ
                      IZT =  - IZP + 1
                      IF (IZP .GT. 0)  IZT = N + IZT

                      TZ = TABI(NINT(ABS(ZNEW-IZP) * FLTB))

                      IF (TZ .NE. 0.0)  THEN
                         DO  LY=-LN2,LN2
                            IYP = IYN + LY
                            IYT = -IYP + 1
                            IF (IYP .GT. 0) IYT = IYT + N

                            TY = TABI(NINT(ABS(YNEW-IYP) * FLTB)) * TZ
                            IF (TY .NE. 0.0)  THEN
                               DO  IXP=IXN-LN2,-1

C                                 GET THE WEIGHT
                                  WG = TABI(NINT(ABS(XNEW-IXP)*FLTB))*TY

                                  IF (WG .NE. 0.0)  THEN
                                     X(-IXP,IYT,IZT) = 
     &                                   X(-IXP,IYT,IZT) + CONJG(BTQ)*WG
                                     W(-IXP,IYT,IZT) =
     &                                   W(-IXP,IYT,IZT) + WG
                                  ENDIF
                               ENDDO
                            ENDIF
                         ENDDO
                      ENDIF
                   ENDDO
                ENDIF
              ENDIF
           ENDIF
C          END J-I LOOP
        ENDDO

        END
