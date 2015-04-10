C++*********************************************************************
C
C RAMP_P.F
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
C  RAMP_P(LUN1,LUN2,NX,NY,NOUT)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

         SUBROUTINE  RAMP_P(LUN1,LUN2,NX,NY,NOUT)

         DIMENSION        X(NX)
         EXTERNAL         BETAI
         DOUBLE PRECISION BETAI
         DOUBLE PRECISION C,D,EPS,B1,B2,A,F,R2,DN1,DN2
         DOUBLE PRECISION Q(6),S(9),QYX1,QYX2,QX1X2
     &                       ,QX1,QX2,QY,SYX1,SYX2,SX1X2,SX1
     &                       ,SX2,SY,SX1Q,SX2Q,SYQ
         EQUIVALENCE (Q(1),QYX1),(Q(2),QYX2),(Q(3),QX1X2),(Q(4),QX1),
     &               (Q(5),QX2),(Q(6),QY)
         EQUIVALENCE (S(1),SYX1),(S(2),SYX2),(S(3),SX1X2),(S(4),SX1),
     &               (S(5),SX2),(S(6),SY),(S(7),SX1Q),
     &               (S(8),SX2Q),(S(9),SYQ)

         DATA  EPS/1.0D-5/

C        ZERO ARRAY S
         S   = 0
         N1  = NX / 2
         N2  = NY / 2
         SX1 = FLOAT(N1) * FLOAT(NX + 1)

         IF (MOD(NX,2) .EQ. 1)   SX1 = SX1 + 1 + N1
         SX2 = FLOAT(N2) * FLOAT(NY + 1)

         IF (MOD(NY,2) .EQ. 1)   SX2 = SX2 + 1 + N2
         SX1   = SX1 * NY
         SX2   = SX2 * NX
         SX1X2 = 0.0D0

         DO  J = 1, NY
           CALL REDLIN(LUN1,X,NX,J)

           DO I = 1, NX
             SYX1 = SYX1 + X(I) * I
             SYX2 = SYX2 + X(I) * J
             SY   = SY   + X(I)
             SX1Q = SX1Q + I * I
             SX2Q = SX2Q + J * J
             SYQ  = SYQ  + X(I) * DBLE(X(I))
           END DO
         END DO

         DN    = FLOAT(NX) * FLOAT(NY)
         QYX1  = SYX1 - SX1 * SY / DN
         QYX2  = SYX2 - SX2 * SY / DN
         QX1X2 = 0.0
         QX1   = SX1Q - SX1 * SX1 / DN
         QX2   = SX2Q - SX2 * SX2 / DN
         QY    = SYQ  - SY  * SY  / DN
         C     = QX1  * QX2 - QX1X2 * QX1X2

         IF (C .GT. EPS) THEN
           B1  = (QYX1 * QX2 - QYX2 * QX1X2) / C
           B2  = (QYX2 * QX1 - QYX1 * QX1X2) / C
           A   = (SY - B1 * SX1 - B2 * SX2)  / DN
           D   = B1 * QYX1 + B2 * QYX2
           R2  = D / QY
           DN1 = 2
           DN2 = DN - 3

           IF (DABS(QY - D) .LT. EPS / 100.0) THEN
              F = 0.0
              P = 0.0
           ELSE
              F = D * (DN - 3.0) / 2 /(QY - D)
              P = 2.0*BETAI(0.5D0 * DN2, 0.5D0 * DN1, DN2 / 
     &             (DN2 + DN1 * F)) 

              IF (P.GT.1.0)  P = 2.0 - P
C +
C     &    (1.0D0-BETAI(0.5D0 * DN1, 0.5D0 * DN2, DN1 / (DN1 + DN2 / F)))
           END IF

           WRITE(NOUT,2020)  A, B1, B2, DSQRT(R2), R2, F, DN2, P
2020       FORMAT(/,
     &     '  Ramp model:    Y = a + b1 * x1 + b2 * x2',/,
     &     '  a  = ',1PD12.5,/,
     &     '  b1 = ',1PD12.5,/,
     &     '  b2 = ',1PD12.5,/,
     &     '  Multiple correlation R = ',0PF10.8,/,
     &     '  R squared              = ',0PF10.8,/,
     &     '  F-statistics  F = ',1PD12.5,
     &     '  with n1=2 and n2=' ,0PF7.0,'  df',/,
     &     '  Significance  p = ',0PF10.8)

           D = A + B1 + B2

           DO I = 1, NY
             QY = D

             CALL REDLIN(LUN1,X,NX,I)

             DO  K = 1, NX
                X(K) = X(K) - QY
                QY   = QY + B1
             END DO

             CALL WRTLIN(LUN2,X,NX,I)

             D = D + B2
           END DO

         ELSE
           WRITE(NOUT,3030)
3030       FORMAT(/,' No solution - image is not modified !')

           DO I = 1,NY
              CALL REDLIN(LUN1,X,NX,I)
              CALL WRTLIN(LUN2,X,NX,I)
           END DO
         END IF

         END


C **********************************************************************
C
C  RAMP_PB(X,NX,NY,SAYSTATS,NOUT)
C
C--*********************************************************************

         SUBROUTINE RAMP_PB(X,NX,NY,SAYSTATS,NOUT)

         IMPLICIT NONE

         INTEGER          :: NX,NY,NOUT
         REAL             :: X(NX,NY)
         LOGICAL          :: SAYSTATS

         EXTERNAL         :: BETAI
         DOUBLE PRECISION :: BETAI
         INTEGER          :: N1,N2,J,I,K,IX,IY
         REAL             :: DN,P

         DOUBLE PRECISION :: C,D,B1,B2,A,F,R2,DN1,DN2
         DOUBLE PRECISION :: Q(6),S(9),QYX1,QYX2,QX1X2
     &                           ,QX1,QX2,QY,SYX1,SYX2,SX1X2,SX1
     &                           ,SX2,SY,SX1Q,SX2Q,SYQ

         EQUIVALENCE (Q(1),QYX1),(Q(2),QYX2),(Q(3),QX1X2),(Q(4),QX1),
     &               (Q(5),QX2),(Q(6),QY)
         EQUIVALENCE (S(1),SYX1),(S(2),SYX2),(S(3),SX1X2),(S(4),SX1),
     &               (S(5),SX2),(S(6),SY),(S(7),SX1Q),
     &               (S(8),SX2Q),(S(9),SYQ)

         DOUBLE PRECISION, PARAMETER  :: EPS = 1.0D-5

C        ZERO ARRAY S
         S   = 0

         N1  = NX / 2
         N2  = NY / 2

         SX1 = FLOAT(N1) * FLOAT(NX + 1)
         IF (MOD(NX,2) ==  1)   SX1 = SX1 + 1 + N1

         SX2 = FLOAT(N2) * FLOAT(NY + 1)
         IF (MOD(NY,2) ==  1)   SX2 = SX2 + 1 + N2

         SX1   = SX1 * NY
         SX2   = SX2 * NX
         SX1X2 = 0.0D0

         DO  J = 1, NY
           DO I = 1, NX
             SYX1 = SYX1 + X(I,J) * I
             SYX2 = SYX2 + X(I,J) * J
             SY   = SY   + X(I,J)
             SX1Q = SX1Q + I * I
             SX2Q = SX2Q + J * J
             SYQ  = SYQ  + X(I,J) * DBLE(X(I,J))
           ENDDO
         ENDDO

         DN    = FLOAT(NX) * FLOAT(NY)
         QYX1  = SYX1 - SX1 * SY / DN
         QYX2  = SYX2 - SX2 * SY / DN
         QX1X2 = 0.0
         QX1   = SX1Q - SX1 * SX1 / DN
         QX2   = SX2Q - SX2 * SX2 / DN
         QY    = SYQ  - SY  * SY  / DN
         C     = QX1  * QX2 - QX1X2 * QX1X2

         IF (C <= EPS) THEN
            WRITE(NOUT,*) ' No ramp solution - image not modified!'
            RETURN
         ENDIF
         
         B1  = (QYX1 * QX2 - QYX2 * QX1X2) / C
         B2  = (QYX2 * QX1 - QYX1 * QX1X2) / C
         A   = (SY - B1 * SX1 - B2 * SX2)  / DN
         D   = B1 * QYX1 + B2 * QYX2
         R2  = D / QY
         DN1 = 2
         DN2 = DN - 3

         IF (DABS(QY - D) < EPS / 100.0) THEN
            F = 0.0
            P = 0.0
         ELSE
            F = D * (DN - 3.0) / 2 /(QY - D)
            P = 2.0*BETAI(0.5D0 * DN2, 0.5D0 * DN1, DN2 / 
     &             (DN2 + DN1 * F)) 
            IF (P > 1.0)  P = 2.0 - P
         ENDIF

         IF (SAYSTATS) THEN
            WRITE(NOUT,2020)  A, B1, B2, DSQRT(R2), R2, F, DN2, P
2020        FORMAT(/,
     &     '  Ramp model:    Y = a + b1 * x1 + b2 * x2',/,
     &     '  a  = ',1PD12.5,/,
     &     '  b1 = ',1PD12.5,/,
     &     '  b2 = ',1PD12.5,/,
     &     '  Multiple correlation R = ',0PF10.8,/,
     &     '  R squared              = ',0PF10.8,/,
     &     '  F-statistics  F = ',1PD12.5,
     &     '  with n1=2 and n2=', 0PF7.0,'  df',/,
     &     '  Significance  p = ',0PF10.8)
         ENDIF

         D = A + B1 + B2

         DO IY = 1, NY
           QY = D

           DO IX = 1,NX
#ifdef NEVER
              if (ix == 1  .and. iy == 1  .or.
     &            ix == 1  .and. iy == ny .or.
     &            ix == nx .and. iy == 1  .or.
     &            ix == nx .and. iy == ny) then
                 write(6,*) 'x,y,val:',ix,iy,x(ix,iy), x(ix,iy)-qy
              endif
#endif
              X(IX,IY) = X(IX,IY) - QY
              QY     = QY + B1
           ENDDO

           D = D + B2
         ENDDO

         END



C **********************************************************************
C
C  RAMP_PB_MASK(X,XMASK,NX,NY,SAYSTATS,NOUT)
C
C--*********************************************************************

         SUBROUTINE RAMP_PB_MASK(X,XMASK,PLIST, 
     &                           NX,NY,SAYSTATS,NOUT)

         IMPLICIT NONE

         INTEGER          :: NX,NY
         REAL             :: X(NX,NY)
         LOGICAL          :: XMASK(NX,NY)
         REAL             :: PLIST(3,NX*NY)
         LOGICAL          :: SAYSTATS
         INTEGER          :: NOUT

         INTEGER          :: IX,IY,NPTS
         REAL             :: QY

         real             :: vnew

         DOUBLE PRECISION :: A, B, C

         NPTS = 0

         DO  IY = 1, NY
           DO IX = 1, NX
             IF (.NOT. XMASK(IX,IY)) THEN
C               USE THIS PIXEL TO FIND NORMALIZE RAMP
                NPTS          = NPTS + 1
                PLIST(1,NPTS) = IX
                PLIST(2,NPTS) = IY
                PLIST(3,NPTS) = X(IX,IY)
              ENDIF
           ENDDO
         ENDDO
      
         CALL LS_PLANE_FIT(PLIST, NPTS, A,B,C)                         

         IF (SAYSTATS) THEN
            WRITE(NOUT,2020) A, B, C
2020        FORMAT(/,
     &     '  Ramp model:   Z = a * x + b * y + c',/,
     &     '  a  = ',1PD12.5,
     &     '    b  = ',1PD12.5,
     &     '    c  = ',1PD12.5)
         ENDIF

         A = -A
         B = -B
         C = -C

C        APPLY RAMP TO ALL PIXELS
         DO IY = 1, NY

           QY  = IY * B

           DO IX = 1,NX

#ifdef NEVER
              if (ix == 1  .and. iy == 1  .or.
     &            ix == 1  .and. iy == ny .or.
     &            ix == nx .and. iy == 1  .or.
     &            ix == nx .and. iy == ny) then

                  vnew =  X(IX,IY) + IX * A + QY + C  
                  write(6,*) 'x,y,val:',ix,iy,x(ix,iy), vnew
              endif
#endif

              X(IX,IY) = X(IX,IY) + IX * A + QY + C

           ENDDO
         ENDDO

         END



C     LEAST-SQUARES-FIT OF PLANE TO ARBITRARY NUMBER OF (X,Y,Z) POINTS
C     PLANE : Z = AX + BY + C                              

      SUBROUTINE LS_PLANE_FIT(PLIST, NPTS, A, B, C)                         

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NPTS
      REAL                :: PLIST(3,NPTS)
      DOUBLE PRECISION    :: A, B, C 

      INTEGER             :: IPT
                                                                   
      DOUBLE PRECISION    :: D = 0                                                     
      DOUBLE PRECISION    :: E = 0                                                     
      DOUBLE PRECISION    :: F = 0                                                     
      DOUBLE PRECISION    :: G = 0                                                     
      DOUBLE PRECISION    :: H = 0                                                     
      DOUBLE PRECISION    :: I = 0                                                     
      DOUBLE PRECISION    :: J = 0                                                     
      DOUBLE PRECISION    :: K = 0                                                     
      DOUBLE PRECISION    :: L = 0                                                     

      DOUBLE PRECISION    :: DENOM                                                  

      DO  IPT = 1, NPTS

         D = D + PLIST(1,IPT) * PLIST(1,IPT)                       
         E = E + PLIST(1,IPT) * PLIST(2,IPT)                       
         F = F + PLIST(1,IPT)
                                  
         G = G + PLIST(2,IPT) * PLIST(2,IPT)                       
         H = H + PLIST(2,IPT)                                  
         I = I + 1 
                                           
         J = J + PLIST(1,IPT) * PLIST(3,IPT)                         
         K = K + PLIST(2,IPT) * PLIST(3,IPT)                         
         L = L + PLIST(3,IPT)                                    

      ENDDO

      DENOM   = F * F * G - 2 * E * F * H + 
     &          D * H * H + E * E * I - D * G * I

      ! X AXIS SLOPE
      A = (H * H * J - G * I * J + E * I * K + 
     &     F * G * L - H * (F * K  + E * L)) / DENOM 
                                                       
      ! Y AXIS SLOPE                                                       
      B = (E * I * J + F * F * K - D * I * K + 
     &     D * H * L - F * (H * J  + E * L)) / DENOM 
                                                       
      ! Z AXIS INTERCEPT
      C = (F * G * J - E * H * J - E * F * K + 
     &     D * H * K + E * E * L - D * G * L) / DENOM
 
      END






