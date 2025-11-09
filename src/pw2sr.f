
C***********************************************************************
C
C PW2SR.F         REFACTORED AND PARALLELIZED      FEB 18 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2018  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C PW2SR(Q,NX,NY,MODE)
C
C PURPOSE: POWER SPECTRUM OF IMAGE, COMPUTED IN-PLACE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

C       -------------------- PW2SR_A ------------------------

        SUBROUTINE  PW2SR_A(Q,NX,LSD,NY,MODE)

C       FOR FOURIER FILE, LSD = NX+2-MOD(NX,2)   eg: 16=18  17=18

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'

        INTEGER      :: NX,NY,LSD !GFORT ERROR IF NOT BEFOR Q,  Sept 2025 al
        REAL         :: Q(LSD,NY)
        CHARACTER    :: MODE

        INTEGER      :: J,I,JN,II,JJ,JB,NXD2,NYD2
        INTEGER      :: LSDM1,LSDD2M1,NYD2PMOD,NXC,NYC,JNL
        REAL         :: SCL,TEMP,SCLSQ,VAL,FVAL,FVALI,FVALR
        REAL         :: FMAXVAL
        LOGICAL      :: EVENY,WANTLOG,WANT2,WANTPH


        REAL, PARAMETER :: THRESH = 1E-9

        REAL, PARAMETER :: PI = 3.141592653589
        REAL, PARAMETER :: DGR_TO_RAD = (PI / 180)

        EVENY    = (MOD(NY,2) == 0)

        LSDM1    = LSD - 1
        LSDD2M1  = LSD/2 - 1

        NXD2     = NX / 2
        NYD2     = NY / 2
        NXC      = NXD2 + 1     ! X CENTER
        NYC      = NYD2 + 1     ! Y CENTER
     
        NYD2PMOD = NYD2 + MOD(NY,2)

C       SCALE FACTOR
        SCL      = 2.0 / FLOAT(NX) / FLOAT(NY)
        SCLSQ    = SCL * SCL

        WANTPH  = (INDEX(MODE,'P') > 0)
        WANTLOG = (INDEX(MODE,'L') > 0)
        WANT2   = (INDEX(MODE,'2') > 0)

C       CONVERT TO AMPLITUDE SPECTRUM

        IF (WANT2)  THEN     
C           SQUARE OF STANDARD POWER SPECTRUM,, UNDOCUMENTED JUN 2011!

C          LOOP OVER ALL ROWS, PUTTING AMPS IN RIGHT HALF OF ROW
c$omp      parallel do private(i,j)
           DO J=1,NY                       ! LOOP OVER ALL ROWS (1:16)
              DO I=LSDM1,1,-2  
                 Q(LSDM1-(LSDM1-I)/2,J) = SCLSQ *
     &             (Q(I+1,J)*Q(I+1,J) + Q(I,J)*Q(I,J))
              ENDDO
           ENDDO

         ELSEIF (WANTLOG) THEN
C          LOG OF STANDARD POWER SPECTRUM

C          LOOP OVER ALL ROWS, PUTTING AMPS IN RIGHT HALF OF ROW
c$omp      parallel do private(i,j)
           DO J=1,NY                       ! LOOP OVER ALL ROWS (1:16)
              DO I=LSDM1,1,-2              ! (17:1 x-2)
                                           ! (17-(17-17)/2) = 17 :17-(17-1)/2=9
                 FVAL = SCL * 
     &                  SQRT(Q(I+1,J)*Q(I+1,J)+Q(I,J)*Q(I,J))

                 IF (FVAL > 0.0) THEN
                    Q(LSDM1-(LSDM1-I)/2,J) = LOG10(FVAL)
                    ! (17:9,1:16) = (17:1x-2, 1:16)
                 ELSE
C                   CAN NOT TAKE LOG!!
                    Q(LSDM1-(LSDM1-I)/2,J) = 0.0
                 ENDIF
              ENDDO
           ENDDO

        ELSEIF (WANTPH) THEN
C          PHASE SPECTRUM
           ! MUST IGNORE NOISE (VERY SMALL NUMBERS) 
           !THRESH = MAX(ABS(Q)) / 10000

#ifdef NEVER
           FMAXVAL = 0.0
           DO J = 1,NY
              DO I = 1,LSD
                 !THRESH = MAX(ABS(Q(I,J)),THRESH)
                 FVAL = ABS(Q(I,J))* SCL
                 IF (FVAL > FMAXVAL) THEN
                   FMAXVAL = FVAL
                   !write(6,'(a,2i7,1Pg12.4)')'  Max value: ',i,j,fmaxval 
                 ENDIF
              ENDDO
           ENDDO
           !write(6,'(a,1Pg12.4)') '  Max value: ',fmaxval
           THRESH = FMAXVAL * 1E-9
#endif

           WRITE(NOUT,'(A,1PG12.4)') '  Zero threshold: ',THRESH
        
c$omp      parallel do private(i,j)
           DO J=1,NY                       
              DO I=LSDM1,1,-2              
                                   
                FVALR = Q(I  ,J) * SCL
                FVALI = Q(I+1,J) * SCL
                IF (FVALR < THRESH) FVALR = 0.0
                IF (FVALI < THRESH) FVALI = 0.0

                Q(LSDM1-(LSDM1-I)/2,J) = ATAN2(FVALI,FVALR) * DGR_TO_RAD
              ENDDO
           ENDDO

         ELSE
C          STANDARD POWER SPECTRUM

C          LOOP OVER ALL ROWS, PUTTING AMPS IN RIGHT HALF OF ROW
c$omp      parallel do private(i,j)
           DO J=1,NY                       ! LOOP OVER ALL ROWS (1:16)
              DO I=LSDM1,1,-2              ! (17:1 x-2)
                                           ! (17-(17-17)/2) = 17 :17-(17-1)/2=9
                 Q(LSDM1-(LSDM1-I)/2,J) =  ! (17:9,1:16)    = (17:1x-2, 1:16)
     &             SCL * SQRT(Q(I+1,J)*Q(I+1,J) + Q(I,J)*Q(I,J))
              ENDDO
           ENDDO
        ENDIF

        !! NYD2PMOD = NYD2 + MOD(NY,2) = 8+0

C       INTERCHANGE TOP AND BOTTEM HALF OF RIGHT HALF COLUMNS          
c$omp   parallel do private(j,jj,i,temp)
        DO J=1,NYD2               ! 1..8           LOOP OVER TOP ROWS
           JJ = J + NYD2PMOD      ! 1+8:8+8 = 9:16 LOOP OVER BOTTEM HALF
           DO I=NXC,LSDM1         ! 9..17          LOOP OVER RIGHT HALF
              TEMP    = Q(I,J)    ! (9:17, 1:8)
              Q(I,J)  = Q(I,JJ)   ! (9:17, 1:8)  = (9:17,9:16)
              Q(I,JJ) = TEMP      ! (9:17, 9:16) = (9:17,1:8)
           ENDDO
c           write(6,'(a, i4,a,  i4,      a, i4,a,2e12.4)') 
c    &              '  ',I, '(',NXC+I/2,',',JN,')',VAL
        ENDDO


        IF (.NOT. EVENY)  THEN    ! ODD HEIGHT

c$omp      parallel do private(i,temp,j)
           DO I=NXC,LSDM1          !  9..16   (for size 15x15)
              TEMP = Q(I,NYD2+1)   ! (9:16,8)
              DO J=NYD2+1,NY-1     !  8..14
                 Q(I,J) = Q(I,J+1) ! (9:16:8:14)
              ENDDO
              Q(I,NY) = TEMP       ! (9:16,15)
          ENDDO
        ENDIF   

C       REVERSE COLS ON RIGHT HALF OF FIRST LINE INTO LEFT HALF
        IF (EVENY) THEN           ! EVEN HEIGHT
           JB = 2                 ! 2
           DO I=1,LSDD2M1         ! 1:8
              II     = LSD - I    ! 18-1:18-8 = 17:10
              Q(I,1) = Q(II,1)    ! (1:8,1)   = (17:10,1) REVERSES 
           ENDDO
        ELSE
           JB = 1
        ENDIF
       
C       REVERSE COLS ON RIGHT HALF OF REMAINING LINES INTO LEFT HALF 

c$omp   parallel do private(j,jj,i,ii)
        DO J=JB,NY                ! 2..16
           JJ = NY-J+JB           ! 16-2+2:16-16+2 = 16:2
           DO I=1,LSDD2M1         ! 1..8
              II     = LSD - I    ! 18-1:18-8  = 17:10
              Q(I,J) = Q(II,JJ)   ! (1:8,2:16) = (17:10,16:2)
           ENDDO
        ENDDO


        !write(6,*)' origin:',NXD2+1,NYD2+1,NXD2,NYD2
        !write(6,*)' origin:',Q(NXD2+1,NYD2+1),Q(NXD2,NYD2)
        !call chkmaxloc('qnyd2-1', q(1:nx,nyd2-1),nx)
        !call chkmaxloc('qnyd2  ', q(1:nx,nyd2),nx)      
        !call chkmaxloc('qnyd2p1', q(1:nx,nyd2+1),nx)
      
C       REPLACE ORIGIN WITH ADJACENT VALUE TO HELP CONTRAST
        Q(NXD2+1, NYD2+1) = Q(NXD2,NYD2)

        !call chkmaxloc('qnyd2  ', q(1:nx,nyd2),nx)      
        !call chkmaxloc('qnyd2p1', q(1:nx,nyd2+1),nx)

        !IF (MODE == 'L') THEN
        !  CALL AAL11(Q,LSD*NY)
        !ENDIF

        END


C       -------------------- PW2SR ------------------------


        SUBROUTINE  PW2SR(Q,NX,NY,MODE)

        IMPLICIT NONE
        INTEGER      :: NX,NY   !GFORT ERROR IF NOT BEFOR Q,  Sept 2025 al
        REAL         :: Q(NX+2-MOD(NX,2), NY)
        CHARACTER*1  :: MODE

        INTEGER      :: NNNN,NSC,J,I,JJ,II,JB
        REAL         :: SCL,TEMP,SCLSQ
        LOGICAL      :: EVENY

        NNNN  = NX + 2 - MOD(NX,2)
        NSC   = NX/2+1
        EVENY = (MOD(NY,2) == 0)

        SCL   = 2.0 / FLOAT(NX) / FLOAT(NY)
        SCLSQ = SCL * SCL

C       CONVERT TO AMPLITUDE SPECTRUM

c$omp   parallel do private(i,j)
        DO J=1,NY
           IF (MODE == '2')  THEN  ! THIS IS NEVER USED JUN 2011!
              DO I=NNNN-1,1,-2
                 Q(NNNN-1-(NNNN-1-I)/2,J)=
     &            SCLSQ * (Q(I+1,J) * Q(I+1,J) + Q(I,J) * Q(I,J))
              ENDDO

           ELSE
              DO I=NNNN-1,1,-2
                 Q(NNNN-1-(NNNN-1-I)/2,J)=
     &            SCL * SQRT(Q(I+1,J) *Q(I+1,J) + Q(I,J) * Q(I,J))
              ENDDO
           ENDIF

C          WANT LOG OF SPECTRUM FOR EASIER VISUALIZATION?
cc         IF (MODE == 'L') CALL AL10(Q(NX/2,J), NX/2) BUGGY IF HERE
        ENDDO


C       MOVE ORIGIN (LOW FREQ) TO CENTER OF SPECTRUM IMAGE

C       INTERCHANGE TOP AND BOTTEM HALF OF RIGHT HALF COLUMNS          
        DO J=1,NY/2
           JJ = J+NY/2+MOD(NY,2)
           DO I=NSC,NNNN-1
              TEMP    = Q(I,J)
              Q(I,J)  = Q(I,JJ)
              Q(I,JJ) = TEMP
           ENDDO
        ENDDO

        IF (.NOT. EVENY) THEN   !ODD ROWS
           DO I=NSC,NNNN-1
              TEMP = Q(I,NY/2+1)
              DO J=NY/2+1,NY-1
                 Q(I,J) = Q(I,J+1)
              ENDDO
              Q(I,NY) = TEMP
          ENDDO
        ENDIF   

C       REVERSE COLS ON RIGHT HALF OF FIRST LINE INTO LEFT HALF
        NSC = NNNN/2-1

        IF (EVENY) THEN       ! EVEN ROWS
           JB = 2
           DO I=1,NSC
              II     = NNNN-I
              Q(I,1) = Q(II,1)
            
           ENDDO
        ELSE
           JB = 1
        ENDIF

C       REVERSE COLS ON RIGHT HALF OF REMAINING LINES INTO LEFT HALF 
        DO J=JB,NY
           JJ = NY-J+JB
           DO I=1,NSC
              II     = NNNN-I
              Q(I,J) = Q(II,JJ)
             
           ENDDO
        ENDDO

        !write(6,'(a,32e10.2)') ' 7:',q(1:17,7)
        !write(6,'(a,32e10.2)') ' 8:',q(1:17,8)

 
        !call chkmaxloc('y=qny/2m1', q(1:nx, ny/2-1),nx)      
        !call chkmaxloc('y=qny/2  ', q(1:nx, ny/2),  nx)      
        !call chkmaxloc('y=qny/2p1', q(1:nx, ny/2+1),nx)
      
        Q(NX/2+1,NY/2+1) = Q(NX/2,NY/2)

        IF (MODE == 'L') THEN
          CALL AAL11(Q,NNNN*NY)
        ENDIF


        END


         SUBROUTINE  AAL11(X,N)

         INTEGER :: N
         REAL    :: X(N)

         DO I=1,N
            IF (X(I) > 0.0) THEN
               X(I) = LOG10(X(I))
            ELSE
C              CAN NOT TAKE LOG!!
               X(I) = 0.0
               !write(6,'(a,i8,e12.4)') 'bad:',i,x(i)
            ENDIF
         ENDDO
         END



