C++*********************************************************************
C
C PW3SR.F                       HALF BUG FIXED   FEB 02  ArDean Leith
C                               PARALLEL         FEB 18  ArDean Leith
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2018  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
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
C PW3SR(Q,NX,NY,NZ,MODE)
C
C PURPOSE: POWER SPECTRUM OF VOLUME
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE PW3SR(Q,NX,NY,NZ,MODE)

        IMPLICIT NONE
        INTEGER      :: NX,NY,NZ !GFORT ERROR IF NOT BEFOR Q,  Sept 2025 al
        REAL         :: Q(NX+2-MOD(NX,2),NY,NZ)!GFORT ERROR,  Sept 2025 al
        CHARACTER    :: MODE

        LOGICAL      :: EVENY,EVENZ,WANTLOG,WANT2,WANTPH
        INTEGER      :: LSD,LSDM1,NSC,K,J,I,JJ,KK,NSL,KB,II,JB
        INTEGER      :: NXD2,NYD2,NZD2
        REAL         :: SCL,SCLSQ,TEMP,FVAL


        LSD   = NX + 2 - MOD(NX,2)
        LSDM1 = LSD -1

        EVENY = MOD(NY,2) == 0
        EVENZ = MOD(NZ,2) == 0

        NSC   = NXD2+1               ! REDEFINED BELOW

        NXD2  = NX/2
        NYD2  = NY/2
        NZD2  = NX/2

        WANTPH  = (INDEX(MODE,'P') > 0)
        WANTLOG = (INDEX(MODE,'L') > 0)
        WANT2   = (INDEX(MODE,'2') > 0)

C       SCALE FACTOR
        SCL   = 2.0 / FLOAT(NX) / FLOAT(NY) / FLOAT(NZ)
        SCLSQ = SCL * SCL

C       LOOP OVER ALL ROWS, PUTTING AMPS IN RIGHT HALF OF ROW

	IF (WANT2)  THEN    ! UNDOCUMENTED
C         SQUARE OF STANDARD POWER SPECTRUM
c$omp     parallel do private(k,j,i)
          DO K=1,NZ
            DO J=1,NY
              DO I=LSDM1,1,-2
                 Q(LSDM1-(LSDM1-I)/2,J,K) = SCLSQ *
     &               (Q(I+1,J,K)*Q(I+1,J,K) + Q(I,J,K)*Q(I,J,K))
              ENDDO
            ENDDO
          ENDDO

        ELSEIF (WANTLOG) THEN
C         LOG OF STANDARD POWER SPECTRUM
c$omp     parallel do private(k,j,i,fval)
          DO K=1,NZ
            DO J=1,NY
              DO I=LSDM1,1,-2
                 FVAL = SCL * 
     &                  SQRT(Q(I+1,J,K)*Q(I+1,J,K)+Q(I,J,K)*Q(I,J,K))

                 IF (FVAL > 0.0) THEN
                    Q(LSDM1-(LSDM1-I)/2,J,K) = LOG10(FVAL)
                 ELSE
C                   CAN NOT TAKE LOG!!
                    Q(LSDM1-(LSDM1-I)/2,J,K) = 0.0
                 ENDIF
              ENDDO
            ENDDO
          ENDDO

        ELSE
C         STANDARD POWER SPECTRUM
c$omp     parallel do private(k,j,i)
          DO K=1,NZ
             DO J=1,NY
               DO I=LSDM1,1,-2
                  Q(LSDM1-(LSDM1-I)/2,J,K) = SCL *
     &                SQRT(Q(I+1,J,K)*Q(I+1,J,K)+Q(I,J,K)*Q(I,J,K))
               ENDDO
             ENDDO
          ENDDO
        ENDIF

       
C       INTERCHANGE TOP AND BOTTEM HALF OF RIGHT HALF COLUMNS
          
c$omp   parallel do private(k,j,jj,i,temp)
        DO K=1,NZ
           DO J=1,NYD2
              JJ = J + NYD2 + MOD(NY,2)
              DO I=NSC,LSDM1
                 TEMP      = Q(I,J,K)
                 Q(I,J,K)  = Q(I,JJ,K)
                 Q(I,JJ,K) = TEMP
              ENDDO
           ENDDO
        ENDDO

        IF (.NOT. EVENY)  THEN
           DO K=1,NZ
              DO I=NSC,LSDM1
                 TEMP = Q(I,NYD2+1,K)
                 DO J=NYD2+1,NY-1
                    Q(I,J,K) = Q(I,J+1,K)
                 ENDDO
                 Q(I,NY,K) = TEMP
              ENDDO
           ENDDO
        ENDIF   

c$omp   parallel do private(k,kk,j,i,temp)
        DO K=1,NZD2
           KK = K + NZD2 + MOD(NZ,2)
           DO J=1,NY
              DO I=NSC,LSDM1
                 TEMP      = Q(I,J,K)
                 Q(I,J,K)  = Q(I,J,KK)
                 Q(I,J,KK) = TEMP
              ENDDO
           ENDDO
        ENDDO

        IF (.NOT. EVENZ)  THEN
           DO J=1,NY
              DO I=NSC,LSDM1
                 TEMP = Q(I,J,NZD2+1)
                 DO K=NZD2+1,NZ-1
                    Q(I,J,K) = Q(I,J,K+1)
                 ENDDO
                 Q(I,J,NZ) = TEMP
              ENDDO
           ENDDO
        ENDIF   

C       REVERSE COLS ON RIGHT HALF OF FIRST LINE INTO LEFT HALF
        NSC = LSD / 2 - 1
        NSL = NZD2
        IF (EVENY)  THEN
           IF (EVENZ)  THEN
              KB = 2
              DO I=1,NSC
                 II       = LSD - I
                 Q(I,1,1) = Q(II,1,1)
              ENDDO
           ELSE
              KB = 1
           ENDIF

           JB = 2
           DO K=KB,NZ
              KK = NZ-K+KB
              DO I=1,NSC
                 II       = LSD - I
                 Q(I,1,K) = Q(II,1,KK)
              ENDDO
           ENDDO
        ELSE
           JB=1
        ENDIF

        IF (EVENZ)  THEN
           KB = 2
           DO J=JB,NY
              JJ = NY-J+JB
              DO I=1,NSC
                 II       = LSD - I
                 Q(I,J,1) = Q(II,JJ,1)     
              ENDDO
           ENDDO
        ELSE
           KB=1
        ENDIF

C       REVERSE COLS ON RIGHT HALF OF REMAINING LINES INTO LEFT HALF 

c$omp parallel do private(k,kk,j,jj,i,ii)
        DO K=KB,NZ
           KK = NZ-K+KB
           DO J=JB,NY
              JJ=NY-J+JB
              DO I=1,NSC
                 II       = LSD - I
                 Q(I,J,K) = Q(II,JJ,KK)
              ENDDO
           ENDDO
        ENDDO

C       REPLACE ORIGIN WITH ADJACENT VALUE TO HELP CONTRAST
	Q(NXD2+1,NYD2+1,NZD2+1) = Q(NXD2,NYD2,NZD2)

        END



#ifdef NEVER
         subroutine  checkit(x,nx,lsd,ny,nz)

         real    :: x(lsd,ny,nz)
         integer :: nx,ny,nz

         !write(6,*)' log(1),log(0):',n,log10(1.0),log10(0.0),log10(-1.0)

         do   k=1,nz
         do   j=1,ny
         do   i=1,nx
            if (x(i,j,k) < 0) then
               write (6,*) 'i,j,k,q(ijk):',i,j,k,x(i,j,k)
               stop
            endif
         enddo
         enddo
         enddo
         end

#endif
