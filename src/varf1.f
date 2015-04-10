C **********************************************************************
C VARF1.F
C              OPFILEC                            FEB  03 ARDEAN LEITH
C              MAXNAM                             JUL  14 ARDEAN LEITH
C
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
C                                      
C VARF1(INC,Y1,WI,ILIST,NANG,NNNN,NSAM,NROW,RMSK)
C       
C23456789012345678901234567890123456789012345678901234567890123456789012
C **********************************************************************

        SUBROUTINE VARF1(INC,Y1,WI,ILIST,NANG,NNNN,NSAM,NROW,RMSK)
     
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        PARAMETER (NDLI=4)

        !CHARACTER*80           FINPAT1,FINPIC1,FINPAT2,FINPIC2
        CHARACTER(LEN=MAXNAM) ::FINPAT1,FINPIC1,FINPAT2,FINPIC2
        COMMON  /F_SPEC/  FINPAT1,FINPIC1,NLET1,FINPAT2,FINPIC2,NLET2

C       Mask may not be same as nlet1 in varf that calls this ?? al

        REAL, ALLOCATABLE, DIMENSION(:,:) :: B1,B2
        DIMENSION DLIST(NDLI),ILIST(NANG),RMSK(NSAM,NROW)

C       LR AND FR ARE AUTOMATIC ARRAYS

        DOUBLE PRECISION FR(INC),LR(INC)
 
        DATA  LUNI/97/

        ALLOCATE (B1(NNNN,NROW), B2(NNNN,NROW), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'VA F; B1 & B2',2*NNNN*NROW)
           RETURN
        ENDIF

C       INITIALIZE THE SUMS
        FR=0.0D0
        LR=0

        DO  K=1,NANG
C          READ ONE PROJECTION
           CALL  FILGET(FINPAT1,FINPIC1,NLET1,ILIST(K),INTFLAG)
           CALL OPFILEC(0,.FALSE.,FINPIC1,LUNI,'O',ITYPE,LSAM,LROW,NSL,
     &             MAXIM,'DUMMY',.FALSE.,IRTFLG1)
           IF (IRTFLG1 .NE. 0) GOTO 2001

           CALL READV(LUNI,B1,NNNN,NROW,NSAM,NROW,1)
           CLOSE(LUNI)

           CALL  FILGET(FINPAT2,FINPIC2,NLET2,ILIST(K),INTFLAG)
           CALL OPFILEC(0,.FALSE.,FINPIC2,LUNI,'O',ITYPE,LSAM,LROW,NSL,
     &             MAXIM,'DUMMY',.FALSE.,IRTFLG2)
           IF (IRTFLG2 .NE. 0)  GOTO 2001

           CALL READV(LUNI,B2,NNNN,NROW,NSAM,NROW,1)
           CLOSE(LUNI)

cc$omp   parallel sections
cc$omp   section
                CALL  SUBACOR(B1,NNNN,NROW)
cc$omp   section
                CALL  SUBACOR(B2,NNNN,NROW)
cc$omp   end parallel sections


c$omp parallel do private(i,j)
              DO    J=1,NROW
                 DO    I=1,NSAM
                    B1(I,J)=B1(I,J)*RMSK(I,J)
                    B2(I,J)=B2(I,J)*RMSK(I,J)
                 ENDDO
              ENDDO

              INV=+1
              CALL  FMRS_2(B1,NSAM,NROW,INV)
              IF(INV.NE.1)  THEN
                 CALL  ERRT(38,'VA F',NE)
                 RETURN
              ENDIF

              CALL  FMRS_2(B2,NSAM,NROW,INV)
              IF(INV.NE.1)  THEN
                 CALL  ERRT(38,'VA F',NE)
                 RETURN
              ENDIF


         DO J=1,NROW
           JJ=J-1
           IF(JJ.GT.NROW/2)  JJ=JJ-NROW
           DO II=1,NSAM,2
            I=II/2+1
            IF(.NOT.(I.EQ.1.AND.JJ.LE.0)) THEN
              PII=0.5*SQRT((FLOAT(JJ)/FLOAT(NROW/2))**2+
     &        (FLOAT(I-1)/FLOAT(NSAM/2))**2)
              IF(PII.LE.0.5)  THEN
                 L=MIN0(MAX0(NINT(PII*Y1/WI)+1,1),INC)
                 LR(L)=LR(L)+1
                 FR(L)=FR(L)+
     &            (B1(I,J)-B2(I,J))**2+(B1(I+1,J)-B2(I+1,J))**2
              ENDIF
            ENDIF
           ENDDO
         ENDDO
        ENDDO

C SAVE RESULTS
        PII=1./REAL(NSAM)/REAL(NROW)/REAL(NSAM)/REAL(NROW)
        DO   L=1,INC
           IF(LR(L).GT.0)  THEN
              DLIST(1)=L
              DLIST(2)=FLOAT(L-1)/FLOAT(INC-1)*0.5
              DLIST(3)=FR(L)*PII/LR(L)
              DLIST(4)=LR(L)
              CALL  SAVD(LUNI,DLIST,NDLI,IRTFLG)
           ENDIF
        ENDDO
        CALL  SAVDC
        CLOSE(LUNI)
C----------------------------------------------
2001    DEALLOCATE (B1)
        DEALLOCATE (B2)
        END 



        SUBROUTINE SUBACOR(PROJ,LSD,N)

        DIMENSION  PROJ(LSD,N)
        DOUBLE     PRECISION QS

        KLP=0
        R=(N/2)**2
        QS=0.0D0

C       ESTIMATE AVERAGE OUTSIDE THE CIRCLE
        NC =N/2+1
c$omp parallel do private(i,j,t,xx),reduction(+:qs,klp)
        DO   J=1,N
           T=J-NC
           XX=T*T
           DO   I=1,N
              T=I-NC
              IF (XX+T*T.GT.R)    THEN
                 QS=QS+DBLE(PROJ(I,J))
                 KLP=KLP+1
              ENDIF
           ENDDO
        ENDDO
        QS = QS/REAL(KLP)
c$omp parallel do private(i,j)
        DO  J=1,N
           DO  I=1,N
              PROJ(I,J)=PROJ(I,J)-QS
           ENDDO
        ENDDO

        END
