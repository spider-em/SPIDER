C ++********************************************************************
C                                                                      *
C  QFACT_P                                                         
C                                                                      *
C                  OPFILEC                         FEB 03 ARDEAN LEITH
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
C                                                                      
C  QFACT_P                                                         
C                                                                      
C  PURPOSE:                                                            
C                                                                      
C  PARAMETERS:                                                         
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE QFACT_P(LUN1,LUN2,LUNQ,INUM,NUMT,NSAM,NROW,HREF,KREF)

        INCLUDE 'CMLIMIT.INC'

        PARAMETER  (NLIST=2)
        DIMENSION  DLIST(NLIST)

        CHARACTER (LEN=MAXNAM) :: FILNAM,FILA,FILQ,FILPAT,DOCFIL

        COMPLEX                              :: QSUM
        DIMENSION                            :: INUM(NUMT)
        INTEGER                              :: HREF,KREF
        LOGICAL                              :: SKIP
        COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: FA,CBUF
        REAL, ALLOCATABLE, DIMENSION(:,:)    :: FABS
        DIMENSION  BUF(NSAM)

        PARAMETER (QMAX=1.0E5)

        DATA   NDOUT/91/

        NROWH=NROW/2
        NSAMH=NSAM/2
        IF ( HREF.EQ.0 . AND. KREF.EQ.0)  THEN
           SKIP=.TRUE.
        ELSE
           IF (HREF.GE.0)  THEN
              HREF=HREF+1
           ELSE
              HREF=-HREF+1
              KREF=-KREF
           ENDIF
           IF (KREF.GE.0)  THEN
              KREF=KREF+1
           ELSE
              KREF=NROW+KREF+1
           ENDIF
           SKIP=.FALSE.
           QSUM=(0.0,0.0)
           SMOD=0.0
        ENDIF

        ALLOCATE (FA(NSAM/2+1,NROW), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'AS F, FA',IER)
           RETURN
        ENDIF

        ALLOCATE (FABS(NSAM/2+1,NROW), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'AS F, FABS',IER)
           DEALLOCATE (FA)
           RETURN
        ENDIF

C       CLEAR BUFFER ARRAYS
        FABS=0.
        FA=(0.,0.)

        ALLOCATE (CBUF(NSAM/2+1,NROW), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'AS F, CBUF',IER)
           DEALLOCATE (FABS)
           DEALLOCATE (FA)
           RETURN
        ENDIF

C       LOOP OVER ALL FILES
        DO  IFIL=1,NUMT
           CALL FILGET(FILPAT,FILNAM,NLET,INUM(IFIL),IRTFLG)
           IF(IRTFLG.NE.0) THEN
              DEALLOCATE (FABS)
              DEALLOCATE (FA)
              RETURN
           ENDIF

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FILNAM,LUN2,'O',ITYPE,
     &                 NSAM,NROW,NSLICE,MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              DEALLOCATE (FABS)
              DEALLOCATE (FA)
              RETURN
           ENDIF
         
           NSLICE =1           
           CALL READV(LUN2,CBUF,NSAM+2,NROW,NSAM,NROW,NSLICE)
           CLOSE(LUN2)
C
           INV=+1
           CALL  FMRS_2(CBUF,NSAM,NROW,INV)
C
c$omp      parallel do private(i,j)
           DO    I=1,NSAMH+1
              DO    J=1,NROW
                 FABS(I,J)=FABS(I,J)+CABS(CBUF(I,J))
                 FA(I,J)=FA(I,J)+CBUF(I,J)
              ENDDO
           ENDDO
C
           IF(.NOT.SKIP)  THEN
              QSUM=QSUM+FA(HREF,KREF)
              SMOD=SMOD+CABS(FA(HREF,KREF))
              DLIST(1)=IFIL
              DLIST(2)=CABS(QSUM)/SMOD
              CALL  SAVD(NDOUT,DLIST,NLIST,IRTFLG)
           ENDIF
C       end of the loop over the images
        ENDDO
        IF(.NOT.SKIP)  THEN
           CALL  SAVDC
           CLOSE(NDOUT)
        ENDIF
C       GENERATE FRIEDEL-RELATED PART OF Q-MAP AND WRITE OUT
        DO  J=1,NROW
           IF(J.LE.NROWH)  THEN
              KJ=J+NROWH
           ELSE
              KJ=J-NROWH
           ENDIF
           IF(J.EQ.1.OR.J.EQ.NROWH+1)  THEN
              MJ=J
           ELSE
              MJ=NROW-J+2
           ENDIF

c$omp      parallel do private(i)
           DO  I=1,NSAMH
              IF(QMAX*FABS(I,J).GT.CABS(FA(I,J)))  THEN
                 BUF(I+NSAMH)=CABS(FA(I,J))/FABS(I,J)
                 BUF(NSAMH-I+1)=CABS(FA(I+1,MJ))/FABS(I+1,MJ)
              ELSE
                 BUF(I+NSAMH)=QMAX
                 BUF(NSAMH-I+1)=QMAX
              ENDIF
           ENDDO
           CALL  WRTLIN(LUNQ,BUF,NSAM,KJ)
        ENDDO
        CLOSE(LUNQ)
C
c$omp   parallel do private(i,j)
        DO    I=1,NSAMH+1
           DO    J=1,NROW
              CBUF(I,J)=FA(I,J)/NUMT
           ENDDO
        ENDDO
        INV=-1
        CALL  FMRS_2(CBUF,NSAM,NROW,INV)

        NSLICE =1
        CALL WRITEV(LUN1,CBUF,NSAM+2,NROW,NSAM,NROW,NSLICE)

        CLOSE (LUN1)

        DEALLOCATE (CBUF)
        DEALLOCATE (FABS)
        DEALLOCATE (FA)

        END



