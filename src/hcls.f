
C++*********************************************************************
C
C  HCLS.F
C               USED ALLOCATE                    JAN 2001 ARDEAN LEITH
C               Q** ALLOC                        MAY 2013 ARDEAN LEITH
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
C  HCLS(LUNF,LUNT,LUND)
C
C  CALL TREE:  HCLS -> HCLP -> DIST_P
C                              CHAVA
C                              DENDRO
C                              DENDRO2
C                              ARBRE
C                              DENLST 
C
C--*********************************************************************

         SUBROUTINE HCLS(LUNF,LUNT,LUND)

         IMPLICIT NONE

         INCLUDE 'CMBLOCK.INC' 
         INCLUDE 'CMLIMIT.INC' 

         INTEGER               :: LUNF,LUNT,LUND

         REAL, ALLOCATABLE     :: Q_K_W(:),    Q_MEMD( :), Q_MEMPK(:)
         REAL, ALLOCATABLE     :: Q_MEMVAL(:), Q_MEMLA(:), Q_MEMLB(:)     &            Q_MEMNT( :),
         REAL, ALLOCATABLE     :: Q_MEMNT(:),  Q_MEMNO(:), Q_MEMIV(:)
         REAL, ALLOCATABLE     :: Q_MEMIW(:),  Q_MEMV(:),  Q_MEMW(:) 

         CHARACTER(LEN=MAXNAM) :: FILNAM
         CHARACTER(LEN=MAXNAM) :: FILNAMT
         CHARACTER(LEN=1)      :: NULL = CHAR(0)

         INTEGER *8            :: MEMTOT_8,MDIM_8

         INTEGER               :: NLET,IRTFLG,IT,ITYPE,NKLA,NFAC
         INTEGER               :: MAXFAC,MAXFA,MINFAC,NE,IER,I
         INTEGER               :: JFIN,NOT_USED,MODE
         REAL                  :: W1

         WRITE(NOUT,*) ' YOU MAY USE A _IMC, _PIX, or _SEQ FILE'
         WRITE(NOUT,*) ' '

         CALL FILERD(FILNAM,NLET,NULL,
     &              'CORAN/PCA FILE (e.g. CORAN_01_IMC)~',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         IT = INDEX(FILNAM,'_SEQ')
         IF (IT > 0) THEN
            ITYPE = 1
            WRITE(NOUT,*)' For _SEQ coordinates file ---'
            CALL OPAUXFILE(.FALSE.,FILNAM,DATEXC,-LUNF,0,
     &                 'O',' ',.TRUE.,IRTFLG)
            IF (IRTFLG .NE. 0) RETURN

            READ(LUNF)  NKLA,NFAC

         ELSEIF(INDEX(FILNAM,'_IMC') > 0) THEN
            CALL OPAUXFILE(.FALSE.,FILNAM,DATEXC,LUNF,0,
     &                 'O',' ',.TRUE.,IRTFLG)
            IF (IRTFLG .NE. 0) RETURN
            ITYPE = 2
            WRITE(NOUT,*)' In image coordinates file ---'
            READ(LUNF,*)  NKLA,NFAC

         ELSEIF(INDEX(FILNAM,'_PIX') > 0) THEN
            CALL OPAUXFILE(.FALSE.,FILNAM,DATEXC,LUNF,0,
     &                 'O',' ',.TRUE.,IRTFLG)
            IF (IRTFLG .NE. 0) RETURN
            ITYPE = 3
            WRITE(NOUT,*)' In pixel coordinates file ---'
            READ(LUNF,*)  NKLA,NFAC

         ELSE
            CALL ERRT(101,'INVALID INPUT FILE TYPE',NE)
            RETURN
         ENDIF

         WRITE(NOUT,*)' Number of objects: ', NKLA
         WRITE(NOUT,*)' Number of factors: ', NFAC

         MAXFAC = NFAC
         MAXFA  = NFAC
         MINFAC = 1
         CALL  RDPRAI(INUMBR,NFAC,MAXFA, MINFAC,MAXFAC,
     &            'FACTOR NUMBERS',NULL,IER)
     
C        DIMENSION  W(MAXFA),COO(NFAC),COB(NFAC),INUM(MAXFA)

         
C        CALCULATE NUMBER OF FACTORS
         MDIM_8    = NKLA                    ! FOR BIG NKLA
         MDIM_8    = MDIM_8 * (NKLA-1) / 2   ! FOR BIG NKLA

         JFIN     = 2 * NKLA - 1

         MEMTOT_8 = NFAC + MDIM_8 + 2*JFIN + 2*NKLA + 2*JFIN +
     &              4*NKLA

         WRITE(NOUT,*) ' DYNAMIC MEMORY ALLOCATION: ',MEMTOT_8

         ALLOCATE(Q_K_W(   NFAC),
     &            Q_MEMD(  MDIM_8),
     &            Q_MEMPK( JFIN),
     &            Q_MEMVAL(JFIN),
     &            Q_MEMLA( NKLA),
     &            Q_MEMLB( NKLA),
     &            Q_MEMNT( JFIN),
     &            Q_MEMNO( JFIN),
     &            Q_MEMIV( NKLA),
     &            Q_MEMIW( NKLA),
     &            Q_MEMV(  NKLA),
     &            Q_MEMW(  NKLA),
     &            STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            WRITE(NOUT,*) '*** UNABLE TO ALLOCATE MEMORY:',MEMTOT_8
            CALL ERRT(46,'HCLS; Q...',MEMTOT_8)
            RETURN
         ENDIF
 
         DO I=1,MAXFA
            Q_K_W(I) = 1.0
         ENDDO

         DO I=1,MAXFA
            CALL  RDPRM(W1,NOT_USED,'FACTOR WEIGHT')
            IF (W1 == 0.0) EXIT
            Q_K_W(I) = W1
         ENDDO

201      WRITE(NOUT,*) ' FACTOR WEIGHTS USED:'
         WRITE(NOUT,23)  (Q_K_W(I),I=1,MAXFA)
23       FORMAT(10(2X,G10.3))

         CALL RDPRI1S(MODE,NOT_USED,
     &                 'CLUSTERING CRITERION (1-5)',IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999

         CALL  HCLP(NKLA,MDIM_8,JFIN,
     &         Q_K_W, INUMBR,MAXFA,NFAC,
     &         Q_MEMD, Q_MEMPK, Q_MEMVAL,
     &         Q_MEMLA,Q_MEMLB, Q_MEMNT,
     &         Q_MEMNO,Q_MEMIV, Q_MEMIW, Q_MEMV,
     &         Q_MEMW, LUNF,LUND,LUNT,MODE,ITYPE)
 

9999     IF (ALLOCATED(Q_K_W))   DEALLOCATE(Q_K_W)
         IF (ALLOCATED(Q_MEMD))  DEALLOCATE(Q_MEMD)
         IF (ALLOCATED(Q_MEMPK)) DEALLOCATE(Q_MEMPK)
         IF (ALLOCATED(Q_MEMVAL))DEALLOCATE(Q_MEMVAL)
         IF (ALLOCATED(Q_MEMLA)) DEALLOCATE(Q_MEMLA)
         IF (ALLOCATED(Q_MEMLB)) DEALLOCATE(Q_MEMLB)
         IF (ALLOCATED(Q_MEMNT)) DEALLOCATE(Q_MEMNT)
         IF (ALLOCATED(Q_MEMNO)) DEALLOCATE(Q_MEMNO)
         IF (ALLOCATED(Q_MEMIV)) DEALLOCATE(Q_MEMIV)
         IF (ALLOCATED(Q_MEMIW)) DEALLOCATE(Q_MEMIW)
         IF (ALLOCATED(Q_MEMV))  DEALLOCATE(Q_MEMV)
         IF (ALLOCATED(Q_MEMW))  DEALLOCATE(Q_MEMW)

         WRITE (NDAT,2600)
2600     FORMAT (/' ',10('-'),'END OF HIERARCHICAL CLUSTERING',10('-')/)

         END
