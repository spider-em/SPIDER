C **********************************************************************
C VARF.F
C              OPFILEC                  FEB  03 ARDEAN LEITH
C              MAXNAM                   JUL  14 ARDEAN LEITH
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
C VARF.F
C 
C PURPOSE:   Calculate the variance in Fourier space
C
C **********************************************************************

        SUBROUTINE  VARF

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        !CHARACTER*80      FINPAT1,FINPIC1,MASK,FINPAT2,FINPIC2
        CHARACTER(LEN=MAXNAM) ::  FINPAT1,FINPIC1,MASK,FINPAT2,FINPIC2
        COMMON  /F_SPEC/  FINPAT1,FINPIC1,NLET1,FINPAT2,FINPIC2,NLET2
 
        REAL, ALLOCATABLE, DIMENSION(:,:)  :: RMSK
        INTEGER, ALLOCATABLE, DIMENSION(:) :: ILIST
 
        DATA  INPIC/98/

        NILMAX = NIMAX
        ALLOCATE (ILIST(NILMAX), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'VA F, ILIST',NILMAX)
           RETURN
        ENDIF

C       -------READ FILE LISTS -------
        CALL FILELIST(.TRUE.,INPIC,FINPAT1,NLET1,ILIST,NILMAX,NANG,
     &                 'FIRST TEMPLATE FOR 2-D IMAGES',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CLOSE(INPIC)

        NIL = 0
        CALL FILELIST(.TRUE.,INPIC,FINPAT2,NLET2,ILIST,NIL,IDUM,
     &                 'SECOND TEMPLATE FOR 2-D IMAGES',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       NANG - TOTAL NUMBER OF IMAGES

C       OPEN THE REAL-SPACE MASK FILE, ALSO YIELDS NSAM AND NROW
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,MASK,INPIC,'O',IFORM,NSAM,NROW,NSLICE,
     &               MAXIM,'MASK',.FALSE.,IRTFLG)

        ALLOCATE (RMSK(NSAM,NROW), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'VA F, RMSK',IER)
           RETURN
        ENDIF

        CALL READV(INPIC,RMSK,NSAM,NROW,NSAM,NROW,1)
        CLOSE(INPIC)

        CALL RDPRM(WI,IDUM,'RING WIDTH')

        Y1     = FLOAT(MAX0(NSAM,NROW))
        INC    = INT(Y1/WI)/2+1
        NNNN   = NSAM+2-MOD(NSAM,2)

        CALL VARF1(INC,Y1,WI,ILIST,NANG,NNNN,NSAM,NROW,RMSK)

        DEALLOCATE (RMSK) 

        END
