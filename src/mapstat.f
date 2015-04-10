
C***********************************************************************
C
C  MAPSTAT.FOR  -- CREATED JAN 91
C
C **********************************************************************
C *  AUTHOR: ArDean Leith 
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
C      MAPSTAT(LUNIN,LUND,NX,NY,NZ,NREC1,NREC2,NVOX,IRTFLG) 
C
C      PURPOSE:  EXAMINES AN IMAGE STACK FOR NUMBERS OF EACH IMAGE VALUE
C                AND CENTER OF MASS OF EACH IMAGE VALUE.  CREATES A DOC.
C                FILE FOR STATISTICS OUTPUT HAVING KEY=IMAGE VALUE + 1, 
C                THEN NUMBER OF CELLS HAVING THAT VALUE, IXCOM, IYCOM, & 
C                IZCOM.  LIMITED TO IMAGEVALUES >= 0 AND <= MAXVOX.
C
C      PARAMETERS:  
C
C      CALLED BY:  IMSTAT 
C
C      CALLS:      REDLIN      WRTLIN   
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

       SUBROUTINE MAPSTAT(LUNIN,LUND,NX,NY,NZ,NREC1,NREC2,
     &            NVOX,IRTFLG)

 
       INCLUDE 'CMLIMIT.INC'
       INCLUDE 'CMBLOCK.INC'


       PARAMETER (MAXBUF = 16000)
       PARAMETER (MAXVOX = 9999)

       LOGICAL :: ONZ
       COMMON BUF(MAXBUF),FVOX(MAXVOX),FXCEN(MAXVOX),FYCEN(MAXVOX),
     &        FZCEN(MAXVOX),ONZ(MAXVOX),NUMSL(MAXVOX)

       REAL                   :: DLIST(7)
       LOGICAL                :: DEBUGGING
       CHARACTER (LEN=MAXNAM) :: DOCNAM
       LOGICAL                :: ADDEXT,GETNAME,ISOLD
       LOGICAL                :: APPEND,MESSAGE,NEWFILE
       CHARACTER(LEN=90)      :: COMMENT                   


       DEBUGGING = .FALSE.

       IOVER = 0
       NVOX  = 0

C      ZERO STATISTICS COUNTERS

       DO I = 1,MAXVOX
          FVOX(I)  = 0.0
          FXCEN(I) = 0.0
          FYCEN(I) = 0.0
          FZCEN(I) = 0.0
          NUMSL(I) = 0
       END DO

       LASTSLICE = -100        
       DO IREC = NREC1,NREC2
          CALL REDLIN(LUNIN,BUF,NX,IREC)
          ISLICE = (IREC / NY) + 1
          IROW   = IREC - (ISLICE - 1) * NY 

          IF (ISLICE .NE. LASTSLICE) THEN
C            STARTING A NEW SLICE, CHANGE ONZ FLAGS
             DO IBLOB = 1,MAXVOX
                ONZ(IBLOB) = .FALSE.
             ENDDO
             LASTSLICE = ISLICE
          ENDIF

          DO ICOL = 1,NX
C            CHECK THIS VOXEL STATISTICS, NOTE THAT TABLE LOCATION IS
C            1 > THAN IMAGE VALUE!!

             IT = BUF(ICOL) + 1

             IF (IT .LE. 0 .OR. IT > MAXVOX) THEN
C               VOXEL VALUE TOO SMALL OR LARGE FOR STATISTICS TABLE
                IOVER = IOVER + 1
                IF (IOVER < 8) THEN
C                   REPORT THE PROBLEM 8 TIMES
                    WRITE(NOUT,*) ' *** VOXEL VALUE: ',IT,
     &                             ' OUTSIDE TABLE SIZE: ',MAXVOX
                ENDIF

             ELSE
C               VOXEL VALUE IS WITHIN TABLE LIMITS
                FVOX(IT)  = FVOX(IT)  + 1
                FXCEN(IT) = FXCEN(IT) + ICOL
                FYCEN(IT) = FYCEN(IT) + IROW
                FZCEN(IT) = FZCEN(IT) + ISLICE

                IF (.NOT. ONZ(IT)) THEN
                   ONZ(IT)   = .TRUE.
                   NUMSL(IT) = NUMSL(IT) + 1
                ENDIF    
                IF (IT > NVOX) NVOX = IT
             ENDIF
          END DO

       END DO

C       OPEN OUTPUT DOC FILE 
        ADDEXT  = .TRUE.                    
        GETNAME = .TRUE.                    
        ISOLD   = .FALSE.                   
        APPEND  = .FALSE.                    
        MESSAGE = .TRUE.                    
        IRTFLG  = -8         ! NO IC USE    

        LUNRET  = LUND
        CALL OPENDOC(DOCNAM,ADDEXT,NLET,LUND,LUNRET,GETNAME,
     &             'OUTPUT STATISTICS DOC',ISOLD,APPEND,MESSAGE,
     &             NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (NEWFILE) THEN
C                     123456789 123456789 123456789 123456789 123456789 1
           COMMENT = 'CONTENTS:   CLUSTER STATISTICS'

           CALL LUNDOCPUTCOM(LUNRET,COMMENT(1:37),IRTFLG)
C                   123456789 123456789 123456789 123456789 123456789
           COMMENT='       #VOXELS        X-COM          Y-COM      '//
     &             '   Z-COM  #OCCUPIED-SLICES '
C                   9 123456789 123456789 123456789 123456789 123456789
           CALL LUNDOCPUTCOM(LUNRET,COMMENT(1:75),IRTFLG)
       ENDIF


C      CALCULATE STATISTICS FOR EACH VOXEL VALUE
       DO I = 1,NVOX
          INUM     = FVOX(I)

          IF (INUM > 0) THEN
C           VOXEL VALUE IS OCCUPIED
            KEY      = I
            DLIST(1) = FVOX(I)
            DLIST(2) = FXCEN(I) / FVOX(I)
            DLIST(3) = FYCEN(I) / FVOX(I)
            DLIST(4) = FZCEN(I) / FVOX(I)
            DLIST(5) = NUMSL(I) 
            IF (KEY == 1) DLIST(5) = DLIST(5) - 1  ! why?? al

C           SAVE STATISTICS IN DOC FILE
            CALL LUNDOCWRTDAT(LUNRET,KEY,DLIST,5,IRTFLG)

            IF (DEBUGGING) THEN
C              WRITE STATISTICS ON STANDARD OUTPUT
               WRITE(NOUT,90) I,INUM,(DLIST(J),J=2,5)
90             FORMAT(2I10,4G12.3)
            ENDIF
          ENDIF
       ENDDO

C      CLOSE THE DOC FILE
9999   CLOSE(LUND)

       END
