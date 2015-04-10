C++*********************************************************************
C
C  GETDOCLIST.F        CREATED                     JUL  03 ArDean Leith
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
C  GETDOCLIST(LUNDOC,PROMPT,ILIST,NANG,IREGO,IREGEND,DOCBUF,IRTFLG)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE GETDOCLIST(PROMPT,LUNDOC,ILIST,NANG,
     &                        IREGO,IREGEND,FLIP,DOCBUF,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'F90ALLOC.INC'

        REAL,    DIMENSION(:,:), POINTER      :: PANG
        INTEGER, DIMENSION(NANG)              :: ILIST
        REAL, DIMENSION(IREGEND-IREGO+1,NANG) :: DOCBUF  
        CHARACTER(LEN=*)                      :: PROMPT
        CHARACTER(LEN=MAXNAM)                 :: DOCNAM 
        LOGICAL                               :: FLIP

C       RETRIEVE ARRAY WITH DOC. FILE DATA IN IT
        MAXXT = 4
        MAXYT = 0                     

C       ALLOCATE PANG IN GETDOCDAT TO RETRIEVE DATA FROM DOC FILE
        CALL GETDOCDAT(PROMPT,.TRUE.,DOCNAM,LUNDOC,
     &                 .TRUE.,MAXXT,MAXYT,PANG,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        DO K=1,NANG
           ITMP    = ILIST(K)  
           ICOUNT = PANG(1,ITMP)
           IF (ICOUNT .LE. 0) THEN
C             MISSING THIS KEY
              IF (FLIP) THEN
                 CALL ERRT(102,'MISSING ANGLES FOR IMAGE',ITMP)
              ELSE
                 CALL ERRT(102,'MISSING DOC. FILE KEY:',ITMP)
              ENDIF
              GOTO 9999
           ENDIF

           IF (FLIP) THEN
C             ORDER IN THE DOCUMENT FILE IS PSI, THETA, PHI  
C             IN ANG ARRAY ORDER SHOULD BE OTHER WAY AROUND
              IT = IREGEND
              DO IREG=IREGO,IREGEND
                 DOCBUF(IT,K) = PANG(IREG+1,ITMP)
                 IT           = IT - 1
              ENDDO
           ELSE
              IT = 1
              DO IREG=IREGO,IREGEND
                 DOCBUF(IT,K) = PANG(IREG+1,ITMP)
                 IT           = IT + 1
              ENDDO
           ENDIF
        ENDDO
        IRTFLG = 0

9999    IF (ASSOCIATED(PANG)) DEALLOCATE(PANG)
        NULLIFY(PANG)

        END

