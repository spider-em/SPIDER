
C++*********************************************************************
C
C    LISTITR.F                FILENAMES DEC 88 ARDEAN LEITH
C                             RENAMED FROM OPNPICP, SEP 96 ARDEAN LEITH
C                             USED REG_,            AUG 00 ARDEAN LEITH
C                             REG(',I5,')           NOV 05 ARDEAN LEITH
C                             GET_NAME              MAR 06 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2013  Health Research Inc.,                         *
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
C    LISTITR(FILNAM,LUN,NX,NY,NZ)
C
C    PURPOSE:   TO PLACE SPIECIFED IMAGE LOCATIONS IN SPECIFIED
C               REGISTERS
C
C    PARAMETERS:
C         FILNAM     FILE NAME 
C         LUN        LOGICAL UNIT NUMBER OF FILE
C         NX,NY,NZ   IMAGE DIMENSIONS
C
C    DESCRIPTION:  
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

	SUBROUTINE LISTITR(FILNAM,LUN,NX,NY,NZ)

        IMPLICIT NONE

	INCLUDE 'CMBLOCK.INC' 
	INCLUDE 'CMLIMIT.INC'
 
        REAL :: BUF(NBUFSIZ)
        COMMON /IOBUF/ BUF

	INTEGER, PARAMETER :: NBUF = 20
 	INTEGER            :: ISAM(NBUF),IROW(NBUF),ILIST(NBUF)
	REAL               :: VALUES(NBUF)
        CHARACTER *(*)     :: FILNAM
        CHARACTER (LEN=22) :: PROMPT
        CHARACTER (LEN=80) :: NAME
        CHARACTER (LEN=1)  :: NULL = CHAR(0)

        INTEGER            :: LUN,NX,NY,NZ,NUM,ITREC,K,NDUM,NLET,IRTFLG
        INTEGER            :: NREG,NGOT,LEN_NAME
        INTEGER            :: LNBLNKN

        CALL REG_GET_USED(NREG)
        
C                 123456789 123456789 1234567890
        PROMPT = 'X-COORDINATES (MAX 20)'

        NUM = MIN(NBUF,NREG)
        CALL RDPRAI(ISAM,NBUF,NUM, 1,NX,PROMPT,NULL,IRTFLG)

C                 123456789 123456789 1234567890
        PROMPT = 'Y-COORDINATES (MAX 20)'
        CALL RDPRAI(IROW,NBUF,NUM, 1,NY,PROMPT,NULL,IRTFLG)

        ITREC = NY * NZ

        DO K =1,NUM
            IF (IROW(K) > ITREC .OR. ISAM(K) > NX) THEN
C              POINT OUTSIDE OF IMAGE 
               CALL ERRT(101,'OUTSIDE IMAGE',NDUM)
               RETURN
            ENDIF

            CALL REDLIN(LUN,BUF,NX,IROW(K))

            VALUES(K) = BUF(ISAM(K))
        ENDDO

        CALL REG_SET_NSELA(NUM,VALUES,.FALSE.,IRTFLG)

        NLET = lnblnkn(FILNAM)
        WRITE(NDAT,91) FILNAM(1:NLET)
91      FORMAT('  FILE: ',A)

        CALL REG_GET_SELS(ILIST,NBUF,NGOT,IRTFLG)

        DO K = 1,NGOT
           CALL REG_GET_NAME(ILIST(K),NAME,LEN_NAME,IRTFLG)

           WRITE(NDAT,90)ISAM(K),IROW(K),NAME(1:LEN_NAME),VALUES(K)
 90        FORMAT('  X = ',I4,'  Y = ',I4,'  ',A,' = ',G12.4)
        ENDDO

 	END


