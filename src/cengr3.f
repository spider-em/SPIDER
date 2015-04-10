
C ++********************************************************************
C                                                                      *
C CENGR3.F                                      6/01/83  Joachim Frank *
C           LONG FILENAMES                       JAN 89  ArDean Leith  *
C           ERROR TRAPS                          FEB 12  ArDean Leith  *
C                                                                      *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C                                                                      *
C  CENGR3(LUN1)                                                        *
C                                                                      *
C  PURPOSE:                                                            *
C    CALCULATES CENTER OF GRAVITY AND RADIUS OF GYRATION OF A DENSITY  *
C    MASS STORED IN A THREE-D FILE. FOR REFERENCE, SEE INTERNATIONAL   *
C    TABLES OF CRYSTALLOGRAPHY, VOL.III P. 327.     --- JF 6/01/83     *
C                                                                      *
C  PARAMETERS:  LUN1     INPUT UNIT                                    *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE CENGR3(LUN1)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        INTEGER                 :: LUN1

        CHARACTER(LEN=MAXNAM)   :: FILNAM
        REAL, ALLOCATABLE       :: AIMG(:)
 
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NX,NY,NZ,MAXIM,
     &               'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
	
        IF (ITYPE .NE. 3 .AND. ITYPE .NE. 1)  THEN
           CLOSE(LUN1)
           CALL  ERRT(2,'CG, CENGR3',NE)
           RETURN
        ENDIF

        CALL RDPRM1S(TH, NOT_USED,'THRESHOLD',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        ALLOCATE (AIMG(NX), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'CENGR3, AIMG',NX)
           RETURN
        ENDIF

C       FIRST CALCULATE CENTER OF GRAVITY
        XM = 0.0
        YM = 0.0
        ZM = 0.0
        BT = 0.0

        DO L=1,NZ
           I1 = (L-1) * NY
           DO  I=1,NY
              CALL REDLIN(LUN1,AIMG,NX,I+I1)
              DO  K=1,NX
                 IF (AIMG(K) > TH) THEN
                    B  = AIMG(K)
                    BT = B  + BT
                    XM = XM + B * K
                    YM = YM + B * I
                    ZM = ZM + B * L
                 ENDIF  
              ENDDO
           ENDDO
        ENDDO

        IF (BT == 0.0) THEN 
           CALL ERRT(101,'NO CG FOUND FOR THIS THRESHOLD',NE)
           GOTO 999
        ENDIF

        ZM = ZM / BT
        YM = YM / BT
        XM = XM / BT
                

        IF (ITYPE == 3) THEN
C          NOW CALCULATE RADIUS OF GYRATION, RG
           BT = 0.0
           RG = 0.0
           DO  L=1,NZ
              I1 = (L-1)*NY
              DO I=1,NY
                 REST = (YM-I)**2 + (ZM-L) **2

                 CALL REDLIN(LUN1,AIMG,NX,I+I1)

                 DO K=1,NX
                    IF (AIMG(K) > TH) THEN
                       RG = RG + AIMG(K) * ((XM-K) ** 2+REST)
                       BT = BT + AIMG(K)
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
           IF (BT == 0.0) THEN 
              CALL ERRT(101,
     &           'NO RADIUS OF GYRATION FOUND FOR THIS THRESHOLD',NE)
              GOTO 999
           ENDIF

           RG = SQRT(RG / BT)
        ENDIF        

        XM = XM - NX/2 -1   
        YM = YM - NY/2 -1 
        ZM = ZM - NZ/2 -1 

        IF (ITYPE == 3) THEN
           WRITE(NOUT,201) XM,YM,ZM, RG
201        FORMAT('  Center Of Gravity:',
     &            '  X= ',F7.2,'  Y= ',F7.2,'  Z= ',F7.2,/,
     &     '  Radius Of Gyration= ',F12.2,/)

           CALL REG_SET_NSEL(1,4,XM,YM,ZM,RG, 0.0,IRTFLG)

        ELSE
           WRITE(NOUT,202) XM,YM
202        FORMAT('  Center Of Gravity:   X = ',F8.2,'  Y = ',F8.2,/)

           CALL REG_SET_NSEL(1,2,XM,YM, 0.0,0.0,0.0,IRTFLG)
        ENDIF

999     IF (ALLOCATED(AIMG)) DEALLOCATE (AIMG)

        END

