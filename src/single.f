
C++*********************************************************************
C
C    SINGLE
C                  MAXNAM                          JUL 14 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C    SINGLE
C
C--*********************************************************************

	SUBROUTINE SINGLE(NSAM,NROW)

 

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        COMMON /SIZE/ISIZ

        DIMENSION               PLIST(7)
        CHARACTER(LEN=MAXNAM):: DOCNAM
        CHARACTER               NULL

        NULL=CHAR(0)

        LUN1   = 10
        LUNDOC = 13

C	GET POS FROM LATTICE OPTION
        WRITE(NOUT,*) ' USE WINDOW DOCUMENT FOR REFLECTION LIST!'
	CALL FILERD(DOCNAM,NLET,NULL,'REFLECTION LIST',IRTFLG4)
	IF (IRTFLG4 .NE. 0)  RETURN

	CALL RDPRMI(IKEY,ID,NOT_USED,
     &     'FIRST KEY FOR REFL LIST (0 TO RETURN)')
	IF (IKEY.EQ.0) RETURN
	CALL RDPRMI(MODE,ID, NOT_USED,
     &       'MODE: MAX (1),CNTR OF DEN (2), NO BKG CORR (-)')
        NOPEN = 0
1	CONTINUE

        IFP = 0
	CALL UNSAV(DOCNAM,NOPEN,LUNDOC,IKEY,PLIST,6,LERR,1)
	NOPEN=1
	IF (LERR.NE.0) RETURN

        IH     = PLIST(1)
        IK     = PLIST(2)
	IXPOS  = PLIST(3)
	IYPOS  = PLIST(4)
        IFDIA  = PLIST(5)
        IFYDIA = PLIST(6)
	WRITE(NDAT,3028) IH,IK
	WRITE(NOUT,3028) IH,IK
3028	FORMAT(' FIRST PASS FOR H =',I5,', AND K = ',I5,' GIVES:')
3002	CALL RDPRMI(IXDIA,IYDIA,NOT_USED,'WINDOW SIZE(ODD #)')
	IF(IXDIA.EQ.0) GO TO 2
	IF(IYDIA.EQ.0)IYDIA=IXDIA
        IF(IFP.EQ.0) THEN
        IXDIA = IFDIA
        IYDIA = IFYDIA
        IFP = 1
        ENDIF
	CALL SPOTWT(MODE,LUN1,IXPOS,IYPOS,IXDIA,IYDIA,PA,PB,PC
     &      ,NEGC,KOX,KOY,CAVG,CMX,BKG,AAVG,AMX,NSAM,NROW,CTOT)
	WRITE(NDAT,3030)IXDIA,IYDIA
	WRITE(NOUT,3030)IXDIA,IYDIA
3030	FORMAT(' PASS FOR WINDOW SIZE ',I5,' BY ',I5,' GIVES:')
	WRITE(NDAT,3029)KOX,KOY,CAVG,CMX,PA,PB,PC,NEGC,BKG,AAVG,AMX,CTOT
3029	FORMAT(' X,Y:',2I6,' CORRECTED AVG,MAX:',2(G12.6,', '),/,
     &    ' PA,PB,PC',
     &    3(G10.3,', '),' NO. OF NEG PTS:',I6,/,' BKG:',G10.3,
     &    ', UNCORRECTED AVG:',G10.3,', UNCORRECTED MAX:',G10.3,/,
     &    ' ,TOTAL CORRECTED INTENSITY: ',G12.4//)

	WRITE(NOUT,3031)KOX,KOY,CAVG,CMX,CTOT
3031	FORMAT(' X,Y: ',2I6,' CORRECTED AVE, MAX, TOTAL: ',3G12.6)
        IXPOS = KOX
        IYPOS = KOY
        GO TO 3002

2	CONTINUE
	CALL RDPRMI(IKEY,ID,NOT_USED,
     &      'NEXT KEY FOR REFL LIST (0 TO RETURN)')
        IF (IKEY.NE.0)GO TO 1
        CLOSE(LUNDOC)

        RETURN
	END	
