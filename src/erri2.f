
C++*********************************************************************
C
C ERRI2.F  -- CREATED JULY 15 1987
C
C **********************************************************************
C *  AUTHOR:  ARDEAN LEITH
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2018  Health Research Inc.,                         *
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
C    ERRI2(NUM1,NUM2,NVAL,ILOW1,IHI1,ILOW2,IHI2)
C
C    PARAMETERS:    NUM1         FIRST INTEGER
C                   NUM2         SECOND INTEGER
C                   NVAL         NUMBER OF VALUES TO BE CHECKED
C                   ILOW1        LOWEST VALUE FOR  NUM1
C                   IHI1         HIGHEST VALUE FOR NUM1
C                   ILOW2        LOWEST VALUE FOR  NUM2
C                   IHI2         HIGHEST VALUE FOR NUM2
C
C    CALLED BY:     
C
C--*******************************************************************

        LOGICAL FUNCTION ERRI2(NUM1,NUM2,NVAL,ILOW1,IHI1,ILOW2,IHI2)

        COMMON /UNITS/ LUNDOC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT
        
        INTEGER  :: NUM1,NUM2,NVAL,ILOW1,IHI1,ILOW2,IHI2
        
        INTEGER  :: IDUM

        ERRI2 = .FALSE.

        IF (NVAL == 1) THEN
          IF (NUM1 < ILOW1 .OR. NUM1 > IHI1) THEN
            WRITE(NOUT,91) ILOW1,IHI1
   91       FORMAT(' ERROR, INPUT RANGE: (',I0,'....',I0,')',/)
            CALL ERRT(100,'ERRI2',IDUM)
            ERRI2 = .TRUE.
            RETURN
          ENDIF

        ELSEIF (NUM1 < ILOW1 .OR. NUM1 > IHI1) THEN
            WRITE(NOUT,92) ILOW1,IHI1
   92       FORMAT(' ERROR, FIRST INPUT RANGE: (',I0,'....',I0,')',/)
            CALL ERRT(100,'ERRI2',IDUM)
            ERRI2 = .TRUE.
        ENDIF

        IF (NVAL >= 2 .AND. (NUM2 < ILOW2 .OR. NUM2 > IHI2)) 
     &      THEN
            WRITE(NOUT,97) ILOW2,IHI2
   97       FORMAT(' ERROR, SECOND INPUT RANGE: (',I0,'...',I0,')',/)
            CALL ERRT(100,'ERRI2',IDUM)
            ERRI2 = .TRUE.
        ENDIF

        END
