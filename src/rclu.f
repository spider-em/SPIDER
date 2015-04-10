C++*********************************************************************
C
C RCLU.F   DERIVED FROM CODE AUTHORED BY: JEAN-PIERRE BRETAUDIERE                          
C         
C **********************************************************************
C * SPIDER - MODULAR IMAGE PROCESSING SYSTEM.    AUTHOR: J.FRANK       *
C * COPYRIGHT 1986 - JEAN-PIERRE BRETAUDIERE                           *                          
C *   THE UNIVERSITY OF TEXAS HEALTH SCIENCE CENTER AT HOUSTON         *                     
C *   MEDICAL SCHOOL - DEPARTMENT OF PATHOLOGY AND LABORATORY MEDICINE *              
C *   P.O. BOX 20708, HOUSTON, TX 77225.                               *
C **********************************************************************
C                                                                     
C RCLU(NGRI,NGUS,MAXMEM)
C PARAMETERS: NGRI   = LOGICAL UNIT NUMBER OF CLUSTER FILE, 
C                      SUPPOSED TO BE OPEN
C PARAMETERS: NGUS   = LOGICAL UNIT NUMBER 
C             MAXMEM = BUFFER SIZE AVAILABLE
C
C MODIFIED VERSION TO INTERFACE WITH SPIDER J.F. 2/5/86
C
C **********************************************************************

	SUBROUTINE RCLU(NGRI,NGUS,MAXMEM)

	INCLUDE 'CMBLOCK.INC' 

	COMMON Q(1)
        COMMON /ENSOR/ LEC, IMP

C       MAXIMUM NUMBER OF FACTORS
	PARAMETER (LFA=16)
        INTEGER   KV(LFA)

	LEC = NIN
	IMP = NDAT

	REWIND NGRI                                                                        
        READ (NGRI)  ICARD, NFAC, NKLA, KFAC, (KV(K), K=1, KFAC)                  
        NK    = 1                                                                 
        NI    = NK + ICARD                                                        
        NPK   = NI + ICARD                                                        
        NGT   = NPK   + NKLA                                                      
        NIV   = NGT   + NKLA * KFAC                                               
        MEM   = NIV   + NKLA                                                      
        CALL RGRI(ICARD, KFAC, NKLA,                                              
     &          Q(NK),Q(NI),Q(NPK),Q(NGT),Q(NIV), NGRI,NGUS)

        CLOSE(NGRI)
	CLOSE(NGUS)

	RETURN
        END                                                                       
