
C ++********************************************************************
C                                                                      *
C                                                                      *
C                                                                      *
C **********************************************************************
C UNIFORM PSEUDO RANDOM NUMBER GENERATOR OVER A RANGE OF 0.0 TO 1.0             
C MULTIPLICATIVE CONGRUENT METHOD U(N+1)=A*U(N) MODULO (2**M)                   
C    WITH  M = 36   A = 64155708247  AND U(0) = A .                             
C NUMBERS ARE SEGMENTED SO THAT NUMERICAL OVERFLOWS ARE CALCULATED              
C WITHOUT HARDWARE OVERFLOW. THUS, THE RESULTS ARE MACHINE INDEPENDENT.         
C THE SPLITTING IN 3 SEGMENTS REQUIRES AT LEAST MACHINE WORDS OF 26 BITS.       
C REFERENCE: K.D. SENNE, J. STOCHASTICS VOL 1, NO 3 (1974), PP 215-238.         
C ADAPTED FROM A FORTRAN PROGRAM BY L. LEBART, A. MORINEAU & J.P. FENELON       
C IN TRAITEMENT DES DONNEES STATISTIQUES, DUNOD PUBL, PARIS, 1979.              
C-------------------------------------------------------------------------      

      FUNCTION  SEN3A ( BIDON )                                                 

C     I DO NOT KNOW IF SAVE IS NEEDED FEB 99 al
      SAVE

C     FIXED DATA                                                           
      DATA M12/ 4096 /                                                          
      DATA F1/2.44140625E-04/,F2/5.96046448E-08/,F3/1.45519152E-11/             
      DATA J1/ 3823 /        ,J2/ 4006 /        ,J3/ 2903 /                     
C     INITIAL DATA THEN CALCULATED VALUES                                  
      DATA I1/ 3823 /        ,I2/ 4006 /        ,I3/ 2903 /                     
C     CONGRUENCE CALCULATION WITH NUMBER SEGMENTATION                      
      K3 = I3*J3                                                                
      L3 = K3/M12                                                               
      K2 = I2*J3 + I3*J2 + L3                                                   
      L2 = K2/M12                                                               
      K1 = I1*J3 + I2*J2 + I3*J1 + L2                                           
      L1 = K1/M12                                                               
      I1 = K1 - L1*M12                                                          
      I2 = K2 - L2*M12                                                          
      I3 = K3 - L3*M12                                                          
      SEN3A = F1*FLOAT(I1)+ F2*FLOAT(I2)+F3*FLOAT(I3)                           
      RETURN                                                                  
      END

