C++*********************************************************************
C
C TSWITCH.F         FILENAMES LENGTHENED           JAN 89 ArDean Leith
C                   REWRITTEN                      MAR 90 ArDean Leith
C                   REWRITTEN                      MAR 93 ArDean Leith
C                   ADDED   'IQ'                   SEP 97 ArDean Leith
C                   ADDED   'AF'                   FEB 98 Pawel Penczek
C                   REMOVED 'PR'                   AUG 98 ArDean Leith
C                   ADDED   'IA'                   SEP 98 ArDean Leith
C                   ADDED   'NE'                   JUN 99 ArDean Leith
C                   ADDED   'FV'                   DEC 99 Pawel Penczek
C                   ADDED   'SO'                   MAR 00 ArDean Leith
C                   ADDED   'ER'                   FEB 01 ArDean Leith
C                   ADDED   'EV'                   APR 01 ArDean Leith
C                   ADDED   'PI'                   JUL 01 ArDean Leith
C                   ADDED   'VA'                   MAY 02 Pawel Penczek
C                   ADDED   'SN'                   MAY 02 Pawel Penczek
C                   ADDED   'LA'                   OCT 02 ArDean Leith
C                   ADDED   'MX'                   MAR 03 Bimal Rath
C                   ADDED   'DIV'  FOR 'MU D'      MAY 03 ArDean Leith
C                   ADDED   'SQRT' FOR 'WU'        MAY 03 ArDean Leith
C                   MOVED   'MD' TO SPIDER         DEC 03 ArDean Leith
C                   ADDED   'PB'                   JAN 04 Pawel Penczek
C                   ADDED   'WA'                   APR 04 ArDean Leith
C                   ADDED   'SY'                   APR 05 ArDean Leith
C                   ADDED   'TS'                   SEP 05 ArDean Leith
C                   ADDED   'DV'                   NOV 05 ArDean Leith
C                   REMOVED REG_READPQ             NOV 05 ArDean Leith
C                   RENAMED 'PB ..'                SEP 06 ArDean Leith
C                   ADDED   'RB ..'                DEC 06 ArDean Leith
C                   ADDED   'BPD ..'               JAN 07 ArDean Leith
C                   ADDED   'RTD ..'               JAN 07 ArDean Leith
C                   'BPD --> BP, BP --> OLD'       JUN 08 ArDean Leith
C                   'LO'                           JUL 08 ArDean Leith
C                   'XM'                           DEC 10 ArDean Leith
C                   'DN'                           FEB 11 ArDean Leith
C                   'ROT'                          SEP 11 ArDean Leith
C                   'IQ VER'                       JAN 12 ArDean Leith
C                   ADDED 'CENT' UTIL_1011         FEB 12 ArDean Leith
C                   ADDED 'FSC','FRC'              FEB 12 ArDean Leith
C                   'VM'='SYS'                     MAR 12 ArDean Leith
C                   ADDED UTIL8 'SH' 23            MAR 12 ArDean Leith
C                   REMOVED 'CTF'                  MAY 12 ArDean Leith
C                   REPLACED 'CTF'                 JUN 12 ArDean Leith
C                   'SPH'                          FEB 13 ArDean Leith
C                   'HIS'                          AUG 13 ArDean Leith
C                   'DIS'                          JAN 15 ArDean Leith
C                   'ML'                           JAN 15 ArDean Leith
C                   'RI'                           APR 16 ArDean Leith
C                   'MRC'                          NOV 19 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2019  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C     TSWITCH(IWHICH,ICOM,MAXDIM,IRTFLG)
C
C     PURPOSE:     SELECTS CORRECT TASK LEVEL CALLING ROUTINE
C
C     PARAMETERS:  IWHICH     TASK LEVEL SELECTION NUMBER    (RETURNED)
C                  ICOM       COMMAND NUMBER                 (RETURNED)
C                  MAXDIM     MAX. SIZE OF COMMON BUFFER     (SENT)
C                  IRTFLG     ERROR FLAG                     (RETURNED)
C
c    NOTE:         TO ADD A NEW COMMAND:
C                  1) INCREMENT NMENU PARAMETER
C                  2) NEW OPERATIONS WHICH DO NOT HAVE ANY SHORT 
C                     OPERATION, SHOULD USE A NUMBER FOR THE SHORT 
C                     OPERATION (TO SAVE REMAINING 2 LETTER SHORT
C                     OPERATIONS FOR THE FUTURE). 
C                  3) PLACE LONG OPERATION IN MENUL(NMENU), SHORT 
C                     OPERATION IN MENU AND THE CORRESPONDING 
C                     SUBROUTINE NUMBER IN IW(NMENU).  
C                  4) UPDATE THE OPERATIONS HANDLED TABLE DIRECTLY BELOW.
C                  5) ADD NEW NUMBER TO GOTO IF YOU HAVE ADDED A NEW
C                     MAIN SUBROUTINE
C                  6) EDIT CORRESPONDING SUBROUTINE TO ADD THE NECESSARY
C                     LISTING TO ITS MENU AND THE NECESSARY CALLS. THE 
C                     OPERATION SENT TO THE SUBROUTINE WILL BE THE 2
C                     CHAR. SHORT OPERATION AT ALL TIMES!  THE LONG
C                     OPERATION IS AVAILABLE IN COMMON /SP_OPER/
C
C     CODE      CALLER          OPERATIONS HANDLED
C	1	UTIL1		DE DU FI HI HD LI MO PK RA RN
C				TT ST TF FS CA GR CG CV CL CTF HIS
C	2	UTIL2		AD BL CP IN IP MU PA PD SQ SU DIV
C				WI CE AR MR DF MA WV PP SZ WU MM CM 
C				PV NK AS MN TH GP RP MX SQRT 14 31
C	3	UTIL3		AF ED RC RT BC CT OR FC SL RO OD MK OP
C                               RTD
C                               DI ER
C	4	CORR1		AC CC CN
C	5	UTIL7		XM DN
C	6	FOUR1		FT PW FF FL FP EF RF GF RD CF FQ FD LO
C	7	PLOT1		CO PF PL CS
C	8	TVWN1          	WT (REMOVED JAN 2012)
C	9	UTIL4		AP HF AT MS IQ NEG VA SN
C	10	UTIL5		RM
C	11	PLOT2		TP
C	12	VTIL2		PS SK PJ BP DC DR MF RB 
C	13	VTIL3		SF SC PH TA SI CR SP ML
C	14	CONF1		CI EP RS PT
C       15      DRIV1           NC VM ME CK TM SR RR PO SA VO AA FR DV
C                               EV PDB   RI
C       16      DRIV2           TV SE WA TS
C       17      DRIV3           UD LD SD SY      
C       18      SPIDER          EN DO LB EX RE IF GO MD OF
C       19      UTIL6		EC IA SO FV TS XM
C       20      MOD1            MY 
C       21      DOCS1           DOC
C       22      MOTIF1          LO
C       23      UTIL_1011       CEN
C       24      UTIL_1110       SH
C       25      UTIL_1010       NORM
C
C  156 DV unused Nov 2005
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE TSWITCH(IWHICH,ICOM,MAXDIM,IRTFLG)

      INCLUDE 'CMBLOCK.INC'

      CHARACTER (LEN=20) :: FCHARS
      COMMON /SP_OPER/      FCHARS

C     NUMBER OF OPERATIONS IN MENU 
      INTEGER, PARAMETER :: NMENU = 173 ! HIGHEST 2 CHAR NUMBER=31

      CHARACTER(LEN=2)   :: MENU(NMENU)
      CHARACTER(LEN=12)  :: MENUL(NMENU)
      INTEGER            :: IW(NMENU)
    
      DATA MENU(124), MENUL(124), IW(124) /'AA','  ',15/
      DATA MENU(1),   MENUL(1),   IW(1)   /'AC','  ',4/
      DATA MENU(2),   MENUL(2),   IW(2)   /'AD','ADD',2/
      DATA MENU(136), MENUL(136), IW(136) /'AF','  ',3/
      DATA MENU(3),   MENUL(3),   IW(3)   /'AP','  ',9/
      DATA MENU(4),   MENUL(4),   IW(4)   /'AR','  ',2/
      DATA MENU(5),   MENUL(5),   IW(5)   /'AS','  ',2/
      DATA MENU(129), MENUL(129), IW(129) /'AT','  ',9/
      DATA MENU(6),   MENUL(6),   IW(6)   /'BC','  ',3/
      DATA MENU(7),   MENUL(7),   IW(7)   /'BL','BLANK',2/
      DATA MENU(8),   MENUL(8),   IW(8)   /'BP','  ',12/
      DATA MENU(9),   MENUL(9),   IW(9)   /'CA','  ',1/
      DATA MENU(10),  MENUL(10),  IW(10)  /'CC','  ',4/
      DATA MENU(11),  MENUL(11),  IW(11)  /'CE','  ',2/
      DATA MENU(163), MENUL(163), IW(163) /'15','CENTER',23/
      DATA MENU(12),  MENUL(12),  IW(12)  /'CF','  ',6/
      DATA MENU(13),  MENUL(13),  IW(13)  /'CG','  ',1/
      DATA MENU(14),  MENUL(14),  IW(14)  /'CI','  ',14/
      DATA MENU(15),  MENUL(15),  IW(15)  /'CK','  ',15/
      DATA MENU(125), MENUL(125), IW(125) /'CL','  ',1/
      DATA MENU(16),  MENUL(16),  IW(16)  /'CM','  ',2/
      DATA MENU(17),  MENUL(17),  IW(17)  /'CN','  ',4/
      DATA MENU(18),  MENUL(18),  IW(18)  /'CO','  ',7/
      DATA MENU(19),  MENUL(19),  IW(19)  /'CP','COPY',2/
      DATA MENU(20),  MENUL(20),  IW(20)  /'CR','  ',13/
      DATA MENU(21),  MENUL(21),  IW(21)  /'CS','  ',7/
      DATA MENU(22),  MENUL(22),  IW(22)  /'CT','  ',3/
      DATA MENU(166), MENUL(166), IW(166) /'18','CTF',1/
      DATA MENU(23),  MENUL(23),  IW(23)  /'CV','  ',1/
      DATA MENU(24),  MENUL(24),  IW(24)  /'DC','  ',12/
      DATA MENU(25),  MENUL(25),  IW(25)  /'DE','  ',1/
      DATA MENU(26),  MENUL(26),  IW(26)  /'DF','  ',2/
      DATA MENU(142), MENUL(142), IW(142) /'DI','  ',3/
      DATA MENU(152), MENUL(152), IW(152) /'12','DIV',2/
      DATA MENU(170), MENUL(170), IW(170) /'30','DIS',2/
      DATA MENU(161), MENUL(161), IW(161) /'DN','  ',5/
      DATA MENU(118), MENUL(118), IW(118) /'DO','  ',18/
      DATA MENU(139), MENUL(139), IW(139) /'11','DOC',21/
      DATA MENU(27),  MENUL(27),  IW(27)  /'DR','  ',12/
      DATA MENU(28),  MENUL(28),  IW(28)  /'DU','  ',1/
      DATA MENU(156), MENUL(156), IW(156) /'DV','  ',15/
      DATA MENU(128), MENUL(128), IW(128) /'EC','  ',19/
      DATA MENU(29),  MENUL(29),  IW(29)  /'ED','  ',3/
      DATA MENU(30),  MENUL(30),  IW(30)  /'EF','  ',6/
      DATA MENU(117), MENUL(117), IW(117) /'EN','END',18/
      DATA MENU(31),  MENUL(31),  IW(31)  /'EP','  ',14/
      DATA MENU(145), MENUL(145), IW(145) /'ER','  ',3/
      DATA MENU(146), MENUL(146), IW(146) /'EV','  ',15/
      DATA MENU(119), MENUL(119), IW(119) /'EX','  ',18/
      DATA MENU(32),  MENUL(32),  IW(32)  /'FC','  ',3/
      DATA MENU(131), MENUL(131), IW(131) /'FD','  ',6/
      DATA MENU(33),  MENUL(33),  IW(33)  /'FF','  ',6/
      DATA MENU(34),  MENUL(34),  IW(34)  /'FI','  ',1/
      DATA MENU(35),  MENUL(35),  IW(35)  /'FL','  ',6/
      DATA MENU(36),  MENUL(36),  IW(36)  /'FP','  ',6/
      DATA MENU(37),  MENUL(37),  IW(37)  /'FQ','  ',6/
      DATA MENU(38),  MENUL(38),  IW(38)  /'FR','  ',15/
      DATA MENU(165), MENUL(165), IW(165) /'17','FRC',6/
      DATA MENU(39),  MENUL(39),  IW(39)  /'FS','  ',1/
      DATA MENU(164), MENUL(164), IW(164) /'16','FSC',6/
      DATA MENU(40),  MENUL(40),  IW(40)  /'FT','  ',6/
      DATA MENU(141), MENUL(141), IW(141) /'FV','  ',19/
      DATA MENU(41),  MENUL(41),  IW(41)  /'GF','  ',6/
      DATA MENU(130), MENUL(130), IW(130) /'GO','  ',18/
      DATA MENU(42),  MENUL(42),  IW(42)  /'GP','  ',2/
      DATA MENU(43),  MENUL(43),  IW(43)  /'GR','  ',1/
      DATA MENU(137), MENUL(137), IW(137) /'HD','  ',1/
      DATA MENU(44),  MENUL(44),  IW(44)  /'HI','  ',1/
      DATA MENU(169), MENUL(169), IW(169) /'20','HIS',1/
      DATA MENU(126), MENUL(126), IW(126) /'HF','  ',9/
      DATA MENU(121), MENUL(121), IW(121) /'IF','  ',18/
      DATA MENU(71),  MENUL(71),  IW(71)  /'IA','  ',19/
      DATA MENU(45),  MENUL(45),  IW(45)  /'IN','  ',2/
      DATA MENU(46),  MENUL(46),  IW(46)  /'IP','INTERP',2/
      DATA MENU(134), MENUL(134), IW(134) /'IQ','  ',9/
      DATA MENU(150), MENUL(150), IW(150) /'LA','  ',19/
      DATA MENU(122), MENUL(122), IW(122) /'LB','  ',18/
      DATA MENU(47),  MENUL(47),  IW(47)  /'LD','  ',17/
      DATA MENU(48),  MENUL(48),  IW(48)  /'LI','LIST',1/
      DATA MENU(159), MENUL(159), IW(159) /'LO','LOC',22/
      DATA MENU(49),  MENUL(49),  IW(49)  /'MA','  ',2/
      DATA MENU(50),  MENUL(50),  IW(50)  /'WA','WARP',16/
      DATA MENU(51),  MENUL(51),  IW(51)  /'ME','MENU',15/
      DATA MENU(52),  MENUL(52),  IW(52)  /'MF','  ',12/
      DATA MENU(123), MENUL(123), IW(123) /'MK','  ',3/
      DATA MENU(171), MENUL(171), IW(171) /'ML','  ',13/
      DATA MENU(53),  MENUL(53),  IW(53)  /'MM','  ',2/
      DATA MENU(54),  MENUL(54),  IW(54)  /'MN','  ',2/
      DATA MENU(55),  MENUL(55),  IW(55)  /'MO','MODEL',1/
      DATA MENU(56),  MENUL(56),  IW(56)  /'MR','  ',2/
      DATA MENU(173), MENUL(173), IW(173) /'31','MRC',2/
      DATA MENU(132), MENUL(132), IW(132) /'MS','  ',9/
      DATA MENU(57),  MENUL(57),  IW(57)  /'MU','MULT',2/
      DATA MENU(151), MENUL(151), IW(151) /'MX','  ',2/
      DATA MENU(133), MENUL(133), IW(133) /'MY','  ',20/
      DATA MENU(58),  MENUL(58),  IW(58)  /'NC','  ',15/
      DATA MENU(140), MENUL(140), IW(140) /'NE','NEG',9/
      DATA MENU(59),  MENUL(59),  IW(59)  /'NK','  ',2/
      DATA MENU(167), MENUL(167), IW(167) /'NO','NORM',25/
      DATA MENU(60),  MENUL(60),  IW(60)  /'OD','  ',3/
      DATA MENU(138), MENUL(138), IW(138) /'OP','  ',3/
      DATA MENU(61),  MENUL(61),  IW(61)  /'OR','  ',3/
      DATA MENU(62),  MENUL(62),  IW(62)  /'PA','PATCH',2/
      DATA MENU(153), MENUL(153), IW(153) /'PB','PDB',15/
      DATA MENU(63),  MENUL(63),  IW(63)  /'PD','PAD',2/
      DATA MENU(64),  MENUL(64),  IW(64)  /'PF','  ',7/
      DATA MENU(65),  MENUL(65),  IW(65)  /'PH','  ',13/
      DATA MENU(147), MENUL(147), IW(147) /'PI','  ',15/
      DATA MENU(66),  MENUL(66),  IW(66)  /'PJ','PROJ',12/
      DATA MENU(67),  MENUL(67),  IW(67)  /'PK','  ',1/
      DATA MENU(68),  MENUL(68),  IW(68)  /'PL','  ',7/
      DATA MENU(69),  MENUL(69),  IW(69)  /'PO','  ',15/
      DATA MENU(70),  MENUL(70),  IW(70)  /'PP','  ',2/
      DATA MENU(72),  MENUL(72),  IW(72)  /'PS','  ',12/
      DATA MENU(73),  MENUL(73),  IW(73)  /'PT','  ',14/
      DATA MENU(74),  MENUL(74),  IW(74)  /'PV','  ',2/
      DATA MENU(75),  MENUL(75),  IW(75)  /'PW','  ',6/
      DATA MENU(76),  MENUL(76),  IW(76)  /'RA','  ',1/
      DATA MENU(157), MENUL(157), IW(157) /'RB','  ',12/
      DATA MENU(77),  MENUL(77),  IW(77)  /'RC','  ',3/
      DATA MENU(78),  MENUL(78),  IW(78)  /'RD','  ',6/
      DATA MENU(120), MENUL(120), IW(120) /'RE','RETURN',18/
      DATA MENU(79),  MENUL(79),  IW(79)  /'RF','  ',6/
      DATA MENU(80),  MENUL(80),  IW(80)  /'RM','  ',10/
      DATA MENU(81),  MENUL(81),  IW(81)  /'RN','  ',1/
      DATA MENU(162), MENUL(162), IW(162) /'14','ROT',3/
      DATA MENU(82),  MENUL(82),  IW(82)  /'RO','  ',3/
      DATA MENU(135), MENUL(135), IW(135) /'RP','  ',2/
      DATA MENU(172), MENUL(172), IW(172) /'RI','  ',15/
      DATA MENU(83),  MENUL(83),  IW(83)  /'RR','  ',15/
      DATA MENU(84),  MENUL(84),  IW(84)  /'RS','  ',14/
      DATA MENU(85),  MENUL(85),  IW(85)  /'RT','  ',3/
      DATA MENU(158), MENUL(158), IW(158) /'13','RTD',3/
      DATA MENU(86),  MENUL(86),  IW(86)  /'SA','  ',15/
      DATA MENU(87),  MENUL(87),  IW(87)  /'SC','  ',13/
      DATA MENU(88),  MENUL(88),  IW(88)  /'SD','SAVE',17/
      DATA MENU(127), MENUL(127), IW(127) /'SE','  ',16/
      DATA MENU(89),  MENUL(89),  IW(89)  /'SF','  ',13/
      DATA MENU(90),  MENUL(90),  IW(90)  /'SH','SHIFT',24/
      DATA MENU(91),  MENUL(91),  IW(91)  /'SI','  ',13/
      DATA MENU(92),  MENUL(92),  IW(92)  /'SK','  ',12/
      DATA MENU(93),  MENUL(93),  IW(93)  /'SL','  ',3/
      DATA MENU(149), MENUL(149), IW(149) /'SN','  ',9/
      DATA MENU(144), MENUL(144), IW(144) /'SO','  ',19/
      DATA MENU(168), MENUL(168), IW(168) /'19','SPH',13/
      DATA MENU(94),  MENUL(94),  IW(94)  /'SP','  ',13/
      DATA MENU(95),  MENUL(95),  IW(95)  /'SQ','  ',2/
      DATA MENU(96),  MENUL(96),  IW(96)  /'SR','  ',15/
      DATA MENU(97),  MENUL(97),  IW(97)  /'ST','SET',1/
      DATA MENU(98),  MENUL(98),  IW(98)  /'SU','SUB',2/
      DATA MENU(154), MENUL(154), IW(154) /'SY','  ',17/
      DATA MENU(99),  MENUL(99),  IW(99)  /'SZ','  ',2/
      DATA MENU(100), MENUL(100), IW(100) /'TA','  ',13/
      DATA MENU(101), MENUL(101), IW(101) /'TF','  ',1/
      DATA MENU(102), MENUL(102), IW(102) /'TH','  ',2/
      DATA MENU(103), MENUL(103), IW(103) /'TM','  ',15/
      DATA MENU(105), MENUL(105), IW(105) /'TP','  ',11/
      DATA MENU(106), MENUL(106), IW(106) /'TR','  ',5/
      DATA MENU(155), MENUL(155), IW(155) /'TS','  ',16/
      DATA MENU(107), MENUL(107), IW(107) /'TT','  ',1/
      DATA MENU(108), MENUL(108), IW(108) /'TV','  ',16/
      DATA MENU(109), MENUL(109), IW(109) /'TW','  ',5/
      DATA MENU(110), MENUL(110), IW(110) /'UD','UNSAVE',17/
      DATA MENU(148), MENUL(148), IW(148) /'VA','  ',9/
      DATA MENU(111), MENUL(111), IW(111) /'VM','SYS',15/
      DATA MENU(112), MENUL(112), IW(112) /'VO','  ',15/
      DATA MENU(113), MENUL(113), IW(113) /'WI','WINDOW',2/
      DATA MENU(114), MENUL(114), IW(114) /'WT','  ',8/
      DATA MENU(115), MENUL(115), IW(115) /'WU','SQRT',2/
      DATA MENU(116), MENUL(116), IW(116) /'WV','  ',2/
      DATA MENU(160), MENUL(160), IW(160) /'XM','  ',5/

      IRTFLG   = 0
      ICOM     = 0
      NLETOP1  = LNBLNKN(FCHAR)
      IF (NLETOP1 <= 0) THEN
         IRTFLG = 1
         RETURN
      ENDIF
      IGO      = INDEX(FCHAR(1:NLETOP1),' ')
      IF (IGO > 1) NLETOP1 = IGO - 1
 
      IF (NLETOP1 == 2) THEN
C        TWO LETTER OPERATION
         DO I = 1,NMENU
            IF (FCHAR(1:2) == MENU(I)(1:2)) THEN
               ICOM   = I
               IWHICH = IW(I)
               FCHARS(1:3) = FCHAR(1:2) // CHAR(0)
            ENDIF
         ENDDO
      ENDIF

      IF (ICOM == 0) THEN
C       OPERATION NOT FOUND IN THE SHORT OPERATION MENU

        LENOP = INDEX(FCHAR,' ') - 1
        IF (LENOP == 0) LENOP = INDEX(FCHAR,CHAR(0)) - 1
        IF (LENOP == 0) LENOP = LEN(FCHAR) - 1
        
        DO I = 1,NMENU
           IF (FCHAR(1:LENOP) == MENUL(I)(1:LENOP)) THEN
              IF (ICOM > 0) THEN
                 ICOM   = 0
                 IWHICH = 0
                 CALL ERRT(101, 
     &            '*** AMBIGUOUS OPERATION - SUPPLY MORE CHARACTERS.',N)
                 IRTFLG = 1
                 RETURN
              ENDIF
              ICOM   = I
              IWHICH = IW(I)
            ENDIF
        ENDDO

        IF (ICOM > 0) THEN
C          SAVE OPERATION IN FCHARS
           FCHARS = FCHAR(1:LENOP) // CHAR(0)

C          SET FCHAR TO SHORT MENU OPERATION
           FCHAR(1:3) = MENU(ICOM)(1:2) // ' '

C          MOVE SUBOPTIONS TO POSITION 4 IN FCHAR
           FCHAR(4:) = FCHAR(LENOP+2:)
        ENDIF
      ENDIF

      IF (ICOM == 0) THEN
         IRTFLG = 1
         RETURN
      ENDIF

      
C     SOME OPERATIONS SHOULD BE SILENT IN BATCH MODE IF .NOT. VERBOSE
      IF ((ICOM == 88 .OR. ICOM == 110) .AND. .NOT. VERBOSE) THEN
         SILENT = .TRUE.
      ELSE
         SILENT = .FALSE.
      ENDIF
    
      SELECT CASE(IWHICH)

      CASE(1)
         IF (FCHAR(1:2) == 'PR') THEN
            CALL ERRT(101,'OBSOLETE OPERATION REMOVED',NE)
         ELSE
            CALL UTIL1(MAXDIM,IRTFLG)
            IF (IRTFLG == 1) THEN
C              used "FI" instead of "FR"
               FCHAR = 'FR'
               RETURN 
            ENDIF
         ENDIF

      CASE(2)
         CALL UTIL2(MAXDIM) 

      CASE(3)
         CALL UTIL3(MAXDIM) 

      CASE(4)
         CALL CORR1

      CASE(5)
         CALL UTIL7()

      CASE(6)
         CALL FOUR1(MAXDIM) 

      CASE(7)
         CALL PLOT1(MAXDIM) 

      CASE(8)
         IF (FCHAR(4:4) == 'T') THEN
            CALL DRIV2(MAXDIM)
         ENDIF

      CASE(9)
         IF (FCHAR(1:6) == 'IQ VER') THEN
            IRTFLG = 2 
         ELSE
            CALL UTIL4(MAXDIM)  
         ENDIF

      CASE(10)
         CALL UTIL5(MAXDIM) 

      CASE(11)
         CALL PLOT2(MAXDIM) 

      CASE(12)
         CALL VTIL2(MAXDIM) 

      CASE(13)
         CALL VTIL3(MAXDIM) 

      CASE(14)
         CALL CONF1(MAXDIM) 

      CASE(15)
         CALL DRIV1(MAXDIM)

      CASE(16)
         CALL DRIV2(MAXDIM)

      CASE(17)
         CALL DRIV3(MAXDIM)

      CASE(18)
C        SPECIAL CASE FOR OPERATIONS IN SPIDER
         IRTFLG = 2

      CASE(19)      
         CALL UTIL6(MAXDIM)

      CASE(20)
         CALL MOD1(MAXDIM)

      CASE(21)
         CALL DOCS1(MAXDIM)

      CASE(22)
         CALL MOTIF1()

      CASE(23)
         CALL UTIL_1011()

      CASE(24)
          CALL UTIL_1110()

      CASE(25)
          CALL UTIL_1010()

      END SELECT

      END
