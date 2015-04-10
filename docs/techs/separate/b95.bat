; MAKES MASKBALL FOR SMALL SUB UNIT 40S; USE JWEB ON output/vol_fq.dat 
; DEREK'S ORIGINAL VOLUME IS  205 X 205 X 205. INTERPOLATED VOLUME OF SIZE 
; 102 X 102 X 102 IS USED

; X = 32 - 93 + 1 = 62
; Y = 10 - 68 + 1 = 59
; Z = 24 - 72 + 1 = 49
; CENTER : X = 32 + 31 = 63; Y = 10 + 30 = 40 ; Z = 24 + 25 = 49
; RADIUS = 62/2 + 1= 32


MO 3
output/maskball
102,102,102
SP
N
1
32,0
63,40,49
0,0,0

; TEST IT BY LOOKING AT junk_test IN JWEB
MU
output/vol_fq
output/maskball
output/junk_test

RE
