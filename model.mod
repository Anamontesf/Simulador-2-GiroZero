param PRICE{i in I, t in 0..T} default PVEH[i]*GAMMA[ceil(t/12),i]*(1+SA[ceil(t/12),i])*(1+DEV)^t;
param BIGM default 10^6;

var Y{v in 1..Veh, i in I, t in 0..T} binary; # purchasing time
var X{v in 1..Veh, i in I, t in 1..T+1} binary; # selling time
var Z{v in 1..Veh, i in I, t in 0..T} binary; # permanence time
var W{v in 1..Veh, i in I, t in 0..T} integer >=0 <=ceil(T/(12*6)); # age range

var FLOW{i in I, t in 0..T+1} >=0;
var KM{v in 1.. Veh, i in I, t in 0..T} >=0;
var ASSET{v in 1.. Veh, i in I, t in 0..T} >=0;
var CO2{i in I, t in 0..T} >=0;
var M1{v in 1.. Veh, i in I, t in 0..T} >=0;
var TC{i in I, t in 0..T} >=0;

var F{i in I, t in 0..T} >=0;
var D{i in I, t in 0..T} >=0;
var A{i in I, t in 0..T} >=0;
var M{i in I, t in 0..T} >=0;
var B{i in I, t in 0..T+1} >=0;

var B1{f in FT, i in I, t in 0..T+1} >=0;
var B2{f in FT, i in I, t in 0..T+1} >=0;
var B3{f in FT, i in I, t in 0..T+1} >=0;
var B4{i in I, t in 0..T+1} >=0;

var BUY{f in FT, v in 1..Veh, i in I, t in 0..T} binary;
var AGE{v in 1..Veh, i in I, t in 0..T} integer >=0;


# Objective function
minimize NPV:
sum{t in 0..T+1, i in I} FLOW[i,t]*(1/((1+WACC)^t))
;

subject to R1a{t in 0..T, i in I}:
FLOW[i,t]=(F[i,t]+A[i,t]+M[i,t]+B[i,t]-D[i,t])
;

subject to R0a{t in 1..T, i in I}:
D[i,t]=EMP*sum{v in 1..Veh}Y[v,i,t]*(PRICE[i,t]*IMP/(DPED*12))*((1-(1+WACC)^(-DPED*12))/WACC)
;

subject to R0b{t in {0}, i in I}:
D[i,t]=0
;

subject to R1b{t in {T+1}, i in I}:
FLOW[i,t]=B[i,t]
;


subject to R5a{i in I, t in 0..T: i in {"HID","ELE"}}:
F[i,t]=(sum{v in 1..Veh}KM[v,i,t]*(1/(EFF[i,TYPE]*FFAC*EVOL[ceil(t/12),i]))*MFAC*(CF_0[i]*(1+SF[ceil(t/12),i])*LAMBDA[ceil(t/12),i]*((1+INF)^t)))/(10^6)
;

subject to R5b{i in I, t in 0..T: i in {"GAS","DIE"}}:
F[i,t]=(sum{v in 1..Veh}KM[v,i,t]*(1/(EFF[i,TYPE]*FFAC*EVOL[ceil(t/12),i]))*MFAC*(CF_0[i]*LAMBDA[ceil(t/12),i]*((1+INF)^t))+TC[i,t])/(10^6);

subject to R5c{i in I, t in 0..T}:
TC[i,t]=TAX[ceil(t/12)]*((1+DEV)^t)*sum{v in 1..Veh}KM[v,i,t]*(1/(EFF[i,TYPE]*FFAC*EVOL[ceil(t/12),i]))*(TTW[i])*MFAC/(10^6);

subject to R9{t in 0..T+1, i in I}:
B[i,t]=sum{f in FT}(B1[f,i,t]+B2[f,i,t]+B3[f,i,t])-B4[i,t]
;

subject to R10a{f in FT, t in 1..T, i in I}:
B1[f,i,t]=sum{v in 1..Veh}BUY[f,v,i,t]*PRICE[i,t]*P[f,1]
;

subject to R10c{f in FT, i in I, t in {0,T}}:
B1[f,i,t]=0
;

subject to 10b{i in I, v in 1..Veh, t in 0..T}:
sum{f in FT}BUY[f,v,i,t]=Y[v,i,t]
;

subject to R12a{f in FT, t in 1..T, i in I: P[f,2]>0}:
B2[f,i,t]=sum{t2 in max(1,(t-n[f]+1))..t, v in 1..Veh}BUY[f,v,i,t2]*PRICE[i,t2]*P[f,2]*(r[f]/(1-(1+r[f])^(-n[f])))
;



subject to R12b{f in FT, i in I, t in {0,T}: P[f,2]>0}:
B2[f,i,t]=0
;


subject to R12c{f in FT, t in 1..T+1, i in I: t-n[f] >=1 && t-n[f]<T+1}:
B3[f,i,t]=sum{v in 1..Veh}BUY[f,v,i,t-n[f]]*PRICE[i,t-n[f]]*P[f,3]
;

subject to R12d{f in FT, i in I, t in {0}}:
B3[f,i,t]=0
;



subject to R12e{f in FT, t in 1..T+1, i in I}:
B4[i,t]=RESALE*PRICE[i,t-1]*sum{v in 1..Veh}X[v,i,t]
;



subject to R12f{i in I, t in {0}}:
B4[i,t]=0
;


subject to 13a{i in I, v in 1..Veh, t in 1..T}:
ASSET[v,i,t]=PRICE[i,0]*Y[v,i,0]*(1-GDEP)^(12*AGE_0[i,v])+sum{t1 in 1..t}PRICE[i,t1]*Y[v,i,t1]*(1-GDEP)^(t-t1)-sum{t2 in 1..t}X[v,i,t2]*max{t3 in 0..T}PRICE[i,t3];


subject to 13C{i in I, v in 1..Veh, t in {0}}:
ASSET[v,i,t]=PRICE[i,t]*Y[v,i,t]*(1-GDEP)^(12*AGE_0[i,v]);


subject to R14{t in 0..T, i in I}:
A[i,t]=sum{v in 1..Veh}ASSET[v,i,t]*OMEGA
;

subject to R31c{i in I, v in 1..Veh, t in 1..T}:
AGE[v,i,t]<=6*W[v,i,t];	

subject to R31a{i in I, v in 1..Veh, t in 1..T}:
M1[v,i,t]=(CM_0*Z[v,i,t]+(W[v,i,t]-1)*DELTA[TYPE])*((1+INF)^t)*MFAC*FM[i];	

subject to R31d{t in 0..T, i in I}:
M[i,t]=KMD*OPD*REL*sum{v in 1..Veh}M1[v,i,t]/(10^6);	

subject to R41b{t in 1..T, i in I, v in 1..Veh}:
(T+1)*(Z[v,i,t])>=AGE[v,i,t];

subject to R41a{t in 1..T, i in I, v in 1..Veh}:
AGE[v,i,t]>=AGE[v,i,t-1]+1-(T+1)*(1-Z[v,i,t])
;


subject to R41c{i in I, v in 1..Veh}:
AGE[v,i,0]=12*AGE_0[i,v]
;


subject to R41d{t in 1..T, i in I, v in 1..Veh}:
12*min(LIM1[i],LIM2)>=AGE[v,i,t]
;




subject to R4{t in 0..T}:
sum{i in I, v in 1..Veh}KM[v,i,t]=KM_0*(1+RG)^t
;


subject to R6{v in 1..Veh, i in I, t in 0..T}:
KM[v,i,t]<=Z[v,i,t]*KMD*OPD;


subject to R11{i in I, v in 1..Veh}:
Y[v,i,0]=INV[i,v];


subject to R2a{i in I, v in 1..Veh}:
sum{t2 in 1..T+1}X[v,i,t2]=sum{t1 in 0..T}Y[v,i,t1]
;


subject to R18{t in 1..ceil(T/12)}:
OPEXBUDGET*(1+RG1)^t>= sum{i in I, t2 in (1+(t-1)*12)..12*t}(F[i,t2]+A[i,t2]+M[i,t2])
;


subject to R19{t in 1..ceil(T/12)}:
CAPEXBUDGET*(1+RG1)^t>= sum{i in I, t2 in (1+(t-1)*12)..12*t}(B[i,t2])
;


subject to R51{t in 1..T,v in 1..Veh}:
Y[v,"ELE",t]<=1-NAC;


subject to R21{v in 1..Veh, t in 1..T, i in I: YEAR[i,1] > ceil(t/12) and ceil(t/12) > YEAR[i,2]}:
Y[v,i,t]=0;

subject to R6b{t in 1..T}:
ceil(KM_0*(1+RG)^t/(KMD*OPD))>=sum{v in 1..Veh, i in I}Z[v,i,t];

subject to R3b{i in I, v in 1..Veh, t1 in 0..T}:
Z[v,i,t1]<=sum{t2 in 0..T:t1>=t2}Y[v,i,t2];
;


subject to R3a{i in I, v in 1..Veh, t1 in 0..T}:
Z[v,i,t1]<=sum{t2 in 1..T+1:t2>t1}X[v,i,t2];
;



subject to R2b{i in I, v in 1..Veh}:
sum{t1 in 1..T+1}t1*X[v,i,t1]>=sum{t2 in 0..T}t2*Y[v,i,t2]
;

subject to R20{t in 1..T, i in I}:
CO2[i,t]=sum{v in 1..Veh}KM[v,i,t]*(1/(EFF[i,TYPE]*FFAC*EVOL[ceil(t/12),i]))*(WTT[i]*BETA[ceil(t/12),i]+TTW[i])/(10^6)
;



subject to R20b{i in I, v in 1..Veh}:
KM[v,i,T]*(1/(EFF[i,TYPE]*FFAC*EVOL[ceil(T/12),i]))*TTW[i]=0
;


subject to R20c:
sum{i in I, v in 1..Veh}KM[v,i,T]*(1/(EFF[i,TYPE]*FFAC*EVOL[ceil(T/12),i]))*TTW[i]=0
;


subject to RNEW{i in {"ELE","HID"}, v in 1..Veh, t in 1..T}:
sum{f in FT}BUY[f,v,i,t]=0;


subject to R22a:
COP21* sum{i1 in I}CO2[i1,0]*(1+RG1)^(96) >= sum{i2 in I}CO2[i2,(96)]
;
