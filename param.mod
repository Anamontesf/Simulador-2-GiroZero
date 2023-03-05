param TYPE integer;
param EMP integer;
param T; # periodos (months)
set I; # Number of technologies to be evaluated
set FT; # Finantial tool available.


# 0. Temporales
param DELTA1{r in 1..3} >=0, <=1;
param GDEP1 >=0, <=1;
param WACC1 >=0, <=1;
param r1{f in FT} >=0, <=1;
param RG1 >=0, <=1;

# 1. Non affected by escenarios
param TTW{i in I} >=0;
param WTT{i in I} >=0;
param FM{i in I} >=0 <=2; # maintenance cost factor change among techs
param DELTA{r in 1..3}; # marginal maintenance cost per additional unit age (%/age (month)). 
param GDEP default ((1+GDEP1)^(1/12))-1; # Periodic depreciation rate of asset Price (%/age (month))
param COP21 >=0, <=1;	
param LIM1{i in I} integer >=0;
param OMEGA >=0, <=1; # Effective monthly rate of price associated to management costs 
param DPED integer >=1;	
param EFF{i in I, r in 1..3} >=0; # units of fuel required per kilometer traveled (unit/km) # Type vehicle
param CM_0>=0; 
param TONNES{r in 1..3} >=0;
param NAC binary;  # si regularmente usa el vehiculo para viajes  mayor a 400 km

#2. Financials (external):
param WACC default ((1+WACC1)^(1/12))-1; # Monthly capital cost of company
param n{f in FT} integer >=0; # Number of periods of each financial tool (months)
param P{f in FT,i in 1..3} >=0, <=1; # Initial, payments proportion and final  amount  of asset Price (%) 
param r{f in FT} default ((1+r1[f])^(1/12))-1;  # Interest rate of each financial tool (monthly effective)	
param RG default ((1+RG1)^(1/12))-1; # Growth market rate reflected in traveled km (monthly %)
param IMP >=0 <=1;	
param CF_0{i in I}  >=0; # cost of a unit of fuel (COP/unit)		

# 3. Affected by scenarios: Optimistic, pesimistic
param LAMBDA{t in 0..ceil(T/12),i in I} >=0, <=3; # Change over time of fuel Price (factor)
param GAMMA{t in 0..ceil(T/12),i in I} >=0, <=3; # Change over time of asset Price (factor)
param RESALE >=0, <=1; # Recovery Price percentage after tech lifetime
param EVOL{t in 0..ceil(T/12),i in I} >=0, <=3; # Expected effcience evolution of fuel (factor)
param BETA{t in 0..ceil(T/12),i in I} >=0, <=3;
param SA{t in 0..ceil(T/12),i in I} >=-1, <=1; # asset  price subsidie (%)
param SF{t in 0..ceil(T/12),i in I} >=-1, <=1; # Fuel price subsidie (%)
param YEAR{i in I,s in 1..2} integer; # market year availability 
param INF ;  # Periodic market inflation rate (%/month)
param DEV ; # Peso devaluation							
param TAX{t in 0..ceil(T/12)} >=0; 

# 4. Operatives (internal): Ligeros Medianos Pesados
param PVEH{i in I} >=0; # Asset Price of a new truck of each tech (COP/truck)
param KM_0>=0;
param KMD integer >=0; # average traveled kilometers per day per truck (KM/truck/day) for that type of vehicle
param OPD integer >=0; # Operative monthly days (days/monthly) 
param Veh default ceil((KM_0*(1+RG)^T)/(OPD*KMD))+1;
param OPEXBUDGET>=0; # Operational budget for that type of vehicle current pesos year 0(COP)
param CAPEXBUDGET >=0; # Capital Budget for that type of vehicle current pesos year 0(COP)
param LIM2;
param INV{i in I,v in 1..Veh} binary; # Number of trucks of each tech in current fleet for that type of vehicle
param AGE_0{i in I, v in 1..Veh} integer >=0; # Average age of each tech in current fleet for that type of vehicle
param KMGAL >=0;
param TONTRIP >=0;

#5. Computed
param FFAC default KMGAL/EFF["DIE",TYPE]; # company fuel efficiency factor
param MFAC default TONTRIP/TONNES[TYPE]; # factor de carga
param REL default KM_0/(KMD*OPD*sum{i in I, v in 1..Veh}INV[i,v]) >=0, <=1 ;


