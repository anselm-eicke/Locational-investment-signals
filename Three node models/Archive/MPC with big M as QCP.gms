Sets
all_t       all hours               /1*10/
t(all_t)    hours                   /1*2/
tec         generators              /base, peak, wind, solar/
con(tec)    conventional generation /base, peak/
all_n       all buses               /west, east, south/
n(all_n)    selected buses          /west, east, south/
;

alias (n,m);
alias (all_n,all_m);

* parameters for supply and demand functions
Parameter elasticity / -0.25 /; 
Parameter p_ref / 65 /;
Parameter capacity_instr / 1 /;
Parameter specific_network_costs /150/;
*Source for network costs: EMMA (3400 EUR/MW/km discontiert mit i = 0.07 ueber 40 Jahre)

Table B(all_n,all_m)        Susceptance of transmission lines
         west   east    south
west        1    250      600
east      250      1      400
south     600    400        1
;

Parameters
* Input Parameters
i_cost(*,*)                 cost data to be loaded from sheet "cost"
i_load(all_t,all_n)         load data to be loaded from sheet "time series" in MWh
i_avail(all_t,tec,all_n)    availability data
* Model Parameters
load(t,n)               hourly load (MWh)
avail(t,tec,n)          availability of wind and solar generation (1)
c_var(tec,n)            variable costs (EUR per MWh)
c_fix(tec,n)            annualized fixed costs (EUR per MW p.a.)
supply_slope
cap_lim(tec,n)          capacity limit of generation in each node
grid_cost
INSTRUMENT(tec,n)
;

*-----------------------------------
* LOAD INPUT PARAMETERS FROM EXCEL


* Load GDX file to GAMS
$GDXIN "in.gdx"
$LOADdc i_cost, i_load, i_avail


*-----------------------------------
* MODEL PARAMETERS
* Time series parameters (loaded)
load(t,n)                = i_load(t,n) / 1000;
avail(t,tec,n)           = i_avail(t,tec,n);
avail(t,con,n)           = 1;

* Variable cost parameters (loaded)
c_var(tec, n)         = i_cost(tec,"cost_var");

* Fixed cost parameters are loaded below in the loop
c_fix(tec, n)         = round(i_cost(tec,"cost_fix") * 1000 * card(t) / 8760);

cap_lim(tec,n)        = 100000;

grid_cost = round(166000 * card(t) / 8760);

supply_slope          = c_fix('wind','east')  * 0.2 / 50;
INSTRUMENT(tec,n)     = 0;
*INSTRUMENT(tec,'west')     = 10;

display c_var, load, avail, c_fix;

Parameters
consumer_surplus
generation_costs
total_network_cost
res_share

o_WF
o_CS
o_GC
o_NC
o_RES_share

o_cap(tec,n)
price(t)
load_deviation(t)
;



Binary variables y1(t,tec,n),y2(t,tec,n),y3(tec,n),y4(tec,n),y5(t,n),y6(t);

Parameter M1 / 1000000/;
Parameter M2 / 1000000/;
Parameter M3 / 1000000/;
Parameter M4 / 1000000/;
Parameter M5 / 1000000/;
Parameter M6 / 1000000/;

Variables
GEN(t,tec,n)
CAP(tec,n)
WF;
;

Positive variables
mu_G_min(t,tec,n)
mu_G_max(t,tec,n)
mu_C_min(tec,n)
mu_C_max(tec,n)
;

alias(tec, ttec);
alias(n, nn);


Equations
gen_min, gen_max,
cap_min, cap_max, 
KKT_GEN, KKT_CAP,
objective,

complementarity1a,
complementarity1b,
complementarity2a,
complementarity2b,
complementarity3a,
complementarity3b,
complementarity4a,
complementarity4b
;

objective..                     WF =e= 1;

gen_min(t,tec,n)..              0 =g= -GEN(t,tec,n);
gen_max(t,tec,n)..              0 =g= GEN(t,tec,n) - CAP(tec,n) * avail(t,tec,n);

cap_min(tec,n)..                0 =g= -CAP(tec,n);
cap_max(tec,n)..                0 =g= CAP(tec,n) - cap_lim(tec,n);
               
*KKT conditions
KKT_GEN(t,tec,n)..              c_var(tec,n) - p_ref * (1-(1/elasticity)) - p_ref * sum((ttec, nn), GEN(t,ttec,nn)) / (elasticity * sum(nn, load(t,nn))) + mu_G_max(t,tec,n) - mu_G_min(t,tec,n) =e= 0;
KKT_CAP(tec,n)..                c_fix(tec,n) + supply_slope * CAP(tec,n) + INSTRUMENT(tec,n) - sum(t,avail(t,tec,n) * mu_G_max(t,tec,n)) + mu_C_max(tec,n) - mu_C_min(tec,n) =e= 0;

complementarity1a(t,tec,n)..    GEN(t,tec,n)        =L= y1(t,tec,n) * M1;
complementarity1b(t,tec,n)..    mu_G_min(t,tec,n)   =L= (1-y1(t,tec,n)) * M1;
complementarity2a(t,tec,n)..    CAP(tec,n) * avail(t,tec,n) - GEN(t,tec,n) =L= y2(t,tec,n) * M2;
complementarity2b(t,tec,n)..    mu_G_max(t,tec,n)   =L= (1-y2(t,tec,n)) * M2;
complementarity3a(tec,n)..      CAP(tec,n)          =L= y3(tec,n) * M3;
complementarity3b(tec,n)..      mu_C_min(tec,n)     =L= (1-y3(tec,n)) * M3;
complementarity4a(tec,n)..      cap_lim(tec,n) - CAP(tec,n) =L= y4(tec,n) * M4;
complementarity4b(tec,n)..      mu_C_max(tec,n)     =L= (1-y4(tec,n)) * M4;

Model LOCI /

objective

gen_min
gen_max
cap_min
cap_max

KKT_GEN
KKT_CAP

complementarity1a
complementarity1b
complementarity2a
complementarity2b
complementarity3a
complementarity3b
complementarity4a
complementarity4b
/;

* Set starting values
*LOAD_real.L(t,n) =load(t,n);


Solve LOCI maximizing WF using MIP;

price(t) = p_ref * (1-(1/elasticity) + (sum((tec,n), GEN.L(t,tec,n)) / sum(n, elasticity * load(t,n))));
load_deviation(t) = sum((tec,n), GEN.L(t,tec,n)) / sum(n,load(t,n));

Display GEN.L, CAP.L, price, load_deviation;
*,WF.L,  GRID.L, LAMBDA.L, LOAD_real.L, CAP_INSTR.L, UP.L, DOWN.L, RESERVE_CAP.L, FLOW.L;
