Sets
all_t       all hours               /1*10/
t(all_t)    hours                   /1*2/
tec         generators              /base, peak, wind, solar/
con(tec)    conventional generation /base, peak/
all_n       all buses               /west, east, south/
n(all_n)    selected buses          /west, east/
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

avail(t,tec,n) $ (avail(t,tec,n) = 0)  =  0.0001 ;

* Variable cost parameters (loaded)
c_var(tec, n)         = i_cost(tec,"cost_var");

* Fixed cost parameters are loaded below in the loop
c_fix(tec, n)         = round(i_cost(tec,"cost_fix") * 1000 * card(t) / 8760);

cap_lim(tec,n)        = 100000;

grid_cost = round(166000 * card(t) / 8760);

supply_slope          = c_fix('wind','east')  * 0.2 / 50;
INSTRUMENT(tec,n)     = 0;
*INSTRUMENT(tec,'west')     = 0.001;

display c_var, load, avail, c_fix;


Parameters
consumer_surplus
generation_costs
total_network_cost
load_deviation(t,n)
res_share

o_WF
o_CS
o_GC
o_NC
o_RES_share

o_cap(tec,n)
o_price(t)
o_cap_instr(tec,n)
load_deviation(t,n)
;


Variables
LOAD_real(t,n)
GEN(t,tec,n)
CAP(tec,n)
RHO(t,tec,n)
WF

profit
;

Positive variables
mu_G_min(t,tec,n)
mu_G_max(t,tec,n)
mu_C_min(tec,n)
mu_C_max(tec,n)
mu_D_min(t,n)

mu_S_min(t,tec,n)
mu_S_max(t,tec,n)
mu_k_min(t,tec)
mu_k_max(t,tec)

UP(t,tec,n)
DOWN(t,tec,n)

CONST(t,tec)
SHARE(t,tec,n)

;


Equations
Objective, gen_min, gen_max, cap_min, cap_max, demand_min,

share_min, share_max, const_min, const_max,

energy_balance, share_def, tie_break_eq;

Objective..             profit =e= sum((t,n), p_ref * LOAD_real(t,n) * (1-1/elasticity +  LOAD_real(t,n) / (2*elasticity* load(t,n))))
                                - sum((tec,n), CAP(tec,n)*(c_fix(tec,n) + CAP(tec,n)* supply_slope + INSTRUMENT(tec,n)))
                                - sum((t,tec,n), GEN(t,tec,n) * c_var(tec,n));
                
gen_max(t,tec,n)..      0 =g= GEN(t,tec,n) - CAP(tec,n) * avail(t,tec,n);
*gen_max(t,tec,n)..      0 =g= SHARE(t,tec,n) - 1;
gen_min(t,tec,n)..      0 =g= -GEN(t,tec,n);
cap_max(tec,n)..        0 =g= CAP(tec,n) - cap_lim(tec,n);
cap_min(tec,n)..        0 =g= -CAP(tec,n);
demand_min(t,n)..       0 =g= -LOAD_real(t,n);

share_min(t,tec,n)..        0 =g= -SHARE(t,tec,n);
share_max(t,tec,n)..        0 =g= SHARE(t,tec,n) - 1;
const_min(t,tec)..          0 =g= -CONST(t,tec);
const_max(t,tec)..          0 =g= CONST(t,tec) - 1;

energy_balance(t)..     0 =e= sum((tec,n),GEN(t,tec,n)) - sum(n,LOAD_real(t,n));

share_def(t,tec,n)..    0 =e= CAP(tec,n) * avail(t,tec,n) * SHARE(t,tec,n) - GEN(t,tec,n);

tie_break_eq(t,tec,n).. 0 =e= SHARE(t,tec,n) - CONST(t,tec);

Model LOCI /

objective
gen_min
gen_max
cap_min
cap_max
demand_min

share_min, share_max, const_min, const_max,

energy_balance
tie_break_eq
share_def
/;

option NLP = Baron;
option QCP = Cplex;

LOCI.optCR = 0.0000001
*LOCI.optca = 0.0000001
Solve LOCI maximizing profit using NLP;



o_price(t) = p_ref * (1 - 1/elasticity + (LOAD_real.L(t,'west')/ load(t,'west') / elasticity));   

load_deviation(t,n) = LOAD_real.L(t,n) / load(t,n);

Display profit.L, SHARE.L, GEN.L, CAP.L, LOAD_real.L, load_deviation;
