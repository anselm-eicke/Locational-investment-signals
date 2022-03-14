Sets
all_t       all hours               /1*10/
t(all_t)    hours                   /1*10/
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
load_deviation(t,n)
res_share

o_WF
o_CS
o_GC
o_NC
o_RES_share

o_cap(tec,n)
o_price(t,n)
o_cap_instr(tec,n)
load_deviation(t,n)
;

Variables
GEN(t,tec,n)
CAP(tec,n)
PRICE(t)
PROFIT
;


Equations
profit_eq,
price_eq,

gen_min, gen_max,
cap_min, cap_max;

** Inner problem
*Primal constraints

alias(tec, ttec);
alias(n, nn);

profit_eq..                 PROFIT =e= sum((t,tec,n), GEN(t,tec,n) * PRICE(t) - c_var(tec,n))
                                        - sum((tec,n), CAP(tec,n) * c_fix(tec,n) (c_fix(tec,n) + CAP(tec,n) * supply_slope));


price_eq(t)..               PRICE(t) =e= p_ref * (1-(1/elasticity) + (sum((tec,n), GEN(t,tec,n)) / sum(n, elasticity * load(t,n))));

gen_min(t,tec,n)..          0 =g= -GEN(t,tec,n);
gen_max(t,tec,n)..          0 =g= GEN(t,tec,n) - CAP(tec,n) * avail(t,tec,n);
cap_min(tec,n)..            0 =g= -CAP(tec,n);
cap_max(tec,n)..            0 =g= CAP(tec,n) - cap_lim(tec,n);

Model LOCI /

gen_min
gen_max
cap_min
cap_max

profit_eq
price_eq
/;

* Set starting values
*LOAD_real.L(t,n) =load(t,n);

LOCI.optCR = 0.00001

option QCP=Cplex;

Solve LOCI max PROFIT using NLP;

display CAP.L, PRICE.L, GEN.L