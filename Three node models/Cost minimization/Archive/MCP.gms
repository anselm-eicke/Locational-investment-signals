Sets
all_t       all hours               /1*48/
t(all_t)    hours                   /1*24/
tec         generators              /base, peak, wind, solar/
con(tec)    conventional generation /base, peak/
all_n       all buses               /north, south, total/
n(all_n)    selected buses          /north, south/
;

alias (n,m);
alias (all_n,all_m);

Parameter specific_network_costs /300/;
*Source for network costs: EMMA (3400 EUR/MW/km discontiert mit i = 0.07 ueber 40 Jahre)

Table B(all_n,all_m)        Susceptance of transmission lines
         north  south   total
north        1     700    250
south      700       1    500
total      250     500      1
;

Parameters
* Input Parameters
i_cost(*,*)                 cost data to be loaded from sheet "cost"
i_load(all_t,all_n)         load data to be loaded from sheet "time series" in MWh
i_avail(all_t,tec,all_n)    availability data

* Model Parameters
load(t,n)                   hourly load in GWh
avail(t,tec,n)              availability of wind and solar generation (1)
c_var(tec,n)                variable costs (EUR per MWh)
c_fix(tec,n)                annualized fixed costs (EUR per MW p.a.)
capacity_slope
cap_lim(tec,n)              capacity limit of generation in each node
grid_cost(n,m)
sc                          scaling factor
    
* Output Parameters
generation_costs
total_network_cost
res_share

o_GC
o_NC
o_RES_share
o_cap(tec,n)
o_gen(t,tec,n)

o_instrument(tec,n)
sum_instrument
network_cost

i_instrument(tec,n)
i_cap(tec,n)
i_lambda(t)
i_gen(t,tec,n)
i_grid_cap(n,m)
i_flow(t,n,m)
i_up(t,tec,n)
i_down(t,tec,n)
i_mu_G_min(t,tec,n)
i_mu_G_max(t,tec,n)
i_mu_C_min(tec,n)
i_mu_C_max(tec,n)
;

* Load data
$GDXIN "in.gdx"
$LOADdc i_cost, i_load, i_avail

* Data assignment
sc = card(t) / 8760;
load(t,n)                   = i_load(t,n) / 1000;
avail(t,tec,n)              = i_avail(t,tec,n);
avail(t,con,n)              = 1;
c_var(tec, n)               = i_cost(tec,"cost_var");
c_fix(tec, n)               = round(i_cost(tec,"cost_fix") * 1000 * sc);
capacity_slope              = 0.15;
cap_lim(tec,n)              = 100;
grid_cost(n,m)              = round(B(n,m) * specific_network_costs * sc);

display c_var, load, avail, c_fix;

Binary variables y1(t,tec,n),y2(t,tec,n),y3(tec,n),y4(tec,n),y5(t);

Parameter M1 / 10000/;
Parameter M2 / 10000/;
Parameter M3 / 10000/;
Parameter M4 / 10000/;
Parameter M5 / 10000/;

Free Variables
GEN(t,tec,n)
CAP(tec,n)
LAMBDA(t)
COSTS
;

Positive variables
mu_G_min(t,tec,n)
mu_G_max(t,tec,n)
mu_C_min(tec,n)
mu_C_max(tec,n)
;


Equations
objective,
gen_min, gen_max,
cap_min, cap_max,
energy_balance,

KKT_GEN, KKT_CAP,

complementarity1a,
complementarity1b,
complementarity2a,
complementarity2b,
complementarity3a,
complementarity3b,
complementarity4a,
complementarity4b,
complementarity5a,
complementarity5b
;

objective..                     COSTS =e= 1;

gen_min(t,tec,n)..              0 =g= -GEN(t,tec,n);
gen_max(t,tec,n)..              0 =g= GEN(t,tec,n) - CAP(tec,n) * avail(t,tec,n);

cap_min(tec,n)..                0 =g= -CAP(tec,n);
cap_max(tec,n)..                0 =g= CAP(tec,n) - cap_lim(tec,n);

energy_balance(t)..             0 =e= sum((tec,n),GEN(t,tec,n)) - sum(n,load(t,n));
               
*KKT conditions
KKT_GEN(t,tec,n)..              c_var(tec,n) + mu_G_max(t,tec,n) - mu_G_min(t,tec,n) - LAMBDA(t) =e= 0;
KKT_CAP(tec,n)..                c_fix(tec,n) + capacity_slope * CAP(tec,n) - sum(t,avail(t,tec,n) * mu_G_max(t,tec,n)) + mu_C_max(tec,n) - mu_C_min(tec,n) =e= 0;

complementarity1a(t,tec,n)..    GEN(t,tec,n)        =L= y1(t,tec,n) * M1;
complementarity1b(t,tec,n)..    mu_G_min(t,tec,n)   =L= (1-y1(t,tec,n)) * M1;
complementarity2a(t,tec,n)..    CAP(tec,n) * avail(t,tec,n) - GEN(t,tec,n) =L= y2(t,tec,n) * M2;
complementarity2b(t,tec,n)..    mu_G_max(t,tec,n)   =L= (1-y2(t,tec,n)) * M2;
complementarity3a(tec,n)..      CAP(tec,n)          =L= y3(tec,n) * M3;
complementarity3b(tec,n)..      mu_C_min(tec,n)     =L= (1-y3(tec,n)) * M3;
complementarity4a(tec,n)..      cap_lim(tec,n) - CAP(tec,n) =L= y4(tec,n) * M4;
complementarity4b(tec,n)..      mu_C_max(tec,n)     =L= (1-y4(tec,n)) * M4;
complementarity5a(t)..          sum(n,load(t,n)) - sum((tec,n),GEN(t,tec,n)) =L= y5(t) * M5;
complementarity5b(t)..          LAMBDA(t)           =L= (1-y5(t)) * M5;

Model LOCI /

objective

gen_min
gen_max
cap_min
cap_max

energy_balance

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
complementarity5a
complementarity5b
/;

GEN.up(t,tec,n) = 100;
GEN.lo(t,tec,n) = 0;

CAP.up(tec,n) = 100;
CAP.lo(tec,n) = 0;

LOCI.nodlim = 25000000;
LOCI.resLim = 40000;

Solve LOCI minimizing COSTS using MIP;

Display GEN.L, CAP.L;
*,WF.L,  GRID.L, LAMBDA.L, LOAD_real.L, CAP_INSTR.L, UP.L, DOWN.L, RESERVE_CAP.L, FLOW.L;
