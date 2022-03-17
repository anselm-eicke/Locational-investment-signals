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

Parameter M1 / 100000/;
Parameter M2 / 100000/;
Parameter M3 / 100000/;
Parameter M4 / 100000/;
Parameter M5 / 100000/;

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


Model LOCI /

objective

gen_min.mu_G_min
gen_max.mu_G_max
cap_min.mu_C_min
cap_max.mu_C_max

energy_balance.LAMBDA

KKT_GEN
KKT_CAP
/;


LOCI.nodlim = 25000000;
LOCI.resLim = 40000;

Solve LOCI using mcp;

Display GEN.L, CAP.L;
*,WF.L,  GRID.L, LAMBDA.L, LOAD_real.L, CAP_INSTR.L, UP.L, DOWN.L, RESERVE_CAP.L, FLOW.L;
