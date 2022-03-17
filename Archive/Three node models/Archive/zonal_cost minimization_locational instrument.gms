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

Binary variables y1(t,tec,n),y2(t,tec,n),y3(tec,n),y4(tec,n),y5(t,n),y6(t);

Parameter M1 / 100000/;
Parameter M2 / 100000/;
Parameter M3 / 100000/;
Parameter M4 / 100000/;
Parameter M5 / 100000/;
Parameter M6 / 100000/;

Free Variables
GEN(t,tec,n)
CAP(tec,n)
COSTS
FLOW(t,n,m)
INSTRUMENT(tec,n)
THETA(t,n)
LAMBDA(t)
;

Positive variables
mu_G_min(t,tec,n)
mu_G_max(t,tec,n)
mu_C_min(tec,n)
mu_C_max(tec,n)
mu_D_min(t,n)

GRID_CAP(n,m)
UP(t,tec,n)
DOWN(t,tec,n)
;

alias(tec, ttec);
alias(n, nn);


Equations

objective, instr_const, 
nodal_energy_balance,
grid_eq1, grid_eq2, grid_eq3, grid_eq4,
redispatch1, redispatch2,

gen_min, gen_max,
cap_min, cap_max,
demand_min,
energy_balance,

KKT_GEN, KKT_CAP, KKT_load,


complementarity1a,
complementarity1b,
complementarity2a,
complementarity2b,
complementarity3a,
complementarity3b,
complementarity4a,
complementarity4b,
complementarity5a,
complementarity5b,
complementarity6a,
complementarity6b
;

objective..                     COSTS =e= sum((tec,n), CAP(tec,n) * c_fix(tec,n) + 0.5 * CAP(tec,n) * CAP(tec,n) * capacity_slope)
                                    + sum((t,tec,n), GEN(t,tec,n) * c_var(tec,n))
                                    + sum((n,m),(GRID_CAP(n,m) * grid_cost(n,m)) / 2)
                                    + sum((t,tec,n), UP(t,tec,n) * c_var(tec,n) - DOWN(t,tec,n) * c_var(tec,n));
                                

nodal_energy_balance(t,n)..     sum(tec,GEN(t,tec,n) - DOWN(t,tec,n) + UP(t,tec,n)) - load(t,n) =E= sum(m,FLOW(t,n,m));

*network constraints
grid_eq1(t,n,m)..               FLOW(t,n,m) =l= GRID_CAP(n,m);
grid_eq2(n,m)..                 GRID_CAP(n,m) =e= GRID_CAP(m,n);
grid_eq3(t,n,m)..               FLOW(t,n,m) =e= B(n,m) *(THETA(t,n) - THETA(t,m));
grid_eq4(t,n)..                 THETA(t,'south') =e= 0;

redispatch1(t,tec,n)..          DOWN(t,tec,n) =L= GEN(t,tec,n);
redispatch2(t,tec,n)..          UP(t,tec,n) =L= CAP(tec,n) * avail(t,tec,n) - GEN(t,tec,n);

instr_const..                   INSTRUMENT('wind','south') =e= 0;

** INNER PROBLEM

gen_min(t,tec,n)..              0 =g= -GEN(t,tec,n);
gen_max(t,tec,n)..              0 =g= GEN(t,tec,n) - CAP(tec,n) * avail(t,tec,n);

cap_min(tec,n)..                0 =g= -CAP(tec,n);
cap_max(tec,n)..                0 =g= CAP(tec,n) - cap_lim(tec,n);

energy_balance(t)..             0 =e= sum((tec,n),GEN(t,tec,n)) - sum(n,load(t,n));
               
*KKT conditions
KKT_GEN(t,tec,n)..              c_var(tec,n) + mu_G_max(t,tec,n) - mu_G_min(t,tec,n) - LAMBDA(t) =e= 0;
KKT_CAP(tec,n)..                c_fix(tec,n) + capacity_slope * CAP(tec,n) - sum(t,avail(t,tec,n) * mu_G_max(t,tec,n)) + mu_C_max(tec,n) - mu_C_min(tec,n) =e= 0;

*+ INSTRUMENT(tec,n)
complementarity1a(t,tec,n)..    GEN(t,tec,n)        =L= y1(t,tec,n) * M1;
complementarity1b(t,tec,n)..    mu_G_min(t,tec,n)   =L= (1-y1(t,tec,n)) * M1;
complementarity2a(t,tec,n)..    CAP(tec,n) * avail(t,tec,n) - GEN(t,tec,n) =L= y2(t,tec,n) * M2;
complementarity2b(t,tec,n)..    mu_G_max(t,tec,n)   =L= (1-y2(t,tec,n)) * M2;
complementarity3a(tec,n)..      CAP(tec,n)          =L= y3(tec,n) * M3;
complementarity3b(tec,n)..      mu_C_min(tec,n)     =L= (1-y3(tec,n)) * M3;
complementarity4a(tec,n)..      cap_lim(tec,n) - CAP(tec,n) =L= y4(tec,n) * M4;
complementarity4b(tec,n)..      mu_C_max(tec,n)     =L= (1-y4(tec,n)) * M4;

complementarity6a(t)..          sum(n,load(t,n)) - sum((tec,n),GEN(t,tec,n)) =L= y6(t) * M6;
complementarity6b(t)..          LAMBDA(t)           =L= (1-y6(t)) * M6;

Model LOCI /


objective
nodal_energy_balance

grid_eq1
grid_eq2
grid_eq3
grid_eq4

instr_const

redispatch1
redispatch2

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
complementarity6a
complementarity6b
/;


GEN.up(t,tec,n) = 100;
GEN.lo(t,tec,n) = 0;

DOWN.up(t,tec,n) = 100;
DOWN.lo(t,tec,n) = 0;

UP.up(t,tec,n) = 100;
UP.lo(t,tec,n) = 0;

LOCI.nodlim = 25000000;
LOCI.resLim = 40000;

*Option optcr = 0.02;

Option MIQCP = Cplex;


Solve LOCI minimizing COSTS using MIQCP;


*price(t) = p_ref * (1-(1/elasticity) + (sum((tec,n), GEN.L(t,tec,n)) / sum(n, elasticity * load(t,n))));
*load_deviation(t) = sum((tec,n), GEN.L(t,tec,n)) / sum(n,load(t,n));
*i_instrument(tec,n) = INSTRUMENT.L(tec,n) / sc / 1000;

network_cost = (sum((n,m),(GRID_CAP.L(n,m) * grid_cost(n,m)) / 2) + sum((t,tec,n), UP.L(t,tec,n) * c_var(tec,n) - DOWN.L(t,tec,n) * c_var(tec,n))) / sc;
*consumer_surplus = sum(t, p_ref * sum((n), LOAD_real.L(t,n)) * (1-1/elasticity + sum((n), LOAD_real.L(t,n)) / (2*elasticity* sum(n,load(t,n))))) / sc;

generation_costs = (sum((tec,n), CAP.L(tec,n) * c_fix(tec,n) + 0.5 * CAP.L(tec,n) * CAP.L(tec,n) * capacity_slope) + sum((t,tec,n), GEN.L(t,tec,n) * c_var(tec,n))) / sc;
                                                                    
*sum_instrument = sum((tec,n), INSTRUMENT.L(tec,n) * CAP.L(tec,n)) / sc;

*load_deviation(t) = sum(n,LOAD_real.L(t,n)) / sum(n,load(t,n));
res_share = 1 - sum((t,con,n), GEN.L(t,con,n)) / sum((t,tec,n), GEN.L(t,tec,n));
o_cap(tec,n) = CAP.L(tec,n);
o_gen(t,tec,n) = GEN.L(t,tec,n);


Display GEN.L, CAP.L, network_cost, GRID_CAP.L;

execute_UNLOAD 'Output/with_instrument.gdx' generation_costs, network_cost, res_share, o_cap, o_gen, price, c_fix;