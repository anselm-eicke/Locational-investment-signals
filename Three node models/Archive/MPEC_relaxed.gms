Sets
all_t       all hours               /1*10/
t(all_t)    hours                   /1*3/
tec         generators              /base, peak, wind, solar/
con(tec)    conventional generation /base, peak/
all_n       all buses               /west, east, south/
n(all_n)    selected buses          /west, south/
;

alias (n,m);
alias (all_n,all_m);

* parameters for supply and demand functions
Parameter elasticity / -0.25 /; 
Parameter p_ref / 65 /;
Parameter capacity_instr / 0 /;
Parameter specific_network_costs /255/;
*Source for network costs: EMMA (3400 EUR/MW/km discontiert mit i = 0.07 ueber 40 Jahre)

Table B(all_n,all_m)        Susceptance of transmission lines
         west   east    south
west        1    300      800
east      300      1      500
south     800    500        1
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
supply_slope
cap_lim(tec,n)              capacity limit of generation in each node
grid_cost(n,m)
sc                          scaling factor
    
* Output Parameters
consumer_surplus
generation_costs
total_network_cost
load_deviation(t,n)
res_share
marginal_fixed_costs(tec,n)  
fixed_linear_costs(tec,n)
max_marginal_price_adder(tec,n)
variable_costs(tec,n)
o_WF
o_CS
o_GC
o_NC
o_RES_share
o_load(t,n)
o_cap(tec,n)
o_price(t,n)
o_cap_instr(tec,n)
o_instrument(tec,n)

i_WF
i_instrument(tec,n)
i_cap(tec,n)
i_lambda(t)
i_load_real(t,n)
i_gen(t,tec,n)
i_grid_cap(n,m)
i_flow(t,n,m)
i_up(t,tec,n)
i_down(t,tec,n)
i_mu_G_min(t,tec,n)
i_mu_G_max(t,tec,n)
i_mu_C_min(tec,n)
i_mu_C_max(tec,n)
i_mu_D_min(t,n)
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
supply_slope                = c_fix('wind','south')  * 0.25 / 50;
cap_lim(tec,n)              = 100;
grid_cost(n,m)              = round(B(n,m) * specific_network_costs * sc);




Binary variables y1(t,tec,n),y2(t,tec,n),y3(tec,n),y4(tec,n),y5(t,n),y6(t);

Parameter M1 / 100000/;
Parameter M2 / 100000/;
Parameter M3 / 100000/;
Parameter M4 / 100000/;
Parameter M5 / 100000/;
Parameter M6 / 100000/;

Free variables
GEN(t,tec,n)
CAP(tec,n)
LAMBDA(t)
WF
FLOW(t,n,m)
INSTRUMENT(tec,n)
THETA(t,n)
;


Positive variables
mu_G_min(t,tec,n)
mu_G_max(t,tec,n)
mu_C_min(tec,n)
mu_C_max(tec,n)
mu_D_min(t,n)
GRID_CAP(n,m)
LOAD_real(t,n)
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

gen_min, gen_max, cap_min, cap_max, demand_min,
energy_balance,
KKT_GEN, KKT_CAP, KKT_load

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

complementarity6a
complementarity6b
;


** outer problem
objective..                 WF =e= sum((t, n), p_ref * LOAD_real(t,n) * ((1-1/elasticity) + 1/2 * LOAD_real(t,n) / elasticity / load(t,n)))
                                    - sum((tec,n), CAP(tec,n) * c_fix(tec,n) + supply_slope * CAP(tec,n) * CAP(tec,n))
                                    - sum((t,tec,n), GEN(t,tec,n) * c_var(tec,n) + supply_slope * GEN(t,tec,n) * GEN(t,tec,n))
                                    - sum((n,m),(GRID_CAP(n,m) * grid_cost(n,m)) / 2)
                                    - sum((t,tec,n), UP(t,tec,n) * c_var(tec,n) - DOWN(t,tec,n) * c_var(tec,n))
                                    ;

nodal_energy_balance(t,n).. sum(tec,GEN(t,tec,n) - DOWN(t,tec,n) + UP(t,tec,n)) - LOAD_real(t,n) =E= sum(m,FLOW(t,n,m));

*network constraints
grid_eq1(t,n,m)..           FLOW(t,n,m) =l= GRID_CAP(n,m);
grid_eq2(n,m)..             GRID_CAP(n,m) =e= GRID_CAP(m,n);
grid_eq3(t,n,m)..           FLOW(t,n,m) =e= B(n,m) *(THETA(t,n) - THETA(t,m));
grid_eq4(t,n)..             THETA(t,'south') =e= 0;

redispatch1(t,tec,n)..      DOWN(t,tec,n) =L= GEN(t,tec,n);
redispatch2(t,tec,n)..      UP(t,tec,n) =L= CAP(tec,n) * avail(t,tec,n) - GEN(t,tec,n);


** Inner problem
*Primal constraints
gen_min(t,tec,n)..          0 =g= -GEN(t,tec,n);
gen_max(t,tec,n)..          0 =g= GEN(t,tec,n) - CAP(tec,n) * avail(t,tec,n);

cap_min(tec,n)..            0 =g= -CAP(tec,n);
cap_max(tec,n)..            0 =g= CAP(tec,n) - cap_lim(tec,n);
               
*KKT conditions
KKT_GEN(t,tec,n)..          c_var(tec,n) - p_ref * (1-(1/elasticity)) - p_ref * sum((ttec, nn), GEN(t,ttec,nn)) / (elasticity * sum(nn, load(t,nn))) + mu_G_max(t,tec,n) - mu_G_min(t,tec,n) =e= 0;
KKT_CAP(tec,n)..            c_fix(tec,n) + supply_slope * CAP(tec,n) + INSTRUMENT(tec,n) - sum(t,avail(t,tec,n) * mu_G_max(t,tec,n)) + mu_C_max(tec,n) - mu_C_min(tec,n) =e= 0;

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
nodal_energy_balance

grid_eq1
grid_eq2
grid_eq3
grid_eq4

redispatch1
redispatch2
gen_min
gen_max
cap_min
cap_max
demand_min.mu_D_min
energy_balance.LAMBDA
KKT_GEN
KKT_CAP
KKT_load

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
complementarity6a
complementarity6b
/;

INSTRUMENT.lo(tec,n) = -250;
INSTRUMENT.up(tec,n) = 750;

LAMBDA.lo(t) = -10;
LAMBDA.up(t) = 250;

LOCI.nodlim = 25000000;
LOCI.resLim = 40000


*option MINLP=Baron;
*Solve LOCI using MINLP max WF;

option MIQCP=Cplex;
Solve LOCI using MIQCP max WF;


consumer_surplus = sum((t, n), p_ref * LOAD_real.L(t,n) * ((1-1/elasticity) + 1/2 * LOAD_real.L(t,n) / elasticity / load(t,n))) / sc * 1000;

marginal_fixed_costs(tec,n)          = c_fix(tec,n) / sc * 1000;
fixed_linear_costs(tec,n)   = supply_slope * CAP.L(tec,n) * CAP.L(tec,n) / sc * 1000;
max_marginal_price_adder(tec,n)     = 0.5 * supply_slope * CAP.L(tec,n) / sc * 1000;
variable_costs(tec,n)       = sum(t, GEN.L(t,tec,n) * c_var(tec,n)) / sc * 1000;

generation_costs = (sum((tec,n), CAP.L(tec,n) * c_fix(tec,n) + supply_slope * CAP.L(tec,n) * CAP.L(tec,n))
                                    + sum((t,tec,n), GEN.L(t,tec,n) * c_var(tec,n))) / sc * 1000;
                                    
total_network_cost = (sum((n,m), GRID_CAP.L(n,m) * grid_cost(n,m) / 2) 
                    + sum((t,tec,n), UP.L(t,tec,n) * c_var(tec,n) + DOWN.L(t,tec,n) * (LAMBDA.L(t) - c_var(tec,n))))/sc * 1000;

load_deviation(t,n) = LOAD_real.L(t,n) / load(t,n);
res_share = 1 - sum((t,con,n), GEN.L(t,con,n)) / sum((t,tec,n), GEN.L(t,tec,n));

o_price(t,n) = p_ref * ((1-1/elasticity) + (LOAD_real.L(t,n)/load(t,n) / elasticity ));
                                   
load_deviation(t,n) = LOAD_real.L(t,n) / load(t,n);

o_instrument(tec,n) = round(INSTRUMENT.L(tec,n) / sc) / 1000;
* conversion into €/kW

o_WF =WF.L;
o_CS = consumer_surplus;
o_GC = generation_costs;
o_NC = total_network_cost;
o_cap(tec,n) = CAP.L(tec,n);
o_load(t,n) = LOAD_real.L(t,n);
o_RES_share = res_share;


i_WF = WF.L;
i_instrument(tec,n) = INSTRUMENT.L(tec,n);
i_cap(tec,n) = CAP.L(tec,n);
i_lambda(t) = LAMBDA.L(t);
i_load_real(t,n) = LOAD_real.L(t,n);
i_gen(t,tec,n) = GEN.L(t,tec,n);
i_grid_cap(n,m) = GRID_CAP.L(n,m); 
i_flow(t,n,m) = FLOW.L(t,n,m);
i_up(t,tec,n) = UP.L(t,tec,n);
i_down(t,tec,n) = DOWN.L(t,tec,n);
i_mu_G_min(t,tec,n) = mu_G_min.L(t,tec,n);
i_mu_G_max(t,tec,n) = mu_G_max.L(t,tec,n);
i_mu_C_min(tec,n) = mu_C_min.L(tec,n);
i_mu_C_max(tec,n) = mu_C_max.L(tec,n);
i_mu_D_min(t,n) = mu_D_min.L(t,n);

Display WF.L, consumer_surplus, generation_costs, total_network_cost, RES_share;
Display CAP.L, GEN.L, o_instrument, INSTRUMENT.L, GRID_CAP.L, FLOW.L, LAMBDA.L, o_price, LOAD_real.L, load, load_deviation, UP.L, DOWN.L;


execute_UNLOAD 'in.gdx' i_cost, i_load, i_avail, i_instrument, i_cap, i_lambda, i_load_real, i_gen, i_grid_cap, i_flow, i_mu_G_min, i_mu_G_max, i_mu_C_min, i_mu_C_max, i_mu_D_min,  i_up, i_down
