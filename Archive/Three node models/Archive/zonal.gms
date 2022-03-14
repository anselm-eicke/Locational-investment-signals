Sets
all_t       all hours               /1*10/
t(all_t)    hours                   /1*5/
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
Parameter capacity_instr / 0 /;
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
load(t,n)                   hourly load in GWh
avail(t,tec,n)              availability of wind and solar generation (1)
c_var(tec,n)                variable costs (EUR per MWh)
c_fix(tec,n)                annualized fixed costs (EUR per MW p.a.)
supply_slope
supply_slope_GEN(n)
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
supply_slope                = c_fix('wind','south')  * 0.2 / 50;
cap_lim(tec,n)              = 100;
grid_cost(n,m)              = round(B(n,m) * specific_network_costs * sc);

supply_slope_GEN(n)         = 0.01;


Parameters
* Additional parameters for bilevel formulation
GEN_FIX(t,tec,n)
CAP_FIX(tec,n)
LOAD_real_FIX(t,n)
LAMBDA(t)
price(t,n)
;

Free variables
GEN(t,tec,n)
CAP(tec,n)
SHARE(t,con)

WF
WF_inner
FLOW(t,n,m)
THETA(t,n)
;

Positive variables
GRID_CAP(n,m)
LOAD_real(t,n)
UP(t,tec,n)
DOWN(t,tec,n)
RESERVE_CAP(tec,n)
;

Equations
inner_objective,
gen_max, gen_min, gen_share,
cap_max, cap_min, demand_min,
zonal_energy_balance,

outer_objective, 
nodal_energy_balance,
grid_eq1, grid_eq2, grid_eq3, grid_eq4,
redispatch1, redispatch2;

** Inner problem                              
inner_objective..           WF_inner =e= sum((t, n), p_ref * LOAD_real(t,n) * ((1-1/elasticity) + 1/2 * LOAD_real(t,n) / elasticity / load(t,n)))
                                    - sum((tec,n), CAP(tec,n) * c_fix(tec,n) + supply_slope * CAP(tec,n) * CAP(tec,n))
                                    - sum((t,tec,n), GEN(t,tec,n) * c_var(tec,n));

gen_max(t,tec,n)..          0 =g= GEN(t,tec,n) - CAP(tec,n) * avail(t,tec,n);
gen_min(t,tec,n)..          0 =g= -GEN(t,tec,n);
gen_share(t,con,n)..        GEN(t,con,n) =e= CAP(con,n) * SHARE(t,con);
cap_max(tec,n)..            0 =g= CAP(tec,n) - cap_lim(tec,n);
cap_min(tec,n)..            0 =g= -CAP(tec,n);
demand_min(t,n)..           0 =g= -LOAD_real(t,n);

zonal_energy_balance(t)..   0 =e= sum((tec,n),GEN(t,tec,n)) - sum(n,LOAD_real(t,n));

** Outer problem
outer_objective..           WF =e= sum((t, n), p_ref * LOAD_real_FIX(t,n) * ((1-1/elasticity) + 1/2 * LOAD_real_FIX(t,n) / elasticity / load(t,n)))
                                    - sum((tec,n), CAP_FIX(tec,n) * c_fix(tec,n) + supply_slope * CAP_FIX(tec,n) * CAP_FIX(tec,n))
                                    - sum((t,tec,n), GEN_FIX(t,tec,n) * c_var(tec,n))                                    
                                    - sum((n,m),(GRID_CAP(n,m) * grid_cost(n,m)) / 2)
                                    - sum((t,tec,n), UP(t,tec,n) * c_var(tec,n) + DOWN(t,tec,n) * (LAMBDA(t) - c_var(tec,n)))
                                    - sum((tec,n), RESERVE_CAP(tec,n) * (c_fix(tec,n)));
                                    
*nodal energy balance
nodal_energy_balance(t,n).. sum(tec,GEN_FIX(t,tec,n) - DOWN(t,tec,n) + UP(t,tec,n)) - LOAD_real_FIX(t,n) =E= sum(m,FLOW(t,n,m));

*network constraints
grid_eq1(t,n,m)..           FLOW(t,n,m) =l= GRID_CAP(n,m);
grid_eq2(n,m)..             GRID_CAP(n,m) =e= GRID_CAP(m,n);
grid_eq3(t,n,m)..           FLOW(t,n,m) =e= B(n,m) *(THETA(t,n) - THETA(t,m));
grid_eq4(t,n)..             THETA(t,'south') =e= 0;

redispatch1(t,tec,n)..      DOWN(t,tec,n) =L= GEN_FIX(t,tec,n);
redispatch2(t,tec,n)..      UP(t,tec,n) =L= CAP_FIX(tec,n) * avail(t,tec,n) - GEN_FIX(t,tec,n) + RESERVE_CAP(tec,n) * avail(t,tec,n);


Model LOCI_inner /
inner_objective
gen_min
gen_max
gen_share
cap_min
cap_max
demand_min
zonal_energy_balance
/;

* Set starting values
Model LOCI_outer /
outer_objective

nodal_energy_balance

grid_eq1
grid_eq2
grid_eq3
grid_eq4

redispatch1
redispatch2
/;

LOAD_real.L(t,n) =load(t,n);
LOAD_real.lo(t,n) = 0.6 * load(t,n);
*CAP.lo('base','south') = 5;

option NLP=Baron;
Solve LOCI_inner using NLP max WF_inner;

price(t,n) = p_ref * ((1-1/elasticity) + LOAD_real.L(t,n) / elasticity / load(t,n));
GEN_FIX(t,tec,n) = GEN.L(t,tec,n);
CAP_FIX(tec,n) = CAP.L(tec,n);
LOAD_real_FIX(t,n) = LOAD_real.L(t,n);
LAMBDA(t) = zonal_energy_balance.M(t);

Display CAP_FIX, GEN_FIX, LOAD_real_FIX, LAMBDA, price;

Solve LOCI_outer using QCP max WF;

consumer_surplus = sum((t,n), p_ref * LOAD_real.L(t,n) * ((1-1/elasticity) + 1/2 * LOAD_real.L(t,n) / elasticity / load(t,n)));

generation_costs = sum((tec,n), CAP.L(tec,n) * c_fix(tec,n) + supply_slope * CAP.L(tec,n) * CAP.L(tec,n))
                                    + sum((t,tec,n), GEN.L(t,tec,n) * c_var(tec,n));
                                    
total_network_cost = sum((n,m), GRID_CAP.L(n,m) * grid_cost(n,m) / 2) 
                    + sum((t,tec,n), UP.L(t,tec,n) * c_var(tec,n) + DOWN.L(t,tec,n) * (LAMBDA(t) - c_var(tec,n)))
                    + sum((tec,n), RESERVE_CAP.L(tec,n) * (c_fix(tec,n)));

load_deviation(t,n) = LOAD_real.L(t,n) / load(t,n);
res_share = 1 - sum((t,con,n), GEN.L(t,con,n)) / sum((t,tec,n), GEN.L(t,tec,n));

                                    
load_deviation(t,n) = LOAD_real.L(t,n) / load(t,n);

Display WF.L, consumer_surplus, generation_costs, total_network_cost, RES_share;

Display CAP.L, GEN.L, GRID_CAP.L, FLOW.L, THETA.L, load_deviation, res_share, LAMBDA,  LOAD_real.L, UP.L, DOWN.L;




o_WF =WF.L;
o_CS = consumer_surplus;
o_GC = generation_costs;
o_NC = total_network_cost;
o_RES_share = res_share;

o_price(t,n) = LAMBDA(t);
o_cap(tec,n) = CAP.L(tec,n);
o_load(t,n) = LOAD_real.L(t,n);


execute_UNLOAD 'Output/Zonal/out.gdx' o_WF, o_CS, o_GC, o_NC, o_RES_share, o_cap, o_load, o_price, GEN;

* Specify how files are organized in Excel 
$onecho > Output/puts.txt
par=o_WF                            rng=overview!b1        rdim=0 cdim=0
par=o_CS                            rng=overview!b2        rdim=0 cdim=0
par=o_GC                            rng=overview!b3        rdim=0 cdim=0
par=o_NC                            rng=overview!b4        rdim=0 cdim=0
par=o_RES_share                     rng=overview!b5        rdim=0 cdim=0 

par=o_price                         rng=price!a1           rdim=1 cdim=0
par=o_cap                           rng=capacity!a1        rdim=1 cdim=1
par=o_load                          rng=load!a1            rdim=1 cdim=1
$offecho
