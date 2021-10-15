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
*, i_instrument, i_cap, i_lambda, i_load_real, i_gen, i_grid_cap, i_flow, i_mu_G_min, i_mu_G_max, i_mu_C_min, i_mu_C_max, i_mu_D_min, i_up, i_down


* Data assignment
sc = card(t) / 8760;
load(t,n)                   = i_load(t,n) / 1000;
avail(t,tec,n)              = i_avail(t,tec,n);
avail(t,con,n)              = 1;
c_var(tec, n)               = i_cost(tec,"cost_var");
c_fix(tec, n)               = round(i_cost(tec,"cost_fix") * 1000 * sc);
supply_slope                = c_fix('wind','east')  * 0.2 / 50;
cap_lim(tec,n)              = 100;
grid_cost(n,m)              = round(B(n,m) * specific_network_costs * sc);




Free variables
GEN(t,tec,n)
CAP(tec,n)
LAMBDA(t)
RHO(t,con,n)
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
SHARE(t,con)
;

Equations
objective, instr_const, 
nodal_energy_balance,
grid_eq1, grid_eq2, grid_eq3, grid_eq4,
redispatch1, redispatch2,
gen_min, gen_max, cap_min, cap_max, demand_min,
zonal_energy_balance,
KKT_GEN, KKT_CAP, KKT_load;

** Outer problem
objective..                 WF =e= sum((t, n), p_ref * LOAD_real(t,n) * ((1-1/elasticity) + 1/2 * LOAD_real(t,n) / elasticity / load(t,n)))
                                    - sum((tec,n), CAP(tec,n) * c_fix(tec,n) + supply_slope * CAP(tec,n) * CAP(tec,n))
                                    - sum((t,tec,n), GEN(t,tec,n) * c_var(tec,n))
                                    - sum((n,m),(GRID_CAP(n,m) * grid_cost(n,m)) / 2)
                                    - sum((t,tec,n), UP(t,tec,n) * c_var(tec,n) + DOWN(t,tec,n) * (LAMBDA(t) - c_var(tec,n)))
                                    ;
                                    
instr_const..               INSTRUMENT('solar','east') =e= 0;

nodal_energy_balance(t,n).. sum(tec,GEN(t,tec,n) - DOWN(t,tec,n) + UP(t,tec,n)) - LOAD_real(t,n) =E= sum(m,FLOW(t,n,m));

*network constraints
grid_eq1(t,n,m)..           FLOW(t,n,m) =l= GRID_CAP(n,m);
grid_eq2(n,m)..             GRID_CAP(n,m) =e= GRID_CAP(m,n);
grid_eq3(t,n,m)..           FLOW(t,n,m) =e= B(n,m) * (THETA(t,n) - THETA(t,m));
grid_eq4(t,n)..             THETA(t,'east') =e= 0;

redispatch1(t,tec,n)..      DOWN(t,tec,n) =L= GEN(t,tec,n);
redispatch2(t,tec,n)..      UP(t,tec,n) =L= CAP(tec,n) * avail(t,tec,n) - GEN(t,tec,n);

** Inner problem
*Primal constraints
gen_max(t,tec,n)..          0 =g= GEN(t,tec,n) - CAP(tec,n) * avail(t,tec,n);
gen_min(t,tec,n)..          0 =g= -GEN(t,tec,n);
cap_max(tec,n)..            0 =g= CAP(tec,n) - cap_lim(tec,n);
cap_min(tec,n)..            0 =g= -CAP(tec,n);
demand_min(t,n)..           0 =g= -LOAD_real(t,n);

zonal_energy_balance(t)..   0 =e= sum((tec,n),GEN(t,tec,n)) - sum(n,LOAD_real(t,n));
               
*KKT conditions
KKT_GEN(t,tec,n)..          c_var(tec,n) + mu_G_max(t,tec,n) - mu_G_min(t,tec,n) - LAMBDA(t) =e= 0;
KKT_CAP(tec,n)..            c_fix(tec,n) + 2 * supply_slope * CAP(tec,n) + INSTRUMENT(tec,n) - sum(t,avail(t,tec,n) * mu_G_max(t,tec,n)) + mu_C_max(tec,n) - mu_C_min(tec,n) =e= 0;
KKT_load(t,n)..             -p_ref* ((1-1/elasticity) + LOAD_real(t,n) / (elasticity * load(t,n))) - mu_D_min(t,n) + LAMBDA(t) =e= 0;                

Model LOCI /
objective
*instr_const
nodal_energy_balance

grid_eq1
grid_eq2
grid_eq3
grid_eq4

redispatch1
redispatch2

gen_min.mu_G_min
gen_max.mu_G_max
cap_min.mu_C_min
cap_max.mu_C_max
demand_min.mu_D_min

zonal_energy_balance.LAMBDA

KKT_GEN
KKT_CAP
KKT_load

/;

* Set starting values

*INSTRUMENT.L(tec,n) = i_instrument(tec,n);
*CAP.L(tec,n) = i_cap(tec,n);
*LAMBDA.L(t) = i_lambda(t);
*LOAD_real.L(t,n) = i_load_real(t,n);
*GEN.L(t,tec,n) = i_gen(t,tec,n);
*GRID_CAP.L(n,m) = i_grid_cap(n,m);
*FLOW.L(t,n,m) = i_flow(t,n,m);

*UP.L(t,tec,n) = i_up(t,tec,n);
*DOWN.L(t,tec,n) = i_down(t,tec,n);

*mu_G_min.L(t,tec,n) = i_mu_G_min(t,tec,n); 
*mu_G_max.L(t,tec,n) = i_mu_G_max(t,tec,n);
*mu_C_min.L(tec,n) = i_mu_C_min(tec,n);
*mu_C_max.L(tec,n) = i_mu_C_max(tec,n);
*mu_D_min.L(t,n) = i_mu_D_min(t,n);

INSTRUMENT.lo(tec,n) = -250;
INSTRUMENT.up(tec,n) = 750;

LAMBDA.lo(t) = -20;
LAMBDA.up(t) = 350;

LOCI.nodlim = 5000000;
LOCI.resLim = 40000

option NLP=Baron;
Solve LOCI using NLP max WF;


consumer_surplus = sum((t, n), p_ref * LOAD_real.L(t,n) * ((1-1/elasticity) + 1/2 * LOAD_real.L(t,n) / elasticity / load(t,n)));

generation_costs = sum((tec,n), CAP.L(tec,n) * c_fix(tec,n) + supply_slope * CAP.L(tec,n) * CAP.L(tec,n))
                                    + sum((t,tec,n), GEN.L(t,tec,n) * c_var(tec,n));
                                    
total_network_cost = sum((n,m), GRID_CAP.L(n,m) * grid_cost(n,m) / 2) 
                    + sum((t,tec,n), UP.L(t,tec,n) * c_var(tec,n) + DOWN.L(t,tec,n) * (LAMBDA.L(t) - c_var(tec,n)));

load_deviation(t,n) = LOAD_real.L(t,n) / load(t,n);
res_share = 1 - sum((t,con,n), GEN.L(t,con,n)) / sum((t,tec,n), GEN.L(t,tec,n));
                                   
load_deviation(t,n) = LOAD_real.L(t,n) / load(t,n);

o_instrument(tec,n) = round(INSTRUMENT.L(tec,n) / card(t) * 8760) / 1000;

o_WF =WF.L;
o_CS = consumer_surplus;
o_GC = generation_costs;
o_NC = total_network_cost;
o_RES_share = res_share;

o_price(t,n) = LAMBDA.L(t);
o_cap(tec,n) = CAP.L(tec,n);
o_load(t,n) = LOAD_real.L(t,n);

Display WF.L, consumer_surplus, generation_costs, total_network_cost, RES_share;
Display CAP.L, GEN.L, o_instrument, INSTRUMENT.L, GRID_CAP.L, FLOW.L, LAMBDA.L, o_price, LOAD_real.L, load, load_deviation, UP.L, DOWN.L;