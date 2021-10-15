Sets
all_t       all hours               /1*10/
t(all_t)    hours                   /1*3/
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







Free variables
LAMBDA(t)
WF
FLOW(t,n,m)
THETA(t,n)

;

Positive variables
LOAD_real(t,n)
GRID_CAP(n,m)
GEN(t,tec,n)
CAP(tec,n)
;

Equations
objective, cap_constraint, cap_limit,

nodal_energy_balance,
grid_eq1, grid_eq2, grid_eq3, grid_eq4;

objective..                  WF =e= sum((t, n), p_ref * LOAD_real(t,n) * ((1-1/elasticity) + 1/2 * LOAD_real(t,n) / elasticity / load(t,n)))
                                    - sum((tec,n), CAP(tec,n) * c_fix(tec,n) + supply_slope * CAP(tec,n) * CAP(tec,n))
                                    - sum((t,tec,n), GEN(t,tec,n) * c_var(tec,n))
                                    - sum((n,m),(GRID_CAP(n,m) * grid_cost(n,m)) / 2);
                                    
cap_constraint(t,tec,n)..   GEN(t,tec,n) =L= CAP(tec,n) * avail(t,tec,n);
cap_limit(tec,n)..          CAP(tec,n) =L= cap_lim(tec,n);

*nodal energy balance
nodal_energy_balance(t,n).. sum(tec,GEN(t,tec,n)) - LOAD_real(t,n) =E= sum(m,FLOW(t,n,m));

*network constraints
grid_eq1(t,n,m)..           FLOW(t,n,m) =l= GRID_CAP(n,m);
grid_eq2(n,m)..             GRID_CAP(n,m) =e= GRID_CAP(m,n);
grid_eq3(t,n,m)..           FLOW(t,n,m) =e= B(n,m) *(THETA(t,n) - THETA(t,m));
grid_eq4(t,n)..             THETA(t,'south') =e= 0;
              
      
Model nodal /
objective

cap_constraint
cap_limit

nodal_energy_balance

grid_eq1
grid_eq2
grid_eq3
grid_eq4
/;

* Set starting point
LOAD_real.L(t,n) =load(t,n);

option QCP=Cplex;
Solve nodal using QCP max WF;


consumer_surplus = sum((t, n), p_ref * LOAD_real.L(t,n) * ((1-1/elasticity) + 1/2 * LOAD_real.L(t,n) / elasticity / load(t,n)));

generation_costs = sum((tec,n), CAP.L(tec,n) * c_fix(tec,n) + supply_slope * CAP.L(tec,n) * CAP.L(tec,n))
                                    + sum((t,tec,n), GEN.L(t,tec,n) * c_var(tec,n));

total_network_cost = sum((n,m),GRID_CAP.L(n,m) * grid_cost(n,m));

load_deviation(t,n) = LOAD_real.L(t,n) / load(t,n);
res_share = 1 - sum((t,con,n), GEN.L(t,con,n)) / sum((t,tec,n), GEN.L(t,tec,n));

o_price(t,n) = p_ref * ((1-1/elasticity) + (LOAD_real.L(t,n)/load(t,n) / elasticity ));

marginal_fixed_costs(tec,n)          = c_fix(tec,n) / sc * 1000;
fixed_linear_costs(tec,n)   = supply_slope * CAP.L(tec,n) * CAP.L(tec,n) / sc * 1000;
max_marginal_price_adder(tec,n)     = 0.5 * supply_slope * CAP.L(tec,n) / sc * 1000;
variable_costs(tec,n)       = sum(t, GEN.L(t,tec,n) * c_var(tec,n)) / sc * 1000;

Display WF.L, consumer_surplus, generation_costs, total_network_cost, RES_share, marginal_fixed_costs, fixed_linear_costs, max_marginal_price_adder, variable_costs;

Display CAP.L, GEN.L, GRID_CAP.L, load_deviation, res_share, o_price,  LOAD_real.L, FLOW.L, THETA.L;

o_WF =WF.L;
o_CS = consumer_surplus;
o_GC = generation_costs;
o_NC = total_network_cost;
o_cap(tec,n) = CAP.L(tec,n);
o_load(t,n) = LOAD_real.L(t,n);
o_RES_share = res_share;

execute_UNLOAD 'Output\Nodal\out.gdx' o_WF, o_CS, o_GC, o_NC, o_RES_share, o_cap, o_load, o_price, GEN;