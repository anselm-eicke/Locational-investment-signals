Sets
all_t       all hours               /1*48/
t(all_t)    hours                   /1*24/
tec         generators              /base, peak, wind, solar/
con(tec)    conventional generation /base, peak/
all_n       all buses               /north, south/
n(all_n)    selected buses          /north, south/
;


alias (n,m);
alias (all_n,all_m);

* parameters for supply and demand functions
Parameter elasticity / -0.25 /; 
Parameter p_ref / 65 /;
Parameter specific_network_costs /300/;
*Source for network costs: EMMA (3400 EUR/MW/km discontiert mit i = 0.07 ueber 40 Jahre)

Table B(all_n,all_m)        Susceptance of transmission lines
         north  south
north        1     700  
south      700       1
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
generation_slope
cap_lim(tec,n)              capacity limit of generation in each node
grid_cost(n,m)
sc                          scaling factor
    
* Output Parameters
consumer_surplus
generation_costs
network_cost
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
o_gen(t,tec,n)
price(t,n)
price_lambda(t,n)
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
capacity_slope              = 0.15;
cap_lim(tec,n)              = 100;
grid_cost(n,m)              = round(B(n,m) * specific_network_costs * sc);

display c_var, load, avail, c_fix;



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


objective..                  WF =e= sum((t,n), p_ref * LOAD_real(t,n) * (1-1/elasticity + LOAD_real(t,n) / (2*elasticity* load(t,n))))
                                    - sum((tec,n), CAP(tec,n) * c_fix(tec,n) + 0.5 * CAP(tec,n) * CAP(tec,n) * capacity_slope)
                                    - sum((t,tec,n), GEN(t,tec,n) * c_var(tec,n))
                                    - sum((n,m),(GRID_CAP(n,m) * grid_cost(n,m)) / 2);
                                    
cap_constraint(t,tec,n)..   GEN(t,tec,n)    =L= CAP(tec,n) * avail(t,tec,n);
cap_limit(tec,n)..          CAP(tec,n)      =L= cap_lim(tec,n);

*nodal energy balance
nodal_energy_balance(t,n).. sum(tec,GEN(t,tec,n)) - LOAD_real(t,n) =E= sum(m,FLOW(t,n,m));

*network constraints
grid_eq1(t,n,m)..           FLOW(t,n,m)     =l= GRID_CAP(n,m);
grid_eq2(n,m)..             GRID_CAP(n,m)   =e= GRID_CAP(m,n);
grid_eq3(t,n,m)..           FLOW(t,n,m)     =e= B(n,m) * (THETA(t,n) - THETA(t,m));
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

Solve nodal using QCP max WF;

price(t,n) = p_ref * (1-1/elasticity + (LOAD_real.L(t,n)) / (elasticity * load(t,n)));

network_cost = sum((n,m),(GRID_CAP.L(n,m) * grid_cost(n,m)) / 2);
consumer_surplus = sum(t, p_ref * sum((n), LOAD_real.L(t,n)) * (1-1/elasticity + sum((n), LOAD_real.L(t,n)) / (2*elasticity* sum(n,load(t,n)))));

generation_costs = (sum((tec,n), CAP.L(tec,n) * c_fix(tec,n) + 0.5 * CAP.L(tec,n) * CAP.L(tec,n) * capacity_slope) + sum((t,tec,n), GEN.L(t,tec,n) * c_var(tec,n)));

load_deviation(t,n) =  LOAD_real.L(t,n) / load(t,n);
res_share = 1 - sum((t,con,n), GEN.L(t,con,n)) / sum((t,tec,n), GEN.L(t,tec,n));
o_cap(tec,n) = CAP.L(tec,n);
o_gen(t,tec,n) = GEN.L(t,tec,n);


Display GEN.L, CAP.L, price, load_deviation, network_cost, GRID_CAP.L;

execute_UNLOAD 'Output/nodal.gdx' consumer_surplus, generation_costs, network_cost, res_share, o_cap, o_gen, price, c_fix;
