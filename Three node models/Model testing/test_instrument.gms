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


* parameters for supply and demand functions
Parameter elasticity / -0.25 /; 
Parameter p_ref / 65 /;
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
generation_slope
cap_lim(tec,n)              capacity limit of generation in each node
grid_cost(n,m)
sc                          scaling factor
    
* Output Parameters
consumer_surplus
generation_costs
total_network_cost
load_deviation(t)
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
price(t)
o_cap_instr(tec,n)
o_instrument(tec,n)
sum_instrument
network_cost

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

Parameters
consumer_surplus
generation_costs
total_network_cost
res_share

o_WF
o_CS
o_GC
o_NC
o_RES_share

o_cap(tec,n)
price(t)
load_deviation(t)
;

Parameters
INSTRUMENT(tec,n)    = 36.12

;

Variables
GEN(t,tec,n)
CAP(tec,n)
LOAD_real(t,n)
LAMBDA(t)
RHO(t,tec,n)
WF
;

Positive variables
mu_G_min(t,tec,n)
mu_G_max(t,tec,n)
mu_C_min(tec,n)
mu_C_max(tec,n)
mu_D_min(t,n)

UP(t,tec,n)
DOWN(t,tec,n)
SHARE(t,tec)
;


Equations
Objective, gen_min, gen_max, cap_min, cap_max;

Objective..             WF =e= sum(t, p_ref * sum((tec,n), GEN(t,tec,n)) * (1-1/elasticity +  sum((tec,n), GEN(t,tec,n)) / (2*elasticity* sum(n,load(t,n)))))
                                - sum((tec,n), CAP(tec,n)* (c_fix(tec,n) + INSTRUMENT(tec,n)) + 0.5 * CAP(tec,n) * CAP(tec,n) * capacity_slope)
                                - sum((t,tec,n), GEN(t,tec,n) * c_var(tec,n));
                
gen_max(t,tec,n)..      0 =g= GEN(t,tec,n) - CAP(tec,n) * avail(t,tec,n);
gen_min(t,tec,n)..      0 =g= -GEN(t,tec,n);
cap_max(tec,n)..        0 =g= CAP(tec,n) - cap_lim(tec,n);
cap_min(tec,n)..        0 =g= -CAP(tec,n);


Model LOCI /

objective
gen_min
gen_max
cap_min
cap_max

/;


Solve LOCI maximizing WF using QCP;

price(t) = p_ref * (1-(1/elasticity) + (sum((tec,n), GEN.L(t,tec,n)) / sum(n, elasticity * load(t,n))));
load_deviation(t) = sum((tec,n), GEN.L(t,tec,n)) / sum(n,load(t,n));

Display GEN.L, CAP.L, price, load_deviation;
*,WF.L,  GRID.L, LAMBDA.L, LOAD_real.L, CAP_INSTR.L, UP.L, DOWN.L, RESERVE_CAP.L, FLOW.L;
