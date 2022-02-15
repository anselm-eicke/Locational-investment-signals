Sets
all_t       all hours               /1*16/
t(all_t)    hours                   /1*12/
tec         generators              /base, peak, wind, solar/
con(tec)    conventional generation /base, peak/
all_n       all buses               /north, south/
n(all_n)    selected buses          /north, south/
;


alias (n,m);
alias (all_n,all_m);

* parameters for supply and demand functions
Parameter elasticity / -0.15 /; 
Parameter p_ref / 70 /;
Parameter specific_network_costs /200/;
Parameter capacity_slope / 0.5 /;
*Source for network costs: EMMA (3400 EUR/MW/km discontiert mit i = 0.07 ueber 40 Jahre)

Table B(all_n,all_m)        Susceptance of transmission lines
         north  south
north        1     700  
south      700       1
;

Parameters
o_instrument;

$GDXIN "Output/without_instrument.gdx"
$LOADdc o_instrument

Parameters
INSTRUMENT(tec,n)
sc          scaling factor
;

sc = card(t) / 8760;
*load instrument
INSTRUMENT(tec,n)           = o_instrument * sc * 1000;


Parameters
* Input Parameters
i_cost(*,*)                 cost data to be loaded from sheet "cost"
i_load(all_t,all_n)         load data to be loaded from sheet "time series" in MWh
i_avail(all_t,tec,all_n)    availability data

* Model Parameters
load_ref(t,n)                   hourly load in GWh
avail(t,tec,n)              availability of wind and solar generation (1)
c_var(tec,n)                variable costs (EUR per MWh)
c_fix(tec,n)                annualized fixed costs (EUR per MW p.a.)
cap_lim(tec,n)              capacity limit of generation in each node
grid_cost(n,m)
a_nodal(t,n)                intercept of inverse nodal demand function
s_nodal(t,n)                slope of inverse nodal demand function

A_zonal(t)                  intercept of inverse zonal demand function
S_zonal(t)                  slope of inverse zonal demand function
    
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
sum_instrument
network_cost


maximum
threshold

;

* Load data
$GDXIN "in.gdx"
$LOADdc i_cost

$GDXIN "load.gdx"
$LOADdc i_load

$GDXIN "avail.gdx"
$LOADdc i_avail
 


* Data assignment
load_ref(t,n)               = i_load(t,n) / 1000;
avail(t,tec,n)              = i_avail(t,tec,n);
avail(t,con,n)              = 1;
c_var(tec, n)               = i_cost(tec,"cost_var");
c_fix(tec, n)               = round(i_cost(tec,"cost_fix") * 1000 * sc);
cap_lim(tec,n)              = 100;
grid_cost(n,m)              = round(B(n,m) * specific_network_costs * sc);

*Inverse demand function at each node
a_nodal(t,n)                = p_ref *(1-1/elasticity);
s_nodal(t,n)                = p_ref *(1/(elasticity*load_ref(t,n)));

* Inverse demand function of the zonal market (only holds if P(t) < a_nodal(t,n) for all t,n)
A_zonal(t)                  = sum(n, a_nodal(t,n) / s_nodal(t,n)) / sum(n, 1/ s_nodal(t,n));
S_zonal(t)                  = 1 / sum(n, 1/ s_nodal(t,n));



display INSTRUMENT;

display c_var, load_ref, avail, c_fix, a_nodal, s_nodal, A_zonal, S_zonal;

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
o_price(t)
o_cap_instr(tec,n)

;


Variables
GEN(t,tec,n)
CAP(tec,n)
LOAD_zonal(t)
LAMBDA(t)
RHO(t,tec,n)
WF
;

Positive variables
mu_G_min(t,tec,n)
mu_G_max(t,tec,n)
mu_C_min(tec,n)
mu_C_max(tec,n)
mu_D_min(t)

UP(t,tec,n)
DOWN(t,tec,n)
SHARE(t,tec)
;

Equations
gen_min, gen_max,
cap_min, cap_max, demand_min,
energy_balance,
KKT_GEN, KKT_CAP, KKT_load;


** Inner problem
*Primal constraints
gen_min(t,tec,n)..          0 =g= -GEN(t,tec,n);
gen_max(t,tec,n)..          0 =g= GEN(t,tec,n) - CAP(tec,n) * avail(t,tec,n);

cap_min(tec,n)..            0 =g= -CAP(tec,n);
cap_max(tec,n)..            0 =g= CAP(tec,n) - cap_lim(tec,n);

demand_min(t)..             0 =g= -LOAD_zonal(t);

energy_balance(t)..         0 =e= sum((tec,n),GEN(t,tec,n)) - LOAD_zonal(t);
           
*KKT conditions
KKT_GEN(t,tec,n)..          c_var(tec,n) + mu_G_max(t,tec,n) - mu_G_min(t,tec,n) - LAMBDA(t) =e= 0;
KKT_CAP(tec,n)..            c_fix(tec,n) + capacity_slope * CAP(tec,n) + INSTRUMENT(tec,n) - sum(t,avail(t,tec,n) * mu_G_max(t,tec,n)) + mu_C_max(tec,n) - mu_C_min(tec,n) =e= 0;
KKT_load(t)..               -(A_zonal(t) + S_zonal(t) * LOAD_zonal(t)) - mu_D_min(t) + LAMBDA(t) =e= 0;              

Model LOCI /

gen_min.mu_G_min
gen_max.mu_G_max
cap_min.mu_C_min
cap_max.mu_C_max
demand_min.mu_D_min
energy_balance.LAMBDA
KKT_GEN
KKT_CAP
KKT_load

/;

* Set starting values
*LOAD_real.L(t,n) =load(t,n);

Option optcr = 0.0005;

Solve LOCI using MCP;

price(t) = A_zonal(t) + S_zonal(t) * LOAD_zonal.L(t);

maximum = smax(t,price(t));
threshold = smin((t,n),a_nodal(t,n));

network_cost_1 = sum((n,m),(GRID_CAP.L(n,m) / 2 * grid_cost(n,m)));
network_cost_2 = sum((t,tec,n), (UP.L(t,tec,n) - DOWN.L(t,tec,n)) * c_var(tec,n));
network_cost_3 = sum((t), A_zonal(t) * LOAD_spot.L(t) + 1/2 * S_zonal(t) * LOAD_spot.L(t) * LOAD_spot.L(t))
                - sum((t,n), a_nodal(t,n) * LOAD_redi.L(t,n) + 1/2 * s_nodal(t,n) * LOAD_redi.L(t,n) * LOAD_redi.L(t,n))
                ;
         
network_cost = network_cost_1 + network_cost_2 + network_cost_3;
       
consumer_surplus = sum((t), A_zonal(t) * LOAD_spot.L(t) + 1/2 * S_zonal(t) * LOAD_spot.L(t) * LOAD_spot.L(t));

generation_costs = (sum((tec,n), CAP.L(tec,n) * c_fix(tec,n) + 0.5 * CAP.L(tec,n) * CAP.L(tec,n) * capacity_slope) + sum((t,tec,n), GEN.L(t,tec,n) * c_var(tec,n)));
                                                                    
sum_instrument = sum((tec,n), INSTRUMENT.L * CAP.L(tec,n));

load_deviation(t,n) = ((SPOT_PRICE.L(t) - a_nodal(t,n)) / s_nodal(t,n)) - load_ref(t,n);
load_shedding(t,n) = LOAD_spot.L(t) - LOAD_redi.L(t,n);
res_share = 1 - sum((t,con,n), GEN.L(t,con,n)) / sum((t,tec,n), GEN.L(t,tec,n));
o_cap(tec,n) = CAP.L(tec,n);
o_gen(t,tec,n) = GEN.L(t,tec,n);
real_generation(t,tec,n) = GEN.L(t,tec,n) + UP.L(t,tec,n) - DOWN.L(t,tec,n);
welfare = WF.L;

display maximum, threshold;
    
Display WF.L, consumer_surplus, generation_costs, network_cost, network_cost_1, network_cost_2, network_cost_3, CAP.L, GEN.L, UP.L, DOWN.L, FLOW.L, price, load_deviation, load_shedding, GRID_CAP.L, LOAD_redi.L, LOAD_spot.L, o_instrument, sum_instrument;

execute_UNLOAD 'Output/without_instrument_MCP.gdx' welfare, consumer_surplus, generation_costs, network_cost, res_share, o_instrument, sum_instrument, o_cap, o_gen, price, c_fix;