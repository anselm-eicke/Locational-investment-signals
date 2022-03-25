Sets
all_t       all hours               /1*10/
t(all_t)    hours                   /1*10/
tec         generators              /base, peak, wind, solar/
con(tec)    conventional generation /base, peak/
all_n       all buses               /north, south/
n(all_n)    selected buses          /north, south/
;


alias (n,m);
alias (all_n,all_m);

* parameters for supply and demand functions
Parameter elasticity / -0.05 /; 
Parameter p_ref / 55 /;
Parameter specific_network_costs /200/;
Parameter capacity_slope / 333 /;
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
load_ref(t,n)               hourly load in GWh
avail(t,tec,n)              availability of wind and solar generation (1)
c_var(tec,n)                variable costs (EUR per MWh)
c_fix(tec,n)                annualized fixed costs (EUR per MW p.a.)
cap_lim(tec,n)              capacity limit of generation in each node
grid_cost(n,m)
sc                          scaling factor
a_nodal(t,n)                intercept of inverse nodal demand function
s_nodal(t,n)                slope of inverse nodal demand function

A_zonal(t)                  intercept of inverse zonal demand function
S_zonal(t)                  slope of inverse zonal demand function
    
* Output Parameters
welfare,
consumer_surplus
generation_costs
load_deviation(t,n)
load_shedding(t,n)
network_cost
network_cost_1
network_cost_2
network_cost_3

res_share
real_generation(t,tec,n)

o_RES_share
o_load(t,n)
o_cap(tec,n)
o_gen(t,tec,n)
price(t)
o_instrument(tec,n)
sum_instrument
redispatch(t,tec,n) 
;

* Load data
$GDXIN "in.gdx"
$LOADdc i_cost

$GDXIN "load.gdx"
$LOADdc i_load

$GDXIN "avail.gdx"
$LOADdc i_avail
 

* Data assignment
sc = card(t) / 8760;
load_ref(t,n)               = i_load(t,n) / 1000;
avail(t,tec,n)              = i_avail(t,tec,n);
avail(t,con,n)              = 1;
c_var(tec, n)               = i_cost(tec,"cost_var");
c_fix(tec, n)               = round(i_cost(tec,"cost_fix") * 1000 * sc);
cap_lim(tec,n)              = 100;
grid_cost(n,m)              = round(B(n,m) * specific_network_costs * sc);
capacity_slope              = capacity_slope * sc;

*Inverse demand function at each node
a_nodal(t,n)                = p_ref *(1-1/elasticity);
s_nodal(t,n)                = p_ref *(1/(elasticity*load_ref(t,n)));

* Inverse demand function of the zonal market (only holds if P(t) < a_nodal(t,n) for all t,n)
A_zonal(t)                  = sum(n, a_nodal(t,n) / s_nodal(t,n)) / sum(n, 1/ s_nodal(t,n));
S_zonal(t)                  = 1 / sum(n, 1/ s_nodal(t,n));

display c_var, load_ref, avail, c_fix, a_nodal, s_nodal, A_zonal, S_zonal;

Binary variables y1(t,tec,n),y2(t,tec,n),y3(tec,n),y4(tec,n),y5(t),y6(t);

Parameter M1 / 100000/;
Parameter M2 / 100000/;
Parameter M3 / 100000/;
Parameter M4 / 100000/;
Parameter M5 / 100000/;
Parameter M6 / 100000/;

Free Variables
GEN(t,tec,n)
CAP(tec,n)
WF
FLOW(t,n,m)
INSTRUMENT(tec,n)
THETA(t,n)
SPOT_PRICE(t)
;

Positive variables
mu_G_min(t,tec,n)
mu_G_max(t,tec,n)
mu_C_min(tec,n)
mu_C_max(tec,n)
mu_D_min(t)

GRID_CAP(n,m)
LOAD_redi(t,n)
LOAD_spot(t)
UP(t,tec,n)
DOWN(t,tec,n)
;

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

objective..                     WF =e= sum((t,n), a_nodal(t,n) * LOAD_redi(t,n) + 1/2 * s_nodal(t,n) * LOAD_redi(t,n) * LOAD_redi(t,n))
                                    - sum((tec,n), CAP(tec,n) * c_fix(tec,n) + 0.5 * CAP(tec,n) * CAP(tec,n) * capacity_slope)
                                    - sum((t,tec,n), GEN(t,tec,n) * c_var(tec,n))
                                    - sum((t,tec,n), UP(t,tec,n) * (c_var(tec,n) + 25) - DOWN(t,tec,n) * (c_var(tec,n) - 25))
                                    - sum((n,m),(GRID_CAP(n,m) * grid_cost(n,m)) / 2)
                                    ;          

nodal_energy_balance(t,n)..     sum(tec,GEN(t,tec,n) - DOWN(t,tec,n) + UP(t,tec,n)) - LOAD_redi(t,n) =E= sum(m,FLOW(t,n,m));

*network constraints
grid_eq1(t,n,m)..               FLOW(t,n,m) =l= GRID_CAP(n,m);
grid_eq2(n,m)..                 GRID_CAP(n,m) =e= GRID_CAP(m,n);
grid_eq3(t,n,m)..               FLOW(t,n,m) =e= B(n,m) *(THETA(t,n) - THETA(t,m));
grid_eq4(t,n)..                 THETA(t,'south') =e= 0;

redispatch1(t,tec,n)..          DOWN(t,tec,n) =L= GEN(t,tec,n);
redispatch2(t,tec,n)..          UP(t,tec,n) =L= CAP(tec,n) * avail(t,tec,n) - GEN(t,tec,n);

** INNER PROBLEM
gen_min(t,tec,n)..              0 =g= -GEN(t,tec,n);
gen_max(t,tec,n)..              0 =g= GEN(t,tec,n) - CAP(tec,n) * avail(t,tec,n);

cap_min(tec,n)..                0 =g= -CAP(tec,n);
cap_max(tec,n)..                0 =g= CAP(tec,n) - cap_lim(tec,n);

demand_min(t)..                 0 =g= -LOAD_spot(t);
energy_balance(t)..             0 =e= sum((tec,n),GEN(t,tec,n)) - LOAD_spot(t);
           
KKT_GEN(t,tec,n)..              c_var(tec,n) + mu_G_max(t,tec,n) - mu_G_min(t,tec,n) - SPOT_PRICE(t) =e= 0;
KKT_CAP(tec,n)..                c_fix(tec,n) + capacity_slope * CAP(tec,n) + INSTRUMENT(tec,n) - sum(t,avail(t,tec,n) * mu_G_max(t,tec,n)) + mu_C_max(tec,n) - mu_C_min(tec,n) =e= 0;
KKT_load(t)..                   -(A_zonal(t) + S_zonal(t) * LOAD_spot(t)) - mu_D_min(t) + SPOT_PRICE(t) =e= 0;              

complementarity1a(t,tec,n)..    GEN(t,tec,n)        =L= y1(t,tec,n) * M1;
complementarity1b(t,tec,n)..    mu_G_min(t,tec,n)   =L= (1-y1(t,tec,n)) * M1;
complementarity2a(t,tec,n)..    CAP(tec,n) * avail(t,tec,n) - GEN(t,tec,n) =L= y2(t,tec,n) * M2;
complementarity2b(t,tec,n)..    mu_G_max(t,tec,n)   =L= (1-y2(t,tec,n)) * M2;
complementarity3a(tec,n)..      CAP(tec,n)          =L= y3(tec,n) * M3;
complementarity3b(tec,n)..      mu_C_min(tec,n)     =L= (1-y3(tec,n)) * M3;
complementarity4a(tec,n)..      cap_lim(tec,n) - CAP(tec,n) =L= y4(tec,n) * M4;
complementarity4b(tec,n)..      mu_C_max(tec,n)     =L= (1-y4(tec,n)) * M4;
complementarity5a(t)..          LOAD_spot(t)        =L= y5(t) * M5;
complementarity5b(t)..          mu_D_min(t)       =L= (1-y5(t)) * M5;
complementarity6a(t)..          LOAD_spot(t) - sum((tec,n),GEN(t,tec,n)) =L= y6(t) * M6;
complementarity6b(t)..          SPOT_PRICE(t)           =L= (1-y6(t)) * M6;


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

demand_min
energy_balance

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


INSTRUMENT.lo(tec,n) = -200;
INSTRUMENT.up(tec,n) = 200;

GEN.up(t,tec,n) = 100;
GEN.lo(t,tec,n) = 0;

DOWN.up(t,tec,n) = 100;
DOWN.lo(t,tec,n) = 0;

UP.up(t,tec,n) = 100;
UP.lo(t,tec,n) = 0;

LOCI.nodlim = 65000000;
LOCI.resLim = 150000;

Option optcr = 0.0000001;

Option MIQCP = Cplex;


Solve LOCI maximizing WF using MIQCP;

price(t) = SPOT_PRICE.L(t);

o_instrument(tec,n) = INSTRUMENT.L(tec,n) / sc / 1000;
                                    
network_cost_1 = sum((n,m),(GRID_CAP.L(n,m) / 2 * grid_cost(n,m)));
network_cost_2 = sum((t,tec,n), (UP.L(t,tec,n) - DOWN.L(t,tec,n)) * c_var(tec,n));
network_cost_3 = sum((t), A_zonal(t) * LOAD_spot.L(t) + 1/2 * S_zonal(t) * LOAD_spot.L(t) * LOAD_spot.L(t))
                - sum((t,n), a_nodal(t,n) * LOAD_redi.L(t,n) + 1/2 * s_nodal(t,n) * LOAD_redi.L(t,n) * LOAD_redi.L(t,n))
                ;
         
network_cost = network_cost_1 + network_cost_2 + network_cost_3;
       
consumer_surplus = sum((t), A_zonal(t) * LOAD_spot.L(t) + 1/2 * S_zonal(t) * LOAD_spot.L(t) * LOAD_spot.L(t));

generation_costs = (sum((tec,n), CAP.L(tec,n) * c_fix(tec,n) + 0.5 * CAP.L(tec,n) * CAP.L(tec,n) * capacity_slope) + sum((t,tec,n), GEN.L(t,tec,n) * c_var(tec,n)));
                                                                    
sum_instrument = sum((tec,n), INSTRUMENT.L(tec,n) * CAP.L(tec,n));

load_deviation(t,n) = LOAD_spot.L(t) - load_ref(t,n);
load_shedding(t,n) = LOAD_spot.L(t) - LOAD_redi.L(t,n);
res_share = 1 - sum((t,con,n), GEN.L(t,con,n)) / sum((t,tec,n), GEN.L(t,tec,n));
o_cap(tec,n) = CAP.L(tec,n);
o_gen(t,tec,n) = GEN.L(t,tec,n);
real_generation(t,tec,n) = GEN.L(t,tec,n) + UP.L(t,tec,n) - DOWN.L(t,tec,n);
welfare = WF.L;

redispatch(t,tec,n) = UP.L(t,tec,n) - DOWN.L(t,tec,n);

Display WF.L, consumer_surplus, generation_costs, network_cost, network_cost_1, network_cost_2, network_cost_3, CAP.L, GEN.L, UP.L, DOWN.L, redispatch, FLOW.L, price, load_deviation, load_shedding, GRID_CAP.L, LOAD_redi.L, LOAD_spot.L, INSTRUMENT.L, o_instrument, sum_instrument;

execute_UNLOAD 'Output/with_instrument_redispatch.gdx' welfare, consumer_surplus, generation_costs, network_cost, network_cost_1, network_cost_2, network_cost_3, res_share, o_instrument, sum_instrument, o_cap, o_gen, price, c_fix;