Sets
all_t       all hours               /1*48/
t(all_t)    hours                   /1*5/
tec         generators              /base, peak, wind, solar/
con(tec)    conventional generation /base, peak/
all_n       all buses               /north, south/
n(all_n)    selected buses          /north, south/
;

alias (n,nn);
alias (n,m);
alias (all_n,all_m);

* parameters for supply and demand functions
Parameter elasticity / -0.25 /; 
Parameter p_ref / 65 /;
Parameter specific_network_costs /0/;
Parameter capacity_slope / 0.5 /;
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
load_ref(t,n)                   hourly load in GWh
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
INSTRUMENT(tec,n)

maximum
threshold

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

INSTRUMENT(tec,n)           = 0;

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

;

Parameter M1 / 10000/;
Parameter M2 / 10000/;
Parameter M3 / 10000/;
Parameter M4 / 10000/;
Parameter M5 / 10000/;
Parameter M6 / 10000/;

Binary variables y1(t,tec,n),y2(t,tec,n),y3(tec,n),y4(tec,n),y5(t),y6(t);



Variables
GEN(t,tec,n)
CAP(tec,n)
LOAD_spot(t)
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
objective

gen_min, gen_max,
cap_min, cap_max, demand_min,
energy_balance,
KKT_GEN, KKT_CAP, KKT_load

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

objective..                     WF =e= 1;

** Inner problem
*Primal constraints
gen_min(t,tec,n)..          0 =g= -GEN(t,tec,n);
gen_max(t,tec,n)..          0 =g= GEN(t,tec,n) - CAP(tec,n) * avail(t,tec,n);

cap_min(tec,n)..            0 =g= -CAP(tec,n);
cap_max(tec,n)..            0 =g= CAP(tec,n) - cap_lim(tec,n);

demand_min(t)..             0 =g= -LOAD_spot(t);

energy_balance(t)..         0 =e= sum((tec,n),GEN(t,tec,n)) - LOAD_spot(t);
           
*KKT conditions
KKT_GEN(t,tec,n)..          c_var(tec,n) + mu_G_max(t,tec,n) - mu_G_min(t,tec,n) - LAMBDA(t) =e= 0;
KKT_CAP(tec,n)..            c_fix(tec,n) + capacity_slope * CAP(tec,n) + INSTRUMENT(tec,n) - sum(t,avail(t,tec,n) * mu_G_max(t,tec,n)) + mu_C_max(tec,n) - mu_C_min(tec,n) =e= 0;
KKT_load(t)..               -(A_zonal(t) + S_zonal(t) * LOAD_spot(t)) - mu_D_min(t) + LAMBDA(t) =e= 0;              

complementarity1a(t,tec,n)..    GEN(t,tec,n)        =L= y1(t,tec,n) * M1;
complementarity1b(t,tec,n)..    mu_G_min(t,tec,n)   =L= (1-y1(t,tec,n)) * M1;
complementarity2a(t,tec,n)..    CAP(tec,n) * avail(t,tec,n) - GEN(t,tec,n) =L= y2(t,tec,n) * M2;
complementarity2b(t,tec,n)..    mu_G_max(t,tec,n)   =L= (1-y2(t,tec,n)) * M2;
complementarity3a(tec,n)..      CAP(tec,n)          =L= y3(tec,n) * M3;
complementarity3b(tec,n)..      mu_C_min(tec,n)     =L= (1-y3(tec,n)) * M3;
complementarity4a(tec,n)..      cap_lim(tec,n) - CAP(tec,n) =L= y4(tec,n) * M4;
complementarity4b(tec,n)..      mu_C_max(tec,n)     =L= (1-y4(tec,n)) * M4;
complementarity5a(t)..          LOAD_spot(t)        =L= y5(t) * M5;
complementarity5b(t)..          mu_D_min(t)         =L= (1-y5(t)) * M5;
complementarity6a(t)..          LOAD_spot(t) - sum((tec,n),GEN(t,tec,n)) =L= y6(t) * M6;
complementarity6b(t)..          LAMBDA(t)           =L= (1-y6(t)) * M6;


Model LOCI /

objective

gen_min
gen_max
cap_min
cap_max

demand_min
energy_balance

KKT_GEN
KKT_CAP
KKT_load

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

/;

* Set starting values
LOAD_spot.L(t) =sum(n,load_ref(t,n));
CAP.L(tec,n) = 50;
GEN.up(t,tec,n) = 100;
GEN.lo(t,tec,n) = 0;

DOWN.up(t,tec,n) = 100;
DOWN.lo(t,tec,n) = 0;

UP.up(t,tec,n) = 100;
UP.lo(t,tec,n) = 0;

*Option optcr = 0.5;

LOCI.nodlim = 25000000;
LOCI.resLim = 40000;

Solve LOCI maximizing WF using MIQCP;

price(t) = A_zonal(t) + S_zonal(t) * LOAD_spot.L(t);

maximum = smax(t,price(t));
threshold = smin((t,n),a_nodal(t,n));

display maximum, threshold;
    
    
load_deviation(t) = LOAD_spot.L(t) / sum(n,load_ref(t,n));

Display GEN.L, CAP.L, LOAD_spot.L, LAMBDA.L, price, load_deviation;
*,WF.L,  GRID.L, LAMBDA.L, LOAD_real.L, CAP_INSTR.L, UP.L, DOWN.L, RESERVE_CAP.L, FLOW.L;