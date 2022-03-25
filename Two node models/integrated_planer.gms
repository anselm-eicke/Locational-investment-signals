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
    
* Output Parameters
welfare
consumer_surplus
generation_costs
network_cost

load_deviation(t,n)
res_share

o_RES_share
o_load(t,n)
o_cap(tec,n)
o_gen(t,tec,n)
price(t,n)
;

* Load data
$GDXIN "Input/in.gdx"
$LOADdc i_cost

$GDXIN "Input/load.gdx"
$LOADdc i_load

$GDXIN "Input/avail.gdx"
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

a_nodal(t,n)                = p_ref *(1-1/elasticity);
s_nodal(t,n)                = p_ref *(1/(elasticity*load_ref(t,n)));

display sc, c_var, load_ref, avail, c_fix, a_nodal, s_nodal;

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

* sum((t,n), p_ref * LOAD_real(t,n) * (1-1/elasticity + LOAD_real(t,n) / (2*elasticity* load_ref(t,n))))
objective..                 WF =e= sum((t,n), a_nodal(t,n) * LOAD_real(t,n) + 1/2 * s_nodal(t,n) * LOAD_real(t,n) * LOAD_real(t,n))
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
LOAD_real.L(t,n) =load_ref(t,n);

Option optcr = 0.00000001;
Solve nodal using QCP max WF;

price(t,n) = p_ref * (1-1/elasticity + (LOAD_real.L(t,n)) / (elasticity * load_ref(t,n)));

network_cost = sum((n,m),(GRID_CAP.L(n,m) * grid_cost(n,m)) / 2);
consumer_surplus = sum((t,n), a_nodal(t,n) * LOAD_real.L(t,n) + 1/2 * s_nodal(t,n) * LOAD_real.L(t,n) * LOAD_real.L(t,n));

generation_costs = (sum((tec,n), CAP.L(tec,n) * c_fix(tec,n) + 0.5 * CAP.L(tec,n) * CAP.L(tec,n) * capacity_slope) + sum((t,tec,n), GEN.L(t,tec,n) * c_var(tec,n)));

load_deviation(t,n) =  LOAD_real.L(t,n) - load_ref(t,n);
res_share = 1 - sum((t,con,n), GEN.L(t,con,n)) / sum((t,tec,n), GEN.L(t,tec,n));
o_cap(tec,n) = CAP.L(tec,n);
o_gen(t,tec,n) = GEN.L(t,tec,n);
welfare = WF.L;

Display WF.L, consumer_surplus, generation_costs, network_cost, CAP.L, GEN.L,  price, load_deviation,  GRID_CAP.L, FLOW.L;

execute_UNLOAD 'Output/nodal.gdx' welfare, consumer_surplus, generation_costs, network_cost, o_gen, o_cap, price, c_fix;
