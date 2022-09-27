function output = PMAP(input)
% 1-D wet vertical plume model based on Plumeria (Mastin, 2007).
% Integrates mass, momentum, and energy equations from BPT using forward
% Euler discretization. Model includes aggregation, particle detrainment,
% and mass coupling effects
%
% Written by Davis Hoffman, 2020 @ Stanford Univeristy

global eps g m_v_amb m_a_amb h_m Z_trop H_trop dTdz_trop ...
    dTdz_strat mw_H2O mw_air Df Nm Nrho Nc Mbin rho_bin Mp rho_p dp rp ...
    ETA T_wi_max T_wi_min hfi hfl hs Tsat R metdata  Zmet Tmet ...
    Tamb_vent Z_vent RHamb rho_m rho_l rho_i emmisivity sigma AggModel Ko ...
    Chi up_t aggregation kB StickModel

% constants/properties
R = 8314; % J/kmol.K
mw_H2O = 18; % molecular weight of water, g/mol
mw_air = 29; % molecular weight of air, g/mol
T_wi_max = 266.65; % max temperature at which liquid water & ice may coexist, K
T_wi_min = 258.15; % min temperature at which liquid water & ice may coexist, K
g = 9.81; % acceleration due to gravity, m/s^2
sigma = 5.67e-8; % Stefan-Boltzmann constant, W/m^2.K^4
kB = 1.381e-23; % Boltzmann constant, J/K
rho_l = 1000; % liquid water density, kg/m^3
rho_i = 920; % ice density, kg/m^3
rho_m = 2500; % magma density, kg/m^3
cp_m = 1000; % specific heat of magma, J/kg.K

% reallocate boundary conditions to new variables
Z_vent = input.Z_vent;
r0 = input.r0;
u0 = input.u0;
Tm = input.Tm+273.15; % convert to Kelvin
no = input.no;
m_ext_l = input.m_ext_l;
T_ext_l = input.T_ext_l+273.15; % convert to Kelvin
m_ext_i = input.m_ext_i;
T_ext_i = input.T_ext_i+273.15; % convert to Kelvin
eps = input.eps;
emmisivity = input.emmisivity;

metdata = input.metdata;
if input.metdata
    [Zmet,Pmet,Tmet,RHmet] = MetReader(input.metfile); % read meteorological data
    Tamb_vent = interp_custom(Zmet,Tmet,Z_vent); % interpolate at vent elevation
    Pamb_vent = interp_custom(Zmet,Pmet,Z_vent);
    RHamb_vent = interp_custom(Zmet,RHmet,Z_vent);
else
    Z_trop = input.Z_trop;
    H_trop = input.H_trop;
    Pamb_SL = input.Pamb_SL;
    Tamb_vent = input.Tamb_vent+273.15; % convert to Kelvin
    RHamb_vent = input.RHamb_vent;
    dTdz_trop = input.dTdz_trop;
    dTdz_strat = input.dTdz_strat;
    Pamb_vent = Pamb_SL*((Tamb_vent-Z_vent*dTdz_trop)/Tamb_vent)^ ...
            (g*mw_air/(R*dTdz_trop)); % ambient pressure at vent, Pa
end

aggregation = input.aggregation;
Df = input.Df;
AggModel = input.AggModel; 
Ko = input.Ko;
StickModel = input.StickModel;
Chi = input.Chi;

numerics = input.numerics;
tol = input.tol;
dz = input.dz;
u_stop = input.u_stop;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% setup volume bins using phi-scale
dphi = input.dphi;
phi_min = input.phi_min;
phi_max = input.phi_max;
phi = phi_max:-dphi:phi_min; % define bin array
dp = 1e-3*2.^-phi; % diameter bin array, m
Vp = pi/6*dp.^3; % equivalent volume of bins, m^3

% setup density bins using linear scale
drho = input.drho;
rho_min= input.rho_min;
rho_max = rho_m;
rho_bin = rho_max:-drho:rho_min; % define bin array
Nrho = numel(rho_bin); % number of rho bins

Mbin = Vp*rho_max; % mass property associated with each bin for parent ash
Nm = numel(Mbin);
Nc = Nm*Nrho+(Nm-1)*(Nrho-1); % total number of particle classes

% load GSD (by mass)
load(input.GSD_filename); Gm = GSD; clear GSD
f = Gm/sum(Gm); % mass fraction of parent ash

% make particle property grid
if Df<3
    [Mp,rho_p,ETA] = gen_secondary_grid;
else % if Df=3, there is no particle density variation (i.e., grid is one-dimensional)
    Mp = zeros(1,Nc); Mp(1:Nm) = Mbin;
    rho_p = zeros(1,Nc); rho_p(1:end) = rho_bin(1);
end

dp = (6/pi*Mp./rho_p).^(1/3); % equivalent diameter property associated with each bin, m
rp = dp/2; % equivalent radius property associated with each bin, m

disp(' ')

% define magma enthalpy function, J/kg
Tref = 273.15; % reference temp for magma, K
h_m = @(T) cp_m*(T-Tref);
clear Tref

% Initialize
Z = Z_vent;
r = r0;
u = u0;
P = Pamb_vent;
Tamb = Tamb_vent;
m_w = (1-m_ext_l-m_ext_i)*no+m_ext_l+m_ext_i; % initial water mass fraction (m_w = m_v+m_l+m_i)
m_a = 0; % initial air mass fraction
m_m = 1-m_w; % initial magma mass fraction
h = (1-m_ext_l-m_ext_i)*no*h_v(Tm)+m_ext_l*h_l(T_ext_l)+m_ext_i*h_i(T_ext_i)+ ...
    m_m*h_m(Tm); % calculate initial enthalpy assuming ideal mixture
Tsat = get_saturation_temp(m_w/mw_H2O/(m_w/mw_H2O+m_a/mw_air)*P);
hs = m_w*h_v(Tsat)+m_a*h_a(Tsat)+m_m*h_m(Tsat);
hfi = get_fi_enthalpy(m_w,m_a,m_m,P);
hfl = get_fl_enthalpy(m_w,m_a,m_m,P);
[T,m_v,m_l,m_i] = solve_temp(m_m,m_a,m_w,h,P); % solve initial mixture temperature
rho_v = P*mw_H2O/R/T; % density of vent water vapor
rho = (m_v/rho_v+m_l/rho_l+m_m/rho_m)^-1; % calculate initial density assuming ideal mixture, kg/m^3
Mtot = pi*r0^2*rho*u0; % initial mass flux, kg/s
Mm = m_m*Mtot; % mass flux of magma at vent, kg/s
Mom = Mtot*u0; % initial momentum flux, kg*m/s^2
E = Mtot*(h+1/2*u0^2+g*Z); % initial energy flux, J/s
n = zeros(1,Nc); % particle flux rate, s^-1
n(1:Nm) = f*m_m*Mtot./Mbin;
n_detrain = zeros(1,Nc); % particle detrainment rate, s^-1
Nd = n/(pi*r0^2*u0); % particle number density, m^-3
Mc = Mp.*n; % particle mass flux rate, kg/s
RH = solve_RH(T,P,m_w,m_a); % initial humidity inside plume
[m_v_amb,m_a_amb] = get_ambient_mass_fractions(Tamb_vent,Pamb_vent,RHamb_vent); % initial vapor/air mass fractions in ambient
h_amb = m_v_amb*h_v(Tamb_vent)+m_a_amb*h_a(Tamb_vent); % initial ambient enthalpy
rho_amb = Pamb_vent/((m_v_amb/mw_H2O+m_a_amb/mw_air)*R*Tamb_vent); % initial ambient density
RHamb = RHamb_vent;
up_t = ones(1,Nc);

% Iterate
iter = 1;
while true
    fprintf('Integrating plume equations through CV #%d (Elevation: %.2f km)\n',iter,Z(iter)/1e3)
   
    % get properties of plume gas phase and calculate particle settling velocities
    mu = 1.458e-6*T(iter).^1.5./(T(iter)+110.4); % dynamic viscosity, kg/m.s
    m_a_g = m_a(iter)/(m_a(iter)+m_v(iter));
    m_v_g = m_v(iter)/(m_a(iter)+m_v(iter));
    rho_g = P(iter)/(m_a_g/mw_air+m_v_g/mw_H2O)/R/T(iter); % plume gas density, kg/m^3
    nu = mu./rho_g; % kinematic viscosity, m^2/s
    up_t = get_terminal_velocity(dp,rho_p,rho_g,nu,g,up_t); % particle settling speeds, m/s
    
    % advance variables to next control volume
    x = Z(iter);
    yarr(1) = Mtot(iter);
    yarr(2) = Mom(iter);
    yarr(3) = E(iter);
    yarr(4) = P(iter);
    yarr(5) = m_m(iter);
    yarr(6) = m_a(iter);
    yarr(7) = m_w(iter);
    yarr(8:8+Nc-1) = n;
    switch numerics
        case '1st-order'
            [x_next,yarr,dz] = advance_forward_Euler(x,yarr,dz);
        case 'Cash-Karp'
            [x_next,yarr,dz] = advance_RKCK_adapt(x,yarr,dz,tol);
    end
    Z(iter+1,1) = x_next;
    Mtot(iter+1,1) = yarr(1);
    Mom(iter+1,1) = yarr(2);
    E(iter+1,1) = yarr(3);
    P(iter+1,1) = yarr(4);
    m_m(iter+1,1) = yarr(5);
    m_a(iter+1,1) = yarr(6);
    m_w(iter+1,1) = yarr(7);
    n = yarr(8:end);
    
    % check the stopping conditions
    u(iter+1,1) = Mom(iter+1)/Mtot(iter+1);
    if u(iter+1)<u_stop % break early if we overshoot the stopping condition
        Z(iter+1) = []; Mtot(iter+1) = []; Mom(iter+1) = []; P(iter+1) = [];
        E(iter+1) = []; m_m(iter+1) = []; m_a(iter+1) = []; m_w(iter+1) = [];
        u(iter+1) = [];
        break
    end
    
    % keep track of particle detrainment and fallout mass
    n_detrain(iter,1:Nc) = Chi*up_t/r(iter)/u(iter).*n*dz;
    Mfall(iter,1) = sum(Mp.*n_detrain(iter,:));

    if metdata
        P(iter+1,1) = interp_custom(Zmet,Pmet,Z(iter+1)); % overwrite pressure variable
        Tamb(iter+1,1) = interp_custom(Zmet,Tmet,Z(iter+1));
        RHamb = interp_custom(Zmet,Tmet,Z(iter+1));
    else
        Tamb(iter+1,1) = solve_temp_amb(Z(iter+1));
    end

    h(iter+1,1) = E(iter+1)/Mtot(iter+1)-1/2*u(iter+1)^2-g*Z(iter+1);
    Tsat = get_saturation_temp(m_w(iter+1)/mw_H2O/(m_w(iter+1)/mw_H2O+m_a(iter+1)/mw_air)*P(iter+1));
    hs = m_w(iter+1)*h_v(Tsat)+m_a(iter+1)*h_a(Tsat)+m_m(iter+1)*h_m(Tsat);
    hfi = get_fi_enthalpy(m_w(iter+1),m_a(iter+1),m_m(iter+1),P(iter+1));
    hfl = get_fl_enthalpy(m_w(iter+1),m_a(iter+1),m_m(iter+1),P(iter+1));
    [T(iter+1,1),m_v(iter+1,1),m_l(iter+1,1),m_i(iter+1,1)] = ...
        solve_temp(m_m(iter+1,1),m_a(iter+1,1),m_w(iter+1,1),h(iter+1,1),P(iter+1,1));
    rho(iter+1,1) = ((m_a(iter+1)/mw_air+m_v(iter+1)/mw_H2O)* ...
        R*T(iter+1)/P(iter+1)+m_l(iter+1)/rho_l+m_i(iter+1)/rho_i+m_m(iter+1)/rho_m)^-1;
    rho_amb(iter+1,1) = P(iter+1)/((m_v_amb/mw_H2O+m_a_amb/mw_air)*R*Tamb(iter+1));
    h_amb(iter+1,1) = m_v_amb*h_v(Tamb(iter+1))+m_a_amb*h_a(Tamb(iter+1));
    r(iter+1,1) = sqrt(Mtot(iter+1)/(pi*rho(iter+1)*u(iter+1)));
    RH(iter+1,1) = solve_RH(T(iter+1),P(iter+1),m_w(iter+1),m_a(iter+1));
    [m_v_amb,m_a_amb] = get_ambient_mass_fractions(Tamb(iter+1),P(iter+1),RHamb);
    Nd(iter+1,:) = n/(pi*r(iter+1)^2*u(iter+1));
    Mc = Mp.*n;
    
    iter = iter+1;
end
% assume all particles remaining in the final control volume are detrained
n_detrain(end,:) = n_detrain(end,:)+n(end,:);

% summarize results into a structure named 'output'
output.Z = Z;
output.P = P;
output.T = T;
output.u = u;
output.r = r;
output.rho = rho;
output.Tamb = Tamb;
output.rho_amb = rho_amb;
output.m_a  = m_a;
output.m_m = m_m;
output.m_w = m_w;
output.m_l = m_l;
output.m_v = m_v;
output.m_i = m_i;
output.Nd = Nd;
output.n_detrain = n_detrain;
output.Mfall = Mfall;
output.Mp = Mp;
output.rho_p = rho_p;
end