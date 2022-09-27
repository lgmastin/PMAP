clear,clc
close all

% Example usage of PMAP - Plume Model for Aggregate Prediction
%
% Written by Davis Hoffman, 2020 @ Stanford Univeristy

%% Model Setup
%%%%%%%%%%%%%%%%%%%%    Input: boundary conditions    %%%%%%%%%%%%%%%%%%%%%
input.Z_vent = 2200;       % vent elevation, m
input.r0 = 110;            % vent radius, m
input.u0 = 100;            % vent velocity, m/s
input.Tm = 910;            % magma temp, C
input.no = 0.03;           % mass fraction of gas in magma-gas mixture
input.m_ext_l = 0;         % mass fraction of external/surface water (liquid) added at vent
input.T_ext_l = 0;         % temperature of external/surface water (liquid) added at vent, C
input.m_ext_i = 0.165;     % mass fraction of external/surface water (ice) added at vent
input.T_ext_i = 0;         % temperature of external/surface water (ice) added at vent, C
input.eps = 0.09;          % entrainment ratio
input.emmisivity = 0.9;    % radiation emmisivity

%%%%%%%%%%%%%%%%%%%%%    Input: ambient conditions    %%%%%%%%%%%%%%%%%%%%%
input.metdata = true;      % read meteorological data from file?
input.metfile = 'input\Redoubt_2009.05.23.12Z_sounding.txt'; % provide source of met data

% if metfile is empty, enter conditions manually
input.Z_trop = 11e3;       % elevation base of tropopause, m
input.H_trop = 9e3;        % height of isothermal tropopause, m
input.Pamb_SL = 101.3e3;   % ambient pressure at sea level, Pa
input.Tamb_vent = -15;      % ambient temp, C
input.RHamb_vent = 0;      % relative humidity
input.dTdz_trop = -0.0065; % C/m
input.dTdz_strat = 0.0016; % C/m

%%%%%%%%%%%%%%%%%%%    Input: aggregation parameters    %%%%%%%%%%%%%%%%%%%
input.aggregation = true;       % run aggregation solver?
input.Df = 2.9;                 % fractal dimension
input.AggModel = [0; ...    % constant kernel
                  0; ...    % Saffman-Turner kernel
                  1; ...    % turbulent-inertial kernel
                  1; ...    % gravitational kernel
                  0];       % Brownian diffusion kernel
input.Ko = 1e-12;               % magnitude for constant kernel(if applicable), m^-3
input.StickModel = 'Hoffman';   % toggle sticking model: 'Hoffman' or 'Costa'
input.Chi = 0.23;               % detrainment parameters
input.GSD_filename = 'input\Redoubt_Event5_GSD.mat'; % filename that stores GSD data (by mass)

%%%%%%%%%%%%%%%%%%%%%%%%%%    Input: numerics    %%%%%%%%%%%%%%%%%%%%%%%%%%
input.numerics = '1st-order';   % '1st-order' or 'Cash-Karp'
input.tol = 2e-2;               % integration tolerance (for Cash-Karp integration only)
input.dz = 25;                  % elevation discritization, m (for 1st-order only)
input.u_stop = 1;               % stopping criteria, m/s

input.dphi = 1; % specify particle bin width for size property, phi units
input.phi_min = -5; % maximum particle size in phi units
input.phi_max = 11; % minimum particle size in phi units

input.drho = 750; % specify particle bin width for density property, kg/m^3
input.rho_min = 1000; % minimum particle density, kg/m^3

savefile = true; % Option: save output to file?
savefilename = 'output\plume_data';

plt = true; % Option: plot results?

%% Run PMAP
% output structure
%     Z - elevation, m
%     P - plume/ambient pressure, Pa
%     T - plume temperature, K
%     u - plume ascent speed, m/s
%     r - plume radius, m
%     rho - plume density, kg/m^3
%     Tamb - ambient temperature, K
%     rho_amb - ambient density, kg/m^3
%     m_a - air mass fraction in plume
%     m_m - magma mass fraction in plume
%     m_w - water mass fraction (all phases) in plume
%     m_l - liquid water mass fraction in plume
%     m_v - water vapor mass fraction in plume
%     m_i - ice mass fraction in plume
%     Nd - particle number density in each property classes, m^-3
%     n_detrain - detrainment rate of particles in each property class, s^-1
%     Mfall - mass of particles that leave each control volume,kg/s
%     Mp - particle mass defining fixed particle classes, kg
%     rho_p - particle density defining fixed particle classes, kg/m^3
output = PMAP(input);

%% Save results
if savefile
    save([savefilename '.mat'],'input','output')
end

%% Plot results
if plt
    post_processor
end
