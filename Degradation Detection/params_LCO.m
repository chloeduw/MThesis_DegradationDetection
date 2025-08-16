%% Params for Electrochemical Model
%   Created May 24, 2012 by Scott Moura
%
%   Most parameters are from DUALFOIL model
%   D_e is from Capiglia et al, for c_e = 1000 mol/m^3
%   Equilibrium potentials are from Bosch Klein TCST 2011
%   Electrode area chosen to correspond to 2.3 Ah cell

% Modified Apr 26, 2016 by Saehong Park
% Added p.epsilon_f_n, p.epsilon_f_p to calculate Bruggeman relationship.

%% Geometric Params
% Thickness of each layer
p.L_n = 100e-6;     % Thickness of negative electrode [m]
p.L_s = 25e-6;     % Thickness of separator [m]
p.L_p = 100e-6;     % Thickness of positive electrode [m]

L_ccn = 25e-6;    % Thickness of negative current collector [m]
L_ccp = 25e-6;    % Thickness of negative current collector [m]

% Particle Radii
p.R_s_n = 10e-6;   % Radius of solid particles in negative electrode [m]
p.R_s_p = 10e-6;   % Radius of solid particles in positive electrode [m]

% Volume fractions
p.epsilon_s_n = 0.6;      % Volume fraction in solid for neg. electrode
p.epsilon_s_p = 0.5;      % Volume fraction in solid for pos. electrode

p.epsilon_e_n = 0.3;   % Volume fraction in electrolyte for neg. electrode
p.epsilon_e_s = 1.0;   % Volume fraction in electrolyte for separator
p.epsilon_e_p = 0.3;   % Volume fraction in electrolyte for pos. electrode

% make element to caclulate phi_{s} by Saehong Park 
p.epsilon_f_n = 0.1;  % Volume fraction of filler in neg. electrode
p.epsilon_f_p = 0.2;  % Volume fraction of filler in pos. electrode

epsilon_f_n = p.epsilon_f_n;  % Volume fraction of filler in neg. electrode
epsilon_f_p = p.epsilon_f_p;  % Volume fraction of filler in pos. electrode


% Specific interfacial surface area
p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;  % Negative electrode [m^2/m^3]
p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;  % Positive electrode [m^2/m^3]

% Mass densities
p.rho_sn = 1800;    % Solid phase in negative electrode [kg/m^3]
p.rho_sp = 5010;    % Solid phase in positive electrode [kg/m^3]
p.rho_e =  1324;    % Electrolyte [kg/m^3]
p.rho_f = 1800;     % Filler [kg/m^3]
p.rho_ccn = 8954;   % Current collector in negative electrode
p.rho_ccp = 2707;   % Current collector in positive electrode

% Compute cell mass [kg/m^2]
p.m_n = p.L_n * (p.rho_e*p.epsilon_e_n + p.rho_sn*p.epsilon_s_n + p.rho_f*epsilon_f_n);
p.m_s = p.L_s * (p.rho_e*p.epsilon_e_n);
p.m_p = p.L_p * (p.rho_e*p.epsilon_e_p + p.rho_sp*p.epsilon_s_p + p.rho_f*epsilon_f_p);
p.m_cc = p.rho_ccn*L_ccn + p.rho_ccp*L_ccp;

% Lumped density [kg/m^2]
p.rho_avg = p.m_n + p.m_s + p.m_p + p.m_cc;

%% Transport Params
% Diffusion coefficient in solid
p.D_s_n0 = 3.9e-14;  % Diffusion coeff for solid in neg. electrode, [m^2/s]
p.D_s_p0 = 1e-13;  % Diffusion coeff for solid in pos. electrode, [m^2/s]

% Diffusional conductivity in electrolyte
p.dactivity = 0;

p.brug = 1.5;       % Bruggeman porosity

% Conductivity of solid
p.sig_n = 100;    % Conductivity of solid in neg. electrode, [1/Ohms*m]
p.sig_p = 10;    % Conductivity of solid in pos. electrode, [1/Ohms*m]

% Miscellaneous
p.t_plus = 0.4;       % Transference number
p.Faraday = 96487;    % Faraday's constant, [Coulumbs/mol]
p.Area = 1;           % Electrode current collector area [m^2]

%% Kinetic Params
p.R = 8.314472;       % Gas constant, [J/mol-K]

p.alph = 0.5;         % Charge transfer coefficients

p.R_f_n = 1e-3;       % Resistivity of SEI layer, [Ohms*m^2]
p.R_f_p = 0;       % Resistivity of SEI layer, [Ohms*m^2]
p.R_c = 0;         % Contact Resistance/Current Collector Resistance, [Ohms-m^2]

% Nominal Reaction rates
p.k_n0 = 1e-5;  % Reaction rate in neg. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]
p.k_p0 = 3e-7; % Reaction rate in pos. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]

%% Thermodynamic Params

% Thermal dynamics
p.C_p = 2000;   % Heat capacity, [J/kg-K]
p.h = 0.36;   % Heat transfer coefficient, [W/K-m^2] 0

% Ambient Temperature
p.T_amb = 298.15; % [K]

% Activation Energies
% Taken from Zhang et al (2014) [Harbin]
% http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
% All units are [J/mol]
p.E.kn = 37.48e3;
p.E.kp = 39.57e3;
p.E.Dsn = 42.77e3;
p.E.Dsp = 18.55e3;
p.E.De = 37.04e3;
p.E.kappa_e = 34.70e3;

% Reference temperature
p.T_ref = 298.15; %[K]

%% Concentrations
% Maxima based on DUALFOIL 
% line 588 in DUALFOIL Fortran code

p.c_s_n_max = 3.6e3 * 372 * 1800 / p.Faraday;   % Max concentration in anode, [mol/m^3]
%p.c_s_n_max = 3.6e3 * 372 * 2260 / p.Faraday;   % Max concentration in anode, [mol/m^3]

%p.c_s_p_max = 3.6e3 * 247 * 5010 / p.Faraday;    % Max concentration in cathode, [mol/m^3]
p.c_s_p_max = 3.6e3 * 274 * 5010 / p.Faraday;    % Max concentration in cathode, [mol/m^3]

p.n_Li_s = 2.5; %2.781;        % Total moles of lithium in solid phase [mol]
p.c_e = 1e3;              % Fixed electrolyte concentration for SPM, [mol/m^3]

%% Cutoff voltages
% p.volt_max = 4.1113; %4.7;
p.volt_min = 3.105; %2.6;
p.volt_max = 4.5;
% p.volt_min = 2.5;

%% Discretization parameters
% Discrete time step
p.delta_t = 1;

% Pade Order
p.PadeOrder = 3;

% Finite difference points along x-coordinate
% p.Nxn = 70;
% p.Nxs = 35;
% p.Nxp = 70;
p.Nxn = 6;
p.Nxs = 3;
p.Nxp = 6;
p.Nx = p.Nxn+p.Nxs+p.Nxp;

p.delta_x_n = 1 / p.Nxn;
p.delta_x_s = 1 / p.Nxs;
p.delta_x_p = 1 / p.Nxp;

%% EHM
p.nu = -0.585340268368373; % SOC+ = nu*SOC- + mu
p.miu = 0.976220713479772;

%Factor between input u (= -3*J/(p.R_s_n * p.c_s_n_max)) and current I [A/m^2]
%pore wall molar flux J = I/(p.Faraday * p.a_s_n * p.L_n) --> Chaturvedi, equation 40 (initial ref on EChM)
p.aa0ne = -3 / (p.R_s_n * p.L_n * p.a_s_n * p.Faraday * p.c_s_n_max); 
p.aa1n = - (2*p.R_s_n) /(7 * p.D_s_n0 * p.L_n * p.a_s_n * p.Faraday * p.c_s_n_max); 
p.bb1n = 1;
p.bb2n = p.R_s_n^2 / (35*p.D_s_n0);
p.Beta = 1 - p.bb2n*p.aa0ne/p.aa1n;
p.gnh  = p.aa0ne/p.aa1n - p.bb2n*p.aa0ne^2/p.aa1n^2; %g

