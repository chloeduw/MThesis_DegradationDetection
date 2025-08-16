clear all; clc; %close all;

%% Gene current
run params_LCO.m;

% Calculate C-rate in terms of [A/m^2] using low/high voltage cutoffs
[cn_low,cp_low] = init_cs(p,p.volt_min); %Initial Solid Concentrations from Voltage
[cn_high,cp_high] = init_cs(p,p.volt_max);
Delta_cn = cn_high-cn_low;
Delta_cp = cp_low-cp_high;
% Current
OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600);%[A/m^2]
capa = 2.3; %[Ah]
A = capa/OneC/(p.a_s_n*p.L_n); % electrode area [m^2]

%%%%%%%%%%%%%%% CONSTANT C-RATE %%%%%%%%%%%%%%%%%%%%%%%%%
p.delta_t = 1;
Ts = p.delta_t;
simTime = 10000; %25500; %2505*19; 10chdisch %2505*8/15 (4/8chdisch), 2450 disch diff, 2505 ch diff, 2605 ch same
Crate = 1;
t = (0:Ts:simTime*Crate-1)'; %2600*10 1800*1  %*100 si C/100
I_ch = -OneC/Crate*A*p.L_n*p.a_s_n; % Current [A] when charging for 1 cell

% For TestEHMCell
I = I_ch*ones(size(t)); %/10 bc 10h for 1 charge, C/10 (/100 if C/100)
I(1:5) = 0;
simin = [t, I]; %I*2

t = t(1):length(I)-1;

Vmax = p.volt_max;
Vmin = p.volt_min;

%% CI new cell
V0 = 3.5; %3.7 for 1800; 3.5 for 2600; 3.608 for 2300

[csn0,csp0] = init_cs(p,V0); %Concentrations initiales dans les électrodes

CI = csn0/p.c_s_n_max; %1- for disch

%% run Simulink model
% SOC |b| 0.2 and 0.9 -> hysteresis min : p.volt_min+0.34
% SOC |b| 0.1 and 0.9 -> hysteresis min : p.volt_min+0.05
% SOC |b| 0.35 and 0.9 -> hysteresis min : p.volt_min+0.435
model = 'TestParallel.slx'; %TestAgeingImpact_CCCV.slx
out = sim(model);