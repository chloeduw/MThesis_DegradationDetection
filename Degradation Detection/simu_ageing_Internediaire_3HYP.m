%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab code for the degradation monitoring test of a fresh cell in
%   parallel to a degraded cell using 3 EKFs.
% Run the code part by part: Gene current, EKF init param, 
%   one of the degraded cells param, run simulink model, run EKFs, 
%   run conditional probabilities calculation
% Author: Chloé Duwaerts, 2024-2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; %close all;

%% Gene current
run params_LCO.m; %p = params of the fresh cell

% Calculate C-rate in terms of [A/m^2] using low/high voltage cutoffs
[cn_low,cp_low] = init_cs(p,p.volt_min); %Initial Solid Concentrations from Voltage
[cn_high,cp_high] = init_cs(p,p.volt_max);
Delta_cn = cn_high-cn_low;
Delta_cp = cp_low-cp_high;
% Current
OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600);% current in [A/m^2]
capa = 2.3; %capacity of the cell [Ah]
A_e = capa/OneC; % electrode area [m^2]

%%%%%%%%%%%%%%% CONSTANT C-RATE %%%%%%%%%%%%%%%%%%%%%%%%%
p.delta_t = 1;
Ts = p.delta_t;
simTime = 3000; 
Crate = 1;
t = (0:Ts:simTime*Crate-1)';
I_ch = -OneC/Crate*A_e; % Current [A] when charging for 1 cell

% For TestEHMCell
I = I_ch*ones(size(t)); 
I(1:5) = 0;
simin = [t, I]; 

t = t(1):length(I)-1;

Vmax = p.volt_max;
Vmin = p.volt_min;

%%CI new cell
V0 = 3.5; %3.7 for charge in ~1800; 3.5 for ~2600; 3.608 for ~2300

[csn0,csp0] = init_cs(p,V0); %Concentrations initiales dans les électrodes

CI = csn0/p.c_s_n_max;

%% Init of the parameters for the EKFs
% 4000-cycle-old cell
p_Aged = p;
p_Aged.epsilon_s_n = (1-0.075)*p.epsilon_s_n; %Delta q_eps^n = L*eps_s^n*c_s,max*A_c
p_Aged.epsilon_s_p = (1-0.049)*p.epsilon_s_p; %Delta q_eps^p = L*eps_s^p*c_s,max*A_c
p_Aged.D_s_n0 = (1-0.28)*p.D_s_n0;            %Delta g_D = 1/tau = D_s/R_s^2
p_Aged.n_Li_s = (1-0.13)*p.n_Li_s;            %Delta n_Li
p_Aged.R_f_n = 1.42*p.R_f_n;                  %Delta r_R = R_f*R_s*A_c*c_s,max
p_Aged = updateDependentParams(p_Aged);

[csn0_Aged,csp0_Aged] = init_cs(p_Aged,V0); %Concentrations initiales dans les électrodes

CI_Aged = csn0_Aged/p_Aged.c_s_n_max;

% 2000-cycle-old cell
p_Int = p;
p_Int.epsilon_s_n = (1-0.04)*p.epsilon_s_n;
p_Int.epsilon_s_p = (1-0.03)*p.epsilon_s_p;
p_Int.D_s_n0 = (1-0.18)*p.D_s_n0;
p_Int.n_Li_s = (1-0.08)*p.n_Li_s;
p_Int.R_f_n = 1.28*p.R_f_n;
p_Int = updateDependentParams(p_Int);

[csn0_Int,csp0_Int] = init_cs(p_Int,V0); %Concentrations initiales dans les électrodes

CI_Int = csn0_Int/p_Int.c_s_n_max; % 1- for disch

%% Init of the parameters of the degraded cell for the diagnosis
% run one of the following degradation levels
cycles = 4000; %degradation level
p_Inc = p;
p_Inc.epsilon_s_n = (1-0.075)*p.epsilon_s_n;
p_Inc.epsilon_s_p = (1-0.049)*p.epsilon_s_p;
p_Inc.D_s_n0 = (1-0.28)*p.D_s_n0;
p_Inc.n_Li_s = (1-0.13)*p.n_Li_s;
p_Inc.R_f_n = 1.42*p.R_f_n;
p_Inc = updateDependentParams(p_Inc);

[csn0_Inc,csp0_Inc] = init_cs(p_Inc,V0); %Concentrations initiales dans les électrodes

CI_Inc = csn0_Inc/p_Inc.c_s_n_max; % 1- for disch

%% Init of the parameters of the degraded cell for the diagnosis
% run one of the following degradation levels
cycles = 3000; %degradation level
p_Inc = p;
p_Inc.epsilon_s_n = (1-0.06)*p.epsilon_s_n;
p_Inc.epsilon_s_p = (1-0.04)*p.epsilon_s_p;
p_Inc.D_s_n0 = (1-0.24)*p.D_s_n0;
p_Inc.n_Li_s = (1-0.11)*p.n_Li_s;
p_Inc.R_f_n = 1.37*p.R_f_n;
p_Inc = updateDependentParams(p_Inc);

[csn0_Inc,csp0_Inc] = init_cs(p_Inc,V0); %Concentrations initiales dans les électrodes

CI_Inc = csn0_Inc/p_Inc.c_s_n_max; % 1- for disch

%% Init of the parameters of the degraded cell for the diagnosis
% run one of the following degradation levels
cycles = 2000; %degradation level
p_Inc = p;
p_Inc.epsilon_s_n = (1-0.04)*p.epsilon_s_n;
p_Inc.epsilon_s_p = (1-0.03)*p.epsilon_s_p;
p_Inc.D_s_n0 = (1-0.18)*p.D_s_n0;
p_Inc.n_Li_s = (1-0.08)*p.n_Li_s;
p_Inc.R_f_n = 1.28*p.R_f_n;
p_Inc = updateDependentParams(p_Inc);

[csn0_Inc,csp0_Inc] = init_cs(p_Inc,V0); %Concentrations initiales dans les électrodes

CI_Inc = csn0_Inc/p_Inc.c_s_n_max; % 1- for disch

%% Init of the parameters of the degraded cell for the diagnosis
% run one of the following degradation levels
cycles = 1400; %degradation level
p_Inc = p;
p_Inc.epsilon_s_n = (1-0.03)*p.epsilon_s_n;
p_Inc.epsilon_s_p = (1-0.025)*p.epsilon_s_p;
p_Inc.D_s_n0 = (1-0.13)*p.D_s_n0;
p_Inc.n_Li_s = (1-0.06)*p.n_Li_s;
p_Inc.R_f_n = 1.19*p.R_f_n;
p_Inc = updateDependentParams(p_Inc);

[csn0_Inc,csp0_Inc] = init_cs(p_Inc,V0); %Concentrations initiales dans les électrodes

CI_Inc = csn0_Inc/p_Inc.c_s_n_max; % 1- for disch

%% Init of the parameters of the degraded cell for the diagnosis
% run one of the following degradation levels
cycles = 1300; %degradation level
p_Inc = p;
p_Inc.epsilon_s_n = (1-0.026)*p.epsilon_s_n;
p_Inc.epsilon_s_p = (1-0.024)*p.epsilon_s_p;
p_Inc.D_s_n0 = (1-0.12)*p.D_s_n0;
p_Inc.n_Li_s = (1-0.055)*p.n_Li_s;
p_Inc.R_f_n = 1.18*p.R_f_n;
p_Inc = updateDependentParams(p_Inc);

[csn0_Inc,csp0_Inc] = init_cs(p_Inc,V0); %Concentrations initiales dans les électrodes

CI_Inc = csn0_Inc/p_Inc.c_s_n_max; % 1- for disch

%% Init of the parameters of the degraded cell for the diagnosis
% run one of the following degradation levels
cycles = 1200; %degradation level
p_Inc = p;
p_Inc.epsilon_s_n = (1-0.024)*p.epsilon_s_n;
p_Inc.epsilon_s_p = (1-0.0225)*p.epsilon_s_p;
p_Inc.D_s_n0 = (1-0.10)*p.D_s_n0;
p_Inc.n_Li_s = (1-0.053)*p.n_Li_s;
p_Inc.R_f_n = 1.165*p.R_f_n;
p_Inc = updateDependentParams(p_Inc);

[csn0_Inc,csp0_Inc] = init_cs(p_Inc,V0); %Concentrations initiales dans les électrodes

CI_Inc = csn0_Inc/p_Inc.c_s_n_max; % 1- for disch

%% Init of the parameters of the degraded cell for the diagnosis
% run one of the following degradation levels
cycles = 1100; %degradation level
p_Inc = p;
p_Inc.epsilon_s_n = (1-0.0232)*p.epsilon_s_n;
p_Inc.epsilon_s_p = (1-0.0213)*p.epsilon_s_p;
p_Inc.D_s_n0 = (1-0.0857)*p.D_s_n0;
p_Inc.n_Li_s = (1-0.05)*p.n_Li_s;
p_Inc.R_f_n = 1.16*p.R_f_n;
p_Inc = updateDependentParams(p_Inc);

[csn0_Inc,csp0_Inc] = init_cs(p_Inc,V0); %Concentrations initiales dans les électrodes

CI_Inc = csn0_Inc/p_Inc.c_s_n_max; % 1- for disch

%% Init of the parameters of the degraded cell for the diagnosis
% run one of the following degradation levels
cycles = 1000; %degradation level
p_Inc = p;
p_Inc.epsilon_s_n = (1-0.02)*p.epsilon_s_n;
p_Inc.epsilon_s_p = (1-0.02)*p.epsilon_s_p;
p_Inc.D_s_n0 = (1-0.07)*p.D_s_n0; 
p_Inc.n_Li_s = (1-0.045)*p.n_Li_s;
p_Inc.R_f_n = 1.15*p.R_f_n;
p_Inc = updateDependentParams(p_Inc);

[csn0_Inc,csp0_Inc] = init_cs(p_Inc,V0); %Concentrations initiales dans les électrodes

CI_Inc = csn0_Inc/p_Inc.c_s_n_max; % 1- for disch

%% Init of the parameters of the degraded cell for the diagnosis
% run one of the following degradation levels
cycles = 800; %degradation level
p_Inc = p;
p_Inc.epsilon_s_n = (1-0.017)*p.epsilon_s_n;
p_Inc.epsilon_s_p = (1-0.018)*p.epsilon_s_p;
p_Inc.D_s_n0 = (1-0.047)*p.D_s_n0; 
p_Inc.n_Li_s = (1-0.0375)*p.n_Li_s;
p_Inc.R_f_n = 1.133*p.R_f_n;
p_Inc = updateDependentParams(p_Inc);

[csn0_Inc,csp0_Inc] = init_cs(p_Inc,V0); %Concentrations initiales dans les électrodes

CI_Inc = csn0_Inc/p_Inc.c_s_n_max; % 1- for disch

%% run Simulink model
% NOW: SOC |b| 0.1 and 0.9 -> hysteresis min : p.volt_min+0.05
% If want to change the SOC limis:
% SOC |b| 0.2 and 0.9 -> hysteresis min : p.volt_min+0.34
% SOC |b| 0.35 and 0.9 -> hysteresis min : p.volt_min+0.435
model = 'TestAgeingImpact_Intermediate_CC_Cycles.slx';
out = sim(model);

%% EKFs
v = sqrt(1e-3);
y = out.Vsim.Data(1:end-1) + randn(simTime,1)*v; %ysim + noise
u1 = out.Isim1.Data(1:end-1)/A_e;% [A/m^2] 
u2 = out.Isim2.Data(1:end-1)/A_e;% [A/m^2]

% Estimation of new cell with params new cell
[xpost_new,y_new,innov_new,S_new,traceSigmapost_new] = ekf(p,simTime,y,u1,0.5*ones(2,1),v,1);
% Estimation of new cell with params 4000-cycle-old cell
[xpost_new_2,y_new_2,innov_new_2,S_new2,traceSigmapost_new2] = ekf(p_Aged,simTime,y,u1,0.5*ones(2,1),v,2);
% Estimation of new cell with params 2000-cycle-old cell
[xpost_new_3,y_new_3,innov_new_3,S_new3,traceSigmapost_new3] = ekf(p_Int,simTime,y,u1,0.5*ones(2,1),v,3);
% Estimation of degraded cell with params new cell
[xpost_Aged_1,y_Aged_1,innov_Aged_1,S_Aged_1,traceSigmapost_Aged_1] = ekf(p,simTime,y,u2,0.5*ones(2,1),v,4);
% Estimation of degraded cell with params 4000-cycle-old cell
[xpost_Aged_2,y_Aged_2,innov_Aged_2,S_Aged_2,traceSigmapost_Aged_2] = ekf(p_Aged,simTime,y,u2,0.5*ones(2,1),v,5);
% Estimation of degraded cell with params 2000-cycle-old cell
[xpost_Aged_3,y_Aged_3,innov_Aged_3,S_Aged_3,traceSigmapost_Aged_3] = ekf(p_Int,simTime,y,u2,0.5*ones(2,1),v,6);

%% Conditional proba calculation 
% 1 = new, 2= 4000, 3 = 2000
f_new = conditionalDensityFct(simTime,innov_new,S_new);
f_new_2 = conditionalDensityFct(simTime,innov_new_2,S_new2);
f_new_3 = conditionalDensityFct(simTime,innov_new_3,S_new3);
f_Aged_1 = conditionalDensityFct(simTime,innov_Aged_1,S_Aged_1);
f_Aged_2 = conditionalDensityFct(simTime,innov_Aged_2,S_Aged_2);
f_Aged_3 = conditionalDensityFct(simTime,innov_Aged_3,S_Aged_3);

start = 2373;
proba_new = zeros(3,simTime); proba_new(:,1:start-1) = 1/2;
proba_Aged = zeros(3,simTime); proba_Aged(:,1:start-1) = 1/2;
denom_new = zeros(1,simTime); denom_Aged = zeros(1,simTime);
for k=start:simTime
    denom_new(k) = f_new(k)*proba_new(1,k-1) + f_new_2(k)*proba_new(2,k-1) + f_new_3(k)*proba_new(3,k-1);
    proba_new(1,k) = f_new(k)*proba_new(1,k-1)/denom_new(k);
    proba_new(2,k) = f_new_2(k)*proba_new(2,k-1)/denom_new(k);
    proba_new(3,k) = f_new_3(k)*proba_new(3,k-1)/denom_new(k);
    denom_Aged(k) = f_Aged_1(k)*proba_Aged(1,k-1) + f_Aged_2(k)*proba_Aged(2,k-1) + f_Aged_3(k)*proba_Aged(3,k-1);
    proba_Aged(1,k) = f_Aged_1(k)*proba_Aged(1,k-1)/denom_Aged(k);
    proba_Aged(2,k) = f_Aged_2(k)*proba_Aged(2,k-1)/denom_Aged(k);
    proba_Aged(3,k) = f_Aged_3(k)*proba_Aged(3,k-1)/denom_Aged(k);
end

figure()
subplot(2,1,1)
plot(t,proba_new(1,:),t,proba_new(2,:),t,proba_new(3,:),"LineWidth",1.5)
legend('EKF new','EKF 4000','EKF 2000')
title('New cell')
xlabel('time [s]')
%xlim([0 10000])
ylabel({'conditional probability','[-]'})
subplot(2,1,2)
plot(t,proba_Aged(1,:),t,proba_Aged(2,:),t,proba_Aged(3,:),"LineWidth",1.5)
legend('EKF new','EKF 4000','EKF 2000')
title('Degraded cell')
xlabel('time [s]')
%xlim([0 10000])
ylabel({'conditional probability','[-]'})
