%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab code for the comparison between the behaviour of a fresh and 
%   degraded (4000-cycle-old) cell
% Author: Chloé Duwaerts, 2024-2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; %close all;

%% Gene current
run params_LCO.m;

% Calculate C-rate in terms of [A/m^2] using low/high voltage cutoffs
[cn_low,cp_low] = init_cs(p,p.volt_min); %Initial Solid Concentrations from Voltage
[cn_high,cp_high] = init_cs(p,p.volt_max);
Delta_cn = cn_high-cn_low;
Delta_cp = cp_low-cp_high;
% Current
OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600); %current in [A/m^2]
capa = 2.3; %capacity of the cell[Ah]
A_e = capa/OneC; % electrode area [m^2]

%%%%%%%%%%%%%%% CONSTANT C-RATE %%%%%%%%%%%%%%%%%%%%%%%%%
p.delta_t = 1;
Ts = p.delta_t;
simTime = 5000;
Crate = 1;
t = (0:Ts:simTime*Crate-1)'; 
I_ch = -OneC/Crate; % Current [A/m^2] when charging for 1 cell

%%CI
V0 = 3.5; %3.7 for charge in ~1800; 3.5 for ~2600; 3.608 for ~2300

[csn0,csp0] = init_cs(p,V0); %Concentrations initiales dans les électrodes

CI = csn0/p.c_s_n_max; 

%% NEW cell
% Simulator 
xsim = zeros(2,simTime);  %SOC- et CSC-
xsim(1:2,1) = [CI; CI];

xp = zeros(2,simTime);
xp(:,1) = EHMlinapprox(p,xsim(1,1)); %Initial condition for SOC+, CSC+
xsim(3:4,1) = xp(1,1);

ysim = zeros(1,simTime-1); %V
[ysim(1,1),Vsplit(:,1)] = output_ehm(p,I_ch,xsim(4,1),xsim(2,1)); %Initial condition for V

I = zeros(1,simTime)';
k = 1;
previous = 0;
while k < simTime-1
    while previous < 0.98 && k < simTime-1
        [xsim(:,k+1), ysim(:,k),Vsplit(:,k)] = EHM_fcn_model(p,I_ch,xsim(:,k));
        I(k)=I_ch;
        previous = xsim(2,k); 
        k = k+1;
    end
    previousk = k;
    while k < simTime-1 
        [xsim(:,k+1), ysim(:,k),Vsplit(:,k)] = EHM_fcn_model(p,-I_ch,xsim(:,k));
        I(k)=0;
        previous = xsim(2,k); 
        k = k+1;
    end
    previousk = k; 
end

simin = [t, I*A_e];
tt = t(1):length(I)-1;

% run Simulink model
model = 'TestEHMCell.slx';
out = sim(model);

%% DEGRADED cell
% Init cell 4000 cycles
p_Aged = p;
p_Aged.epsilon_s_n = (1-0.075)*p.epsilon_s_n; %Delta q_eps^n = L*eps_s^n*c_s,max*A_c
p_Aged.epsilon_s_p = (1-0.049)*p.epsilon_s_p; %Delta q_eps^p = L*eps_s^p*c_s,max*A_c
p_Aged.D_s_n0 = (1-0.28)*p.D_s_n0;            %Delta g_D = 1/tau = D_s/R_s^2
p_Aged.n_Li_s = (1-0.13)*p.n_Li_s;            %Delta n_Li
p_Aged.R_f_n = 1.42*p.R_f_n;                  %Delta r_R = R_f*R_s*A_c*c_s,max
p_Aged = updateDependentParams(p_Aged);

[csn0_Aged,csp0_Aged] = init_cs(p_Aged,V0); %Concentrations initiales dans les électrodes

CI_Aged = csn0_Aged/p_Aged.c_s_n_max; % 1- for disch

% Simulateur 
xsim = zeros(2,simTime);  %SOC- et CSC-
xsim(1:2,1) = [CI_Aged; CI_Aged];

xp = zeros(2,simTime);
xp(:,1) = EHMlinapprox(p_Aged,xsim(1,1)); %Initial condition for SOC+, CSC+
xsim(3:4,1) = xp(1,1);

ysim = zeros(1,simTime-1); %V
[ysim(1,1),Vsplit(:,1)] = output_ehm(p_Aged,I_ch,xsim(4,1),xsim(2,1)); %Initial condition for V

I_Aged = zeros(1,simTime)';
k = 1;
previous = 0;
while k < simTime-1
    while previous < 0.98 && k < simTime-1%0.8
        [xsim(:,k+1), ysim(:,k),Vsplit(:,k)] = EHM_fcn_model(p_Aged,I_ch,xsim(:,k));
        I_Aged(k)=I_ch; %want to apply same current as new cell to see different behaviour
        previous = xsim(2,k); %2
        k = k+1;
    end
    previousk = k;
    while k < simTime-1 %0.045 %k < (previousk+500) &&
        [xsim(:,k+1), ysim(:,k),Vsplit(:,k)] = EHM_fcn_model(p_Aged,-I_ch,xsim(:,k));
        I_Aged(k)=0;
        previous = xsim(2,k); %2
        k = k+1;
    end
    previousk = k; 
end

simin_Aged = [t, I_Aged*A_e];
tt_Aged = t(1):length(I)-1;

% run Simulink model
model = 'TestEHMCellDegraded.slx';
out_Aged = sim(model);

%% Figures
figure()
subplot(211)
plot(t(1:3000),out.Vsim.Data(1:3000),t(1:3000),out_Aged.Vsim.Data(1:3000),LineWidth=1.5)
title('Output voltage')
legend('New','Degraded','Location','northwest')
xlabel('time [s]')
ylabel('[V]')
axis tight;
subplot(212)
plot(t(1:3000),I(1:3000)*A_e,t(1:3000),I_Aged(1:3000)*A_e,LineWidth=1.5)
legend('New','Degraded','Location','northwest')
title('Input current')
xlabel('time [s]')
ylabel('[A]')
axis tight;

figure()
subplot(211)
plot(t(1:3000),out.SOCn.Data(1:3000),t(1:3000),out_Aged.SOCn.Data(1:3000),LineWidth=1.5)
legend('New','Degraded','Location','northwest')
xlabel('time [s]')
ylabel('[-]')
title('SOC')
axis tight;
subplot(212)
plot(t(1:3000),out.CSCn.Data(1:3000),t(1:3000),out_Aged.CSCn.Data(1:3000),LineWidth=1.5)
legend('New','Degraded','Location','northwest')
xlabel('time [s]')
ylabel('[-]')
title('CSC')
axis tight;

%%
% e_SOC_Aged = -out.SOCn.Data(1:end-3)+out_Aged.SOCn.Data(1:end-3);
% e_CSC_Aged = -out.CSCn.Data(1:end-3)+out_Aged.CSCn.Data(1:end-3);
% e_y_Aged = -out.Vsim.Data(1:end-3)+out_Aged.Vsim.Data(1:end-3);
% 
% figure()
% subplot(311)
% plot(t(1:end-2),e_SOC_Aged,LineWidth=1.5)
% xlabel('time [s]')
% ylabel('[-]')
% title('\Delta SOC')
% axis tight;
% subplot(312)
% plot(t(1:end-2),e_CSC_Aged,LineWidth=1.5)
% xlabel('time [s]')
% ylabel('[-]')
% title('\Delta CSC')
% axis tight;
% subplot(313)
% plot(t(1:end-2),e_y_Aged,LineWidth=1.5)
% xlabel('time [s]')
% ylabel('[V]')
% title('\Delta y')
% axis tight;