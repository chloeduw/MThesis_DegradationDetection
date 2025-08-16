%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab code for the validation of the EKF on a fresh cell
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

% EHM Simulator 
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
        previous = xsim(2,k); %control charge/discharge cycle with CSC
        k = k+1;
    end
    while previous > 0.1 && k < simTime-1 
        [xsim(:,k+1), ysim(:,k),Vsplit(:,k)] = EHM_fcn_model(p,-I_ch,xsim(:,k));
        I(k)=-I_ch;
        previous = xsim(2,k); 
        k = k+1;
    end
end

simin = [t, I*A_e];
tt = t(1):length(I)-1;

% run Simulink model
model = 'TestEHMCell.slx';
out = sim(model);

%% EKF validation
v = sqrt(1e-3);
y = out.Vsim.Data(1:end-1) + randn(simTime,1)*v; %ysim + noise

% Estimation of new cell with params new cell
[xpost_new,y_new,innov_new,S_new,traceSigmapost_new] = ekf(p,simTime,y,I,0.5*ones(2,1),v,1);

%% Figures
figure()
subplot(211)
plot(t(1:end-1),out.SOCn.Data(1:end-2),LineWidth=1.5)
hold on
plot(t(1:end-1),xpost_new(1,1:end-1),LineWidth=1.5)
xlabel('time [s]')
ylabel('[-]')
legend('Real','EKF')
title('SOC estimation')
axis tight;
subplot(212)
plot(t(1:end-1),out.CSCn.Data(1:end-2),LineWidth=1.5)
hold on
plot(t(1:end-1),xpost_new(2,1:end-1),LineWidth=1.5)
xlabel('time [s]')
ylabel('[-]')
legend('Real','EKF')
title('CSC estimation')
axis tight;

figure()
subplot(211)
plot(t(1:end-2),y(1:end-2),t(1:end-2),out.Vsim.Data(1:end-3),t(1:end-2),y_new(1:end-2),LineWidth=1.5)
legend('Measured','Real','EKF')
xlabel('time [s]')
ylabel('[V]')
title('y')
axis tight;
subplot(212)
plot(t(1:end-2),I(1:end-2)*A_e,LineWidth=1.5)
xlabel('time [s]')
ylabel('[A]')
title('I')
axis tight;

Delta_SOC_new = xpost_new(1,1:end-1) - out.SOCn.Data(1:end-2)';
Delta_CSC_new = xpost_new(2,1:end-1) - out.CSCn.Data(1:end-2)';
Delta_y_new = y_new(1,1:end-1) - out.Vsim.Data(1:end-2)';

figure()
subplot(311)
plot(t(1000:end-1),Delta_SOC_new(1000:end),LineWidth=1.5)
xlabel('time [s]')
ylabel('[-]')
ylim([-0.1 0.1])
title('\Delta SOC')
axis tight;
subplot(312)
plot(t(1000:end-1),Delta_CSC_new(1000:end),LineWidth=1.5)
xlabel('time [s]')
ylabel('[-]')
ylim([-0.1 0.1])
title('\Delta CSC')
axis tight;
subplot(313)
plot(t(1000:end-1),Delta_y_new(1000:end),LineWidth=1.5)
xlabel('time [s]')
ylabel('[V]')
ylim([-0.1 0.1])
title('\Delta y')
axis tight;

%% histogram
start=1000;
[var_new1,mean_new1] = var(innov_new(start:end));
var_th = mean(S_new(1000:end));

fprintf(1,'Innovation new cell (new EKF): variance = %2.2e, mean = %2.2e',var_new1,mean_new1);
disp(".");

figure()
histogram(innov_new(start:end),'Normalization','probability',BinWidth=1e-2)
title({'Innovation distribution',['(\sigma^2 = ',num2str(var_new1,'%.3e'),', \mu = ',num2str(mean_new1,'%.3e'),')']})