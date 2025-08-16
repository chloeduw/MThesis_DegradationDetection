clear all
%clc
%close all

%load P_V0_Luis_Daniel.mat
%load P_Scott_Moura.mat
%load current_profile.mat
run params_LCO.m

% Calculate C-rate in terms of [A/m^2] using low/high voltage cutoffs
[cn_low,cp_low] = init_cs(p,p.volt_min); %Initial Solid Concentrations from Voltage
[cn_high,cp_high] = init_cs(p,p.volt_max);
Delta_cn = cn_high-cn_low;
Delta_cp = cp_low-cp_high;
OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600);

%%%%%%%%%%%%%%% CONSTANT C-RATE %%%%%%%%%%%%%%%%%%%%%%%%%
p.delta_t = 1;
% t = 0:p.delta_t:(120);
t = 0:p.delta_t:2600*1; %*100 si C/100
I = -1*OneC*ones(size(t))/1; %/100 si C/100

I = [I zeros(1,1)];
t = t(1):length(I)-1;
% break

V0 = 3.5;
p.nu = -0.585340268368373;
p.miu = 0.976220713479772;
p.R_f_p = 0;
p.R_f_n = p.R_f_p;

NT = length(t);

[csn0,csp0] = init_cs(p,V0); %Concentrations initiales dans les électrodes

xsim = zeros(2,NT);  %SOC- et CSC-
xsim(1:2,1) = [csn0/p.c_s_n_max; csn0/p.c_s_n_max]; %Initial conditions

xp = zeros(2,NT);
xp(:,1) = EHMlinapprox(p,xsim(1,1)); %Initial condition for SOC+
xsim(3:4,1) = xp(1,1);

ysim = zeros(1,NT-1); %V
[ysim(1,1),Vsplit(:,1)] = output_ehm(p,I(1),xsim(4,1),xsim(2,1)); %Initial condition for V

for k = 1:(NT-1)
    
    [xsim(:,k+1), ysim(:,k),Vsplit(:,k)] = EHM_fcn_model(p,I(k),xsim(:,k));
    
    fprintf(1,'t: %3.0f sec | Cur: %2.2f C-rate | Volt: %2.3fV | [%1.3f/%1.3f/%1.3f/%1.3f] \n',...
        t(k),I(k)/OneC, ysim(1,k), xsim(1,k),xsim(2,k),xsim(3,k),xsim(4,k));

end

%% Figures
figure(21)
subplot(211)
plot(t,I/OneC,LineWidth=1.5)
title('current')
axis tight;
subplot(212)
plot(t(1:length(ysim)),ysim(:),LineWidth=1.5)
title('voltage')
axis tight;

figure()
subplot(211)
plot(t,xsim(1,:),t,xsim(3,:),'--',LineWidth=1.5)
legend('SOC-','SOC+')
title('SOCs')
axis tight;
subplot(212)
plot(t,xsim(2,:),t,xsim(4,:),'--',LineWidth=1.5)
legend('CSC-','CSC+')
title('CSCs')
axis tight;

%% KF
j = -3 / (p.R_s_n * p.L_n * p.a_s_n * p.Faraday * p.c_s_n_max)*I(2:end);
u = [j; Vsplit(4,:)];
Beta = 7/10;
gnh = 147/20*p.D_s_n0/p.R_s_n^2;
A = [1 0; gnh/(Beta*(1-Beta)) 1-gnh/(Beta*(1-Beta))];
B = 1/(p.c_s_n_max*p.Faraday*p.epsilon_e_s*p.L_n)*[1 1/(1-Beta)];
C = [0 1];
Rinit = 1e-6; %cov v -> y -> V
Qinit = [1e-3 0; 0 1e-4]; %cov w -> x -> SOC and CSC
[xpred,xpost,deltax] = kf(NT,A,B,C,Rinit,Qinit,ysim,u);

%% Fig KF
figure()
subplot(211)
plot(t,xpred(1,:),t,xpost(1,:),'--',LineWidth=1.5)
legend('SOC_{pred}','SOC_{post}')
title('SOCs')
axis tight;
subplot(212)
plot(t,xpred(2,:),t,xpost(2,:),'--',LineWidth=1.5)
legend('CSC_{pred}','CSC_{post}')
title('CSCs')
axis tight;

figure()
subplot(211)
plot(t,xsim(1,:),t,xpost(1,:),'--',LineWidth=1.5)
legend('SOC_{EHM-}','SOC_{KF}')
title('SOCs')
axis tight;
subplot(212)
plot(t,xsim(2,:),t,xpost(2,:),'--',LineWidth=1.5)
legend('CSC_{EHM-}','CSC_{post}')
title('CSCs')
axis tight;