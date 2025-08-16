function [x_out, y_out,Vsplit] = EHM_fcn_model(p,I,x0)%, y0)
Ts = 1;

aa0ne = -3 / (p.R_s_n * p.L_n * p.a_s_n * p.Faraday * p.c_s_n_max); %Facteur entre l'entrée u et le courant I
%molar flux J = I/(p.Faraday * p.a_s_n * p.L_n) --> equation 40 de Chaturvedi (ref initiale sur modele electrochimique)
%entrée u = 3*J/(p.R_s_n * p.c_s_n_max)

SOCnh = x0(1,1); %états à l'instant k
CSCnh = x0(2,1);

%aa0n = -3 / (p.R_s_n * p.L_n * p.a_s_n * p.Faraday * p.c_s_n_max);
aa1n = - (2*p.R_s_n) /(7 * p.D_s_n0 * p.L_n * p.a_s_n * p.Faraday * p.c_s_n_max); 
bb1n = 1;
bb2n = p.R_s_n^2 / (35*p.D_s_n0);
Ben = 1 - bb2n*aa0ne/aa1n; %Beta

%Cen1h = x0(3);
%Cen2h = x0(4);

%T1h = x0(5);
%T2h = x0(6);

%volt = y0(1);

%% EKF implementation
%clc
%count=count+1

    % EHM MODEL FOR ELECTRODE - TIME VARYING
    gnh  = aa0ne/aa1n - bb2n*aa0ne^2/aa1n^2; %g
    
    
    % EHM MODEL FOR ELECTROLYTE - TIME VARYING
    % Params/Matrices
    %[Ae_raw,Be_raw,~] = ce_comp(p,T1h);   %Arrh
    %Ae_ct = Ae_raw;
    %Be_ct = Be_raw;
    
    SOCn_prv = SOCnh;
    xsp_prv = EHMlinapprox(p,SOCn_prv);
    %SOCp_prv = xsp_prv(1,:);
    %UpT = refPotentialCathode(p,SOCp_prv);
    %UnT = refPotentialAnode(p,SOCn_prv);        

    % State Model
    f1 = SOCnh + Ts*aa0ne*I;
    f2 = Ts*gnh/(Ben*(1-Ben))*SOCnh + (1-Ts*gnh/(Ben*(1-Ben)))*CSCnh + Ts*aa0ne/(1-Ben)*I;

    f = [f1;f2];
    % EKF estimation
    x_out = f;

    % Output Model
    xp = EHMlinapprox(p,x_out(1,1)); %SOC+,CSC+
    x_out = [x_out;xp];
    
    [yhc,Vsplit] = output_ehm(p,I,xp(2),x_out(2));
    y_out = [x_out(2)]; % [yhc];

    
    
    