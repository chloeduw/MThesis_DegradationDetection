function [x_out, y_out,Vsplit] = EHM_fcn_model(p,I,x0)%, y0)
    Ts = 1; %sampling time
    
    SOCnh = x0(1,1); %états à l'instant k
    CSCnh = x0(2,1);

%% EKF implementation

    % EHM MODEL FOR ELECTRODE - TIME VARYING
    
    SOCn_prv = SOCnh;
    xsp_prv = EHMlinapprox(p,SOCn_prv);      

    % State Model
    f1 = SOCnh + Ts*p.aa0ne*I;
    f2 = Ts*p.gnh/(p.Beta*(1-p.Beta))*SOCnh + (1-Ts*p.gnh/(p.Beta*(1-p.Beta)))*CSCnh + Ts*p.aa0ne/(1-p.Beta)*I;

    f = [f1;f2];
    % EKF estimation
    x_out = f;

    % Output Model
    xp = EHMlinapprox(p,x_out(1,1)); %SOC+,CSC+
    x_out = [x_out;xp];
    
    [yhc,Vsplit] = output_ehm(p,I,xp(2),x_out(2));
    y_out = yhc; 

    
    
    