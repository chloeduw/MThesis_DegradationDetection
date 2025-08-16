function [h,H] = JacobianMeasMat(p)

    syms SOCn CSCn Cur
    CSCp = p.nu*SOCn + p.miu; %=SOCp

    T = p.T_ref; 

    a = p.R*T/(p.alph*p.Faraday);
    b_p = -Cur/(2*p.a_s_p*p.L_p);
    c_p = p.k_p0*p.c_s_p_max*sqrt(p.c_e);
    i0_p = c_p*(CSCp*(1-CSCp)).^p.alph;
    eta_p = a*asinh(b_p/i0_p);
    deta_p = p.nu*deta(a,b_p,c_p,CSCp); 

    b_n = Cur/(2*p.a_s_n*p.L_n);
    c_n = p.k_n0*p.c_s_n_max*sqrt(p.c_e);
    i0_n = c_n*(CSCn*(1-CSCn)).^p.alph;
    eta_n = a*asinh(b_n/i0_n);
    deta_n = deta(a,b_n,c_n,CSCn); 
    
    [Up, dUp] = refPotentialCathode(p,CSCp); %dUp = dUp/dCSCp*1/c_s_n_max
    [Un, dUn] = refPotentialAnode(p,CSCn); 

    Rf = p.R_f_p/(p.L_p*p.a_s_p)-p.R_f_n/(p.L_n*p.a_s_n);

    h = eta_p - eta_n + Up - Un + Rf*Cur;

    Dhx1 = deta_p + p.nu*dUp*p.c_s_p_max; %*c_s_p_max bc see comment line 21
    Dhx2 = -deta_n - dUn*p.c_s_n_max; %*c_s_p_max bc see comment line 21
    H = [Dhx1 Dhx2];
end