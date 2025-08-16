function p = updateDependentParams(p)
    % Specific interfacial surface area
    p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;
    p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;

    % Cell mass [kg/m^2]
    p.m_n = p.L_n * (p.rho_e*p.epsilon_e_n + p.rho_sn*p.epsilon_s_n + p.rho_f*p.epsilon_f_n);
    p.m_p = p.L_p * (p.rho_e*p.epsilon_e_p + p.rho_sp*p.epsilon_s_p + p.rho_f*p.epsilon_f_p);
    
    % Lumped density [kg/m^2]
    p.rho_avg = p.m_n + p.m_s + p.m_p + p.m_cc;

    % Other
    p.aa0ne = -3 / (p.R_s_n * p.L_n * p.a_s_n * p.Faraday * p.c_s_n_max); 
    p.aa1n = - (2*p.R_s_n) /(7 * p.D_s_n0 * p.L_n * p.a_s_n * p.Faraday * p.c_s_n_max);
    p.bb2n = p.R_s_n^2 / (35*p.D_s_n0);
    p.Beta = 1 - p.bb2n*p.aa0ne/p.aa1n;
    p.gnh  = p.aa0ne/p.aa1n - p.bb2n*p.aa0ne^2/p.aa1n^2; 
end