function [V,Vsplit] = output_ehm(p,I,CSCp,CSCn)
                            %p,ysim,xp(2),x(2)

Cen = p.c_e; %Avant =10
Cep = p.c_e;% avant = 10;
T1 = p.T_ref;
k_n = p.k_n0;% * exp(p.E_kn/p.R*(1/p.T_ref - 1./T1)); reaction rate
k_p = p.k_p0;% * exp(p.E_kp/p.R*(1/p.T_ref - 1./T1)); reaction rate

Gan = 2*p.L_n*p.a_s_n*k_n*p.c_s_n_max*sqrt(Cen);
Gap = 2*p.L_p*p.a_s_p*k_p*p.c_s_p_max*sqrt(Cep);
beta0 = p.R*T1 / (p.alph*p.Faraday);

%Surface overpotential
etap = beta0*asinh((1./(Gap .* ((1 - CSCp ) .* CSCp ).^p.alph)).*-I); %denom in asin = j_n,0
%eq 14 de EHM, avec J = I/(F*a*L) eq 40 Chaturvedi
etan = beta0*asinh((1./(Gan .* ((1 - CSCn ) .* CSCn ).^p.alph)).* I);

%SEI
Rf   =  p.R_f_n/(p.a_s_n * p.L_n) + p.R_f_p/(p.a_s_p * p.L_p); %Solid-electrolyte interphase film resistance
%Pas certains de la justesse de l'équation, pour moi cela devrait être 
%Rp = Rf*u = p.R_f_n*3*I/(p.R_s_n * p.L_n * p.a_s_n * p.F * p.c_s_n_max) + la même
%pour l'électrode positive
%Car eq12 EHM : Rp = Rf*u; u = 3J/R*Cmax et J = I/F*a*L
%Pas de problème ici car p.R_f_n = p.R_f_p = 0;
%Rp = Rf*I; 
Rp = p.R_f_n*3*I/(p.R_s_n * p.L_n * p.a_s_n * p.Faraday * p.c_s_n_max)+p.R_f_p*3*I/(p.R_s_p * p.L_p * p.a_s_p * p.Faraday * p.c_s_p_max);

%Equilibrium potential (OCV)
Up = refPotentialCathode(p,CSCp );
Un = refPotentialAnode(p,CSCn );

V = etap - etan + Up - Un - Rp; %Equation 12 EHM
Vsplit = [etap; etan; Up; Un];

end