function [V,Vsplit] = output_ehm(p,I,CSCp,CSCn)

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
Rf   =  p.R_f_p/(p.a_s_p * p.L_p) - p.R_f_n/(p.a_s_n * p.L_n); %Solid-electrolyte interphase film resistance

%Equilibrium potential (OCV)
Up = refPotentialCathode(p,CSCp );
Un = refPotentialAnode(p,CSCn );

V = etap - etan + Up - Un + Rf*I; %Equation 12 EHM
Vsplit = [etap; etan; Up; Un];

end