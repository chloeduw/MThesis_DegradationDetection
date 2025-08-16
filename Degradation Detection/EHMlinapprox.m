function [xp_det] = EHMlinapprox(p,SOCn)
%% EHM MODEL
% Positive electrode linear approximation
SOCp = (p.nu*SOCn + p.miu);
CSCp = SOCp;

xp_det = [SOCp; CSCp];



