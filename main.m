%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Create phantoms      %
%  9 susceptibility volumes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gen_phantom();  % Save phantom in .data folder
load('data\chi_phantom.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Forward                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STI_forward();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Inverse                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STI_inverse(STIParams, maxit, tol, alpha, beta);