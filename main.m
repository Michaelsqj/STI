%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Create phantoms      %
%  9 susceptibility volumes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gen_phantom('chi_phantom');  % Save phantom in .data folder
% Content: 'Params' (B0, gamma, fov, sizeVol, voxSize, pathname, filename_chi)
%          'chi', 'BrainMask', 'chiavg', 'chiani'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Forward                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STI_forward('chi_phantom');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Inverse                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STI_inverse(STIParams, maxit, tol, alpha, beta);