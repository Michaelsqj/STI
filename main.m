%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Create phantoms      %
%  9 susceptibility volumes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfolder('data')
    mkdir('data');
end
gen_phantom('data/chi_phantom');  % Save phantom in .data folder
% Content: 'Params' (B0, gamma, fov, sizeVol, voxSize, pathname, filename_chi)
%          'chi', 'BrainMask', 'chiavg', 'chiani'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    Forward     %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STI_forward('data/chi_phantom');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    Inverse     %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       STIParams includes fields:
%           OriNum: Orientation Number
%           fov, sizeVol, VoxSize: field of view, size of Volume, voxel size
%           filenameFreqMap: cell array containing the filenames of frequency maps
%           maskMSA: White Matter Mask
%           filenameQSM: 
%
%       maxit: maximum number of iteration
%       tol: convergence tolerence for LSQR
%       alpha: regularization parameter
%       beta: regularization parameter 
alpha = 10;
beta = 0.1;
tol = 1e-4;
maxit = 50;
STIParams_filename = 'data/STIParams.mat';
[chi_r,flag,relres,iter,resvec,lsvec]=STI_inverse(STIParams_filename, maxit, tol, alpha);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    Evaluate    %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chi_t = load('data/chi_phantom.mat').chi;
[RE_MMS, RE_MSA, AE] = evaluate(chi_r, chi_t);