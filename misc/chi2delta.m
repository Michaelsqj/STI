function delta = chi2delta(chi11, chi12, chi13, chi22, chi23, chi33, Params)
% calculate the relative field change from symmetric susceptibility tensor
% at certain orientation as defined in Params
% INPUT:
% 
% chi_ij:   tensor components in the subject frame
% Params:   data parameter structure
%           should have Params.sizeVol, Params.fov, Params.TAng
%
% OUTPUT:
% delta: relative magnetic field change

% delta(k) = (1/3)*H0*FT{chi}*H0 - H0*k*(k*FT{chi}*H0/k^2)

warning off all

% k-space data of tensor components
chi11k = fftn(chi11);                     
chi12k = fftn(chi12);
chi13k = fftn(chi13);
chi22k = fftn(chi22);
chi23k = fftn(chi23);
chi33k = fftn(chi33);

[a11, a12, a13, a22, a23, a33] = conv_kernel_tensor(Params);

delta = a11.*chi11k + a12.*chi12k + a13.*chi13k + ...
            a22.*chi22k + a23.*chi23k + a33.*chi33k ;

delta = real(ifftn(delta));

