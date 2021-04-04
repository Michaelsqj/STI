function gen_phantom(filename)
    Params.B0 = 7;              % Tesla
    Params.gamma = 42.577e6;      % Hz/T
    Params.fov = [128, 128, 128];        % mm     
    Params.sizeVol = [128, 128, 128];    
    Params.voxSize = Params.fov./Params.sizeVol;    % mm
    Params.pathname = 'data/';
    Params.filename_chi = ['data/', filename];
    Params.chemshift = 0.2;     %chemical shift = 0.7ppm

    chi = cell(3,3);
    % chi_distribution of numerical phantoms
    % head phantom, with anisotropic susceptibility tensor, in subject frame
    %      A,       a ,     b,      c,      x0,     y0,     z0, phi, theta, psi
    E11 = [0.05,    0.6,    0.85,   0.7,    0,      0,      0, 0, 0, 0; ...         % bg isotropic
           0.05,    0.15,   0.15,   0.3,    0,      0.5,    0, 0, 0, 0; ...         % isotropic
           -0.138,  0.1,    0.35,   0.4,    -0.25,  -0.15,  0, 0, 0, -15; ...       % x pointing tensor 
           -0.156,  0.15,   0.4,    0.4,    0.25,   -0.15,  0, 0, 0, 12; ...        % y pointing tensor 
           -0.156,  0.1,    0.1,    0.5,    0,      -0.5,   0, 0, 0, 0; ...         % z pointing tensor     
         ];
    chi{1,1} = phantom3d(E11, 128);

    E22 = [0.05,    0.6,    0.85,   0.7,    0,      0,      0, 0, 0, 0; ...         % bg isotropic
           0.05,    0.15,   0.15,   0.3,    0,      0.5,    0, 0, 0, 0; ...         % isotropic
           -0.156,  0.1,    0.35,   0.4,    -0.25,  -0.15,  0, 0, 0, -15; ...       % x pointing tensor 
           -0.138,  0.15,   0.4,    0.4,    0.25,   -0.15,  0, 0, 0, 12; ...        % y pointing tensor 
           -0.156,  0.1,    0.1,    0.5,    0,      -0.5,   0, 0, 0, 0; ...         % z pointing tensor     
         ];
    chi{2,2} = phantom3d(E22, 128);

    E33 = [0.05,    0.6,    0.85,   0.7,    0,      0,      0, 0, 0, 0; ...         % bg isotropic
           0.05,    0.15,   0.15,   0.3,    0,      0.5,    0, 0, 0, 0; ...         % isotropic
           -0.156,  0.1,    0.35,   0.4,    -0.25,  -0.15,  0, 0, 0, -15; ...       % x pointing tensor 
           -0.156,  0.15,   0.4,    0.4,    0.25,   -0.15,  0, 0, 0, 12; ...        % y pointing tensor 
           -0.138,  0.1,    0.1,    0.5,    0,      -0.5,   0, 0, 0, 0; ...         % z pointing tensor     
         ];
    chi{3,3} = phantom3d(E33, 128);

    for i = 1:3
        for j = 1:3
            if i~=j
                chi{i,j}=zeros(size(chi{1,1}));
            end
        end
    end
    
    BrainMask = (abs(chi{1,1}) > 0);    
    
    % generate chemical shift also using phantom3d
    e = [Params.chemshift,  0.1,    0.35,   0.4,    -0.25,  -0.15,  0, 0, 0, -15];
    chemshift_map = phantom3d(e, 128);

    % Compute MMS (chiavg), MSA (chiani), eigvs (eigen vectors), 
    %         PEV (principle eig vector), CMAP (?)
    [chi1,chi2, chi3, ~, ~, ~] = tensor2eig(chi);
    
    chiavg = (chi1+chi2+chi3)/3;
    chiani = chi1 - (chi2+chi3)/2;
    % save phantom data
    save(Params.filename_chi, 'Params', 'chi', 'BrainMask', 'chiavg', 'chiani', 'chemshift_map')
    disp([Params.filename_chi, ' saved.'])
end
