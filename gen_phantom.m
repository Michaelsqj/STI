function gen_phantom()
    Params.B0 = 7;              % Tesla
    Params.gamma = 42.577e6;      % Hz/T
    Params.fov = [128, 128, 128];        % mm     
    Params.sizeVol = [128, 128, 128];    
    Params.voxSize = Params.fov./Params.sizeVol;    % mm
    Params.pathname = 'data\';
    Params.filename_chi = 'data\chi_phantom.mat';

    % chi_distribution of numerical phantoms
    % head phantom, with anisotropic susceptibility tensor, in subject frame
    %     A, a , b, c, x0, y0, z0, phi, theta, psi
    E11 = [0.05, 0.6, 0.85, 0.7, 0, 0, 0, 0, 0, 0; ...        % bg isotropic
         0.05, 0.15, 0.15, 0.3, 0, 0.5, 0, 0, 0, 0; ...       % isotropic
         -0.138, 0.1, 0.35, 0.4, -0.25, -0.15, 0, 0, 0, -15; ...   % x pointing tensor 
         -0.156, 0.15, 0.4, 0.4, 0.25, -0.15, 0, 0, 0, 12; ...     % y pointing tensor 
         -0.156, 0.1, 0.1, 0.5, 0, -0.5, 0, 0, 0, 0; ...     % z pointing tensor     
         ];
    chi11 = phantom3d(E11, 128);

    E22 = [0.05, 0.6, 0.85, 0.7, 0, 0, 0, 0, 0, 0; ...        % bg isotropic
         0.05, 0.15, 0.15, 0.3, 0, 0.5, 0, 0, 0, 0; ...       % isotropic
         -0.156, 0.1, 0.35, 0.4, -0.25, -0.15, 0, 0, 0, -15; ...   % x pointing tensor 
         -0.138, 0.15, 0.4, 0.4, 0.25, -0.15, 0, 0, 0, 12; ...     % y pointing tensor 
         -0.156, 0.1, 0.1, 0.5, 0, -0.5, 0, 0, 0, 0; ...     % z pointing tensor     
         ];
    chi22 = phantom3d(E22, 128);

    E33 = [0.05, 0.6, 0.85, 0.7, 0, 0, 0, 0, 0, 0; ...        % bg isotropic
         0.05, 0.15, 0.15, 0.3, 0, 0.5, 0, 0, 0, 0; ...       % isotropic
         -0.156, 0.1, 0.35, 0.4, -0.25, -0.15, 0, 0, 0, -15; ...   % x pointing tensor 
         -0.156, 0.15, 0.4, 0.4, 0.25, -0.15, 0, 0, 0, 12; ...     % y pointing tensor 
         -0.138, 0.1, 0.1, 0.5, 0, -0.5, 0, 0, 0, 0; ...     % z pointing tensor     
         ];
    chi33 = phantom3d(E33, 128);

    chi12 = zeros(size(chi11));                 % off diagonal tensor components
    chi13 = zeros(size(chi11));
    chi23 = zeros(size(chi11));    

    BrainMask = (abs(chi11) > 0);    

    % Compute MMS (chiavg), MSA (chiani), eigvs (eigen vectors), 
    %         PEV (principle eig vector), CMAP (?)

    
    % save phantom data
    save(Params.filename_chi, 'Params', 'chi11', 'chi12', 'chi13', 'chi22', 'chi23', 'chi33', 'BrainMask', ...
            'chiavg', 'chiani', 'eigvs', 'PEV', 'CMAP')
    disp([Params.filename_chi, ' saved.'])
end