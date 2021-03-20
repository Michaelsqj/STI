% Params.phantomID = 2
% 
%         A, a , b, c, x0, y0, z0, phi, theta, psi
% E11 = [0.05, 0.6, 0.85, 0.7, 0, 0, 0, 0, 0, 0; ...        % bg isotropic
%      -0.138, 0.15, 0.15, 0.3, 0, 0.5, 0, 0, 0, 0; ...       % x pointing tensor 
% 
%      -0.156, 0.05,0.3, 0.3, -0.125, 0.05, 0, 0, 0, 15; ...   % z pointing tensor
%      -0.138, 0.05,0.3, 0.3, -0.125, 0.05, 0, 0, 0, 15; ...   % x pointing tensor
%      -0.097-0.05, 0.05,0.3, 0.3, -0.125, 0.05, 0, 0, 0, 15; ... % 45 degree  
%      -0.1015-0.05, 0.05,0.3, 0.3, -0.125, 0.05, 0, 0, 0, 15; ... % 30 degree
%       -0.0925-0.05, 0.05,0.3, 0.3, -0.125, 0.05, 0, 0, 0, 15; ... % 60 degree
% 
%       0.05, 0.1, 0.25, 0.4, -0.25, 0.05, 0, 0, 0, 15; ...    % neighbouring isotropic        
% 
%       -0.156, 0.05, 0.3, 0.3, 0.125, 0.05, 0, 0, 0, -12; ...     % z pointing tensor
%       -0.138, 0.05, 0.3, 0.3, 0.125, 0.05, 0, 0, 0, -12; ...     % x pointing tensor 
%       -0.097-0.05, 0.05, 0.3, 0.3, 0.125, 0.05, 0, 0, 0, -12; ...  % 45 degree
%       -0.1015-0.05, 0.05, 0.3, 0.3, 0.125, 0.05, 0, 0, 0, -12; ...  % 30 degree          
%       -0.0925-0.05, 0.05, 0.3, 0.3, 0.125, 0.05, 0, 0, 0, -12; ...  % 60 degree       
% 
%       0.05, 0.1, 0.25, 0.4, 0.25, 0.05, 0, 0, 0, -12; ...   % neighbouring isotropic 
% 
%      -0.156, 0.1, 0.1, 0.5, 0, -0.5, 0, 0, 0, 0; ...     % y pointing tensor            
%      ];
% 
% chi11 = phantom3d(E11, 128);
% 
% E22 = [0.05, 0.6, 0.85, 0.7, 0, 0, 0, 0, 0, 0; ...        % bg isotropic
%      -0.156, 0.15, 0.15, 0.3, 0, 0.5, 0, 0, 0, 0; ...       % x pointing tensor 
% 
%      -0.156, 0.05,0.3, 0.3, -0.125, 0.05, 0, 0, 0, 15; ...   % z pointing tensor 
%       0.05, 0.1, 0.25, 0.4, -0.25, 0.05, 0, 0, 0, 15; ...    % neighbouring isotropic        
% 
%      -0.156, 0.05, 0.3, 0.3, 0.125, 0.05, 0, 0, 0, -12; ...     % z pointing tensor 
%       0.05, 0.1, 0.25, 0.4, 0.25, 0.05, 0, 0, 0, -12; ...   % neighbouring isotropic
% 
%      -0.138, 0.1, 0.1, 0.5, 0, -0.5, 0, 0, 0, 0; ...     % y pointing tensor            
%      ];
% 
% chi22 = phantom3d(E22, 128);
% 
% E33 = [0.05, 0.6, 0.85, 0.7, 0, 0, 0, 0, 0, 0; ...        % bg isotropic
%      -0.156, 0.15, 0.15, 0.3, 0, 0.5, 0, 0, 0, 0; ...       % x pointing tensor 
% 
%      -0.138, 0.05,0.3, 0.3, -0.125, 0.05, 0, 0, 0, 15; ...   % z pointing tensor 
%      -0.156, 0.05,0.3, 0.3, -0.125, 0.05, 0, 0, 0, 15; ...   % x pointing tensor
%      -0.097-0.05, 0.05,0.3, 0.3, -0.125, 0.05, 0, 0, 0, 15; ... % 45 degree  
%      -0.0925-0.05, 0.05,0.3, 0.3, -0.125, 0.05, 0, 0, 0, 15; ... % 30 degree
%      -0.1015-0.05, 0.05,0.3, 0.3, -0.125, 0.05, 0, 0, 0, 15; ... % 60 degree
% 
%       0.05, 0.1, 0.25, 0.4, -0.25, 0.05, 0, 0, 0, 15; ...    % neighbouring isotropic        
% 
%      -0.138, 0.05, 0.3, 0.3, 0.125, 0.05, 0, 0, 0, -12; ...     % z pointing tensor
%      -0.156, 0.05, 0.3, 0.3, 0.125, 0.05, 0, 0, 0, -12; ...     % z pointing tensor
%      -0.097-0.05, 0.05, 0.3, 0.3, 0.125, 0.05, 0, 0, 0, -12; ...  % 45 degree
%      -0.0925-0.05, 0.05, 0.3, 0.3, 0.125, 0.05, 0, 0, 0, -12; ...  % 30 degree
%      -0.1015-0.05, 0.05, 0.3, 0.3, 0.125, 0.05, 0, 0, 0, -12; ...  % 60 degree
% 
%       0.05, 0.1, 0.25, 0.4, 0.25, 0.05, 0, 0, 0, -12; ...   % neighbouring isotropic        
% 
%      -0.156, 0.1, 0.1, 0.5, 0, -0.5, 0, 0, 0, 0; ...     % y pointing tensor            
%      ];
% 
% chi33 = phantom3d(E33, 128);
% 
% chi12 = zeros(size(chi11));                 % off diagonal tensor components
% chi13 = zeros(size(chi11));
% chi23 = zeros(size(chi11));
% 
% E12 =[
%     0.009, 0.05,0.3, 0.3, -0.125, 0.05, 0, 0, 0, 15; ... % 45 degree 
%     0.009, 0.05, 0.3, 0.3, 0.125, 0.05, 0, 0, 0, -12; ...  % 45 degree        
%     0.0078, 0.05,0.3, 0.3, -0.125, 0.05, 0, 0, 0, 15; ... % 30 degree 
%      0.0078, 0.05, 0.3, 0.3, 0.125, 0.05, 0, 0, 0, -12; ...  % 30 degree                 
%     ];
% chi13 = phantom3d(E12, 128);

