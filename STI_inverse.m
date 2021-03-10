function [chi11, chi12, chi13, chi22, chi23, chi33, flag, relres, iter, resvec] = STI_inverse(STIParams, maxit, tol, alpha, beta)
    % Input:  
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
    %
    % Output:
    %     chi{i}{j}, 3x3 9 susceptibility tensor components
    %     flag: convergence flag: 0 means converged (similar as lsqr)
    %     relres: final relative residual
    %     iter: total iteration number
    %     resvec: relative residual history

    TV = TVOP;

    % load MO dataset
    OriNum = STIParams.OriNum;
    sizeVol = STIParams.sizeVol;

    deltaBArray = zeros(sizeVol(1), sizeVol(2), sizeVol(3), OriNum, 'single');        % in real space
    ParamsArray = cell(OriNum, 1);

    % filenameFreqMap contains: Params, deltaB, BrainMask

    for OriInd = 1:OriNum
        disp(['loading data ', num2str(OriInd), ' ...'])    
        S = load(STIParams.filenameFreqMap{OriInd});
        Params = S.Params;
        if OriInd == 1
            BrainMask = S.BrainMask;
        end
        deltaBArray(:,:,:,OriInd) = S.deltaB.*BrainMask;                    % deltaB in real space, masked
        Params.Weight = BrainMask;
        ParamsArray{OriInd} = Params;
    end

    clear S
    disp('Done.')

    Params = ParamsArray{1};

    % --------------- loading the isotropic chi from QSM fitting
    disp('loading QSM data...')
    S = load(STIParams.filenameQSM);
    chiavg_qsm = S.chiavg;      % simulation
    clear S    
    
    wG = gradient_mask_all(chiavg_qsm, BrainMask, STIParams.EdgePer, 0, 0.1);  % wG is 4D

    % using least square method for solving STI
    % 
    % setting up b for solving Ax = b;
    VoxNum = prod(sizeVol);
    b = zeros((OriNum+12)*VoxNum, 1, 'single');    % 6 for normal STI, 9 for L2, 12 for L2+maskChiani

    for OriInd = 1:OriNum    
        temp = deltaBArray(:,:,:,OriInd).*ParamsArray{OriInd}.Weight;    
        b(((OriInd - 1)*VoxNum+1) : OriInd*VoxNum) = temp(:);
    end
    b((OriNum*VoxNum+1):end) = 0;           % For regularization

    clear temp deltaBArray
    aii_array = conv_kernel_tensor_arr(Params, OriNum);

    disp('solving Susceptibility Tensor Imaging (STI) with MMSR-STI ...')

    tic
    [chi_tensor, flag, relres, iter, resvec] = lsqr(@afun,b,tol,maxit);
    toc

    % ------------------------------------------------------------------
    % change solution format
    chi11 = BrainMask.*reshape(chi_tensor(1:VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
    chi12 = BrainMask.*reshape(chi_tensor((1*VoxNum+1):2*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
    chi13 = BrainMask.*reshape(chi_tensor((2*VoxNum+1):3*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
    chi22 = BrainMask.*reshape(chi_tensor((3*VoxNum+1):4*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
    chi23 = BrainMask.*reshape(chi_tensor((4*VoxNum+1):5*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
    chi33 = BrainMask.*reshape(chi_tensor((5*VoxNum+1):6*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));

    resvec = resvec./norm(b);                       % convert to relative residual vector

    % internal function including afun
    function y = afun(x,transp_flag)
            
        if strcmp(transp_flag,'transp')              % y = A'*x
            y = zeros(6*VoxNum, 1, 'single');          % chitensor (6Nv X 1)
            x = single(x);                             % change to single format ((Nori*Nv x 1))      
            
            for orient_i = 1:OriNum
                    Params = ParamsArray{orient_i};
                    delta = x(((orient_i - 1)*VoxNum+1) : orient_i*VoxNum);            
                    delta = reshape(delta, sizeVol(1), sizeVol(2), sizeVol(3));
                    delta = fftn(delta.*Params.Weight);
                    
                    % from the 6N X N matrix
                    chi11_temp = real(ifftn(a11_array(:,:,:,orient_i).*delta));
                    chi12_temp = real(ifftn(a12_array(:,:,:,orient_i).*delta));
                    chi13_temp = real(ifftn(a13_array(:,:,:,orient_i).*delta));
                    chi22_temp = real(ifftn(a22_array(:,:,:,orient_i).*delta));
                    chi23_temp = real(ifftn(a23_array(:,:,:,orient_i).*delta));
                    chi33_temp = real(ifftn(a33_array(:,:,:,orient_i).*delta));

                    y = y + [chi11_temp(:); chi12_temp(:); chi13_temp(:); chi22_temp(:); chi23_temp(:); chi33_temp(:)];                                    
            end
            
            % regularize x11, x22, x33 on ~BrainMask;
            % regularize x12, x13, x23 on ~STIParams.maskMSA(:)
            y = y + alpha*x((OriNum*VoxNum+1):(OriNum+6)*VoxNum).*cat(1, ~BrainMask(:), ~STIParams.maskMSA(:), ...
                        ~STIParams.maskMSA(:), ~BrainMask(:), ~STIParams.maskMSA(:), ~BrainMask(:));
        
            temp = alpha*x(((OriNum+6)*VoxNum+1):(OriNum+9)*VoxNum).*repmat(~STIParams.maskMSA(:), [3,1]); 
            y(1:VoxNum)              = y(1:VoxNum) + temp(1:VoxNum) - temp(2*VoxNum+1:3*VoxNum);                         % x11
            y((3*VoxNum+1):4*VoxNum) = y((3*VoxNum+1):4*VoxNum) - temp(1:VoxNum) + temp(VoxNum+1:2*VoxNum);              % -x22              
            y((5*VoxNum+1):6*VoxNum) = y((5*VoxNum+1):6*VoxNum) - temp(VoxNum+1:2*VoxNum) + temp(2*VoxNum+1:3*VoxNum);   % -x33       
            
            temp = beta*(TV'*(reshape(wG(:).*x(((OriNum+9)*VoxNum+1):(OriNum+12)*VoxNum), size(wG))));
            y(1:VoxNum) = y(1:VoxNum) + temp(:);                                 % x11
            y((3*VoxNum+1):4*VoxNum) = y((3*VoxNum+1):4*VoxNum) + temp(:);       % x22
            y((5*VoxNum+1):6*VoxNum) = y((5*VoxNum+1):6*VoxNum) + temp(:);       % x33

            disp('Iteration of transpose A ...')

        elseif strcmp(transp_flag,'notransp')            % y = A*x
                y = zeros((OriNum+12)*VoxNum, 1, 'single');        % x is chi_tensor, y is N deltaB and 12 regularization
                x = single(x);                                     % change to single format
            
                x11 = reshape(x(1:VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
                x22 = reshape(x((3*VoxNum+1):4*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
                x33 = reshape(x((5*VoxNum+1):6*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
                        
                chi11k = fftn(x11);
                chi12k = fftn(reshape(x((1*VoxNum+1):2*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3)));
                chi13k = fftn(reshape(x((2*VoxNum+1):3*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3)));
                chi22k = fftn(x22);
                chi23k = fftn(reshape(x((4*VoxNum+1):5*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3)));
                chi33k = fftn(x33);
            
                % data fidelity terms
            for orient_i = 1:OriNum
                    Params = ParamsArray{orient_i};    
                    
                    delta = a11_array(:,:,:,orient_i).*chi11k + a12_array(:,:,:,orient_i).*chi12k + a13_array(:,:,:,orient_i).*chi13k + ...
                                a22_array(:,:,:,orient_i).*chi22k + a23_array(:,:,:,orient_i).*chi23k + a33_array(:,:,:,orient_i).*chi33k ;

                    delta = Params.Weight.*real(ifftn(delta));
                
                    y(((orient_i - 1)*VoxNum+1) : orient_i*VoxNum) = delta(:);          
            end

            % regularize x11, x22, x33 on ~BrainMask; 
            % regularize x12, x13, x23 on ~STIParams.maskMSA(:)
            y(OriNum*VoxNum+1:(OriNum+6)*VoxNum) = alpha*x(:).*cat(1, ~BrainMask(:), ~STIParams.maskMSA(:), ...
                        ~STIParams.maskMSA(:), ~BrainMask(:), ~STIParams.maskMSA(:), ~BrainMask(:));

            y(((OriNum+6)*VoxNum+1):((OriNum+9)*VoxNum)) = alpha*repmat(~STIParams.maskMSA(:), [3,1]).*cat(1, ...
                                x11(:) - x22(:), x22(:) - x33(:), x33(:) - x11(:));             
                    
            % L2 regularization on edge prior
            TVreg = beta*wG.*(TV*(x11 + x22 + x33));
            y(((OriNum+9)*VoxNum+1):((OriNum+12)*VoxNum)) = TVreg(:);

            disp('Iteration of nontranspose A ...')
                
        end
    end
end
