function [chi, flag, relres, iter, resvec] = STI_inverse(STIParams, maxit, tol, alpha, beta)
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

    % --------------- using least square method for solving STI
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
    
    tic
    [chi_tensor, flag, relres, iter, resvec] = lsqr(@afun,b,tol,maxit);
    toc

    chi = cell(3,3);
    for ii = 1:3
        for jj = 1:3
            chi{ii,jj} = BrainMask.*reshape(chi_tensor(((ii-1)*3+jj-1)*VoxNum+1:((ii-1)*3+jj)*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
        end
    end

    resvec = resvec./norm(b);                       % convert to relative residual vector

    % internal function including afun
    function y = afun(x,transp_flag)
        % A * chi = delta
        % chi = A' * delta
        % size(A) = OriNum+9+3+3, 9*VoxNum
        if strcmp(transp_flag,'transp')              % y = A'*x
            y = zeros(9*VoxNum, 1, 'single');          % chitensor (9VoxNum X 1)
            x = single(x);                             % change to single format ((OriNum+9+3+3)*VoxNum)      
            
            % N orientation
            delta = reshape(x(1:OriNum*VoxNum), VoxNum, OriNum);
            for i = 1:3
                for j = 1:3
                    aii = reshape(aii_array{i,j}, VoxNum, OriNum);
                    y(((i-1)*3+j-1)*VoxNum+1:((i-1)*3+j)*VoxNum) = y(((i-1)*3+j-1)*VoxNum+1:((i-1)*3+j)*VoxNum) + ...
                                                                   sum(aii.*delta,2);
                end
            end
            
            % Regularize on BrainMask and maskMSA
            delta = reshape(x(OriNum*VoxNum+1: (OriNum+9)*VoxNum), VoxNum, 9);
            mask_arr = cat(2, ~BrainMask(:), ~STIParams.maskMSA(:), ~STIParams.maskMSA(:),...
                              ~STIParams.maskMSA(:), ~BrainMask(:), ~STIParams.maskMSA(:),...
                              ~STIParams.maskMSA(:), ~STIParams.maskMSA(:), ~BrainMask(:));
            y = y + alpha*delta(:).*mask_arr(:);

            % Regularize of chi11=chi22=chi33
            delta = reshape(x((OriNum+9)*VoxNum+1: (OriNum+12)*VoxNum), VoxNum, 3);
            y(1:VoxNum) = y(1:VoxNum) + (~STIParams.maskMSA(:)).* (delta(:,1)+delta(:,3));
            y(4*VoxNum+1:5*VoxNum) = y(4*VoxNum+1:5*VoxNum) + (~STIParams.maskMSA(:)).* (-delta(:,1)+delta(:,2));
            y(8*VoxNum+1:9*VoxNum) = y(8*VoxNum+1:9*VoxNum) + (~STIParams.maskMSA(:)).* (-delta(:,2)-delta(:,3));
            
            % Wg, divergence is transpose of gradient
            delta = reshape(x((OriNum+12)*VoxNum+1: (OriNum+15)*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3), 3);
            temp = wG.*divg(delta);
            y(1:VoxNum) = y(1:VoxNum) + temp(:);
            y(4*VoxNum+1:5*VoxNum) = y(4*VoxNum+1:5*VoxNum) + temp(:);
            y(8*VoxNum+1:9*VoxNum) = y(8*VoxNum+1:9*VoxNum) + temp(:);
            
            disp('Iteration of transpose A ...')

        elseif strcmp(transp_flag,'notransp')            % y = A*x
                y = zeros((OriNum+12)*VoxNum, 1, 'single');        % x is chi_tensor, y is N deltaB and 12 regularization
                x = single(x);                                     % change to single format
            
                x11 = reshape(x(1:VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
                x22 = reshape(x((3*VoxNum+1):4*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
                x33 = reshape(x((5*VoxNum+1):6*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3));
                
                chi_k = cell(3, 3);
                for i = 1:3
                    for j = 1:3
                        chi_k{i, j} = fftn(reshape(x((3*(i-1)+j-1)*VoxNum+1: (3*(i-1)+j)*VoxNum), sizeVol(1), sizeVol(2), sizeVol(3)));
                    end
                end
            
                % OriNum A*chi=delta
                for orient_i = 1:OriNum
                        Params = ParamsArray{orient_i};    
                        delta = zeros(VoxNum, 'single');
                        for i = 1:3
                            for j = 1:3
                                delta = delta + aii_array{i,j}(:,:,:,orient_i).*chi_k{i,j};
                            end
                        end
                        delta = Params.Weight.*real(ifftn(delta));
                        y(((orient_i - 1)*VoxNum+1) : orient_i*VoxNum) = delta(:);          
                end

                % regularize x11, x22, x33 on ~BrainMask; 
                % BrainMask is the region of brain on the NxNxN grid, 
                % out of the brain, chi should 0
                % regularize x12, x13, x21, x23, x31, x32 on ~STIParams.maskMSA(:)
                y(OriNum*VoxNum+1:(OriNum+9)*VoxNum) = alpha*x(:).*cat(1, ~BrainMask(:), ~STIParams.maskMSA(:), ~STIParams.maskMSA(:),...
                                                                          ~STIParams.maskMSA(:), ~BrainMask(:), ~STIParams.maskMSA(:),...
                                                                          ~STIParams.maskMSA(:), ~STIParams.maskMSA(:), ~BrainMask(:));

                y(((OriNum+9)*VoxNum+1):((OriNum+12)*VoxNum)) = alpha*repmat(~STIParams.maskMSA(:), [3,1]).*cat(1, ...
                                    x11(:) - x22(:), x22(:) - x33(:), x33(:) - x11(:));             

                % L2 regularization on edge prior
                temp = beta*wG.*(grad(x11 + x22 + x33));
                y(((OriNum+12)*VoxNum+1):((OriNum+15)*VoxNum)) = temp(:);
                
                disp('Iteration of nontranspose A ...')
        end
    end
end
