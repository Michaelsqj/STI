function [chi,flag,relres,iter,resvec,lsvec] = STI_inverse(STIParams_filename, maxit, tol, alpha)
    % Input:  
    %       STIParams includes fields:
    %           OriNum: Orientation Number
    %           fov, sizeVol, VoxSize: field of view, size of Volume, voxel size
    %           filename_deltaB: cell array containing the filenames of deltaB maps
    %           maskMSA: White Matter Mask
    %           QSM: 
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

    % --------------- Load Data ---------------
    STIParams = load(['data/', STIParams_filename]).STIParams;
    OriNum = STIParams.OriNum;
    sizeVol = STIParams.sizeVol;
    deltaBArray = zeros(sizeVol(1), sizeVol(2), sizeVol(3), OriNum, 'single');        % in real space
    ParamsArray = cell(OriNum, 1);
    for OriInd = 1:OriNum
        disp(['loading data ', num2str(OriInd), ' ...'])    
        temp = load(STIParams.filename_deltaB{OriInd});
        Params = temp.Params;
        if OriInd == 1
            BrainMask = temp.BrainMask;
        end
        deltaBArray(:,:,:,OriInd) = temp.deltaB.*BrainMask;                    % deltaB in real space, masked
        Params.BrainMask = BrainMask;
        ParamsArray{OriInd} = Params;
    end
    disp('Finished Loading Data.')    
    
    % --------------- using least square method for solving STI
    % setting up b for solving Ax = b;
    VoxNum = prod(sizeVol);
    b = zeros((OriNum+12)*VoxNum, 1, 'single');
    for OriInd = 1:OriNum    
        temp = deltaBArray(:,:,:,OriInd).*BrainMask;    
        b(((OriInd - 1)*VoxNum+1) : OriInd*VoxNum) = temp(:);
    end
    aii_array = conv_kernel_tensor_arr(ParamsArray, OriNum);
    
    tic
    [chi_tensor,flag,relres,iter,resvec,lsvec] = lsqr(@afun,b,tol,maxit);
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
            mask_arr = cat(1, ~BrainMask(:), ~STIParams.maskMSA(:), ~STIParams.maskMSA(:),...
                              ~STIParams.maskMSA(:), ~BrainMask(:), ~STIParams.maskMSA(:),...
                              ~STIParams.maskMSA(:), ~STIParams.maskMSA(:), ~BrainMask(:));
            y = y + alpha*delta(:).*mask_arr(:);

            % Regularize of chi11=chi22=chi33
            delta = reshape(x((OriNum+9)*VoxNum+1: (OriNum+12)*VoxNum), VoxNum, 3);
            y(1:VoxNum) = y(1:VoxNum) + (~STIParams.maskMSA(:)).* (delta(:,1)+delta(:,3));
            y(4*VoxNum+1:5*VoxNum) = y(4*VoxNum+1:5*VoxNum) + (~STIParams.maskMSA(:)).* (-delta(:,1)+delta(:,2));
            y(8*VoxNum+1:9*VoxNum) = y(8*VoxNum+1:9*VoxNum) + (~STIParams.maskMSA(:)).* (-delta(:,2)-delta(:,3));
            
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
                        delta = zeros(sizeVol, 'double');
                        for i = 1:3
                            for j = 1:3
                                delta = delta + aii_array{i,j}(:,:,:,orient_i).*chi_k{i,j};
                            end
                        end
                        delta = Params.BrainMask.*real(ifftn(delta));
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
                
                disp('Iteration of nontranspose A ...')
        end
    end
end
