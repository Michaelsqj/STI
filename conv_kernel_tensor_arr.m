function aii_array = conv_kernel_tensor_arr(P, OriNum)
    % Return: aii_array: 3x3 cell, 
    %           each element is a [Nx, Ny, Nz, OriNum] array 
    %           if OriNum=1, each element is [Nx, Ny, Nz] array
    if isa(P, 'struct')
        Params = P;
        ParamsArray = {Params};
        OriNum = 1;
    elseif isa(P, 'cell')
        ParamsArray = P;
        Params = ParamsArray{1};
    else
        disp('conv_kernel_tensor_arr: wrong input type');
    end
    
    Nx = Params.sizeVol(1);
    Ny = Params.sizeVol(2);
    Nz = Params.sizeVol(3);

    dkx = 1/Params.fov(1);
    dky = 1/Params.fov(2);
    dkz = 1/Params.fov(3);

    kx = linspace(-Nx/2+1, Nx/2, Nx).*dkx;
    ky = linspace(-Ny/2+1, Ny/2, Ny).*dky;
    kz = linspace(-Nz/2+1, Nz/2, Nz).*dkz;

    [KX_Grid, KY_Grid, KZ_Grid] = meshgrid(kx, ky, kz);  
    % shift the kernel instead of the k-space data
    KX_Grid = ifftshift(KX_Grid);
    KY_Grid = ifftshift(KY_Grid);
    KZ_Grid = ifftshift(KZ_Grid);
    KSq = KX_Grid.^2 + KY_Grid.^2 + KZ_Grid.^2;          % k^2
    K_Grid = {KX_Grid, KY_Grid, KZ_Grid};
    clear kx ky kz dkx dky dkz Nx Ny Nz

    H0 = [0, 0, 1]';                                                    % H0 in lab frame
    H0subArray = zeros(3, OriNum);
    for OriInd = 1:OriNum
        Params = ParamsArray{OriInd};   
        H0subArray(:,OriInd) = Params.TAng'*H0;                          % H0 in image/subject space     
    end

    % precalculated coefficient matrix
    aii_array = cell(3, 3);
    for i = 1:3
        for j = 1:3
            aii_array{i, j} = zeros([Params.sizeVol, OriNum]);
        end
    end

    for OriInd = 1:OriNum
        h = [H0subArray(1,OriInd), H0subArray(2,OriInd), H0subArray(3,OriInd)];
        KHdKSq = (K_Grid{1}*h(1) + K_Grid{2}*h(2) + K_Grid{3}*h(3))./KSq; % nan at the center
        for i = 1:3
            for j = 1:3
                aii_array{i, j}(:,:,:,OriInd) = h(i).*h(j)/3 - KHdKSq.*K_Grid{i}*h(j);
            end
        end
    end

    % remove nan
    for i = 1:3
        for j = 1:3
            aii_array{i, j}(isnan(aii_array{i, j})) = 0;
            aii_array{i, j} = squeeze(aii_array{i, j});
        end
    end
    
end