function [a11, a12, a13, a22, a23, a33] = conv_kernel_tensor(Params)
    Nx = Params.sizeVol(1);
    Ny = Params.sizeVol(2);
    Nz = Params.sizeVol(3);

    dkx = 1/Params.fov(1);
    dky = 1/Params.fov(2);
    dkz = 1/Params.fov(3);

    % convolution kernel 
    kx = linspace(-Nx/2+1, Nx/2, Nx).*dkx;
    ky = linspace(-Ny/2+1, Ny/2, Ny).*dky;
    kz = linspace(-Nz/2+1, Nz/2, Nz).*dkz;

    [KX_Grid, KY_Grid, KZ_Grid] = meshgrid(ky, kx, kz);  % mesh in k space
    % shift the kernel instead of the k-space data
    KX_Grid = ifftshift(KX_Grid);
    KY_Grid = ifftshift(KY_Grid);
    KZ_Grid = ifftshift(KZ_Grid);
    KSq = KX_Grid.^2 + KY_Grid.^2 + KZ_Grid.^2;          % k^2

    H0 = [0, 0, 1]';                                     % H0 in lab frame
    H0 = Params.TAng'*H0;                                % H0 in image/subject space
    h1 = H0(1);                                          % Params.TAng' is equivalent to Params.TAnginv
    h2 = H0(2);
    h3 = H0(3);

    KHdKSq = (KX_Grid*h1 + KY_Grid*h2 + KZ_Grid*h3)./KSq; % nan at the center

    a11 = h1.^2/3 - KHdKSq.*KX_Grid*h1;
    a22 = h2.^2/3 - KHdKSq.*KY_Grid*h2;
    a33 = h3.^2/3 - KHdKSq.*KZ_Grid*h3;

    a12 = 2*h1.*h2/3 - KHdKSq.*(KX_Grid*h2 + KY_Grid*h1);
    a13 = 2*h1.*h3/3 - KHdKSq.*(KX_Grid*h3 + KZ_Grid*h1);
    a23 = 2*h2.*h3/3 - KHdKSq.*(KY_Grid*h3 + KZ_Grid*h2);              

    % remove nan
    a11(isnan(a11)) = 0;  
    a12(isnan(a12)) = 0;   
    a13(isnan(a13)) = 0;
    a22(isnan(a22)) = 0;  
    a23(isnan(a23)) = 0;   
    a33(isnan(a33)) = 0;

end