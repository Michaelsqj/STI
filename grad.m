function g = grad(image)
    % Input volume: dim=3
    % Output the gradient of the input volume [Nx, Ny, Nz, 3]
    
    Dx = image([2:end,end],:,:) - image;
    Dy = image(:,[2:end,end],:) - image;
    Dz = image(:,:,[2:end,end]) - image;
    g = cat(4, Dx, Dy, Dz);
end