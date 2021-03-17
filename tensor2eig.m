function [chi1, chi2, chi3, eig0, eig1, eig2] = tensor2eig(chi)
    % calculate the the eigenvalue of the full susceptibility tensor
    % OUTPUT:
    % chi1, chi2, chi3: three eigen values of the full susceptibility tensor
    % 
    % INPUT:
    % chi_ij:   full susceptibility tensor components (symmetric tensor)
    % eigvs:    structure defining the eigen vectors 
    %           should have eigvs.eig0, eigvs.eig1, eigvs.eig2, which are
    %           filenames referenced to eigenvector data files exported from DTI


    [Nx, Ny, Nz] = size(chi{1,1});

    chi1 = zeros(size(chi{1,1}));
    chi2 = zeros(size(chi{1,1}));
    chi3 = zeros(size(chi{1,1}));

    eig0 = zeros(3, Nx, Ny, Nz);
    eig1 = eig0;
    eig2 = eig0;

    % loop through each voxel
    h = waitbar(0, 'doing diagonalization...');
    for iz = 1:Nz
        for iy = 1:Ny
            for ix = 1:Nx      
                chiT = zeros(3,3);
                for i = 1:3
                    for j = 1:3
                        chiT(i,j)=chi{i,j}(ix,iy,iz);
                    end
                end
                [V, D] = eig(chiT);
                chidiag = diag(D);
                [B, Ind] = sort(chidiag);

                chi1(iy, ix, iz) = B(end);
                chi2(iy, ix, iz) = B(end-1);
                chi3(iy, ix, iz) = B(end-2);

                eig0(:, iy, ix, iz) = V(:, Ind(end));   % priciple eigen vector
                eig1(:, iy, ix, iz) = V(:, Ind(end-1));
                eig2(:, iy, ix, iz) = V(:, Ind(end-2));

            end
        end
        waitbar(iz/Nz);
    end
    close(h);
end
