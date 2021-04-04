function delta = chitensor2delta(chi, Params)
    % Compute field shift (delta) from chi and Params 
    % Input : 
    %           chi: 3x3 cell
    %           Params: 
    delta = zeros(size(chi{1,1}));
    aii = conv_kernel_tensor_arr(Params);
    chik = cell(3,3);
    for i=1:3
        for j=1:3
            chik{i,j} = fftn(chi{i,j});
            delta = delta + aii{i,j}.*chik{i,j};
        end
    end
    delta = real(ifftn(delta));
end