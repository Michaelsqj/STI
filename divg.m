function div = divg(image)
    % Compute divergence of the input
    % Input: dim=4, 
    div = adjDx(image(:,:,:,1)) + adjDy(image(:,:,:,2)) + adjDz(image(:,:,:,3));
    
end

function res = adjDx(x)
    res = x([1,1:end-1],:,:) - x;
    res(1,:,:) = -x(1,:,:);
    res(end,:,:) = x(end-1,:,:);
end

function res = adjDz(x)
    res = x(:,:,[1,1:end-1]) - x;
    res(:,:,1) = -x(:,:,1);
    res(:,:,end) = x(:,:,end-1);
end

function res = adjDy(x)
    res = x(:,[1,1:end-1],:) - x;
    res(:,1,:) = -x(:,1,:);
    res(:,end,:) = x(:,end-1,:);
end
