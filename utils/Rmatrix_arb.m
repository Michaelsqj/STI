function R=Rmatrix_arb(rl, ap, fh, orderR)
    % convert the rotation angle into rotation matrix
    if nargin < 4
        orderR = [1, 2, 3];
    end

    T(:,:,1) = [ 1, 0, 0; ...
            0, cosd(rl), -sind(rl); ...
            0, sind(rl), cosd(rl)];       % rotation about x, counter clockwise

    T(:,:,2) = [ cosd(ap), 0, sind(ap); ...
             0, 1, 0; ...
            -sind(ap), 0, cosd(ap)];      % rotation about y, counter clockwise

    T(:,:,3) = [cosd(fh), -sind(fh), 0;...
            sind(fh), cosd(fh), 0;...
            0, 0, 1];                   % rotation about z, counter clockwise


    R = T(:,:,orderR(1))*T(:,:,orderR(2))*T(:,:,orderR(3));

end