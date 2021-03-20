function wG = gradient_mask_all(rawImage, Mask)
    % Generate weighting using gradient of rawImage
    % with percentage of voxels in Mask to be edges
    % edge is 0, non-edge is 1 for regularization, note wG.^2 = wG
    
    % ---------------- Compute Gradient ----------------
    rawImage = rawImage.*Mask;
    wG = abs(grad(rawImage)); % gradient map
    
    % ---------------- Compute Threshold ----------------
    % initial guess of the threshold
    wGtemp = sqrt(sum(wG.^2, 4));       % average over dimentions
    wG = wGtemp<0.05;
end