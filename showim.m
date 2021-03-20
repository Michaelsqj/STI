function showim(image, z)
    if nargin<2
        shape = size(image);
        z = round(shape(3)/2);
    end
    mini = min(min(min(image)));
    maxi = max(max(max(image)));
    imshow(squeeze(image(:,:,z)),[mini,maxi]);
    colorbar();
end