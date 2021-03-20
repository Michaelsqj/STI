function showim(image, z)
    if nargin<2
        shape = size(image);
        z = round(shape(3)/2);
    end
    mini = double(min(image, [], 'all'));
    maxi = double(max(image, [], 'all'));
    imshow(squeeze(image(:,:,z)),[mini,maxi]);
    colorbar();
end