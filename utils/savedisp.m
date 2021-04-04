function flag = savedisp(pathname, filename, var, clim, slice2disp, unit)
% savedisp(plotpath, 'chi11', chi11, [-0.11, 0.11]);

if nargin < 5
    slice2disp = [];
    unit = 'ppm';
elseif nargin < 6
    unit = 'ppm';
end

flag = 0;

N = ndims(var);

if N > 2 && isempty(slice2disp)
    mimage(var); 
    if N < 4
        colormap(gray); caxis(clim); colorbar; 
    end
    title([inputname(3), ' (', unit, ')'])
else
    % slice of 3D volume
    figure;
    if N > 2
        imagesc(var(:,:,slice2disp));
    else
        imagesc(var);
    end
    axis equal; axis tight;
    caxis(clim); colormap(gray); colorbar; title([inputname(3), '\_s', num2str(slice2disp), ' (', unit, ')']);
    axis off;
    set(gca, 'fontsize', 20)
end

print(gcf, fullfile(pathname, filename), '-dpng'); close(gcf)

end

