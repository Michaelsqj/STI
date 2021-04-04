function p = phantom3d(ellipse, n)
   %
   %   N is a scalar that specifies the grid size of P.
   %   If you omit the argument, N defaults to 64.
   % 
   %   P = PHANTOM3D(E,N) generates a user-defined phantom, where each row
   %   of the matrix E specifies an ellipsoid in the image.  E has ten columns,
   %   with each column containing a different parameter for the ellipsoids:
   %   
   %     Column 1:  A      the additive intensity value of the ellipsoid
   %     Column 2:  a      the length of the x semi-axis of the ellipsoid 
   %     Column 3:  b      the length of the y semi-axis of the ellipsoid
   %     Column 4:  c      the length of the z semi-axis of the ellipsoid
   %     Column 5:  x0     the x-coordinate of the center of the ellipsoid
   %     Column 6:  y0     the y-coordinate of the center of the ellipsoid
   %     Column 7:  z0     the z-coordinate of the center of the ellipsoid
   %     Column 8:  phi    phi Euler angle (in degrees) (rotation about z-axis)
   %     Column 9:  theta  theta Euler angle (in degrees) (rotation about x-axis)
   %     Column 10: psi    psi Euler angle (in degrees) (rotation about z-axis)
   %
   %   For purposes of generating the phantom, the domains for the x-, y-, and 
   %   z-axes span [-1,1].  Columns 2 through 7 must be specified in terms
   %   of this range.
   %   Example
   %   -------
   %        ph = phantom3d(128);
   %        figure, imshow(squeeze(ph(64,:,:)))

   p = zeros([n n n]);

   rng =  ( (0:n-1)-(n-1)/2 ) / ((n-1)/2); 

   [x,y,z] = meshgrid(rng,rng,rng);

   coord = [flatten(x); flatten(y); flatten(z)];

   p = flatten(p);

   for k = 1:size(ellipse,1)    
      A     = ellipse(k,1);            % Amplitude change for this ellipsoid
      asq   = ellipse(k,2)^2;        % a^2
      bsq   = ellipse(k,3)^2;        % b^2
      csq   = ellipse(k,4)^2;        % c^2
      x0    = ellipse(k,5);           % x offset
      y0    = ellipse(k,6);           % y offset
      z0    = ellipse(k,7);           % z offset
      phi   = ellipse(k,8)*pi/180;   % first Euler angle in radians
      theta = ellipse(k,9)*pi/180; % second Euler angle in radians
      psi   = ellipse(k,10)*pi/180;  % third Euler angle in radians
      
      cphi     = cos(phi);
      sphi     = sin(phi);
      ctheta   = cos(theta);
      stheta   = sin(theta);
      cpsi     = cos(psi);
      spsi     = sin(psi);
      
      % Euler rotation matrix
      alpha = [cpsi*cphi-ctheta*sphi*spsi   cpsi*sphi+ctheta*cphi*spsi  spsi*stheta;
               -spsi*cphi-ctheta*sphi*cpsi  -spsi*sphi+ctheta*cphi*cpsi cpsi*stheta;
               stheta*sphi                  -stheta*cphi                ctheta];        
      
      % rotated ellipsoid coordinates
      coordp = alpha*coord;
      
      idx = find((coordp(1,:)-x0).^2./asq + (coordp(2,:)-y0).^2./bsq + (coordp(3,:)-z0).^2./csq <= 1);
      p(idx) = p(idx) + A;
   end

   p = reshape(p,[n n n]);

return;


function out = flatten(in)

   out = reshape(in,[1 numel(in)]);

return;  
