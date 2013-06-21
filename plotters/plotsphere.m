function plotsphere( points, radius )
%plotsphere This function plots samples from noisy spheres in either 2 or 3
%dimensions.

% For 2 dimensions, it plots all the points.
% For 3 dimensions, it plots half the sphere to make for a better
% visualization.

dim=size(points,2);
rad=ceil(max(pdist(points))*2)/2;   % rounds radius upward to nearest 0.5

switch dim
    case 2
       plot(points(:,1),points(:,2),'.');
       
       hold on;
       axis([-rad/2,rad/2,-rad/2,rad/2]);
       axis square;
       
       t = 0:0.1:2*pi;
       xcirc = radius*cos(t);
       ycirc = radius*sin(t);
       plot(xcirc,ycirc,'-r') 
       
       hold off;
       
    case 3
       z = points(:,3); ind = find(z>0);
       filtpts = points(ind,:);     % These two lines extract elements with x>0
       
       plot3(filtpts(:,1),filtpts(:,2),filtpts(:,3),'.');
       
       hold on;
       xlim([-rad/2,rad/2]); ylim([-rad/2,rad/2]); zlim([-rad/2,rad/2]);
       
       [x,y,z] = sphere(100);  % Makes a 101-by-101 point sphere
       x = x(51:end,:);       % Keep top 51 x points
       y = y(51:end,:);       % Keep top 51 y points
       z = z(51:end,:);       % Keep top 51 z points
       
       mesh(radius.*x,radius.*y,radius.*z);  % Plot the surface
    
       axis square;
       hold off;
       
    otherwise
       exception = MException('VerifyOutput:OutOfBounds', ...
       'Cannot plot anything but 2- and 3-D spheres!'); 
       throw(exception);
end
       
end

