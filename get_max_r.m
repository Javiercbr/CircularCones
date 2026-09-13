function [rMax] = get_max_r(theta, x, y, z, R)


[ theta_, x_, y_, z_  ]= ndgrid( theta, x, y, z);
radio = sqrt((x_-R*cos(theta_)).^2 + (y_ -R*sin(theta_)).^2 + z_.^2);
rMax = max(max(max(max(radio))));

 
end