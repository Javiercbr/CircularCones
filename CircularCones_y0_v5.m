function [go,  varphi, s] = CircularCones_y0_v5(f, theta, Nvarphi, Ns, Nr, Npsi, R, x, y, z, rMax, J, show_legend)


% Calculates the 3D Radon transform of function
% Javier Cebeiro and Marcela Morvidone
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                   VERSION 0
% March 2026
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% INPUT
% f(x,y,z): 3D phantom dim=N x N x N
% Nphi
% Npsi
% Ns
% Nvarphi
% Nr
% symmetry: if 1 values for s<0 are obtained by symmetry property,if 0 values
% for s<0 are obtained by straigthforward calculation


% OUTPUT
% g
% psi : polar angle of the sphere  ( 0<=psi<2pi  )
% phi : azimuthual angle of the sphere  ( 0<=phi<=pi  )
% s
% x ,y , z



% to call the program

% [g, theta, varphi] =


% DISCRETIZATION of the transform variables



dvarphi = pi/(Nvarphi-1);
varphi  = -pi/2:  dvarphi  :pi/2;





dpsi = 2*pi/Npsi ;
psi  = 0:  dpsi   :2*pi-dpsi;

L = floor(Npsi/J);



Ns = 2*floor(Ns/2) + 1 ; % array length must be odd (to have an exact s=0)
ds = 2/(Ns-1);
%s = - 1 + ds : ds : 1 - ds; % full array
s = - 1  : ds : 1 ; % full array

dr = rMax/(Nr-1);
r  = 0 : dr : rMax; % full array



Ntheta = length(theta);

go = zeros(Ntheta, Nvarphi, Ns);

barra = waitbar(0,'Calculation of the cone projection');

for k=1: Ntheta



        waitbar(gather(k/Ntheta), barra,sprintf('Calculation of the cone projection z %d/%d', k, Ntheta));
    
    if show_legend
        fprintf('Calculation of the cone projection z %d%%\n', ...
            round(100*k/Ntheta));
    end
   % waitbar(gather(k/Ntheta));
    g = 0;
    for m=1:Nr                % first integral (phi)



        g0 = 0;

        for j =1:J



            [   varphi_, s_, psi_ ] = ndgrid(  varphi, s, psi(1 + L*(j-1):L*j));








            xi = (-r(m)*s_.*cos(varphi_) + r(m)*sqrt(1-s_.^2).*sin(psi_).*sin(varphi_)   + R).*cos(-theta(k)) + (r(m)*s_.*sin(varphi_) + r(m)*sqrt(1-s_.^2).*sin(psi_).*cos(varphi_)).*sin(-theta(k));
            yi = r(m)*sqrt(1-s_.^2).*cos(psi_);
            zi = -(-r(m)*s_.*cos(varphi_) + r(m)*sqrt(1-s_.^2).*sin(psi_).*sin(varphi_)  + R).*sin(-theta(k)) + (r(m)*s_.*sin(varphi_) + r(m)*sqrt(1-s_.^2).*sin(psi_).*cos(varphi_)).*cos(-theta(k));
            a = interpn(x,y,z,f,xi,yi,zi);  % interpolation
            a( isnan(a))=0;                 % for points outside the limits(returned as NaN)
            %a = sqrt(1-s_.^2).*a;

            g0= g0 + sum(a,3);


        end

        [  varphi_, s_, psi_ ] = ndgrid(  varphi, s, psi(1 + L*J:end));
        xi = (-r(m)*s_.*cos(varphi_) + r(m)*sqrt(1-s_.^2).*sin(psi_).*sin(varphi_)   + R).*cos(-theta(k)) + (r(m)*s_.*sin(varphi_) + r(m)*sqrt(1-s_.^2).*sin(psi_).*cos(varphi_)).*sin(-theta(k));
        yi = r(m)*sqrt(1-s_.^2).*cos(psi_);
        zi = -(-r(m)*s_.*cos(varphi_) + r(m)*sqrt(1-s_.^2).*sin(psi_).*sin(varphi_)  + R).*sin(-theta(k)) + (r(m)*s_.*sin(varphi_) + r(m)*sqrt(1-s_.^2).*sin(psi_).*cos(varphi_)).*cos(-theta(k));
        a = interpn(x,y,z,f,xi,yi,zi);  % interpolation
        a( isnan(a))=0;                 % for points outside the limits(returned as NaN)
        %a = sqrt(1-s_.^2).*a;
        g0= g0 + sum(a,3);




        g = g + r(m)*g0;                          % acumulation


    end








%g=  sin(omega_).*g;

% second integral angle
g=dpsi*dr*squeeze(g);                                      % tu adjust dimension

go(k,:,:) = g;

end





close(barra)

end


