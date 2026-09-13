clear all, clc, close all, reset(gpuDevice)

 
load('phSL3D32.mat')
 
M = 2*32; 
Ntheta  = 2*M; %8*M;
R = 1;

theta = make_array(0, 2*pi, Ntheta, true);

% Discretization of axis: x, y, z


N=length(p); 

xM = R/(2*sqrt(2) + 0.01);
x = make_array(-xM, xM, N);
%y = make_array( 0.02, 2*xM, N);
y = x;
z = y;



[rMax] = get_max_r(theta, x, y, z, R);



% variables para el W
Nphi = 4* M-1; % 2* M+1;
Npsi = 6*M+2; %  azimuth integration variable


phi = make_array(   0, 2*pi, Nphi, true);
psi = make_array(0, pi, Npsi, true);
 
 

%%
clc
Nr = 500; 

[p, theta,  phi, psi, R, Nr, x, y, z, rMax] = type_and_device2(p, theta,  phi, psi, R, Nr, x, y, z, rMax);


[Pfz0] = CalcularP_z0(p, theta,  psi, phi, R, Nr, x, y, z, rMax);

 


%%
clc
% nuevas dimensiones
Nalpha = 2*M+1;
Npsi0  = 67;
Nphi0  = 101;


phi0 = make_array(   pi,   2*pi, Nphi0, false);
psi0 = make_array( pi/4, 3*pi/4, Npsi0, false);

[Fz1]  = full_plane_integration_polar_z0( x, y, z, psi, phi,  Pfz0, theta,  Nalpha, R,  1,  psi0, phi0);
Fz1 = gather(Fz1);
 

phi0 = make_array(   0, pi, Nphi0, false);
[Fz2]  = full_plane_integration_polar_z0( x, y, z, psi, phi,  Pfz0, theta,  Nalpha, R, -1,  psi0, phi0);
Fz2 = gather(Fz2);
 

clear Pfz0
%%
clc
%Nphi = 4*M-1; % 2* M+1;
%Npsi = 4*M+2; %  azimuth integration variable


phi = make_array(   0, pi, Nphi, true);
psi = make_array(0, pi, Npsi, true);
 

[Pfy0] = CalcularP_y0(p, theta,  psi, phi, R, Nr, x, y, z, rMax);



phi0 = make_array(   0, 2*pi, Nphi0, false);
psi0 = make_array(   0, pi/4, Npsi0, false);
[Fy1 ] = full_plane_integration_polar_y0( x, y, z, psi, phi,  Pfy0, theta,  Nalpha, R,  1 ,  psi0, phi0);
Fy1=gather(Fy1);


psi0 = make_array(3*pi/4, pi, Npsi0, false);
[Fy2 ] = full_plane_integration_polar_y0( x, y, z, psi, phi,  Pfy0, theta,  Nalpha, R,  1 ,  psi0, phi0);
Fy2=gather(Fy2);
clear Pfy0

%%


phi = make_array(  pi, 2*pi, Nphi, true);
psi = make_array(0, pi, Npsi, true);


[Pfy0] = CalcularP_y0(p, theta,  psi, phi, R, Nr, x, y, z, rMax);



phi0 = make_array(   0, 2*pi, Nphi0, false);
psi0 = make_array(   0, pi/4, Npsi0, false);
[Fy3 ] = full_plane_integration_polar_y0( x, y, z, psi, phi,  Pfy0, theta,  Nalpha, R,  1 ,  psi0, phi0);
Fy3=gather(Fy3);


psi0 = make_array(3*pi/4, pi, Npsi0, false);
[Fy4 ] = full_plane_integration_polar_y0( x, y, z, psi, phi,  Pfy0, theta,  Nalpha, R,  1 ,  psi0, phi0);
Fy4=gather(Fy4);

%%
clc

F =  Fz1 + Fz2  + Fy1+ Fy2 + Fy3 + Fy4;

dx = x(2) - x(1);
dy = y(2) - y(1);
dz = z(2) - z(1);


f_ = -6*del2( F  , dx, dy, dz)/(8*pi^2); 



mostrar_muchos('gray','3',f_, p);



