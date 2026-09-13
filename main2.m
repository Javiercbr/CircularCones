clear all, clc, close all, reset(gpuDevice)

 
load('phSL3D32.mat')
 
M = 2*32; 
Ntheta  = 2*M +1;
R = 1;


theta = make_array(0, 2*pi, Ntheta, true);
%theta = theta(54:64 );

% Discretization of axis: x, y, z


N=length(p); 
% 
xM = R/(2*sqrt(2) + 0.01);
x = make_array(-xM, xM, N);
y = make_array( 0.02, 2*xM, N);
z = y;


% zM = R/(sqrt(6) + 0.01);
% z = make_array( 0 ,  zM, N);
% y = z;
% x = make_array( -zM/2 ,  zM/2, N);
%  




[rMax] = get_max_r(theta, x, y, z, R);

 
Nr = 100; % 200; 
J =64; % 32; %16;
Nvarphi =round(8 *M) +1;
Ns =16*M +1;
Npsi = 8*M +1;
[p, theta,  M, R, Nr, x, y, z, rMax, Ns, Nvarphi, Npsi] = type_and_device2(p, theta,  M, R, Nr, x, y, z, rMax, Ns, Nvarphi, Npsi);

 
[gz0, varphi, s] = CircularCones_z0_v5(p,theta, Nvarphi, Ns, Nr, Npsi, R, x, y, z, rMax, J );
save  'gz0.mat' 'gz0' 'theta' 'varphi' 's' 'R' 'x' 'y' 'z'
[gy0, varphi, s] = CircularCones_y0_v5(p,theta, Nvarphi, Ns, Nr, Npsi, R, x, y, z, rMax, J );
save  'gy0.mat' 'gy0' 'theta' 'varphi' 's' 'R' 'x' 'y' 'z'
 
